"""
PhyloGenie v3 Pipeline — Real Tools
Stages use real APIs: NCBI Entrez, NCBI BLAST, Augustus web, Biopython.

Improvements over v2:
- Retries with backoff around flaky network calls (NCBI, Augustus)
- Bounded concurrency (per-accession fetch/ORF/structure run in parallel;
  a global semaphore caps how many jobs run at once)
- Cooperative job cancellation
- Bio.Align.PairwiseAligner instead of the deprecated Bio.pairwise2
- Structured logging instead of print/silent excepts
"""

import asyncio
import logging
import re
import time
from io import StringIO
from pathlib import Path
from typing import Any, Optional

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import Phylo as BioPhylo

import requests

import config
import db
from pipeline.validators import OPTIONAL_TOOLS

log = logging.getLogger("phylogenie.pipeline")

# ── Entrez setup ──
Entrez.email = config.NCBI_EMAIL or "anonymous@example.com"
if config.NCBI_API_KEY:
    Entrez.api_key = config.NCBI_API_KEY

# ── Job store (in-memory cache, persisted to SQLite) ──
JOB_STORE: dict[str, dict] = {}

Path(config.OUTPUT_DIR).mkdir(exist_ok=True, parents=True)

# Caps how many jobs run their pipeline concurrently, and how many
# per-accession network calls happen at once within a job.
_JOB_SEMAPHORE = asyncio.Semaphore(max(1, config.MAX_CONCURRENT_JOBS))
_FETCH_SEMAPHORE = asyncio.Semaphore(4 if config.NCBI_API_KEY else 2)


class JobCancelled(Exception):
    pass


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def _load_jobs_from_disk():
    for job in db.load_all_jobs():
        # A job still marked "running" or "queued" means the process was
        # killed/restarted mid-flight — there's no in-flight asyncio task to
        # resume, so it would otherwise sit "running" forever. Mark it failed.
        if job.get("status") in ("running", "queued"):
            job["status"] = "error"
            job["error"] = "Interrupted by server restart"
            job["current_step"] = "error"
            job["completed_at"] = job.get("completed_at") or time.time()
            db.save_job(job)
        JOB_STORE[job["id"]] = job


_load_jobs_from_disk()


def _save_job(job: dict):
    try:
        db.save_job(job)
    except Exception:
        log.exception("Failed to persist job %s", job.get("id"))


def job_update(job: dict, step: str, status: str, progress: int, message: str, data: Optional[dict] = None):
    job["steps"][step] = {
        "status": status,
        "message": message,
        "data": data or {},
        "updated_at": time.time(),
    }
    job["progress"] = progress
    job["current_step"] = step
    _save_job(job)


def _check_cancelled(job: dict):
    if job.get("_cancel_requested"):
        raise JobCancelled(f"Job {job['id']} was cancelled by user")


async def _retry(fn, *args, attempts: int = 3, base_delay: float = 1.5, **kwargs):
    """Run a blocking-style async call with exponential backoff retries."""
    last_exc = None
    for attempt in range(1, attempts + 1):
        try:
            return await fn(*args, **kwargs)
        except Exception as e:
            last_exc = e
            if attempt < attempts:
                delay = base_delay * (2 ** (attempt - 1))
                log.warning("Attempt %d/%d failed for %s: %s — retrying in %.1fs",
                            attempt, attempts, getattr(fn, "__name__", fn), e, delay)
                await asyncio.sleep(delay)
    raise last_exc


# ─────────────────────────────────────────────
# Step 1 — Sequence Retrieval (NCBI Entrez)
# ─────────────────────────────────────────────

async def fetch_sequence(accession: str) -> dict:
    """Fetch FASTA + GenBank record from NCBI Entrez."""
    loop = asyncio.get_event_loop()

    def _fetch():
        results = {}

        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            fasta_text = handle.read()
            handle.close()
            records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
            if not records:
                raise RuntimeError(f"No FASTA record returned for '{accession}' — check the accession is correct")
            rec = records[0]
            results["fasta"] = fasta_text
            results["sequence"] = str(rec.seq)
            results["seq_id"] = rec.id
            results["description"] = rec.description
            results["length"] = len(rec.seq)
        except Exception as e:
            raise RuntimeError(f"FASTA fetch failed for {accession}: {e}") from e

        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            gb_text = handle.read()
            handle.close()
            gb_records = list(SeqIO.parse(StringIO(gb_text), "genbank"))
            if gb_records:
                gb = gb_records[0]
                features = []
                for feat in gb.features[:20]:
                    q = {k: (v[0] if isinstance(v, list) and v else v) for k, v in feat.qualifiers.items()}
                    features.append({"type": feat.type, "location": str(feat.location), "qualifiers": q})
                results["genbank"] = {
                    "accession": accession,
                    "name": gb.name,
                    "organism": gb.annotations.get("organism", "Unknown"),
                    "taxonomy": gb.annotations.get("taxonomy", []),
                    "definition": gb.description,
                    "keywords": gb.annotations.get("keywords", []),
                    "references": [str(r.title) for r in gb.annotations.get("references", [])[:3]],
                    "features": features,
                    "molecule_type": gb.annotations.get("molecule_type", ""),
                    "topology": gb.annotations.get("topology", ""),
                    "date": gb.annotations.get("date", ""),
                    "source": gb.annotations.get("source", ""),
                }
        except Exception as e:
            log.warning("GenBank fetch failed for %s: %s", accession, e)
            results["genbank"] = {"error": str(e)}

        results["accession"] = accession
        return results

    async with _FETCH_SEMAPHORE:
        result = await _retry(lambda: loop.run_in_executor(None, _fetch), attempts=3, base_delay=2)
        await asyncio.sleep(0.35 if config.NCBI_API_KEY else 0.5)  # respect NCBI rate limits
        return result


# ─────────────────────────────────────────────
# Step 2 — ORF Finder
# ─────────────────────────────────────────────

def find_orfs(sequence: str, min_len: int = 100) -> list[dict]:
    """Find all ORFs in 6 reading frames."""
    seq = Seq(sequence.upper())
    orfs = []
    frames = [
        (seq, "+"),
        (seq[1:], "+"),
        (seq[2:], "+"),
        (seq.reverse_complement(), "-"),
        (seq.reverse_complement()[1:], "-"),
        (seq.reverse_complement()[2:], "-"),
    ]
    for frame_idx, (s, strand) in enumerate(frames):
        offset = frame_idx % 3
        protein = str(s.translate())
        start = 0
        while True:
            m_pos = protein.find("M", start)
            if m_pos == -1:
                break
            stop_pos = protein.find("*", m_pos)
            if stop_pos == -1:
                stop_pos = len(protein)
            orf_protein = protein[m_pos:stop_pos]
            orf_len = len(orf_protein)
            if orf_len >= min_len // 3:
                if strand == "+":
                    nuc_start = m_pos * 3 + offset
                    nuc_end = stop_pos * 3 + offset
                else:
                    nuc_start = len(sequence) - (stop_pos * 3 + offset)
                    nuc_end = len(sequence) - (m_pos * 3 + offset)
                orfs.append({
                    "frame": f"{strand}{frame_idx % 3 + 1}",
                    "strand": strand,
                    "start": nuc_start,
                    "end": nuc_end,
                    "length_nt": (stop_pos - m_pos) * 3,
                    "length_aa": orf_len,
                    "protein": orf_protein[:200],
                    "has_stop": stop_pos < len(protein),
                    "gc_content": round(
                        (str(s[m_pos * 3:stop_pos * 3]).count("G") +
                         str(s[m_pos * 3:stop_pos * 3]).count("C")) /
                        max((stop_pos - m_pos) * 3, 1) * 100, 1)
                })
            start = m_pos + 1

    orfs.sort(key=lambda x: x["length_aa"], reverse=True)
    return orfs[:20]


async def run_orf_finder(sequence: str) -> dict:
    loop = asyncio.get_event_loop()
    orfs = await loop.run_in_executor(None, find_orfs, sequence)
    return {
        "total_orfs": len(orfs),
        "orfs": orfs,
        "longest_orf": orfs[0] if orfs else None,
        "sequence_length": len(sequence),
    }


# ─────────────────────────────────────────────
# Step 3 — Augustus Gene Prediction (Web API, with local fallback)
# ─────────────────────────────────────────────

async def run_augustus(sequence: str, organism: str = "human") -> dict:
    loop = asyncio.get_event_loop()

    def _submit_and_poll():
        species_map = {
            "human": "human", "mouse": "mouse", "arabidopsis": "arabidopsis",
            "drosophila": "fly", "zebrafish": "zebrafish", "rice": "rice",
        }
        species = species_map.get(organism.lower(), "human")
        seq_to_submit = sequence[:50000]
        fasta = f">query_sequence\n{seq_to_submit}\n"

        try:
            resp = requests.post(
                config.AUGUSTUS_URL,
                data={"species": species, "sequence": fasta, "format": "plain", "report": "gene"},
                timeout=30,
            )
            resp.raise_for_status()

            job_id_match = re.search(r"jobid=(\w+)", resp.url + resp.text)
            if job_id_match:
                aug_job_id = job_id_match.group(1)
                for _ in range(config.AUGUSTUS_POLL_ATTEMPTS):
                    time.sleep(config.AUGUSTUS_POLL_INTERVAL)
                    result_url = f"https://bioinf.uni-greifswald.de/augustus/results/{aug_job_id}"
                    r = requests.get(result_url, timeout=15)
                    if r.status_code == 200 and "# start gene" in r.text:
                        return _parse_augustus_output(r.text, organism)
            return _parse_augustus_output(resp.text, organism)
        except Exception as e:
            log.warning("Augustus web submission failed: %s", e)
            return {"error": str(e), "genes": [], "genes_predicted": 0}

    result = await loop.run_in_executor(None, _submit_and_poll)

    if result.get("error") or not result.get("genes"):
        orfs = find_orfs(sequence, min_len=300)
        genes = []
        for i, orf in enumerate(orfs[:6]):
            genes.append({
                "id": f"gene{i + 1}",
                "start": orf["start"], "end": orf["end"], "strand": orf["strand"],
                "transcript": f"transcript{i + 1}",
                "cds_start": orf["start"], "cds_end": orf["end"],
                "length_aa": orf["length_aa"], "protein": orf["protein"],
                "source": "ORF-based fallback",
            })
        return {
            "genes_predicted": len(genes),
            "genes": genes,
            "organism_model": organism,
            "source": "ORF fallback (Augustus web unavailable)",
            "raw_output": "",
        }

    return result


def _parse_augustus_output(text: str, organism: str) -> dict:
    genes = []
    current_gene = None

    for line in text.splitlines():
        line = line.strip()
        if line.startswith("# start gene"):
            gene_id = line.split()[-1]
            current_gene = {"id": gene_id, "exons": [], "source": "Augustus"}
        elif line.startswith("# end gene") and current_gene:
            genes.append(current_gene)
            current_gene = None
        elif line and not line.startswith("#") and current_gene is not None:
            parts = line.split("\t")
            if len(parts) >= 8:
                feature = parts[2]
                try:
                    start, end, strand = int(parts[3]), int(parts[4]), parts[6]
                    current_gene["strand"] = strand
                    if feature == "gene":
                        current_gene["start"] = start
                        current_gene["end"] = end
                        current_gene["length_nt"] = end - start
                    elif feature == "CDS":
                        current_gene["exons"].append({"start": start, "end": end, "length": end - start})
                except (ValueError, IndexError):
                    pass
        elif line.startswith("# ") and current_gene and "protein" in text:
            prot_match = re.search(r"\[(.+?)\]", line)
            if prot_match:
                current_gene["protein"] = prot_match.group(1)[:200]

    for g in genes:
        g.setdefault("length_aa", g.get("length_nt", 0) // 3)
        g.setdefault("protein", "")

    return {
        "genes_predicted": len(genes),
        "genes": genes,
        "organism_model": organism,
        "source": "Augustus Web Server",
        "raw_output": text[:3000],
    }


# ─────────────────────────────────────────────
# Step 4 — BLASTN (nucleotide vs nt)
# ─────────────────────────────────────────────

async def run_blastn(sequence: str, accession: str) -> dict:
    loop = asyncio.get_event_loop()

    def _blast():
        try:
            query_seq = sequence[:2000]
            result_handle = NCBIWWW.qblast(
                "blastn", "nt", query_seq,
                hitlist_size=config.BLAST_HITLIST_SIZE,
                expect=config.BLAST_EXPECT,
                format_type="XML",
            )
            raw = result_handle.read()
            result_handle.close()
            blast_records = list(NCBIXML.parse(StringIO(raw)))
            hits = []
            for record in blast_records:
                for alignment in record.alignments[:config.BLAST_HITLIST_SIZE]:
                    hsp = alignment.hsps[0]
                    hits.append({
                        "title": alignment.title[:120],
                        "accession": alignment.accession,
                        "length": alignment.length,
                        "score": hsp.score,
                        "bits": round(hsp.bits, 1),
                        "e_value": hsp.expect,
                        "identities": hsp.identities,
                        "gaps": hsp.gaps,
                        "align_length": hsp.align_length,
                        "identity_pct": round(hsp.identities / hsp.align_length * 100, 1) if hsp.align_length else 0,
                        "query_cover": round(hsp.align_length / len(query_seq) * 100, 1),
                        "query_seq": hsp.query[:80],
                        "subject_seq": hsp.sbjct[:80],
                        "midline": hsp.match[:80],
                    })
            hits.sort(key=lambda x: x["bits"], reverse=True)
            return {"hits": hits, "query_length": len(query_seq), "query_accession": accession, "program": "blastn", "database": "nt"}
        except Exception as e:
            log.warning("BLASTN failed for %s: %s", accession, e)
            return {"error": str(e), "hits": [], "query_length": len(sequence), "program": "blastn"}

    return await loop.run_in_executor(None, _blast)


# ─────────────────────────────────────────────
# Step 5 — BLASTP (protein vs nr)
# ─────────────────────────────────────────────

async def run_blastp(protein_seq: str, gene_id: str = "query") -> dict:
    loop = asyncio.get_event_loop()

    def _blast():
        if not protein_seq or len(protein_seq) < 10:
            return {"error": "No protein sequence", "hits": []}
        try:
            query = protein_seq[:500]
            result_handle = NCBIWWW.qblast("blastp", "nr", query, hitlist_size=15, expect=0.001, format_type="XML")
            raw = result_handle.read()
            result_handle.close()
            blast_records = list(NCBIXML.parse(StringIO(raw)))
            hits = []
            for record in blast_records:
                for alignment in record.alignments[:15]:
                    hsp = alignment.hsps[0]
                    hits.append({
                        "title": alignment.title[:120],
                        "accession": alignment.accession,
                        "length": alignment.length,
                        "score": hsp.score,
                        "bits": round(hsp.bits, 1),
                        "e_value": hsp.expect,
                        "identity_pct": round(hsp.identities / hsp.align_length * 100, 1) if hsp.align_length else 0,
                        "query_cover": round(hsp.align_length / len(query) * 100, 1),
                        "positives_pct": round(hsp.positives / hsp.align_length * 100, 1) if hsp.align_length else 0,
                    })
            hits.sort(key=lambda x: x["bits"], reverse=True)
            return {"hits": hits, "query_length": len(query), "gene_id": gene_id, "program": "blastp", "database": "nr"}
        except Exception as e:
            log.warning("BLASTP failed for %s: %s", gene_id, e)
            return {"error": str(e), "hits": [], "program": "blastp"}

    return await loop.run_in_executor(None, _blast)


# ─────────────────────────────────────────────
# Step 6 — MSA (only if multiple sequences)
# ─────────────────────────────────────────────

async def run_msa(sequences: list[dict]) -> dict:
    """Pairwise-align every sequence to the first, using Bio.Align.PairwiseAligner
    (the modern, non-deprecated replacement for Bio.pairwise2)."""
    if len(sequences) < 2:
        return {"skipped": True, "reason": "MSA requires 2+ sequences"}

    loop = asyncio.get_event_loop()

    def _msa():
        try:
            aligner = PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5

            seqs = [s["sequence"][:1000] for s in sequences]
            ids = [s["accession"] for s in sequences]

            ref = seqs[0]
            aligned_seqs = [ref]
            scores = [100.0]

            for seq in seqs[1:]:
                alignments = aligner.align(ref, seq)
                best = alignments[0]
                a_ref, a_seq = str(best[0]), str(best[1])
                aligned_seqs.append(a_seq)
                identity = sum(a == b for a, b in zip(a_ref, a_seq) if a != "-" and b != "-")
                total = max(len(a_ref), 1)
                scores.append(round(identity / total * 100, 1))

            max_len = max(len(s) for s in aligned_seqs)
            aligned_seqs = [s.ljust(max_len, "-") for s in aligned_seqs]

            conservation = []
            for col_i in range(min(max_len, 200)):
                col = [s[col_i] for s in aligned_seqs if col_i < len(s)]
                if all(c == col[0] and c != "-" for c in col):
                    conservation.append("*")
                elif len(set(c for c in col if c != "-")) <= 2:
                    conservation.append(":")
                else:
                    conservation.append(" ")

            return {
                "num_sequences": len(sequences),
                "alignment_length": max_len,
                "method": "Biopython PairwiseAligner (global)",
                "aligned": [
                    {"id": ids[i], "sequence": aligned_seqs[i][:500], "identity_to_ref": scores[i]}
                    for i in range(len(ids))
                ],
                "conservation": "".join(conservation[:200]),
            }
        except Exception as e:
            log.exception("MSA failed")
            return {"error": str(e), "num_sequences": len(sequences)}

    return await loop.run_in_executor(None, _msa)


# ─────────────────────────────────────────────
# Step 7 — Phylogenetic Tree (Biopython NJ)
# ─────────────────────────────────────────────

async def run_phylo_tree(sequences: list[dict]) -> dict:
    if len(sequences) < 2:
        return {"skipped": True, "reason": "Phylogeny requires 2+ sequences"}

    loop = asyncio.get_event_loop()

    def _phylo():
        try:
            ids = [s["accession"] for s in sequences]
            seqs_raw = [s["sequence"][:600] for s in sequences]
            min_len = min(len(s) for s in seqs_raw)
            seqs_raw = [s[:min_len] for s in seqs_raw]

            aln = MultipleSeqAlignment([
                SeqRecord(Seq(s), id=ids[i], description="") for i, s in enumerate(seqs_raw)
            ])

            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(aln)

            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)

            def tree_to_dict(clade):
                node = {"name": clade.name or "", "branch_length": round(float(clade.branch_length or 0), 5)}
                if clade.clades:
                    node["children"] = [tree_to_dict(c) for c in clade.clades]
                return node

            tree_dict = tree_to_dict(tree.root)

            buf = StringIO()
            BioPhylo.write(tree, buf, "newick")
            newick = buf.getvalue().strip()

            return {
                "method": "Neighbor-Joining (Biopython)",
                "distance_model": "identity",
                "taxa": len(sequences),
                "tree": tree_dict,
                "newick": newick,
            }
        except Exception as e:
            log.exception("Phylo tree construction failed")
            return {"error": str(e), "taxa": len(sequences)}

    return await loop.run_in_executor(None, _phylo)


# ─────────────────────────────────────────────
# Step 8 — Protein Secondary Structure (Biopython ProtParam)
# ─────────────────────────────────────────────

async def run_protein_analysis(protein_seq: str, gene_id: str = "protein") -> dict:
    loop = asyncio.get_event_loop()

    def _analyse():
        clean = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", protein_seq.upper())
        if len(clean) < 5:
            return {"error": "Protein too short or invalid", "gene_id": gene_id}

        try:
            analysis = ProteinAnalysis(clean)

            aa_comp = analysis.get_amino_acids_percent()
            mw = analysis.molecular_weight()
            ip = analysis.isoelectric_point()
            instability = analysis.instability_index()
            gravy = analysis.gravy()
            aromaticity = analysis.aromaticity()

            try:
                helix, turn, sheet = analysis.secondary_structure_fraction()
            except Exception:
                helix, turn, sheet = 0.0, 0.0, 0.0
            coil = max(0.0, 1.0 - helix - turn - sheet)

            try:
                ec_cys, ec_no_cys = analysis.molar_extinction_coefficient()
            except Exception:
                ec_cys, ec_no_cys = 0, 0

            try:
                flex = analysis.flexibility()
                avg_flex = round(sum(flex) / len(flex), 4) if flex else 0
            except Exception:
                avg_flex = 0

            seq_len = len(clean)
            n_helix, n_turn, n_sheet = int(helix * seq_len), int(turn * seq_len), int(sheet * seq_len)
            n_coil = seq_len - n_helix - n_turn - n_sheet
            ss_string = "H" * n_helix + "T" * n_turn + "E" * n_sheet + "C" * max(0, n_coil)

            return {
                "gene_id": gene_id,
                "protein_length": len(clean),
                "molecular_weight": round(mw, 2),
                "isoelectric_point": round(ip, 2),
                "instability_index": round(instability, 2),
                "is_stable": instability < 40,
                "gravy": round(gravy, 4),
                "aromaticity": round(aromaticity, 4),
                "avg_flexibility": avg_flex,
                "extinction_coeff_with_cys": ec_cys,
                "extinction_coeff_no_cys": ec_no_cys,
                "secondary_structure": {
                    "helix": round(helix * 100, 1), "turn": round(turn * 100, 1),
                    "sheet": round(sheet * 100, 1), "coil": round(coil * 100, 1),
                },
                "ss_string": ss_string[:200],
                "amino_acid_composition": {
                    k: round(v * 100, 2) for k, v in sorted(aa_comp.items(), key=lambda x: -x[1])[:10]
                },
                "protein_sequence": clean[:300],
            }
        except Exception as e:
            log.exception("Protein analysis failed for %s", gene_id)
            return {"error": str(e), "gene_id": gene_id, "protein_sequence": clean[:100]}

    return await loop.run_in_executor(None, _analyse)


# ─────────────────────────────────────────────
# Full Pipeline Orchestrator
# ─────────────────────────────────────────────

async def run_full_pipeline(job_id: str, accessions: list[str], options: dict):
    job = JOB_STORE[job_id]

    job_update(job, "retrieval", "queued", 0, "Waiting for a free execution slot...")
    async with _JOB_SEMAPHORE:
        await _run_pipeline_body(job, accessions, options)


async def _run_pipeline_body(job: dict, accessions: list[str], options: dict):
    job_id = job["id"]
    multiple = len(accessions) > 1
    tools_opt = options.get("tools") or {}

    def wants(tool: str) -> bool:
        return bool(tools_opt.get(tool, True))

    try:
        _check_cancelled(job)

        # ── Step 1: Retrieve sequences (parallel, bounded) ──
        job_update(job, "retrieval", "running", 5, f"Fetching {len(accessions)} sequence(s) from NCBI...")
        try:
            sequences = await asyncio.gather(*(fetch_sequence(acc) for acc in accessions))
        except Exception as e:
            job_update(job, "retrieval", "error", 5, f"Sequence retrieval failed: {e}")
            raise
        sequences = list(sequences)

        job_update(job, "retrieval", "complete", 18,
                   f"Retrieved {len(sequences)} sequence(s) from NCBI GenBank",
                   {"sequences": [{k: v for k, v in s.items() if k != "fasta"} for s in sequences]})
        _check_cancelled(job)

        # ── Step 2: ORF Finder (parallel) ──
        job_update(job, "orf", "running", 20, "Finding open reading frames (6 frames)...")
        orf_data_list = await asyncio.gather(*(run_orf_finder(s["sequence"]) for s in sequences))
        orf_results = [{"accession": sequences[i]["accession"], "result": orf_data_list[i]} for i in range(len(sequences))]
        job_update(job, "orf", "complete", 32,
                   f"Found {sum(r['result']['total_orfs'] for r in orf_results)} ORFs total",
                   {"orf_results": orf_results})
        _check_cancelled(job)

        # ── Step 3: Augustus gene prediction (parallel, optional) ──
        organism_model = options.get("organism_model", "human")
        if wants("augustus"):
            job_update(job, "augustus", "running", 34, "Submitting to Augustus gene predictor (web)...")
            aug_data_list = await asyncio.gather(*(run_augustus(s["sequence"], organism_model) for s in sequences))
            augustus_results = [{"accession": sequences[i]["accession"], "result": aug_data_list[i]} for i in range(len(sequences))]
            total_genes = sum(r["result"].get("genes_predicted", 0) for r in augustus_results)
            job_update(job, "augustus", "complete", 46, f"Augustus predicted {total_genes} gene(s)",
                       {"augustus_results": augustus_results})
        else:
            augustus_results = [
                {"accession": s["accession"], "result": {
                    "genes_predicted": 0, "genes": [], "organism_model": organism_model,
                    "source": "Skipped by user selection", "raw_output": "",
                }} for s in sequences
            ]
            job_update(job, "augustus", "skipped", 46, "Skipped — disabled for this job", {"skipped": True})
        _check_cancelled(job)

        # ── Step 4: BLASTN (sequential — NCBI rate-limits web BLAST heavily; optional) ──
        if wants("blastn"):
            job_update(job, "blastn", "running", 48, "Running BLASTN against NCBI nt database...")
            blastn_results = []
            for seq in sequences:
                _check_cancelled(job)
                bn = await run_blastn(seq["sequence"], seq["accession"])
                blastn_results.append({"accession": seq["accession"], "result": bn})
                await asyncio.sleep(1)
            job_update(job, "blastn", "complete", 60,
                       f"BLASTN complete — {sum(len(r['result'].get('hits', [])) for r in blastn_results)} hits",
                       {"blastn_results": blastn_results})
        else:
            job_update(job, "blastn", "skipped", 60, "Skipped — disabled for this job", {"skipped": True})

        # ── Step 5: BLASTP (use longest ORF/Augustus protein; optional) ──
        if wants("blastp"):
            job_update(job, "blastp", "running", 62, "Running BLASTP on predicted proteins...")
            blastp_results = []
            for i, seq in enumerate(sequences):
                _check_cancelled(job)
                protein = ""
                orfs = orf_results[i]["result"].get("orfs", [])
                if orfs:
                    protein = orfs[0]["protein"]
                aug_genes = augustus_results[i]["result"].get("genes", [])
                if aug_genes and aug_genes[0].get("protein"):
                    protein = aug_genes[0]["protein"]

                if protein:
                    bp = await run_blastp(protein, seq["accession"])
                    blastp_results.append({"accession": seq["accession"], "result": bp})
                    await asyncio.sleep(1)
            job_update(job, "blastp", "complete", 74,
                       f"BLASTP complete — {sum(len(r['result'].get('hits', [])) for r in blastp_results)} hits",
                       {"blastp_results": blastp_results})
        else:
            job_update(job, "blastp", "skipped", 74, "Skipped — disabled for this job", {"skipped": True})
        _check_cancelled(job)

        # ── Step 6: MSA (optional, requires 2+ sequences) ──
        if multiple and wants("msa"):
            job_update(job, "msa", "running", 76, "Running multiple sequence alignment...")
            msa_result = await run_msa(sequences)
            job_update(job, "msa", "complete", 82, "MSA complete", {"msa": msa_result})
        elif not multiple:
            job_update(job, "msa", "skipped", 82, "Skipped — MSA requires 2+ sequences",
                       {"skipped": True, "reason": "Single sequence"})
        else:
            job_update(job, "msa", "skipped", 82, "Skipped — disabled for this job", {"skipped": True})

        # ── Step 7: Phylogenetic tree (optional, requires 2+ sequences) ──
        if multiple and wants("phylo"):
            job_update(job, "phylo", "running", 83, "Building phylogenetic tree (Neighbor-Joining)...")
            phylo_result = await run_phylo_tree(sequences)
            job_update(job, "phylo", "complete", 88, "Phylogenetic tree built", {"tree": phylo_result})
        elif not multiple:
            job_update(job, "phylo", "skipped", 88, "Skipped — requires 2+ sequences",
                       {"skipped": True, "reason": "Single sequence"})
        else:
            job_update(job, "phylo", "skipped", 88, "Skipped — disabled for this job", {"skipped": True})
        _check_cancelled(job)

        # ── Step 8: Protein structure analysis (parallel, optional) ──
        if wants("structure"):
            job_update(job, "structure", "running", 89, "Analysing protein physicochemical properties...")
            struct_inputs = []
            for i in range(len(sequences)):
                orfs = orf_results[i]["result"].get("orfs", [])
                if orfs:
                    struct_inputs.append((orfs[0]["protein"], sequences[i]["accession"]))
            structure_results = list(await asyncio.gather(
                *(run_protein_analysis(protein, gene_id) for protein, gene_id in struct_inputs)
            ))
            job_update(job, "structure", "complete", 97,
                       f"Protein analysis complete for {len(structure_results)} sequence(s)",
                       {"structures": structure_results})
        else:
            job_update(job, "structure", "skipped", 97, "Skipped — disabled for this job", {"skipped": True})

        job["status"] = "complete"
        job["progress"] = 100
        job["current_step"] = "done"
        job["completed_at"] = time.time()
        _save_job(job)
        log.info("Job %s completed successfully", job_id)

    except JobCancelled:
        job["status"] = "cancelled"
        job["current_step"] = "cancelled"
        job["completed_at"] = time.time()
        _save_job(job)
        log.info("Job %s was cancelled", job_id)

    except Exception as e:
        job["status"] = "error"
        job["error"] = str(e)
        job["current_step"] = "error"
        job["completed_at"] = time.time()
        _save_job(job)
        log.exception("Job %s failed", job_id)
