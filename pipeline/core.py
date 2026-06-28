"""
PhyloGenie Pipeline — Real Tools
All stages use real APIs: NCBI Entrez, NCBI BLAST, Augustus web, Biopython
"""

import asyncio
import time
import json
import os
import re
import math
import hashlib
from pathlib import Path
from typing import Optional, List, Dict, Any
from io import StringIO

# Biopython
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import requests
import config

# ── Entrez setup ──
Entrez.email = config.NCBI_EMAIL
if config.NCBI_API_KEY:
    Entrez.api_key = config.NCBI_API_KEY

# ── Job store ──
JOB_STORE: Dict[str, Dict] = {}

Path(config.OUTPUT_DIR).mkdir(exist_ok=True)


# ─────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────

def job_update(job: Dict, step: str, status: str, progress: int, message: str, data: Dict = None):
    job["steps"][step] = {
        "status": status,
        "message": message,
        "data": data or {},
        "updated_at": time.time()
    }
    job["progress"] = progress
    job["current_step"] = step
    _save_job(job)


def _save_job(job: Dict):
    path = Path(config.OUTPUT_DIR) / f"job_{job['id']}.json"
    try:
        with open(path, "w") as f:
            json.dump(job, f, indent=2, default=str)
    except Exception:
        pass


def _load_jobs_from_disk():
    """Load persisted jobs on startup."""
    for f in Path(config.OUTPUT_DIR).glob("job_*.json"):
        try:
            with open(f) as fh:
                job = json.load(fh)
                JOB_STORE[job["id"]] = job
        except Exception:
            pass


_load_jobs_from_disk()


# ─────────────────────────────────────────────
# Step 1 — Sequence Retrieval (NCBI Entrez)
# ─────────────────────────────────────────────

async def fetch_sequence(accession: str) -> Dict:
    """Fetch FASTA + GenBank record from NCBI Entrez."""
    loop = asyncio.get_event_loop()

    def _fetch():
        results = {}

        # Fetch FASTA
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="fasta", retmode="text")
            fasta_text = handle.read()
            handle.close()
            records = list(SeqIO.parse(StringIO(fasta_text), "fasta"))
            if records:
                rec = records[0]
                results["fasta"] = fasta_text
                results["sequence"] = str(rec.seq)
                results["seq_id"] = rec.id
                results["description"] = rec.description
                results["length"] = len(rec.seq)
        except Exception as e:
            raise RuntimeError(f"FASTA fetch failed for {accession}: {e}")

        # Fetch GenBank
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gb", retmode="text")
            gb_text = handle.read()
            handle.close()
            gb_records = list(SeqIO.parse(StringIO(gb_text), "genbank"))
            if gb_records:
                gb = gb_records[0]
                features = []
                for feat in gb.features[:20]:
                    q = {k: (v[0] if isinstance(v, list) and v else v)
                         for k, v in feat.qualifiers.items()}
                    features.append({
                        "type": feat.type,
                        "location": str(feat.location),
                        "qualifiers": q
                    })
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
            results["genbank"] = {"error": str(e)}

        results["accession"] = accession
        return results

    return await loop.run_in_executor(None, _fetch)


# ─────────────────────────────────────────────
# Step 2 — ORF Finder
# ─────────────────────────────────────────────

def find_orfs(sequence: str, min_len: int = 100) -> List[Dict]:
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
                    "protein": orf_protein[:200],  # cap for display
                    "has_stop": stop_pos < len(protein),
                    "gc_content": round(
                        (str(s[m_pos*3:stop_pos*3]).count("G") +
                         str(s[m_pos*3:stop_pos*3]).count("C")) /
                        max((stop_pos - m_pos) * 3, 1) * 100, 1)
                })
            start = m_pos + 1

    orfs.sort(key=lambda x: x["length_aa"], reverse=True)
    return orfs[:20]


async def run_orf_finder(sequence: str) -> Dict:
    loop = asyncio.get_event_loop()
    orfs = await loop.run_in_executor(None, find_orfs, sequence)
    return {
        "total_orfs": len(orfs),
        "orfs": orfs,
        "longest_orf": orfs[0] if orfs else None,
        "sequence_length": len(sequence)
    }


# ─────────────────────────────────────────────
# Step 3 — Augustus Gene Prediction (Web API)
# ─────────────────────────────────────────────

async def run_augustus(sequence: str, organism: str = "human") -> Dict:
    """Submit to Augustus web server and poll for results."""
    loop = asyncio.get_event_loop()

    def _submit_and_poll():
        # Map organism names to Augustus species codes
        species_map = {
            "human": "human", "mouse": "mouse", "arabidopsis": "arabidopsis",
            "drosophila": "fly", "zebrafish": "zebrafish", "rice": "rice",
        }
        species = species_map.get(organism.lower(), "human")

        # Truncate sequence to avoid timeout (Augustus web limits ~100kb)
        seq_to_submit = sequence[:50000] if len(sequence) > 50000 else sequence
        fasta = f">query_sequence\n{seq_to_submit}\n"

        try:
            resp = requests.post(
                "https://bioinf.uni-greifswald.de/augustus/submission.php",
                data={
                    "species": species,
                    "sequence": fasta,
                    "format": "plain",
                    "report": "gene",
                },
                timeout=30
            )
            resp.raise_for_status()

            # Parse job ID from response
            job_id_match = re.search(r'jobid=(\w+)', resp.url + resp.text)
            if job_id_match:
                job_id = job_id_match.group(1)
                # Poll for results
                for _ in range(30):
                    time.sleep(10)
                    result_url = f"https://bioinf.uni-greifswald.de/augustus/results/{job_id}"
                    r = requests.get(result_url, timeout=15)
                    if r.status_code == 200 and "# start gene" in r.text:
                        return _parse_augustus_output(r.text, organism)
            # If can't poll, parse inline result
            return _parse_augustus_output(resp.text, organism)
        except Exception as e:
            return {"error": str(e), "genes": [], "genes_predicted": 0}

    result = await loop.run_in_executor(None, _submit_and_poll)

    # Fallback: if Augustus web failed, do local ORF-based gene prediction
    if result.get("error") or not result.get("genes"):
        orfs = find_orfs(sequence, min_len=300)
        genes = []
        for i, orf in enumerate(orfs[:6]):
            genes.append({
                "id": f"gene{i+1}",
                "start": orf["start"],
                "end": orf["end"],
                "strand": orf["strand"],
                "transcript": f"transcript{i+1}",
                "cds_start": orf["start"],
                "cds_end": orf["end"],
                "length_aa": orf["length_aa"],
                "protein": orf["protein"],
                "source": "ORF-based fallback"
            })
        return {
            "genes_predicted": len(genes),
            "genes": genes,
            "organism_model": organism,
            "source": "ORF fallback (Augustus web unavailable)",
            "raw_output": ""
        }

    return result


def _parse_augustus_output(text: str, organism: str) -> Dict:
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
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    current_gene["strand"] = strand
                    if feature == "gene":
                        current_gene["start"] = start
                        current_gene["end"] = end
                        current_gene["length_nt"] = end - start
                    elif feature == "CDS":
                        current_gene["exons"].append({"start": start, "end": end, "length": end - start})
                except Exception:
                    pass
        elif "# protein sequence" in line and current_gene:
            pass
        elif line.startswith("# ") and current_gene and "protein" in text:
            prot_match = re.search(r'\[(.+?)\]', line)
            if prot_match:
                current_gene["protein"] = prot_match.group(1)[:200]

    for g in genes:
        if "length_aa" not in g:
            g["length_aa"] = g.get("length_nt", 0) // 3
        if "protein" not in g:
            g["protein"] = ""

    return {
        "genes_predicted": len(genes),
        "genes": genes,
        "organism_model": organism,
        "source": "Augustus Web Server",
        "raw_output": text[:3000]
    }


# ─────────────────────────────────────────────
# Step 4 — BLASTN (nucleotide vs nt)
# ─────────────────────────────────────────────

async def run_blastn(sequence: str, accession: str) -> Dict:
    loop = asyncio.get_event_loop()

    def _blast():
        try:
            # Use first 2000bp for speed; NCBI web BLAST has limits
            query_seq = sequence[:2000]
            result_handle = NCBIWWW.qblast(
                "blastn", "nt", query_seq,
                hitlist_size=config.BLAST_HITLIST_SIZE,
                expect=config.BLAST_EXPECT,
                format_type="XML"
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
            return {"error": str(e), "hits": [], "query_length": len(sequence), "program": "blastn"}

    return await loop.run_in_executor(None, _blast)


# ─────────────────────────────────────────────
# Step 5 — BLASTP (protein vs nr)
# ─────────────────────────────────────────────

async def run_blastp(protein_seq: str, gene_id: str = "query") -> Dict:
    loop = asyncio.get_event_loop()

    def _blast():
        if not protein_seq or len(protein_seq) < 10:
            return {"error": "No protein sequence", "hits": []}
        try:
            query = protein_seq[:500]  # Cap for speed
            result_handle = NCBIWWW.qblast(
                "blastp", "nr", query,
                hitlist_size=15,
                expect=0.001,
                format_type="XML"
            )
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
            return {"error": str(e), "hits": [], "program": "blastp"}

    return await loop.run_in_executor(None, _blast)


# ─────────────────────────────────────────────
# Step 6 — MSA (only if multiple sequences)
# ─────────────────────────────────────────────

async def run_msa(sequences: List[Dict]) -> Dict:
    """Run MSA using NCBI BLAST pairwise + Biopython ClustalW-style alignment."""
    if len(sequences) < 2:
        return {"skipped": True, "reason": "MSA requires 2+ sequences"}

    loop = asyncio.get_event_loop()

    def _msa():
        try:
            from Bio import pairwise2
            from Bio.pairwise2 import format_alignment

            seqs = [s["sequence"][:1000] for s in sequences]  # cap for speed
            ids = [s["accession"] for s in sequences]

            # Build aligned records using global alignment
            records = []
            # Use the first sequence as reference, align all others to it
            ref = seqs[0]
            aligned_seqs = [ref]
            scores = [100.0]

            for i, seq in enumerate(seqs[1:], 1):
                alns = pairwise2.align.globalms(ref, seq, 2, -1, -2, -0.5)
                if alns:
                    aligned_seqs.append(alns[0].seqB)
                    score = alns[0].score
                    identity = sum(a == b for a, b in zip(alns[0].seqA, alns[0].seqB)) / len(alns[0].seqA) * 100
                    scores.append(round(identity, 1))
                else:
                    aligned_seqs.append(seq.ljust(len(ref), "-"))
                    scores.append(0.0)

            max_len = max(len(s) for s in aligned_seqs)
            aligned_seqs = [s.ljust(max_len, "-") for s in aligned_seqs]

            # Conservation string
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
                "method": "Biopython global pairwise alignment",
                "aligned": [
                    {"id": ids[i], "sequence": aligned_seqs[i][:500], "identity_to_ref": scores[i]}
                    for i in range(len(ids))
                ],
                "conservation": "".join(conservation[:200])
            }
        except Exception as e:
            return {"error": str(e), "num_sequences": len(sequences)}

    return await loop.run_in_executor(None, _msa)


# ─────────────────────────────────────────────
# Step 7 — Phylogenetic Tree (Biopython NJ)
# ─────────────────────────────────────────────

async def run_phylo_tree(sequences: List[Dict]) -> Dict:
    if len(sequences) < 2:
        return {"skipped": True, "reason": "Phylogeny requires 2+ sequences"}

    loop = asyncio.get_event_loop()

    def _phylo():
        try:
            from Bio.Align import MultipleSeqAlignment
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
            from io import StringIO
            from Bio import Phylo as BioPhylo

            ids = [s["accession"] for s in sequences]
            seqs_raw = [s["sequence"][:600] for s in sequences]
            min_len = min(len(s) for s in seqs_raw)
            seqs_raw = [s[:min_len] for s in seqs_raw]

            aln = MultipleSeqAlignment([
                SeqRecord(Seq(s), id=ids[i], description="")
                for i, s in enumerate(seqs_raw)
            ])

            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(aln)

            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)

            # Convert to JSON-serialisable nested dict
            def tree_to_dict(clade):
                node = {
                    "name": clade.name or "",
                    "branch_length": round(float(clade.branch_length or 0), 5)
                }
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
                "newick": newick
            }
        except Exception as e:
            return {"error": str(e), "taxa": len(sequences)}

    return await loop.run_in_executor(None, _phylo)


# ─────────────────────────────────────────────
# Step 8 — Protein Secondary Structure (Biopython ProtParam)
# ─────────────────────────────────────────────

async def run_protein_analysis(protein_seq: str, gene_id: str = "protein") -> Dict:
    loop = asyncio.get_event_loop()

    def _analyse():
        # Clean protein seq — remove stop codons and non-standard AAs
        clean = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', protein_seq.upper())
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

            # Secondary structure fraction (helix, turn, sheet)
            try:
                helix, turn, sheet = analysis.secondary_structure_fraction()
            except Exception:
                helix, turn, sheet = 0.0, 0.0, 0.0
            coil = max(0.0, 1.0 - helix - turn - sheet)

            # Extinction coefficient
            try:
                ec_cys, ec_no_cys = analysis.molar_extinction_coefficient()
            except Exception:
                ec_cys, ec_no_cys = 0, 0

            # Flexibility
            try:
                flex = analysis.flexibility()
                avg_flex = round(sum(flex) / len(flex), 4) if flex else 0
            except Exception:
                avg_flex = 0

            # Build SS string from fractions for display
            seq_len = len(clean)
            n_helix = int(helix * seq_len)
            n_turn = int(turn * seq_len)
            n_sheet = int(sheet * seq_len)
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
                    "helix": round(helix * 100, 1),
                    "turn": round(turn * 100, 1),
                    "sheet": round(sheet * 100, 1),
                    "coil": round(coil * 100, 1),
                },
                "ss_string": ss_string[:200],
                "amino_acid_composition": {
                    k: round(v * 100, 2) for k, v in sorted(aa_comp.items(), key=lambda x: -x[1])[:10]
                },
                "protein_sequence": clean[:300],
            }
        except Exception as e:
            return {"error": str(e), "gene_id": gene_id, "protein_sequence": clean[:100]}

    return await loop.run_in_executor(None, _analyse)


# ─────────────────────────────────────────────
# Full Pipeline Orchestrator
# ─────────────────────────────────────────────

async def run_full_pipeline(job_id: str, accessions: List[str], options: Dict):
    job = JOB_STORE[job_id]
    multiple = len(accessions) > 1

    try:
        # ── Step 1: Retrieve sequences ──
        job_update(job, "retrieval", "running", 5, f"Fetching {len(accessions)} sequence(s) from NCBI...")
        sequences = []
        for acc in accessions:
            try:
                seq_data = await fetch_sequence(acc.strip())
                sequences.append(seq_data)
                await asyncio.sleep(0.4)  # NCBI rate limit (3/sec without key, 10/sec with key)
            except Exception as e:
                job_update(job, "retrieval", "error", 5, f"Failed to fetch {acc}: {e}")
                raise

        job_update(job, "retrieval", "complete", 18,
                   f"Retrieved {len(sequences)} sequence(s) from NCBI GenBank",
                   {"sequences": [
                       {k: v for k, v in s.items() if k not in ("fasta",)} for s in sequences
                   ]})

        # ── Step 2: ORF Finder ──
        job_update(job, "orf", "running", 20, "Finding open reading frames (6 frames)...")
        orf_results = []
        for seq in sequences:
            orf_data = await run_orf_finder(seq["sequence"])
            orf_results.append({"accession": seq["accession"], "result": orf_data})
        job_update(job, "orf", "complete", 32,
                   f"Found {sum(r['result']['total_orfs'] for r in orf_results)} ORFs total",
                   {"orf_results": orf_results})

        # ── Step 3: Augustus gene prediction ──
        job_update(job, "augustus", "running", 34, "Submitting to Augustus gene predictor (web)...")
        augustus_results = []
        for seq in sequences:
            aug = await run_augustus(seq["sequence"], options.get("organism_model", "human"))
            augustus_results.append({"accession": seq["accession"], "result": aug})
        total_genes = sum(r["result"].get("genes_predicted", 0) for r in augustus_results)
        job_update(job, "augustus", "complete", 46,
                   f"Augustus predicted {total_genes} gene(s)",
                   {"augustus_results": augustus_results})

        # ── Step 4: BLASTN ──
        job_update(job, "blastn", "running", 48, "Running BLASTN against NCBI nt database...")
        blastn_results = []
        for seq in sequences:
            bn = await run_blastn(seq["sequence"], seq["accession"])
            blastn_results.append({"accession": seq["accession"], "result": bn})
            await asyncio.sleep(1)  # NCBI rate limit between BLAST calls
        job_update(job, "blastn", "complete", 60,
                   f"BLASTN complete — {sum(len(r['result'].get('hits',[])) for r in blastn_results)} hits",
                   {"blastn_results": blastn_results})

        # ── Step 5: BLASTP (use longest ORF protein) ──
        job_update(job, "blastp", "running", 62, "Running BLASTP on predicted proteins...")
        blastp_results = []
        for i, seq in enumerate(sequences):
            # Get best protein from ORF or Augustus
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
                   f"BLASTP complete — {sum(len(r['result'].get('hits',[])) for r in blastp_results)} hits",
                   {"blastp_results": blastp_results})

        # ── Step 6: MSA (only if multiple sequences) ──
        if multiple:
            job_update(job, "msa", "running", 76, "Running multiple sequence alignment...")
            msa_result = await run_msa(sequences)
            job_update(job, "msa", "complete", 82, "MSA complete", {"msa": msa_result})
        else:
            job_update(job, "msa", "skipped", 82, "Skipped — MSA requires 2+ sequences",
                       {"skipped": True, "reason": "Single sequence"})

        # ── Step 7: Phylogenetic tree (only if multiple) ──
        if multiple:
            job_update(job, "phylo", "running", 83, "Building phylogenetic tree (Neighbor-Joining)...")
            phylo_result = await run_phylo_tree(sequences)
            job_update(job, "phylo", "complete", 88, "Phylogenetic tree built", {"tree": phylo_result})
        else:
            job_update(job, "phylo", "skipped", 88, "Skipped — requires 2+ sequences",
                       {"skipped": True, "reason": "Single sequence"})

        # ── Step 8: Protein structure analysis ──
        job_update(job, "structure", "running", 89, "Analysing protein physicochemical properties...")
        structure_results = []
        for i, seq in enumerate(sequences):
            orfs = orf_results[i]["result"].get("orfs", [])
            if orfs:
                protein = orfs[0]["protein"]
                struct = await run_protein_analysis(protein, seq["accession"])
                structure_results.append(struct)
        job_update(job, "structure", "complete", 97,
                   f"Protein analysis complete for {len(structure_results)} sequence(s)",
                   {"structures": structure_results})

        job["status"] = "complete"
        job["progress"] = 100
        job["current_step"] = "done"
        job["completed_at"] = time.time()
        _save_job(job)

    except Exception as e:
        job["status"] = "error"
        job["error"] = str(e)
        job["current_step"] = "error"
        _save_job(job)
