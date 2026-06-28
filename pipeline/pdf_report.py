"""
PDF Report Generator — PhyloGenie
Generates a full analysis PDF using ReportLab
"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, PageBreak, KeepTogether
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT
from datetime import datetime
import time
from pathlib import Path
from typing import Dict
import config


# ── Colour palette ──
NEON_BLUE  = colors.HexColor("#00AAFF")
DARK_BG    = colors.HexColor("#0A0F1A")
MID_BG     = colors.HexColor("#111827")
BORDER     = colors.HexColor("#1E2D45")
MUTED      = colors.HexColor("#64748B")
WHITE      = colors.white
GREEN      = colors.HexColor("#22C55E")
YELLOW     = colors.HexColor("#F59E0B")
RED        = colors.HexColor("#EF4444")
PURPLE     = colors.HexColor("#8B5CF6")


def _styles():
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle("title", fontSize=22, textColor=NEON_BLUE,
                                 fontName="Helvetica-Bold", spaceAfter=6, alignment=TA_LEFT),
        "h1": ParagraphStyle("h1", fontSize=14, textColor=NEON_BLUE,
                               fontName="Helvetica-Bold", spaceBefore=14, spaceAfter=6),
        "h2": ParagraphStyle("h2", fontSize=11, textColor=WHITE,
                               fontName="Helvetica-Bold", spaceBefore=8, spaceAfter=4),
        "body": ParagraphStyle("body", fontSize=9, textColor=colors.HexColor("#CBD5E1"),
                                fontName="Helvetica", spaceAfter=4, leading=13),
        "mono": ParagraphStyle("mono", fontSize=8, textColor=colors.HexColor("#7DD3FC"),
                                fontName="Courier", spaceAfter=4, leading=12),
        "muted": ParagraphStyle("muted", fontSize=8, textColor=MUTED,
                                 fontName="Helvetica", spaceAfter=3),
        "label": ParagraphStyle("label", fontSize=8, textColor=MUTED,
                                  fontName="Helvetica-Bold", spaceAfter=2),
        "value": ParagraphStyle("value", fontSize=9, textColor=WHITE,
                                  fontName="Helvetica", spaceAfter=4),
    }


def _table_style(header_color=NEON_BLUE):
    return TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), MID_BG),
        ("TEXTCOLOR",  (0, 0), (-1, 0), header_color),
        ("FONTNAME",   (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",   (0, 0), (-1, 0), 8),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
        ("TOPPADDING",    (0, 0), (-1, 0), 6),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#0F1825")),
        ("TEXTCOLOR",  (0, 1), (-1, -1), colors.HexColor("#CBD5E1")),
        ("FONTNAME",   (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE",   (0, 1), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
         [colors.HexColor("#0F1825"), colors.HexColor("#111827")]),
        ("GRID",       (0, 0), (-1, -1), 0.3, BORDER),
        ("TOPPADDING",    (0, 1), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 5),
        ("LEFTPADDING",   (0, 0), (-1, -1), 8),
        ("RIGHTPADDING",  (0, 0), (-1, -1), 8),
    ])


def _hr():
    return HRFlowable(width="100%", thickness=0.5, color=BORDER, spaceAfter=8, spaceBefore=8)


def generate_pdf(job: Dict, output_path: str) -> str:
    S = _styles()
    steps = job.get("steps", {})
    accessions = job.get("accessions", [])
    created = datetime.fromtimestamp(job.get("created_at", time.time())).strftime("%Y-%m-%d %H:%M:%S")
    completed = datetime.fromtimestamp(job.get("completed_at", time.time())).strftime("%Y-%m-%d %H:%M:%S") \
                if job.get("completed_at") else "N/A"

    doc = SimpleDocTemplate(
        output_path,
        pagesize=A4,
        leftMargin=2*cm, rightMargin=2*cm,
        topMargin=2*cm, bottomMargin=2*cm,
        title=f"PhyloGenie Report — Job {job['id']}",
    )

    story = []

    # ── Cover ──
    story.append(Paragraph("🧬 PhyloGenie", S["title"]))
    story.append(Paragraph("Genomic Sequence Analysis Report", ParagraphStyle(
        "sub", fontSize=12, textColor=MUTED, fontName="Helvetica", spaceAfter=4)))
    story.append(_hr())

    meta = [
        ["Job ID", job["id"]],
        ["Accessions", ", ".join(accessions)],
        ["Status", job.get("status", "").upper()],
        ["Created", created],
        ["Completed", completed],
        ["Sequences", str(len(accessions))],
    ]
    t = Table(meta, colWidths=[4*cm, 12*cm])
    t.setStyle(_table_style())
    story.append(t)
    story.append(Spacer(1, 0.4*cm))

    # ── 1. Sequence Retrieval ──
    ret = steps.get("retrieval", {})
    if ret.get("status") == "complete":
        story.append(Paragraph("1. Sequence Retrieval", S["h1"]))
        story.append(_hr())
        for seq in ret.get("data", {}).get("sequences", []):
            story.append(Paragraph(f"Accession: {seq.get('accession','')}", S["h2"]))
            gb = seq.get("genbank", {})
            info = [
                ["Description", seq.get("description", "")[:100]],
                ["Organism", gb.get("organism", "")],
                ["Length", f"{seq.get('length', 0):,} bp"],
                ["Molecule Type", gb.get("molecule_type", "")],
                ["Topology", gb.get("topology", "")],
                ["Source", gb.get("source", "")[:80]],
                ["Date", gb.get("date", "")],
            ]
            t = Table(info, colWidths=[4*cm, 12*cm])
            t.setStyle(_table_style())
            story.append(t)
            story.append(Spacer(1, 0.2*cm))

            # First 300bp of sequence
            sequence = seq.get("sequence", "")
            if sequence:
                story.append(Paragraph("Nucleotide Sequence (first 300 bp):", S["label"]))
                story.append(Paragraph(sequence[:300] + ("..." if len(sequence) > 300 else ""), S["mono"]))

            # GenBank features
            features = gb.get("features", [])
            if features:
                story.append(Paragraph("GenBank Features:", S["label"]))
                rows = [["Type", "Location", "Key Qualifier"]]
                for feat in features[:10]:
                    q = feat.get("qualifiers", {})
                    key_q = next(iter(q.values()), "") if q else ""
                    if isinstance(key_q, list): key_q = key_q[0]
                    rows.append([feat.get("type",""), feat.get("location","")[:30], str(key_q)[:50]])
                t = Table(rows, colWidths=[3*cm, 5*cm, 8*cm])
                t.setStyle(_table_style())
                story.append(t)

            story.append(Spacer(1, 0.3*cm))

    # ── 2. ORF Finder ──
    orf_step = steps.get("orf", {})
    if orf_step.get("status") == "complete":
        story.append(Paragraph("2. ORF Finder", S["h1"]))
        story.append(_hr())
        for orf_data in orf_step.get("data", {}).get("orf_results", []):
            story.append(Paragraph(f"Accession: {orf_data['accession']}", S["h2"]))
            result = orf_data["result"]
            story.append(Paragraph(
                f"Total ORFs found: {result.get('total_orfs', 0)} | Sequence length: {result.get('sequence_length', 0):,} bp",
                S["body"]))
            orfs = result.get("orfs", [])
            if orfs:
                rows = [["Frame", "Strand", "Start", "End", "Length (aa)", "GC %", "Protein (first 60aa)"]]
                for orf in orfs[:10]:
                    rows.append([
                        orf.get("frame",""), orf.get("strand",""),
                        str(orf.get("start","")), str(orf.get("end","")),
                        str(orf.get("length_aa","")), f"{orf.get('gc_content','')}%",
                        orf.get("protein","")[:60]
                    ])
                t = Table(rows, colWidths=[1.5*cm, 1.5*cm, 2*cm, 2*cm, 2.5*cm, 1.5*cm, 5*cm])
                t.setStyle(_table_style())
                story.append(t)
            story.append(Spacer(1, 0.3*cm))

    # ── 3. Augustus ──
    aug_step = steps.get("augustus", {})
    if aug_step.get("status") == "complete":
        story.append(Paragraph("3. Augustus Gene Prediction", S["h1"]))
        story.append(_hr())
        for aug_data in aug_step.get("data", {}).get("augustus_results", []):
            story.append(Paragraph(f"Accession: {aug_data['accession']}", S["h2"]))
            r = aug_data["result"]
            story.append(Paragraph(
                f"Source: {r.get('source','')} | Organism model: {r.get('organism_model','')} | Genes predicted: {r.get('genes_predicted',0)}",
                S["body"]))
            genes = r.get("genes", [])
            if genes:
                rows = [["Gene ID", "Start", "End", "Strand", "Length (aa)", "Protein (first 60aa)"]]
                for g in genes[:8]:
                    rows.append([
                        g.get("id",""), str(g.get("start","")), str(g.get("end","")),
                        g.get("strand",""), str(g.get("length_aa","")),
                        g.get("protein","")[:60]
                    ])
                t = Table(rows, colWidths=[2.5*cm, 2.5*cm, 2.5*cm, 1.5*cm, 2.5*cm, 4.5*cm])
                t.setStyle(_table_style())
                story.append(t)
            story.append(Spacer(1, 0.3*cm))

    # ── 4. BLASTN ──
    bn_step = steps.get("blastn", {})
    if bn_step.get("status") == "complete":
        story.append(Paragraph("4. BLASTN Results (Nucleotide vs nt)", S["h1"]))
        story.append(_hr())
        for bn_data in bn_step.get("data", {}).get("blastn_results", []):
            story.append(Paragraph(f"Query: {bn_data['accession']}", S["h2"]))
            hits = bn_data["result"].get("hits", [])
            if hits:
                rows = [["Accession", "Description", "Score", "E-value", "Identity %", "Qry Cover %"]]
                for h in hits[:15]:
                    rows.append([
                        h.get("accession","")[:15],
                        h.get("title","")[:55],
                        str(h.get("bits","")),
                        f"{h.get('e_value',''):.2e}" if isinstance(h.get("e_value"), float) else str(h.get("e_value","")),
                        f"{h.get('identity_pct','')}%",
                        f"{h.get('query_cover','')}%",
                    ])
                t = Table(rows, colWidths=[2.5*cm, 6*cm, 1.5*cm, 2*cm, 2*cm, 2*cm])
                t.setStyle(_table_style())
                story.append(t)
            else:
                story.append(Paragraph("No significant hits found.", S["muted"]))
            story.append(Spacer(1, 0.3*cm))

    # ── 5. BLASTP ──
    bp_step = steps.get("blastp", {})
    if bp_step.get("status") == "complete":
        story.append(Paragraph("5. BLASTP Results (Protein vs nr)", S["h1"]))
        story.append(_hr())
        for bp_data in bp_step.get("data", {}).get("blastp_results", []):
            story.append(Paragraph(f"Query: {bp_data['accession']}", S["h2"]))
            hits = bp_data["result"].get("hits", [])
            if hits:
                rows = [["Accession", "Description", "Score", "E-value", "Identity %", "Positives %"]]
                for h in hits[:15]:
                    rows.append([
                        h.get("accession","")[:15],
                        h.get("title","")[:55],
                        str(h.get("bits","")),
                        f"{h.get('e_value',''):.2e}" if isinstance(h.get("e_value"), float) else str(h.get("e_value","")),
                        f"{h.get('identity_pct','')}%",
                        f"{h.get('positives_pct','')}%",
                    ])
                t = Table(rows, colWidths=[2.5*cm, 6*cm, 1.5*cm, 2*cm, 2*cm, 2*cm])
                t.setStyle(_table_style())
                story.append(t)
            else:
                story.append(Paragraph("No significant hits found.", S["muted"]))
            story.append(Spacer(1, 0.3*cm))

    # ── 6. MSA ──
    msa_step = steps.get("msa", {})
    if msa_step.get("status") == "complete":
        story.append(Paragraph("6. Multiple Sequence Alignment", S["h1"]))
        story.append(_hr())
        msa = msa_step.get("data", {}).get("msa", {})
        story.append(Paragraph(
            f"Method: {msa.get('method','')} | Sequences: {msa.get('num_sequences','')} | Alignment length: {msa.get('alignment_length','')} bp",
            S["body"]))
        for aln in msa.get("aligned", []):
            story.append(Paragraph(f"{aln['id']}  ({aln.get('identity_to_ref','')}% identity to ref)", S["label"]))
            story.append(Paragraph(aln.get("sequence","")[:120] + "...", S["mono"]))
        story.append(Spacer(1, 0.3*cm))

    # ── 7. Phylo ──
    phylo_step = steps.get("phylo", {})
    if phylo_step.get("status") == "complete":
        story.append(Paragraph("7. Phylogenetic Tree", S["h1"]))
        story.append(_hr())
        tree = phylo_step.get("data", {}).get("tree", {})
        story.append(Paragraph(
            f"Method: {tree.get('method','')} | Model: {tree.get('distance_model','')} | Taxa: {tree.get('taxa','')}",
            S["body"]))
        newick = tree.get("newick","")
        if newick:
            story.append(Paragraph("Newick Format:", S["label"]))
            story.append(Paragraph(newick[:300], S["mono"]))
        story.append(Spacer(1, 0.3*cm))

    # ── 8. Protein Structure ──
    struct_step = steps.get("structure", {})
    if struct_step.get("status") == "complete":
        story.append(Paragraph("8. Protein Physicochemical Analysis", S["h1"]))
        story.append(_hr())
        for st in struct_step.get("data", {}).get("structures", []):
            if st.get("error"):
                story.append(Paragraph(f"{st.get('gene_id','')} — Error: {st['error']}", S["muted"]))
                continue
            story.append(Paragraph(f"Gene: {st.get('gene_id','')}", S["h2"]))
            ss = st.get("secondary_structure", {})
            rows = [
                ["Property", "Value"],
                ["Protein Length", f"{st.get('protein_length','')} aa"],
                ["Molecular Weight", f"{st.get('molecular_weight','')} Da"],
                ["Isoelectric Point (pI)", str(st.get("isoelectric_point",""))],
                ["Instability Index", f"{st.get('instability_index','')} ({'stable' if st.get('is_stable') else 'unstable'})"],
                ["GRAVY Score", str(st.get("gravy",""))],
                ["Aromaticity", str(st.get("aromaticity",""))],
                ["Avg Flexibility", str(st.get("avg_flexibility",""))],
                ["Extinction Coeff (with Cys)", str(st.get("extinction_coeff_with_cys",""))],
                ["Alpha Helix %", f"{ss.get('helix','')}%"],
                ["Beta Sheet %",  f"{ss.get('sheet','')}%"],
                ["Turns %",       f"{ss.get('turn','')}%"],
                ["Coil %",        f"{ss.get('coil','')}%"],
            ]
            t = Table(rows, colWidths=[7*cm, 9*cm])
            t.setStyle(_table_style())
            story.append(t)

            prot_seq = st.get("protein_sequence","")
            if prot_seq:
                story.append(Paragraph("Protein Sequence (first 200 aa):", S["label"]))
                story.append(Paragraph(prot_seq[:200] + ("..." if len(prot_seq) > 200 else ""), S["mono"]))

            aa_comp = st.get("amino_acid_composition", {})
            if aa_comp:
                story.append(Paragraph("Top Amino Acid Composition:", S["label"]))
                aa_rows = [["AA", "Percentage"]] + [[k, f"{v}%"] for k, v in aa_comp.items()]
                t = Table(aa_rows, colWidths=[3*cm, 4*cm])
                t.setStyle(_table_style())
                story.append(t)

            story.append(Spacer(1, 0.3*cm))

    # ── Footer ──
    story.append(_hr())
    story.append(Paragraph(
        f"Generated by PhyloGenie · {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} · "
        f"NCBI Entrez · NCBI BLAST · Augustus · Biopython",
        ParagraphStyle("footer", fontSize=7, textColor=MUTED, fontName="Helvetica", alignment=TA_CENTER)
    ))

    doc.build(story)
    return output_path
