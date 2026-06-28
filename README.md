# 🧬 PhyloGenie v2 — Real Tools Edition

A full genomic analysis platform using **real** NCBI, Augustus, BLAST and Biopython tools.

---

## ⚡ Windows Setup (Step by Step)

### Step 1 — Install Python
1. Go to https://www.python.org/downloads/
2. Download Python **3.11** or newer
3. Run the installer — ✅ check **"Add Python to PATH"**
4. Click Install Now

### Step 2 — Get your NCBI email & API key (5 minutes, free)
1. Go to https://www.ncbi.nlm.nih.gov/account/
2. Click **Sign In** → **Create account** (free)
3. After signing in, go to: https://www.ncbi.nlm.nih.gov/account/settings/
4. Scroll to **API Key Management** → click **Create an API Key**
5. Copy your API key

### Step 3 — Configure the app
Open `config.py` in Notepad and fill in:
```python
NCBI_EMAIL = "your_actual_email@gmail.com"
NCBI_API_KEY = "your_api_key_here"
```

### Step 4 — Open Command Prompt in the phylogenie2 folder
- Open the `phylogenie2` folder in File Explorer
- Click the address bar at the top
- Type `cmd` and press Enter

### Step 5 — Install dependencies
```cmd
pip install -r requirements.txt
```
This downloads: FastAPI, Uvicorn, Biopython, ReportLab, requests
Takes 1–3 minutes.

### Step 6 — Run the app
```cmd
python main.py
```
You'll see:
```
🧬 PhyloGenie starting...
📖 Open: http://localhost:8000
```

### Step 7 — Open your browser
Go to → **http://localhost:8000**

To stop the server: press `Ctrl + C` in the command prompt.

---

## 🔬 How to Run an Analysis

1. Click **New Analysis** in the sidebar
2. Type an accession number (e.g. `NM_001301717`) and click **Add**, OR
   click one of the **sample buttons** (BRCA1, TP53, BRAF, EGFR)
3. For MSA and Phylogenetic tree — add **2 or more** accessions
4. Click **Run Full Pipeline**
5. Watch each step run in real time (~5–15 minutes for BLAST steps)
6. When complete — click **View Results** to explore each section
7. Click **⬇ Download PDF** for the full analysis report

---

## ⏱ Expected Runtime

| Step | Time |
|------|------|
| Sequence Retrieval (NCBI) | 3–10 sec |
| ORF Finder | < 1 sec |
| Augustus (web) | 30–120 sec |
| BLASTN (NCBI web) | 2–8 min |
| BLASTP (NCBI web) | 2–8 min |
| MSA (Biopython) | 2–10 sec |
| Phylo Tree | < 1 sec |
| Protein Analysis | < 1 sec |
| **Total (single seq)** | **~10–20 min** |
| **Total (4 sequences)** | **~30–50 min** |

> BLAST is slow because NCBI rate-limits free web BLAST. With an API key it's ~3× faster.

---

## 📁 Project Structure

```
phylogenie2/
├── main.py              ← Start the server (run this)
├── config.py            ← Put your NCBI email + API key here
├── requirements.txt     ← pip install -r requirements.txt
├── api/
│   └── routes.py        ← REST API endpoints
├── pipeline/
│   ├── core.py          ← All real tool integrations
│   └── pdf_report.py    ← PDF generator (ReportLab)
├── templates/
│   └── index.html       ← Full frontend
└── outputs/             ← Auto-created: job JSONs + PDFs saved here
```

---

## 🛠 Tools Used

| Tool | What it does | How |
|------|-------------|-----|
| **Biopython Entrez** | Fetches FASTA + GenBank from NCBI | HTTP API |
| **ORF Finder** | Finds ORFs in all 6 reading frames | Pure Python |
| **Augustus** | Predicts genes (exons/introns) | Web API |
| **NCBI BLASTN** | Nucleotide homology search | Biopython `NCBIWWW` |
| **NCBI BLASTP** | Protein homology search | Biopython `NCBIWWW` |
| **Biopython MSA** | Multiple sequence alignment | `pairwise2` |
| **Biopython Phylo NJ** | Neighbor-Joining phylogenetic tree | `DistanceTreeConstructor` |
| **Biopython ProtParam** | MW, pI, instability, 2° structure | `ProteinAnalysis` |
| **ReportLab** | PDF generation | Local |

---

## ❓ Common Issues

**`ModuleNotFoundError: No module named 'Bio'`**
→ Run: `pip install biopython`

**`BLAST is taking forever`**
→ Normal. NCBI free BLAST can take 5–10 min. Add your API key to `config.py`.

**`Augustus returned no genes`**
→ Augustus web server can be slow or down. The pipeline falls back to ORF-based prediction automatically.

**`Connection error on NCBI`**
→ Check your internet. NCBI may be temporarily rate-limiting you (wait 30 min).
