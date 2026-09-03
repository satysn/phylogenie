# 🧬 PhyloGenie v3 — Real Tools, Production-Ready

A genomic sequence analysis platform built on **real** NCBI, Augustus and Biopython tools —
sequence retrieval, ORF finding, gene prediction, BLAST, multiple sequence alignment,
phylogenetics and protein structure analysis, all from one web UI.

v3 is a hardening pass over the v2 "Real Tools Edition": input validation, retries and
bounded concurrency around flaky network calls, SQLite-backed job persistence, job
cancellation, rate limiting, structured logging, tests, and Docker/one-click deploy configs.

---

## What's new in v3

| Area | v2 | v3 |
|---|---|---|
| Secrets | NCBI API key hardcoded in `config.py` (and committed to git!) | Read from environment / `.env`, never committed |
| Job storage | One JSON file per job on disk | SQLite (`phylogenie.db`), single file, queryable, still human-inspectable |
| Input validation | None — bad accessions fail deep in the pipeline | Accessions validated up front with a clear 400 error |
| Concurrency | Sequences processed one at a time | Retrieval/ORF/structure run in parallel per job (bounded); jobs queue past `MAX_CONCURRENT_JOBS` |
| Reliability | A flaky NCBI call fails the whole job | Automatic retries with backoff on Entrez calls |
| Cancellation | Not possible once launched | `POST /api/jobs/{id}/cancel` stops the job between steps |
| MSA | Deprecated `Bio.pairwise2` | `Bio.Align.PairwiseAligner` (the modern replacement) |
| Rate limiting | None — one client can hammer the API | Per-IP sliding-window limiter |
| Deployability | `python main.py` only | Dockerfile, docker-compose, Procfile, render.yaml |
| Tests | None | pytest suite for the pure/pipeline logic |
| Path handling | Relative paths break if launched from another directory | All paths resolved relative to the app, not the CWD |

---

## Quick start (local)

```bash
cd "phylogenie v3"
python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # macOS/Linux

pip install -r requirements.txt
cp .env.example .env          # then edit .env — see below
python main.py
```

Open **http://localhost:8000**.

### Configure `.env`

Only `NCBI_EMAIL` is required — NCBI asks for a contact email on every Entrez/BLAST request.

```bash
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=                 # optional, but ~3-10x faster — free at ncbi.nlm.nih.gov/account/settings
```

> ⚠️ **If you used v1/v2 of this project**: the old `config.py` had a real NCBI API key
> committed directly to git history on GitHub. Treat that key as burned — [rotate/delete it
> in your NCBI account settings](https://www.ncbi.nlm.nih.gov/account/settings/) and generate
> a fresh one for `.env`. v3 never commits secrets — `.env` is gitignored.

---

## Running an analysis

1. **New Analysis** → add an accession (e.g. `NM_001301717`) or click a sample gene button
2. Add **2+ accessions** to also get MSA and a phylogenetic tree
3. **Run Full Pipeline** and watch each of the 8 stages complete in real time
4. **View Results**, then download the **PDF report** or raw **JSON**

### Expected runtime

| Step | Time |
|---|---|
| Sequence Retrieval (NCBI) | 3–10 sec |
| ORF Finder | < 1 sec |
| Augustus (web, falls back to ORF-based prediction if unavailable) | 30–120 sec |
| BLASTN / BLASTP (NCBI web) | 2–8 min each |
| MSA / Phylo Tree / Protein Analysis | < 5 sec |
| **Total (1 sequence)** | **~10–20 min** |
| **Total (4 sequences)** | **~30–50 min** — retrieval/ORF/structure run in parallel, BLAST stays sequential to respect NCBI's rate limits |

Jobs beyond `MAX_CONCURRENT_JOBS` (default 3) queue automatically rather than all hammering
NCBI at once. You can cancel a running job at any time from its status card.

### Quick single-sequence lookup

For a fast check without running the full pipeline:

```bash
curl -X POST http://localhost:8000/api/sequence/retrieve \
  -H "Content-Type: application/json" -d '{"accession": "NM_000546"}'
```

---

## API

Interactive docs at `/docs` (Swagger) and `/redoc`.

| Method | Endpoint | Description |
|---|---|---|
| `GET`  | `/api/health` | Service status, config warnings, active job count |
| `POST` | `/api/jobs` | Launch a new pipeline job |
| `GET`  | `/api/jobs` | List all jobs (paginated: `?limit=&offset=`) |
| `GET`  | `/api/jobs/{id}` | Job status (lightweight — no result payloads) |
| `GET`  | `/api/jobs/{id}/full` | Job status + full results |
| `POST` | `/api/jobs/{id}/cancel` | Request cancellation of a running job |
| `DELETE` | `/api/jobs/{id}` | Delete a job and its cached PDF |
| `GET`  | `/api/jobs/{id}/pdf` | Download the PDF report (generated on first request) |
| `GET`  | `/api/jobs/{id}/json` | Download the full job as JSON |
| `POST` | `/api/sequence/retrieve` | Quick single-sequence lookup (seconds) |
| `POST` | `/api/blast` | Standalone BLAST on a pasted sequence (no job) |
| `GET`  | `/api/sample-accessions` | Demo accession numbers |

---

## Deployment

### Docker (recommended)

```bash
cp .env.example .env   # fill in NCBI_EMAIL
docker compose up --build
```

Job data and the SQLite DB persist in named volumes across restarts.

### Render / Railway / Fly.io / Heroku-style platforms

- `render.yaml` is included for one-click Render deploys — set `NCBI_EMAIL` (and optionally
  `NCBI_API_KEY`) as environment variables in the dashboard (marked `sync: false` so they're
  never stored in the repo).
- `Procfile` works for any Heroku-style buildpack platform: `web: uvicorn main:app --host 0.0.0.0 --port $PORT`

### Bare metal / VM

```bash
pip install -r requirements.txt
NCBI_EMAIL=you@example.com PORT=8000 python main.py
```

Put it behind a reverse proxy (nginx/Caddy) for TLS in production. The app is a single
process — for real concurrent load, run multiple instances behind a load balancer with
`DB_PATH` pointed at a shared volume, or swap the SQLite job store for Postgres.

---

## Running tests

```bash
pip install -r requirements-dev.txt
pytest -v
```

Tests cover the pure pipeline logic (ORF finding, protein analysis, MSA, accession
validation, the SQLite job store) — no network calls, so they run in seconds and in CI.

---

## Project structure

```
phylogenie v3/
├── main.py                 # FastAPI app, lifespan, rate limiting, static/template serving
├── config.py                # Env-driven configuration (no secrets committed)
├── db.py                    # SQLite job store
├── api/
│   └── routes.py            # REST API endpoints
├── pipeline/
│   ├── core.py               # The 8-stage pipeline (retrieval → structure analysis)
│   ├── validators.py         # Accession/input validation
│   └── pdf_report.py         # PDF report generator (ReportLab)
├── templates/
│   └── index.html            # Full single-page frontend (vanilla JS + D3.js)
├── static/                   # Static assets
├── tests/                    # pytest suite (no network required)
├── Dockerfile / docker-compose.yml / .dockerignore
├── Procfile / render.yaml
├── .env.example / .gitignore
└── requirements.txt / requirements-dev.txt
```

## Tools used

| Tool | What it does | How |
|---|---|---|
| Biopython Entrez | Fetches FASTA + GenBank from NCBI | HTTP API, retried with backoff |
| ORF Finder | Finds ORFs in all 6 reading frames | Pure Python |
| Augustus | Predicts genes (exons/introns) | Web API, falls back to ORF-based prediction if unreachable |
| NCBI BLASTN / BLASTP | Nucleotide / protein homology search | Biopython `NCBIWWW` |
| Biopython `PairwiseAligner` | Multiple sequence alignment | Global pairwise alignment to a reference sequence |
| Biopython Phylo NJ | Neighbor-Joining phylogenetic tree | `DistanceTreeConstructor` |
| Biopython ProtParam | MW, pI, instability, secondary structure | `ProteinAnalysis` |
| ReportLab | PDF report generation | Local |

## Common issues

**`ModuleNotFoundError: No module named 'Bio'`**
→ `pip install -r requirements.txt`

**`BLAST is taking forever`**
→ Normal — NCBI's free web BLAST queues jobs and can take 2–8 min per sequence. Add an
`NCBI_API_KEY` in `.env` for higher throughput.

**`Augustus returned no genes`**
→ The Augustus web server can be slow or down; the pipeline automatically falls back to
ORF-based gene prediction so the job still completes.

**`429 Rate limit exceeded`**
→ You're hitting `/api/*` faster than `RATE_LIMIT_PER_MINUTE` (default 30/min per IP).
Raise it in `.env` if you're running trusted internal load.

**Health check shows `"status": "degraded"`**
→ `NCBI_EMAIL` isn't set. Real jobs will still attempt to run, but NCBI may throttle or
reject requests without a contact email.
