"""
PhyloGenie v3 — Configuration

All settings are read from environment variables (optionally via a local
.env file). Copy .env.example to .env and fill in your values — never
commit real secrets to git.
"""

import os
from pathlib import Path

# Load a local .env file if python-dotenv is available. This is optional —
# in a real deployment (Docker, Render, Railway, ...) env vars are injected
# by the platform and no .env file is needed.
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass


def _bool(name: str, default: bool) -> bool:
    val = os.getenv(name)
    if val is None:
        return default
    return val.strip().lower() in ("1", "true", "yes", "on")


def _int(name: str, default: int) -> int:
    try:
        return int(os.getenv(name, default))
    except (TypeError, ValueError):
        return default


# ── NCBI Entrez / BLAST ──
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "").strip()
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()

# ── Storage ──
# Paths are always resolved relative to this file's directory, not the
# process's current working directory, so the app behaves the same
# whether it's launched as `python main.py`, from a Docker CMD, or from
# a process manager that starts it from a different cwd.
BASE_DIR = Path(__file__).resolve().parent


def _resolve(env_var: str, default_name: str) -> str:
    raw = os.getenv(env_var, "").strip()
    if not raw:
        return str(BASE_DIR / default_name)
    p = Path(raw)
    return str(p if p.is_absolute() else BASE_DIR / p)


OUTPUT_DIR = _resolve("OUTPUT_DIR", "outputs")
DB_PATH = _resolve("DB_PATH", "phylogenie.db")

# ── Job limits ──
MAX_SEQUENCES = _int("MAX_SEQUENCES", 10)
JOB_RETENTION_DAYS = _int("JOB_RETENTION_DAYS", 30)
MAX_CONCURRENT_JOBS = _int("MAX_CONCURRENT_JOBS", 3)

# ── BLAST settings ──
BLAST_HITLIST_SIZE = _int("BLAST_HITLIST_SIZE", 20)
BLAST_EXPECT = float(os.getenv("BLAST_EXPECT", "0.001"))
BLAST_TIMEOUT_SECONDS = _int("BLAST_TIMEOUT_SECONDS", 180)

# ── Augustus web API ──
AUGUSTUS_URL = os.getenv("AUGUSTUS_URL", "https://bioinf.uni-greifswald.de/augustus/submission.php")
AUGUSTUS_POLL_ATTEMPTS = _int("AUGUSTUS_POLL_ATTEMPTS", 12)
AUGUSTUS_POLL_INTERVAL = _int("AUGUSTUS_POLL_INTERVAL", 10)

# ── Server ──
HOST = os.getenv("HOST", "0.0.0.0")
PORT = _int("PORT", 8000)
RELOAD = _bool("RELOAD", False)
CORS_ORIGINS = [o.strip() for o in os.getenv("CORS_ORIGINS", "*").split(",") if o.strip()]
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()

# ── Rate limiting (requests per minute per client IP) ──
# The frontend polls /api/jobs/{id} every 2s (30 req/min) while a job runs,
# plus occasional history/health calls on top — 30 was too tight and caused
# the poll's own 429s to visibly reset the progress panel. This still caps
# obvious abuse/scraping while giving normal single-job polling headroom,
# including several teammates polling their own jobs behind the same IP.
RATE_LIMIT_PER_MINUTE = _int("RATE_LIMIT_PER_MINUTE", 120)


def validate() -> list[str]:
    """Return a list of human-readable configuration problems (empty = OK)."""
    problems = []
    if not NCBI_EMAIL:
        problems.append(
            "NCBI_EMAIL is not set. NCBI requires a contact email for Entrez/BLAST "
            "requests — set it in your .env file or environment before running real jobs."
        )
    return problems
