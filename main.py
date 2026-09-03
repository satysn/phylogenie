"""
PhyloGenie v3 — Real Tools Edition
FastAPI + Biopython + NCBI + Augustus
"""

import asyncio
import logging
import sys
import time
from collections import defaultdict, deque
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

# Some Windows terminals default to a legacy console codepage (cp1252) that
# can't encode emoji — reconfigure to UTF-8 where possible so startup
# messages and logs never crash the process.
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        try:
            _stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass

import config
import db
from api.routes import router
from pipeline import JOB_STORE

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL, logging.INFO),
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
)
log = logging.getLogger("phylogenie")


async def _cleanup_old_jobs():
    """Periodically purge jobs older than JOB_RETENTION_DAYS."""
    while True:
        try:
            cutoff = time.time() - config.JOB_RETENTION_DAYS * 86400
            stale_ids = [jid for jid, j in JOB_STORE.items() if j.get("created_at", 0) < cutoff]
            for jid in stale_ids:
                JOB_STORE.pop(jid, None)
                pdf_path = Path(config.OUTPUT_DIR) / f"job_{jid}.pdf"
                if pdf_path.exists():
                    pdf_path.unlink()
            removed = db.delete_jobs_older_than(cutoff)
            if removed:
                log.info("Cleanup: purged %d job(s) older than %d days", removed, config.JOB_RETENTION_DAYS)
        except Exception:
            log.exception("Cleanup task failed")
        await asyncio.sleep(6 * 3600)  # every 6 hours


@asynccontextmanager
async def lifespan(app: FastAPI):
    problems = config.validate()
    for p in problems:
        log.warning(p)
    log.info("PhyloGenie v3 starting — %d job(s) loaded from disk", len(JOB_STORE))
    cleanup_task = asyncio.create_task(_cleanup_old_jobs())
    yield
    cleanup_task.cancel()


app = FastAPI(title="PhyloGenie", version="3.0.0", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_methods=["*"],
    allow_headers=["*"],
)
app.add_middleware(GZipMiddleware, minimum_size=1000)

# ── Simple in-memory rate limiter (per client IP, sliding 60s window) ──
_request_log: dict[str, deque] = defaultdict(deque)


@app.middleware("http")
async def rate_limit(request: Request, call_next):
    if request.url.path.startswith("/api/"):
        client_ip = request.client.host if request.client else "unknown"
        now = time.time()
        window = _request_log[client_ip]
        while window and now - window[0] > 60:
            window.popleft()
        if len(window) >= config.RATE_LIMIT_PER_MINUTE:
            return JSONResponse(
                status_code=429,
                content={"detail": "Rate limit exceeded — please slow down."},
            )
        window.append(now)
        # Bound memory growth: occasionally drop IPs with no recent activity
        # instead of letting _request_log accumulate one entry per client
        # forever on a long-running public instance.
        if len(_request_log) > 10_000:
            for ip in [k for k, w in _request_log.items() if not w or now - w[-1] > 60]:
                del _request_log[ip]
    return await call_next(request)


app.include_router(router, prefix="/api")

Path(config.OUTPUT_DIR).mkdir(exist_ok=True, parents=True)
_STATIC_DIR = config.BASE_DIR / "static"
_STATIC_DIR.mkdir(exist_ok=True)
app.mount("/static", StaticFiles(directory=str(_STATIC_DIR)), name="static")

_INDEX_PATH = config.BASE_DIR / "templates" / "index.html"
_INDEX_HTML = None


@app.get("/", response_class=HTMLResponse)
async def root():
    global _INDEX_HTML
    if _INDEX_HTML is None or config.RELOAD:
        _INDEX_HTML = _INDEX_PATH.read_text(encoding="utf-8")
    return _INDEX_HTML


if __name__ == "__main__":
    import uvicorn
    print("\n🧬 PhyloGenie v3 starting...")
    print(f"📖 Open: http://localhost:{config.PORT}")
    print(f"📚 API docs: http://localhost:{config.PORT}/docs\n")
    uvicorn.run("main:app", host=config.HOST, port=config.PORT, reload=config.RELOAD)
