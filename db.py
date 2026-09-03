"""
SQLite-backed job store.

Jobs are kept in an in-memory dict for fast access during a run, and
mirrored to SQLite (single file, no external DB needed) so they survive
restarts and can be listed/queried without re-reading every JSON blob
from disk. The full job document (steps, results, etc.) is stored as a
JSON column — simple and adequate for a single-instance deployment.
"""

import json
import sqlite3
import threading
from pathlib import Path
from typing import Any, Optional

import config

_lock = threading.Lock()
_conn: Optional[sqlite3.Connection] = None


def _get_conn() -> sqlite3.Connection:
    global _conn
    if _conn is None:
        Path(config.DB_PATH).parent.mkdir(parents=True, exist_ok=True)
        _conn = sqlite3.connect(config.DB_PATH, check_same_thread=False)
        _conn.execute("""
            CREATE TABLE IF NOT EXISTS jobs (
                id TEXT PRIMARY KEY,
                status TEXT NOT NULL,
                created_at REAL NOT NULL,
                completed_at REAL,
                accessions TEXT NOT NULL,
                data TEXT NOT NULL
            )
        """)
        _conn.execute("CREATE INDEX IF NOT EXISTS idx_jobs_created_at ON jobs(created_at DESC)")
        _conn.commit()
    return _conn


def save_job(job: dict[str, Any]) -> None:
    conn = _get_conn()
    with _lock:
        conn.execute(
            """INSERT INTO jobs (id, status, created_at, completed_at, accessions, data)
               VALUES (?, ?, ?, ?, ?, ?)
               ON CONFLICT(id) DO UPDATE SET
                 status=excluded.status,
                 completed_at=excluded.completed_at,
                 data=excluded.data""",
            (
                job["id"],
                job.get("status", "unknown"),
                job.get("created_at", 0),
                job.get("completed_at"),
                json.dumps(job.get("accessions", [])),
                json.dumps(job, default=str),
            ),
        )
        conn.commit()


def load_job(job_id: str) -> Optional[dict[str, Any]]:
    conn = _get_conn()
    row = conn.execute("SELECT data FROM jobs WHERE id = ?", (job_id,)).fetchone()
    return json.loads(row[0]) if row else None


def load_all_jobs() -> list[dict[str, Any]]:
    conn = _get_conn()
    rows = conn.execute("SELECT data FROM jobs ORDER BY created_at DESC").fetchall()
    return [json.loads(r[0]) for r in rows]


def delete_job(job_id: str) -> bool:
    conn = _get_conn()
    with _lock:
        cur = conn.execute("DELETE FROM jobs WHERE id = ?", (job_id,))
        conn.commit()
        return cur.rowcount > 0


def delete_jobs_older_than(cutoff_timestamp: float) -> int:
    conn = _get_conn()
    with _lock:
        cur = conn.execute("DELETE FROM jobs WHERE created_at < ?", (cutoff_timestamp,))
        conn.commit()
        return cur.rowcount
