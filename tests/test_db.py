import importlib

import config
import db


def test_save_load_delete_roundtrip(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DB_PATH", str(tmp_path / "test.db"))
    monkeypatch.setattr(db, "_conn", None)  # force a fresh connection against the temp db

    job = {"id": "ABC123", "status": "complete", "created_at": 1.0, "completed_at": 2.0, "accessions": ["NM_000546"]}
    db.save_job(job)

    loaded = db.load_job("ABC123")
    assert loaded["id"] == "ABC123"
    assert loaded["status"] == "complete"

    all_jobs = db.load_all_jobs()
    assert any(j["id"] == "ABC123" for j in all_jobs)

    assert db.delete_job("ABC123") is True
    assert db.load_job("ABC123") is None
    assert db.delete_job("ABC123") is False


def test_delete_jobs_older_than(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DB_PATH", str(tmp_path / "test2.db"))
    monkeypatch.setattr(db, "_conn", None)

    db.save_job({"id": "OLD", "status": "complete", "created_at": 100.0, "accessions": []})
    db.save_job({"id": "NEW", "status": "complete", "created_at": 10_000_000.0, "accessions": []})

    removed = db.delete_jobs_older_than(1_000_000.0)
    assert removed == 1
    remaining_ids = {j["id"] for j in db.load_all_jobs()}
    assert remaining_ids == {"NEW"}
