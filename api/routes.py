"""
PhyloGenie API Routes
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import List, Optional, Dict
import uuid, time, os
from pathlib import Path

from pipeline import run_full_pipeline, JOB_STORE, generate_pdf, _save_job
import config

router = APIRouter()


class JobRequest(BaseModel):
    accessions: List[str]
    options: Optional[Dict] = {}


@router.get("/health")
async def health():
    return {"status": "ok", "version": "2.0.0", "tools": ["NCBI Entrez", "ORF Finder", "Augustus", "BLASTN", "BLASTP", "MSA", "Phylo NJ", "ProtParam"]}


@router.post("/jobs")
async def create_job(req: JobRequest, background_tasks: BackgroundTasks):
    accessions = [a.strip().upper() for a in req.accessions if a.strip()]
    if not accessions:
        raise HTTPException(400, "At least one accession required")
    if len(accessions) > config.MAX_SEQUENCES:
        raise HTTPException(400, f"Maximum {config.MAX_SEQUENCES} sequences per job")

    job_id = str(uuid.uuid4())[:8].upper()
    multiple = len(accessions) > 1

    job = {
        "id": job_id,
        "status": "running",
        "progress": 0,
        "current_step": "retrieval",
        "accessions": accessions,
        "options": req.options or {},
        "created_at": time.time(),
        "completed_at": None,
        "error": None,
        "multiple_sequences": multiple,
        "steps": {
            "retrieval": {"status": "pending", "message": "Queued", "data": {}},
            "orf":       {"status": "pending", "message": "Queued", "data": {}},
            "augustus":  {"status": "pending", "message": "Queued", "data": {}},
            "blastn":    {"status": "pending", "message": "Queued", "data": {}},
            "blastp":    {"status": "pending", "message": "Queued", "data": {}},
            "msa":       {"status": "pending" if multiple else "skipped",
                          "message": "Queued" if multiple else "Skipped — single sequence", "data": {}},
            "phylo":     {"status": "pending" if multiple else "skipped",
                          "message": "Queued" if multiple else "Skipped — single sequence", "data": {}},
            "structure": {"status": "pending", "message": "Queued", "data": {}},
        }
    }
    JOB_STORE[job_id] = job
    _save_job(job)
    background_tasks.add_task(run_full_pipeline, job_id, accessions, req.options or {})
    return {"job_id": job_id, "status": "running"}


@router.get("/jobs")
async def list_jobs():
    jobs = sorted(JOB_STORE.values(), key=lambda j: j["created_at"], reverse=True)
    return {"jobs": [
        {
            "id": j["id"], "status": j["status"], "progress": j["progress"],
            "accessions": j["accessions"], "created_at": j["created_at"],
            "completed_at": j.get("completed_at"), "error": j.get("error"),
            "multiple": j.get("multiple_sequences", False)
        } for j in jobs
    ]}


@router.get("/jobs/{job_id}")
async def get_job(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, f"Job {job_id} not found")
    # Return without large data blobs for status polling
    return {k: v for k, v in job.items() if k != "steps"} | {
        "steps": {
            k: {sk: sv for sk, sv in v.items() if sk != "data"}
            for k, v in job["steps"].items()
        }
    }


@router.get("/jobs/{job_id}/full")
async def get_job_full(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, f"Job {job_id} not found")
    return job


@router.delete("/jobs/{job_id}")
async def delete_job(job_id: str):
    if job_id not in JOB_STORE:
        raise HTTPException(404, f"Job {job_id} not found")
    del JOB_STORE[job_id]
    # Remove disk files
    for ext in ["json", "pdf"]:
        p = Path(config.OUTPUT_DIR) / f"job_{job_id}.{ext}"
        if p.exists():
            p.unlink()
    return {"deleted": job_id}


@router.get("/jobs/{job_id}/pdf")
async def download_pdf(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    if job["status"] != "complete":
        raise HTTPException(400, "Job not complete yet")

    pdf_path = str(Path(config.OUTPUT_DIR) / f"job_{job_id}.pdf")

    # Generate if not cached
    if not Path(pdf_path).exists():
        try:
            generate_pdf(job, pdf_path)
        except Exception as e:
            raise HTTPException(500, f"PDF generation failed: {e}")

    return FileResponse(
        pdf_path,
        media_type="application/pdf",
        filename=f"PhyloGenie_{job_id}.pdf"
    )


@router.get("/sample-accessions")
async def samples():
    return {"samples": [
        {"accession": "NM_001301717", "gene": "BRCA1", "organism": "Homo sapiens", "desc": "DNA repair"},
        {"accession": "NM_000546",    "gene": "TP53",  "organism": "Homo sapiens", "desc": "Tumor suppressor"},
        {"accession": "NM_004333",    "gene": "BRAF",  "organism": "Homo sapiens", "desc": "Proto-oncogene"},
        {"accession": "NM_005228",    "gene": "EGFR",  "organism": "Homo sapiens", "desc": "Receptor kinase"},
        {"accession": "NM_000059",    "gene": "BRCA2", "organism": "Homo sapiens", "desc": "DNA repair"},
    ]}
