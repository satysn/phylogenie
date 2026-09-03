"""
PhyloGenie API Routes
"""

import logging
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional

from fastapi import APIRouter, BackgroundTasks, HTTPException
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field, field_validator

import config
import db
from pipeline import JOB_STORE, _save_job, generate_pdf, run_full_pipeline
from pipeline.core import fetch_sequence, run_blastn, run_blastp
from pipeline.validators import ValidationError, validate_accession, validate_accessions, validate_organism_model, validate_tools

log = logging.getLogger("phylogenie.api")

router = APIRouter()


class JobRequest(BaseModel):
    accessions: List[str] = Field(..., min_length=1)
    options: Optional[Dict] = {}

    @field_validator("accessions")
    @classmethod
    def _non_empty_strings(cls, v):
        if not any(a.strip() for a in v):
            raise ValueError("At least one non-empty accession is required")
        return v


@router.get("/health")
async def health():
    problems = config.validate()
    return {
        "status": "ok" if not problems else "degraded",
        "version": "3.0.0",
        "tools": ["NCBI Entrez", "ORF Finder", "Augustus", "BLASTN", "BLASTP", "MSA", "Phylo NJ", "ProtParam"],
        "warnings": problems,
        "active_jobs": sum(1 for j in JOB_STORE.values() if j.get("status") == "running"),
        "total_jobs": len(JOB_STORE),
    }


@router.post("/jobs")
async def create_job(req: JobRequest, background_tasks: BackgroundTasks):
    try:
        accessions = validate_accessions(req.accessions, config.MAX_SEQUENCES)
    except ValidationError as e:
        raise HTTPException(400, str(e))

    options = dict(req.options or {})
    if "organism_model" in options:
        options["organism_model"] = validate_organism_model(options["organism_model"])
    tools = validate_tools(options.get("tools"))
    options["tools"] = tools

    def initial_step(name: str, requires_multiple: bool = False):
        if requires_multiple and not (len(accessions) > 1):
            return {"status": "skipped", "message": "Skipped — requires 2+ sequences", "data": {}}
        if not tools.get(name, True):
            return {"status": "skipped", "message": "Skipped — disabled for this job", "data": {}}
        return {"status": "pending", "message": "Queued", "data": {}}

    job_id = str(uuid.uuid4())[:8].upper()
    multiple = len(accessions) > 1

    job = {
        "id": job_id,
        "status": "running",
        "progress": 0,
        "current_step": "retrieval",
        "accessions": accessions,
        "options": options,
        "created_at": time.time(),
        "completed_at": None,
        "error": None,
        "multiple_sequences": multiple,
        "_cancel_requested": False,
        "steps": {
            "retrieval": {"status": "pending", "message": "Queued", "data": {}},
            "orf":       {"status": "pending", "message": "Queued", "data": {}},
            "augustus":  initial_step("augustus"),
            "blastn":    initial_step("blastn"),
            "blastp":    initial_step("blastp"),
            "msa":       initial_step("msa", requires_multiple=True),
            "phylo":     initial_step("phylo", requires_multiple=True),
            "structure": initial_step("structure"),
        },
    }
    JOB_STORE[job_id] = job
    _save_job(job)
    background_tasks.add_task(run_full_pipeline, job_id, accessions, options)
    log.info("Created job %s for accessions=%s", job_id, accessions)
    return {"job_id": job_id, "status": "running"}


@router.get("/jobs")
async def list_jobs(limit: int = 100, offset: int = 0):
    jobs = sorted(JOB_STORE.values(), key=lambda j: j["created_at"], reverse=True)
    page = jobs[offset: offset + max(1, min(limit, 500))]
    return {
        "total": len(jobs),
        "jobs": [
            {
                "id": j["id"], "status": j["status"], "progress": j["progress"],
                "accessions": j["accessions"], "created_at": j["created_at"],
                "completed_at": j.get("completed_at"), "error": j.get("error"),
                "multiple": j.get("multiple_sequences", False),
            } for j in page
        ],
    }


@router.get("/jobs/{job_id}")
async def get_job(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, f"Job {job_id} not found")
    return {k: v for k, v in job.items() if k not in ("steps", "_cancel_requested")} | {
        "steps": {k: {sk: sv for sk, sv in v.items() if sk != "data"} for k, v in job["steps"].items()}
    }


@router.get("/jobs/{job_id}/full")
async def get_job_full(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, f"Job {job_id} not found")
    return {k: v for k, v in job.items() if k != "_cancel_requested"}


@router.post("/jobs/{job_id}/cancel")
async def cancel_job(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, f"Job {job_id} not found")
    if job["status"] not in ("running",):
        raise HTTPException(400, f"Job is '{job['status']}' and cannot be cancelled")
    job["_cancel_requested"] = True
    return {"job_id": job_id, "cancel_requested": True}


@router.delete("/jobs/{job_id}")
async def delete_job(job_id: str):
    if job_id not in JOB_STORE:
        raise HTTPException(404, f"Job {job_id} not found")
    del JOB_STORE[job_id]
    db.delete_job(job_id)
    for ext in ["pdf"]:
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

    if not Path(pdf_path).exists():
        try:
            generate_pdf(job, pdf_path)
        except Exception as e:
            log.exception("PDF generation failed for job %s", job_id)
            raise HTTPException(500, f"PDF generation failed: {e}")

    return FileResponse(pdf_path, media_type="application/pdf", filename=f"PhyloGenie_{job_id}.pdf")


@router.get("/jobs/{job_id}/json")
async def download_json(job_id: str):
    job = JOB_STORE.get(job_id)
    if not job:
        raise HTTPException(404, "Job not found")
    clean = {k: v for k, v in job.items() if k != "_cancel_requested"}
    return JSONResponse(
        content=clean,
        headers={"Content-Disposition": f'attachment; filename="PhyloGenie_{job_id}.json"'},
    )


class SequenceRequest(BaseModel):
    accession: str


@router.post("/sequence/retrieve")
async def retrieve_sequence(req: SequenceRequest):
    """Quick single-sequence lookup — bypasses the full pipeline (seconds, not minutes)."""
    try:
        accession = validate_accession(req.accession)
    except ValidationError as e:
        raise HTTPException(400, str(e))
    try:
        result = await fetch_sequence(accession)
    except Exception as e:
        raise HTTPException(502, f"NCBI retrieval failed: {e}")
    return {k: v for k, v in result.items() if k != "fasta"} | {"fasta": result.get("fasta", "")}


class BlastRequest(BaseModel):
    sequence: str
    program: str = "blastn"  # blastn | blastp


@router.post("/blast")
async def standalone_blast(req: BlastRequest):
    """Run a standalone BLAST search against a pasted sequence (no job/pipeline)."""
    seq = req.sequence.strip().upper()
    if not seq:
        raise HTTPException(400, "sequence is required")
    if len(seq) > 5000:
        raise HTTPException(400, "sequence too long for standalone BLAST (max 5000 chars) — use a full pipeline job instead")

    if req.program == "blastp":
        result = await run_blastp(seq, gene_id="standalone")
    elif req.program == "blastn":
        result = await run_blastn(seq, accession="standalone")
    else:
        raise HTTPException(400, "program must be 'blastn' or 'blastp'")
    return result


@router.get("/sample-accessions")
async def samples():
    return {"samples": [
        {"accession": "NM_001301717", "gene": "BRCA1", "organism": "Homo sapiens", "desc": "DNA repair"},
        {"accession": "NM_000546",    "gene": "TP53",  "organism": "Homo sapiens", "desc": "Tumor suppressor"},
        {"accession": "NM_004333",    "gene": "BRAF",  "organism": "Homo sapiens", "desc": "Proto-oncogene"},
        {"accession": "NM_005228",    "gene": "EGFR",  "organism": "Homo sapiens", "desc": "Receptor kinase"},
        {"accession": "NM_000059",    "gene": "BRCA2", "organism": "Homo sapiens", "desc": "DNA repair"},
    ]}
