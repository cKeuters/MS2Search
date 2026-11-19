import os, tempfile, time, uuid
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from fastapi import FastAPI, HTTPException, Response, UploadFile
from dreams.api import dreams_embeddings

app = FastAPI()
executor = ThreadPoolExecutor(max_workers=1)
jobs: dict[str, dict] = {}  
JOB_TTL = 3600   


def _cleanup():
    now = time.time()
    for jid in [j for j, v in jobs.items() if now - v["created"] > JOB_TTL]:
        jobs.pop(jid, None)


def _get(job_id: str) -> dict:
    job = jobs.get(job_id)
    if job is None:
        raise HTTPException(404, "unknown job")
    return job


def _run(job: dict, path: str):
    job["status"], job["started"] = "running", time.time()
    try:
        job["result"] = np.asarray(dreams_embeddings(path), dtype="float32")
        job["status"] = "done"
    except Exception as e:
        job["status"], job["error"] = "error", repr(e)
    finally:
        os.remove(path)


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/jobs")
async def create_job(file: UploadFile):
    _cleanup()
    data = await file.read()
    with tempfile.NamedTemporaryFile(suffix=".mgf", delete=False) as tmp:
        tmp.write(data)
    job_id = uuid.uuid4().hex
    job = {"status": "queued", "created": time.time()}
    jobs[job_id] = job
    executor.submit(_run, job, tmp.name)
    return {"job_id": job_id}


@app.get("/jobs/{job_id}")
def job_status(job_id: str):
    job = _get(job_id)
    ahead = sum(1 for j in jobs.values()
                if j["status"] in ("queued", "running") and j["created"] < job["created"])
    elapsed = int(time.time() - job["started"]) if job["status"] == "running" else None
    return {"status": job["status"], "error": job.get("error"),
            "ahead": ahead, "elapsed": elapsed}


@app.get("/jobs/{job_id}/result")
def job_result(job_id: str):
    job = _get(job_id)
    if job["status"] != "done":
        raise HTTPException(409, "not finished")
    arr = job["result"]
    return Response(content=arr.tobytes(), media_type="application/octet-stream",
                    headers={"X-Rows": str(arr.shape[0]), "X-Cols": str(arr.shape[1])})



