# app.py
import tempfile, numpy as np
from fastapi import FastAPI, UploadFile
from fastapi.responses import Response
from dreams.api import dreams_embeddings
import os

app = FastAPI()
CKPT = os.environ["DREAMS_CKPT"]

@app.get("/health")
def health():
    return {"status": "ok"}

@app.post("/embed")
async def embed(file: UploadFile):
    with tempfile.NamedTemporaryFile(suffix=".mgf") as tmp:
        tmp.write(await file.read())
        tmp.flush()
        embs =  dreams_embeddings(tmp.name,  model_pth =CKPT)
    arr = np.asarray(embs, dtype="float32")
    return Response(
        content=arr.tobytes(),
        media_type="application/octet-stream",
        headers={"X-Rows": str(arr.shape[0]), "X-Cols": str(arr.shape[1])},
    )