FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
      git build-essential curl ca-certificates \
      libxrender1 libxext6 libsm6 libx11-6 libglib2.0-0 \
 && rm -rf /var/lib/apt/lists/*

RUN curl -L -o /opt/embedding_model.ckpt \
  "https://zenodo.org/records/10997887/files/embedding_model.ckpt?download=1" \
&& echo "4f6913c242b219b75e667dd1ed5add59  /opt/embedding_model.ckpt" | md5sum -c -

RUN pip install --no-cache-dir torch==2.2.1 --index-url https://download.pytorch.org/whl/cpu
RUN git clone https://github.com/pluskal-lab/DreaMS.git /opt/DreaMS \
 && pip install --no-cache-dir -e /opt/DreaMS \
 && pip install --no-cache-dir fastapi uvicorn python-multipart

ENV DREAMS_CKPT=/opt/embedding_model.ckpt

COPY app.py /srv/app.py
WORKDIR /srv
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]