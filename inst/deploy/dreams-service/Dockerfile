FROM python:3.11-slim
RUN apt-get update && apt-get install -y --no-install-recommends git build-essential \
 && rm -rf /var/lib/apt/lists/*
 
RUN git clone https://github.com/pluskal-lab/DreaMS.git /opt/DreaMS \
 && pip install --no-cache-dir -e /opt/DreaMS \
 && pip install --no-cache-dir fastapi uvicorn python-multipart

COPY app.py /srv/app.py
WORKDIR /srv
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "1"]