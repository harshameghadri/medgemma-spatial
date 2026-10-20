# MedGemma Spatial Transcriptomics - production image.
# Build from the repository root:  docker build -f app/Dockerfile .
# The same file is used by CI and by the HuggingFace Spaces deploy.

FROM python:3.12-slim AS base

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

RUN apt-get update && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*

FROM base AS builder
COPY requirements.txt /tmp/requirements.txt
RUN pip install --prefix=/install -r /tmp/requirements.txt

FROM base
RUN useradd -m -u 1000 medgemma
COPY --from=builder /install /usr/local

WORKDIR /app
COPY --chown=medgemma:medgemma src/ /app/src/
COPY --chown=medgemma:medgemma app/streamlit_app.py /app/streamlit_app.py
COPY --chown=medgemma:medgemma .streamlit/ /app/.streamlit/
COPY --chown=medgemma:medgemma data/PanglaoDB_markers_27_Mar_2020.tsv /app/data/

USER medgemma
ENV HOME=/home/medgemma \
    PYTHONPATH=/app \
    HF_HOME=/home/medgemma/.cache/huggingface

HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

EXPOSE 8501
CMD ["streamlit", "run", "streamlit_app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--server.fileWatcherType=none", \
     "--browser.gatherUsageStats=false", \
     "--server.maxUploadSize=1024"]
