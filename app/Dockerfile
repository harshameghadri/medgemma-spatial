# MedGemma Spatial Transcriptomics - Production Dockerfile
# Multi-stage build for efficient image size

# Stage 1: Base dependencies
FROM python:3.10-slim as base

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Stage 2: Python dependencies
FROM base as builder

WORKDIR /tmp

# Copy requirements
COPY app/requirements.txt /tmp/requirements.txt

# Install Python packages
RUN pip install --user -r requirements.txt

# Stage 3: Runtime
FROM base

# Create non-root user for security
RUN useradd -m -u 1000 medgemma && \
    mkdir -p /app /data /output && \
    chown -R medgemma:medgemma /app /data /output

# Copy installed packages from builder
COPY --from=builder /root/.local /home/medgemma/.local

# Set working directory
WORKDIR /app

# Copy application code
COPY --chown=medgemma:medgemma src/ /app/src/
COPY --chown=medgemma:medgemma notebooks/ /app/notebooks/
COPY --chown=medgemma:medgemma app/streamlit_app.py /app/
COPY --chown=medgemma:medgemma data/PanglaoDB_markers_27_Mar_2020.tsv /app/data/

# Switch to non-root user
USER medgemma

# Add local packages to PATH
ENV PATH=/home/medgemma/.local/bin:$PATH \
    PYTHONPATH=/app:$PYTHONPATH

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Expose Streamlit port
EXPOSE 8501

# Run Streamlit app
CMD ["streamlit", "run", "streamlit_app.py", \
     "--server.port=8501", \
     "--server.address=0.0.0.0", \
     "--server.headless=true", \
     "--server.fileWatcherType=none", \
     "--browser.gatherUsageStats=false"]
