# DeepCatch Bioinformatics Validation Suite
#
# Python 3.11-slim base image with all dependencies pinned to exact versions
# for full computational reproducibility.
#
# Build:
#   docker build -t deepcatch-validation:latest .
#
# Run:
#   docker run --rm -v $(pwd)/results:/app/results \
#       deepcatch-validation:latest python run_bioinfo_validation.py
#
# Entry point runs the full 10-module bioinformatics validation suite.

FROM python:3.11-slim

LABEL maintainer="DeepCatch Project"
LABEL description="Bioinformatics validation suite for deep-learning cancer screening"
LABEL version="1.0.0"

# ── System dependencies ──────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# ── Python environment ───────────────────────────────────────────────────
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    OMP_NUM_THREADS=2 \
    OPENBLAS_NUM_THREADS=2 \
    MKL_NUM_THREADS=2

WORKDIR /app

# ── Dependencies (all pinned) ────────────────────────────────────────────
# SciPy ecosystem
RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    scipy==1.11.4 \
    scikit-learn==1.4.2 \
    pandas==2.2.2

# Additional utilities
RUN pip install --no-cache-dir \
    matplotlib==3.8.4 \
    seaborn==0.13.2 \
    joblib==1.4.2 \
    tqdm==4.66.4 \
    pyyaml==6.0.1

# ── Copy source code ─────────────────────────────────────────────────────
COPY validation_framework.py /app/
COPY validation/ /app/validation/
COPY reproducibility/ /app/reproducibility/
COPY run_bioinfo_validation.py /app/

# ── Create output directory ──────────────────────────────────────────────
RUN mkdir -p /app/results

# ── Verify imports ───────────────────────────────────────────────────────
RUN python -c "import validation_framework; print('validation_framework OK')" && \
    python -c "from validation import *; print('validation modules OK')"

# ── Entry point ──────────────────────────────────────────────────────────
# Default: run ALL validation modules
CMD ["python", "run_bioinfo_validation.py"]
