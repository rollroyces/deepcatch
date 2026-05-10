"""
DeepCatch Fragmentomics API — FastAPI Application
==================================================

Local REST API for cfDNA fragment end motif cancer detection.

.. rubric:: Endpoints

- ``GET  /health`` — Health check with model/device info
- ``POST /predict`` — Predict cancer from DNA sequence or FASTA file
- ``GET  /docs`` — Interactive Swagger/OpenAPI documentation

.. rubric:: Quick Start

.. code-block:: bash

    pip install fastapi uvicorn torch
    cd deepcatch
    python -m api.main

    # Or with uvicorn directly:
    uvicorn api.main:app --host 0.0.0.0 --port 8000 --reload

.. rubric:: Design Notes for Apple Silicon

- The model and tokenizer are loaded once at startup into shared memory.
- Requests are processed synchronously (model is <1 MB — inference is fast).
- MPS device is auto-detected and used when available.
- All inputs are validated via Pydantic schemas.
- OpenAPI docs are auto-generated at ``/docs``.

.. rubric:: Integration With Agentic Frameworks

The OpenAPI spec at ``/openapi.json`` enables zero-code integration
with LangChain tools, OpenAI function calling, and other agent frameworks.

Example LangChain tool:

.. code-block:: python

    from langchain.tools import StructuredTool
    import requests

    def predict_cfDNA(sequence: str) -> dict:
        resp = requests.post(
            "http://localhost:8000/predict",
            json={"sequence": sequence, "method": "ensemble"}
        )
        return resp.json()

    tool = StructuredTool.from_function(predict_cfDNA)
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path
from typing import Optional, Union

# Ensure deepcatch is importable when run as module
_workspace = Path(__file__).resolve().parent.parent
if str(_workspace) not in sys.path:
    sys.path.insert(0, str(_workspace))

# ── FastAPI imports ──────────────────────────────────────────────
try:
    from fastapi import FastAPI, HTTPException
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.openapi.utils import get_openapi
except ImportError:
    raise ImportError(
        "fastapi is required. Install: pip install fastapi uvicorn"
    )

import numpy as np

from .schemas import (
    ClassicalResult,
    ExplanationResult,
    FileInput,
    HealthResponse,
    MotifAttribution,
    NeuralResult,
    PredictResponse,
    SequenceInput,
    TokenizationInfo,
)

# ── Model imports (lazy) ────────────────────────────────────────
_model: Optional["MotifDiversityModel"] = None
_tokenizer: Optional["DNATokenizer"] = None
_device: str = "cpu"
_torch_available: bool = False


def _init_model():
    """Initialize model and tokenizer (called once at startup)."""
    global _model, _tokenizer, _device, _torch_available

    from model.tokenizer import get_or_train_tokenizer, DNATokenizer
    from model.motif_model import (
        MotifDiversityModel,
        _get_device,
        _HAS_TORCH,
    )

    _torch_available = _HAS_TORCH
    _tokenizer = get_or_train_tokenizer()

    if _HAS_TORCH:
        import torch
        _device = _get_device()
        _model = MotifDiversityModel(
            vocab_size=_tokenizer.vocab_size,
            d_model=64,
            n_heads=4,
            n_layers=2,
            dropout=0.1,
        )

        # Try loading trained checkpoint
        ckpt_path = Path(__file__).resolve().parent.parent / "model" / "motif_model_checkpoint.pt"
        if ckpt_path.exists():
            try:
                ckpt = torch.load(ckpt_path, map_location=_device, weights_only=True)
                if ckpt.get("config", {}).get("vocab_size") == _tokenizer.vocab_size:
                    _model.load_state_dict(ckpt["model_state_dict"])
                    print(f"  Loaded trained checkpoint: {ckpt_path}")
                else:
                    print(f"  Checkpoint vocab mismatch, using untrained model")
            except Exception as e:
                print(f"  Checkpoint load failed: {e}, using untrained model")

        _model.to(_device)
        _model.eval()


def _read_fasta(path: str) -> str:
    """Read a single sequence from a FASTA file."""
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


# ── FastAPI App ─────────────────────────────────────────────────

app = FastAPI(
    title="DeepCatch Fragmentomics API",
    description=(
        "Local cancer detection from cfDNA fragment end motif patterns. "
        "Provides classical MDS, neural self-attention model, and "
        "Integrated Gradients explanations."
    ),
    version="2.1.0",
    docs_url="/docs",
    redoc_url="/redoc",
    openapi_url="/openapi.json",
)

# CORS: allow local agents and development tools
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
async def startup_event():
    """Load model and tokenizer into shared memory."""
    _init_model()


# ── Endpoints ───────────────────────────────────────────────────

@app.get("/health", response_model=HealthResponse)
async def health_check():
    """
    Health check.

    Returns model load status, compute device, vocab size,
    and PyTorch availability.
    """
    return HealthResponse(
        status="healthy",
        model_loaded=_model is not None,
        device=_device,
        vocab_size=_tokenizer.vocab_size if _tokenizer else 0,
        torch_available=_torch_available,
    )


@app.post("/predict", response_model=PredictResponse)
async def predict(
    request: Union[SequenceInput, FileInput],
):
    """
    Predict cancer from cfDNA fragment end motif frequencies.

    Accepts either a raw nucleotide sequence or a FASTA file path.
    The pipeline:

    1. Tokenize the DNA sequence with BPE subword tokenizer.
    2. Compute classical MDS (Simpson + Shannon diversity).
    3. Run the neural self-attention model.
    4. Optionally compute Integrated Gradients explanations.
    5. Return structured JSON with all results.

    **Example request (sequence):**

    .. code-block:: json

        {
            "sequence": "ACGTACGTNNACGTACGTACGT...",
            "method": "ensemble",
            "threshold": 0.5,
            "include_explanation": true
        }

    **Example request (file):**

    .. code-block:: json

        {
            "fasta_path": "/data/samples/sample_001.fasta",
            "method": "ensemble"
        }
    """
    t0 = time.perf_counter()

    # ── Determine input source ──
    if isinstance(request, SequenceInput):
        sequence = request.sequence.upper()
        method = request.method.value
        threshold = request.threshold
        include_explanation = request.include_explanation
        include_attention = request.include_attention
    else:
        sequence = _read_fasta(request.fasta_path)
        method = request.method.value
        threshold = request.threshold
        include_explanation = request.include_explanation
        include_attention = request.include_attention

    if not sequence:
        raise HTTPException(
            status_code=400, detail="Empty sequence after parsing."
        )

    # ── Tokenize ──
    if _tokenizer is None:
        raise HTTPException(
            status_code=503,
            detail="Tokenizer not initialized. Call /health first.",
        )

    frequencies = _tokenizer.count_frequencies(sequence)
    num_tokens = int(np.sum(frequencies * _tokenizer.vocab_size))

    tokenization = TokenizationInfo(
        vocab_size=_tokenizer.vocab_size,
        num_tokens=len(_tokenizer.encode(sequence)),
        sequence_length=len(sequence),
        method="bpe",
    )

    # ── Classical MDS ──
    classical: Optional[ClassicalResult] = None
    if method in ("classical", "ensemble"):
        from model.motif_model import compute_classical_mds

        # Reconstruct raw counts for MDS
        total_tokens = len(_tokenizer.encode(sequence))
        if total_tokens > 0:
            counts = (frequencies * total_tokens).astype(np.int64)
        else:
            counts = np.zeros(_tokenizer.vocab_size, dtype=np.int64)

        mds_simpson = compute_classical_mds(counts, method="simpson")
        mds_shannon = compute_classical_mds(counts, method="shannon")

        # Simple threshold for classical MDS (data-driven, default 0.5)
        classical_pred = "cancer" if mds_simpson > 0.5 else "healthy"

        classical = ClassicalResult(
            mds_simpson=round(mds_simpson, 6),
            mds_shannon=round(mds_shannon, 6),
            prediction=classical_pred,
        )

    # ── Neural Model ──
    neural: Optional[NeuralResult] = None
    explanation: Optional[ExplanationResult] = None
    ensemble_verdict: Optional[str] = None

    if method in ("neural", "ensemble"):
        if _model is None or not _torch_available:
            neural = NeuralResult(
                probability=0.0,
                prediction="unknown",
                threshold=threshold,
            )
        else:
            from model.motif_model import MotifPredictor

            predictor = MotifPredictor(
                _model, device=_device, threshold=threshold
            )
            result = predictor.predict(frequencies)

            neural = NeuralResult(
                probability=round(result["probability"], 6),
                prediction=result["prediction"],
                threshold=threshold,
            )

            # ── Integrated Gradients ──
            if include_explanation:
                from model.interpret import explain_prediction

                expl = explain_prediction(
                    frequencies,
                    _model,
                    tokenizer=_tokenizer,
                    top_k=10,
                    n_steps=50,
                )
                explanation = ExplanationResult(
                    top_motifs=[
                        MotifAttribution(**m) for m in expl["top_motifs"]
                    ],
                    bottom_motifs=[
                        MotifAttribution(**m) for m in expl["bottom_motifs"]
                    ],
                    convergence_delta=expl.get("convergence_delta"),
                )

    # ── Ensemble Verdict ──
    if method == "ensemble" and classical is not None and neural is not None:
        # Ensemble: neural model is primary, classical is tiebreaker
        if neural.probability >= threshold:
            ensemble_verdict = "cancer"
        elif neural.probability < 0.3 and classical.prediction == "healthy":
            ensemble_verdict = "healthy"
        else:
            # Uncertain zone: use log-odds average
            import math

            def _logit(p: float) -> float:
                eps = 1e-8
                return math.log(max(p, eps) / max(1 - p, eps))

            def _sigmoid(x: float) -> float:
                return 1.0 / (1.0 + math.exp(-x))

            neural_logit = _logit(neural.probability)
            classical_logit = (
                1.0 if classical.mds_simpson > 0.5 else -1.0
            )
            avg_logit = (neural_logit + classical_logit) / 2
            ensemble_prob = _sigmoid(avg_logit)
            ensemble_verdict = (
                "cancer" if ensemble_prob >= threshold else "healthy"
            )

    elapsed_ms = (time.perf_counter() - t0) * 1000

    return PredictResponse(
        status="ok",
        method=method,
        classical=classical,
        neural=neural,
        explanation=explanation,
        tokenization=tokenization,
        ensemble_verdict=ensemble_verdict,
        compute_time_ms=round(elapsed_ms, 2),
    )


@app.get("/")
async def root():
    """Redirect to docs."""
    return {
        "service": "DeepCatch Fragmentomics API",
        "version": "2.1.0",
        "docs": "/docs",
        "health": "/health",
        "predict": "/predict",
    }


# ── Entry Point ─────────────────────────────────────────────────

def main():
    """Entry point for ``python -m api.main``."""
    import uvicorn

    host = os.environ.get("DEEPCATCH_HOST", "0.0.0.0")
    port = int(os.environ.get("DEEPCATCH_PORT", "8000"))

    print(f"🦾 DeepCatch API starting on http://{host}:{port}")
    print(f"   Device:     {_device}")
    print(f"   Model:      {'loaded' if _model else 'not loaded'}")
    print(f"   Tokenizer:  vocab_size={_tokenizer.vocab_size if _tokenizer else 'N/A'}")
    print(f"   Docs:       http://{host}:{port}/docs")

    uvicorn.run(
        "api.main:app",
        host=host,
        port=port,
        reload=False,
        log_level="info",
    )


if __name__ == "__main__":
    main()
