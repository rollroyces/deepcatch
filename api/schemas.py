"""
Pydantic schemas for the DeepCatch Fragmentomics API.

All request/response models live here to keep ``main.py`` lean.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import List, Optional

from pydantic import BaseModel, Field, field_validator


# ── Request Models ──────────────────────────────────────────────

class PredictionMethod(str, Enum):
    """Available prediction methods."""

    CLASSICAL = "classical"     # Statistical MDS only
    NEURAL = "neural"           # Neural model only
    ENSEMBLE = "ensemble"       # Both (default)


class SequenceInput(BaseModel):
    """
    Input for /predict when the sequence is provided directly.

    Example
    -------
    .. code-block:: json

        {
            "sequence": "ACGTACGTNNACGT...",
            "method": "ensemble",
            "threshold": 0.5,
            "include_explanation": true,
            "token_vocab_size": 256
        }
    """

    sequence: str = Field(
        ...,
        min_length=4,
        description="Raw nucleotide sequence (A, C, G, T, N).",
        examples=["ACGTACGTACGTNNNNACGTACGT"],
    )
    method: PredictionMethod = Field(
        PredictionMethod.ENSEMBLE,
        description="Prediction method.",
    )
    threshold: float = Field(
        0.5,
        ge=0.0,
        le=1.0,
        description="Decision threshold for binary prediction.",
    )
    include_explanation: bool = Field(
        True,
        description="Include Integrated Gradients feature attribution.",
    )
    include_attention: bool = Field(
        False,
        description="Include self-attention weight matrix.",
    )
    token_vocab_size: int = Field(
        256,
        ge=32,
        le=512,
        description="BPE vocabulary size for tokenization.",
    )

    @field_validator("sequence")
    @classmethod
    def sequence_must_be_dna(cls, v: str) -> str:
        v = v.upper()
        invalid = set(v) - {"A", "C", "G", "T", "N"}
        if invalid:
            raise ValueError(
                f"Sequence contains invalid characters: {invalid}. "
                f"Only A, C, G, T, N are allowed."
            )
        return v


class FileInput(BaseModel):
    """
    Input for /predict when a FASTA file path is provided.

    Example
    -------
    .. code-block:: json

        {
            "fasta_path": "/data/sample_001.fasta",
            "method": "ensemble",
            "threshold": 0.5,
            "include_explanation": true
        }
    """

    fasta_path: str = Field(
        ...,
        description="Path to a .fasta or .fa file.",
    )
    method: PredictionMethod = Field(PredictionMethod.ENSEMBLE)
    threshold: float = Field(0.5, ge=0.0, le=1.0)
    include_explanation: bool = Field(True)
    include_attention: bool = Field(False)
    token_vocab_size: int = Field(256, ge=32, le=512)

    @field_validator("fasta_path")
    @classmethod
    def validate_fasta_extension(cls, v: str) -> str:
        p = Path(v)
        if p.suffix.lower() not in (".fasta", ".fa", ".fna", ".ffn"):
            raise ValueError(
                f"Expected a FASTA file extension (.fasta, .fa, .fna, .ffn), "
                f"got: {p.suffix}"
            )
        if not p.exists():
            raise FileNotFoundError(f"FASTA file not found: {v}")
        return str(p.resolve())


# ── Response Models ──────────────────────────────────────────────

class MotifAttribution(BaseModel):
    """A single motif attribution entry."""

    token: str = Field(..., description="DNA subword token string.")
    id: int = Field(..., description="Token ID in vocabulary.")
    attribution: float = Field(
        ..., description="Integrated Gradients attribution score."
    )


class ClassicalResult(BaseModel):
    """Classical MDS baseline result."""

    mds_simpson: float = Field(
        ..., description="Normalized Simpson diversity MDS ∈ [0, 1]."
    )
    mds_shannon: float = Field(
        ..., description="Normalized Shannon entropy MDS ∈ [0, 1]."
    )
    prediction: str = Field(
        ..., description="Binary prediction: 'cancer' or 'healthy'."
    )


class NeuralResult(BaseModel):
    """Neural model prediction result."""

    probability: float = Field(
        ..., description="Cancer probability ∈ [0, 1]."
    )
    prediction: str = Field(
        ..., description="Binary prediction: 'cancer' or 'healthy'."
    )
    threshold: float = Field(
        ..., description="Decision threshold used."
    )


class ExplanationResult(BaseModel):
    """Integrated Gradients explanation."""

    top_motifs: List[MotifAttribution] = Field(
        default_factory=list,
        description="Motifs most strongly increasing cancer probability.",
    )
    bottom_motifs: List[MotifAttribution] = Field(
        default_factory=list,
        description="Motifs most strongly decreasing cancer probability.",
    )
    convergence_delta: Optional[float] = Field(
        None,
        description="IG approximation error (lower is better).",
    )


class TokenizationInfo(BaseModel):
    """Information about tokenization."""

    vocab_size: int
    num_tokens: int
    sequence_length: int
    method: str = "bpe"


class PredictResponse(BaseModel):
    """
    Response schema for the /predict endpoint.

    Example
    -------
    .. code-block:: json

        {
            "status": "ok",
            "method": "ensemble",
            "classical": {
                "mds_simpson": 0.723,
                "mds_shannon": 0.691,
                "prediction": "cancer"
            },
            "neural": {
                "probability": 0.847,
                "prediction": "cancer",
                "threshold": 0.5
            },
            "explanation": {
                "top_motifs": [
                    {"token": "CG", "id": 12, "attribution": 0.034}
                ],
                "bottom_motifs": [],
                "convergence_delta": 0.0012
            },
            "tokenization": {
                "vocab_size": 256,
                "num_tokens": 15320,
                "sequence_length": 98304,
                "method": "bpe"
            },
            "ensemble_verdict": "cancer",
            "compute_time_ms": 42.7
        }
    """

    status: str = Field("ok", description="Request status.")
    method: str = Field(..., description="Prediction method used.")
    classical: Optional[ClassicalResult] = None
    neural: Optional[NeuralResult] = None
    explanation: Optional[ExplanationResult] = None
    tokenization: Optional[TokenizationInfo] = None
    ensemble_verdict: Optional[str] = Field(
        None, description="Final ensemble verdict."
    )
    compute_time_ms: float = Field(
        ..., description="Total compute time in milliseconds."
    )


class HealthResponse(BaseModel):
    """Health-check response."""

    status: str = "healthy"
    model_loaded: bool = Field(..., description="Whether model is loaded.")
    device: str = Field(..., description="Compute device (mps/cpu).")
    vocab_size: int = Field(0, description="Tokenizer vocabulary size.")
    torch_available: bool = Field(..., description="PyTorch availability.")
