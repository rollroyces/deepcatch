"""
DeepCatch Model — Lightweight neural fragmentomics for Apple Silicon.

Provides a memory-efficient BPE tokenizer, attention-based motif model,
and integrated-gradients interpretability for cfDNA cancer detection.

.. rubric:: Modules

- ``tokenizer`` — DNA-optimized BPE subword tokenizer
- ``motif_model`` — Lightweight CNN + self-attention model
- ``interpret`` — Integrated Gradients feature attribution

Designed for Apple Silicon (M1–M3, 16 GB unified memory).
Uses MPS-accelerated PyTorch when available; pure NumPy fallback for inference.
"""

# Tokenizer always works (pure Python)
from .tokenizer import (
    DNATokenizer,
    train_bpe_tokenizer,
    get_or_train_tokenizer,
    tokenize_sequence,
    DEFAULT_VOCAB_SIZE,
    MAX_VOCAB_SIZE,
)

# Model — works with or without PyTorch
from .motif_model import (
    MotifDiversityModel,
    MotifSelfAttention,
    MotifPredictor,
    compute_classical_mds,
    _get_device,
    _HAS_TORCH,
)

# Integrated Gradients — works only with PyTorch, graceful fallback
from .interpret import (
    IntegratedGradientsExplainer,
    explain_prediction,
)

__all__ = [
    # Tokenizer
    "DNATokenizer",
    "train_bpe_tokenizer",
    "get_or_train_tokenizer",
    "tokenize_sequence",
    "DEFAULT_VOCAB_SIZE",
    "MAX_VOCAB_SIZE",
    # Model
    "MotifDiversityModel",
    "MotifSelfAttention",
    "MotifPredictor",
    "compute_classical_mds",
    "_get_device",
    "_HAS_TORCH",
    # Interpret
    "IntegratedGradientsExplainer",
    "explain_prediction",
]
