#!/usr/bin/env python3
"""
Lightweight Fragmentomics Neural Model
=======================================

A compact neural network that predicts cancer from cfDNA fragment end
motif frequency histograms.  Designed to run comfortably on Apple Silicon
with 16 GB of unified memory.

.. rubric:: Architecture

::

    Frequency Histogram (vocab_size,)
         │
         ▼
    Embedding (vocab_size → d_model)
         │
         ▼
    Positional Encoding (learned, per-motif)
         │
         ▼
    ┌─────────────────────────────────┐
    │   Multi-Head Self-Attention     │  ← motifs attend to each other
    │   (2 heads, d_k=d_model//2)     │     to learn co-occurrence
    ├─────────────────────────────────┤
    │   LayerNorm + FFN (d_model→4×)  │
    └──────────────┬──────────────────┘
         │
         ▼
    Global Average Pooling
         │
         ▼
    MLP Head:  d_model → 32 → 1  (sigmoid)
         │
         ▼
    Cancer Probability ∈ [0, 1]

.. rubric:: Memory Budget (vocab_size=256, d_model=64)

- Embedding:      256 × 64  = 16,384  floats  (65.5 KB)
- Pos encoding:   256 × 64  = 16,384  floats  (65.5 KB)
- Q/K/V/O proj:   4 × 64²  = 16,384  floats  (65.5 KB)
- FFN:            64×256 + 256×64 = 32,768 floats (131 KB)
- MLP head:       64×32 + 32×1 = 2,080 floats   (8.3 KB)
- **Total**: ~86K params (~344 KB FP32)

Even with vocab_size=512 and d_model=128: ~350K params (~1.4 MB).

.. rubric:: Device Selection

On Apple Silicon, :func:`_get_device` returns ``mps`` when PyTorch ≥1.12
and ``torch.backends.mps.is_available()`` is True.  Falls back to CPU
otherwise (no penalty—the model is tiny).

.. rubric:: Baseline Comparison

The original MDS (normalized Simpson diversity) is computed alongside
the model output for every prediction.  The model can learn interactions
that pure statistical diversity cannot capture — e.g., the diagnostic
significance of a *specific* motif being depleted while another is enriched.
"""

from __future__ import annotations

import math
import warnings
from typing import Dict, Optional, Tuple

import numpy as np

# ── PyTorch imports (lazy, with informative errors) ──────────────
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F

    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False
    torch = None  # type: ignore[assignment]
    nn = None      # type: ignore[assignment]
    F = None       # type: ignore[assignment]


def _get_device() -> str:
    """Return the best available torch device ('mps', 'cuda', or 'cpu')."""
    if not _HAS_TORCH:
        return "cpu"
    # Use a try/except for MPS (may not be compiled)
    try:
        if torch.backends.mps.is_available():
            return "mps"
    except Exception:
        pass
    if torch.cuda.is_available():
        return "cuda"
    return "cpu"


# ── Classical MDS (kept as baseline) ────────────────────────────

def compute_classical_mds(
    counts: np.ndarray, method: str = "simpson"
) -> float:
    """
    Compute the classical Motif Diversity Score from motif counts.

    Parameters
    ----------
    counts : np.ndarray
        1-D array of motif counts.
    method : str
        ``"simpson"`` — Normalized Simpson diversity (original DeepCatch).
        ``"shannon"`` — Normalized Shannon entropy.

    Returns
    -------
    float
        MDS in [0, 1].

    Notes
    -----
    Shannon MDS:

    .. math::

        \\text{MDS}_\\text{Shannon} =
        -\\frac{1}{\\log_2(n)} \\sum_{i=1}^{n} p_i \\log_2(p_i)

    Simpson MDS:

    .. math::

        \\text{MDS}_\\text{Simpson} =
        \\frac{1 - \\sum p_i^2}{1 - 1/n}
    """
    n = len(counts)
    total = counts.sum()
    if total == 0:
        return 0.0

    p = counts / total
    p = p[p > 0]  # avoid log(0)

    if method == "shannon":
        entropy = -np.sum(p * np.log2(p))
        return float(entropy / np.log2(n))
    else:  # simpson (default)
        simpson = np.sum((counts / total) ** 2)
        return float((1 - simpson) / (1 - 1.0 / n))


# ── PyTorch Model Classes (conditionally defined) ────────────────

_STUB_CLASSES: Dict[str, type] = {}

if _HAS_TORCH:
    # ------------------------------------------------------------------
    # Real PyTorch implementations
    # ------------------------------------------------------------------

    class MotifSelfAttention(nn.Module):
        """
        Multi-head self-attention over motif embeddings.

        Unlike standard NLP attention (which operates over sequence
        positions), this attention operates over the *vocabulary axis*.
        Each motif is a "token," and attention weights capture which
        motifs co-occur or mutually exclude each other in cfDNA samples.

        Parameters
        ----------
        d_model : int
            Embedding dimension.
        n_heads : int
            Number of attention heads (default 4).
        dropout : float
            Dropout rate (default 0.1).
        """

        def __init__(
            self, d_model: int, n_heads: int = 4, dropout: float = 0.1
        ):
            super().__init__()
            assert d_model % n_heads == 0, "d_model must be divisible by n_heads"
            self.d_model = d_model
            self.n_heads = n_heads
            self.d_k = d_model // n_heads
            self.q_proj = nn.Linear(d_model, d_model)
            self.k_proj = nn.Linear(d_model, d_model)
            self.v_proj = nn.Linear(d_model, d_model)
            self.out_proj = nn.Linear(d_model, d_model)
            self.dropout = nn.Dropout(dropout)

        def forward(self, x):
            B, N, D = x.shape
            q = self.q_proj(x).view(B, N, self.n_heads, self.d_k).transpose(1, 2)
            k = self.k_proj(x).view(B, N, self.n_heads, self.d_k).transpose(1, 2)
            v = self.v_proj(x).view(B, N, self.n_heads, self.d_k).transpose(1, 2)
            scale = math.sqrt(self.d_k)
            attn_weights = torch.matmul(q, k.transpose(-2, -1)) / scale
            attn_weights = F.softmax(attn_weights, dim=-1)
            attn_weights = self.dropout(attn_weights)
            attn_out = torch.matmul(attn_weights, v)
            attn_out = attn_out.transpose(1, 2).contiguous().view(B, N, D)
            return self.out_proj(attn_out)


    class MotifDiversityModel(nn.Module):
        """
        Lightweight neural motif diversity model for cfDNA cancer detection.

        Accepts a frequency histogram over BPE token IDs and predicts
        cancer probability.  Self-attention over the vocabulary axis
        captures dependencies between motifs that classical MDS cannot
        express.

        Parameters
        ----------
        vocab_size : int
            Number of tokens in vocabulary (default 256).
        d_model : int
            Internal embedding dimension (default 64).
        n_heads : int
            Attention heads (default 4).
        n_layers : int
            Number of attention + FFN blocks (default 2).
        dropout : float
            Dropout rate (default 0.1).
        use_attention : bool
            If False, model is a plain MLP — useful as baseline.
        """

        def __init__(
            self,
            vocab_size: int = 256,
            d_model: int = 64,
            n_heads: int = 4,
            n_layers: int = 2,
            dropout: float = 0.1,
            use_attention: bool = True,
        ):
            super().__init__()
            self.vocab_size = vocab_size
            self.d_model = d_model
            self.use_attention = use_attention

            self.freq_embed = nn.Sequential(
                nn.Linear(1, d_model),
                nn.GELU(),
            )
            self.pos_embed = nn.Parameter(
                torch.randn(1, vocab_size, d_model) * 0.02
            )

            self.layers = nn.ModuleList()
            for _ in range(n_layers):
                self.layers.append(
                    nn.ModuleDict(
                        {
                            "attn": (
                                MotifSelfAttention(d_model, n_heads, dropout)
                                if use_attention
                                else nn.Identity()
                            ),
                            "norm1": nn.LayerNorm(d_model),
                            "ffn": nn.Sequential(
                                nn.Linear(d_model, d_model * 4),
                                nn.GELU(),
                                nn.Dropout(dropout),
                                nn.Linear(d_model * 4, d_model),
                                nn.Dropout(dropout),
                            ),
                            "norm2": nn.LayerNorm(d_model),
                        }
                    )
                )

            self.head = nn.Sequential(
                nn.Linear(d_model, 32),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(32, 1),
                nn.Sigmoid(),
            )
            self._init_weights()

        def _init_weights(self):
            for m in self.modules():
                if isinstance(m, nn.Linear):
                    nn.init.xavier_uniform_(m.weight, gain=0.5)
                    if m.bias is not None:
                        nn.init.zeros_(m.bias)

        def forward(self, frequencies, return_attention=False):
            B, V = frequencies.shape
            assert V == self.vocab_size, (
                f"Expected vocab_size={self.vocab_size}, got {V}"
            )
            x = frequencies.unsqueeze(-1)          # (B, V, 1)
            x = self.freq_embed(x)                 # (B, V, d_model)
            x = x + self.pos_embed                 # (B, V, d_model)

            last_attn = None
            for layer in self.layers:
                if self.use_attention:
                    attn_out = layer["attn"](x)
                    if return_attention:
                        last_attn = self._capture_attention(layer["attn"], x)
                    x = layer["norm1"](x + attn_out)
                ffn_out = layer["ffn"](x)
                x = layer["norm2"](x + ffn_out)

            x = x.mean(dim=1)                      # (B, d_model)
            probs = self.head(x)                   # (B, 1)

            if return_attention:
                return probs, last_attn
            return probs

        @staticmethod
        def _capture_attention(attn_module, x):
            B, N, D = x.shape
            q = attn_module.q_proj(x).view(
                B, N, attn_module.n_heads, attn_module.d_k
            ).transpose(1, 2)
            k = attn_module.k_proj(x).view(
                B, N, attn_module.n_heads, attn_module.d_k
            ).transpose(1, 2)
            scale = math.sqrt(attn_module.d_k)
            weights = torch.matmul(q, k.transpose(-2, -1)) / scale
            return F.softmax(weights, dim=-1)

        def predict_proba(self, frequencies):
            self.eval()
            with torch.no_grad():
                return self.forward(frequencies)

        @property
        def num_parameters(self) -> int:
            return sum(
                p.numel() for p in self.parameters() if p.requires_grad
            )


    _STUB_CLASSES["MotifSelfAttention"] = MotifSelfAttention
    _STUB_CLASSES["MotifDiversityModel"] = MotifDiversityModel

else:
    # ------------------------------------------------------------------
    # Stub classes when PyTorch is not installed
    # (classical MDS + tokenizer still work)
    # ------------------------------------------------------------------

    class _TorchStub:
        """Base stub for torch-dependent classes."""

        def __init__(self, *args, **kwargs):
            self.vocab_size = kwargs.get("vocab_size", 256)
            self.d_model = kwargs.get("d_model", 64)
            warnings.warn(
                "PyTorch is not installed. Neural model is unavailable. "
                "Classical MDS fallback works.  Install: pip install torch"
            )

        def __call__(self, *args, **kwargs):
            raise RuntimeError(
                "PyTorch not installed — neural model unavailable. "
                "Use classical MDS via compute_classical_mds()."
            )

        def to(self, *args, **kwargs):
            return self

        def eval(self):
            return self

        def train(self, mode=True):
            return self

        def parameters(self):
            return []

        def modules(self):
            return []

        def state_dict(self):
            return {}

        def load_state_dict(self, state_dict, strict=True):
            pass

        @property
        def num_parameters(self) -> int:
            return 0

        def predict_proba(self, frequencies):
            raise RuntimeError("PyTorch not installed — model unavailable.")


    class MotifSelfAttention(_TorchStub):
        """Stub for MotifSelfAttention when PyTorch is absent."""
        pass


    class MotifDiversityModel(_TorchStub):
        """Stub for MotifDiversityModel when PyTorch is absent."""
        pass

    _STUB_CLASSES["MotifSelfAttention"] = MotifSelfAttention
    _STUB_CLASSES["MotifDiversityModel"] = MotifDiversityModel


# ── High-level Predictor Class ──────────────────────────────────

class MotifPredictor:
    """
    End-to-end predictor wrapping the BPE tokenizer, neural model,
    and classical MDS baseline.

    Parameters
    ----------
    model : MotifDiversityModel
        Trained PyTorch model.
    device : str
        Torch device string (``"mps"``, ``"cuda"``, ``"cpu"``).
    threshold : float
        Decision threshold for binary prediction (default 0.5).
    """

    def __init__(
        self,
        model: "MotifDiversityModel",
        device: str = "cpu",
        threshold: float = 0.5,
    ):
        self.model = model.to(device) if hasattr(model, "to") else model
        self.device = device
        self.threshold = threshold
        if hasattr(self.model, "eval"):
            self.model.eval()

    def predict(
        self, frequencies: np.ndarray
    ) -> Dict[str, object]:
        """
        Run a full prediction on a frequency histogram.

        Parameters
        ----------
        frequencies : np.ndarray of shape ``(vocab_size,)``
            Normalized motif frequencies.

        Returns
        -------
        dict with keys:
            - ``probability``: cancer probability (float)
            - ``prediction``: ``"cancer"`` or ``"healthy"``
            - ``classical_mds``: MDS from statistical method (float)
            - ``threshold``: decision threshold used
        """
        classical_mds = compute_classical_mds(
            frequencies, method="simpson",
        )

        if not _HAS_TORCH:
            return {
                "probability": 0.5,
                "prediction": (
                    "cancer" if classical_mds > self.threshold else "healthy"
                ),
                "classical_mds": classical_mds,
                "threshold": self.threshold,
            }

        freq_tensor = torch.from_numpy(frequencies).float().unsqueeze(0)
        freq_tensor = freq_tensor.to(self.device)

        with torch.no_grad():
            prob = self.model(freq_tensor)
            prob = float(prob.cpu().flatten()[0])

        return {
            "probability": prob,
            "prediction": "cancer" if prob >= self.threshold else "healthy",
            "classical_mds": classical_mds,
            "threshold": self.threshold,
        }

    @classmethod
    def from_checkpoint(
        cls,
        checkpoint_path: str,
        device: Optional[str] = None,
        **kwargs,
    ) -> "MotifPredictor":
        """
        Load a predictor from a saved model checkpoint.

        Parameters
        ----------
        checkpoint_path : str
            Path to ``.pt`` or ``.pth`` file.
        device : str, optional
            Device override.

        Returns
        -------
        MotifPredictor
        """
        if not _HAS_TORCH:
            raise ImportError(
                "PyTorch is required to load model checkpoints. "
                "Install: pip install torch"
            )

        if device is None:
            device = _get_device()

        checkpoint = torch.load(
            checkpoint_path, map_location=device, weights_only=True
        )
        model = MotifDiversityModel(**checkpoint.get("config", {}))
        model.load_state_dict(checkpoint["model_state_dict"])
        return cls(model, device=device, **kwargs)


__all__ = [
    "MotifDiversityModel",
    "MotifSelfAttention",
    "MotifPredictor",
    "compute_classical_mds",
    "_get_device",
    "_HAS_TORCH",
]
