"""Real implementation for temporal_transformer.py.

Original location: src/agent3-longitudinal/temporal_transformer.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Lightweight temporal model — a linear-attention surrogate for the
Transformer described in README §"Foundation model + longitudinal Bayes
Kalman".

This implementation uses numpy-only positional + value mixing, suitable
as a CPU-friendly smoke-test before the production Transformer checkpoint
is loaded.
"""
import numpy as np


class TemporalSurrogate:
    def __init__(self, n_features: int, n_heads: int = 2):
        self.n_features = n_features
        self.n_heads = n_heads
        rng = np.random.default_rng(0)
        self.W_q = rng.normal(0, 0.1, size=(n_features, n_features))
        self.W_k = rng.normal(0, 0.1, size=(n_features, n_features))
        self.W_v = rng.normal(0, 0.1, size=(n_features, n_features))

    def forward(self, X: np.ndarray) -> np.ndarray:
        # X shape: (T, F)
        Q = X @ self.W_q
        K = X @ self.W_k
        V = X @ self.W_v
        scores = Q @ K.T / np.sqrt(self.n_features)
        weights = np.exp(scores - scores.max(axis=-1, keepdims=True))
        weights /= weights.sum(axis=-1, keepdims=True) + 1e-12
        return weights @ V


__all__ = ["TemporalSurrogate"]
