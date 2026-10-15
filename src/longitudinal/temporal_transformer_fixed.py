"""Real implementation for temporal_transformer_fixed.py.

Original location: src/agent3-longitudinal/temporal_transformer_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version of :mod:`temporal_transformer` — uses scaled dot-product
with a numerical stability fix for very long sequences.
"""
from .temporal_transformer import TemporalSurrogate
import numpy as np


class FixedTemporalSurrogate(TemporalSurrogate):
    def forward(self, X: np.ndarray) -> np.ndarray:
        Q = X @ self.W_q
        K = X @ self.W_k
        V = X @ self.W_v
        scale = np.sqrt(self.n_features) + 1e-6
        scores = Q @ K.T / scale
        # Numerical stability: subtract per-row max before softmax.
        scores -= scores.max(axis=-1, keepdims=True)
        weights = np.exp(scores)
        weights /= weights.sum(axis=-1, keepdims=True) + 1e-12
        return weights @ V


__all__ = ["TemporalSurrogate", "FixedTemporalSurrogate"]
