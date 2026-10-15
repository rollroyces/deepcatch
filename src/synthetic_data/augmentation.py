"""Real implementation for augmentation.py.

Original location: src/agent5-synthetic-data/augmentation.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""Data augmentation for cfDNA fragmentomics samples.

Implements GC-bias shuffle, motif jitter, and multiplicative size noise.
"""
import numpy as np


def gc_bias(sizes: np.ndarray, gc_target: float = 0.5, strength: float = 0.1,
            seed: int | None = None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    bias = 1.0 + strength * (rng.normal(gc_target, 0.05, size=len(sizes)) - gc_target)
    return (sizes * bias).astype(int).clip(50, 500)


def motif_jitter(motifs: np.ndarray, p: float = 0.01,
                 seed: int | None = None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = motifs.copy()
    flip = rng.random(out.shape) < p
    perm = rng.integers(0, out.shape[1], size=out.shape)
    out[flip] = out[flip, perm[flip]]
    return out


def size_noise(sizes: np.ndarray, sigma_bp: float = 2.0,
               seed: int | None = None) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return (sizes + rng.normal(0, sigma_bp, size=len(sizes))).astype(int).clip(50, 500)


__all__ = ["gc_bias", "motif_jitter", "size_noise"]
