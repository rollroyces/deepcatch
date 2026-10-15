"""Real implementation for diffusion_model.py.

Original location: src/agent5-synthetic-data/diffusion_model.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""Lightweight forward/reverse diffusion on (size, motif) embeddings.

This is *not* a full DDPM — it is a Gaussian forward noising + linear-score
reverse process used for synthetic cohort enrichment (see
METHYLATION_PROJECT.md §5 for the diffusion extension plan).
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class DiffusionConfig:
    n_steps: int = 100
    beta_start: float = 1e-4
    beta_end: float = 0.02


class GaussianDiffusion1D:
    def __init__(self, cfg: DiffusionConfig | None = None):
        self.cfg = cfg or DiffusionConfig()
        self.betas = np.linspace(self.cfg.beta_start, self.cfg.beta_end,
                                 self.cfg.n_steps)
        self.alphas = 1.0 - self.betas
        self.alpha_bars = np.cumprod(self.alphas)

    def forward(self, x0: np.ndarray, t: int,
                noise: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
        rng = np.random.default_rng()
        if noise is None:
            noise = rng.normal(0, 1, size=x0.shape)
        a = self.alpha_bars[t]
        return np.sqrt(a) * x0 + np.sqrt(1 - a) * noise, noise

    def reverse_step(self, x_t: np.ndarray, t: int,
                     predicted_noise: np.ndarray) -> np.ndarray:
        a = self.alphas[t]
        ab = self.alpha_bars[t]
        mean = (x_t - self.betas[t] / np.sqrt(1 - ab) * predicted_noise) / np.sqrt(a)
        if t > 0:
            mean += np.sqrt(self.betas[t]) * np.random.normal(0, 1, size=x_t.shape)
        return mean


__all__ = ["DiffusionConfig", "GaussianDiffusion1D"]
