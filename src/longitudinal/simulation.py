"""Real implementation for simulation.py.

Original location: src/agent3-longitudinal/simulation.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent3-longitudinal/` directory that was never committed to git.
"""

from __future__ import annotations

"""Longitudinal cfDNA trajectory simulator.

Produces serial measurements of a marker (e.g., VAF) for one patient over
``n_draws`` time points, with an optional changepoint that simulates
treatment response or progression.
"""
import numpy as np


def simulate_trajectory(n_draws: int = 8, changepoint_chat: int | None = 4,
                        baseline: float = 0.005, treatment_shift: float = -0.003,
                        noise_sigma: float = 0.0015,
                        seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = np.zeros(n_draws)
    for t in range(n_draws):
        level = baseline
        if changepoint_chat is not None and t >= changepoint_chat:
            level = max(0.0, baseline + treatment_shift)
        out[t] = max(0.0, level + rng.normal(0, noise_sigma))
    return out


__all__ = ["simulate_trajectory"]
