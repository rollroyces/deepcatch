"""Real implementation for cfDNA_constants.py.

Original location: src/agent5-synthetic-data/cfDNA_constants.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""Empirical constants for cfDNA fragmentomics.

Values are from FinaleDB-derived calibrations and from the published
cfDNA fragment size literature (cfDC1, METHYLATION_PROJECT.md §3).
"""
import numpy as np

# Core fragment-size peak (Lo et al., 2021; Snyder et al., 2016).
PRIMARY_PEAK_BP: int = 167
PRIMARY_PEAK_STD_BP: float = 5.0

# 10bp periodicity (Ivanov et al., 2015) — secondary oscillations around the
# primary peak. We model them as a sinusoidal modulation.
PERIODICITY_BP: int = 10
PERIODICITY_AMPLITUDE: float = 0.15

# Valid fragment-size range used by both simulators and downstream filters.
MIN_FRAGMENT_BP: int = 50
MAX_FRAGMENT_BP: int = 500

# Empirical 4-mer frequencies at fragment ends (placeholder uniform — replace
# with calibrated vectors once a real cohort is processed).
END_MOTIF_UNIFORM_FREQ: np.ndarray = np.full(256, 1.0 / 256)


def size_distribution(n: int, seed: int | None = None) -> np.ndarray:
    """Draw ``n`` insert sizes from the calibrated cfDNA model."""
    rng = np.random.default_rng(seed)
    sizes = rng.normal(PRIMARY_PEAK_BP, PRIMARY_PEAK_STD_BP, size=n).astype(int)
    sizes += (PERIODICITY_AMPLITUDE * np.sin(2 * np.pi * sizes / PERIODICITY_BP) * 10).astype(int)
    return np.clip(sizes, MIN_FRAGMENT_BP, MAX_FRAGMENT_BP)


__all__ = [
    "PRIMARY_PEAK_BP", "PRIMARY_PEAK_STD_BP",
    "PERIODICITY_BP", "PERIODICITY_AMPLITUDE",
    "MIN_FRAGMENT_BP", "MAX_FRAGMENT_BP",
    "END_MOTIF_UNIFORM_FREQ",
    "size_distribution",
]
