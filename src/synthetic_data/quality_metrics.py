"""Real implementation for quality_metrics.py.

Original location: src/agent5-synthetic-data/quality_metrics.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent5-synthetic-data/` directory that was never committed to git.
"""

from __future__ import annotations

"""Quality metrics for synthetic cfDNA cohorts.

Used by the smoke validation pipeline (PHASE0_PLAN.md) to confirm the
synthetic generator matches empirical cfDNA fragment-size moments before
training downstream classifiers.
"""
import numpy as np

from .cfDNA_constants import PRIMARY_PEAK_BP, PRIMARY_PEAK_STD_BP


def fragment_size_moments(sizes: np.ndarray) -> dict[str, float]:
    mean = float(sizes.mean())
    std = float(sizes.std())
    skew = float(((sizes - mean) ** 3).mean() / (std ** 3 + 1e-12))
    return {
        "mean_bp": mean,
        "std_bp": std,
        "skewness": skew,
        "delta_mean_from_167": mean - PRIMARY_PEAK_BP,
    }


def fft_periodicity_score(sizes: np.ndarray, target_period_bp: int = 10) -> float:
    """Peak FFT magnitude at 1/``target_period_bp`` Hz relative to mean."""
    counts, _ = np.histogram(sizes, bins=range(50, 500))
    yf = np.abs(np.fft.rfft(counts - counts.mean()))
    freqs = np.fft.rfftfreq(len(counts), 1.0)
    target_hz = 1.0 / target_period_bp
    idx = int(np.argmin(np.abs(freqs - target_hz)))
    return float(yf[idx] / (yf.mean() + 1e-12))


__all__ = ["fragment_size_moments", "fft_periodicity_score"]
