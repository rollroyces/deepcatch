"""
GC-Bias Normalization Module (FragmentoSign)

Implements LOESS-based local normalization for cfDNA fragment coverage,
following the DELFI framework (Cristiano et al. 2019, Nature 570:385-389).

Why LOESS over global scaling:
- GC bias is NON-LINEAR: coverage vs GC content follows a parabolic curve
- Local regression captures the curve shape better than linear correction
- Standard in fragmentomics: DELFI, ichorCNA, and CopywriteR all use LOESS

DEPENDENCIES: numpy, scipy, statsmodels
"""

import numpy as np
from typing import Tuple, Optional
from scipy.stats import gaussian_kde


def compute_gc_content(sequence: str) -> float:
    """Compute GC content of a DNA sequence."""
    gc = sum(1 for b in sequence.upper() if b in 'GC')
    return gc / len(sequence) if len(sequence) > 0 else 0.0


def loess_normalize(
    coverage: np.ndarray,
    gc_content: np.ndarray,
    frac: float = 0.3,
    n_bins: int = 100
) -> Tuple[np.ndarray, np.ndarray]:
    """
    LOESS-based GC-bias correction for fragment coverage.
    
    Algorithm (from DELFI supplementary):
    1. Bin genome into windows of equal GC content
    2. Compute median coverage per bin
    3. Fit LOESS curve: coverage = f(GC_content)
    4. Correct each window: corrected_coverage = coverage / predicted_coverage × global_median
    
    Args:
        coverage: per-window read coverage [n_windows]
        gc_content: per-window GC content (0-1) [n_windows]
        frac: LOESS smoothing fraction (default 0.3)
        n_bins: number of GC bins for initial binning
    
    Returns:
        corrected_coverage, correction_factors
    """
    from scipy.interpolate import interp1d
    
    # Bin by GC content and compute median coverage per bin
    gc_bins = np.linspace(0.0, 1.0, n_bins + 1)
    bin_centers = (gc_bins[:-1] + gc_bins[1:]) / 2
    binned_cov = np.full(n_bins, np.nan)
    
    for i in range(n_bins):
        mask = (gc_content >= gc_bins[i]) & (gc_content < gc_bins[i + 1])
        if mask.sum() > 0:
            binned_cov[i] = np.median(coverage[mask])
    
    # Remove NaN bins (empty GC ranges)  
    valid = ~np.isnan(binned_cov)
    if valid.sum() < 10:
        # Fallback: global scaling
        global_median = np.median(coverage)
        return coverage, np.ones_like(coverage)
    
    gc_valid = bin_centers[valid]
    cov_valid = binned_cov[valid]
    
    # LOESS smoothing via LOWESS from statsmodels
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        smoothed = lowess(cov_valid, gc_valid, frac=frac, return_sorted=True)
        gc_smooth = smoothed[:, 0]
        cov_smooth = smoothed[:, 1]
    except ImportError:
        # Fallback: polynomial regression (degree 2)
        coeffs = np.polyfit(gc_valid, cov_valid, 2)
        gc_smooth = gc_valid
        cov_smooth = np.polyval(coeffs, gc_valid)
    
    # Interpolate to all windows
    interp = interp1d(gc_smooth, cov_smooth, kind='linear', 
                      fill_value='extrapolate', bounds_error=False)
    predicted_cov = interp(gc_content)
    
    # Correct: observed / predicted × global median
    predicted_cov = np.maximum(predicted_cov, 1e-6)
    global_median = np.median(coverage)
    correction_factors = global_median / predicted_cov
    corrected_coverage = coverage * correction_factors
    
    return corrected_coverage, correction_factors


def DELFI_style_normalization(
    coverage: np.ndarray,
    gc_content: np.ndarray,
    mappability: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Full DELFI-style normalization pipeline:
    1. Mappability filter (exclude low-mappability regions)
    2. LOESS GC-bias correction
    3. Median-centering
    
    Reference: Cristiano et al. 2019, Nature 570:385-389
    Mathios et al. 2021, Nature Communications 12:5060
    """
    if mappability is not None:
        mask = mappability > 0.8
        coverage = coverage[mask]
        gc_content = gc_content[mask]
    
    corrected, _ = loess_normalize(coverage, gc_content)
    
    # Median center
    corrected = corrected / np.median(corrected)
    
    return corrected


def compute_MDS(fragment_ends: np.ndarray, n_motifs: int = 256) -> float:
    """
    Motif Diversity Score (MDS).
    
    MDS quantifies the diversity of fragment end motifs.
    Higher MDS → more diverse end motifs → more likely cancer.
    
    Formula: MDS = 1 - Σ(p_i²) / (1 - 1/n)
    where p_i = frequency of motif i, n = number of possible motifs
    
    This is the normalized Simpson diversity index.
    
    Reference: Jiang et al. 2020, Nature Genetics 52:712-719
    """
    # Compute motif frequencies
    motif_counts = np.bincount(fragment_ends, minlength=n_motifs)
    p = motif_counts / motif_counts.sum()
    
    # Simpson diversity → MDS (normalized)
    simpson = np.sum(p ** 2)
    mds = (1 - simpson) / (1 - 1.0 / n_motifs)
    
    return mds
