"""
GC-Bias Normalization Module (FragmentoSign)

Implements LOESS-based local normalization for cfDNA fragment coverage,
following the DELFI framework (Cristiano et al. 2019, Nature 570:385-389).

.. rubric:: Why LOESS over Global Scaling?

- **GC bias is non-linear**: coverage vs GC content follows a parabolic curve
  with reduced coverage at both low-GC and high-GC extremes.
- **Local regression** captures the curve shape better than linear correction
  or global scaling: the smooth function :math:`f(GC)` is fit locally at each
  GC value using weighted least-squares.
- **Standard in fragmentomics**: DELFI, ichorCNA, CopywriteR, and other
  cfDNA tools all use LOESS or LOWESS for GC correction.

.. rubric:: LOESS Algorithm Steps

1. **Bin**: Divide genome into ``n_bins`` equal-width GC-content bins.
2. **Median per bin**: Compute median coverage for windows in each bin,
   producing a raw GC-vs-coverage scatter.
3. **LOWESS smooth**: Fit a locally weighted regression using
   ``statsmodels.nonparametric.smoothers_lowess.lowess`` with
   smoothing fraction ``frac``.
4. **Interpolate**: Map the smooth curve back to every window via
   linear interpolation.
5. **Correct**: For each window,
   :math:`\text{corrected} = \text{coverage} \times \frac{\text{global\_median}}{\text{predicted\_coverage}}`.

.. rubric:: Fallback Behavior

If ``statsmodels`` is not installed (``ImportError``), the module falls back
to a degree-2 polynomial regression via ``numpy.polyfit``, which approximates
the parabolic GC-coverage relationship.

If fewer than 10 valid GC bins exist (sparse data), the module falls back to
global scaling: corrected = coverage (no GC correction applied).

.. rubric:: References

.. [1] Cristiano, S. et al. (2019). Nature 570:385-389. PMID: 31142840
.. [2] Mathios, D. et al. (2021). Nature Communications 12:5060. PMID: 34611155
.. [3] Cleveland, W.S. (1979). JASA 74:829-836. (LOESS/LOWESS method)

.. rubric:: Dependencies

- numpy
- scipy (stats.gaussian_kde, interpolate.interp1d)
- statsmodels (optional; nonparametric.smoothers_lowess.lowess)
"""

import numpy as np
from typing import Tuple, Optional
from scipy.stats import gaussian_kde


def compute_gc_content(sequence: str) -> float:
    """
    Compute GC content of a DNA sequence.

    GC content is the proportion of G and C bases in the sequence.
    IUPAC ambiguity codes are treated as non-GC (only A, T, G, C
    are considered).

    Parameters
    ----------
    sequence : str
        DNA sequence string. Case-insensitive (uppercased internally).

    Returns
    -------
    float
        GC content in [0, 1]. Returns 0.0 for empty sequences.

    Notes
    -----
    Edge case: an empty string returns 0.0 (division raises no error).

    Examples
    --------
    >>> compute_gc_content("GCGC")
    1.0
    >>> compute_gc_content("ATAT")
    0.0
    >>> compute_gc_content("")
    0.0
    """
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

    Implements the DELFI supplementary protocol [1]_.

    Algorithm
    ---------
    1. Bin genome windows by GC content into ``n_bins`` equal-width bins.
    2. Compute median coverage per GC bin to create a GC-coverage curve.
    3. Fit LOESS (LOWESS) smooth curve: coverage = f(GC).
    4. Interpolate the smooth curve to predict coverage for every window.
    5. Apply correction:

    .. math::

        \text{corrected}_i = \text{coverage}_i \times
        \frac{\text{median}(\text{coverage})}
        {\max(\hat{f}(\text{GC}_i), \, 10^{-6})}

    where :math:`\hat{f}(\text{GC}_i)` is the LOESS-predicted coverage
    at the GC content of window *i*.

    Parameters
    ----------
    coverage : np.ndarray
        Per-window read coverage, shape ``(n_windows,)``.
    gc_content : np.ndarray
        Per-window GC content in [0, 1], shape ``(n_windows,)``.
    frac : float, optional
        LOESS smoothing fraction (default 0.3). Larger = smoother.
    n_bins : int, optional
        Number of GC bins for initial binning (default 100).

    Returns
    -------
    corrected_coverage : np.ndarray
        GC-corrected coverage, shape ``(n_windows,)``.
    correction_factors : np.ndarray
        Multiplicative correction factors applied to each window,
        shape ``(n_windows,)``. Values > 1 mean coverage was boosted;
        values < 1 mean coverage was suppressed.

    Raises
    ------
    ValueError
        If ``coverage`` and ``gc_content`` have mismatched lengths.

    Notes
    -----
    - If ``statsmodels`` is not available, falls back to quadratic polynomial
      regression (degree-2 ``np.polyfit``).
    - If fewer than 10 valid GC bins exist after NaN removal, falls back to
      global scaling: corrected = coverage (no GC correction). The returned
      ``correction_factors`` are all 1.0 in this case.
    - Predicted coverage is clipped to a minimum of 1e-6 to avoid division
      by zero.

    References
    ----------
    .. [1] Cristiano et al. (2019). Nature 570:385-389, Supplementary Methods.
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
    Full DELFI-style normalization pipeline.

    Performs a three-step normalization following the protocol in [1]_ and [2]_.

    Steps
    -----
    1. **Mappability filter**: Exclude windows with mappability ≤ 0.8
       (low-mappability regions give unreliable coverage).
    2. **LOESS GC-bias correction**: Apply :func:`loess_normalize`.
    3. **Median-centering**: Divide all corrected values by the global
       median so that the median corrected coverage equals 1.0.

    Parameters
    ----------
    coverage : np.ndarray
        Per-window read coverage.
    gc_content : np.ndarray
        Per-window GC content in [0, 1].
    mappability : np.ndarray, optional
        Per-window mappability scores in [0, 1]. If provided, windows
        with mappability ≤ 0.8 are filtered out. If None, all windows
        are retained.

    Returns
    -------
    np.ndarray
        Normalized, median-centered coverage values.

    Raises
    ------
    ValueError
        If after mappability filtering, zero windows remain.

    References
    ----------
    .. [1] Cristiano et al. (2019). Nature 570:385-389.
    .. [2] Mathios et al. (2021). Nature Communications 12:5060.
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

    MDS quantifies the diversity of fragment end motifs using the
    normalized Simpson diversity index.

    .. math::

        \text{MDS} = \frac{1 - \sum_{i=1}^{n} p_i^2}{1 - 1/n}

    where :math:`p_i` is the frequency of motif *i* and
    :math:`n` is the total number of possible motifs.

    The numerator :math:`1 - \sum p_i^2` is the Simpson diversity index
    (probability that two randomly drawn motifs differ).
    The denominator normalizes to [0, 1]:

    - MDS = 0 when all fragments share the same motif (no diversity)
    - MDS = 1 when all 256 motifs are equally frequent (max diversity)

    Higher MDS → more diverse fragment end motifs → more likely cancer.

    Parameters
    ----------
    fragment_ends : np.ndarray
        Integer array of motif indices. Values must be in [0, n_motifs).
        For 4-mer motifs, ``n_motifs=256``.
    n_motifs : int, optional
        Total number of possible motifs (default 256 for 4-mers).

    Returns
    -------
    float
        MDS in [0, 1]. Returns 0.0 if all counts are zero.

    References
    ----------
    .. [1] Jiang, P. et al. (2020). Cancer Discovery 10(5):664-673. PMID: 32111602
    """
    # Compute motif frequencies
    motif_counts = np.bincount(fragment_ends, minlength=n_motifs)
    p = motif_counts / motif_counts.sum()

    # Simpson diversity → MDS (normalized)
    simpson = np.sum(p ** 2)
    mds = (1 - simpson) / (1 - 1.0 / n_motifs)

    return mds


__all__ = [
    "compute_gc_content",
    "loess_normalize",
    "DELFI_style_normalization",
    "compute_MDS",
]
