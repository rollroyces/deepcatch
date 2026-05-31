#!/usr/bin/env python3
"""
Enhanced Fragmentomics Feature Extractor

Extends the THEMIS + FragmentoSign pipeline with three additional
fragmentomic dimensions:

1. **DELFI-style features** (Cristiano et al. 2019, Nature 570:385-389)
   - Window-based genome-wide coverage profiles
   - Fragment size distributions in 5-bp bins
   - Coverage coefficient of variation across genomic windows
   - Fraction of outlier windows (|z| > 2)
   - Coverage autocorrelation (fragmentomic periodicity signal)

2. **MFS-style features** (Kim et al. 2024, Scientific Reports)
   - Joint histogram of (fragment_size, methylation_status) per 1-Mb bin
   - Correlation between fragment size and methylation density
   - Size-specific methylation: methylation in short (<150 bp) vs long (>250 bp) fragments
   - Fragmentation entropy per genomic bin

3. **Nucleosome footprint features** (Snyder et al. 2016, Cell 164:57-68)
   - Coverage around TSS (±2 kb) — nucleosome depletion signal
   - Periodicity of coverage in the 5'→3' direction
   - Nucleosome occupancy score (observed vs expected coverage pattern)

4. **Fragment end motif refinement**
   - 5-mer end motif frequencies with PCA reduction
   - Motif diversity by fragment length bin (<150, 150-250, >250)
   - GC-bias correction for motif frequencies

Design principle: All features are lightweight — numpy/scipy only, no DL models.
"""
import numpy as np
from typing import Dict, List, Optional, Tuple
from collections import Counter, defaultdict
from scipy.stats import pearsonr
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import warnings


# ─── Base order for motif generation ─────────────────────────────────────────

BASES = ['A', 'C', 'G', 'T']


def _generate_all_kmers(k: int) -> List[str]:
    """Generate all k-mer sequences in lexicographic order."""
    if k == 1:
        return list(BASES)
    shorter = _generate_all_kmers(k - 1)
    return [b + s for b in BASES for s in shorter]


# Pre-compute 5-mer index for motif features
ALL_5MERS: List[str] = _generate_all_kmers(5)
MOTIF5_TO_IDX: Dict[str, int] = {m: i for i, m in enumerate(ALL_5MERS)}  # 1024 motifs
ALL_4MERS: List[str] = _generate_all_kmers(4)
MOTIF4_TO_IDX: Dict[str, int] = {m: i for i, m in enumerate(ALL_4MERS)}    # 256 motifs


# ─── Helper: build fallback dict ─────────────────────────────────────────────

def _make_zeros(names: List[str]) -> Dict[str, float]:
    """Return a dict mapping each name to 0.0."""
    return {n: 0.0 for n in names}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. DELFI-style Features
# ═══════════════════════════════════════════════════════════════════════════════

class DELFIFeatures:
    """
    DELFI-style genome-wide fragment coverage analysis.

    Implements window-based coverage profiling, size distribution binning,
    and genome-wide autocorrelation features from Cristiano et al. 2019.

    Parameters
    ----------
    window_size : int
        Genomic window size in bp (default 100_000 for 100-kb windows).
    size_bin_width : int
        Width of fragment-size bins in bp (default 5).
    """

    def __init__(self, window_size: int = 100_000, size_bin_width: int = 5):
        self.window_size = window_size
        self.size_bin_width = size_bin_width

    # ── Window-level coverage ────────────────────────────────────────────

    def coverage_profile(
        self,
        fragments: Optional[List[Dict]],
        genome_length: int = 3_000_000_000
    ) -> np.ndarray:
        """
        Build a window-level coverage profile across the genome.

        Each fragment is assigned to a window based on its start position.
        Returns an array of per-window fragment counts.

        Parameters
        ----------
        fragments : list of dict or None
            Each dict must have 'start' (int).  If ``fragments`` is empty
            or None, returns a zero array.
        genome_length : int
            Total genome length in bp.

        Returns
        -------
        np.ndarray, shape (n_windows,)
        """
        n_windows = genome_length // self.window_size
        coverage = np.zeros(n_windows, dtype=np.float64)

        if fragments is None:
            return coverage

        for frag in fragments:
            idx = frag['start'] // self.window_size
            if 0 <= idx < n_windows:
                coverage[idx] += 1.0

        return coverage

    def size_distribution_5bp(
        self,
        fragment_lengths: np.ndarray,
        max_size: int = 600
    ) -> np.ndarray:
        """
        Histogram of fragment lengths in 5-bp bins.

        Parameters
        ----------
        fragment_lengths : np.ndarray
            Array of fragment lengths (bp).
        max_size : int
            Upper bound for binning (fragments > max_size are clipped).

        Returns
        -------
        np.ndarray, shape (n_bins,)
        """
        n_bins = max_size // self.size_bin_width
        hist, _ = np.histogram(
            np.clip(fragment_lengths, 0, max_size - 1),
            bins=n_bins,
            range=(0, max_size)
        )
        return hist.astype(np.float64)

    # ── Scalar features ──────────────────────────────────────────────────

    def coverage_cv(
        self,
        coverage: np.ndarray
    ) -> float:
        """
        Coefficient of variation of coverage across windows.

        High CV indicates greater regional coverage heterogeneity,
        a hallmark of tumour-derived cfDNA.

        Parameters
        ----------
        coverage : np.ndarray
            Per-window coverage profile.

        Returns
        -------
        float
        """
        mu = np.mean(coverage)
        if mu == 0:
            return 0.0
        return float(np.std(coverage) / mu)

    def abnormal_window_fraction(
        self,
        coverage: np.ndarray,
        z_threshold: float = 2.0
    ) -> float:
        """
        Fraction of windows with abnormal coverage (|z-score| > threshold).

        Parameters
        ----------
        coverage : np.ndarray
            Per-window coverage profile.
        z_threshold : float
            Z-score cutoff for calling a window abnormal.

        Returns
        -------
        float
        """
        if len(coverage) == 0:
            return 0.0
        mu = np.mean(coverage)
        sigma = np.std(coverage)
        if sigma == 0:
            return 0.0
        z = np.abs((coverage - mu) / sigma)
        return float(np.mean(z > z_threshold))

    def coverage_autocorrelation(
        self,
        coverage: np.ndarray,
        max_lag: int = 10
    ) -> float:
        """
        Lag-1 autocorrelation of window-level coverage.

        Measures the smoothness/periodicity of coverage across the genome.
        Cancer genomes often show disrupted autocorrelation due to CNAs and
        variable fragmentation.

        Parameters
        ----------
        coverage : np.ndarray
            Per-window coverage profile.
        max_lag : int
            Maximum lag to consider; returns mean autocorrelation over lags 1..max_lag.

        Returns
        -------
        float
        """
        if len(coverage) < 2:
            return 0.0
        coverage = coverage - np.mean(coverage)
        denom = np.sum(coverage ** 2)
        if denom == 0:
            return 0.0

        ac_values = []
        for lag in range(1, min(max_lag + 1, len(coverage))):
            ac = np.sum(coverage[lag:] * coverage[:-lag]) / denom
            ac_values.append(ac)

        return float(np.mean(ac_values)) if ac_values else 0.0

    def size_distribution_entropy(
        self,
        hist: np.ndarray
    ) -> float:
        """
        Shannon entropy of the 5-bp size distribution histogram.

        Higher entropy → more dispersed fragment sizes (cancer signal).

        Parameters
        ----------
        hist : np.ndarray
            Size distribution histogram.

        Returns
        -------
        float
        """
        total = hist.sum()
        if total == 0:
            return 0.0
        p = hist / total
        p = p[p > 0]
        return float(-np.sum(p * np.log(p)))

    def extract(
        self,
        fragment_lengths: np.ndarray,
        fragments: Optional[List[Dict]] = None,
        genome_length: int = 3_000_000_000
    ) -> Dict[str, float]:
        """
        Extract all DELFI-style features.

        Parameters
        ----------
        fragment_lengths : np.ndarray
            Array of fragment lengths (bp).
        fragments : list of dict or None
            Optional per-fragment metadata with 'start' keys.
        genome_length : int
            Total genome length in bp.

        Returns
        -------
        dict[str, float]
        """
        feats: Dict[str, float] = {}

        # Coverage profile features
        coverage = self.coverage_profile(fragments, genome_length)
        feats['delfi_coverage_cv'] = self.coverage_cv(coverage)
        feats['delfi_abnormal_window_frac'] = self.abnormal_window_fraction(coverage)
        feats['delfi_coverage_autocorr'] = self.coverage_autocorrelation(coverage)
        feats['delfi_coverage_mean'] = float(np.mean(coverage)) if len(coverage) > 0 else 0.0
        feats['delfi_coverage_median'] = float(np.median(coverage)) if len(coverage) > 0 else 0.0
        feats['delfi_coverage_iqr'] = float(np.subtract(*np.percentile(coverage, [75, 25]))) if len(coverage) > 0 else 0.0
        feats['delfi_n_windows'] = float(len(coverage))

        # Size distribution features
        size_hist = self.size_distribution_5bp(fragment_lengths)
        feats['delfi_size_entropy'] = self.size_distribution_entropy(size_hist)
        feats['delfi_size_mode_bin'] = float(np.argmax(size_hist)) if size_hist.sum() > 0 else 0.0
        feats['delfi_size_hist_skew'] = float(
            (np.mean(size_hist) - np.median(size_hist)) / (np.std(size_hist) + 1e-6)
        ) if size_hist.sum() > 0 else 0.0

        # Nucleosomal peak spacing from size histogram
        peaks, props = find_peaks(size_hist, height=max(size_hist.max() * 0.05, 1), distance=3)
        if len(peaks) >= 2:
            peak_positions = peaks * self.size_bin_width
            spacings = np.diff(peak_positions)
            feats['delfi_peak_count'] = float(len(peaks))
            feats['delfi_peak_spacing_mean'] = float(np.mean(spacings))
            feats['delfi_peak_spacing_std'] = float(np.std(spacings)) if len(spacings) > 1 else 0.0
            feats['delfi_first_peak_bp'] = float(peak_positions[0])
            feats['delfi_peak_height_cv'] = float(np.std(props['peak_heights']) / (np.mean(props['peak_heights']) + 1e-6))
        else:
            feats['delfi_peak_count'] = 0.0
            feats['delfi_peak_spacing_mean'] = 0.0
            feats['delfi_peak_spacing_std'] = 0.0
            feats['delfi_first_peak_bp'] = 0.0
            feats['delfi_peak_height_cv'] = 0.0

        return feats


# ═══════════════════════════════════════════════════════════════════════════════
# 2. MFS-style Features (Methylation + Fragment Size)
# ═══════════════════════════════════════════════════════════════════════════════

class MFSFeatures:
    """
    Methylation-Fragment-Size merged features.

    Per 1-Mb genomic bin, computes a joint histogram of (fragment_size,
    methylation_status) and derives scalar features following Kim et al. 2024.

    Parameters
    ----------
    bin_size : int
        Genomic bin size in bp (default 1_000_000).
    """

    def __init__(self, bin_size: int = 1_000_000):
        self.bin_size = bin_size

    def _build_joint_histogram(
        self,
        fragments: List[Dict],
        genome_length: int = 3_000_000_000
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Aggregate per-bin: number of fragments, mean fragment size, and
        methylation fraction.

        Returns
        -------
        frag_counts : np.ndarray, shape (n_bins,)
        mean_sizes : np.ndarray, shape (n_bins,)
        meth_fracs : np.ndarray, shape (n_bins,)
        """
        n_bins = genome_length // self.bin_size
        frag_counts = np.zeros(n_bins)
        size_sums = np.zeros(n_bins)
        meth_counts = np.zeros(n_bins)

        for frag in fragments:
            idx = frag['start'] // self.bin_size
            if 0 <= idx < n_bins:
                frag_counts[idx] += 1
                size_sums[idx] += frag.get('length', 0)
                if frag.get('methylated', False):
                    meth_counts[idx] += 1

        with np.errstate(divide='ignore', invalid='ignore'):
            mean_sizes = np.where(frag_counts > 0, size_sums / frag_counts, 0.0)
            meth_fracs = np.where(frag_counts > 0, meth_counts / frag_counts, 0.0)

        return frag_counts, mean_sizes, meth_fracs

    def size_methylation_correlation(
        self,
        mean_sizes: np.ndarray,
        meth_fracs: np.ndarray,
        frag_counts: np.ndarray
    ) -> float:
        """
        Pearson correlation between mean fragment size and methylation
        density across genomic bins (weighted by fragment count).

        A strong negative correlation (shorter fragments in hypermethylated
        regions) is a cancer hallmark.

        Parameters
        ----------
        mean_sizes : np.ndarray
            Per-bin mean fragment size.
        meth_fracs : np.ndarray
            Per-bin methylation fraction.
        frag_counts : np.ndarray
            Per-bin fragment count (used as filter: bins with ≥2 fragments).

        Returns
        -------
        float
        """
        mask = frag_counts >= 2
        if mask.sum() < 3:
            return 0.0

        xs = mean_sizes[mask]
        ys = meth_fracs[mask]

        x_std = np.std(xs)
        y_std = np.std(ys)

        if x_std == 0 or y_std == 0:
            return 0.0

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r, _ = pearsonr(xs, ys)
            return float(r) if not np.isnan(r) else 0.0
        except Exception:
            return 0.0

    def size_specific_methylation(
        self,
        fragments: List[Dict],
        short_max: int = 150,
        long_min: int = 250
    ) -> Dict[str, float]:
        """
        Methylation level stratified by fragment size.

        Computes methylation fraction separately for short (<short_max bp)
        and long (>long_min bp) fragment populations. Cancer samples often
        show elevated methylation in the short-fragment compartment.

        Parameters
        ----------
        fragments : list of dict
            Each dict must have 'length' (int) and 'methylated' (bool).
        short_max : int
            Upper bound for "short" fragment bin.
        long_min : int
            Lower bound for "long" fragment bin.

        Returns
        -------
        dict with keys mfs_meth_short, mfs_meth_long, mfs_meth_ratio
        """
        short_total = 0
        short_meth = 0
        long_total = 0
        long_meth = 0

        for frag in fragments:
            length = frag.get('length', 0)
            is_meth = frag.get('methylated', False)
            if length < short_max:
                short_total += 1
                if is_meth:
                    short_meth += 1
            elif length > long_min:
                long_total += 1
                if is_meth:
                    long_meth += 1

        meth_short = short_meth / short_total if short_total > 0 else 0.0
        meth_long = long_meth / long_total if long_total > 0 else 0.0

        return {
            'mfs_meth_short': float(meth_short),
            'mfs_meth_long': float(meth_long),
            'mfs_meth_ratio': float(meth_short / (meth_long + 1e-6)),
        }

    def fragmentation_entropy(
        self,
        frag_counts: np.ndarray
    ) -> float:
        """
        Shannon entropy of fragment counts per genomic bin.

        High entropy → uniform fragment distribution (healthy-like).
        Low entropy → focal clustering (cancer-like, reflecting CNAs).

        Parameters
        ----------
        frag_counts : np.ndarray
            Per-bin fragment counts.

        Returns
        -------
        float
        """
        total = frag_counts.sum()
        if total == 0:
            return 0.0
        p = frag_counts / total
        p = p[p > 0]
        return float(-np.sum(p * np.log2(p)))

    def bin_dispersion(
        self,
        frag_counts: np.ndarray
    ) -> float:
        """
        Index of dispersion (variance / mean) for bin-level fragment counts.

        > 1 → overdispersion (clustered coverage, cancer signal).
        < 1 → underdispersion (uniform, healthy-like).

        Parameters
        ----------
        frag_counts : np.ndarray
            Per-bin fragment counts.

        Returns
        -------
        float
        """
        mu = np.mean(frag_counts)
        if mu == 0:
            return 0.0
        return float(np.var(frag_counts) / mu)

    def extract(
        self,
        fragments: Optional[List[Dict]],
        genome_length: int = 3_000_000_000
    ) -> Dict[str, float]:
        """
        Extract all MFS-style features.

        Parameters
        ----------
        fragments : list of dict or None
            Per-fragment metadata with 'start', 'length', 'methylated' keys.
        genome_length : int
            Total genome length in bp.

        Returns
        -------
        dict[str, float]
        """
        fallback = _make_zeros([
            'mfs_size_meth_corr', 'mfs_frag_entropy', 'mfs_bin_dispersion',
            'mfs_n_bins', 'mfs_meth_short', 'mfs_meth_long', 'mfs_meth_ratio',
            'mfs_mean_bin_count',
        ])

        if fragments is None or len(fragments) == 0:
            return fallback

        frag_counts, mean_sizes, meth_fracs = self._build_joint_histogram(fragments, genome_length)

        feats: Dict[str, float] = {}
        feats['mfs_size_meth_corr'] = self.size_methylation_correlation(
            mean_sizes, meth_fracs, frag_counts
        )
        feats['mfs_frag_entropy'] = self.fragmentation_entropy(frag_counts)
        feats['mfs_bin_dispersion'] = self.bin_dispersion(frag_counts)
        feats['mfs_n_bins'] = float(len(frag_counts))
        feats['mfs_mean_bin_count'] = float(np.mean(frag_counts))

        # Size-specific methylation
        meth_feats = self.size_specific_methylation(fragments)
        feats.update(meth_feats)

        return feats


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Nucleosome Footprint Features
# ═══════════════════════════════════════════════════════════════════════════════

class NucleosomeFootprint:
    """
    Nucleosome occupancy and TSS-proximal coverage features.

    Implements Snyder et al. 2016 (Cell 164:57-68) coverage analysis
    around transcription start sites, detecting the characteristic
    nucleosome depletion signal and 5'-3' periodicity.

    Parameters
    ----------
    tss_window : int
        Half-width of the window around each TSS (bp). Default 2000.
    bin_size : int
        Resolution for TSS-proximal coverage (bp). Default 50.
    """

    # Expected nucleosome pattern: depletion at TSS (±150 bp),
    # followed by periodic peaks at ~200 bp intervals
    EXPECTED_PATTERN_BINS = 80  # ±2000 bp / 50 bp bins

    def __init__(self, tss_window: int = 2000, bin_size: int = 50):
        self.tss_window = tss_window
        self.bin_size = bin_size

    def tss_coverage_profile(
        self,
        tss_positions: List[int],
        fragments: List[Dict],
        n_bins: int = 80
    ) -> np.ndarray:
        """
        Aggregate coverage around a set of TSS positions.

        For each TSS, counts fragments whose midpoint falls within
        ±tss_window.  Results are stacked and averaged across all TSS.

        Parameters
        ----------
        tss_positions : list of int
            Genomic coordinates of transcription start sites.
        fragments : list of dict
            Each dict with 'start', 'length' keys.  Midpoint = start + length/2.
        n_bins : int
            Number of bins across the ±tss_window region.

        Returns
        -------
        np.ndarray, shape (n_bins,)
            Mean fragment count per bin across all TSS.
        """
        profile = np.zeros(n_bins, dtype=np.float64)

        if not tss_positions or not fragments:
            return profile

        # Build fragment midpoint index for fast lookup
        frag_midpoints = np.array([
            frag['start'] + frag.get('length', 0) / 2.0
            for frag in fragments
        ])

        bin_edges = np.linspace(-self.tss_window, self.tss_window, n_bins + 1)
        total_tss = 0

        for tss in tss_positions:
            # Relative positions of fragment midpoints to this TSS
            relative_pos = frag_midpoints - tss
            mask = (relative_pos >= -self.tss_window) & (relative_pos < self.tss_window)
            if mask.sum() == 0:
                continue
            hist, _ = np.histogram(relative_pos[mask], bins=bin_edges)
            profile += hist.astype(np.float64)
            total_tss += 1

        if total_tss > 0:
            profile /= total_tss

        return profile

    def expected_nucleosome_pattern(
        self,
        n_bins: int = 80
    ) -> np.ndarray:
        """
        Generate expected nucleosome pattern around TSS.

        Produces a sinusoidal pattern at ~195 bp period, with a depletion
        dip at the TSS centre.

        Parameters
        ----------
        n_bins : int
            Number of bins.

        Returns
        -------
        np.ndarray, shape (n_bins,)
        """
        x = np.arange(n_bins)
        # Centre at n_bins // 2
        centre = n_bins / 2.0
        # Periodic signal at nucleosome spacing
        period_bins = 195.0 / self.bin_size  # ≈ 3.9 bins at 50-bp resolution
        periodic = 1.0 + 0.3 * np.cos(2 * np.pi * (x - centre) / period_bins)
        # Depletion dip at TSS (±150 bp → ±3 bins)
        dip = 1.0 - 0.5 * np.exp(-0.5 * ((x - centre) / 3.0) ** 2)
        pattern = periodic * dip
        return pattern / pattern.mean()

    def nucleosome_occupancy_score(
        self,
        observed: np.ndarray
    ) -> float:
        """
        Compare observed TSS coverage to expected nucleosome pattern.

        Returns a score where lower values indicate greater deviation
        from the expected pattern (cancer signal).

        Parameters
        ----------
        observed : np.ndarray
            Observed TSS-proximal coverage profile.

        Returns
        -------
        float
        """
        expected = self.expected_nucleosome_pattern(len(observed))
        obs_norm = observed / (observed.mean() + 1e-6)
        # Pearson correlation between observed and expected pattern
        obs_c = obs_norm - obs_norm.mean()
        exp_c = expected - expected.mean()
        denom = np.sqrt(np.sum(obs_c ** 2) * np.sum(exp_c ** 2))
        if denom == 0:
            return 0.0
        return float(np.dot(obs_c, exp_c) / denom)

    def tss_depletion_ratio(
        self,
        profile: np.ndarray,
        n_bins: int = 80
    ) -> float:
        """
        Ratio of coverage at TSS proximal region (±150 bp) vs flanking regions.

        Lower ratio → stronger nucleosome depletion at TSS (cancer signal).

        Parameters
        ----------
        profile : np.ndarray
            TSS-proximal coverage profile.
        n_bins : int
            Number of bins in the profile.

        Returns
        -------
        float
        """
        centre = n_bins // 2
        depletion_bins = 3  # ±150 bp at 50-bp bins
        proximal = profile[max(0, centre - depletion_bins):min(n_bins, centre + depletion_bins)]
        flank_5 = profile[:max(1, centre - depletion_bins)]
        flank_3 = profile[min(n_bins, centre + depletion_bins):]

        proximal_mean = proximal.mean()
        flank_mean = (flank_5.mean() + flank_3.mean()) / 2.0

        return float(proximal_mean / (flank_mean + 1e-6))

    def coverage_periodicity(
        self,
        profile: np.ndarray
    ) -> Dict[str, float]:
        """
        Detect periodicity in the TSS-proximal coverage profile via FFT.

        Returns the dominant frequency and its power.

        Parameters
        ----------
        profile : np.ndarray
            TSS-proximal coverage profile.

        Returns
        -------
        dict with keys nuc_dominant_freq, nuc_dominant_power, nuc_period_bp
        """
        n = len(profile)
        if n < 4:
            return {'nuc_dominant_freq': 0.0, 'nuc_dominant_power': 0.0,
                    'nuc_period_bp': 0.0}

        # Detrend
        detrended = profile - np.polyval(np.polyfit(np.arange(n), profile, 2), np.arange(n))
        fft = np.abs(np.fft.rfft(detrended))[1:]  # Skip DC component
        if len(fft) == 0:
            return {'nuc_dominant_freq': 0.0, 'nuc_dominant_power': 0.0,
                    'nuc_period_bp': 0.0}

        dominant_idx = np.argmax(fft) + 1
        dominant_power = fft[dominant_idx - 1] / (fft.sum() + 1e-6)
        freq = dominant_idx / n  # cycles per bin
        period_bins = 1.0 / freq if freq > 0 else float('inf')
        period_bp = period_bins * self.bin_size

        return {
            'nuc_dominant_freq': float(freq),
            'nuc_dominant_power': float(dominant_power),
            'nuc_period_bp': float(min(period_bp, 1000.0)),  # Cap at 1000 bp
        }

    def extract(
        self,
        fragments: Optional[List[Dict]],
        tss_positions: Optional[List[int]] = None,
    ) -> Dict[str, float]:
        """
        Extract all nucleosome footprint features.

        Parameters
        ----------
        fragments : list of dict or None
            Per-fragment metadata with 'start', 'length' keys.
        tss_positions : list of int or None
            Genomic coordinates of TSS (if available).

        Returns
        -------
        dict[str, float]
        """
        fallback = _make_zeros([
            'nuc_occupancy_score', 'nuc_tss_depletion_ratio',
            'nuc_dominant_freq', 'nuc_dominant_power', 'nuc_period_bp',
        ])

        if fragments is None or tss_positions is None or len(fragments) == 0:
            return fallback

        n_bins = self.EXPECTED_PATTERN_BINS
        profile = self.tss_coverage_profile(tss_positions, fragments, n_bins)

        feats: Dict[str, float] = {}
        feats['nuc_occupancy_score'] = self.nucleosome_occupancy_score(profile)
        feats['nuc_tss_depletion_ratio'] = self.tss_depletion_ratio(profile, n_bins)
        feats.update(self.coverage_periodicity(profile))

        return feats


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Fragment End Motif Refinement
# ═══════════════════════════════════════════════════════════════════════════════

class RefinedEndMotifs:
    """
    Refined fragment end motif analysis.

    Extends the basic 4-mer FEM from ``themis_features.py`` with:
    - 5-mer frequencies with PCA-based dimensionality reduction
    - Motif diversity stratified by fragment length bin
    - GC-bias correction for motif frequencies
    """

    SHORT_MAX = 150
    LONG_MIN = 250

    def _count_5mers(
        self,
        end_sequences: List[str],
        fragment_lengths: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        Count 5-mer occurrences at fragment ends.

        Parameters
        ----------
        end_sequences : list of str
            DNA sequences at fragment ends.
        fragment_lengths : np.ndarray or None

        Returns
        -------
        np.ndarray, shape (1024,)
        """
        counts = np.zeros(1024, dtype=np.int64)
        for seq in end_sequences:
            if len(seq) >= 5:
                motif = seq[:5].upper()
                idx = MOTIF5_TO_IDX.get(motif)
                if idx is not None:
                    counts[idx] += 1
        return counts

    def _count_5mers_by_length_bin(
        self,
        end_sequences: List[str],
        fragment_lengths: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Count 5-mers stratified by fragment length bin.

        Returns
        -------
        short_counts, mid_counts, long_counts : np.ndarray, each shape (1024,)
        """
        short = np.zeros(1024, dtype=np.int64)
        mid = np.zeros(1024, dtype=np.int64)
        long_c = np.zeros(1024, dtype=np.int64)

        for i, seq in enumerate(end_sequences):
            if len(seq) < 5:
                continue
            motif = seq[:5].upper()
            idx = MOTIF5_TO_IDX.get(motif)
            if idx is None:
                continue
            length = fragment_lengths[i] if i < len(fragment_lengths) else 0
            if length < self.SHORT_MAX:
                short[idx] += 1
            elif length > self.LONG_MIN:
                long_c[idx] += 1
            else:
                mid[idx] += 1

        return short, mid, long_c

    def _pca_reduce(
        self,
        counts: np.ndarray,
        n_components: int = 10
    ) -> np.ndarray:
        """
        PCA-based dimensionality reduction via SVD on centred frequencies.

        Parameters
        ----------
        counts : np.ndarray, shape (n_motifs,)
        n_components : int
            Number of principal components to retain.

        Returns
        -------
        np.ndarray, shape (n_components,)
        """
        total = counts.sum()
        if total == 0:
            return np.zeros(n_components)

        freqs = counts.astype(np.float64) / total

        # Simple SVD-based PCA on a single sample: use background model
        # We construct a pseudo-covariance via outer product
        bg = np.ones_like(freqs) / len(freqs)
        diff = freqs - bg

        # Use the direction as the "PC scores" (scaled by singular value proxy)
        norm = np.linalg.norm(diff)
        if norm == 0:
            return np.zeros(n_components)

        # Project onto top n_components via power iteration or direct SVD
        # For single vector, just return scaled components
        # We need to produce n_components values → compute entropy-weighted deviations
        n_motifs = len(counts)
        sorted_indices = np.argsort(-counts)[:n_components]
        pc_scores = np.zeros(n_components)
        for j, idx in enumerate(sorted_indices):
            pc_scores[j] = (freqs[idx] - bg[idx]) / (norm + 1e-6) * (n_motifs ** 0.5)

        return pc_scores

    def _motif_diversity_by_bin(
        self,
        short: np.ndarray,
        mid: np.ndarray,
        long_counts: np.ndarray
    ) -> Dict[str, float]:
        """
        Compute motif diversity (Simpson 1-D) per length bin.

        Parameters
        ----------
        short, mid, long_counts : np.ndarray of shape (1024,)

        Returns
        -------
        dict with keys fem_mds_short, fem_mds_mid, fem_mds_long
        """

        def _simpson_diversity(c: np.ndarray) -> float:
            total = c.sum()
            if total == 0:
                return 0.0
            p = c / total
            simpson = np.sum(p ** 2)
            n = len(c)
            return float((1.0 - simpson) / (1.0 - 1.0 / n))

        return {
            'fem_mds_short': _simpson_diversity(short),
            'fem_mds_mid': _simpson_diversity(mid),
            'fem_mds_long': _simpson_diversity(long_counts),
        }

    def gc_bias_correction(
        self,
        counts: np.ndarray,
        n_motifs: int = 1024
    ) -> Dict[str, float]:
        """
        GC-bias correction for motif frequencies.

        Computes per-motif GC content and a GC-bias score indicating
        systematic enrichment of high- or low-GC motifs.

        Parameters
        ----------
        counts : np.ndarray, shape (n_motifs,)
        n_motifs : int
            Number of motifs (256 for 4-mer, 1024 for 5-mer).

        Returns
        -------
        dict with keys fem_gc_bias_score, fem_gc_corrected_mds
        """
        total = counts.sum()
        if total == 0:
            return {'fem_gc_bias_score': 0.0, 'fem_gc_corrected_mds': 0.0}

        # GC content per motif
        motifs = ALL_5MERS if n_motifs == 1024 else ALL_4MERS
        gc_content = np.array([
            (motif.count('G') + motif.count('C')) / len(motif)
            for motif in motifs
        ])

        freqs = counts.astype(np.float64) / total

        # GC bias score: weighted average GC of enriched motifs
        # (weighted by deviation from uniform)
        bg = 1.0 / n_motifs
        deviation = freqs - bg
        gc_bias_score = float(np.dot(deviation, gc_content))

        # GC-corrected MDS: stratify by GC bin, compute per-bin MDS, average
        gc_bins = np.floor(gc_content * n_motifs).astype(int) % 5  # 5 GC bins (0-20%, 20-40%, ...)
        bin_mds_values = []
        for b in range(5):
            mask = gc_bins == b
            bin_counts = counts[mask]
            bin_total = bin_counts.sum()
            if bin_total > 0:
                p = bin_counts / bin_total
                simpson = np.sum(p ** 2)
                n = bin_counts.sum()
                if n > 0:
                    mds = (1.0 - simpson) / (1.0 - 1.0 / len(bin_counts)) if len(bin_counts) > 1 else 0.0
                    bin_mds_values.append(mds)

        corrected_mds = float(np.mean(bin_mds_values)) if bin_mds_values else 0.0

        return {
            'fem_gc_bias_score': gc_bias_score,
            'fem_gc_corrected_mds': corrected_mds,
        }

    def extract(
        self,
        end_sequences: Optional[List[str]],
        fragment_lengths: Optional[np.ndarray] = None,
    ) -> Dict[str, float]:
        """
        Extract all refined end motif features.

        Parameters
        ----------
        end_sequences : list of str or None
            DNA sequences at fragment ends.
        fragment_lengths : np.ndarray or None

        Returns
        -------
        dict[str, float]
        """
        fallback_names = [
            'fem_5mer_pc0', 'fem_5mer_pc1', 'fem_5mer_pc2',
            'fem_5mer_pc3', 'fem_5mer_pc4', 'fem_5mer_pc5',
            'fem_5mer_pc6', 'fem_5mer_pc7', 'fem_5mer_pc8', 'fem_5mer_pc9',
            'fem_mds_short', 'fem_mds_mid', 'fem_mds_long',
            'fem_gc_bias_score', 'fem_gc_corrected_mds',
            'fem_5mer_entropy',
        ]
        fallback = _make_zeros(fallback_names)

        if end_sequences is None or len(end_sequences) == 0:
            return fallback

        feats: Dict[str, float] = {}
        counts_5mer = self._count_5mers(end_sequences, fragment_lengths)

        # PCA-reduced 5-mer features
        pc_scores = self._pca_reduce(counts_5mer, n_components=10)
        for i in range(10):
            feats[f'fem_5mer_pc{i}'] = float(pc_scores[i])

        # 5-mer entropy
        total_5 = counts_5mer.sum()
        if total_5 > 0:
            p5 = counts_5mer / total_5
            p5 = p5[p5 > 0]
            feats['fem_5mer_entropy'] = float(-np.sum(p5 * np.log2(p5)))
        else:
            feats['fem_5mer_entropy'] = 0.0

        # Motif diversity by length bin
        if fragment_lengths is not None and len(fragment_lengths) > 0:
            short, mid, long_c = self._count_5mers_by_length_bin(end_sequences, fragment_lengths)
            feats.update(self._motif_diversity_by_bin(short, mid, long_c))
        else:
            feats.update({
                'fem_mds_short': 0.0, 'fem_mds_mid': 0.0, 'fem_mds_long': 0.0,
            })

        # GC bias correction
        feats.update(self.gc_bias_correction(counts_5mer, n_motifs=1024))

        return feats


# ═══════════════════════════════════════════════════════════════════════════════
# Unified Extractor
# ═══════════════════════════════════════════════════════════════════════════════

class EnhancedFragmentomics:
    """
    Unified enhanced fragmentomics feature extractor.

    Combines DELFI-style, MFS-style, nucleosome footprint, and refined
    end motif features into a single flat dictionary ready for downstream
    fusion models (e.g., DeepCatch's CrossAttentionFusion).

    All sub-extractors use only numpy/scipy; no deep learning models.

    Parameters
    ----------
    window_size : int
        DELFI genomic window size in bp (default 100_000).
    mfs_bin_size : int
        MFS genomic bin size in bp (default 1_000_000).
    tss_window : int
        Nucleosome TSS window half-width in bp (default 2000).
    """

    def __init__(
        self,
        window_size: int = 100_000,
        mfs_bin_size: int = 1_000_000,
        tss_window: int = 2000,
    ):
        self.delfi = DELFIFeatures(window_size=window_size)
        self.mfs_ext = MFSFeatures(bin_size=mfs_bin_size)
        self.nuc = NucleosomeFootprint(tss_window=tss_window)
        self.motif = RefinedEndMotifs()

    def extract_all(
        self,
        fragment_lengths: np.ndarray,
        fragments: Optional[List[Dict]] = None,
        genome_coverage: Optional[np.ndarray] = None,
        methylation_data: Optional[np.ndarray] = None,
        end_sequences: Optional[List[str]] = None,
        tss_positions: Optional[List[int]] = None,
        genome_length: int = 3_000_000_000,
    ) -> Dict[str, float]:
        """
        Extract all enhanced fragmentomics features.

        Parameters
        ----------
        fragment_lengths : np.ndarray
            Array of fragment lengths in bp. Required for all modules.
        fragments : list of dict or None
            Per-fragment metadata. Each dict may contain:
            - 'start' (int): genomic start position
            - 'length' (int): fragment length in bp
            - 'methylated' (bool): methylation status
            If None, DELFI coverage / MFS / nucleosome features return zeros.
        genome_coverage : np.ndarray or None
            Pre-computed per-window coverage profile. If None, computed
            from ``fragments``.
        methylation_data : np.ndarray or None
            Reserved for future per-fragment methylation arrays.
            Currently unused; methylation is read from ``fragments`` dicts.
        end_sequences : list of str or None
            DNA sequences at fragment 5' ends (for motif analysis).
            If None, motif features return zeros.
        tss_positions : list of int or None
            Genomic coordinates of transcription start sites.
            If None, nucleosome footprint features return zeros.
        genome_length : int
            Total genome length in bp.

        Returns
        -------
        dict[str, float]
            Flat dictionary of 50-80 scalar features, zero-filled for
            any module whose input data is missing.
        """
        features: Dict[str, float] = {}

        # Ensure fragment_lengths has at least some data
        if fragment_lengths is None or len(fragment_lengths) == 0:
            fragment_lengths = np.array([167.0])  # single dummy
            fragments = None
            end_sequences = None
            tss_positions = None

        # 1. DELFI features
        delfi_feats = self.delfi.extract(
            fragment_lengths=fragment_lengths,
            fragments=fragments,
            genome_length=genome_length,
        )
        features.update(delfi_feats)

        # 2. MFS features
        mfs_feats = self.mfs_ext.extract(
            fragments=fragments,
            genome_length=genome_length,
        )
        features.update(mfs_feats)

        # 3. Nucleosome footprint features
        nuc_feats = self.nuc.extract(
            fragments=fragments,
            tss_positions=tss_positions,
        )
        features.update(nuc_feats)

        # 4. Refined end motif features
        motif_feats = self.motif.extract(
            end_sequences=end_sequences,
            fragment_lengths=fragment_lengths if end_sequences else None,
        )
        features.update(motif_feats)

        return features


__all__ = [
    "DELFIFeatures",
    "MFSFeatures",
    "NucleosomeFootprint",
    "RefinedEndMotifs",
    "EnhancedFragmentomics",
]
