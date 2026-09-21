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

        # Be paranoid about zero-variance inputs. Some numpy/scipy
        # versions return NaN from pearsonr() on constant input,
        # others return 0.0 — and either way the function shouldn't
        # be returning NaN silently. Force 0.0 here.
        with np.errstate(invalid="ignore", divide="ignore"):
            x_std = float(np.std(xs))
            y_std = float(np.std(ys))

        if not np.isfinite(x_std) or not np.isfinite(y_std) or x_std < 1e-12 or y_std < 1e-12:
            return 0.0

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r, _ = pearsonr(xs.astype(float), ys.astype(float))
            if not np.isfinite(r):
                return 0.0
            return float(r)
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

    # Default nucleosome parameters. These are the canonical values
    # for human cfDNA reported in Snyder et al. 2016 (Cell 164:57-68)
    # and Jiang et al. 2020 (Cancer Discovery 10:664-673):
    #   - Nucleosome repeat length ~195 bp (mean nucleosome + linker)
    #   - TSS depletion dip ~150 bp half-width
    #   - Dip amplitude 0.5 (50% coverage reduction at the TSS)
    #   - Sinusoidal amplitude 0.3 (modest periodicity signal)
    # These can be overridden per-instance for non-canonical cfDNA
    # sources (e.g. yeast nucleosomes ~165 bp, mouse ES ~190 bp).
    DEFAULT_NUCLEOSOME_PERIOD_BP = 195.0
    DEFAULT_TSS_DIP_HALFWIDTH_BP = 150.0
    DEFAULT_TSS_DIP_AMPLITUDE = 0.5
    DEFAULT_PERIODIC_AMPLITUDE = 0.3

    def __init__(
        self,
        tss_window: int = 2000,
        bin_size: int = 50,
        nucleosome_period_bp: float = DEFAULT_NUCLEOSOME_PERIOD_BP,
        tss_dip_halfwidth_bp: float = DEFAULT_TSS_DIP_HALFWIDTH_BP,
        tss_dip_amplitude: float = DEFAULT_TSS_DIP_AMPLITUDE,
        periodic_amplitude: float = DEFAULT_PERIODIC_AMPLITUDE,
    ):
        self.tss_window = tss_window
        self.bin_size = bin_size
        self.nucleosome_period_bp = nucleosome_period_bp
        self.tss_dip_halfwidth_bp = tss_dip_halfwidth_bp
        self.tss_dip_amplitude = tss_dip_amplitude
        self.periodic_amplitude = periodic_amplitude

    def tss_coverage_profile(
        self,
        tss_positions: List,
        fragments: List[Dict],
        n_bins: int = 80
    ) -> np.ndarray:
        """Aggregate coverage around a set of TSS positions.

        Each TSS must be a ``(chrom, pos)`` tuple. Fragments must have
        a ``chrom`` field; fragments on chromosomes with no TSS in the
        provided list are skipped. This per-chromosome matching
        prevents the previous bug where a fragment on chr5 was
        incorrectly counted against a chr1 TSS, inflating apparent
        coverage of unrelated promoters.

        For backwards compatibility, a ``TSS`` passed as a bare int is
        treated as chromosome "unknown"; fragments must also carry
        ``chrom="unknown"`` to match. The recommended path is to pass
        ``[(chrom, pos), ...]`` from a UCSC refFlat TSS table.

        Parameters
        ----------
        tss_positions : list of (chrom, pos) tuples or list of int
            Genomic TSS coordinates. Bare ints are accepted for legacy
            callers but the chrom-agnostic mode is unsafe for whole-
            genome cfDNA — prefer the (chrom, pos) form.
        fragments : list of dict
            Each dict must have ``start``, ``length``, and ``chrom``
            keys.  Midpoint = start + length/2.
        n_bins : int
            Number of bins across the ±tss_window region.

        Returns
        -------
        np.ndarray, shape (n_bins,)
            Mean fragment count per bin across all TSS, normalised by
            the number of TSS that received at least one matching
            fragment. Empty / no-match returns an all-zero profile.
        """
        profile = np.zeros(n_bins, dtype=np.float64)

        if not tss_positions or not fragments:
            return profile

        # Normalize TSS to (chrom, pos) tuples.
        tss_by_chrom: Dict[str, List[int]] = {}
        for tss in tss_positions:
            if isinstance(tss, (tuple, list)) and len(tss) == 2:
                chrom, pos = tss[0], tss[1]
            else:
                # Legacy bare-int path — only matches fragments with
                # chrom == "unknown". Documented as unsafe for WGS data.
                chrom, pos = "unknown", int(tss)
            tss_by_chrom.setdefault(chrom, []).append(int(pos))

        # Pre-group fragment midpoints per chromosome for O(N+M) work.
        frag_by_chrom_lists: Dict[str, List[float]] = {}
        for frag in fragments:
            chrom = frag.get("chrom", "unknown")
            mid = frag["start"] + frag.get("length", 0) / 2.0
            frag_by_chrom_lists.setdefault(chrom, []).append(mid)
        frag_by_chrom: Dict[str, np.ndarray] = {
            chrom: np.asarray(mids, dtype=np.float64)
            for chrom, mids in frag_by_chrom_lists.items()
        }

        bin_edges = np.linspace(-self.tss_window, self.tss_window, n_bins + 1)
        total_tss = 0

        for chrom, tss_list in tss_by_chrom.items():
            mids = frag_by_chrom.get(chrom)
            if mids is None or len(mids) == 0:
                # No fragments on this chromosome → skip all TSS.
                continue
            tss_arr = np.asarray(tss_list, dtype=np.float64)
            # Vectorize: relative position of every fragment midpoint
            # vs every TSS on this chromosome → (n_mids, n_tss).
            relative_pos = mids[:, None] - tss_arr[None, :]
            mask = (relative_pos >= -self.tss_window) & (
                relative_pos < self.tss_window
            )
            for j in range(tss_arr.shape[0]):
                col = relative_pos[:, j]
                col_mask = mask[:, j]
                if not col_mask.any():
                    continue
                hist, _ = np.histogram(col[col_mask], bins=bin_edges)
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

        Produces a sinusoidal pattern at the configured nucleosome
        repeat length (default ~195 bp from Snyder 2016 / Jiang 2020),
        with a Gaussian depletion dip at the TSS centre. The hard-
        coded 195 bp period and 150 bp dip half-width from the
        original implementation are now class-level defaults and can be
        overridden per-instance for non-canonical sources.

        Parameters
        ----------
        n_bins : int
            Number of bins.

        Returns
        -------
        np.ndarray, shape (n_bins,)
        """
        x = np.arange(n_bins)
        centre = n_bins / 2.0
        # Periodic signal at the configured nucleosome repeat length.
        period_bins = self.nucleosome_period_bp / self.bin_size
        periodic = 1.0 + self.periodic_amplitude * np.cos(
            2 * np.pi * (x - centre) / period_bins
        )
        # Gaussian depletion dip at TSS.
        dip_halfwidth_bins = self.tss_dip_halfwidth_bp / self.bin_size
        dip = 1.0 - self.tss_dip_amplitude * np.exp(
            -0.5 * ((x - centre) / max(dip_halfwidth_bins, 1e-6)) ** 2
        )
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
    - Top-N per-motif deviation scores (NOT a true PCA — see
      ``_top_motif_deviations`` for the algorithm description).
    - Motif diversity stratified by fragment length bin.
    - GC-bias correction for motif frequencies.
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

    def _top_motif_deviations(
        self,
        counts: np.ndarray,
        n_components: int = 10
    ) -> np.ndarray:
        """Top-N per-motif deviation scores (NOT a true PCA projection).

        This is a per-sample rank-weighted deviation vector, NOT a
        principal-component projection. There is no inter-sample
        covariance structure; each sample's motifs are scored against
        the uniform-background null independently. The function was
        previously named ``_pca_reduce`` and emitted keys prefixed
        ``fem_5mer_pc*``, which misleadingly suggested a true PCA. The
        keys are kept as ``fem_5mer_pc*`` for backwards compatibility
        with callers/tests, but the implementation is honest about what
        it computes.

        Algorithm
        ---------
        1. Convert raw counts to frequencies ``freqs``.
        2. Subtract the uniform-background null ``bg = 1/n_motifs``.
        3. Sort motif indices by count (descending); take the top N.
        4. For each of the top-N motifs, emit the per-motif deviation
           ``(freqs[i] - bg[i]) / norm * sqrt(n_motifs)`` where
           ``norm = ||freqs - bg||``.

        For a true cohort-level PCA reduction (motif × sample matrix
        → first N principal components per sample) use a separate
        external pipeline; this function is the per-sample deviation
        summary used by ``extract()``.

        Parameters
        ----------
        counts : np.ndarray, shape (n_motifs,)
            Per-sample motif counts (1024 for 5-mers, 256 for 4-mers).
        n_components : int
            Number of top-motif deviation scores to return.

        Returns
        -------
        np.ndarray, shape (n_components,)
            Top-N per-motif deviation scores. Zero-filled when counts
            are empty or all-zero.
        """
        total = counts.sum()
        if total == 0:
            return np.zeros(n_components)

        freqs = counts.astype(np.float64) / total

        # Uniform-background null. This is an OK null for 5-mers in
        # cfDNA (Jiang 2020 reports near-uniform background over
        # non-cancer controls) but a more accurate null would condition
        # on the local GC content — left for a future revision.
        bg = np.ones_like(freqs) / len(freqs)
        diff = freqs - bg
        norm = np.linalg.norm(diff)
        if norm == 0:
            return np.zeros(n_components)

        # Sort by count descending; emit the top-N per-motif deviations.
        # We sort by count (not by |deviation|) so the ordering is
        # stable across samples with similar but not identical motifs.
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

        For each length bin, the Simpson diversity is normalized
        against the *effective* motif alphabet in that bin (the number
        of motifs that received at least one count), not the universe
        size of 1024. The previous version used the universe
        denominator for every bin, which inflated MDS for length bins
        with sparse motif coverage (e.g. short fragments <150 bp
        typically only exercise a subset of all 1024 5-mers because
        DNASE1L3 cutting preferences are biased). Normalizing against
        the effective alphabet makes MDS values comparable across bins
        regardless of coverage sparsity.

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
            simpson = float(np.sum(p ** 2))
            # Effective alphabet: count motifs with at least one
            # observation in this bin. For sparse bins this can be
            # <<1024 and the normalizer (1 - 1/n_eff) prevents MDS
            # from being artificially inflated by zero-count motifs
            # that could never contribute to the numerator.
            n_eff = int((c > 0).sum())
            if n_eff <= 1:
                return 0.0
            return (1.0 - simpson) / (1.0 - 1.0 / n_eff)

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

        # Per-motif deviation scores (formerly mis-named 'PCA'; see
        # ``_top_motif_deviations`` for the honest algorithm description).
        dev_scores = self._top_motif_deviations(counts_5mer, n_components=10)
        for i in range(10):
            feats[f'fem_5mer_pc{i}'] = float(dev_scores[i])

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
            Optional per-fragment methylation array (bool or 0/1).
            Length must match ``fragments`` when provided. If set and
            ``fragments`` dicts lack a ``methylated`` field, this
            array is attached to each fragment so MFS-style features
            get a real methylation signal. When ``methylation_data``
            is None and ``fragments`` has no ``methylated`` field,
            MFS methylation-dependent features return zero (the
            documented fallback).
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

        # Backfill per-fragment methylation status from the standalone
        # ``methylation_data`` array if the fragments don't already
        # carry a ``methylated`` field. This is the path real cfDNA
        # methylation pipelines take (e.g. FinaleMe outputs a
        # per-fragment β-value thresholded at 0.5) and previously the
        # array was silently ignored.
        if (
            fragments is not None
            and methylation_data is not None
            and len(methylation_data) == len(fragments)
        ):
            for frag, meth in zip(fragments, methylation_data):
                if "methylated" not in frag:
                    frag["methylated"] = bool(meth)

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
