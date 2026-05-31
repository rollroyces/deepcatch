#!/usr/bin/env python3
"""
Full test suite for Enhanced Fragmentomics Feature Extractors.

Tests DELFIFeatures, MFSFeatures, NucleosomeFootprint, RefinedEndMotifs,
and EnhancedFragmentomics (unified extractor).

Coverage: 30+ tests covering normal operation, edge cases, missing data
fallbacks, deterministic reproducibility, and numerical range sanity.
"""

import sys
import os
import unittest
import numpy as np

# ── Path setup ────────────────────────────────────────────────────────────────
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)
)))
sys.path.insert(0, PROJECT_ROOT)

from src.fragmentomics.enhanced_features import (
    DELFIFeatures,
    MFSFeatures,
    NucleosomeFootprint,
    RefinedEndMotifs,
    EnhancedFragmentomics,
)

# ── Test Configuration ────────────────────────────────────────────────────────
TEST_SEED = 42
N_FRAGMENTS = 500
GENOME_LENGTH = 3_000_000_000
WINDOW_SIZE = 100_000
MFS_BIN_SIZE = 1_000_000

rng = np.random.RandomState(TEST_SEED)

BASES = ['A', 'C', 'G', 'T']


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_fragment_lengths(n=N_FRAGMENTS, seed=TEST_SEED) -> np.ndarray:
    """Generate realistic fragment length distribution (~167 bp mean)."""
    loc_rng = np.random.RandomState(seed)
    lengths = loc_rng.gamma(shape=12, scale=14, size=n)
    return np.clip(lengths, 20, 600).astype(np.float64)


def _make_fragments(n=N_FRAGMENTS, seed=TEST_SEED):
    """Generate mock fragment dicts with start, length, methylated."""
    loc_rng = np.random.RandomState(seed)
    fragments = []
    for _ in range(n):
        start = loc_rng.randint(0, GENOME_LENGTH - 500)
        length = int(loc_rng.gamma(shape=12, scale=14))
        length = max(20, min(600, length))
        fragments.append({
            'start': start,
            'length': length,
            'methylated': loc_rng.random() < 0.3,
        })
    return fragments


def _make_end_sequences(n=N_FRAGMENTS, seed=TEST_SEED):
    """Generate random 5'-end sequences (10 bp each)."""
    loc_rng = np.random.RandomState(seed)
    seqs = []
    for _ in range(n):
        seq = ''.join(BASES[loc_rng.randint(0, 4)] for _ in range(10))
        seqs.append(seq)
    return seqs


def _make_tss_positions(n=50, seed=TEST_SEED):
    """Generate mock TSS positions spread across the genome."""
    loc_rng = np.random.RandomState(seed)
    return list(loc_rng.randint(0, GENOME_LENGTH, size=n))


# ═══════════════════════════════════════════════════════════════════════════════
# A. DELFIFeatures Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestDELFIFeatures(unittest.TestCase):
    """Tests for DELFI-style genome-wide fragment coverage features."""

    def setUp(self):
        self.delfi = DELFIFeatures(window_size=100_000, size_bin_width=5)
        self.lengths = _make_fragment_lengths()
        self.fragments = _make_fragments()

    def test_extract_all_keys_present(self):
        """extract() returns all 15 expected feature keys."""
        result = self.delfi.extract(self.lengths, self.fragments)
        expected_keys = [
            'delfi_coverage_cv', 'delfi_abnormal_window_frac',
            'delfi_coverage_autocorr', 'delfi_coverage_mean',
            'delfi_coverage_median', 'delfi_coverage_iqr',
            'delfi_n_windows', 'delfi_size_entropy',
            'delfi_size_mode_bin', 'delfi_size_hist_skew',
            'delfi_peak_count', 'delfi_peak_spacing_mean',
            'delfi_peak_spacing_std', 'delfi_first_peak_bp',
            'delfi_peak_height_cv',
        ]
        for key in expected_keys:
            self.assertIn(key, result, f"Missing key: {key}")
        self.assertEqual(len(result), 15, "Should have exactly 15 keys")

    def test_extract_empty_array(self):
        """extract() with empty array returns most features as zeros gracefully."""
        result = self.delfi.extract(np.array([]))
        # Most features should be zero for empty input
        zero_keys = [
            'delfi_coverage_cv', 'delfi_abnormal_window_frac',
            'delfi_coverage_autocorr', 'delfi_coverage_mean',
            'delfi_coverage_median', 'delfi_coverage_iqr',
            'delfi_size_entropy', 'delfi_size_mode_bin',
            'delfi_size_hist_skew', 'delfi_peak_count',
            'delfi_peak_spacing_mean', 'delfi_peak_spacing_std',
            'delfi_first_peak_bp', 'delfi_peak_height_cv',
        ]
        for key in zero_keys:
            self.assertEqual(result[key], 0.0,
                             f"Key {key} should be 0.0 for empty input, got {result[key]}")
        # delfi_n_windows is structural (genome_length / window_size), not data-dependent
        self.assertEqual(result['delfi_n_windows'], float(GENOME_LENGTH // 100_000))

    def test_extract_all_identical_lengths(self):
        """extract() with all identical fragment lengths is an edge case – no crash."""
        identical = np.full(100, 167.0)
        result = self.delfi.extract(identical)
        self.assertGreaterEqual(len(result), 15)
        # Histogram skew should be near-zero for uniform distribution
        self.assertIsInstance(result['delfi_size_hist_skew'], float)

    def test_coverage_profile_basic(self):
        """coverage_profile() returns correct window counts."""
        frags = [
            {'start': 50_000},
            {'start': 150_000},
            {'start': 50_000},  # Same window as first
            {'start': 250_000},
        ]
        coverage = self.delfi.coverage_profile(frags, genome_length=500_000)
        self.assertEqual(len(coverage), 5)  # 500k / 100k = 5 windows
        self.assertEqual(coverage[0], 2.0)  # Two fragments in window 0
        self.assertEqual(coverage[1], 1.0)
        self.assertEqual(coverage[2], 1.0)

    def test_coverage_profile_none_fragments(self):
        """coverage_profile() with None fragments returns zeros."""
        coverage = self.delfi.coverage_profile(None)
        n_windows = GENOME_LENGTH // 100_000
        self.assertEqual(len(coverage), n_windows)
        self.assertTrue(np.all(coverage == 0.0))

    def test_size_distribution_sums_correctly(self):
        """size_distribution_5bp() histogram sums to total fragment count."""
        hist = self.delfi.size_distribution_5bp(self.lengths, max_size=600)
        self.assertAlmostEqual(hist.sum(), len(self.lengths), delta=1)

    def test_coverage_cv_positive(self):
        """coverage_cv() is non-negative."""
        coverage = self.delfi.coverage_profile(self.fragments)
        cv = self.delfi.coverage_cv(coverage)
        self.assertGreaterEqual(cv, 0.0)

    def test_autocorrelation_with_few_windows(self):
        """coverage_autocorrelation() handles < 2 windows."""
        ac = self.delfi.coverage_autocorrelation(np.array([1.0]))
        self.assertEqual(ac, 0.0)

    def test_abnormal_window_fraction_edge(self):
        """abnormal_window_fraction() with constant coverage returns 0."""
        const = np.ones(100)
        frac = self.delfi.abnormal_window_fraction(const, z_threshold=2.0)
        self.assertEqual(frac, 0.0)

    def test_extract_without_fragments_still_has_size_features(self):
        """extract() without fragments dict still computes size-distribution features."""
        result = self.delfi.extract(self.lengths, fragments=None)
        # Size features should be computed from lengths always
        self.assertGreater(result['delfi_size_entropy'], 0.0,
                           "Size entropy should be >0 for real distribution")
        self.assertTrue(np.isfinite(result['delfi_size_entropy']))


# ═══════════════════════════════════════════════════════════════════════════════
# B. MFSFeatures Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestMFSFeatures(unittest.TestCase):
    """Tests for Methylation-Fragment-Size merged features."""

    def setUp(self):
        self.mfs = MFSFeatures(bin_size=1_000_000)
        self.fragments = _make_fragments()

    def test_extract_all_keys_present(self):
        """extract() returns all 8 expected feature keys."""
        result = self.mfs.extract(self.fragments)
        expected_keys = [
            'mfs_size_meth_corr', 'mfs_frag_entropy', 'mfs_bin_dispersion',
            'mfs_n_bins', 'mfs_meth_short', 'mfs_meth_long', 'mfs_meth_ratio',
            'mfs_mean_bin_count',
        ]
        for key in expected_keys:
            self.assertIn(key, result, f"Missing key: {key}")
        self.assertEqual(len(result), len(expected_keys))

    def test_extract_none_fragments_returns_zeros(self):
        """extract() with None returns all zero fallback."""
        result = self.mfs.extract(None)
        for key, value in result.items():
            self.assertEqual(value, 0.0,
                             f"Key {key} should be 0.0 for None input, got {value}")

    def test_extract_empty_fragments_returns_zeros(self):
        """extract() with empty list returns all zero fallback."""
        result = self.mfs.extract([])
        for key, value in result.items():
            self.assertEqual(value, 0.0,
                             f"Key {key} should be 0.0 for empty input, got {value}")

    def test_constant_methylation_yields_zero_corr(self):
        """size_methylation_correlation() with constant methylation returns 0."""
        frags = []
        loc_rng = np.random.RandomState(TEST_SEED)
        for _ in range(200):
            start = loc_rng.randint(0, GENOME_LENGTH)
            length = int(loc_rng.gamma(shape=12, scale=14))
            frags.append({
                'start': start,
                'length': max(20, min(600, length)),
                'methylated': True,  # All methylated = constant
            })
        result = self.mfs.extract(frags)
        # With all methylated, there's no variance in meth_fracs across bins
        self.assertAlmostEqual(result['mfs_size_meth_corr'], 0.0, delta=0.05)

    def test_perfect_correlation_size_and_meth(self):
        """When methylation perfectly correlates with size, corr should be high."""
        frags = []
        for i in range(300):
            # Spread across bins
            start = (i % 100) * 1_000_000
            length = 100 + i  # Increases linearly
            methylated = i > 150  # Threshold at 150
            frags.append({
                'start': start,
                'length': length,
                'methylated': methylated,
            })
        result = self.mfs.extract(frags)
        # Should have strong positive correlation
        self.assertGreater(abs(result['mfs_size_meth_corr']), 0.3)

    def test_single_fragment_edge_case(self):
        """extract() with a single fragment works (not enough for correlation)."""
        single = [{'start': 1_000_000, 'length': 167, 'methylated': False}]
        result = self.mfs.extract(single)
        self.assertIn('mfs_size_meth_corr', result)
        self.assertEqual(result['mfs_size_meth_corr'], 0.0)

    def test_size_specific_methylation_values(self):
        """size_specific_methylation() returns fractions in [0,1]."""
        feats = self.mfs.size_specific_methylation(self.fragments)
        self.assertGreaterEqual(feats['mfs_meth_short'], 0.0)
        self.assertLessEqual(feats['mfs_meth_short'], 1.0)
        self.assertGreaterEqual(feats['mfs_meth_long'], 0.0)
        self.assertLessEqual(feats['mfs_meth_long'], 1.0)

    def test_bin_dispersion_non_negative(self):
        """bin_dispersion() is non-negative."""
        frag_counts, _, _ = self.mfs._build_joint_histogram(self.fragments)
        disp = self.mfs.bin_dispersion(frag_counts)
        self.assertGreaterEqual(disp, 0.0)


# ═══════════════════════════════════════════════════════════════════════════════
# C. NucleosomeFootprint Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestNucleosomeFootprint(unittest.TestCase):
    """Tests for nucleosome footprint and TSS-proximal coverage features."""

    def setUp(self):
        self.nuc = NucleosomeFootprint(tss_window=2000, bin_size=50)
        self.fragments = _make_fragments()
        self.tss_positions = _make_tss_positions()

    def test_extract_all_keys_present(self):
        """extract() returns all 5 expected feature keys."""
        result = self.nuc.extract(self.fragments, self.tss_positions)
        expected_keys = [
            'nuc_occupancy_score', 'nuc_tss_depletion_ratio',
            'nuc_dominant_freq', 'nuc_dominant_power', 'nuc_period_bp',
        ]
        for key in expected_keys:
            self.assertIn(key, result, f"Missing key: {key}")
        self.assertEqual(len(result), len(expected_keys))

    def test_extract_without_tss_returns_zeros(self):
        """extract() without tss_positions returns all zeros."""
        result = self.nuc.extract(self.fragments, tss_positions=None)
        for key, value in result.items():
            self.assertEqual(value, 0.0,
                             f"Key {key} should be 0.0 without TSS, got {value}")

    def test_extract_with_none_fragments_returns_zeros(self):
        """extract() with None fragments returns all zeros."""
        result = self.nuc.extract(None, self.tss_positions)
        for key, value in result.items():
            self.assertEqual(value, 0.0,
                             f"Key {key} should be 0.0 without fragments, got {value}")

    def test_tss_coverage_profile_shape(self):
        """tss_coverage_profile() returns array of correct shape."""
        profile = self.nuc.tss_coverage_profile(self.tss_positions, self.fragments, n_bins=80)
        self.assertEqual(len(profile), 80)
        self.assertTrue(np.all(profile >= 0))

    def test_periodicity_known_signal(self):
        """coverage_periodicity() detects known periodic signal correctly."""
        # Generate a clean sinusoid at ~195 bp period (3.9 bins × 50 bp)
        n = 80
        x = np.arange(n)
        period_bins = 195.0 / 50.0  # ≈ 3.9
        signal = 10.0 + 3.0 * np.cos(2 * np.pi * x / period_bins)
        result = self.nuc.coverage_periodicity(signal)
        # Dominant period should be close to 195 bp
        self.assertGreater(result['nuc_dominant_power'], 0.0,
                           "Should detect non-zero dominant power for periodic signal")
        # The period in bp should be close to 195 (within reasonable tolerance)
        self.assertGreater(result['nuc_period_bp'], 100.0)
        self.assertLess(result['nuc_period_bp'], 300.0)

    def test_expected_nucleosome_pattern_normalized(self):
        """expected_nucleosome_pattern() has mean ≈ 1."""
        pattern = self.nuc.expected_nucleosome_pattern(n_bins=80)
        self.assertAlmostEqual(pattern.mean(), 1.0, delta=0.01)

    def test_occupancy_score_range(self):
        """nucleosome_occupancy_score() is in [-1, 1] (correlation)."""
        profile = self.nuc.tss_coverage_profile(self.tss_positions, self.fragments, n_bins=80)
        score = self.nuc.nucleosome_occupancy_score(profile)
        self.assertGreaterEqual(score, -1.1)
        self.assertLessEqual(score, 1.1)

    def test_depletion_ratio_non_negative(self):
        """tss_depletion_ratio() is non-negative."""
        profile = self.nuc.tss_coverage_profile(self.tss_positions, self.fragments, n_bins=80)
        ratio = self.nuc.tss_depletion_ratio(profile, n_bins=80)
        self.assertGreaterEqual(ratio, 0.0)


# ═══════════════════════════════════════════════════════════════════════════════
# D. RefinedEndMotifs Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestRefinedEndMotifs(unittest.TestCase):
    """Tests for refined fragment end motif features."""

    def setUp(self):
        self.motif = RefinedEndMotifs()
        self.end_seqs = _make_end_sequences()
        self.lengths = _make_fragment_lengths()

    def test_extract_all_keys_present(self):
        """extract() returns all 16 expected feature keys."""
        result = self.motif.extract(self.end_seqs, self.lengths)
        expected_keys = [
            'fem_5mer_pc0', 'fem_5mer_pc1', 'fem_5mer_pc2',
            'fem_5mer_pc3', 'fem_5mer_pc4', 'fem_5mer_pc5',
            'fem_5mer_pc6', 'fem_5mer_pc7', 'fem_5mer_pc8', 'fem_5mer_pc9',
            'fem_mds_short', 'fem_mds_mid', 'fem_mds_long',
            'fem_gc_bias_score', 'fem_gc_corrected_mds',
            'fem_5mer_entropy',
        ]
        for key in expected_keys:
            self.assertIn(key, result, f"Missing key: {key}")
        self.assertEqual(len(result), len(expected_keys))

    def test_extract_without_end_sequences_returns_zeros(self):
        """extract() with None end_sequences returns all zeros."""
        result = self.motif.extract(None)
        for key, value in result.items():
            self.assertEqual(value, 0.0,
                             f"Key {key} should be 0.0 without end_sequences, got {value}")

    def test_extract_empty_sequences_returns_zeros(self):
        """extract() with empty sequence list returns all zeros."""
        result = self.motif.extract([])
        for key, value in result.items():
            self.assertEqual(value, 0.0,
                             f"Key {key} should be 0.0 for empty input, got {value}")

    def test_uniform_motif_distribution(self):
        """Uniform motif distribution produces near-zero bias score."""
        # Generate sequences that are uniformly distributed
        n = 4096  # 4× per motif to get reasonable coverage
        loc_rng = np.random.RandomState(TEST_SEED)
        all_5mers = [a+b+c+d+e for a in BASES for b in BASES
                     for c in BASES for d in BASES for e in BASES]
        seqs = [m + 'XXXXX' for m in all_5mers * 4]
        result = self.motif.extract(seqs)
        # With uniform distribution, GC bias should be near 0
        self.assertAlmostEqual(result['fem_gc_bias_score'], 0.0, delta=0.01)

    def test_gc_bias_correction_biased_case(self):
        """GC bias correction: known GC-rich bias produces positive score."""
        # Generate only GC-rich motifs (all G/C)
        n = 500
        loc_rng = np.random.RandomState(TEST_SEED)
        gc_seqs = []
        for _ in range(n):
            seq = ''.join(loc_rng.choice(['G', 'C'], size=5)) + 'XXXXX'
            gc_seqs.append(seq)
        # Also add some AT-rich to avoid zero counts causing NaN
        for _ in range(n // 2):
            seq = ''.join(loc_rng.choice(['A', 'T'], size=5)) + 'XXXXX'
            gc_seqs.append(seq)
        result = self.motif.extract(gc_seqs)
        # GC bias should be positive (enrichment of high-GC motifs)
        self.assertGreater(result['fem_gc_bias_score'], 0.0,
                           "GC-biased input should produce positive GC bias score")

    def test_length_bin_mds_values(self):
        """Motif diversity scores by length bin are in [0, 1]."""
        result = self.motif.extract(self.end_seqs, self.lengths)
        for key in ['fem_mds_short', 'fem_mds_mid', 'fem_mds_long']:
            self.assertGreaterEqual(result[key], 0.0,
                                    f"{key} should be >= 0")
            self.assertLessEqual(result[key], 1.0,
                                 f"{key} should be <= 1")

    def test_single_fragment_edge(self):
        """extract() with single fragment does not crash."""
        single_seq = ['ACGTACGTAC']
        result = self.motif.extract(single_seq)
        self.assertIn('fem_5mer_entropy', result)
        self.assertIsInstance(result['fem_5mer_pc0'], float)

    def test_5mer_entropy_positive_for_varied_motifs(self):
        """5-mer entropy is positive for varied motif distribution."""
        result = self.motif.extract(self.end_seqs)
        self.assertGreater(result['fem_5mer_entropy'], 0.0,
                           "Entropy should be >0 for varied motifs")


# ═══════════════════════════════════════════════════════════════════════════════
# E. EnhancedFragmentomics Integration Tests
# ═══════════════════════════════════════════════════════════════════════════════

class TestEnhancedFragmentomics(unittest.TestCase):
    """Integration tests for the unified EnhancedFragmentomics extractor."""

    def setUp(self):
        self.enhanced = EnhancedFragmentomics(
            window_size=100_000,
            mfs_bin_size=1_000_000,
            tss_window=2000,
        )
        self.lengths = _make_fragment_lengths()
        self.fragments = _make_fragments()
        self.end_seqs = _make_end_sequences()
        self.tss_positions = _make_tss_positions()

    def test_extract_all_with_only_lengths(self):
        """extract_all() with only fragment_lengths returns features with expected keys."""
        result = self.enhanced.extract_all(self.lengths)
        # Should have DELFI features (size-based ones at least)
        self.assertIn('delfi_size_entropy', result)
        self.assertIn('delfi_coverage_cv', result)
        # MFS should be zeros (no fragments dict)
        self.assertEqual(result['mfs_size_meth_corr'], 0.0)
        # Motif should be zeros (no end sequences)
        self.assertEqual(result['fem_5mer_entropy'], 0.0)
        # Nucleosome should be zeros (no TSS)
        self.assertEqual(result['nuc_occupancy_score'], 0.0)

    def test_extract_all_with_full_data(self):
        """extract_all() with all data returns non-zero features where expected."""
        result = self.enhanced.extract_all(
            fragment_lengths=self.lengths,
            fragments=self.fragments,
            end_sequences=self.end_seqs,
            tss_positions=self.tss_positions,
        )
        # Features that should be non-zero with real data
        self.assertGreater(result['delfi_size_entropy'], 0.0)
        self.assertGreater(result['delfi_n_windows'], 0.0)
        self.assertGreater(result['fem_5mer_entropy'], 0.0)
        # With methylation data present, correlation may be non-zero
        self.assertIsInstance(result['mfs_size_meth_corr'], float)

    def test_extract_all_empty_array_no_crash(self):
        """extract_all() with empty array returns all zeros, no crash."""
        result = self.enhanced.extract_all(np.array([]))
        self.assertIsInstance(result, dict)
        # Check a few representative keys
        for key in ['delfi_size_entropy', 'mfs_size_meth_corr',
                     'nuc_occupancy_score', 'fem_5mer_entropy']:
            self.assertIn(key, result)
            self.assertEqual(result[key], 0.0,
                             f"{key} should be 0.0 for empty input")

    def test_extract_all_deterministic(self):
        """extract_all() is deterministic: same input → same output."""
        kwargs = dict(
            fragment_lengths=self.lengths,
            fragments=self.fragments,
            end_sequences=self.end_seqs,
            tss_positions=self.tss_positions,
        )
        result1 = self.enhanced.extract_all(**kwargs)
        result2 = self.enhanced.extract_all(**kwargs)
        for key in result1:
            self.assertEqual(result1[key], result2[key],
                             f"Key {key} differs between runs")

    def test_feature_names_union(self):
        """All feature keys from sub-extractors are present in extract_all()."""
        result = self.enhanced.extract_all(
            fragment_lengths=self.lengths,
            fragments=self.fragments,
            end_sequences=self.end_seqs,
            tss_positions=self.tss_positions,
        )
        # DELFI keys
        self.assertIn('delfi_coverage_cv', result)
        self.assertIn('delfi_size_entropy', result)
        # MFS keys
        self.assertIn('mfs_size_meth_corr', result)
        self.assertIn('mfs_meth_ratio', result)
        # Nucleosome keys
        self.assertIn('nuc_occupancy_score', result)
        self.assertIn('nuc_tss_depletion_ratio', result)
        # Motif keys
        self.assertIn('fem_5mer_pc0', result)
        self.assertIn('fem_gc_bias_score', result)

    def test_output_values_in_reasonable_ranges(self):
        """extract_all() output values are within reasonable ranges."""
        result = self.enhanced.extract_all(
            fragment_lengths=self.lengths,
            fragments=self.fragments,
            end_sequences=self.end_seqs,
            tss_positions=self.tss_positions,
        )
        # Fractions should be in [0, 1]
        for key in ['delfi_abnormal_window_frac', 'mfs_meth_short', 'mfs_meth_long']:
            self.assertGreaterEqual(result[key], 0.0,
                                    f"{key}={result[key]} should be >= 0")
        # Entropies should be ≥ 0
        for key in ['delfi_size_entropy', 'mfs_frag_entropy', 'fem_5mer_entropy']:
            self.assertGreaterEqual(result[key], 0.0,
                                    f"{key}={result[key]} should be >= 0")
        # Occupancy score (correlation) ∈ [-1, 1]
        self.assertGreaterEqual(result['nuc_occupancy_score'], -1.1)
        self.assertLessEqual(result['nuc_occupancy_score'], 1.1)
        # GC bias score should be finite
        self.assertTrue(np.isfinite(result['fem_gc_bias_score']),
                        f"fem_gc_bias_score={result['fem_gc_bias_score']} should be finite")
        # MDS values ∈ [0, 1]
        for key in ['fem_mds_short', 'fem_mds_mid', 'fem_mds_long', 'fem_gc_corrected_mds']:
            self.assertGreaterEqual(result[key], 0.0,
                                    f"{key}={result[key]} should be >= 0")
            self.assertLessEqual(result[key], 1.0,
                                 f"{key}={result[key]} should be <= 1")

    def test_total_feature_count(self):
        """extract_all() total feature count matches sum of sub-extractors."""
        result = self.enhanced.extract_all(
            fragment_lengths=self.lengths,
            fragments=self.fragments,
            end_sequences=self.end_seqs,
            tss_positions=self.tss_positions,
        )
        total = len(result)
        # DELFI: 15 + MFS: 8 + Nucleosome: 5 + Motif: 16 = 44
        self.assertEqual(total, 44,
                         f"Expected 44 features, got {total}")

    def test_extract_all_none_lengths_does_not_crash(self):
        """extract_all() with fragment_lengths=None falls back gracefully."""
        result = self.enhanced.extract_all(None)
        self.assertIsInstance(result, dict)
        self.assertGreater(len(result), 0)


# ═══════════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # Run tests with verbosity
    suite = unittest.TestLoader().loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Print summary
    print("\n" + "=" * 60)
    print("  Enhanced Fragmentomics Test Suite Summary")
    print("=" * 60)
    print(f"  Tests run:    {result.testsRun}")
    print(f"  Failures:     {len(result.failures)}")
    print(f"  Errors:       {len(result.errors)}")
    print(f"  Skipped:      {len(result.skipped)}")
    if result.wasSuccessful():
        print("  Status:       ✅ ALL TESTS PASSED")
    else:
        print("  Status:       ❌ SOME TESTS FAILED")
        for test, traceback in result.failures + result.errors:
            print(f"\n  FAIL: {test}")
            print(f"  {traceback[:200]}...")
    print("=" * 60)
    sys.exit(0 if result.wasSuccessful() else 1)
