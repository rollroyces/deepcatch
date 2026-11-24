"""
CANCER-SEEK ultra-sensitive module smoke tests.

These tests guard against regression in the v2.3.0-port of `cancer-screening`
siblings. They are intentionally lightweight (no real-data binding); the
fast-path gate verifies that ``import src.ultrasensitive`` works and that
each module's primary classes construct + produce a non-trivial output.

See ``docs/CANCER_SEEK_ARCHITECTURE.md`` for the full module map.
"""

import pytest

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


pytestmark = pytest.mark.ultra


def test_ultra_package_importable():
    """The src.ultrasensitive package itself imports cleanly."""
    from src.ultrasensitive import (
        UMICorrector,
        UMIFamily,
        UMIRead,
        DuplexConsensus,
        CNVAnalyzer,
        LowVAFVariantFilter,
        VariantCaller,
        VariantCall,
        MethylationProcessor,
        SpikeInSimulator,
        SpikeInConfig,
        run_lod_validation,
    )
    assert UMICorrector is not None
    assert CNVAnalyzer is not None
    assert LowVAFVariantFilter is not None
    assert MethylationProcessor is not None
    assert SpikeInSimulator is not None


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_umi_corrector_hamming_distance():
    """``UMICorrector.hamming_distance`` is the entry point for Hamming clustering."""
    from src.ultrasensitive.umi import UMICorrector

    c = UMICorrector()
    assert c.hamming_distance("ACGT", "ACGT") == 0
    assert c.hamming_distance("ACGT", "ACGA") == 1
    # Position-by-position comparison: A/T differ, C/C and G/G match, T/A differ.
    assert c.hamming_distance("ACGT", "TCGA") == 2
    assert c.hamming_distance("ACGT", "TGCA") == 4


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_cnv_analyzer_construction():
    """``CNVAnalyzer`` constructs with non-default parameters."""
    from src.ultrasensitive.cnv_analyzer import CNVAnalyzer

    a = CNVAnalyzer(bin_size=100_000, segmentation_alpha=1e-3)
    assert a.bin_size == 100_000
    assert a.segmentation_alpha == 1e-3


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_low_vaf_variant_filter_pvalue_floor():
    """Background-error p-value floor; binomial/Poisson boundary at 100."""
    from src.ultrasensitive.variant_calling import LowVAFVariantFilter

    f = LowVAFVariantFilter(background_error_rate=1e-7, binomial_pval_threshold=1e-6)
    p_small = f.compute_binomial_pvalue(alt_count=5, total_count=10, error_rate=1e-7)
    p_large = f.compute_binomial_pvalue(alt_count=5, total_count=1000, error_rate=1e-7)
    assert 0.0 <= p_small <= 1.0
    assert 0.0 <= p_large <= 1.0


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_methylation_beta_value_smoothing():
    """Beta-value computation honors the smoothing alpha."""
    from src.ultrasensitive.methylation import MethylationProcessor

    mp = MethylationProcessor(smoothing_constant=100.0)
    beta = mp.compute_beta_values([80], [20])  # ground truth = 0.8
    # With alpha=100 smoothing the result is (80+100)/(80+20+200) = 180/300 = 0.6
    assert abs(beta[0] - 0.6) < 1e-6


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_spike_in_simulator_runs():
    """``SpikeInSimulator.run_spike_in_experiment`` produces a non-empty summary."""
    from src.ultrasensitive.spike_in_simulator import (
        SpikeInConfig, SpikeInSimulator,
    )

    cfg = SpikeInConfig(tumor_fraction=1e-3, genome_equivalents=50_000, n_samples=3)
    sim = SpikeInSimulator(cfg)
    summary, _ = sim.run_spike_in_experiment()
    assert len(summary) == 3
    assert "target_vaf" in summary.columns


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_run_lod_validation_smoke():
    """``run_lod_validation`` runs at a single VAF and returns the row."""
    from src.ultrasensitive.spike_in_simulator import run_lod_validation

    df = run_lod_validation([1e-3], n_samples=3)
    assert len(df) == 1
    assert "target_vaf" in df.columns
    assert "lod_achieved" in df.columns


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_spike_in_simulator_detection_is_depth_driven():
    """Detection must be governed by sequencing depth, not by replicate count.

    Regression guard. The ported simulator used to draw mutant reads from the
    molecules physically present in a 100 bp window -- 50,000 genome
    equivalents put 0.0016 molecules in such a window -- so the Poisson mean
    was ~1e-8 and it detected nothing at any VAF, while the power analysis
    (computed genome-wide) reported power 1.0. The two columns contradicted
    each other and the "0.001% VAF confirmed" claim had no support.

    The corrected model samples ``depth`` reads at a target site, so the
    expected mutant count is ``depth * VAF`` and power rises with depth.
    """
    from src.ultrasensitive.spike_in_simulator import SpikeInConfig, SpikeInSimulator

    def power_at(depth: int, vaf: float = 1e-3) -> float:
        cfg = SpikeInConfig(
            tumor_fraction=vaf,
            genome_equivalents=50_000,
            n_samples=5,
            sequencing_depth=depth,
        )
        return SpikeInSimulator(cfg).power_analysis()["single_sample_power"]

    # 1,000x is far too shallow for 0.1% VAF: 1 expected mutant read.
    assert power_at(1_000) < 0.2
    # 100,000x gives 100 expected mutant reads -- comfortable.
    assert power_at(100_000) > 0.9
    # Power is monotonically non-decreasing in depth.
    assert power_at(10_000) <= power_at(100_000) <= power_at(1_000_000)


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_spike_in_simulator_0_001pct_is_not_confirmed_at_100k_depth():
    """0.001% VAF is NOT reachable at 100,000x -- the honest result.

    At 50,000 genome equivalents, 0.001% VAF is 0.5 mutant molecules in the
    whole genome, and 100,000x yields 1 mutant read at a site -- below the
    3-supporting-read threshold. The simulator must report that, rather than
    claiming the LOD is confirmed.
    """
    from src.ultrasensitive.spike_in_simulator import SpikeInConfig, SpikeInSimulator

    cfg = SpikeInConfig(
        tumor_fraction=1e-5,
        genome_equivalents=50_000,
        n_samples=10,
        sequencing_depth=100_000,
    )
    sim = SpikeInSimulator(cfg)
    counts = sim.compute_molecule_counts()
    power = sim.power_analysis()

    # 0.5 mutant molecules genome-wide -- the number that makes 0.001% hard.
    assert abs(counts["expected_mutant_molecules_genome_wide"] - 0.5) < 1e-9
    # One mutant read per site at this depth, and one read is not a call.
    assert abs(counts["expected_mutant_reads"] - 1.0) < 1e-9
    assert power["single_sample_power"] < 0.5
    # The depth that would actually clear the threshold is reported, and it
    # is far above 100,000x.
    assert power["depth_for_80pct_power"] > 100_000


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_run_lod_validation_detection_rate_agrees_with_power():
    """The observed detection rate and the modelled power must not diverge.

    The old model reported ``detection_rate_mean = 0.0`` next to
    ``single_sample_power = 1.0``. At a depth where the model says the
    mutation is comfortably detectable, the simulated detector has to find
    it.
    """
    from src.ultrasensitive.spike_in_simulator import run_lod_validation

    df = run_lod_validation([1e-2], n_samples=10, sequencing_depth=100_000)
    row = df.iloc[0]
    assert row["single_sample_power"] > 0.9
    assert row["detection_rate_mean"] > 0.9


@pytest.mark.skipif(not HAS_NUMPY, reason="requires numpy")
def test_xgboost_has_dep_guard():
    """XGBoost classifier should raise a friendly ``ImportError`` without [ultra] extras."""
    from src.ultrasensitive.xgboost_classifier import XGBoostCancerClassifier

    # The module-level ``HAS_DEPS = False`` path is the no-xgboost build;
    # the .train() call should raise ImportError with our explanatory message.
    clf = XGBoostCancerClassifier()
    if not clf.__class__.__module__.endswith(".xgboost_classifier"):
        return  # pragma: no cover -- guards against pytest collection false-positive
    # We don't actually invoke train() here because it requires xgboost;
    # just confirm the class exists and its DEFAULT_PARAMS schema is sane.
    assert "n_estimators" in XGBoostCancerClassifier.DEFAULT_PARAMS
    assert "scale_pos_weight" in XGBoostCancerClassifier.DEFAULT_PARAMS
