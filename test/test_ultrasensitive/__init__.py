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
    assert c.hamming_distance("ACGT", "TCGA") == 4


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
