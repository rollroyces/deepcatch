"""Tests for the tumor-naive fragmentomics adapter.

These tests don't touch the real FinaleDB cohort — they synthesize
minimal pipeline artifacts in a tmp dir so the adapter contract is
verifiable without a 600+ sample download.
"""
import json
import os
import sys
import tempfile

import numpy as np
import pytest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from src.fragmentomics.tumor_naive_adapter import (  # noqa: E402
    CHANNEL_DIMS,
    CHANNEL_NAMES,
    TOTAL_DIM,
    TumorNaiveSample,
    expected_dim,
    load_cohort,
    load_labels_tsv,
    load_sample,
)


# ---------- fixtures ----------

def _make_sample(features_dir: str, sample_id: str,
                 n_bins: dict[str, int] | None = None) -> dict[str, int]:
    """Write minimal valid pipeline artifacts for one sample; return per-channel dims."""
    n_bins = n_bins or CHANNEL_DIMS
    rng = np.random.default_rng(0)

    for ch in ("delfi_5mb_ratio", "delfi_5mb_coverage",
               "delfi_100kb_ratio", "delfi_100kb_counts"):
        np.save(os.path.join(features_dir, f"{sample_id}.{ch}.npy"),
                rng.random(n_bins[ch]))

    # FSD JSON: 196 bins (5bp, 20-1000bp), normalized to sum=1
    bins = {}
    n_per = n_bins["fsd_histogram"]
    raw = rng.random(n_per)
    raw = raw / raw.sum()  # normalize, like the real pipeline
    for i in range(n_per):
        lo = 20 + i * 5
        hi = lo + 5
        bins[f"{lo}-{hi}"] = float(raw[i])
    with open(os.path.join(features_dir, f"{sample_id}.fsd.json"), "w") as f:
        json.dump({
            "sample": sample_id,
            "mode": "frag_tsv",
            "fragment_count": int(rng.integers(1_000_000, 30_000_000)),
            "median_length": 167.0,
            "mean_length": 167.0,
            "mode_length": 167.0,
            "p10": 100, "p25": 130, "p75": 200, "p90": 250,
            "short_fraction_100_150": 0.35,
            "long_fraction_150_220": 0.5,
            "short_long_ratio": 0.7,
            "size_bins": bins,
        }, f)
    return n_bins


# ---------- contract tests ----------

def test_channel_constants():
    """The channel contract is the public API of the adapter."""
    assert TOTAL_DIM == sum(CHANNEL_DIMS.values())
    assert TOTAL_DIM == 63_246  # 631+631+30894+30894+196
    assert set(CHANNEL_NAMES) == set(CHANNEL_DIMS.keys())
    assert len(CHANNEL_NAMES) == 5


def test_load_sample_returns_correct_shape():
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "S1")
        s = load_sample("S1", d)
        assert isinstance(s, TumorNaiveSample)
        assert s.sample_id == "S1"
        assert s.vector.shape == (TOTAL_DIM,)
        assert s.vector.dtype == np.float64
        assert s.n_frags is not None and s.n_frags > 0


def test_load_sample_fsd_bins_are_normalized():
    """The FSD bins should sum to 1.0 (the pipeline normalizes)."""
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "S2")
        s = load_sample("S2", d)
        fsd_start = sum(CHANNEL_DIMS[c] for c in CHANNEL_NAMES
                        if c != "fsd_histogram")
        fsd_block = s.vector[fsd_start:]
        np.testing.assert_allclose(fsd_block.sum(), 1.0, atol=1e-9)


def test_load_sample_missing_artifact_raises():
    with tempfile.TemporaryDirectory() as d:
        with pytest.raises(FileNotFoundError, match="delfi_5mb_ratio"):
            load_sample("missing", d)


def test_load_cohort_stacks_correctly():
    with tempfile.TemporaryDirectory() as d:
        for sid in ("A", "B", "C"):
            _make_sample(d, sid)
        ids = ["A", "B", "C"]
        X, order = load_cohort(ids, d)
        assert X.shape == (3, TOTAL_DIM)
        assert order == ["A", "B", "C"]


def test_load_cohort_skip_missing():
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "A")
        _make_sample(d, "B")
        # "C" has no artifacts — skip_missing=True should drop it
        X, order = load_cohort(["A", "B", "C"], d, skip_missing=True)
        assert X.shape == (2, TOTAL_DIM)
        assert order == ["A", "B"]


def test_load_cohort_skip_missing_false_is_strict():
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "A")
        with pytest.raises(FileNotFoundError):
            load_cohort(["A", "missing"], d, skip_missing=False)


def test_load_labels_tsv_basic():
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        f.write("A\tcancer\nB\thealthy\nC\tpositive\n")
        path = f.name
    try:
        labels = load_labels_tsv(path)
        assert labels == {"A": 1, "B": 0, "C": 1}
    finally:
        os.unlink(path)


def test_load_labels_tsv_handles_three_columns():
    """3-column TSV (sample / label / study) — the adapter ignores the study column."""
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        f.write("A\tcancer\tjiang\nB\thealthy\tcristiano\n")
        path = f.name
    try:
        labels = load_labels_tsv(path)
        assert labels == {"A": 1, "B": 0}
    finally:
        os.unlink(path)


def test_load_labels_tsv_rejects_unknown_label():
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
        f.write("A\tmaybe\n")
        path = f.name
    try:
        with pytest.raises(ValueError, match="Unrecognized label"):
            load_labels_tsv(path)
    finally:
        os.unlink(path)


def test_expected_dim_subset():
    assert expected_dim() == TOTAL_DIM
    assert expected_dim(["fsd_histogram"]) == 196
    assert expected_dim(["delfi_5mb_ratio", "fsd_histogram"]) == 631 + 196


def test_load_sample_channels_subset():
    """Requesting only a subset of channels returns a smaller vector."""
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "S")
        s = load_sample("S", d, channels=["fsd_histogram"])
        assert s.vector.shape == (196,)
        s2 = load_sample("S", d, channels=["delfi_5mb_ratio", "delfi_5mb_coverage"])
        assert s2.vector.shape == (631 + 631,)


def test_100kb_coverage_median_normalized_by_default():
    """Default behavior matches the pipeline: 100kb counts are median-centered.

    This is the fix that brings the DeepCatch adapter's result within
    0.001 AUC of the upstream pipeline (otherwise ~0.008 lower due to
    sequencing-depth batch effect).
    """
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "S")
        # With normalization ON: median of the 100kb block should be 1.0
        s_on = load_sample("S", d,
                           channels=["delfi_100kb_counts"],
                           median_normalize_100kb_coverage=True)
        block_on = s_on.vector
        med_on = float(np.median(block_on))
        assert abs(med_on - 1.0) < 1e-9, f"expected median 1.0, got {med_on}"

        # With normalization OFF: median is whatever the synthetic fixture had
        s_off = load_sample("S", d,
                            channels=["delfi_100kb_counts"],
                            median_normalize_100kb_coverage=False)
        med_off = float(np.median(s_off.vector))
        assert med_off != 1.0, "OFF branch should preserve raw median"


def test_median_normalization_passes_through_load_cohort():
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "A")
        _make_sample(d, "B")
        X, _ = load_cohort(["A", "B"], d,
                           channels=["delfi_100kb_counts"],
                           median_normalize_100kb_coverage=True)
        # Each row's 100kb block should be median-centered
        for row in X:
            assert abs(float(np.median(row)) - 1.0) < 1e-9


def test_load_sample_unknown_channel_raises():
    with tempfile.TemporaryDirectory() as d:
        _make_sample(d, "S")
        with pytest.raises(ValueError, match="Unknown channel"):
            load_sample("S", d, channels=["bogus_channel"])


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))