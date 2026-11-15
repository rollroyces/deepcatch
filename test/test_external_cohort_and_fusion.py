"""Tests for the external-cohort feature extractor and the fusion estimators.

Two modules are covered here:

- ``src.fragmentomics.finale_features`` -- the self-contained FinaleDB
  fragment-file extractor. The important property is that it accepts
  BOTH FinaleDB fragment schemas (5-column ``chrom start end mapq
  strand`` and 6-column bedtools ``chrom start end name mapq strand``)
  and produces the same DELFI/FSD features as the sibling pipeline's
  extractor, which the benchmark's cached training features came from.
- ``src.foundation.regularized_fusion`` -- the late-fusion estimators.
  The property that matters is that they are honest: they must not
  manufacture signal where there is none, and the BMA weights must be a
  proper distribution.
"""
from __future__ import annotations

import gzip
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.foundation.regularized_fusion import (  # noqa: E402
    bayesian_model_average,
    log_loss,
    regularized_late_fusion,
    zscore,
)
from src.fragmentomics.finale_features import (  # noqa: E402
    CHROM_SIZES,
    FSD_BIN_START,
    N_BINS,
    extract_delfi,
    extract_fsd,
)


# --------------------------------------------------------------------------- #
# Feature extraction
# --------------------------------------------------------------------------- #
def _write_frags(path: Path, n_rows: int = 20000, six_columns: bool = False,
                 seed: int = 0) -> None:
    """Write a synthetic fragment file in one of the two schemas."""
    rng = np.random.default_rng(seed)
    chroms = list(CHROM_SIZES)
    with gzip.open(path, "wt") as f:
        for _ in range(n_rows):
            c = chroms[rng.integers(len(chroms))]
            size = CHROM_SIZES[c]
            start = int(rng.integers(0, max(size - 1000, 1)))
            length = int(rng.integers(100, 221))
            end = min(start + length, size)
            mapq = int(rng.integers(0, 60))
            if six_columns:
                f.write(f"{c}\t{start}\t{end}\tfrag1\t{mapq}\t"
                        f"{'+' if rng.random() > 0.5 else '-'}\n")
            else:
                f.write(f"{c}\t{start}\t{end}\t{mapq}\t"
                        f"{'+' if rng.random() > 0.5 else '-'}\n")


def test_both_frag_schemas_give_identical_features(tmp_path):
    """5-column and 6-column files must yield the same numbers.

    The sibling pipeline's FSD extractor required 6 columns and silently
    dropped every fragment of the 5-column files FinaleDB currently
    serves; this module has to accept both.
    """
    five = tmp_path / "five.frag.tsv.bgz"
    six = tmp_path / "six.frag.tsv.bgz"
    _write_frags(five, six_columns=False)
    _write_frags(six, six_columns=True)

    d5, d6 = extract_delfi(str(five)), extract_delfi(str(six))
    assert d5["5mb_ratio"].shape == (631,)
    assert d5["100kb_ratio"].shape == (30894,)
    assert np.allclose(d5["5mb_ratio"], d6["5mb_ratio"])
    assert np.allclose(d5["100kb_ratio"], d6["100kb_ratio"])
    assert np.allclose(d5["100kb_counts"], d6["100kb_counts"])
    assert d5["genome_short_long_ratio"] == pytest.approx(
        d6["genome_short_long_ratio"])

    f5, f6 = extract_fsd(str(five)), extract_fsd(str(six))
    assert f5["fragment_count"] == f6["fragment_count"]
    assert f5["fragment_count"] > 0
    assert len(f5["size_bins"]) == N_BINS == 196
    assert list(f5["size_bins"])[0] == f"{FSD_BIN_START}-{FSD_BIN_START + 5}"
    assert sum(f5["size_bins"].values()) == pytest.approx(1.0)


def test_fsd_rejects_an_empty_file(tmp_path):
    empty = tmp_path / "empty.frag.tsv.bgz"
    with gzip.open(empty, "wt") as f:
        f.write("")
    with pytest.raises(ValueError):
        extract_fsd(str(empty))


def test_short_and_long_fragments_are_counted_separately(tmp_path):
    path = tmp_path / "mixed.frag.tsv.bgz"
    rng = np.random.default_rng(3)
    with gzip.open(path, "wt") as f:
        for _ in range(5000):
            start = int(rng.integers(0, 100_000))
            length = int(rng.choice([120, 200]))  # short vs long
            f.write(f"chr1\t{start}\t{start + length}\t30\t+\n")
    d = extract_delfi(str(path))
    # every counted fragment is 100-220bp and classified as short or long
    assert d["mean_short_frac"] == pytest.approx(0.5, abs=0.02)
    assert d["n_windows_100kb"] == 30894


# --------------------------------------------------------------------------- #
# Fusion estimators
# --------------------------------------------------------------------------- #
def _toy_cohort(n_patients=40, seed=0):
    """A paired cohort with a real signal plus patient-invariant noise."""
    rng = np.random.default_rng(seed)
    base = rng.normal(0, 3, n_patients)
    jitter = rng.normal(0, 1, n_patients)
    panel = np.concatenate([base + 3 + jitter, base + jitter])
    frag = np.concatenate([base + 1 + jitter, base + jitter])
    y = np.concatenate([np.ones(n_patients, int), np.zeros(n_patients, int)])
    groups = np.concatenate([np.arange(n_patients), np.arange(n_patients)])
    return panel, frag, y, groups


def test_zscore_uses_only_the_fitting_rows():
    x = np.array([1.0, 2.0, 3.0, 100.0])
    z = zscore(x, np.array([0, 1, 2]))
    assert abs(float(z[:3].mean())) < 1e-9
    # the out-of-fit row is shifted, not centred
    assert z[3] > 5


def test_log_loss_is_minimised_by_a_perfect_prediction():
    y = np.array([0, 0, 1, 1])
    assert log_loss(y, np.array([0.0, 0.0, 1.0, 1.0])) < log_loss(
        y, np.array([0.5, 0.5, 0.5, 0.5]))


def test_regularized_late_fusion_beats_a_random_guesser():
    panel, frag, y, groups = _toy_cohort(seed=1)
    rng = np.random.default_rng(0)
    idx = np.arange(len(y))
    tr, te = idx[:60], idx[60:]
    out = regularized_late_fusion(panel, frag, y, tr, te, groups)
    assert out["score"].shape == te.shape
    assert np.all(np.isfinite(out["score"]))
    assert 0.0 <= out["score"].min() and out["score"].max() <= 1.0
    assert out["selected_C"] is not None
    assert 0.0 <= out["selected_shrink"] <= 1.0


def test_bayesian_model_average_weights_form_a_distribution():
    panel, frag, y, groups = _toy_cohort(seed=2)
    idx = np.arange(len(y))
    tr, te = idx[:60], idx[60:]
    out = bayesian_model_average(panel, frag, y, tr, te, groups)
    w = out["weights"]
    assert sum(w.values()) == pytest.approx(1.0, abs=1e-3)
    assert all(v >= 0 for v in w.values())
    assert out["score"].shape == te.shape
    assert np.all(np.isfinite(out["score"]))


def test_fusion_estimators_are_deterministic():
    panel, frag, y, groups = _toy_cohort(seed=3)
    idx = np.arange(len(y))
    tr, te = idx[:60], idx[60:]
    a = regularized_late_fusion(panel, frag, y, tr, te, groups)
    b = regularized_late_fusion(panel, frag, y, tr, te, groups)
    assert np.allclose(a["score"], b["score"])
    assert a["selected_C"] == b["selected_C"]


def test_fusion_estimators_never_see_the_test_labels():
    """Flipping the scored rows' labels must not change their scores."""
    panel, frag, y, groups = _toy_cohort(seed=4)
    idx = np.arange(len(y))
    tr, te = idx[:60], idx[60:]
    y2 = y.copy()
    y2[te] = 1 - y2[te]
    a = regularized_late_fusion(panel, frag, y, tr, te, groups)
    b = regularized_late_fusion(panel, frag, y2, tr, te, groups)
    assert np.allclose(a["score"], b["score"])
