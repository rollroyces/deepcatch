"""The external validation must harmonise each publication on its own.

What went wrong
---------------
``scripts/external_validation_finaledb.py`` documents the protocol as
"per-study z-score harmonization, each publication normalized by its own
feature statistics" and claims to be "identical to the LOSO benchmark" --
but it z-scored the held-out cohort as a SINGLE pool. With 238 case-only
samples from publication 9 in that pool, publications 1 and 7 were scored
against publication 9's mean and variance, which is neither what the
protocol says nor what the LOSO benchmark does.

The difference is not small. Measured on the held-out cohort:

    train global  / test global          pub 1 = 0.4615, pooled = 0.6228
    train global  / test per-publication  pub 1 = 0.6205, pooled = 0.7407
    train per-study / test per-publication pub 1 = 0.6154, pooled = 0.7037

These tests pin the harmonisation rule, so the pool-wide version cannot
come back silently.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "scripts"))
sys.path.insert(0, str(_REPO_ROOT))


def _load():
    spec = importlib.util.spec_from_file_location(
        "extval", _REPO_ROOT / "scripts" / "external_validation_finaledb.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ev = _load()


def test_zscore_by_group_centres_each_group_on_its_own():
    """Two groups with different offsets: each must be centred on its own
    mean, not on the pooled mean."""
    rng = np.random.default_rng(0)
    a = rng.normal(10.0, 1.0, size=(50, 4))
    b = rng.normal(-5.0, 2.0, size=(30, 4))
    X = np.vstack([a, b])
    groups = np.array(["a"] * 50 + ["b"] * 30)
    out = ev._zscore_by_group(X, groups)

    assert np.allclose(out[:50].mean(axis=0), 0.0, atol=1e-9)
    assert np.allclose(out[50:].mean(axis=0), 0.0, atol=1e-9)
    assert np.allclose(out[:50].std(axis=0), 1.0, atol=1e-9)
    assert np.allclose(out[50:].std(axis=0), 1.0, atol=1e-9)


def test_pooled_zscoring_would_centre_the_pool_not_the_groups():
    """The control: the same matrix z-scored as ONE pool does NOT centre
    each group, which is what the old code did."""
    rng = np.random.default_rng(0)
    a = rng.normal(10.0, 1.0, size=(50, 4))
    b = rng.normal(-5.0, 2.0, size=(30, 4))
    X = np.vstack([a, b])
    pooled = ev._zscore(X)
    assert not np.allclose(pooled[:50].mean(axis=0), 0.0, atol=1e-6)


def test_the_script_harmonises_per_publication():
    """The held-out cohort must be normalised per publication, and the
    training cohort per study."""
    source = (_REPO_ROOT / "scripts"
              / "external_validation_finaledb.py").read_text()
    assert "_zscore_by_group(Xte, stte)" in source, (
        "the held-out cohort must be z-scored per publication, not as one "
        "pool")
    assert "_zscore_by_group(Xtr, sttr)" in source, (
        "the training cohort must be z-scored per study, matching the "
        "in-cohort benchmark")
    assert "Xte_z = _zscore(Xte)" not in source, (
        "the held-out cohort must not be z-scored as a single pool")


def test_the_artifact_records_the_harmonisation_sensitivity():
    """The artifact must carry the other two harmonisations, so a reader
    can see how much the choice matters."""
    source = (_REPO_ROOT / "scripts"
              / "external_validation_finaledb.py").read_text()
    assert "harmonization_sensitivity" in source
    for tag in ("train_global_test_global",
                "train_global_test_per_publication",
                "train_per_study_test_per_publication"):
        assert tag in source, f"the sensitivity block is missing {tag}"


@pytest.mark.parametrize("groups", [[], ["a"]])
def test_zscore_by_group_handles_degenerate_groups(groups):
    X = np.zeros((0, 3))
    if groups:
        X = np.ones((1, 3))
    g = np.array(groups * 1 if groups else [])
    out = ev._zscore_by_group(X, g)
    assert out.shape == X.shape
