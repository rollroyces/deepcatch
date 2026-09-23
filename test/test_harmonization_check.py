"""Tests for scripts/harmonization_check.py.

Covers:
- per-study z-score contract: returns zero mean / unit std on TRAIN data
- the train-only fit (regression test for the data-leak pitfall the
  cfdna-fragmentomics skill warns about)
- unseen-study fallback in `apply_per_study_zscore` is identity
- end-to-end: synthetic multi-cohort fixture produces a non-zero AUC
  delta when a clear batch effect is injected
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPTS))

# Skip everything if real MAF cache is missing — we need it for the
# end-to-end test that calls build_fixture (which uses load_tcga_cohort).
# The unit tests for the z-score functions don't need any data.
TCGA_CACHE = ROOT / "validation" / "tcga" / "tcga_cache"
_HAS_MAF_CACHE = TCGA_CACHE.exists() and any(TCGA_CACHE.glob("*.maf.gz"))


# ── 1. z-score fit returns zero mean / unit std on TRAIN data ─────────

def test_fit_per_study_zscore_returns_zero_mean_unit_std_on_train():
    from harmonization_check import fit_per_study_zscore
    rng = np.random.default_rng(0)
    # Two studies with very different means/stds
    X = np.vstack([
        rng.normal(loc=5.0, scale=2.0, size=(40, 2)),
        rng.normal(loc=-2.0, scale=0.5, size=(40, 2)),
    ])
    study = np.array([0] * 40 + [1] * 40)

    params = fit_per_study_zscore(X, study)

    # Re-derive the per-study mean/std on the SAME training matrix — must
    # be (≈0, ≈1) per study after applying the fit.
    for s, (mu, sd) in params.items():
        mask = study == s
        transformed = (X[mask] - mu) / sd
        assert np.allclose(transformed.mean(axis=0), 0.0, atol=1e-10), (
            f"study {s} mean should be 0 after fit, got {transformed.mean(axis=0)}"
        )
        assert np.allclose(transformed.std(axis=0), 1.0, atol=1e-6), (
            f"study {s} std should be 1 after fit, got {transformed.std(axis=0)}"
        )


def test_fit_per_study_zscore_does_not_fit_on_full_matrix():
    """Regression test: the z-score must be fit on TRAIN ONLY.

    If `fit_per_study_zscore` ever got refactored to accept the full
    matrix, the test AUC would inflate (test data leaks through the
    transform). This test pins the function signature to the train-only
    contract by checking that re-fitting on (train, test) — separately —
    gives DIFFERENT params than using only the train matrix.
    """
    from harmonization_check import fit_per_study_zscore
    rng = np.random.default_rng(1)
    # Two clusters: train and test have very different means for the same
    # study. The TRAIN-only fit must encode only the train distribution.
    X_train = np.vstack([
        rng.normal(loc=0.0, scale=1.0, size=(30, 1)),
        rng.normal(loc=10.0, scale=1.0, size=(30, 1)),
    ])
    study_train = np.array([0] * 30 + [1] * 30)

    X_test = np.vstack([
        rng.normal(loc=3.0, scale=1.0, size=(30, 1)),
        rng.normal(loc=12.0, scale=1.0, size=(30, 1)),
    ])
    study_test = np.array([0] * 30 + [1] * 30)

    params_train_only = fit_per_study_zscore(X_train, study_train)
    params_full = fit_per_study_zscore(
        np.vstack([X_train, X_test]),
        np.concatenate([study_train, study_test]),
    )

    # The TRAIN-only mean for each study must match the empirical mean of
    # the TRAIN rows (not the pooled rows).
    for s, (mu_tr, _) in params_train_only.items():
        train_mean = X_train[study_train == s].mean(axis=0)
        assert np.allclose(mu_tr, train_mean, atol=1e-6), (
            f"study {s}: fit-param mean={mu_tr}, expected train mean={train_mean}"
        )
    # And the full-matrix params must differ from train-only params for
    # at least one study (the test rows shifted the mean).
    diffs = [
        float(np.abs(params_full[s][0] - params_train_only[s][0]).max())
        for s in params_train_only
    ]
    assert max(diffs) > 1e-3, (
        "train-only and full-matrix params should differ for at least "
        "one study; got diffs = " + str(diffs)
    )


def test_apply_per_study_zscore_unseen_study_falls_back_to_identity():
    """If a study appears in the matrix but not in the fit params, the
    transform must leave it unchanged (no NaN, no zero, no random)."""
    from harmonization_check import apply_per_study_zscore
    X = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]])
    study = np.array([0, 0, 1, 1])
    # Only fit params for study 0 — study 1 is unseen.
    params = {0: (np.array([2.0, 3.0]), np.array([1.0, 1.0]))}
    out = apply_per_study_zscore(X, study, params)
    # Study 0 should be centered exactly to (-1, -1) and (1, 1)
    assert np.allclose(out[study == 0], np.array([[-1, -1], [1, 1]]))
    # Study 1 should be UNCHANGED (no params → identity)
    assert np.allclose(out[study == 1], X[study == 1])


# ── 4. end-to-end AUC delta with clear batch effect ───────────────────

@pytest.mark.skipif(
    not _HAS_MAF_CACHE,
    reason="real TCGA MAF cache missing — fixture build needs it",
)
def test_synthetic_fixture_auc_delta_when_batch_effect_is_dominant():
    """If study label perfectly predicts the cancer label, harmonization
    should DRAMATICALLY reduce AUC. This is the study-as-classifier
    trap (cfdna-fragmentomics skill: 0.999 → 0.497 on only-Jiang-cancer
    vs only-Cristiano-healthy).

    We inject a synthetic version of this: study A = 100% cancer,
    study D = 100% control. The raw pipeline would learn "high score
    means study_A" and inflate AUC. Per-study z-scoring removes the
    study signal and AUC should drop toward chance.
    """
    from harmonization_check import (
        build_feature_matrix,
        cv_evaluate,
        fit_per_study_zscore,
        apply_per_study_zscore,
        STUDIES,
    )
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold

    rng = np.random.default_rng(7)
    n_per_study = 40
    # Cancer signal: study_A is 100% cancer, study_D is 100% control.
    # Studies B/C are mixed so the model can't trivially learn study.
    # Strong raw panel signal: high for A, low for D, middling for B/C.
    panel = np.concatenate([
        rng.normal(loc=+2.0, scale=0.5, size=n_per_study),  # study 0: A — cancer
        rng.normal(loc=+0.5, scale=0.5, size=n_per_study),  # study 1: B
        rng.normal(loc=+0.0, scale=0.5, size=n_per_study),  # study 2: C
        rng.normal(loc=-2.0, scale=0.5, size=n_per_study),  # study 3: D — control
    ])
    frag = rng.normal(0, 1, size=4 * n_per_study)
    # Cancer label — INDEPENDENTLY drawn within each study
    # (so the per-study class distribution is what creates the confound).
    # For the test we make study 0 100% cancer and study 3 100% control:
    y = np.array([1] * n_per_study + [0] * n_per_study + [0] * n_per_study + [0] * n_per_study)
    # Re-cast so studies 0 and 3 actually have a y-induced batch effect.
    # Study 0 = all cancer (y=1) — strong cancer label bias.
    # Study 3 = all control (y=0) — strong control label bias.
    # Studies 1 and 2 are mixed 50/50 so they don't add information.
    y = np.concatenate([
        np.ones(n_per_study, dtype=int),                        # study 0
        rng.binomial(1, 0.5, n_per_study),                       # study 1
        rng.binomial(1, 0.5, n_per_study),                       # study 2
        np.zeros(n_per_study, dtype=int),                        # study 3
    ])
    study = np.array([0] * n_per_study + [1] * n_per_study + [2] * n_per_study + [3] * n_per_study)
    assert len(y) == len(panel) == 4 * n_per_study

    X = build_feature_matrix(panel, frag)

    # Raw evaluation — should be HIGH AUC because study_A scores are
    # systematically higher and study_A is all cancer.
    raw_result = cv_evaluate(X, y, study, harmonize=False, seed=42)
    # Harmonized — should DROP because per-study z-scoring removes the
    # study-mean signal that was acting as a class proxy.
    harm_result = cv_evaluate(X, y, study, harmonize=True, seed=42)

    raw_auc = raw_result["pooled_auc"]
    harm_auc = harm_result["pooled_auc"]
    assert raw_auc > 0.85, (
        f"raw AUC should be high (>0.85) under the study-as-classifier "
        f"confound; got {raw_auc:.4f}"
    )
    assert harm_auc < raw_auc, (
        f"harmonized AUC ({harm_auc:.4f}) should drop below raw AUC "
        f"({raw_auc:.4f}) when study perfectly predicts class"
    )
    assert harm_auc < 0.7, (
        f"harmonized AUC should drop toward chance (<0.7) when study "
        f"is a perfect class proxy; got {harm_auc:.4f}"
    )


@pytest.mark.skipif(
    not _HAS_MAF_CACHE,
    reason="real TCGA MAF cache missing — fixture build needs it",
)
def test_synthetic_fixture_with_and_without_harmonization_produce_different_aucs():
    """The end-to-end test from the task spec: build the real fixture
    via `build_fixture`, run the pipeline with and without
    harmonization, and assert they produce DIFFERENT AUCs (i.e. the
    harmonization code path is wired up correctly and the bias is
    non-trivial).
    """
    from harmonization_check import build_fixture, build_feature_matrix, cv_evaluate

    fix = build_fixture(seed=42)
    X = build_feature_matrix(fix["panel_scores"], fix["frag_scores"])
    y = fix["y"]
    study = fix["study"]

    raw = cv_evaluate(X, y, study, harmonize=False, seed=42)
    harm = cv_evaluate(X, y, study, harmonize=True, seed=42)

    # The two AUCs should differ — the bias is injected, so harmonization
    # has to do SOMETHING. The magnitude is documented in the JSON; we
    # just assert they're not bit-identical.
    assert raw["pooled_auc"] != harm["pooled_auc"], (
        "raw and harmonized pipelines produced identical AUC — "
        "harmonization code path appears not wired up correctly"
    )
    # Both should still produce a valid AUC > chance
    assert raw["pooled_auc"] > 0.5
    assert harm["pooled_auc"] > 0.5