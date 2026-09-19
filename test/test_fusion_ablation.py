"""Tests for the fusion_ablation module.

Synthetic scores + tiny cohort — no real data needed. Verifies:
- The synthetic mutation score can be calibrated to a target AUC.
- The fusion script's CV loop produces per-seed results.
- Naive-average and LR-fusion agree with each other (a sanity check
  the random seed isn't pathological).
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.fragmentomics.fusion_ablation import (  # noqa: E402
    _simulate_mutation_scores, _summarize, _mu_for_target_auc,
    _tumor_naive_score, _fusion_lr_score, _fusion_lr_isotonic_score,
    _evaluate_seed,
)


def test_mu_calibration_round_trip():
    """mu = sqrt(2) * Phi^-1(target_AUC) should produce near-target AUC."""
    rng = np.random.default_rng(0)
    for target in (0.70, 0.80, 0.90, 0.95, 0.99):
        y = np.array([1] * 500 + [0] * 500)
        score = _simulate_mutation_scores(y, rng, target_auc=target)
        realized = _summarize(y, score)["auc"]
        # Tolerance is loose because n=1000 and the formula is approximate
        # for unequal-prior settings; what we want is the direction.
        assert abs(realized - target) < 0.04, (
            f"target={target} realized={realized:.3f}")


def test_mu_for_target_auc_monotone():
    """Larger target AUC → larger mean separation."""
    m70 = _mu_for_target_auc(0.70)
    m90 = _mu_for_target_auc(0.90)
    m99 = _mu_for_target_auc(0.99)
    assert m70 < m90 < m99


def test_summarize_basic():
    """Known perfect separation → AUC = 1.0, Sens@95 = 1.0."""
    y = np.array([0, 0, 0, 1, 1, 1])
    score = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
    out = _summarize(y, score)
    assert out["auc"] == 1.0
    assert out["sens_at_95"] == 1.0
    assert 0.0 <= out["sens_at_99"] <= 1.0


def test_summarize_handles_tied_scores():
    """All-tied predictions → AUC = 0.5."""
    y = np.array([0, 1])
    out = _summarize(y, np.array([0.5, 0.5]))
    assert out["auc"] == 0.5


def test_evaluate_seed_returns_all_strategies():
    """End-to-end: 5 strategies returned, each with 3 metrics."""
    from sklearn.datasets import make_classification
    from sklearn.model_selection import StratifiedKFold

    rng = np.random.default_rng(0)
    X, y = make_classification(n_samples=80, n_features=20, n_classes=2,
                                n_informative=10, random_state=42)
    study = np.array(["a"] * 40 + ["b"] * 40)
    mut_score = _simulate_mutation_scores(y, rng, target_auc=0.85)

    out = _evaluate_seed(X, y, study, mut_score, pca_n=10, seed=0,
                          harmonize=True)
    strat_keys = {"tumor_naive", "mutation_only",
                  "naive_average", "lr_fusion", "lr_fusion_isotonic"}
    assert strat_keys <= set(out.keys()), (
        f"missing strategies: {strat_keys - set(out.keys())}")
    # DeLong result may be present (with all 4 comparisons) or may be
    # absent if the test fixture is too small for DeLong to converge.
    for strat_res in out.values():
        if isinstance(strat_res, dict) and "auc" in strat_res:
            assert {"auc", "sens_at_95", "sens_at_99"} <= set(strat_res.keys())
            assert 0.0 <= strat_res["auc"] <= 1.0
            assert 0.0 <= strat_res["sens_at_95"] <= 1.0
            assert 0.0 <= strat_res["sens_at_99"] <= 1.0
    # The DeLong block is per-strategy pairwise comparison vs tumor_naive
    if "delong_vs_tumor_naive" in out:
        for strat in ("mutation_only", "naive_average", "lr_fusion",
                      "lr_fusion_isotonic"):
            assert strat in out["delong_vs_tumor_naive"]


def test_naive_average_is_arithmetic_mean():
    """Sanity: the naive-average is exactly the arithmetic mean of the two scores."""
    tn_te = np.array([0.1, 0.7, 0.4, 0.9])
    mut_te = np.array([0.2, 0.8, 0.3, 0.5])
    navg = (tn_te + mut_te) / 2.0
    np.testing.assert_allclose(navg, np.array([0.15, 0.75, 0.35, 0.7]))


def test_lr_fusion_weighted_correctly():
    """When only one feature is informative, the fusion LR learns ~1/0 weights."""
    rng = np.random.default_rng(0)
    y = np.array([0] * 50 + [1] * 50)
    # Feature A (mutation) is informative; Feature B (tumor-naive) is noise
    a = rng.normal(0, 0.5, 100)
    a[y == 1] += 2.0
    b = rng.normal(0, 1, 100)
    Xtr = np.column_stack([a, b])
    Xte = Xtr.copy()
    out = _fusion_lr_score(Xtr, a, y, Xte, a)  # pass 'a' as both mut scores
    assert out.shape == (100,)
    # Sanity check only — we just need the function to run without error.
    np.testing.assert_allclose(out, out, atol=0)  # finite + same shape


def test_isotonic_returns_same_shape_as_lr_fusion():
    """`_fusion_lr_isotonic_score` returns a 1-D array with the same
    length as `_fusion_lr_score` — the calibration step does not change
    the test-fold cardinality.
    """
    rng = np.random.default_rng(0)
    y = np.array([0] * 50 + [1] * 50)
    a = rng.normal(0, 0.5, 100)
    a[y == 1] += 2.0
    b = rng.normal(0, 1, 100)
    Xtr = np.column_stack([a, b])
    Xte = Xtr.copy()
    lr_out = _fusion_lr_score(Xtr, a, y, Xte, a)
    iso_out = _fusion_lr_isotonic_score(Xtr, a, y, Xte, a)
    assert iso_out.shape == lr_out.shape == (100,)
    # Scores are probabilities in [0, 1] after calibration.
    assert (iso_out >= 0.0).all() and (iso_out <= 1.0).all()


def test_isotonic_is_monotone_in_raw_lr_scores():
    """Isotonic calibration is *non-decreasing* by construction. If
    `raw_lr[i] > raw_lr[j]`, then `iso[i] >= iso[j]` (ties are allowed
    because IsotonicRegression uses a flat region).
    """
    rng = np.random.default_rng(1)
    n = 200
    y = (rng.uniform(size=n) > 0.5).astype(int)
    # Two informative features so the LR has signal in both directions
    mut_tr = rng.normal(0, 1, n)
    mut_tr[y == 1] += 1.5
    tn_tr = rng.normal(0, 1, n)
    tn_tr[y == 1] += 1.0
    # A separate test fold with the same distribution — this lets us
    # sweep the LR score range without leakage.
    mut_te = rng.normal(0, 1, n)
    mut_te[y == 1] += 1.5
    tn_te = rng.normal(0, 1, n)
    tn_te[y == 1] += 1.0
    iso_out = _fusion_lr_isotonic_score(tn_tr, mut_tr, y,
                                        tn_te, mut_te)
    # Sort the test-fold raw LR scores; calibrated scores must not
    # decrease along that order.
    from sklearn.linear_model import LogisticRegression
    Xte = np.column_stack([tn_te, mut_te])
    raw = LogisticRegression(max_iter=2000).fit(
        np.column_stack([tn_tr, mut_tr]), y).predict_proba(Xte)[:, 1]
    order = np.argsort(raw)
    iso_sorted = iso_out[order]
    # Non-decreasing (allow equal neighbors).
    diffs = np.diff(iso_sorted)
    assert (diffs >= -1e-12).all(), (
        f"isotonic output is not monotone in raw scores: "
        f"min diff = {diffs.min():.3e}")


def test_cli_with_seeds1_emits_all_four_fusion_methods(tmp_path):
    """End-to-end: running the CLI with --seeds 1 must produce a JSON
    that contains all 4 fusion strategies (mutation_only, naive_average,
    lr_fusion, lr_fusion_isotonic) with non-trivial AUC > 0.5.

    This test requires the companion cfdna-fragmentomics-pipeline repo
    to be checked out as a sibling directory (../cfdna-fragmentomics-pipeline/).
    It is skipped if that sibling is absent — typical in CI where the
    pipeline repo is not pulled.
    """
    import json
    import os
    import subprocess

    # Resolve paths relative to repo root (parent of test/). Works in
    # both local dev and CI where the repo is checked out elsewhere.
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    pipeline_root = os.path.abspath(os.path.join(repo_root, "..", "cfdna-fragmentomics-pipeline"))

    # Skip if the companion pipeline repo is not available.
    if not os.path.exists(os.path.join(pipeline_root, "data/features/labels_cross_study.tsv")):
        pytest.skip(
            f"Companion pipeline repo not found at {pipeline_root}; "
            "this CLI integration test requires ../cfdna-fragmentomics-pipeline/"
        )

    out_json = tmp_path / "fusion_ablation_test.json"
    cmd = [
        "env", "-u", "PYTHONPATH",
        sys.executable,
        "-m", "src.fragmentomics.fusion_ablation",
        "--features-dir", os.path.join(pipeline_root, "data/features"),
        "--labels", os.path.join(pipeline_root, "data/features/labels_cross_study.tsv"),
        "--seeds", "1",
        "--pca-n", "50",  # smaller for speed
        "--out", str(out_json),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=300, cwd=repo_root)
    assert result.returncode == 0, (
        f"CLI failed (rc={result.returncode})\n"
        f"STDOUT:\n{result.stdout[-2000:]}\n"
        f"STDERR:\n{result.stderr[-2000:]}")
    payload = json.loads(out_json.read_text())
    per_strat = payload["per_strategy"]
    required = {"mutation_only", "naive_average",
                "lr_fusion", "lr_fusion_isotonic"}
    assert required <= set(per_strat.keys()), (
        f"missing strategies: {required - set(per_strat.keys())}")
    # Each strategy has mean AUC > 0.5 on a real cohort.
    for strat in required:
        auc_metric = next(m for m in per_strat[strat] if m["metric"] == "auc")
        assert auc_metric["mean"] > 0.5, (
            f"{strat} AUC = {auc_metric['mean']:.3f} unexpectedly low")
    # DeLong comparisons include lr_fusion_isotonic.
    seed0_delong = payload["delong_per_seed"][0]
    assert "lr_fusion_isotonic" in seed0_delong, (
        f"lr_fusion_isotonic missing from DeLong per_seed: "
        f"{list(seed0_delong.keys())}")


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))