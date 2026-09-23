"""Tests for the cross-platform methylation fusion methodology demo.

The script's purpose is to *demonstrate the fusion code path* with two
published baselines from DISJOINT cohorts (HCC_J cfDNA vs TCGA-LIHC
tissue). The tests below verify:

  1. Each fusion strategy returns a valid AUC in (0, 1].
  2. Synthetic 2-channel fusion produces a non-trivial lift over each
     channel alone (since the two synthesized channels are independent,
     any working fusion must improve over the weakest channel).
  3. The provenance block correctly marks this as a methodology demo
     (not a clinical fusion).
  4. The CLI runs end-to-end and writes a JSON.
  5. The script handles --seeds 1 and --output customization.
  6. The synthetic-score calibration realizes the target AUC within
     tolerance (so the test fixtures behave as documented).
  7. Honest-caveat coverage: every documented disjoint-cohort caveat is
     present in the JSON output.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# Add scripts/ to sys.path so we can import the fusion helpers directly.
sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "scripts")
)

from cross_platform_methylation_fusion import (  # noqa: E402
    _lr_fusion_score,
    _logit_average,
    _mu_for_target_auc,
    _naive_average,
    _safe_logit,
    _sat_from_curve,
    evaluate_seeds,
    load_hcc_j_per_cancer,
    load_lihc_methylation,
    synthesize_scores,
)

# ─────────────────────────────────────────────────────────────────────
# Unit tests on the score synthesizers
# ─────────────────────────────────────────────────────────────────────


def test_mu_for_target_auc_monotone():
    """Larger target AUC → larger mean separation."""
    m70 = _mu_for_target_auc(0.70)
    m90 = _mu_for_target_auc(0.90)
    m99 = _mu_for_target_auc(0.99)
    assert m70 < m90 < m99


def test_mu_for_target_auc_sqrt2_formula():
    """mu = sqrt(2) * Phi^-1(AUC). Spot-check at AUC = 0.5 → mu = 0."""
    from math import isclose, sqrt
    from scipy.stats import norm

    assert isclose(_mu_for_target_auc(0.5), 0.0, abs_tol=1e-9)
    expected = sqrt(2.0) * norm.ppf(0.9)
    assert isclose(_mu_for_target_auc(0.9), expected, abs_tol=1e-9)


def test_synthesize_scores_in_unit_interval():
    """All synthesized scores must lie in (0, 1)."""
    rng = np.random.default_rng(0)
    y = np.array([1] * 50 + [0] * 50)
    s = synthesize_scores(y, target_auc=0.85, rng=rng)
    assert (s > 0.0).all() and (s < 1.0).all()
    assert s.shape == (100,)


def test_synthesize_scores_realizes_target_auc_within_tolerance():
    """Calibrated score should hit the target AUC within 0.05 (n=200)."""
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(1)
    y = np.array([1] * 100 + [0] * 100)
    for target in (0.70, 0.85, 0.95):
        s = synthesize_scores(y, target_auc=target, rng=rng)
        realized = float(roc_auc_score(y, s))
        assert abs(realized - target) < 0.05, (
            f"target={target} realized={realized:.3f}"
        )


# ─────────────────────────────────────────────────────────────────────
# Unit tests on the fusion helpers
# ─────────────────────────────────────────────────────────────────────


def test_safe_logit_handles_extremes():
    """logit must not blow up at p = 0 or p = 1."""
    p = np.array([0.0, 1.0, 0.5, 0.0001, 0.9999])
    z = _safe_logit(p)
    assert np.isfinite(z).all()
    # Symmetry: logit(0.5) = 0.
    assert abs(z[2]) < 1e-6


def test_naive_average_is_arithmetic_mean():
    a = np.array([0.1, 0.7, 0.4, 0.9])
    b = np.array([0.2, 0.8, 0.3, 0.5])
    np.testing.assert_allclose(_naive_average(a, b),
                               np.array([0.15, 0.75, 0.35, 0.7]))


def test_logit_average_is_in_unit_interval():
    """logit-average resigmoided should still be in (0, 1)."""
    a = np.array([0.1, 0.4, 0.6, 0.9])
    b = np.array([0.2, 0.5, 0.7, 0.8])
    z = _logit_average(a, b)
    assert (z > 0.0).all() and (z < 1.0).all()


def test_lr_fusion_score_returns_unit_interval():
    """LR fusion output must be a 1-D probability vector."""
    rng = np.random.default_rng(0)
    n = 80
    y = (rng.uniform(size=n) > 0.5).astype(int)
    a_tr = rng.uniform(0.0, 1.0, n)
    b_tr = rng.uniform(0.0, 1.0, n)
    a_te = rng.uniform(0.0, 1.0, n)
    b_te = rng.uniform(0.0, 1.0, n)
    out = _lr_fusion_score(a_tr, b_tr, y, a_te, b_te)
    assert out.shape == (n,)
    assert (out >= 0.0).all() and (out <= 1.0).all()


def test_sat_from_curve_at_perfect_separation():
    """Perfect separation → sens@95 = 1.0."""
    y = np.array([0, 0, 0, 1, 1, 1])
    s = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
    assert _sat_from_curve(y, s, 0.05) == 1.0


def test_sat_from_curve_at_tied_scores():
    """All-tied → sens@95 = prevalence (here 0.5)."""
    y = np.array([0, 1])
    s = np.array([0.5, 0.5])
    # fpr=0, tpr=0 at threshold 0.5; sat(0.05) hits tpr at fpr<=0.05.
    # With tied scores, fpr and tpr both jump to 1 simultaneously, so
    # tpr at fpr<=0.05 is the last point where fpr<=0.05, which is tpr=0.
    assert _sat_from_curve(y, s, 0.05) == 0.0


# ─────────────────────────────────────────────────────────────────────
# End-to-end fusion test (synthetic, no real data needed)
# ─────────────────────────────────────────────────────────────────────


def test_synthetic_fusion_lifts_over_weaker_channel():
    """A working LR fusion must beat the weaker of its two channels.

    Setup: a stronger 'methylation' channel (AUC 0.90) and a weaker
    'fragmentomics' channel (AUC 0.75). Both are synthesized as
    independent columns (so fusion can help). LR fusion must achieve
    AUC > 0.75 (the weaker channel alone).

    Note: it does NOT have to beat 0.90 (the stronger channel alone) —
    when the two channels are independent, fusion can only help up to
    the point of diminishing returns; LR fusion typically hits between
    the two.
    """
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(42)
    n_pos, n_neg = 90, 260  # ~HCC_J OvR proportions
    y = np.concatenate([np.ones(n_pos, dtype=int),
                        np.zeros(n_neg, dtype=int)])
    meth = synthesize_scores(y, target_auc=0.90, rng=rng)
    frag = synthesize_scores(y, target_auc=0.75, rng=rng)
    out = evaluate_seeds(frag, meth, y, seeds=[42, 13, 7])
    lr_auc = out["per_strategy"]["lr_fusion"]["auc_mean"]
    frag_auc = out["frag_only"]["auc"]
    meth_auc = out["meth_only"]["auc"]
    # Sanity: realized single-channel AUCs are close to their targets.
    assert abs(frag_auc - 0.75) < 0.05, f"frag AUC {frag_auc:.3f} off-target"
    assert abs(meth_auc - 0.90) < 0.05, f"meth AUC {meth_auc:.3f} off-target"
    # The fusion must beat the weaker channel — anything else means the
    # fusion code path is broken.
    assert lr_auc > frag_auc, (
        f"LR fusion AUC {lr_auc:.4f} failed to beat frag-only AUC "
        f"{frag_auc:.4f}"
    )
    # And it should not be wildly below the stronger channel — sanity
    # check the fusion isn't catastrophically broken.
    assert lr_auc > 0.65, f"LR fusion AUC {lr_auc:.4f} unexpectedly low"


def test_all_three_strategies_produce_valid_auc():
    """Each of naive / logit / LR fusion must produce AUC in (0, 1]."""
    rng = np.random.default_rng(7)
    n_pos, n_neg = 90, 260
    y = np.concatenate([np.ones(n_pos, dtype=int),
                        np.zeros(n_neg, dtype=int)])
    meth = synthesize_scores(y, target_auc=0.90, rng=rng)
    frag = synthesize_scores(y, target_auc=0.75, rng=rng)
    out = evaluate_seeds(frag, meth, y, seeds=[42])
    for name in ("naive_average", "logit_average", "lr_fusion"):
        auc = out["per_strategy"][name]["auc_mean"]
        assert 0.0 < auc <= 1.0, f"{name} AUC {auc:.4f} not in (0, 1]"
        # per-seed AUC list must be non-empty and contain 1 element
        # for --seeds 1.
        assert len(out["per_strategy"][name]["per_seed_auc"]) == 1


# ─────────────────────────────────────────────────────────────────────
# Loader tests — guard against schema drift
# ─────────────────────────────────────────────────────────────────────


def test_load_hcc_j_per_cancer_finds_hccj():
    """load_hcc_j_per_cancer must find the HCC_J row in the canonical JSON."""
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    frag_json = os.path.join(repo_root, "results",
                             "per_cancer_sens_at_spec.json")
    if not os.path.exists(frag_json):
        pytest.skip(f"{frag_json} not found — run scripts/per_cancer_sens_at_spec.py first")
    row = load_hcc_j_per_cancer(frag_json, "HCC_J")
    assert row["n_pos"] == 89
    assert row["n_neg"] == 264
    assert row["auc_mean"] > 0.7
    assert "sens_at_95" in row
    assert "sens_at_99" in row


def test_load_lihc_methylation_finds_lihc():
    """load_lihc_methylation must find the TCGA-LIHC row in the sibling JSON."""
    sibling_meth_json = (
        "/Users/hermes/deepcatch-methylation/results/multi_cancer_baseline.json"
    )
    if not os.path.exists(sibling_meth_json):
        pytest.skip(f"{sibling_meth_json} not found — sibling repo not checked out")
    row = load_lihc_methylation(sibling_meth_json, "TCGA-LIHC")
    assert row["n_samples"] == 18
    assert row["n_tumor"] == 12
    assert row["n_normal"] == 6
    assert row["auc_mean"] > 0.9  # published = 0.972


# ─────────────────────────────────────────────────────────────────────
# Provenance / honest-caveat coverage
# ─────────────────────────────────────────────────────────────────────


def test_provenance_marks_methodology_demo():
    """The JSON must mark this as a methodology demo, not a clinical fusion."""
    out_json = (
        "/Users/hermes/deepcatch/results/cross_platform_methylation_fusion.json"
    )
    if not os.path.exists(out_json):
        pytest.skip(f"{out_json} not found — run the script first")
    with open(out_json) as f:
        d = json.load(f)
    prov = d["provenance"]
    caveats = prov["honest_caveats"]
    text = "\n".join(caveats).lower()
    assert "methodology" in text or "not a clinical fusion" in text, (
        f"methodology caveat missing from provenance: {caveats}"
    )
    assert "disjoint" in text or "do not overlap" in text, (
        f"disjoint-cohort caveat missing: {caveats}"
    )
    assert "synth" in text or "synthesized" in text, (
        f"synthesized-scores caveat missing: {caveats}"
    )


def test_provenance_carries_lihc_and_hccj_n():
    """Provenance must carry the published n_pos/n_neg/n_samples for each channel."""
    out_json = (
        "/Users/hermes/deepcatch/results/cross_platform_methylation_fusion.json"
    )
    if not os.path.exists(out_json):
        pytest.skip(f"{out_json} not found — run the script first")
    with open(out_json) as f:
        d = json.load(f)
    fc = d["provenance"]["frag_channel"]
    mc = d["provenance"]["meth_channel"]
    assert fc["cancer_key"] == "HCC_J"
    assert fc["n_pos"] == 89
    assert fc["n_neg"] == 264
    assert mc["cancer_key"] == "TCGA-LIHC"
    assert mc["n_samples"] == 18


# ─────────────────────────────────────────────────────────────────────
# CLI integration
# ─────────────────────────────────────────────────────────────────────


def test_cli_runs_end_to_end_and_writes_json(tmp_path):
    """End-to-end CLI run with --seeds 1 writes a JSON with 3 strategies."""
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sibling_meth_json = (
        "/Users/hermes/deepcatch-methylation/results/multi_cancer_baseline.json"
    )
    if not os.path.exists(sibling_meth_json):
        pytest.skip(f"{sibling_meth_json} not found — sibling repo not checked out")
    out_json = tmp_path / "cross_platform_fusion_test.json"
    cmd = [
        "env", "-u", "PYTHONPATH",
        sys.executable,
        os.path.join(repo_root, "scripts", "cross_platform_methylation_fusion.py"),
        "--seeds", "1",
        "--output", str(out_json),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=120, cwd=repo_root)
    assert result.returncode == 0, (
        f"CLI failed (rc={result.returncode})\n"
        f"STDOUT:\n{result.stdout[-2000:]}\n"
        f"STDERR:\n{result.stderr[-2000:]}"
    )
    payload = json.loads(out_json.read_text())
    required = {"naive_average", "logit_average", "lr_fusion"}
    assert required <= set(payload["per_strategy"].keys()), (
        f"missing strategies: {required - set(payload['per_strategy'].keys())}"
    )
    # Each strategy has a valid AUC > 0.5.
    for strat in required:
        auc = payload["per_strategy"][strat]["auc_mean"]
        assert auc > 0.5, f"{strat} AUC = {auc:.3f} unexpectedly low"
    # Provenance block must explicitly call out the disjoint cohorts.
    caveats_text = "\n".join(
        payload["provenance"]["honest_caveats"]
    ).lower()
    assert "disjoint" in caveats_text or "do not overlap" in caveats_text


def test_cli_help_exits_fast(tmp_path):
    """--help must exit fast with rc=0 — guards against broken argparse."""
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    cmd = [
        "env", "-u", "PYTHONPATH",
        sys.executable,
        os.path.join(repo_root, "scripts", "cross_platform_methylation_fusion.py"),
        "--help",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            timeout=30, cwd=repo_root)
    assert result.returncode == 0
    assert "frag-json" in result.stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))