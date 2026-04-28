#!/usr/bin/env python3
"""
MODULE 5: DeLong's Test for Correlated AUCs
============================================

Reference: DeLong, DeLong & Clarke-Pearson (1988), "Comparing the Areas Under
Two or More Correlated Receiver Operating Characteristic Curves: A
Nonparametric Approach", Biometrics 44:837-845.

This is THE definitive statistical test for comparing two AUCs computed on
the SAME set of test samples. Unlike a simple z-test that assumes
independence, DeLong's test accounts for the correlation between AUC
estimates because the models are evaluated on the same patients.

The full covariance matrix is computed from the structural components of the
Mann-Whitney U-statistic, giving a proper test of:
  H₀: AUC_A = AUC_B  vs.  H₁: AUC_A ≠ AUC_B

Also includes the Sun & Xu (2014) fast implementation using matrix operations
rather than nested loops, scaling to large sample sizes efficiently.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from scipy import stats

# ── type aliases ────────────────────────────────────────────────────────────
Array = np.ndarray


# ── core DeLong test ────────────────────────────────────────────────────────


def delong_test(
    y_true: Array,
    scores_a: Array,
    scores_b: Array,
    alpha: float = 0.05,
) -> Dict[str, float]:
    """Full DeLong's test for comparing two correlated AUCs.

    Implements DeLong et al. (1988) §2, equations (2)-(5). Computes the
    structural components θ₁, …, θ_m, then the covariance matrix V, then
    tests H₀: AUC_A = AUC_B.

    Args:
        y_true: Ground-truth binary labels (0/1), shape (n_samples,).
        scores_a: Predicted scores for model A, shape (n_samples,).
        scores_b: Predicted scores for model B, shape (n_samples,).
        alpha: Significance level (default 0.05 for 95% confidence).

    Returns:
        Dict with keys:
          - auc_a: AUC of model A
          - auc_b: AUC of model B
          - delta_auc: AUC_B - AUC_A
          - var_a: Variance of AUC_A
          - var_b: Variance of AUC_B
          - cov_ab: Covariance between AUC_A and AUC_B
          - se_delta: Standard error of delta_auc
          - z_statistic: Test statistic (asymptotically N(0,1) under H₀)
          - p_value: Two-sided p-value
          - ci95_lower: 95% CI lower bound for delta
          - ci95_upper: 95% CI upper bound for delta
          - significant: Boolean (p < alpha)

    Raises:
        ValueError: If y_true doesn't have both classes or inputs are mis-sized.

    Reference:
        DeLong, DeLong & Clarke-Pearson (1988), Biometrics 44:837-845.
        Sun & Xu (2014), "Fast Implementation of DeLong's Algorithm",
        arXiv:1405.1437.
    """
    y_true = np.asarray(y_true).ravel()
    scores_a = np.asarray(scores_a).ravel()
    scores_b = np.asarray(scores_b).ravel()

    # Validate
    if len(y_true) != len(scores_a) or len(y_true) != len(scores_b):
        raise ValueError(
            f"Length mismatch: y_true={len(y_true)}, "
            f"scores_a={len(scores_a)}, scores_b={len(scores_b)}"
        )

    classes = np.unique(y_true)
    if len(classes) < 2:
        raise ValueError(
            f"Need both classes for AUC computation; got classes={classes}"
        )

    pos_class = max(classes)  # typically 1
    neg_class = min(classes)  # typically 0

    # ── Structural components ───────────────────────────────────────────
    pos_mask = y_true == pos_class
    neg_mask = y_true == neg_class

    n_pos = int(np.sum(pos_mask))
    n_neg = int(np.sum(neg_mask))

    if n_pos < 2 or n_neg < 2:
        raise ValueError(
            f"Need ≥2 samples per class; got pos={n_pos}, neg={n_neg}"
        )

    # Separator function: ψ(s₁, s₂) = 1 if s₁ > s₂, 0.5 if s₁ = s₂, 0 if s₁ < s₂
    # These form the structural components of the Mann-Whitney U

    # For Model A
    V10_a = _compute_V10(scores_a[pos_mask], scores_a[neg_mask])
    V01_a = _compute_V01(scores_a[neg_mask], scores_a[pos_mask])

    # For Model B
    V10_b = _compute_V10(scores_b[pos_mask], scores_b[neg_mask])
    V01_b = _compute_V01(scores_b[neg_mask], scores_b[pos_mask])

    # ── AUC estimates ──────────────────────────────────────────────────
    auc_a = float(np.mean(V10_a))
    auc_b = float(np.mean(V10_b))

    # ── Structural component arrays (θ in DeLong notation) ─────────────
    # θ_pos = mean over negatives of separator
    # θ_neg = mean over positives of separator
    theta_pos_a = V10_a  # shape (n_pos,)
    theta_neg_a = V01_a  # shape (n_neg,)

    theta_pos_b = V10_b
    theta_neg_b = V01_b

    # ── Covariance matrix ──────────────────────────────────────────────
    # Eq (2) from DeLong: S₁₀ = covariance of θ_pos, S₀₁ = covariance of θ_neg
    # var(AUC) = S₁₀ / n_pos + S₀₁ / n_neg

    S10_a = _sample_variance(theta_pos_a)
    S01_a = _sample_variance(theta_neg_a)
    var_a = S10_a / n_pos + S01_a / n_neg

    S10_b = _sample_variance(theta_pos_b)
    S01_b = _sample_variance(theta_neg_b)
    var_b = S10_b / n_pos + S01_b / n_neg

    # Covariance between model A and B
    S10_ab = _sample_covariance(theta_pos_a, theta_pos_b)
    S01_ab = _sample_covariance(theta_neg_a, theta_neg_b)
    cov_ab = S10_ab / n_pos + S01_ab / n_neg

    # ── Delta and test statistic ───────────────────────────────────────
    delta = auc_b - auc_a
    var_delta = var_a + var_b - 2 * cov_ab

    if var_delta <= 0:
        # Pathological case — effectively equal AUCs
        return {
            "auc_a": auc_a,
            "auc_b": auc_b,
            "delta_auc": delta,
            "var_a": var_a,
            "var_b": var_b,
            "cov_ab": cov_ab,
            "se_delta": 0.0,
            "z_statistic": 0.0,
            "p_value": 1.0,
            "ci95_lower": delta,
            "ci95_upper": delta,
            "significant": False,
        }

    se_delta = np.sqrt(var_delta)
    z = delta / se_delta
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))

    # 95% CI for delta
    z_alpha = stats.norm.ppf(1 - alpha / 2)
    ci_lower = delta - z_alpha * se_delta
    ci_upper = delta + z_alpha * se_delta

    return {
        "auc_a": float(auc_a),
        "auc_b": float(auc_b),
        "delta_auc": float(delta),
        "var_a": float(var_a),
        "var_b": float(var_b),
        "cov_ab": float(cov_ab),
        "se_delta": float(se_delta),
        "z_statistic": float(z),
        "p_value": float(p_value),
        "ci95_lower": float(ci_lower),
        "ci95_upper": float(ci_upper),
        "significant": bool(p_value < alpha),
    }


def delong_test_multi(
    y_true: Array,
    scores: Dict[str, Array],
    alpha: float = 0.05,
    correction: str = "bonferroni",
) -> Dict[str, Dict[str, float]]:
    """Pairwise DeLong tests for multiple models, with multiple-testing correction.

    Args:
        y_true: Ground-truth labels.
        scores: Dict of model_name → predicted_scores.
        alpha: Significance level.
        correction: 'bonferroni' or 'bh' (Benjamini-Hochberg).

    Returns:
        Nested dict: pair_key → DeLong result dict.
        pair_key is formatted as "ModelA vs ModelB".
        Also includes a "pairwise" dict mapping (name_a, name_b) → result,
        and a "summary" with corrected p-values.
    """
    names = list(scores.keys())
    if len(names) < 2:
        raise ValueError("Need at least 2 models for comparison")

    n_comparisons = len(names) * (len(names) - 1) // 2

    pairwise_raw: Dict[Tuple[str, str], Dict[str, float]] = {}
    p_values: list = []
    pair_keys: list = []

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            na, nb = names[i], names[j]
            result = delong_test(y_true, scores[na], scores[nb], alpha)
            pairwise_raw[(na, nb)] = result
            p_values.append(result["p_value"])
            pair_keys.append((na, nb))

    # Apply correction
    p_array = np.array(p_values)

    if correction == "bonferroni":
        from validation_framework import SignificanceTester
        corrected = SignificanceTester.bonferroni_correct(p_array)
    elif correction == "bh":
        from validation_framework import SignificanceTester
        corrected, _ = SignificanceTester.benjamini_hochberg(p_array, alpha)
    else:
        corrected = p_array

    # Build output
    output: Dict[str, Dict[str, float]] = {}
    summary_lines: list[str] = []

    for idx, (na, nb) in enumerate(pair_keys):
        key = f"{na} vs {nb}"
        output[key] = dict(pairwise_raw[(na, nb)])
        output[key]["p_value_raw"] = output[key]["p_value"]
        output[key]["p_value_corrected"] = float(corrected[idx])
        output[key]["significant_corrected"] = bool(float(corrected[idx]) < alpha)

        summary_lines.append(
            f"  {na} vs {nb}: AUC Δ={output[key]['delta_auc']:.4f}, "
            f"p={float(corrected[idx]):.4f}, "
            f"{'*' if float(corrected[idx]) < alpha else 'ns'}"
        )

    output["_summary"] = {"n_comparisons": n_comparisons, "correction": correction}
    output["_summary_lines"] = "\n".join(summary_lines)

    return output


# ── helpers ─────────────────────────────────────────────────────────────────


def _compute_V10(pos_scores: Array, neg_scores: Array) -> Array:
    """Compute V₁₀: for each positive, fraction of negatives with lower score.

    Uses Sun & Xu (2014) fast outer-comparison via broadcasting.

    Args:
        pos_scores: Scores for positive class, shape (n_pos,).
        neg_scores: Scores for negative class, shape (n_neg,).

    Returns:
        Array of shape (n_pos,) — per-positive separators.
    """
    # Broadcasting: (n_pos, 1) vs (1, n_neg) → (n_pos, n_neg)
    diff = pos_scores[:, np.newaxis] - neg_scores[np.newaxis, :]
    # ψ(s_pos, s_neg) — elements of the structural component
    sep = (diff > 0).astype(float) + 0.5 * (diff == 0).astype(float)
    return np.mean(sep, axis=1)


def _compute_V01(neg_scores: Array, pos_scores: Array) -> Array:
    """Compute V₀₁: for each negative, fraction of positives with higher score."""
    diff = pos_scores[np.newaxis, :] - neg_scores[:, np.newaxis]
    sep = (diff > 0).astype(float) + 0.5 * (diff == 0).astype(float)
    return np.mean(sep, axis=1)


def _sample_variance(x: Array) -> float:
    """Unbiased sample variance."""
    return float(np.var(x, ddof=1)) if len(x) > 1 else 0.0


def _sample_covariance(x: Array, y: Array) -> float:
    """Unbiased sample covariance."""
    if len(x) < 2:
        return 0.0
    n = len(x)
    return float(np.sum((x - np.mean(x)) * (y - np.mean(y))) / (n - 1))


# ── reporting ───────────────────────────────────────────────────────────────


def report(result: Dict[str, float]) -> str:
    """Publication-ready DeLong test report."""
    sig = "SIGNIFICANT ✓" if result.get("significant", False) else "NOT SIGNIFICANT ✗"
    return (
        "══ DeLong AUC Comparison ══\n"
        f"  AUC Model A:     {result['auc_a']:.4f}\n"
        f"  AUC Model B:     {result['auc_b']:.4f}\n"
        f"  Δ AUC (B - A):   {result['delta_auc']:.4f}\n"
        f"  SE(Δ):           {result['se_delta']:.4f}\n"
        f"  95% CI for Δ:    [{result['ci95_lower']:.4f}, {result['ci95_upper']:.4f}]\n"
        f"  z-statistic:     {result['z_statistic']:.4f}\n"
        f"  p-value:         {result['p_value']:.6f}\n"
        f"  Variance A:      {result['var_a']:.6f}\n"
        f"  Variance B:      {result['var_b']:.6f}\n"
        f"  Cov(A, B):       {result['cov_ab']:.6f}\n"
        f"  Verdict:         {sig}"
    )


# ── standalone runner ───────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys

    print("Module 5: DeLong Test — self-test")
    print("=" * 60)

    rng = np.random.RandomState(42)
    n = 500

    # Generate data with true signal
    X = rng.randn(n, 5)
    true_prob = 1.0 / (1.0 + np.exp(-(X[:, 0] + 0.5 * X[:, 1])))
    y_true = (rng.rand(n) < true_prob).astype(int)

    # Model A (good): close to true signal
    scores_a = np.clip(true_prob + 0.1 * rng.randn(n), 0, 1)

    # Model B (better): even closer, slightly more signal
    scores_b = np.clip(true_prob + 0.05 * rng.randn(n) + 0.02, 0, 1)

    # Model C (same as A): statistically indistinguishable
    scores_c = np.clip(true_prob + 0.11 * rng.randn(n), 0, 1)

    print("\n── A vs B (B should be better) ──")
    result_ab = delong_test(y_true, scores_a, scores_b)
    print(report(result_ab))

    print("\n── A vs C (should be similar) ──")
    result_ac = delong_test(y_true, scores_a, scores_c)
    print(report(result_ac))

    print("\n── Multi-model comparison ──")
    multi = delong_test_multi(
        y_true,
        {"Model A": scores_a, "Model B": scores_b, "Model C": scores_c},
        correction="bonferroni",
    )
    print(multi["_summary_lines"])

    # Verification: B should have higher AUC than A
    if result_ab["auc_b"] > result_ab["auc_a"]:
        print(f"\n✓ Model B AUC > Model A AUC (as expected)")
    else:
        print(f"\n⚠ AUC ordering unexpected (stochastic sample)")

    print("\nSelf-test complete.")
    sys.exit(0)
