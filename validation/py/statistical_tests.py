"""
Statistical test implementations used throughout the DeepCatch pipeline.

All tests match the Node.js implementations in:
  validation/node/realHeadToHead.js  (AUC, DeLong, bootstrap)
  validation/node/realCET.js         (AUC computation)

References:
  DeLong ER, DeLong DM, Clarke-Pearson DL (1988) Biometrics 44:837-845
  Efron B, Tibshirani RJ (1993) "An Introduction to the Bootstrap"
"""

import numpy as np
from scipy.stats import norm
from typing import Tuple, List, Optional, Dict


def compute_auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """
    Compute AUC via trapezoidal integration (identical to JS implementation).

    Args:
        scores: Array of prediction scores (higher = more likely positive).
        labels: Binary array (1 = positive, 0 = negative).

    Returns:
        AUC (float, 0.5–1.0). Returns 0.5 if either class is empty.

    References:
        Fawcett T (2006) Pattern Recognit Lett 27:861-874
    """
    n = len(scores)
    if n == 0:
        return 0.5

    n_pos = int(np.sum(labels))
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5

    # Sort by descending score
    order = np.argsort(-scores)
    sorted_labels = labels[order]
    sorted_scores = scores[order]

    auc = 0.0
    prev_fpr = 0.0
    prev_tpr = 0.0
    tp = 0
    fp = 0

    for i in range(n):
        if sorted_labels[i]:
            tp += 1
        else:
            fp += 1

        # Emit point at score boundary
        if i == n - 1 or sorted_scores[i] != sorted_scores[i + 1]:
            tpr = tp / n_pos
            fpr = fp / n_neg
            auc += (fpr - prev_fpr) * (tpr + prev_tpr) / 2.0
            prev_fpr = fpr
            prev_tpr = tpr

    return float(auc)


def bootstrap_auc(y_true: np.ndarray, y_score: np.ndarray,
                  n_bootstrap: int = 2000, alpha: float = 0.05,
                  rng: Optional[np.random.RandomState] = None) -> Dict:
    """
    Stratified bootstrap confidence intervals for AUC.

    Args:
        y_true: Binary labels.
        y_score: Prediction scores.
        n_bootstrap: Number of bootstrap iterations.
        alpha: Significance level (default 0.05 → 95% CI).
        rng: Seeded RandomState for reproducibility.

    Returns:
        dict with keys: point, lo, hi, se, n_boot
    """
    if rng is None:
        rng = np.random.RandomState(42)

    n = len(y_true)
    aucs = np.empty(n_bootstrap)

    for b in range(n_bootstrap):
        idx = rng.randint(0, n, size=n)
        # Reject degenerate bootstraps (all-same-class)
        bs_labels = y_true[idx]
        bs_scores = y_score[idx]
        if bs_labels.sum() == 0 or bs_labels.sum() == len(bs_labels):
            aucs[b] = 0.5
        else:
            aucs[b] = compute_auc(bs_scores, bs_labels)

    point = compute_auc(y_score, y_true)
    aucs.sort()
    lo_idx = int(np.floor(n_bootstrap * alpha / 2))
    hi_idx = int(np.floor(n_bootstrap * (1 - alpha / 2)))
    lo = aucs[max(0, lo_idx)]
    hi = aucs[min(n_bootstrap - 1, hi_idx)]
    se = np.std(aucs, ddof=1)

    return {'point': point, 'lo': lo, 'hi': hi, 'se': se, 'n_boot': n_bootstrap}


# Alias for backward compatibility
bootstrap_ci = bootstrap_auc


def delong_test(y_true: np.ndarray, scores_a: np.ndarray,
                scores_b: np.ndarray, z_threshold: float = 5.0) -> Dict:
    """
    DeLong's test for comparing two correlated AUCs.

    Tests H0: AUC_A = AUC_B  vs  H1: AUC_A ≠ AUC_B

    Implementation follows DeLong et al. (1988) exactly as in
    realHeadToHead.js, computing V10 and V01 structural components.

    Args:
        y_true: Binary ground-truth labels.
        scores_a: Scores from classifier A.
        scores_b: Scores from classifier B.
        z_threshold: Clamp |z| to avoid numerical overflow.

    Returns:
        dict with: auc_a, auc_b, delta_auc, z, p_value, significant, se
    """
    n = len(y_true)
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    n_pos = len(pos_idx)
    n_neg = len(neg_idx)

    if n_pos == 0 or n_neg == 0:
        return {
            'auc_a': 0.5, 'auc_b': 0.5, 'delta_auc': 0.0,
            'z': 0.0, 'p_value': 1.0, 'significant': False, 'se': float('inf'),
        }

    auc_a = compute_auc(scores_a, y_true)
    auc_b = compute_auc(scores_b, y_true)

    # V10: For each negative, proportion of positives with higher score
    v10_a = np.zeros(n_neg)
    v10_b = np.zeros(n_neg)
    for j, idx in enumerate(neg_idx):
        v10_a[j] = np.mean(scores_a[pos_idx] > scores_a[idx]) + \
                   0.5 * np.mean(scores_a[pos_idx] == scores_a[idx])
        v10_b[j] = np.mean(scores_b[pos_idx] > scores_b[idx]) + \
                   0.5 * np.mean(scores_b[pos_idx] == scores_b[idx])

    # V01: For each positive, proportion of negatives with lower score
    v01_a = np.zeros(n_pos)
    v01_b = np.zeros(n_pos)
    for j, idx in enumerate(pos_idx):
        v01_a[j] = np.mean(scores_a[idx] > scores_a[neg_idx]) + \
                   0.5 * np.mean(scores_a[idx] == scores_a[neg_idx])
        v01_b[j] = np.mean(scores_b[idx] > scores_b[neg_idx]) + \
                   0.5 * np.mean(scores_b[idx] == scores_b[neg_idx])

    # Covariance
    if n_neg > 1:
        s10_12 = np.sum((v10_a - auc_a) * (v10_b - auc_b)) / (n_neg - 1)
        var_a_neg = np.sum((v10_a - auc_a) ** 2) / ((n_neg - 1) * n_neg) if n_neg > 1 else 0
        var_b_neg = np.sum((v10_b - auc_b) ** 2) / ((n_neg - 1) * n_neg) if n_neg > 1 else 0
    else:
        s10_12 = 0.0
        var_a_neg = 0.0
        var_b_neg = 0.0

    if n_pos > 1:
        s01_12 = np.sum((v01_a - auc_a) * (v01_b - auc_b)) / (n_pos - 1)
        var_a_pos = np.sum((v01_a - auc_a) ** 2) / ((n_pos - 1) * n_pos) if n_pos > 1 else 0
        var_b_pos = np.sum((v01_b - auc_b) ** 2) / ((n_pos - 1) * n_pos) if n_pos > 1 else 0
    else:
        s01_12 = 0.0
        var_a_pos = 0.0
        var_b_pos = 0.0

    var_a = var_a_neg + var_a_pos
    var_b = var_b_neg + var_b_pos
    cov_ab = s10_12 / max(1, n_neg) + s01_12 / max(1, n_pos)

    se_diff = np.sqrt(max(0, var_a + var_b - 2 * cov_ab))
    if se_diff > 0:
        z = (auc_a - auc_b) / se_diff
    else:
        z = 0.0

    # Clamp z to avoid overflow
    z = np.clip(z, -z_threshold, z_threshold)
    p_value = 2.0 * (1.0 - norm.cdf(abs(z)))

    return {
        'auc_a': auc_a, 'auc_b': auc_b, 'delta_auc': auc_a - auc_b,
        'z': float(z), 'p_value': float(p_value),
        'significant': bool(p_value < 0.05), 'se': float(se_diff),
    }


def _normal_cdf(x: float) -> float:
    """Accurate normal CDF using Abramowitz & Stegun approximation."""
    return float(norm.cdf(x))


def sensitivity_at_specificity(scores: np.ndarray, labels: np.ndarray,
                              target_spec: float = 0.99) -> float:
    """
    Compute sensitivity at a given specificity threshold.

    Identical algorithm to realHeadToHead.js sensitivityAtSpecificity().
    """
    n_pos = int(np.sum(labels))
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.0

    order = np.argsort(-scores)
    sorted_labels = labels[order]
    sorted_scores = scores[order]

    best_sens = 0.0
    tp = 0
    fp = 0

    for i in range(len(scores)):
        if sorted_labels[i]:
            tp += 1
        else:
            fp += 1
        if i == len(scores) - 1 or sorted_scores[i] != sorted_scores[i + 1]:
            spec = 1.0 - fp / n_neg
            sens = tp / n_pos
            if spec >= target_spec:
                best_sens = max(best_sens, sens)

    return best_sens


def bonferroni_correct(p_values: List[float]) -> List[float]:
    """
    Bonferroni correction: p_adj = min(1, p_i * n).

    Safest (most conservative) multiple-testing correction.
    """
    n = len(p_values)
    return [min(1.0, p * n) for p in p_values]


def benjamini_hochberg(p_values: List[float], fdr: float = 0.05) -> List[bool]:
    """
    Benjamini-Hochberg FDR correction.

    Returns: List[bool] indicating rejected null hypotheses.

    References:
        Benjamini & Hochberg (1995) JRSSB 57:289-300
    """
    n = len(p_values)
    if n == 0:
        return []

    idx_order = np.argsort(p_values)
    sorted_p = np.array(p_values)[idx_order]
    thresholds = (np.arange(1, n + 1) / n) * fdr

    # Find max k such that p(k) <= k/n * FDR
    reject_mask = sorted_p <= thresholds
    if np.any(reject_mask):
        k_star = np.max(np.where(reject_mask)[0])
        reject_sorted = np.zeros(n, dtype=bool)
        reject_sorted[:k_star + 1] = True
    else:
        reject_sorted = np.zeros(n, dtype=bool)

    # Re-order back to original
    reject = np.zeros(n, dtype=bool)
    reject[idx_order] = reject_sorted
    return reject.tolist()


# ── Demo ──────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    rng = np.random.RandomState(42)
    n = 200
    y_true = np.concatenate([np.ones(100), np.zeros(100)])

    # Model A: signal ~ cancer, some noise
    scores_a = np.where(y_true, 0.6 + rng.normal(0, 0.2, n), 0.3 + rng.normal(0, 0.2, n))
    scores_a = np.clip(scores_a, 0, 1)

    # Model B: slightly better signal
    scores_b = np.where(y_true, 0.65 + rng.normal(0, 0.2, n), 0.25 + rng.normal(0, 0.2, n))
    scores_b = np.clip(scores_b, 0, 1)

    auc_a = compute_auc(scores_a, y_true)
    auc_b = compute_auc(scores_b, y_true)
    print(f"AUC A = {auc_a:.4f}, AUC B = {auc_b:.4f}")

    boot = bootstrap_auc(y_true, scores_b, n_bootstrap=1000, rng=rng)
    print(f"Bootstrap CI (B): {boot['point']:.4f} [{boot['lo']:.4f}, {boot['hi']:.4f}]")

    delong = delong_test(y_true, scores_b, scores_a)
    print(f"DeLong: ΔAUC={delong['delta_auc']:.4f}, z={delong['z']:.2f}, "
          f"p={delong['p_value']:.4f}, sig={delong['significant']}")

    sens99 = sensitivity_at_specificity(scores_b, y_true, 0.99)
    print(f"Sensitivity @ 99% spec: {sens99:.3f}")

    # Multiple testing demo
    p_vals = [0.001, 0.04, 0.0001, 0.5, 0.02]
    print(f"Bonferroni: {bonferroni_correct(p_vals)}")
    print(f"BH rejected: {benjamini_hochberg(p_vals)}")
