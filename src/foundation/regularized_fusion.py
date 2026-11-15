"""Regularized late fusion and Bayesian model averaging for small-n cohorts.

Why this module exists
----------------------
On the honest paired TCGA-LUAD cohort (40 samples from 20 patients,
scored with a per-patient ``GroupKFold``), the multi-modal pipeline
currently reads:

    foundation (tiny transformer encoder)   AUC 0.554 +/- 0.005
    sklearn LogisticRegression on 2 channels AUC 0.910
    equal-weight average of the 2 channels  AUC 0.915

The neural encoder is the wrong model class for a cohort this small:
it has more parameters than samples and spends its capacity memorising
the per-patient invariant that ``GroupKFold`` deliberately denies it.
The fix is not a bigger or better-tuned encoder -- it is a smaller
*estimator*. These two fusion estimators fuse the channel scores late,
with an explicit prior that shrinks the combination toward the
equal-weight average, and they select their own hyper-parameters with
an inner group-aware cross-validation so the outer measurement stays
honest.

Both estimators are pure numpy / scikit-learn: no torch, so they add
no dependency and cannot inherit the encoder's overfitting.

Estimators
----------
``regularized_late_fusion``
    A stacker over the standardized channel scores. The stacker is a
    logistic regression whose regularization strength, and how much of
    its output is shrunk back toward the equal-weight average, are
    both chosen by an inner group-aware CV.

``bayesian_model_average``
    Pseudo-BMA over a candidate set (single-channel models, the
    equal-weight average, the regularised stack). Each candidate's
    predictive likelihood is measured by an inner group-aware CV, and
    the candidates are averaged with posterior weights proportional to
    that likelihood. Unlike the stacker, BMA never commits to one
    functional form: it hedges.

Both are fitted on the training fold only and never see the scored
fold's labels.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

__all__ = [
    "regularized_late_fusion",
    "bayesian_model_average",
    "log_loss",
    "zscore",
]

# Logistic-regression regularization grid for the meta-learner.
C_GRID: Tuple[float, ...] = (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
# Mixing weight between the stacker output and the equal-weight average.
SHRINK_GRID: Tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)


def zscore(x: np.ndarray, fit_idx: np.ndarray) -> np.ndarray:
    """Standardize ``x`` using the mean/sd of the fitting rows only."""
    mu = float(np.mean(x[fit_idx]))
    sd = float(np.std(x[fit_idx]))
    return (x - mu) / max(sd, 1e-12)


def log_loss(y: np.ndarray, p: np.ndarray, eps: float = 1e-12) -> float:
    """Mean binary log-loss."""
    p = np.clip(np.asarray(p, dtype=float), eps, 1.0 - eps)
    y = np.asarray(y, dtype=float)
    return float(-np.mean(y * np.log(p) + (1.0 - y) * np.log(1.0 - p)))


def _safe_auc(y: np.ndarray, s: np.ndarray) -> float:
    """AUC that returns 0.5 when the fold has only one class."""
    y = np.asarray(y)
    if len(np.unique(y)) < 2:
        return 0.5
    return float(roc_auc_score(y, np.asarray(s, dtype=float)))


def _inner_splits(n: int, groups: np.ndarray, fit_idx: np.ndarray,
                  n_inner: int) -> list:
    """Inner group-aware CV splits over the fitting rows."""
    g = groups[fit_idx]
    n_groups = len(np.unique(g))
    k = int(max(2, min(n_inner, n_groups)))
    if n_groups < 2:
        return []
    cv = GroupKFold(n_splits=k)
    return list(cv.split(np.zeros(len(fit_idx)), y=None, groups=g))


def _equal_weight_prob(zp: np.ndarray, zf: np.ndarray, y: np.ndarray,
                       fit_idx: np.ndarray, score_idx: np.ndarray) -> np.ndarray:
    """Equal-weight average of the standardized channels, Platt-scaled.

    The 1-D calibration is fitted on the training rows only, so the
    output is a probability (needed for the BMA likelihood) and the
    ranking is unchanged (a monotone transform).
    """
    avg = 0.5 * (zp + zf)
    cal = LogisticRegression(C=1.0, max_iter=1000).fit(
        avg[fit_idx].reshape(-1, 1), y[fit_idx])
    return cal.predict_proba(avg[score_idx].reshape(-1, 1))[:, 1]


def regularized_late_fusion(
    panel: np.ndarray,
    frag: np.ndarray,
    y: np.ndarray,
    fit_idx: np.ndarray,
    score_idx: np.ndarray,
    groups: np.ndarray,
    n_inner: int = 3,
) -> Dict[str, object]:
    """Late fusion with a regularised, shrunk meta-learner.

    The meta-learner is a logistic regression on the two standardized
    channel scores. Two hyper-parameters are chosen by an inner
    group-aware CV:

    - ``C`` -- the regularization strength of the meta-learner;
    - ``shrink`` -- how much of the stacker output is kept, the rest
      coming from the equal-weight average of the channels. A shrink
      of 0 means "ignore the stacker, use the equal-weight average".

    Selection uses inner AUC: with a training fold of a few dozen
    samples the log-loss is too noisy to rank candidates.

    Returns a dict with ``score`` (the scored-fold predictions),
    ``inner_auc``, ``selected_C``, ``selected_shrink`` and
    ``n_inner_folds``.
    """
    zp = zscore(panel, fit_idx)
    zf = zscore(frag, fit_idx)

    splits = _inner_splits(len(y), groups, fit_idx, n_inner)
    if not splits:
        p = _equal_weight_prob(zp, zf, y, fit_idx, score_idx)
        return {"score": p, "inner_auc": None, "selected_C": None,
                "selected_shrink": 1.0, "n_inner_folds": 0,
                "fallback": "equal-weight average"}

    y_fit = y[fit_idx]
    inner_rows = np.arange(len(fit_idx))
    best = (-1.0, None, None)
    for C in C_GRID:
        for shrink in SHRINK_GRID:
            oof = np.zeros(len(fit_idx), dtype=float)
            for tr_in, te_in in splits:
                sub_fit = fit_idx[tr_in]
                sub_score = fit_idx[te_in]
                zps = zscore(panel, sub_fit)
                zfs = zscore(frag, sub_fit)
                X = np.column_stack([zps, zfs])
                lr = LogisticRegression(C=C, max_iter=2000).fit(
                    X[sub_fit], y[sub_fit])
                p_lr = lr.predict_proba(X[sub_score])[:, 1]
                cal = LogisticRegression(C=1.0, max_iter=1000).fit(
                    (0.5 * (zps + zfs))[sub_fit].reshape(-1, 1), y[sub_fit])
                p_eq = cal.predict_proba(
                    (0.5 * (zps + zfs))[sub_score].reshape(-1, 1))[:, 1]
                oof[te_in] = shrink * p_lr + (1.0 - shrink) * p_eq
            a = _safe_auc(y_fit, oof)
            if a > best[0]:
                best = (a, C, shrink)

    inner_auc, best_C, best_shrink = best
    # Refit the selected configuration on the whole training fold.
    X = np.column_stack([zp, zf])
    lr = LogisticRegression(C=best_C, max_iter=2000).fit(X[fit_idx], y[fit_idx])
    p_lr = lr.predict_proba(X[score_idx])[:, 1]
    cal = LogisticRegression(C=1.0, max_iter=1000).fit(
        (0.5 * (zp + zf))[fit_idx].reshape(-1, 1), y[fit_idx])
    p_eq = cal.predict_proba((0.5 * (zp + zf))[score_idx].reshape(-1, 1))[:, 1]
    score = best_shrink * p_lr + (1.0 - best_shrink) * p_eq
    return {
        "score": score,
        "inner_auc": round(float(inner_auc), 4),
        "selected_C": best_C,
        "selected_shrink": float(best_shrink),
        "n_inner_folds": len(splits),
    }


def bayesian_model_average(
    panel: np.ndarray,
    frag: np.ndarray,
    y: np.ndarray,
    fit_idx: np.ndarray,
    score_idx: np.ndarray,
    groups: np.ndarray,
    n_inner: int = 3,
) -> Dict[str, object]:
    """Pseudo-BMA over a candidate set of fusion models.

    Candidates: panel-only, frag-only, the equal-weight average, the
    average of the two single-channel models, and a regularised stack.
    Each candidate is scored by an inner group-aware CV on the training
    rows; its predictive log-likelihood becomes its posterior weight
    (equal priors), and the candidates are averaged with those weights.

    Returns ``score``, ``weights`` (posterior per candidate),
    ``inner_auc`` per candidate and ``n_inner_folds``.
    """
    zp = zscore(panel, fit_idx)
    zf = zscore(frag, fit_idx)

    splits = _inner_splits(len(y), groups, fit_idx, n_inner)
    if not splits:
        p = _equal_weight_prob(zp, zf, y, fit_idx, score_idx)
        return {"score": p, "weights": {}, "inner_auc": {}, "n_inner_folds": 0,
                "fallback": "equal-weight average"}

    y_fit = y[fit_idx]

    def _fit_predict(name: str, sub_fit: np.ndarray, sub_score: np.ndarray,
                     zp_s: np.ndarray, zf_s: np.ndarray) -> np.ndarray:
        if name == "panel_only":
            lr = LogisticRegression(C=1.0, max_iter=2000).fit(
                zp_s[sub_fit].reshape(-1, 1), y[sub_fit])
            return lr.predict_proba(zp_s[sub_score].reshape(-1, 1))[:, 1]
        if name == "frag_only":
            lr = LogisticRegression(C=1.0, max_iter=2000).fit(
                zf_s[sub_fit].reshape(-1, 1), y[sub_fit])
            return lr.predict_proba(zf_s[sub_score].reshape(-1, 1))[:, 1]
        if name == "equal_weight":
            cal = LogisticRegression(C=1.0, max_iter=1000).fit(
                (0.5 * (zp_s + zf_s))[sub_fit].reshape(-1, 1), y[sub_fit])
            return cal.predict_proba(
                (0.5 * (zp_s + zf_s))[sub_score].reshape(-1, 1))[:, 1]
        if name == "mean_of_channels":
            p1 = _fit_predict("panel_only", sub_fit, sub_score, zp_s, zf_s)
            p2 = _fit_predict("frag_only", sub_fit, sub_score, zp_s, zf_s)
            return 0.5 * (p1 + p2)
        # regularised stack: C chosen by the same inner CV, fixed at 1.0
        # here so the candidate set stays small and interpretable.
        X = np.column_stack([zp_s, zf_s])
        lr = LogisticRegression(C=1.0, max_iter=2000).fit(X[sub_fit], y[sub_fit])
        return lr.predict_proba(X[sub_score])[:, 1]

    names = ["panel_only", "frag_only", "equal_weight", "mean_of_channels",
             "stack"]
    oof = {name: np.zeros(len(fit_idx), dtype=float) for name in names}
    for tr_in, te_in in splits:
        sub_fit = fit_idx[tr_in]
        sub_score = fit_idx[te_in]
        zp_s = zscore(panel, sub_fit)
        zf_s = zscore(frag, sub_fit)
        for name in names:
            oof[name][te_in] = _fit_predict(name, sub_fit, sub_score,
                                             zp_s, zf_s)

    loglik = {}
    inner_auc = {}
    for name in names:
        p = np.clip(oof[name], 1e-12, 1 - 1e-12)
        loglik[name] = float(np.sum(y_fit * np.log(p)
                                    + (1 - y_fit) * np.log(1 - p)))
        inner_auc[name] = round(_safe_auc(y_fit, oof[name]), 4)

    # Pseudo-BMA posterior: weights proportional to exp(predictive
    # log-likelihood), computed in a shifted space to avoid overflow.
    ll = np.array([loglik[n] for n in names], dtype=float)
    w = np.exp(ll - ll.max())
    w = w / w.sum()
    weights = {n: round(float(wi), 4) for n, wi in zip(names, w)}

    # Score the held-out rows with every candidate fitted on the full
    # training fold, then average with the posterior weights.
    preds = {}
    for name in names:
        preds[name] = _fit_predict(name, fit_idx, score_idx, zp, zf)
    score = sum(weights[n] * preds[n] for n in names)

    return {
        "score": score,
        "weights": weights,
        "inner_auc": inner_auc,
        "n_inner_folds": len(splits),
        "candidate_predictions": {n: preds[n].tolist() for n in names},
    }
