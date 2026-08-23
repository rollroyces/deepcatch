"""
Decision curve analysis + per-specificity operating table.

Reference: Vickers & Elkin (2006), "Decision Curve Analysis: A Novel Method
for Evaluating Prediction Models", Medical Decision Making 26:565-574.

Net benefit at threshold p_t:
    NB(p_t) = (TP/N) - (FP/N) * (p_t / (1 - p_t))
            = sensitivity * prevalence - (1 - specificity) * (1 - prevalence) * p_t / (1 - p_t)

Treat-all baseline:
    NB_all(p_t) = prevalence - (1 - prevalence) * p_t / (1 - p_t)
Treat-none baseline:
    NB_none = 0

A model has clinical value at threshold p_t iff NB_model > max(NB_all, NB_none).
The output is a table of (threshold, NB_model, NB_all, NB_none) and the
range of thresholds where the model beats both baselines.

This module is self-contained — no framework dependencies — so it's
trivial to unit-test and safe to call from the CLI without setup.
"""
from __future__ import annotations

import json
from typing import Optional

import numpy as np


def net_benefit(y_true: np.ndarray, y_score: np.ndarray,
                threshold: float) -> tuple[float, float, float]:
    """Net benefit at threshold p_t for (model, treat-all, treat-none).

    Treats y_score as the model's *predicted probability of being positive*.
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    n = len(y_true)
    if n == 0:
        return 0.0, 0.0, 0.0
    prev = float(y_true.mean())
    if prev <= 0 or prev >= 1:
        # Degenerate cohort (no positives or no negatives): undefined.
        return 0.0, 0.0, 0.0
    pred_pos = y_score >= threshold
    tp = float(((pred_pos == 1) & (y_true == 1)).sum())
    fp = float(((pred_pos == 1) & (y_true == 0)).sum())
    # NB = TP/N - FP/N * w,  where w = p_t/(1 - p_t).
    w = threshold / (1.0 - threshold) if threshold < 1.0 else float("inf")
    nb_model = tp / n - (fp / n) * w
    nb_all = prev - (1.0 - prev) * w
    nb_none = 0.0
    return nb_model, nb_all, nb_none


def decision_curve(y_true: np.ndarray, y_score: np.ndarray,
                   thresholds: Optional[np.ndarray] = None) -> dict:
    """Compute the full decision curve.

    Returns dict with arrays `thresholds`, `nb_model`, `nb_all`, `nb_none`,
    and a `clinical_value_range` field listing the (lo, hi) threshold range
    where the model beats both baselines.
    """
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    if thresholds is None:
        # Fine grid from 0.01 to 0.50 (above 0.5 the NB weight explodes
        # and the model rarely helps; below 0.01 the treat-all baseline
        # dominates).
        thresholds = np.linspace(0.01, 0.50, 50)
    nb_model, nb_all, nb_none = [], [], []
    for t in thresholds:
        a, b, c = net_benefit(y_true, y_score, float(t))
        nb_model.append(a); nb_all.append(b); nb_none.append(c)
    nb_model = np.asarray(nb_model)
    nb_all = np.asarray(nb_all)
    nb_none = np.asarray(nb_none)
    # Clinical-value range: thresholds where model beats BOTH baselines.
    beats = (nb_model > nb_all) & (nb_model > nb_none)
    if beats.any():
        idx = np.where(beats)[0]
        clinical_value_range = (float(thresholds[idx[0]]),
                                float(thresholds[idx[-1]]))
    else:
        clinical_value_range = None
    return {
        "thresholds": thresholds.tolist(),
        "nb_model": nb_model.tolist(),
        "nb_all": nb_all.tolist(),
        "nb_none": nb_none.tolist(),
        "clinical_value_range": clinical_value_range,
    }


def per_specificity_table(y_true: np.ndarray, y_score: np.ndarray,
                          specificities: Optional[list[float]] = None
                          ) -> list[dict]:
    """Sensitivity at a fixed list of specificities.

    A clinician picks an operating point by setting the specificity they
    require (e.g. 95% for screening, 99% for confirmatory). This is the
    table they'd look at to decide.

    Returns list of {specificity, sensitivity, threshold} dicts.
    """
    from sklearn.metrics import roc_curve
    if specificities is None:
        specificities = [0.80, 0.85, 0.90, 0.95, 0.98, 0.99]
    fpr, tpr, thr = roc_curve(y_true, y_score)
    # fpr = 1 - specificity; we want sens at >= specified specificity,
    # so read where fpr <= 1 - specificity.
    out: list[dict] = []
    for sp in specificities:
        target_fpr = 1.0 - sp
        idx = np.where(fpr <= target_fpr)[0]
        if len(idx):
            sens = float(tpr[idx[-1]])
            op_thr = float(thr[idx[-1]])
        else:
            sens = 0.0
            op_thr = float("nan")
        out.append({"specificity": sp, "sensitivity": sens,
                    "operating_threshold": op_thr})
    return out


def write_json_summary(out_path: str, y_true: np.ndarray,
                      score_name_to_arr: dict[str, np.ndarray]) -> dict:
    """One-call helper: compute decision curve + specificity table for each
    named score array, dump JSON.

    Returns the JSON-serializable dict.
    """
    payload: dict = {"per_strategy": {}}
    for name, scores in score_name_to_arr.items():
        dc = decision_curve(y_true, scores)
        ps = per_specificity_table(y_true, scores)
        payload["per_strategy"][name] = {
            "decision_curve": dc,
            "per_specificity": ps,
        }
    if out_path:
        with open(out_path, "w") as f:
            json.dump(payload, f, indent=2)
    return payload