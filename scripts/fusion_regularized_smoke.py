#!/usr/bin/env python3
"""Honest evaluation of regularized late fusion and Bayesian model averaging.

This is the measurement the multi-modal layer needs: the foundation
model currently scores AUC 0.554 on the paired TCGA-LUAD cohort while a
plain logistic regression on the same two channels scores 0.910. This
script runs the *same* protocol as ``scripts/foundation_real_smoke.py``
-- real TCGA-derived panel scores, per-patient ``GroupKFold``,
pair-broken shuffled-label control -- and adds two late-fusion
estimators that are appropriate for a cohort this small:

- ``regularized_late_fusion``: a regularised, shrunk stacker over the
  standardized channel scores;
- ``bayesian_model_average``: pseudo-BMA over a candidate set of
  fusion models, weighted by inner-CV predictive likelihood.

The script also re-computes the LR baseline on the identical folds, so
the new numbers are directly comparable with the published baseline
rather than with a number from a different run.

Usage:
    python scripts/fusion_regularized_smoke.py --seeds 3 --out-json \\
        results/fusion_regularized_smoke.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.model_selection import GroupKFold

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))

from src.foundation.label_control import (  # noqa: E402
    assert_not_degenerate,
    foldwise_permuted_labels,
)
from src.foundation.regularized_fusion import (  # noqa: E402
    bayesian_model_average,
    regularized_late_fusion,
)

VARIANTS = (
    "lr_baseline",
    "naive_average",
    "regularized_late_fusion",
    "bayesian_model_average",
)


def _auc(y, s) -> float:
    y = np.asarray(y)
    if len(np.unique(y)) < 2:
        return 0.5
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(y, s))


def _sens_at_spec(y, s, spec: float) -> float:
    fpr, tpr, _ = roc_curve(y, s)
    idx = np.where(fpr <= (1.0 - spec))[0]
    return float(tpr[idx[-1]]) if len(idx) else 0.0


# The label-permutation controls live in ``src/foundation/label_control.py``
# so that every script shares ONE implementation. That matters here
# because this file previously carried its own copy of the permutation,
# while ``scripts/foundation_real_smoke.py`` carried a DIFFERENT one
# whose construction was a no-op: it permuted the positive index set
# among itself and then wrote 1s back to those same positions, so
# ``y_shuf == y_true`` for every seed and the "shuffled" AUCs it
# reported were the honest AUCs. A control implemented twice drifts;
# a control implemented once cannot.
#
# MEASURED NULL ON THIS COHORT: 0.399 +/- 0.118 (200 permutations),
# NOT 0.5. The reason is structural: the paired cohort gives every
# patient one positive and one negative arm whose panel scores differ
# by +127 on average (100% of patients have the positive arm higher).
# A model fitted to random labels still recovers that within-pair axis,
# because it is the direction along which the paired points differ,
# and the resulting predictions are anti-correlated with the permuted
# labels they are scored against.
#
# So a control below 0.5 is NOT on its own evidence of leakage on this
# design; it has to be compared with the measured null.


def _arm_swap_labels(y: np.ndarray, seed: int) -> np.ndarray:
    """Swap the two arms of a random half of the patients.

    The paired cohort has one positive and one negative arm per
    patient. This control keeps the class balance but re-pairs half
    the patients, so a model that leans on the per-patient invariant
    (the mutation-derived score that is constant within a pair) scores
    *below* 0.5 instead of above it. It is the sharper test of the
    exact leakage Audit-2 was worried about.
    """
    rng = np.random.default_rng(seed + 3131)
    n_pat = len(y) // 2
    out = y.copy()
    for i in range(n_pat):
        if rng.random() < 0.5:
            out[i], out[n_pat + i] = out[n_pat + i], out[i]
    return out


def _measure_null(panel, frag, y, groups, folds, seed, trials) -> dict:
    """Empirical null of the permuted-label control, measured not assumed.

    The paired cohort does NOT give a null centred at 0.5: every patient
    contributes one positive and one negative arm whose panel scores
    differ by ~+127 (the positive arm is higher for 100% of patients),
    and a model fitted to random labels still recovers that within-pair
    axis. The resulting predictions are anti-correlated with the
    permuted labels, so the null sits below 0.5. Measuring it is the
    only honest way to read the control.
    """
    X = np.column_stack([panel, frag])
    aucs = []
    for t in range(trials):
        yp = foldwise_permuted_labels(y, folds, seed + 10000 * (t + 1))
        oof = np.zeros(len(y))
        for tr, te in folds:
            m = LogisticRegression(C=1.0, max_iter=1000).fit(X[tr], yp[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
        aucs.append(_auc(yp, oof))
    a = np.asarray(aucs)
    return {
        "trials": int(trials),
        "mean": round(float(a.mean()), 4),
        "std": round(float(a.std()), 4),
        "q05": round(float(np.percentile(a, 5)), 4),
        "q95": round(float(np.percentile(a, 95)), 4),
        "fraction_below_0.5": round(float((a < 0.5).mean()), 3),
    }


def _run_seed(panel, frag, y, groups, seed, n_folds, n_inner):
    """Outer GroupKFold over the two channel scores."""
    n = len(y)
    n_groups = len(np.unique(groups))
    cv = GroupKFold(n_splits=int(min(n_folds, n_groups)))
    oof = {v: np.zeros(n, dtype=float) for v in VARIANTS}
    diag = {"selected_C": [], "selected_shrink": [], "bma_weights": []}
    folds = list(cv.split(np.zeros(n), y, groups=groups))
    y_perm = foldwise_permuted_labels(y, folds, seed)
    # A control that reproduces the real labels is not a control: it
    # reports the honest AUC back as if it were a null. Raise rather
    # than publish a meaningless number.
    assert_not_degenerate(y, y_perm)
    y_swap = _arm_swap_labels(y, seed)
    oof_perm = {v: np.zeros(n, dtype=float) for v in VARIANTS}
    oof_swap = {v: np.zeros(n, dtype=float) for v in VARIANTS}

    for tr, te in folds:
        # --- honest labels ---
        oof["naive_average"][te] = (panel[te] + frag[te]) / 2.0
        Xtr = np.column_stack([panel[tr], frag[tr]])
        Xte = np.column_stack([panel[te], frag[te]])
        lr = LogisticRegression(C=1.0, max_iter=1000).fit(Xtr, y[tr])
        oof["lr_baseline"][te] = lr.predict_proba(Xte)[:, 1]
        rlf = regularized_late_fusion(panel, frag, y, tr, te, groups,
                                      n_inner=n_inner)
        oof["regularized_late_fusion"][te] = rlf["score"]
        diag["selected_C"].append(rlf["selected_C"])
        diag["selected_shrink"].append(rlf["selected_shrink"])
        bma = bayesian_model_average(panel, frag, y, tr, te, groups,
                                      n_inner=n_inner)
        oof["bayesian_model_average"][te] = bma["score"]
        diag["bma_weights"].append(bma["weights"])

        # --- control 1: globally permuted labels (null) ---
        oof_perm["naive_average"][te] = (panel[te] + frag[te]) / 2.0
        lr_p = LogisticRegression(C=1.0, max_iter=1000).fit(Xtr, y_perm[tr])
        oof_perm["lr_baseline"][te] = lr_p.predict_proba(Xte)[:, 1]
        oof_perm["regularized_late_fusion"][te] = regularized_late_fusion(
            panel, frag, y_perm, tr, te, groups, n_inner=n_inner)["score"]
        oof_perm["bayesian_model_average"][te] = bayesian_model_average(
            panel, frag, y_perm, tr, te, groups, n_inner=n_inner)["score"]

        # --- control 2: arms swapped within half the patients ---
        oof_swap["naive_average"][te] = (panel[te] + frag[te]) / 2.0
        lr_s = LogisticRegression(C=1.0, max_iter=1000).fit(Xtr, y_swap[tr])
        oof_swap["lr_baseline"][te] = lr_s.predict_proba(Xte)[:, 1]
        oof_swap["regularized_late_fusion"][te] = regularized_late_fusion(
            panel, frag, y_swap, tr, te, groups, n_inner=n_inner)["score"]
        oof_swap["bayesian_model_average"][te] = bayesian_model_average(
            panel, frag, y_swap, tr, te, groups, n_inner=n_inner)["score"]

    out = {}
    for v in VARIANTS:
        out[v] = {
            "auc": round(_auc(y, oof[v]), 4),
            "sens_at_95": round(_sens_at_spec(y, oof[v], 0.95), 4),
            "sens_at_99": round(_sens_at_spec(y, oof[v], 0.99), 4),
            "permuted_label_auc": round(_auc(y_perm, oof_perm[v]), 4),
            "arm_swap_auc": round(_auc(y_swap, oof_swap[v]), 4),
        }
    return out, diag


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--n-folds", type=int, default=5)
    ap.add_argument("--n-inner", type=int, default=3)
    ap.add_argument("--null-trials", type=int, default=50,
                    help="permutations used to measure the control's "
                         "empirical null (0 = skip)")
    ap.add_argument("--n-patients", type=int, default=20)
    ap.add_argument("--tumor-fraction", type=float, default=0.001)
    ap.add_argument("--cfdna-depth", type=int, default=5000)
    ap.add_argument("--bg-error-rate", type=float, default=0.002)
    ap.add_argument("--out-json", default="results/fusion_regularized_smoke.json")
    args = ap.parse_args()

    from foundation_real_smoke import _try_load_real_panel_scores

    seeds = list(range(args.seeds))
    print(f"[1/3] building real TCGA-derived panel scores "
          f"({args.n_patients} patients, {len(seeds)} seeds) ...", flush=True)
    per_seed = _try_load_real_panel_scores(
        n_patients=args.n_patients, seeds=seeds,
        tumor_fraction=args.tumor_fraction, cfdna_depth=args.cfdna_depth,
        bg_error_rate=args.bg_error_rate)
    if per_seed is None:
        print("ERROR: could not build the real panel scores "
              "(TCGA cache missing?)")
        return 1

    print("[2/3] running the honest protocol ...", flush=True)
    results = {}
    for seed in seeds:
        d = per_seed[seed]
        res, diag = _run_seed(d["panel_scores"], d["frag_scores"], d["y_true"],
                               d["patient_groups"], seed, args.n_folds,
                               args.n_inner)
        if args.null_trials:
            folds = list(GroupKFold(n_splits=min(args.n_folds, len(
                np.unique(d["patient_groups"])))).split(
                    np.zeros(len(d["y_true"])), d["y_true"],
                    groups=d["patient_groups"]))
            diag["empirical_null"] = _measure_null(
                d["panel_scores"], d["frag_scores"], d["y_true"],
                d["patient_groups"], folds, seed, args.null_trials)
            print(f"      seed {seed}  measured null (LR baseline, "
                  f"{args.null_trials} permutations): "
                  f"{diag['empirical_null']['mean']:.4f} "
                  f"+/- {diag['empirical_null']['std']:.4f}")
        results[str(seed)] = {"metrics": res, "diagnostics": diag}
        for v in VARIANTS:
            print(f"      seed {seed}  {v:<26} AUC={res[v]['auc']:.4f}  "
                  f"perm={res[v]['permuted_label_auc']:.4f} "
                  f"swap={res[v]['arm_swap_auc']:.4f}")

    summary = {}
    for v in VARIANTS:
        aucs = [results[str(s)]["metrics"][v]["auc"] for s in seeds]
        shuf = [results[str(s)]["metrics"][v]["permuted_label_auc"]
                for s in seeds]
        swap = [results[str(s)]["metrics"][v]["arm_swap_auc"] for s in seeds]
        summary[v] = {
            "auc_mean": round(float(np.mean(aucs)), 4),
            "auc_std": round(float(np.std(aucs)), 4),
            "permuted_label_auc_mean": round(float(np.mean(shuf)), 4),
            "arm_swap_auc_mean": round(float(np.mean(swap)), 4),
            "per_seed_auc": aucs,
        }

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "cohort": "TCGA-LUAD paired (20 patients x 2 arms), real MAF-derived "
                      "panel scores and real mutation-derived fragmentomics "
                      "channel; simulated plasma at TF="
                      f"{args.tumor_fraction}",
            "cv": f"GroupKFold(n_splits=min({args.n_folds}, n_patients)) over "
                  "patient groups -- both arms of every patient stay together",
            "controls": {
                "permuted_label_auc": "labels permuted within each fold's "
                                      "train and test blocks independently; "
                                      "on this paired design the measured "
                                      "null is ~0.42, not 0.5 (see "
                                      "per_seed[*].diagnostics.empirical_null)",
                "arm_swap_auc": "arms swapped within half the patients; a "
                                "model leaning on the per-patient invariant "
                                "lands below 0.5",
            },
            "seeds": seeds,
            "n_inner_folds": args.n_inner,
        },
        "summary": summary,
        "per_seed": results,
        "baseline_reference": {
            "source": "scripts/foundation_real_smoke.py (Audit-2 honest numbers)",
            "foundation_auc": 0.554,
            "lr_baseline_auc": 0.910,
            "naive_average_auc": 0.915,
        },
        "honest_framing": "Same data, same folds, same control as the foundation "
                          "smoke: the LR baseline is re-computed here so the "
                          "comparison is within-run. A fusion estimator that "
                          "beats the LR baseline on this cohort is evidence "
                          "about the estimator class, not about clinical "
                          "performance -- n=40 samples from 20 patients, "
                          "simulated plasma, one cancer type.",
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"[3/3] wrote {args.out_json}")
    print()
    for v in VARIANTS:
        s = summary[v]
        print(f"  {v:<26} AUC {s['auc_mean']:.4f} +/- {s['auc_std']:.4f}   "
              f"permuted {s['permuted_label_auc_mean']:.4f}   "
              f"arm-swap {s['arm_swap_auc_mean']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
