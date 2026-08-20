#!/usr/bin/env python3
"""Compute decision curve + per-specificity operating table for the 627-sample
cross-study cohort, using the same per-seed OOF predictions as the fusion
ablation. Self-contained: re-runs the tumor-naive LR-on-PCA inside CV folds,
then writes a JSON summary a clinician can read.

Outputs results/decision_curve_627.json with:
    {
      "per_strategy": {
        "tumor_naive": {
          "decision_curve": {"thresholds": [...], "nb_model": [...],
                              "nb_all": [...], "nb_none": [...],
                              "clinical_value_range": [lo, hi]},
          "per_specificity": [
            {"specificity": 0.80, "sensitivity": 0.95, "operating_threshold": 0.42},
            ...
          ]
        },
        "naive_average": { ... },
        ...
      },
      "metadata": {n_samples, n_cancer, n_healthy, seeds}
    }
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Optional

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

if __package__ in (None, ""):
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    from fragmentomics.tumor_naive_adapter import (
        load_cohort, load_labels_tsv, CHANNEL_NAMES,
    )
    from fragmentomics.decision_curve import (
        decision_curve, per_specificity_table,
    )
    from fragmentomics.fusion_ablation import (
        _simulate_mutation_scores, _summarize, _harmonize,
        _MUT_CALIB,
    )
else:
    from .tumor_naive_adapter import load_cohort, load_labels_tsv, CHANNEL_NAMES
    from .decision_curve import decision_curve, per_specificity_table
    from .fusion_ablation import (
        _simulate_mutation_scores, _summarize, _harmonize, _MUT_CALIB,
    )


def _tumor_naive_score_oof(X: np.ndarray, y: np.ndarray,
                             study: np.ndarray, pca_n: int, seeds: int,
                             ) -> np.ndarray:
    """OOF tumor-naive scores across `seeds` runs, concatenated per seed.

    Returns shape (seeds * n, ) — each seed gives a full OOF prediction.
    Caller averages / pools as needed.
    """
    out = []
    for s in range(seeds):
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=s)
        per_seed = np.zeros(len(y))
        for tr, te in cv.split(X, y):
            Xtr, sc = _harmonize(X[tr], study[tr], None)
            Xte, _ = _harmonize(X[te], study[te], sc)
            max_pca = min(Xtr.shape[0], Xtr.shape[1])
            m = LogisticRegression(max_iter=2000).fit(
                PCA(n_components=min(pca_n, max_pca)).fit(Xtr).transform(Xtr), y[tr])
            per_seed[te] = m.predict_proba(
                PCA(n_components=min(pca_n, max_pca)).fit(Xtr).transform(Xte)
            )[:, 1]
        out.append(per_seed)
    return np.concatenate(out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features-dir", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--out", required=True,
                    help="Output JSON path")
    ap.add_argument("--pca-n", type=int, default=200)
    ap.add_argument("--seeds", type=int, default=5,
                    help="Number of CV seeds (default 5; pooled across seeds "
                         "for stable estimates)")
    ap.add_argument("--target-auc", type=float,
                    default=_MUT_CALIB["target_auc"])
    ap.add_argument("--thresholds", default="0.01,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.40,0.50",
                    help="Comma-separated decision-curve thresholds")
    args = ap.parse_args()

    labels = load_labels_tsv(args.labels)
    study_of: dict[str, str] = {}
    with open(args.labels) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                study_of[parts[0]] = parts[2]
            else:
                study_of[parts[0]] = "default"
    sample_ids = sorted(labels)
    X, order = load_cohort(sample_ids, args.features_dir,
                           channels=CHANNEL_NAMES, skip_missing=True)
    y = np.asarray([labels[s] for s in order], dtype=int)
    study = np.asarray([study_of.get(s, "default") for s in order])
    print(f"[decision_curve] loaded {X.shape[0]} samples "
          f"({(y==1).sum()} cancer / {(y==0).sum()} healthy)",
          file=sys.stderr)

    # OOF tumor-naive scores (5 seeds × n samples, then average).
    tn_per_seed = _tumor_naive_score_oof(X, y, study, args.pca_n, args.seeds)
    tn_score = tn_per_seed.reshape(args.seeds, -1).mean(axis=0)
    print(f"[decision_curve] pooled tumor-naive OOF AUC: "
          f"{_summarize(y, tn_score)['auc']:.4f}", file=sys.stderr)

    # Synthetic mutation score (matches DeepCatch calibration).
    rng = np.random.default_rng(int(_MUT_CALIB["seed"]))
    mut_score = _simulate_mutation_scores(y, rng,
                                          target_auc=args.target_auc)
    navg = (tn_score + mut_score) / 2.0

    thresholds = np.asarray([float(x) for x in args.thresholds.split(",")])
    payload = {
        "metadata": {
            "n_samples": int(X.shape[0]),
            "n_cancer": int((y == 1).sum()),
            "n_healthy": int((y == 0).sum()),
            "seeds": args.seeds,
            "pca_n": args.pca_n,
            "thresholds": thresholds.tolist(),
            "mutation_score_target_auc": args.target_auc,
        },
        "per_strategy": {},
    }
    for name, scores in [("tumor_naive", tn_score),
                          ("mutation_only", mut_score),
                          ("naive_average", navg)]:
        payload["per_strategy"][name] = {
            "auc": _summarize(y, scores)["auc"],
            "decision_curve": decision_curve(y, scores, thresholds=thresholds),
            "per_specificity": per_specificity_table(
                y, scores, specificities=[0.80, 0.85, 0.90, 0.95, 0.98, 0.99]),
        }
        print(f"[decision_curve] {name}: AUC="
              f"{payload['per_strategy'][name]['auc']:.4f}", file=sys.stderr)

    with open(args.out, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"[decision_curve] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())