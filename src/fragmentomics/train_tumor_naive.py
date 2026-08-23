#!/usr/bin/env python3
"""
Tumor-naive cfDNA classifier — DeepCatch wrapper around the
cfdna-fragmentomics-pipeline pre-computed artifacts.

Reads the pipeline's per-sample `.npy` and `.fsd.json` outputs, runs
the same harmonized-PCA + Logistic Regression pipeline used by
`scripts/train_classifier.py` in the standalone repo, and reports the
multi-seed honest result.

Mirrors the pipeline's honest benchmark methodology:
- 5 random seeds for stratified 5-fold CV (configurable)
- Per-study z-score harmonization inside each fold (no leakage)
- PCA inside each fold (no leakage)
- Pooled out-of-fold predictions → single ROC / AUC / fixed-spec sens
- Mean +/- std across seeds reported for AUC

The "tumor-naive" label means: no mutation information is used. This
is the channel that complements DeepCatch's mutation-informed detection
(see real_tcga_validation.py).

Usage:
    python -m deepcatch.fragmentomics.train_tumor_naive \\
        --features-dir ../cfdna-fragmentomics-pipeline/data/features \\
        --labels ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \\
        --pca-n 200 --seeds 5

Outputs JSON to --out.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List, Optional

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

# Allow running as `python src/fragmentomics/train_tumor_naive.py` too
if __package__ in (None, ""):
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
    from fragmentomics.tumor_naive_adapter import (
        load_cohort, load_labels_tsv, CHANNEL_NAMES, CHANNEL_DIMS,
    )
else:
    from .tumor_naive_adapter import (
        load_cohort, load_labels_tsv, CHANNEL_NAMES, CHANNEL_DIMS,
    )


# ---------- study harmonization (fit-on-train, apply-to-test) ----------

def _harmonize(X: np.ndarray, study: np.ndarray,
               scalers: Optional[dict] = None) -> tuple[np.ndarray, dict]:
    """Per-study z-score: fit StandardScaler on train, apply to test.

    For the test fold we pass the *train* scalers dict, so test-set
    statistics never leak into the transformation.
    """
    if scalers is None:
        scalers = {}
        for st in np.unique(study):
            mask = study == st
            if mask.sum() > 1:
                scalers[st] = StandardScaler().fit(X[mask])
    out = np.empty_like(X, dtype=float)
    for st, sc in scalers.items():
        mask = study == st
        if mask.any():
            out[mask] = sc.transform(X[mask])
    return out, scalers


# ---------- core CV ----------

def _evaluate_seed(X: np.ndarray, y: np.ndarray, study: np.ndarray,
                   pca_n: int, seed: int, harmonize: bool = True) -> dict:
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    yt: list[int] = []
    ys: list[float] = []
    for tr, te in cv.split(X, y):
        if harmonize:
            Xtr, scalers = _harmonize(X[tr], study[tr], None)
            Xte, _ = _harmonize(X[te], study[te], scalers)
        else:
            sc = StandardScaler().fit(X[tr])
            Xtr = sc.transform(X[tr]); Xte = sc.transform(X[te])
        # PCA inside fold → no test leakage
        max_pca = min(Xtr.shape[0], Xtr.shape[1])
        pca = PCA(n_components=min(pca_n, max_pca)).fit(Xtr)
        model = LogisticRegression(max_iter=2000).fit(
            pca.transform(Xtr), y[tr])
        ys.extend(model.predict_proba(pca.transform(Xte))[:, 1].tolist())
        yt.extend(y[te].tolist())
    yt_arr = np.asarray(yt); ys_arr = np.asarray(ys)
    auc = float(roc_auc_score(yt_arr, ys_arr))
    fpr, tpr, _ = roc_curve(yt_arr, ys_arr)
    def _sens_at(target: float) -> float:
        idx = np.where(fpr <= target)[0]
        return float(tpr[idx[-1]]) if len(idx) else 0.0
    return {
        "seed": seed,
        "auc": auc,
        "sens_at_95": _sens_at(0.05),
        "sens_at_99": _sens_at(0.01),
    }


# ---------- CLI ----------

def _parse_study(labels_tsv: str) -> dict[str, str]:
    """Read `{sample}\\t{label}[\\t{study}]` → sample→study string.

    Tolerates both 2-column and 3-column label files.
    """
    out: dict[str, str] = {}
    with open(labels_tsv) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                out[parts[0]] = parts[2]
            else:
                # default to a single study
                out[parts[0]] = "default"
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--features-dir", required=True,
                    help="Directory containing pipeline artifacts "
                         "(.delfi_*.npy and .fsd.json)")
    ap.add_argument("--labels", required=True,
                    help="TSV `{sample}\\t{label}` or `{sample}\\t{label}\\t{study}`")
    ap.add_argument("--out", default=None,
                    help="Output JSON path (default: stdout)")
    ap.add_argument("--pca-n", type=int, default=200,
                    help="PCA components inside each fold (default 200)")
    ap.add_argument("--seeds", type=int, default=5,
                    help="Number of random seeds for stratified K-fold "
                         "(default 5; honest benchmark methodology)")
    ap.add_argument("--channels", default=",".join(CHANNEL_NAMES),
                    help="Comma-separated subset of channels to use; "
                         f"default: {','.join(CHANNEL_NAMES)}")
    ap.add_argument("--skip-missing", action="store_true", default=True,
                    help="Skip samples missing pipeline artifacts instead "
                         "of raising (default: True — labels files often "
                         "contain a few orphans after data-quality filters)")
    ap.add_argument("--strict", dest="skip_missing", action="store_false",
                    help="Disable --skip-missing; raise on the first "
                         "missing artifact instead")
    ap.add_argument("--no-harmonize", action="store_true",
                    help="Skip per-study z-score harmonization (only valid "
                         "for single-study cohorts)")
    args = ap.parse_args()

    channels = [c.strip() for c in args.channels.split(",") if c.strip()]
    for c in channels:
        if c not in CHANNEL_DIMS:
            ap.error(f"unknown channel {c!r}; valid: {list(CHANNEL_DIMS)}")

    labels = load_labels_tsv(args.labels)
    study_of = _parse_study(args.labels)
    sample_ids = sorted(labels)
    print(f"[adapter] loading {len(sample_ids)} samples from {args.features_dir}",
          file=sys.stderr)
    X, order = load_cohort(sample_ids, args.features_dir,
                           channels=channels,
                           skip_missing=args.skip_missing)
    if X.shape[0] == 0:
        print("[adapter] no samples loaded — check --features-dir",
              file=sys.stderr)
        return 1
    print(f"[adapter] loaded X.shape={X.shape}; "
          f"{len(order)}/{len(sample_ids)} samples had artifacts",
          file=sys.stderr)
    if X.shape[0] != len(order):
        print("[adapter] FATAL: matrix length vs order length mismatch",
              file=sys.stderr)
        return 2
    y = np.asarray([labels[s] for s in order], dtype=int)
    study = np.asarray([study_of.get(s, "default") for s in order])

    print(f"[adapter] class balance: "
          f"{(y == 1).sum()} cancer / {(y == 0).sum()} healthy; "
          f"studies={sorted(set(study))}", file=sys.stderr)

    # Multi-seed CV
    seeds = list(range(args.seeds))
    per_seed = [_evaluate_seed(X, y, study, args.pca_n, s,
                               harmonize=not args.no_harmonize)
                for s in seeds]
    aucs = np.asarray([r["auc"] for r in per_seed])
    s95s = np.asarray([r["sens_at_95"] for r in per_seed])
    s99s = np.asarray([r["sens_at_99"] for r in per_seed])

    result = {
        "n_samples": int(X.shape[0]),
        "n_cancer": int((y == 1).sum()),
        "n_healthy": int((y == 0).sum()),
        "studies": sorted(set(study)),
        "channels": channels,
        "feature_dim": int(X.shape[1]),
        "pca_n": args.pca_n,
        "n_seeds": args.seeds,
        "harmonized": not args.no_harmonize,
        "auc_mean": float(aucs.mean()),
        "auc_std": float(aucs.std()),
        "sens_at_95_mean": float(s95s.mean()),
        "sens_at_95_std": float(s95s.std()),
        "sens_at_99_mean": float(s99s.mean()),
        "sens_at_99_std": float(s99s.std()),
        "per_seed": per_seed,
    }
    payload = json.dumps(result, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(payload + "\n")
        print(f"[adapter] wrote {args.out}", file=sys.stderr)
    print(payload)
    return 0


if __name__ == "__main__":
    sys.exit(main())