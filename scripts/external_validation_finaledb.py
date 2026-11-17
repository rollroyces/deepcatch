#!/usr/bin/env python3
"""Pure external validation: train on FinaleDB 6+8, test on held-out studies.

This is the external-validation counterpart of
``scripts/cross_study_loso.py``. The LOSO benchmark holds out a
publication *inside* the cached cohort; this script goes further and
evaluates on samples DeepCatch has never touched in any form, fetched
from the public FinaleDB API and S3 bucket by
``scripts/fetch_external_cohort.py``:

  train : FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019),
          the 627-sample cohort the cross-study benchmark uses
  test  : publications 1 (Snyder 2016), 7 (Sun 2019),
          9 (Adalsteinsson 2017) -- blood plasma only, never trained on

Protocol
--------
Identical to the LOSO benchmark, so the numbers are comparable:

- per-study z-score harmonization, each publication normalized by its
  own feature statistics (the "new site normalizes itself" assumption);
- PCA (n_components = min(200, n_train)) fit on the training fold only;
- LogisticRegression(max_iter=2000) fit on the training fold only;
- the held-out publication is scored once, never fitted on.

The harmonization is applied on BOTH sides, and the two sides are not
interchangeable. On the training side it is per study, because the
training cohort pools two studies (6 and 8) and the in-cohort benchmark
harmonises per study inside every fold -- the true-confound control
shows what happens otherwise (without it the classifier learns "which
study is this from" at AUC 0.999). On the held-out side it is per
publication, because the held-out cohort pools several studies and
z-scoring it as one pool lets the largest of them anchor the others.

This code used to z-score the held-out cohort as a SINGLE pool, which
contradicted the protocol above and the LOSO benchmark it claims to
mirror. It is recorded in the artifact's ``harmonization_sensitivity``
block, together with the two other harmonisations, because the
difference is not small.

Reported: pooled external AUC with a DeLong 95% CI, sensitivity at the
screening-grade specificity grid, per-publication AUC, and a
per-cancer/disease breakdown.

Usage:
    python scripts/external_validation_finaledb.py \
        --external-dir data/external_cohort \
        --out-json results/external_validation_finaledb.json \
        --out-md   docs/EXTERNAL_VALIDATION.md
"""
from __future__ import annotations
import sys as _sys
from pathlib import Path as _Path

# The repository's own path layout, resolved from this file: the
# sibling pipeline, the feature cache and the labels are looked up
# next to the repository, then in the home directory, with the
# environment variables overriding both.
_sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "scripts"))
from _paths import SIBLING_PIPELINE  # noqa: E402


import argparse
import json
import os
import sys
from datetime import datetime, timezone

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)
sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))

from src.per_cancer_sens_at_spec import (  # noqa: E402
    DEFAULT_SPECIFICITIES,
    SCREENING_SPECIFICITIES,
    delong_auc_ci,
)

SISTER = os.environ.get(
    "DEEPCATCH_SISTER_PIPELINE", str(SIBLING_PIPELINE))
TRAIN_PUBLICATIONS = {"6", "8"}

# The canonical label vocabulary lives in the package, so every reader of
# a labels file in this repository resolves labels the same way.
_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from fragmentomics.labels import (  # noqa: E402
    looks_like_header, resolve_columns, resolve_label,
)


def _read_labels(labels_path: str, publications: set[str]):
    """Read the labels file into (labels, studies, disease_class, publication).

    Split out of `_load_cohort` so the label reading can be tested
    WITHOUT the sibling pipeline's checkout: the cohort loader imports
    `honest_benchmark.load5`, which is absent in CI and in a fresh
    clone, while the label vocabulary is not.
    """
    labels, studies, disease_class, publication = {}, {}, {}, {}
    with open(labels_path) as f:
        headered = None
        label_col = study_col = None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if headered is None:
                # The header is decided ONCE, from the first row, and the
                # column indices come from it.
                headered = looks_like_header(parts)
                label_col, study_col = resolve_columns(parts, headered=headered)
                if headered:
                    continue
            if label_col is None or label_col >= len(parts):
                continue
            study = parts[study_col] if study_col is not None and study_col < len(parts) else ""
            if not study or study in ("study", ""):
                continue
            s = parts[0]
            pub = parts[4] if len(parts) >= 5 and parts[4] not in ("publication",) \
                else {"jiang": "6", "cristiano": "8"}.get(study, "")
            if pub not in publications:
                continue
            labels[s] = resolve_label(parts[label_col], s)
            studies[s] = study
            disease_class[s] = parts[1] if len(parts) > 1 else ""
            publication[s] = pub
    return labels, studies, disease_class, publication


def _load_cohort(features_dir: str, labels_path: str, publications: set[str]):
    """Load the 5-channel features + labels for the given publications.

    Uses the sibling pipeline's ``load5`` so the training and external
    features are assembled by exactly the same code path.
    """
    sys.path.insert(0, SISTER)
    sys.path.insert(0, os.path.join(SISTER, "scripts"))
    from honest_benchmark import load5

    labels, studies, disease_class, publication = _read_labels(
        labels_path, publications)
    X, y, st = load5(labels, studies, feat_dir=features_dir)
    order = [s for s in sorted(labels)
             if all(os.path.exists(os.path.join(
                 features_dir, f"{s}.{c}")) for c in
                 ("delfi_5mb_ratio.npy", "delfi_5mb_coverage.npy",
                  "delfi_100kb_ratio.npy", "delfi_100kb_counts.npy"))
             and os.path.exists(os.path.join(features_dir, f"{s}.fsd.json"))]
    if len(order) != X.shape[0]:
        raise RuntimeError(f"row alignment failed: {len(order)} vs {X.shape[0]}")
    disease = np.array([disease_class.get(s, "UNKNOWN") for s in order])
    return X, np.asarray(y, dtype=int), np.asarray(st), disease, order


def _load_external(external_dir: str):
    """Load the externally fetched cohort from its manifest + features."""
    manifest_path = os.path.join(external_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"{manifest_path} not found -- run scripts/fetch_external_cohort.py "
            "first")
    with open(manifest_path) as f:
        manifest = json.load(f)
    samples = manifest.get("samples", [])
    feat_dir = os.path.join(external_dir, "features")
    keep = [s for s in samples
            if os.path.exists(os.path.join(feat_dir, f"{s['sample']}.fsd.json"))]
    labels = {s["sample"]: int(s["label"]) for s in keep}
    studies = {s["sample"]: str(s["publication"]) for s in keep}
    X, y, st = _load5_from(feat_dir, labels, studies)
    # ``load5`` emits rows in sorted-sample order, so every per-sample
    # array must be built in that same order -- the manifest is not
    # sorted, and mixing the two silently mislabels the rows.
    order = sorted(labels)
    disease_by_sample = {s["sample"]: (s["disease"] or "UNKNOWN") for s in keep}
    disease = np.array([disease_by_sample[s] for s in order])
    return X, y, st, disease, keep, manifest


def _load5_from(feat_dir: str, labels: dict, studies: dict):
    sys.path.insert(0, SISTER)
    sys.path.insert(0, os.path.join(SISTER, "scripts"))
    from honest_benchmark import load5
    return load5(labels, studies, feat_dir=feat_dir)


def _zscore(X: np.ndarray) -> np.ndarray:
    """Z-score each feature column of one matrix."""
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    return (X - mu) / np.maximum(sd, 1e-12)


def _zscore_by_group(X: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Z-score each matrix within its own group.

    The harmonisation assumption is "a new site normalises itself": each
    publication is centred and scaled by ITS OWN feature statistics, so
    no publication's distribution anchors another's. Z-scoring the whole
    held-out cohort as one pool instead lets the largest publication
    anchor the others -- with 238 case-only samples from publication 9 in
    the pool, publications 1 and 7 are scored against publication 9's
    mean and variance, which is not what the protocol says and not what
    the LOSO benchmark does.

    The same applies to the training cohort, which pools two studies:
    the in-cohort benchmark harmonises per study inside every fold, and
    the true-confound control shows why (without it, the classifier
    learns "which study is this from" at AUC 0.999).
    """
    out = np.empty_like(X, dtype=float)
    for g in np.unique(groups):
        m = groups == g
        out[m] = _zscore(X[m])
    return out


def _sens_at_spec(y, s, spec: float) -> float:
    fpr, tpr, _ = roc_curve(y, s)
    idx = np.where(fpr <= (1.0 - spec))[0]
    return float(tpr[idx[-1]]) if len(idx) else 0.0


def _auc_ci(y, s) -> dict:
    d = delong_auc_ci(np.asarray(y), np.asarray(s))
    return {
        "auc": round(float(d["auc"]), 4),
        "ci_lo": round(float(d["ci_lo"]), 4),
        "ci_hi": round(float(d["ci_hi"]), 4),
        "n_pos": int(d["n_pos"]),
        "n_neg": int(d["n_neg"]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features-dir",
                    default=os.path.join(SISTER, "data/features"))
    ap.add_argument("--labels-multiclass",
                    default=os.path.join(SISTER, "labels_multiclass.tsv"))
    ap.add_argument("--external-dir", default="data/external_cohort")
    ap.add_argument("--pca", type=int, default=200)
    ap.add_argument("--out-json", default="results/external_validation_finaledb.json")
    ap.add_argument("--out-md", default="docs/EXTERNAL_VALIDATION.md")
    args = ap.parse_args()

    print("[1/4] loading the training cohort (FinaleDB 6+8) ...", flush=True)
    Xtr, ytr, sttr, dtr, _ = _load_cohort(
        args.features_dir, args.labels_multiclass, TRAIN_PUBLICATIONS)
    print(f"      n={len(ytr)} (cancer {int((ytr==1).sum())}, "
          f"healthy {int((ytr==0).sum())}), features={Xtr.shape[1]}")

    print("[2/4] loading the external cohort ...", flush=True)
    Xte, yte, stte, dte, samples, manifest = _load_external(args.external_dir)
    print(f"      n={len(yte)} (cancer {int((yte==1).sum())}, "
          f"non-cancer {int((yte==0).sum())}), from {sorted(set(stte.tolist()))}")

    # ---- harmonize each publication by its own statistics ----
    # Per study on the training side (two studies pooled), per publication
    # on the held-out side. See _zscore_by_group.
    Xtr_z = _zscore_by_group(Xtr, sttr)
    Xte_z = _zscore_by_group(Xte, stte)

    print("[3/4] fitting on the training cohort and scoring the external one ...",
          flush=True)
    n_components = int(min(args.pca, Xtr_z.shape[0], Xtr_z.shape[1]))
    pca = PCA(n_components=n_components, random_state=0).fit(Xtr_z)
    model = LogisticRegression(max_iter=2000).fit(pca.transform(Xtr_z), ytr)
    score = model.predict_proba(pca.transform(Xte_z))[:, 1]

    pooled = _auc_ci(yte, score)
    sens = {f"{s:g}": round(_sens_at_spec(yte, score, s), 4)
            for s in tuple(DEFAULT_SPECIFICITIES) + tuple(SCREENING_SPECIFICITIES)}

    # ---- threshold transferred from the TRAINING cohort ----
    # A held-out publication that carries only cancer samples (no
    # controls) cannot contribute to an AUC, but it can still be scored
    # at a threshold fixed in training -- which is how a case-only
    # external cohort is normally reported.
    from cross_study_finallydb import pooled_oof
    ytr_oof, score_tr_oof, _ = pooled_oof(Xtr, ytr, sttr, list(range(5)),
                                           args.pca, harmonize=True)
    thr_spec = 0.99
    fpr, tpr, thr = roc_curve(ytr, score_tr_oof)
    idx = np.where(fpr <= (1.0 - thr_spec))[0]
    train_threshold = float(thr[idx[-1]]) if len(idx) else float("nan")
    train_spec_at_thr = round(float(1.0 - fpr[idx[-1]]), 4) if len(idx) else None

    # ---- publications with both classes -> AUC; case-only -> sensitivity ----
    per_pub = {}
    matched_mask = np.zeros(len(yte), dtype=bool)
    for pub in sorted(set(stte.tolist())):
        m = stte == pub
        classes = set(yte[m].tolist())
        if len(classes) < 2:
            pos = int((yte[m] == 1).sum())
            det = round(float((score[m] >= train_threshold).mean()), 4) if pos else None
            per_pub[pub] = {
                "n": int(m.sum()), "n_cancer": pos, "n_non_cancer": 0,
                "auc": None,
                "reason": "case-only publication (no non-cancer controls in "
                          "this study); reported as sensitivity at the "
                          "training-derived threshold",
                "sensitivity_at_train_threshold_99spec": det,
            }
            continue
        matched_mask |= m
        per_pub[pub] = {"n": int(m.sum()), **_auc_ci(yte[m], score[m])}

    # Pooled AUC over the publications that carry both classes (the only
    # combination where an AUC is meaningful).
    if matched_mask.any():
        pooled_matched = _auc_ci(yte[matched_mask], score[matched_mask])
    else:
        pooled_matched = None

    # ---- harmonisation sensitivity -------------------------------------
    # The protocol says each publication is normalised by its own
    # statistics. It is worth showing what the other two harmonisations
    # give, because the difference is not small: z-scoring the held-out
    # cohort as ONE pool lets publication 9 (238 case-only samples)
    # anchor publications 1 and 7.
    harmonization_sensitivity = {}
    for tag, tr_z, te_z in (
        ("train_global_test_global",
         _zscore(Xtr), _zscore(Xte)),
        ("train_global_test_per_publication",
         _zscore(Xtr), _zscore_by_group(Xte, stte)),
        ("train_per_study_test_per_publication",
         _zscore_by_group(Xtr, sttr), _zscore_by_group(Xte, stte)),
    ):
        pca_s = PCA(n_components=n_components, random_state=0).fit(tr_z)
        lr_s = LogisticRegression(max_iter=2000).fit(pca_s.transform(tr_z), ytr)
        sc_s = lr_s.predict_proba(pca_s.transform(te_z))[:, 1]
        mm = np.zeros(len(yte), dtype=bool)
        per_pub_s = {}
        for pub in sorted(set(stte.tolist())):
            m = stte == pub
            if len(set(yte[m].tolist())) > 1:
                per_pub_s[pub] = _auc_ci(yte[m], sc_s[m])
                mm |= m
        harmonization_sensitivity[tag] = {
            "per_publication": per_pub_s,
            "pooled_matched": (_auc_ci(yte[mm], sc_s[mm]) if mm.any()
                               else None),
        }

    # ---- per-disease breakdown (vs the external non-cancer controls) ----
    per_dis = {}
    healthy = yte == 0
    for dis in sorted({d for d in dte if d != "HEALTHY"}):
        idx = dte == dis
        if idx.sum() < 5 or healthy.sum() < 5:
            per_dis[dis] = {"skipped": True, "reason": f"n={int(idx.sum())}"}
            continue
        yy = np.concatenate([np.ones(idx.sum()), np.zeros(healthy.sum())])
        ss = np.concatenate([score[idx], score[healthy]])
        per_dis[dis] = {"n": int(idx.sum()), **_auc_ci(yy, ss)}

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "design": {
            "train": "FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019) "
                     f"-- n={len(ytr)}",
            "test": "FinaleDB publications " + ", ".join(
                f"{p} ({manifest.get('publications', [])})" for p in
                sorted(set(stte.tolist()))) + " -- blood plasma only, never "
            "trained on",
            "data_source": "public FinaleDB REST API + public S3 fragment "
                           "objects; no IRB, no data-use agreement",
            "harmonization": "per-study z-score on the training side and "
                              "per-publication z-score on the held-out side "
                              "(each publication normalized by its own "
                              "feature statistics -- the 'new site normalises "
                              "itself' assumption, as in the LOSO benchmark)",
            "classifier": f"PCA(n_components={n_components}) + "
                          "LogisticRegression(max_iter=2000), fitted on the "
                          "training cohort only",
            "labels": "FinaleDB disease field: a diagnosis ending in 'cancer' "
                      "-> cancer; Healthy / Hepatitis B / Cirrhosis / "
                      "inflammatory bowel disease / SLE / liver transplant -> "
                      "non-cancer control",
        },
        "cohorts": {
            "train": {"n": int(len(ytr)),
                      "n_cancer": int((ytr == 1).sum()),
                      "n_healthy": int((ytr == 0).sum())},
            "external": {"n": int(len(yte)),
                         "n_cancer": int((yte == 1).sum()),
                         "n_non_cancer": int((yte == 0).sum()),
                         "per_publication": {
                             p: int((stte == p).sum())
                             for p in sorted(set(stte.tolist()))}},
        },
        "pooled_external_auc": pooled,
        "pooled_external_auc_matched_publications_only": pooled_matched,
        "harmonization_sensitivity": harmonization_sensitivity,
        "training_derived_threshold": {
            "specificity_target": thr_spec,
            "threshold": round(train_threshold, 6),
            "train_specificity_at_threshold": train_spec_at_thr,
            "note": "threshold fitted on the training cohort's out-of-fold "
                    "scores; used to score case-only external publications",
        },
        "sens_at_spec": sens,
        "per_publication_auc": per_pub,
        "per_disease_auc": per_dis,
        "honest_framing": (
            "The test publication was never used for training, feature "
            "extraction statistics, or hyper-parameter selection. This is an "
            "external validation on public data, not a clinical validation: "
            "the non-cancer arm is small and comes from different "
            "publications than the cancer arm in places, so the pooled AUC "
            "mixes a study shift with a cohort-composition shift. Read the "
            "per-publication and per-disease rows before quoting the pooled "
            "number."),
    }

    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"[4/4] wrote {args.out_json}")

    md = [
        "# External Validation on Held-Out FinaleDB Publications\n",
        "> Train: FinaleDB 6 + 8 (the cross-study cohort). Test: publications "
        + ", ".join(str(p) for p in sorted(set(stte.tolist())))
        + " -- never trained on. Public data, no IRB required.\n",
        f"- Generated: `{payload['generated_at']}`\n",
        f"- Train n={len(ytr)} (cancer {int((ytr==1).sum())}) | "
        f"External n={len(yte)} (cancer {int((yte==1).sum())}, "
        f"non-cancer {int((yte==0).sum())})\n",
        f"- Classifier: {payload['design']['classifier']}\n",
        # The threshold the case-only publications are scored at is a
        # quoted number (the methods paper names it), so it has to appear
        # in this narrative too -- otherwise a reader cannot connect the
        # sensitivity below to an operating point.
        f"- Training-derived threshold "
        f"(99% specificity on the training OOF): "
        f"{payload['training_derived_threshold']['threshold']} "
        f"(train specificity at that threshold: "
        f"{payload['training_derived_threshold']['train_specificity_at_threshold']})\n",
        "\n## Pooled external AUC\n",
        f"**{pooled['auc']}** (DeLong 95% CI {pooled['ci_lo']}-{pooled['ci_hi']}), "
        f"n_pos={pooled['n_pos']}, n_neg={pooled['n_neg']}\n",
        # The number the manuscript quotes: the publications that carry BOTH
        # classes. Pooling the case-only publication into an AUC denominator
        # produces a figure that means nothing, so it is reported separately.
        ("\n**Pooled over the publications carrying both classes**: "
         f"**{pooled_matched['auc']}** (DeLong 95% CI {pooled_matched['ci_lo']}"
         f"-{pooled_matched['ci_hi']}), n_pos={pooled_matched['n_pos']}, "
         f"n_neg={pooled_matched['n_neg']}\n"
         if pooled_matched else "\n_Only one held-out publication carries "
         "both classes, so no matched pooled AUC is reported._\n"),
        "\n| Sens@spec | " + " | ".join(sens) + " |",
        "|---|" + "---:|" * len(sens),
        "| Sensitivity | " + " | ".join(str(sens[s]) for s in sens) + " |",
        "\n## Per-publication AUC\n",
        "| Publication | n | AUC [DeLong 95% CI] |",
        "|---|---:|---|",
    ]
    for pub, row in per_pub.items():
        if row.get("auc") is None:
            md.append(f"| {pub} | {row['n']} (case-only) | sensitivity at the "
                      f"training threshold = "
                      f"{row.get('sensitivity_at_train_threshold_99spec')} |")
        else:
            md.append(f"| {pub} | {row['n']} | {row['auc']} "
                      f"[{row['ci_lo']}, {row['ci_hi']}] |")
    md += ["\n## Per-disease AUC (vs the external non-cancer controls)\n",
           "| Disease | n | AUC [DeLong 95% CI] |", "|---|---:|---|"]
    for dis, row in per_dis.items():
        if row.get("skipped"):
            md.append(f"| {dis} | - | skipped ({row['reason']}) |")
        else:
            md.append(f"| {dis} | {row['n']} | {row['auc']} "
                      f"[{row['ci_lo']}, {row['ci_hi']}] |")
    md += ["\n## Honest framing\n", payload["honest_framing"] + "\n"]
    os.makedirs(os.path.dirname(args.out_md) or ".", exist_ok=True)
    with open(args.out_md, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"      wrote {args.out_md}")

    print()
    print(f"  POOLED EXTERNAL AUC (all publications): {pooled['auc']} "
          f"[{pooled['ci_lo']}, {pooled['ci_hi']}]")
    if pooled_matched:
        print(f"  POOLED EXTERNAL AUC (publications with both classes): "
              f"{pooled_matched['auc']} [{pooled_matched['ci_lo']}, "
              f"{pooled_matched['ci_hi']}]")
    print(f"  training-derived threshold (99% spec): {train_threshold:.4f}")
    for pub, row in per_pub.items():
        if row.get("auc") is None:
            print(f"    pub {pub}: n={row['n']} case-only, sensitivity at the "
                  f"training threshold = "
                  f"{row.get('sensitivity_at_train_threshold_99spec')}")
        else:
            print(f"    pub {pub}: n={row['n']} AUC={row['auc']} "
                  f"[{row['ci_lo']}, {row['ci_hi']}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
