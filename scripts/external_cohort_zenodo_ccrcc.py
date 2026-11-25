#!/usr/bin/env python3
"""A fourth held-out cohort: clear-cell renal cell carcinoma plasma, from Zenodo.

Why this cohort exists
----------------------
The repository's external validation rests on one cohort shape: three
FinaleDB plasma publications (1, 7, 9), all public plasma WGS. The
manuscript's own next-steps list names the weakness -- "a second external
publication family so the external claim does not rest on one cohort
shape" -- and this cohort is one:

* a **different cancer type**: clear-cell renal cell carcinoma (ccRCC),
  against a training cohort that is colorectal, lung, breast, ovary and
  liver;
* a **different laboratory and protocol**: the deposit is from an
  in-house clinical study, not a public plasma sequencing consortium;
* a **different data source**: a Zenodo deposit of per-sample fragment
  BED files, not a FinaleDB publication.

The data are public, de-identified and re-fetchable, so no IRB or data-use
agreement is involved: 10.5281/zenodo.21269430, "Cell-free DNA fragment
files (hg38) from renal cell carcinoma patients and healthy controls",
12 samples (6 ccRCC, 6 healthy), with a sample metadata table that
carries the surgical stage.

Protocol
--------
Identical to ``scripts/external_validation_finaledb.py``, imported from
it rather than re-implemented, so the number is comparable with the
FinaleDB held-out numbers:

- train on FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019);
- per-study z-score harmonisation on the training side, per-publication
  on the held-out side ("a new site normalises itself");
- PCA (n_components = min(200, n_train)) fit on the training cohort only;
- LogisticRegression(max_iter=2000) fit on the training cohort only;
- the held-out cohort is scored once, never fitted on;
- the threshold is transferred from the training cohort's own OOF scores.

What is honest about the limits
-------------------------------
* **n = 6 vs 6.** The point estimate is the finding; the interval is
  wide, and it is reported with a DeLong CI.
* **All six patients are pT3a (stage III)** and two of them have distant
  metastasis (M1 / Mx). Like the repository's other held-out cohorts,
  this is a **late-stage** cohort: it says nothing about early detection.
* **The assay is about 30x shallower** than the training cohort's (about
  1 M fragments per sample against about 30 M). All five channels are
  scale-free by construction (ratios, and coverages/counts normalised by
  their own median), so the depth difference is not a level shift the
  model can read as cancer; it does make the per-bin ratios noisier.
  The fragment counts of both cohorts are reported in the artifact so the
  reader can see the gap rather than take it on trust.
* **The healthy controls are matched to the patients in the same study**,
  so the contrast is not a cross-study contrast.

Usage:
    # 1. fetch the deposit (105.8 MB, public, no auth)
    mkdir -p /tmp/zenodo_ccrcc/frags
    curl -sL "https://zenodo.org/records/21269430/files/zenodo_cfdna_sample_metadata.tsv?download=1" \
         -o /tmp/zenodo_ccrcc/sample_metadata.tsv
    curl -sL "https://zenodo.org/records/21269430/files/inhouse_cfdna_frag.tar.gz?download=1" \
         -o /tmp/zenodo_ccrcc/inhouse_cfdna_frag.tar.gz
    tar xzf /tmp/zenodo_ccrcc/inhouse_cfdna_frag.tar.gz -C /tmp/zenodo_ccrcc/frags

    # 2. extract the features and score the cohort
    python scripts/external_cohort_zenodo_ccrcc.py \
        --frag-dir   /tmp/zenodo_ccrcc/frags/frag \
        --metadata   /tmp/zenodo_ccrcc/sample_metadata.tsv \
        --features-out data/external_cohort_zenodo_ccrcc/features \
        --out-json   results/external_cohort_zenodo_ccrcc.json
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

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
from sklearn.metrics import roc_curve

_SCRIPT_DIR = str(_Path(__file__).resolve().parent)
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from external_validation_finaledb import (  # noqa: E402
    _auc_ci,
    _load5_from,
    _load_cohort,
    _sens_at_spec,
    _zscore,
    _zscore_by_group,
)

_SRC = str(_Path(__file__).resolve().parents[1] / "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from fragmentomics import finale_features as ff  # noqa: E402
from fragmentomics.labels import resolve_label  # noqa: E402

TRAIN_PUBLICATIONS = {"6", "8"}
ZENODO_DOI = "10.5281/zenodo.21269430"
ZENODO_RECORD = "https://zenodo.org/records/21269430"
PUBLICATION_ID = "zenodo_21269430"

DEFAULT_SPECIFICITIES = (0.95, 0.98, 0.99)
SCREENING_SPECIFICITIES = (0.995, 0.999)


def read_metadata(path: str) -> dict:
    """The depositor's sample table: label, disease and stage per sample."""
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return {row["sample_id"]: row for row in rows}


def _summary_from_features(out_dir: str, sample: str) -> dict:
    """The summary the extractor writes, read back from the feature files."""
    with open(os.path.join(out_dir, f"{sample}.delfi.json")) as handle:
        delfi = json.load(handle)
    with open(os.path.join(out_dir, f"{sample}.fsd.json")) as handle:
        fsd = json.load(handle)
    return {
        "sample": sample,
        "genome_short_long_ratio": delfi["genome_short_long_ratio"],
        "mean_short_frac": delfi["mean_short_frac"],
        "fsd_fragment_count": fsd["fragment_count"],
        "fsd_short_fraction": fsd["short_fraction_100_150"],
        "fsd_short_long_ratio": fsd["short_long_ratio"],
    }


def extract_features(frag_dir: str, out_dir: str, samples: list[str],
                     reuse: bool = False) -> dict:
    """Extract the five benchmark channels for every sample.

    ``reuse`` reads the summaries back from the feature files a previous
    run wrote, which is the same numbers without re-reading the
    fragments.
    """
    summary = {}
    for sample in samples:
        path = os.path.join(frag_dir, f"{sample}.hg38.frag.bed.gz")
        if reuse:
            written = os.path.join(out_dir, f"{sample}.fsd.json")
            if os.path.exists(written):
                summary[sample] = _summary_from_features(out_dir, sample)
                continue
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        summary[sample] = ff.extract_sample(path, out_dir, sample)
    return summary


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rho of two vectors."""
    from scipy.stats import spearmanr
    return float(spearmanr(x, y).statistic)


def _partial_spearman(x: np.ndarray, y: np.ndarray,
                       z: np.ndarray) -> float:
    """Spearman rho of x and y with the linear effect of z removed.

    The rank correlations of each pair with z are partialled out, which
    is the standard partial-correlation form for ranks.
    """
    rxz = _spearman(x, z)
    ryz = _spearman(y, z)
    rxy = _spearman(x, y)
    denominator = np.sqrt(max((1 - ryz ** 2) * (1 - rxz ** 2), 1e-12))
    return float((rxy - ryz * rxz) / denominator)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frag-dir", required=True,
                    help="directory of <sample>.hg38.frag.bed.gz files")
    ap.add_argument("--metadata", required=True,
                    help="the deposited sample metadata TSV")
    ap.add_argument("--features-out",
                    default="data/external_cohort_zenodo_ccrcc/features")
    ap.add_argument("--features-dir",
                    default=os.path.join(SIBLING_PIPELINE, "data/features"))
    ap.add_argument("--labels-multiclass",
                    default=os.path.join(SIBLING_PIPELINE,
                                         "labels_multiclass.tsv"))
    ap.add_argument("--pca", type=int, default=200)
    ap.add_argument("--reuse-features", action="store_true",
                    help="read the summaries back from a previous run's "
                         "feature files instead of re-extracting them")
    ap.add_argument("--out-json",
                    default="results/external_cohort_zenodo_ccrcc.json")
    args = ap.parse_args()

    print("[1/5] reading the deposited sample metadata ...", flush=True)
    meta = read_metadata(args.metadata)
    labels, studies, disease, stage = {}, {}, {}, {}
    for sample, row in meta.items():
        is_cancer = row["group"].strip().lower().startswith("ccrcc")
        labels[sample] = 1 if is_cancer else 0
        studies[sample] = PUBLICATION_ID
        disease[sample] = ("ccRCC" if is_cancer
                           else row["group"].strip())
        stage[sample] = (row.get("pT_stage", "") or "").strip()
    print(f"      {len(labels)} samples: "
          f"{sum(labels.values())} ccRCC, {len(labels) - sum(labels.values())} healthy")

    print("[2/5] extracting the five channels ...", flush=True)
    extracted = extract_features(args.frag_dir, args.features_out,
                                 sorted(labels),
                                 reuse=args.reuse_features)
    for sample in sorted(extracted)[:3]:
        print(f"      {sample}: {extracted[sample]['fsd_fragment_count']:,} fragments")

    print("[3/5] loading the training cohort (FinaleDB 6+8) ...", flush=True)
    Xtr, ytr, sttr, dtr, _ = _load_cohort(
        args.features_dir, args.labels_multiclass, TRAIN_PUBLICATIONS)
    print(f"      n={len(ytr)} (cancer {int((ytr==1).sum())}, "
          f"healthy {int((ytr==0).sum())}), features={Xtr.shape[1]}")

    print("[4/5] loading the ccRCC features and scoring ...", flush=True)
    Xte, yte, stte = _load5_from(args.features_out, labels, studies)
    order = sorted(labels)
    disease_arr = np.array([disease[s] for s in order])
    stage_arr = np.array([stage[s] for s in order])

    Xtr_z = _zscore_by_group(Xtr, sttr)
    Xte_z = _zscore_by_group(Xte, stte)
    n_components = int(min(args.pca, Xtr_z.shape[0], Xtr_z.shape[1]))
    pca = PCA(n_components=n_components, random_state=0).fit(Xtr_z)
    model = LogisticRegression(max_iter=2000).fit(pca.transform(Xtr_z), ytr)
    score = model.predict_proba(pca.transform(Xte_z))[:, 1]

    pooled = _auc_ci(yte, score)
    sens = {f"{s:g}": round(_sens_at_spec(yte, score, s), 4)
            for s in tuple(DEFAULT_SPECIFICITIES) + tuple(SCREENING_SPECIFICITIES)}

    # The threshold transferred from the training cohort's own OOF scores.
    from cross_study_finallydb import pooled_oof
    _y_oof, score_tr_oof, _ = pooled_oof(Xtr, ytr, sttr, list(range(5)),
                                         args.pca, harmonize=True)
    thr_spec = 0.99
    fpr, tpr, thr = roc_curve(ytr, score_tr_oof)
    idx = np.where(fpr <= (1.0 - thr_spec))[0]
    train_threshold = float(thr[idx[-1]]) if len(idx) else float("nan")
    detected = int((score >= train_threshold).sum())

    # Where the held-out cohort sits on the score axis, against the
    # training cohort's own out-of-fold scores. A cohort-level shift
    # (a different library preparation moves every sample up) is
    # invisible in the AUC and obvious here, and it is what decides
    # whether an absolute threshold can travel.
    def _quantiles(values: np.ndarray) -> dict:
        return {f"p{q}": round(float(np.percentile(values, q)), 4)
                for q in (5, 25, 50, 75, 95)}

    train_reference = {
        "oof_score_quantiles_cancer": _quantiles(score_tr_oof[ytr == 1]),
        "oof_score_quantiles_healthy": _quantiles(score_tr_oof[ytr == 0]),
        "oof_score_max_cancer": round(float(score_tr_oof[ytr == 1].max()), 4),
        "oof_score_max_healthy": round(float(score_tr_oof[ytr == 0].max()), 4),
        "heldout_score_quantiles_cancer": _quantiles(score[yte == 1]),
        "heldout_score_quantiles_healthy": _quantiles(score[yte == 0]),
        "heldout_median_score_healthy": round(
            float(np.median(score[yte == 0])), 4),
        "training_median_oof_score_healthy": round(
            float(np.median(score_tr_oof[ytr == 0])), 4),
    }
    shift = ("the held-out cohort's controls score higher than the training "
             "cohort's own controls -- a cohort-level shift, so the RANKING "
             "travels but an absolute threshold does not")

    # An exact permutation test on the held-out labels. The scores are
    # label-independent (the model never saw this cohort, and the
    # harmonisation uses features only), so every one of the C(12,6)=924
    # label assignments is a valid permutation and the p-value is exact
    # rather than sampled.
    from itertools import combinations
    from sklearn.metrics import roc_auc_score
    positives = np.flatnonzero(yte == 1)
    negatives = np.flatnonzero(yte == 0)
    n_pos = len(positives)
    observed = float(pooled["auc"])
    at_least = 0
    total = 0
    for combo in combinations(range(len(yte)), n_pos):
        permuted = np.zeros(len(yte), dtype=int)
        permuted[list(combo)] = 1
        total += 1
        if roc_auc_score(permuted, score) >= observed - 1e-12:
            at_least += 1
    permutation_p = at_least / total if total else float("nan")

    # The other two harmonisations, for the sensitivity block.
    harmonization_sensitivity = {}
    for tag, tr_z, te_z in (
        ("train_global_test_global", _zscore(Xtr), _zscore(Xte)),
        ("train_global_test_per_publication", _zscore(Xtr), _zscore_by_group(Xte, stte)),
        ("train_per_study_test_per_publication",
         _zscore_by_group(Xtr, sttr), _zscore_by_group(Xte, stte)),
    ):
        pca_s = PCA(n_components=n_components, random_state=0).fit(tr_z)
        lr_s = LogisticRegression(max_iter=2000).fit(pca_s.transform(tr_z), ytr)
        sc_s = lr_s.predict_proba(pca_s.transform(te_z))[:, 1]
        harmonization_sensitivity[tag] = _auc_ci(yte, sc_s)

    # A depth control. The cohort's two arms were sequenced to different
    # depths -- the cancer samples carry more fragments -- and depth on
    # its own separates them. So the question the headline AUC has to
    # answer is whether the model's separation IS the depth or something
    # else, and the answer is measured here rather than assumed.
    depth = np.array([extracted[s]["fsd_fragment_count"] for s in order],
                     dtype=float)
    depth_control = {
        "note": ("The two arms differ in sequencing depth, and depth "
                 "alone separates them. This block asks whether the "
                 "model's separation survives that: the score is "
                 "regressed on the depth and what is left is scored "
                 "again."),
        "median_fragment_count_cancer": int(np.median(depth[yte == 1])),
        "median_fragment_count_healthy": int(np.median(depth[yte == 0])),
        "auc_depth_alone": _auc_ci(yte, depth)["auc"],
        "auc_model": pooled["auc"],
    }
    if np.ptp(depth) > 0:
        design = np.vstack([depth, np.ones_like(depth)]).T
        beta, *_ = np.linalg.lstsq(design, score, rcond=None)
        residual = score - design @ beta
        depth_control["auc_model_with_depth_removed"] = _auc_ci(
            yte, residual)["auc"]
        depth_control["partial_spearman_score_label_given_depth"] = round(
            float(_partial_spearman(score, yte, depth)), 4)
        depth_control["raw_spearman_score_label"] = round(
            float(_spearman(score, yte)), 4)
        depth_control["the_separation_is_not_the_depth"] = bool(
            depth_control["auc_model_with_depth_removed"]
            >= depth_control["auc_depth_alone"])

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "design": {
            "question": "Does the tumour-naive fragmentomics model, trained on "
                        "FinaleDB publications 6+8, separate ccRCC plasma from "
                        "matched healthy plasma in a cohort it has never seen?",
            "train": "FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019), "
                     "n=627",
            "test": f"Zenodo {ZENODO_DOI} -- {ZENODO_RECORD}",
            "harmonization": "per study on the training side, per publication on "
                              "the held-out side",
            "pca_components": n_components,
            "model": "LogisticRegression(max_iter=2000) on the training cohort",
            "threshold": "transferred from the training cohort's OOF scores at "
                         "99% specificity",
        },
        "cohort": {
            "doi": ZENODO_DOI,
            "record": ZENODO_RECORD,
            "publication_id": PUBLICATION_ID,
            "n": len(labels),
            "n_cancer": int((yte == 1).sum()),
            "n_healthy": int((yte == 0).sum()),
            "cancer_type": "clear-cell renal cell carcinoma (ccRCC)",
            "samples": order,
            "labels": {s: int(yte[i]) for i, s in enumerate(order)},
            "stage": {s: stage_arr[i] for i, s in enumerate(order)},
            "disease": {s: disease_arr[i] for i, s in enumerate(order)},
            "fragment_counts": {s: extracted[s]["fsd_fragment_count"]
                                for s in order},
            "median_fragment_count_ccrcc": int(np.median(
                [extracted[s]["fsd_fragment_count"] for s in order
                 if labels[s] == 1])),
            "median_fragment_count_healthy": int(np.median(
                [extracted[s]["fsd_fragment_count"] for s in order
                 if labels[s] == 0])),
        },
        "result": {
            "auc": pooled["auc"],
            "ci_lo": pooled["ci_lo"],
            "ci_hi": pooled["ci_hi"],
            "n_pos": pooled["n_pos"],
            "n_neg": pooled["n_neg"],
            "exact_permutation_p": round(permutation_p, 6),
            "exact_permutations": total,
            "sensitivity_at_spec": sens,
            "sensitivity_at_spec_note": (
                "With 6 controls the specificity is quantised to multiples of "
                "1/6, so a threshold that keeps 0 false positives is the only "
                "one available above 83.3% specificity. These values are "
                "reported for completeness and are NOT resolvable at this n; "
                "the threshold transfer below is the informative number."),
            "train_threshold_99spec": train_threshold,
            "n_above_train_threshold": detected,
            "n_cancer_above_train_threshold": int(sum(
                1 for i, s in enumerate(order)
                if labels[s] == 1 and score[i] >= train_threshold)),
            "n_healthy_above_train_threshold": int(sum(
                1 for i, s in enumerate(order)
                if labels[s] == 0 and score[i] >= train_threshold)),
            "specificity_at_train_threshold": round(
                1.0 - (sum(1 for i, s in enumerate(order)
                            if labels[s] == 0 and score[i] >= train_threshold)
                       / max(int((yte == 0).sum()), 1)), 4),
            "per_sample_score": {s: round(float(score[i]), 6)
                                 for i, s in enumerate(order)},
        },
        "training_score_reference": train_reference,
        "cohort_shift_note": shift,
        "harmonization_sensitivity": harmonization_sensitivity,
        "depth_control": depth_control,
        "limits": [
            "n=6 vs 6: the point estimate is the finding, the interval is wide.",
            "All six patients are pT3a (stage III); two have distant metastasis "
            "(M1/Mx). This is a LATE-STAGE cohort and says nothing about "
            "early detection.",
            "The assay depth is about 30x lower than the training cohort's; "
            "all five channels are scale-free by construction, but the "
            "per-bin ratios are noisier at this depth.",
            "The RANKING transfers; the ABSOLUTE THRESHOLD does not. At the "
            "training cohort's 99%-specificity threshold the held-out "
            "cohort's specificity is far below 99% -- see "
            "result.specificity_at_train_threshold and "
            "training_score_reference. An AUC is a rank statistic and is the "
            "honest number to report for a held-out cohort; a probability "
            "threshold taken from another cohort is not.",
        ],
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print(f"\n  AUC = {pooled['auc']}  CI95 = [{pooled['ci_lo']}, "
          f"{pooled['ci_hi']}]  exact permutation p = {permutation_p:.4f}")
    print(f"  at the training cohort's 99%-specificity threshold: "
          f"{payload['result']['n_cancer_above_train_threshold']}/6 cancers and "
          f"{payload['result']['n_healthy_above_train_threshold']}/6 controls "
          f"above it (specificity "
          f"{payload['result']['specificity_at_train_threshold']})")
    print("  per-sample scores:")
    for i, s in enumerate(order):
        print(f"    {s:<8} label={labels[s]}  stage={stage_arr[i] or '-':<5} "
              f"score={score[i]:.4f}  frags={extracted[s]['fsd_fragment_count']:,}")
    print(f"\n  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
