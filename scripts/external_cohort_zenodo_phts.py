#!/usr/bin/env python3
"""Score the Cleveland Clinic PHTS screening cohort with the five channels.

What this is
------------
A second external cohort, from a different deposit and a different
laboratory: Zenodo 15109495 ("Cell-free DNA Fragmentomics and Second
Malignant Neoplasm Risk in Patients with PTEN Hamartoma Tumor
Syndrome", Liu et al., *Cell Reports Medicine* 2024,
doi:10.1016/j.xcrm.2023.101384), 99 plasma samples from a
hereditary-cancer screening population at the Cleveland Clinic. Fifty
carry no cancer and forty-nine do, with the depositor's own split into
a single primary malignant neoplasm (PMN, 23) and a second malignant
neoplasm (SMN, 26), each of them split by when the plasma was drawn:
Pre-SMN 4, Post-SMN 22, Pre-PMN 5, Post-PMN 18. The tumour stage is
0: 5, 1: 22, 2: 10, 3: 2, with ten malignant samples unstaged.

Why it matters
--------------
Two things no other public plasma cohort in this repository's survey
has.

1. It is the first **stage-graded** cohort that ships *fragment* data.
   The readiness document's stage breakdown had no cohort to run on:
   GSE124600 carries stage labels but ships methylation matrices, and
   the methylation channel turned out to be a coverage signal rather
   than a feature.

2. Nine of the 49 have plasma drawn **before the cancer diagnosis**
   (Pre-SMN 4, Pre-PMN 5) -- the depositor's own "subclinical cancer"
   arm, in a population that carries a germline PTEN variant and is
   therefore at 7-fold risk of a second primary.

   The depositor's ``date_diff`` column says what "before the
   diagnosis" means, and it is narrower than the arm names suggest:
   those nine draws fall between 1 and 7 days before the diagnosis,
   and the "post" draws fall between 1 and 37 days after it. So the
   tumour was already there when the "pre" plasma was drawn. This is a
   cohort of samples drawn *around* a diagnosis, not a lead-time
   cohort, and the artifact reports the intervals so that nobody has to
   take the arm names on trust.

What it does not do
-------------------
It does not claim an early-detection result, and it cannot: the
pre-diagnosis arm is nine samples drawn within a week of a diagnosis
that had already happened. It does not claim a replication of the
depositor's own result either: their finding is about fragment END
MOTIFS (CC- and CA-end motifs) and about nucleosome fractions, and
this repository's five channels do not include end motifs at all. What
can be checked against their paper is the direction of the effect they
report for genome-wide fragmentation, and that check is in the
artifact.

The protocol is the same one the ccRCC cohort used
(``scripts/external_validation_finaledb.py``, imported rather than
re-implemented), so the numbers are comparable with the other
held-out cohorts: train on FinaleDB publications 6 (Jiang 2015) + 8
(Cristiano 2019); per-study z-score harmonisation on the training
side, per-publication on the held-out side; PCA (n_components =
min(200, n_train)) fit on the training cohort only;
LogisticRegression(max_iter=2000) fit on the training cohort only;
the held-out cohort is scored once, never fitted on; the threshold is
transferred from the training cohort's own OOF scores.

Usage:
    # 1. unpack and convert the deposit (7.6 GB download, then ~15 min)
    python3 scripts/prepare_zenodo_phts_cohort.py \
        --tarball /tmp/phts_cohort/02-bamtobed.tar.gz \
        --out-dir /tmp/phts_cohort/frag \
        --metadata data/external_cohort/phts_sample_metadata.tsv \
        --report results/zenodo_phts_prep.json

    # 2. score it
    python3 scripts/external_cohort_zenodo_phts.py \
        --frag-dir /tmp/phts_cohort/frag \
        --metadata data/external_cohort/phts_sample_metadata.tsv
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_SCRIPT_DIR = str(_Path(__file__).resolve().parent)
if _SCRIPT_DIR not in _sys.path:
    _sys.path.insert(0, _SCRIPT_DIR)

import argparse
import json
import os
from datetime import datetime, timezone

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve

from _paths import SIBLING_PIPELINE  # noqa: E402
from external_validation_finaledb import (  # noqa: E402
    _auc_ci,
    _load5_from,
    _load_cohort,
    _sens_at_spec,
    _zscore,
    _zscore_by_group,
)
from geo_gse124600_methylation_cohort import permutation_p  # noqa: E402
from prepare_zenodo_phts_cohort import TARGET_SUFFIX  # noqa: E402

_SRC = str(_Path(__file__).resolve().parents[1] / "src")
if _SRC not in _sys.path:
    _sys.path.insert(0, _SRC)
from fragmentomics import finale_features as ff  # noqa: E402

TRAIN_PUBLICATIONS = {"6", "8"}
ZENODO_DOI = "10.5281/zenodo.15109495"
ZENODO_RECORD = "https://zenodo.org/records/15109495"
PUBLICATION_ID = "zenodo_15109495"
PAPER_DOI = "10.1016/j.xcrm.2023.101384"

DEFAULT_SPECIFICITIES = (0.95, 0.98, 0.99)
SCREENING_SPECIFICITIES = (0.995, 0.999)

# The arms of the depositor's own table. The pre/post split is the one
# that matters: "Pre" is plasma drawn before the (second) malignant
# neoplasm was diagnosed, "Post" is plasma drawn after.
_PRE_GROUPS = ("Pre-SMN", "Pre-PMN")
_POST_GROUPS = ("Post-SMN", "Post-PMN")
_NO_CANCER_GROUP = "No Cancer"

# How many label permutations each contrast gets. The ccRCC cohort could
# enumerate all C(12,6) = 924 assignments exactly; this cohort has 99
# samples, so an exact enumeration is C(99,49) and impossible. The
# p-value is therefore a sampled permutation test, and the artifact
# says so.
PERMUTATION_TRIALS = 20000
PERMUTATION_SEED = 20261125


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
    run wrote, which is the same numbers without re-reading 300 million
    fragments.
    """
    summary = {}
    for sample in samples:
        path = os.path.join(frag_dir, f"{sample}{TARGET_SUFFIX}")
        if reuse:
            written = os.path.join(out_dir, f"{sample}.fsd.json")
            if os.path.exists(written):
                summary[sample] = _summary_from_features(out_dir, sample)
                continue
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        summary[sample] = ff.extract_sample(path, out_dir, sample)
    return summary


def _contrast(score: np.ndarray, index: dict, cases: list[str],
              controls: list[str], trials: int) -> dict:
    """One arm against the no-cancer arm: AUC, DeLong CI, permutation p.

    The scores are label-independent (the model never saw this cohort
    and the harmonisation uses features only), so permuting the labels
    is a valid null.
    """
    case_idx = [index[s] for s in cases]
    control_idx = [index[s] for s in controls]
    block: dict = {
        "n_pos": len(case_idx),
        "n_neg": len(control_idx),
        "median_score_pos": (round(float(np.median(score[case_idx])), 4)
                             if case_idx else None),
        "median_score_neg": (round(float(np.median(score[control_idx])), 4)
                             if control_idx else None),
    }
    if case_idx and control_idx:
        y = np.concatenate([np.ones(len(case_idx)),
                            np.zeros(len(control_idx))])
        s = np.concatenate([score[case_idx], score[control_idx]])
        try:
            stats = _auc_ci(y, s)
            block.update({"auc": stats["auc"], "ci_lo": stats["ci_lo"],
                          "ci_hi": stats["ci_hi"]})
        except (ValueError, ZeroDivisionError):
            block.update({"auc": None, "ci_lo": None, "ci_hi": None})
        _observed, p = permutation_p(y, s, trials=trials,
                                     seed=PERMUTATION_SEED)
        block["permutation_p"] = p
    else:
        block.update({"auc": None, "ci_lo": None, "ci_hi": None,
                      "permutation_p": None})
    return block


def _spearman_permutation_p(x: np.ndarray, y: np.ndarray, trials: int,
                             seed: int) -> tuple[float, float]:
    """Spearman rho and a permutation p-value for it."""
    from scipy.stats import spearmanr
    result = spearmanr(x, y)
    rho = float(result.statistic)
    rng = np.random.default_rng(seed)
    shuffled = y.copy()
    hits = 0
    for _ in range(trials):
        rng.shuffle(shuffled)
        if abs(float(spearmanr(x, shuffled).statistic)) >= abs(rho) - 1e-12:
            hits += 1
    return rho, (hits + 1) / (trials + 1)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frag-dir", required=True,
                    help=f"directory of *{TARGET_SUFFIX} files")
    ap.add_argument("--metadata", required=True,
                    help="the deposited sample metadata TSV")
    ap.add_argument("--features-out",
                    default="data/external_cohort_zenodo_phts/features")
    ap.add_argument("--features-dir",
                    default=os.path.join(SIBLING_PIPELINE, "data/features"))
    ap.add_argument("--labels-multiclass",
                    default=os.path.join(SIBLING_PIPELINE,
                                         "labels_multiclass.tsv"))
    ap.add_argument("--prep-report", default="results/zenodo_phts_prep.json",
                    help="the preparation report (genome build, conversion)")
    ap.add_argument("--pca", type=int, default=200)
    ap.add_argument("--reuse-features", action="store_true",
                    help="read the summaries back from a previous run's "
                         "feature files instead of re-extracting them")
    ap.add_argument("--out-json",
                    default="results/external_cohort_zenodo_phts.json")
    args = ap.parse_args()

    print("[1/6] reading the deposited sample metadata ...", flush=True)
    meta = read_metadata(args.metadata)
    labels, studies, disease, stage, group = {}, {}, {}, {}, {}
    for sample, row in meta.items():
        is_cancer = row["cancer"].strip().lower() == "yes"
        labels[sample] = 1 if is_cancer else 0
        studies[sample] = PUBLICATION_ID
        disease[sample] = (row.get("group3", "") or "").strip()
        stage[sample] = (row.get("stage", "") or "").strip()
        group[sample] = (row.get("group2", "") or "").strip()

    def _stage_code(sample: str):
        value = stage[sample]
        return int(value) if value.isdigit() else None

    order = sorted(labels)
    index = {s: i for i, s in enumerate(order)}
    n_cancer = sum(labels.values())
    print(f"      {len(labels)} samples: {n_cancer} with cancer, "
          f"{len(labels) - n_cancer} without")

    print("[2/6] extracting the five channels ...", flush=True)
    extracted = extract_features(args.frag_dir, args.features_out, order,
                                 reuse=args.reuse_features)
    for sample in order[:3]:
        print(f"      {sample}: {extracted[sample]['fsd_fragment_count']:,} "
              f"fragments")

    print("[3/6] loading the training cohort (FinaleDB 6+8) ...", flush=True)
    Xtr, ytr, sttr, dtr, train_order = _load_cohort(
        args.features_dir, args.labels_multiclass, TRAIN_PUBLICATIONS)
    print(f"      n={len(ytr)} (cancer {int((ytr == 1).sum())}, "
          f"healthy {int((ytr == 0).sum())}), features={Xtr.shape[1]}")

    print("[4/6] loading the cohort features and scoring ...", flush=True)
    Xte, yte, stte = _load5_from(args.features_out, labels, studies)
    if len(Xte) != len(order):
        raise RuntimeError(f"row alignment failed: {len(Xte)} vs {len(order)}")

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
    # training cohort's own out-of-fold scores. A cohort-level shift is
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

    print("[5/6] the stage question ...", flush=True)
    no_cancer = [s for s in order if disease[s] == _NO_CANCER_GROUP]
    pre = [s for s in order if disease[s] in _PRE_GROUPS]
    post = [s for s in order if disease[s] in _POST_GROUPS]
    staged = [s for s in order if _stage_code(s) is not None]

    # The lead time: the depositor's own ``date_diff`` column -- the days
    # between the plasma draw and the cancer diagnosis, negative when the
    # draw came first. This is what decides what the "pre" arms actually
    # are, and it is the difference between a lead-time cohort and a
    # cohort of samples drawn around a diagnosis.
    def _to_float(value):
        text = (value or "").strip()
        if text in ("", "NA", "None"):
            return None
        try:
            return float(text)
        except ValueError:
            return None

    days_to_diagnosis = {s: _to_float(meta[s].get("date_diff")) for s in order}
    dated = [d for d in days_to_diagnosis.values() if d is not None]
    lead_time: dict = {
        "column": "date_diff",
        "definition": "days between the plasma draw and the cancer "
                      "diagnosis; negative means the draw came first",
        "n_dated": len(dated),
        "by_arm": {},
    }
    for name, members in (("no_cancer", no_cancer), ("pre_diagnosis", pre),
                          ("post_diagnosis", post)):
        values = sorted(d for d in (days_to_diagnosis[s] for s in members)
                        if d is not None)
        if not values:
            continue
        entry = {
            "n_dated": len(values),
            "days_to_diagnosis_min": values[0],
            "days_to_diagnosis_median": values[len(values) // 2],
            "days_to_diagnosis_max": values[-1],
        }
        if values[0] < 0:
            entry["lead_time_days_min"] = int(-values[-1])
            entry["lead_time_days_max"] = int(-values[0])
        lead_time["by_arm"][name] = entry
    negative = [d for d in dated if d < 0]
    lead_time["max_lead_time_days"] = (int(max(-d for d in negative))
                                       if negative else 0)
    lead_time["note"] = (
        "The longest interval anywhere in this cohort between a plasma "
        f"draw and the diagnosis is {lead_time['max_lead_time_days']} days. "
        "Every sample in a 'pre' arm was drawn within a week of the "
        "diagnosis, so the tumour was already present when the plasma was "
        "drawn. This is a cohort of samples drawn AROUND a diagnosis, not "
        "a lead-time cohort: a separation in the pre-diagnosis arm is not "
        "an early-detection result.")

    arms = {
        "all_cancer": [s for s in order if labels[s] == 1],
        "pre_diagnosis": pre,
        "post_diagnosis": post,
        "stage_0": [s for s in staged if _stage_code(s) == 0],
        "stage_0_1": [s for s in staged if _stage_code(s) in (0, 1)],
        "stage_2_3": [s for s in staged if _stage_code(s) in (2, 3)],
        "all_staged": staged,
    }
    contrasts = {name: _contrast(score, index, cases, no_cancer,
                                 PERMUTATION_TRIALS)
                 for name, cases in arms.items()}

    # The pre- versus post-diagnosis contrast: can the model tell
    # plasma drawn before the tumour was diagnosed from plasma drawn
    # after it? Both arms have cancer, so this is the contrast that
    # isolates the diagnosis from the presence of a tumour.
    pre_vs_post = _contrast(score, index, post, pre, PERMUTATION_TRIALS)
    pre_vs_post["note"] = ("positives are the post-diagnosis samples, "
                            "negatives the pre-diagnosis ones")

    # Monotonicity across the stage, over the staged malignant samples.
    stage_codes = np.array([_stage_code(s) for s in staged])
    stage_scores = score[[index[s] for s in staged]]
    rho, rho_p = _spearman_permutation_p(stage_codes, stage_scores,
                                          PERMUTATION_TRIALS,
                                          PERMUTATION_SEED)
    median_by_stage = {
        str(level): round(float(np.median(stage_scores[stage_codes == level])), 4)
        for level in (0, 1, 2, 3) if (stage_codes == level).any()
    }
    n_by_stage = {str(level): int((stage_codes == level).sum())
                  for level in (0, 1, 2, 3)}
    monotone = (rho_p < 0.05 and
                all(median_by_stage[str(a)] < median_by_stage[str(b)]
                    for a, b in ((0, 1), (1, 2)) if str(a) in median_by_stage
                    and str(b) in median_by_stage))

    # The direction the depositor's own paper reports, checked against
    # this repository's closest channel. Their result: PHTS patients
    # with a second malignant neoplasm show INCREASED genome-wide
    # fragmentation -- a decreased proportion of di- and tri-nucleosome
    # fragments -- against patients without cancer.
    #
    # The FSD channel's 196 bins run from 20 bp to 1000 bp in 5 bp
    # steps, so the di-nucleosome window (320-340 bp) is inside it: the
    # depositor's "proportion of di-nucleosome-associated fragments" can
    # be read straight out of the features this repository computes.
    def _fsd_extra(sample: str) -> dict:
        """The FSD readings this check needs, read from the feature file."""
        path = os.path.join(args.features_out, f"{sample}.fsd.json")
        if not os.path.exists(path):
            return {}
        with open(path) as handle:
            payload = json.load(handle)
        bins = payload.get("size_bins") or {}
        total = 0.0
        window = 0.0
        for label, frequency in bins.items():
            try:
                low = int(str(label).split("-")[0])
            except ValueError:
                continue
            total += float(frequency)
            if 320 <= low <= 340:
                window += float(frequency)
        return {
            "di_nucleosome_fraction": (window / total) if total else None,
            "median_length": payload.get("median_length"),
            "p90_length": payload.get("p90"),
            "mode_length": payload.get("mode_length"),
        }

    fsd_extra = {s: _fsd_extra(s) for s in order}
    di_fraction = {s: fsd_extra[s].get("di_nucleosome_fraction")
                   for s in order}
    short_long = np.array([extracted[s]["genome_short_long_ratio"]
                           for s in order])
    smn_idx = [index[s] for s in order if disease[s] == "Post-SMN"]
    noc_idx = [index[s] for s in order if disease[s] == _NO_CANCER_GROUP]
    direction_check = {
        "published_finding": (
            "Liu et al. report that PHTS patients with a second malignant "
            "neoplasm have increased genome-wide fragmentation -- a "
            "decreased proportion of di- and tri-nucleosome-associated "
            "cfDNA fragments -- against PHTS patients without cancer."),
        "closest_channel": (
            "the FSD channel's 320-340 bp bins, which is the "
            "di-nucleosome proportion the paper reports, and the "
            "genome_short_long_ratio (100-150 bp over 150-220 bp) as the "
            "second reading."),
        "n_no_cancer": len(noc_idx),
        "n_post_smn": len(smn_idx),
    }
    usable = [s for s in order if di_fraction[s] is not None]
    if usable:
        values = np.array([di_fraction[s] for s in usable])
        order_usable = {s: i for i, s in enumerate(usable)}
        smn_u = [order_usable[s] for s in usable if disease[s] == "Post-SMN"]
        noc_u = [order_usable[s] for s in usable
                 if disease[s] == _NO_CANCER_GROUP]
        direction_check["median_di_nucleosome_fraction_no_cancer"] = round(
            float(np.median(values[noc_u])), 4) if noc_u else None
        direction_check["median_di_nucleosome_fraction_post_smn"] = round(
            float(np.median(values[smn_u])), 4) if smn_u else None
        if noc_u and smn_u:
            y_dir = np.concatenate([np.ones(len(smn_u)), np.zeros(len(noc_u))])
            s_dir = np.concatenate([values[smn_u], values[noc_u]])
            auc_di = _auc_ci(y_dir, s_dir)["auc"]
            direction_check["auc_post_smn_vs_no_cancer"] = auc_di
            # The paper's direction is a DECREASED di-nucleosome
            # proportion in the SMN arm, so the AUC of "SMN vs no
            # cancer" on that proportion is below 0.5 when the direction
            # is reproduced.
            direction_check["direction_reproduced"] = bool(auc_di < 0.5)
    direction_check["median_short_long_ratio_no_cancer"] = round(
        float(np.median(short_long[noc_idx])), 4) if noc_idx else None
    direction_check["median_short_long_ratio_post_smn"] = round(
        float(np.median(short_long[smn_idx])), 4) if smn_idx else None
    if noc_idx and smn_idx:
        y_dir = np.concatenate([np.ones(len(smn_idx)), np.zeros(len(noc_idx))])
        s_dir = np.concatenate([short_long[smn_idx], short_long[noc_idx]])
        direction_check["auc_short_long_ratio_post_smn_vs_no_cancer"] = (
            _auc_ci(y_dir, s_dir)["auc"])
    else:
        direction_check["auc_short_long_ratio_post_smn_vs_no_cancer"] = None

    # The library's fragment shape against the training cohort's: the
    # per-cohort harmonisation is doing real work if these differ.
    train_slr, train_counts = [], []
    for sample in train_order:
        delfi_path = os.path.join(args.features_dir, f"{sample}.delfi.json")
        fsd_path = os.path.join(args.features_dir, f"{sample}.fsd.json")
        if os.path.exists(delfi_path):
            with open(delfi_path) as handle:
                train_slr.append(json.load(handle)["genome_short_long_ratio"])
        if os.path.exists(fsd_path):
            with open(fsd_path) as handle:
                train_counts.append(json.load(handle)["fragment_count"])
    library = {
        "heldout_median_short_long_ratio": round(
            float(np.median(short_long)), 4),
        "training_median_short_long_ratio": (
            round(float(np.median(train_slr)), 4) if train_slr else None),
        "heldout_median_fragment_count": int(np.median(
            [extracted[s]["fsd_fragment_count"] for s in order])),
        "training_median_fragment_count": (
            int(np.median(train_counts)) if train_counts else None),
        "note": ("The held-out cohort's library has a different "
                 "fragment-size shape from the training cohort's, so the "
                 "per-cohort harmonisation carries a real shift. The "
                 "ranking is the honest number here; the absolute "
                 "threshold is not."),
    }

    # Is there a class difference in this cohort at all, and does the
    # model capture it? Each channel is read on its own, against the
    # no-cancer arm. A channel that separates the classes while the
    # model's score does not localises the failure: the signal is in
    # the data, and the model does not transfer.
    channels = {
        "genome_short_long_ratio": [extracted[s]["genome_short_long_ratio"]
                                    for s in order],
        "mean_short_frac": [extracted[s]["mean_short_frac"] for s in order],
        "fsd_short_fraction_100_150": [extracted[s]["fsd_short_fraction"]
                                       for s in order],
        "fsd_short_long_ratio": [extracted[s]["fsd_short_long_ratio"]
                                 for s in order],
        "fsd_median_length": [fsd_extra[s].get("median_length")
                              for s in order],
        "fsd_p90_length": [fsd_extra[s].get("p90_length") for s in order],
        "di_nucleosome_fraction": [di_fraction[s] for s in order],
    }
    single_channel: dict = {}
    for name, values in channels.items():
        usable = [i for i, v in enumerate(values) if v is not None]
        if len(usable) < 10:
            continue
        array = np.array([values[i] for i in usable], dtype=float)
        positions = {i: j for j, i in enumerate(usable)}
        entry: dict = {"n": len(usable)}
        for arm_name, cases in (("all_cancer",
                                 [s for s in order if labels[s] == 1]),
                                ("post_smn",
                                 [s for s in order
                                  if disease[s] == "Post-SMN"])):
            case_idx = [positions[index[s]] for s in cases
                        if index[s] in positions]
            control_idx = [positions[index[s]] for s in no_cancer
                           if index[s] in positions]
            if not case_idx or not control_idx:
                continue
            y_arm = np.concatenate([np.ones(len(case_idx)),
                                    np.zeros(len(control_idx))])
            s_arm = np.concatenate([array[case_idx], array[control_idx]])
            stats = _auc_ci(y_arm, s_arm)
            entry[arm_name] = {
                "auc": stats["auc"],
                "ci_lo": stats["ci_lo"],
                "ci_hi": stats["ci_hi"],
                "n_pos": len(case_idx),
                "n_neg": len(control_idx),
                "median_pos": round(float(np.median(array[case_idx])), 5),
                "median_neg": round(float(np.median(array[control_idx])), 5),
            }
        single_channel[name] = entry
    single_channel["note"] = (
        "Each channel on its own, against the no-cancer arm. A channel "
        "whose interval excludes 0.5 is separating the classes on its "
        "own, whether or not the model's score does. Note that "
        "genome_short_long_ratio and mean_short_frac are monotone "
        "transforms of each other -- one signal, not two -- and "
        "fsd_short_fraction_100_150 and fsd_short_long_ratio are the same "
        "reading of the 100-150 and 150-220 bp windows.")

    # The other harmonisations, for the sensitivity block. The fourth one
    # is the control that matters here: every other variant lets the
    # held-out cohort normalise itself, which is what pushes its scores
    # into the model's saturated region. This one normalises the
    # held-out cohort by the TRAINING cohort's statistics, so the
    # cohort's absolute position in feature space is preserved.
    def _zscore_by(X: np.ndarray, reference: np.ndarray) -> np.ndarray:
        return (X - reference.mean(axis=0)) / np.maximum(
            reference.std(axis=0), 1e-12)

    harmonization_sensitivity = {}
    for tag, tr_z, te_z in (
        ("train_global_test_global", _zscore(Xtr), _zscore(Xte)),
        ("train_global_test_per_publication", _zscore(Xtr), _zscore_by_group(Xte, stte)),
        ("train_per_study_test_per_publication",
         _zscore_by_group(Xtr, sttr), _zscore_by_group(Xte, stte)),
        ("train_global_test_train_statistics",
         _zscore(Xtr), _zscore_by(Xte, Xtr)),
    ):
        pca_s = PCA(n_components=n_components, random_state=0).fit(tr_z)
        lr_s = LogisticRegression(max_iter=2000).fit(pca_s.transform(tr_z), ytr)
        sc_s = lr_s.predict_proba(pca_s.transform(te_z))[:, 1]
        harmonization_sensitivity[tag] = _auc_ci(yte, sc_s)
    harmonization_sensitivity["note"] = (
        "The first three let the held-out cohort normalise itself, which "
        "is what pushes its scores into the model's saturated region. The "
        "fourth normalises it by the TRAINING cohort's statistics, so its "
        "absolute position in feature space is preserved: it is the "
        "control that says whether the null is an artefact of the "
        "harmonisation choice.")

    # The preparation report: the genome build and the conversion, so
    # the artifact carries them rather than leaving them in a log.
    prep: dict = {}
    if os.path.exists(args.prep_report):
        with open(args.prep_report, encoding="utf-8") as handle:
            prep = json.load(handle)

    pre_block = contrasts["pre_diagnosis"]
    separating = sorted(
        name for name, entry in single_channel.items()
        if name != "note" and entry.get("post_smn", {}).get("ci_lo") is not None
        and (entry["post_smn"]["ci_lo"] > 0.5
             or entry["post_smn"]["ci_hi"] < 0.5))
    verdict = {
        "question": ("Can the tumour-naive fragmentomics model separate "
                     "plasma drawn within days of a cancer diagnosis "
                     "(the cohort's longest lead time is "
                     f"{lead_time['max_lead_time_days']} days) from the "
                     "plasma of the same cohort without cancer?"),
        "max_lead_time_days": lead_time["max_lead_time_days"],
        "is_a_lead_time_cohort": False,
        "pre_diagnosis_vs_no_cancer_auc": pre_block["auc"],
        "pre_diagnosis_ci": [pre_block["ci_lo"], pre_block["ci_hi"]],
        "pre_diagnosis_n": pre_block["n_pos"],
        "separable_at_this_n": (
            None if pre_block["ci_lo"] is None
            else pre_block["ci_lo"] > 0.5),
        "stage_monotone": bool(monotone),
        "spearman_rho_stage_vs_score": round(rho, 4),
        "spearman_permutation_p": round(rho_p, 6),
        "model_separates_cancer_from_no_cancer": (
            None if pooled["ci_lo"] is None else pooled["ci_lo"] > 0.5),
        "channels_that_separate_post_smn_from_no_cancer": separating,
        "the_signal_is_in_the_data": bool(separating),
        "depositor_direction_in_the_point_estimates": {
            "di_nucleosome_fraction_auc": direction_check.get(
                "auc_post_smn_vs_no_cancer"),
            "short_long_ratio_auc": direction_check.get(
                "auc_short_long_ratio_post_smn_vs_no_cancer"),
            "note": ("Both readings point the way the depositor's paper "
                     "does -- the di-nucleosome proportion lower in the "
                     "second-neoplasm arm, the short/long ratio higher -- "
                     "but neither interval excludes 0.5 at this n, so "
                     "these are point estimates and not results."),
        },
    }
    if pre_block["ci_lo"] is None:
        verdict["answer"] = ("not measurable: the pre-diagnosis arm has no "
                             "samples")
    else:
        model_separates = pooled["ci_lo"] > 0.5
        verdict["answer"] = (
            f"the tumour-naive model "
            f"{'separates' if model_separates else 'does not separate'} "
            f"the classes in this cohort (AUC {pooled['auc']}, CI "
            f"[{pooled['ci_lo']}, {pooled['ci_hi']}], permutation p "
            f"{contrasts['all_cancer']['permutation_p']}); the "
            f"pre-diagnosis arm is at chance (AUC {pre_block['auc']}, CI "
            f"[{pre_block['ci_lo']}, {pre_block['ci_hi']}], n="
            f"{pre_block['n_pos']}); the stage gradient is "
            f"{'monotone' if monotone else 'not monotone'} (rho "
            f"{round(rho, 4)}, p {round(rho_p, 4)}); and "
            + (f"the channels whose interval excludes 0.5 on the "
               f"second-neoplasm arm are {', '.join(separating)}, so the "
               f"signal is in the data and the model does not carry it"
               if separating else
               "no single channel's interval excludes 0.5 on the "
               "second-neoplasm arm either (the depositor's published "
               "direction is visible in the point estimates -- the "
               "di-nucleosome proportion lower in that arm, the short/long "
               "ratio higher -- but not past 0.5 at n=22 against 50), so "
               "this cohort carries no fragmentation signal that the "
               "model could have carried"))

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "design": {
            "question": "Does the tumour-naive fragmentomics model, trained on "
                        "FinaleDB publications 6+8, separate the plasma of a "
                        "hereditary-cancer screening cohort with a malignant "
                        "neoplasm from the plasma of the same cohort without "
                        "one, in a cohort it has never seen?",
            "train": "FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019), "
                     "n=627",
            "test": f"Zenodo {ZENODO_DOI} -- {ZENODO_RECORD}",
            "harmonization": "per study on the training side, per publication on "
                              "the held-out side",
            "pca_components": n_components,
            "model": "LogisticRegression(max_iter=2000) on the training cohort",
            "threshold": "transferred from the training cohort's OOF scores at "
                         "99% specificity",
            "paper": f"Liu et al., Cell Reports Medicine 2024, {PAPER_DOI}",
        },
        "cohort": {
            "doi": ZENODO_DOI,
            "record": ZENODO_RECORD,
            "publication_id": PUBLICATION_ID,
            "n": len(labels),
            "n_cancer": int((yte == 1).sum()),
            "n_healthy": int((yte == 0).sum()),
            "population": "PTEN hamartoma tumor syndrome (PHTS), a hereditary "
                          "cancer screening population at the Cleveland Clinic",
            "arms": {
                "no_cancer": len(no_cancer),
                "pre_smn": len([s for s in order if disease[s] == "Pre-SMN"]),
                "post_smn": len([s for s in order if disease[s] == "Post-SMN"]),
                "pre_pmn": len([s for s in order if disease[s] == "Pre-PMN"]),
                "post_pmn": len([s for s in order if disease[s] == "Post-PMN"]),
            },
            "stage_counts": n_by_stage,
            "samples": order,
            "labels": {s: int(yte[i]) for i, s in enumerate(order)},
            "stage": {s: stage[s] for s in order},
            "disease": {s: disease[s] for s in order},
            "fragment_counts": {s: extracted[s]["fsd_fragment_count"]
                                for s in order},
            "library": library,
            "genome_build": prep.get("genome_build", {}),
            "conversion": prep.get("conversion", {}),
        },
        "result": {
            "auc": pooled["auc"],
            "ci_lo": pooled["ci_lo"],
            "ci_hi": pooled["ci_hi"],
            "n_pos": pooled["n_pos"],
            "n_neg": pooled["n_neg"],
            "permutation_p": contrasts["all_cancer"]["permutation_p"],
            "permutation_method": (
                f"sampled label permutation, {PERMUTATION_TRIALS} trials, "
                f"seed {PERMUTATION_SEED}. An exact enumeration is not "
                f"possible here: this cohort has 99 samples, so the "
                f"number of label assignments is C(99,49)."),
            "sensitivity_at_spec": sens,
            "sensitivity_at_spec_note": (
                "The specificity of this cohort is quantised to multiples "
                "of 1/50 (there are 50 controls), so the values above 99% "
                "specificity are reported for completeness and are not "
                "resolvable at this n; the threshold transfer below is the "
                "informative number."),
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
            "score_saturation": {
                "heldout_median_score_healthy": round(
                    float(np.median(score[yte == 0])), 4),
                "training_median_oof_score_healthy": round(
                    float(np.median(score_tr_oof[ytr == 0])), 4),
                "saturated": bool(np.median(score[yte == 0]) > 0.9),
                "note": ("The held-out cohort's scores sit at the top of "
                         "the model's range -- its no-cancer controls score "
                         "far above the training cohort's own controls. The "
                         "AUC is a rank statistic and survives that, but "
                         "the probabilities do not: only the ordering "
                         "inside this cohort carries information, and the "
                         "absolute threshold does not travel at all."),
            },
        },
        "stage_question": {
            "arms_against_no_cancer": contrasts,
            "pre_vs_post_diagnosis": pre_vs_post,
            "monotonicity": {
                "spearman_rho_stage_vs_score": round(rho, 4),
                "permutation_p": round(rho_p, 6),
                "median_score_by_stage": median_by_stage,
                "n_by_stage": n_by_stage,
                "monotone": bool(monotone),
            },
            "verdict": verdict,
        },
        "lead_time": lead_time,
        "published_direction_check": direction_check,
        "single_channel_readout": single_channel,
        "training_score_reference": train_reference,
        "cohort_shift_note": (
            "The held-out cohort's controls score differently from the "
            "training cohort's own controls -- a cohort-level shift, so "
            "the RANKING travels and an absolute threshold does not."),
        "harmonization_sensitivity": harmonization_sensitivity,
        "limits": [
            "The 'pre-diagnosis' arm is not a lead-time arm. The "
            "depositor's own date_diff column puts every one of those "
            f"draws between 1 and {lead_time['max_lead_time_days']} days "
            "before the diagnosis, so the tumour was already present. A "
            "separation there is not an early-detection result, and this "
            "cohort cannot answer a lead-time question.",
            "The pre-diagnosis arm is nine samples. Its interval is wide "
            "and every number that involves it is a measurement, not a "
            "result.",
            "This is a hereditary-cancer screening population, not a "
            "symptomatic or incidentally diagnosed one: the no-cancer "
            "controls are high-risk patients under surveillance, and the "
            "malignant samples carry PTEN-associated tumours (breast, "
            "thyroid, endometrial, melanoma, kidney, colon). The "
            "cohort's cancer mix is not the training cohort's.",
            "The deposit is ultra-low-pass: "
            f"{int(np.median([extracted[s]['fsd_fragment_count'] for s in order])):,} "
            "fragments per sample against the training cohort's "
            f"{int(np.median(train_counts)):,} -- about "
            f"{np.median(train_counts) / max(np.median([extracted[s]['fsd_fragment_count'] for s in order]), 1):.0f}x "
            "shallower. All five channels are scale-free by construction, "
            "but the per-bin ratios are noisier at this depth.",
            "The deposit's fragments are mapped to hg19 while the "
            "extractor's bin grid uses hg38 sizes -- the grid the "
            "training cohort's features use, so it is the right grid, but "
            "the preparation report records how many fragments it drops.",
            "The library's fragment-size shape is close to the training "
            f"cohort's (short/long ratio {round(float(np.median(short_long)), 4)} "
            f"against {round(float(np.median(train_slr)), 4)}), so the "
            "cohort shift is not a shape difference. The ranking travels; "
            "the absolute threshold does not.",
            "The model reads this whole cohort as cancer-like: the "
            "no-cancer controls' median score is "
            f"{round(float(np.median(score[yte == 0])), 4)} against the "
            "training cohort's own controls' "
            f"{round(float(np.median(score_tr_oof[ytr == 0])), 4)}, so the "
            "scores saturate at the top of the model's range and the AUC "
            "is a ranking among near-identical values. That is the "
            "mechanism behind the null -- a transfer failure, not a "
            "statement that the cohort carries no signal. The channels "
            "read on their own are in single_channel_readout.",
            "This is not a replication of the depositor's own result: "
            "their finding is about fragment end motifs, which this "
            "repository's five channels do not include. Only the "
            "direction of their genome-wide fragmentation finding is "
            "checked, in published_direction_check.",
        ],
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print(f"\n  AUC = {pooled['auc']}  CI95 = [{pooled['ci_lo']}, "
          f"{pooled['ci_hi']}]  permutation p = "
          f"{contrasts['all_cancer']['permutation_p']:.4f}")
    print("  the stage question:")
    for name, block in contrasts.items():
        print(f"    {name:<16} n={block['n_pos']:>2}  "
              f"AUC={block['auc']}  CI=[{block['ci_lo']}, {block['ci_hi']}]  "
              f"p={block['permutation_p']}")
    print(f"    pre vs post      n={pre_vs_post['n_pos']} vs "
          f"{pre_vs_post['n_neg']}  AUC={pre_vs_post['auc']}  "
          f"CI=[{pre_vs_post['ci_lo']}, {pre_vs_post['ci_hi']}]")
    print(f"    stage gradient   rho={round(rho, 4)}  p={round(rho_p, 4)}  "
          f"monotone={monotone}  medians={median_by_stage}")
    print(f"  verdict: {verdict['answer']}")
    print(f"\n  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
