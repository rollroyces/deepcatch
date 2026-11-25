#!/usr/bin/env python3
"""Score the SAUNA cohort: a tumour-fraction dose-response readout.

What this is
------------
The SAUNA deposit (Zenodo 15097509, "SAUNA -- Simulated Annealing for
Unique Nucleosome Arrangements"; Balazs, Zsolt) is nine cfDNA fragment
files from colorectal-cancer patients and benign controls, each mapped to
**hg38** -- the same build the repository's feature grid uses -- at
14-33x coverage, which is the training cohort's own depth. The files
are already in the fragment form the extractor reads, so nothing has to
be converted.

Why it matters
--------------
It is the only public cfDNA fragment cohort this repository knows of
that carries a **graded tumour fraction** in the same cancer type as
the training cohort: 0.0049, 0.0071, 0.0174, 0.0249, 0.2903, 0.3271
and 0.4029 across the seven CRC samples, plus two benign samples at
0.0062 and 0.0072. One patient (SeCT-26) was sampled three times, with
the tumour fraction rising 0.2903 -> 0.3271 -> 0.4029.

That makes it a **dose-response** cohort rather than a cancer/control
cohort, and dose response is the sharper question: does the score rise
with the tumour burden, over an eighty-fold range of it? The PHTS
cohort's null (the tumour-naive model does not transfer to a screening
population) is ambiguous between "the channel does not work" and "that
cohort was shifted". This cohort is the control for that ambiguity: same
build, same depth, same cancer type, and one patient measured against
himself.

What it does not do
-------------------
It does not claim a cancer-detection result. Nine samples, and the
cancer/control contrast among them is close to meaningless on its own:
the two benign samples carry a tumour fraction of 0.0062 and 0.0072,
and three of the CRC samples carry 0.0049 to 0.0249 -- the same order of
magnitude. What the cohort can carry is the gradient and the
longitudinal series, and both are reported as measurements with their
intervals.
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
from sklearn.metrics import roc_auc_score

from _paths import SIBLING_PIPELINE  # noqa: E402
from external_validation_finaledb import (  # noqa: E402
    _auc_ci,
    _load5_from,
    _load_cohort,
    _zscore,
    _zscore_by_group,
)

_SRC = str(_Path(__file__).resolve().parents[1] / "src")
if _SRC not in _sys.path:
    _sys.path.insert(0, _SRC)
from fragmentomics import finale_features as ff  # noqa: E402

TRAIN_PUBLICATIONS = {"6", "8"}
ZENODO_DOI = "10.5281/zenodo.15097509"
ZENODO_RECORD = "https://zenodo.org/records/15097509"
PUBLICATION_ID = "zenodo_15097509"

PERMUTATION_TRIALS = 20000
PERMUTATION_SEED = 20261125


def read_metadata(path: str) -> dict:
    """The depositor's sample table: group, coverage and tumour fraction."""
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return {row["sample"]: row for row in rows}


def _spearman_permutation_p(x: np.ndarray, y: np.ndarray, trials: int,
                             seed: int) -> tuple[float, float]:
    """Spearman rho and a permutation p-value for it."""
    from scipy.stats import spearmanr
    result = spearmanr(x, y)
    rho = float(result.statistic)
    rng = np.random.default_rng(seed)
    shuffled = np.array(y, dtype=float)
    hits = 0
    for _ in range(trials):
        rng.shuffle(shuffled)
        if abs(float(spearmanr(x, shuffled).statistic)) >= abs(rho) - 1e-12:
            hits += 1
    return rho, (hits + 1) / (trials + 1)


def _zscore_by(X: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return (X - reference.mean(axis=0)) / np.maximum(
        reference.std(axis=0), 1e-12)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frag-dir", required=True,
                    help="directory of the deposit's <sample>_hg38.frag.tsv.bgz")
    ap.add_argument("--metadata", required=True,
                    help="the deposited metadata.tsv")
    ap.add_argument("--features-out",
                    default="data/external_cohort_zenodo_sauna/features")
    ap.add_argument("--features-dir",
                    default=os.path.join(SIBLING_PIPELINE, "data/features"))
    ap.add_argument("--labels-multiclass",
                    default=os.path.join(SIBLING_PIPELINE,
                                         "labels_multiclass.tsv"))
    ap.add_argument("--prep-report", default="results/zenodo_sauna_prep.json")
    ap.add_argument("--pca", type=int, default=200)
    ap.add_argument("--out-json",
                    default="results/external_cohort_zenodo_sauna.json")
    args = ap.parse_args()

    print("[1/6] reading the deposited metadata ...", flush=True)
    meta = read_metadata(args.metadata)
    order = sorted(meta)
    index = {s: i for i, s in enumerate(order)}
    tumour_fraction = np.array([float(meta[s]["tf"]) for s in order])
    is_cancer = np.array([meta[s]["group"].strip().lower() == "crc"
                          for s in order])
    print(f"      {len(order)} samples: {int(is_cancer.sum())} CRC, "
          f"{int((~is_cancer).sum())} benign; tumour fraction "
          f"{tumour_fraction.min():.4f}-{tumour_fraction.max():.4f}")

    print("[2/6] extracting the five channels ...", flush=True)
    extracted = {}
    for sample in order:
        path = os.path.join(args.frag_dir, f"{sample}_hg38.frag.tsv.bgz")
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        extracted[sample] = ff.extract_sample(path, args.features_out, sample)
        print(f"      {sample}: {extracted[sample]['fsd_fragment_count']:,} "
              f"fragments", flush=True)

    print("[3/6] loading the training cohort (FinaleDB 6+8) ...", flush=True)
    Xtr, ytr, sttr, dtr, train_order = _load_cohort(
        args.features_dir, args.labels_multiclass, TRAIN_PUBLICATIONS)
    print(f"      n={len(ytr)} (cancer {int((ytr == 1).sum())}, "
          f"healthy {int((ytr == 0).sum())}), features={Xtr.shape[1]}")

    print("[4/6] scoring the cohort ...", flush=True)
    labels = {s: 1 if is_cancer[i] else 0 for i, s in enumerate(order)}
    studies = {s: PUBLICATION_ID for s in order}
    Xte, yte, stte = _load5_from(args.features_out, labels, studies)
    if len(Xte) != len(order):
        raise RuntimeError(f"row alignment failed: {len(Xte)} vs {len(order)}")
    # The deposit's files are sorted by chromosome and position, so a
    # truncated file leaves most of the 100 kb bins empty -- the loader
    # then divides the counts by a zero median and hands back a matrix
    # full of NaN, which scikit-learn reports as a puzzling "Input X
    # contains NaN". The preparation script checks every file's size
    # first; this is the second line of defence, and it says what the
    # real cause is.
    if not np.isfinite(Xte).all():
        bad = int((~np.isfinite(Xte)).sum())
        raise SystemExit(
            f"the loaded features hold {bad} non-finite values. The most "
            "likely cause is a truncated fragment file: these files are "
            "sorted by position, so a partial file leaves most 100 kb bins "
            "empty and the loader divides by a zero median. Check the "
            "files' sizes against the survey's record and re-run "
            "scripts/prepare_zenodo_sauna_cohort.py before scoring.")

    Xtr_z = _zscore_by_group(Xtr, sttr)
    Xte_z = _zscore_by_group(Xte, stte)
    n_components = int(min(args.pca, Xtr_z.shape[0], Xtr_z.shape[1]))
    pca = PCA(n_components=n_components, random_state=0).fit(Xtr_z)
    model = LogisticRegression(max_iter=2000).fit(pca.transform(Xtr_z), ytr)
    score = model.predict_proba(pca.transform(Xte_z))[:, 1]

    # The dose response: does the score rise with the tumour fraction?
    rho_model, rho_model_p = _spearman_permutation_p(
        tumour_fraction, score, PERMUTATION_TRIALS, PERMUTATION_SEED)
    crc = np.flatnonzero(is_cancer)
    rho_crc, rho_crc_p = _spearman_permutation_p(
        tumour_fraction[crc], score[crc], PERMUTATION_TRIALS, PERMUTATION_SEED)

    # The same question for each channel on its own. If the channels
    # track the tumour fraction and the model's score does not, the
    # model is the thing that fails; if neither does, the channel does
    # not work at this cohort's depth.
    def _channel_values(key: str) -> np.ndarray:
        return np.array([extracted[s][key] for s in order], dtype=float)

    channels = {
        "genome_short_long_ratio": _channel_values("genome_short_long_ratio"),
        "mean_short_frac": _channel_values("mean_short_frac"),
        "fsd_short_fraction_100_150": _channel_values("fsd_short_fraction"),
        "fsd_short_long_ratio": _channel_values("fsd_short_long_ratio"),
    }
    single_channel = {}
    for name, values in channels.items():
        rho, p = _spearman_permutation_p(tumour_fraction, values,
                                          PERMUTATION_TRIALS,
                                          PERMUTATION_SEED)
        rho_c, p_c = _spearman_permutation_p(tumour_fraction[crc],
                                              values[crc],
                                              PERMUTATION_TRIALS,
                                              PERMUTATION_SEED)
        single_channel[name] = {
            "rho_vs_tumour_fraction": round(rho, 4),
            "permutation_p": round(p, 6),
            "rho_vs_tumour_fraction_crc_only": round(rho_c, 4),
            "permutation_p_crc_only": round(p_c, 6),
            "median_crc": round(float(np.median(values[crc])), 5),
            "median_benign": round(float(np.median(values[~is_cancer])), 5),
        }

    # The longitudinal series: one patient, three draws, a rising tumour
    # fraction. This is the comparison that cannot be moved by a cohort
    # shift, because all three samples are the same patient in the same
    # laboratory.
    draws_per_patient: dict = {}
    for s in order:
        patient = meta[s].get("patient")
        draws_per_patient[patient] = draws_per_patient.get(patient, 0) + 1
    top_patient = max(draws_per_patient, key=lambda p: draws_per_patient[p])
    series = sorted([s for s in order if meta[s].get("patient") == top_patient],
                    key=lambda s: float(meta[s]["tf"]))
    def _median_length(sample: str):
        path = os.path.join(args.features_out, f"{sample}.fsd.json")
        if not os.path.exists(path):
            return None
        with open(path) as handle:
            return json.load(handle).get("median_length")

    longitudinal = {
        "patient": meta[series[0]]["patient"] if series else None,
        "n_draws": len(series),
        "note": ("Same patient, same laboratory, rising tumour fraction. "
                 "A cohort shift cannot move this comparison."),
        "draws": [
            {
                "sample": s,
                "tumour_fraction": float(meta[s]["tf"]),
                "model_score": round(float(score[index[s]]), 6),
                "short_long_ratio": round(
                    float(extracted[s]["genome_short_long_ratio"]), 5),
                "fsd_median_length": _median_length(s),
                "fragments": extracted[s]["fsd_fragment_count"],
            }
            for s in series
        ],
    }
    if len(series) >= 3:
        values = [longitudinal["draws"][i]["model_score"]
                  for i in range(len(series))]
        longitudinal["model_score_monotone"] = all(
            values[i] < values[i + 1] for i in range(len(values) - 1))
        slr = [longitudinal["draws"][i]["short_long_ratio"]
               for i in range(len(series))]
        longitudinal["short_long_ratio_monotone"] = all(
            slr[i] < slr[i + 1] for i in range(len(slr) - 1))

    # The cancer/control contrast, reported with the caveat that it is
    # close to meaningless here: the benign samples carry a tumour
    # fraction of the same order as the lowest CRC samples.
    if is_cancer.any() and (~is_cancer).any():
        y_bin = is_cancer.astype(int)
        contrast = _auc_ci(y_bin, score)
    else:
        contrast = {"auc": None, "ci_lo": None, "ci_hi": None}
    benign_tf = [float(meta[s]["tf"]) for s in order if not is_cancer[index[s]]]
    crc_tf = [float(meta[s]["tf"]) for s in order if is_cancer[index[s]]]

    # The same four harmonisations as the PHTS cohort, so the two
    # cohorts' numbers are read the same way.
    harmonization_sensitivity = {}
    for tag, tr_z, te_z in (
        ("train_global_test_global", _zscore(Xtr), _zscore(Xte)),
        ("train_global_test_per_publication", _zscore(Xtr),
         _zscore_by_group(Xte, stte)),
        ("train_per_study_test_per_publication", _zscore_by_group(Xtr, sttr),
         _zscore_by_group(Xte, stte)),
        ("train_global_test_train_statistics", _zscore(Xtr),
         _zscore_by(Xte, Xtr)),
    ):
        pca_s = PCA(n_components=n_components, random_state=0).fit(tr_z)
        lr_s = LogisticRegression(max_iter=2000).fit(pca_s.transform(tr_z), ytr)
        sc_s = lr_s.predict_proba(pca_s.transform(te_z))[:, 1]
        block = _auc_ci(y_bin, sc_s)
        block["rho_vs_tumour_fraction"] = round(float(
            _spearman_permutation_p(tumour_fraction, sc_s, 2000,
                                    PERMUTATION_SEED)[0]), 4)
        harmonization_sensitivity[tag] = block

    prep: dict = {}
    if os.path.exists(args.prep_report):
        with open(args.prep_report, encoding="utf-8") as handle:
            prep = json.load(handle)

    # The training cohort's own library numbers, so the cohort shift is
    # visible rather than implied.
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
        "heldout_median_short_long_ratio": round(float(np.median(
            channels["genome_short_long_ratio"])), 4),
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

    verdict = {
        "question": ("Does the tumour-naive fragmentomics model track the "
                     "tumour fraction in colorectal plasma, over an "
                     "eighty-fold range of it, at the training cohort's own "
                     "depth and genome build?"),
        "rho_vs_tumour_fraction": round(rho_model, 4),
        "permutation_p": round(rho_model_p, 6),
        "rho_vs_tumour_fraction_crc_only": round(rho_crc, 4),
        "permutation_p_crc_only": round(rho_crc_p, 6),
        "tracks_the_tumour_fraction": bool(rho_model_p < 0.05 and rho_model > 0),
        "longitudinal_series_monotone": longitudinal.get(
            "model_score_monotone"),
        "channels_that_track_it": sorted(
            name for name, entry in single_channel.items()
            if entry["permutation_p"] < 0.05
            and entry["rho_vs_tumour_fraction"] > 0),
    }
    verdict["answer"] = (
        f"the score's Spearman rho against the tumour fraction is "
        f"{verdict['rho_vs_tumour_fraction']} (permutation p "
        f"{verdict['permutation_p']}), "
        + ("so it tracks it" if verdict["tracks_the_tumour_fraction"]
           else "so it does not track it at this n")
        + f"; within the CRC samples alone rho is "
        f"{verdict['rho_vs_tumour_fraction_crc_only']} (p "
        f"{verdict['permutation_p_crc_only']}); and the longitudinal series "
        f"is {'monotone' if verdict['longitudinal_series_monotone'] else 'not monotone'}")

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "design": {
            "question": "Does the tumour-naive fragmentomics model, trained "
                        "on FinaleDB publications 6+8, track the tumour "
                        "fraction in a colorectal cohort it has never seen, "
                        "and does it separate CRC from benign plasma there?",
            "train": "FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019), "
                     "n=627",
            "test": f"Zenodo {ZENODO_DOI} -- {ZENODO_RECORD}",
            "harmonization": "per study on the training side, per publication on "
                              "the held-out side",
            "pca_components": n_components,
            "model": "LogisticRegression(max_iter=2000) on the training cohort",
            "primary_readout": "the tumour-fraction gradient, not a "
                                "cancer/control AUC",
        },
        "cohort": {
            "doi": ZENODO_DOI,
            "record": ZENODO_RECORD,
            "publication_id": PUBLICATION_ID,
            "n": len(order),
            "n_crc": int(is_cancer.sum()),
            "n_benign": int((~is_cancer).sum()),
            "build": "hg38 -- the build the feature grid uses",
            "coverage": {s: meta[s].get("coverage") for s in order},
            "tumour_fraction": {s: float(meta[s]["tf"]) for s in order},
            "patient": {s: meta[s].get("patient") for s in order},
            "samples": order,
            "labels": {s: int(is_cancer[i]) for i, s in enumerate(order)},
            "fragment_counts": {s: extracted[s]["fsd_fragment_count"]
                                for s in order},
            "library": library,
            "format": prep.get("format", {}),
        },
        "result": {
            "crc_vs_benign_auc": contrast["auc"],
            "crc_vs_benign_ci": [contrast["ci_lo"], contrast["ci_hi"]],
            "crc_vs_benign_note": (
                "Reported for completeness and close to meaningless on its "
                "own: the two benign samples carry a tumour fraction of "
                f"{min(benign_tf):.4f}-{max(benign_tf):.4f} and three CRC "
                f"samples carry {min(crc_tf):.4f}-{max(crc_tf):.4f}, so the "
                "contrast is between a low tumour fraction and a slightly "
                "different low tumour fraction. The gradient and the "
                "longitudinal series are what this cohort carries."),
            "dose_response": {
                "rho_model_vs_tumour_fraction": round(rho_model, 4),
                "permutation_p": round(rho_model_p, 6),
                "rho_crc_only": round(rho_crc, 4),
                "permutation_p_crc_only": round(rho_crc_p, 6),
                "permutation_method": (
                    f"sampled permutation, {PERMUTATION_TRIALS} trials, seed "
                    f"{PERMUTATION_SEED}"),
            },
            "per_sample_score": {s: round(float(score[i]), 6)
                                 for i, s in enumerate(order)},
            "score_saturation": {
                "heldout_score_min": round(float(score.min()), 4),
                "heldout_score_max": round(float(score.max()), 4),
                "note": ("Whether the model's scores sit at the top of its "
                         "range in this cohort, as they did in the PHTS "
                         "cohort. See harmonization_sensitivity for the "
                         "control that does not let the cohort normalise "
                         "itself."),
            },
        },
        "single_channel_readout": single_channel,
        "longitudinal_series": longitudinal,
        "harmonization_sensitivity": harmonization_sensitivity,
        "verdict": verdict,
        "limits": [
            "Nine samples. The gradient's Spearman rho over nine points "
            "has a wide interval, and the longitudinal series is three "
            "draws of one patient. Both are measurements.",
            "The cancer/control contrast in this cohort is close to "
            "meaningless on its own -- the benign samples carry a tumour "
            "fraction of the same order as the lowest CRC samples -- so "
            "the AUC must not be read as a detection result.",
            "One patient supplies the longitudinal series, so the "
            "comparison cannot be generalised beyond him; what it can do "
            "is show whether the score moves with his own tumour burden.",
            "The tumour fraction is the depositor's, and it is a fraction "
            "of the reads assigned to the tumour; it is used here as the "
            "dose axis, not as a label the model is allowed to see.",
        ],
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print(f"\n  CRC vs benign AUC = {contrast['auc']} "
          f"(CI [{contrast['ci_lo']}, {contrast['ci_hi']}]) -- close to "
          f"meaningless here; the gradient is the readout")
    print(f"  dose response: rho(score, tumour fraction) = "
          f"{rho_model:.4f} (p {rho_model_p:.4f}); CRC only "
          f"{rho_crc:.4f} (p {rho_crc_p:.4f})")
    print("  per sample:")
    for i, s in enumerate(order):
        print(f"    {s:<14} tf={tumour_fraction[i]:.4f}  "
              f"label={'CRC' if is_cancer[i] else 'benign'}  "
              f"score={score[i]:.6f}  "
              f"slr={extracted[s]['genome_short_long_ratio']:.4f}")
    print(f"  verdict: {verdict['answer']}")
    print(f"\n  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
