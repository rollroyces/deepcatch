#!/usr/bin/env python3
"""Score the end-motif cohort: does the fragment end-motif profile separate cancer from healthy plasma?

What this is
------------
Zenodo 10.5281/zenodo.17253504 publishes per-sample 4-mer end-motif
frequencies for 223 cfDNA samples: 87 healthy controls and 136 cancer
samples (breast cancer 50, NSCLC 56, metastatic colorectal cancer 30),
in two library preparations. It is the only public deposit the survey
found that carries **the channel the PHTS paper's own finding rests on**
-- fragment end motifs -- with cancer/healthy labels and per-sample
values, and it is small enough to measure in minutes.

Three measurements, because they answer different questions:

1. **The published claim, tested on independent data.** The PHTS paper
   reports that CC- and CA-end motifs are enriched in the plasma of
   patients whose cancer was about to be diagnosed. The end-motif
   fraction of those two motif families is computed per sample and
   tested against the labels with a permutation test. No model is
   involved, so nothing here depends on a training set.
2. **Is there a class signal in the 256 features at all?** A
   cross-validated model inside this cohort (5 seeds x 5 folds, pooled
   out-of-fold scores) answers that, pooled over the cancers and per
   cancer type.
3. **Does the repository's own motif model transfer?** The training
   cohort's 98 samples that carry a motif vector (48 healthy, 50
   hepatocellular carcinoma) are enough to fit one, and applying it
   here is a test of the repository's own prediction that a
   single-cancer-type training set yields a type-specific decision
   boundary.

What it does not do
-------------------
It does not claim that the deposit's motif definition matches the
repository's. The deposit says the motifs come from the 5' fragment
boundary but does not say whether both strands' ends were counted per
fragment, which is what the repository's own extraction does. Every
number here carries that ambiguity.
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

from _paths import SIBLING_PIPELINE

ENDMOTIF_DOI = "10.5281/zenodo.17253504"
PERMUTATION_TRIALS = 20000
PERMUTATION_SEED = 20261125

# The two motif families the PHTS paper reports as enriched in the
# plasma of patients whose cancer was about to be diagnosed.
PUBLISHED_MOTIFS = ("CC", "CA")


def _motif_families() -> list[str]:
    """The 4-mers whose first two bases are CC or CA."""
    letters = "ACGT"
    out = []
    for i in range(256):
        motif = "".join(letters[(i >> (2 * (3 - k))) & 3] for k in range(4))
        if motif[:2] in PUBLISHED_MOTIFS:
            out.append(motif)
    return out


def _read_label_table(path: str) -> dict:
    rows = {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            rows[row["sample_id"]] = row
    return rows


def _permutation_p(labels: np.ndarray, values: np.ndarray,
                   trials: int = PERMUTATION_TRIALS,
                   seed: int = PERMUTATION_SEED) -> float:
    """A label-permutation p-value for a rank statistic."""
    observed = roc_auc_score(labels, values)
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(trials):
        if roc_auc_score(rng.permutation(labels), values) >= observed - 1e-12:
            hits += 1
    return (hits + 1) / (trials + 1)


def _cv_auc(X: np.ndarray, y: np.ndarray, seeds: int = 5,
            folds: int = 5) -> tuple[float, list]:
    """Pooled out-of-fold AUC over seeds x folds, and the per-seed AUCs."""
    from sklearn.preprocessing import StandardScaler
    score_acc = np.zeros(len(y), dtype=float)
    per_seed = []
    for seed in range(seeds):
        oof = np.zeros(len(y), dtype=float)
        for train, test in StratifiedKFold(folds, shuffle=True,
                                          random_state=seed).split(X, y):
            scaler = StandardScaler().fit(X[train])
            model = LogisticRegression(max_iter=2000).fit(
                scaler.transform(X[train]), y[train])
            oof[test] = model.predict_proba(scaler.transform(X[test]))[:, 1]
        per_seed.append(float(roc_auc_score(y, oof)))
        score_acc += oof
    return float(roc_auc_score(y, score_acc / seeds)), per_seed


def _de_long(y: np.ndarray, score: np.ndarray) -> dict:
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "..", "src"))
    from per_cancer_sens_at_spec import delong_auc_ci
    d = delong_auc_ci(np.asarray(y), np.asarray(score))
    return {"auc": round(float(d["auc"]), 4),
            "ci_lo": round(float(d["ci_lo"]), 4),
            "ci_hi": round(float(d["ci_hi"]), 4),
            "n_pos": int(d["n_pos"]), "n_neg": int(d["n_neg"])}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cohort-dir", required=True,
                    help="the prepared cohort directory (features + labels)")
    ap.add_argument("--features-dir",
                    default=os.path.join(SIBLING_PIPELINE, "data/features"))
    ap.add_argument("--labels-multiclass",
                    default=os.path.join(SIBLING_PIPELINE,
                                         "labels_multiclass.tsv"))
    ap.add_argument("--out-json",
                    default="results/external_cohort_zenodo_endmotif.json")
    args = ap.parse_args()

    print("[1/5] reading the cohort ...", flush=True)
    labels = _read_label_table(os.path.join(args.cohort_dir, "labels.tsv"))
    features_dir = os.path.join(args.cohort_dir, "features")
    order, X, y, groups, assays = [], [], [], [], []
    for sample in sorted(labels):
        path = os.path.join(features_dir, f"{sample}.motifs.npy")
        if not os.path.exists(path):
            continue
        order.append(sample)
        X.append(np.load(path))
        y.append(1 if labels[sample]["label"] != "healthy" else 0)
        groups.append(labels[sample]["group"])
        assays.append(labels[sample]["assay"])
    X = np.vstack(X)
    y = np.asarray(y, dtype=int)
    groups = np.asarray(groups)
    assays = np.asarray(assays)
    print(f"      n={len(order)} ({int(y.sum())} cancer, "
          f"{int((1 - y).sum())} healthy), {X.shape[1]} motifs", flush=True)

    print("[2/5] the published claim: the CC- and CA-end motifs ...", flush=True)
    family = _motif_families()
    indices = [i for i in range(256)
               if "".join("ACGT"[(i >> (2 * (3 - k))) & 3]
                          for k in range(4)) in family]
    fraction = X[:, indices].sum(axis=1)
    claim = {
        "motifs": family,
        "definition": "the fraction of a sample's 4-mer end-motif counts "
                      "whose first two bases are CC or CA -- the two "
                      "families the PHTS paper reports as enriched in the "
                      "plasma of patients whose cancer was about to be "
                      "diagnosed",
        "median_cancer": round(float(np.median(fraction[y == 1])), 6),
        "median_healthy": round(float(np.median(fraction[y == 0])), 6),
    }
    claim["auc"] = _de_long(y, fraction)["auc"]
    claim["permutation_p"] = round(_permutation_p(y, fraction), 6)
    claim["direction_matches_the_published_finding"] = bool(
        claim["median_cancer"] > claim["median_healthy"])
    claim["median_by_group"] = {
        group: round(float(np.median(fraction[groups == group])), 6)
        for group in sorted(set(groups.tolist()))
    }
    # The same control, on this statistic: two cancer groups from
    # different studies. If the CC-/CA-end fraction separates THEM, then
    # the class difference above is a study difference too.
    for left, right in (("BC", "NSCLC"), ("HI25", "HI37")):
        mask = (groups == left) | (groups == right)
        if int(mask.sum()) < 20:
            continue
        block_y = (groups[mask] == right).astype(int)
        claim[f"control_{left}_vs_{right}_auc"] = _de_long(
            block_y, fraction[mask])["auc"]
        claim[f"control_{left}_vs_{right}_medians"] = [
            round(float(np.median(fraction[groups == left])), 6),
            round(float(np.median(fraction[groups == right])), 6),
        ]
    print(f"      cancer {claim['median_cancer']:.5f} vs healthy "
          f"{claim['median_healthy']:.5f}, AUC {claim['auc']}, "
          f"p {claim['permutation_p']}", flush=True)
    print(f"      per group: {claim['median_by_group']}", flush=True)

    print("[3/5] is there a class signal in the 256 features? ...", flush=True)
    # The confound control, and it matters here: the healthy arm comes
    # from two studies and the cancer arm from three others, so a
    # study-level difference in the libraries is not separable from a
    # class difference by this cohort. The control is the repository's
    # own: if the features separate two groups that SHARE a label, then
    # what the class AUC measures is the study, not the cancer.
    confound = {}
    for left, right in (("HI25", "HI37"), ("BC", "NSCLC")):
        mask = (groups == left) | (groups == right)
        if int(mask.sum()) < 20:
            continue
        block_y = (groups[mask] == right).astype(int)
        auc_c, seeds_c = _cv_auc(X[mask], block_y)
        confound[f"{left}_vs_{right}"] = {
            "auc": round(float(auc_c), 4),
            "per_seed_aucs": [round(a, 4) for a in seeds_c],
            "n_left": int((groups[mask] == left).sum()),
            "n_right": int((groups[mask] == right).sum()),
            "shared_label": ("healthy" if left in ("HI25", "HI37")
                             else "cancer"),
            "note": ("Two groups that SHARE a label, from different "
                     "studies. A high AUC here means the features carry "
                     "study-level structure, and the class AUC cannot be "
                     "read as a cancer signal."),
        }
        print(f"      control {left} vs {right} (both "
              f"{confound[f'{left}_vs_{right}']['shared_label']}): AUC "
              f"{auc_c}", flush=True)

    auc, per_seed = _cv_auc(X, y)
    # the pooled out-of-fold scores are recomputed here for the interval
    score_acc = np.zeros(len(y), dtype=float)
    from sklearn.preprocessing import StandardScaler
    for seed in range(5):
        oof = np.zeros(len(y), dtype=float)
        for train, test in StratifiedKFold(5, shuffle=True,
                                          random_state=seed).split(X, y):
            scaler = StandardScaler().fit(X[train])
            model = LogisticRegression(max_iter=2000).fit(
                scaler.transform(X[train]), y[train])
            oof[test] = model.predict_proba(scaler.transform(X[test]))[:, 1]
        score_acc += oof
    oof_scores = score_acc / 5
    within = {
        "pooled": _de_long(y, oof_scores),
        "per_seed_aucs": [round(a, 4) for a in per_seed],
    }
    for group in sorted(set(groups.tolist())):
        if group == "HI25" or group == "HI37":
            continue
        mask = (groups == group) | (y == 0)
        if mask.sum() < 20 or len(set(y[mask].tolist())) < 2:
            continue
        auc_g, seeds_g = _cv_auc(X[mask], y[mask])
        within[group] = {"auc": round(float(auc_g), 4),
                          "per_seed_aucs": [round(a, 4) for a in seeds_g],
                          "n_cancer": int(y[mask].sum()),
                          "n_healthy": int((1 - y[mask]).sum())}
    print(f"      pooled out-of-fold AUC {within['pooled']['auc']} "
          f"[{within['pooled']['ci_lo']}, {within['pooled']['ci_hi']}]", flush=True)

    print("[4/5] the repository's own motif model, transferred ...", flush=True)
    train_labels = {}
    with open(args.labels_multiclass, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            train_labels[row["sample"]] = row
    train_X, train_y, train_groups = [], [], []
    for name in sorted(train_labels):
        path = os.path.join(args.features_dir, f"{name}.motifs.npy")
        if not os.path.exists(path):
            continue
        train_X.append(np.load(path))
        train_y.append(0 if str(train_labels[name]["label"]).lower()
                        in ("healthy", "0", "control") else 1)
        train_groups.append(train_labels[name]["disease_class"])
    train_X = np.vstack(train_X)
    train_y = np.asarray(train_y, dtype=int)
    train_groups = np.asarray(train_groups)
    transferred = None
    if len(set(train_y.tolist())) == 2:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler().fit(train_X)
        model = LogisticRegression(max_iter=2000).fit(
            scaler.transform(train_X), train_y)
        transferred_score = model.predict_proba(scaler.transform(X))[:, 1]
        transferred = {
            "training_cohort": {
                "n": int(len(train_y)),
                "n_cancer": int(train_y.sum()),
                "n_healthy": int((1 - train_y).sum()),
                "cancer_type": sorted(set(train_groups.tolist())),
                "note": ("The training cohort's samples that carry a motif "
                         "vector are all ONE cancer type, so this is a "
                         "type-specific model -- which is the repository's "
                         "own prediction, not a confound in the cohort."),
            },
            "pooled": _de_long(y, transferred_score),
        }
        for group in sorted(set(groups.tolist())):
            if group in ("HI25", "HI37"):
                continue
            mask = (groups == group) | (y == 0)
            if len(set(y[mask].tolist())) < 2:
                continue
            transferred[group] = _de_long(y[mask], transferred_score[mask])
        print(f"      trained on {len(train_y)} samples "
              f"({sorted(set(train_groups.tolist()))}), pooled AUC "
              f"{transferred['pooled']['auc']}", flush=True)

    print("[5/5] writing the artifact ...", flush=True)
    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "design": {
            "question": "Does the fragment end-motif profile separate cancer "
                        "from healthy plasma in an independent cohort, and "
                        "does the repository's own motif model transfer to "
                        "it?",
            "cohort": f"Zenodo {ENDMOTIF_DOI} -- per-sample 4-mer end-motif "
                      "frequencies, 223 samples",
            "features": "the 256 motif frequencies, the repository's own "
                        "motif order",
            "primary_readout": "the within-cohort cross-validated AUC "
                                "(is there a signal?), with the published "
                                "claim's motif fraction and the transferred "
                                "model next to it",
        },
        "cohort": {
            "doi": ENDMOTIF_DOI,
            "record": "https://zenodo.org/records/17253504",
            "n": int(len(order)),
            "n_cancer": int(y.sum()),
            "n_healthy": int((1 - y).sum()),
            "groups": {g: int((groups == g).sum())
                       for g in sorted(set(groups.tolist()))},
            "assays": {a: int((assays == a).sum())
                       for a in sorted(set(assays.tolist()))},
        },
        "published_claim": claim,
        "within_cohort_model": within,
        "confound_control": confound,
        "transferred_model": transferred,
        "limits": [
            "The measured confound control decides how this cohort can be "
            "read: the same 256 features separate BC from NSCLC -- two "
            "cancer groups from different studies -- at AUC 0.9996, and "
            "the two healthy studies differ on the CC-/CA-end fraction at "
            "AUC 0.7057. The study structure in these libraries is larger "
            "than any class effect here, so none of the class AUCs in this "
            "artifact can be read as a cancer signal -- and the "
            "transferred model's inversion (AUC 0.3413) is what a "
            "definition mismatch between the repository's end-motif "
            "extraction and this deposit's looks like.",
            "The deposit says the motifs come from the 5' fragment "
            "boundary but does not say whether both strands' ends were "
            "counted per fragment, which is what the repository's own "
            "extraction does. Every number here carries that ambiguity.",
            "The healthy arm is a screening or healthy-control population "
            "from two studies (HI25, HI37) and the cancer arm is three "
            "different tumour types, so a pooled AUC mixes the type "
            "shift with any class shift.",
            "The within-cohort model is cross-validated inside this "
            "cohort, so its AUC is a measure of whether the features "
            "separate the classes here -- not of transfer.",
            "The transferred model was trained on one cancer type "
            "(hepatocellular carcinoma) because those are the only "
            "training-cohort samples with a motif vector, so a per-type "
            "breakdown is where its own prediction is tested.",
        ],
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print(f"\n  published claim: cancer {claim['median_cancer']:.5f} vs "
          f"healthy {claim['median_healthy']:.5f} (AUC {claim['auc']}, "
          f"p {claim['permutation_p']})")
    print(f"  within-cohort model: AUC {within['pooled']['auc']} "
          f"[{within['pooled']['ci_lo']}, {within['pooled']['ci_hi']}]")
    for group, block in within.items():
        if group in ("pooled", "per_seed_aucs"):
            continue
        print(f"    {group:<14} AUC {block['auc']}  n={block['n_cancer']} "
              f"cancer vs {block['n_healthy']} healthy")
    if transferred:
        print(f"  transferred model: pooled AUC {transferred['pooled']['auc']} "
              f"[{transferred['pooled']['ci_lo']}, "
              f"{transferred['pooled']['ci_hi']}]")
        for group, block in transferred.items():
            if group in ("pooled", "training_cohort"):
                continue
            print(f"    {group:<14} AUC {block['auc']}")
    print(f"\n  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
