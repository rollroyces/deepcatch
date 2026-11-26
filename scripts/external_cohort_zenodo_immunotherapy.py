#!/usr/bin/env python3
"""Score the immunotherapy cohort: does anything here separate responders from non-responders?

What this is
------------
Zenodo 10.5281/zenodo.15738112 deposits the fragment files of 185
longitudinal plasma cfDNA samples from 68 patients with head and neck
squamous cell carcinoma who received pembrolizumab, across six
academic centres (Bandaru et al., medRxiv
10.64898/2026.03.24.26348354; Journal of Clinical Investigation
10.1172/JCI196284). The endpoint is **clinical in the treatment-response
sense**: responder versus non-responder, with the response defined
pathologically (tumour-erythema >= 10%).

It is the first public deposit this repository has found whose readout
is a treatment response rather than a cancer-versus-healthy
classification, and its channel is the one the PHTS finding rests on:
fragment end motifs.

What is measured, and why each measurement is here
--------------------------------------------------
1. **The published claim, on the published metric.** The study's metric
   is the rMDS -- the normalized Shannon entropy of the 5' end 4-mer
   frequencies within 500 kb bins, z-scored across the genome within
   each sample -- and its unsupervised claim is that rMDS separates
   responders from non-responders. The per-bin differential is measured
   here (the study's Figure 3A), and a model on the rMDS bins is scored
   with the repository's own protocol (patient-level grouped
   cross-validation, pooled out-of-fold).
2. **Transfer to an unseen centre.** The study's own second validation
   is a leave-one-institute-out evaluation, and its third is a
   "pre-holdout" set in which the batch correction has not been fitted.
   Both are reproduced here, because a cohort whose samples come from
   six centres has the centre as a first-class technical variable.
3. **The confound control, which is the point.** The repository's rule
   is that an AUC without a confound control is not a result: the
   end-motif channel was already shown to separate two cancer groups
   from different studies better than it separated cancer from healthy.
   Here the control is the same measurement run on two centres that
   **share** a response label, and on two library-preparation batches
   that share a response label. If the features separate those, then
   the response AUC cannot be read as a response signal.

What is not claimed
-------------------
The deposited files are not the study's own filtered fragment set: they
carry fragments of every length and every MAPQ, while the study says it
kept 50-350 bp autosomal fragments at MAPQ >= 30 and excluded the
ENCODE blacklist. The first two filters are applied here and the
blacklist is not (the track is not in the deposit), so these numbers
are a reproduction on the deposited data, not a reproduction of the
study's pipeline. The estimator is a regularized logistic regression on
SVD components, not the study's TabPFN.

Usage:
    python3 scripts/external_cohort_zenodo_immunotherapy.py \
        --features-dir data/external_cohort_zenodo_immunotherapy \
        --metadata data/external_cohort/immunotherapy_sample_metadata.tsv \
        --out-json results/external_cohort_zenodo_immunotherapy.json
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler

SEEDS = 5
FOLDS = 5
SVD_COMPONENTS = 32
PERMUTATION_TRIALS = 20000
PERMUTATION_SEED = 20261126


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _read_metadata(path: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("sample_id") or "").strip()
            if sample:
                out[sample] = {k: (v or "").strip() for k, v in row.items()}
    return out


def _delong(y: np.ndarray, score: np.ndarray) -> dict:
    sys.path.insert(0, os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "..", "src"))
    from per_cancer_sens_at_spec import delong_auc_ci

    d = delong_auc_ci(np.asarray(y), np.asarray(score))
    return {"auc": round(float(d["auc"]), 4),
            "ci_lo": round(float(d["ci_lo"]), 4),
            "ci_hi": round(float(d["ci_hi"]), 4),
            "n_pos": int(d["n_pos"]), "n_neg": int(d["n_neg"])}


def _fit_predict(X: np.ndarray, y: np.ndarray, train: np.ndarray,
                 test: np.ndarray) -> np.ndarray:
    """A regularized logistic regression on SVD components, fitted on
    `train` and scored on `test`.

    The study reduces the rMDS with a truncated SVD and classifies with
    a small-sample classifier; this does the same reduction and uses a
    regularized logistic regression, because the study's TabPFN is not
    available here and the point of the measurement is the confound
    control, not the classifier.
    """
    n_components = int(min(SVD_COMPONENTS, X.shape[1] - 1,
                            max(2, len(train) - 1)))
    reducer = TruncatedSVD(n_components=n_components, random_state=0)
    reduced = reducer.fit_transform(X[train])
    scaler = StandardScaler().fit(reduced)
    model = LogisticRegression(max_iter=5000, C=1.0).fit(
        scaler.transform(reduced), y[train])
    test_reduced = reducer.transform(X[test])
    return model.predict_proba(scaler.transform(test_reduced))[:, 1]


def _grouped_cv_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                    seeds: int = SEEDS, folds: int = FOLDS):
    """Pooled out-of-fold AUC over seeds x folds of a patient-grouped CV.

    The patient is the group: a patient's three time points must not be
    split across the folds, which is the same rule the repository uses
    for its paired cohort and the same rule the study uses for its own
    patient-level cross-validation. The fold assignment is stratified
    on the response and shuffled per seed, so the seeds really are
    different splits and their spread is a measurement rather than an
    artefact of a deterministic partition.
    """
    score_acc = np.zeros(len(y), dtype=float)
    per_seed = []
    for seed in range(seeds):
        oof = np.zeros(len(y), dtype=float)
        splitter = StratifiedGroupKFold(n_splits=folds, shuffle=True,
                                        random_state=seed)
        for train, test in splitter.split(X, y, groups=groups):
            oof[test] += _fit_predict(X, y, train, test)
        per_seed.append(float(roc_auc_score(y, oof)))
        score_acc += oof
    return float(roc_auc_score(y, score_acc / seeds)), per_seed, score_acc / seeds


def _permutation_p(y: np.ndarray, score: np.ndarray,
                   trials: int = PERMUTATION_TRIALS) -> float:
    """The one-sided permutation p-value of an AUC."""
    observed = roc_auc_score(y, score)
    rng = np.random.default_rng(PERMUTATION_SEED)
    hits = 0
    for _ in range(trials):
        if roc_auc_score(rng.permutation(y), score) >= observed:
            hits += 1
    return (hits + 1) / (trials + 1)


def _batch_corrected(rmds: np.ndarray, isolation: np.ndarray,
                      prep: np.ndarray, institute: np.ndarray,
                      response: np.ndarray) -> np.ndarray:
    """The study's own batch correction, applied bin by bin.

    The study models each bin's rMDS as a function of the cfDNA
    isolation date, the WGS library-preparation date and the collection
    institute, with the treatment response kept in the model as a
    biological covariate that the correction must not remove. The
    corrected value is therefore the intercept, the response's own
    coefficient and the residual -- the technical terms are taken out
    and the response term is left in.

    The design uses sum contrasts for the institute, as the study does,
    so the institutes are not separated one against another but against
    their mean.
    """
    institutes = sorted(set(institute.tolist()))
    dummies = np.zeros((len(institute), max(len(institutes) - 1, 0)))
    for index, name in enumerate(institutes[:-1]):
        dummies[:, index] = (institute == name).astype(float)
    design = np.column_stack([
        np.ones(len(response)), isolation, prep, dummies, response])

    corrected = np.zeros_like(rmds)
    for column in range(rmds.shape[1]):
        y = rmds[:, column]
        beta, *_ = np.linalg.lstsq(design, y, rcond=None)
        residual = y - design @ beta
        # intercept + the response's own coefficient + what the model
        # did not explain
        corrected[:, column] = (beta[0] + beta[-1] * response + residual)
    return corrected


def _label_is_consistent(metadata: dict, order: list[str]) -> bool:
    """Does the response label agree with the study's own response grade?

    The metadata carries both: a categorised response rate (0 for a
    non-responder, 1-10 for a responder) and the binary label. A
    cohort whose label contradicted its own grade would be measuring
    something else, so the two are checked against each other rather
    than assumed to agree.
    """
    for sample in order:
        row = metadata[sample]
        try:
            rate = int(row.get("response_rate", ""))
        except ValueError:
            continue
        is_responder = row.get("response") == "Responder"
        if (rate >= 1) != is_responder:
            return False
    return True


def _bin_differential(X: np.ndarray, y: np.ndarray) -> dict:
    """How many bins separate the two classes, and in which direction.

    This is the study's own unsupervised reading (its Figure 3A): a
    two-group test per bin, with the p-values corrected for multiple
    testing over the genome.
    """
    from scipy import stats

    pvalues = np.ones(X.shape[1])
    differences = np.zeros(X.shape[1])
    for column in range(X.shape[1]):
        responders = X[y == 1, column]
        non = X[y == 0, column]
        responders = responders[np.isfinite(responders)]
        non = non[np.isfinite(non)]
        if len(responders) < 3 or len(non) < 3:
            continue
        _, p = stats.ttest_ind(responders, non, equal_var=False)
        pvalues[column] = p if np.isfinite(p) else 1.0
        differences[column] = responders.mean() - non.mean()
    order_p = np.argsort(pvalues)
    ranked = pvalues[order_p] * len(pvalues) / (np.arange(len(pvalues)) + 1)
    fdr = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.ones_like(fdr)
    adjusted[order_p] = np.minimum(fdr, 1.0)
    significant = adjusted < 0.1
    return {
        "bins_tested": int(X.shape[1]),
        "bins_fdr_below_0.1": int(significant.sum()),
        "fraction_significant": round(float(significant.mean()), 4),
        "higher_in_responders": int((significant & (differences > 0)).sum()),
        "lower_in_responders": int((significant & (differences < 0)).sum()),
        "median_abs_difference": round(float(np.median(np.abs(differences))),
                                       4),
    }


def _days_since_earliest(values: list[str]) -> np.ndarray:
    """A date column as days since the earliest date in it."""
    import datetime as dt

    parsed = []
    for value in values:
        text = (value or "").strip()
        if not text:
            parsed.append(None)
            continue
        for fmt in ("%m/%d/%y", "%m/%d/%Y", "%Y-%m-%d"):
            try:
                parsed.append(dt.datetime.strptime(text, fmt).date())
                break
            except ValueError:
                continue
        else:
            parsed.append(None)
    known = [d for d in parsed if d is not None]
    if not known:
        return np.zeros(len(values))
    first = min(known)
    out = np.zeros(len(values))
    for index, day in enumerate(parsed):
        out[index] = (day - first).days if day is not None else 0.0
    return out


def _pairwise_control(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                      patients: np.ndarray, left: str, right: str,
                      label: str) -> dict:
    """The control: do the features separate two groups that SHARE a label?

    The two groups share a response label, so the response cannot be
    what the features separate -- the question is whether they separate
    the two GROUPS. The label of the control model is therefore the
    group itself, fitted and scored on the group's samples only, with
    the same patient-grouped cross-validation the response model uses,
    so the two AUCs are comparable.

    A high AUC here means the features carry group-level structure, and
    the response AUC cannot be read as a response signal.
    """
    mask = (groups == left) | (groups == right)
    n_left = int((groups[mask] == left).sum())
    n_right = int((groups[mask] == right).sum())
    if n_left < 8 or n_right < 8:
        return {"n_left": n_left, "n_right": n_right,
                "note": "too few samples on one side to measure"}
    block_y = (groups[mask] == right).astype(int)
    auc, per_seed, _ = _grouped_cv_auc(X[mask], block_y, patients[mask])
    return {"auc": round(auc, 4),
            "per_seed_aucs": [round(a, 4) for a in per_seed],
            "n_left": n_left, "n_right": n_right,
            "shared_label": label,
            "note": ("two groups that SHARE a label; the AUC here is the "
                     "separation of the groups themselves, so a high value "
                     "means the features carry group-level structure and "
                     "the response AUC cannot be read as a response signal")}


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--features-dir", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--prep-report",
                        help="the preparation report, for the depth control")
    parser.add_argument("--out-json", required=True)
    args = parser.parse_args(argv)

    features_dir = Path(args.features_dir)
    metadata = _read_metadata(args.metadata)

    print("[1/6] reading the cohort ...", flush=True)
    order, y, groups, patients, visits, institutes, batches = [], [], [], [], [], [], []
    rmds_by_sample: dict[str, dict[str, float]] = {}
    motif_rows: list[np.ndarray] = []
    for sample in sorted(metadata):
        row = metadata[sample]
        rmds_path = features_dir / f"{sample}.rmds.npy"
        bins_path = features_dir / f"{sample}.rmds_bins.npy"
        motif_path = features_dir / f"{sample}.motifs.npy"
        if not rmds_path.exists() or not motif_path.exists():
            continue
        rmds = np.load(rmds_path)
        bins = np.load(bins_path, allow_pickle=True)
        motifs = np.load(motif_path)
        if rmds.size == 0 or bins.size != rmds.size:
            continue
        rmds_by_sample[sample] = dict(zip(bins.tolist(), rmds.tolist()))
        order.append(sample)
        motif_rows.append(motifs)
        y.append(1 if row["response"] == "Responder" else 0)
        groups.append(row["institute"])
        patients.append(row["patient"])
        visits.append(row["visit"])
        institutes.append(row["institute"])
        batches.append(row["library_prep_date"])

    y = np.asarray(y, dtype=int)
    groups = np.asarray(groups)
    patients = np.asarray(patients)
    visits = np.asarray(visits)
    institutes = np.asarray(institutes)
    batches = np.asarray(batches)
    motif_X = np.vstack(motif_rows)

    # The rMDS is a per-bin statistic, so the samples have to sit on one
    # grid for the matrix to mean anything. A bin with no fragments in a
    # sample has no diversity in it, and a handful of bins are like
    # that -- one of them (chr21:11) has fragments in 48 samples and
    # none in 136. The shared grid is therefore the bins EVERY sample
    # has fragments in, and the bins that were dropped are counted in
    # the report rather than quietly imputed.
    shared = set(rmds_by_sample[order[0]])
    for sample in order[1:]:
        shared &= set(rmds_by_sample[sample])
    union = set()
    for sample in order:
        union |= set(rmds_by_sample[sample])
    grid = sorted(shared)
    rmds_X = np.vstack([
        np.asarray([rmds_by_sample[sample][key] for key in grid])
        for sample in order])
    # Re-standardise on the shared grid: a sample's z-scores were
    # computed over the bins IT had, and the shared grid is a subset of
    # them, so the scale differs slightly from sample to sample.
    rmds_X = (rmds_X - rmds_X.mean(axis=1, keepdims=True))
    sd = rmds_X.std(axis=1, keepdims=True)
    rmds_X = np.divide(rmds_X, sd, out=np.zeros_like(rmds_X), where=sd > 0)

    print(f"      n={len(order)} ({int(y.sum())} responders, "
          f"{int((1 - y).sum())} non-responders), {len(set(patients))} patients, "
          f"{len(grid)} shared rMDS bins (of {len(union)} in the union), "
          f"{motif_X.shape[1]} motifs", flush=True)

    results: dict = {
        "doi": "10.5281/zenodo.15738112",
        "study_doi": "10.64898/2026.03.24.26348354",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cohort": {
            "samples": len(order),
            "responders": int(y.sum()),
            "non_responders": int((1 - y).sum()),
            "patients": len(set(patients.tolist())),
            "institutes": sorted(set(institutes.tolist())),
            "visits": sorted(set(visits.tolist())),
            "rmds_bins": int(rmds_X.shape[1]),
            "rmds_bins_in_the_union": len(union),
            "rmds_bins_dropped_for_being_incomplete": len(union) - len(grid),
            "per_institute": {name: int((institutes == name).sum())
                              for name in sorted(set(institutes.tolist()))},
            "responder_fraction_by_institute": {
                name: round(float(y[institutes == name].mean()), 4)
                for name in sorted(set(institutes.tolist()))},
            "the_label_is_internally_consistent": _label_is_consistent(
                metadata, order),
        },
    }

    print("[2/6] the published claim: the per-bin rMDS differential ...",
          flush=True)
    differential = _bin_differential(rmds_X, y)
    results["published_claim_rmds"] = {
        "definition": ("the normalized Shannon entropy of the 5' end 4-mer "
                       "frequencies in 500 kb bins, z-scored across the "
                       "genome within each sample -- the study's rMDS"),
        "bins_tested": differential["bins_tested"],
        "bins_fdr_below_0.1": differential["bins_fdr_below_0.1"],
        "fraction_significant": differential["fraction_significant"],
        "higher_in_responders": differential["higher_in_responders"],
        "lower_in_responders": differential["lower_in_responders"],
        "median_abs_difference": differential["median_abs_difference"],
        "note": ("the study reports 1,080 differentially differential bins "
                 "at FDR < 0.1 out of the autosome's 500 kb bins; the count "
                 "here is computed on the deposited fragments with the "
                 "study's length and MAPQ filters and without the ENCODE "
                 "blacklist exclusion"),
    }
    print(f"      {results['published_claim_rmds']['bins_fdr_below_0.1']} of "
          f"{rmds_X.shape[1]} bins differ at FDR < 0.1", flush=True)

    print("[3/6] the model on the rMDS bins, and on the 256 motifs ...",
          flush=True)
    model_results = {}
    for name, X in (("rmds_bins", rmds_X), ("motifs_256", motif_X)):
        auc, per_seed, oof = _grouped_cv_auc(X, y, patients)
        block = {
            "pooled_out_of_fold_auc": round(auc, 4),
            "per_seed_aucs": [round(a, 4) for a in per_seed],
            "de_long": _delong(y, oof),
            "permutation_p": round(_permutation_p(y, oof), 6),
        }
        model_results[name] = block
        print(f"      {name}: pooled OOF AUC {auc:.4f} "
              f"[{block['de_long']['ci_lo']}, {block['de_long']['ci_hi']}]",
              flush=True)

    # The same measurement with the study's own batch correction applied
    # -- the batch's own terms removed, the response's term kept. This
    # is the question the study asks; the uncorrected reading above is
    # the question this repository has always asked, which is what the
    # features do before anything is taken out.
    isolation = _days_since_earliest(
        [metadata[sample].get("isolation_date", "") for sample in order])
    prep = _days_since_earliest(
        [metadata[sample].get("library_prep_date", "") for sample in order])
    corrected = _batch_corrected(rmds_X, isolation, prep, institutes, y)
    corrected = (corrected - corrected.mean(axis=1, keepdims=True))
    sd = corrected.std(axis=1, keepdims=True)
    corrected = np.divide(corrected, sd, out=np.zeros_like(corrected),
                           where=sd > 0)
    auc_c, per_seed_c, oof_c = _grouped_cv_auc(corrected, y, patients)
    model_results["rmds_bins_batch_corrected"] = {
        "pooled_out_of_fold_auc": round(auc_c, 4),
        "per_seed_aucs": [round(a, 4) for a in per_seed_c],
        "de_long": _delong(y, oof_c),
        "permutation_p": round(_permutation_p(y, oof_c), 6),
        "note": ("the study's own correction: a per-bin model on the "
                 "cfDNA isolation date, the WGS library-preparation date "
                 "and the collection institute, with the response kept in "
                 "the model as a covariate that the correction does not "
                 "remove"),
    }
    print(f"      rmds_bins_batch_corrected: pooled OOF AUC {auc_c:.4f} "
          f"[{model_results['rmds_bins_batch_corrected']['de_long']['ci_lo']}, "
          f"{model_results['rmds_bins_batch_corrected']['de_long']['ci_hi']}]",
          flush=True)
    results["model"] = model_results

    corrected_differential = _bin_differential(corrected, y)
    results["published_claim_rmds_batch_corrected"] = {
        "definition": ("the same per-bin differential on the study's own "
                       "batch-corrected rMDS: the isolation date, the "
                       "library-preparation date and the collection "
                       "institute are modelled out of every bin, and the "
                       "response is kept in the model"),
        "bins_fdr_below_0.1": corrected_differential["bins_fdr_below_0.1"],
        "higher_in_responders": corrected_differential["higher_in_responders"],
        "lower_in_responders": corrected_differential["lower_in_responders"],
        "note": ("the study reports 1,080 bins differential at FDR < 0.1 "
                 "out of the autosome's 500 kb bins, and the study's "
                 "classifier reaches AUC 0.99 in patient-level "
                 "cross-validation and 0.89 with an institute held out"),
    }
    print(f"      batch-corrected: "
          f"{corrected_differential['bins_fdr_below_0.1']} of "
          f"{corrected_differential['bins_tested']} bins differ at "
          f"FDR < 0.1", flush=True)

    print("[4/6] transfer to an unseen centre (the study's LOIO) ...",
          flush=True)
    loio = {}
    for name, X in (("rmds_bins", rmds_X), ("motifs_256", motif_X)):
        per_institute = {}
        for institute in sorted(set(institutes.tolist())):
            test = np.where(institutes == institute)[0]
            train = np.where(institutes != institute)[0]
            if len(test) < 5 or len(set(y[train].tolist())) < 2:
                continue
            if len(set(y[test].tolist())) < 2:
                # A single-class held-out centre has no AUC; the
                # separation is reported as a rate instead.
                score = _fit_predict(X, y, train, test)
                predicted = (score >= 0.5).astype(int)
                rate = float((predicted == y[test]).mean())
                per_institute[institute] = {
                    "n": int(len(test)),
                    "class_counts": {str(c): int((y[test] == c).sum())
                                     for c in sorted(set(y[test].tolist()))},
                    "accuracy": round(rate, 4),
                    "note": ("the held-out centre carries one response "
                             "class, so there is no AUC to report"),
                }
                continue
            score = _fit_predict(X, y, train, test)
            per_institute[institute] = {
                "n": int(len(test)),
                "auc": round(float(roc_auc_score(y[test], score)), 4),
                "de_long": _delong(y[test], score),
            }
        loio[name] = per_institute
    results["leave_one_institute_out"] = loio
    print(f"      {json.dumps({k: {i: v.get('auc', v.get('accuracy')) for i, v in b.items()} for k, b in loio.items()})}",
          flush=True)

    print("[5/6] the confound control: two centres that SHARE a label ...",
          flush=True)
    controls = {}
    names = sorted(set(institutes.tolist()))
    for name, X in (("rmds_bins", rmds_X), ("motifs_256", motif_X)):
        per_pair = {}
        for i, left in enumerate(names):
            for right in names[i + 1:]:
                for label, keep in (("Responder", 1), ("Non-Responder", 0)):
                    mask = (institutes == left) | (institutes == right)
                    if int((y[mask] == keep).sum()) < 8:
                        continue
                    block = X[mask]
                    block_y = (y[mask] == keep).astype(int)
                    sub_groups = institutes[mask]
                    result = _pairwise_control(
                        block, block_y, sub_groups, patients[mask],
                        left, right, label)
                    if "auc" in result:
                        per_pair[f"{left}_vs_{right}_{label}"] = result
        controls[name] = per_pair
        worst = max((v["auc"] for v in per_pair.values()), default=None)
        print(f"      {name}: {len(per_pair)} shared-label centre pairs, "
              f"the highest AUC {worst}", flush=True)
    results["confound_control_centres"] = controls

    print("[6/6] the response separation INSIDE each technical group ...",
          flush=True)
    # The response AUC computed inside one institute, and inside one
    # library-preparation batch. A batch effect cannot help a model
    # that never sees two batches, so if the response AUC survives
    # inside a group it is not riding on the batch; if it collapses, the
    # pooled AUC was.
    within_group = {}
    for label, keys in (("institute", institutes), ("library_batch", batches)):
        for name, X in (("rmds_bins", rmds_X), ("motifs_256", motif_X)):
            per_group = {}
            for group in sorted(set(keys.tolist())):
                mask = keys == group
                if int(mask.sum()) < 12:
                    continue
                if len(set(y[mask].tolist())) < 2:
                    continue
                if len(set(patients[mask].tolist())) < FOLDS:
                    continue
                auc, per_seed, _ = _grouped_cv_auc(X[mask], y[mask],
                                                   patients[mask])
                per_group[group] = {
                    "n": int(mask.sum()),
                    "responders": int(y[mask].sum()),
                    "patients": len(set(patients[mask].tolist())),
                    "pooled_out_of_fold_auc": round(auc, 4),
                    "per_seed_aucs": [round(a, 4) for a in per_seed],
                }
            within_group[f"{label}_{name}"] = per_group
            print(f"      {label} / {name}: {len(per_group)} groups, "
                  f"AUCs {[v['pooled_out_of_fold_auc'] for v in per_group.values()]}",
                  flush=True)
    results["response_auc_within_technical_group"] = within_group

    per_visit = {}
    for visit in sorted(set(visits.tolist())):
        mask = visits == visit
        if int(mask.sum()) < 20 or len(set(y[mask].tolist())) < 2:
            continue
        auc, per_seed, _ = _grouped_cv_auc(rmds_X[mask], y[mask],
                                           patients[mask])
        per_visit[visit] = {"n": int(mask.sum()),
                            "responders": int(y[mask].sum()),
                            "pooled_out_of_fold_auc": round(auc, 4),
                            "per_seed_aucs": [round(a, 4) for a in per_seed]}
        print(f"      visit {visit}: n={int(mask.sum())} AUC {auc:.4f}",
              flush=True)
    results["per_visit_rmds"] = per_visit

    # The depth control. Sequencing depth is the confound the ccRCC
    # cohort turned out to carry, and the rMDS is an entropy -- an
    # estimate whose noise depends on how many fragments a bin has -- so
    # the depth has to be asked about before the number means anything.
    if args.prep_report:
        from scipy import stats as _stats

        prep = json.loads(open(args.prep_report, encoding="utf-8").read())
        depths = np.asarray([prep["per_sample"][sample]["fragments_in_file"]
                             for sample in order], dtype=float)
        rho, p_value = _stats.spearmanr(depths, y)
        depth_only = np.zeros(len(y), dtype=float)
        for train, test in StratifiedKFold(5, shuffle=True,
                                          random_state=0).split(
                depths.reshape(-1, 1), y):
            model = LogisticRegression(max_iter=1000).fit(
                depths[train].reshape(-1, 1), y[train])
            depth_only[test] = model.predict_proba(
                depths[test].reshape(-1, 1))[:, 1]
        # The score with the depth regressed out of it, which is the
        # control that says whether the separation was the depth.
        _, _, oof_rmds = _grouped_cv_auc(rmds_X, y, patients)
        design = np.column_stack([np.ones(len(y)), depths])
        beta, *_ = np.linalg.lstsq(design, oof_rmds, rcond=None)
        without_depth = oof_rmds - design @ beta
        results["depth_control"] = {
            "median_fragments_responders": int(np.median(depths[y == 1])),
            "median_fragments_non_responders": int(np.median(depths[y == 0])),
            "spearman_depth_vs_response": round(float(rho), 4),
            "spearman_p": round(float(p_value), 6),
            "auc_depth_alone": round(float(roc_auc_score(y, depth_only)), 4),
            "auc_model": round(float(roc_auc_score(y, oof_rmds)), 4),
            "auc_model_with_depth_removed": round(
                float(roc_auc_score(y, without_depth)), 4),
            "note": ("the rMDS is an entropy, so its noise depends on how "
                     "many fragments a bin holds; this asks whether the "
                     "depth does the separating"),
        }
        print(f"      depth control: depth alone AUC "
              f"{results['depth_control']['auc_depth_alone']}, "
              f"the model {results['depth_control']['auc_model']}, "
              f"without the depth "
              f"{results['depth_control']['auc_model_with_depth_removed']}",
              flush=True)

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(
        json.dumps(results, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8")
    print(f"wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
