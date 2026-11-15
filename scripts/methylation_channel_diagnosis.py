#!/usr/bin/env python3
"""Diagnose whether a cfDNA methylation channel carries signal or coverage.

Why this script exists
----------------------
``scripts/methylation_cohort_channel.py`` scores a methylome channel: it
trains FinaleMe's HMM on a cohort and reports a leave-one-out AUC for the
decoded beta summary. On a cohort whose fragment libraries do not carry
the digestion signature FinaleMe reads, that AUC can be produced entirely
by library complexity, and the AUC alone cannot tell the two apart.

This script is the diagnostic that makes the difference visible. It works
from the same per-sample matrices the decoder consumes and adds:

* **observation depth** -- how many fragment observations land on each CpG.
  A methylation state is only estimable where a CpG is observed
  repeatedly; below ~2 observations per CpG the state is not identified.
* **a model-free methylation estimator** -- the fraction of a CpG's
  observations whose fragment ends *at* that CpG. FinaleMe's feature file
  stores ``Dist_frag_end`` = distance from the CpG to the nearest
  fragment end. A CpG that lies exactly at a fragment end was cut out of
  the library, which is the digestion signature of an unmethylated CpG; a
  CpG that is interior to every fragment was not cut, which is the
  signature of a methylated one. So the boundary fraction is a direct
  estimate of the unmethylated probability that needs no HMM, no prior
  and no training, and it can be checked against the HMM output.
* **decoder degeneracy** -- how much of the HMM output is a hard 0 or 1
  rather than a graded value.
* **the coverage baseline next to every score**, so a channel result is
  never reported without the number it has to beat.

The script reports, for every feature it scores, the AUC that the number
of covered CpGs alone achieves. A channel feature that does not beat its
own coverage baseline is not evidence of methylation.

Usage
-----
    python scripts/methylation_channel_diagnosis.py \
        --matrices-dir  /path/to/matrices \
        --labels        /path/to/cohort_labels.tsv \
        --decoded-dir   /path/to/decoded      # optional: *.decoded.txt.gz \
        --out-json      results/methylation_channel_diagnosis.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from collections import defaultdict
from datetime import datetime, timezone

import numpy as np

# A CpG needs repeated observations before its methylation state is
# identified; below this the estimator is reported but flagged.
DEFAULT_MIN_OBS = 2


# --------------------------------------------------------------------------
# matrix parsing
# --------------------------------------------------------------------------
def read_matrix(path: str) -> dict:
    """Summarise one ``CpgMultiMetricsStatsV2`` matrix.

    Columns: 0 chr, 1 start, 2 end, 3 readName, 4 FragLen, 5 Frag_strand,
    6 methy_stat, 7 Norm_Frag_cov, 8 baseQ, 9 Offset_frag,
    10 Dist_frag_end, 11 kmer (replaced by the prior downstream).
    """
    n_obs = 0
    per_cpg: dict[tuple[str, int], list[int]] = defaultdict(list)
    methy_stat_counts: dict[str, int] = defaultdict(int)
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#") or line.startswith("chr\tstart"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 11:
                    continue
                try:
                    key = (parts[0], int(parts[1]))
                    dist_to_end = int(parts[10])
                except ValueError:
                    continue
                per_cpg[key].append(dist_to_end)
                methy_stat_counts[parts[6]] += 1
                n_obs += 1
    except (OSError, EOFError) as exc:
        print(f"skipping unreadable matrix {path}: {exc}", file=sys.stderr)
        return {"n_observations": 0, "n_cpgs": 0}

    n_boundary = sum(1 for v in per_cpg.values() for d in v if d == 0)
    obs_counts = np.array([len(v) for v in per_cpg.values()], dtype=float)
    return {
        "n_observations": int(n_obs),
        "n_cpgs": int(len(per_cpg)),
        "mean_obs_per_cpg": round(float(obs_counts.mean()), 3) if obs_counts.size else 0.0,
        "median_obs_per_cpg": int(np.median(obs_counts)) if obs_counts.size else 0,
        "max_obs_per_cpg": int(obs_counts.max()) if obs_counts.size else 0,
        "n_observations_at_fragment_end": int(n_boundary),
        "boundary_observation_rate": (round(n_boundary / n_obs, 5)
                                      if n_obs else None),
        "methy_stat_counts": dict(methy_stat_counts),
        "_per_cpg": per_cpg,
    }


def boundary_estimates(per_cpg: dict, min_obs: int) -> dict:
    """Per-sample methylation estimate from fragment-boundary evidence.

    For every CpG the fraction of its observations that sit at a fragment
    end estimates the probability that the CpG was unmethylated (it was
    cut out of the library). CpGs seen fewer than ``min_obs`` times are
    excluded: a single observation cannot distinguish "interior" from
    "at the end".
    """
    vals = []
    n_used = 0
    for dists in per_cpg.values():
        if len(dists) < min_obs:
            continue
        n_used += 1
        vals.append(sum(1 for d in dists if d == 0) / len(dists))
    if not vals:
        return {"n_cpgs_used": 0, "mean_unmethylated_fraction": None,
                "median_unmethylated_fraction": None,
                "std_unmethylated_fraction": None,
                "fraction_cpgs_fully_unmethylated": None,
                "fraction_cpgs_fully_methylated": None}
    a = np.asarray(vals, dtype=float)
    return {
        "n_cpgs_used": int(a.size),
        "mean_unmethylated_fraction": round(float(a.mean()), 5),
        "median_unmethylated_fraction": round(float(np.median(a)), 5),
        "std_unmethylated_fraction": round(float(a.std()), 5),
        "fraction_cpgs_fully_unmethylated": round(float((a >= 0.999).mean()), 5),
        "fraction_cpgs_fully_methylated": round(float((a <= 0.001).mean()), 5),
    }


def read_decoded(path: str) -> dict:
    """Per-CpG decoded values plus a degeneracy profile."""
    betas: dict[tuple[str, int], float] = {}
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    betas[(parts[0], int(parts[1]))] = float(parts[3]) / 100.0
                except ValueError:
                    continue
    except (OSError, EOFError):
        return {"n_cpgs": 0, "_betas": {}}
    if not betas:
        return {"n_cpgs": 0, "_betas": {}}
    a = np.asarray(list(betas.values()), dtype=float)
    return {
        "n_cpgs": int(a.size),
        "beta_mean": round(float(a.mean()), 5),
        "beta_std": round(float(a.std()), 5),
        "fraction_exactly_1": round(float((a >= 0.999).mean()), 5),
        "fraction_exactly_0": round(float((a <= 0.001).mean()), 5),
        "fraction_intermediate": round(
            float(((a > 0.001) & (a < 0.999)).mean()), 5),
        "_betas": betas,
    }


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def loo_auc(names: list[str], y: np.ndarray, X: np.ndarray) -> float | None:
    """Leave-one-out AUC of a logistic model on X."""
    if len(names) < 5 or len(set(y.tolist())) < 2:
        return None
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import LeaveOneOut
        from sklearn.metrics import roc_auc_score
    except ImportError:
        return None
    oof = np.zeros(len(names))
    for tr, te in LeaveOneOut().split(X):
        if len(set(y[tr].tolist())) < 2:
            continue
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        clf = LogisticRegression(max_iter=1000).fit((X[tr] - mu) / sd, y[tr])
        oof[te] = clf.predict_proba((X[te] - mu) / sd)[:, 1]
    if not np.any(oof):
        return None
    return round(float(roc_auc_score(y, oof)), 4)


def simple_auc(y: np.ndarray, x: np.ndarray) -> float | None:
    """Rank AUC of a raw feature, in the feature's own direction.

    A rank AUC is direction-sensitive: a feature that is *lower* in the
    positive class scores below 0.5 even though it separates perfectly.
    Every score is therefore reported with ``auc_directionless`` as well,
    so a feature that separates in the opposite direction is not mistaken
    for a feature that does not separate at all.
    """
    if len(set(y.tolist())) < 2:
        return None
    try:
        from sklearn.metrics import roc_auc_score
    except ImportError:
        return None
    return round(float(roc_auc_score(y, x)), 4)


def separation(x: np.ndarray, y: np.ndarray) -> dict:
    """Direction-agnostic separation of one feature, with its direction."""
    auc = simple_auc(y, x)
    if auc is None:
        return {"auc": None, "auc_directionless": None, "direction": None}
    dl = max(auc, 1.0 - auc)
    direction = ("higher_feature_means_positive_class" if auc >= 0.5
                 else "lower_feature_means_positive_class")
    return {"auc": auc, "auc_directionless": round(dl, 4),
            "direction": direction}


def matched_locus_test(per_cpg_values: dict, labels: dict,
                       n_permutations: int = 200) -> dict:
    """Per-CpG cancer-minus-control difference, with a permutation null.

    Every CpG is tested on the samples that actually cover it, because the
    two arms do not cover identical loci. The reported statistic is the
    MEDIAN of the per-CpG differences; the null comes from shuffling the
    sample labels.
    """
    vals_by_cpg: dict[tuple[str, int], dict[str, float]] = defaultdict(dict)
    for name, per_cpg in per_cpg_values.items():
        for key, v in per_cpg.items():
            vals_by_cpg[key][name] = v
    rows = []
    for key, d in vals_by_cpg.items():
        ys = [labels.get(n, 0) for n in d]
        if len(set(ys)) < 2:
            continue
        pos = [d[n] for n in d if labels.get(n, 0) == 1]
        neg = [d[n] for n in d if labels.get(n, 0) == 0]
        if not pos or not neg:
            continue
        rows.append(float(np.mean(pos)) - float(np.mean(neg)))
    if len(rows) < 50:
        return {"n_cpgs_tested": len(rows),
                "note": "too few CpGs covered by both arms"}
    diffs = np.asarray(rows, dtype=float)
    obs = float(np.median(diffs))
    rng = np.random.default_rng(0)
    names = sorted(per_cpg_values)
    null = []
    for _ in range(n_permutations):
        perm = {n: int(rng.integers(0, 2)) for n in names}
        d2 = []
        for key, d in vals_by_cpg.items():
            ys = [perm.get(n, 0) for n in d]
            if len(set(ys)) < 2:
                continue
            pos = [d[n] for n in d if perm.get(n, 0) == 1]
            neg = [d[n] for n in d if perm.get(n, 0) == 0]
            if pos and neg:
                d2.append(float(np.mean(pos)) - float(np.mean(neg)))
        if d2:
            null.append(float(np.median(d2)))
    null = np.asarray(null)
    return {
        "n_cpgs_tested": int(len(rows)),
        "median_difference_cancer_minus_control": round(obs, 6),
        "fraction_cpgs_higher_in_cancer": round(float((diffs > 0).mean()), 5),
        "permutation_p_value": (round(float((np.abs(null) >= abs(obs)).mean()), 5)
                                if null.size else None),
        "null_median_sd": (round(float(null.std()), 8) if null.size else None),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrices-dir", required=True,
                    help="directory of per-sample *.matrix.txt.gz files")
    ap.add_argument("--labels", required=True,
                    help="TSV: sample<TAB>label(0/1)<TAB>disease")
    ap.add_argument("--decoded-dir", default=None,
                    help="optional directory of *.decoded.txt.gz to compare "
                         "the HMM output with the model-free estimator")
    ap.add_argument("--min-obs", type=int, default=DEFAULT_MIN_OBS,
                    help="observations a CpG needs before its methylation "
                         f"state is used (default {DEFAULT_MIN_OBS})")
    ap.add_argument("--out-json", default="results/methylation_channel_diagnosis.json")
    args = ap.parse_args()

    labels: dict[str, int] = {}
    diseases: dict[str, str] = {}
    with open(args.labels) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] != "sample":
                labels[parts[0]] = int(parts[1])
                diseases[parts[0]] = parts[2] if len(parts) > 2 else "unknown"

    matrices = sorted(os.path.join(args.matrices_dir, n)
                      for n in os.listdir(args.matrices_dir)
                      if n.endswith("matrix.txt.gz"))
    if not matrices:
        print(f"no matrices in {args.matrices_dir}", file=sys.stderr)
        return 1

    per_sample: dict[str, dict] = {}
    decoded_values: dict[str, dict] = {}
    per_cpg_boundary: dict[str, dict] = {}
    for path in matrices:
        name = os.path.basename(path).replace(".chr22_matrix.txt.gz", "")
        info = read_matrix(path)
        if not info["n_observations"]:
            continue
        per_cpg = info.pop("_per_cpg")
        est = boundary_estimates(per_cpg, args.min_obs)
        info.update(est)
        # The same estimator over every covered CpG, including the ones
        # seen only once. At a library depth where almost no CpG is seen
        # twice this is the only variant with enough loci to score, so it
        # is reported next to the identified-state version.
        info["all_covered_cpgs"] = boundary_estimates(per_cpg, 1)
        info["fraction_cpgs_identified"] = round(
            info["n_cpgs_used"] / info["n_cpgs"], 5) if info["n_cpgs"] else None
        info["label"] = labels.get(name)
        info["disease"] = diseases.get(name, "unknown")
        per_sample[name] = info
        # Keep the per-CpG boundary estimates for the matched-locus test.
        per_cpg_boundary[name] = {
            key: sum(1 for d in dists if d == 0) / len(dists)
            for key, dists in per_cpg.items() if len(dists) >= args.min_obs
        }

        if args.decoded_dir:
            cand = os.path.join(args.decoded_dir, f"{name}.decoded.txt.gz")
            if os.path.exists(cand):
                dec = read_decoded(cand)
                if dec.get("n_cpgs"):
                    info["decoded"] = {k: v for k, v in dec.items()
                                       if not k.startswith("_")}
                    decoded_values[name] = dec["_betas"]

    if not per_sample:
        print("no readable matrices", file=sys.stderr)
        return 1

    names = sorted(per_sample)
    y = np.array([per_sample[n]["label"] if per_sample[n]["label"] is not None
                  else 0 for n in names], dtype=int)

    def col(key: str) -> np.ndarray:
        return np.array([per_sample[n].get(key) if per_sample[n].get(key)
                         is not None else np.nan for n in names], dtype=float)

    n_cpgs = col("n_cpgs")
    boundary = col("mean_unmethylated_fraction")
    depth = col("mean_obs_per_cpg")
    scores: dict[str, object] = {}

    def safe_sep(x: np.ndarray) -> dict:
        if np.isnan(x).any() or not np.isfinite(x).any():
            return {"auc": None, "auc_directionless": None, "direction": None}
        return separation(x, y)

    scores["coverage_baseline_n_cpgs"] = safe_sep(n_cpgs)
    scores["coverage_baseline_obs_per_cpg"] = safe_sep(depth)
    scores["boundary_estimator"] = safe_sep(boundary)

    if decoded_values and len(decoded_values) >= 5:
        dnames = sorted(decoded_values)
        dy = np.array([per_sample[n]["label"] if per_sample[n]["label"] is not None
                       else 0 for n in dnames], dtype=int)
        dbeta = np.array([np.mean(list(decoded_values[n].values()))
                          for n in dnames], dtype=float)
        scores["hmm_beta_mean"] = separation(dbeta, dy)
        inter = np.array([per_sample[n]["decoded"]["fraction_intermediate"]
                          for n in dnames], dtype=float)
        scores["hmm_fraction_intermediate_mean"] = round(float(inter.mean()), 5)
        scores["hmm_fraction_intermediate_min"] = round(float(inter.min()), 5)
        scores["hmm_fraction_intermediate_max"] = round(float(inter.max()), 5)

    # Matched-locus comparison: the same statistic computed only on CpGs
    # both arms cover, where composition cannot manufacture a difference.
    matched = matched_locus_test(per_cpg_boundary, labels)
    matched_hmm = (matched_locus_test(decoded_values, labels)
                   if decoded_values else None)

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cohort": {
            "n_samples": len(names),
            "n_samples_with_decoded_output": len(decoded_values),
            "min_observations_per_cpg": args.min_obs,
        },
        "depth": {
            "mean_obs_per_cpg_across_samples": (round(float(np.nanmean(depth)), 3)
                                                if np.isfinite(depth).any() else None),
            "min_obs_per_cpg_across_samples": (round(float(np.nanmin(depth)), 3)
                                              if np.isfinite(depth).any() else None),
            "max_obs_per_cpg_across_samples": (round(float(np.nanmax(depth)), 3)
                                              if np.isfinite(depth).any() else None),
            "note": "a methylation state is only identified where a CpG is "
                    "observed repeatedly; below ~2 observations per CpG the "
                    "per-CpG state carries no information",
        },
        "scores": scores,
        "matched_locus_boundary_estimator": matched,
        "matched_locus_hmm_beta": matched_hmm,
        "per_sample": per_sample,
        "interpretation": {
            "coverage_baseline": "coverage_baseline_n_cpgs_auc is the number "
                                 "of covered CpGs alone; a methylation feature "
                                 "must beat it to mean anything",
            "boundary_estimator": "mean_unmethylated_fraction is the share of "
                                   "a CpG's observations that sit at a fragment "
                                   "end; it estimates the unmethylated state "
                                   "without a model, a prior or training",
            "matched_locus": "the matched-locus tests restrict both arms to "
                             "the CpGs they share, so a difference there "
                             "cannot be produced by library composition",
            "decoder_degeneracy": "fraction_intermediate near zero means the "
                                  "decoder returned hard 0/1 calls, so the "
                                  "beta summary is a count of called states "
                                  "rather than a graded methylation level",
        },
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps({k: v for k, v in payload.items()
                      if k not in ("per_sample",)}, indent=2))
    print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
