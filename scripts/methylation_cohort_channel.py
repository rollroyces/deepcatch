#!/usr/bin/env python3
"""Train one HMM across a cohort and decode per-sample beta values.

The per-sample beta profile is what a methylation *channel* needs: one
number set per sample that can be scored against the cancer / control
label. A single-sample HMM cannot provide it (it has no healthy
reference and its states collapse), so this script:

1. concatenates the per-sample per-CpG matrices produced by
   ``CpgMultiMetricsStatsV2`` into one cohort matrix;
2. attaches the WGBS methylation prior to every row (the column the
   decoder expects at index 11);
3. trains the HMM ONCE on the cohort -- the decoder trains when
   ``-decodeModeOnly`` is off -- and keeps the trained model;
4. decodes each sample separately with that trained model
   (``-decodeModeOnly``), which yields per-sample beta values;
5. summarises each sample's beta profile and scores the summary
   against the cancer / non-cancer label with leave-one-out CV.

Usage:
    python scripts/methylation_cohort_channel.py \
        --matrices-dir  /path/to/matrices \
        --labels        /path/to/cohort_labels.tsv \
        --out-json      results/methylation_cohort_channel.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO_ROOT)

FINALEME_HOME = os.environ.get(
    "FINALEME_HOME", "/Users/hermes/.hermes/.local/finaleme")
JAR = os.path.join(FINALEME_HOME, "FinaleMe-0.58-jar-with-dependencies.jar")
CLASS_PATH = (
    os.path.join(FINALEME_HOME, "fm_src/FinaleMe-main/target/classes")
    + ":" + os.path.join(FINALEME_HOME, "fm_src/FinaleMe-main/target/lib/*")
)
PRIOR = os.path.join(FINALEME_HOME, "wgbs_buffyCoat_jensen2015GB.methy.hg19.bw")
# Fraction of samples that must cover a CpG for it to enter the common
# set used by the composition control.
COVER_FRACTION = 0.8


def attach_prior(matrix_path: str, out_path: str) -> dict:
    """Concatenate one or more matrices and attach the WGBS prior.

    Keeps columns 0-10 (the decoder reads index 4, 7, 9, 10, 11) and
    replaces the trailing k-mer column with the prior, which is what the
    decoder parses at index 11.
    """
    import pyBigWig

    bw = pyBigWig.open(PRIOR)
    n_rows = n_prior = n_nan = 0
    with gzip.open(out_path, "wt", compresslevel=6) as fo:
        fo.write("\t".join(["chr", "start", "end", "readName", "FragLen",
                             "Frag_strand", "methy_stat", "Norm_Frag_cov",
                             "baseQ", "Offset_frag", "Dist_frag_end",
                             "methyPrior"]) + "\n")
        for matrix_path_i in matrix_path:
            # Read the whole matrix first: a file that is still being
            # written raises EOFError mid-iteration, which a per-line
            # guard would not catch.
            try:
                with gzip.open(matrix_path_i, "rt") as fi:
                    matrix_lines = fi.readlines()
            except (OSError, EOFError) as exc:
                print(f"skipping unreadable matrix {matrix_path_i}: {exc}",
                      file=sys.stderr)
                continue
            if True:
                for line in matrix_lines[1:]:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 11:
                        continue
                    prior = "NaN"
                    try:
                        vals = bw.values(parts[0], int(parts[1]), int(parts[2]))
                        if vals and vals[0] is not None:
                            val = float(vals[0])
                            if math.isnan(val):
                                n_nan += 1
                            else:
                                prior = f"{val:.4f}"
                                n_prior += 1
                    except Exception:
                        n_nan += 1
                    fo.write("\t".join(parts[:11] + [prior]) + "\n")
                    n_rows += 1
    bw.close()
    return {"rows": n_rows, "with_prior": n_prior, "nan_prior": n_nan}


def decode(model: str, matrix: str, out: str, decode_only: bool) -> int:
    cmd = ["java", "-Xmx8G", "-cp", CLASS_PATH,
           "org.cchmc.epifluidlab.finaleme.hmm.FinaleMe"]
    if decode_only:
        cmd.append("-decodeModeOnly")
    cmd += [model, matrix, out]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"decode failed for {matrix}:\n{proc.stderr[-1500:]}",
              file=sys.stderr)
    return proc.returncode



def load_decoded_cpg_betas(decoded_path: str) -> dict:
    """Map (chrom, start) -> beta for one decoded sample."""
    out = {}
    try:
        with gzip.open(decoded_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    out[(parts[0], int(parts[1]))] = float(parts[3]) / 100.0
                except ValueError:
                    continue
    except (OSError, EOFError):
        return {}
    return out


def summarise(decoded_path: str, max_cpgs: int | None = None) -> dict:
    """Summarise one sample's decoded beta values.

    ``max_cpgs`` depth-matches the sample: if it has more CpGs than the
    cap, every k-th CpG is kept, so the retained set still spans the
    genome. This is how the coverage confound is removed -- two samples
    with 10,000 and 50,000 covered CpGs are not comparable otherwise.
    """
    betas = []
    try:
        with gzip.open(decoded_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                try:
                    betas.append(float(parts[3]) / 100.0)
                except ValueError:
                    continue
    except (OSError, EOFError) as exc:
        # a truncated gz means the file is still being written
        print(f"skipping unreadable decode {decoded_path}: {exc}",
              file=sys.stderr)
        return {"n_cpgs": 0}
    a = np.asarray(betas, dtype=float)
    if a.size == 0:
        return {"n_cpgs": 0}
    n_total = int(a.size)
    if max_cpgs and a.size > max_cpgs:
        stride = -(-a.size // max_cpgs)
        a = a[::stride][:max_cpgs]
    return {
        "n_cpgs": int(a.size),
        "n_cpgs_before_match": n_total,
        "beta_mean": round(float(a.mean()), 4),
        "beta_median": round(float(np.median(a)), 4),
        "beta_std": round(float(a.std()), 4),
        "fraction_methylated": round(float((a > 0.5).mean()), 4),
        "fraction_fully_methylated": round(float((a >= 0.999).mean()), 4),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrices-dir", required=True,
                    help="directory of per-sample *.matrix.txt.gz files")
    ap.add_argument("--labels", required=True,
                    help="TSV: sample<TAB>label(0/1)<TAB>disease")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--out-json", default="results/methylation_cohort_channel.json")
    args = ap.parse_args()

    work = args.work_dir or os.path.join(os.path.dirname(args.out_json),
                                          "_methylation_cohort_work")
    os.makedirs(work, exist_ok=True)

    labels = {}
    with open(args.labels) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] != "sample":
                labels[parts[0]] = (int(parts[1]),
                                     parts[2] if len(parts) > 2 else "unknown")

    matrices = []
    for name in sorted(os.listdir(args.matrices_dir)):
        if name.endswith("matrix.txt.gz"):
            matrices.append(os.path.join(args.matrices_dir, name))
    if not matrices:
        print(f"no matrices in {args.matrices_dir}", file=sys.stderr)
        return 1
    print(f"[1/4] cohort matrix from {len(matrices)} samples", flush=True)
    cohort_matrix = os.path.join(work, "cohort.decoder_input.txt.gz")
    stats = attach_prior(matrices, cohort_matrix)
    print(f"      {stats}")

    model = os.path.join(work, "cohort.hmm_model")
    cohort_decoded = os.path.join(work, "cohort.decoded.txt.gz")
    print("[2/4] training the HMM on the cohort ...", flush=True)
    if decode(model, cohort_matrix, cohort_decoded, decode_only=False) != 0:
        return 1
    print(f"      trained model: {model}")

    print("[3/4] decoding each sample with the trained model ...", flush=True)
    per_sample = {}
    for m in matrices:
        name = os.path.basename(m).replace(".chr22_matrix.txt.gz", "")
        one = os.path.join(work, f"{name}.decoder_input.txt.gz")
        out = os.path.join(work, f"{name}.decoded.txt.gz")
        attach_prior([m], one)
        if decode(model, one, out, decode_only=True) != 0:
            continue
        s = summarise(out)
        if s.get("n_cpgs"):
            lab, dis = labels.get(name, (None, "unknown"))
            s["label"] = lab
            s["disease"] = dis
            per_sample[name] = s

    # Depth-matched pass: cap every sample at the smallest CpG count
    # seen, so the two arms are compared on the same number of loci.
    counts = [v["n_cpgs"] for v in per_sample.values()]
    cap = min(counts) if counts else 0
    matched = {}
    if cap >= 100:
        for name in list(per_sample):
            out = os.path.join(work, f"{name}.decoded.txt.gz")
            if not os.path.exists(out):
                continue
            m = summarise(out, max_cpgs=cap)
            m["label"] = per_sample[name].get("label")
            m["disease"] = per_sample[name].get("disease")
            matched[name] = m

    print("[4/4] scoring the beta summary against the labels ...", flush=True)
    names = sorted(per_sample)
    y = np.array([per_sample[n]["label"] if per_sample[n]["label"] is not None
                  else 0 for n in names])
    feats = np.array([[per_sample[n]["beta_mean"], per_sample[n]["beta_median"],
                        per_sample[n]["fraction_methylated"],
                        per_sample[n]["fraction_fully_methylated"]]
                       for n in names], dtype=float)
    n_cpgs = np.array([per_sample[n]["n_cpgs"] for n in names], dtype=float)

    def _loo_auc(X: np.ndarray) -> float | None:
        """Leave-one-out AUC of a small logistic model on X."""
        if len(names) < 5 or len(set(y.tolist())) < 2:
            return None
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import LeaveOneOut
        from sklearn.metrics import roc_auc_score
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

    auc = _loo_auc(feats)
    # Confound baseline: does the CpG-coverage count alone separate the
    # classes? If it does, the beta summary is not adding anything over
    # library complexity.
    auc_depth_only = _loo_auc(n_cpgs.reshape(-1, 1))
    # Single-feature AUC of beta_mean alone (no model, no shrinkage).
    try:
        from sklearn.metrics import roc_auc_score as _rauc
        auc_beta_mean_only = (round(float(_rauc(y, feats[:, 0])), 4)
                              if len(set(y.tolist())) > 1 else None)
    except Exception:
        auc_beta_mean_only = None

    auc_matched = None
    if matched:
        ym = np.array([matched[n]["label"] if matched[n]["label"] is not None
                        else 0 for n in sorted(matched)])
        fm = np.array([[matched[n]["beta_mean"], matched[n]["beta_median"],
                        matched[n]["fraction_methylated"],
                        matched[n]["fraction_fully_methylated"]]
                       for n in sorted(matched)], dtype=float)
        nm = np.array([matched[n]["n_cpgs"] for n in sorted(matched)],
                       dtype=float)
        auc_matched = _loo_auc(fm)
        auc_matched_depth = _loo_auc(nm.reshape(-1, 1))

    # Third control: restrict every sample to the SAME CpG set (the
    # intersection), which removes compositional differences as well as
    # depth. Depth matching caps the count but each sample still keeps a
    # different subset; a common set makes the comparison exact.
    common_summary = {}
    common_size = 0
    decoded_maps = {}
    for name in names:
        out = os.path.join(work, f"{name}.decoded.txt.gz")
        if os.path.exists(out):
            decoded_maps[name] = load_decoded_cpg_betas(out)
    if decoded_maps:
        # A strict intersection collapses on this cohort -- one sample
        # covers a different region set, so the intersection over all
        # samples can be empty even though pairs overlap in thousands of
        # CpGs. Use a coverage-frequency threshold instead: the CpGs that
        # at least COVERAGE_FRACTION of the samples see.
        from collections import Counter
        cover = Counter()
        for m in decoded_maps.values():
            cover.update(m.keys())
        n_samples = len(decoded_maps)
        need = max(2, int(np.ceil(COVER_FRACTION * n_samples)))
        common = {k for k, c in cover.items() if c >= need}
        if len(common) < 50:
            # Fall back to the most widely covered loci.
            common = {k for k, _ in cover.most_common(max(50, n_samples * 50))}
        common_size = len(common)
        if common_size >= 50:
            for name, m in decoded_maps.items():
                vals = [m[k] for k in common if k in m]
                if vals:
                    a = np.asarray(vals, dtype=float)
                    common_summary[name] = {
                        "n_cpgs": int(a.size),
                        "beta_mean": round(float(a.mean()), 4),
                        "beta_median": round(float(np.median(a)), 4),
                        "fraction_methylated": round(float((a > 0.5).mean()), 4),
                        "fraction_fully_methylated": round(
                            float((a >= 0.999).mean()), 4),
                        "label": per_sample[name].get("label"),
                        "disease": per_sample[name].get("disease"),
                    }
    auc_common = None
    if common_summary:
        cn = sorted(common_summary)
        yc = np.array([common_summary[n]["label"] if common_summary[n]["label"]
                        is not None else 0 for n in cn])
        fc = np.array([[common_summary[n]["beta_mean"],
                        common_summary[n]["beta_median"],
                        common_summary[n]["fraction_methylated"],
                        common_summary[n]["fraction_fully_methylated"]]
                       for n in cn], dtype=float)
        saved_names, saved_y = names, y

        names, y = cn, yc
        auc_common = _loo_auc(fc)
        names, y = saved_names, saved_y

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "cohort": {
            "n_samples_decoded": len(names),
            "n_samples_total": len(matrices),
            "per_sample": per_sample,
        },
        "cohort_matrix": stats,
        "beta_summary_auc_loo": auc,
        "beta_mean_only_auc": auc_beta_mean_only,
        "cpgs_covered_only_auc_loo": auc_depth_only,
        "common_cpg_set": {
            "n_shared_cpgs": int(common_size),
            "beta_summary_auc_loo": auc_common,
            "per_sample": common_summary,
        },
        "depth_matched": {
            "cpgs_cap": int(cap),
            "beta_summary_auc_loo": auc_matched,
            "cpgs_covered_only_auc_loo": auc_matched_depth,
            "per_sample": matched,
        },
        "confound_note": (
            "cpgs_covered_only_auc_loo is the baseline: if the number of "
            "CpGs with coverage alone separates the classes as well as the "
            "beta summary does, then the beta summary is not adding "
            "anything over library complexity, and the two must be "
            "compared before the channel is claimed."),
        "honest_framing": (
            "The HMM was trained once across the cohort and each sample was "
            "decoded separately, so every beta value is derived from that "
            "sample's own fragments with a shared model. The score above is "
            "a leave-one-out AUC on a four-number beta summary over a cohort "
            "of tens of samples: it is a feasibility measurement of the "
            "channel, not a validated feature. It will not survive a "
            "leave-one-study-out check, it uses no external cohort, and "
            "the GNN that would consume it is still untrained. It also "
            "has n=31 at most, so the leave-one-out AUC is unstable; read "
            "it next to the CpG-coverage baseline in the same artifact. "
            "CONFOUND TO CHECK BEFORE BELIEVING IT: the number of CpGs with "
            "coverage differs between samples (library complexity and depth), "
            "and the beta summary is not depth-normalised -- a difference in "
            "the summary can therefore reflect depth rather than biology. "
            "The per-sample n_cpgs column is the thing to look at first."),
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"wrote {args.out_json}")
    if auc is not None:
        print(f"  beta-summary AUC (LOO)      : {auc}")
    if auc_beta_mean_only is not None:
        print(f"  beta_mean alone (no model)  : {auc_beta_mean_only}")
    if auc_depth_only is not None:
        print(f"  CpGs-covered only (confound): {auc_depth_only}")
    if auc_common is not None:
        print(f"  common-CpG-set beta summary: {auc_common} "
              f"({common_size} shared CpGs)")
    if auc_matched is not None:
        print(f"  depth-matched beta summary   : {auc_matched} "
              f"(cap {cap} CpGs)")
        print(f"  depth-matched CpGs-only      : {auc_matched_depth}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
