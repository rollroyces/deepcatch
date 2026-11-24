#!/usr/bin/env python3
"""GSE124600's plasma methylation matrix, measured: is it a methylome or a call matrix?

Why this exists
---------------
`docs/PUBLIC_DATA_SOURCES.md` §4 records GSE124600 as the one public
plasma cohort with per-sample stage labels (142 cancer, 132 healthy,
`crc grade: Stage_I/II/III/IV`) and per-sample methylation matrices.
`ULTRA_EARLY_READINESS.md` §1.3 calls the stage breakdown the one
readiness gap, and says that closing it would need "a methylation
feature extractor and its own model".

The matrix is the data, so it was measured rather than assumed. Two
questions, both answered by the deposited files:

1. **What is in the matrix?** The values are 0 for 91.1% of the
   205,250 x 274 cells and reach 6,949.25 elsewhere. That is a
   methylated-CpG CALL matrix with normalised counts, not a fractional
   methylome: the zero means "this CpG was not called methylated in this
   sample", not "this CpG is 0% methylated".
2. **Is anything separable, and is any of it methylation?** The
   per-locus summaries are ANTI-correlated with cancer (AUC 0.27-0.41),
   and the number of CpGs called methylated per sample separates the
   classes at AUC 0.7688. That is a burden/coverage signal, which is
   exactly what this repository's own methylation channel found on
   FinaleMe -- where the decode emits hard 0/1 state assignments and the
   "methylation" the channel consumed was the fraction of covered loci
   called methylated.

The stage breakdown of that coverage signal is reported per stage and
tested for monotonicity, because the readiness document's question is
whether a stage-labelled public cohort delivers a stage-stratified
result. It does not: the called-CpG count does not rise with stage
(Spearman rho = 0.096, p = 0.26 over the 142 cancer samples), and the
per-stage AUCs are noisy at these arm sizes.

What this does not do
---------------------
It does not claim a methylation result, and it does not replace the
repository's model with a new one. It measures what the cohort contains
and what can be read from it, so the readiness document's gap is closed
with a number instead of a plan.

Usage:
    # 1. fetch the deposited matrices (23.6 MB, public, no auth)
    mkdir -p /tmp/gse124600
    cd /tmp/gse124600
    base=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124600/suppl
    for f in GSE124600_sumP_20181112CRC_uploadGEO_nodup1.5_6_12_C_umepm.txt.gz \\
             GSE124600_sumP_20181112CRC_uploadGEO_dup1.5_6_12_C_umepm.txt.gz; do
      curl -sL "$base/$f" -o "$f"
    done
    curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124600&targ=gsm&form=text&view=brief" \
         -o gse124600_gsm.txt

    # 2. measure it
    python scripts/geo_gse124600_methylation_cohort.py \\
        --matrix-dir /tmp/gse124600 \\
        --geo-records /tmp/gse124600/gse124600_gsm.txt \\
        --out-json results/methylation_cohort_geo_gse124600.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
from datetime import datetime, timezone

import numpy as np

PLASMA_MATRICES = (
    "GSE124600_sumP_20181112CRC_uploadGEO_nodup1.5_6_12_C_umepm.txt.gz",
    "GSE124600_sumP_20181112CRC_uploadGEO_dup1.5_6_12_C_umepm.txt.gz",
)

STAGE_CODES = {"Stage_I": 1, "Stage_II": 2, "Stage_III": 3, "Stage_IV": 4}


def read_matrix(path: str) -> tuple[list[str], np.ndarray]:
    """The matrix: rows are loci, columns are samples.

    The plasma headers carry their sample block space-separated inside a
    single tab field, while the data rows are tab-separated -- so the
    header has to be split on whitespace and the data on tabs.
    """
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        names: list[str] = []
        for field in header[7:]:
            names.extend(field.split())
        rows = []
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rows.append([float(x) for x in parts[7:]])
    return names, np.asarray(rows, dtype=np.float32)


def read_geo_records(path: str) -> dict[str, tuple[str, str]]:
    """sample title -> (disease state, crc grade) from the GEO records."""
    records: dict[str, tuple[str, str]] = {}
    text = open(path, encoding="utf-8", errors="replace").read()
    for block in text.split("^SAMPLE = ")[1:]:
        title = re.search(r"!Sample_title = (.+)", block)
        if not title:
            continue
        disease = re.search(
            r"!Sample_characteristics_ch1 = disease state: (.+)", block)
        stage = re.search(r"!Sample_characteristics_ch1 = crc grade: (\S+)", block)
        records[title.group(1).strip()] = (
            disease.group(1).strip() if disease else "",
            stage.group(1).strip() if stage else "",
        )
    return records


def match_sample(name: str, records: dict[str, tuple[str, str]]) -> tuple[str, str]:
    """The matrix's column names carry the depositor's suffixes."""
    candidates = [name]
    for suffix in ("_dup1.5", "_m"):
        if candidates[-1].endswith(suffix):
            candidates.append(candidates[-1][: -len(suffix)])
    for candidate in candidates:
        if candidate in records:
            return records[candidate]
    return ("", "")


def auc(y_true, y_score) -> float:
    """The rank AUC (ties count half)."""
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score, dtype=float)
    pos, neg = y_score[y_true == 1], y_score[y_true == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    order = np.argsort(y_score, kind="mergesort")
    ranks = np.empty(len(y_score), dtype=float)
    sorted_scores = y_score[order]
    index = 0
    while index < len(sorted_scores):
        last = index
        while last + 1 < len(sorted_scores) and sorted_scores[last + 1] == sorted_scores[index]:
            last += 1
        ranks[order[index:last + 1]] = (index + last) / 2.0 + 1.0
        index = last + 1
    return float((ranks[y_true == 1].sum() - len(pos) * (len(pos) + 1) / 2.0)
                 / (len(pos) * len(neg)))


def permutation_p(y_true, y_score, trials: int = 20000, seed: int = 20261124):
    """A label-permutation p-value for the rank statistic."""
    rng = np.random.default_rng(seed)
    observed = auc(y_true, y_score)
    index = np.arange(len(y_true))
    hits = 0
    for _ in range(trials):
        rng.shuffle(index)
        if auc(y_true[index], y_score) >= observed - 1e-12:
            hits += 1
    return round(observed, 4), round((hits + 1) / (trials + 1), 6)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix-dir", required=True)
    ap.add_argument("--geo-records", required=True)
    ap.add_argument("--out-json",
                    default="results/methylation_cohort_geo_gse124600.json")
    args = ap.parse_args()

    names: list[str] = []
    blocks: list[np.ndarray] = []
    for name in PLASMA_MATRICES:
        sample_names, matrix = read_matrix(os.path.join(args.matrix_dir, name))
        names += sample_names
        blocks.append(matrix)
    matrix = np.concatenate(blocks, axis=1)
    print(f"matrix: {matrix.shape[0]} loci x {matrix.shape[1]} samples")

    records = read_geo_records(args.geo_records)
    labels, stages = [], []
    unmatched = 0
    for name in names:
        disease, stage = match_sample(name, records)
        if not disease:
            unmatched += 1
        labels.append(1 if "cancer" in disease.lower() else 0)
        stages.append(stage)
    labels = np.asarray(labels)
    stages = np.asarray(stages)
    if unmatched:
        print(f"  WARNING: {unmatched} columns did not match a GEO record")

    zero_fraction = float((matrix == 0).mean())
    nonzero = (matrix > 0).sum(axis=0).astype(float)
    sums = matrix.sum(axis=0).astype(float)
    coverage = (matrix > 0).sum(axis=1)

    summaries = {}
    for fraction in (0.5, 0.9):
        keep = coverage >= fraction * matrix.shape[1]
        shared_mean = matrix[keep].mean(axis=0)
        mean_nonzero = np.array([
            matrix[:, j][matrix[:, j] > 0].mean()
            if (matrix[:, j] > 0).any() else np.nan
            for j in range(matrix.shape[1])])
        block = {"n_shared_loci": int(keep.sum())}
        for label, summary in (
            ("mean_over_shared_loci", shared_mean),
            ("mean_of_nonzero_values", mean_nonzero),
            ("called_loci_count", nonzero),
            ("total_signal", sums),
        ):
            ok = ~np.isnan(summary)
            observed, p = permutation_p(labels[ok], summary[ok])
            block[label] = {"auc": observed, "permutation_p": p}
        summaries[f"shared_loci_at_{int(fraction * 100)}pct"] = block

    # The stage breakdown of the only class signal, and whether it is
    # monotone in stage.
    healthy = nonzero[labels == 0]
    per_stage = {}
    for stage in sorted(set(stages[labels == 1].tolist())):
        mask = (labels == 1) & (stages == stage)
        if not mask.any():
            continue
        y = np.concatenate([np.ones(int(mask.sum())), np.zeros(len(healthy))])
        s = np.concatenate([nonzero[mask], healthy])
        per_stage[stage] = {
            "n": int(mask.sum()),
            "median_called_loci": float(np.median(nonzero[mask])),
            "auc_vs_healthy": round(auc(y, s), 4),
        }

    codes = np.array([STAGE_CODES.get(s, 0) for s in stages])
    cancer = labels == 1
    rho, rho_p = None, None
    if len(np.unique(codes[cancer])) > 1:
        from scipy.stats import spearmanr
        result = spearmanr(codes[cancer], nonzero[cancer])
        rho, rho_p = round(float(result.statistic), 4), round(float(result.pvalue), 4)

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "accession": "GSE124600",
        "question": "What does the deposited GSE124600 plasma methylation matrix "
                    "contain, and what can be read from it -- including the "
                    "stage breakdown this repository's readiness document was "
                    "waiting for?",
        "cohort": {
            "n_samples": int(matrix.shape[1]),
            "n_cancer": int((labels == 1).sum()),
            "n_healthy": int((labels == 0).sum()),
            "n_loci": int(matrix.shape[0]),
            "stage_distribution": {
                stage: int(((labels == 1) & (stages == stage)).sum())
                for stage in sorted(set(stages[labels == 1].tolist()))
            },
            "median_called_loci_healthy": float(np.median(nonzero[labels == 0])),
            "median_called_loci_cancer": float(np.median(nonzero[labels == 1])),
        },
        "matrix_semantics": {
            "fraction_of_cells_exactly_zero": round(zero_fraction, 4),
            "max_value": float(matrix.max()),
            "min_nonzero_value": float(matrix[matrix > 0].min()),
            "reading": "A methylated-CpG CALL matrix with normalised counts, "
                       "not a fractional methylome: 91% of the cells are 0, "
                       "which means 'this CpG was not called methylated in "
                       "this sample', not 'this CpG is 0% methylated'. The "
                       "per-locus summaries are anti-correlated with cancer "
                       "because the signal is spread over more loci, which "
                       "is what a normalised count matrix does.",
        },
        "summaries": summaries,
        "per_stage_called_loci": per_stage,
        "stage_monotonicity": {
            "spearman_rho_stage_vs_called_loci": rho,
            "p": rho_p,
            "n_cancer": int(cancer.sum()),
            "reading": "The called-CpG count does NOT rise with stage, so the "
                       "stage-labelled cohort does not deliver a "
                       "stage-stratified early-detection result. The "
                       "per-stage AUCs are noisy at these arm sizes.",
        },
        "conclusion": "The cohort reproduces this repository's methylation "
                      "finding on an independent, larger, stage-labelled "
                      "deposit: the matrix yields a coverage/burden signal, "
                      "not a methylation-level signal. It therefore closes "
                      "the readiness gap with a measurement rather than a "
                      "plan, and the negative answer is the answer.",
    }
    os.makedirs(os.path.dirname(args.out_json) or ".", exist_ok=True)
    with open(args.out_json, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    print(f"  zero cells: {zero_fraction:.3f}; max value {matrix.max():.2f}")
    print(f"  called loci: healthy median {np.median(nonzero[labels == 0]):,.0f}, "
          f"cancer median {np.median(nonzero[labels == 1]):,.0f}")
    for key, block in summaries.items():
        print(f"  {key}: " + "  ".join(
            f"{name}={value['auc']:.4f}" for name, value in block.items()
            if isinstance(value, dict) and "auc" in value))
    for stage, block in per_stage.items():
        print(f"  {stage}: n={block['n']} AUC vs healthy={block['auc_vs_healthy']:.4f}")
    print(f"  Spearman(stage, called loci) = {rho} (p={rho_p})")
    print(f"\n  wrote {args.out_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
