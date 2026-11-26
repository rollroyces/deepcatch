#!/usr/bin/env python3
"""Turn the end-motif deposit's tables into the repository's motif features.

What this is
------------
Zenodo 10.5281/zenodo.17253504 ("cirDNA End-Motif Frequencies (2- and
4-mers) from sWGS") publishes, per sample, the counts of each 4-mer at
the 5' end of cfDNA fragments, for healthy controls and for breast
cancer, NSCLC and metastatic colorectal cancer, in two library
preparations (SSP and DSP). The tables are spreadsheets: one sheet per
group and assay, the 256 motifs as rows and the samples as columns.

Why it is worth converting
--------------------------
This is the only public deposit the survey found that carries **the
channel the PHTS paper's own finding rests on** -- fragment end
motifs -- for a cohort with cancer/healthy labels and per-sample
values, and it is 2.3 MB rather than gigabytes. The repository's
motif channel has so far only been measured inside the training
cohort.

The layout, which the parser has to get right
----------------------------------------------
Each sheet is a 256-row table of motifs (AAAA ... TTTT, the
lexicographic ACGT order -- the same order as the repository's
``_kmer_index``, which is checked here) with the samples in two column
blocks: the counts on the left, the fractions on the right, separated
by an empty column and with a trailing "Mean" column. The two blocks
must describe the same samples, and the parser verifies that rather
than assuming it.

What it writes
--------------
``<sample>.motifs.npy`` per sample, in the format the repository's own
loader expects (the 256 frequencies in the repository's motif order),
plus a label table with the group, the assay and the class, and a
report of what was found.

Usage:
    python3 scripts/prepare_zenodo_endmotif_cohort.py \
        --dir /tmp/endmotif \
        --out-dir data/external_cohort_zenodo_endmotif \
        --report results/zenodo_endmotif_prep.json
"""
from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone

import numpy as np

ENDMOTIF_DOI = "10.5281/zenodo.17253504"

# The groups the deposit's file names carry, and what they mean.
#   HI  -- healthy controls
#   BC  -- breast cancer
#   NSCLC -- non-small-cell lung cancer
#   mCRC -- metastatic colorectal cancer
_GROUPS = {
    "HI25": ("healthy", "healthy control"),
    "HI37": ("healthy", "healthy control"),
    "BC": ("breast_cancer", "breast cancer"),
    "NSCLC": ("nsclc", "non-small-cell lung cancer"),
    "mCRC": ("mcr_colorectal", "metastatic colorectal cancer"),
}


def _kmer_index(kmer: str) -> int:
    """The repository's motif index (A=0, C=1, G=2, T=3, lexicographic)."""
    m = {"A": 0, "C": 1, "G": 2, "T": 3}
    if any(c not in m for c in kmer):
        return -1
    return sum(m[c] << (2 * (3 - i)) for i, c in enumerate(kmer))


def _read_sheet(path: str) -> tuple[list[str], dict, list[str]]:
    """Read one sheet into (sample ids, counts per sample, motif order)."""
    import openpyxl
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    sheet = workbook[workbook.sheetnames[0]]
    rows = list(sheet.iter_rows(values_only=True))
    header = list(rows[0])
    body = rows[1:]
    # The last row is a SUM formula; the motifs are the rows above it.
    body = [r for r in body if isinstance(r[0], str) and len(r[0]) == 4]
    motifs = [str(r[0]) for r in body]

    separator = [i for i, h in enumerate(header) if h in (None, "")]
    if not separator:
        raise SystemExit(f"{path}: no separator column between the two blocks")
    cut = separator[0]
    left = [(i, header[i]) for i in range(1, cut) if header[i] not in (None, "")]
    right = [(i, header[i]) for i in range(cut + 1, len(header))
             if header[i] not in (None, "", "Mean")]

    left_ids = [str(name) for _i, name in left]
    right_ids = [str(name) for _i, name in right]
    if left_ids != right_ids:
        raise SystemExit(
            f"{path}: the two blocks describe different samples "
            f"({len(left_ids)} vs {len(right_ids)})")

    counts: dict = {}
    for _i, name in left:
        column = header.index(name, 1)
        values = np.array([float(r[column]) for r in body], dtype=float)
        total = values.sum()
        counts[str(name)] = values / total if total > 0 else values
    return left_ids, counts, motifs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", required=True,
                    help="directory holding the deposit's .xlsx tables")
    ap.add_argument("--out-dir", default="data/external_cohort_zenodo_endmotif")
    ap.add_argument("--report", default="results/zenodo_endmotif_prep.json")
    args = ap.parse_args()

    files = sorted(name for name in os.listdir(args.dir)
                   if name.endswith("_4bp_5pr.xlsx"))
    if not files:
        raise SystemExit(f"no 4-mer tables in {args.dir}")

    features_dir = os.path.join(args.out_dir, "features")
    os.makedirs(features_dir, exist_ok=True)

    per_sample: dict = {}
    groups: dict = {}
    problems: list[str] = []
    for name in files:
        stem = name[:-len(".xlsx")]
        group, _assay = stem.split("_", 1)
        assay = stem.split("_")[1]
        if group not in _GROUPS:
            problems.append(f"{name}: unknown group {group}")
            continue
        label, description = _GROUPS[group]
        sample_ids, counts, motifs = _read_sheet(
            os.path.join(args.dir, name))

        # The motif order has to be the repository's own.
        expected = [m for i in range(256)
                    for m in [_index_to_kmer(i)] if m]
        if motifs != expected:
            problems.append(
                f"{name}: the motif rows are not in the repository's order "
                f"(first difference at "
                f"{next((i for i, (a, b) in enumerate(zip(motifs, expected)) if a != b), 0)})")
            continue

        for sample in sample_ids:
            key = f"{sample}__{group}__{assay}"
            np.save(os.path.join(features_dir, f"{key}.motifs.npy"),
                    counts[sample].astype(np.float64))
            per_sample[key] = {
                "sample": sample,
                "group": group,
                "group_label": label,
                "group_description": description,
                "assay": assay,
                "n_motifs": int(counts[sample].size),
                "total_count": int(round(float(counts[sample].sum()))),
            }
        groups.setdefault(group, {"label": label, "assays": [], "n_samples": 0})
        groups[group]["assays"].append(assay)
        groups[group]["n_samples"] += len(sample_ids)
        print(f"      {stem}: {len(sample_ids)} samples ({label})", flush=True)

    if problems:
        print("[FAILED] the cohort is not ready:")
        for problem in problems:
            print(f"      - {problem}")
        return 1

    labels_path = os.path.join(args.out_dir, "labels.tsv")
    with open(labels_path, "w", encoding="utf-8") as handle:
        handle.write("sample_id\tlabel\tgroup\tassay\tdescription\n")
        for key, row in sorted(per_sample.items()):
            handle.write(f"{key}\t{row['group_label']}\t{row['group']}\t"
                         f"{row['assay']}\t{row['group_description']}\n")

    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {
            "record": "https://zenodo.org/records/17253504",
            "doi": ENDMOTIF_DOI,
            "paper": "Pisareva & Richaud -- cirDNA End-Motif Frequencies "
                     "(2- and 4-mers) from sWGS; preprint "
                     "10.1101/2025.09.17.676909",
        },
        "format": {
            "table": "one sheet per group and assay; 256 motif rows, the "
                     "samples in two column blocks (counts on the left, "
                     "fractions on the right) separated by an empty column, "
                     "with a trailing Mean column",
            "motif_order": "lexicographic ACGT -- the repository's own "
                           "_kmer_index order, verified against it here",
            "features": "<sample>.motifs.npy -- the 256 frequencies in the "
                        "repository's motif order, the format the "
                        "repository's own loader reads",
            "library": "shallow WGS, primary alignments, proper pairs, "
                       "MAPQ >= 30, GRCh38; motifs at the 5' fragment "
                       "boundary",
            "note": ("The deposit states that the motifs come from the 5' "
                     "boundary but does not say whether both strands' ends "
                     "were counted per fragment, which is what the "
                     "repository's own motif extraction does. That "
                     "ambiguity is carried with every number computed "
                     "from these features."),
        },
        "groups": groups,
        "samples": {
            "n": len(per_sample),
            "n_healthy": sum(1 for r in per_sample.values()
                             if r["group_label"] == "healthy"),
            "n_cancer": sum(1 for r in per_sample.values()
                            if r["group_label"] != "healthy"),
        },
        "per_sample": dict(sorted(per_sample.items())),
    }
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(f"[ok] {report['samples']['n']} samples "
          f"({report['samples']['n_cancer']} cancer, "
          f"{report['samples']['n_healthy']} healthy); wrote {args.report}")
    return 0


def _index_to_kmer(index: int) -> str:
    letters = "ACGT"
    out = []
    for i in range(4):
        out.append(letters[(index >> (2 * (3 - i))) & 3])
    return "".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
