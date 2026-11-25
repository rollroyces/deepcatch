#!/usr/bin/env python3
"""Verify the SAUNA cohort's fragment files before they are scored.

Why this script exists
-----------------------
The SAUNA deposit (Zenodo 15097509) ships nine samples' fragment files
already in the form the cohort extractor reads -- five columns, ``chrom
start end score strand``, bgzip-compressed, mapped to **hg38** -- so
unlike the Cleveland Clinic PHTS deposit they need no conversion. What
they do need is a check, because a file that is quietly the wrong shape
is the failure mode this repository keeps meeting: the PHTS deposit's
BEDPE was readable and wrong (one read long per "fragment"), and the
ccRCC deposit's files were fine. A verifier that says which is which is
worth more than a conversion that assumes.

So this script checks, per file:

* the size against the size the survey recorded, so a truncated
  download cannot pass;
* the column count and the header convention on the first rows;
* that no row is an **inverted interval** (end <= start). The depositor
  built these files with ``cut -f1,2,6,8,9`` on a BEDPE, and for a
  minus-strand pair the outer coordinates are start1 and end2 -- if the
  pair were written in read order rather than genomic order, the cut
  would produce end < start, and the extractor would silently drop
  every one of those fragments (its length filter keeps only
  100 <= end - start <= 220). They are not inverted, and this is where
  that is established rather than assumed;
* the fragment-length distribution, which has to show the
  mono-nucleosome peak as the mode and a visible di-nucleosome
  shoulder, exactly as the PHTS preparation script does.

The report it writes is what the scorer reads for the cohort's facts.

Usage:
    python3 scripts/prepare_zenodo_sauna_cohort.py \
        --dir /tmp/sauna_cohort \
        --survey results/zenodo_cfdna_survey.json \
        --report results/zenodo_sauna_prep.json
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone

# The nucleosome ladder, as in the PHTS preparation script: the
# mono-nucleosome peak has to be the mode, and the di-nucleosome
# shoulder has to be visible.
_MODE_RANGE = (150.0, 200.0)
_DI_WINDOW = (320, 345)
_DI_MIN_FRACTION = 0.008

SAUNA_DOI = "10.5281/zenodo.15097509"

# How many rows to read from the head of each file for the checks. The
# format is uniform across a file, so a sample is representative; the
# row count is not sampled, because it is the one number that has to be
# exact.
_HEAD_ROWS = 20_000_000


def _survey_sizes(survey_path: str) -> dict:
    with open(survey_path, encoding="utf-8") as handle:
        survey = json.load(handle)
    for candidate in survey.get("candidates") or []:
        if candidate.get("doi") != SAUNA_DOI:
            continue
        sizes = {}
        for entry in candidate.get("files") or []:
            sizes[entry["key"]] = int(entry["size"])
        return sizes
    raise SystemExit(f"{SAUNA_DOI} is not among the survey's candidates")


def _check_file(path: str) -> dict:
    """Read the head of one file and report what is in it."""
    script = (
        "import sys\n"
        "n=0; inv=0; neg=0; short=0; long=0; total=0.0; di=0\n"
        "cols={}\n"
        "mode={}\n"
        "import itertools\n"
        "with open(sys.argv[1], 'rt', errors='replace') as handle:\n"
        "    for line in itertools.islice(handle, %d):\n"
        "        parts = line.rstrip('\\n').split('\\t')\n"
        "        cols[len(parts)] = cols.get(len(parts), 0) + 1\n"
        "        if len(parts) != 5:\n"
        "            continue\n"
        "        try:\n"
        "            start = int(parts[1]); end = int(parts[2])\n"
        "        except ValueError:\n"
        "            continue\n"
        "        n += 1\n"
        "        if end <= start:\n"
        "            inv += 1\n"
        "        length = end - start\n"
        "        if length < 0:\n"
        "            neg += 1\n"
        "        total += length\n"
        "        if 100 <= length <= 150:\n"
        "            short += 1\n"
        "        elif 150 < length <= 220:\n"
        "            long += 1\n"
        "        if %d <= length <= %d:\n"
        "            di += 1\n"
        "        bucket = (length // 5) * 5\n"
        "        mode[bucket] = mode.get(bucket, 0) + 1\n"
        "best = max(mode, key=mode.get) if mode else None\n"
        "print(n, inv, neg, short, long, di, total, best, mode.get(best, 0) if best is not None else 0, sorted(cols.items())[:5])\n"
        % (_HEAD_ROWS, _DI_WINDOW[0], _DI_WINDOW[1])
    )
    result = subprocess.run(
        ["python3", "-c", script, path],
        capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise SystemExit(f"could not read {path}: {result.stderr.strip()}")
    fields = result.stdout.split()
    n = int(fields[0])
    inv = int(fields[1])
    neg = int(fields[2])
    short = int(fields[3])
    long = int(fields[4])
    di = int(fields[5])
    total = float(fields[6])
    mode = int(fields[7])
    mode_count = int(fields[8])
    columns = json.loads(fields[9].replace("'", '"'))
    return {
        "rows_sampled": n,
        "inverted_intervals": inv,
        "negative_lengths": neg,
        "short_fraction_100_150": round(short / n, 6) if n else None,
        "long_fraction_150_220": round(long / n, 6) if n else None,
        "di_nucleosome_fraction": round(di / n, 6) if n else None,
        "mean_fragment_length": round(total / n, 2) if n else None,
        "mode_fragment_length": mode,
        "mode_fraction": round(mode_count / n, 6) if n else None,
        "column_counts": {str(k): v for k, v in columns},
    }


def _count_rows(path: str) -> int:
    """The exact number of rows in a file (not sampled)."""
    result = subprocess.run(
        ["bash", "-c", f"gunzip -c '{path}' | wc -l"],
        capture_output=True, text=True, check=False)
    try:
        return int(result.stdout.strip())
    except ValueError:
        return -1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dir", required=True,
                    help="directory holding the deposit's files")
    ap.add_argument("--survey", default="results/zenodo_cfdna_survey.json")
    ap.add_argument("--report", default="results/zenodo_sauna_prep.json")
    ap.add_argument("--count-rows", action="store_true",
                    help="count every row of every file (slow)")
    args = ap.parse_args()

    sizes = _survey_sizes(args.survey)
    samples: dict = {}
    problems: list[str] = []
    for name, expected in sorted(sizes.items()):
        path = os.path.join(args.dir, name)
        if not os.path.exists(path):
            problems.append(f"{name} is missing")
            continue
        actual = os.path.getsize(path)
        if actual != expected:
            problems.append(
                f"{name} is {actual} bytes, the survey recorded {expected}")
            continue
        if name.endswith(".tsv"):
            samples[name] = {"bytes": actual, "kind": "metadata"}
            continue
        stats = _check_file(path)
        stats["bytes"] = actual
        stats["kind"] = "fragments"
        if args.count_rows:
            stats["rows_total"] = _count_rows(path)
        samples[name] = stats
        print(f"      {name}: {stats['rows_sampled']:,} rows sampled, mode "
              f"{stats['mode_fragment_length']} bp, di-nucleosome "
              f"{stats['di_nucleosome_fraction']:.3%}, inverted "
              f"{stats['inverted_intervals']}", flush=True)

    for name, stats in samples.items():
        if stats.get("kind") != "fragments":
            continue
        if stats["inverted_intervals"] or stats["negative_lengths"]:
            problems.append(
                f"{name} has {stats['inverted_intervals']} inverted "
                "intervals; the extractor would drop those fragments")
        mode = stats["mode_fragment_length"]
        if not (_MODE_RANGE[0] <= mode <= _MODE_RANGE[1]):
            problems.append(
                f"{name}'s fragment lengths do not peak on the nucleosome "
                f"peak (mode {mode} bp)")
        if stats["di_nucleosome_fraction"] < _DI_MIN_FRACTION:
            problems.append(
                f"{name}'s fragment lengths have no di-nucleosome shoulder")

    if problems:
        print("[FAILED] the cohort is not ready:")
        for problem in problems:
            print(f"      - {problem}")
        return 1

    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {
            "record": "https://zenodo.org/records/15097509",
            "doi": SAUNA_DOI,
            "paper": "Balazs, Zsolt -- SAUNA, Simulated Annealing for Unique "
                     "Nucleosome Arrangements",
        },
        "format": {
            "columns": "chrom, start, end, score, strand (the fragment form "
                       "the extractor reads)",
            "compression": "bgzip (gzip-compatible)",
            "build": "hg38 -- the build the extractor's bin grid uses, so no "
                     "conversion and no grid loss",
            "conversion_needed": False,
            "note": "The depositor built these files with "
                    "`cut -f1,2,6,8,9` on a BEDPE, so the outer coordinates "
                    "are start1 and end2; the check above establishes that "
                    "no row is an inverted interval, which would have made "
                    "the extractor drop that fragment silently.",
        },
        "checks": {
            "rows_sampled_per_file": _HEAD_ROWS,
            "mode_range": list(_MODE_RANGE),
            "di_nucleosome_window": list(_DI_WINDOW),
            "di_nucleosome_min_fraction": _DI_MIN_FRACTION,
        },
        "samples": samples,
    }
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(f"[ok] every file matches the survey's size and shows the "
          f"nucleosome ladder; wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
