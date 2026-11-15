#!/usr/bin/env python3
"""Check whether a cfDNA library carries the digestion signature FinaleMe needs.

Why
---
FinaleMe predicts DNA methylation from the *fragmentation pattern* of
cfDNA: in a MspI-digested library an unmethylated CpG is cut out of the
molecule, so it sits at a fragment end, while a methylated CpG is not cut
and stays in the interior of a fragment. The per-CpG methylation state is
therefore only recoverable where the fragment ends still mark the CpG
sites.

Two things destroy that information, and both are visible in the fragment
file alone, without running the decoder:

1. **Size selection.** If the library was trimmed to a maximum length,
   the fragment ends are cut positions rather than digestion sites, so no
   CpG sits at a fragment end any more. This shows up as a hard ceiling
   in the fragment length distribution.
2. **Depth.** A methylation state is only identified where a CpG is
   observed repeatedly. The number of observations a CpG receives is
   roughly (fragments covering it) / (CpGs in the region), so a shallow
   library cannot identify any state, however well prepared it is.

This script reports both, plus the direct measurement of the digestion
signature: the share of fragment ends that land on a CpG, compared with
the share expected by chance.

Usage
-----
    python scripts/methylation_library_check.py \
        --frag-files /path/a.frag.tsv.bgz /path/b.frag.tsv.bgz \
        [--chrom 22] [--cpg-bedgraph CG_motif.hg19.chr22.bedgraph.gz] \
        --out-json results/methylation_library_check.json
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import os
import sys
from datetime import datetime, timezone

# Fragment length bins (bp).
LENGTH_BINS = [(0, 150), (150, 200), (200, 250), (250, 300), (300, 400),
               (400, 500), (500, 700), (700, 1000), (1000, 2000),
               (2000, 10 ** 9)]
LENGTH_LABELS = ["<150", "150-200", "200-250", "250-300", "300-400",
                 "400-500", "500-700", "700-1k", "1k-2k", ">2k"]
# A library whose longest fragments are bounded near this value has been
# size-selected; the digestion signature is truncated with it.
SIZE_SELECTION_CEILING = 300


def normalise_chrom(name: str) -> str:
    """Compare chromosome names regardless of a ``chr`` prefix.

    The FinaleDB and CRAG fragment files use bare names (``22``) while the
    reference CpG bedgraphs use UCSC names (``chr22``).
    """
    return name[3:] if name.startswith("chr") else name


def load_cpg_positions(path: str, chrom: str) -> list[int]:
    """0-based CpG start positions from a bedgraph of 1-bp intervals."""
    positions = []
    target = normalise_chrom(chrom)
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if chrom and normalise_chrom(parts[0]) != target:
                continue
            try:
                positions.append(int(parts[1]))
            except ValueError:
                continue
    return sorted(positions)


def nearest_distance(positions: list[int], p: int) -> int:
    i = bisect.bisect_left(positions, p)
    best = 10 ** 9
    for j in (i - 1, i):
        if 0 <= j < len(positions):
            best = min(best, abs(positions[j] - p))
    return best


def check_file(path: str, chrom: str, cpg_positions: list[int]) -> dict:
    n = 0
    at_cpg_start = 0
    at_cpg_end = 0
    length_counts = [0] * len(LENGTH_BINS)
    lengths = []
    target = normalise_chrom(chrom) if chrom else None
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if target and normalise_chrom(parts[0]) != target:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            n += 1
            length = end - start
            lengths.append(length)
            for i, (lo, hi) in enumerate(LENGTH_BINS):
                if lo <= length < hi:
                    length_counts[i] += 1
                    break
            if cpg_positions:
                if nearest_distance(cpg_positions, start) == 0:
                    at_cpg_start += 1
                if nearest_distance(cpg_positions, end) == 0:
                    at_cpg_end += 1
    if not n:
        return {"file": os.path.basename(path), "n_fragments": 0}
    lengths.sort()
    out = {
        "file": os.path.basename(path),
        "n_fragments": n,
        "fragment_length_median": int(lengths[len(lengths) // 2]),
        "fragment_length_mean": round(sum(lengths) / len(lengths), 2),
        "fragment_length_max": int(lengths[-1]),
        "fraction_at_or_above_300bp": round(
            sum(length_counts[i] for i, (lo, hi) in enumerate(LENGTH_BINS)
                if lo >= SIZE_SELECTION_CEILING) / n, 6),
        "length_profile": {LENGTH_LABELS[i]: round(length_counts[i] / n, 6)
                            for i in range(len(LENGTH_BINS))},
    }
    if cpg_positions:
        # Chance rate: how often a random position sits on a CpG, taken
        # from the CpG list's own density over its span.
        span = cpg_positions[-1] - cpg_positions[0]
        density = len(cpg_positions) / max(span, 1)
        out["fragment_ends_on_a_cpg"] = {
            "start_rate": round(at_cpg_start / n, 6),
            "end_rate": round(at_cpg_end / n, 6),
            "chance_rate_per_position": round(density, 6),
            "start_enrichment_over_chance": round(
                (at_cpg_start / n) / density, 3) if density else None,
            "end_enrichment_over_chance": round(
                (at_cpg_end / n) / density, 3) if density else None,
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frag-files", nargs="+", required=True,
                    help="fragment files (*.frag.tsv.bgz or *.frag.bed.gz)")
    ap.add_argument("--chrom", default="22",
                    help="chromosome to inspect (default 22; FinaleDB uses "
                         "bare names, the CRAG files too)")
    ap.add_argument("--cpg-bedgraph", default=None,
                    help="optional bedgraph of CpG positions; enables the "
                         "digestion-signature measurement")
    ap.add_argument("--out-json", default="results/methylation_library_check.json")
    args = ap.parse_args()

    paths = []
    for p in args.frag_files:
        if os.path.isdir(p):
            paths.extend(os.path.join(p, n) for n in sorted(os.listdir(p))
                         if n.endswith(".gz"))
        else:
            paths.append(p)
    if not paths:
        print("no fragment files given", file=sys.stderr)
        return 1

    cpg_positions = []
    if args.cpg_bedgraph:
        cpg_positions = load_cpg_positions(args.cpg_bedgraph, args.chrom)
        print(f"loaded {len(cpg_positions)} CpG positions on chr{args.chrom}")

    results = []
    for path in paths:
        try:
            results.append(check_file(path, args.chrom, cpg_positions))
        except (OSError, EOFError) as exc:
            print(f"skipping {path}: {exc}", file=sys.stderr)
            continue

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "chromosome": args.chrom,
        "libraries": results,
        "interpretation": {
            "size_selection": "fraction_at_or_above_300bp near zero means the "
                              "library was size-selected; fragment ends are "
                              "then truncation points, not digestion sites, "
                              "and no CpG marks an unmethylated site",
            "digestion_signature": "start_enrichment_over_chance near 1 means "
                                    "fragment ends are unrelated to CpGs, so "
                                    "the methylation state is not identifiable "
                                    "from the fragment geometry",
            "depth": "fragments per CpG is the number of observations a CpG "
                     "receives; a state needs more than one",
        },
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out_json)), exist_ok=True)
    with open(args.out_json, "w") as f:
        json.dump(payload, f, indent=2)
    print(json.dumps(payload["libraries"], indent=2))
    print(f"\nwrote {args.out_json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
