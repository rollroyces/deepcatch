#!/usr/bin/env python3
"""Prepare the Cleveland Clinic PHTS cohort for the cohort scorer.

What this does
--------------
The deposit ships one tarball of paired-end fragment files
(``02-bamtobed.tar.gz``, 7.6 GB) whose members are **BEDPE** files: one
row per read *pair*, with the two mates' coordinates in separate
columns. The cohort extractor (``src/fragmentomics/finale_features.py``)
reads a *fragment* BED -- one row per fragment, with the fragment's
outer coordinates in columns 2 and 3 -- and the scorer looks for
``<sample_id>.<build>.frag.bed.gz``, the name the ccRCC deposit shipped.

So this script unpacks the tarball and rewrites each member into the
fragment form the extractor expects: the fragment's outer coordinates
(min start, max end), the pair's MAPQ and the strand, gzipped on the
way through. Nothing is parsed or altered beyond that -- the bases are
the deposit's.

Why the conversion is not optional
----------------------------------
Feeding a BEDPE to the extractor as-is does not fail. It reads quietly,
and it is wrong: the extractor's first three columns are the FIRST
mate's coordinates, so every "fragment" would be one read long (about
146 bp instead of about 167 bp) and the fragment-size distribution would
collapse onto the read-length peak instead of the nucleosome ladder.
The check at the end of this script is the guard against exactly that:
every converted sample has to show the nucleosome ladder -- the
mono-nucleosome peak as the mode, and a visible di-nucleosome
shoulder -- and the script refuses to report success if it does not.

The genome build
----------------
The deposit's fragments are mapped to **hg19**, while the extractor's
bin grid uses hg38 chromosome sizes (the grid the training cohort's
features were computed on, so the held-out cohort must be computed on
it too). A bin index is a coordinate division, so the two builds agree
on every bin except the tail of each chromosome that hg19 has and hg38
does not. This script measures that tail instead of assuming it: it
records the observed maximum coordinate of every chromosome, checks
them against both builds' lengths, and counts the fragments the hg38
grid would drop. The report it writes says which build the data is and
how much of it the grid discards.

Usage
-----
    python3 scripts/prepare_zenodo_phts_cohort.py \
        --tarball /tmp/phts_cohort/02-bamtobed.tar.gz \
        --out-dir /tmp/phts_cohort/frag \
        --metadata data/external_cohort/phts_sample_metadata.tsv \
        --report results/zenodo_phts_prep.json
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import tarfile
from datetime import datetime, timezone

import numpy as np
import pandas as pd

# The extractor reads a gzipped fragment BED and the scorer looks for
# this name; the deposit ships "<sample_id>.bedpe".
TARGET_SUFFIX = ".hg19.frag.bed.gz"

# hg38 chromosome sizes -- the extractor's grid (src/fragmentomics/
# finale_features.py). Kept here so the report can say exactly how much
# of an hg19 deposit the grid drops, instead of leaving it implicit.
_HG38 = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

# hg19 chromosome sizes -- the build the deposit is mapped to.
_HG19 = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
}

# The nucleosome ladder. A converted sample's fragment lengths must show
# the mono-nucleosome peak as the MODE and a visible di-nucleosome
# shoulder. A sample read as if the BEDPE's first-mate columns were
# fragment columns would show neither: its lengths would pile up at the
# read length and the di-nucleosome window would be empty. These are the
# numbers the guard below checks.
_MODE_RANGE = (150.0, 200.0)
_DI_WINDOW = (320, 345)
_DI_MIN_FRACTION = 0.008

_CHUNK = 1_000_000


def _bedpe_columns() -> list[str]:
    return ["chrom1", "start1", "end1", "chrom2", "start2", "end2",
            "name", "mapq", "strand1", "strand2", "tlen"]


def convert_member(source, target: str, stats: dict) -> None:
    """Rewrite one BEDPE member as a gzipped fragment BED.

    The fragment's outer coordinates are the mates' outer coordinates:
    ``min(start1, start2)`` to ``max(end1, end2)``. A pair on two
    different chromosomes is not a fragment and is dropped, as is a row
    whose coordinates do not increase.
    """
    reader = pd.read_csv(source, sep="\t", header=None,
                         names=_bedpe_columns(),
                         usecols=[0, 1, 2, 3, 4, 5, 7, 8],
                         chunksize=_CHUNK)
    rows = 0
    dropped = 0
    lengths: list[np.ndarray] = []
    max_coord: dict[str, int] = {}
    beyond_hg38 = 0
    with gzip.open(target, "wt") as out:
        for chunk in reader:
            chrom1 = chunk["chrom1"].to_numpy()
            start = np.minimum(chunk["start1"].to_numpy(),
                               chunk["start2"].to_numpy())
            end = np.maximum(chunk["end1"].to_numpy(), chunk["end2"].to_numpy())
            keep = (chrom1 == chunk["chrom2"].to_numpy()) & (end > start)
            dropped += int((~keep).sum())
            chrom = chrom1[keep]
            start = start[keep]
            end = end[keep]
            if start.size == 0:
                continue
            lengths.append(end - start)
            rows += int(start.size)

            # The build evidence: the largest coordinate seen on each
            # chromosome, and how many fragments start past the hg38
            # grid's last bin on that chromosome.
            uniq, inverse = np.unique(chrom, return_inverse=True)
            first = np.flatnonzero(np.append(True, inverse[1:] != inverse[:-1]))
            for name, value in zip(uniq, np.maximum.reduceat(end, first)):
                key = str(name)
                if int(value) > max_coord.get(key, 0):
                    max_coord[key] = int(value)
            bounds = np.array([_HG38.get(str(c), 0) for c in chrom])
            beyond_hg38 += int((start >= bounds).sum())

            pd.DataFrame({
                "chrom": chrom,
                "start": start,
                "end": end,
                "mapq": chunk["mapq"].to_numpy()[keep],
                "strand": chunk["strand1"].to_numpy()[keep],
            }).to_csv(out, sep="\t", header=False, index=False,
                      lineterminator="\n")
    all_lengths = np.concatenate(lengths) if lengths else np.zeros(0, dtype=int)
    stats["rows"] = rows
    stats["rows_dropped"] = dropped
    stats["median_fragment_length"] = (float(np.median(all_lengths))
                                        if all_lengths.size else None)
    stats["mean_fragment_length"] = (float(all_lengths.mean())
                                     if all_lengths.size else None)
    stats["p10_fragment_length"] = (float(np.percentile(all_lengths, 10))
                                    if all_lengths.size else None)
    stats["p90_fragment_length"] = (float(np.percentile(all_lengths, 90))
                                    if all_lengths.size else None)
    if all_lengths.size:
        bins = all_lengths // 5
        counts = np.bincount(bins)
        stats["mode_fragment_length"] = float(int(np.argmax(counts)) * 5)
        in_di = ((all_lengths >= _DI_WINDOW[0])
                 & (all_lengths <= _DI_WINDOW[1]))
        stats["di_nucleosome_fraction"] = float(in_di.mean())
    else:
        stats["mode_fragment_length"] = None
        stats["di_nucleosome_fraction"] = None
    stats["max_coordinate_by_chromosome"] = dict(sorted(max_coord.items()))
    stats["fragments_beyond_hg38_grid"] = beyond_hg38
    stats["fraction_beyond_hg38_grid"] = (
        round(beyond_hg38 / rows, 6) if rows else None)


def _build_verdict(max_coord: dict[str, int]) -> dict:
    """Say which genome build the observed coordinates imply.

    A chromosome whose observed maximum runs past hg38's length but
    stays inside hg19's is hg19. A chromosome where the two builds
    agree cannot tell the two apart, so it is recorded as agreeing.
    """
    beyond_hg38 = []
    inside_hg19 = []
    disagree = []
    for chrom, observed in max_coord.items():
        if chrom not in _HG38 or chrom not in _HG19:
            continue
        if observed > _HG38[chrom]:
            beyond_hg38.append(chrom)
        if observed <= _HG19[chrom]:
            inside_hg19.append(chrom)
        if _HG38[chrom] != _HG19[chrom]:
            disagree.append(chrom)
    verdict = "hg19" if beyond_hg38 else "hg38"
    return {
        "verdict": verdict,
        "chromosomes_past_the_hg38_length": sorted(beyond_hg38),
        "chromosomes_within_the_hg19_length": sorted(inside_hg19),
        "n_chromosomes_where_the_builds_differ": len(disagree),
        "max_observed_chr1": max_coord.get("chr1"),
        "hg38_chr1_length": _HG38["chr1"],
        "hg19_chr1_length": _HG19["chr1"],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tarball", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--metadata", required=True,
                    help="the cohort's label table (tab separated, "
                         "keyed by sample_id)")
    ap.add_argument("--report", default="results/zenodo_phts_prep.json")
    args = ap.parse_args()

    with open(args.metadata, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = [line.split("\t")[0] for line in handle if line.strip()]
    print(f"[1/4] {len(samples)} labelled samples in the metadata")

    os.makedirs(args.out_dir, exist_ok=True)
    per_sample: dict[str, dict] = {}
    seen = 0
    truncated = False
    with tarfile.open(args.tarball, "r:gz") as tar:
        try:
            for member in tar:
                if not member.isfile():
                    continue
                sample = os.path.basename(member.name).removesuffix(".bedpe")
                if sample not in samples:
                    continue
                source = tar.extractfile(member)
                assert source is not None
                target = os.path.join(args.out_dir, f"{sample}{TARGET_SUFFIX}")
                stats: dict = {}
                convert_member(source, target, stats)
                per_sample[sample] = stats
                seen += 1
                print(f"      {sample}: {stats['rows']:,} fragments, mode "
                      f"{stats['mode_fragment_length']:.0f} bp, di-nucleosome "
                      f"{stats['di_nucleosome_fraction']:.3%}", flush=True)
        except (EOFError, tarfile.ReadError, OSError) as error:
            # A download that is still running ends mid-member. That is
            # not a crash: the completeness check below is what decides
            # whether the cohort is complete.
            truncated = True
            print(f"      the archive ends early ({type(error).__name__}: "
                  f"{error}) -- it is a partial download")
    print(f"[2/4] converted {seen}/{len(samples)}")
    if truncated:
        print("      the archive is truncated, so the completeness check "
              "below will fail -- download the whole 7.6 GB first")

    missing = [s for s in samples
               if not os.path.exists(
                   os.path.join(args.out_dir, f"{s}{TARGET_SUFFIX}"))]
    if missing:
        print(f"[3/4] FAILED: {len(missing)} labelled samples have no file in "
              f"the archive: {missing[:5]}"
              + (" ..." if len(missing) > 5 else ""))
        print("      (a partial download looks exactly like this)")
        return 1
    # The guard: a converted sample's fragment lengths must show the
    # nucleosome ladder -- the mono-nucleosome peak as the mode, and a
    # visible di-nucleosome shoulder. If they do not, the coordinates
    # were read wrong and every channel computed from them is garbage.
    bad: dict[str, dict] = {}
    for sample, st in per_sample.items():
        mode = st.get("mode_fragment_length")
        di = st.get("di_nucleosome_fraction")
        if (mode is None or not (_MODE_RANGE[0] <= mode <= _MODE_RANGE[1])
                or di is None or di < _DI_MIN_FRACTION):
            bad[sample] = {"mode": mode, "di_nucleosome_fraction": di}
    if bad:
        print(f"[3/4] FAILED: {len(bad)} samples do not show the nucleosome "
              f"ladder (mode in {_MODE_RANGE}, di-nucleosome fraction >= "
              f"{_DI_MIN_FRACTION}): "
              + ", ".join(f"{s}(mode={v['mode']}, di={v['di_nucleosome_fraction']})"
                         for s, v in list(bad.items())[:5]))
        print("      (a fragment-size distribution without the nucleosome "
              "ladder means the BEDPE columns were read as if they were "
              "fragment columns)")
        return 1
    print(f"[3/4] every converted sample shows the nucleosome ladder "
          f"(mode in {_MODE_RANGE}, di-nucleosome fraction >= "
          f"{_DI_MIN_FRACTION})")

    max_coord: dict[str, int] = {}
    for stats in per_sample.values():
        for chrom, value in stats["max_coordinate_by_chromosome"].items():
            if value > max_coord.get(chrom, 0):
                max_coord[chrom] = value
    verdict = _build_verdict(max_coord)
    total_rows = sum(st["rows"] for st in per_sample.values())
    total_beyond = sum(st["fragments_beyond_hg38_grid"]
                       for st in per_sample.values())
    report = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": {
            "record": "https://zenodo.org/records/15109495",
            "doi": "10.5281/zenodo.15109495",
            "archive_member": "02-bamtobed.tar.gz",
        },
        "conversion": {
            "from": "BEDPE (one row per read pair, 11 columns)",
            "to": "fragment BED (chrom, start, end, mapq, strand), gzipped",
            "rule": "outer coordinates = min(start1, start2) .. max(end1, end2); "
                    "pairs on two chromosomes are dropped",
            "target_suffix": TARGET_SUFFIX,
        },
        "genome_build": {
            **verdict,
            "note": "The extractor's bin grid uses hg38 chromosome sizes, "
                    "so the hg38 grid is what the features are computed on "
                    "-- the same grid the training cohort's features use. "
                    "Fragments starting past an hg38 chromosome end are "
                    "dropped by that grid.",
            "hg38_grid_fragments_dropped": total_beyond,
            "hg38_grid_fraction_dropped": (
                round(total_beyond / total_rows, 6) if total_rows else None),
        },
        "samples": {
            "n_labelled": len(samples),
            "n_converted": len(per_sample),
            "total_fragments": total_rows,
        },
        "per_sample": dict(sorted(per_sample.items())),
    }
    os.makedirs(os.path.dirname(args.report) or ".", exist_ok=True)
    with open(args.report, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
    print(f"[4/4] the data is {verdict['verdict']}; the hg38 grid drops "
          f"{total_beyond:,}/{total_rows:,} fragments "
          f"({report['genome_build']['hg38_grid_fraction_dropped']:.4%})")
    print(f"      wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
