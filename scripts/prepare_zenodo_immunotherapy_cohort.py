#!/usr/bin/env python3
"""Turn the immunotherapy deposit's fragment files into the repository's features.

What this is
------------
Zenodo 10.5281/zenodo.15738112 ("Genome-Wide Variations of End Motif in
Cell-Free DNA Fragments Distinguish Immunotherapy Responders from
Non-Responders in Head and Neck Cancer: A Multi-Institute Prospective
Study", Bandaru et al., medRxiv 10.64898/2026.03.24.26348354; published
in the Journal of Clinical Investigation) deposits the post-processed
fragment files of **185 longitudinal plasma cfDNA samples from 68
patients** with locally advanced HNSCC who received neoadjuvant and
adjuvant pembrolizumab, across six academic centres.

It is the first deposit the survey found whose endpoint is **clinical
in the treatment-response sense** rather than a cancer-versus-healthy
classification: the readout is *responder versus non-responder*, and
the study is a prospective multi-institutional phase II trial
(NCT02641093).

The channel it publishes is the one this repository's PHTS finding
rests on -- fragment end motifs -- and the study's own metric is a
regional version of it, the **rMDS**: the normalized Shannon entropy of
the 5' end 4-mer frequencies inside 500 kb bins, z-scored across the
genome within each sample. The repository computes that metric here
from the deposited fragments, so the published claim can be tested on
the published data.

Three things are computed per sample
------------------------------------
1. **The 256 4-mer end-motif frequencies** -- the repository's motif
   channel, in the repository's own motif order, with both fragment
   ends counted and each oriented to the + strand.
2. **The per-500 kb-bin rMDS vector** -- the study's metric, computed
   from the same fragment set, so the published claim is tested on the
   published definition rather than on a paraphrase of it.
3. **The DELFI 5-channel features** -- the repository's standard
   external-cohort channel, so this cohort is comparable with the
   FinaleDB cohorts.

The filters, and where they come from
-------------------------------------
The study says its fragments are "only autosomal fragments of 50-350 bp
in length" with "a mapping quality of 30 or greater", and that reads
overlapping the ENCODE blacklist were excluded. **The deposited files
are not filtered that way** -- they carry fragments of every length and
every MAPQ -- so the first two filters are applied here and the
blacklist exclusion is not (the blacklist track is not in the
deposit). How much the deposited file differs from the study's own
filter is measured and written to the report rather than assumed.

The DELFI features are computed by the repository's own extractor
without those extra filters, because that extractor is the repository's
definition of the channel; the two channels are therefore computed
from slightly different fragment sets, and the report says so.

Usage:
    python3 scripts/prepare_zenodo_immunotherapy_cohort.py \
        --frag-dir /Volumes/Download/projects/deepcatch/cohorts/immunotherapy/Frag_Files \
        --reference /Volumes/Download/projects/deepcatch/inputs/hg38/hg38.fa \
        --metadata data/external_cohort/immunotherapy_sample_metadata.tsv \
        --out-dir data/external_cohort_zenodo_immunotherapy \
        --report results/zenodo_immunotherapy_prep.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "src"))

IMMUNOTHERAPY_DOI = "10.5281/zenodo.15738112"
STUDY_DOI = "10.64898/2026.03.24.26348354"

# The study's own filters (JCI 196284, Methods).
MIN_FRAGMENT_LENGTH = 50
MAX_FRAGMENT_LENGTH = 350
MIN_MAPQ = 30
BIN_SIZE = 500_000

# The chromosomes the study says it kept: the autosome.
AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))

# The base codes, in the repository's motif order: A=0, C=1, G=2, T=3.
_BASE_CODE = np.full(256, 255, dtype=np.uint8)
_BASE_CODE[ord("A")] = 0
_BASE_CODE[ord("C")] = 1
_BASE_CODE[ord("G")] = 2
_BASE_CODE[ord("T")] = 3


def _chromosome_codes(fasta: str, chrom: str) -> np.ndarray:
    """The chromosome's sequence as base codes (255 for anything else)."""
    import pysam

    with pysam.FastaFile(fasta) as handle:
        sequence = handle.fetch(chrom).upper()
    raw = np.frombuffer(sequence.encode("ascii", "ignore"), dtype=np.uint8)
    return _BASE_CODE[raw]


def _read_metadata(path: str) -> dict[str, dict[str, str]]:
    """The label table: {sample_id: {column: value}}."""
    import csv

    out: dict[str, dict[str, str]] = {}
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("sample_id") or "").strip()
            if sample:
                out[sample] = {k: (v or "").strip() for k, v in row.items()}
    return out


def _sample_paths(frag_dir: str) -> dict[str, str]:
    """{sample_id: fragment file} for every deposited sample."""
    out: dict[str, str] = {}
    for path in sorted(Path(frag_dir).glob("*.hg38.frag.bed.gz")):
        sample = path.name.split(".")[0]
        out[sample] = str(path)
    return out


def _fragments_for(path: str, chrom: str):
    """The fragments on one chromosome, as (start, end, mapq) arrays.

    The deposited files are bgzipped and tabix-indexed, so one
    chromosome is read without decompressing the whole file. The file is
    a bedtools BEDPE with a header: chrom, start, end, name, mapq, and
    the fragment's strand.
    """
    import pysam

    starts: list[int] = []
    ends: list[int] = []
    mapqs: list[int] = []
    with pysam.TabixFile(path) as tabix:
        try:
            rows = tabix.fetch(chrom)
        except ValueError:
            return None
        for line in rows:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                starts.append(int(fields[1]))
                ends.append(int(fields[2]))
                mapqs.append(int(fields[4]))
            except ValueError:
                continue
    if not starts:
        return None
    return (np.asarray(starts, dtype=np.int64),
            np.asarray(ends, dtype=np.int64),
            np.asarray(mapqs, dtype=np.int64))


def _motif_codes(codes: np.ndarray, start: np.ndarray, end: np.ndarray):
    """The 4-mer code of each fragment's two 5' ends, and the valid mask.

    A fragment's 5' end is its left end when the fragment is on the +
    strand and its right end when it is on the - strand; the second
    end's bases are reverse-complemented so that every motif is read
    on the + strand. The repository's motif index is the two-bit
    encoding of the four bases, A=0, C=1, G=2, T=3.
    """
    length = codes.shape[0]
    left_ok = (start + 4 <= length) & (start >= 0)
    right_ok = (end - 4 >= 0) & (end <= length)
    valid = left_ok & right_ok

    left = np.zeros(start.shape, dtype=np.int64)
    right = np.zeros(start.shape, dtype=np.int64)

    if np.any(left_ok):
        idx = np.where(left_ok)[0]
        window = codes[start[idx][:, None] + np.arange(4)[None, :]]
        valid_bases = window < 255
        left[idx] = np.where(
            valid_bases.all(axis=1),
            (window.astype(np.int64) << (2 * np.arange(3, -1, -1))).sum(axis=1),
            -1)
    if np.any(right_ok):
        idx = np.where(right_ok)[0]
        window = codes[end[idx][:, None] - 4 + np.arange(4)[None, :]]
        # The complement of a base code is 3 - code, and the order is
        # reversed: the reverse complement of the last four bases.
        complemented = np.where(window < 255, 3 - window.astype(np.int64), 255)
        complemented = complemented[:, ::-1]
        valid_bases = complemented < 255
        right[idx] = np.where(
            valid_bases.all(axis=1),
            (complemented << (2 * np.arange(3, -1, -1))).sum(axis=1),
            -1)

    valid &= (left >= 0) & (right >= 0)
    return left, right, valid


def _shannon_entropy(counts: np.ndarray) -> np.ndarray:
    """The normalized Shannon entropy of each row of counts.

    Normalization is by log(256), the entropy of a uniform distribution
    over the 256 4-mers, which is what the published rMDS does: the
    score runs from 0 (one motif only) to 1 (all 256 equally frequent).
    """
    total = counts.sum(axis=1, keepdims=True).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(total > 0, counts / np.maximum(total, 1), 0.0)
        entropy = -(np.where(p > 0, p * np.log(p), 0.0)).sum(axis=1)
    return entropy / np.log(256.0)


def _chromosome_worker(job):
    """Everything one chromosome contributes, for every sample.

    The per-bin motif COUNTS are not kept: a bin's rMDS is a function of
    its counts alone, so the entropy is computed as soon as the
    chromosome is done and only the entropy is carried forward. Keeping
    the counts instead would mean 185 samples x 5,500 bins x 256 counts
    -- two gigabytes of counters that say the same thing as 185 x 5,500
    floats.
    """
    (reference, chrom, paths, samples) = job
    codes = _chromosome_codes(reference, chrom)
    motif_counts = {sample: np.zeros(256, dtype=np.int64) for sample in samples}
    bin_counts: dict[str, dict[int, np.ndarray]] = {s: {} for s in samples}
    fragment_counts = {sample: 0 for sample in samples}
    filtered_counts = {sample: 0 for sample in samples}
    length_histogram = {sample: np.zeros(64, dtype=np.int64) for sample in samples}
    mapq_histogram = {sample: np.zeros(256, dtype=np.int64) for sample in samples}

    for sample in samples:
        found = _fragments_for(paths[sample], chrom)
        if found is None:
            continue
        start, end, mapq = found
        fragment_counts[sample] += start.shape[0]
        lengths = end - start
        length_histogram[sample] += np.bincount(
            np.clip(lengths, 0, 63), minlength=64)
        mapq_histogram[sample] += np.bincount(
            np.clip(mapq, 0, 255), minlength=256)

        keep = ((lengths >= MIN_FRAGMENT_LENGTH)
                & (lengths <= MAX_FRAGMENT_LENGTH)
                & (mapq >= MIN_MAPQ))
        if not np.any(keep):
            continue
        filtered_counts[sample] += int(keep.sum())
        left, right, valid = _motif_codes(codes, start[keep], end[keep])
        valid &= (start[keep] >= 0) & (end[keep] <= codes.shape[0])
        if not np.any(valid):
            continue
        motif_counts[sample] += np.bincount(left[valid], minlength=256)
        motif_counts[sample] += np.bincount(right[valid], minlength=256)

        bins = (start[keep][valid] // BIN_SIZE).astype(np.int64)
        codes_left = left[valid]
        codes_right = right[valid]
        per_bin = bin_counts[sample]
        for bin_index, code in zip(bins, codes_left):
            row = per_bin.get(bin_index)
            if row is None:
                row = np.zeros(256, dtype=np.int64)
                per_bin[bin_index] = row
            row[code] += 1
        for bin_index, code in zip(bins, codes_right):
            row = per_bin.get(bin_index)
            if row is None:
                row = np.zeros(256, dtype=np.int64)
                per_bin[bin_index] = row
            row[code] += 1

    # The entropy of each bin, computed now and the counts dropped.
    entropies: dict[str, dict[int, float]] = {}
    for sample, per_bin in bin_counts.items():
        if not per_bin:
            entropies[sample] = {}
            continue
        keys = sorted(per_bin)
        counts = np.vstack([per_bin[key] for key in keys])
        entropies[sample] = {int(key): float(value) for key, value
                             in zip(keys, _shannon_entropy(counts))}
    return (chrom, motif_counts, entropies, fragment_counts,
            filtered_counts, length_histogram, mapq_histogram)


def _delfi_worker(job):
    """The repository's DELFI features for one sample."""
    (sample, path, out_dir) = job
    from fragmentomics.finale_features import extract_delfi

    delfi = extract_delfi(path)
    for key, value in delfi.items():
        np.save(Path(out_dir) / f"{sample}.{key}.npy", value)
    return sample, {k: int(np.asarray(v).size) for k, v in delfi.items()}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--frag-dir", required=True,
                        help="the directory of deposited fragment files")
    parser.add_argument("--reference", required=True,
                        help="the hg38 reference FASTA (with its .fai)")
    parser.add_argument("--metadata", required=True,
                        help="the sample label table")
    parser.add_argument("--out-dir", required=True,
                        help="where the per-sample features are written")
    parser.add_argument("--report", required=True,
                        help="the JSON report to write")
    parser.add_argument("--max-samples", type=int, default=0,
                        help="stop after this many samples (0: all)")
    parser.add_argument("--workers", type=int, default=0,
                        help="parallel workers (0: one per core, capped)")
    args = parser.parse_args(argv)

    started = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata = _read_metadata(args.metadata)
    paths = _sample_paths(args.frag_dir)
    samples = sorted(set(metadata) & set(paths))
    if args.max_samples:
        samples = samples[:args.max_samples]
    missing = sorted(set(metadata) - set(paths))
    extra = sorted(set(paths) - set(metadata))

    workers = args.workers or max(1, min((os.cpu_count() or 2) - 1, 6))
    print(f"{len(samples)} samples have both a label and a fragment file; "
          f"{workers} workers", flush=True)
    if missing:
        print(f"  {len(missing)} labelled samples have no file: "
              f"{missing[:5]}{'...' if len(missing) > 5 else ''}", flush=True)
    if extra:
        print(f"  {len(extra)} deposited files have no label: "
              f"{extra[:5]}{'...' if len(extra) > 5 else ''}", flush=True)

    # The reference is read once per chromosome, and every sample is
    # computed against it in the same pass: reading the reference per
    # sample would read 24 chromosomes x 185 samples of a 3 GB genome.
    motif_counts = {sample: np.zeros(256, dtype=np.int64) for sample in samples}
    bin_entropy: dict[str, dict[tuple[str, int], float]] = {
        sample: {} for sample in samples}
    fragment_counts = {sample: 0 for sample in samples}
    filtered_counts = {sample: 0 for sample in samples}
    length_histogram = {sample: np.zeros(64, dtype=np.int64) for sample in samples}
    mapq_histogram = {sample: np.zeros(256, dtype=np.int64) for sample in samples}

    jobs = [(args.reference, chrom, paths, samples) for chrom in AUTOSOMES]
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(workers) as pool:
            results = []
            for result in pool.imap_unordered(_chromosome_worker, jobs):
                results.append(result)
                print(f"  {result[0]}: done", flush=True)
    else:
        results = []
        for job in jobs:
            result = _chromosome_worker(job)
            results.append(result)
            print(f"  {result[0]}: done", flush=True)

    for (_chrom, motifs, entropies, fragments, filtered, lengths,
         mapqs) in results:
        for sample in samples:
            motif_counts[sample] += motifs[sample]
            fragment_counts[sample] += fragments[sample]
            filtered_counts[sample] += filtered[sample]
            length_histogram[sample] += lengths[sample]
            mapq_histogram[sample] += mapqs[sample]
            for bin_index, value in entropies[sample].items():
                bin_entropy[sample][(_chrom, bin_index)] = value

    print("computing the per-sample features ...", flush=True)
    written = []
    for sample in samples:
        motifs = motif_counts[sample]
        total = motifs.sum()
        frequencies = motifs / total if total else motifs.astype(np.float64)
        np.save(out_dir / f"{sample}.motifs.npy", frequencies)

        bins = bin_entropy[sample]
        if bins:
            keys = sorted(bins)
            entropy = np.asarray([bins[key] for key in keys])
            z = (entropy - entropy.mean()) / (entropy.std() or 1.0)
            np.save(out_dir / f"{sample}.rmds.npy", z.astype(np.float64))
            np.save(out_dir / f"{sample}.rmds_bins.npy",
                    np.asarray([f"{key[0]}:{key[1]}" for key in keys]))
        else:
            np.save(out_dir / f"{sample}.rmds.npy", np.zeros(0))
            np.save(out_dir / f"{sample}.rmds_bins.npy", np.zeros(0))
        written.append(sample)
        print(f"  {sample}: {total} end motifs, "
              f"{len(bins)} bins, {fragment_counts[sample]} fragments",
              flush=True)

    # The DELFI channel, from the repository's own extractor, one sample
    # per worker: it reads the whole fragment file, so it is the part
    # that parallelises.
    delfi_jobs = [(sample, paths[sample], str(out_dir)) for sample in written]
    if workers > 1 and delfi_jobs:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(workers) as pool:
            for sample, _sizes in pool.imap_unordered(_delfi_worker,
                                                     delfi_jobs):
                print(f"  {sample}: DELFI features done", flush=True)
    else:
        for job in delfi_jobs:
            sample, _sizes = _delfi_worker(job)
            print(f"  {sample}: DELFI features done", flush=True)

    report = {
        "doi": IMMUNOTHERAPY_DOI,
        "study_doi": STUDY_DOI,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "reference": os.path.abspath(args.reference),
        "filters": {
            "min_fragment_length": MIN_FRAGMENT_LENGTH,
            "max_fragment_length": MAX_FRAGMENT_LENGTH,
            "min_mapq": MIN_MAPQ,
            "chromosomes": list(AUTOSOMES),
            "blacklist_excluded": False,
        },
        "samples": {
            "labelled": len(metadata),
            "with_file": len(paths),
            "prepared": len(written),
            "labelled_without_file": missing,
            "file_without_label": extra,
        },
        "per_sample": {
            sample: {
                "fragments_in_file": int(fragment_counts[sample]),
                "fragments_after_filters": int(filtered_counts[sample]),
                "end_motifs": int(motif_counts[sample].sum()),
                "bins": len(bin_entropy[sample]),
                "length_histogram": length_histogram[sample].tolist(),
                "mapq_histogram": mapq_histogram[sample].tolist(),
            }
            for sample in written
        },
        "elapsed_seconds": round(time.time() - started, 1),
    }
    Path(args.report).parent.mkdir(parents=True, exist_ok=True)
    Path(args.report).write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
