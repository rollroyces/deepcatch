"""FinaleDB fragment-file feature extraction (self-contained).

This module reproduces the exact 5-channel DELFI feature set that the
cross-study benchmark was built on, directly from a FinaleDB
``*.frag.tsv.bgz`` file, with no dependency on the sibling
``cfdna-fragmentomics-pipeline`` checkout.

Why it lives here
-----------------
The headline benchmark's features were extracted by the sibling repo's
``scripts/extract_delfi.py`` + ``scripts/extract_fsd.py``. The FSD
extractor there expects the 6-column bedtools ``bamtobed -bedpe``
schema (``chrom start end name mapq strand``), but FinaleDB currently
serves a 5-column schema (``chrom start end mapq strand``). Its
``len(parts) < 6`` guard therefore drops **every** fragment and the
extraction dies with "no fragments in input". This module accepts both
schemas.

The two schemas are equivalent for our purposes: ``start``/``end`` sit
at indices 1/2 in both, and with the default ``mapq_threshold=0`` every
fragment is kept regardless of which index holds the MAPQ. So a
histogram computed here matches the one the sibling extractor produced
from the 6-column files, and the two feature sets are directly
comparable.

Reference
---------
- DELFI short/long binning: Shen et al. 2018 (PNAS) / DELFI protocol.
- Fragment-size distribution: Mouliere et al. 2018 (Nature).
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pandas as pd

# hg38 chromosome sizes (canonical chromosomes 1-22, X, Y) -- identical to
# the sibling repo's extract_delfi.py so the bin grids match exactly.
CHROM_SIZES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
BIN_100KB = 100_000
BIN_5MB = 5_000_000

# DELFI fragment classes
SHORT_MAX = 150      # 100-150 bp
LONG_MAX = 220       # 150-220 bp
MIN_LEN = 100

# Fragment-size distribution bins (196 x 5 bp from 20 bp)
FSD_BIN_START = 20
FSD_BIN_END = 1000
FSD_BIN_STRIDE = 5
N_BINS = (FSD_BIN_END - FSD_BIN_START) // FSD_BIN_STRIDE  # 196

DELFI_FEATURES = (
    "delfi_5mb_ratio",
    "delfi_5mb_coverage",
    "delfi_100kb_ratio",
    "delfi_100kb_counts",
)


def _bin_map(chrom: str, bin_size: int) -> list[tuple[int, int]]:
    """chrom -> [(start, end), ...] windows covering the chromosome."""
    n = (CHROM_SIZES[chrom] + bin_size - 1) // bin_size
    return [(i * bin_size, min((i + 1) * bin_size, CHROM_SIZES[chrom]))
            for i in range(n)]


def _mapq_index(n_fields: int) -> int:
    """Index of the MAPQ column for either FinaleDB schema.

    5 columns: chrom start end mapq strand
    6 columns: chrom start end name mapq strand  (bedtools bamtobed -bedpe)
    """
    return 3 if n_fields == 5 else 4


# Chromosome order and per-chromosome bin counts (identical grid to the
# sibling repo's extract_delfi.py, so features are directly comparable).
_CHROMS = list(CHROM_SIZES)
_CHROM_CODE = {c: i for i, c in enumerate(_CHROMS)}
_N_BINS_100KB = np.array([len(_bin_map(c, BIN_100KB)) for c in _CHROMS])
_N_BINS_5MB = np.array([len(_bin_map(c, BIN_5MB)) for c in _CHROMS])
_OFFSET_100KB = np.concatenate([[0], np.cumsum(_N_BINS_100KB)[:-1]])
_OFFSET_5MB = np.concatenate([[0], np.cumsum(_N_BINS_5MB)[:-1]])
_TOTAL_100KB = int(_N_BINS_100KB.sum())
_TOTAL_5MB = int(_N_BINS_5MB.sum())


def read_fragments(path: str, mapq_threshold: int = 0):
    """Read a FinaleDB fragment file into arrays.

    Returns (chrom_codes, start, end, mapq) where ``chrom_codes`` is the
    index of the chromosome in ``_CHROMS`` (-1 for contigs outside the
    canonical set, which are dropped by the callers).

    Vectorised with pandas: the training-cohort files hold ~37 M
    fragments each, and a per-line Python loop takes minutes per sample.
    """
    header = pd.read_csv(path, sep="\t", header=None, nrows=2,
                         compression="gzip", engine="c")
    n_fields = header.shape[1]
    usecols = [0, 1, 2, _mapq_index(n_fields)]
    df = pd.read_csv(path, sep="\t", header=None, usecols=usecols,
                     compression="gzip", engine="c",
                     names=["chrom", "start", "end", "mapq"])
    chrom = df["chrom"].astype(str).to_numpy()
    start = df["start"].to_numpy(dtype=np.int64)
    end = df["end"].to_numpy(dtype=np.int64)
    mapq = df["mapq"].to_numpy(dtype=np.int64)
    keep = mapq >= mapq_threshold
    codes = np.array([_CHROM_CODE.get(c, -1) for c in chrom], dtype=np.int64)
    return codes[keep], start[keep], end[keep], mapq[keep]


def extract_delfi(path: str) -> dict:
    """DELFI 5 Mb / 100 kb regional fragmentation features.

    Returns a dict with the 4 arrays the cross-study benchmark consumes:
    ``5mb_ratio``, ``5mb_coverage`` (median-normalised), ``100kb_ratio``
    and ``100kb_counts`` -- the same quantities the sibling repo writes
    to ``<sample>.delfi_*.npy``.
    """
    codes, start, end, _mapq = read_fragments(path)
    flen = end - start

    keep = (codes >= 0) & (flen >= MIN_LEN) & (flen <= LONG_MAX)
    codes = codes[keep]
    start = start[keep]
    flen = flen[keep]

    is_short = flen <= SHORT_MAX
    total_short = int(is_short.sum())
    total_long = int((~is_short).sum())

    # 100 kb bins
    bin100 = start // BIN_100KB
    valid100 = bin100 < _N_BINS_100KB[codes]
    flat100 = _OFFSET_100KB[codes] + bin100
    short_100 = np.bincount(flat100[valid100 & is_short],
                            minlength=_TOTAL_100KB).astype(np.int64)
    long_100 = np.bincount(flat100[valid100 & ~is_short],
                           minlength=_TOTAL_100KB).astype(np.int64)

    # 5 Mb bins
    bin5 = start // BIN_5MB
    valid5 = bin5 < _N_BINS_5MB[codes]
    flat5 = _OFFSET_5MB[codes] + bin5
    short_5 = np.bincount(flat5[valid5 & is_short],
                           minlength=_TOTAL_5MB).astype(np.int64)
    long_5 = np.bincount(flat5[valid5 & ~is_short],
                          minlength=_TOTAL_5MB).astype(np.int64)

    ratios_100kb = short_100 / np.maximum(long_100, 1)
    counts_100kb = (short_100 + long_100).astype(np.int64)
    ratio_5mb = short_5 / np.maximum(long_5, 1)
    cover_5mb = (short_5 + long_5).astype(np.int64)
    med = float(np.median(cover_5mb[cover_5mb > 0]))
    cover_norm = cover_5mb.astype(float) / max(med, 1.0)

    return {
        "5mb_ratio": ratio_5mb,
        "5mb_coverage": cover_norm,
        "100kb_ratio": ratios_100kb,
        "100kb_counts": counts_100kb,
        "genome_short_long_ratio": float(total_short / max(total_long, 1)),
        "mean_short_frac": float(total_short / max(total_short + total_long, 1)),
        "n_windows_100kb": int(len(ratios_100kb)),
    }


def extract_fsd(path: str, mapq_threshold: int = 0) -> dict:
    """Fragment-size distribution histogram (196 x 5 bp bins from 20 bp).

    Returns the same ``size_bins`` mapping the sibling repo writes to
    ``<sample>.fsd.json`` (normalised frequencies), so the benchmark's
    ``load5`` can read the result unchanged.
    """
    _codes, start, end, _mapq = read_fragments(path, mapq_threshold)
    arr = (end - start)
    arr = arr[(arr >= FSD_BIN_START) & (arr <= FSD_BIN_END)]
    n = int(arr.size)
    if n == 0:
        raise ValueError("no fragments in input")
    hist, edges = np.histogram(
        arr, bins=range(FSD_BIN_START, FSD_BIN_END + 1, FSD_BIN_STRIDE))
    labels = [f"{edges[i]}-{edges[i + 1]}" for i in range(len(edges) - 1)]
    size_bins = {labels[i]: float(hist[i] / n) for i in range(len(hist))}
    short = int(((arr >= 100) & (arr <= 150)).sum())
    long_ = int(((arr >= 150) & (arr <= 220)).sum())
    return {
        "fragment_count": n,
        "median_length": float(np.median(arr)),
        "mean_length": float(np.mean(arr)),
        "mode_length": float(edges[int(np.argmax(hist))] + 2.5),
        "p10": float(np.percentile(arr, 10)),
        "p25": float(np.percentile(arr, 25)),
        "p75": float(np.percentile(arr, 75)),
        "p90": float(np.percentile(arr, 90)),
        "short_fraction_100_150": float(short / n),
        "long_fraction_150_220": float(long_ / n),
        "short_long_ratio": float(short / max(long_, 1)),
        "size_bins": size_bins,
    }


def extract_sample(frag_path: str, out_dir: str, sample: str) -> dict:
    """Extract all 5 benchmark channels for one sample and write them.

    Writes ``<sample>.delfi_5mb_ratio.npy``, ``.delfi_5mb_coverage.npy``,
    ``.delfi_100kb_ratio.npy``, ``.delfi_100kb_counts.npy``,
    ``.delfi.json`` and ``.fsd.json`` -- the exact file set and naming
    the cross-study benchmark's loader expects.
    """
    os.makedirs(out_dir, exist_ok=True)
    delfi = extract_delfi(frag_path)
    fsd = extract_fsd(frag_path)

    np.save(os.path.join(out_dir, f"{sample}.delfi_5mb_ratio.npy"),
            delfi["5mb_ratio"])
    np.save(os.path.join(out_dir, f"{sample}.delfi_5mb_coverage.npy"),
            delfi["5mb_coverage"])
    np.save(os.path.join(out_dir, f"{sample}.delfi_100kb_ratio.npy"),
            delfi["100kb_ratio"])
    np.save(os.path.join(out_dir, f"{sample}.delfi_100kb_counts.npy"),
            delfi["100kb_counts"])

    with open(os.path.join(out_dir, f"{sample}.delfi.json"), "w") as f:
        json.dump({
            "sample": sample,
            "genome_short_long_ratio": delfi["genome_short_long_ratio"],
            "mean_short_frac": delfi["mean_short_frac"],
            "n_windows_100kb": delfi["n_windows_100kb"],
        }, f)

    with open(os.path.join(out_dir, f"{sample}.fsd.json"), "w") as f:
        json.dump({"sample": sample, "mode": "frag", **fsd}, f)

    return {
        "sample": sample,
        "genome_short_long_ratio": delfi["genome_short_long_ratio"],
        "mean_short_frac": delfi["mean_short_frac"],
        "fsd_fragment_count": fsd["fragment_count"],
        "fsd_short_fraction": fsd["short_fraction_100_150"],
        "fsd_short_long_ratio": fsd["short_long_ratio"],
    }
