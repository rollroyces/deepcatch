"""
Tumor-naive fragmentomics channel — DeepCatch adapter for the
cfdna-fragmentomics-pipeline pre-computed artifacts.

This module reads the pipeline's per-sample `.npy` / `.fsd.json` outputs
and assembles the same 5-channel tumor-naive profile used by the
standalone pipeline's classifier (AUC 0.9753 ± 0.002 on 627 cross-study
pan-cancer samples, 5-seed CV).

Channels assembled per sample:
    1. 5Mb short/long ratio (631 bins)            — DELFI ratio profile
    2. 5Mb coverage (median-normalized, 631 bins) — coarse CNA
    3. 100kb short/long ratio (30,894 bins)       — finer DELFI ratio
    4. 100kb coverage (median-normalized)         — finer CNA
    5. FSD size histogram (196 bins, 5bp bins)    — fragment-length shape

Expected pipeline artifacts (per sample):
    {features_dir}/{sample}.delfi_5mb_ratio.npy
    {features_dir}/{sample}.delfi_5mb_coverage.npy
    {features_dir}/{sample}.delfi_100kb_ratio.npy
    {features_dir}/{sample}.delfi_100kb_counts.npy
    {features_dir}/{sample}.fsd.json

The pipeline repo is standalone and not a pip dependency — DeepCatch only
reads the artifacts. This keeps the adapter zero-dependency w.r.t. the
pipeline and prevents tight coupling.

Design choice — why a *reader*, not a *re-implementer*?
    The pipeline's 5-channel profile is the empirically tuned, honestly
    ablated input (motifs and mean-length both within noise). Re-computing
    it inside DeepCatch would silently drift from the upstream result.
    Reading the artifacts makes the two repos automatically stay in sync
    on the data side; the model code is the only thing that lives in
    DeepCatch.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np


# Channel order is the contract between this adapter and any downstream
# classifier. Don't reorder without updating every consumer.
CHANNEL_NAMES = [
    "delfi_5mb_ratio",
    "delfi_5mb_coverage",
    "delfi_100kb_ratio",
    "delfi_100kb_counts",
    "fsd_histogram",
]
CHANNEL_DIMS = {
    "delfi_5mb_ratio": 631,
    "delfi_5mb_coverage": 631,
    "delfi_100kb_ratio": 30894,
    "delfi_100kb_counts": 30894,
    "fsd_histogram": 196,
}
TOTAL_DIM = sum(CHANNEL_DIMS.values())  # 63,246


@dataclass
class TumorNaiveSample:
    """Per-sample tumor-naive feature vector with provenance."""
    sample_id: str
    vector: np.ndarray  # shape (TOTAL_DIM,), float64
    n_frags: Optional[int] = None  # from FSD JSON metadata


def _load_fsd_histogram(fsd_json_path: str) -> np.ndarray:
    """Read the 196-bin (5bp, 20-1000bp) fragment-length histogram.

    The pipeline writes bins as keys like "100-105" → count. We sort
    them numerically so the output is deterministic regardless of
    dict insertion order.
    """
    with open(fsd_json_path) as f:
        d = json.load(f)
    bins = d["size_bins"]
    keys = sorted(bins, key=lambda k: int(k.split("-")[0]))
    return np.asarray([bins[k] for k in keys], dtype=np.float64)


def _load_or_none(path: str) -> Optional[np.ndarray]:
    return np.load(path) if os.path.exists(path) else None


def load_sample(
    sample_id: str,
    features_dir: str,
    channels: Optional[List[str]] = None,
    median_normalize_100kb_coverage: bool = True,
) -> TumorNaiveSample:
    """Load the tumor-naive feature vector for one sample.

    Returns a TumorNaiveSample; raises FileNotFoundError with the missing
    path if any of the requested channels is absent (strict by default
    — silent zeros would mask pipeline bugs).

    Parameters
    ----------
    sample_id : str
        Sample identifier (matches `{sample}.delfi_*.npy` filenames).
    features_dir : str
        Directory containing the pipeline artifacts.
    channels : list of str or None
        Subset of CHANNEL_NAMES to load; None = all 5.
    median_normalize_100kb_coverage : bool
        If True (default), divide the 100kb coverage counts by their
        per-sample median. This matches the pipeline's `load_full_profile`
        and removes sequencing-depth batch effects — without it the
        AUC drops by ~0.008 on cross-study pan-cancer (validated).
    """
    if channels is None:
        channels = CHANNEL_NAMES

    pieces: list[np.ndarray] = []
    n_frags: Optional[int] = None

    for ch in channels:
        if ch in ("delfi_5mb_ratio", "delfi_5mb_coverage",
                 "delfi_100kb_ratio", "delfi_100kb_counts"):
            p = os.path.join(features_dir, f"{sample_id}.{ch}.npy")
            arr = _load_or_none(p)
            if arr is None:
                raise FileNotFoundError(
                    f"Missing pipeline artifact for {sample_id!r} "
                    f"channel {ch!r}: {p}"
                )
            # Per-sample median-normalize 100kb coverage. Matches the
            # pipeline's load_full_profile — without it the cross-study
            # AUC drops ~0.008 (sequencing-depth batch effect).
            if ch == "delfi_100kb_counts" and median_normalize_100kb_coverage:
                med = float(np.median(arr))
                if med > 0:
                    arr = arr.astype(np.float64) / med
            pieces.append(arr)
        elif ch == "fsd_histogram":
            p = os.path.join(features_dir, f"{sample_id}.fsd.json")
            if not os.path.exists(p):
                raise FileNotFoundError(
                    f"Missing pipeline artifact for {sample_id!r} "
                    f"channel 'fsd_histogram': {p}"
                )
            with open(p) as f:
                d = json.load(f)
            n_frags = d.get("fragment_count")
            pieces.append(_load_fsd_histogram(p))
        else:
            raise ValueError(f"Unknown channel: {ch!r}; "
                             f"valid: {CHANNEL_NAMES}")

    return TumorNaiveSample(
        sample_id=sample_id,
        vector=np.concatenate(pieces).astype(np.float64),
        n_frags=n_frags,
    )


def load_cohort(
    sample_ids: List[str],
    features_dir: str,
    channels: Optional[List[str]] = None,
    skip_missing: bool = False,
    median_normalize_100kb_coverage: bool = True,
) -> tuple[np.ndarray, list[str]]:
    """Load a cohort as a stacked matrix.

    Returns (X, order) where X is shape (n_loaded, TOTAL_DIM or sum of
    requested channel dims) and `order[i]` is the sample_id for row i.

    With skip_missing=False (default), any missing artifact is fatal.
    With skip_missing=True, missing samples are silently dropped (and
    reported via shorter `order`).
    """
    rows: list[np.ndarray] = []
    order: list[str] = []
    for s in sample_ids:
        try:
            sample = load_sample(
                s, features_dir,
                channels=channels,
                median_normalize_100kb_coverage=median_normalize_100kb_coverage,
            )
        except FileNotFoundError:
            if not skip_missing:
                raise
            continue
        rows.append(sample.vector)
        order.append(s)
    if not rows:
        return np.empty((0, TOTAL_DIM)), order
    return np.stack(rows, axis=0), order


def load_labels_tsv(labels_tsv: str) -> Dict[str, int]:
    """Read `{sample}\\t{label}` TSV. Returns dict sample_id → int (0/1).

    Recognized labels (case-insensitive): cancer/positive/tumor/y/1/true
    → 1; healthy/control/normal/n/0/false → 0. Anything else raises.
    """
    POS = {"cancer", "positive", "tumor", "y", "1", "true", "case"}
    NEG = {"healthy", "control", "normal", "n", "0", "false"}
    out: Dict[str, int] = {}
    with open(labels_tsv) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            s = parts[0]
            lab = parts[1].strip().lower()
            if lab in POS:
                out[s] = 1
            elif lab in NEG:
                out[s] = 0
            else:
                raise ValueError(
                    f"Unrecognized label {lab!r} for sample {s!r}; "
                    f"expected one of {sorted(POS | NEG)}"
                )
    return out


def expected_dim(channels: Optional[List[str]] = None) -> int:
    """Return the total feature dimension for a channel subset."""
    if channels is None:
        return TOTAL_DIM
    return sum(CHANNEL_DIMS[c] for c in channels)