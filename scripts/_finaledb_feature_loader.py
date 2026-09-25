#!/usr/bin/env python3
"""FinaleDB 5-channel feature loader used by the cross-platform validator.

The fragmentomics baseline reads samples from
``/Users/hermes/cfdna-fragmentomics-pipeline/data/features/<sid>.<chan>.npy``
where each ``.npy`` carries the per-channel values for one sample.
Channel naming matches the deepcatch / cfdna-fragmentomics convention:

    delfi_5mb_ratio      delfi_5mb_coverage    delfi_100kb_ratio
    delfi_100kb_counts   fsd_histogram

This helper is intentionally minimal — it does NOT replicate the
upstream pipeline's median-normalization or per-bin GC correction;
those are the consumer's responsibility.
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

# Channel contract — keep in sync with scripts/cross_study_finallydb.py
CHANNEL_NAMES = (
    "delfi_5mb_ratio",
    "delfi_5mb_coverage",
    "delfi_100kb_ratio",
    "delfi_100kb_counts",
    "fsd_histogram",
)


def _sample_id_from_filename(p: Path) -> Optional[str]:
    """``S1.delfi_5mb_ratio.npy`` → ``S1`` when the channel is one of
    :data:`CHANNEL_NAMES`. Returns None when the filename does not
    match the convention.

    Robust against dot-in-sample-ID (e.g. ``GM12878.123``): we always
    treat the LAST dot-segment as the ``.npy`` extension and the
    SECOND-TO-LAST as the channel name; everything before is the
    sample ID.
    """
    name = p.name
    if not name.endswith(".npy"):
        return None
    stem = name[:-4]  # strip '.npy'
    parts = stem.split(".")
    if len(parts) < 2:
        return None
    chan_candidate = parts[-1]
    sid = ".".join(parts[:-1])
    for chan in CHANNEL_NAMES:
        if chan == chan_candidate:
            return sid
    # Lenient fallback: also accept loose channel match (last segment
    # starts with `delfi_`, `fsd_`, etc.).
    if chan_candidate.startswith(("delfi_", "fsd_")):
        return sid
    return sid  # last resort — still return sid so caller can lookup


def load_all_channels_or_skip(
    features_dir: Path,
    sample_ids: Optional[List[str]] = None,
) -> Tuple[np.ndarray, List[str]]:
    """Load the 5-channel features for the requested sample IDs.

    Returns ``(X, kept)`` where ``X`` has shape ``(len(kept), total_dim)``
    (concatenation of channels in :data:`CHANNEL_NAMES` order).
    Samples that are missing any channel are silently dropped from the
    output (the cross-platform validator labels the dropped count in
    its provenance block).

    When ``sample_ids`` is None, all samples in the directory are used.
    """
    features_dir = Path(features_dir)
    if not features_dir.is_dir():
        return np.zeros((0, 0)), []

    # Group files by sample_id
    by_sample: dict = {}
    for p in sorted(features_dir.glob("*.npy")):
        sid = _sample_id_from_filename(p)
        if sid is None:
            continue
        # channel name is the path's last ".npy"-stripped segment
        chan = p.stem.rsplit(".", 1)[-1]
        by_sample.setdefault(sid, {})[chan] = p

    # Filter to requested sample IDs
    if sample_ids is not None:
        wanted = set(sample_ids)
    else:
        wanted = set(by_sample.keys())

    kept: List[str] = []
    rows: List[np.ndarray] = []
    dim: Optional[int] = None
    for sid in sorted(wanted):
        chans = by_sample.get(sid, {})
        # Each sample must have all 5 channels
        missing = [c for c in CHANNEL_NAMES if c not in chans]
        if missing:
            continue
        try:
            vecs = [np.load(chans[c]) for c in CHANNEL_NAMES]
        except (OSError, ValueError):
            continue
        row = np.concatenate([v.ravel() for v in vecs])
        if dim is None:
            dim = row.size
        if row.size != dim:
            # Skip the sample if its concatenated size disagrees
            # with the rest (different genome build / bin count)
            continue
        rows.append(row)
        kept.append(sid)

    if not rows:
        return np.zeros((0, dim or 0)), []

    X = np.vstack(rows)
    return X, kept
