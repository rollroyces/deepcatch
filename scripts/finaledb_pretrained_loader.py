#!/usr/bin/env python3
"""Loader for the FinaleDB pretrained checkpoint.

The pretraining script (scripts/pretrain_real_finaledb.py) saves a
flat ``X`` matrix of shape (n_samples, 2256) where the first 83 columns
are the modality summary features and the remaining 2173 are raw DELFI
profile (ratio, coverage, meanlen, motifs, wps). This loader splits
the flat X back into the per-modality dict that
``FoundationDownstream._validate_modalities`` requires.

The layout, in order:

    [frag_basic (4), frag_enhanced (44), cnv (6), sero (4),
     gnn (1), tissue (24), ratio (631), cov (631), meanlen (631),
     motifs (256), wps (24)]
       total = 4+44+6+4+1+24 + 631+631+631+256+24 = 83 + 2173 = 2256

Usage
-----

    from scripts.finaledb_pretrained_loader import (
        load_real_cohort,
        modalities_from_flat_X,
    )

    cohort = load_real_cohort(
        "/Users/hermes/deepcatch/data/finaledb_pretrain_cohort.npz"
    )
    modalities = modalities_from_flat_X(cohort["X"])
    # modalities is dict[str, ndarray] keyed by MODALITY_NAMES
    # each value is (n_samples, expected_dim_i)

    from src.foundation.downstream import FoundationDownstream
    fd = FoundationDownstream(
        pretrained=True,
        checkpoint_path="/Users/hermes/deepcatch/checkpoints/foundation_pretrained_finaledb.pt",
    )
    proba = fd.predict_proba(modalities)
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np

# Layout constants — must match scripts/pretrain_real_finaledb.py
_SPLIT_OFFSETS = {
    # modality summaries (first 83 dims)
    "frag_basic": (0, 4),
    "frag_enhanced": (4, 48),
    "cnv": (48, 54),
    "sero": (54, 58),
    "gnn": (58, 59),
    "tissue": (59, 83),
    # raw DELFI features (last 2173 dims) — exposed as aux channels
    # but NOT consumed by FoundationDownstream (which only takes the
    # 6 MODALITY_NAMES). Useful for downstream ad-hoc analyses.
    "_ratio_5mb": (83, 714),
    "_cov_5mb": (714, 1345),
    "_meanlen_5mb": (1345, 1976),
    "_motifs": (1976, 2232),
    "_wps_100kb": (2232, 2256),
}

# Sanity check at import time
_TOTAL = sum(end - start for start, end in _SPLIT_OFFSETS.values())
assert _TOTAL == 2256, (
    f"Layout constants sum to {_TOTAL}, expected 2256. "
    "Update both scripts/pretrain_real_finaledb.py and this loader."
)


def load_real_cohort(npz_path: str | Path) -> Dict[str, object]:
    """Load the pretrained-cohort npz. Returns a dict of arrays."""
    data = np.load(npz_path, allow_pickle=True)
    return {
        "X": data["X"],
        "y": data["y"],
        "sample_ids": data["sample_ids"],
        "studies": data["studies"],
        "feature_dim": int(data["feature_dim"]),
        "n_healthy": int(data["n_healthy"]),
        "n_cancer": int(data["n_cancer"]),
    }


def modalities_from_flat_X(
    X: np.ndarray,
    keys: tuple = ("frag_basic", "frag_enhanced", "cnv", "sero", "gnn", "tissue"),
) -> Dict[str, np.ndarray]:
    """Split a flat (n_samples, 2256) matrix into per-modality arrays.

    Parameters
    ----------
    X : (n_samples, 2256) ndarray
        Flat feature matrix from the pretrained-cohort npz.
    keys : tuple of str
        Which modality keys to extract. Default is the 6 keys
        ``FoundationDownstream._validate_modalities`` requires.

    Returns
    -------
    modalities : dict[str, (n_samples, dim_i) ndarray]
        Each value is the slice of X corresponding to that modality.

    Raises
    ------
    ValueError
        If X.shape[1] != 2256.
    """
    if X.ndim != 2 or X.shape[1] != 2256:
        raise ValueError(
            f"Expected X shape (n, 2256), got {X.shape}. "
            "The flat matrix comes from scripts/pretrain_real_finaledb.py."
        )
    out: Dict[str, np.ndarray] = {}
    for k in keys:
        if k not in _SPLIT_OFFSETS:
            raise KeyError(
                f"Unknown modality {k!r}. Valid keys: {list(_SPLIT_OFFSETS)}"
            )
        start, end = _SPLIT_OFFSETS[k]
        out[k] = X[:, start:end].astype(np.float32, copy=False)
    return out


__all__ = ["load_real_cohort", "modalities_from_flat_X"]
