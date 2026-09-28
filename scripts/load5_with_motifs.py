#!/usr/bin/env python3
"""Extended load5 loader: original 5-channel DELFI features PLUS the 256-dim
4-mer end-motif distribution (Jiang 2020 [12], Snyder 2016 [1], Cristiano 2019
[7]).

The original load5 (cfdna-fragmentomics-pipeline/scripts/honest_benchmark.py
load5) concatenates:
    r5 (5mb_ratio) + c5 (5mb_coverage) + r100 (100kb_ratio)
    + cn (median-normalized 100kb_counts) + sb (FSD-196)
to produce a 63,246-dim feature vector.

This module produces:
    load5_vec + motifs (256,)
The motif file is optional — when missing for a sample, the loader
**silently skips** that sample (and emits a one-line warning) rather than
failing the whole benchmark. This matches the graceful-missing-artifact
behaviour already established by load5 itself (samples without all 5
channels are dropped, not errored).

Motif file: ``{sample_id}.motifs.npy`` (4-mer frequencies, shape (256,),
sum to ~1.0 per sample; see honest_benchmark.py load8 for the convention).

This is the source-of-truth implementation referenced by both:
- ``scripts/cross_study_finallydb.py`` (the cross-study benchmark)
- ``test/test_load5_with_motifs.py`` (the regression suite)

The reviewer-flagged inconsistency (METHODS_PAPER §1/§2.1 claims motifs
are part of the headline 5-channel benchmark, but load5 does not actually
load them) is fixed by wiring this loader through the --include-motifs
flag on the cross-study benchmark. When the flag is OFF (default), load5
behaviour is preserved bit-identical.
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from typing import Optional, Tuple

import numpy as np

# Default features dir matches the upstream cfdna-fragmentomics-pipeline layout.
DEFAULT_FEATURES_DIR = "/Users/hermes/cfdna-fragmentomics-pipeline/data/features"

# Motif file shape and channel-set contract.
MOTIF_DIM = 256  # 4^4 possible 4-mers over {A,C,G,T}

# Path glob pattern for the 5 DELFI channels + the FSD histogram file.
# The 5th 'FSD' entry is a JSON, not a .npy — handled separately.
_DELFI_NPY_FILES = (
    "{s}.delfi_5mb_ratio.npy",
    "{s}.delfi_5mb_coverage.npy",
    "{s}.delfi_100kb_ratio.npy",
    "{s}.delfi_100kb_counts.npy",
)
_FSD_FILE = "{s}.fsd.json"
_MOTIF_FILE = "{s}.motifs.npy"


def _fsd_vec(s: str, feat_dir: str) -> Optional[np.ndarray]:
    """196-bin fragment-length histogram from the FSD JSON. None if absent."""
    p = os.path.join(feat_dir, _FSD_FILE.format(s=s))
    if not os.path.exists(p):
        return None
    with open(p) as f:
        sb = json.load(f)["size_bins"]
    keys = sorted(sb, key=lambda k: int(k.split("-")[0]))
    return np.array([sb[k] for k in keys], dtype=float)


def load5_with_motifs(labels, studies, feat_dir: str = DEFAULT_FEATURES_DIR):
    """Load the 5-channel DELFI vector PLUS 256-dim 4-mer end-motif counts.

    Concatenated feature vector per sample:
        np.concatenate([
            np.load(r5),          # 5mb_ratio
            np.load(c5),          # 5mb_coverage
            np.load(r100),        # 100kb_ratio
            c100 / median(c100),  # median-normalized 100kb_counts (cn)
            fsd_hist(s),          # FSD-196 fragment-size histogram
            np.load(mot),         # 256-dim 4-mer end-motif distribution
        ])

    Missing-artifact policy:
    - Samples missing ANY of the 5 DELFI channels or the FSD JSON are
      silently dropped (matches load5 in honest_benchmark.py).
    - Samples missing ONLY the motif file (motifs.npy) are silently
      dropped from THIS loader (we are loading motifs — so a sample
      without motifs is excluded). The caller may instead use
      ``load5_with_motifs_optional`` if it wants to keep those samples
      with motifs zero-filled; this loader is the strict contract.

    Returns:
        (X, y, study) tuple matching load5's signature.
        X.shape = (n_kept, 63_246 + 256) = (n_kept, 63_502)
    """
    rows, order, y = [], [], []
    n_skipped_no_motif = 0
    for s in sorted(labels):
        r5 = os.path.join(feat_dir, _DELFI_NPY_FILES[0].format(s=s))
        c5 = os.path.join(feat_dir, _DELFI_NPY_FILES[1].format(s=s))
        r100 = os.path.join(feat_dir, _DELFI_NPY_FILES[2].format(s=s))
        c100 = os.path.join(feat_dir, _DELFI_NPY_FILES[3].format(s=s))
        mot = os.path.join(feat_dir, _MOTIF_FILE.format(s=s))
        if not all(os.path.exists(p) for p in (r5, c5, r100, c100, mot)):
            if not os.path.exists(mot):
                n_skipped_no_motif += 1
            continue
        sb = _fsd_vec(s, feat_dir)
        if sb is None:
            continue
        cn = np.load(c100) / np.median(np.load(c100))
        v = np.concatenate([
            np.load(r5), np.load(c5), np.load(r100), cn, sb,
            np.load(mot),
        ])
        rows.append(v)
        order.append(s)
        y.append(labels[s])
    if n_skipped_no_motif:
        warnings.warn(
            f"load5_with_motifs: dropped {n_skipped_no_motif} samples whose "
            f"motifs.npy was missing in {feat_dir}. This is expected on the "
            f"local FinaleDB cache (only 98 of 627 samples have motifs). "
            f"Use load5 (without motifs) or load5_with_motifs_optional "
            f"(motifs zero-filled) if you need the full cohort.",
            stacklevel=2,
        )
    return (
        np.asarray(rows),
        np.asarray(y),
        np.array([studies.get(s, "6") for s in order]),
    )


def load5_with_motifs_optional(
    labels, studies, feat_dir: str = DEFAULT_FEATURES_DIR,
):
    """Variant of :func:`load5_with_motifs` that KEEPS samples missing the
    motif file by zero-filling the motif block instead of dropping them.

    Use this for head-to-head comparisons on the SAME cohort: load5 returns
    n=627 (the cross-study cohort), load5_with_motifs drops to ~n=98 (only
    the motifs-cache subset), and load5_with_motifs_optional keeps n=627
    with zero-filled motif blocks for the no-motif samples.

    Zero-filling motif frequencies is NOT biologically correct (a sample
    with zero motif frequencies is a flat distribution, not "no data"),
    but it provides a head-to-head comparison on the full cohort so the
    AUC delta reflects only the new feature, not the cohort filter.
    """
    rows, order, y = [], [], []
    for s in sorted(labels):
        r5 = os.path.join(feat_dir, _DELFI_NPY_FILES[0].format(s=s))
        c5 = os.path.join(feat_dir, _DELFI_NPY_FILES[1].format(s=s))
        r100 = os.path.join(feat_dir, _DELFI_NPY_FILES[2].format(s=s))
        c100 = os.path.join(feat_dir, _DELFI_NPY_FILES[3].format(s=s))
        if not all(os.path.exists(p) for p in (r5, c5, r100, c100)):
            continue
        sb = _fsd_vec(s, feat_dir)
        if sb is None:
            continue
        cn = np.load(c100) / np.median(np.load(c100))
        mot_path = os.path.join(feat_dir, _MOTIF_FILE.format(s=s))
        if os.path.exists(mot_path):
            mot_vec = np.load(mot_path)
        else:
            mot_vec = np.zeros(MOTIF_DIM, dtype=float)
        v = np.concatenate([
            np.load(r5), np.load(c5), np.load(r100), cn, sb, mot_vec,
        ])
        rows.append(v)
        order.append(s)
        y.append(labels[s])
    return (
        np.asarray(rows),
        np.asarray(y),
        np.array([studies.get(s, "6") for s in order]),
    )


def load5_only(labels, studies, feat_dir: str = DEFAULT_FEATURES_DIR):
    """Reproduce the upstream load5 contract exactly (63,246-dim, no motifs).

    Re-implemented here to avoid the cross-pipeline ``from honest_benchmark
    import load5`` import (which would couple deepcatch to
    cfdna-fragmentomics-pipeline at import time). The implementation is
    identical to honest_benchmark.py:load5 line-for-line so bit-for-bit
    identity is preserved on the same cohort.
    """
    rows, order, y = [], [], []
    for s in sorted(labels):
        r5 = os.path.join(feat_dir, _DELFI_NPY_FILES[0].format(s=s))
        c5 = os.path.join(feat_dir, _DELFI_NPY_FILES[1].format(s=s))
        r100 = os.path.join(feat_dir, _DELFI_NPY_FILES[2].format(s=s))
        c100 = os.path.join(feat_dir, _DELFI_NPY_FILES[3].format(s=s))
        if not all(os.path.exists(p) for p in (r5, c5, r100, c100)):
            continue
        sb = _fsd_vec(s, feat_dir)
        if sb is None:
            continue
        cn = np.load(c100) / np.median(np.load(c100))
        v = np.concatenate([np.load(r5), np.load(c5), np.load(r100), cn, sb])
        rows.append(v)
        order.append(s)
        y.append(labels[s])
    return (
        np.asarray(rows),
        np.asarray(y),
        np.array([studies.get(s, "6") for s in order]),
    )


# --------------------------------------------------------------------------- #
# CLI: sanity-check the loader on the local cache and print the dimension delta
# --------------------------------------------------------------------------- #
def _main():
    """Standalone CLI: report the dimension delta on the local features cache.

    Run as:
        python scripts/load5_with_motifs.py --features-dir <DIR>
    """
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--features-dir", default=DEFAULT_FEATURES_DIR,
        help=f"Directory of {{sample}}.{{chan}}.npy + .fsd.json + .motifs.npy "
             f"(default: {DEFAULT_FEATURES_DIR})",
    )
    ap.add_argument(
        "--labels-tsv", default=None,
        help="Optional labels TSV (sample<TAB>label<TAB>study). "
             "If omitted, scans the features dir and treats any sample "
             "with a motifs.npy as 'present'.",
    )
    args = ap.parse_args()

    if args.labels_tsv:
        labels, studies = {}, {}
        with open(args.labels_tsv) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 2 or p[0] == "sample":
                    continue
                labels[p[0]] = 1 if p[1] == "cancer" else 0
                studies[p[0]] = p[2] if len(p) >= 3 else "unknown"
        X5, y5, st5 = load5_only(labels, studies, args.features_dir)
        Xm, ym, stm = load5_with_motifs(labels, studies, args.features_dir)
        Xo, _, sto = load5_with_motifs_optional(labels, studies, args.features_dir)
        print(f"load5_only              : X.shape={X5.shape}, "
              f"n_cancer={int((y5==1).sum())}, n_healthy={int((y5==0).sum())}")
        print(f"load5_with_motifs       : X.shape={Xm.shape}, "
              f"n_cancer={int((ym==1).sum())}, n_healthy={int((ym==0).sum())}")
        print(f"load5_with_motifs_optional: X.shape={Xo.shape}")
        # Validation: motif block is exactly the last 256 dims of Xm.
        if Xm.shape[0] > 0:
            motif_sum = float(Xm[0, -MOTIF_DIM:].sum())
            print(f"sample 0 motif block sum: {motif_sum:.6f} (should be ~1.0)")
    else:
        # No labels file — just count motif vs DELFI file presence.
        fdir = args.features_dir
        n_delfi, n_motif, n_both = 0, 0, 0
        sample_ids = set()
        for fname in os.listdir(fdir):
            if fname.endswith(".motifs.npy"):
                n_motif += 1
                sample_ids.add(fname[:-len(".motifs.npy")])
        for fname in os.listdir(fdir):
            if fname.endswith(".delfi_5mb_ratio.npy"):
                n_delfi += 1
        for s in sample_ids:
            if all(os.path.exists(os.path.join(fdir, t.format(s=s)))
                   for t in _DELFI_NPY_FILES) and \
               os.path.exists(os.path.join(fdir, _FSD_FILE.format(s=s))):
                n_both += 1
        print(f"features dir          : {fdir}")
        print(f"samples with delfi    : {n_delfi}")
        print(f"samples with motifs   : {n_motif}")
        print(f"samples with both     : {n_both}")


if __name__ == "__main__":
    sys.exit(_main() or 0)
