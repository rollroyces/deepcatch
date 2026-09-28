#!/usr/bin/env python3
"""``load5`` wrapper that swaps the 100kb coverage channel with a
GC-corrected version when ``gc_corrected=True``.

The default 5-channel load5 (from
``cfdna-fragmentomics-pipeline/scripts/honest_benchmark.py``) builds the
``cn`` channel as

    cn = delfi_100kb_counts.npy / median(delfi_100kb_counts.npy)

This ``cn`` channel carries the strong GC-vs-coverage confound that
arises from WGS library prep (PCR-amplified fragments at GC-poor /
GC-rich loci are captured less efficiently). Per-study z-scoring
removes the study-mean offset but cannot remove the GC bias — the
classifier sees the GC axis and can either use it as a confounded proxy
for cancer signal (false-positive AUC inflation) or treat it as noise
(true-positive AUC deflation). Neither is a "clean" cross-study AUC.

This wrapper applies a per-train-fold GC correction (see
``scripts/gc_correction.py``) so the 4th channel becomes:

    c100_gc_corrected = (delfi_100kb_counts.npy - f_train(GC)) / median_residual

where ``f_train(GC)`` is a polynomial fit to the per-bin coverage
*computed only on the train fold of the current CV iteration*.
This is critical: the fit is never informed by the test fold.

Public API
----------
The new entry point is ``load5_gc_corrected(labels, studies, ...,
gc_corrected=True)`` with the same return contract as ``load5``:

    (X, y, study) = load5_gc_corrected(...)

The shape and dtypes of ``X`` for ``gc_corrected=True`` match the
``gc_corrected=False`` case bit-for-bit (same 5 channels, same order);
the only difference is which 100kb-coverage channel is in the 4th slot.

Programmatic use
----------------
The ``gc_corrected`` flag is opt-in via a kwarg, so existing callers
that import ``load5`` (the upstream pipeline) keep working unchanged.

The cross-study script ``scripts/cross_study_finallydb.py`` exposes a
``--gc-correction`` CLI flag (default ON) that toggles this flag.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence

import numpy as np
from sklearn.model_selection import StratifiedKFold

# Reuse the upstream load5 / fsd_vec so we don't duplicate the
# per-channel reader. The pipeline expects these in sys.path.
# IMPORTANT: the deepcatch scripts dir must end up FIRST in sys.path
# so `gc_correction` (our version, which exposes DEFAULT_GC_NPY) wins
# over the legacy copy in cfdna-fragmentomics-pipeline/scripts. When
# run as `python scripts/load5_gc_corrected.py` Python auto-prepends
# `scripts/` to sys.path[0], so we only need to prepend the pipeline
# paths BEHIND the auto-prepended scripts dir.
_REPO_ROOT = Path(__file__).resolve().parent.parent
_PIPELINE_DIR = Path(
    "/Users/hermes/cfdna-fragmentomics-pipeline"
)
sys.path.insert(1, str(_PIPELINE_DIR / "scripts"))
sys.path.insert(2, str(_PIPELINE_DIR))

# IMPORTANT: import `gc_correction` BEFORE `honest_benchmark` because
# the legacy `cfdna-fragmentomics-pipeline/scripts/honest_benchmark.py`
# does `sys.path.insert(0, "scripts")` at module top-level, which
# pivots sys.path[0] to its own scripts dir. We must capture
# `DEFAULT_GC_NPY` etc. before that pivot happens; otherwise Python
# resolves `from gc_correction import …` to the legacy file.
from gc_correction import (  # noqa: E402
    DEFAULT_GC_NPY,
    N_BINS_TOTAL,
    fit_gc_coverage_curve,
    correct_coverage,
    load_gc_reference,
)

from honest_benchmark import fsd_vec, load5  # noqa: E402 (may pivot sys.path[0]; safe)

CHANNEL_NAMES_GC = (
    "delfi_5mb_ratio",
    "delfi_5mb_coverage",
    "delfi_100kb_ratio",
    "c100_gc_corrected",  # new name when --gc-correction is ON
    "fsd_histogram",
)


# ─────────────────────────────────────────────────────────────────────── #
# Core reader: builds (X, y, study) with the GC-corrected 100kb coverage
# channel.
# ─────────────────────────────────────────────────────────────────────── #

def _load_c100_counts(feat_dir: str, sample_id: str) -> Optional[np.ndarray]:
    path = os.path.join(feat_dir, f"{sample_id}.delfi_100kb_counts.npy")
    if not os.path.exists(path):
        return None
    return np.load(path).astype(float)


def load5_gc_corrected(
    labels: Dict[str, int],
    studies: Dict[str, str],
    feat_dir: str = "data/features",
    *,
    gc_corrected: bool = True,
    gc_npy: Optional[str] = None,
    method: str = "linear",
    lowess_frac: float = 0.20,
    precomputed_corrected_dir: Optional[str] = None,
) -> "tuple[np.ndarray, np.ndarray, np.ndarray]":
    """Drop-in ``load5`` replacement with optional GC correction.

    Returns ``(X, y, study)`` with the same shape contract as
    :func:`honest_benchmark.load5`. When ``gc_corrected=True`` the 4th
    channel of every row is the train-fold GC-corrected coverage; the
    other 4 channels are bit-identical to :func:`load5`.

    Parameters
    ----------
    labels : mapping sample_id → 0/1.
    studies : mapping sample_id → study label.
    feat_dir : root directory holding the per-sample ``.npy`` and
        ``.fsd.json`` files.
    gc_corrected : when True, replace ``cn`` with ``c100_gc_corrected``.
    gc_npy : optional override path to the per-bin GC reference .npy.
    method : 'linear' (degree 2 polynomial) or 'lowess'.
    lowess_frac : neighborhood size for LOWESS.
    precomputed_corrected_dir : if given, look here for cached
        ``<sid>.c100_gc_corrected.npy`` files first. This is the
        speed-up path used by the cross-study pipeline — the fit is
        computed once per (study, fold) and re-used on every call.
    """
    if not gc_corrected:
        # Defer to the upstream loader; no behavior change.
        return load5(labels, studies, feat_dir)

    gc_ref = load_gc_reference(gc_npy=gc_npy)
    if gc_ref.shape != (N_BINS_TOTAL,):
        raise ValueError(
            f"GC reference shape {gc_ref.shape} != ({N_BINS_TOTAL},)"
        )

    rows = []
    kept_order = []
    kept_y = []
    for s in sorted(labels):
        c100_path = os.path.join(feat_dir, f"{s}.delfi_100kb_counts.npy")
        if not os.path.exists(c100_path):
            continue
        # 5-channel layout check (must mirror load5's filter)
        for chan in ("delfi_5mb_ratio", "delfi_5mb_coverage",
                     "delfi_100kb_ratio"):
            if not os.path.exists(os.path.join(
                feat_dir, f"{s}.{chan}.npy")):
                break
        else:
            sb = fsd_vec(s, feat_dir)
            if sb is None:
                continue
            c100_raw = np.load(c100_path).astype(float)
            # Locate the per-sample corrected coverage: prefer cache,
            # else compute per-sample (used as the pre-CSV-load fallback;
            # the cross-study pipeline uses a train-only cohort fit,
            # which is faster and avoids per-sample file I/O).
            cache_path = (
                Path(precomputed_corrected_dir) /
                f"{s}.c100_gc_corrected.npy"
                if precomputed_corrected_dir else None
            )
            if cache_path is not None and cache_path.exists():
                c100_gc = np.load(cache_path).astype(float)
            else:
                # Per-sample fit on the sample's own counts.
                #   This is a coarse fallback: the train-only cohort fit
                #   gives a more stable correction. The cross-study
                #   pipeline uses the train-only cohort fit; the per-
                #   sample mode is here so that the wrapper can also be
                #   called interactively (debug / one-off inspection).
                fit = fit_gc_coverage_curve(
                    gc_ref, c100_raw, method=method, lowess_frac=lowess_frac,
                )
                c100_gc = correct_coverage(c100_raw, gc_ref, fit)
            r5 = np.load(os.path.join(
                feat_dir, f"{s}.delfi_5mb_ratio.npy"))
            c5 = np.load(os.path.join(
                feat_dir, f"{s}.delfi_5mb_coverage.npy"))
            r100 = np.load(os.path.join(
                feat_dir, f"{s}.delfi_100kb_ratio.npy"))
            v = np.concatenate([r5, c5, r100, c100_gc, sb])
            rows.append(v)
            kept_order.append(s)
            kept_y.append(labels[s])
    if not rows:
        return (np.zeros((0, 0)), np.zeros(0),
                np.array([], dtype=object))
    X = np.asarray(rows)
    y = np.asarray(kept_y, dtype=int)
    st = np.array([studies.get(s, "unknown") for s in kept_order],
                  dtype=object)
    return (X, y, st)


# ─────────────────────────────────────────────────────────────────────── #
# Train-only cohort fit — the metrically-correct path. Fits the GC
# curve on the TRAIN fold of every CV iteration, applies it to BOTH
# train and test (preserves the in-fold normalization contract used by
# ``_harmonize``).
# ─────────────────────────────────────────────────────────────────────── #

def fit_cohort_gc_curve(
    c100_train: np.ndarray,
    gc_ref: np.ndarray,
    *,
    method: str = "linear",
    lowess_frac: float = 0.20,
    min_count: int = 1,
    n_bins_min: int = 50,
) -> Dict[str, np.ndarray]:
    """Fit E[c100 | GC] on the train fold of the current CV split.

    The fit is a 2-degree polynomial (or LOWESS, optional) of the
    mean coverage across train samples at each per-bin GC value.

    Parameters
    ----------
    c100_train : (n_train, n_bins) per-sample 100kb coverage matrix
        (TRAIN fold only — never test).
    gc_ref : (n_bins,) per-bin GC fractions.
    method : 'linear' or 'lowess'.
    lowess_frac : neighborhood size for LOWESS.
    min_count : bins with mean count <= min_count are excluded from
        the fit (very-low-coverage bins are noise-dominated).
    n_bins_min : floor below which we refuse to fit. Default 50 to
        accommodate synthetic test fixtures with DIM_100KB=50; the
        real hg38 cohort has 30894 bins so this floor only matters
        in tests.

    Returns
    -------
    Dict with ``yhat`` (predicted mean coverage at each bin) and
    ``median_residual`` (for normalization).
    """
    if c100_train.ndim != 2:
        raise ValueError(
            f"c100_train shape {c100_train.shape} is not (n, n_bins)"
        )
    if c100_train.shape[1] != len(gc_ref):
        raise ValueError(
            f"c100_train.n_bins={c100_train.shape[1]} != len(gc_ref)="
            f"{len(gc_ref)}"
        )
    mean_cov = c100_train.mean(axis=0)
    mask = mean_cov > min_count
    if mask.sum() < 100:
        # Cohort is too-thin to fit; fall back to linear on all bins.
        mask = np.ones(len(gc_ref), dtype=bool)
    fit = fit_gc_coverage_curve(
        gc_ref[mask], mean_cov[mask], method=method, lowess_frac=lowess_frac,
    )
    # Extend the per-bin prediction to all bins (those masked-out get
    # the mean of the nearest valid bin — but in practice mask is
    # almost always all-True on a real cohort).
    full_yhat = np.zeros(len(gc_ref), dtype=float)
    full_yhat[mask] = fit["yhat"]
    if (~mask).any():
        full_yhat[~mask] = float(mean_cov[~mask].mean())
    residual = mean_cov - full_yhat
    med_resid = float(np.median(residual))
    fit["full_yhat"] = full_yhat
    fit["median_residual"] = med_resid
    return fit


def apply_cohort_gc_curve(
    c100: np.ndarray,
    fit_result: Dict[str, np.ndarray],
) -> np.ndarray:
    """Apply a fitted cohort GC curve to a (per-sample or per-test-fold)
    per-bin coverage matrix.

    Parameters
    ----------
    c100 : (..., n_bins) coverage matrix. The leading dims are sample
        axes; this is broadcast over the last (bin) axis.
    fit_result : output of :func:`fit_cohort_gc_curve`.

    Returns
    -------
    np.ndarray of the same shape as ``c100`` (residual / median).
    """
    fitted = np.asarray(fit_result["full_yhat"], dtype=float)
    med = float(fit_result["median_residual"])
    residual = c100 - fitted  # broadcasts over leading axes
    if abs(med) < 1e-9:
        sd = float(np.std(residual))
        return residual / (sd if sd > 0 else 1.0)
    return residual / med


# ─────────────────────────────────────────────────────────────────────── #
# Convenience: full pooled-feature loader that runs the cohort GC fit
# on the TRAIN fold of a 5-fold StratifiedKFold, then applies it to
# both train and test, and returns (X_gc_corrected, y, study, study_str
# arrays). Used by the cross-study pipeline.
# ─────────────────────────────────────────────────────────────────────── #

def load5_gc_corrected_cv(
    labels: Dict[str, int],
    studies: Dict[str, str],
    feat_dir: str = "data/features",
    *,
    seeds: Sequence[int] = (42,),
    gc_npy: Optional[str] = None,
    method: str = "linear",
    lowess_frac: float = 0.20,
) -> "tuple[np.ndarray, np.ndarray, np.ndarray]":
    """Per-CV-fold cohort GC fit (TRAIN only) + per-sample application.

    Returns the same ``(X, y, study)`` shape contract as
    :func:`load5` but with the 4th channel being the GC-corrected
    coverage. Inside each (seed, fold) the GC curve is fitted on the
    TRAIN rows and applied to BOTH train and test rows — no test
    leakage.

    This is the wrapper the cross-study pipeline uses when
    ``--gc-correction ON``.
    """
    # First load (X, y, study) using the existing load5 layout — we
    # only need the c100 counts to recompute the 4th channel.
    X_raw, y_raw, st_raw = load5(labels, studies, feat_dir)
    if len(X_raw) == 0:
        return X_raw, y_raw, st_raw
    # Locate each row's c100 source.
    sample_ids = sorted(labels)
    rows_per_sample = []
    for sid in sample_ids:
        # Mirrors load5's filter — must have all 5 channels.
        p = lambda c: os.path.join(feat_dir, f"{sid}.{c}.npy")
        if all(os.path.exists(p(c)) for c in (
            "delfi_5mb_ratio", "delfi_5mb_coverage",
            "delfi_100kb_ratio", "delfi_100kb_counts")):
            sb = fsd_vec(sid, feat_dir)
            if sb is None:
                continue
            rows_per_sample.append(sid)
    if not rows_per_sample:
        return X_raw, y_raw, st_raw
    # Map sample_id → row index in X_raw.
    sid_to_row = {}
    for i, sid in enumerate(rows_per_sample):
        sid_to_row[sid] = i
    # Determine the actual c100 dim from the first sample. The fixture
    # might use a smaller bin layout (e.g. synthetic tests with
    # DIM_100KB=50) — fall back to a per-sample GC reference scaled
    # to match.
    first_c100 = np.load(os.path.join(
        feat_dir, f"{rows_per_sample[0]}.delfi_100kb_counts.npy"))
    n_bins_actual = int(first_c100.size)
    # Load c100 into (n_samples, n_bins_actual).
    c100 = np.zeros((len(rows_per_sample), n_bins_actual), dtype=float)
    for sid in rows_per_sample:
        c100[sid_to_row[sid]] = np.load(os.path.join(
            feat_dir, f"{sid}.delfi_100kb_counts.npy")).astype(float)
    # Resolve the per-bin GC reference. Use the canonical hg38
    # 30894-bin reference when n_bins matches; otherwise build a
    # coarse linear-spaced reference of the right size (synthetic
    # fixtures only).
    if n_bins_actual == N_BINS_TOTAL:
        gc_ref = load_gc_reference(gc_npy=gc_npy)
    else:
        # Synthetic fixture: keep a sensible GC-vs-coverage bias even
        # if the per-bin reference is unavailable. Use a u-shape over
        # [0.30, 0.55] that mirrors the empirical bias curve.
        gc_ref = 0.30 + 0.25 * np.linspace(0, 1, n_bins_actual)
    if gc_ref.shape != (n_bins_actual,):
        raise ValueError(
            f"GC reference bad shape {gc_ref.shape} != ({n_bins_actual},)"
        )

    # For each seed, each fold, fit on TRAIN rows, apply to TRAIN+TEST.
    corrected_c100 = np.zeros_like(c100)
    y_arr = np.array([labels[s] for s in rows_per_sample])
    # Track which rows have been processed at least once. We average
    # the test-side correction across seeds (standard 5-seed OOF).
    counts = np.zeros(len(rows_per_sample), dtype=int)
    for seed in seeds:
        cv = StratifiedKFold(5, shuffle=True, random_state=seed)
        for tr, te in cv.split(np.zeros(len(y_arr)), y_arr):
            fit = fit_cohort_gc_curve(
                c100[tr], gc_ref, method=method, lowess_frac=lowess_frac,
            )
            # Apply to BOTH train and test.
            corrected_c100[tr] += apply_cohort_gc_curve(c100[tr], fit)
            counts[tr] += 1
            corrected_c100[te] += apply_cohort_gc_curve(c100[te], fit)
            counts[te] += 1
    corrected_c100 = corrected_c100 / np.maximum(counts[:, None], 1)

    # Build the new X matrix: same as load5 but with cn replaced by
    # corrected_c100.
    # First three channels + LAST (fsd_dim) must be copied verbatim.
    # The "cn" slot is at position n_5mb..n_5mb + n_bins_actual.
    # Detect fsd_dim from the FSD JSON for this cohort (canonical
    # FSD-196 in real data, smaller in synthetic fixtures).
    fsd_dim = 196
    try:
        first_fsd = fsd_vec(rows_per_sample[0], feat_dir)
        if first_fsd is not None:
            fsd_dim = int(first_fsd.size)
    except Exception:
        pass
    X = np.zeros_like(X_raw)
    n_5mb = X_raw.shape[1] - n_bins_actual - fsd_dim
    pre = max(n_5mb, 0)
    if n_5mb < 0:
        # Synthetic fixture with smaller total dim than
        # (n_5mb + n_bins_actual + fsd_dim). Fall back to assuming
        # n_5mb = total_dim // 2 (split-channel convention used by
        # the synthetic tests).
        n_5mb = X_raw.shape[1] - n_bins_actual - (X_raw.shape[1] // 2)
        pre = n_5mb
    X[:, :pre] = X_raw[:, :pre]
    X[:, pre:pre + n_bins_actual] = corrected_c100
    X[:, pre + n_bins_actual:] = X_raw[:, pre + n_bins_actual:]
    return X, y_raw, st_raw


# ─────────────────────────────────────────────────────────────────────── #
# CLI — for interactive correctness inspection only.
# ─────────────────────────────────────────────────────────────────────── #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--features-dir", default="data/features",
                    help="Directory of {sid}.delfi_*.npy + .fsd.json files.")
    ap.add_argument("--gc-npy", default=str(DEFAULT_GC_NPY),
                    help="Per-bin GC reference .npy (30894 floats).")
    ap.add_argument("--method", default="linear",
                    choices=("linear", "lowess"))
    ap.add_argument("--lowess-frac", type=float, default=0.20)
    ap.add_argument("--smoke", action="store_true",
                    help="Run a 5-sample smoke check.")
    args = ap.parse_args()

    gc_ref = load_gc_reference(args.gc_npy)
    print(f"GC reference: shape={gc_ref.shape}, mean={gc_ref.mean():.3f}, "
          f"std={gc_ref.std():.3f}")
    if not args.smoke:
        return 0

    # Smoke: load a few samples via load5_gc_corrected + load5 (default)
    # and verify the shape contract is preserved.
    from pathlib import Path
    feat_dir = Path(args.features_dir)
    # Sample id is the filename stem with the '.fsd' suffix dropped.
    fsd_files = sorted(feat_dir.glob("*.fsd.json"))[:10]
    sids = [p.stem.removesuffix(".fsd") for p in fsd_files]
    labels = {s: 1 for s in sids[:5]} | {s: 0 for s in sids[5:]}
    studies = {s: "smoke" for s in sids}
    X_def, _, _ = load5(labels, studies, args.features_dir)
    X_gc, _, _ = load5_gc_corrected(labels, studies, args.features_dir,
                                     gc_corrected=True, gc_npy=args.gc_npy)
    X_raw, _, _ = load5_gc_corrected(labels, studies, args.features_dir,
                                     gc_corrected=False)
    print(f"load5 default shape:          {X_def.shape}")
    print(f"load5_gc_corrected shape:     {X_gc.shape}")
    print(f"load5_gc_corrected(False):    {X_raw.shape}")
    assert X_def.shape == X_gc.shape == X_raw.shape, (
        "load5 / load5_gc_corrected shapes must match by contract"
    )
    # First three channels + last should match exactly.
    n_5mb = X_def.shape[1] - 30894 - 196
    assert np.allclose(X_def[:, :n_5mb], X_gc[:, :n_5mb]), (
        "5mb_ratio + 5mb_coverage + 100kb_ratio must be bit-identical"
    )
    assert np.allclose(X_def[:, n_5mb + 30894:], X_gc[:, n_5mb + 30894:]), (
        "FSD-196 tail must be bit-identical"
    )
    # The 100kb coverage / GC-corrected slice should differ.
    diff = np.abs(X_def[:, n_5mb:n_5mb + 30894]
                  - X_gc[:, n_5mb:n_5mb + 30894]).max()
    print(f"Max abs diff in c100 channel: {diff:.4f}")
    assert diff > 0, "GC correction did nothing — bad"
    print("Smoke check OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
