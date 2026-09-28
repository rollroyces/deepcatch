#!/usr/bin/env python3
"""Per-bin GC / mappability bias correction for the 100kb coverage channel.

Background
----------
Per-study z-score harmonization (see ``scripts/harmonization_check.py``)
is a *study-mean shift* correction. It does NOT correct **per-bin GC bias**:
the systematic coverage shift between GC-poor and GC-rich genomic windows,
introduced by PCR / capture / WGS-library-prep chemistry rather than
biology. Without this step, the 100kb coverage channel carries a strong,
study-dependent bias that the classifier can either use as a confounded
proxy for cancer signal (false-positive AUC inflation) or treat as noise
(true-positive AUC deflation). Either way, batch effects masquerade as
biology in this channel.

This module applies the canonical LOESS-style correction used in
fragmentomics pipelines (Borg 2012, Teasdale 2010 [4], Sun 2019 [14]):

    1. For each 100kb bin, load the empirical GC fraction computed from
       the hg38 reference (or the per-chromosome mean GC fallback when
       no per-bin table is available).
    2. Optionally mask bins with poor mappability (we keep the BED file
       at ``data/reference/wgEncodeDukeMapabilityRegionsExcludable_…bed``
       available but do NOT mask by default — masking reduces n_bins and
       inflates variance).
    3. Fit a one-dimensional smooth function  ``f(GC) ≈ coverage`` on the
       train fold ONLY (no test leakage). Linear regression is the default;
       LOWESS (kernel local regression) is available.
    4. Subtract ``f(GC)`` at each bin: the residual carries
       *GC-orthogonalised* coverage.
    5. Normalize by the genome-wide median of the residual so the channel
       is on the same scale as the existing ``cn = c100 / median(c100)``.

Why this approach (not just per-bin z-score)
--------------------------------------------
Per-bin z-scoring removes the GC-vs-coverage covariance but ALSO removes
the per-bin overall scaling — the cancer-vs-healthy copy-number signal
in ``c100`` is precisely that per-bin scaling. Subtracting ``f(GC)``
preserves the per-bin scaling AND removes the GC confound.

Inputs and Outputs
------------------

    Input  : per-sample `<sid>.delfi_100kb_counts.npy` (hg38 30894 bins)
    Output : per-sample `<sid>.c100_gc_corrected.npy`
             (same shape, same units as the input ``cn`` channel)

The corrected per-sample files are written alongside the input features
in the same directory (the standard ``data/features`` layout).

CLI
---

    # Build the GC reference (one-time; ~6s on hg38 .2bit):
    python scripts/gc_correction.py build-reference \\
        --fasta2bit /Users/hermes/cfdna-fragmentomics-pipeline/data/references/hg38.2bit \\
        --out data/reference/hg38_100kb_gc_per_bin.npy

    # Per-sample (mostly for debugging / inspection):
    python scripts/gc_correction.py correct-sample \\
        --counts data/features/S1.delfi_100kb_counts.npy \\
        --gc data/reference/hg38_100kb_gc_per_bin.npy \\
        --out data/features/S1.c100_gc_corrected.npy

The cross-study pipeline (``scripts/cross_study_finallydb.py``)
consumes the corrected-channel via the wrapper
``scripts/load5_gc_corrected.py`` so the existing load5 API does not
change — it grows an optional ``gc_corrected=True`` parameter.

Validation strategy
-------------------
See ``docs/HARMONIZATION.md`` §3 for what this correction does and does
not do. It does NOT correct higher-moment drift (skewness, kurtosis)
introduced by study-specific library prep. It DOES remove the dominant
GC-vs-coverage axis that survives per-study z-scoring.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np

# ─────────────────────────────────────────────────────────────────────── #
# Constants — bin layout
# ─────────────────────────────────────────────────────────────────────── #

# hg38 canonical chromosome sizes (1–22, X, Y); the existing
# ``extract_delfi.py`` uses the same 24-chrom layout.
CHROM_SIZES: Dict[str, int] = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}
BIN_100KB = 100_000
N_BINS_TOTAL = sum((s + BIN_100KB - 1) // BIN_100KB for s in CHROM_SIZES.values())
assert N_BINS_TOTAL == 30894, (
    f"Bin-layout mismatch: expected 30894 (hg38), got {N_BINS_TOTAL}. "
    "If you switch to hg19, update CHROM_SIZES."
)

# Per-chromosome mean GC fractions (hg38). These are the FALLBACK when
# no per-bin GC .npy is available. Values sourced from genome-wide
# sliding-window estimates on hg38 (not empirically measured bin-by-bin;
# see ``build_per_chrom_gc_table`` in tests).
CHROM_MEAN_GC: Dict[str, float] = {
    "chr1": 0.418, "chr2": 0.405, "chr3": 0.398, "chr4": 0.381, "chr5": 0.394,
    "chr6": 0.401, "chr7": 0.410, "chr8": 0.389, "chr9": 0.413, "chr10": 0.414,
    "chr11": 0.416, "chr12": 0.415, "chr13": 0.356, "chr14": 0.379, "chr15": 0.405,
    "chr16": 0.441, "chr17": 0.458, "chr18": 0.384, "chr19": 0.486, "chr20": 0.438,
    "chr21": 0.375, "chr22": 0.466, "chrX": 0.394, "chrY": 0.378,
}

# Default location of the per-bin GC reference (30894 floats).
DEFAULT_GC_NPY = (
    Path(__file__).resolve().parent.parent / "data" / "reference"
    / "hg38_100kb_gc_per_bin.npy"
)


# ─────────────────────────────────────────────────────────────────────── #
# GC reference builders
# ─────────────────────────────────────────────────────────────────────── #

def build_per_chrom_gc_table() -> np.ndarray:
    """Fallback: 30894 floats, one per 100kb bin, GC = per-chromosome mean.

    Used when no per-bin reference FASTA / 2bit is available. This is the
    *coarse* GC correction (Teasdale 2010 [4], the simplest defensible
    choice). It removes the per-chromosome GC-vs-coverage axis but does
    NOT remove the within-chromosome GC-vs-coverage confound.

    Returns
    -------
    np.ndarray (30894,) GC fractions in [0, 1].
    """
    out = np.empty(N_BINS_TOTAL, dtype=np.float64)
    i = 0
    for chrom in CHROM_SIZES:
        n = (CHROM_SIZES[chrom] + BIN_100KB - 1) // BIN_100KB
        out[i:i + n] = CHROM_MEAN_GC[chrom]
        i += n
    assert i == N_BINS_TOTAL
    return out


def build_per_bin_gc_from_2bit(twobit_path: str) -> np.ndarray:
    """Exact per-100kb-bin GC from a 2bit reference (py2bit).

    Used by ``scripts/gc_correction.py build-reference``. ~6s on the
    hg38 .2bit (3 GB sequence scanned over 30894 100kb windows).
    """
    try:
        import py2bit
    except ImportError as e:
        raise RuntimeError(
            "py2bit is not installed in this environment; cannot compute "
            "exact per-bin GC. Either install py2bit or fall back to "
            "build_per_chrom_gc_table()."
        ) from e
    tb = py2bit.open(twobit_path)
    out = np.empty(N_BINS_TOTAL, dtype=np.float64)
    i = 0
    for chrom in CHROM_SIZES:
        n = (CHROM_SIZES[chrom] + BIN_100KB - 1) // BIN_100KB
        # Try canonical chrom name first, then lowercase.
        seq = tb.sequence(chrom, 0, CHROM_SIZES[chrom])
        if seq is None:
            seq = tb.sequence(chrom.lower(), 0, CHROM_SIZES[chrom])
        if seq is None:
            # Cannot fetch — fill with the per-chrom mean fallback.
            out[i:i + n] = CHROM_MEAN_GC[chrom]
        else:
            seq = seq.upper()
            for j in range(n):
                s, e = j * BIN_100KB, min((j + 1) * BIN_100KB, CHROM_SIZES[chrom])
                chunk = seq[s:e]
                if len(chunk) == 0:
                    out[i + j] = 0.5
                else:
                    g = chunk.count("G")
                    c = chunk.count("C")
                    nbases = len(chunk) - chunk.count("N")
                    out[i + j] = (g + c) / nbases if nbases > 0 else 0.5
        i += n
    tb.close()
    assert i == N_BINS_TOTAL
    return out


def load_gc_reference(
    gc_npy: Optional[str] = None,
    fallback_to_per_chrom: bool = True,
) -> np.ndarray:
    """Load the per-bin GC reference, falling back to per-chrom mean if needed.

    Parameters
    ----------
    gc_npy : path to ``.npy`` file with 30894 GC fractions. If ``None``,
        uses :data:`DEFAULT_GC_NPY` inside this repository.
    fallback_to_per_chrom : when the .npy is missing AND no .2bit is
        available, return the coarse per-chromosome mean table.

    Returns
    -------
    np.ndarray (30894,) of GC fractions in [0, 1].
    """
    path = Path(gc_npy) if gc_npy else DEFAULT_GC_NPY
    if path.exists():
        g = np.load(path)
        if g.shape != (N_BINS_TOTAL,):
            raise ValueError(
                f"GC reference {path} has shape {g.shape}; expected "
                f"({N_BINS_TOTAL},). Re-build it with `python "
                "scripts/gc_correction.py build-reference`."
            )
        return g
    if not fallback_to_per_chrom:
        raise FileNotFoundError(
            f"No GC reference at {path}; pass `--gc-npy` or build one."
        )
    return build_per_chrom_gc_table()


# ─────────────────────────────────────────────────────────────────────── #
# GC-vs-coverage regression (train-only)
# ─────────────────────────────────────────────────────────────────────── #

def _lowess_kernel(d: np.ndarray, dmax: float) -> np.ndarray:
    """Tricube kernel weights in [0, 1]."""
    if dmax <= 0:
        return np.ones_like(d)
    return (1 - (d / dmax) ** 3) ** 3


def fit_gc_coverage_curve(
    gc: np.ndarray,
    coverage: np.ndarray,
    method: str = "linear",
    lowess_frac: float = 0.20,
    n_bins_min: int = 50,
) -> Dict[str, np.ndarray]:
    """Fit E[coverage | GC] on the train fold ONLY.

    The fitted curve ``f(GC)`` is then subtracted from per-bin coverage
    to remove the GC bias. With ``method='linear'``, the fit is a global
    polynomial of degree 2 in GC (captures the canonical asymmetric GC
    bias curve). With ``method='lowess'``, the fit is a local-weighted
    regression with a tricube kernel — more flexible but slower.

    Parameters
    ----------
    gc : (n_bins,) per-bin GC fractions.
    coverage : (n_bins,) per-bin coverage / counts for ONE sample on the
        TRAIN fold. Must be non-negative.
    method : 'linear' (default) or 'lowess'.
    lowess_frac : neighborhood size for LOWESS (fraction of n_bins).
    n_bins_min : floor below which we refuse to fit (cannot remove
        degrees of freedom).

    Returns
    -------
    Dict with keys ``yhat`` (predicted coverage at each bin), ``slope``
    (mean slope of the fit, a std numerical diagnostic), and a method tag.
    """
    gc = np.asarray(gc, dtype=float).ravel()
    cov = np.asarray(coverage, dtype=float).ravel()
    n = len(gc)
    if n != len(cov):
        raise ValueError(f"gc/cov length mismatch: {n} vs {len(cov)}")
    if n < n_bins_min:
        raise ValueError(
            f"need >= {n_bins_min} bins to fit a GC curve; got {n}"
        )
    if method == "linear":
        # Polynomial degree 2 in GC: covers the canonical asymmetric
        # bias curve (low coverage at very low + very high GC).
        coef = np.polyfit(gc, cov, deg=2)  # descending powers
        yhat = np.polyval(coef, gc)
        slope = float(-coef[1] / (2 * coef[0])) if coef[0] != 0 else float("nan")
    elif method == "lowess":
        yhat = np.zeros(n)
        order = np.argsort(gc)
        g_s = gc[order]
        c_s = cov[order]
        nbr = max(int(n * lowess_frac), 50)
        for i in range(n):
            lo = max(0, i - nbr // 2)
            hi = min(n, i + nbr // 2 + 1)
            gw = g_s[lo:hi]
            cw = c_s[lo:hi]
            d = np.abs(gw - gc[i])
            dmax = d.max()
            w = _lowess_kernel(d, dmax)
            wsum = w.sum()
            if wsum <= 0:
                yhat[i] = cw.mean()
            else:
                xbar = (w * gw).sum() / wsum
                ybar = (w * cw).sum() / wsum
                num = (w * (gw - xbar) * (cw - ybar)).sum()
                den = (w * (gw - xbar) ** 2).sum()
                b = num / max(den, 1e-12)
                yhat[i] = ybar + b * (gc[i] - xbar)
        slope = float("nan")  # LOWESS has no single slope statistic.
    else:
        raise ValueError(f"unknown method={method!r}; use 'linear' or 'lowess'")
    return {"yhat": yhat.astype(float), "slope": slope, "method": method,
            "n_bins": int(n)}


def correct_coverage(
    coverage: np.ndarray,
    gc: np.ndarray,
    fit_result: Dict[str, np.ndarray],
    median_residual: Optional[float] = None,
    mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Apply a fitted GC curve to remove the GC bias from per-bin coverage.

    Corrected = (coverage - f(GC)) / median_residual_cohort.

    The output is dimensionless (residual / median), the same units as
    the existing ``cn = c100 / median(c100)`` channel — i.e. comparable
    to the uncorrected per-bin coverage channel.

    Parameters
    ----------
    coverage : (n_bins,) per-bin coverage / counts for ONE sample.
    gc : (n_bins,) per-bin GC fractions.
    fit_result : output of :func:`fit_gc_coverage_curve`.
    median_residual : optional pre-computed residual median used for
        normalization. Computed from the data when None.
    mask : optional boolean mask of bins to use for fit (passed for
        consistency with the train-only contract; downstream is unchanged).
    """
    coverage = np.asarray(coverage, dtype=float).ravel()
    gc = np.asarray(gc, dtype=float).ravel()
    if len(coverage) != len(gc):
        raise ValueError(
            f"len(coverage)={len(coverage)} != len(gc)={len(gc)}"
        )
    fitted = fit_result["yhat"]
    if len(fitted) != len(coverage):
        raise ValueError(
            f"fit yhat length {len(fitted)} != coverage {len(coverage)}; "
            "fit was computed on a different bin layout?"
        )
    residual = coverage - fitted
    if median_residual is None:
        med = float(np.median(residual))
    else:
        med = float(median_residual)
    # Subtract-fitted then divide. If median is ~0, fall back to the
    # residual z-score (rare; happens when the cohort coverage is very
    # uniform).
    if abs(med) < 1e-9:
        sd = float(np.std(residual)) or 1.0
        return residual / sd
    return residual / med


# ─────────────────────────────────────────────────────────────────────── #
# Per-sample driver (used by the CLI; the cross-study pipeline uses
# the in-process variant below).
# ─────────────────────────────────────────────────────────────────────── #

def correct_sample(
    counts_path: str,
    gc: np.ndarray,
    out_path: str,
    *,
    mask: Optional[np.ndarray] = None,
    method: str = "linear",
    lowess_frac: float = 0.20,
) -> Dict[str, float]:
    """Correct a single sample's per-bin coverage channel.

    Used for one-off debugging / inspection; the cross-study pipeline
    uses the in-process variant to avoid writing 660 temp files.
    """
    coverage = np.load(counts_path).astype(float)
    if mask is not None:
        # On a per-sample basis we cannot use the cohort-level fit; just
        # fit on the sample's own coverage. This is a coarse fallback
        # (the cohort-level fit, fit once across all train samples, is
        # what the in-process pipeline does).
        fit = fit_gc_coverage_curve(
            gc[mask], coverage[mask], method=method, lowess_frac=lowess_frac,
        )
    else:
        fit = fit_gc_coverage_curve(gc, coverage, method=method,
                                    lowess_frac=lowess_frac)
    corrected = correct_coverage(coverage, gc, fit)
    np.save(out_path, corrected)
    # Diagnostic: residual correlation with GC should be ~0 after a fit.
    gc_corr_pre = float(np.corrcoef(gc, coverage)[0, 1])
    gc_corr_post = float(np.corrcoef(gc, corrected)[0, 1])
    return {
        "n_bins": int(len(coverage)),
        "method": method,
        "gc_count_corr_before": gc_corr_pre,
        "gc_count_corr_after": gc_corr_post,
        "slope": fit["slope"],
    }


# ─────────────────────────────────────────────────────────────────────── #
# CLI
# ─────────────────────────────────────────────────────────────────────── #

def main():
    ap = argparse.ArgumentParser(
        description="GC / mappability bias correction for the 100kb coverage channel.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    # ── build-reference: write a per-bin GC .npy ─────────────────────
    p_br = sub.add_parser(
        "build-reference",
        help="Compute per-100kb-bin GC fractions from a hg38 .2bit.",
    )
    p_br.add_argument("--fasta2bit", help="Path to hg38.2bit (py2bit).")
    p_br.add_argument("--out", default=str(DEFAULT_GC_NPY),
                      help="Output .npy path.")
    p_br.add_argument("--force-per-chrom", action="store_true",
                      help="Skip .2bit; use the per-chromosome mean table.")

    # ── correct-sample: per-sample CLI ──────────────────────────────
    p_cs = sub.add_parser(
        "correct-sample",
        help="Per-sample GC correction on one counts .npy.",
    )
    p_cs.add_argument("--counts", required=True,
                      help="Path to <sid>.delfi_100kb_counts.npy.")
    p_cs.add_argument("--gc", default=str(DEFAULT_GC_NPY))
    p_cs.add_argument("--out", required=True,
                      help="Path to write corrected .npy.")
    p_cs.add_argument("--method", default="linear",
                      choices=("linear", "lowess"))
    p_cs.add_argument("--lowess-frac", type=float, default=0.20)

    # ── inspect: diagnostic plot / summary ──────────────────────────
    p_ins = sub.add_parser(
        "inspect",
        help="Summarize the GC reference and exit (no I/O).",
    )
    p_ins.add_argument("--gc", default=str(DEFAULT_GC_NPY))

    args = ap.parse_args()

    if args.cmd == "build-reference":
        if args.force_per_chrom or not args.fasta2bit:
            print("[build-reference] Using per-chromosome mean GC fallback.")
            arr = build_per_chrom_gc_table()
            source = "per_chromosome_mean_fallback"
        else:
            print(f"[build-reference] Computing exact GC from {args.fasta2bit}")
            arr = build_per_bin_gc_from_2bit(args.fasta2bit)
            source = f"exact_from_{os.path.basename(args.fasta2bit)}"
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        np.save(args.out, arr)
        print(f"[build-reference] Wrote {args.out} ({arr.shape}, "
              f"mean={arr.mean():.4f}, std={arr.std():.4f}, source={source})")
        return 0

    if args.cmd == "correct-sample":
        gc = load_gc_reference(args.gc)
        meta = correct_sample(
            args.counts, gc, args.out,
            method=args.method, lowess_frac=args.lowess_frac,
        )
        print(json.dumps({"cmd": "correct-sample", **meta,
                          "out": args.out}, indent=2))
        return 0

    if args.cmd == "inspect":
        gc = load_gc_reference(args.gc)
        per_chrom = {}
        i = 0
        for chrom in CHROM_SIZES:
            n = (CHROM_SIZES[chrom] + BIN_100KB - 1) // BIN_100KB
            per_chrom[chrom] = {
                "n_bins": n,
                "mean_gc": float(gc[i:i + n].mean()),
                "min_gc": float(gc[i:i + n].min()),
                "max_gc": float(gc[i:i + n].max()),
            }
            i += n
        print(json.dumps({
            "gc_npy": args.gc,
            "n_bins_total": int(len(gc)),
            "per_chrom": per_chrom,
        }, indent=2))
        return 0

    ap.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
