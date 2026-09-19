"""
AlphaGenome Atlas AVI-weighted panel-LLR comparison.

Runs the real-TCGA-LUAD panel-LLR benchmark with four aggregation
strategies, in parallel on the same cohort and seeds:

  1. panel_llr_uniform     — Σ LLR_i (existing baseline; AUC 0.921 @ 0.1%)
  2. panel_llr_avi         — Σ w_i · LLR_i   (AVI-weighted, full panel)
  3. panel_llr_topk_500    — Σ w_i · LLR_i   restricted to top-500 AVI
  4. panel_llr_topk_1000   — Σ w_i · LLR_i   restricted to top-1000 AVI
  5. panel_llr_topk_2000   — Σ w_i · LLR_i   restricted to top-2000 AVI

AVI weights are loaded once from the AlphaGenome Atlas API if an API key
is available, otherwise from a Tabix cache, otherwise from a deterministic
proxy that mirrors the published AVI training distribution.  All runs
share the same seeds and simulated read counts, so per-patient deltas are
attributable to the weighting change alone.

Writes:
  results/alphagenome_weighted_llr.json
  docs/ALPHAGENOME_WEIGHTED_LLR.md  (handled in a separate step)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import real_tcga_validation as rtv  # noqa: E402
import alphagenome_weights as aw    # noqa: E402

LOG = logging.getLogger("alphagenome_panel_run")


# ---------------------------------------------------------------------------
# Per-position LLR with optional weight
# ---------------------------------------------------------------------------

def _pos_llr(depths, obs_alt, errors):
    return rtv.compute_llr_scores(depths, obs_alt, errors)


def _panel_score_with_weights(per_pos_llr: np.ndarray,
                              weights: Optional[np.ndarray],
                              top_k: Optional[int] = None) -> float:
    """Deterministic weighted-LLR panel score.

    If weights is None: uniform sum (baseline).
    If top_k is given: only the top-K weighted positions contribute.

    Top-K here is selected by `w * LLR` ranking — a joint criterion that
    conflates AVI prior with LLR magnitude. See
    `_panel_score_with_weights_by_avi` for the pure-prior selector that
    decouples the two.
    """
    if weights is None:
        contrib = per_pos_llr
    else:
        contrib = weights * per_pos_llr
    if top_k is not None and top_k < len(contrib):
        # Use the top-K LARGEST contributions (positive LLR, ranked)
        idx = np.argsort(contrib)[::-1][:top_k]
        return float(contrib[idx].sum())
    return float(contrib.sum())


def _panel_score_with_weights_by_avi(per_pos_llr: np.ndarray,
                                     avi_norm: np.ndarray,
                                     weights: np.ndarray,
                                     top_k: int) -> float:
    """Top-K-by-raw-AVI panel score: select positions by AVI priority only,
    then sum their weighted LLR (weights * LLR).

    Decouples AVI prior (`avi_norm`) from LLR magnitude: positions are
    chosen purely by their AVI priority, regardless of whether the LLR at
    that position is large or small. The aggregation is still weighted,
    so the contribution of a chosen position scales with its AVI weight.

    Unlike `_panel_score_with_weights(..., top_k=...)` which ranks by the
    joint `w * LLR` criterion, this method ranks by AVI alone. The two
    selectors are equivalent only when LLR ordering happens to match AVI
    ordering rank-by-rank, which is rare.
    """
    if len(avi_norm) != len(per_pos_llr):
        raise ValueError(
            f"avi_norm length ({len(avi_norm)}) != per_pos_llr length "
            f"({len(per_pos_llr)})"
        )
    # Select top-K by raw AVI priority (largest first), independent of LLR.
    idx = np.argsort(avi_norm)[::-1][:top_k]
    return float((weights[idx] * per_pos_llr[idx]).sum())


# ---------------------------------------------------------------------------
# Main comparison runner
# ---------------------------------------------------------------------------

def run_comparison(
    cohort: Dict[str, Any],
    mutations_with_ref: Sequence[Dict],
    avi_weights: Dict[str, aw.AVIWeight],
    *,
    tumor_fractions: Sequence[float] = (0.1, 0.05, 0.01, 0.005, 0.001),
    seeds: Sequence[int] = (42, 123, 456, 789, 1024),
    cfdna_depth: int = 5000,
    bg_error_rate: float = 0.002,
    topk_values: Sequence[int] = (500, 1000, 2000),
) -> Dict[str, Any]:
    """Run the 5-strategy comparison and aggregate across seeds.

    Per (seed, tumor_fraction, strategy):
      - sim per-patient → (pos_score, neg_score)
      - ROC AUC across patients for that seed
    Per (tumor_fraction, strategy):
      - mean ± std of per-seed AUC
    """
    methods = ["panel_llr_uniform", "panel_llr_avi"]
    methods += [f"panel_llr_topk_{k}" for k in topk_values]      # top-K by w*LLR (joint)
    methods += [f"panel_llr_topk_by_avi_{k}" for k in topk_values]  # top-K by raw AVI (pure prior)

    # Pre-build per-patient weights (sorted by AVI desc within patient)
    patients = list(cohort["patients"].keys())
    patient_weights: Dict[str, np.ndarray] = {}
    patient_avi: Dict[str, np.ndarray] = {}      # raw AVI norm (for top-K-by-AVI selector)
    patient_keys: Dict[str, List[str]] = {}
    for p in patients:
        mlist = cohort["patients"][p]
        # Build AVI vector in mutation-list order. Missing weights -> 0 (treated
        # as "no signal" — the LLR contribution is also 0 for error-only reads
        # so this is safe).
        ws: List[float] = []
        avis: List[float] = []
        ks: List[str] = []
        for m in mlist:
            # Look up by mutation identity: the cohort mutation dicts lack
            # ref/alt alleles, so we match by (sample, chrom, pos, gene)
            # against the mutations_with_ref list.
            vkey = None
            for mr in mutations_with_ref:
                if (mr["sample"] == m.get("sample")
                        and mr["chrom"] == m.get("chrom")
                        and mr["pos"] == m.get("pos")
                        and mr["gene"] == m.get("gene")):
                    vkey = aw.variant_key(mr["chrom"], mr["pos"], mr["ref"], mr["alt"])
                    break
            if vkey and vkey in avi_weights:
                ws.append(avi_weights[vkey].avi_norm)
                avis.append(avi_weights[vkey].avi_norm)
            else:
                ws.append(0.0)
                avis.append(0.0)
            ks.append(vkey or "?")
        patient_weights[p] = np.asarray(ws, dtype=float)
        patient_avi[p] = np.asarray(avis, dtype=float)
        patient_keys[p] = ks

    # Aggregation container
    results: Dict[str, List[Dict]] = {m: [] for m in methods}

    METRIC_FNS = {
        "auc": lambda y, s: float(__import__("sklearn.metrics", fromlist=["roc_auc_score"]).roc_auc_score(y, s)),
        "sens_at_95_spec": lambda y, s: rtv.sensitivity_at_specificity(y, s, 0.95),
        "paired_win_rate": lambda y, s: float(np.mean(s[: len(s) // 2] > s[len(s) // 2 :])),
    }

    t0 = time.time()
    for tf in tumor_fractions:
        print(f"\n=== Tumor fraction {tf*100:.2f}% ===")
        # Per-method per-seed AUC arrays
        method_auc_by_seed: Dict[str, Dict[int, float]] = {m: {} for m in methods}

        for seed in seeds:
            # Per-patient simulated data
            sim_p: Dict[str, Dict] = {}
            sim_n: Dict[str, Dict] = {}
            llr_p: Dict[str, np.ndarray] = {}
            llr_n: Dict[str, np.ndarray] = {}
            for p in patients:
                muts = cohort["patients"][p]
                dp = rtv.simulate_cfdna_from_real(
                    muts, tumor_fraction=tf, cfdna_depth=cfdna_depth,
                    seed=seed, bg_error_rate=bg_error_rate,
                )
                dn = rtv.simulate_cfdna_from_real(
                    muts, tumor_fraction=0.0, cfdna_depth=cfdna_depth,
                    seed=seed, bg_error_rate=bg_error_rate,
                )
                nv_p, nv_n = dp["n_variants"], dn["n_variants"]
                panel_size = min(nv_p, nv_n, len(muts))
                # Per-position LLR (truncated to panel_size for both)
                lp = _pos_llr(dp["depths"][:panel_size],
                              dp["X"][:, 1].astype(int)[:panel_size],
                              dp["X"][:, 3][:panel_size])
                ln = _pos_llr(dn["depths"][:panel_size],
                              dn["X"][:, 1].astype(int)[:panel_size],
                              dn["X"][:, 3][:panel_size])
                llr_p[p] = lp
                llr_n[p] = ln

            for method in methods:
                pos_scores, neg_scores = [], []
                for p in patients:
                    lp = llr_p[p]
                    ln = llr_n[p]
                    w = patient_weights[p]
                    if method == "panel_llr_uniform":
                        sp = _panel_score_with_weights(lp, None)
                        sn = _panel_score_with_weights(ln, None)
                    elif method == "panel_llr_avi":
                        sp = _panel_score_with_weights(lp, w)
                        sn = _panel_score_with_weights(ln, w)
                    elif method.startswith("panel_llr_topk_by_avi_"):
                        # Top-K by raw AVI priority (decoupled from LLR)
                        k = int(method.split("_")[-1])
                        sp = _panel_score_with_weights_by_avi(lp, patient_avi[p], w, top_k=k)
                        sn = _panel_score_with_weights_by_avi(ln, patient_avi[p], w, top_k=k)
                    else:
                        # panel_llr_topk_{K} — top-K by w*LLR (joint criterion)
                        k = int(method.split("_")[-1])
                        sp = _panel_score_with_weights(lp, w, top_k=k)
                        sn = _panel_score_with_weights(ln, w, top_k=k)
                    pos_scores.append(sp)
                    neg_scores.append(sn)

                y = np.array([1] * len(pos_scores) + [0] * len(neg_scores))
                s = np.array(pos_scores + neg_scores)
                auc = float(__import__("sklearn.metrics", fromlist=["roc_auc_score"]).roc_auc_score(y, s))
                method_auc_by_seed[method][seed] = auc

            print(f"  seed {seed:>4}:  "
                  + "  ".join(f"{m}={method_auc_by_seed[m][seed]:.4f}" for m in methods))

        # Aggregate across seeds
        for method in methods:
            auc_vals = list(method_auc_by_seed[method].values())
            results[method].append({
                "tumor_fraction": float(tf),
                "metric": "auc",
                "mean": float(np.mean(auc_vals)),
                "std": float(np.std(auc_vals, ddof=1)) if len(auc_vals) > 1 else 0.0,
                "per_seed": {str(s): float(v) for s, v in method_auc_by_seed[method].items()},
            })

    elapsed = time.time() - t0
    return {
        "methods": methods,
        "results": results,
        "elapsed_seconds": elapsed,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="AVI-weighted panel LLR comparison")
    parser.add_argument("--cache-dir",
                        default="validation/tcga/tcga_cache",
                        help="TCGA MAF cache directory")
    parser.add_argument("--output",
                        default="results/alphagenome_weighted_llr.json")
    parser.add_argument("--n-patients", type=int, default=20)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--cfdna-depth", type=int, default=5000)
    parser.add_argument("--bg-error-rate", type=float, default=0.002)
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--tabix", default=None)
    parser.add_argument("--weights-cache",
                        default="results/alphagenome_avi_weights.json")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    print("=" * 72)
    print("  DeepCatch — AVI-weighted Panel-LLR Comparison (v2.2)")
    print("=" * 72)

    # 1. Load cohort + ref/alt-augmented mutations
    cohort = rtv.load_tcga_cohort(
        args.cache_dir,
        n_patients=args.n_patients,
        cancer_types=["LUAD"],
        allow_download=True,
    )
    mutations_with_ref = aw.load_tcga_mutations_with_ref_alt(args.cache_dir)
    print(f"  Cohort: {cohort['n_patients']} patients, {cohort['n_mutations']} mutations")
    print(f"  Ref/alt-augmented pool: {len(mutations_with_ref)} mutations "
          f"({sum(1 for m in mutations_with_ref if len(m['ref'])==1 and len(m['alt'])==1)} SNVs)")

    # 2. Weight cohort
    weights_cache = Path(args.weights_cache)
    avi_weights, primary_source = aw.weight_cohort(
        mutations_with_ref,
        api_key=args.api_key,
        tabix_path=Path(args.tabix) if args.tabix else None,
        cache_path=weights_cache,
    )
    n_real = sum(1 for w in avi_weights.values() if w.source in ("atlas_api", "tabix_local"))
    n_proxy = sum(1 for w in avi_weights.values() if w.source == "proxy")
    n_missing = sum(1 for w in avi_weights.values() if w.source == "missing")
    print(f"  AVI source: {primary_source}  "
          f"(real={n_real}  proxy={n_proxy}  missing={n_missing})")

    # 3. Normalise weights
    raw_scores = np.array([w.avi_score for w in avi_weights.values()], dtype=float)
    norm = aw.normalize_avi(raw_scores)
    for vk, norm_val in zip(avi_weights.keys(), norm):
        avi_weights[vk].avi_norm = float(norm_val)
    print(f"  AVI raw range:  {raw_scores.min():.2f} .. {raw_scores.max():.2f}")
    print(f"  AVI norm range: {norm.min():.3f} .. {norm.max():.3f}")

    # 4. Run comparison
    seeds = [42, 123, 456, 789, 1024][: args.seeds]
    tumor_fractions = [0.1, 0.05, 0.01, 0.005, 0.001]

    t0 = time.time()
    out = run_comparison(
        cohort,
        mutations_with_ref,
        avi_weights,
        tumor_fractions=tumor_fractions,
        seeds=seeds,
        cfdna_depth=args.cfdna_depth,
        bg_error_rate=args.bg_error_rate,
        topk_values=(500, 1000, 2000),
    )
    elapsed = time.time() - t0

    # 5. Build output JSON
    payload = {
        "metadata": {
            "runner": "alphagenome_panel_run.py",
            "date": time.strftime("%Y-%m-%d %H:%M:%S"),
            "data_source": cohort["source"],
            "cancer_types": ["LUAD"],
            "n_patients": cohort["n_patients"],
            "n_mutations_in_cohort": cohort["n_mutations"],
            "n_mutations_with_ref_alt": len(mutations_with_ref),
            "cfdna_depth": args.cfdna_depth,
            "bg_error_rate": args.bg_error_rate,
            "seeds_used": list(seeds),
            "tumor_fractions": list(tumor_fractions),
            "avi_primary_source": primary_source,
            "avi_n_real": n_real,
            "avi_n_proxy": n_proxy,
            "avi_n_missing": n_missing,
            "avi_raw_range": [float(raw_scores.min()), float(raw_scores.max())],
            "avi_norm_range": [float(norm.min()), float(norm.max())],
            "methods_compared": out["methods"],
            "elapsed_seconds": elapsed,
            # Honest framing — identical to the existing pipeline
            "pipeline_type": "REAL_MUTATIONS_+_SIMULATED_PLASMA_READS",
            "note": (
                "Ground-truth variants come from real TCGA-LUAD MAF data; "
                "plasma reads are simulated by Poisson sampling. Panel LLR "
                "is computed with AVI weights from AlphaGenome Atlas when "
                "available; otherwise a deterministic per-variant-class proxy "
                "is used (see avi_primary_source). Avi weights are applied "
                "as fixed scalars in a deterministic aggregation — not as "
                "features to a learned model — which complies with the "
                "AlphaGenome Terms of Use non-training clause."
            ),
        },
        "results": out["results"],
        "elapsed_seconds": elapsed,
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"\n  📁 Results saved to {args.output}")
    print(f"  ⏱  Wall time: {elapsed:.1f}s")

    # 6. Print compact comparison table
    print("\n" + "=" * 72)
    print("  PANEL-LLR COMPARISON  (mean AUC ± std across seeds)")
    print("=" * 72)
    hdr = f"  {'TF':<6}" + "".join(f"{m:>22}" for m in out["methods"])
    print(hdr)
    for tf in tumor_fractions:
        line = f"  {tf*100:5.2f}%"
        for m in out["methods"]:
            row = next((r for r in out["results"][m]
                        if r["tumor_fraction"] == tf and r["metric"] == "auc"), None)
            if row is None:
                line += f"    {'n/a':>18}"
            else:
                line += f"  {row['mean']:.4f}±{row['std']:.4f}".rjust(22)
        print(line)
    print("=" * 72)

    # 7. Bottom-line deltas vs baseline (uniform)
    print("\n  Bottom-line Δ vs panel_llr_uniform (AUC @ 0.1%):")
    tf_key = 0.001
    base = next((r["mean"] for r in out["results"]["panel_llr_uniform"]
                 if r["tumor_fraction"] == tf_key), None)
    if base is not None:
        for m in out["methods"]:
            if m == "panel_llr_uniform":
                continue
            row = next((r for r in out["results"][m]
                        if r["tumor_fraction"] == tf_key), None)
            if row is not None:
                d = row["mean"] - base
                flag = "✓" if d >= 0.02 else ("≈" if d >= 0 else "✗")
                print(f"    {flag}  {m:>22}: AUC {row['mean']:.4f}  Δ={d:+.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())