"""
AlphaMissense-weighted panel-LLR comparison for DeepCatch v2.2.

Runs the real-TCGA-LUAD panel-LLR benchmark with four aggregation
strategies, on the same cohort and seeds:

  1. panel_llr_uniform     — Σ LLR_i  (existing baseline; AUC 0.921 @ 0.1%)
  2. panel_llr_am          — Σ w_i · LLR_i   (AM-pathogenicity-weighted)
  3. panel_llr_topk_500    — Σ w_i · LLR_i   restricted to top-500 AM
  4. panel_llr_topk_1000   — Σ w_i · LLR_i   restricted to top-1000 AM
  5. panel_llr_topk_2000   — Σ w_i · LLR_i   restricted to top-2000 AM
  6. panel_llr_pathogenic  — Σ w_i · LLR_i   restricted to AM-pathogenic
                              (score >= 0.564 per Cheng 2023 thresholds)

AM weights are loaded once from a local per-Uniprot pickle if present,
else from the raw AlphaMissense_hg38.tsv if present, else from a
deterministic per-variant-class published-prior PROXY (clearly flagged
in the output). All runs share the same seeds and simulated read counts
so per-patient deltas are attributable to the weighting change alone.

Writes:
  results/alphamissense_weighted_llr.json

The companion `docs/ALPHAMISSENSE_WEIGHTED_LLR.md` is written by a
separate step and contains the honest findings table.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import real_tcga_validation as rtv  # noqa: E402
import alphamissense_weights as amw  # noqa: E402

LOG = logging.getLogger("alphamissense_panel_run")


# ---------------------------------------------------------------------------
# Per-position LLR with optional weight
# ---------------------------------------------------------------------------

def _pos_llr(depths, obs_alt, errors):
    return rtv.compute_llr_scores(depths, obs_alt, errors)


def _panel_score_with_weights(per_pos_llr: np.ndarray,
                              weights: Optional[np.ndarray],
                              top_k: Optional[int] = None,
                              pathogenic_only: bool = False,
                              am_norm: Optional[np.ndarray] = None) -> float:
    """Deterministic weighted-LLR panel score.

    - weights is None: uniform sum (baseline).
    - top_k given: only the top-K contributions count (largest positives).
    - pathogenic_only: only loci with AM score >= 0.564 count.
    - am_norm: parallel vector of normalised AM scores (0..1).
    """
    if weights is None:
        contrib = per_pos_llr
    else:
        contrib = weights * per_pos_llr
    if pathogenic_only and am_norm is not None:
        # Use raw AM score (not normalised) for the threshold check
        # The threshold is on the probability scale [0,1].
        # The proxy uses values that include >=0.564 for LoF classes.
        # For pathogenic_only we re-mask using raw scores.
        raise RuntimeError("pathogenic_only requires am_raw via _score_pathogenic")
    if top_k is not None and top_k < len(contrib):
        idx = np.argsort(contrib)[::-1][:top_k]
        return float(contrib[idx].sum())
    return float(contrib.sum())


def _panel_score_pathogenic(per_pos_llr: np.ndarray,
                            weights: np.ndarray,
                            am_raw: np.ndarray,
                            threshold: float = amw.THRESHOLD_PATHOGENIC,
                            top_k: Optional[int] = None) -> float:
    """Weighted LLR summed over only AM-pathogenic loci (score >= threshold).

    am_raw is the parallel vector of raw AM pathogenicity probabilities.
    Loci with am_raw < threshold contribute 0 (treated as not prioritised).
    """
    mask = am_raw >= threshold
    contrib = weights * per_pos_llr * mask.astype(float)
    if top_k is not None and top_k < len(contrib):
        idx = np.argsort(contrib)[::-1][:top_k]
        return float(contrib[idx].sum())
    return float(contrib.sum())


# ---------------------------------------------------------------------------
# Main comparison runner
# ---------------------------------------------------------------------------

def run_comparison(
    cohort: Dict[str, Any],
    cohort_mutations: Sequence[Dict],
    am_weights: Dict[str, amw.AMWeight],
    am_norm: np.ndarray,
    am_raw: np.ndarray,
    *,
    tumor_fractions: Sequence[float] = (0.1, 0.05, 0.01, 0.005, 0.001),
    seeds: Sequence[int] = (42, 123, 456, 789, 1024),
    cfdna_depth: int = 5000,
    bg_error_rate: float = 0.002,
    topk_values: Sequence[int] = (500, 1000, 2000),
) -> Dict[str, Any]:
    """Run the multi-strategy comparison and aggregate across seeds."""
    methods = ["panel_llr_uniform", "panel_llr_am"]
    methods += [f"panel_llr_topk_{k}" for k in topk_values]
    methods += ["panel_llr_pathogenic"]

    # Build per-patient AM weight vectors aligned to cohort['patients'] order.
    # Each patient has its own missense SNVs; we keep them in the MAF
    # order (lexicographic) and emit a weight vector of the same length.
    # Mutations without an AM lookup key get weight=0 (i.e. contribute
    # nothing in the weighted sum).
    patients = list(cohort["patients"].keys())
    patient_w: Dict[str, np.ndarray] = {}
    patient_am_raw: Dict[str, np.ndarray] = {}
    patient_keys: Dict[str, List[str]] = {}
    for p in patients:
        mlist = cohort["patients"][p]
        # For each cohort mutation, find the matching AM weight.
        # Cohort mutations don't carry uniprot/prot_pos — match via
        # (sample, gene, chrom, pos, ref, alt) against cohort_mutations.
        ws: List[float] = []
        rs: List[float] = []
        ks: List[str] = []
        missense_seen: List[Dict] = []
        for m in mlist:
            ref = m.get("ref", "")
            alt = m.get("alt", "")
            if not ref or not alt:
                # Match from the augmented pool
                for mr in cohort_mutations:
                    if (mr["sample"] == m.get("sample")
                            and mr["chrom"] == m.get("chrom")
                            and mr["pos"] == m.get("pos")
                            and mr["gene"] == m.get("gene")):
                        ref, alt = mr["ref"], mr["alt"]
                        break
            vc = m.get("variant_class", "")
            if vc != "Missense_Mutation" or not ref or not alt \
                    or len(ref) != 1 or len(alt) != 1:
                ws.append(0.0)
                rs.append(0.0)
                ks.append("")
                continue
            # Look up via cohort_mutations for the Uniprot / prot_pos.
            uniprot = None; prot_pos = None
            for mr in cohort_mutations:
                if (mr["sample"] == m.get("sample")
                        and mr["chrom"] == m.get("chrom")
                        and mr["pos"] == m.get("pos")
                        and mr["gene"] == m.get("gene")
                        and mr["ref"] == ref and mr["alt"] == alt):
                    uniprot = mr["uniprot"]
                    prot_pos = mr["prot_pos"]
                    break
            if uniprot is None or prot_pos is None:
                ws.append(0.0)
                rs.append(0.0)
                ks.append("")
                continue
            key = amw.am_key(uniprot, prot_pos, ref, alt)
            w = am_weights.get(key)
            if w is None:
                ws.append(0.0)
                rs.append(0.0)
                ks.append(key)
            else:
                ws.append(w.norm)
                rs.append(w.score)
                ks.append(key)
        patient_w[p] = np.asarray(ws, dtype=float)
        patient_am_raw[p] = np.asarray(rs, dtype=float)
        patient_keys[p] = ks

    results: Dict[str, List[Dict]] = {m: [] for m in methods}

    t0 = time.time()
    for tf in tumor_fractions:
        print(f"\n=== Tumor fraction {tf*100:.2f}% ===")
        method_auc_by_seed: Dict[str, Dict[int, float]] = {m: {} for m in methods}

        for seed in seeds:
            # Per-patient simulated data — IDENTICAL across methods for
            # this seed, so per-patient deltas are pure weighting deltas.
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
                panel_size = min(len(muts), len(patient_w[p]))
                lp = _pos_llr(dp["depths"][:panel_size],
                              dp["X"][:, 1].astype(int)[:panel_size],
                              dp["X"][:, 3][:panel_size])
                ln = _pos_llr(dn["depths"][:panel_size],
                              dn["X"][:, 1].astype(int)[:panel_size],
                              dn["X"][:, 3][:panel_size])
                llr_p[p] = lp
                llr_n[p] = ln

            for method in methods:
                pos_scores: List[float] = []
                neg_scores: List[float] = []
                for p in patients:
                    lp = llr_p[p]
                    ln = llr_n[p]
                    w = patient_w[p]
                    raw = patient_am_raw[p]
                    if method == "panel_llr_uniform":
                        sp = _panel_score_with_weights(lp, None)
                        sn = _panel_score_with_weights(ln, None)
                    elif method == "panel_llr_am":
                        sp = _panel_score_with_weights(lp, w)
                        sn = _panel_score_with_weights(ln, w)
                    elif method.startswith("panel_llr_topk_"):
                        k = int(method.split("_")[-1])
                        sp = _panel_score_with_weights(lp, w, top_k=k)
                        sn = _panel_score_with_weights(ln, w, top_k=k)
                    elif method == "panel_llr_pathogenic":
                        sp = _panel_score_pathogenic(lp, w, raw)
                        sn = _panel_score_pathogenic(ln, w, raw)
                    else:
                        raise RuntimeError(f"unknown method: {method}")
                    pos_scores.append(sp)
                    neg_scores.append(sn)

                y = np.array([1] * len(pos_scores) + [0] * len(neg_scores))
                s = np.array(pos_scores + neg_scores)
                auc = float(roc_auc_score(y, s))
                method_auc_by_seed[method][seed] = auc

            print(f"  seed {seed:>4}:  "
                  + "  ".join(f"{m}={method_auc_by_seed[m][seed]:.4f}" for m in methods))

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
    parser = argparse.ArgumentParser(description="AlphaMissense-weighted panel LLR comparison")
    parser.add_argument("--cache-dir",
                        default="validation/tcga/tcga_cache",
                        help="TCGA MAF cache directory")
    parser.add_argument("--output",
                        default="results/alphamissense_weighted_llr.json")
    parser.add_argument("--index", default=None,
                        help="Path to a per-Uniprot AM pickle (optional)")
    parser.add_argument("--tsv", default=None,
                        help="Path to AlphaMissense_hg38.tsv[.gz] (optional)")
    parser.add_argument("--n-patients", type=int, default=20)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--cfdna-depth", type=int, default=5000)
    parser.add_argument("--bg-error-rate", type=float, default=0.002)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    print("=" * 72)
    print("  DeepCatch — AlphaMissense-weighted Panel-LLR Comparison (v2.2)")
    print("=" * 72)
    print("  Data: AlphaMissense (Cheng et al., Science 2023)")
    print("  License: CC BY-NC-SA 4.0  (non-commercial; see docs)")
    print("=" * 72)

    # 1) Load cohort
    cohort = rtv.load_tcga_cohort(
        args.cache_dir,
        n_patients=args.n_patients,
        cancer_types=["LUAD"],
        allow_download=True,
    )

    # 2) Build the augmented pool (carries Uniprot + protein position)
    cohort_mutations = amw.load_tcga_mutations_with_protein(args.cache_dir)
    print(f"  Cohort: {cohort['n_patients']} patients, {cohort['n_mutations']} mutations")
    print(f"  Augmented pool (Uniprot + protein-position): {len(cohort_mutations)} "
          f"missense SNVs (across {len(set(m['sample'] for m in cohort_mutations))} patients)")

    # 3) Weight the cohort
    am_weights, info = amw.weight_cohort(
        cohort_mutations,
        index_path=Path(args.index) if args.index else None,
        tsv_path=Path(args.tsv) if args.tsv else None,
    )
    print(f"  AM primary source: {info['primary_source']}")
    print(f"  AM real={info['n_real']}  proxy={info['n_proxy']}  missing={info['n_missing']}")
    print(f"  Missense total: {info['missense_total']}  "
          f"with full key: {info['missense_with_full_key']}")

    # 4) Normalise
    raw = np.array([w.score for w in am_weights.values()], dtype=float)
    norm = amw.normalize_am(raw)
    for k, v in zip(am_weights.keys(), norm):
        am_weights[k].norm = float(v)
    print(f"  AM raw range:   {raw.min():.3f} .. {raw.max():.3f}")
    print(f"  AM norm range:  {norm.min():.3f} .. {norm.max():.3f}")
    n_pathogenic = sum(1 for w in am_weights.values()
                       if w.score >= amw.THRESHOLD_PATHOGENIC)
    print(f"  Pathogenic (>= {amw.THRESHOLD_PATHOGENIC}): {n_pathogenic}/{len(am_weights)}")

    # 5) Run comparison
    seeds = [42, 123, 456, 789, 1024][: args.seeds]
    tumor_fractions = [0.1, 0.05, 0.01, 0.005, 0.001]

    am_raw_arr = np.array([w.score for w in am_weights.values()], dtype=float)
    am_norm_arr = np.array([w.norm for w in am_weights.values()], dtype=float)

    t0 = time.time()
    out = run_comparison(
        cohort, cohort_mutations, am_weights, am_norm_arr, am_raw_arr,
        tumor_fractions=tumor_fractions,
        seeds=seeds,
        cfdna_depth=args.cfdna_depth,
        bg_error_rate=args.bg_error_rate,
        topk_values=(500, 1000, 2000),
    )
    elapsed = time.time() - t0

    # 6) Build output JSON
    payload = {
        "metadata": {
            "runner": "alphamissense_panel_run.py",
            "date": time.strftime("%Y-%m-%d %H:%M:%S"),
            "data_source": cohort["source"],
            "cancer_types": ["LUAD"],
            "n_patients": cohort["n_patients"],
            "n_mutations_in_cohort": cohort["n_mutations"],
            "n_missense_in_cohort": info["missense_total"],
            "n_missense_with_full_am_key": info["missense_with_full_key"],
            "match_rate_missense_to_am": (
                info["missense_with_full_key"] / info["missense_total"]
                if info["missense_total"] else 0.0
            ),
            "cfdna_depth": args.cfdna_depth,
            "bg_error_rate": args.bg_error_rate,
            "seeds_used": list(seeds),
            "tumor_fractions": list(tumor_fractions),
            "am_primary_source": info["primary_source"],
            "am_n_real": info["n_real"],
            "am_n_proxy": info["n_proxy"],
            "am_n_missing": info["n_missing"],
            "am_raw_range": [float(raw.min()), float(raw.max())],
            "am_norm_range": [float(norm.min()), float(norm.max())],
            "am_n_pathogenic": n_pathogenic,
            "am_index_path": info["index_path"],
            "am_tsv_path": info["tsv_path"],
            "methods_compared": out["methods"],
            "elapsed_seconds": elapsed,
            "pipeline_type": "REAL_MUTATIONS_+_SIMULATED_PLASMA_READS",
            "license": "CC BY-NC-SA 4.0 (AlphaMissense — non-commercial)",
            "citation": ("Cheng J, Novati G, Pan J, et al. "
                         "'Accurate proteome-wide missense variant effect "
                         "prediction with AlphaMissense.' Science 381, "
                         "eadg7492 (2023). https://doi.org/10.1126/science.adg7492"),
            "note": (
                "Ground-truth variants come from real TCGA-LUAD MAF data; "
                "plasma reads are simulated by Poisson sampling. Panel LLR is "
                "computed with AM pathogenicity weights loaded from a local "
                "index if present, otherwise a per-variant-class prior "
                "(see am_primary_source). AM weights are applied as fixed "
                "scalars in a deterministic aggregation — not as features to "
                "a learned model."
            ),
        },
        "weight_info": info,
        "results": out["results"],
        "elapsed_seconds": elapsed,
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"\n  Results saved to {args.output}")
    print(f"  Wall time: {elapsed:.1f}s")

    # Compact comparison table
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

    # Bottom-line deltas
    print("\n  Bottom-line Delta vs panel_llr_uniform (AUC @ 0.1%):")
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
                flag = "OK" if d >= 0.02 else ("~" if d >= 0 else "X")
                print(f"    {flag}  {m:>22}: AUC {row['mean']:.4f}  Delta={d:+.4f}")

    return 0


if __name__ == "__main__":
    sys.exit(main())