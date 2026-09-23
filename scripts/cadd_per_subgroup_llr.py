"""
Per-subgroup CADD panel LLR for DeepCatch v2.2.

Honest framing:
  - The 5,738 TCGA mutations are all from 20 TCGA-LUAD patients. There is no
    OV / PAAD data here, so "per-cancer" panels are biologically meaningless
    for this cohort.
  - The actual subgroups available are:
      * TP53-mutant vs TP53-wildtype (12 vs 8 patients)
      * KRAS-mutant vs KRAS-wildtype (4 vs 16 patients)
      * STK11-mutant vs STK11-wildtype (5 vs 15 patients)
      * Top-half mutation burden vs bottom-half (10 vs 10)
  - For each subgroup we measure:
      * Baseline uniform LLR per-subgroup AUC + Sens@99% (at 0.1% ctDNA)
      * CADD Top-K=500 per-subgroup LLR (selecting high-CADD mutations
        within that subgroup's patient set only)
  - We also test a *driver-only* panel — restricted to mutations in the LUAD
    driver set (TP53, KRAS, EGFR, STK11, KEAP1, CDKN2A, SMARCA4, NKX2-1) —
    on the full 20-patient cohort, to compare against the existing Top-K=500
    finding.

The honest answer this script aims to surface: does a *per-subgroup* CADD
panel lift Sens@99% (at 0.1% ctDNA, AUC ≥ 0.921) compared to the published
whole-cohort Top-K=500 result (Sens@99% = 0.64)?

If subgroup panel selection beats the whole-cohort 0.64 within its subgroup,
that's a real per-subgroup lift (and the per-subgroup uniform LLR is the
proper baseline, not the whole-cohort number).

Outputs:
  - /Users/hermes/deepcatch/results/cadd_per_subgroup_llr.json
  - /Users/hermes/deepcatch/docs/CADD_PER_SUBGROUP_LLR.md
"""
import json
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from sklearn.metrics import roc_auc_score

ROOT = Path("/Users/hermes/deepcatch")
sys.path.insert(0, str(ROOT))

from real_tcga_validation import (  # noqa: E402
    compute_llr_scores,
    simulate_cfdna_from_real,
    _fisher_scores,
    _panel_metrics,
    sensitivity_at_specificity,
)
from cadd_weighted_llr import (  # noqa: E402
    load_or_match_caddings,
    build_per_patient_weights,
    build_topk_per_patient_weights,
)

COHORT_MUTATIONS = ROOT / "results" / "cadd_mutations_full.json"
RESULTS_PATH = ROOT / "results" / "cadd_per_subgroup_llr.json"
DOCS_PATH = ROOT / "docs" / "CADD_PER_SUBGROUP_LLR.md"

LUAD_DRIVERS = [
    "TP53", "KRAS", "EGFR", "STK11", "KEAP1",
    "CDKN2A", "SMARCA4", "NKX2-1",
]

SEEDS = [42, 123, 456, 789, 1024]
TUMOR_FRACTIONS = [0.1, 0.05, 0.01, 0.005, 0.001]


# ─────────────────────────────────────────────────────────────────────────────
# Subgroup definitions
# ─────────────────────────────────────────────────────────────────────────────
def define_subgroups(cohort: Dict[str, List[Dict]]) -> Dict[str, List[str]]:
    """Return subgroup_name -> list of patient IDs."""
    patients = list(cohort.keys())
    mut_burden = [(p, len(cohort[p])) for p in patients]
    mut_burden_sorted = sorted(mut_burden, key=lambda x: x[1])
    n = len(patients)
    half = n // 2
    low_half = sorted([p for p, _ in mut_burden_sorted[:half]])
    high_half = sorted([p for p, _ in mut_burden_sorted[half:]])

    def has_gene(p: str, gene: str) -> bool:
        return any(m["gene"] == gene for m in cohort[p])

    return {
        "all_20_patients": patients,
        "TP53_mutant":  sorted([p for p in patients if has_gene(p, "TP53")]),
        "TP53_wildtype": sorted([p for p in patients if not has_gene(p, "TP53")]),
        "KRAS_mutant":  sorted([p for p in patients if has_gene(p, "KRAS")]),
        "KRAS_wildtype": sorted([p for p in patients if not has_gene(p, "KRAS")]),
        "STK11_mutant": sorted([p for p in patients if has_gene(p, "STK11")]),
        "STK11_wildtype": sorted([p for p in patients if not has_gene(p, "STK11")]),
        "high_burden_top10": high_half,
        "low_burden_bottom10": low_half,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Per-patient weight builders restricted to a subgroup
# ─────────────────────────────────────────────────────────────────────────────
def build_uniform_weights_for_subgroup(
    subgroup_patients: List[str],
    cohort: Dict[str, List[Dict]],
) -> Dict[str, np.ndarray]:
    """Standard uniform LLR: weight = 1 for every mutation of every patient in
    the subgroup, 0 (drop) for patients outside the subgroup."""
    weights = {}
    for p in subgroup_patients:
        weights[p] = np.ones(len(cohort[p]))
    return weights


def build_topk_weights_for_subgroup(
    subgroup_patients: List[str],
    cohort: Dict[str, List[Dict]],
    matches: Dict[str, Any],
    top_k: int,
    min_phred: float = 0.0,
) -> Dict[str, np.ndarray]:
    """Top-K CADD weights *restricted to the subgroup's patients*. Patients
    outside the subgroup are absent from the dict (drop them from the ROC).
    """
    weights_all = build_topk_per_patient_weights(
        cohort, matches, top_k=top_k, min_phred=min_phred
    )
    return {p: weights_all[p] for p in subgroup_patients if p in weights_all}


def build_driver_only_weights(
    cohort: Dict[str, List[Dict]],
    drivers: List[str],
) -> Dict[str, np.ndarray]:
    """Weight = 1 for mutations in driver genes, 0 otherwise."""
    driver_set = set(drivers)
    weights = {}
    for p, muts in cohort.items():
        w = np.array([1.0 if m["gene"] in driver_set else 0.0 for m in muts])
        weights[p] = w
    return weights


# ─────────────────────────────────────────────────────────────────────────────
# Subgroup panel detection — uses the same per-seed ROC pipeline as
# cadd_weighted_llr.run_weighted_panel_detection but operates on the supplied
# weight dict (so we can restrict to a patient subset).
# ─────────────────────────────────────────────────────────────────────────────
def run_subgroup_panel_detection(
    cohort: Dict[str, List[Dict]],
    weights_per_patient: Dict[str, np.ndarray],
    subgroup_patients: List[str],
    tumor_fractions: Optional[List[float]] = None,
    seeds: Optional[List[int]] = None,
    cfdna_depth: int = 5000,
    bg_error_rate: float = 0.002,
) -> Dict[str, List[Dict]]:
    """Per-subgroup panel detection. The ROC is computed *within* the
    subgroup (positives and negatives come from the same patients), so the
    sens_at_99_spec is a within-subgroup statistic.

    Drops any patient whose weights sum to 0 (no panel loci left after
    filtering) — flagged in the returned metadata.
    """
    if tumor_fractions is None:
        tumor_fractions = TUMOR_FRACTIONS
    if seeds is None:
        seeds = SEEDS

    # Drop zero-weight patients up-front and record
    used = [p for p in subgroup_patients
            if p in weights_per_patient and weights_per_patient[p].sum() > 0]
    dropped = [p for p in subgroup_patients if p not in used]

    results = {"panel_llr_weighted": []}

    for tf in tumor_fractions:
        print(f"      TF={tf*100:.2f}%  ({len(used)} patients × {len(seeds)} seeds"
              + (f", dropped {dropped}" if dropped else "") + ")")
        per_seed = {}
        for seed in seeds:
            pos_w, neg_w = [], []
            for patient in used:
                muts = cohort[patient]
                weights = weights_per_patient[patient]
                dp = simulate_cfdna_from_real(muts, tumor_fraction=tf,
                                              cfdna_depth=cfdna_depth, seed=seed,
                                              bg_error_rate=bg_error_rate)
                dn = simulate_cfdna_from_real(muts, tumor_fraction=0.0,
                                              cfdna_depth=cfdna_depth, seed=seed,
                                              bg_error_rate=bg_error_rate)
                lp = compute_llr_scores(dp["depths"], dp["X"][:, 1].astype(int), dp["X"][:, 3])
                ln = compute_llr_scores(dn["depths"], dn["X"][:, 1].astype(int), dn["X"][:, 3])
                nv_p, nv_n = dp["n_variants"], dn["n_variants"]
                panel_size = min(nv_p, nv_n)
                wp = weights[:panel_size]
                wn = weights[:nv_n]
                pos_w.append(float((lp[:panel_size] * wp).sum()))
                neg_w.append(float((ln[:panel_size] * wn).sum()))
            y = np.array([1] * len(pos_w) + [0] * len(neg_w))
            per_seed[seed] = _panel_metrics(y, np.array(pos_w + neg_w))
        for m in ("auc", "sens_at_95_spec", "sens_at_99_spec", "paired_win_rate"):
            vals = [per_seed[s][m] for s in seeds]
            results["panel_llr_weighted"].append({
                "tumor_fraction": tf,
                "metric": m,
                "mean": float(np.mean(vals)),
                "std": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                "per_seed": {str(s): per_seed[s][m] for s in seeds},
                "n_patients": len(used),
                "dropped_patients": dropped,
            })
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Helpers for reporting
# ─────────────────────────────────────────────────────────────────────────────
def sens99_at(results: Dict[str, List[Dict]], tf: float = 0.001) -> Dict[str, float]:
    rows = [r for r in results["panel_llr_weighted"]
            if r["tumor_fraction"] == tf and r["metric"] == "sens_at_99_spec"]
    if not rows:
        return {"mean": float("nan"), "std": float("nan")}
    return {"mean": rows[0]["mean"], "std": rows[0]["std"]}


def auc_at(results: Dict[str, List[Dict]], tf: float = 0.001) -> Dict[str, float]:
    rows = [r for r in results["panel_llr_weighted"]
            if r["tumor_fraction"] == tf and r["metric"] == "auc"]
    if not rows:
        return {"mean": float("nan"), "std": float("nan")}
    return {"mean": rows[0]["mean"], "std": rows[0]["std"]}


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("=" * 78)
    print("Per-subgroup CADD panel LLR for DeepCatch v2.2")
    print("=" * 78)

    print("\n[1] Loading cohort and CADD matches...")
    cohort_data = json.load(open(COHORT_MUTATIONS))
    cohort = cohort_data["cohort_20_patients"]
    n_muts = sum(len(v) for v in cohort.values())
    print(f"  Cohort: {len(cohort)} patients, {n_muts} mutations")

    matches = load_or_match_caddings()
    n_matched = len(matches["snv_matched"]) + len(matches["indel_matched"])
    print(f"  CADD matches: {n_matched}/{n_muts} ({n_matched/n_muts*100:.1f}%)")

    subgroups = define_subgroups(cohort)
    print("\n[2] Subgroups defined:")
    for k, v in subgroups.items():
        print(f"  {k}: n={len(v)}")

    # Driver-mutation-only panel (whole cohort, no per-subgroup split)
    print("\n[3] Driver-only panel (LUAD driver gene set, whole cohort)...")
    weights_driver = build_driver_only_weights(cohort, LUAD_DRIVERS)
    n_driver_muts = sum(int(w.sum()) for w in weights_driver.values())
    print(f"  Driver-only mutations kept: {n_driver_muts}/{n_muts} "
          f"({n_driver_muts/n_muts*100:.1f}%)")
    driver_genes_seen = sorted({
        m["gene"]
        for p, muts in cohort.items() for m in muts
        if m["gene"] in set(LUAD_DRIVERS)
    })
    print(f"  Driver genes with ≥1 mutation in cohort: {driver_genes_seen}")
    res_driver = run_subgroup_panel_detection(
        cohort, weights_driver, list(cohort.keys()),
        tumor_fractions=[0.001, 0.005, 0.01, 0.05, 0.1],
    )
    sens99_driver = sens99_at(res_driver, 0.001)
    auc_driver = auc_at(res_driver, 0.001)
    print(f"  Driver-only @ 0.1% ctDNA: AUC={auc_driver['mean']:.4f}±{auc_driver['std']:.4f} "
          f"  Sens@99%={sens99_driver['mean']:.3f}±{sens99_driver['std']:.3f}")

    # Per-subgroup: uniform vs CADD Top-K=500 (within subgroup)
    out = {
        "research_use_only": True,
        "experiment": "Per-subgroup CADD panel LLR for DeepCatch v2.2",
        "cohort": {"n_patients": len(cohort), "n_mutations": n_muts},
        "cadd_match_rate": n_matched / n_muts,
        "seeds": SEEDS,
        "tumor_fractions": TUMOR_FRACTIONS,
        "luad_drivers_defined": LUAD_DRIVERS,
        "luad_drivers_observed": driver_genes_seen,
        "driver_only_panel_results": res_driver,
        "subgroups": {},
    }

    print("\n[4] Per-subgroup panel comparisons at 0.1% ctDNA...")
    for sg_name, sg_patients in subgroups.items():
        if sg_name == "all_20_patients":
            continue
        if len(sg_patients) < 2:
            print(f"  -- {sg_name}: too few patients (n={len(sg_patients)}), skip")
            continue
        print(f"\n  -- {sg_name} (n={len(sg_patients)})")

        weights_uniform = build_uniform_weights_for_subgroup(sg_patients, cohort)
        res_uniform = run_subgroup_panel_detection(cohort, weights_uniform, sg_patients)

        weights_top500 = build_topk_weights_for_subgroup(
            sg_patients, cohort, matches, top_k=500
        )
        res_top500 = run_subgroup_panel_detection(cohort, weights_top500, sg_patients)

        # Top-K=200 within subgroup (driver context — tighter sub-panel)
        weights_top200 = build_topk_weights_for_subgroup(
            sg_patients, cohort, matches, top_k=200
        )
        res_top200 = run_subgroup_panel_detection(cohort, weights_top200, sg_patients)

        s_u = sens99_at(res_uniform, 0.001)
        s_t5 = sens99_at(res_top500, 0.001)
        s_t2 = sens99_at(res_top200, 0.001)
        a_u = auc_at(res_uniform, 0.001)
        a_t5 = auc_at(res_top500, 0.001)
        a_t2 = auc_at(res_top200, 0.001)
        print(f"     uniform LLR       AUC={a_u['mean']:.4f}±{a_u['std']:.4f}  "
              f"Sens@99%={s_u['mean']:.3f}±{s_u['std']:.3f}")
        print(f"     CADD Top-K=500    AUC={a_t5['mean']:.4f}±{a_t5['std']:.4f}  "
              f"Sens@99%={s_t5['mean']:.3f}±{s_t5['std']:.3f}")
        print(f"     CADD Top-K=200    AUC={a_t2['mean']:.4f}±{a_t2['std']:.4f}  "
              f"Sens@99%={s_t2['mean']:.3f}±{s_t2['std']:.3f}")

        out["subgroups"][sg_name] = {
            "n_patients": len(sg_patients),
            "patients": sg_patients,
            "uniform_llr": res_uniform,
            "cadd_topk_500": res_top500,
            "cadd_topk_200": res_top200,
        }

    # Also run uniform vs Top-K=500 on the WHOLE cohort as the head-to-head
    # comparison anchor (the published number is 0.64 Sens@99%).
    print("\n[5] Whole-cohort anchor: uniform LLR vs CADD Top-K=500...")
    weights_uniform_all = build_uniform_weights_for_subgroup(list(cohort.keys()), cohort)
    res_uniform_all = run_subgroup_panel_detection(cohort, weights_uniform_all, list(cohort.keys()))
    weights_top500_all = build_topk_weights_for_subgroup(list(cohort.keys()), cohort, matches, top_k=500)
    res_top500_all = run_subgroup_panel_detection(cohort, weights_top500_all, list(cohort.keys()))
    s_u_all = sens99_at(res_uniform_all, 0.001)
    s_t5_all = sens99_at(res_top500_all, 0.001)
    a_u_all = auc_at(res_uniform_all, 0.001)
    a_t5_all = auc_at(res_top500_all, 0.001)
    print(f"  whole-cohort uniform LLR: AUC={a_u_all['mean']:.4f}±{a_u_all['std']:.4f}  "
          f"Sens@99%={s_u_all['mean']:.3f}±{s_u_all['std']:.3f}")
    print(f"  whole-cohort CADD TopK500: AUC={a_t5_all['mean']:.4f}±{a_t5_all['std']:.4f}  "
          f"Sens@99%={s_t5_all['mean']:.3f}±{s_t5_all['std']:.3f}")

    out["whole_cohort_anchor"] = {
        "uniform_llr": res_uniform_all,
        "cadd_topk_500": res_top500_all,
    }

    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
