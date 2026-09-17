"""CADD Top-K per-subgroup validation on the GDC TCGA-LUAD cohort (382 patients).

Honest framing:
- New GDC cohort has 382 patients × 124k mutations vs original 20-patient
  cohort's 5,738 mutations (19x patient count, 22x mutations).
- CADD tabix match rate: 10,922/123,162 SNVs (8.9%) — much lower than the
  original cohort's 86% because GDC bulk WXS includes many passenger
  mutations not in gnomAD r3.0 (CADD only covers observed variants).
- Per-patient median matched mutations: 21. Top-K=200 is infeasible for
  most patients (only 1 has >=200 matches). We test Top-K=20 (median) and
  Top-K=50, plus a driver-restricted panel as alternatives.

Subgroups (defined on the 382-patient cohort):
- TP53 mutant/wildtype (195 / 187)
- KRAS mutant/wildtype (112 / 270)
- STK11 mutant/wildtype (51 / 331)
- top/bottom half mutation burden (191 / 191)
"""
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path("/Users/hermes/deepcatch")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from real_tcga_validation import (  # noqa: E402
    compute_llr_scores, simulate_cfdna_from_real, _panel_metrics,
)
from cadd_weighted_llr import build_topk_per_patient_weights  # noqa: E402

COHORT_PATH = ROOT / "results" / "cadd_mutations_gdc_validation.json"
MATCHES_PATH = ROOT / "results" / "cadd_matches_gdc_validation.json"
RESULTS_OUT = ROOT / "results" / "cadd_per_subgroup_llr_GDC_VALIDATION.json"
LUAD_DRIVERS = ["TP53", "KRAS", "EGFR", "STK11", "KEAP1",
                "CDKN2A", "SMARCA4", "NKX2-1"]
SEEDS = [42, 123, 456]  # 3 seeds for tractable runtime
TUMOR_FRACTIONS = [0.001]  # 0.1% ctDNA — the headline metric


def define_subgroups(cohort):
    patients = sorted(cohort.keys())
    mb = sorted([(p, len(cohort[p])) for p in patients], key=lambda x: x[1])
    n = len(patients)
    half = n // 2
    low_half = sorted([p for p, _ in mb[:half]])
    high_half = sorted([p for p, _ in mb[half:]])

    def has_gene(p, gene):
        return any(m["gene"] == gene for m in cohort[p])

    return {
        "all_patients": patients,
        "TP53_mutant":  sorted([p for p in patients if has_gene(p, "TP53")]),
        "TP53_wildtype": sorted([p for p in patients if not has_gene(p, "TP53")]),
        "KRAS_mutant":  sorted([p for p in patients if has_gene(p, "KRAS")]),
        "KRAS_wildtype": sorted([p for p in patients if not has_gene(p, "KRAS")]),
        "STK11_mutant": sorted([p for p in patients if has_gene(p, "STK11")]),
        "STK11_wildtype": sorted([p for p in patients if not has_gene(p, "STK11")]),
        "high_burden_top_half": high_half,
        "low_burden_bottom_half": low_half,
    }


def build_uniform_weights(patients, cohort):
    return {p: np.ones(len(cohort[p])) for p in patients}


def build_topk_weights(patients, cohort, matches, top_k):
    weights_all = build_topk_per_patient_weights(cohort, matches, top_k=top_k)
    return {p: weights_all[p] for p in patients if p in weights_all}


def build_driver_only_weights(cohort, drivers):
    s = set(drivers)
    w = {}
    for p, muts in cohort.items():
        w[p] = np.array([1.0 if m["gene"] in s else 0.0 for m in muts])
    return w


def _simulate_patient(args):
    """Helper for parallel execution: simulate one patient, return (pos_score, neg_score)."""
    patient, cohort, weights, tf, cfdna_depth, seed, bg_error_rate = args
    muts = cohort[patient]
    weights_arr = weights[patient]
    dp = simulate_cfdna_from_real(
        muts, tumor_fraction=tf, cfdna_depth=cfdna_depth,
        seed=seed, bg_error_rate=bg_error_rate,
    )
    dn = simulate_cfdna_from_real(
        muts, tumor_fraction=0.0, cfdna_depth=cfdna_depth,
        seed=seed, bg_error_rate=bg_error_rate,
    )
    lp = compute_llr_scores(dp["depths"], dp["X"][:, 1].astype(int), dp["X"][:, 3])
    ln = compute_llr_scores(dn["depths"], dn["X"][:, 1].astype(int), dn["X"][:, 3])
    nv_p, nv_n = dp["n_variants"], dn["n_variants"]
    panel_size = min(nv_p, nv_n)
    wp = weights_arr[:panel_size]
    wn = weights_arr[:nv_n]
    return (
        float((lp[:panel_size] * wp).sum()),
        float((ln[:panel_size] * wn).sum()),
    )


def run_panel(cohort, weights_per_patient, subgroup_patients,
              tumor_fractions=TUMOR_FRACTIONS, seeds=SEEDS,
              cfdna_depth=5000, bg_error_rate=0.002, n_workers=8):
    from concurrent.futures import ProcessPoolExecutor
    used = [p for p in subgroup_patients
            if p in weights_per_patient and weights_per_patient[p].sum() > 0]
    dropped = [p for p in subgroup_patients if p not in used]
    results = {"panel_llr_weighted": []}
    for tf in tumor_fractions:
        per_seed = {}
        for seed in seeds:
            args_list = [(p, cohort, weights_per_patient, tf, cfdna_depth, seed, bg_error_rate)
                         for p in used]
            pos_w, neg_w = [], []
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                for ps, ns in ex.map(_simulate_patient, args_list, chunksize=4):
                    pos_w.append(ps)
                    neg_w.append(ns)
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


def sens99_at(results, tf=0.001):
    rows = [r for r in results["panel_llr_weighted"]
            if r["tumor_fraction"] == tf and r["metric"] == "sens_at_99_spec"]
    if not rows:
        return {"mean": float("nan"), "std": float("nan")}
    return {"mean": rows[0]["mean"], "std": rows[0]["std"]}


def auc_at(results, tf=0.001):
    rows = [r for r in results["panel_llr_weighted"]
            if r["tumor_fraction"] == tf and r["metric"] == "auc"]
    if not rows:
        return {"mean": float("nan"), "std": float("nan")}
    return {"mean": rows[0]["mean"], "std": rows[0]["std"]}


def main():
    print("=" * 78)
    print("CADD per-subgroup VALIDATION on GDC TCGA-LUAD (~382 patients)")
    print("=" * 78)

    print("\n[1] Loading GDC cohort and CADD matches...")
    cohort_full = json.load(open(COHORT_PATH))['cohort']
    matches = json.load(open(MATCHES_PATH))
    n_patients_full = len(cohort_full)
    n_muts_full = sum(len(v) for v in cohort_full.values())

    # For tractable runtime, sample a subset if cohort is very large.
    # The task target is "5x the original 20" = 100+ patients.
    TARGET_N_PATIENTS = 150
    if n_patients_full > TARGET_N_PATIENTS:
        import random
        random.seed(42)
        sampled_patients = sorted(random.sample(list(cohort_full.keys()),
                                                 TARGET_N_PATIENTS))
        cohort = {p: cohort_full[p] for p in sampled_patients}
        print(f"  Subsampled {TARGET_N_PATIENTS}/{n_patients_full} patients for runtime")
    else:
        cohort = cohort_full
    n_patients = len(cohort)
    n_muts = sum(len(v) for v in cohort.values())
    n_matched = len(matches['snv_matched'])
    n_total_snv = matches['n_total_snv']
    print(f"  Cohort: {n_patients} patients, {n_muts} mutations")
    print(f"  CADD matched: {n_matched}/{n_total_snv} ({n_matched/n_total_snv*100:.1f}%)")
    print(f"  PHRED stats: {matches['phred_stats']}")

    per_patient_matched = Counter(m['sample'] for m in matches['snv_matched'])
    counts = sorted(per_patient_matched.values())
    n = len(counts)
    print(f"  Per-patient matched mutations:")
    print(f"    min={min(counts)} median={counts[n//2]} mean={sum(counts)/n:.1f} max={max(counts)}")
    print(f"    patients with >=20 matches: {sum(1 for c in counts if c >= 20)}")
    print(f"    patients with >=100 matches: {sum(1 for c in counts if c >= 100)}")

    subgroups = define_subgroups(cohort)
    print("\n[2] Subgroups:")
    for k, v in subgroups.items():
        print(f"  {k}: n={len(v)}")

    out = {
        "experiment": "CADD per-subgroup VALIDATION on GDC TCGA-LUAD",
        "cohort": {"n_patients": n_patients, "n_mutations": n_muts},
        "cadd_match_rate_snv": n_matched / n_total_snv,
        "phred_stats_matched": matches['phred_stats'],
        "per_patient_matched_distribution": {
            "min": min(counts), "median": counts[n//2], "mean": sum(counts)/n,
            "max": max(counts),
            "n_with_ge_20": sum(1 for c in counts if c >= 20),
            "n_with_ge_100": sum(1 for c in counts if c >= 100),
        },
        "seeds": SEEDS,
        "tumor_fractions": TUMOR_FRACTIONS,
        "luad_drivers": LUAD_DRIVERS,
        "subgroups": {},
    }

    # Driver-restricted panel (whole cohort)
    print("\n[3] Driver-only panel (LUAD driver gene set, whole cohort)...")
    weights_driver = build_driver_only_weights(cohort, LUAD_DRIVERS)
    n_driver_muts = sum(int(w.sum()) for w in weights_driver.values())
    n_drivers_in_data = sorted({
        m['gene'] for p, muts in cohort.items() for m in muts
        if m['gene'] in set(LUAD_DRIVERS)
    })
    print(f"  Driver-only mutations: {n_driver_muts}/{n_muts} "
          f"({n_driver_muts/n_muts*100:.1f}%)")
    print(f"  Driver genes seen: {n_drivers_in_data}")
    res_driver = run_panel(cohort, weights_driver, list(cohort.keys()))
    s_d = sens99_at(res_driver, 0.001)
    a_d = auc_at(res_driver, 0.001)
    print(f"  Driver-only @ 0.1% ctDNA: AUC={a_d['mean']:.4f}±{a_d['std']:.4f}  "
          f"Sens@99%={s_d['mean']:.3f}±{s_d['std']:.3f}")
    out["driver_only_whole_cohort"] = {
        "n_mutations_kept": n_driver_muts,
        "n_genes_observed": len(n_drivers_in_data),
        "genes_observed": n_drivers_in_data,
        "results": res_driver,
    }

    print("\n[4] Per-subgroup panel comparisons...")
    for sg_name, sg_patients in subgroups.items():
        if sg_name == "all_patients":
            continue
        if len(sg_patients) < 2:
            print(f"  -- {sg_name}: too few (n={len(sg_patients)}), skip")
            continue
        print(f"\n  -- {sg_name} (n={len(sg_patients)})")

        weights_uniform = build_uniform_weights(sg_patients, cohort)
        res_uniform = run_panel(cohort, weights_uniform, sg_patients)

        weights_k20 = build_topk_weights(sg_patients, cohort, matches, top_k=20)
        res_k20 = run_panel(cohort, weights_k20, sg_patients)

        weights_k50 = build_topk_weights(sg_patients, cohort, matches, top_k=50)
        res_k50 = run_panel(cohort, weights_k50, sg_patients)

        s_u = sens99_at(res_uniform, 0.001)
        s_k20 = sens99_at(res_k20, 0.001)
        s_k50 = sens99_at(res_k50, 0.001)
        a_u = auc_at(res_uniform, 0.001)
        a_k20 = auc_at(res_k20, 0.001)
        a_k50 = auc_at(res_k50, 0.001)
        print(f"     uniform LLR    AUC={a_u['mean']:.4f}  Sens@99%={s_u['mean']:.3f}")
        print(f"     CADD TopK=20   AUC={a_k20['mean']:.4f}  Sens@99%={s_k20['mean']:.3f}")
        print(f"     CADD TopK=50   AUC={a_k50['mean']:.4f}  Sens@99%={s_k50['mean']:.3f}")

        out["subgroups"][sg_name] = {
            "n_patients": len(sg_patients),
            "patients": sg_patients,
            "uniform_llr": res_uniform,
            "cadd_topk_20": res_k20,
            "cadd_topk_50": res_k50,
        }

    print("\n[5] Whole-cohort anchor...")
    weights_uniform_all = build_uniform_weights(list(cohort.keys()), cohort)
    res_uniform_all = run_panel(cohort, weights_uniform_all, list(cohort.keys()))
    weights_k20_all = build_topk_weights(list(cohort.keys()), cohort, matches, top_k=20)
    res_k20_all = run_panel(cohort, weights_k20_all, list(cohort.keys()))
    weights_k50_all = build_topk_weights(list(cohort.keys()), cohort, matches, top_k=50)
    res_k50_all = run_panel(cohort, weights_k50_all, list(cohort.keys()))
    s_u_all = sens99_at(res_uniform_all, 0.001)
    s_k20_all = sens99_at(res_k20_all, 0.001)
    s_k50_all = sens99_at(res_k50_all, 0.001)
    a_u_all = auc_at(res_uniform_all, 0.001)
    a_k20_all = auc_at(res_k20_all, 0.001)
    a_k50_all = auc_at(res_k50_all, 0.001)
    print(f"  whole-cohort uniform: AUC={a_u_all['mean']:.4f}  "
          f"Sens@99%={s_u_all['mean']:.3f}")
    print(f"  whole-cohort TopK=20: AUC={a_k20_all['mean']:.4f}  "
          f"Sens@99%={s_k20_all['mean']:.3f}")
    print(f"  whole-cohort TopK=50: AUC={a_k50_all['mean']:.4f}  "
          f"Sens@99%={s_k50_all['mean']:.3f}")

    out["whole_cohort_anchor"] = {
        "uniform_llr": res_uniform_all,
        "cadd_topk_20": res_k20_all,
        "cadd_topk_50": res_k50_all,
    }

    RESULTS_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_OUT, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults saved to {RESULTS_OUT}")


if __name__ == "__main__":
    main()