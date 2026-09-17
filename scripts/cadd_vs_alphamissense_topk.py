"""
Direct apples-to-apples CADD vs AlphaMissense Top-K panel comparison
on the SAME 20-patient TCGA-LUAD cohort.

Hypothesis: AlphaMissense (Cheng 2023, Science) provides ~2x the
coverage of CADD on private/rare somatic mutations (CADD only scores
gnomAD-observed variants; AM scores the full human missome via the
canonical transcript), and the per-mutation AM pathogenicity score is
a stronger per-locus weight than CADD PHRED for panel selection.

Both runners (cadd_weighted_llr.build_topk_per_patient_weights and
alphamissense_panel_run._panel_score_topk_select) use the SAME
apples-to-apples operation:
   1. Per patient, rank all mutations by the score.
   2. Keep the top-K.
   3. Sum uniform LLR over those K.

This script wires both score sources into a single comparison harness
so the deltas are attributable only to the scoring function, not the
selection logic.

Writes:
    results/cadd_vs_alphamissense_topk_20patient.json
"""

from __future__ import annotations

import argparse
import json
import logging
import statistics
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/Users/hermes/deepcatch")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from real_tcga_validation import (  # noqa: E402
    compute_llr_scores,
    simulate_cfdna_from_real,
    _panel_metrics,
)
from cadd_weighted_llr import build_topk_per_patient_weights  # noqa: E402
import alphamissense_weights as amw  # noqa: E402


LOG = logging.getLogger("cadd_vs_am_topk")


# ─────────────────────────────────────────────────────────────────────────────
# Per-score source selection → {patient_id → weights array}
# ─────────────────────────────────────────────────────────────────────────────
def build_cadd_topk_weights(cohort, cadd_matches, top_k):
    """Reuse the CADD runner's per-patient Top-K selector."""
    return build_topk_per_patient_weights(cohort, cadd_matches, top_k=top_k)


def build_am_topk_select_weights(cohort, cohort_mutations, am_weights, top_k):
    """Build per-patient AM Top-K selector aligned to the cohort's mutation order."""
    patients = list(cohort["patients"].keys())
    weights_per_patient = {}
    matched_per_patient = {}
    for p in patients:
        mlist = cohort["patients"][p]
        scored = []
        for i, m in enumerate(mlist):
            if m.get("variant_class") != "Missense_Mutation":
                continue
            ref = m.get("ref", "")
            alt = m.get("alt", "")
            if not ref or not alt:
                for mr in cohort_mutations:
                    if (mr["sample"] == m.get("sample")
                            and mr["chrom"] == m.get("chrom")
                            and mr["pos"] == m.get("pos")
                            and mr["gene"] == m.get("gene")):
                        ref, alt = mr["ref"], mr["alt"]
                        break
            if len(ref) != 1 or len(alt) != 1:
                continue
            # Find the AM-keyed lookup
            for mr in cohort_mutations:
                if (mr["sample"] == m.get("sample")
                        and mr["chrom"] == m.get("chrom")
                        and mr["pos"] == m.get("pos")
                        and mr["gene"] == m.get("gene")
                        and mr["ref"] == ref and mr["alt"] == alt):
                    ref_aa = mr.get("ref_aa") or ref
                    alt_aa = mr.get("alt_aa") or alt
                    key = amw.am_key(mr["uniprot"], mr["prot_pos"], ref_aa, alt_aa)
                    break
            else:
                continue
            w = am_weights.get(key)
            if w is None:
                continue
            scored.append((i, w.score))
        scored.sort(key=lambda x: -x[1])
        kept = {i for i, _ in scored[:top_k]}
        w = np.array([1.0 if i in kept else 0.0 for i in range(len(mlist))])
        weights_per_patient[p] = w
        matched_per_patient[p] = len(kept)
    return weights_per_patient, matched_per_patient


def build_am_topk_select_weights_from_dict(cohort_dict, cohort_mutations, am_weights, top_k):
    """Same as `build_am_topk_select_weights` but takes a {patient: [mut_dict]} dict.

    The CADD runner's `build_topk_per_patient_weights` uses this dict
    format directly, so we use the same here to guarantee the per-patient
    mutation ordering is identical between CADD and AM Top-K runs.
    """
    weights_per_patient = {}
    matched_per_patient = {}
    for p, mlist in cohort_dict.items():
        scored = []
        for i, m in enumerate(mlist):
            if m.get("variant_class") != "Missense_Mutation":
                continue
            ref = m.get("ref", "")
            alt = m.get("alt", "")
            if len(ref) != 1 or len(alt) != 1:
                continue
            # Find the AM-keyed lookup
            for mr in cohort_mutations:
                if (mr["sample"] == m.get("sample")
                        and mr["chrom"] == m.get("chrom")
                        and mr["pos"] == m.get("pos")
                        and mr["gene"] == m.get("gene")
                        and mr["ref"] == ref and mr["alt"] == alt):
                    ref_aa = mr.get("ref_aa") or ref
                    alt_aa = mr.get("alt_aa") or alt
                    key = amw.am_key(mr["uniprot"], mr["prot_pos"], ref_aa, alt_aa)
                    break
            else:
                continue
            w = am_weights.get(key)
            if w is None:
                continue
            scored.append((i, w.score))
        scored.sort(key=lambda x: -x[1])
        kept = {i for i, _ in scored[:top_k]}
        w = np.array([1.0 if i in kept else 0.0 for i in range(len(mlist))])
        weights_per_patient[p] = w
        matched_per_patient[p] = len(kept)
    return weights_per_patient, matched_per_patient


# ─────────────────────────────────────────────────────────────────────────────
# Panel detection runner (mirrors cadd_weighted_llr.run_weighted_panel_detection)
# ─────────────────────────────────────────────────────────────────────────────
def run_panel_detection(cohort, weights_per_patient, tumor_fractions, seeds,
                        cfdna_depth=5000, bg_error_rate=0.002):
    """`cohort` may be a {patient: [mut]} dict (CADD-style) or a real_tcga_validation-style
    cohort object with a 'patients' key. Both forms are accepted.
    """
    if isinstance(cohort, dict) and "patients" in cohort:
        patients = list(cohort["patients"].keys())
        cohort_muts = cohort["patients"]
    else:
        patients = list(cohort.keys())
        cohort_muts = cohort
    results = {"panel_llr_weighted": []}
    for tf in tumor_fractions:
        per_seed = {}
        for seed in seeds:
            pos_w, neg_w = [], []
            for patient in patients:
                muts = cohort_muts[patient]
                weights = weights_per_patient[patient]
                dp = simulate_cfdna_from_real(
                    muts, tumor_fraction=tf, cfdna_depth=cfdna_depth,
                    seed=seed, bg_error_rate=bg_error_rate,
                )
                dn = simulate_cfdna_from_real(
                    muts, tumor_fraction=0.0, cfdna_depth=cfdna_depth,
                    seed=seed, bg_error_rate=bg_error_rate,
                )
                lp = compute_llr_scores(dp['depths'], dp['X'][:, 1].astype(int), dp['X'][:, 3])
                ln = compute_llr_scores(dn['depths'], dn['X'][:, 1].astype(int), dn['X'][:, 3])
                nv_p, nv_n = dp['n_variants'], dn['n_variants']
                panel_size = min(nv_p, nv_n)
                wp = weights[:panel_size]
                wn = weights[:nv_n]
                pos_w.append(float((lp[:panel_size] * wp).sum()))
                neg_w.append(float((ln[:panel_size] * wn).sum()))
            y = np.array([1] * len(pos_w) + [0] * len(neg_w))
            per_seed[seed] = _panel_metrics(y, np.array(pos_w + neg_w))
        for m in ('auc', 'sens_at_95_spec', 'sens_at_99_spec', 'paired_win_rate'):
            vals = [per_seed[s][m] for s in seeds]
            results["panel_llr_weighted"].append({
                'tumor_fraction': tf,
                'metric': m,
                'mean': float(np.mean(vals)),
                'std': float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                'per_seed': {str(s): per_seed[s][m] for s in seeds},
            })
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--am-index",
                        default=str(ROOT / "results/alphamissense_real_index_full.pkl"))
    parser.add_argument("--cadd-matches",
                        default=str(ROOT / "results/cadd_matches_augmented.json"))
    parser.add_argument("--cohort-mutations",
                        default=str(ROOT / "results/cadd_mutations_full.json"))
    parser.add_argument("--output",
                        default=str(ROOT / "results/cadd_vs_alphamissense_topk_20patient.json"))
    parser.add_argument("--n-patients", type=int, default=20)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--topk-values", default="20,200,500")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    print("=" * 80)
    print("  Apples-to-apples CADD vs AlphaMissense Top-K panel LLR on 20-pt cohort")
    print("=" * 80)

    # Use cadd_mutations_full.json as the single source of truth for the
    # 20-patient cohort (matches CADD matches, has ref/alt/variant_class).
    cohort_data = json.load(open(args.cohort_mutations))
    cohort_dict_for_cadd = {
        p: [{**m, "chrom": m.get("chrom", m.get("chromosome", ""))}
            for m in cohort_data["cohort_20_patients"][p]]
        for p in cohort_data["cohort_20_patients"]
    }
    n_patients = len(cohort_dict_for_cadd)
    n_muts = sum(len(v) for v in cohort_dict_for_cadd.values())
    print(f"  Cohort: {n_patients} patients, {n_muts} mutations")

    # CADD matches
    cadd_matches = json.load(open(args.cadd_matches))
    n_cadd = len(cadd_matches["snv_matched"]) + len(cadd_matches["indel_matched"])
    print(f"  CADD matched: {n_cadd}/{n_muts} ({n_cadd/n_muts*100:.1f}%)")

    # Augmented pool + AM weights
    aug = amw.load_tcga_mutations_with_protein(str(ROOT / "validation/tcga/tcga_cache"))
    am_w, info = amw.weight_cohort(aug, index_path=Path(args.am_index))
    n_am = sum(1 for v in am_w.values() if v.source in ("picklelocal", "tsvlocal"))
    print(f"  AM primary source: {info['primary_source']}, n_real={info['n_real']}, n_proxy={info['n_proxy']}")

    tumor_fractions = [0.1, 0.05, 0.01, 0.005, 0.001]
    seeds = [42, 123, 456, 789, 1024][: args.seeds]
    topk_values = [int(x) for x in args.topk_values.split(",")]

    # Uniform baseline
    print("\n[1] Uniform baseline...")
    weights_uniform = {p: np.ones(len(cohort_dict_for_cadd[p])) for p in cohort_dict_for_cadd}
    res_uniform = run_panel_detection(cohort_dict_for_cadd, weights_uniform, tumor_fractions, seeds)

    schemes = {"uniform": res_uniform}

    # CADD Top-K per K
    for k in topk_values:
        print(f"\n[2] CADD Top-K={k} selection...")
        w = build_cadd_topk_weights(cohort_dict_for_cadd, cadd_matches, top_k=k)
        kept = sum(int((ww > 0).sum()) for ww in w.values())
        print(f"  total kept: {kept} (avg per patient {kept/n_patients:.1f})")
        schemes[f"cadd_topk_{k}"] = run_panel_detection(
            cohort_dict_for_cadd, w, tumor_fractions, seeds)

    # AM Top-K per K
    # Use the same cohort_dict_for_cadd so the mutation ordering is identical
    # between CADD and AM Top-K runs.
    for k in topk_values:
        print(f"\n[3] AM Top-K={k} selection...")
        w, kept_per_patient = build_am_topk_select_weights_from_dict(
            cohort_dict_for_cadd, aug, am_w, top_k=k,
        )
        kept = sum(kept_per_patient.values())
        print(f"  total kept: {kept} (avg per patient {kept/n_patients:.1f})")
        schemes[f"am_topk_{k}"] = run_panel_detection(
            cohort_dict_for_cadd, w, tumor_fractions, seeds)

    # Headline at 0.1%
    print("\n" + "=" * 80)
    print("  Headline @ 0.1% ctDNA  (mean AUC ± std)")
    print("=" * 80)
    headline = {}
    for name, res in schemes.items():
        row = next((r for r in res["panel_llr_weighted"]
                    if r["tumor_fraction"] == 0.001 and r["metric"] == "auc"), None)
        if row:
            print(f"    {name:>20}: {row['mean']:.4f} ± {row['std']:.4f}")
            headline[name] = {"mean": row["mean"], "std": row["std"]}

    out = {
        "experiment": "CADD vs AlphaMissense Top-K panel selection on 20-patient TCGA-LUAD cohort",
        "cohort": {"n_patients": n_patients, "n_mutations": n_muts},
        "match_rates": {
            "cadd": {"n_matched": n_cadd, "match_rate": n_cadd / n_muts},
            "alphamissense": {
                "primary_source": info["primary_source"],
                "n_real": info["n_real"],
                "n_proxy": info["n_proxy"],
                "augmented_pool_n_missense": info["missense_total"],
            },
        },
        "tumor_fractions": tumor_fractions,
        "seeds": seeds,
        "topk_values": topk_values,
        "results": {name: res for name, res in schemes.items()},
        "headline_at_0_1pct": headline,
        "license": (
            "CADD: CC BY-NC-SA 4.0 (non-commercial). "
            "AlphaMissense: CC BY-NC-SA 4.0 (non-commercial)."
        ),
        "references": [
            "Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. "
            "A general framework for estimating the relative pathogenicity of human "
            "genetic variants. Nat Genet. 2014;46(3):310-315. doi:10.1038/ng.2892",
            "Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense "
            "variant effect prediction with AlphaMissense. Science 381, eadg7492 "
            "(2023). doi:10.1126/science.adg7492",
        ],
        "note": (
            "All schemes use the SAME per-patient Top-K SELECTION logic "
            "(rank by score, keep top-K, sum uniform LLR over them). The only "
            "difference is the scoring function. CADD matches 4,882/19,421 mutations "
            "via the gnomAD-only TSV; AlphaMissense matches 15,462/15,903 missense "
            "SNVs (~97%) via the canonical-transcript AlphaMissense TSV. The "
            "asymmetric match rates are a key limitation of CADD on private/rare "
            "somatic mutations and motivate AlphaMissense-based panel selection."
        ),
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults saved to {args.output}")


if __name__ == "__main__":
    main()