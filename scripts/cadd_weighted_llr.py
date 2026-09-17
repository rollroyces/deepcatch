"""
CADD-weighted panel LLR for DeepCatch v2.2.

Hypothesis: weighting each mutation's LLR by its CADD score improves ultra-low
ctDNA detection. CADD (Combined Annotation Dependent Depletion, Kircher et al.
2014) is a public variant impact score (CC BY-NC-SA 4.0, non-commercial).

Citation:
    Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J.
    A general framework for estimating the relative pathogenicity of human
    genetic variants. Nat Genet. 2014;46(3):310-315.
    https://doi.org/10.1038/ng.2892

Data sources (both downloaded once into /Users/hermes/deepcatch/data/cadd/):
  - CADD v1.6 GRCh38 gnomAD genomes r3.0 SNVs (5.9 GB advertised, ~6.35 GB actual):
      https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.snv.tsv.gz
  - CADD v1.7 GRCh38 gnomAD genomes r4.0 InDels (1.2 GB):
      https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz

Match-rate note. The whole-genome SNV TSV (81 GB GRCh38, 78 GB GRCh37) was
beyond the local disk budget (23 GiB free on the bench path), and per-chromosome
files are not provided for v1.7 GRCh38. We therefore used the gnomAD-only TSVs
(observed variants) which are far smaller. Many TCGA-LUAD somatic mutations are
rare/private and do not appear in gnomAD; match rate is honestly reported below.

License: CADD data and scores are CC BY-NC-SA 4.0 (non-commercial). This work
is research use.
"""
import gzip
import json
import os
import statistics
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from sklearn.metrics import roc_auc_score

# Re-use the existing real_tcga_validation infrastructure
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from real_tcga_validation import (  # noqa: E402
    compute_llr_scores,
    simulate_cfdna_from_real,
    _panel_metrics,
    sensitivity_at_specificity,
    _fisher_scores,
)


# ─────────────────────────────────────────────────────────────────────────────
# Paths
# ─────────────────────────────────────────────────────────────────────────────
ROOT = Path("/Users/hermes/deepcatch")
DATA_CADD = ROOT / "data" / "cadd"
SNV_TSV = DATA_CADD / "cadd_v1.6_gnomad_r3_snv.tsv.gz"
SNV_TBI = DATA_CADD / "cadd_v1.6_gnomad_r3_snv.tsv.gz.tbi"
INDEL_TSV = DATA_CADD / "cadd_v1.7_gnomad_r4_indel.tsv.gz"
INDEL_TBI = DATA_CADD / "cadd_v1.7_gnomad_r4_indel.tsv.gz.tbi"
COHORT_MUTATIONS = ROOT / "results" / "cadd_mutations_full.json"
MATCHES_PATH = ROOT / "results" / "cadd_matches.json"
RESULTS_PATH = ROOT / "results" / "cadd_weighted_llr.json"
DOCS_PATH = ROOT / "docs" / "CADD_WEIGHTED_LLR.md"


# ─────────────────────────────────────────────────────────────────────────────
# CADD lookup helpers
# ─────────────────────────────────────────────────────────────────────────────
def _cadd_lookup_one(chrom: str, pos: int, ref: str, alt: str,
                     tsv_path: Path) -> Optional[Tuple[float, float]]:
    """Return (phred, raw) or None if (chrom, pos, ref, alt) is not in the TSV."""
    chrom_num = chrom.replace("chr", "")
    try:
        result = subprocess.run(
            ["tabix", str(tsv_path), f"{chrom_num}:{pos}-{pos}"],
            capture_output=True, text=True, timeout=10,
        )
    except subprocess.TimeoutExpired:
        return None
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 6:
            continue
        if fields[2] == ref and fields[3] == alt:
            try:
                return float(fields[5]), float(fields[4])
            except ValueError:
                continue
    return None


def load_or_match_caddings(force_rematch: bool = False) -> Dict[str, Any]:
    """Load matched CADD scores from cache, or run tabix matching if not cached.

    Uses the augmented file (tabix + REST API) when present, falling back to
    the tabix-only file. The augmented file contains CADD v1.6 scores for the
    3 possible alt alleles at each SNV position (matched by ref+alt).
    """
    augmented = MATCHES_PATH.parent / "cadd_matches_augmented.json"
    if augmented.exists() and not force_rematch:
        return json.load(open(augmented))
    if MATCHES_PATH.exists() and not force_rematch:
        return json.load(open(MATCHES_PATH))

    data = json.load(open(COHORT_MUTATIONS))
    cohort = data["cohort_20_patients"]
    all_muts = [m for v in cohort.values() for m in v]
    snv = [m for m in all_muts if len(m["ref"]) == 1 and len(m["alt"]) == 1]
    indel = [m for m in all_muts if not (len(m["ref"]) == 1 and len(m["alt"]) == 1)]

    print(f"Matching {len(snv)} SNVs against {SNV_TSV.name}...")
    s = time.time()
    matched_snv, unmatched_snv = [], []
    for i, m in enumerate(snv):
        hit = _cadd_lookup_one(m["chrom"], m["pos"], m["ref"], m["alt"], SNV_TSV)
        if hit is not None:
            phred, raw = hit
            matched_snv.append({**m, "cadd_phred": phred, "cadd_raw": raw})
        else:
            unmatched_snv.append(m)
        if (i + 1) % 1000 == 0:
            print(f"  ...{i+1}/{len(snv)} elapsed={time.time()-s:.1f}s")
    print(f"  SNV matched: {len(matched_snv)}/{len(snv)} in {time.time()-s:.1f}s")

    print(f"Matching {len(indel)} indels against {INDEL_TSV.name}...")
    s = time.time()
    matched_indel, unmatched_indel = [], []
    for m in indel:
        hit = _cadd_lookup_one(m["chrom"], m["pos"], m["ref"], m["alt"], INDEL_TSV)
        if hit is not None:
            phred, raw = hit
            matched_indel.append({**m, "cadd_phred": phred, "cadd_raw": raw})
        else:
            unmatched_indel.append(m)
    print(f"  Indel matched: {len(matched_indel)}/{len(indel)} in {time.time()-s:.1f}s")

    phreds = [m["cadd_phred"] for m in matched_snv]
    out = {
        "snv_matched": matched_snv,
        "snv_unmatched": unmatched_snv,
        "indel_matched": matched_indel,
        "indel_unmatched": unmatched_indel,
        "n_total_snv": len(snv),
        "n_matched_snv": len(matched_snv),
        "n_total_indel": len(indel),
        "n_matched_indel": len(matched_indel),
        "snv_phred_stats": {
            "min": min(phreds) if phreds else None,
            "median": statistics.median(phreds) if phreds else None,
            "mean": statistics.mean(phreds) if phreds else None,
            "max": max(phreds) if phreds else None,
        },
        "data_sources": {
            "snv_file": str(SNV_TSV),
            "snv_size_bytes": SNV_TSV.stat().st_size if SNV_TSV.exists() else None,
            "indel_file": str(INDEL_TSV),
            "indel_size_bytes": INDEL_TSV.stat().st_size if INDEL_TSV.exists() else None,
            "snv_source_url": "https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.snv.tsv.gz",
            "snv_description": "CADD v1.6 GRCh38 SNV scores for gnomAD genomes r3.0 observed variants (advertised 5.9 GB)",
            "indel_source_url": "https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz",
            "indel_description": "CADD v1.7 GRCh38 InDel scores for gnomAD genomes r4.0 observed variants (1.2 GB)",
        },
        "cadd_version": "v1.6 (SNV) / v1.7 (indel)",
        "reference": "Kircher et al. 2014, Nature Genetics, doi:10.1038/ng.2892",
        "license": "CC BY-NC-SA 4.0 (non-commercial)",
    }
    with open(MATCHES_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved matches to {MATCHES_PATH}")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Build per-patient weight vectors
# ─────────────────────────────────────────────────────────────────────────────
def build_per_patient_weights(cohort: Dict[str, List[Dict]],
                              matches: Dict[str, Any],
                              imputation: str = "median"
                              ) -> Dict[str, np.ndarray]:
    """Return {patient_id: weights_array} aligned with the cohort's mutation list.

    imputation:
      - 'zero'     : unmatched weights = 0 (effectively drops those loci)
      - 'median'   : unmatched weights = median PHRED of matched
      - 'one'      : unmatched weights = 1.0 (same as uniform, sanity check)
    """
    matched_by_key = {}
    for m in matches["snv_matched"] + matches["indel_matched"]:
        key = (m["sample"], m["chrom"], m["pos"], m["ref"], m["alt"])
        matched_by_key[key] = m["cadd_phred"]

    phreds = list(matched_by_key.values())
    if imputation == "median":
        default = statistics.median(phreds) if phreds else 1.0
    elif imputation == "zero":
        default = 0.0
    elif imputation == "one":
        default = 1.0
    else:
        raise ValueError(f"unknown imputation: {imputation}")

    weights_per_patient = {}
    for patient, muts in cohort.items():
        w = np.zeros(len(muts))
        for i, m in enumerate(muts):
            key = (patient, m["chrom"], m["pos"], m["ref"], m["alt"])
            w[i] = matched_by_key.get(key, default)
        weights_per_patient[patient] = w
    return weights_per_patient


def build_topk_per_patient_weights(cohort: Dict[str, List[Dict]],
                                   matches: Dict[str, Any],
                                   top_k: int,
                                   min_phred: float = 0.0) -> Dict[str, np.ndarray]:
    """Weight vector that is 1 for the top-K CADD mutations (by PHRED) per
    patient, 0 otherwise. Optionally restrict to PHRED >= min_phred.
    """
    matched_by_key = {}
    for m in matches["snv_matched"] + matches["indel_matched"]:
        key = (m["sample"], m["chrom"], m["pos"], m["ref"], m["alt"])
        matched_by_key[key] = m["cadd_phred"]

    weights_per_patient = {}
    for patient, muts in cohort.items():
        # build (idx, phred) for matched
        scored = []
        for i, m in enumerate(muts):
            key = (patient, m["chrom"], m["pos"], m["ref"], m["alt"])
            phred = matched_by_key.get(key)
            if phred is not None and phred >= min_phred:
                scored.append((i, phred))
        # sort by PHRED desc, take top K
        scored.sort(key=lambda x: -x[1])
        kept_indices = {i for i, _ in scored[:top_k]}
        w = np.array([1.0 if i in kept_indices else 0.0 for i in range(len(muts))])
        weights_per_patient[patient] = w
    return weights_per_patient


# ─────────────────────────────────────────────────────────────────────────────
# Weighted panel detection
# ─────────────────────────────────────────────────────────────────────────────
def run_weighted_panel_detection(
    cohort: Dict[str, List[Dict]],
    weights_per_patient: Dict[str, np.ndarray],
    tumor_fractions: Optional[List[float]] = None,
    seeds: Optional[List[int]] = None,
    cfdna_depth: int = 5000,
    bg_error_rate: float = 0.002,
) -> Dict[str, List[Dict]]:
    """Per-sample detection using CADD-weighted LLR aggregation.

    Mirrors real_tcga_validation.run_panel_detection but uses per-patient
    weight vectors: score = Σ w_i * LLR_i.
    """
    if tumor_fractions is None:
        tumor_fractions = [0.1, 0.05, 0.01, 0.005, 0.001]
    if seeds is None:
        seeds = [42, 123, 456, 789, 1024]

    patients = list(cohort.keys())
    results = {'panel_llr_weighted': [], 'panel_llr_uniform': [],
               'panel_fisher': [], 'panel_strand': []}

    for tf in tumor_fractions:
        print(f"\n  Panel detection @ TF={tf*100:.2f}% ({len(patients)} patients × {len(seeds)} seeds)")
        llr_w_by_seed = {}
        llr_u_by_seed = {}
        fisher_by_seed = {}
        strand_by_seed = {}
        for seed in seeds:
            pos_w, neg_w = [], []
            pos_u, neg_u = [], []
            pos_fish, neg_fish = [], []
            pos_strand, neg_strand = [], []
            for patient in patients:
                muts = cohort[patient]
                weights = weights_per_patient[patient]

                dp = simulate_cfdna_from_real(muts, tumor_fraction=tf,
                                              cfdna_depth=cfdna_depth, seed=seed,
                                              bg_error_rate=bg_error_rate)
                dn = simulate_cfdna_from_real(muts, tumor_fraction=0.0,
                                              cfdna_depth=cfdna_depth, seed=seed,
                                              bg_error_rate=bg_error_rate)
                lp = compute_llr_scores(dp['depths'], dp['X'][:, 1].astype(int), dp['X'][:, 3])
                ln = compute_llr_scores(dn['depths'], dn['X'][:, 1].astype(int), dn['X'][:, 3])

                nv_p, nv_n = dp['n_variants'], dn['n_variants']
                panel_size = min(nv_p, nv_n)
                # In simulate_cfdna_from_real, mutations are placed in the first
                # n_variants slots of dp['X'], so weights for muts[i] line up
                # with dp['X'][i].
                wp = weights[:panel_size]
                wn = weights[:nv_n]

                pos_u.append(float(lp[:panel_size].sum()))
                neg_u.append(float(ln[:panel_size].sum()))
                pos_w.append(float((lp[:panel_size] * wp).sum()))
                neg_w.append(float((ln[:panel_size] * wn).sum()))

                fp = _fisher_scores(dp['depths'][:panel_size],
                                    dp['X'][:, 1].astype(int)[:panel_size],
                                    dp['X'][:, 3][:panel_size])
                fn = _fisher_scores(dn['depths'][:nv_n],
                                    dn['X'][:, 1].astype(int)[:nv_n],
                                    dn['X'][:, 3][:nv_n])
                pos_fish.append(float(fp.sum()))
                neg_fish.append(float(fn.sum()))
                sc_p = dp['strand_conc'][:panel_size]
                sc_n = dn['strand_conc'][:nv_n]
                pos_strand.append(float((fp * sc_p).sum()))
                neg_strand.append(float((fn * sc_n).sum()))

            y = np.array([1] * len(pos_w) + [0] * len(neg_w))
            llr_w_by_seed[seed] = _panel_metrics(y, np.array(pos_w + neg_w))
            llr_u_by_seed[seed] = _panel_metrics(y, np.array(pos_u + neg_u))
            fisher_by_seed[seed] = _panel_metrics(y, np.array(pos_fish + neg_fish))
            strand_by_seed[seed] = _panel_metrics(y, np.array(pos_strand + neg_strand))
            print(f"    seed {seed}: "
                  f"weighted-AUC={llr_w_by_seed[seed]['auc']:.4f}  "
                  f"uniform-AUC={llr_u_by_seed[seed]['auc']:.4f}  "
                  f"Fisher-AUC={fisher_by_seed[seed]['auc']:.4f}")

        for key, by_seed in (
            ('panel_llr_weighted', llr_w_by_seed),
            ('panel_llr_uniform', llr_u_by_seed),
            ('panel_fisher', fisher_by_seed),
            ('panel_strand', strand_by_seed),
        ):
            for m in ('auc', 'sens_at_95_spec', 'sens_at_99_spec', 'paired_win_rate'):
                if m not in by_seed[seeds[0]]:
                    continue
                vals = [by_seed[s][m] for s in seeds]
                results[key].append({
                    'tumor_fraction': tf,
                    'metric': m,
                    'mean': float(np.mean(vals)),
                    'std': float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                    'per_seed': {str(s): by_seed[s][m] for s in seeds},
                })
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Main entry: run all comparisons and save results
# ─────────────────────────────────────────────────────────────────────────────
def main():
    print("=" * 80)
    print("CADD-weighted panel LLR for DeepCatch v2.2")
    print("=" * 80)

    print("\n[1] Loading TCGA-LUAD cohort and CADD matches...")
    cohort_data = json.load(open(COHORT_MUTATIONS))
    cohort = cohort_data["cohort_20_patients"]
    n_muts = sum(len(v) for v in cohort.values())
    n_patients = len(cohort)
    print(f"  Cohort: {n_patients} patients, {n_muts} mutations")

    matches = load_or_match_caddings()
    n_matched = len(matches["snv_matched"]) + len(matches["indel_matched"])
    print(f"  CADD matches: SNV={len(matches['snv_matched'])}/{matches['n_total_snv']} "
          f"({len(matches['snv_matched'])/matches['n_total_snv']*100:.1f}%) "
          f"indel={len(matches['indel_matched'])}/{matches['n_total_indel']} "
          f"({len(matches['indel_matched'])/matches['n_total_indel']*100:.1f}%) "
          f"overall={n_matched}/{n_muts} ({n_matched/n_muts*100:.1f}%)")

    tumor_fractions = [0.1, 0.05, 0.01, 0.005, 0.001]
    seeds = [42, 123, 456, 789, 1024]

    print("\n[2] CADD-weighted LLR (median imputation for unmatched)...")
    weights_median = build_per_patient_weights(cohort, matches, imputation="median")
    res_median = run_weighted_panel_detection(cohort, weights_median,
                                              tumor_fractions, seeds)

    print("\n[3] CADD-weighted LLR (zero weight for unmatched)...")
    weights_zero = build_per_patient_weights(cohort, matches, imputation="zero")
    res_zero = run_weighted_panel_detection(cohort, weights_zero,
                                            tumor_fractions, seeds)

    print("\n[4] Top-K=500 by CADD PHRED (1 if top-500 else 0)...")
    weights_top500 = build_topk_per_patient_weights(cohort, matches, top_k=500)
    res_top500 = run_weighted_panel_detection(cohort, weights_top500,
                                              tumor_fractions, seeds)

    print("\n[5] Top-K=1000 by CADD PHRED (1 if top-1000 else 0)...")
    weights_top1000 = build_topk_per_patient_weights(cohort, matches, top_k=1000)
    res_top1000 = run_weighted_panel_detection(cohort, weights_top1000,
                                               tumor_fractions, seeds)

    print("\n[6] Top-K=2000 by CADD PHRED (1 if top-2000 else 0)...")
    weights_top2000 = build_topk_per_patient_weights(cohort, matches, top_k=2000)
    res_top2000 = run_weighted_panel_detection(cohort, weights_top2000,
                                               tumor_fractions, seeds)

    print("\n[7] Top-K=2000 AND PHRED >= 20 (top 1% deleterious)...")
    weights_top2000_del = build_topk_per_patient_weights(cohort, matches,
                                                         top_k=2000, min_phred=20.0)
    res_top2000_del = run_weighted_panel_detection(cohort, weights_top2000_del,
                                                   tumor_fractions, seeds)

    print("\n[8] Top-K=1000 AND PHRED >= 20 (top 1% deleterious)...")
    weights_top1000_del = build_topk_per_patient_weights(cohort, matches,
                                                         top_k=1000, min_phred=20.0)
    res_top1000_del = run_weighted_panel_detection(cohort, weights_top1000_del,
                                                   tumor_fractions, seeds)

    # Compile results
    out = {
        "experiment": "CADD-weighted panel LLR for DeepCatch v2.2",
        "cohort": {"n_patients": n_patients, "n_mutations": n_muts},
        "cadd_matches": {
            "n_total_snv": matches["n_total_snv"],
            "n_matched_snv": len(matches["snv_matched"]),
            "n_total_indel": matches["n_total_indel"],
            "n_matched_indel": len(matches["indel_matched"]),
            "match_rate_snv": len(matches["snv_matched"]) / matches["n_total_snv"],
            "match_rate_overall": (len(matches["snv_matched"]) +
                                   len(matches["indel_matched"])) / n_muts,
            "phred_stats_matched": matches["snv_phred_stats"],
            "data_sources": matches["data_sources"],
            "reference": matches["reference"],
            "license": matches["license"],
            "cadd_version": matches["cadd_version"],
        },
        "tumor_fractions": tumor_fractions,
        "seeds": seeds,
        "results": {
            "cadd_weighted_median_imputation": res_median,
            "cadd_weighted_zero_imputation": res_zero,
            "topk_500": res_top500,
            "topk_1000": res_top1000,
            "topk_2000": res_top2000,
            "topk_2000_phred_ge_20": res_top2000_del,
            "topk_1000_phred_ge_20": res_top1000_del,
        },
    }
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()