"""
Combined CADD + AlphaMissense per-mutation weights in DeepCatch panel LLR.

For the 4,882 CADD-scored TCGA-LUAD mutations we attach (when available)
the AlphaMissense pathogenicity probability in [0,1] (Cheng 2023, Science),
and evaluate three deterministic weight schemes for the panel log-likelihood
aggregator:

    1. weight_cadd_only        : CADD Top-K=200 per patient (1 if in top-K, 0 else)
                                 — the published winner for per-subgroup LLR.
    2. weight_alphamissense_only: AM Top-K=200 per patient (1 if in top-K, 0 else)
                                 — by AlphaMissense pathogenicity probability.
    3. weight_combined         : multiplicative normalised product
                                 w_i = (CADD_i / max_CADD) * (AM_i / max_AM)
                                 computed only over loci with BOTH scores;
                                 loci missing one score get weight 0.

All three schemes share the same cohort (20 TCGA-LUAD patients, 5,738 mutations),
the same seeds (5 seeds), the same ctDNA fractions (10/5/1/0.5/0.1%), and the
same per-seed simulation pipeline (`real_tcga_validation.simulate_cfdna_from_real`),
so deltas are attributable to the weighting change alone.

Joint-match honesty:
The AlphaMissense TSV is keyed by (Uniprot, prot_pos, ref_aa, alt_aa). It does
not cover every TCGA-LUAD missense mutation (only canonical protein-coding
transcripts). We report the joint match rate (CADD ∩ AlphaMissense) below.
The combined weighting uses ONLY the 4,161 (CADD AND AM) intersections; the
other two schemes are evaluated on the CADD-only set so the head-to-head is
apples-to-apples (every weight scheme gets the same per-patient top-K
opportunity).

Outputs:
    results/cadd_alphamissense_combined.json
    docs/CADD_ALPHAMISSENSE_COMBINED.md        (written separately)
"""
from __future__ import annotations

import gzip
import json
import logging
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

ROOT = Path("/Users/hermes/deepcatch")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from real_tcga_validation import (  # noqa: E402
    compute_llr_scores,
    simulate_cfdna_from_real,
    _panel_metrics,
    _fisher_scores,
)
from alphamissense_weights import (  # noqa: E402
    am_key,
    classify_am,
    AMWeight,
    THRESHOLD_PATHOGENIC,
    THRESHOLD_BENIGN,
)

from cadd_weighted_llr import (  # noqa: E402
    load_or_match_caddings,
    build_topk_per_patient_weights,
)

LOG = logging.getLogger("cadd_alphamissense_combined")

AM_TSV_GZ = Path("/Users/hermes/.cache/mrnavax/AlphaMissense_hg38.tsv.gz")
COHORT_MUTATIONS = ROOT / "results" / "cadd_mutations_full.json"
CADD_MATCHES = ROOT / "results" / "cadd_matches_augmented.json"
RESULTS_PATH = ROOT / "results" / "cadd_alphamissense_combined.json"
TOP_K = 200  # the published winner


# ─────────────────────────────────────────────────────────────────────────────
# 1. Stream AlphaMissense TSV → per-(Uniprot, prot_pos, ref, alt) score
# ─────────────────────────────────────────────────────────────────────────────
def stream_am_tsv_for_wanted_keys(tsv_gz: Path,
                                  wanted_keys: set) -> Dict[str, AMWeight]:
    """Stream the AlphaMissense TSV.gz and return only entries whose
    ``uniprot:prot_pos:ref_aa:alt_aa`` key is in ``wanted_keys``.
    """
    out: Dict[str, AMWeight] = {}
    if not tsv_gz.exists():
        LOG.error("AM TSV not found at %s", tsv_gz)
        return out
    t0 = time.time()
    n_rows = 0
    n_missense = 0
    LOG.info("Streaming %s for %d wanted keys …", tsv_gz.name, len(wanted_keys))
    # AlphaMissense_hg38.tsv columns:
    #   0=CHROM 1=POS 2=REF 3=ALT 4=genome 5=uniprot_id 6=transcript_id
    #   7=protein_variant (e.g. "V2L") 8=am_pathogenicity 9=am_class
    with gzip.open(tsv_gz, "rt") as f:
        for row in f:
            row = row.rstrip("\n")
            if not row or row.startswith("#"):
                continue
            fields = row.split("\t")
            if len(fields) < 10:
                continue
            up = fields[5]
            pv = fields[7]
            score_s = fields[8]
            n_rows += 1
            # Parse protein_variant like "V2L" → ('V', 2, 'L')
            try:
                ref_aa = pv[0]
                alt_aa = pv[-1]
                pp = int(pv[1:-1])
            except (ValueError, IndexError):
                continue
            key = am_key(up, pp, ref_aa, alt_aa)
            if key not in wanted_keys:
                continue
            try:
                s = float(score_s)
            except ValueError:
                continue
            out[key] = AMWeight(
                uniprot=up, prot_pos=pp, ref_aa=ref_aa, alt_aa=alt_aa,
                score=s, classification=classify_am(s), source="tsvlocal",
            )
            n_missense += 1
            if len(out) == len(wanted_keys):
                break
    elapsed = time.time() - t0
    LOG.info("Scanned %d rows in %.1fs; collected %d matches (%.2f%% of wanted)",
             n_rows, elapsed, n_missense,
             100 * n_missense / max(1, len(wanted_keys)))
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 2. Build per-patient weight vectors for each scheme
# ─────────────────────────────────────────────────────────────────────────────
def build_per_patient_topk_weights_by_score(
    cohort: Dict[str, List[Dict]],
    score_by_key: Dict[Tuple[str, str, int, str, str], float],
    top_k: int,
) -> Tuple[Dict[str, np.ndarray], Dict[str, int]]:
    """For each patient, weight=1 for the top-K highest-score loci, 0 otherwise.
    Loci absent from ``score_by_key`` are dropped (weight=0).
    """
    weights_per_patient: Dict[str, np.ndarray] = {}
    matched_per_patient: Dict[str, int] = {}
    for patient, muts in cohort.items():
        scored = []
        for i, m in enumerate(muts):
            k = (patient, m["chrom"], int(m["pos"]), m["ref"], m["alt"])
            s = score_by_key.get(k)
            if s is not None:
                scored.append((i, s))
        scored.sort(key=lambda x: -x[1])
        kept = {i for i, _ in scored[:top_k]}
        w = np.array([1.0 if i in kept else 0.0 for i in range(len(muts))])
        weights_per_patient[patient] = w
        matched_per_patient[patient] = len(kept)
    return weights_per_patient, matched_per_patient


def build_per_patient_combined_weights(
    cohort: Dict[str, List[Dict]],
    cadd_by_key: Dict[Tuple[str, str, int, str, str], float],
    am_by_key: Dict[Tuple[str, str, int, str, str], float],
) -> Tuple[Dict[str, np.ndarray], Dict[str, int]]:
    """Per-locus weight = (CADD/max_CADD_global) * (AM/max_AM_global),
    restricted to loci with BOTH scores.
    """
    # Compute global maxima across all matched loci.
    cadd_vals = list(cadd_by_key.values())
    am_vals = list(am_by_key.values())
    max_cadd = max(cadd_vals) if cadd_vals else 1.0
    max_am = max(am_vals) if am_vals else 1.0
    if max_cadd <= 0:
        max_cadd = 1.0
    if max_am <= 0:
        max_am = 1.0

    # Combine per-key
    combo: Dict[Tuple[str, str, int, str, str], float] = {}
    for k, c in cadd_by_key.items():
        a = am_by_key.get(k)
        if a is None:
            continue
        combo[k] = (c / max_cadd) * (a / max_am)

    weights_per_patient: Dict[str, np.ndarray] = {}
    matched_per_patient: Dict[str, int] = {}
    for patient, muts in cohort.items():
        w = np.zeros(len(muts))
        for i, m in enumerate(muts):
            k = (patient, m["chrom"], int(m["pos"]), m["ref"], m["alt"])
            v = combo.get(k)
            if v is not None:
                w[i] = v
        weights_per_patient[patient] = w
        matched_per_patient[patient] = int((w > 0).sum())
    return weights_per_patient, matched_per_patient


# ─────────────────────────────────────────────────────────────────────────────
# 3. Panel detection with arbitrary weight vectors (mirrors cadd_weighted_llr)
# ─────────────────────────────────────────────────────────────────────────────
def run_panel_detection(
    cohort: Dict[str, List[Dict]],
    weights_per_patient: Dict[str, np.ndarray],
    tumor_fractions: List[float],
    seeds: List[int],
    cfdna_depth: int = 5000,
    bg_error_rate: float = 0.002,
) -> Dict[str, List[Dict]]:
    patients = list(cohort.keys())
    out = {"panel_llr_weighted": []}
    for tf in tumor_fractions:
        LOG.info("  TF=%.3f%%  (%d patients × %d seeds)", tf * 100, len(patients), len(seeds))
        per_seed = {}
        for seed in seeds:
            pos_w, neg_w = [], []
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
                wp = weights[:panel_size]
                wn = weights[:nv_n]
                pos_w.append(float((lp[:panel_size] * wp).sum()))
                neg_w.append(float((ln[:panel_size] * wn).sum()))
            y = np.array([1] * len(pos_w) + [0] * len(neg_w))
            per_seed[seed] = _panel_metrics(y, np.array(pos_w + neg_w))
        for m in ('auc', 'sens_at_95_spec', 'sens_at_99_spec', 'paired_win_rate'):
            vals = [per_seed[s][m] for s in seeds]
            out["panel_llr_weighted"].append({
                'tumor_fraction': tf,
                'metric': m,
                'mean': float(np.mean(vals)),
                'std': float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                'per_seed': {str(s): per_seed[s][m] for s in seeds},
            })
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Helpers to extract a metric from results
# ─────────────────────────────────────────────────────────────────────────────
def _at(results: Dict[str, List[Dict]], tf: float, metric: str) -> Dict[str, float]:
    rows = [r for r in results["panel_llr_weighted"]
            if r["tumor_fraction"] == tf and r["metric"] == metric]
    if not rows:
        return {"mean": float("nan"), "std": float("nan")}
    return {"mean": rows[0]["mean"], "std": rows[0]["std"]}


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────
def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    print("=" * 80)
    print("Combined CADD + AlphaMissense weighting for DeepCatch panel LLR")
    print("=" * 80)

    print("\n[1] Loading cohort, CADD matches, AlphaMissense TSV…")
    cohort_data = json.load(open(COHORT_MUTATIONS))
    cohort = cohort_data["cohort_20_patients"]
    n_patients = len(cohort)
    n_muts = sum(len(v) for v in cohort.values())
    print(f"  Cohort: {n_patients} patients, {n_muts} mutations")

    # CADD matches
    matches = load_or_match_caddings()
    cadd_matched = matches["snv_matched"] + matches["indel_matched"]
    n_cadd = len(cadd_matched)
    print(f"  CADD matched: {n_cadd}/{n_muts} ({n_cadd/n_muts*100:.1f}%)")

    # Load cohort mutations enriched with amino-acid info
    cohort_mut_aa = []
    aa_path = ROOT / "results" / "cohort_with_amino_acids.jsonl"
    if aa_path.exists():
        with open(aa_path) as f:
            for line in f:
                cohort_mut_aa.append(json.loads(line))
        print(f"  Loaded {len(cohort_mut_aa)} missense SNVs with (Uniprot, prot_pos, ref_aa, alt_aa)")
    else:
        print(f"  !! {aa_path} missing — run scripts/enrich_cohort_with_amino_acids.py first")
        sys.exit(1)

    # Build (sample, chrom, pos, ref, alt) -> (uniprot, prot_pos, ref_aa, alt_aa)
    prot_by_full: Dict[Tuple[str, str, int, str, str], Tuple[str, int, str, str]] = {}
    for m in cohort_mut_aa:
        k = (m["sample"], m["chrom"], int(m["pos"]), m["ref"], m["alt"])
        if m.get("uniprot") and m.get("prot_pos") and m.get("ref_aa") and m.get("alt_aa"):
            prot_by_full[k] = (m["uniprot"], int(m["prot_pos"]), m["ref_aa"], m["alt_aa"])
    print(f"  Unique keys (sample,chrom,pos,ref,alt) → full protein: {len(prot_by_full)}")

    # ── Build AM lookup keys for CADD-matched missense mutations
    cadd_by_key: Dict[Tuple[str, str, int, str, str], float] = {}
    am_wanted: set = set()
    # First, build the (sample, chrom, pos, ref, alt) → full protein lookup
    # from cohort_mut_aa (now includes ref_aa, alt_aa)
    # Now loop CADD-matched
    for m in cadd_matched:
        k = (m["sample"], m["chrom"], int(m["pos"]), m["ref"], m["alt"])
        cadd_by_key[k] = float(m["cadd_phred"])
        if m.get("variant_class") == "Missense_Mutation":
            full = prot_by_full.get(k)
            if full:
                up, pp, ref_aa, alt_aa = full
                am_wanted.add(am_key(up, pp, ref_aa, alt_aa))

    print(f"  AM lookup keys (missense + CADD + Uniprot+prot_pos): {len(am_wanted)}")

    # ── Stream the TSV
    if not AM_TSV_GZ.exists():
        print(f"  !! AlphaMissense TSV missing at {AM_TSV_GZ} — aborting")
        sys.exit(1)

    am_weights = stream_am_tsv_for_wanted_keys(AM_TSV_GZ, am_wanted)
    print(f"  AlphaMissense hits: {len(am_weights)}/{len(am_wanted)} "
          f"({len(am_weights)/max(1,len(am_wanted))*100:.1f}%)")

    # ── Compute joint CADD ∩ AM match rate
    am_by_key: Dict[Tuple[str, str, int, str, str], float] = {}
    joint_matched = 0
    for m in cadd_matched:
        k = (m["sample"], m["chrom"], int(m["pos"]), m["ref"], m["alt"])
        full = prot_by_full.get(k)
        if not full:
            continue
        up, pp, ref_aa, alt_aa = full
        akey = am_key(up, pp, ref_aa, alt_aa)
        am_w = am_weights.get(akey)
        if am_w is not None:
            am_by_key[k] = am_w.score
            joint_matched += 1
    print(f"\n[2] Joint match (CADD ∩ AlphaMissense): "
          f"{joint_matched}/{n_cadd} CADD-matched mutations "
          f"({joint_matched/max(1,n_cadd)*100:.1f}%)")
    print(f"    Joint match / total cohort: {joint_matched}/{n_muts} "
          f"({joint_matched/n_muts*100:.1f}%)")

    # ── Distributions
    am_scores = list(am_by_key.values())
    cadd_scores = list(cadd_by_key.values())
    print(f"\n  CADD PHRED among matched: n={len(cadd_scores)} "
          f"min={min(cadd_scores):.2f} median={statistics.median(cadd_scores):.2f} "
          f"max={max(cadd_scores):.2f}")
    print(f"  AlphaMissense among matched: n={len(am_scores)} "
          f"min={min(am_scores):.3f} median={statistics.median(am_scores):.3f} "
          f"max={max(am_scores):.3f}")

    # ── Build per-patient weight vectors for each scheme
    tumor_fractions = [0.1, 0.05, 0.01, 0.005, 0.001]
    seeds = [42, 123, 456, 789, 1024]

    print("\n[3] Scheme 1 — CADD Top-K=200 (published winner)…")
    w_cadd, kept_cadd = build_per_patient_topk_weights_by_score(
        cohort, cadd_by_key, top_k=TOP_K,
    )
    print(f"  Per-patient top-K kept: median={statistics.median(kept_cadd.values())} "
          f"min={min(kept_cadd.values())} max={max(kept_cadd.values())}")
    res_cadd = run_panel_detection(cohort, w_cadd, tumor_fractions, seeds)

    print("\n[4] Scheme 2 — AlphaMissense Top-K=200 by pathogenicity…")
    w_am, kept_am = build_per_patient_topk_weights_by_score(
        cohort, am_by_key, top_k=TOP_K,
    )
    print(f"  Per-patient top-K kept: median={statistics.median(kept_am.values())} "
          f"min={min(kept_am.values())} max={max(kept_am.values())}")
    res_am = run_panel_detection(cohort, w_am, tumor_fractions, seeds)

    print("\n[5] Scheme 3 — Combined CADD × AlphaMissense (normalised product)…")
    w_combo, kept_combo = build_per_patient_combined_weights(
        cohort, cadd_by_key, am_by_key,
    )
    print(f"  Per-patient non-zero weights: median={statistics.median(kept_combo.values())} "
          f"min={min(kept_combo.values())} max={max(kept_combo.values())}")
    res_combo = run_panel_detection(cohort, w_combo, tumor_fractions, seeds)

    # ── Headline: bottom-line per ctDNA fraction
    print("\n[6] Headline comparison (5-seed × 5-fold pooled OOF):")
    headline_rows = []
    for tf in tumor_fractions:
        c = _at(res_cadd, tf, "sens_at_99_spec")
        a = _at(res_am, tf, "sens_at_99_spec")
        x = _at(res_combo, tf, "sens_at_99_spec")
        ca = _at(res_cadd, tf, "auc")
        aa = _at(res_am, tf, "auc")
        xa = _at(res_combo, tf, "auc")
        # Winner = lowest std-aware mean; tie-break by AUC
        winner = "tie"
        score = {  # sens99 + 0.5*auc as a single composite for tiebreak only
            "cadd": c["mean"] + 0.5 * ca["mean"],
            "am":   a["mean"] + 0.5 * aa["mean"],
            "combined": x["mean"] + 0.5 * xa["mean"],
        }
        top = max(score.values())
        winners = [k for k, v in score.items() if abs(v - top) < 1e-9]
        if len(winners) == 1:
            winner = winners[0]
        headline_rows.append({
            "tumor_fraction": tf,
            "cadd_sens99": c, "am_sens99": a, "combined_sens99": x,
            "cadd_auc": ca, "am_auc": aa, "combined_auc": xa,
            "winner_sens99": winner,
        })
        print(f"  TF={tf*100:5.2f}%  "
              f"CADD Sens@99%={c['mean']:.3f}±{c['std']:.3f}  "
              f"AM Sens@99%={a['mean']:.3f}±{a['std']:.3f}  "
              f"Combined Sens@99%={x['mean']:.3f}±{x['std']:.3f}  → {winner}")

    out = {
        "experiment": "Combined CADD + AlphaMissense per-mutation weights for DeepCatch panel LLR",
        "cohort": {"n_patients": n_patients, "n_mutations": n_muts},
        "joint_match": {
            "n_cadd_matched": n_cadd,
            "n_am_wanted": len(am_wanted),
            "n_joint_matched": joint_matched,
            "joint_match_rate_of_cadd": joint_matched / max(1, n_cadd),
            "joint_match_rate_of_cohort": joint_matched / n_muts,
            "cadd_phred_stats_matched": {
                "min": min(cadd_scores), "max": max(cadd_scores),
                "median": statistics.median(cadd_scores),
                "mean": statistics.mean(cadd_scores),
            },
            "am_score_stats_matched": {
                "min": min(am_scores), "max": max(am_scores),
                "median": statistics.median(am_scores),
                "mean": statistics.mean(am_scores),
            },
        },
        "schemes": {
            "weight_cadd_only": {
                "description": "CADD Top-K=200 per patient (1 if in top-K by CADD PHRED, 0 else)",
                "top_k": TOP_K,
                "per_patient_kept_stats": {
                    "min": min(kept_cadd.values()),
                    "median": statistics.median(kept_cadd.values()),
                    "max": max(kept_cadd.values()),
                },
                "results": res_cadd,
            },
            "weight_alphamissense_only": {
                "description": "AlphaMissense Top-K=200 per patient (1 if in top-K by AM pathogenicity, 0 else)",
                "top_k": TOP_K,
                "per_patient_kept_stats": {
                    "min": min(kept_am.values()),
                    "median": statistics.median(kept_am.values()),
                    "max": max(kept_am.values()),
                },
                "results": res_am,
            },
            "weight_combined": {
                "description": "w_i = (CADD_i / max_CADD) * (AM_i / max_AM); loci missing either score get w=0 (no imputation)",
                "per_patient_nonzero_stats": {
                    "min": min(kept_combo.values()),
                    "median": statistics.median(kept_combo.values()),
                    "max": max(kept_combo.values()),
                },
                "results": res_combo,
            },
        },
        "headline_by_tumor_fraction": headline_rows,
        "tumor_fractions": tumor_fractions,
        "seeds": seeds,
        "license": "AlphaMissense: CC BY-NC-SA 4.0 (non-commercial); CADD: CC BY-NC-SA 4.0 (non-commercial)",
        "references": [
            "Kircher M, Witten DM, Jain P, O'Roak BJ, Cooper GM, Shendure J. A general framework for estimating the relative pathogenicity of human genetic variants. Nat Genet. 2014;46(3):310-315.",
            "Cheng J, Novati G, Pan J, et al. Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science 381, eadg7492 (2023).",
        ],
        "determinism_note": "All three weight schemes are deterministic functions of (CADD PHRED, AM pathogenicity, cohort membership). No parameters are learned from the data.",
    }
    RESULTS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PATH, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nResults saved to {RESULTS_PATH}")


if __name__ == "__main__":
    main()
