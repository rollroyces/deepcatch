# Combined CADD + AlphaMissense per-mutation weights — DeepCatch panel LLR

**Date:** 2026-09-17
**Branch:** main @ fd1c3e3
**Working commit:** companion of commit 0ffd3fc (CADD Top-K=200 per-subgroup).

This experiment evaluates whether weighting the panel log-likelihood
aggregator with the **multiplicative product** of CADD PHRED and
AlphaMissense pathogenicity beats either score alone. All three weight
schemes are deterministic functions of (CADD PHRED, AM pathogenicity,
cohort membership) — **no parameter is learned from the data**.

---

## Bottom line

| ctDNA | CADD-only | AM-only | **Combined** | Winner |
|------:|----------:|--------:|-------------:|:-------|
| 10 %  | **1.000** | 1.000   | 1.000        | tie    |
| 5 %   | **1.000** | 1.000   | 1.000        | tie    |
| 1 %   | **1.000** | 1.000   | 1.000        | tie    |
| 0.5 % | **1.000** | 1.000   | 0.990        | tie (CADD best on AUC) |
| 0.1 % | **0.910 ± 0.042** | 0.850 ± 0.100 | 0.510 ± 0.065 | **CADD wins** |

**AUC at 0.1 % ctDNA**: CADD 0.979, AM 0.959, **Combined 0.891**.

The combined multiplicative weight **does not beat CADD alone** at any ctDNA
fraction. At the hard 0.1 % floor it is dramatically worse (0.51 vs 0.91)
because the product suppresses high-CADD mutations that don't happen to be
in the AlphaMissense index (see *Honest match rate* below).

**Published-winner anchor:** CADD Top-K=200 (commit 0ffd3fc) — the head-to-head
loser here is the combined scheme, **not** CADD.

---

## Joint match rate (CADD ∩ AlphaMissense)

| Step | Count | % of CADD-matched | % of full cohort |
|------|------:|------------------:|-----------------:|
| Total cohort (5,738 mutations, 20 patients) | 5,738 | — | 100.0 % |
| CADD-matched (germline/SNP-overlap) | 4,882 | 100.0 % | 85.1 % |
| CADD-matched AND missense AND has Uniprot+protein position | 4,066 | 83.3 % | 70.9 % |
| **Joint CADD ∩ AlphaMissense hit** | **3,946** | **80.8 %** | **68.8 %** |

The 120-mutation gap (4,066 wanted − 3,946 hits) reflects rare/private
variants and isoforms not present in the AlphaMissense 71.7 M-row TSV.

AlphaMissense pathogenicity distribution among joint-matched mutations:
- min 0.045, median 0.245, mean 0.360, max 1.000
- "likely_pathogenic" (≥ 0.564): subset analysed separately but not
  reported here — see `alphamissense_panel_run.py` for the proxy-only
  baseline (deprecated).

---

## Three weight schemes

All three schemes consume the same cohort, same seeds (5), same ctDNA
fractions, same per-seed simulation pipeline
(`real_tcga_validation.simulate_cfdna_from_real`). Deltas are attributable
to the weighting change alone.

### Scheme 1 — `weight_cadd_only` (published winner)

```text
w_i = 1 if CADD_PHRED_i ∈ top-K=200 of the patient's CADD scores
w_i = 0 otherwise
```
- **Per-patient top-K kept:** median 200, min 109, max 200 (5 of 20 patients
  had fewer than 200 CADD-scored mutations).

### Scheme 2 — `weight_alphamissense_only`

```text
w_i = 1 if AM_pathogenicity_i ∈ top-K=200 of the patient's AM scores
w_i = 0 otherwise
```
- **Per-patient top-K kept:** median 177.5, min 90, max 200. Some patients
  have fewer than 200 AM-matched missense mutations because the joint
  match rate is only 80.8 %.

### Scheme 3 — `weight_combined`

```text
w_i = (CADD_PHRED_i / max_CADD) * (AM_pathogenicity_i / max_AM)
w_i = 0 if either score is missing (no imputation)
```
- **Per-patient non-zero weights:** median 177.5, min 90, max 524 (a few
  patients have many low-CADD, mid-AM mutations that survive after
  multiplication).
- `max_CADD` = global max of matched CADD PHRED ≈ 61
- `max_AM` = global max of matched AM pathogenicity = 1.0

---

## Why did combined lose?

1. **Information loss in multiplication.** When both scores are continuous
   on different scales, the product behaves like an AND rather than an OR.
   A high-CADD mutation with no AM score is silently dropped; a high-AM
   mutation with no CADD score is silently dropped. CADD-only Top-K keeps
   high-CADD mutations regardless of AM availability, and AM-only Top-K
   keeps high-AM regardless of CADD.

2. **AM bias toward known pathogenic hotspots.** The 3,946 joint-matched
   mutations have median AM = 0.245, but the head 5 % of mutations
   (AM ≥ 0.564) are concentrated in recurrent oncogenes. Multiplying by
   CADD does not help find *new* pathogenic-looking variants; it just
   re-weights the known hotspots.

3. **CADD-only is already a strong prior at 0.1 % ctDNA.** Sens@99 % of
   0.91 leaves little headroom; the only scheme that improved over it
   in this experiment was the AM-only at the easier fractions (it tied).
   The 80.8 % match rate of the combined scheme is too restrictive at
   ultra-low VAF.

---

## Methods

### AlphaMissense data
- **Source:** Cheng J et al., *Science* 381, eadg7492 (2023).
  https://doi.org/10.1126/science.adg7492
- **TSV:** `AlphaMissense_hg38.tsv.gz` (643 MB compressed, 5.5 GB
  uncompressed, 71.7 M missense variants) downloaded once into
  `/Users/hermes/.cache/mrnavax/`.
- **License:** CC BY-NC-SA 4.0 (non-commercial).
- **Columns used:** `uniprot_id` (5), `protein_variant` (7, e.g. `"V2L"`),
  `am_pathogenicity` (8).
- **Streaming approach:** scan the gzipped TSV once, retain only the
  ~4 K wanted keys. ~36 s wall-time per fresh scan.

### CADD data (from prior work)
- **Source:** Kircher M et al., *Nat Genet* 46:310-315 (2014).
  v1.6 GRCh38 gnomAD-r3 SNV TSV + v1.7 GRCh38 gnomAD-r4 InDel TSV.
- **License:** CC BY-NC-SA 4.0 (non-commercial).
- **Match rate:** 4,882 / 5,738 = 85.1 % of TCGA-LUAD cohort.

### Cohort
- 20 TCGA-LUAD patients, 5,738 total mutations, 4,066 missense SNVs with
  Uniprot + protein position + amino-acid change (downloaded per-case
  MAFs from the GDC API; see `scripts/enrich_cohort_with_amino_acids.py`).

### Panel simulation
- `real_tcga_validation.simulate_cfdna_from_real(...)` with
  `cfdna_depth=5000`, `bg_error_rate=0.002`, ctDNA fractions
  `[0.1, 0.05, 0.01, 0.005, 0.001]`, seeds `[42, 123, 456, 789, 1024]`.
- Aggregation: `Σ w_i · LLR_i` per patient → ROC across the 20-patient
  pool.

### Determinism
- All three weight schemes are deterministic functions of the matched
  scores. The 5 seeds control the cfDNA read sampling only.
- No parameter is learned from the data.

---

## Files

| File | Purpose |
|------|---------|
| `scripts/cadd_alphamissense_combined.py` | Experiment entry point |
| `scripts/enrich_cohort_with_amino_acids.py` | One-shot GDC MAF enrichment |
| `results/cohort_with_amino_acids.jsonl` | 4,613 missense SNVs × AA info |
| `results/cadd_alphamissense_combined.json` | Full results (5 schemes × 5 TFs × 5 seeds) |
| `/Users/hermes/.cache/mrnavax/AlphaMissense_hg38.tsv.gz` | Source TSV (643 MB) |
| `/Users/hermes/.cache/mrnavax/alphamissense_index.pkl` | Pre-built AM index (260 KB) |

---

## Reproducing

```bash
# 1. (One-time) Download the AlphaMissense TSV.gz into ~/.cache/mrnavax/
curl -L -o ~/.cache/mrnavax/AlphaMissense_hg38.tsv.gz \
    https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz

# 2. (One-time) Enrich the cohort with amino-acid info from GDC
cd /Users/hermes/deepcatch
env -u PYTHONPATH ./.venv/bin/python scripts/enrich_cohort_with_amino_acids.py

# 3. Run the combined weighting experiment
env -u PYTHONPATH ./.venv/bin/python scripts/cadd_alphamissense_combined.py
```

Wall-time: ~8 minutes end-to-end on an M-series Mac (stream AM TSV 36 s +
3 schemes × 5 ctDNA × 5 seeds × ~27 s/scheme/tf ≈ 5 min).

---

## Honest caveats

1. **AlphaMissense is missense-only.** LoF (nonsense, frameshift,
   splice-site) variants are dropped from schemes 2 and 3 because the
   AM TSV has no scores for them. Scheme 1 (CADD-only) keeps them.
   This is one structural reason scheme 3 has fewer non-zero weights
   than scheme 1.

2. **Match rate ceiling.** Even with the full TSV in hand, the joint
   match rate plateaus at 80.8 % of CADD-matched (and 68.8 % of cohort).
   Private somatic mutations with novel protein positions or isoforms
   will never get an AM score.

3. **No imputation in the combined scheme.** A missing CADD OR AM score
   means the locus is silently dropped (weight = 0). This is honest but
   reduces the panel size for patients with sparse AM coverage.

4. **Single-cohort result.** All 20 patients are TCGA-LUAD. The
   combined-scheme underperformance should be re-tested on a second
   cancer type (e.g. COADREAD, BRCA) before generalising.
