# Cross-Platform Methylation + Fragmentomics Fusion (methodology demo)

This document describes the cross-platform fusion methodology demo: combining
a fragmentomics channel from this repo with a methylation channel from the
sibling [`deepcatch-methylation`](https://github.com/deepcatch/deepcatch-methylation)
repo.

## TL;DR

| Channel | Source | Cohort | Published AUC | n |
|---------|--------|--------|---------------|---|
| Fragmentomics | this repo, `results/per_cancer_sens_at_spec.json` HCC_J OvR row | Jiang 2015 cfDNA WGS (HCC plasma) | **0.7532** (Delong CI 0.688–0.819) | 89 cancer + 264 healthy |
| Methylation | sibling repo, `results/multi_cancer_baseline.json` TCGA-LIHC row | GDC TCGA-LIHC tissue HM450 β-values | **0.9722 ± 0.0088** (5-seed × 5-fold CV) | 12 tumor + 6 normal (n=18) |

**The two channels are from disjoint cohorts. This is a methodology demo, not
a clinical fusion.** Per-sample scores for both channels are synthesized from
the published aggregate metrics (calibrated Gaussian + sigmoid). The demo
verifies that the three fusion strategies run end-to-end and produce valid
AUCs in the expected range.

## Why this demo exists

The sibling repo's methylation pipeline is built around TCGA-LIHC tissue
methylation arrays (HM450), while the fragmentomics pipeline in this repo is
built around cfDNA WGS (Cristiano 2019 + Jiang 2015). These are the only
real public cfDNA-vs-methylation cohorts available at scale, but they do not
overlap in patients, in tissue type (plasma cfDNA vs solid-tissue β-values),
or in assay chemistry (WGS vs HM450 array).

A *real* cross-platform fusion would require paired fragmentomics +
methylation on the same patients (the so-called "DELFI + methylation" cohort
that several biotechs are building but have not yet published at scale).
Until that exists, any cross-platform fusion claim must be marked as a
methodology demonstration, not a clinical result.

This script (`scripts/cross_platform_methylation_fusion.py`) shows what the
fusion code path would look like if such a cohort existed, using the
currently-available published numbers as inputs.

## What the script does

1. Loads the HCC_J per-cancer OvR row from
   `results/per_cancer_sens_at_spec.json` (n=353, target AUC=0.7532).
2. Loads the TCGA-LIHC methylation baseline row from
   `../deepcatch-methylation/results/multi_cancer_baseline.json` (n=18,
   target AUC=0.9722).
3. Synthesizes per-sample scores for a 353-sample HCC_J cohort calibrated to
   the published fragmentomics AUC (calibrated Gaussian, mean separation
   mu = sqrt(2) · Φ⁻¹(AUC), then sigmoid). The "methylation column" is
   synthesized with the same calibration targeting the LIHC AUC (0.972).
   This represents the per-sample methylation score that a paired cohort
   would produce — i.e. what a real-world DELFI + methylation assay would
   return for those 353 samples.
4. Runs three fusion strategies with 5 seeds × 5-fold CV:
   - **`naive_average`** — arithmetic mean of the two scores.
   - **`logit_average`** — mean of logit(p), resigmoided.
   - **`lr_fusion`** — L2-LR stacking the two scores (C=1.0).
5. Writes `results/cross_platform_methylation_fusion.json` with per-strategy
   AUCs, sens@95, sens@99, and an explicit provenance block calling out the
   disjoint-cohort caveat.

## Honest caveats (also embedded in the JSON output)

- **Disjoint cohorts.** The fragmentomics channel is from Jiang 2015 cfDNA
  WGS (HCC plasma). The methylation channel is from TCGA-LIHC tissue HM450.
  These cohorts do not overlap in patients, tissue type, or assay chemistry.
- **Synthetic per-sample scores.** Both channels' per-sample scores are
  synthesized from the published aggregate metrics. This demonstrates the
  *fusion code path*, not a clinical result.
- **No paired cohort exists publicly as of 2026-Q1.** A true cross-platform
  fusion requires paired fragmentomics + methylation on the same patients.
- The methylation-proxy features in the sibling repo are NOT true
  methylation; they are fragment-length/coverage-derived proxies (see
  `deepcatch-methylation/docs/METHYLATION_PROXY_RESULTS.md`). The TCGA-LIHC
  numbers in `multi_cancer_baseline.json` ARE true tissue HM450 β-values
  and are the more credible methylation reference for this demo.
- The HCC_J per-cancer OvR AUC (0.753) is harder than the pooled
  cross-study AUC (0.9747, see `RESULTS.md`). OvR is per-cancer-vs-all-healthy.

## How to run

```bash
# Full run (5 seeds × 5-fold CV, writes the JSON)
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/cross_platform_methylation_fusion.py

# Smoke run (1 seed only)
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/cross_platform_methylation_fusion.py --seeds 1

# Help / arg list
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    scripts/cross_platform_methylation_fusion.py --help
```

## Tests

`test/test_cross_platform_fusion.py` covers:

- Score synthesizer calibration (`mu = sqrt(2) · Φ⁻¹(AUC)` correctness).
- All three fusion helpers (naive / logit / LR) return unit-interval
  probabilities.
- Synthetic 2-channel fusion beats the weaker of its two channels (a
  working fusion code path test).
- All three strategies produce valid AUCs in (0, 1].
- Provenance block correctly marks the JSON as a methodology demo and
  carries both cohorts' published n_pos / n_neg.
- End-to-end CLI integration runs and writes the JSON.

```bash
env -u PYTHONPATH /Users/hermes/deepcatch/.venv/bin/python \
    -m pytest test/test_cross_platform_fusion.py -v
```

## Related references

- `docs/CROSS_STUDY_BENCHMARK.md` — pooled cross-study fragmentomics numbers
  (n=627, AUC 0.9747).
- `docs/FUSION_ISOTONIC.md` — tumor-naive + mutation-score fusion ablation
  (sibling fusion experiment on the same cohort, real data).
- `deepcatch-methylation/docs/METHYLATION_PROXY_RESULTS.md` — sibling repo's
  methylation-proxy head-to-head numbers and the detailed honest caveats.
- `deepcatch-methylation/BIORXIV_PAPER_FRAGMENTOMICS.md` — sibling paper
  draft with the cross-platform narrative framing.