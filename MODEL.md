# Model Card — DeepCatch v2.2.0

## Overview

DeepCatch is an open-source pipeline for detecting cancer from cell-free
DNA (cfDNA) sequencing data. It provides **mutation-informed** detection
(variant calling + panel-LLR + strand-score + Fisher's exact) and,
through its cross-repo adapter, a **tumor-naive** detection channel
(fragmentomics: 5Mb + 100 kb short/long ratio, median-normalized
coverage, and full fragment-size-distribution histogram).

**This is research-grade software, not a clinical diagnostic tool.**
It has not been validated on a held-out clinical cohort. See
[Limitations](#limitations) below.

## Model details

| | |
|---|---|
| **Model name** | DeepCatch v2.2.0 |
| **Version** | 2.2.0 (release tag, commit `0842600`) |
| **Date** | 2026-08-11 |
| **Authors** | Royce (Independent Researcher, ORCID [pending]) |
| **License** | MIT |
| **Repository** | https://github.com/rollroyces/deepcatch |
| **Companion pipeline** | [cfdna-fragmentomics-pipeline](https://github.com/rollroyces/cfdna-fragmentomics-pipeline) (the source of the tumor-naive channel) |
| **Companion data** | 627 cross-study FinaleDB samples (Jiang 2015 + Cristiano 2019) |
| **Reference paper** | `paper/PAPER.md` (242 lines) and `paper/paper.tex` (177 lines) |

## Intended use

- **Primary**: research-grade benchmarking of mutation-informed ctDNA
  detection methods against simulated plasma at sub-1% VAF.
- **Secondary**: demonstration of cross-repo integration with a
  fragmentomics pipeline (pairing via file-format contract, not
  import-time dependency).
- **Tertiary**: an extensible framework for adding new ctDNA
  detection modalities (the multimodal_fusion module is the
  integration point).

## Out-of-scope use

- **Clinical diagnostics.** This model has not been validated on
  held-out clinical cohorts. The TCGA-LUAD headline number is from
  Poisson-sampled simulated plasma reads, not real patient samples.
- **Population screening.** The 99%-specificity operating point is
  in-sample and not externally validated.
- **Non-LUAD cancer types.** The headline number is for TCGA-LUAD
  only. Other cancer types would require retraining.

## Training data

| Channel | Source | Cohort | Type |
|---|---|---|---|
| Mutation-informed (panel LLR, Fisher, strand) | Real TCGA MAFs (GDC) + Poisson-sampled reads | 20 TCGA-LUAD patients, 5,738 mutations | Real mutations, simulated plasma |
| Tumor-naive (5-channel fragmentomics) | FinaleDB (live cfDNA WGS) | 363 cancer + 264 healthy across Jiang 2015 + Cristiano 2019 | Real cfDNA (low-pass WGS) |

Mutation channel data: fetched from the GDC open-access API on first
run; cached in `validation/tcga/tcga_cache/`. Synthetic data is
**deliberately refused** — see `real_tcga_validation.py` for the
fail-loud guard.

Tumor-naive channel data: pre-computed by the companion pipeline
(per-study z-score-harmonized 5-channel profile). The on-disk
format is documented in `src/fragmentomics/tumor_naive_adapter.py`.

## Performance

### Mutation-informed (5-seed, real TCGA mutations + simulated plasma)

| VAF | Panel LLR AUC | Fisher AUC | Strand AUC | Sens@95% | Sens@99% |
|---|---|---|---|---|---|
| 5% | 1.000 | 0.999 | 0.999 | 1.000 | 1.000 |
| 1% | 0.998 | 0.998 | 0.997 | 1.000 | 1.000 |
| 0.5% | 0.999 | 0.999 | 0.996 | 0.999 | 1.000 |
| **0.1%** | **0.921** | **0.834** | **0.831** | **0.770** | **0.460** |

Source: `results/real_tcga_validation.json` (commit `0842600`).

### Tumor-naive (5-seed, real FinalDB cross-study)

| Cohort | Headline AUC (5-seed) | Sens@95% | Sens@99% |
|---|---|---|---|
| Single-study Jiang 2015 (121 samples) | 0.9716 ± 0.003 | 0.894 | 0.811 |
| Cross-study pan-cancer (627 samples, harmonized) | 0.9745 ± 0.002 | 0.888 | 0.774 |

Source: `cfdna-fragmentomics-pipeline/scripts/honest_benchmark.py`.

### Fusion (mutation-informed + tumor-naive, 10-seed cross-study)

| Strategy | AUC | Sens@95% | Sens@99% |
|---|---|---|---|
| Tumor-naive only | 0.9743 ± 0.002 | 0.883 | 0.760 |
| Mutation-only (synthetic, AUC 0.92) | 0.9242 | 0.656 | 0.336 |
| **Naive average** | **0.9886** | **0.927** | **0.859** |
| **LR fusion** | **0.9887** | **0.937** | **0.845** |

**Paired t-test (10 seeds)**: LR-fusion AUC − tumor-naive AUC = **+0.0143**
(t = 31.96, p < 0.0001, 95% bootstrap CI [0.0135, 0.0152]).
**DeLong p < 0.0015 on every one of 5 individual seeds**.

Source: `src/fragmentomics/fusion_ablation.py` (DeepCatch side).

## Limitations

1. **No held-out clinical validation.** Every headline number is
   from in-sample 5-fold CV. The mutation-informed channel is
   real-TCGA mutations + simulated plasma reads. The tumor-naive
   channel is real FinaleDB samples but the train/test split is
   pooled OOF on the same cohort, not a held-out second study.
2. **Synthetic mutation channel in the fusion experiment.** The
   fusion result is calibrated against a synthetic mutation channel
   matching the TCGA-LUAD headline. A real mutation channel on real
   plasma (with all the noise of real variant calling) might be
   weaker; the calibration-sensitivity sweep in BENCHMARK.md shows
   that below mutation AUC 0.80, fusion provides no benefit.
3. **Cohort size.** The 627-sample cross-study cohort is an
   improvement over DELFI's 481 but still small relative to clinical
   validation cohorts (typically 1000+).
4. **Healthy controls.** 264 controls vs DELFI's 245 — the 99%-spec
   operating point is statistically stable but not as tight as a
   prospective study would require.
5. **RNA edits, methylation, and other modalities are not included.**
   The fusion experiment is mutation + fragmentomics only.
6. **The TCGA-LUAD headline is simulation-based.** "TCGA-LUAD"
   means real TCGA tumor mutations from the GDC API; the "plasma"
   reads are Poisson-sampled at the stated tumor fraction. This is
   a dilution/spike-in benchmark, not a clinical plasma validation.
   The README flags this explicitly.

## Ethical considerations

- False-positive ctDNA detection could lead to unnecessary follow-up
  procedures. The Sens@95% / Sens@99% operating points are
  appropriate for screening workflows, not for confirmatory
  diagnosis.
- The training mutation set is from TCGA, which has known
  demographic bias (predominantly European-ancestry patients). The
  model may not generalize to all populations.
- The FinaleDB cohorts (Jiang 2015, Cristiano 2019) have their own
  cohort biases; transferability to other populations is not
  validated.

## How to use

```bash
# Install (no GPU required)
pip install -e .

# Real-TCGA mutation-informed validator (downloads GDC data on first run)
python real_tcga_validation.py --help

# Cross-repo tumor-naive adapter (requires the pipeline's features/)
deepcatch-tumornaive \
  --features-dir ../cfdna-fragmentomics-pipeline/data/features \
  --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
  --seeds 5 --pca-n 200

# Mutation-informed + tumor-naive fusion ablation
deepcatch-fusion \
  --features-dir ../cfdna-fragmentomics-pipeline/data/features \
  --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
  --seeds 5 --pca-n 200

# Decision curve + per-specificity operating table (clinician-facing)
deepcatch-decisioncurve \
  --features-dir ../cfdna-fragmentomics-pipeline/data/features \
  --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
  --out decision_curve.json

# Synthetic-cohort AUC gate (no network, <30 s)
python scripts/adapter_auc_gate.py
```

## Citation

```bibtex
@software{deepcatch_v2_2_0,
  author = {Royce},
  title  = {DeepCatch: Mutation-informed Ultra-Sensitive cfDNA Detection},
  year   = {2026},
  url    = {https://github.com/rollroyces/deepcatch},
  version = {2.2.0},
}
```

## Versioning

This model card applies to **DeepCatch v2.2.0** (commit `0842600`).
Subsequent versions may update the headline numbers; check the
`results/` directory and the README for the latest validated metrics.
