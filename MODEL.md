# Model Card — DeepCatch v2.3.0

> **⚠️ STATUS: research-use-only software benchmark.**
> NOT clinically validated. NOT a medical device. NOT FDA-approved.
> NOT intended for clinical diagnosis, treatment selection, screening,
> or patient stratification.
> For methods research only.

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
| **Model name** | DeepCatch v2.3.0 |
| **Version** | 2.3.0 (commit `7f9613b` on `main`) |
| **Date** | 2026-09-23 |
| **Authors** | Yu Ching Lam (Independent Researcher, ORCID [0009-0008-9113-769X](https://orcid.org/0009-0008-9113-769X)) |
| **License** | MIT |
| **Repository** | https://github.com/rollroyces/deepcatch |
| **Companion pipeline** | [cfdna-fragmentomics-pipeline](https://github.com/rollroyces/cfdna-fragmentomics-pipeline) (the source of the tumor-naive channel) |
| **Companion data** | 627 cross-study FinaleDB samples (Jiang 2015 + Cristiano 2019) |
| **Reference paper** | `paper/paper.tex` — the manuscript (compiled to `paper/biorxiv_submission_v2.3.0.pdf`); `paper/PAPER.md` — the earlier panel-only revision, kept for reference |

## Intended use

- **Primary**: research-grade benchmarking of mutation-informed ctDNA
  detection methods against simulated plasma at sub-1% VAF.
- **Secondary**: demonstration of cross-repo integration with a
  fragmentomics pipeline (pairing via file-format contract, not
  import-time dependency).
- **Tertiary**: an extensible framework for adding new ctDNA
  detection modalities (the multimodal_fusion module is the
  integration point).

## Intended users

- ML researchers working on cfDNA fragmentomics or ctDNA detection.
- Fragmentomics / epigenomics methods developers.
- Bioinformatics scientists evaluating open-data benchmarks.

## NOT intended for

- **Clinical diagnosis** of any cancer type.
- **Treatment selection** (e.g. therapy assignment based on ctDNA
  detection).
- **Population screening** (no validated sensitivity / specificity
  at the prevalence typical of a screening cohort, ~0.4%).
- **Patient stratification** (e.g. MRD-positive vs MRD-negative
  treatment escalation decisions).
- **Any clinical decision-making** — DeepCatch is research-use-only
  software, not a medical device, and is NOT FDA-approved.

## Training data

| Channel | Source | Cohort | Type |
|---|---|---|---|
| Mutation-informed (panel LLR, Fisher, strand) | Real TCGA MAFs (GDC) + Poisson-sampled reads | 20 TCGA-LUAD patients, 5,738 mutations | Real mutations, simulated plasma |
| Tumor-naive (5-channel fragmentomics) | FinaleDB (live cfDNA WGS) | 363 cancer + 264 healthy across Jiang 2015 + Cristiano 2019 | Real cfDNA (low-pass WGS) |
| Foundation-model pretraining | Open FinaleDB cache (publications 6 + 8) | 200 samples, `data/finaledb_pretrain_cohort_PRODUCTION.npz` | Real cfDNA WGS (single open cohort) |

Mutation channel data: fetched from the GDC open-access API on first
run; cached in `validation/tcga/tcga_cache/`. Synthetic data is
**deliberately refused** — see `real_tcga_validation.py` for the
fail-loud guard.

Tumor-naive channel data: pre-computed by the companion pipeline
(per-study z-score-harmonized 5-channel profile). The on-disk
format is documented in `src/fragmentomics/tumor_naive_adapter.py`.

Foundation pretraining data: 200 samples from `data/finaledb_pretrain_cohort_PRODUCTION.npz`
(real FinaleDB WGS, see [`docs/PRETRAINING.md`](docs/PRETRAINING.md)). The
pre-fix `*SYNTHETIC_v0*` artifacts in `checkpoints/` and `data/` are
**NOT real-data-trained** — they bypassed the real cohort assembly
under the synthetic-bypass bug documented in
[`docs/PRETRAIN_BUG.md`](docs/PRETRAIN_BUG.md). Do not load them as
"real-data-derived" checkpoints.

## Validation data

Same as training data: **internal cross-study on the same open
cohort, not external.** The cross-study benchmark pools two
published FinaleDB WGS studies (Jiang 2015 + Cristiano 2019) into a
single n=627 cohort with per-study z-score harmonization. There is
**no held-out second cohort**; pooled OOF on the same cohort is the
honest validation surface.

The healthy controls in the cross-study benchmark are **real
healthy-donor plasma**: FinaleDB records every training-cohort control
as `disease: Healthy`, `tissue: blood plasma` (262 Cristiano + 32
Jiang donors). What the harmonisation does *not* fix is the class mix:
the per-study z-score uses a class-mixed per-study mean, and the two
training studies have very different class proportions (Jiang is ~70%
HCC, Cristiano is mostly healthy).

## Performance metrics

Detailed per-cancer sens@spec with DeLong 95% CIs, pooled AUCs, and
true-confound controls are in
[`docs/CROSS_STUDY_BENCHMARK.md`](docs/CROSS_STUDY_BENCHMARK.md). Headline
numbers below for quick reference.

### Mutation-informed (5-seed, real TCGA mutations + simulated plasma)

| VAF | Panel LLR AUC | Fisher AUC | Strand AUC | Sens@95% | Sens@99% |
|---|---|---|---|---|---|
| 5% | 1.000 | 0.999 | 0.999 | 1.000 | 1.000 |
| 1% | 0.998 | 0.998 | 0.997 | 1.000 | 1.000 |
| 0.5% | 0.999 | 0.999 | 0.996 | 0.999 | 1.000 |
| **0.1%** | **0.921** | **0.834** | **0.831** | **0.770** | **0.460** |

Source: `results/real_tcga_validation.json`.
**Honest framing**: real TCGA mutations + Poisson-sampled reads at
the stated tumor fraction — a spike-in/dilution benchmark, **not**
a clinical plasma validation.

### Tumor-naive (5-seed, real FinaleDB cross-study)

| Cohort | Headline AUC (5-seed) | Sens@95% | Sens@99% |
|---|---|---|---|
| Single-study Jiang 2015 (121 samples) | 0.9661 ± 0.0031 | — | — |
| Cross-study pan-cancer (627 samples, harmonized) | 0.9684 ± 0.0027 | 0.865 | 0.700 |
| Cross-study pan-cancer, legacy no-GC correction | 0.9770 ± 0.0013 | 0.912 | 0.785 |

Source: `results/cross_study_finallydb_gc_corrected.json` (rows 1 and 2) and
`results/cross_study_finallydb.json` (row 3). The per-publication row carries
no sensitivity in the artifact, so those cells are left empty rather than
filled from another run. The two cross-study rows are **different
benchmarks** — the GC-corrected one is the headline, the legacy one is what
the same pipeline gives without the GC correction — and they are not
interchangeable.

(The numbers this table used to carry, 0.9716 and 0.9746, came from the
sibling repository's `scripts/honest_benchmark.py` before its PCA was
seeded; 0.9746 is one seed's value out of the legacy artifact's five, not
the pooled mean. The sibling's benchmark moves once it is seeded, so the
table is now written from this repository's own artifacts.)

### Fusion (mutation-informed + tumor-naive, 10-seed cross-study)

| Strategy | AUC | Sens@95% | Sens@99% |
|---|---|---|---|
| Tumor-naive only | 0.9740 ± 0.0013 | 0.880 | 0.742 |
| Mutation-only (synthetic, AUC 0.92) | 0.9020 | 0.595 | 0.267 |
| **Naive average** | **0.9884** | **0.927** | **0.841** |
| **LR fusion** | **0.9884** | **0.932** | **0.833** |

**Paired t-test (10 seeds)**: LR-fusion AUC − tumor-naive AUC = **+0.0143**
(t = 49.22, p < 0.0001, 95% CI [0.0138, 0.0149]).
**DeLong p < 0.0029 on every one of the 10 individual seeds**.

Source: `src/fragmentomics/fusion_ablation.py` (DeepCatch side); artifact
`results/fusion_ablation_isotonic.json`; narrative
[`docs/FUSION_ISOTONIC.md`](docs/FUSION_ISOTONIC.md).
**Honest framing**: the "mutation channel" in this fusion is a
**synthetic** score calibrated to AUC 0.92, not a measurement from
the same plasma as the fragmentomics channel. This is a what-if
experiment, not an end-to-end plasma measurement.

### Foundation-model smoke (3-seed, 20 paired, Audit-2 honest framing)

| Variant | AUC | Honest framing |
|---|---|---|
| Foundation (3-ensemble tiny transformer) | 0.55 ± 0.01 | Honest per-patient GroupKFold, real TCGA panel-LLR + real per-patient frag channel, pair-broken shuffled-label control |
| sklearn LR baseline (same channels) | 0.91 | Foundation is honestly out-performed by LR on this cohort |
| Shuffled-label foundation | 0.77 | Pair-broken shuffled control |
| delta_auc_normalized | −3.04 | foundation < shuffled (negative result, reported honestly) |

Source: `results/foundation_real_smoke.json`, commit `1d97f78` (Audit-2).
The pre-fix headline "foundation AUC 0.93 ± 0.03" was inflated by
three P0 bugs (CV patient-identity leak, synthetic-Gaussian
substitution of the real frag channel, no-op shuffled-label
control). See `AUDIT_2_FINDINGS.md` for the full bug writeup.

### Late fusion — the estimator that does work at this sample size

The neural encoder is the wrong model class for n=40. Two late-fusion
estimators were measured on the identical protocol (real TCGA-derived
channels, per-patient GroupKFold, LR baseline re-computed in the same
run):

| Estimator | AUC (3 seeds) | permuted-label control | arm-swap control |
|---|---|---|---|
| foundation encoder | 0.554 ± 0.005 | — | — |
| `lr_baseline` (re-computed) | 0.9042 ± 0.0174 | 0.569 | 0.383 |
| `naive_average` | 0.9150 ± 0.0178 | 0.614 | 0.513 |
| **`regularized_late_fusion`** | **0.9825 ± 0.0094** | 0.573 | 0.436 |
| **`bayesian_model_average`** | **0.9650 ± 0.0089** | 0.511 | 0.413 |

Source: `results/fusion_regularized_smoke.json`
(`scripts/fusion_regularized_smoke.py`); implementation in
`src/foundation/regularized_fusion.py`. The re-computed LR baseline
(0.9042) and equal-weight average (0.9150) reproduce the published
0.910 / 0.915, so the comparison is within-run.

**The control had to be repaired to mean anything.** The audit-2
"pair-broken shuffle" permuted the positive and negative index sets
separately and reassigned 1s and 0s to those same sets, which
reproduces `y_true` exactly — it was a no-op, and the "shuffled"
numbers in the audit-2 table are the honest numbers. The corrected
control permutes labels inside each fold; its null on this paired
design is **0.4063 ± 0.1143, not 0.5** (measured over 150 permutations
seed), because every patient contributes one positive and one negative
arm whose panel scores differ by ~+127 with the positive arm higher in
100% of patients, so a model fitted to random labels still recovers
that within-pair axis and lands below 0.5. A control below 0.5 is
therefore not on its own evidence of leakage; it must be read against
the measured null.

**Scope**: 40 samples from 20 patients, simulated plasma at TF=0.001,
one cancer type. This is an estimator-class result, not a clinical
performance claim.

### External validation on held-out public publications

Training uses FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano
2019), n=627. The external cohort is publications 1, 7 and 9, fetched
from the public FinaleDB API and S3 (no IRB, no data-use agreement);
278 of 283 samples were extracted. Source:
`results/external_validation_finaledb.json`; narrative
[`docs/EXTERNAL_VALIDATION.md`](docs/EXTERNAL_VALIDATION.md).

| Publication | Assay | n (ca/non-ca) | External AUC [DeLong 95% CI] |
|---|---|---|---|
| 1 — Snyder 2016 | WGS | 52 (43/9) | 0.708 [0.5148, 0.9012] |
| 7 — Sun 2019 | WGS | 31 (18/13) | **0.9444** [0.8731, 1.0000] |
| 9 — Adalsteinsson 2017 | WES | 239 (239/0) | case-only: sensitivity 0.6569 at the training-derived 99%-specificity threshold |
| **pooled over publications carrying both classes** | — | 83 (61/22) | **0.7385** [0.6281, 0.8488] |

Transfer between whole-genome studies is real but far from uniform: the
held-out colorectal publication transfers well (0.9444; colorectal
one-vs-rest 0.9182), a second whole-genome publication covering thirteen
tumour types transfers less well (0.708), and the whole-exome one only
partly (breast 0.6938, prostate 0.7127 — above chance but well below the
whole-genome studies). Per-publication z-score harmonization removes a
study's location/scale shift but not a change of assay platform. Publication 9 ships 1,508 cancer samples
and zero controls, so it is reported as sensitivity at a
training-fixed threshold rather than folded into an AUC denominator.

### Leave-one-study-out transfer

Holding a publication out of the *training* cohort is a different
measurement: `scripts/cross_study_loso.py` gives 0.6805
[0.6344, 0.7267] for train-Jiang → test-Cristiano against in-domain
0.9789 / 0.969. The asymmetry is informative — Jiang is HCC-only, so
training on it yields an HCC-specific boundary, while a pan-cancer
training set yields a generic one (train-Cristiano → test-Jiang
0.9361). The pooled 0.9684 does **not** survive a publication
holdout and must not be read as a generalisation claim.

### Honest negative results preserved

- **Continuous per-cancer weighting** (CADD-style linear / sigmoid
  weights on channels) — regresses vs hard top-K.
- **CADD + AlphaMissense multiplicative weight** — underperforms CADD
  alone at all ctDNA fractions.
- **Driver-only panel** (TP53 + KRAS + EGFR + …) — Sens@99% = 0.19
  (catastrophic).
- **Isotonic post-hoc calibration of LR fusion** — regresses on every
  seed (mean Δ AUC = −0.0035). See [`docs/FUSION_ISOTONIC.md`](docs/FUSION_ISOTONIC.md).
- **SensAtSpecLoss ablation** — no detectable effect on the n=20
  foundation cohort. See [`docs/SENS_AT_SPEC_ABLATION.md`](docs/SENS_AT_SPEC_ABLATION.md).
- **Sparse-aware projection ablation** — null effect on this cohort.
  See [`docs/SPARSE_AWARE_ABLATION.md`](docs/SPARSE_AWARE_ABLATION.md).
- **Foundation vs LR** — foundation honestly loses to LR on the
  n=20 cohort. See `AUDIT_2_FINDINGS.md`.
- **Methylation (V3) channel** — measured on 31 real FinaleDB samples and
  reported as a **null**: the decoded β summary reaches LOO AUC 0.9444,
  which is exactly the AUC of the covered-CpG count and of the
  observations-per-CpG count alone, and at the 103,863 loci where both arms
  have coverage the median cancer-minus-control difference is exactly 0.0
  (each locus tested on the samples that cover it; permutation p = 1.0). The libraries are size-selected at ~300 bp and
  deliver 1.09 observations per CpG, so the digestion signature FinaleMe
  reads is absent. The extraction path runs on real public data; the
  channel is **not** a validated feature. See
  [`docs/METHYLATION_PIPELINE.md`](docs/METHYLATION_PIPELINE.md).

## Limitations

1. **Cohort sizes are small.**
   - n=20 paired synthetic for the foundation smoke (TCGA-LUAD).
   - n=627 for the open-data cross-study benchmark (FinaleDB pan-cancer).
   - n=129 processed frequency vectors for the Jiang motif benchmark.
   These are appropriate for methods research but are far below the
   cohort sizes typical of clinical validation (1,000+).

2. **The harmonisation is not class-balanced.** The per-study z-score
   uses a class-mixed per-study mean, and the two training studies have
   very different class proportions (Jiang is ~70% HCC, Cristiano is
   mostly healthy), so a study's mean partly encodes its case mix. A
   per-study-per-class scaler would be more rigorous. The control arm
   itself is sound: FinaleDB records every training-cohort control as
   `disease: Healthy`, `tissue: blood plasma` (262 Cristiano + 32 Jiang
   donors), i.e. real healthy-donor plasma.

3. **No IRB approval, no prospective collection, no clinical
   outcomes.** The validation is on open-data WGS fragments from
   public repositories.

4. **External validation exists, but it is public-data external
   validation, not clinical validation.** Publications 1/7/9 were never
   trained on and give AUC 0.9444 (pub 7, WGS), 0.708 (pub 1, WGS),
   0.7385 pooled over the publications carrying both classes, and 0.6569
   sensitivity on the whole-exome publication. What it still lacks: prospective follow-up,
   stage labels, and a large control arm — the held-out non-cancer arm
   is 22 public samples. A leave-one-study-out run inside the training
   cohort gives 0.6805, so the pooled 0.9684 is a within-cohort
   separation number, not a generalisation claim.

5. **Mutation-informed channel in the fusion experiment is
   synthetic.** The fusion AUC 0.9884 is a what-if pairing of a
   real fragmentomics channel with a synthetic mutation channel
   calibrated to AUC 0.92 — not a real measurement of both channels
   on the same plasma. Pairing on real plasma is the next phase.

6. **The TCGA-LUAD headline is simulation-based.** "TCGA-LUAD"
   means real TCGA tumor mutations from the GDC API; the "plasma"
   reads are Poisson-sampled at the stated tumor fraction. This is
   a dilution/spike-in benchmark, not a clinical plasma validation.

7. **Batch effects are not exhaustively characterized.** The
   harmonization check ([`docs/HARMONIZATION.md`](docs/HARMONIZATION.md))
   is a synthetic fixture with NEUTRAL verdict; on a real
   multi-study plasma pool, larger batch effects are expected.

8. **Pretraining is on 200 samples from a single open cohort.**
   Foundation-model pretraining uses 200 samples from
   `data/finaledb_pretrain_cohort_PRODUCTION.npz`. The pre-fix
   `*SYNTHETIC_v0*` artifacts are not real-data-trained (see
   [`docs/PRETRAIN_BUG.md`](docs/PRETRAIN_BUG.md)).

9. **The 99%-specificity operating point is in-sample.** Not
   externally validated, not calibrated against a screening cohort
   with prevalence ~0.4%.

10. **No FDA pathway, no regulatory submission.** DeepCatch is
    research-use-only software.

## Honest evaluation (target venue)

The intended submission venue is a **methods-focused** journal such
as *PLOS Computational Biology* or *Bioinformatics Advances* — NOT
a clinical-validation venue such as *Nature Medicine* or *Cancer
Discovery*. DeepCatch does not have the prospective-cohort data,
IRB approval, or external validation required for a clinical-
validation paper.

## Ethical considerations

- False-positive ctDNA detection could lead to unnecessary follow-up
  procedures. The Sens@95% / Sens@99% operating points are
  appropriate for **methods research only**, not for confirmatory
  diagnosis.
- The training mutation set is from TCGA, which has known
  demographic bias (predominantly European-ancestry patients). The
  model may not generalize to all populations.
- The FinaleDB cohorts (Jiang 2015, Cristiano 2019) have their own
  cohort biases; transferability to other populations is not
  validated.
- **Do not use this software to make medical decisions.** It is
  research-use-only.

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

# Decision curve + per-specificity operating table (clinician-facing — research use only)
deepcatch-decisioncurve \
  --features-dir ../cfdna-fragmentomics-pipeline/data/features \
  --labels      ../cfdna-fragmentomics-pipeline/data/features/labels_cross_study.tsv \
  --out decision_curve.json

# Synthetic-cohort AUC gate (no network, <30 s)
python scripts/adapter_auc_gate.py
```

## Citation

```bibtex
@software{deepcatch_v2_3_0,
  author = {Yu Ching Lam},
  title  = {DeepCatch: Mutation-informed Ultra-Sensitive cfDNA Detection
            (research-use-only open-data methods benchmark)},
  year   = {2026},
  url    = {https://github.com/rollroyces/deepcatch},
  version = {2.3.0},
  note   = {Research-use-only open-data benchmark. NOT clinically validated.
            NOT FDA-approved. NOT a medical device. For methods research only.}
}
```

## Versioning

This model card applies to **DeepCatch v2.3.0** (commit `7f9613b` on
`main`). Subsequent versions may update the headline numbers; check
the `results/` directory and the README for the latest validated
metrics.
