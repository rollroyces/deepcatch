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
| **Model name** | DeepCatch v2.2.0 |
| **Version** | 2.2.0 (commit `fb1e523` on `main`) |
| **Date** | 2026-09-23 |
| **Authors** | Yu Ching Lam (Independent Researcher, ORCID [0009-0008-9113-769X](https://orcid.org/0009-0008-9113-769X)) |
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
(real FinaleDB WGS, see [`docs/PRETRAINING.md`](PRETRAINING.md)). The
pre-fix `*SYNTHETIC_v0*` artifacts in `checkpoints/` and `data/` are
**NOT real-data-trained** — they bypassed the real cohort assembly
under the synthetic-bypass bug documented in
[`docs/PRETRAIN_BUG.md`](PRETRAIN_BUG.md). Do not load them as
"real-data-derived" checkpoints.

## Validation data

Same as training data: **internal cross-study on the same open
cohort, not external.** The cross-study benchmark pools two
published FinaleDB WGS studies (Jiang 2015 + Cristiano 2019) into a
single n=627 cohort with per-study z-score harmonization. There is
**no held-out second cohort**; pooled OOF on the same cohort is the
honest validation surface.

The healthy controls in the cross-study benchmark are FinaleDB
TF=0 arms (modeled cancer-at-TF=0 read depth), **not real
healthy-donor plasma**.

## Performance metrics

Detailed per-cancer sens@spec with DeLong 95% CIs, pooled AUCs, and
true-confound controls are in
[`docs/CROSS_STUDY_BENCHMARK.md`](CROSS_STUDY_BENCHMARK.md). Headline
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
| Single-study Jiang 2015 (121 samples) | 0.9716 ± 0.003 | 0.894 | 0.811 |
| Cross-study pan-cancer (627 samples, harmonized) | 0.9746 ± 0.002 | 0.888 | 0.774 |

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
design is **0.42 ± 0.12, not 0.5** (measured over 50 permutations per
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
| 1 — Snyder 2016 | WGS | 9 (8/1) | 0.875 [0.630, 1.000] |
| 7 — Sun 2019 | WGS | 31 (18/13) | **0.9915** [0.9712, 1.0000] |
| 9 — Adalsteinsson 2017 | WES | 238 (238/0) | case-only: sensitivity 0.403 at the training-derived 99%-specificity threshold |
| **pooled over publications carrying both classes** | — | 40 (26/14) | **0.8819** [0.7701, 0.9936] |

The assay platform, not the study, is the binding constraint: the
whole-genome held-out publication transfers essentially intact
(colorectal one-vs-rest 0.9474) while the whole-exome one does not
(breast 0.4119, prostate 0.3193 — below chance). Per-publication
z-score harmonization removes a study's location/scale shift but not
a change of assay platform. Publication 9 ships 1,508 cancer samples
and zero controls, so it is reported as sensitivity at a
training-fixed threshold rather than folded into an AUC denominator.

### Leave-one-study-out transfer

Holding a publication out of the *training* cohort is a different
measurement: `scripts/cross_study_loso.py` gives 0.6805
[0.6344, 0.7267] for train-Jiang → test-Cristiano against in-domain
0.9789 / 0.9701. The asymmetry is informative — Jiang is HCC-only, so
training on it yields an HCC-specific boundary, while a pan-cancer
training set yields a generic one (train-Cristiano → test-Jiang
0.9361). The pooled 0.9670 does **not** survive a publication
holdout and must not be read as a generalisation claim.

### Honest negative results preserved

- **Continuous per-cancer weighting** (CADD-style linear / sigmoid
  weights on channels) — regresses vs hard top-K.
- **CADD + AlphaMissense multiplicative weight** — underperforms CADD
  alone at all ctDNA fractions.
- **Driver-only panel** (TP53 + KRAS + EGFR + …) — Sens@99% = 0.19
  (catastrophic).
- **Isotonic post-hoc calibration of LR fusion** — regresses on every
  seed (mean Δ AUC = −0.0035). See [`docs/FUSION_ISOTONIC.md`](FUSION_ISOTONIC.md).
- **SensAtSpecLoss ablation** — no detectable effect on the n=20
  foundation cohort. See [`docs/SENS_AT_SPEC_ABLATION.md`](SENS_AT_SPEC_ABLATION.md).
- **Sparse-aware projection ablation** — null effect on this cohort.
  See [`docs/SPARSE_AWARE_ABLATION.md`](SPARSE_AWARE_ABLATION.md).
- **Foundation vs LR** — foundation honestly loses to LR on the
  n=20 cohort. See `AUDIT_2_FINDINGS.md`.

## Limitations

1. **Cohort sizes are small.**
   - n=20 paired synthetic for the foundation smoke (TCGA-LUAD).
   - n=627 for the open-data cross-study benchmark (FinaleDB pan-cancer).
   - n=129 processed frequency vectors for the Jiang motif benchmark.
   These are appropriate for methods research but are far below the
   cohort sizes typical of clinical validation (1,000+).

2. **Healthy controls in the cross-study benchmark are FinaleDB TF=0
   arms (modeled cancer-at-TF=0 read depth), not real healthy-donor
   plasma.**

3. **No IRB approval, no prospective collection, no clinical
   outcomes.** The validation is on open-data WGS fragments from
   public repositories.

4. **External validation exists, but it is public-data external
   validation, not clinical validation.** Publications 1/7/9 were never
   trained on and give AUC 0.9915 (pub 7, WGS), 0.8819 pooled over the
   publications carrying both classes, and 0.403 sensitivity on the
   whole-exome publication. What it still lacks: prospective follow-up,
   stage labels, and real healthy-donor controls — the non-cancer arm is
   14 public samples. A leave-one-study-out run inside the training
   cohort gives 0.6805, so the pooled 0.9670 is a within-cohort
   separation number, not a generalisation claim.

5. **Mutation-informed channel in the fusion experiment is
   synthetic.** The fusion AUC 0.9886 is a what-if pairing of a
   real fragmentomics channel with a synthetic mutation channel
   calibrated to AUC 0.92 — not a real measurement of both channels
   on the same plasma. Pairing on real plasma is the next phase.

6. **The TCGA-LUAD headline is simulation-based.** "TCGA-LUAD"
   means real TCGA tumor mutations from the GDC API; the "plasma"
   reads are Poisson-sampled at the stated tumor fraction. This is
   a dilution/spike-in benchmark, not a clinical plasma validation.

7. **Batch effects are not exhaustively characterized.** The
   harmonization check ([`docs/HARMONIZATION.md`](HARMONIZATION.md))
   is a synthetic fixture with NEUTRAL verdict; on a real
   multi-study plasma pool, larger batch effects are expected.

8. **Pretraining is on 200 samples from a single open cohort.**
   Foundation-model pretraining uses 200 samples from
   `data/finaledb_pretrain_cohort_PRODUCTION.npz`. The pre-fix
   `*SYNTHETIC_v0*` artifacts are not real-data-trained (see
   [`docs/PRETRAIN_BUG.md`](PRETRAIN_BUG.md)).

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
@software{deepcatch_v2_2_0,
  author = {Yu Ching Lam},
  title  = {DeepCatch: Mutation-informed Ultra-Sensitive cfDNA Detection
            (research-use-only open-data methods benchmark)},
  year   = {2026},
  url    = {https://github.com/rollroyces/deepcatch},
  version = {2.2.0},
  note   = {Research-use-only open-data benchmark. NOT clinically validated.
            NOT FDA-approved. NOT a medical device. For methods research only.}
}
```

## Versioning

This model card applies to **DeepCatch v2.2.0** (commit `fb1e523` on
`main`). Subsequent versions may update the headline numbers; check
the `results/` directory and the README for the latest validated
metrics.
