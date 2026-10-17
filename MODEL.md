# Model Card — DeepCatch v2.3.0

> **⚠️ STATUS: research-use-only software benchmark.**
> NOT clinically validated. NOT a medical device. NOT FDA-approved.
> NOT intended for clinical diagnosis, treatment selection, screening,
> or patient stratification.
> For methods research only.

> **v2.3.0 (2026-10-16) — CANCER-SEEK integration.** This release folds the
> CANCER-SEEK ultra-sensitive screening modules (`src/ultrasensitive/`),
> the Adaptive Screening Policy RL agent (`src/screening_policy/`), and
> the cfRNA persistence coupling (`src/methylation/`) into the DeepCatch
> umbrella. See [`docs/CANCER_SEEK_ARCHITECTURE.md`](CANCER_SEEK_ARCHITECTURE.md)
> and [`docs/ROADMAP_CANCER_SEEK.md`](ROADMAP_CANCER_SEEK.md). All previous
> model-card claims are unchanged; the modules below are reference-port
> integrations and have **not** been validated against a new cohort.
> Validate `src/ultrasensitive/spike_in_simulator.py` against the existing
> `ULTRA_EARLY_READINESS.md` calibration before publishing any LOD numbers.

## Overview

DeepCatch is an open-source pipeline for detecting cancer from cell-free
DNA (cfDNA) sequencing data. It provides **mutation-informed** detection
(variant calling + panel-LLR + strand-score + Fisher's exact) and,
through its cross-repo adapter, a **tumor-naive** detection channel
(fragmentomics: 5Mb + 100 kb short/long ratio, median-normalized
coverage, and full fragment-size-distribution histogram).

The optional **ultra-sensitive CANCER-SEEK modules** (`src/ultrasensitive/`,
requires `pip install 'deepcatch[ultra]'`) extend DeepCatch toward
0.001% VAF detection: UMI consensus, low-VAF variant calling with CHIP
filter, CNV CBS segmentation, beta-value methylation processing, XGBoost
triage + SHAP, β-VAE anomaly detection, multi-modal FusionNet, and a
spike-in LOD simulator. The optional **Adaptive Screening Policy**
(`src/screening_policy/`) ships a SAC-like RL agent over a
10,000-patient simulated cohort, calibrating screening cadence
against QALY reward. The optional **cfRNA coupling**
(`src/methylation/`) is a structural companion to `mrnavax-rare` for
companion-diagnostic read-side.

**This is research-grade software, not a clinical diagnostic tool.**
It has not been validated on a held-out clinical cohort. See
[Limitations](#limitations) below.

## Model details

| | |
|---|---|
| **Model name** | DeepCatch v2.3.0 |
| **Version** | 2.3.0 (commit on `feat/cancer-seek-integration` branch) |
| **Date** | 2026-10-16 |
| **Authors** | Yu Ching Lam (Independent Researcher, ORCID [0009-0008-9113-769X](https://orcid.org/0009-0008-9113-769X)) |
| **License** | MIT |
| **Repository** | https://github.com/rollroyces/deepcatch |
| **Companion pipeline** | [cfdna-fragmentomics-pipeline](https://github.com/rollroyces/cfdna-fragmentomics-pipeline) (the source of the tumor-naive channel) |
| **Companion data** | 627 cross-study FinaleDB samples (Jiang 2015 + Cristiano 2019) |
| **Reference paper** | `paper/PAPER.md` (242 lines) and `paper/paper.tex` (177 lines) |

## v2.3.0 — CANCER-SEEK integration (merged from sibling repos)

On 2026-10-16 the following sibling repositories were consolidated into the
DeepCatch umbrella:

| Source repo                                                          | Status                | Destination                                       |
|----------------------------------------------------------------------|-----------------------|--------------------------------------------------|
| [`rollroyces/cancer-screening`](https://github.com/rollroyces/cancer-screening) | Archived              | `src/ultrasensitive/`, `src/screening_policy/`   |
| [`rollroyces/deepcatch-methylation`](https://github.com/rollroyces/deepcatch-methylation) | Archived              | `src/methylation/`                                |
| [`rollroyces/cfDNA-Fragmentomics`](https://github.com/rollroyces/cfDNA-Fragmentomics) | Archived              | (superseded by `cfdna-fragmentomics-pipeline` + `src/fragmentomics/`) |

Activation:

```bash
# Core (always available)
pip install -e .

# Ultra-sensitive modules (XGBoost + SHAP)
pip install -e ".[ultra]"

# Adaptive Screening Policy (SAC upgrade)
pip install -e ".[ultra-rl]"
```

### Why a separate `src/ultrasensitive/` package?

`src/ultrasensitive/` is **isolated** from the DeepCatch core
(`src/fragmentomics/`, `src/foundation/`, etc.). The 189-test fast-path gate
(189 passed, 12 deselected) covers the existing DeepCatch modules and does
**not** import `src/ultrasensitive/`; the optional `[ultra]` extras gate
keeps the CANCER-SEEK modules out of the default install. Internal tests for
the merged modules are added under `test/test_ultrasensitive/` and
`test/test_methylation/`.

### Inserted modules (ports, not new designs)

| Module                                | Source path                                          | Notes                                                                                            |
|---------------------------------------|------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| `src/ultrasensitive/umi.py`           | `cancer-screening/src/core/umi.py` (SHA c3a933601)   | Bug fix applied: added `min_variants` guard to `estimate_tumor_fraction`.                        |
| `src/ultrasensitive/cnv_analyzer.py`  | `cancer-screening/src/core/cnv_analyzer.py` (SHA 7dd7d48f) | Bug fix applied: `logging.getLogger(__Logger(__name__))` -> `logging.getLogger(__name__)`.        |
| `src/ultrasensitive/variant_calling.py` | `cancer-screening/src/core/variant_calling.py` (SHA 01dae958) | Identical port.                                                                                  |
| `src/ultrasensitive/methylation.py`   | `cancer-screening/src/core/methylation.py` (SHA 4b880edb) | Identical port.                                                                                  |
| `src/ultrasensitive/spike_in_simulator.py` | `cancer-screening/src/utils/spike_in_simulator.py` (SHA 72f878c3) | Added `out_csv=...` for committed `results/` artifacts.                                          |
| `src/ultrasensitive/xgboost_classifier.py` | `cancer-screening/src/ml_models/xgboost_classifier.py` (SHA f4b13827) | Identical port; gated on `[ultra]` extras.                                                       |
| `src/ultrasensitive/methylation_encoder.py` | `cancer-screening/src/ml_models/methylation_encoder.py` (SHA ac70301) | Identical port; architecture-only. Real-data validation is Phase 4 of `METHYLATION_PROJECT.md`.    |
| `src/ultrasensitive/fusion_net.py`    | `cancer-screening/src/ml_models/fusion_net.py` (SHA 2fe61590) | Identical port of 4-modality (g+m+p+c) FusionNet. Reference; routes through `src/multimodal_fusion/` post-v2.4. |
| `src/ultrasensitive/vae_anomaly.py`   | `cancer-screening/src/ml_models/vae_anomaly.py` (SHA 697712be) | Identical port.                                                                                  |
| `src/screening_policy/rl_agent.py`    | `cancer-screening/src/screening_policy/rl_agent.py` (SHA f24ce6f3) | Identical port; SAC-like linear Q-function. Phase 5 = upgrade to NN Q head.                       |
| `src/methylation/cfrna_coupling.py`   | `deepcatch-methylation/deepcatch_methylation/cfrna_coupling.py` (v0.10.0) | Identical port.                                                                                  |
| `src/methylation/methylation_reader.py` | `deepcatch-methylation/deepcatch_methylation/methylation_reader.py` (v0.11.0) | Identical port.                                                                                  |

## Intended use

- **Primary**: research-grade benchmarking of mutation-informed ctDNA
  detection methods against simulated plasma at sub-1% VAF.
- **Secondary**: demonstration of cross-repo integration with a
  fragmentomics pipeline (pairing via file-format contract, not
  import-time dependency).
- **Tertiary**: an extensible framework for adding new ctDNA
  detection modalities (the multimodal_fusion module is the
  integration point).
- **Quaternary (v2.3.0+)**: optional CANCER-SEEK ultra-sensitive pipeline
  scaffold for UMI-based 0.001% VAF detection, CNV segmentation, and
  SHAP-aware XGBoost triage. Requires `[ultra]` extras; not part of the
  default install.

## Intended users

- ML researchers working on cfDNA fragmentomics or ctDNA detection.
- Fragmentomics / epigenomics methods developers.
- Bioinformatics scientists evaluating open-data benchmarks.
- (v2.3.0+) Researchers studying ultra-low-VAF detection or adaptive
  screening policy.

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
| Ultra-sensitive simulator (v2.3.0+) | Synthetic spike-in (no real data) | configurable VAF 1e-6 to 1e-2 | Synthetic; for LOD validation only |

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

(v2.3.0+) Ultra-sensitive modules: the ports reference synthetic
data and synthetic cohorts generated by `MultiModalDataGenerator` and
`SpikeInSimulator`. The CANCER-SEEK `UMI` and `variant_calling`
modules target the GDC + FinaleDB ingestion path; real-data
calibration is a Phase 5 follow-up.

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

### Ultra-sensitive CANCER-SEEK (v2.3.0+) — new

(v2.3.0+) **No new validation numbers.** The CANCER-SEEK ports are
**reference implementations only** at this version. They include:

- `src/ultrasensitive/spike_in_simulator.py` — synthetic LOD validation
  across VAF 1e-6 to 1e-2. **Required** to validate before publishing
  any 0.001% VAF headline number; binds to `ULTRA_EARLY_READINESS.md`.
- `src/ultrasensitive/xgboost_classifier.py` — XGBoost + SHAP baseline
  on synthetic multi-omics features. **Not yet validated against
  TCGA-Frozen or FinaleDB cohort** — required before any claim.
- `src/methylation/cfrna_coupling.py` — read-side companion to
  `mrnavax-rare`. v0.10.0 thresholds are **research-use defaults**;
  no patient data processed.
- `src/screening_policy/rl_agent.py` — SAC-like RL agent on a
  10,000-patient simulated cohort. No clinical screening data;
  detect-rate / FPR are honest simulation outputs only.

**Honest disclosure**: the CANCER-SEEK ULOD (ultra-low detection)
modules do not yet have a real-plasma or real-cohort benchmark; treat
any future number with the same "synthetic fixture" caveat as the
existing `Foundation` smoke from Audit-2.

## Limitations

1. **Cohort sizes are small.**
   - n=20 paired synthetic for the foundation smoke (TCGA-LUAD).
   - n=627 for the open-data cross-study benchmark (FinaleDB pan-cancer).
   - n=129 processed frequency vectors for the Jiang motif benchmark.
   These are appropriate for methods research but are far below the
   cohort sizes typical of clinical validation (1,000+).
   - (v2.3.0+) n=0 for the CANCER-SEEK ultra-sensitive modules.

2. **Healthy controls in the cross-study benchmark are FinaleDB TF=0
   arms (modeled cancer-at-TF=0 read depth), not real healthy-donor
   plasma.**

3. **No IRB approval, no prospective collection, no clinical
   outcomes.** The validation is on open-data WGS fragments from
   public repositories.

4. **No external held-out validation.** Every headline number is
   pooled OOF on the same cohort the model was trained on. The
   cross-study benchmark pools two studies with per-study
   z-score harmonization but does not hold one study out.

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

11. **(v2.3.0+) CANCER-SEEK modules are reference ports.** All ultra-
    sensitive modules require `[ultra]` extras and have **not been
    benchmarked against real plasma**. Treat any LOD number as a
    synthetic-fixture simulation output.

12. **(v2.3.0+) The Adaptive Screening Policy uses a simplified
    linear Q-function** — not a full neural-network SAC. Real-world
    deployment requires a Phase 5 upgrade.

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

# Optional: enable CANCER-SEEK ultra-sensitive modules
pip install -e ".[ultra]"

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

# CANCER-SEEK LOD validation (with [ultra] extras)
deepcatch-lod --vaf_range "[1e-2, 1e-3, 1e-4, 1e-5]" --n_samples 20
```

## Citation

```bibtex
@software{deepcatch_v2_3_0,
  author = {Yu Ching Lam},
  title  = {DeepCatch: Mutation-informed Ultra-Sensitive cfDNA Detection
            (research-use-only open-data methods benchmark,
            with CANCER-SEEK integration)},
  year   = {2026},
  url    = {https://github.com/rollroyces/deepcatch},
  version = {2.3.0},
  note   = {Research-use-only open-data benchmark. NOT clinically validated.
            NOT FDA-approved. NOT a medical device. For methods research only.
            CANCER-SEEK ultra-sensitive modules are reference ports; no
            real-cohort validation yet.}
}
```

## Versioning

This model card applies to **DeepCatch v2.3.0** (on the
`feat/cancer-seek-integration` branch at merge). Subsequent versions
may update the headline numbers; check the `results/` directory and
the README for the latest validated metrics.
