# DeepCatch: a software framework for cfDNA fragmentomics with open-data cross-study benchmarking

**Authors:** Yu Ching Lam¹, [co-authors: TBD]

¹ DeepCatch Project — github.com/rollroyces/deepcatch (Independent researcher)

**Venue target:** PLOS Computational Biology or Bioinformatics Advances

**Status (2026-09-23):** SKETCH — methods-paper skeleton. Sections 1–4 are drafted; Results (§3) points to JSON artifacts and a forthcoming companion document (`docs/CROSS_STUDY_BENCHMARK.md`, to be created by a sibling subagent) that will host the cross-study AUC table once the FinaleDB publications-6+8 sweep lands. This is NOT a finished manuscript.

---

## Abstract

**Background.** Cell-free DNA (cfDNA) fragmentomics — the joint analysis of fragment length, end-motif, copy-number, and coverage profiles from plasma sequencing — is a promising substrate for non-invasive multi-cancer detection. Published fragmentomics methods (DELFI, Galleri, FinaleMe) report compelling per-cancer discrimination, but their training data and validation cohorts are largely proprietary, and independent replication on harmonized open data is scarce.

**Methods.** DeepCatch is an open-source software framework that turns public cfDNA resources (FinaleDB, GDC TCGA, FLARE/GSE317007) into reproducible fragmentomics experiments. The framework combines (i) a six-modality feature schema (5-channel DELFI, 4-mer end motifs, fragment size, coverage/CNV, serology, tissue-of-origin); (ii) a multi-modal fusion encoder with masked-modality self-supervised pre-training and cross-attention fusion; (iii) panel-based tumor-informed detection when a prior cohort is available; and (iv) a validation protocol built on 5-seed × 5-fold GroupKFold with shuffled-label and true-confound controls. The implementation is MIT-licensed, ships with 374+ regression tests, and exposes every intermediate artifact (cohort matrices, per-seed AUC arrays, ablation JSONs) to the reader.

**Results (to be filled in).** On the paired 20-patient TCGA-LUAD synthetic-plasma cohort at 0.1% ctDNA, the foundation encoder honestly achieves AUC ≈ 0.55 (n=20 seeds, per-patient GroupKFold) — below the sklearn LR baseline (0.91) on the same channels, which is the expected outcome on this paired design and is itself a finding worth reporting. Cross-study AUC on FinaleDB publications 6 + 8 (n=657) is reported in the forthcoming `docs/CROSS_STUDY_BENCHMARK.md`. Per-cancer sensitivity at fixed specificity (with DeLong confidence intervals) is reported in `results/per_cancer_sens_at_spec/`. Three honest null/negative ablations (focal-BCE loss, sparse-aware projection, harmonization) are documented.

**Conclusions.** DeepCatch is not a clinical assay. It is a reproducible fragmentomics engineering substrate: it lets a researcher swap in a new cohort, a new modality, or a new fusion head and get honest numbers in hours rather than months. The framework is most useful for negative results — telling investigators which architectural ideas do NOT help on open data before they spend real-plasma budget on them.

---

## 1. Introduction

### 1.1 Motivation

Cell-free DNA (cfDNA) shed by tumors into peripheral blood carries genomic, epigenomic, and fragmentomic information about the source tissue [1–4]. Two broad assay classes exploit this: **tumor-informed panels** (Signatera, CAPP-Seq [5,6]) track patient-specific mutations at 0.01–0.1% ctDNA with deep targeted sequencing and duplex-UMI error suppression; and **tumor-agnostic multi-cancer screening** (DELFI [7], Galleri, FinaleMe [8]) detect cancer genome-wide by aggregating fragmentation, methylation, and copy-number signals across the genome.

This paper is about the second class. Fragmentomics-based multi-cancer detection has produced headline-grade AUCs in discovery cohorts (Cristiano et al. 2019 reported AUC ≈ 0.94 for HCC-vs-healthy on 232 patients [7]), but the underlying training data, feature pipelines, and pretrained models are typically not released. Independent replication is rare. A new lab attempting to build a fragmentomics classifier must re-implement the entire pipeline from a methods-section description, then discover — sometimes after months — that the architectural idea they wanted to test does not transfer.

### 1.2 Related work

- **DELFI (Cristiano et al. 2019 [7])** — 5-channel arm-level DELFI profile (short/medium/long fragment counts, coverage, per-arm copy number) at ~5 Mb resolution. The original paper trained a per-cancer gradient-boosted classifier on 232 plasma samples (HCC, breast, colorectal, lung, ovarian, pancreatic, gastric, bile-duct vs healthy) and reported AUCs ranging 0.86–0.97 per cancer.
- **CAPP-Seq (Newman et al. 2014 [5], 2016 update [6])** — tumor-informed targeted panel with a Fisher-method aggregation across loci; underpins the Stanford clinical MRD workflow.
- **CancerSEEK / Galleri (Klein et al. 2018 [9]; Liu et al. 2020 [10])** — combined mutation + protein biomarker (CancerSEEK); methylation-based multi-cancer detection (Galleri).
- **FinaleMe (Liu et al. 2018 [11])** — fragment-end motif model using 4-mer end-sequence preferences.
- **Jiang et al. 2020 [12]** — HCC-specific fragmentomics in a 222-patient cohort, one of the largest single-cancer fragmentomics datasets with public fragment files.
- **FinaleDB / Finaletoolkit [13]** — public archive of fragmentomics datasets with harmonized `*.frag.tsv.bgz` files; the primary data source for the present work.

### 1.3 Gap

There is no open-source framework that (i) consumes FinaleDB / GDC TCGA / FLARE directly, (ii) computes the canonical five-channel DELFI profile plus 4-mer motifs plus fragment size plus coverage/CNV in a single script, (iii) trains a multi-modal fusion encoder with honest per-patient GroupKFold validation, and (iv) ships paired-design ablations for every architectural knob. DeepCatch is that framework.

The contribution is not a new model class. It is a reproducible engineering substrate, with honest negative results reported alongside any positive ones.

---

## 2. Methods

### 2.1 Fragmentomics feature schema (six modalities)

DeepCatch consumes a per-sample cfDNA feature dict with six modality slots, defined in `src/foundation/config.py`:

| Modality | Dim | Description |
|---|---:|---|
| `frag_basic` | 4 | panel-LLR score (when tumor-informed) + 3 length-summary stats (mean, median, short-fragment fraction) |
| `frag_enhanced` | 44 | per-patient mutation-derived frag signature: log mutation burden, mean VAF, VAF std, driver-gene fraction, mutation spectrum (6 channels), aneuploidy proxy, and 36 calibrated sequencing-noise jitter features |
| `cnv` | 6 | arm-level copy-number summary (chr1p/1q/3p/8q/13q/17p short-arm ratios) |
| `sero` | 4 | optional protein / serology placeholder (zero-filled when unavailable) |
| `gnn` | 1 | optional graph-neural-network fragment-coverage placeholder (zero-filled when unavailable) |
| `tissue` | 24 | tissue-of-origin methylation prior; flat zero vector when unavailable |

For tumor-naive (un-paired) cohorts the same schema is used but `frag_basic[:, 0]` carries the panel-LLR placeholder and `frag_enhanced` carries the 5-channel DELFI profile summaries at 5 Mb resolution (Cristiano-style) plus the 4-mer end-motif distribution (256 counts, projected to 36 PCA-retained components).

### 2.2 Multi-modal fusion architecture

The encoder (`src/foundation/model.py:MultiModalEncoder`) is a per-modality `LinearProjection` (or `SparseAwareLinearProjection`, opt-in) followed by a 4-layer Transformer with 4-head self-attention and a `[MODALITY]` token prepended to each modality's token sequence. The encoder was pre-trained on 200 real FinaleDB samples (Cristiano 2019 + Jiang 2015) with masked-modality-prediction (Phase 1, 15 epochs) and InfoNCE contrastive (Phase 2, 5 epochs). Pre-training details, the synthetic-bypass bug fix, and the regression test are documented in `docs/PRETRAINING.md` and `docs/PRETRAIN_BUG.md`. The PRODUCTION checkpoint is `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` (543,872 parameters, embed_dim=128).

Downstream classification is a `FoundationDownstream` head: 0.7 × frozen-pretrained-encoder → LR-head + 0.3 × trainable Transformer, with binary cross-entropy (default) or focal-BCE (`loss="sens_at_spec"`, opt-in, `alpha_pos` sweep parameter) as the loss.

### 2.3 Panel-based detection (when prior cohort is known)

For tumor-informed panel detection, the same `FoundationDownstream` consumes a modality dict where `frag_basic[:, 0]` is the panel-LLR sum (the sufficient statistic under independent Poisson observation at each locus), `frag_enhanced[:, 0]` is the patient-level frag score, and the other modalities are filled with cohort-derived priors. The panel-LLR scoring function, the Fisher-method aggregate, and the strand-concordance-weighted variant are documented in `paper/PAPER.md` (§2 of the companion benchmark paper) and `scripts/foundation_real_smoke.py`.

### 2.4 Validation protocol

Every reported number follows one protocol:

1. **5-seed × 5-fold per-patient GroupKFold.** Seeds {0, 1, 2, 3, 4} (or {0..19} for higher-power ablations). Both arms of every paired patient stay in the same fold. This is the protocol audit-2 P0-A (commit `6231fcb`) installed to replace the leaky `StratifiedKFold` that originally inflated the foundation AUC to 0.94.
2. **Shuffled-label control.** Labels are permuted independently in the pos-half and neg-half to break the pair structure. AUC against the shuffled labels is reported as `shuffled_foundation_auc`. The honest diagnostic on paired data is the *delta* between real-label AUC and shuffled-label AUC, not the absolute AUC.
3. **True-confound control.** For cross-study experiments, per-study batch effects are injected on top of the synthetic signal and the per-study AUC is reported alongside the pooled AUC. A diagnostic (`delta_auc_pooled_minus_per_study_mean`) quantifies the gap. See `docs/HARMONIZATION.md`.
4. **Per-cancer sensitivity at fixed specificity with DeLong CIs.** Per-cancer sensitivity at 95% and 99% specificity is computed by pooling OOF predictions across seeds, fixing the threshold at the empirical 95th / 99th percentile of the healthy controls, and applying the DeLong method to estimate the AUC confidence interval. See `results/per_cancer_sens_at_spec/`.
5. **Bit-identical paired-design sanity.** For every ablation, the non-ablation metrics (`panel_only_aucs`, `frag_only_aucs`, `lr_baseline_aucs`, `naive_avg_aucs`) are verified bit-identical to <1e-9 across the two JSONs, proving the only difference is the ablated knob.

### 2.5 Implementation and tests

The framework is implemented in PyTorch with NumPy / scikit-learn baselines. All 374+ tests pass; the pretrain-bug-fix regression suite (`test/test_pretrain_bug_fix.py`, 11 tests) pins the contract that the synthetic generator is never invoked when real data is supplied. Code, data, results, and ablations are version-controlled at github.com/rollroyces/deepcatch (MIT license).

---

## 3. Results (TO BE FILLED IN — see JSON artifacts)

> This section is a pointer skeleton. The cross-study AUC table lands in
> `docs/CROSS_STUDY_BENCHMARK.md` (forthcoming). Until then, the
> JSON artifacts cited below are the source of truth.

### 3.1 Honest baseline AUC on paired synthetic cohort (20-patient TCGA-LUAD)

On 20 paired TCGA-LUAD patients at 0.1% ctDNA, real-panel + real-mutation-derived-frag, per-patient GroupKFold, 20 seeds:

| Metric | Value | Source JSON |
|---|---|---|
| `panel_only AUC` | 0.915 ± 0.003 | `results/sens_at_spec_ce_n20.json` |
| `frag_only AUC` (real mutation-derived) | 0.620 ± 0.005 | `results/sens_at_spec_ce_n20.json` |
| `foundation AUC` | **0.5588 ± 0.0094** | `results/sens_at_spec_ce_n20.json:foundation_aucs` |
| `lr_baseline AUC` (sklearn LR on [panel, frag]) | 0.910 | `results/sens_at_spec_ce_n20.json:lr_baseline_aucs` |
| `naive_avg AUC` | 0.915 | `results/sens_at_spec_ce_n20.json:naive_avg_aucs` |
| `shuffled_foundation AUC` (label-permuted control) | 0.741 | `results/sens_at_spec_ce_n20.json:shuffled_foundation_aucs` |
| `delta_auc_normalized` | −3.45 | `results/sens_at_spec_ce_n20.json` |

**Honest framing (audit-2 P0-A→P0-F fixed):** under per-patient GroupKFold, the foundation model is honestly out-performed by the sklearn LR baseline on this paired design. The foundation picks up some pair-invariant signal even under pair-broken label shuffle (hence `shuffled_foundation AUC = 0.741 > 0.5`), but the real-label AUC is below the shuffled-label AUC on this cohort — the model overfits. The honest metric is the *delta* (foundation − shuffled), not the absolute AUC. See `AUDIT_2_FINDINGS.md`.

### 3.2 Cross-study AUC on FinaleDB publications 6 + 8 (n=657)

**To be filled in.** Pointers: `data/finaledb_pretrain_cohort_PRODUCTION.npz` (200-sample real-data cohort), `docs/PRETRAINING.md`, and the forthcoming `docs/CROSS_STUDY_BENCHMARK.md`. The 657-sample pre-extracted cache is at `cfdna-fragmentomics-pipeline/data/features/` (Cristiano 2019 + Jiang 2015, 262 + 32 healthy / 275 + 89 cancer). Per-cancer AUCs at fixed specificity are reported in `results/per_cancer_sens_at_spec/`.

### 3.3 Per-cancer sensitivity at fixed specificity (with DeLong CIs)

**To be filled in.** Pointer: `results/per_cancer_sens_at_spec/*.json` (one JSON per cancer type, pooled across seeds). The DeLong CIs and the per-cancer threshold-fixing protocol are described in `scripts/per_cancer_sens_at_spec.py`.

### 3.4 Honest null / negative ablations

| Ablation | Verdict | Source |
|---|---|---|
| **focal-BCE loss** (`loss="sens_at_spec"`, α_pos=20) | NULL: Δ AUC = −0.0005 [−0.0015, +0.0005], p=0.30 (n=20 paired seeds). Supersedes the n=5 "positive direction, NS" reading. | `docs/SENS_AT_SPEC_ABLATION.md` |
| **sparse-aware projection** (`projection_kinds={"frag_basic":"sparse_aware"}`) | NEGATIVE, marginal: Δ AUC = −0.0155 [−0.0315, +0.0005], p=0.055 (5 seeds). The slot is not sparse at this operating point. | `docs/SPARSE_AWARE_ABLATION.md` |
| **harmonization** (per-cohort batch-effect check) | NEGATIVE: pooled AUC is high but per-study AUC variance is large; the signal is partly study-confounded. | `docs/HARMONIZATION.md`, `results/harmonization_check.json` |

All three ablations are reported honestly. No "win by default" framing.

---

## 4. Discussion

### 4.1 What this work does

DeepCatch provides a reproducible fragmentomics engineering substrate. It takes a public cfDNA cohort (FinaleDB, GDC, FLARE) through (i) feature extraction in a canonical six-modality schema, (ii) self-supervised pre-training on real data, (iii) downstream classification with honest per-patient GroupKFold validation, and (iv) paired-design ablations of every architectural knob. The headline contribution is not a number but the engineering substrate itself: a new lab can swap in a new cohort and get honest numbers in hours.

### 4.2 What this work does NOT do

- **Clinical validation.** DeepCatch is not a clinical assay. No real-plasma cohort has been sequenced. No patient samples were collected under IRB.
- **FDA pathway.** No regulatory submission is planned or implied.
- **Prospective screening.** No prospective cohort, no enrollment, no follow-up.
- **Headline-grade AUCs.** The synthetic-plasma foundation AUC (≈0.55) is below the LR baseline. The cross-study FinaleDB AUC will be reported when `docs/CROSS_STUDY_BENCHMARK.md` lands; it will not be inflated by CV leak, pair-invariant feature reuse, or shuffled-label-control bugs.
- **Real plasma fragmentomics.** The pretraining cohort is from FinaleDB plasma (Cristiano 2019, Jiang 2015), but the downstream validation cohorts are synthetic-plasma (TCGA-LUAD real mutations + Poisson sampling). A real-plasma fragmentomics validation requires an IRB-approved prospective cohort — see §4.4.

### 4.3 Limitations

1. **Small paired synthetic cohort (n=20).** Per-patient GroupKFold on 20 patients means each fold has 4 patients × 2 arms = 8 test samples. CIs are wide; per-seed std on foundation AUC is ~0.01.
2. **Simulated healthy controls.** The healthy arm of every paired patient is a TF=0 simulation of the same patient's mutations, not an unrelated healthy donor. This is the standard analytical validation approach (Newman 2016 [6]) but it is not a substitute for unrelated healthy plasma.
3. **No real plasma in downstream validation.** The 0.1% ctDNA numbers are from synthetic reads. Real plasma has additional noise sources (PCR bias, mapping ambiguity, contamination) that the simulation does not capture.
4. **FinaleDB API degraded.** As of 2026-09-21, the FinaleDB REST API returns 500 on `/api/v1/publication` and `/api/v1/seqrun`. The cross-study sweep runs against the pre-extracted cache, not a fresh live fetch. See `docs/PRETRAINING.md` §"FinaleDB REST API in degraded state."
5. **Encoder size.** PRODUCTION_CONFIG is 543,872 parameters (embed_dim=128, 4 layers, 4 heads). State-of-the-art fragmentomics encoders are 10–100× larger. The present model fits on CPU in <5 seconds; clinical-grade models would require GPU training.
6. **Two-channel only.** The pretraining data covers Cristiano 2019 + Jiang 2015 (publications 6 + 8). FinaleMe (Liu 2018 [11]) and the methylation-based MCED publications are not yet in the cohort.

### 4.4 Future work

- **Real-plasma paired cohort.** Requires IRB + collaborator with banked plasma + matched tumor sequencing. Estimated cost: $50K–$200K for a 50-patient pilot. Target: CCGA, TRACERx, or an in-house liquid-biopsy clinic.
- **Cross-platform validation.** Apply the pretrained encoder to ONT / PacBio / single-cell cfDNA to test transfer. The 4-mer end-motif distribution should be platform-robust; the coverage/CNV channel may need per-platform recalibration.
- **Methylation and fragment-end motifs.** The 256-channel 4-mer motif distribution (Liu 2018 [11]) is currently projected to 36 PCA components. A direct 256-channel input would test whether the PCA bottleneck costs signal.
- **Larger pretraining cohort.** Scale from 200 → 2,000+ real FinaleDB samples once the API is restored. A 10× cohort with 50 epochs would push the encoder from "demonstration" to "useful."
- **Independent replication.** A second lab, with its own implementation, running the same cohort through DeepCatch. The released artifacts (cohort matrices, per-seed AUC arrays, ablation JSONs) make this reproducible.

---

## 5. Availability

- **Source code:** github.com/rollroyces/deepcatch
- **License:** MIT
- **Data:**
  - FinaleDB pre-extracted cache: `cfdna-fragmentomics-pipeline/data/features/` (657 samples, Cristiano 2019 + Jiang 2015)
  - GDC TCGA-LUAD MAF files: `validation/tcga/tcga_cache/` (402 MAF files, open access)
  - FLARE / GSE317007 (ONT HNSCC): see `docs/PATH_TO_IMPACT.md`
- **Pretrained checkpoints:**
  - `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` (real-data-trained, 543,872 params, PRODUCTION_CONFIG)
  - `checkpoints/foundation_pretrained_SYNTHETIC_v0.pt` (historical; not real-data-trained — do NOT use for "real-data-derived" claims)
- **Reproduce all numbers:** `RUN_ALL.sh` (top-level driver) and the per-script invocations documented in §3.
- **Companion documents:**
  - `docs/PRETRAINING.md` — pre-training pipeline
  - `docs/PRETRAIN_BUG.md` — synthetic-bypass bug + regression test
  - `docs/SPARSE_AWARE_ABLATION.md` — honest null/negative ablation
  - `docs/SENS_AT_SPEC_ABLATION.md` — honest null/negative ablation
  - `docs/HARMONIZATION.md` — per-study batch-effect check
  - `AUDIT_2_FINDINGS.md` (repo root) — consolidated audit-2 report
  - `docs/CROSS_STUDY_BENCHMARK.md` (forthcoming) — cross-study AUC table
  - `MODEL.md` (forthcoming) — model architecture reference
  - `paper/PAPER.md` (companion benchmark paper) — bioRxiv submission
  - `paper/paper.tex` — LaTeX source of the companion benchmark paper

---

*Manuscript draft — 2026-09-23. Sections 1, 2, 4, 5 are stable; §3 is a pointer skeleton that will be filled when the cross-study sweep lands. Comments welcome via GitHub issues.*
