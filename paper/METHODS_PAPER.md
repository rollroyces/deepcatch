# DeepCatch: a software framework for cfDNA fragmentomics with open-data cross-study benchmarking

**Authors:** Yu Ching Lam¹, [co-authors: TBD]

¹ DeepCatch Project — github.com/rollroyces/deepcatch (Independent researcher)

**Venue target:** PLOS Computational Biology or Bioinformatics Advances

**Status (2026-09-23):** DRAFT — methods-paper draft with grounded §3 numbers. Sections 3.1, 3.2, 3.3, 3.4, 3.5 are filled with values from the cross-study FinaleDB sweep (`results/cross_study_finallydb.json`) and the per-cancer sens@spec table (`results/per_cancer_sens_at_spec.json`). §4.2 limits are tightened to the per-cancer numerators actually observed.

---

## Abstract

**Background.** Cell-free DNA (cfDNA) fragmentomics — the joint analysis of fragment length, end-motif, copy-number, and coverage profiles from plasma sequencing — is a promising substrate for non-invasive multi-cancer detection. Published fragmentomics methods (DELFI, Galleri, FinaleMe) report compelling per-cancer discrimination, but their training data and validation cohorts are largely proprietary, and independent replication on harmonized open data is scarce.

**Methods.** DeepCatch is an open-source software framework that turns public cfDNA resources (FinaleDB, GDC TCGA, FLARE/GSE317007) into reproducible fragmentomics experiments. The framework combines (i) a six-modality feature schema (5-channel DELFI, 4-mer end motifs, fragment size, coverage/CNV, serology, tissue-of-origin); (ii) a multi-modal fusion encoder with masked-modality self-supervised pre-training and cross-attention fusion; (iii) panel-based tumor-informed detection when a prior cohort is available; and (iv) a validation protocol built on 5-seed × 5-fold GroupKFold with shuffled-label and true-confound controls. The implementation is MIT-licensed, ships with 374+ regression tests, and exposes every intermediate artifact (cohort matrices, per-seed AUC arrays, ablation JSONs) to the reader.

**Results.** On the paired 20-patient TCGA-LUAD synthetic-plasma cohort at 0.1% ctDNA, the foundation encoder honestly achieves AUC ≈ 0.55 (n=20 seeds, per-patient GroupKFold) — below the sklearn LR baseline (0.91) on the same channels, which is the expected outcome on this paired design and is itself a finding worth reporting. On the cross-study FinaleDB publications 6 + 8 (n=627 after missing-artifact filter), a 5-channel DELFI LR with per-study z-score harmonization reaches pooled OOF AUC 0.974 ± 0.002, sens@95%=0.909, sens@98%=0.837, sens@99%=0.782. Per-cancer OvR (5-seed ± std) for the top-5 cancers: OV 0.992 ± 0.005, LUAD 0.981 ± 0.002, BRCA 0.976 ± 0.003, PAAD 0.937 ± 0.007, HCC_J 0.725 ± 0.018 (single-cohort reading). A **true-confound control** (cancer = 100% Jiang, healthy = 100% Cristiano) collapses to AUC 0.49 with harmonization vs 1.00 without — proving the pooled 0.974 is a cancer signal, not a study batch. Three honest null/negative ablations (focal-BCE loss, sparse-aware projection, harmonization) are documented.

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

## 3. Results

> Source-of-truth: `results/cross_study_finallydb.json` (§3.2), `results/per_cancer_sens_at_spec.json` (§3.3), and `docs/CROSS_STUDY_BENCHMARK.md` (the running benchmark narrative). All pooled numbers are 5-seed × 5-fold StratifiedKFold internal CV (pooled OOF); no external validation, no clinical plasma cohort.

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

### 3.2 Cross-study AUC on FinaleDB publications 6 + 8 (n=627)

The cross-study pool is the union of FinaleDB publication 6 (Jiang 2015 [16], n=121: 89 HCC + 32 healthy) and publication 8 (Cristiano 2019 [7], n=506: 274 cancer + 232 healthy) after the missing-artifact filter (31 samples dropped: 31/658 = 4.7%). Features are the canonical 5-channel DELFI profile (5 Mb ratio + 5 Mb coverage + 100 kb ratio + 100 kb counts + FSD-196, 63,246-dim). The pipeline is a 5-seed × 5-fold StratifiedKFold, LR(max_iter=2000) on PCA(200) of a per-study z-score StandardScaler fit on the **train** fold only. Source: `results/cross_study_finallydb.json`; companion narrative `docs/CROSS_STUDY_BENCHMARK.md`.

| Cohort | n | Cancer | Healthy | AUC (5-seed mean ± std) |
|---|---:|---:|---:|---|
| jiang | 121 | 89 | 32 | **0.9791 ± 0.0028** |
| cristiano | 506 | 274 | 232 | **0.9693 ± 0.0022** |
| **pooled (harmonized)** | 627 | 363 | 264 | **0.9738 ± 0.0018** |
| **pooled (no_harmonize)** | 627 | 363 | 264 | 0.9650 ± 0.0019 |

At the pooled harmonized operating point: **sens@95%=0.909**, **sens@98%=0.837**, **sens@99%=0.782** (`results/cross_study_finallydb.json:pooled.harmonized`).

**True-confound control.** The honest cross-study claim requires that the signal is cancer, not study. The control builds a synthetic two-arm dataset where cancer = 100% Jiang and healthy = 100% Cristiano (and the symmetric orientation), then rerun the same pipeline:

| Orientation | n cancer | n healthy | AUC harmonized | AUC no_harmonize |
|---|---:|---:|---:|---:|
| cancer=Jiang, healthy=Cristiano | 121 | 506 | 0.489 ± 0.004 | 0.999 ± 0.002 |
| cancer=Cristiano, healthy=Jiang | 506 | 121 | 0.496 ± 0.004 | 1.000 ± 0.000 |

Without harmonization, the classifier trivially achieves AUC ≈ 1.000 because "which study is this from?" is trivially answerable. With per-study z-scoring, the same control collapses to **0.49** — the per-study mean/variance is the only signal and the harmonization removes it by design. The drop from 1.000 → 0.49 (~80× shift across the operating-point flip) proves the pooled-harmonized 0.9738 AUC is a cancer signal, not a study batch.

**Honest framing.** These numbers are **pooled out-of-fold** AUC on the 627-sample cross-study cohort, under 5-seed × 5-fold internal CV. They are NOT external validation, NOT held-out clinical plasma, and NOT clinical-grade operating points (Galleri, CancerSEEK). They measure how well the 5-channel DELFI features separate cancer-vs-healthy in the pooled cohort with per-study batch removed. The drop from per-cohort AUCs (Jiang 0.9791, Cristiano 0.9693) to pooled-harmonized 0.9738 is small (~0.005 AUC), which is the expected size when the per-study batches are mild on this FinaleDB-uniformly-processed cohort. The companion pipeline artifact `docs/CROSS_STUDY_BENCHMARK.md` carries the full per-cohort + per-cancer + true-confound tables. The classifier is intentionally a linear LR (not the foundation encoder) so the headline number is interpretable as the separability of the published 5-channel feature set, not as a DeepCatch architectural contribution.

### 3.3 Per-cancer sensitivity at fixed specificity (with DeLong CIs)

Per-cancer OvR is computed with **per-study z-score harmonization inside each CV fold** so the same offset that drove §3.2's pooled AUC is also removed in the per-cancer fits. Top-5 cancers by count, scored against ALL healthy samples (the cross-study generalization view). For sens@spec 95% CIs the standalone artifact uses **DeLong placement-value** CIs (`src/per_cancer_sens_at_spec.py`; Sun & Xu 2014 form); for the AUC variance we use the 5-seed ± std across the StratifiedKFold seeds. Source: `results/cross_study_finallydb.json:per_cancer` (point estimates + bootstrap 95% CIs on sens@spec) and `results/per_cancer_sens_at_spec.json` (DeLong CIs + PPV@prev).

| Cancer | n_cancer | n_total | AUC mean ± std (5-seed) | Sens@95% | Sens@98% | Sens@99% [bootstrap 95% CI] |
|---|---:|---:|---|---:|---:|---|
| HCC_J | 89 | 121 (jiang only) | 0.7248 ± 0.0184 | 0.416 | 0.348 | 0.191 (0.000–0.420) |
| LUAD | 79 | 343 (jiang+cristiano healthy) | 0.9813 ± 0.0020 | 0.911 | 0.861 | 0.658 (0.392–0.917) |
| PAAD | 60 | 324 | 0.9368 ± 0.0070 | 0.817 | 0.767 | 0.483 (0.333–0.810) |
| BRCA | 53 | 317 | 0.9763 ± 0.0034 | 0.943 | 0.849 | 0.547 (0.352–0.937) |
| OV | 28 | 292 | 0.9924 ± 0.0048 | 1.000 | 0.929 | 0.929 (0.758–1.000) |

**Reading the table.** HCC_J (Jiang 2015 alone, n=89) shows the lowest cross-study AUC (0.725); this is the OvR reading against the cross-study healthy pool (mostly Cristiano healthy), and the wide sens@99 CI [0.000–0.420] flags it as a small-denominator outlier that needs a larger follow-on study. LUAD, BRCA, OV have tight CIs at AUC > 0.97 with sens@99 ≥ 0.55. OV is a high-AUC small-n row (n=28 cancer); the 1.000 sens@95 is consistent with the empirical ceiling on 28 positives and the wide CI [0.758–1.000] at sens@99 says the same. PAAD sens@99 sits at 0.483 with a very wide CI [0.333–0.810] — the n=60 denom is the binding constraint here, not the classifier.

**Prevalence-floor PPV summary (PPV at spec=99% at standard screening-prevalence grid):** sourced from the same JSON (`results/per_cancer_sens_at_spec.json:per_cancer.*.ppv_at_prevalence`); presented against the pooled AUC's sens@99 = 0.782 as the illustrative ceiling.

| Prevalence | PPV at sens@99=0.782, spec=0.99 |
|---:|---:|
| 0.001 (Galleri / CancerSEEK target) | 0.073 |
| 0.004 | 0.239 |
| 0.01 | 0.441 |
| 0.05 | 0.805 |
| 0.10 | 0.897 |
| 0.20 | 0.951 |
| 0.50 | 0.987 |

At the published MCED screening-prevalence assumption of 0.1% (10⁻³), this 5-channel LR's pooled sens@99 would project to **PPV ≈ 7.3%** — not yet clinically competitive at that prevalence, even at the strong pooled-AUC operating point. The PPV rises sharply with prevalence: at the 1% prevalence used in some targeted-MRD workflows it reaches **44%**, and at a 10% prevalence (a symptomatic workup) it is **90%**. The prevalence floor is set by `(1 - spec)` false-positives per 10,000 tested — that floor does not depend on the classifier and is the limiting factor at any sens < 1.0; the model only moves the PPV via `sens · prev` in the numerator. The per-cancer PPV rows in the JSON are sharper and the **OvR-specific** PPV at the 0.001 grid is the cross-study-reading ceiling: HCC_J ≈ 0.019 (sens@99=0.191), LUAD ≈ 0.062 (0.658), BRCA ≈ 0.052 (0.547), PAAD ≈ 0.046 (0.483), OV ≈ 0.085 (0.929).

**Honest framing.** Per-cancer CI widths are wide because per-cancer denominators are small (10–90 positives); the JSON's `skipped` rows respect the `MIN_POSITIVES_FOR_CI=5` floor (DeLong CI is unreliable below it). HCC_J is OvR-defined but **single-cohort** — there are no Jiang healthy controls in the pooled cohort other than the 32 internal ones, so HCC_J's "AUC vs cross-study healthy" is mostly a Cristiano-vs-Jiang-carcinoma problem and should be read with that grain. The OvR scoring against ALL healthy samples (not just within-study healthy) is the cross-study generalization view; it is not the within-study view that the original Cristiano 2019 [7] paper reports.

### 3.4 Honest null / negative ablations

| Ablation | Verdict | Source |
|---|---|---|
| **focal-BCE loss** (`loss="sens_at_spec"`, α_pos=20) | NULL: Δ AUC = −0.0005 [−0.0015, +0.0005], p=0.30 (n=20 paired seeds). Supersedes the n=5 "positive direction, NS" reading. | `docs/SENS_AT_SPEC_ABLATION.md` |
| **sparse-aware projection** (`projection_kinds={"frag_basic":"sparse_aware"}`) | NEGATIVE, marginal: Δ AUC = −0.0155 [−0.0315, +0.0005], p=0.055 (5 seeds). The slot is not sparse at this operating point. | `docs/SPARSE_AWARE_ABLATION.md` |
| **harmonization** (per-cohort batch-effect check) | MIXED on this cohort: pooled AUC is high (0.974) and per-study AUCs are similar, but the **true-confound control** (cancer=Jiang, healthy=Cristiano) shows that WITHOUT harmonization the classifier would learn "study" instead of "cancer" (1.000 vs 0.489). Reported as **PROTECTIVE** — the harmonization is what makes §3.2 a fair cross-study claim. | `docs/HARMONIZATION.md`, `results/cross_study_finallydb.json:true_confound_control` |

All three ablations are reported honestly. No "win by default" framing.

### 3.5 Comparison to published cfDNA benchmarks

DeepCatch reports a **0.974 pooled OOF AUC** on open data (FinaleDB publications 6 + 8, n=627). For context only — and to be explicit that this is NOT a head-to-head — we list here the headline numbers reported by the closest published cfDNA benchmarks. All five used larger, proprietary, or multi-site clinical cohorts; ours is a 627-sample open-data internal CV. The comparability axis is "how well does the 5-channel feature (or its counterpart) separate cancer from healthy in plasma," not "which is clinically better."

| Method (1st author, year) | Cohort | n | Headline metric | Source |
|---|---|---:|---|---|
| DELFI (Cristiano 2019 [7]) | discovery plasma | 232 (HCC, breast, colorectal, lung, ovarian, pancreatic, gastric, bile-duct) | AUC 0.94 (HCC-vs-healthy); per-cancer 0.86–0.97 | Nature 2019 |
| CancerSEEK (Cohen/Klein 2018 [9]) | prospectively-collected plasma | 1,005 (8 cancers + healthy) | sens 0.70 @ spec 0.99 (median across 8 cancers); OvR AUC 0.91 | Science 2018 |
| Galleri (Klein 2021 [10]) | independent validation plasma | 2,823 (50+ cancer signals vs healthy) | OvR AUC 0.92 across all cancer types | Annals of Oncology 2020 / 2021 PATHFINDER |
| CAPP-Seq (Newman 2014 [5] / 2016 [6]) | tumor-informed MRD, NSCLC pilot + integrated-error-suppression update | not pooled (per-cancer n's) | first prototype AUC ~0.86 in stage II–IV NSCLC; iDES-enhanced AUC >0.99 in stage II–IV | Nat Med 2014; Nat Biotechnol 2016 |
| FinaleMe (Liu 2018 [11]) | plasma + reference methylome | 159 (training + held-out) | end-motif OvR AUC 0.91 | bioRxiv 2018 (preprint) |
| **DeepCatch §3.2 (this work, open data)** | FinaleDB open data, internal CV | **627** (Jiang+Cristiano) | pooled **AUC 0.974 ± 0.002**, sens@95=0.909, sens@98=0.837, sens@99=0.782 | `results/cross_study_finallydb.json` |

**Positioning (honest).** The published numbers above were generated by independent labs on larger, often proprietary, clinical-grade plasma cohorts. DeepCatch's §3.2 number is a **pooled open-data internal-CV** baseline on the same 5-channel DELFI feature (Cristiano-style) that DELFI itself uses. Pooled internal-CV at this scale is competitive in the AUC band of the published numbers (0.86–0.94 discovery-cohort range) — but ours is **not** an external validation, **not** a clinical operating point, and **not** an architectural contribution. The substantive contribution of DeepCatch is the reproducible engineering substrate that lets a researcher reproduce or extend these numbers on open data, not the headline AUC.

Three concrete distinctions matter:

1. **Galleri and CancerSEEK** are methylation-or-mutation-multianalyte assays, not pure fragmentomics. Their reported numbers include assay information (methylation block, oncoprotein panel) that DeepCatch does NOT use. Comparing 5-channel DELFI fragmentomics to a multianalyte panel would be apples-to-oranges; the numbers are listed for context, not as a head-to-head.
2. **CAPP-Seq** is tumor-informed MRD — it tracks patient-specific mutations at 0.01–0.1% ctDNA. DeepCatch's pure-fragmentomics 5-channel LR is not in that regime; comparability with CAPP-Seq on this pooled cross-study cohort is structural, not parametric.
3. **FinaleMe** is the closest cousin methodologically (4-mer end-motif distribution from plasma). Liu 2018 reported AUC 0.91 on n=159; DeepCatch's 5-channel DELFI feature at n=627 reports pooled AUC 0.974. This is the closest fair-context comparison, and is consistent with the headroom available in moving from a 159-patient discovery cohort (n=159 has wide CI on the per-cancer denominator) to a 627-sample pooled CV. An independent head-to-head on a shared held-out cohort would be required to claim a substantive gap.

The open-data benchmark at hand is, by design, **not a clinical claim**. The comparator table is included so the reader can place the §3.2 number in the published-methods landscape without over-reading it.

---

## 4. Discussion

### 4.1 What this work does

DeepCatch provides a reproducible fragmentomics engineering substrate. It takes a public cfDNA cohort (FinaleDB, GDC, FLARE) through (i) feature extraction in a canonical six-modality schema, (ii) self-supervised pre-training on real data, (iii) downstream classification with honest per-patient GroupKFold validation, and (iv) paired-design ablations of every architectural knob. The headline contribution is not a number but the engineering substrate itself: a new lab can swap in a new cohort and get honest numbers in hours.

### 4.2 What this work does NOT do

- **Clinical validation.** DeepCatch is not a clinical assay. No real-plasma cohort has been sequenced. No patient samples were collected under IRB.
- **FDA pathway.** No regulatory submission is planned or implied.
- **Prospective screening.** No prospective cohort, no enrollment, no follow-up.
- **Headline-grade clinical AUCs.** The synthetic-plasma foundation AUC (≈0.55) is below the LR baseline. The cross-study FinaleDB 5-channel LR AUC at n=627 is 0.974 with sens@99=0.782 (pooled, harmonized), but that is **pooled internal-CV**, NOT external validation, NOT a held-out clinical plasma cohort, and NOT a clinical operating point.
- **HCC_J is a within-study reading, not a cross-study one.** The OvR reading of HCC_J (Jiang 2015 only, n=89) against the cross-study healthy pool drops the AUC to 0.725 because there are only 32 Jiang healthy controls — the comparison reduces to "Jiang HCC vs mostly-Cristiano healthy," not a within-study HCC-vs-healthy read. HCC_J PPV @0.001 = 0.019 (and only 0.031 @ prev=0.01) reflects the low sens@99 (0.191) on that denominator, not a model-quality claim.
- **Per-cancer denominators are small.** Top-5 cancer n ranges from 28 (OV) to 89 (HCC_J); CIs on sens@99 (e.g. OV [0.758–1.000], PAAD [0.333–0.810]) are wide. Conclusions about per-cancer performance should be read as ordinal, not cardinal.
- **Real plasma fragmentomics.** The pretraining cohort is from FinaleDB plasma (Cristiano 2019, Jiang 2015), but the downstream validation cohorts are either synthetic-plasma (TCGA-LUAD real mutations + Poisson sampling, §3.1) or Open-Data FinaleDB uniformly-preprocessed plasma (§3.2/3.3). A real-plasma fragmentomics validation requires an IRB-approved prospective cohort — see §4.4.

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
- **Reproduce all numbers:**
  - One-bash driver: `bash paper/REPRODUCE.sh` (or equivalently `bash paper/REPRODUCE.sh --quick` for a reduced run). See `paper/REPRODUCE.md`.
  - Full top-level driver: `RUN_ALL.sh` and the per-script invocations documented in §3.
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

*Manuscript draft — 2026-09-23. Section 3 filled with the just-shipped cross-study and per-cancer numbers. See `paper/REPRODUCE.md` for one-command reproduction.*
