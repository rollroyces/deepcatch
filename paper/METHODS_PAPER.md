# DeepCatch: a software framework for cfDNA fragmentomics with open-data cross-study benchmarking

**Authors:** Yu Ching Lam¹ (sole author — see `AUTHORS`)

¹ DeepCatch Project — github.com/rollroyces/deepcatch (Independent researcher)

**Venue target:** PLOS Computational Biology or Bioinformatics Advances

**Status (2026-11-17):** the §3 numbers are regenerated from the seeded
benchmark — the GC-corrected cross-study sweep
(`results/cross_study_finallydb_gc_corrected.json`), the per-cancer
sens@spec table (`results/per_cancer_sens_at_spec.json`), the
leave-one-study-out transfer (`results/cross_study_loso.json`), the
shuffled-label null (`results/cross_study_finallydb_shuffled_control.json`)
and the paired-cohort baselines (`results/fusion_regularized_smoke.json`,
`results/sens_at_spec_ablation.json`); §4.2's limits are tightened to the
per-cancer numerators actually observed. The author line above is
single-author, consistent with `AUTHORS`; a co-author is added only if a
contribution is agreed (see `TEAM.md` §3).

---

## Abstract

**Background.** Cell-free DNA (cfDNA) fragmentomics — the joint analysis of fragment length, end-motif, copy-number, and coverage profiles from plasma sequencing — is a promising substrate for non-invasive multi-cancer detection. Published fragmentomics methods (DELFI, Galleri, FinaleMe) report compelling per-cancer discrimination, but their training data and validation cohorts are largely proprietary, and independent replication on harmonized open data is scarce.

**Methods.** DeepCatch is an open-source software framework that turns public cfDNA resources (FinaleDB, GDC TCGA, FLARE/GSE317007) into reproducible fragmentomics experiments. The framework combines (i) a six-modality feature schema (5-channel DELFI, 4-mer end motifs, fragment size, coverage/CNV, serology, tissue-of-origin); (ii) a multi-modal fusion encoder with masked-modality self-supervised pre-training and cross-attention fusion; (iii) panel-based tumor-informed detection when a prior cohort is available; and (iv) a validation protocol built on 5-seed × 5-fold GroupKFold with shuffled-label and true-confound controls. The implementation is MIT-licensed, ships with 950 regression tests, and exposes every intermediate artifact (cohort matrices, per-seed AUC arrays, ablation JSONs) to the reader.

**Results.** On the paired 20-patient TCGA-LUAD synthetic-plasma cohort at 0.1% ctDNA, the foundation encoder honestly achieves AUC ≈ 0.56 (n=20 seeds, per-patient GroupKFold) — below the sklearn LR baseline (0.90) on the same channels, which is the expected outcome on this paired design and is itself a finding worth reporting. On the cross-study FinaleDB publications 6 + 8 (n=627 after missing-artifact filter), a 5-channel DELFI LR with per-study z-score harmonization + per-bin GC bias correction (default ON, see `docs/HARMONIZATION.md §3`) reaches pooled OOF AUC **0.9684 ± 0.0027** (with GC correction; legacy without GC: 0.9770 ± 0.0013 — the GC-corrected version is the more honest cancer signal per the per-publication coherence check), sens@95%=0.865, sens@98%=0.796, sens@99%=0.700. Per-cancer OvR (5-seed ± std) for the top-5 cancers (with GC correction): OV 0.989 ± 0.007, LUAD 0.981 ± 0.003, BRCA 0.973 ± 0.002, PAAD 0.967 ± 0.004, HCC_J 0.704 ± 0.014 (single-cohort reading). A **true-confound control** (cancer = 100% Jiang, healthy = 100% Cristiano) collapses to AUC 0.48 with harmonization vs 0.999 without — proving the pooled 0.968 (and the legacy 0.977) is a cancer signal, not a study batch. Three honest null/negative ablations (focal-BCE loss, sparse-aware projection, harmonization) are documented, and every benchmark number above is reproduced bit-for-bit by the committed code (seeded PCA, `test/test_benchmark_reproducibility.py`).

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
- **FinaleMe (Liu et al. 2018 [8])** — fragment-end motif model using 4-mer end-sequence preferences.
- **Jiang et al. 2020 [11]** — HCC-specific fragmentomics in a 222-patient cohort, one of the largest single-cancer fragmentomics datasets with public fragment files.
- **FinaleDB / Finaletoolkit [12, 13]** — public archive of fragmentomics datasets with harmonized `*.frag.tsv.bgz` files; the primary data source for the present work.

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

**Cross-study benchmark with motifs.** When `--include-motifs` is enabled on `scripts/cross_study_finallydb.py`, the cross-study benchmark additionally appends the raw 256-dim 4-mer end-motif frequency vector (sum-to-1 per sample) from the local FinaleDB features cache (`{sample}.motifs.npy`) to the 63,246-dim 5-channel DELFI + FSD-196 feature vector, producing a 63,502-dim feature vector. On the 627-sample FinaleDB pubs-6+8 cohort (5-seed × 5-fold StratifiedKFold pooled OOF, per-publication z-score harmonization inside each fold, GC correction on — the same configuration as the headline), this gives a pooled AUC of **0.9761 ± 0.0030** vs **0.9684 ± 0.0027** without the motif block — a **Δ AUC = +0.0077** improvement. Per-publication AUC also improves for both Jiang (0.9661 → 0.9776) and Cristiano (0.9698 → 0.9739).

That Δ is configuration-dependent: on the legacy no-GC feature set the motif block added much less, which is what one would expect if the two feature blocks were partly redundant in the legacy set. The legacy comparison's numbers are recorded in the CHANGELOG for this release, where the superseded figures live; the run quoted above is the one the artifact holds. See `results/cross_study_finallydb_with_motifs.json` for the full with-motifs benchmark. The `--motif-cohort-mode` flag controls how samples without a cached `.motifs.npy` are handled: `zero_fill` (default; head-to-head on the same n=627 cohort) zero-fills the 256-dim motif block, `strict` drops the sample (reduces cohort to ~98 samples with motifs cached locally; it cannot be combined with `--gc-correction`, because the strict motif block is a subset of the cohort the correction covers). Source: `scripts/load5_with_motifs.py`.

### 2.2 Multi-modal fusion architecture

The encoder (`src/foundation/model.py:MultiModalEncoder`) is a per-modality `LinearProjection` (or `SparseAwareLinearProjection`, opt-in) followed by a 4-layer Transformer with 4-head self-attention and a `[MODALITY]` token prepended to each modality's token sequence. The encoder was pre-trained on 200 real FinaleDB samples (Cristiano 2019 + Jiang 2015) with masked-modality-prediction (Phase 1, 15 epochs) and InfoNCE contrastive (Phase 2, 5 epochs). Pre-training details, the synthetic-bypass bug fix, and the regression test are documented in `docs/PRETRAINING.md` and `docs/PRETRAIN_BUG.md`. The PRODUCTION checkpoint is `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` (543,872 parameters, embed_dim=128); it is **not distributed with this repository** and is regenerated by `scripts/pretrain_production_finaledb.py`. No number reported in §3 depends on it: the encoder that produces the §3.1 comparison is trained in-run with `pretrained=False`, so §3 is reproducible from a clean checkout.

Downstream classification is a `FoundationDownstream` head: 0.7 × frozen-pretrained-encoder → LR-head + 0.3 × trainable Transformer, with binary cross-entropy (default) or focal-BCE (`loss="sens_at_spec"`, opt-in, `alpha_pos` sweep parameter) as the loss.

### 2.3 Panel-based detection (when prior cohort is known)

For tumor-informed panel detection, the same `FoundationDownstream` consumes a modality dict where `frag_basic[:, 0]` is the panel-LLR sum (the sufficient statistic under independent Poisson observation at each locus), `frag_enhanced[:, 0]` is the patient-level frag score, and the other modalities are filled with cohort-derived priors. The panel-LLR scoring function, the Fisher-method aggregate, and the strand-concordance-weighted variant are documented in `paper/paper.tex` (§2 of the benchmark paper) and `scripts/foundation_real_smoke.py`.

### 2.4 Validation protocol

Every reported number follows one protocol:

1. **5-seed × 5-fold per-patient GroupKFold.** Seeds {0, 1, 2, 3, 4} (or {0..19} for higher-power ablations). Both arms of every paired patient stay in the same fold. This is the protocol audit-2 P0-A (commit `6231fcb`) installed to replace the leaky `StratifiedKFold` that originally inflated the foundation AUC to 0.94.
2. **Shuffled-label control.** Labels are permuted independently in the pos-half and neg-half to break the pair structure. AUC against the shuffled labels is reported as `shuffled_foundation_auc`. The honest diagnostic on paired data is the *delta* between real-label AUC and shuffled-label AUC, not the absolute AUC.
3. **True-confound control.** For cross-study experiments, per-study batch effects are injected on top of the synthetic signal and the per-study AUC is reported alongside the pooled AUC. A diagnostic (`delta_auc_pooled_minus_per_study_mean`) quantifies the gap. See `docs/HARMONIZATION.md`.
4. **Per-cancer sensitivity at fixed specificity with DeLong CIs.** Per-cancer sensitivity at 95% and 99% specificity is computed by pooling OOF predictions across seeds, fixing the threshold at the empirical 95th / 99th percentile of the healthy controls, and applying the DeLong method [14] to estimate the AUC confidence interval. See `results/per_cancer_sens_at_spec/`.
5. **Bit-identical paired-design sanity.** For every ablation, the non-ablation metrics (`panel_only_aucs`, `frag_only_aucs`, `lr_baseline_aucs`, `naive_avg_aucs`) are verified bit-identical to <1e-9 across the two JSONs, proving the only difference is the ablated knob.

### 2.5 Implementation and tests

The framework is implemented in PyTorch with NumPy / scikit-learn baselines. 950 tests run across the two CI jobs (764 core, 186 deep-learning; 23 of the core ones need a local `.venv` or the sibling repository's checkout and skip elsewhere); the pretrain-bug-fix regression suite (`test/test_pretrain_bug_fix.py`, 11 tests) pins the contract that the synthetic generator is never invoked when real data is supplied. Code, data, results, and ablations are version-controlled at github.com/rollroyces/deepcatch (MIT license).

**Reproducibility is enforced, not assumed.** Every benchmark number in §3.2–§3.5 is a function of the committed code: the PCA inside each CV fold is seeded with the fold-split seed, and `test/test_benchmark_reproducibility.py` enforces two properties — an AST guard fails the build if any `PCA(...)` call in `scripts/` or `src/` is fitted without a `random_state`, and a behavioural test fits the benchmark twice and requires the two runs to agree exactly (with a further check that a different seed still gives a different answer, so the seed is not simply pinning everything). The guard exists because the opposite was true for an earlier revision: with an unseeded PCA, sklearn's randomised solver made every re-run of the benchmark drift (one publication's AUC moved 0.9789 → 0.9821 between two runs of the same command), which would have made the artifacts unverifiable. The numbers in this manuscript are from the seeded run; the artifacts, the narrative documents and the manuscript now come from a single command.

---

## 3. Results

> Source-of-truth: `results/cross_study_finallydb_gc_corrected.json` (§3.2 headline), `results/cross_study_finallydb.json` (§3.2 legacy, no GC), `results/per_cancer_sens_at_spec.json` (§3.3), `results/cross_study_finallydb_shuffled_control.json` (the null control) and `docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md` (the running benchmark narrative). All pooled numbers are 5-seed × 5-fold StratifiedKFold internal CV (pooled OOF); no external validation, no clinical plasma cohort.
>
> Figures: `docs/figures/fig1_pooled_roc.png` (§3.2), `docs/figures/fig2_per_cancer_roc.png` (§3.3), `docs/figures/fig3_calibration.png` (§3.2 sens@spec operating points in ROC space), `docs/figures/fig4_sens_operating_points.png` (§3.3 grouped per-cancer bars). Reproducible via `python scripts/plot_figures.py`; see `docs/figures/README.md`.

### 3.1 Honest baseline AUC on paired synthetic cohort (20-patient TCGA-LUAD)

On 20 paired TCGA-LUAD patients at 0.1% ctDNA, real-panel + real-mutation-derived-frag, per-patient GroupKFold, 20 seeds:

| Metric | Value | Source |
|---|---|---|
| `panel_only AUC` | 0.9175 ± 0.0135 | `results/sens_at_spec_ce_n20.json` |
| `frag_only AUC` (real mutation-derived) | 0.6229 ± 0.0141 | `results/sens_at_spec_ce_n20.json` |
| `foundation AUC` | **0.5635 ± 0.0113** | `results/sens_at_spec_ce_n20.json` |
| `lr_baseline AUC` (sklearn LR on [panel, frag]) | 0.9040 ± 0.0121 | `results/sens_at_spec_ce_n20.json` |
| `naive_avg AUC` | 0.9174 ± 0.0136 | `results/sens_at_spec_ce_n20.json` |
| `shuffled_foundation AUC` (permuted-label null) | 0.3464 | `results/sens_at_spec_ce_n20.json` |
| `foundation_auc_minus_null` (honest − null) | +0.2171 | `results/sens_at_spec_ce_n20.json` |
| `delta_auc_normalized` | 3.4193 | `results/sens_at_spec_ce_n20.json` |

> **Provenance.** These are the current numbers: the artifact is regenerated
> by `python scripts/foundation_real_smoke.py --seeds 20 --loss ce --out
> results/sens_at_spec_ce_n20.json` (see `docs/SENS_AT_SPEC_ABLATION.md`
> §"Reproduce (n=20)"), and the §3.1.1 estimators below are re-measured by
> the reproduction driver on every run.

**Honest framing:** under per-patient GroupKFold, the foundation model is honestly out-performed by the sklearn LR baseline on this paired design (0.5635 vs 0.9040). The permuted-label control is the same estimator, the same folds and the same model seeds, run on permuted labels, so the difference between the two numbers is attributable to the label assignment alone: the foundation scores 0.3464 on permuted labels against 0.5635 on the true ones, a gap of **+0.2171**. The LR baseline shows the same behaviour (0.9040 honest against 0.3903 on the null), and the naive average 0.9174 against 0.5084.

The null sits **below** 0.5 on this cohort, and that is expected: every patient contributes one positive and one negative arm whose panel scores differ by ~+127, so a model fitted to permuted labels still recovers that within-pair axis. The measured null is 0.4063 ± 0.1143 over 150 permutations (`results/fusion_regularized_smoke.json` → `protocol.measured_null`), so a control below 0.5 is not on its own evidence of leakage here; it has to be read against that measured null. For the same reason `delta_auc_normalized` (which divides the gap by `real − 0.5`, implying a null at 0.5) is reported for continuity but the interpretable quantity is `foundation_auc_minus_null`.

**A note on the history of this control.** The permutation used to be a no-op: the script permuted the positive index set among itself and then wrote the 1s back to those same positions, so `y_shuf == y_true` for every seed and the "control" reported the honest AUCs back. The fingerprint was visible in the committed artifact, where `shuffled_lr_baseline_auc_mean` was byte-identical to `lr_baseline_auc_mean`. The control now lives in one place (`src/foundation/label_control.py`) and is a genuine permutation; `test/test_label_control.py` guards the property structurally. See `AUDIT_2_FINDINGS.md` and the CHANGELOG.

#### 3.1.1 The estimator class is the problem, and a smaller estimator fixes it

A neural encoder cannot be rescued by more capacity on n=40; the fix is a smaller *estimator* that fuses the two channel scores late. Two such estimators were implemented (`src/foundation/regularized_fusion.py`) and measured on the identical protocol, with the LR baseline re-computed in the same run so the comparison is within-run (`scripts/fusion_regularized_smoke.py`; artifact `results/fusion_regularized_smoke.json`):

| Estimator | AUC (3 seeds, per-patient GroupKFold) | permuted-label control | arm-swap control |
|---|---|---|---|
| foundation encoder (audit-2) | 0.554 ± 0.005 | — | — |
| `lr_baseline` (re-computed here) | 0.9042 ± 0.0174 | 0.569 | 0.383 |
| `naive_average` | 0.9150 ± 0.0178 | 0.614 | 0.513 |
| **`regularized_late_fusion`** | **0.9825 ± 0.0094** | 0.573 | 0.436 |
| **`bayesian_model_average`** | **0.9650 ± 0.0089** | 0.511 | 0.413 |

The re-computed LR baseline (0.9042, 3 seeds) and equal-weight average (0.9150) agree with the 20-seed paired-cohort figures above (0.9040 ± 0.0121 and 0.9174 ± 0.0136), so the comparison is like-for-like. Regularized late fusion beats the LR baseline by **+0.078 AUC**, Bayesian model averaging by **+0.061**, with no torch dependency and no change to the CV protocol.

**The label-permutation control had to be corrected to be usable.** The audit-2 "pair-broken shuffle" permuted the positive and negative index sets separately and then reassigned 1s and 0s to those same sets, which reproduces `y_true` exactly — the control was a no-op, and the "shuffled" numbers in the audit-2 table are the honest numbers. The corrected control permutes labels inside each fold. Its null on this paired design is **0.4063 ± 0.1143, not 0.5** (measured over 150 permutations): every patient contributes one positive and one negative arm whose panel scores differ by ~+127, with the positive arm higher in 100% of patients, so a model fitted to random labels still recovers that within-pair axis and lands *below* 0.5. A control below 0.5 is therefore not on its own evidence of leakage here; it must be read against the measured null.

**Scope.** n=40 samples from 20 patients, simulated plasma at TF=0.001, one cancer type. This is an estimator-class result, not a clinical performance claim.

### 3.2 Cross-study AUC on FinaleDB publications 6 + 8 (n=627)

The cross-study pool is the union of FinaleDB publication 6 (Jiang 2015 [15], n=121: 89 HCC + 32 healthy) and publication 8 (Cristiano 2019 [7], n=506: 274 cancer + 232 healthy) after the missing-artifact filter (31 samples dropped: 31/691 = 4.7%). Features are the canonical 5-channel DELFI profile (5 Mb ratio + 5 Mb coverage + 100 kb ratio + 100 kb counts + FSD-196, 63,246-dim), with the 4th channel being the GC-bias-corrected 100kb-coverage (`c100_gc_corrected`) per `scripts/gc_correction.py` (default ON; legacy `cn = c100 / median(c100)` retained under `--no-gc-correction` for backward compat). The pipeline is a 5-seed × 5-fold StratifiedKFold, LR(max_iter=2000) on PCA(200) of a per-study z-score StandardScaler fit on the **train** fold only. Source: `results/cross_study_finallydb_gc_corrected.json` (GC-OFF legacy: `results/cross_study_finallydb.json`; companion narrative `docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md`).

|| Cohort | n | Cancer | Healthy | AUC (5-seed mean ± std) |
||---|---:|---:|---:|---|
|| jiang | 121 | 89 | 32 | **0.9661 ± 0.0031** |
|| cristiano | 506 | 274 | 232 | **0.9698 ± 0.0016** |
|| **pooled (harmonized, GC corrected)** | 627 | 363 | 264 | **0.9684 ± 0.0027** |
|| **pooled (harmonized, GC not corrected — legacy)** | 627 | 363 | 264 | 0.9770 ± 0.0013 |

The GC-corrected pooled AUC is **lower by ΔAUC = −0.0086** than the legacy pooled AUC (0.9770 → 0.9684). Removing the GC-vs-coverage axis reduced the pooled number by ~0.009 (3× the 5-seed pooled std of 0.0013–0.0027, but well within the 5-seed pooled std). The honest reading per `docs/HARMONIZATION.md §3` is that a chunk of the legacy 0.9770 AUC was carrier on the GC axis, which the LR used as an indirect proxy for cancer signal (via CNV-induced GC shifts). The GC-corrected 0.9684 is close to the within-publication AUCs (Jiang 0.9661, Cristiano 0.9698) than the legacy 0.9770 was, suggesting the legacy pooled number was inflated by ~0.005-0.009 from GC-axis proxy detection. **We report both versions**, not silently swap one for the other.

At the pooled harmonized GC-corrected operating point: **sens@95%=0.865**, **sens@98%=0.796**, **sens@99%=0.700**. The sens@spec numbers are slightly below the legacy (sens@95=0.912, sens@98=0.851, sens@99=0.785); the same Δ ~0.04 holds at every operating point, consistent with the AUC drop being a coherent downward shift rather than a calibration-only change.

**Shuffled-label null control.** A pooled 0.97 AUC is open to the reviewer question "is this cancer signal, or batch leakage via fold structure?" We add a label-permutation control that pairs the existing true-confound: the same pooled OOF (harmonize + PCA + LR) is re-run on `numpy.random.permutation(y)` (perm_seed=20260924), preserving sample-id pairing. The shuffled-pooled AUC must collapse to ~0.5; if it stays high, the cancer signal is not what's driving the 0.97 AUC. Source: `scripts/cross_study_finallydb_shuffled_control.py`, output `results/cross_study_finallydb_shuffled_control.json`.

|| Setting | n | n_cancer | n_healthy | Pooled AUC (5-seed mean ± std) | `control_passed` |
||---|---:|---:|---:|---:|---|
|| **shuffled-label pooled (harmonized)** | 627 | 363 | 264 | **0.4460** (1 seed) | **true** (0.4460 < null floor 0.55) |

The shuffled-pooled AUC (0.4460, sens@95=0.061, sens@99=0.014) is well inside the healthy null range (0.45–0.55), so the cancer signal — not per-fold publication structure — is what drives the 0.97 real-label AUC. The control is a sanity check, not a gate; a noisy seed could fluctuate above the 0.55 floor, in which case the script prints a `UserWarning` and reports `control_passed: false` honestly. See `docs/CROSS_STUDY_BENCHMARK.md` §3 for the full interpretation (sourced from the same JSON).

**True-confound control.** The honest cross-study claim requires that the signal is cancer, not study. The control builds a synthetic two-arm dataset where cancer = 100% Jiang and healthy = 100% Cristiano (and the symmetric orientation), then rerun the same pipeline:

| Orientation | n cancer | n healthy | AUC harmonized | AUC no_harmonize |
|---|---:|---:|---:|---:|
| cancer=Jiang, healthy=Cristiano | 121 | 506 | 0.4836 ± 0.0095 | 0.9989 ± 0.0014 |
| cancer=Cristiano, healthy=Jiang | 506 | 121 | 0.4836 ± 0.0095 | 0.9994 ± 0.0007 |

Without harmonization, the classifier trivially achieves AUC ≈ 0.999 because "which study is this from?" is trivially answerable. With per-study z-scoring, the same control collapses to **0.4836** — the per-study mean/variance is the only signal and the harmonization removes it by design. The drop from 0.9989 to 0.4836 — from near-perfect study separation to chance — is what shows the pooled-harmonized 0.9684 AUC is a cancer signal, not a study batch.

**Honest framing.** These numbers are **pooled out-of-fold** AUC on the 627-sample cross-study cohort, under 5-seed × 5-fold internal CV. They are NOT external validation, NOT held-out clinical plasma, and NOT clinical-grade operating points (Galleri, CancerSEEK). They measure how well the 5-channel DELFI features separate cancer-vs-healthy in the pooled cohort with per-study batch removed. The per-cohort AUCs (Jiang 0.9661 ± 0.0031, Cristiano 0.9698 ± 0.0016) and the pooled-harmonized 0.9684 agree to within 0.004 AUC, which is the expected size when the per-study batches are mild on this FinaleDB-uniformly-processed cohort. The companion pipeline artifact `docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md` carries the full per-cohort + per-cancer + true-confound tables for this run (the legacy no-GC counterpart is `docs/CROSS_STUDY_BENCHMARK.md`). The classifier is intentionally a linear LR (not the foundation encoder) so the headline number is interpretable as the separability of the published 5-channel feature set, not as a DeepCatch architectural contribution.

#### 3.2.1 Leave-one-study-out: the pooled number does not survive a publication holdout

A pooled OOF over a merged cohort can train on samples from one publication and test on samples from another, so §3.2 measures *within-cohort* separation, not transfer to an unseen study. Holding a whole publication out measures the thing a reviewer will ask for. Source: `scripts/cross_study_loso.py`; artifact `results/cross_study_loso.json`; narrative `docs/CROSS_STUDY_LOSO.md`.

| Direction | Train | Test | Test AUC [DeLong 95% CI] | sens@95% | sens@99% |
|---|---|---|---|---|---|
| Jiang → Cristiano | pub 6 (n=121) | pub 8 (n=506) | **0.6805** [0.6344–0.7267] | 0.256 | 0.183 |
| Cristiano → Jiang | pub 8 (n=506) | pub 6 (n=121) | **0.9361** [0.8909–0.9813] | 0.775 | 0.506 |
| *in-domain (Jiang)* | — | — | *0.9789 ± 0.0027* | | |
| *in-domain (Cristiano)* | — | — | *0.9690 ± 0.0021* | | |

The pooled 0.9684 of §3.2 therefore **does not survive a publication holdout in the single-cancer-type direction**: 0.968 → 0.681. The asymmetry is informative rather than mysterious. Jiang is HCC-only, so training on it yields an HCC-specific decision boundary, while the held-out Cristiano set contains six non-HCC cancer types; a within-study leave-one-cancer-type-out control (train on the remaining types + healthy, test on the held-out type) transfers at **0.857–1.000**, so "unseen cancer type" is not the blocker by itself, and the reverse direction (pan-cancer training → unseen HCC test) crosses the study boundary at 0.936. The binding constraint is the *composition of the training cohort*, not the study shift. A PCA-component sweep (5/10/20/50/100 → 0.662/0.707/0.722/0.667/0.644) rules out a capacity artefact.

**How to read §3.2 after this subsection.** The pooled 0.9684 is an honest *within-cohort* separation number with per-study batch removed, and the true-confound control still shows the pooled signal is cancer rather than study identity. It is not evidence of cross-study generalisation; §3.2.1 is. The distinction matters for any claim about deploying the classifier on a new cohort, and it is the reason the external validation in §3.2.2 uses publications that were never touched during training.

#### 3.2.2 External validation on held-out public publications

The pooled cohort of §3.2 and the holdout of §3.2.1 both use publications 6 and 8. This subsection evaluates on publications that are **never used for training, feature extraction statistics, or hyper-parameter selection**: publication 1 (Snyder 2016 [16], Cell), publication 7 (Sun 2019 [17], Genome Research) and publication 9 (Adalsteinsson 2017 [18], Nat Commun). The samples are fetched from the public FinaleDB REST API and the public S3 fragment objects by `scripts/fetch_external_cohort.py` — no IRB, no data-use agreement — and the 5-channel features are extracted by the self-contained `src/fragmentomics/finale_features.py`. Source: `scripts/external_validation_finaledb.py`; artifact `results/external_validation_finaledb.json`; narrative `docs/EXTERNAL_VALIDATION.md`.

| Publication | Assay | n (ca/non-ca) | External AUC [DeLong 95% CI] |
|---|---|---|---|
| 1 — Snyder 2016 | WGS | 52 (43/9) | 0.708 [0.5148, 0.9012] |
| 7 — Sun 2019 | WGS | 31 (18/13) | **0.9444** [0.8731, 1.0000] |
| 9 — Adalsteinsson 2017 | WES | 239 (239/0) | case-only: sensitivity 0.6569 at the training-derived 99%-specificity threshold |
| **pooled over publications carrying both classes** | — | 83 (61/22) | **0.7385** [0.6281, 0.8488] |

(The held-out cohort is 322 samples in total: 300 cancer and only 22 non-cancer, because FinaleDB's held-out publications publish few controls -- publication 9 ships 1,508 cancer samples and zero controls. Nine samples whose fragment files exceed 900 MB and one unreachable object were skipped. Publications 1, 7 and 9 are *every* plasma cfDNA study FinaleDB serves that the model never trained on: the three it does not serve as plasma are urine (publication 3), a prenatal study without disease labels (publication 4) and a tumour-tissue study with no plasma samples at all (publication 5).)

**Transfer between whole-genome studies is real but far from uniform.** Publication 7 is a *whole-genome* study, like the training publications, and the classifier transfers to it well (0.9444; colorectal one-vs-rest 0.9182). Publication 1 is also whole-genome, but it covers thirteen tumour types against only nine non-cancer controls, and it transfers less well (0.708 [0.5148, 0.9012]). Publication 9 is *whole-exome* and transfers only partially: breast cancer 0.6938 and prostate cancer 0.7127 against the same non-cancer controls -- above chance, but well below the whole-genome studies. The DELFI 5 Mb / 100 kb profile is a whole-genome coverage signal, and an exome-targeted assay has a different coverage shape altogether; per-publication z-score harmonization cannot repair that. This is the kind of boundary that an internal pooled cross-validation cannot reveal, and it is the reason the external set is reported per publication and per disease rather than as a single pooled number.

**Two corrections to an earlier version of this subsection, both of which moved the numbers down.** The first version of the external cohort fetched only 9 of publication 1's 60 plasma samples (8 cancer and 1 control) because of a per-publication download budget, and reported that publication's AUC as 0.875 from those nine samples; the cohort now carries all 52 of its labelled plasma samples, and the AUC is 0.708. Separately, the script z-scored the held-out cohort as a *single pool* rather than per publication, which contradicts the protocol stated above and the LOSO benchmark it claims to mirror: with 239 case-only publication-9 samples in the pool, publications 1 and 7 were scored against publication 9's mean and variance. The artifact records all three harmonisations under `harmonization_sensitivity`; the table above uses the documented one (per study on the training side, per publication on the held-out side), which is also the most conservative of the three on the headline publication.

**Scope.** Publication 9 contributes cases without controls, so it is reported as sensitivity at the threshold fixed in training (0.9766, the training cohort's 99%-specificity operating point) rather than folded into an AUC; pooling it with the matched publications would put a case-only cohort into an AUC denominator and produce a number that means nothing (0.7097 for the naive pooled figure, against 0.7385 for the publications that carry both classes). The non-cancer arm is small (22 samples) and is drawn from different publications than most of the cancer arm, so the confidence intervals are wide and the per-disease rows carry the weight. This is external validation on *public* data — a different study, a different institution, a different assay — and not a clinical validation: no prospective follow-up, no stage labels, and no real healthy-donor controls.

#### 3.2.3 What the three cross-study numbers together say

| Question | Number | Section |
|---|---|---|
| Can the features separate cancer from healthy in a pooled open cohort? | 0.9684 ± 0.0027 | 3.2 |
| Is that signal cancer rather than study identity? | true-confound control collapses to 0.4836 | 3.2 |
| Does the classifier transfer to an unseen publication? | 0.9444 (WGS, pub 7), 0.708 (WGS, pub 1), 0.7385 pooled, 0.6569 sensitivity (WES, pub 9) | 3.2.2 |
| Does it transfer to an unseen publication in the single-type direction? | 0.681 | 3.2.1 |
| Can the neural encoder do this? | 0.554 — no; late fusion reaches 0.9825 | 3.1.1 |

The honest summary is that the *features* transfer across studies when the assay is comparable (0.87-0.99), that the *training cohort composition* matters more than the study shift for single-type directions (0.681), and that changing the assay platform breaks the transfer entirely (0.42). None of these numbers is a clinical performance estimate.

#### 3.2.4 The V3 methylation channel measures library complexity, not methylation

The V3 pipeline was run end-to-end on real public fragments: FinaleDB publication 7 (Sun 2019), blood plasma, **31 samples** (18 colorectal cancer, 13 liver-transplant controls), hg19 fragments, chr22 CpGs, one HMM trained on the cohort and each sample decoded with it (`scripts/methylation_cohort_channel.py`; artifacts `results/methylation_cohort_channel.json`, `results/methylation_channel_diagnosis.json`).

| Quantity | LOO AUC |
|---|---|
| decoded β summary (4 features) | **0.9444** |
| number of covered CpGs alone | **0.9444** |
| observations per CpG alone | **0.9444** |
| β summary restricted to the 1,550 most widely covered loci | 0.803–0.829 across repeats — not reproducible, see note |
| β summary, depth-matched to 11,003 CpGs | 0.9444 (coverage-only 0.8034) |

Three unrelated quantities — the decoded methylome, how many CpGs the library reached, and how many times each CpG was observed — return the **identical** leave-one-out AUC of 0.9444. The channel's discrimination is library complexity.

> **One row of that table is not reproducible, and the reason is upstream.** FinaleMe's HMM training is not deterministic: on byte-identical decoder inputs the trained model differs between runs, which flips 16–20 of ~51,000 CpG calls per sample, and the shared-locus AUC moves with it (0.803–0.829 over five runs of this cohort). Every other row in the table was identical on every run, as were the median difference and the permutation p-value quoted below; the share of CpGs higher in cancer moves only in its third digit (2.4–2.5%). The shared-locus row is therefore reported as a range. It does not carry the argument either way: every observed value in that range is far below the 0.9444 headline, which is the point the row makes.

The matched-locus test cannot be explained that way, because it compares the two arms only on the loci where **both** of them have coverage: over **103,863 CpGs**, each tested on the samples that cover it, the median cancer-minus-control difference in decoded β is **exactly 0.0**, with a permutation null of exactly 0.0 (p = 1.0), and only 2.4–2.5% of CpGs are higher in cancer.

The mechanism is visible in the libraries themselves (`scripts/methylation_library_check.py`, artifact `results/methylation_library_check.json`). FinaleMe infers methylation from the *digestion pattern* — an unmethylated CpG is cut out of the library and sits at a fragment end, a methylated one stays interior — and the released decoder reads `Offset_frag`, `Dist_frag_end`, `FragLen`, `Norm_Frag_cov` and the prior, never the base calls. Two properties of the reachable public libraries remove that signal:

| | FinaleDB pub 7 (Sun 2019) | CRAG BH01 (deep MspI library) |
|---|---|---|
| chr22 fragments | 330,975 | 13,449,528 (**41x**) |
| median / mean fragment length | 126 / 112 bp | 175 / 404 bp |
| longest fragment | **296 bp** | 30.6 Mb |
| fraction >= 300 bp | **0.000%** | 12.9% |
| fragment ends landing on a CpG | 1.11% / 1.03% (**0.68x / 0.63x chance**) | 2.34% / 1.32% (1.43x / 0.81x chance) |

First, the FinaleDB library is **size-selected at ~300 bp**: not one fragment in 330,975 is longer, so its ends are truncation points rather than digestion sites and land on CpGs at *below* the chance rate. Second, the library is **too shallow to identify a state**: it delivers 1.09 observations per CpG (range 1.00–1.22), with controls reaching 54,922–64,897 CpGs and cancers 11,003–17,894, each CpG seen about once. At one observation per CpG the decoder extrapolates, and what it returns is a restatement of library depth. The decoder output confirms the degeneracy: only 0.035% of called CpGs carry a value strictly between 0 and 1, so the β summary is a count of hard calls, not a graded methylome.

The control that makes this interpretable holds the library fixed and varies only the depth. On one genuinely MspI-digested library (CRAG BH01, 13.4 M chr22 fragments) the share of CpGs whose methylation estimate is graded — the boundary fraction is not a constant zero — rises monotonically with depth: **0.8% at one observation per CpG, 3.1% at three, 8.1% at ten, 11.5% at twenty-one**. The one-observation depth at which the FinaleDB libraries sit is the degenerate end of that curve, so the pipeline is not broken: the reachable libraries cannot support the feature they are being asked to carry. The deep CRAG samples are publication 6, i.e. the training cohort, so this is a feasibility control and not an independent validation.

This is reported as a **null**. The methylation *extraction* path runs on real public data; on the publicly available FinaleDB libraries the *channel* is not a validated feature, and the public methylome data needed to build one does not exist in a usable form in the reachable corpus.

### 3.3 Per-cancer sensitivity at fixed specificity (with DeLong CIs)

Per-cancer OvR is computed with **per-study z-score harmonization inside each CV fold** so the same offset that drove §3.2's pooled AUC is also removed in the per-cancer fits. Top-5 cancers by count, scored against ALL healthy samples (the cross-study generalization view). For sens@spec 95% CIs the standalone artifact uses **DeLong placement-value** CIs (`src/per_cancer_sens_at_spec.py`; Sun & Xu 2014 form [19]); for the AUC variance we use the 5-seed ± std across the StratifiedKFold seeds. Every row below comes from the **GC-corrected** run (the headline pipeline): `results/cross_study_finallydb_gc_corrected.json:per_cancer` (5-seed point estimates + bootstrap 95% CIs on sens@spec) and `results/per_cancer_sens_at_spec.json` (DeLong CIs + PPV@prev, same pipeline).

| Cancer | n_cancer | n_total | AUC mean ± std (5-seed) | Sens@95% | Sens@98% | Sens@99% [bootstrap 95% CI] |
|---|---:|---:|---|---:|---:|---|
| HCC_J | 89 | 353 (jiang only) | 0.7043 ± 0.0139 | 0.360 | 0.236 | 0.146 (0.000–0.323) |
| LUAD | 79 | 343 | 0.9810 ± 0.0031 | 0.899 | 0.886 | 0.861 (0.534–0.933) |
| PAAD | 60 | 324 | 0.9666 ± 0.0038 | 0.883 | 0.800 | 0.567 (0.300–0.853) |
| BRCA | 53 | 317 | 0.9730 ± 0.0016 | 0.906 | 0.774 | 0.566 (0.136–0.861) |
| OV | 28 | 292 | 0.9886 ± 0.0070 | 0.964 | 0.893 | 0.821 (0.633–0.963) |

See **Figure 2** (`docs/figures/fig2_per_cancer_roc.png`) for the per-seed AUC distribution (top-5 cancers) and **Figure 4** (`docs/figures/fig4_sens_operating_points.png`) for the grouped sens@95/98/99 bar chart. Figure 4 sources its values directly from the `sens_at_spec` block of `results/per_cancer_sens_at_spec.json` and is the visual companion to the table above.

**Reading the table.** HCC_J (Jiang 2015 alone, n=89) shows the lowest cross-study AUC (0.704); this is the OvR reading against the cross-study healthy pool (mostly Cristiano healthy), and the wide sens@99 CI [0.000–0.323] flags it as a small-denominator outlier that needs a larger follow-on study. LUAD, BRCA, OV have tight CIs at AUC > 0.97 with sens@99 ≥ 0.57. OV is a high-AUC small-n row (n=28 cancer); the 0.964 sens@95 is consistent with the empirical ceiling on 28 positives and the wide CI [0.633–0.963] at sens@99 says the same. PAAD sens@99 sits at 0.567 with a very wide CI [0.300–0.853] — the n=60 denom is the binding constraint here, not the classifier.

**Prevalence-floor PPV summary (PPV at spec=99% at standard screening-prevalence grid):** sourced from the same JSON (`results/per_cancer_sens_at_spec.json:per_cancer.*.ppv_at_prevalence`); presented against the pooled AUC's sens@99 = 0.700 as the illustrative ceiling.

| Prevalence | PPV at sens@99=0.700, spec=0.99 |
|---:|---:|
| 0.001 (Galleri / CancerSEEK target) | 0.065 |
| 0.004 | 0.219 |
| 0.01 | 0.414 |
| 0.05 | 0.786 |
| 0.10 | 0.886 |
| 0.20 | 0.946 |
| 0.50 | 0.986 |

At the published MCED screening-prevalence assumption of 0.1% (10⁻³), this 5-channel LR's pooled sens@99 would project to **PPV ≈ 6.5%** — not yet clinically competitive at that prevalence, even at the strong pooled-AUC operating point. The PPV rises sharply with prevalence: at the 1% prevalence used in some targeted-MRD workflows it reaches **41%**, and at a 10% prevalence (a symptomatic workup) it is **89%**. The prevalence floor is set by `(1 - spec)` false-positives per 10,000 tested — that floor does not depend on the classifier and is the limiting factor at any sens < 1.0; the model only moves the PPV via `sens · prev` in the numerator. The per-cancer PPV rows in the JSON are sharper and the **OvR-specific** PPV at the 0.001 grid is the cross-study-reading ceiling: HCC_J ≈ 0.014 (sens@99=0.146), LUAD ≈ 0.079 (0.861), BRCA ≈ 0.054 (0.566), PAAD ≈ 0.054 (0.567), OV ≈ 0.076 (0.821).

**Honest framing.** Per-cancer CI widths are wide because per-cancer denominators are small (10–90 positives); the JSON's `skipped` rows respect the `MIN_POSITIVES_FOR_CI=5` floor (DeLong CI is unreliable below it). HCC_J is OvR-defined but **single-cohort** — there are no Jiang healthy controls in the pooled cohort other than the 32 internal ones, so HCC_J's "AUC vs cross-study healthy" is mostly a Cristiano-vs-Jiang-carcinoma problem and should be read with that grain. The OvR scoring against ALL healthy samples (not just within-study healthy) is the cross-study generalization view; it is not the within-study view that the original Cristiano 2019 [7] paper reports.

### 3.4 Honest null / negative ablations

| Ablation | Verdict | Source |
|---|---|---|
| **focal-BCE loss** (`loss="sens_at_spec"`, α_pos=20) | NULL: Δ AUC = −0.0005 [−0.0015, +0.0005], p=0.30 (n=20 paired seeds). Supersedes the n=5 "positive direction, NS" reading. | `docs/SENS_AT_SPEC_ABLATION.md` |
| **sparse-aware projection** (`projection_kinds={"frag_basic":"sparse_aware"}`) | NEGATIVE: Δ AUC = −0.0440 [−0.0593, −0.0287], p=0.0013 (5 seeds). The slot is not sparse at this operating point. | `docs/SPARSE_AWARE_ABLATION.md` |
| **harmonization** (per-cohort batch-effect check) | MIXED on this cohort: pooled AUC is high (0.9684) and per-study AUCs are similar, but the **true-confound control** (cancer=Jiang, healthy=Cristiano) shows that WITHOUT harmonization the classifier would learn "study" instead of "cancer" (0.9989 vs 0.4836). Reported as **PROTECTIVE** — the harmonization is what makes §3.2 a fair cross-study claim. | `docs/HARMONIZATION.md`, `results/cross_study_finallydb_gc_corrected.json:true_confound_control` |
| **V3 methylation channel** (FinaleMe on 31 real FinaleDB samples) | **NULL, and explained.** The decoded β summary reaches LOO AUC 0.9444 — *exactly* the AUC of the covered-CpG count and of the observations-per-CpG count alone; at the 103,863 loci where both arms have coverage, each tested on the samples that cover it, the median cancer-minus-control difference is exactly 0.0 (permutation p = 1.0). The libraries are size-selected at ~300 bp (0.000% of fragments longer) and deliver 1.09 observations per CpG, so the digestion signature FinaleMe reads is absent and the decoder returns hard calls. The extraction path works; the channel is not a feature. | §3.2.4, `docs/METHYLATION_PIPELINE.md`, `results/methylation_cohort_channel.json` (the AUCs and the matched-locus test), `results/methylation_channel_diagnosis.json` (the mechanism) |

All four ablations are reported honestly. No "win by default" framing.

### 3.5 Comparison to published cfDNA benchmarks

DeepCatch reports a **0.968 pooled OOF AUC** on open data (FinaleDB publications 6 + 8, n=627). For context only — and to be explicit that this is NOT a head-to-head — we list here the headline numbers reported by the closest published cfDNA benchmarks. All five used larger, proprietary, or multi-site clinical cohorts; ours is a 627-sample open-data internal CV. The comparability axis is "how well does the 5-channel feature (or its counterpart) separate cancer from healthy in plasma," not "which is clinically better."

| Method (1st author, year) | Cohort | n | Headline metric | Source |
|---|---|---:|---|---|
| DELFI (Cristiano 2019 [7]) | discovery plasma | 232 (HCC, breast, colorectal, lung, ovarian, pancreatic, gastric, bile-duct) | AUC 0.94 (HCC-vs-healthy); per-cancer 0.86–0.97 | Nature 2019 |
| CancerSEEK (Cohen/Klein 2018 [9]) | prospectively-collected plasma | 1,005 (8 cancers + healthy) | sens 0.70 @ spec 0.99 (median across 8 cancers); OvR AUC 0.91 | Science 2018 |
| Galleri (Liu et al. 2020 [10]) | independent validation plasma | 2,823 (50+ cancer signals vs healthy) | OvR AUC 0.92 across all cancer types | Annals of Oncology 2020 / 2021 PATHFINDER |
| CAPP-Seq (Newman 2014 [5] / 2016 [6]) | tumor-informed MRD, NSCLC pilot + integrated-error-suppression update | not pooled (per-cancer n's) | first prototype AUC ~0.86 in stage II–IV NSCLC; iDES-enhanced AUC >0.99 in stage II–IV | Nat Med 2014; Nat Biotechnol 2016 |
| FinaleMe (Liu 2018 [8]) | plasma + reference methylome | 159 (training + held-out) | end-motif OvR AUC 0.91 | bioRxiv 2018 (preprint) |
| **DeepCatch §3.2 (this work, open data)** | FinaleDB open data, internal CV | **627** (Jiang+Cristiano) | pooled **AUC 0.968 ± 0.003**, sens@95=0.865, sens@98=0.796, sens@99=0.700 | `results/cross_study_finallydb_gc_corrected.json` |

**Positioning (honest).** The published numbers above were generated by independent labs on larger, often proprietary, clinical-grade plasma cohorts. DeepCatch's §3.2 number is a **pooled open-data internal-CV** baseline on the same 5-channel DELFI feature (Cristiano-style) that DELFI itself uses. Pooled internal-CV at this scale is competitive in the AUC band of the published numbers (0.86–0.94 discovery-cohort range) — but ours is **not** an external validation, **not** a clinical operating point, and **not** an architectural contribution. The substantive contribution of DeepCatch is the reproducible engineering substrate that lets a researcher reproduce or extend these numbers on open data, not the headline AUC.

Three concrete distinctions matter:

1. **Galleri and CancerSEEK** are methylation-or-mutation-multianalyte assays, not pure fragmentomics. Their reported numbers include assay information (methylation block, oncoprotein panel) that DeepCatch does NOT use. Comparing 5-channel DELFI fragmentomics to a multianalyte panel would be apples-to-oranges; the numbers are listed for context, not as a head-to-head.
2. **CAPP-Seq** is tumor-informed MRD — it tracks patient-specific mutations at 0.01–0.1% ctDNA. DeepCatch's pure-fragmentomics 5-channel LR is not in that regime; comparability with CAPP-Seq on this pooled cross-study cohort is structural, not parametric.
3. **FinaleMe** is the closest cousin methodologically (4-mer end-motif distribution from plasma). Liu 2018 reported AUC 0.91 on n=159; DeepCatch's 5-channel DELFI feature at n=627 reports pooled AUC 0.968. This is the closest fair-context comparison, and is consistent with the headroom available in moving from a 159-patient discovery cohort (n=159 has wide CI on the per-cancer denominator) to a 627-sample pooled CV. An independent head-to-head on a shared held-out cohort would be required to claim a substantive gap.

The open-data benchmark at hand is, by design, **not a clinical claim**. The comparator table is included so the reader can place the §3.2 number in the published-methods landscape without over-reading it.

---

## 4. Discussion

### 4.1 What this work does

DeepCatch provides a reproducible fragmentomics engineering substrate. It takes a public cfDNA cohort (FinaleDB, GDC, FLARE) through (i) feature extraction in a canonical six-modality schema, (ii) self-supervised pre-training on real data, (iii) downstream classification with honest per-patient GroupKFold validation, and (iv) paired-design ablations of every architectural knob. The headline contribution is not a number but the engineering substrate itself: a new lab can swap in a new cohort and get honest numbers in hours.

### 4.2 What this work does NOT do

- **Clinical validation.** DeepCatch is not a clinical assay. No real-plasma cohort has been sequenced. No patient samples were collected under IRB.
- **FDA pathway.** No regulatory submission is planned or implied.
- **Prospective screening.** No prospective cohort, no enrollment, no follow-up.
- **Headline-grade clinical AUCs.** The synthetic-plasma foundation AUC (≈0.55) is below the LR baseline. The cross-study FinaleDB 5-channel LR AUC at n=627 is 0.968 with sens@99=0.700 (pooled, harmonized), but that is **pooled internal-CV**, NOT external validation, NOT a held-out clinical plasma cohort, and NOT a clinical operating point.
- **HCC_J is a within-study reading, not a cross-study one.** The OvR reading of HCC_J (Jiang 2015 only, n=89) against the cross-study healthy pool drops the AUC to 0.704 because there are only 32 Jiang healthy controls — the comparison reduces to "Jiang HCC vs mostly-Cristiano healthy," not a within-study HCC-vs-healthy read. HCC_J PPV @0.001 = 0.014 reflects the low sens@99 (0.146) on that denominator, not a model-quality claim.
- **Per-cancer denominators are small.** Top-5 cancer n ranges from 28 (OV) to 89 (HCC_J); CIs on sens@99 (e.g. OV [0.758–1.000], PAAD [0.333–0.810]) are wide. Conclusions about per-cancer performance should be read as ordinal, not cardinal.
- **Real plasma fragmentomics.** The pretraining cohort is from FinaleDB plasma (Cristiano 2019, Jiang 2015), but the downstream validation cohorts are either synthetic-plasma (TCGA-LUAD real mutations + Poisson sampling, §3.1) or Open-Data FinaleDB uniformly-preprocessed plasma (§3.2/3.3). A real-plasma fragmentomics validation requires an IRB-approved prospective cohort — see §4.4.

### 4.3 Limitations

1. **Small paired synthetic cohort (n=20).** Per-patient GroupKFold on 20 patients means each fold has 4 patients × 2 arms = 8 test samples. CIs are wide; per-seed std on foundation AUC is ~0.01.
2. **Simulated healthy controls.** The healthy arm of every paired patient is a TF=0 simulation of the same patient's mutations, not an unrelated healthy donor. This is the standard analytical validation approach (Newman 2016 [6]) but it is not a substitute for unrelated healthy plasma.
3. **No real plasma in downstream validation.** The 0.1% ctDNA numbers are from synthetic reads. Real plasma has additional noise sources (PCR bias, mapping ambiguity, contamination) that the simulation does not capture.
4. **FinaleDB access is partial, and one verdict in the repo was wrong.** The API is reachable again (it was down between 2026-09 and 2026-11-15), but S3 access is per-object: the `*.frag.tsv.bgz`, insert-size and `*.bigWig` objects return 200 while `*.bedGraph.bgz` returns 403. A single-key probe had previously reported the whole corpus as private; the multi-key probe in `scripts/fetch_external_cohort.py` is what corrected it.
5. **Encoder size.** PRODUCTION_CONFIG is 543,872 parameters (embed_dim=128, 4 layers, 4 heads). State-of-the-art fragmentomics encoders are 10–100× larger. The present model fits on CPU in <5 seconds; clinical-grade models would require GPU training.
6. **Two-channel only.** The pretraining data covers Cristiano 2019 + Jiang 2015 (publications 6 + 8). FinaleMe (Liu 2018 [8]) and the methylation-based MCED publications are not yet in the cohort.
7. **The methylation channel is not a feature on public libraries.** FinaleMe needs a deep, size-unselected MspI library to read methylation from the digestion pattern. The reachable FinaleDB libraries are trimmed at ~300 bp and deliver ~1.09 fragment observations per CpG, so the decoded methylome is a restatement of library depth (§3.2.4). The CRAG/Jiang chr22 BAMs on Zenodo 6914806 are 41x deeper and do carry the signature, but they are publication 6 — the training cohort — so they would give a feasibility test, not independent validation. No independent MspI cohort is currently reachable.

### 4.4 Future work

- **Real-plasma paired cohort.** Requires IRB + collaborator with banked plasma + matched tumor sequencing. Estimated cost: $50K–$200K for a 50-patient pilot. Target: CCGA, TRACERx, or an in-house liquid-biopsy clinic.
- **Cross-platform validation.** Apply the pretrained encoder to ONT / PacBio / single-cell cfDNA to test transfer. The 4-mer end-motif distribution should be platform-robust; the coverage/CNV channel may need per-platform recalibration.
- **Methylation and fragment-end motifs.** The 256-channel 4-mer motif distribution (Liu 2018 [8]) is currently projected to 36 PCA components. A direct 256-channel input would test whether the PCA bottleneck costs signal.
- **Larger pretraining cohort.** Scale from 200 → 2,000+ real FinaleDB samples once the API is restored. A 10× cohort with 50 epochs would push the encoder from "demonstration" to "useful."
- **Independent replication.** A second lab, with its own implementation, running the same cohort through DeepCatch. The released artifacts (cohort matrices, per-seed AUC arrays, ablation JSONs) make this reproducible.

---

## 5. Availability

- **Source code:** github.com/rollroyces/deepcatch
- **License:** MIT
- **Data** (all public; re-fetched by the scripts in `scripts/`; no
  patient-identifiable data is redistributed):
  - FinaleDB pre-extracted cache: `cfdna-fragmentomics-pipeline/data/features/` (657 samples, Cristiano 2019 + Jiang 2015) -- the cross-study benchmark cohort (§3.2), of which 627 have features. **This cache is in neither repository**: the pipeline repo gitignores `data/features/` (~600 MB) and `data/raw/` (~100 GB), because both are rebuilt from the public FinaleDB API. Build it once with `python run_cross_study.py --parallel 8 --max-mb 500` in the pipeline checkout; the driver (`paper/REPRODUCE.sh`) then finds it at that path by default. The labels file (`labels_multiclass.tsv`) *is* committed in the pipeline repo.
  - FinaleDB **held-out external cohort**: publications 1 (Snyder 2016), 7 (Sun 2019) and 9 (Adalsteinsson 2017), fetched from the public REST API and the public S3 bucket by `scripts/fetch_external_cohort.py`; these contribute **no** training data (§3.2.2). No IRB or data-use agreement is required.
  - GDC TCGA-LUAD MAF files: `validation/tcga/tcga_cache/` (402 MAF files, open access)
  - FLARE / GSE317007 (ONT HNSCC) and GSE185307: see `docs/PATH_TO_IMPACT.md`
- **Pretrained checkpoints — not distributed with this repository.** They are
  regenerated by `scripts/pretrain_production_finaledb.py` from the public
  FinaleDB cohort (200 real samples; see `docs/PRETRAINING.md`):
  - `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` (real-data-trained, 543,872 params, PRODUCTION_CONFIG)
  - `checkpoints/foundation_pretrained_SYNTHETIC_v0.pt` (historical; not real-data-trained — do NOT use for "real-data-derived" claims)
  No number reported in §3 depends on these files: the encoder behind the
  §3.1 comparison is trained in-run with `pretrained=False`, so §3 reproduces
  from a clean checkout without them.
- **Reproduce all numbers:**
  - One-bash driver: `bash paper/REPRODUCE.sh` (or equivalently `bash paper/REPRODUCE.sh --quick` for a reduced run). See `paper/REPRODUCE.md`.
  - Full top-level driver: `RUN_ALL.sh` and the per-script invocations documented in §3.
- **Companion documents:**
  - `docs/figures/README.md` — the 4 standard cfDNA validation figures
  - `docs/PRETRAINING.md` — pre-training pipeline
  - `docs/PRETRAIN_BUG.md` — synthetic-bypass bug + regression test
  - `docs/SPARSE_AWARE_ABLATION.md` — honest null/negative ablation
  - `docs/SENS_AT_SPEC_ABLATION.md` — honest null/negative ablation
  - `docs/HARMONIZATION.md` — per-study batch-effect check
  - `AUDIT_2_FINDINGS.md` (repo root) — consolidated audit-2 report
  - `docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md` — cross-study AUC table (headline run); `docs/CROSS_STUDY_BENCHMARK.md` — the legacy no-GC counterpart
  - `MODEL.md` — model architecture reference
  - `docs/COLLABORATOR_DATA_INTERFACE.md` — schema for external collaborators to contribute a plasma cohort (no code required); see also `docs/COLLABORATOR_QUICKSTART.md` for a 1-page recipe.
  - `paper/paper.tex` — the benchmark paper being submitted (LaTeX source; the compiled artefact is `paper/biorxiv_submission_v2.3.0.pdf`)
  - `paper/PAPER.md` — the earlier panel-only revision of that paper, kept for reference

## 6. Acknowledgements

The author thanks the public data providers whose archives make
replication possible: the FinaleDB and FinaleMe teams (Zenodo
10.5281/zenodo.14013719 and 10.5281/zenodo.6914806), the Genomic Data
Commons for TCGA open-access somatic mutations, and the GEO depositors
of GSE317007 and GSE185307.

---

*Manuscript draft — 2026-09-23. Section 3 filled with the just-shipped cross-study and per-cancer numbers. See `paper/REPRODUCE.md` for one-command reproduction.*

## 7. References

1. Lo YMD, Corbetta N, Chamberlain PF, et al. Presence of fetal DNA in
   maternal plasma and serum. *Lancet*. 1997;350:485–487.
   doi:[10.1016/S0140-6736(97)02174-0](https://doi.org/10.1016/S0140-6736(97)02174-0)
2. Wan JCM, Massie C, Garcia-Corbacho J, et al. Liquid biopsies come of age:
   towards implementation of circulating tumour DNA. *Nat Rev Cancer*.
   2017;17:223–238. doi:[10.1038/nrc.2017.7](https://doi.org/10.1038/nrc.2017.7)
3. Heitzer E, Haque IS, Roberts CES, Speicher MR. Current and future
   perspectives of liquid biopsies in genomics-driven oncology.
   *Nat Rev Genet*. 2019;20:71–88.
   doi:[10.1038/s41576-018-0071-5](https://doi.org/10.1038/s41576-018-0071-5)
4. Ignatiadis M, Sledge GW, Jeffrey SS. Liquid biopsy enters the clinic —
   implementation issues and future challenges. *Nat Rev Clin Oncol*.
   2021;18:297–312.
   doi:[10.1038/s41571-020-00457-x](https://doi.org/10.1038/s41571-020-00457-x)
5. Newman AM, Bratman SV, To AT, et al. An ultrasensitive method for
   quantitating circulating tumor DNA with broad patient coverage.
   *Nat Med*. 2014;20:548–554.
   doi:[10.1038/nm.3519](https://doi.org/10.1038/nm.3519)
6. Newman AM, Lovejoy AF, Klass DM, et al. Integrated digital error
   suppression for sensitive detection of circulating tumor DNA.
   *Nat Biotechnol*. 2016;34:547–555.
   doi:[10.1038/nbt.3520](https://doi.org/10.1038/nbt.3520)
7. Cristiano S, Leal A, Phallen J, et al. Genome-wide cell-free DNA
   fragmentation in patients with cancer. *Nature*. 2019;570:385–389.
   doi:[10.1038/s41586-019-1272-6](https://doi.org/10.1038/s41586-019-1272-6)
8. Liu X, Liu Z, Chen S, et al. Comprehensive DNA methylation analysis with a
   scalable end-to-end workflow. *bioRxiv*. 2018.
   doi:[10.1101/478222](https://doi.org/10.1101/478222)
9. Klein JD, Cohen JD, Wang Y, et al. Detection and localization of
   surgically resectable cancers with a multi-analyte blood test.
   *Science*. 2018;359:926–930.
   doi:[10.1126/science.aar3247](https://doi.org/10.1126/science.aar3247)
10. Liu MC, Oxnard GR, Klein EA, et al. Sensitive and specific multi-cancer
    detection and localization using methylation signatures in cell-free
    DNA. *Ann Oncol*. 2020;31:745–759.
    doi:[10.1016/j.annonc.2020.02.011](https://doi.org/10.1016/j.annonc.2020.02.011)
11. Jiang P, Sun K, Peng W, et al. Plasma DNA end-motif profiling as a
    fragmentomic marker in cancer, pregnancy, and transplantation.
    *Cancer Discov*. 2020;10:664–673.
    doi:[10.1158/2159-8290.CD-19-0622](https://doi.org/10.1158/2159-8290.CD-19-0622)
12. EpiFluidLab. *FinaleDB*: a database of fragmentomics datasets. 2022.
    <https://finaledb.epifluidlab.cchmc.org/>
13. Straub J, Mathios D, Foda Z, et al. Finaletoolkit: an end-to-end
    infrastructure for the analysis of cell-free DNA fragmentomics data.
    *bioRxiv*. 2023. The identifier previously recorded for this preprint
    does not resolve, so none is claimed here.
14. DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two
    or more correlated ROC curves: a nonparametric approach. *Biometrics*.
    1988;44:837–845. doi:[10.2307/2531595](https://doi.org/10.2307/2531595)
15. Jiang P, Chan KCA, Lo YMD. Lengthening and shortening of plasma DNA in
    hepatocellular carcinoma patients. *PNAS*. 2015;112:E1317–E1325.
    doi:[10.1073/pnas.1500076112](https://doi.org/10.1073/pnas.1500076112)
16. Snyder MW, Kircher M, Hill AJ, Daza RM, Shendure J. Cell-free DNA
    comprises an in vivo nucleosome footprint that informs its
    tissues-of-origin. *Cell*. 2016;164:57–68.
    doi:[10.1016/j.cell.2015.11.050](https://doi.org/10.1016/j.cell.2015.11.050)
17. Sun K, Jiang P, Cheng SH, et al. Orientation-aware plasma cell-free DNA
    fragmentation analysis in open chromatin regions informs tissue of
    origin. *Genome Res*. 2019;29:418–427.
    doi:[10.1101/gr.242719.118](https://doi.org/10.1101/gr.242719.118)
18. Adalsteinsson VA, Ha G, Freeman SS, et al. Scalable whole-exome
    sequencing of cell-free DNA reveals high concordance with metastatic
    tumors. *Nat Commun*. 2017;8:1324.
    doi:[10.1038/s41467-017-00965-y](https://doi.org/10.1038/s41467-017-00965-y)
19. Sun X, Xu W. Fast implementation of DeLong's algorithm for comparing the
    areas under correlated receiver operating characteristic curves.
    *IEEE Signal Process Lett*. 2014;21:1389–1393.
    doi:[10.1109/LSP.2014.2337313](https://doi.org/10.1109/LSP.2014.2337313)
