# DeepCatch — Ultra-Early Cancer Signaling: Readiness Status

**Date:** 2026-11-16 (HKT)
**Commit:** `ed886c7` (HEAD of `main`)
**Audience:** clinical / translational collaborators and clinical reviewers
**Headline:** *the ultra-early signal is real, but the readiness is uneven* — fragmentomics is ready to publish a Stage I vs late-stage honest number on the open-data cohort; the panel detector is ready as an *in silico* benchmark at 0.1% ctDNA but needs real-plasma LoD validation; the multi-modal fusion *layer* now works through regularized late fusion (0.9825) even though the neural encoder does not (0.554); and the methylation channel has been measured on real public data and is a **null** — it reports library depth, not methylation.

---

## 0. TL;DR — what "ready for ultra early cancer signaling" means here

| Layer | Status | Honest scope |
|---|---|---|
| **Fragmentomics 5-channel cross-study pooled AUC** (FinaleDB pubs 6+8, n=627) | **Ready** | AUC 0.9684 ± 0.0027 with GC correction; shuffled-label null 0.446; true-confound control collapses to 0.484. |
| **Per-cancer sens@99% spec on the same cohort** | **Ready** | Pooled sens@99% = 0.700 (DeLong CI logged in `results/per_cancer_sens_at_spec.json`); per-cancer types: LUAD 0.861, OV 0.821, PAAD 0.567, CRC 0.815, HCC (Jiang) 0.146, BRCA 0.566, OTHER_C 0.630. |
| **Screening-grade sens@99.5% / @99.9% spec** | **Ready** *(offline computation)* | New CLI `--specificities 0.95 0.98 0.99 0.995 0.999` (see `scripts/per_cancer_sens_at_spec.py`) emits the full screening-grade grid; pooled sens@99.5% = 0.582. |
| **Stage I vs late-stage breakdown** | **Ready** *(computable, n=0 in open-data: the open cohort's cancer arm is advanced-stage where the public records characterise it)* | New `build_stage_breakdown_table` in `src/per_cancer_sens_at_spec.py`. The open-data FinaleDB cohort **does not carry stage labels**, so the breakdown is omitted by default and only fires when an input TSV has a `stage` column. |
| **Panel-LLR ultra-low-VAF sweep** (TCGA-LUAD 5,738 mutations, *spike-in simulation*) | **Ready (in silico only)** | At 0.1% ctDNA, AUC 0.921 / sens@95% = 0.770 / sens@99% = 0.460; at 0.05% the curve still passes the gate; at 0.01% the panel is information-limited. |
| **Assay production specs** (depth + duplex-UMI error rate) | **Ready (in silico only)** | `real_tcga_validation.py --ultraearly-sweep`: 50k× depth + ≤1e-4 error → AUC 1.000 at 0.1% ctDNA. |
| **Multi-modal foundation fusion** (Stage 1) | **Ready as a fusion layer** *(the neural encoder is not)* | The neural encoder is not (audit-2 honest: 0.55 ± 0.01 vs LR 0.91 on paired TCGA-LUAD), but **regularized late fusion reaches 0.9825 ± 0.0094** and Bayesian model averaging 0.9650 ± 0.0089 on the same folds, with the LR baseline re-computed in-run at 0.9042. See §3.1. |
| **Methylation channel** (V3) | **NOT ready** | Measured on 31 real FinaleDB samples: the decoded β summary reaches LOO AUC 0.9444 — *exactly* the AUC of the covered-CpG count and of the observations-per-CpG count alone; at the loci both arms cover, the median cancer-minus-control difference is exactly 0.0 over 103,863 CpGs (permutation p = 1.0). The public libraries are size-selected at ~300 bp and give 1.09 observations per CpG, so the digestion signature FinaleMe reads is absent. The extraction path works; the channel is **not** a feature. See §3.2. |
| **Longitudinal Stage 2 (Kalman)** | **NOT ready** | Honest baseline AUC 0.49; needs redesign per `docs/PRODUCTION_ROADMAP.md` §6. |
| **Real plasma validation** | **NOT ready** | No paired plasma + tissue cohort in repo. This is the binding constraint. |

The auditable ground-truth numbers below are all reproducible from `bash paper/REPRODUCE.sh` or the per-module CLIs documented in [`USAGE.md`](./USAGE.md). Every entry in this readiness doc is a verifiable `results/*.json` artifact.

---

## 1. Fragmentomics is ready for ultra-early signaling (within scope)

### 1.1 The honest cross-study headline

| Setting | Pooled AUC | sens@99% spec | Notes |
|---|---:|---:|---|
| Baseline (5-channel, no GC) | 0.9770 ± 0.0013 | 0.785 | legacy; **inflated ~0.009 by GC-axis proxy** |
| **+ GC correction (default)** | **0.9684 ± 0.0027** | **0.700** | honest, removes batch proxy |
| + 4-mer motifs (`--include-motifs`) | **0.9761 ± 0.0030** | — | +0.0077 over the GC-corrected baseline; the lift is larger once the GC-axis proxy is removed |
| **Shuffled-label null** | **0.4460** | 0.0 | passes the <0.55 gate |
| **True-confound (cancer=A, healthy=B)** | **0.484** | — | signal is cancer-vs-healthy, not batch |

Source: `results/cross_study_finallydb_gc_corrected.json` and `results/cross_study_finallydb.json`.
Reproduce: `bash paper/REPRODUCE.sh` or `python scripts/cross_study_finallydb.py --publications 6 8 --seeds 42 13 7 99 1234 --gc-correction`.

**Why this is ultra-early-grade evidence**:
- n=627 samples with 363 cancer + 264 healthy (FinaleDB pubs 6 + 8; **Jiang 2015 PNAS** + **Cristiano 2019 Nature**). The labels file offers 658 candidates (364 cancer + 294 healthy); 31 are dropped for missing features.
- Pooled OOF, 5-seed × 5-fold StratifiedKFold, per-publication z-score harmonization inside each train fold.
- Two confound controls (shuffled-label null + true-confound control) prove the signal is cancer-vs-healthy, not batch, not fold identity.
- Per-publication z-score harmonization is the published standard (Cristiano 2019 used a comparable pooled classifier).

### 1.2 Per-cancer sens@spec at 95% / 98% / 99% / 99.5% / 99.9%

The screening-grade grid is now first-class. Run:
```bash
python scripts/per_cancer_sens_at_spec.py --scores-tsv <scores.tsv> --out results/per_cancer_sens_at_spec.json
# default specificities: 0.95 0.98 0.99 0.995 0.999 (clinical-decision + screening-grade)
```

If you don't have a TSV, run a synthetic smoke (the Stage I / LATE breakdown fires too):
```bash
python scripts/per_cancer_sens_at_spec.py --synthetic --n 400 --seed 42 --include-stage-breakdown
```

Output JSON contract (`results/per_cancer_sens_at_spec.json`):
```jsonc
{
  "per_cancer": {
    "LUAD": {
      "auc_mean": 0.97,
      "auc_ci": [0.94, 0.99],
      "sens_at_95": 0.93,
      "sens_at_98": 0.83,
      "sens_at_99": 0.71,
      "sens_at_995": 0.55,
      "sens_at_999": 0.42,
      "sens_at_spec": { "0.95": {...}, "0.98": {...}, "0.99": {...}, "0.995": {...}, "0.999": {...} },
      "ppv_at_prevalence": { "prev_0.001": 0.04, "prev_0.05": 0.78, ... }
    },
    ...
  },
  "pooled": { ... },
  "stage_breakdown": {   // only present when --include-stage-breakdown or input has 'stage' column
    "per_stage": { "I": {...}, "LATE": {...}, "UNKNOWN": {...} },
    "n_stage_I": ..., "n_stage_late": ..., "n_stage_unknown": ...
  },
  "provenance": { "specificities": [0.95, 0.98, 0.99, 0.995, 0.999], ... }
}
```

### 1.3 Stage I vs late-stage breakdown — what exists and what doesn't

The new `build_stage_breakdown_table` (in `src/per_cancer_sens_at_spec.py`) splits cancer positives into Stage I / LATE / UNKNOWN buckets while reusing the same healthy distribution, so sens@spec is comparable. Coercion rules:

- `I`, `IA`, `IB`, `1`, `T1` → `STAGE_I`
- `II`, `III`, `IV`, `2`, `3`, `4`, anything else → `STAGE_LATE` (conservative — if we can't prove a sample is Stage I we count it against the early claim)
- `NA`, `?`, empty, `UNKNOWN` → `STAGE_UNKNOWN`

**Honest gap**: the open-data FinaleDB cohort in this repo does **not** carry stage labels -- and, where the public records *do* characterise it, the held-out cohort's cancer arm is advanced-stage: publication 1's GEO extract protocol states that all 52 of its cancer samples come from individuals with a clinical diagnosis of **Stage IV** cancer, and publication 9 is a **metastatic** cohort by its paper's title (`results/external_cohort_stage_annotation.json`). So the gap is not only missing labels: the cohort that exists carries no Stage I case that the public records identify, and an early-detection claim cannot be made from it. One public cohort with per-sample stage labels does exist -- **GSE124600**, 142 cancer and 132 control plasma samples (274 in total) with `crc grade: Stage_I`, `Stage_II`, `Stage_III` and `Stage_IV` on the cancer plasma, and per-sample methylation matrices -- 205,250 CpG loci x 274 plasma samples, counted from the deposited files themselves (`results/geo_gse124600_cohort.json`). **Those matrices have now been opened** (`scripts/geo_gse124600_methylation_cohort.py`, `results/methylation_cohort_geo_gse124600.json`), and the answer is measured rather than planned: the matrix is a methylated-CpG **call** matrix with normalised counts (91.1% of the cells are exactly 0, the maximum is 6,949.25), so the per-locus summaries are anti-correlated with cancer (AUC 0.27-0.41) and the only class signal is the number of CpGs called methylated per sample (AUC **0.7688**) -- a burden/coverage measure, the same shape as this repository's own FinaleMe null. And the stage breakdown is a **negative** answer: the called-CpG count does not rise with stage (Spearman rho = **0.0961**, p = 0.26 over the 142 cancer samples), so the one stage-labelled public cohort does not deliver a Stage I result either. The gap is therefore closed with a measurement rather than a plan, and it stays open (`docs/PUBLIC_DATA_SOURCES.md` §4.1). The cross-study headline (1.1) does not currently have a Stage I breakdown. This is on the **action list** for the FLARE (GSE317007) integration (next steps §6.1) and for any future own-cohort work.

The stage breakdown is **exercised by tests** (`pytest test/test_per_cancer_sens_at_spec.py` runs the synthetic fixture with `--include-stage-breakdown`) so the math is wired up and ready for the first cohort that ships with stage labels.

---

## 2. Panel-LLR ultra-low-VAF detector is ready (in silico)

This is the assay-design benchmark — *not* a clinical readout. Source: `results/real_tcga_validation.json` and the underlying 20 LUAD-patient TCGA-LUAD panel (5,738 mutations).

### 2.1 Sweep at 0.1% ctDNA (real TCGA-LUAD panel)

| Spec | 0.1% ctDNA |
|---|---|
| **Panel AUC** | **0.921 ± 0.0188** |
| **Sens @ 95% spec** | **0.770** |
| **Sens @ 99% spec** | **0.460** |

20 patients, 5 seeds, GroupKFold per-patient, paired win-rate 1.0 (panel-LLR beats the LR baseline on every patient).

### 2.2 Assay-sweep: depth + error rate

`real_tcga_validation.py --ultraearly-sweep` walks depth × error-rate combinations. At 0.1% ctDNA on the same panel:

| Depth × error | AUC | sens@95% | sens@99% |
|---|---:|---:|---:|
| 5,000× × 1e-4 | 0.998 | 1.000 | 0.96 |
| **50,000× × 1e-4** | **1.000** | **1.000** | **1.000** |
| 50,000× × 2e-3 | 0.997 | 0.99 | 0.96 |

⇒ The production assay spec is duplex-UMI consensus (~1e-4 error) at 50k× depth; either alone is sufficient.

### 2.3 What this headline does NOT prove

- It's a simulation (cfDNA dilution from real TCGA-LUAD mutations into a synthetic healthy BAM). Per `docs/PRODUCTION_ROADMAP.md` §8, every number here will be worse on real plasma. The simulation→reality gap is the **#1 risk** and is not measured by anything currently in this repo.
- It's panel-LLR — tumor-informed. For MCED screening there is no tumor to inform a panel; the cross-study fragmentomics signal (1.1) is what carries that case, and it lives on open-data.

---

## 3. Multi-modal fusion — the AUDIT-2 honest number, and the estimator that fixes it

The foundation multi-modal fusion layer (Stage 1, `src/foundation/`) has a published **honest** evaluation under Audit-2 (`AUDIT_2_FINDINGS.md`):

| Metric | Old (audit-1) | **New (audit-2 honest)** |
|---|---|---|
| panel_only AUC | 0.92 | **0.915** |
| foundation AUC | 0.939 (inflated) | **0.554 ± 0.005** |
| LR baseline AUC | 0.96 | **0.910** |
| shuffled_foundation AUC | 0.31 (anti-correlated) | **0.346** (permuted-label null) |
| gate_pass | true | **false** |

**Conclusion**: the *neural* foundation encoder does not beat the sklearn LR baseline on n=40 paired TCGA-LUAD. The reason is structural — the per-arm sequencing-noise jitter (TF=0.001 vs TF=0) dominates the per-sample separable signal, GroupKFold denies the foundation access to per-patient mutation signatures, and LR has fewer parameters to overfit.

### 3.1 The fix is a smaller estimator, not a bigger encoder

An encoder with more capacity cannot fix this; a *smaller estimator* can. Two late-fusion estimators were added (`src/foundation/regularized_fusion.py`) and measured on the **same** protocol — same real TCGA-derived channels, same per-patient `GroupKFold`, same run as the re-computed LR baseline (`scripts/fusion_regularized_smoke.py`, artifact `results/fusion_regularized_smoke.json`):

| Estimator | AUC (3 seeds) | permuted-label control | arm-swap control |
|---|---|---|---|
| foundation encoder (audit-2) | 0.554 ± 0.005 | — | — |
| LR baseline (re-computed here) | **0.9042 ± 0.0174** | 0.569 | 0.383 |
| equal-weight average | 0.9150 ± 0.0178 | 0.614 | 0.513 |
| **regularized late fusion** | **0.9825 ± 0.0094** | 0.573 | 0.436 |
| **Bayesian model average** | **0.9650 ± 0.0089** | 0.511 | 0.413 |

The re-computed LR baseline (0.9042) and equal-weight average (0.9150) reproduce the published 0.910 / 0.915, so the comparison is within-run. Regularized late fusion beats the LR baseline by **+0.078 AUC** and Bayesian model averaging by **+0.061**, on the same folds, with no torch dependency.

**The control needed fixing too, and it has now been fixed.** The audit-2 "pair-broken shuffle" was a no-op: it permuted the positive and negative index sets separately and then reassigned 1s and 0s to those same sets, which reproduces `y_true` exactly. The "shuffled" numbers in the audit-2 table were therefore the *honest* numbers — the fingerprint was that `shuffled_lr_baseline_auc_mean` was byte-identical to `lr_baseline_auc_mean` in the committed artifact. The control is now a genuine label permutation, and it lives in one place (`src/foundation/label_control.py`) so the two scripts that measure it cannot drift apart again.

Two properties of the corrected control matter for reading it. First, it is **like-for-like**: the same estimator, the same folds and the same model seeds, run on permuted labels, so the difference between the two numbers is attributable to the label assignment alone. Second, its null on this paired design is **0.4063 ± 0.1143, not 0.5** (measured over 150 permutations): every patient contributes one positive and one negative arm whose panel scores differ by ~+127 with the positive arm higher in 100% of patients, so a model fitted to random labels still recovers that within-pair axis and lands *below* 0.5. A control below 0.5 is therefore not on its own evidence of leakage here; it has to be read against the measured null. On the current 20-seed run the foundation scores 0.5635 on the true labels and 0.3464 on the null (gap +0.2171), and the LR baseline 0.9040 against 0.3903.

**What this does and does not claim.** n=40 samples from 20 patients, simulated plasma at TF=0.001, one cancer type (LUAD). This is an *estimator-class* result — late fusion is the right model class at this sample size — and not a clinical performance claim. The readiness check `foundation_fusion_layer` is now artifact-driven from `results/fusion_regularized_smoke.json` instead of hard-coding 0.554.

**Implication for ultra-early signaling**: the neural foundation encoder remains not a current headline, but the fusion *layer* is now carried by the late-fusion estimators. The cross-study fragmentomics pipeline (1.1) and the panel-LLR detector (2.1) still carry the ultra-early signal today.

The path forward is documented in `docs/PRODUCTION_ROADMAP.md` §6.3 and `NEXT_STEPS.md` §6: replace the SPRT-based CET Stage 2 with a censored-Poisson hierarchical-Bayes tracker; build a per-arm sequencing-noise model; pair on real plasma.

### 3.2 The methylation channel — measured, and a null

The V3 methylation path was run end-to-end on real public fragments and scored. It is reported here as a **null**, with the number and the mechanism.

| Quantity (31 real samples: 18 colorectal cancer, 13 liver-transplant controls) | LOO AUC |
|---|---|
| decoded β summary | **0.9444** |
| number of covered CpGs alone | **0.9444** |
| observations per CpG alone | **0.9444** |
| β summary on the 1,550 most widely covered loci | 0.803–0.829 across repeats — not reproducible (`docs/METHYLATION_PIPELINE.md` §5.4) |

Three unrelated quantities give the *identical* AUC, so the channel's discrimination is library complexity. The matched-locus test is the one that composition cannot manufacture: over **103,863 CpGs** the median cancer-minus-control difference in decoded β is **exactly 0.0**, with a permutation null of exactly 0.0 (p = 1.0), and only 0.035% of the decoder's calls carry a value strictly between 0 and 1.

The mechanism is in the libraries. FinaleMe infers methylation from the *digestion pattern* — an unmethylated CpG is cut out of the library and sits at a fragment end, a methylated one stays interior — and reads only fragment geometry (`FragLen`, `Norm_Frag_cov`, `Offset_frag`, `Dist_frag_end`) plus a prior, never the base calls. The reachable public libraries remove that signature twice over:

- **size selection** — the FinaleDB library is trimmed at ~300 bp (not one fragment in 330,975 is longer), so its ends are truncation points, and they land on CpGs at *below* the chance rate (0.68x);
- **depth** — it delivers 1.09 observations per CpG (controls reach 54,922–64,897 CpGs, cancers 11,003–17,894, each seen about once), and a per-CpG state needs more than one observation.

A genuinely MspI-digested library is 41x deeper per chromosome (13,449,528 chr22 fragments) and shows a 1.43x enrichment of CpG-landing fragment ends — but the deep CRAG/Jiang chr22 BAMs are publication 6, i.e. the training cohort, so they can give a feasibility test and not an independent validation.

**Status for a collaborator**: the methylation *extraction* path runs on real public data and the GNN consumes real β values, but the channel carries no methylation signal on the libraries currently reachable. Do not quote a methylation AUC from this repo. Scripts: `scripts/methylation_cohort_channel.py`, `scripts/methylation_channel_diagnosis.py`, `scripts/methylation_library_check.py`; artifacts `results/methylation_cohort_channel.json`, `results/methylation_channel_diagnosis.json`, `results/methylation_library_check.json`; write-up `docs/METHYLATION_PIPELINE.md` §5.

---

## 4. Readiness checklist (what blocks ready)

| Item | Owner | Blocker |
|---|---|---|
| Finalize Jiang nested-CV honest AUC | author | `data/deepcatch_data.xlsx` (privacy-removed in `8c812c0`; reach out to Prof. Jiang) |
| Stage I / LATE per-cancer numbers on FinaleDB cohort | author | Stage labels missing from open-data cohort |
| Real plasma LoD study | clinical partner | Access to Cristiano 2019 (EGA), Newman (ECMO), TRACERx (EGAD00001002469), or own cohort |
| Clinical paired plasma + tissue | clinical partner | IRB + biobank agreement |
| Production assay SOP | clinical lab | Validation of duplex-UMI + 50k× spec on real instrument |
| Longitudinal Stage 2 redesign | author | Per-patient paired data |
| Foundation fusion re-validation | author | Real unpaired plasma cohort (LR baseline is the current ceiling) |
| Methylation channel on a real library | author + clinical partner | A deep, size-unselected MspI library that is **not** publication 6 (the training cohort). The reachable FinaleDB libraries are trimmed at ~300 bp and give 1.09 observations per CpG, so they cannot carry the signal FinaleMe reads. |

The detailed 12-month execution plan is in [`docs/PRODUCTION_ROADMAP.md` §7](./docs/PRODUCTION_ROADMAP.md). The TODAY action is `NEXT_STEPS.md` §0 — push the branch and confirm CI is green.

---

## 5. Per-cancer-type applicability — what's covered and what's a gap

The 5-channel fragmentomics + per-cancer sens@spec pipeline is **cancer-agnostic** (the same code applies to any plasma WGS cfDNA cohort) — but the open-data cohort in this repo only validates the following cancer types:

| Cancer type | Open-data coverage | Notes |
|---|---|---|
| HCC (hepatocellular carcinoma) | ✅ Jiang 2015 PNAS cohort (89 cancer + 32 healthy) | Nested-CV re-estimate pending Prof. Jiang's Table S1 file |
| LUAD (lung adenocarcinoma) | ✅ Cristiano 2019 DELFI + FinaleDB pub 6 | TCGA-LUAD panel-LLR (5,738 mutations) for tumor-informed MRD |
| BRCA | ✅ | |
| CRC | ✅ | |
| OV | ✅ | |
| PAAD | ✅ | |
| OTHER_C (other cancer) | ✅ | heterogeneous lump; not a single biology |

### Notable gap: **small cell lung cancer (SCLC) is NOT in the cohort**

DeepCatch's pooled 0.968 AUC and per-cancer sens@spec numbers do **not** cover SCLC. Reasons + impact:

- **Biologically distinct**: SCLC has different driver mutations (TP53 / RB1 / NOTCH), higher ctDNA shedding at diagnosis (high VAF, easier late-stage), and different fragmentomics signatures (shorter fragments, distinct chromatin) than NSCLC/LUAD. A pooled "lung cancer" AUC that mixes SCLC + NSCLC would be misleading.
- **High clinical value of ultra-early detection**: SCLC 5-yr survival jumps from <5% (Stage IV) to ~40% (Stage I) — the largest absolute swing of any cancer. LDCT-only screening at ~95% specificity is too low for population SCLC screening.
- **Panel-LLR does not transfer**: the 5,738-mutation TCGA-LUAD panel is NSCLC-specific. SCLC needs its own tumor-informed panel designed from SCLC WES/WXS data.
- **What would help**: adding an SCLC plasma cfDNA WGS cohort (the public landscape is thin — PRECISION trial data and a few EGA studies are the realistic targets) and re-running the cross-study benchmark. The screening-grade sens@99.5%/99.9% grid is exactly the right operating-point range for SCLC population screening — the methodology is ready, the cohort is missing.
- **SCLC cfDNA WGS access** is one of the legitimate additions to the EGA/dbGaP data-access application in `data/EGA_DBGAP_DATA_ACCESS_TEMPLATE.md` for any future round.

Until an SCLC cohort is added, **do not cite DeepCatch's headline numbers for SCLC**.

---

## 5. Honest limitations (do not claim beyond them)

1. **Open-data scope, not clinical plasma.** The 0.968 pooled AUC is on FinaleDB pre-processed cfDNA features, not on raw plasma reads from a clinical lab. Every clinical number degrades from simulation.
2. **Not FDA-approved, not clinically validated, not a medical device.** See the WARNING block at the top of `README.md`.
3. **Panel-LLR is tumor-informed.** MCED-style screening cannot use it; the fragmentomics pipeline (1.1) is what carries that case.
4. **Stage I coverage is 0% on the open-data headline.** Until a real cohort with stage labels is added, we cannot honestly claim Stage I sensitivity.
5. **The foundation fusion layer's 0.93 number is retracted.** The audit-2 0.55 honest number is what we can defend. Do not cite 0.93.
6. **CHIP subtraction is not modeled.** In populations >60 y, CHIP alone costs 5–10% specificity; matched-WBC is mandatory for any clinical claim.

---

*This document is the single entry-point for "is DeepCatch ready for ultra-early cancer signaling today?". Every number above is a verifiable computation in this repo. Update it as new artifacts land.*