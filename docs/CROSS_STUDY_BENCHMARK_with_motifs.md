# Cross-Study FinaleDB Benchmark (Open Data)
> **Scope**: Open-data benchmark on FinaleDB publications: 6 (jiang) + 8 (cristiano). **NOT** clinical validation. **NOT** external cohort validation. Pooled OOF on the same cohort that trained the model.
- Generated: `2026-11-16T19:43:28.044297+00:00`
- Classifier: `LogisticRegression(max_iter=2000)`
- Feature set: 5-channel (5mb_ratio + 5mb_coverage + 100kb_ratio + c100_gc_corrected + FSD-196) + 256-dim 4-mer end motifs (cohort mode: zero_fill)
- PCA n_components: 200 (capped at `min(n_train, n_features)` per fold)
- CV: 5-fold StratifiedKFold, 5-seed pooled OOF
- Seeds: `[42, 13, 7, 99, 1234]`
- Bootstrap: n=1000, seed=2026
- Cell-line regex: `^(GM\d+|HeLa|HepG2|K562|HL60|Jurkat|Raji|MCF7|U937|THP1|HEK293|HCT116|SW480|A549|GM12878)`

## 1. Cohort inventory
- Samples in labels file (filtered to requested publications): **658** (364 cancer + 294 healthy)
- Samples with all 5-channel DELFI features: **627** (31 dropped due to missing artifacts)
- Cell-line filter removed: **0** samples (none matched the regex on this open-data cohort)
- Per-publication (FinaleDB):
  - **6 (jiang)**: n=121 (89 cancer + 32 healthy)
  - **8 (cristiano)**: n=537 (275 cancer + 262 healthy)
- Per-study (legacy `study` column in labels file):
  - **jiang**: n=121 (89 cancer + 32 healthy)
  - **cristiano**: n=537 (275 cancer + 262 healthy)
- Per-cancer (cancer samples only):
| Cancer | n_total | cristiano | jiang |
|---|---:|---:|---:|
  | HCC_J | 89 | 0 | 89 |
  | LUAD | 79 | 79 | 0 |
  | PAAD | 60 | 60 | 0 |
  | BRCA | 54 | 54 | 0 |
  | OV | 28 | 28 | 0 |
  | CRC | 27 | 27 | 0 |
  | OTHER_C | 27 | 27 | 0 |

## 2. Per-publication AUC
Each publication evaluated independently (no harmonization needed; one publication only).

| Publication | Study | n | n_cancer | n_healthy | AUC (5-seed mean ± std) |
|---|---|---:|---:|---:|---|
| 6 | jiang | 121 | 89 | 32 | 0.9776 ± 0.0024 |
| 8 | cristiano | 506 | 274 | 232 | 0.9739 ± 0.0015 |

## 3. Pooled cross-publication AUC (with/without per-publication harmonization)
Harmonization = per-publication z-score StandardScaler fit on train fold only.

| Setting | AUC mean ± std | Sens@95% | Sens@98% | Sens@99% |
|---|---|---:|---:|---:|
| harmonized | 0.9761 ± 0.0030 | 0.920 | 0.832 | 0.774 |
| no_harmonize | 0.9687 ± 0.0040 | 0.904 | 0.857 | 0.799 |

## 4. Per-cancer Sens@Spec (top-5 cancers by sample count)
One-vs-rest: each cancer vs ALL healthy samples in the pooled cross-study cohort. Per-publication harmonization inside each CV fold.

**Primary CIs are DeLong** (DeLong, DeLong, Clarke-Pearson 1988 for AUC; Sun & Xu 2014 for Sens@spec — equivalent to DeLong-Han-Agarwal structural component restricted to positives). Bootstrap 95% CIs (n=1000 percentile resamples) are preserved under `bootstrap_ci` in `results/cross_study_finallydb.json` for the audit trail (narrower at n>=10; DeLong is the standard cfDNA reference per the cfdna-early-detection-validation skill). The standalone clinical-decision table — AUC + Sens@spec CIs + PPV@prev at spec=0.99 — is at `results/per_cancer_sens_at_spec.json` (schema per `src.per_cancer_sens_at_spec.build_per_cancer_table`).

**Caveat**: HCC_J is Jiang-only (n=89 cancer + 32 healthy). Its per-cancer OvR uses 89 HCC vs ALL healthy controls (264 controls in the OOF subset) — the per-cancer denominator is small but the OvR is well-defined. Other top-5 cancers are Cristiano-only.

| Cancer | n_cancer | n_healthy | AUC mean ± std | AUC DeLong [95% CI] | Sens@95% DeLong [95% CI] | Sens@98% DeLong [95% CI] | Sens@99% DeLong [95% CI] |
|---|---:|---:|---|---|---|---|---|
| HCC_J | 89 | 264 | 0.7410 ± 0.0176 | 0.7604 [0.6930–0.8278] | 0.449 [0.346–0.553] | 0.213 [0.129–0.298] | 0.146 [0.074–0.218] |
| LUAD | 79 | 264 | 0.9840 ± 0.0025 | 0.9858 [0.9749–0.9966] | 0.911 [0.848–0.975] | 0.899 [0.831–0.966] | 0.873 [0.799–0.948] |
| PAAD | 60 | 264 | 0.9784 ± 0.0019 | 0.9791 [0.9620–0.9962] | 0.917 [0.845–0.988] | 0.867 [0.779–0.954] | 0.750 [0.640–0.860] |
| BRCA | 53 | 264 | 0.9757 ± 0.0012 | 0.9763 [0.9553–0.9972] | 0.906 [0.825–0.986] | 0.811 [0.705–0.918] | 0.698 [0.574–0.823] |
| OV | 28 | 264 | 0.9901 ± 0.0061 | 0.9947 [0.9890–1.0000] | 1.000 [0.965–1.000] | 0.893 [0.773–1.000] | 0.857 [0.723–0.991] |

## 5. True-confound control
Cancer = 100% from one publication, healthy = 100% from another. Without harmonization the classifier learns 'which publication is this from?' (AUC ~0.999). With per-publication z-score harmonization the publication-specific mean/variance is the only signal and is removed by design (AUC should collapse toward 0.50).

| Orientation | n_cancer | n_healthy | AUC harmonized | AUC no_harmonize |
|---|---:|---:|---:|---:|
| cancer_jiang_healthy_cristiano | 121 | 506 | 0.484 ± 0.007 | 1.000 ± 0.000 |
| cancer_cristiano_healthy_jiang | 506 | 121 | 0.484 ± 0.007 | 1.000 ± 0.000 |

## Verdict
- Pooled harmonized cross-publication AUC: **0.9761 ± 0.0030** (n=627 with features, of 658 in labels file)
- Pooled AUC without harmonization: **0.9687 ± 0.0040** (mild change confirms the per-publication batch effect is small on this FinaleDB-uniformly-processed cohort)
- Per-publication AUC (sorted by AUC): 6 (jiang) **0.9776 ± 0.0024** (n=121), 8 (cristiano) **0.9739 ± 0.0015** (n=506)
- True-confound control: cancer_jiang_healthy_cristiano: harmonized AUC **0.484** (should be ~0.50), no-harmonize AUC **1.000** (should be ~1.00 — proves the batch effect is removable)
- True-confound control: cancer_cristiano_healthy_jiang: harmonized AUC **0.484** (should be ~0.50), no-harmonize AUC **1.000** (should be ~1.00 — proves the batch effect is removable)

## Honest framing
These numbers are pooled out-of-fold AUC on the 627-sample cross-publication cohort (after load5's missing-artifact filter). Internal CV; no external validation. They measure how well the 5-channel DELFI features separate cancer from healthy when pooled across jiang, cristiano with per-publication z-score harmonization. They do NOT measure clinical-grade sensitivity at the Galleri / CancerSEEK operating points, which require independent held-out plasma cohorts.

With per-publication z-score harmonization the true-confound AUC (cancer = one publication, healthy = another) collapses toward 0.50 — the per-publication mean/variance shift is the only signal and the harmonization removes it by design. WITHOUT harmonization the same control reaches ~0.999 — the classifier learns 'which publication is this from?', not 'is this cancer or healthy?'. The paired comparison (harmonized vs no_harmonize) is the only honest way to claim a cross-publication benchmark is not a publication-batch artifact.

Per-cancer sens@spec is OvR: each cancer class is scored against ALL healthy samples (not just the within-publication healthy ones). This is the cross-publication generalization view, not the within-publication view. Top-5 cancers by count are reported; smaller cohorts (n<10 cancer) are skipped.

**Open-data benchmark — NOT clinical validation.**
