# Cross-Study FinaleDB Benchmark (Open Data)
> **Scope**: Open-data benchmark on FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019). **NOT** clinical validation. **NOT** external cohort validation. Pooled OOF on the same cohort that trained the model.
- Generated: `2026-09-23T08:11:52.390185+00:00`
- Classifier: `LogisticRegression(max_iter=2000)`
- Feature set: 5-channel (5mb_ratio + 5mb_coverage + 100kb_ratio + 100kb_counts + FSD-196)
- PCA n_components: 200 (capped at `min(n_train, n_features)` per fold)
- CV: 5-fold StratifiedKFold, 5-seed pooled OOF
- Seeds: `[42, 13, 7, 99, 1234]`
- Bootstrap: n=1000, seed=2026
- Cell-line regex: `^(GM\d+|HeLa|HepG2|K562|HL60|Jurkat|Raji|MCF7|U937|THP1|HEK293|HCT116|SW480|A549|GM12878)`

## 1. Cohort inventory
- Samples in labels_multiclass.tsv: **658** (364 cancer + 294 healthy)
- Samples with all 5-channel DELFI features: **627** (31 dropped due to missing artifacts)
- Cell-line filter removed: **0** samples (none matched the regex on this open-data cohort)
- Per-study:
  - **jiang**: n=121 (89 cancer + 32 healthy)
  - **cristiano**: n=537 (275 cancer + 262 healthy)
- Per-cancer (cancer samples only):
  | Cancer | n_total | n_jiang | n_cristiano |
  |---|---:|---:|---:|
  | HCC_J | 89 | 89 | 0 |
  | LUAD | 79 | 0 | 79 |
  | PAAD | 60 | 0 | 60 |
  | BRCA | 54 | 0 | 54 |
  | OV | 28 | 0 | 28 |
  | CRC | 27 | 0 | 27 |
  | OTHER_C | 27 | 0 | 27 |

## 2. Per-cohort AUC
Each study evaluated independently (no harmonization needed; one study only).

| Study | n | n_cancer | n_healthy | AUC (5-seed mean ± std) |
|---|---:|---:|---:|---|
| jiang | 121 | 89 | 32 | 0.9791 ± 0.0028 |
| cristiano | 506 | 274 | 232 | 0.9704 ± 0.0013 |

## 3. Pooled cross-study AUC (with/without per-study harmonization)
Harmonization = per-study z-score StandardScaler fit on train fold only.

| Setting | AUC mean ± std | Sens@95% | Sens@98% | Sens@99% |
|---|---|---:|---:|---:|
| harmonized | 0.9747 ± 0.0012 | 0.909 | 0.832 | 0.793 |
| no_harmonize | 0.9661 ± 0.0019 | 0.890 | 0.818 | 0.755 |

## 4. Per-cancer Sens@Spec (top-5 cancers by sample count)
One-vs-rest: each cancer vs ALL healthy samples in the pooled cross-study cohort. Per-study harmonization inside each CV fold.

**Primary CIs are DeLong** (DeLong, DeLong, Clarke-Pearson 1988 for AUC; Sun & Xu 2014 for Sens@spec — equivalent to DeLong-Han-Agarwal structural component restricted to positives). Bootstrap 95% CIs (n=1000 percentile resamples) are preserved under `bootstrap_ci` in `results/cross_study_finallydb.json` for the audit trail (narrower at n>=10; DeLong is the standard cfDNA reference per the cfdna-early-detection-validation skill). The standalone clinical-decision table — AUC + Sens@spec CIs + PPV@prev at spec=0.99 — is at `results/per_cancer_sens_at_spec.json` (schema per `src.per_cancer_sens_at_spec.build_per_cancer_table`).

**Caveat**: HCC_J is Jiang-only (n=89 cancer + 32 healthy). Its per-cancer OvR uses 89 HCC vs ALL healthy controls (264 controls in the OOF subset) — the per-cancer denominator is small but the OvR is well-defined. Other top-5 cancers are Cristiano-only.

| Cancer | n_cancer | n_healthy | AUC mean ± std | AUC DeLong [95% CI] | Sens@95% DeLong [95% CI] | Sens@98% DeLong [95% CI] | Sens@99% DeLong [95% CI] |
|---|---:|---:|---|---|---|---|---|
| HCC_J | 89 | 264 | 0.7292 ± 0.0205 | 0.7532 [0.6879–0.8186] | 0.382 [0.281–0.483] | 0.292 [0.198–0.386] | 0.124 [0.057–0.190] |
| LUAD | 79 | 264 | 0.9806 ± 0.0021 | 0.9830 [0.9715–0.9946] | 0.911 [0.848–0.975] | 0.848 [0.768–0.928] | 0.633 [0.526–0.740] |
| PAAD | 60 | 264 | 0.9402 ± 0.0056 | 0.9416 [0.9023–0.9809] | 0.833 [0.738–0.929] | 0.700 [0.583–0.817] | 0.467 [0.341–0.593] |
| BRCA | 53 | 264 | 0.9747 ± 0.0023 | 0.9751 [0.9518–0.9983] | 0.943 [0.879–1.000] | 0.849 [0.751–0.947] | 0.491 [0.356–0.625] |
| OV | 28 | 264 | 0.9924 ± 0.0048 | 0.9969 [0.9927–1.0000] | 1.000 [0.965–1.000] | 0.929 [0.827–1.000] | 0.929 [0.827–1.000] |

## 5. True-confound control
Cancer = 100% from one study, healthy = 100% from the other. Without harmonization the classifier learns 'which study is this from?' (AUC ~0.999). With per-study z-score harmonization the study-specific mean/variance is the only signal and is removed by design (AUC should collapse toward 0.50).

| Orientation | n_cancer | n_healthy | AUC harmonized | AUC no_harmonize |
|---|---:|---:|---:|---:|
| cancer_jiang_healthy_cristiano | 121 | 506 | 0.499 ± 0.013 | 0.999 ± 0.001 |
| cancer_cristiano_healthy_jiang | 506 | 121 | 0.494 ± 0.012 | 0.999 ± 0.001 |

## Verdict
- Pooled harmonized cross-study AUC: **0.9747 ± 0.0012** (n=627 with features, of 658 in labels file)
- Pooled AUC without harmonization: **0.9661 ± 0.0019** (mild change confirms the per-study batch effect is small on this FinaleDB-uniformly-processed cohort)
- Per-cohort AUC: Jiang **0.9791 ± 0.0028** (n=121), Cristiano **0.9704 ± 0.0013** (n=506)
- True-confound control: cancer=Jiang + healthy=Cristiano: harmonized AUC **0.499** (should be ~0.50), no-harmonize AUC **0.999** (should be ~1.00 — proves the batch effect is removable)
- True-confound control: cancer=Cristiano + healthy=Jiang: harmonized AUC **0.494** (should be ~0.50), no-harmonize AUC **0.999**

## Honest framing
These numbers are pooled out-of-fold AUC on the 627-sample cross-study cohort (after load5's missing-artifact filter). Internal CV; no external validation. They measure how well the 5-channel DELFI features separate cancer from healthy when pooled across Jiang 2015 and Cristiano 2019 with per-study z-score harmonization. They do NOT measure clinical-grade sensitivity at the Galleri / CancerSEEK operating points, which require independent held-out plasma cohorts.

With per-study z-score harmonization the true-confound AUC (cancer = one study, healthy = the other) collapses toward 0.50 — the per-study mean/variance shift is the only signal and the harmonization removes it by design. WITHOUT harmonization the same control reaches ~0.999 — the classifier learns 'which study is this from?', not 'is this cancer or healthy?'. The paired comparison (harmonized vs no_harmonize) is the only honest way to claim a cross-study benchmark is not a study-batch artifact.

Per-cancer sens@spec is OvR: each cancer class is scored against ALL healthy samples (not just the within-study healthy ones). This is the cross-study generalization view, not the within-study view. Top-5 cancers by count are reported; smaller cohorts (n<10 cancer) are skipped.

**Open-data benchmark — NOT clinical validation.**
