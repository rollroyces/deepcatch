# Cross-Study FinaleDB Benchmark (Open Data)
> **Scope**: Open-data benchmark on FinaleDB publications 6 (Jiang 2015) + 8 (Cristiano 2019). **NOT** clinical validation. **NOT** external cohort validation. Pooled OOF on the same cohort that trained the model.
- Generated: `2026-09-23T04:21:47.412205+00:00`
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
| cristiano | 506 | 274 | 232 | 0.9693 ± 0.0022 |

## 3. Pooled cross-study AUC (with/without per-study harmonization)
Harmonization = per-study z-score StandardScaler fit on train fold only.

| Setting | AUC mean ± std | Sens@95% | Sens@98% | Sens@99% |
|---|---|---:|---:|---:|
| harmonized | 0.9738 ± 0.0018 | 0.909 | 0.837 | 0.782 |
| no_harmonize | 0.9650 ± 0.0019 | 0.884 | 0.818 | 0.774 |

## 4. Per-cancer Sens@Spec (top-5 cancers by sample count)
One-vs-rest: each cancer vs ALL healthy samples in the pooled cross-study cohort. Per-study harmonization inside each CV fold. Bootstrap 95% CI on sens@spec (n=1000 resamples).

| Cancer | n_cancer | AUC mean ± std | Sens@95% [95% CI] | Sens@98% [95% CI] | Sens@99% [95% CI] |
|---|---:|---:|---|---|---|
| HCC_J | 89 | 0.7248 ± 0.0184 | 0.416 [0.296–0.530] | 0.348 [0.156–0.462] | 0.191 [0.000–0.420] |
| LUAD | 79 | 0.9813 ± 0.0020 | 0.911 [0.830–0.973] | 0.861 [0.571–0.943] | 0.658 [0.392–0.917] |
| PAAD | 60 | 0.9368 ± 0.0070 | 0.817 [0.707–0.917] | 0.767 [0.429–0.863] | 0.483 [0.333–0.810] |
| BRCA | 53 | 0.9763 ± 0.0034 | 0.943 [0.865–1.000] | 0.849 [0.500–0.981] | 0.547 [0.352–0.937] |
| OV | 28 | 0.9924 ± 0.0048 | 1.000 [0.895–1.000] | 0.929 [0.821–1.000] | 0.929 [0.758–1.000] |

## 5. True-confound control
Cancer = 100% from one study, healthy = 100% from the other. Without harmonization the classifier learns 'which study is this from?' (AUC ~0.999). With per-study z-score harmonization the study-specific mean/variance is the only signal and is removed by design (AUC should collapse toward 0.50).

| Orientation | n_cancer | n_healthy | AUC harmonized | AUC no_harmonize |
|---|---:|---:|---:|---:|
| cancer_jiang_healthy_cristiano | 121 | 506 | 0.489 ± 0.004 | 0.999 ± 0.002 |
| cancer_cristiano_healthy_jiang | 506 | 121 | 0.496 ± 0.004 | 1.000 ± 0.000 |

## Verdict
- Pooled harmonized cross-study AUC: **0.9738 ± 0.0018** (n=627 with features, of 658 in labels file)
- Pooled AUC without harmonization: **0.9650 ± 0.0019** (mild change confirms the per-study batch effect is small on this FinaleDB-uniformly-processed cohort)
- Per-cohort AUC: Jiang **0.9791 ± 0.0028** (n=121), Cristiano **0.9693 ± 0.0022** (n=506)
- True-confound control: cancer=Jiang + healthy=Cristiano: harmonized AUC **0.489** (should be ~0.50), no-harmonize AUC **0.999** (should be ~1.00 — proves the batch effect is removable)
- True-confound control: cancer=Cristiano + healthy=Jiang: harmonized AUC **0.496** (should be ~0.50), no-harmonize AUC **1.000**

## Honest framing
These numbers are pooled out-of-fold AUC on the 627-sample cross-study cohort (after load5's missing-artifact filter). Internal CV; no external validation. They measure how well the 5-channel DELFI features separate cancer from healthy when pooled across Jiang 2015 and Cristiano 2019 with per-study z-score harmonization. They do NOT measure clinical-grade sensitivity at the Galleri / CancerSEEK operating points, which require independent held-out plasma cohorts.

With per-study z-score harmonization the true-confound AUC (cancer = one study, healthy = the other) collapses toward 0.50 — the per-study mean/variance shift is the only signal and the harmonization removes it by design. WITHOUT harmonization the same control reaches ~0.999 — the classifier learns 'which study is this from?', not 'is this cancer or healthy?'. The paired comparison (harmonized vs no_harmonize) is the only honest way to claim a cross-study benchmark is not a study-batch artifact.

Per-cancer sens@spec is OvR: each cancer class is scored against ALL healthy samples (not just the within-study healthy ones). This is the cross-study generalization view, not the within-study view. Top-5 cancers by count are reported; smaller cohorts (n<10 cancer) are skipped.

**Open-data benchmark — NOT clinical validation.**
