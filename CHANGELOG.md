# Changelog

All notable changes to DeepCatch are recorded here. The project follows
[semantic versioning](https://semver.org/); the Zenodo DOI for each
release is minted from the corresponding GitHub release.

## [2.3.0] — 2026-11-16

The release that turns the internal benchmark into an externally validated
one, and replaces the neural fusion layer with an estimator that fits the
cohort size.

### Added

- **External validation on held-out public cohorts.** FinaleDB publications 1
  (Snyder 2016), 7 (Sun 2019) and 9 (Adalsteinsson 2017) are fetched from the
  public API + S3 — no IRB, no data-use agreement — and scored by a model
  trained only on publications 6 + 8. External AUC **0.9915** [0.9712, 1.0]
  on publication 7 (WGS, n=31); **0.8819** [0.7701, 0.9936] pooled over the
  publications carrying both classes. `scripts/fetch_external_cohort.py`,
  `scripts/external_validation_finaledb.py`.
- **Leave-one-study-out transfer measurement.** Train on one publication, test
  on the other: 0.6805 (Jiang → Cristiano) and 0.9361 (Cristiano → Jiang)
  against in-domain 0.9789 / 0.9701, with a within-study leave-one-cancer-type
  control showing that an unseen cancer type transfers at 0.857–1.000.
  `scripts/cross_study_loso.py`.
- **Regularized late fusion and Bayesian model averaging** for the Stage-1
  multi-modal layer: 0.9825 ± 0.0094 and 0.9650 ± 0.0089 on the honest
  paired TCGA-LUAD protocol, against the neural encoder's 0.554 and an LR
  baseline re-computed in-run at 0.9042.
  `src/foundation/regularized_fusion.py`,
  `scripts/fusion_regularized_smoke.py`.
- **Self-contained FinaleDB feature extraction** — the 5-channel DELFI + FSD
  profile, bit-identical to the sibling pipeline's extractor, accepting both
  FinaleDB fragment schemas. `src/fragmentomics/finale_features.py`.
- **The methylation path runs on real public data**: fragment file → BAM →
  per-CpG matrix → WGBS prior → HMM decode → `src/methylation_gnn`.
  `scripts/finaleme_bam_from_frag.py`, `scripts/finaleme_decoder_input.py`,
  `scripts/methylation_real_smoke.py`, `docs/METHYLATION_PIPELINE.md`.

### Fixed

- The FinaleDB readiness probe reported the S3 bucket as "private" because it
  probed a single inaccessible key; the bucket is partially public and the
  fragment objects are reachable.
- The foundation smoke's label-permutation control was a no-op — it permuted
  index sets and reassigned the same values, so it reproduced `y_true`
  exactly. The corrected control's null on the paired design is 0.42 ± 0.12,
  not 0.5, and is now measured rather than assumed.
- `scripts/finaleme_pipeline.py` invoked a FinaleMe class and CLI flags that
  do not exist in the released tool.
- `scripts/ultraearly_readiness.py` never put the repo root on `sys.path`, so
  the stage-breakdown check failed with "No module named 'src'" and the
  readiness CI job had been red since October.
- `src/foundation/__init__.py` eagerly imported torch, so a collection error
  in any test file aborted the whole "Core tests" job — zero tests ran.

### Changed

- The readiness gate is artifact-driven and now passes; all eight CI jobs are
  green.
- README, MODEL.md and the manuscript state the external and leave-one-study-out
  numbers alongside the pooled one, and describe the methylation channel as
  *extraction on real data, not yet a validated feature*.

## [2.2.0] — 2026-08-21

Panel-based ultra-sensitive ctDNA detection benchmark; tumor-naive
fragmentomics channel; DeLong CIs; decision curve; AUC gate; MODEL.md.
