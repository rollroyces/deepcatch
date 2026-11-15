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
- **A methylation-channel diagnosis that reports a null.** The channel was
  measured on 31 real samples (FinaleDB publication 7; 17 colorectal cancer,
  14 liver-transplant controls) and it separates the arms with leave-one-out
  AUC **0.9444** — exactly the AUC of the covered-CpG count alone and of the
  observations-per-CpG count alone. Restricted to the 1,550 loci both arms
  cover, the median cancer-minus-control difference in decoded β is exactly
  **0.0** over 103,863 CpGs (permutation p = 1.0), and only 0.035% of the
  decoder's calls carry a value strictly between 0 and 1. The cause is in the
  libraries: they are size-selected at ~300 bp (0.000% of fragments longer)
  and deliver 1.09 observations per CpG, so the digestion signature FinaleMe
  reads is absent. `scripts/methylation_channel_diagnosis.py`,
  `scripts/methylation_library_check.py`,
  `scripts/methylation_cohort_channel.py` (depth instrumentation),
  `docs/METHYLATION_PIPELINE.md` §5.

### Fixed

- **The benchmark is now reproducible.** The cross-study benchmark fitted PCA
  without `random_state`. sklearn picks the randomised solver when
  `n_components < 0.8 * min(n_samples, n_features)` — true for every fold
  here — and that solver is seeded from `random_state`, so re-running the
  benchmark produced a slightly different number every time (one
  publication's AUC moved 0.9789 → 0.9821 between two runs of the same
  command). Four call sites are now seeded:
  `scripts/cross_study_finallydb.py`, `src/fragmentomics/decision_curve_cli.py`
  (which also fitted the PCA twice per fold), `src/fragmentomics/train_tumor_naive.py`
  and `src/fragmentomics/fusion_ablation.py`.
  `test/test_benchmark_reproducibility.py` installs an AST guard that fails
  on any unseeded `PCA(...)` in `scripts/` or `src/`, plus a behavioural test
  that requires two runs of the fold fit to agree exactly and a different
  seed to still give a different answer.
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

- **The benchmark artifacts were regenerated from the seeded pipeline**, and every
  document now quotes that one run. The committed set was internally
  inconsistent before this: `results/cross_study_finallydb.json` came from a
  PCA-200 run (pooled 0.9747) while `docs/CROSS_STUDY_BENCHMARK.md` and
  `results/per_cancer_sens_at_spec.json` came from a PCA-20 run (pooled
  0.9225), and the paper quoted numbers from both. The canonical numbers are
  now: **pooled harmonized 0.9684 ± 0.0027** with GC correction (the
  headline) and **0.9770 ± 0.0013** legacy without it — a GC-axis effect of
  −0.0086, the same finding as before, now measured on one reproducible
  footing. Per-cancer, sens@spec, the true-confound control and the
  shuffled-label null were regenerated from the same runs
  (`scripts/benchmark_numbers.py` prints all of them).
- **The fusion ablation was re-run from this repository, with the labels/study
  columns resolved from the header.** It had never run here before: the labels
  loader raised on the canonical 4-column `labels_multiclass.tsv`, and behind
  that lay a second bug — the study column was read positionally, so every
  sample was assigned to a "study" called cancer or control, the per-study
  harmonisation removed the signal it was meant to isolate, and the pooled
  tumor-naive AUC came out at 0.4915 (chance). The corrected, seeded run
  reproduces the published numbers: tumor-naive **0.9740 ± 0.0013**, naive
  average **0.9884 ± 0.0010** (paired over 10 seeds: +0.0144, t = 54.50,
  positive on every seed), LR fusion 0.9884 ± 0.0009, isotonic variant
  0.9849 ± 0.0014 (it regresses on every seed, as before). The preprint
  abstract now quotes these and states that the fused mutation channel is a
  calibrated synthetic proxy rather than the real panel LLR.
  Artifact: `results/fusion_ablation_isotonic.json`.
- **The leave-one-study-out artifact was regenerated and reproduced itself.**
  Re-running `scripts/cross_study_loso.py` against the seeded benchmark
  changed nothing but the timestamp and a stale note (the note said the
  headline benchmark's PCA "is unseeded"; it was not any more) -- every
  transfer number came back identical: jiang -> cristiano 0.6805,
  cristiano -> jiang 0.9361, in-domain 0.9789 / 0.9690. That is the
  cleanest available evidence that the pipeline is now deterministic.
- **A test now keeps the manuscript and the artifacts in step.**
  `test/test_paper_numbers_match_artifacts.py` reads the headline
  numbers out of the committed JSON and requires the documents that cite
  them to contain them, so regenerating an artifact without updating the
  prose fails the build. It caught one gap immediately: the preprint
  quoted the single-type transfer direction (0.6805) but not the reverse
  (0.9361), which the abstract now states.
- The readiness gate is artifact-driven and now passes; all eight CI jobs are
  green.
- README, MODEL.md and the manuscript state the external and leave-one-study-out
  numbers alongside the pooled one, and describe the methylation channel as
  *extraction on real data, not yet a validated feature*.
- `USAGE.md` pointed at `scripts/honest_benchmark.py`, which does not exist in
  this repository; it now gives this repo's benchmark command and its expected
  output.

## [2.2.0] — 2026-08-21

Panel-based ultra-sensitive ctDNA detection benchmark; tumor-naive
fragmentomics channel; DeLong CIs; decision curve; AUC gate; MODEL.md.
