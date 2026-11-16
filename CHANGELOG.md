# Changelog

All notable changes to DeepCatch are recorded here. The project follows
[semantic versioning](https://semver.org/); the Zenodo DOI for each
release is minted from the corresponding GitHub release.

## [2.3.0] — 2026-11-16

The release that turns the internal benchmark into an externally validated
one, and replaces the neural fusion layer with an estimator that fits the
cohort size.

### Added

- **The methylation GNN can run on real public data.**
  `src/methylation_gnn/` was architecture-only: the graph builder filled
  the node set with randomly generated regions when none were loaded, and
  the data module filled the beta matrix with random values, so the
  "methylation" the GNN consumed had no relationship to any real sample.
  `src/methylation_gnn/real_data.py` reads the same FinaleMe decoded
  output the channel was measured on and returns the structure the graph
  builder expects, with the **real CpG coordinates** as the graph's nodes.
  On the published cohort it loads the 31 real samples over the same
  1,550 most widely covered loci the channel used, applying the same
  selection rule and recording which rule it used. Region-level
  annotations the FinaleMe output does not carry (GC content, CpG
  density) are left absent rather than filled with random numbers, and a
  run pointed at real loci can no longer fall back to synthetic regions.
  `src/methylation_gnn/real_data.py`,
  `test/test_methylation_gnn_real_data.py`.
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
  observations-per-CpG count alone. At the 103,863 loci where both arms have
  coverage, each tested on the samples that cover it, the median
  cancer-minus-control difference in decoded β is exactly **0.0** (permutation
  p = 1.0), and only 0.035% of the
  decoder's calls carry a value strictly between 0 and 1. The cause is in the
  libraries: they are size-selected at ~300 bp (0.000% of fragments longer)
  and deliver 1.09 observations per CpG, so the digestion signature FinaleMe
  reads is absent. `scripts/methylation_channel_diagnosis.py`,
  `scripts/methylation_library_check.py`,
  `scripts/methylation_cohort_channel.py` (depth instrumentation),
  `docs/METHYLATION_PIPELINE.md` §5.
- **A reproduction driver that checks itself.** `paper/REPRODUCE.sh` now
  regenerates *every* artifact the manuscript quotes — the two cross-study
  sweeps, the shuffled-label null, the LOSO transfer, the external
  validation, the fusion ablation, the methylation channel and the figures —
  not just the two cross-study sweeps it covered before. It then compares
  the regenerated numbers against the committed ones, number by number
  (`scripts/reproduce_check.py`), so a run that succeeds but writes a
  different number fails instead of passing quietly. `--quick` runs the
  same pipeline with fewer seeds into a scratch directory and leaves the
  repository untouched. The data paths are environment-overridable
  (`DEEPCATCH_FEATURES_DIR`, `DEEPCATCH_LABELS`, `DEEPCATCH_EXTERNAL_DIR`,
  `DEEPCATCH_METHYLATION_DIR`). `paper/REPRODUCE.sh`,
  `scripts/reproduce_check.py`, `paper/REPRODUCE.md`.
  The driver's pytest gate now runs on macOS's bash 3.2 as well (expanding
  an empty array under `set -u` raised "unbound variable" there), and the
  `--timeout` flag is only passed when the `pytest-timeout` plugin is
  actually installed — it was never declared as a dependency, so on a
  clean checkout pytest rejected the unknown argument and the gate
  failed.

### Fixed

- **`import src.methylation_gnn` needed torch, which silenced the whole
  "Core tests" CI job.** The package imported the model, the trainer and
  the inference module eagerly, and those subclass `torch.nn.Module` at
  import time — so the package raised in any environment without torch.
  The "Core tests" job installs only numpy/scipy/scikit-learn/pandas/
  matplotlib/seaborn/pytest, so a test file that imported the package
  failed at **collection**, and pytest aborted the run: zero tests ran.
  This is the same class of defect `src/foundation` had. The
  torch-dependent names are now resolved lazily (PEP 562), the
  torch-free modules (`graph_builder`, `config`, `integration`, `data`,
  `real_data`) import eagerly as before, and asking for a torch-dependent
  name without torch raises an ImportError that names the dependency
  instead of "'NoneType' object has no attribute 'Module'".
  `src/methylation_gnn/__init__.py`, `src/methylation_gnn/integration.py`,
  `src/methylation_gnn/gnn_model.py`,
  `test/test_methylation_gnn_real_data.py`.
- **The manuscript conflated two different locus sets.** "At the 1,550
  CpGs that both arms cover ... over 103,863 CpGs" mixed two numbers that
  come from two different tests. The 1,550 are the *most widely covered*
  loci — the pool the sample-level shared-locus AUC is computed on, and
  the pool the channel script falls back to when too few loci reach its
  80% coverage threshold (only 16 do). The 103,863 are the loci where
  *both* arms have coverage, each tested on the samples that cover it,
  and they are what the per-CpG median difference of exactly 0.0 is
  measured over. The two sets are not the same, and the 1,550 loci are
  not all covered by both arms. Corrected in the preprint, the methods
  paper, the markdown twin, README, MODEL, REVIEWERS, the readiness doc,
  the landing page and the changelog.
- **The submission PDF's metadata was stale.** The PDF is 8 pages, not
  the 7 three documents claimed, and 264 KB, not 261 KB — the
  dataset-naming edit pushed the references onto an eighth page. The
  PDF has been rebuilt from the corrected `paper.tex`.
- **The bioRxiv upload form quoted a stale abstract length** (2905
  against the manuscript's 2897). `scripts/check_biorxiv_form.py` now
  fails if the form and the manuscript disagree, and a test asserts the
  changelog states the current length.
- **The methylation channel's shared-locus number is not reproducible, and
  the cause is upstream.** FinaleMe's HMM training is not deterministic: on
  byte-identical decoder inputs the trained model differs between runs,
  which flips 16–20 of ~51,000 CpG calls per sample, and the shared-locus
  ("common CpG set") AUC moves with it — 0.8034 / 0.8120 / 0.812 /
  0.8162 / 0.8291 over five runs of the same cohort. Every headline
  number was unaffected: the decoded-β AUC (0.9444), the covered-CpG
  baseline (0.9444), the observations-per-CpG baseline (0.9444), the
  depth-matched pair (0.9444 / 0.8034) and the per-CpG median difference
  (exactly 0.0 over 103,863 CpGs, permutation p = 1.0) were identical on
  every run. The shared-locus row is now reported as a range in the
  documents that quote it, and the artifact carries a
  `reproducibility_note` saying so. A secondary order-dependence in the
  scoring code — the shared loci were read straight out of a `set`, so the
  float summation behind `np.mean` depended on `PYTHONHASHSEED` — was
  fixed at the same time: the loci are now read in sorted order, and two
  hash seeds on the same decoded files return the same value.
  `scripts/methylation_cohort_channel.py`,
  `test/test_methylation_hashseed_stability.py`.
- **The external-validation artifact predated the PCA-seeding fix.**
  `results/external_validation_finaledb.json` was generated before the
  PCA was seeded and was never regenerated, so the training-derived
  99%-specificity threshold it recorded (0.9887) came from an unseeded
  run. The seeded pipeline gives **0.9766**, and the sensitivity of the
  case-only publication 9 at that threshold moves from 0.4034 to
  **0.4454**. The external AUCs the paper quotes (publication 7 at
  **0.9915**, pooled over the publications carrying both classes at
  **0.8819**) are unchanged. The artifact, the narrative and the two
  numbers in `paper/METHODS_PAPER.md` §3.2.2 have been corrected, the
  threshold is now quoted in `docs/EXTERNAL_VALIDATION.md`, and both
  numbers are on the quotable-number list so the drift guard covers
  them.
- **`--include-screening` crashed the per-cancer writer.** `SCREENING_SPEC_GRID`
  was defined with a trailing comma, which made it a one-element tuple
  *containing* the grid rather than the grid. The per-cancer standalone
  writer then passed that tuple where a float specificity was expected and
  died with `TypeError: unsupported operand type(s) for -: 'float' and
  'tuple'` — an error naming neither the flag nor the grid. It broke exactly
  the combination the reproduction driver needs (`--include-screening` +
  `--out-per-cancer-json`). The grid is flat again, the normaliser now
  rejects a nested grid with a message that names both, and the artifact's
  own `provenance.specificities` reports the grid that was actually used
  instead of always claiming the legacy three.
  `scripts/cross_study_finallydb.py`, `test/test_screening_spec_grid.py`.
- **The reproduction driver could overwrite a published artifact.**
  `cross_study_finallydb.py` defaults `--out-per-cancer-json` to
  `results/per_cancer_sens_at_spec.json`, which belongs to the GC-corrected
  run. The driver ran the legacy (no-GC) sweep without naming that flag, so
  the legacy run overwrote the committed GC-corrected artifact with legacy
  numbers (BRCA 0.9738 → 0.9794) — a silent drift of a number the paper
  quotes. Every published-path output flag is now named explicitly on every
  invocation, and a test walks each script's argparse defaults and fails if
  the driver leaves one unset.
  `paper/REPRODUCE.sh`, `test/test_paper_reproduce_driver.py`.
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
  (0.9361), which the abstract now states. A second check fails if a
  *superseded* figure reappears in a current-facing document. The
  quotable-number list lives once, in `scripts/benchmark_numbers.py`,
  which prints every quoted figure with its artifact and which the test
  imports, so the test and the printer cannot disagree.
- **The figures were showing the legacy benchmark.** `scripts/plot_figures.py`
  read `results/cross_study_finallydb.json` (the no-GC artifact, pooled AUC
  0.9770) while the manuscript quotes the GC-corrected 0.9684 -- so Figure 1
  in the submission PDF plotted a different curve from the one the text
  describes. The script now reads the headline artifact and falls back to
  the legacy one only when it is absent; the three affected figures were
  regenerated and the PDF rebuilt. The paper-number guard now also
  asserts which artifact the figure script reads, so the same drift
  cannot come back through the images.
- **The abstract fitted bioRxiv's character limit again.** It grew with the
  v2.3.0 results to 3119 plain-text characters against bioRxiv's 3000
  limit, so the upload would have been refused; it is now 2897 with every
  number intact. A test fails the build if the abstract exceeds 3000
  characters, or creeps within 50 of the limit, and
  `scripts/check_biorxiv_form.py` fails if the upload form quotes a
  different length than the manuscript actually has.
- **Eight tests that no workflow ever ran now run in CI.**
  `test/test_foundation_smoke.py` needs torch, so the core job skipped
  every test in it, and the deep-learning job does not collect `test/`
  -- so its checks (the LR baseline must beat the naive average, the
  foundation must stay inside the LR band, and the pipeline must not
  fall back to synthetic data when no cache is present) were never
  executed anywhere. The deep-learning job now runs that file, and the
  seven of its tests that use the default invocation share one smoke
  run through a module-scoped fixture instead of each repeating the
  ~27 s training run -- so the file costs 37 s rather than 3.5 minutes.
- **The drift guard now covers the collaborator-facing documents.** It
  checked the manuscript but not README.md, MODEL.md,
  ULTRA_EARLY_READINESS.md or the per-experiment write-ups -- the files a
  wet-lab collaborator reads before deciding whether to send us data.
  Those now have to quote the committed artifacts too, and the
  superseded-value sweep covers them as well.
- **The preprint's data-availability statement named the wrong datasets.**
  It listed TCGA-LUAD and FLARE but not FinaleDB -- which is where the
  cross-study benchmark (n=627) and the held-out external cohort
  (publications 1, 7, 9) come from, i.e. the paper's headline numbers.
  It now names every source, states that no patient-identifiable data is
  redistributed, and records that the external cohort needs no IRB. The
  manuscript also gains the standard author-contributions, funding and
  competing-interests sections that life-sciences journals require.
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
