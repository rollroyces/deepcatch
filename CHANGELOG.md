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
  trained only on publications 6 + 8. External AUC **0.9444** [0.8731, 1.0]
  on publication 7 (WGS, n=31); **0.708** [0.5148, 0.9012] on publication 1
  (WGS, n=52); **0.7385** [0.6281, 0.8488] pooled over the publications
  carrying both classes. `scripts/fetch_external_cohort.py`,
  `scripts/external_validation_finaledb.py`.
- **Leave-one-study-out transfer measurement.** Train on one publication, test
  on the other: 0.6805 (Jiang → Cristiano) and 0.9361 (Cristiano → Jiang)
  against in-domain 0.9789 / 0.969, with a within-study leave-one-cancer-type
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

- **The drift guard covered fewer documents than quote the numbers.**
  An audit that walks every quoted number and finds every document
  quoting it showed that the legacy (no-GC) benchmark's pooled AUC was
  quoted by nine documents while the guard named four, and that the
  leave-one-study-out in-domain pair and both transfer directions, and
  the true-confound control, were quoted by the team page, the model
  card, the results index and the Simplified-Chinese README without
  being named at all — so those copies could drift with nothing
  failing. The legacy and GC-corrected pooled AUCs, the true-confound
  control and the LOSO in-domain pair are now required of every
  document that quotes them, and the per-cancer sens@99 row is required
  of the readiness documents that print it.
  `scripts/benchmark_numbers.py`.
- **A cohort count mixed two different numbers.** The readiness
  document said "n=627 samples with 364 cancer + 263 healthy"; 627 is
  363 + 264, and 364 + 294 is the *candidate* count the labels file
  holds before 31 samples are dropped for missing features. Both numbers
  are now stated, and each in its place.
  `ULTRA_EARLY_READINESS.md`.
- **The results directory's index described the v1 simulation
  benchmark.** `results/README.md` was last updated on 2026-04-28: its
  headline table was the old simulation sweep (AUC 0.9610 at 1% ctDNA
  against re-implemented Bie/CAPP-Seq/iDES), its file index named three
  files that are not in the repository, and its honest-limitations list
  said there was no tissue-of-origin prediction — the module exists and
  is wired into the fusion, though its 29-tissue reference atlas is
  synthetic because the real cfSort atlas is not here, and no
  tissue-of-origin artifact is committed. It is now a v2.3.0 index: the
  real-data headline results, each figure naming the artifact it comes
  from, the honest limitations as they now stand, and the file index
  that matches what is actually committed.
  `results/README.md`.
- **The model card's tumor-naive table quoted numbers no artifact
  holds.** Its two rows carried 0.9716 (Jiang, 121 samples) and 0.9746
  (cross-study, 627 samples), taken from the sibling repository's
  `scripts/honest_benchmark.py` before its PCA was seeded — and 0.9746
  is one seed's value out of the legacy artifact's five, not the pooled
  mean, so it did not describe the run it was labelled with. The table
  is now written from this repository's own artifacts (0.9661 ± 0.0031
  for Jiang, 0.9684 ± 0.0027 GC-corrected pooled, 0.9770 ± 0.0013
  legacy pooled), the two cross-study rows are kept apart because they
  are different benchmarks, and the per-publication row's sensitivity
  cells are left empty rather than filled from another run. Each row is
  pinned to its artifact by the drift guard.
  `MODEL.md`, `scripts/benchmark_numbers.py`.
- **A validation document named three artifacts that are not in the
  repository.** `docs/CADD_GDC_VALIDATION.md` is a dated one-off report
  (2026-09-17) whose per-subgroup table came from a 150-patient GDC run
  whose cohort, CADD-match and results JSONs are gitignored and
  uncommitted, and whose builder resolves its cache directory from an
  absolute path — so the run cannot be repeated from a clean checkout
  and a reader cannot verify a single figure in it. The document now
  says so at the top, and points at the committed, reproducible
  20-patient counterpart (`results/cadd_per_subgroup_llr.json`).
  `docs/CADD_GDC_VALIDATION.md`.
- **The decision-curve CLI read the study column positionally.** It
  resolved the *label* with the header-aware resolver but took the
  *study* from column 2 of each line. In the canonical FinaleDB layout
  (`sample  disease_class  label  study`) column 2 is the label, so the
  study became "cancer"/"healthy", the per-study harmonisation
  z-scored each class separately, and the run reported a chance-level
  AUC (0.4841 against 0.9773) instead of failing. The adapter's own
  docstring warns about exactly this, and the bug stayed invisible
  because the only labels file the CLI was ever run against was the
  headerless `labels_cross_study.tsv`, where the positional read
  happens to be right. Both the label and the study now come from the
  header-aware resolver, and a regression test builds a synthetic
  cohort with a class signal, runs the CLI against the canonical
  layout, and requires the signal to survive.
  `src/fragmentomics/decision_curve_cli.py`,
  `test/test_decision_curve_labels.py`.
- **The decision-curve table in RESULTS.md predated the seeded PCA.**
  Section 4's operating points were computed on 2026-08-28, before the
  PCA was seeded and before the CLI stopped refitting the PCA on the
  test fold (commit `96b38da`). They could not be reproduced by the
  committed code, and the method they came from is no longer the method
  in the repository. The table is now read out of a committed
  artifact, `results/decision_curve.json`, which the reproduction
  driver regenerates; the run is deterministic (two runs are
  byte-identical) and the operating points are pinned by the drift
  guard.
  `results/decision_curve.json`, `RESULTS.md`, `paper/REPRODUCE.sh`.
- **RESULTS.md's TL;DR quoted a superseded generation.** Its
  tumor-naive row carried the legacy 0.978 rather than the GC-corrected
  headline, its fusion row carried 0.989 where the fusion artifact holds
  0.9884, and its harmonization row carried a "0.499 (control) → 0.500
  (harmonized on true confound)" pair that no artifact backs — the
  true-confound control is 0.9989 without harmonization and 0.4836
  with it. The precision note also told readers that ±0.002 AUC
  between runs is convergence noise, which was true while the PCA was
  unseeded and is no longer true at all. All four are corrected, and
  RESULTS.md — the repository's headline results document — was covered
  by no check whatsoever, so it is now in the drift guard's document
  list, with its TL;DR rows pinned to one run each.
  `RESULTS.md`, `scripts/benchmark_numbers.py`.
- **Two PPV figures did not reproduce.** The 3.5%-prevalence column of
  the screening table read 38.7% and 75.6%; recomputed from the
  operating points the sibling repository's `results/ppv_screening.json`
  actually holds (sens 0.91 at spec 0.95, and sens 0.82 at spec 0.99),
  they are 39.8% and 74.8%. The table now names the repository the
  script and artifact live in, and says which columns are the
  artifact's own values and which are computed at the stated
  prevalences.
  `RESULTS.md`.
- **The project-team state table was a generation behind.** `TEAM.md` §7
  — the document a wet-lab lab reads first — still listed the
  per-cancer AUC, the methylation channel and the mutation panel as not
  done, quoted a fusion AUC the fusion artifact no longer holds, and
  said the hold-out validation was missing (the public hold-out
  validation exists; the *clinical* one does not). It is now written
  from the committed artifacts, and the honest list says what is
  genuinely not done.
  `TEAM.md`.
- **A committed generated document could drift from its generator.** The
  same hole that let the readiness report call the held-out cohort
  unusable was open for every other generated document, because nothing
  compared a document with the generator that writes it. A new test
  regenerates each offline generator's document *and* its artifact into
  a temporary directory and asserts the committed copies are exactly what
  the generator produces. The audit that motivated it found the others
  clean: the leave-one-study-out report, the external-validation report
  and the readiness report all regenerate byte-identically (modulo the
  generation timestamp and the commit stamp). A hand edit to a generated
  file now fails the build.
  `test/test_generated_docs_match_their_generators.py`.
- **The ORCID/affiliation flag now states the verified public record.**
  `docs/JOURNAL_SUBMISSION.md` flagged that the declared "Independent
  Researcher" affiliation and the ORCID record disagree, without saying
  what the record contains. It does now: one public employment
  (University of Michigan), no start date, no end date, no department and
  no role title — so ORCID renders it as a *current* employment — and no
  other employment, education, work or affiliation of any kind. Both
  resolutions are spelled out, with the files to change for each; which
  one is true is the author's fact, not the repository's.
  `docs/JOURNAL_SUBMISSION.md`.
- **The readiness report called the external cohort unusable.** The
  publication-readiness diagnostic only sees the *training* labels file
  and the *training* feature cache, so the held-out publications come
  out as `labels-missing` — and the document's "Honest reading" then
  called them "**BLOCKED** ... blocked on **FinaleDB infrastructure** ...
  the API/S3 are not reachable to fetch them", while the same document's
  own API verdict read `ok` and its S3 verdict read `partial`. They are
  not unusable: publications 1, 7 and 9 **are** the external test set,
  fetched from the public API/S3 and scored by
  `scripts/external_validation_finaledb.py`. The correcting note had also
  been hand-edited into a *generated* file, so the next run of the
  generator would have deleted it. The generator now reads the held-out
  cohort from the external-validation artifact, emits the note itself,
  separates "held out" from "blocked", and a test writes a deliberately
  different cohort into a fixture artifact and asserts the document
  follows it.
  `scripts/_verify_publication_readiness.py`,
  `docs/PUBLICATION_READINESS.md`.
- **The Chinese README stated the fast-path gate from an older run.** Its
  table carried the gate's collected/passed/skipped/deselected split and
  its file tree the fast-path test count, while the English original
  carried different ones — the translation was a test behind and
  nothing compared the two. Both READMEs now state what the gate
  actually produces, and the suite-table guard checks the English one
  against a real collection.
  `README.md`, `README.zh-CN.md`.

- **A stated limitation was not true.** The repository said twice that the
  cross-study benchmark's healthy arm was "FinaleDB TF=0 arms (modeled
  cancer-at-TF=0 read depth), not real healthy-donor plasma". Checked
  against FinaleDB's own metadata (2026-11-16), every one of the 294
  training-cohort control samples is recorded as `disease: Healthy`,
  `tissue: blood plasma` — 262 Cristiano donors and 32 Jiang donors. The
  control arm is real donor plasma, and the limitation has been replaced
  with the one that IS true: the per-study z-score uses a class-mixed
  per-study mean, and the two training studies have very different case
  mixes (Jiang is ~70% HCC, Cristiano mostly healthy), so a study's mean
  partly encodes its case mix.
  `MODEL.md`, `docs/PUBLIC_DATA_SOURCES.md`.
- **Three more stale numbers in the methods paper, caught by widening the
  drift guard.** The §3.2 honest-framing paragraph still quoted the
  pooled harmonised AUC as 0.9738 and the per-cohort rows as Jiang 0.9791
  / Cristiano 0.9693 — all from an unseeded run; the current values are
  0.9684, 0.9661 and 0.9698. The true-confound table also MIXED two
  artifacts: one row from the GC-corrected artifact (0.484) and one from
  the legacy (0.496). All now come from the GC-corrected artifact, at the
  artifact's own precision, and the guard covers the per-cohort rows, the
  true-confound value and the top-5 per-cancer AUCs and sens@99 as well,
  so they cannot drift again. `paper/METHODS_PAPER.md`,
  `docs/HARMONIZATION.md`, `scripts/benchmark_numbers.py`.
- **The readiness report carried a number no artifact contained.**
  `scripts/ultraearly_readiness.py` wrote the literal `0.499` into the
  true-confound note of the generated report, while the artifact it
  describes says 0.4836 — the 0.499 was the value of an older,
  unseeded run. So `docs/ULTRA_EARLY_READINESS_AUTO.md` and
  `results/ultraearly_readiness.json` contradicted the artifact they
  report on, and nothing noticed, because the generator wrote the
  number itself rather than reading it. It is now read from the
  artifact, and a test writes an artifact with a distinctive value and
  asserts the report follows it.
- **The methods paper's summary table still quoted the external
  validation as it was BEFORE the v2.3.0 correction** — 0.986 (WGS),
  0.87 pooled, 0.42 sensitivity (WES) — the figures from the 9-sample
  publication-1 cohort and the single-pooled harmonisation. The row now
  quotes the corrected 0.9444 / 0.708 / 0.7385 / 0.6569. The README
  carried the old publication-7 AUC 0.9915 in the same way.
- **A "~80× shift" claim with no artifact behind it.** The methods
  paper described the true-confound collapse as a drop "from 1.000 to
  0.49 (~80× shift across the operating-point flip)". No committed
  artifact contains a sensitivity for the true-confound control — it
  has only ever carried AUC fields — so the ratio was not derived from
  anything. The sentence now states the collapse the artifact does
  support: 0.9989 → 0.4836, near-perfect study separation to chance.
- **The drift guard could not see a table row that mixed two runs.**
  Every number in a row can exist somewhere in the repository while the
  row itself quotes two different artifacts, and the document-level
  check stays green because each number is present somewhere. The
  methods paper's harmonization row did exactly that — the legacy
  pooled AUC beside the GC-corrected true-confound value, citing the
  legacy artifact. The guard now pins the rows that carry the headline
  numbers: each must quote the right values, must not quote the same
  quantity from another run, and must name that run's artifact. It
  matches numbers at the precision they are written at, so "0.970"
  no longer stands in for both the GC-corrected 0.9698 and the legacy
  0.9704.
- **The fusion table had drifted to a mix of two runs.** MODEL.md and
  RESULTS.md carried a pre-regeneration copy of it: a stale tumor-naive
  AUC (0.9743), a stale mutation-only AUC (0.9242) and a stale LR-fusion
  AUC (0.9887) beside the current naive-average figure, with the paired
  t-statistic from the older run (t = 31.96 instead of 49.22), the CI
  from it, and "5 individual seeds" for a 10-seed benchmark. Every place
  that quotes the table now quotes the artifact's own numbers, and each
  row of it is pinned by the drift guard. `docs/FUSION_ISOTONIC.md`
  already had the right ones, which is how the drift was spotted.
- **The readiness document's per-cancer sensitivity row was a whole
  generation out of date.** It quoted pooled sens@99% = 0.669 with LUAD
  0.71, OV 0.90, PAAD 0.69, CRC 0.83, HCC 0.78, BRCA 0.30 and OTHER_C
  0.61 — figures from an older run of the per-cancer table, some of them
  (HCC 0.78, BRCA 0.30) far from the current values. It now quotes
  0.700 and the artifact's per-cancer figures.
- **The citation metadata was not checked at all.** `CITATION.cff` is
  what GitHub's "Cite this repository" button hands a reader, and no
  guard covered a `.cff` file. It still quoted the superseded pooled
  AUC (0.9670) and a title the manuscript does not use. It now quotes
  the manuscript's own title, the current pooled AUC, and the
  affiliation the manuscript declares, and a test holds the three in
  agreement.
- **The Simplified-Chinese README carried the same stale generation as
  the English one.** Its fusion table still quoted a pre-regeneration
  run (tumor-naive 0.9734, naive average 0.9882, LR fusion 0.9883, and
  the older t-statistic 31.96), its test badge said "70/70 passing"
  when the suite has 720 tests, and it still described the neural
  foundation encoder as the fusion layer — it is a v2.2 document that
  predates the whole v2.3.0 release, and nothing checked it. The fusion
  table and the badge now carry the current figures, and a banner says
  plainly which parts of the translation lag the English original
  rather than leaving the reader to find the discrepancy.
- **A test-count update corrupted a benchmark figure, twice.** The stated
  test count is updated by replacing a string across the documents, and
  a raw replacement rewrites any longer number that contains it: the
  external-validation figure 0.7097 became 0.7207 in one update and
  0.7237 in the next, so the methods paper said "0.7237 for the naive
  pooled figure" where the artifact says **0.7097**. The count updates
  now match only standalone integers, the figure is quoted explicitly
  in the guard's list so its exact rendering is checked, and a new test
  compares the count the documents state against the count the suite
  actually collects — so a count that drifts from reality fails the
  build instead of being agreed upon unanimously.
- **The English README's own test-coverage table was a generation out of
  date.** It said the biomedical-review suite had 24 tests (it has 37)
  and that the fast-path gate collects 233 (it collects 540, running 504
  with 23 skipped and 13 deselected). A front-page claim that
  contradicts the suite it describes. The table now carries the
  measured numbers, and a test collects each suite and the gate and
  compares — so the table cannot drift from the repository again.
- **The Simplified-Chinese README was retranslated.** It had been a v2.2
  document for several releases: it described the neural foundation
  encoder as the fusion layer, omitted the external validation, the
  leave-one-study-out measurement and the methylation null entirely, and
  carried a pre-regeneration fusion table. It is now a faithful
  translation of the current English README — the same figures, the same
  tables, the same links — and a test requires it to carry every figure
  the English original states, so a figure added to one cannot be
  forgotten in the other. The English original remains authoritative.
- **The English README contradicted itself on one figure.** The
  "What it does" bullet gave the tumor-naive + mutation fusion AUC as
  0.989 while the quick-start section gave 0.9884 ± 0.0010 for the
  same quantity; the artifact's mean rounds to 0.988, so the bullet was
  quoting the older run's 0.9886. Both are now 0.988.
- **The external validation was measured on a cohort that was too small,
  and on a harmonisation the code did not implement.** Two independent
  problems, both of which moved the headline numbers down, and both now
  corrected:
  - *The cohort.* The first fetch applied a per-publication download
    budget that took only 9 of publication 1's 60 plasma samples (8 cancer
    and 1 control) and reported that publication's external AUC as
    **0.875** from those nine. The cohort now carries all 52 of
    publication 1's labelled plasma samples, and its AUC is **0.708**
    [0.5148, 0.9012]. The held-out cohort is 322 samples (300 cancer, 22
    controls), and the fetcher gained a `--skip-existing` flag so a
    re-run only pays for the samples it is missing.
  - *The harmonisation.* `scripts/external_validation_finaledb.py`
    documents "per-study z-score harmonization, each publication
    normalized by its own feature statistics" and claims to be
    "identical to the LOSO benchmark" -- but it z-scored the held-out
    cohort as a **single pool**. With 239 case-only publication-9 samples
    in that pool, publications 1 and 7 were scored against publication 9's
    mean and variance. The harmonisation is now per study on the training
    side and per publication on the held-out side, as documented, and the
    artifact records all three harmonisations under
    `harmonization_sensitivity` because the difference is not small.

  The corrected numbers: publication 7 (Sun 2019, whole-genome) **0.9444**
  [0.8731, 1.0000], publication 1 (Snyder 2016, whole-genome) **0.708**
  [0.5148, 0.9012], pooled over the publications carrying both classes
  **0.7385** [0.6281, 0.8488] on n=83, and the case-only publication 9
  as sensitivity **0.6569** at the training-derived threshold 0.9766. The
  previous 0.9915 / 0.8819 are gone, and with them the claim that the
  assay platform was the only binding constraint: transfer between
  whole-genome studies is real but far from uniform (0.9444 against
  0.708), and the whole-exome publication transfers only partly rather
  than not at all. The manuscript, the abstract, the README, MODEL.md,
  CITATION.cff, the landing page and the readiness documents have been
  updated to the corrected numbers, and the drift guard now covers the
  per-publication and pooled external figures.

  The held-out cohort is also now demonstrably the **complete set**: of
  the eight studies FinaleDB publishes, the three that are not plasma
  cancer studies are excluded for properties of the data (publication 3 is
  urine, publication 4 has no disease labels, publication 5 has no
  plasma seqruns at all), so publications 1, 7 and 9 are every plasma
  cfDNA study FinaleDB serves that the model never trained on.
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
  run. The seeded pipeline gives **0.9766**. The artifact, the narrative
  and the two numbers in `paper/METHODS_PAPER.md` §3.2.2 have been
  corrected, the threshold is now quoted in `docs/EXTERNAL_VALIDATION.md`,
  and both numbers are on the quotable-number list so the drift guard
  covers them. (The external AUCs themselves were superseded later the
  same day, when the held-out cohort was enlarged and the harmonisation
  was corrected -- see the entry above.)
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
  limit, so the upload would have been refused; it is now 2937 with every
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
