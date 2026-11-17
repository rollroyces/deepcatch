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
  measured on 31 real samples (FinaleDB publication 7; 18 colorectal cancer,
  13 liver-transplant controls) and it separates the arms with leave-one-out
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

### Added

- **The files a software-engineering journal asks for.** The JOSS
  checklist wants a contribution guide, an author list and a code of
  conduct; none of the three was in the repository. `CONTRIBUTING.md`
  leads with the rule that is unusual here — every number must be
  traceable to a committed artifact, and the artifact must be the one
  the code actually produced — and then gives the setup, the two CI
  test commands, the dataset-intake procedure and the release policy.
  `AUTHORS` records the lead author with the ORCID iD already public
  in `CITATION.cff`; `CODE_OF_CONDUCT.md` is plain-language and scoped
  to the contribution (be accurate, keep patient data out, do not
  present a simulation result as a clinical one). None of the three
  adds a personal contact detail to the repository.

### Fixed

- **The pooled null the documents quote was the wrong seed's, at a
  precision the guard could not see.** Six documents stated the
  paired-cohort null as "0.42 ± 0.12". Those are the HIGHEST seed's
  mean and standard deviation (0.4176 and 0.1187) rounded to two
  decimals; the cohort's pooled null over all 150 permutations is
  **0.4063 ± 0.1143**. The per-seed nulls (0.3997, 0.4017, 0.4176) sit
  within a rounding unit of each other, so nothing flagged it — and at
  two decimals the number guard could not see it at all, because it
  requires at least three decimals to identify a value. A figure
  quoted at two decimals can drift freely while the guard stays green.
  The fusion artifact now carries the pooled figure at
  `protocol.measured_null` (mean, std, and the per-seed values it
  pools), the documents quote it at four decimals, and two new pins
  hold them to it. The pins are red-green checked in both directions:
  writing the old two-decimal rendering back into a document now
  FAILS the guard, where before it passed.
- **The shuffled-label control was a no-op, and now it is a control.**
  `scripts/foundation_real_smoke.py` built its "pair-broken shuffle" by
  permuting the positive index set among itself and the negative index
  set among itself, then writing 1s back to the permuted positive
  positions and 0s to the permuted negative positions. A permutation of
  a set has the same set as its image, so that reproduces `y_true`
  exactly: `y_shuf == y_true` for every seed, and the "shuffled" AUCs
  the script reported were the **honest** AUCs rather than a null. The
  fingerprint was in the committed artifact, where
  `shuffled_lr_baseline_auc_mean` was byte-identical to
  `lr_baseline_auc_mean` (both 0.904) and
  `shuffled_naive_avg_auc_mean` was byte-identical to
  `naive_avg_auc_mean` (both 0.9174) — a control that reproduces the
  real number exactly is not a control.
  The control is now a genuine label permutation, and it lives in one
  place: `src/foundation/label_control.py`. That module is the single
  implementation used by `foundation_real_smoke.py`,
  `fusion_regularized_smoke.py` and `cross_study_finallydb.py`,
  because a control implemented three times is a control that drifts —
  and this one had already drifted into a no-op in one of them.
  The honest run and the control run also went through **different**
  code paths (the control trained the frozen-encoder variant alone
  while the honest score was the 0.7/0.3 blend), so the reported
  difference mixed a change of estimator with a change of labels; both
  now go through the same function, with the same folds and the same
  model seeds, so the two numbers differ only in the label assignment.
  On the current 20-seed paired cohort the foundation scores **0.5635**
  on the true labels against **0.3464** on the permuted-label null (gap
  **+0.2171**), and the LR baseline **0.9040** against **0.3903**. The
  honest numbers are unchanged to the last digit — the defect was
  confined to the control. A new artifact field
  `foundation_auc_minus_null` carries the interpretable gap, because
  `delta_auc_normalized` divides by `real − 0.5` and so assumes a null
  at 0.5, whereas on this paired design the measured null is **0.42 ±
  0.12** (below chance, because a model fitted to permuted labels still
  recovers the within-pair axis) — a control below 0.5 has to be read
  against that measured null, not against 0.5.
  `test/test_label_control.py` guards the property structurally rather
  than by number: the control must differ from the honest number, the
  permutation must preserve the class balance, and no script may
  re-implement the permutation or restore the no-op construction.
- **The sparse-aware ablation's artifacts were four months stale, and
  the control fix exposed it.** `results/sparse_aware_linear.json`,
  `results/sparse_aware_panel.json` and `results/sparse_aware_ablation.json`
  were generated by the **pre-audit-2** code (leaky `StratifiedKFold`)
  and had never been regenerated since; their `lr_baseline_auc_mean` was
  0.9725, while the post-audit-2 code gives **0.9100** on the same
  cohort. Their per-seed arrays are the first five entries of the
  pre-audit-2 n=20 array, which is how the staleness can be seen
  directly. Because the control was broken in exactly those artifacts
  (`shuffled_lr_baseline_auc_mean` byte-identical to
  `lr_baseline_auc_mean`), they could not simply be left alone.
  Regenerating them moves the sparse-aware projection from a *marginal*
  negative to a clearly significant one: Δ AUC **−0.0440**
  [−0.0593, −0.0287], p=**0.0013**, 5/5 seeds negative, against the
  previously reported −0.0155, p=0.055. The direction is unchanged and
  the recommendation (keep `LinearProjection` as the default) is
  unchanged, but the margin is larger than the document claimed.
  `docs/SPARSE_AWARE_ABLATION.md` now carries the regenerated figures,
  and its per-seed tables were updated with them.
- **What the audit of the paired-cohort work found.** An independent
  audit of the four commits around it found three real defects, all now
  fixed. (1) The control the table called a "label-permuted control" is
  not a permutation: `scripts/foundation_real_smoke.py` writes 1 back to
  the same indices, so `y_shuf == y_true` — as
  `scripts/fusion_regularized_smoke.py` already recorded. The row is
  renamed the pair-broken control it is, and the paragraph that leaned
  on it now says what the control actually breaks and that the LR and
  naive-average channels reproduce their real AUC under it. (2) The
  driver's new fourteenth step regenerated an artifact that nothing
  compared: `scripts/reproduce_check.py` collected the quotable list and
  the extra numbers, and the table-row pins were in neither. The
  paired-cohort figures are now in the comparison (69 numbers instead of
  57), and a perturbed artifact fails it. (3) `paper/REPRODUCE.md` still
  quoted the pre-regeneration value. The audit also found the methods
  paper's abstract, the ablation document's prose and the README's
  test-file line describing the old run; they now carry the measured
  figures. The two artifacts disagreed on the standard deviation
  (population versus sample) for the same per-seed array, so the smoke
  script now uses ddof=1 and the two documents quote one number.
  One correction to the previous entry: the coverage test caught the
  missing driver step only once the table's rows were pinned to the
  artifact — it was the audit that found the table did not match the
  artifact at all.
- **The methods paper's paired-cohort table was a record, not a
  reproduction — now it is one.** Section 3.1's seven rows cited
  `results/sens_at_spec_ce_n20.json` and six of them did not match it:
  that artifact was an older run of the same protocol, and the numbers
  in the table came from an audit-2 run whose output was never
  committed. Re-running the documented commands against the current code
  (`scripts/foundation_real_smoke.py --seeds 20 --loss ce`, the
  `sens_at_spec` variant, and `scripts/sens_at_spec_ablation.py`) gives
  the current figures — panel-only 0.9175 ± 0.0135, frag-only 0.6229 ±
  0.0141, foundation 0.5635 ± 0.0111, LR baseline 0.9040 ± 0.0121,
  naive average 0.9174 ± 0.0136, permuted-label control 0.7421, delta
  −2.813 — and the three artifacts are committed in that state. (Those two
  control figures were not a null: the permutation was a no-op, so they
  were the honest numbers; the control has since been fixed and
  re-measured — see the first entry in this section.) The
  table carries them and names the artifact, seven row pins hold it
  there, and `docs/SENS_AT_SPEC_ABLATION.md`'s n=20 comparison
  follows. The driver has a fourteenth step that re-measures the
  paired-cohort baseline with the documented command, so §3.1's figures
  are compared like the rest (it skips itself, with the reason, when
  the GDC cache is absent). The guard's number parser also learned to
  read a typographic minus (U+2212): a negative figure written the way
  the manuscripts write it parsed as its positive magnitude, so it
  could never match the artifact it came from.
- **The readiness document's panel table still carried the wrong
  sensitivity, and the guard could not see it.** The document's table
  said **0.490** where the artifact says **0.460** — while another
  table in the same file already said 0.460, so the document
  contradicted itself and the document-level guard stayed green (it
  asks whether the value appears *somewhere* in the file). Two row
  pins now ask the ROW to quote the artifact's value; the guard was
  red-green checked, and the panel AUC's standard deviation was
  truncated (0.018) where the artifact says 0.0188.
- **The methods paper's §3.1 cited an artifact that does not hold its
  numbers.** Seven of its rows named `results/sens_at_spec_ce_n20.json`;
  six of them did not match it (0.915 vs 0.9869, 0.620 vs 0.7238,
  0.910 vs 0.9663, 0.915 vs 0.9869, 0.741 vs 0.6452, −3.45 vs −1.47).
  That artifact is an EARLIER run of the same protocol: its
  permuted-label entries for the LR and naive-average channels still
  repeat the unpermuted values, which the audit-2 fix in
  `scripts/foundation_real_smoke.py` says they must not. The table now
  names the artifact behind each figure — the foundation, LR-baseline
  and naive-average figures are in committed artifacts — and says
  plainly that the panel-only, frag-only and permuted-label figures
  come from an audit-2 run whose output artifact is not committed, so
  they are a record rather than a reproduction.
- **The README named a workflow that does not exist.** It said the 13
  deselected slow tests "run in dedicated GitHub Actions workflows
  (`foundation-real-smoke`, `validate.yml`)"; no such workflow has
  ever existed in this repository. Eleven of the thirteen live in
  `src/foundation/test_integration.py` and do run in the deep-learning
  job of `validate.yml`, which does not deselect them; the other two
  need the companion repository's feature cache and run only where it
  is present. The PDF guard's figure list is now read out of
  `paper.tex` rather than written out, so a new `\includegraphics` is
  pixel-compared instead of merely counted.
- **The last pointers to the earlier manuscript are gone, and the
  acknowledgements the team page promises do not exist.** The project
  metadata's "Paper" URL, the README's research-paper link, the methods
  paper's BibTeX note and the CI trigger list still pointed at
  `paper/PAPER.md` — the earlier panel-only revision — and the
  compilation workflow rebuilt the submission PDF whenever that file
  changed. They now point at `paper/paper.tex`, the manuscript the PDF
  is compiled from. `TEAM.md` claimed the data providers and the
  reviewers "are acknowledged in paper/PAPER.md and BENCHMARK.md
  citations"; neither document has an acknowledgements section, so the
  claim was false. It now says what is true — the data are named in the
  manuscripts' data-availability statements — and the journal checklist
  carries the missing section as an open item, because whom to name
  there, and how to disclose the AI-assisted audits, is the author's
  call.
- **The readiness document quoted the panel's sensitivity at two
  decimals, and the guard could not see it.**
  `ULTRA_EARLY_READINESS.md` said "AUC 0.921 / sens@99% = 0.49" where
  the artifact says **0.460** — the figure every other document
  quotes. The number guard requires at least three decimals, because
  "0.49" cannot identify a value, so the figure was invisible to it.
  It is written as 0.460 now and named in the quotable list, so a later
  rewrite of that document fails the build; the readiness document had
  never been named for the panel's AUC either, and now it is.
- **The long-form manuscript still carried a September status line.**
  `paper/METHODS_PAPER.md` described itself as a 2026-09-23 draft whose
  §3 came from the legacy sweep (`results/cross_study_finallydb.json`),
  while its numbers come from the GC-corrected one; it now names the
  artifacts its figures actually come from, and records that the author
  line is a placeholder to settle before submission. `paper/README.md`
  called that manuscript "the JOSS / PLOS-style article" while the
  venue map sends JOSS to the root `paper.md`; it also quoted the
  submission PDF's byte size, which changes on every rebuild.
- **Nothing compared the submission PDF with the manuscript.**
  `paper/biorxiv_submission_v2.3.0.pdf` is the file bioRxiv receives,
  and it is a committed binary: a PDF that was not rebuilt after the
  manuscript changed would carry every figure the manuscript has since
  changed, and the drift guard — which checks the .tex against the
  artifacts — would not see it. This repository has already had a
  submission PDF whose Figure 1 plotted a different curve from the one
  the text described (v2.2.0). The new guard reads the PDF's own text
  and requires every figure *and* count the manuscript states to be in
  it, and compares the embedded images against the committed figure
  files pixel for pixel, so a figure regenerated without rebuilding the
  PDF fails the build. The list is derived from `paper/paper.tex`, so
  the guard cannot drift from the manuscript. `pypdf` is a test-only
  dependency; the CI job installs it.
  `test/test_submission_pdf_matches_the_manuscript.py`.
- **The upload package listed a supplementary that belongs to a different
  paper.** `paper/supplementary.tex` is the deprecated v1.0 draft's
  supplementary: it describes figures S1–S10 and tables S1–S5 for the
  multi-modal product framing the project moved away from (CTC/miRNA
  detectors, contrastive variant calling, MAML, longitudinal CET
  monitoring), and it ships no image files. The bioRxiv upload list put
  it next to `paper.tex`, so the upload would have described experiments
  the manuscript does not contain. The form no longer lists it, and
  `paper/README.md` says what it actually is.
  `paper/README.md`, `paper/BIORXIV_SUBMISSION.md`.
- **Five documents still called the earlier manuscript the current one.**
  The labelling pass covered `docs/README.md` and the docs site, but
  `docs/PORTFOLIO_ROADMAP.md`, `paper/METHODS_PAPER.md` (twice),
  `docs/JOURNAL_SUBMISSION.md` and `MODEL.md` still sent a reader to
  `paper/PAPER.md` — the earlier panel-only revision — as the bioRxiv
  manuscript, and `MODEL.md` quoted line counts ("242 lines", "177
  lines") that had already drifted and that no test guards. They now
  name `paper/paper.tex` as the manuscript, and the unguarded line
  counts are gone.
- **The sibling pipeline's drifting headline was measured rather than
  assumed.** `ZENODO_SETUP.md` told the author to seed the PCA in
  `cfdna-fragmentomics-pipeline` before tagging it, because that
  repository's headline moved on every run. Two runs on an isolated
  worktree at its HEAD (its uncommitted work untouched) put the
  mechanism to it: scikit-learn's PCA picks the randomised SVD when
  `n_components < 0.8 * min(n_samples, n_features)`, and with PCA(200)
  on the 627-sample cohort that is 200 < 400, so the subspace was
  random. The unseeded run gave 0.9744 ± 0.0032 and raised an lbfgs
  "failed to converge" warning; the seeded run gives 0.9770 ± 0.0013
  with no warning — identical to DeepCatch's own legacy pooled number.
  The one-line fix and the sibling's document updates are prepared on
  its local branch `fix/seed-pca` (not pushed), and the seven further
  unseeded PCA calls in that repository were seeded too: the sens@spec
  and calibration figures now give a pooled OOF AUC of 0.977018 ±
  0.0013, bit-identical to the cross-study headline, while the control
  that fits no PCA does not move at all. Two findings came out of it —
  channel reweighting is now exactly the baseline (Δ +0.0000), so the
  "-0.0331" the sibling's report used to print was unseeded noise, and
  that report had been carrying hardcoded numbers that had already
  drifted out of step with the table above them.
  `ZENODO_SETUP.md`.
- **The panel-LLR headline — the abstract's first figure — was not
  regenerated by anything.** The readiness verdict depends on it and
  the driver never re-ran it, so a regression there would have been
  invisible. Re-running it against the committed inputs reproduces the
  artifact exactly (820 numeric fields; only the elapsed time
  differs), so it is now step 11 of the driver, compared like the rest,
  and the panel AUC is part of the quotable-number list.
  `paper/REPRODUCE.sh`, `scripts/benchmark_numbers.py`.
- **The drift guard named fewer documents than quote the numbers.** The
  audit walked every quoted figure and every current-facing document and
  found the copies the guard did not name: `CITATION.cff` — the file
  GitHub's "Cite this repository" button hands a reader — quoted seven
  headline figures and was named by none of them, and the LOSO transfer
  directions, the in-domain pair and the per-cancer rows were quoted by
  a dozen documents the guard never looked at. Every document that
  quotes a guarded figure is now named, so a copy of one cannot drift
  unnoticed.
  `scripts/benchmark_numbers.py`.
- **The with-motifs benchmark never applied the GC correction, and the
  paper compared it against a GC-corrected headline anyway.** The
  loader chain took the `--include-motifs` branch first, so
  `--gc-correction` was silently unapplied on that path — the run was
  legacy-footed underneath while its artifact recorded
  `gc_correction: true` and a `c100_gc_corrected` feature set. The
  methods paper therefore reported "+0.0101 over the GC-corrected
  baseline", which was mostly the GC correction and not the motif
  block: the like-for-like lift on the legacy set was +0.0015. The GC
  correction and the motif block are independent feature blocks, so the
  corrected matrix with the motif block appended IS that configuration;
  it is now what the code runs, and the measured lift on the
  GC-corrected set is **+0.0077** (0.9684 → 0.9761). The artifact
  records the EFFECTIVE configuration (`gc_correction_applied`), and
  `--motif-cohort-mode strict` now refuses to run with the correction
  instead of silently ignoring it.
  `scripts/load5_with_motifs.py`, `scripts/cross_study_finallydb.py`,
  `results/cross_study_finallydb_with_motifs.json`, `paper/METHODS_PAPER.md`.
- **Four readers of a labels file did not agree.** The vocabulary now
  lives in one module (`src/fragmentomics/labels.py`) and every reader
  resolves the label AND the column layout the same way: by the header
  when there is one, by the width when there is not. They used to
  disagree on the same file — the package's reader took the disease
  class for a headerless four-column file, the two script readers
  skipped a three-column file and returned an empty cohort, and the
  package kept a private copy of the vocabulary the commit claimed was
  shared. A test now reads the same file through all four readers and
  requires one answer, and another test walks the collaborator's merge
  path end to end.
  `src/fragmentomics/labels.py`, `src/fragmentomics/tumor_naive_adapter.py`,
  `scripts/cross_study_finallydb.py`,
  `scripts/external_validation_finaledb.py`,
  `scripts/adapter_local_cohort.py`, `test/test_label_vocabulary.py`.
- **A superseded figure was banned globally although a current figure
  renders the same digits.** 0.9915 is the old publication-7 external
  AUC *and* the decision curve's naive-average AUC; 0.9670 is the OV
  per-cancer DeLong CI upper bound; 0.9701 one of the pooled per-seed
  AUCs; 0.9753 the OTHER_C per-cancer AUC. A global literal ban on
  those would forbid a correct document, so each is now banned only
  where it would be a regression — the documents that quote the
  quantity it replaced.
  `scripts/benchmark_numbers.py`,
  `test/test_paper_numbers_match_artifacts.py`.
- **Four stale figures in a status document, and a false requirement
  that held it in place.** `paper/METHODS_PAPER_STATUS.md` quoted
  Jiang 0.9791, Cristiano 0.9693, pooled 0.9738 and no-harmonize 0.9650
  where the artifact holds 0.9789, 0.9704, 0.9770 and 0.9678. It also
  satisfied a guard requirement for the Bayesian model average (0.9650)
  only by that same stale number, so correcting it turned the guard red;
  the requirement was dropped instead. `NEXT_STEPS.md`'s panel AUC
  0.9215 became 0.921, and that headline — the first figure in the
  preprint's abstract — is now pinned to every document that quotes it.
  `paper/METHODS_PAPER_STATUS.md`, `NEXT_STEPS.md`,
  `scripts/benchmark_numbers.py`.
- **The with-motifs step could not be reproduced off this machine.** The
  cached 4-mer motif features are tracked in neither repository, so
  without them the step writes a *different* benchmark; it is now
  skipped with an explanation, `paper/REPRODUCE.md` carries the twelve
  steps and names the two that need inputs the repository does not
  contain, and the readiness verdict is compared rather than merely
  regenerated.
  `paper/REPRODUCE.sh`, `paper/REPRODUCE.md`,
  `scripts/reproduce_check.py`.
- **The collaborator adapter could read a header row as data.** It
  detected a header only when the first column was spelled `sample_id`,
  so the canonical FinaleDB header (`sample`, `disease_class`, `label`,
  `study`) fell through to the data branch and became a phantom sample:
  `sample_id` "sample", `disease_class` "disease_class", `label`
  "label" — which the cross-study loader then counts as a non-cancer
  control. The header is now recognised by its column names rather than
  by one exact first-column spelling, and a test covers both layouts
  plus the headerless one.
  `scripts/adapter_local_cohort.py`, `test/test_adapter_local_cohort.py`.
- **The preprint mis-stated the methylation cohort split.** Six
  documents said the 31 real samples were "17 colorectal cancer and 14
  liver-transplant controls"; the artifact holds **18 colorectal cancer
  and 13 liver-transplant controls** (and the manuscript contradicted
  itself two paragraphs earlier, where it had the right numbers). The
  composition was never pinned, so no check could catch it.
  `paper/paper.tex`, `paper/METHODS_PAPER.md`, `README.md`,
  `README.zh-CN.md`, `ULTRA_EARLY_READINESS.md`,
  `docs/METHYLATION_PIPELINE.md`.
- **Two headline readers turned any label that was not literally
  "cancer" into a healthy control.** `scripts/cross_study_finallydb.py`
  and `scripts/external_validation_finaledb.py` resolved the label cell
  by exact string equality, so "Cancer", "tumor" or "control" became
  controls — silently, in the scripts that write the headline
  cross-study AUC and the external validation. The package's own
  resolver read the same file differently, which is what exposed it. The
  vocabulary now lives in one place (`src/fragmentomics/labels.py`), an
  unrecognised label is an error rather than an invented control, and
  the collaborator adapter raises the schema error its interface
  promises. 24 tests cover it. The fix is behaviour-preserving on the
  real data: re-running the GC-corrected sweep reproduces all 483
  numeric fields of the committed artifact (pooled AUC
  0.9684259954921112) and all 329 of the per-cancer artifact, so only a
  label that WAS mis-filed changes behaviour.
  `src/fragmentomics/labels.py`, `scripts/cross_study_finallydb.py`,
  `scripts/external_validation_finaledb.py`,
  `scripts/adapter_local_cohort.py`, `test/test_label_vocabulary.py`.
- **The bioRxiv form checker did not check anything.** It reported
  which artifact-backed numbers the form contained, against pairs that
  still named the superseded external AUCs (0.9915, 0.8819) and the
  wrong artifact for the methylation null — so a stale form would have
  passed. Every number the form quotes is now resolved against the
  artifact it claims to come from, and a figure the v2.3.0 regeneration
  superseded is refused outright. The list of superseded figures moved
  to `scripts/benchmark_numbers.py` so the form and the drift guard
  cannot disagree about what is stale.
  `scripts/check_biorxiv_form.py`, `scripts/benchmark_numbers.py`.
- **The release's central claim was unguarded and unreproducible.**
  The regularized late fusion (0.9825) and the Bayesian model average
  (0.9650) are the goal-2 headline, quoted by ten documents, pinned by
  none, and regenerated by no driver step. They are now pinned, and the
  driver regenerates the artifact. The with-motifs benchmark and the
  shuffled-label null (0.4460) are pinned too, and the driver runs the
  with-motifs sweep; `quoted_artifacts()` now counts the table-row pins,
  so an artifact pinned only by a row is no longer invisible to the
  driver's coverage test. Both new steps were executed against the
  committed artifacts: the with-motifs sweep reproduces all 483 numeric
  fields of `cross_study_finallydb_with_motifs.json` exactly, and the
  late-fusion smoke reproduces all 218 of
  `fusion_regularized_smoke.json`.
  `scripts/benchmark_numbers.py`, `paper/REPRODUCE.sh`.
- **The readiness verdict was regenerated after the proof step and
  never compared.** It is derived from the artifacts the run just
  regenerated, and ~20 documents quote it. It now runs before the
  proof step and its numbers are compared like the rest.
  `paper/REPRODUCE.sh`, `scripts/reproduce_check.py`.
- **The driver could not run on a clean machine.** The cross-study,
  LOSO and external-validation steps import the sibling pipeline's
  loader, so a reproducer without that checkout died with "No module
  named 'honest_benchmark'". The prerequisite is now documented and
  pre-flighted with a message that says what to clone.
  `paper/REPRODUCE.sh`.
- **Stale claims in source docstrings and one wrong citation.** The
  adapter's docstring still quoted the sibling pipeline's superseded
  0.9753, two modules named output files they do not write, and
  RESULTS.md attributed the methylation null to an artifact that does
  not contain it. `RESULTS.md` §2 also stated the true-confound pair
  from the August snapshot without saying so; it now carries the
  current pair as well.
- **A committed generator could not run from a clean checkout.** The
  CADD GDC cohort builder resolved its paths from a hard-coded
  `/Users/hermes/deepcatch`, so it only worked on the machine that
  produced the report. Its paths are now derived from the repository
  and overridable, and its report says plainly what is still missing
  (the 400-file MAF cache and the CADD matching TSV).
  `scripts/build_gdc_validation_cohort.py`,
  `scripts/cadd_alphamissense_combined.py`,
  `docs/CADD_GDC_VALIDATION.md`.
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
  exactly. The corrected control's null on the paired design is 0.4063 ± 0.1143,
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
- **The text submitted to bioRxiv was not the paper.**
  `docs/PREPRINT_ABSTRACT.md` — the document `docs/PATH_TO_IMPACT.md`
  tells the author to submit — was written separately and had drifted:
  it quoted an older abstract (panel AUC 0.922 instead of 0.921,
  Sens@95% 0.600 instead of 0.770, the Fisher score 0.849 instead of
  0.834), carried none of the external-validation or
  leave-one-study-out results, had an older title, and named
  "Royce et al." on a single-author paper. The abstract and the PDF
  would have been different papers. The document is now generated from
  `paper/paper.tex` by `scripts/generate_preprint_abstract.py`, with the
  LaTeX rendered to the plain text bioRxiv's field accepts, and the
  generator's vocabulary is closed: an unknown command raises rather
  than silently dropping a symbol. A test fails the build if the
  committed document is not what the generator produces.
  `scripts/generate_preprint_abstract.py`,
  `test/test_preprint_abstract_generator.py`,
  `docs/PREPRINT_ABSTRACT.md`.
- **The abstract fitted bioRxiv's character limit again.** It grew with the
  v2.3.0 results to 3119 plain-text characters against bioRxiv's 3000
  limit, so the upload would have been refused; it is now 2940 with every
  number intact. A test fails the build if the abstract exceeds 3000
  characters, or creeps within 50 of the limit, and
  `scripts/check_biorxiv_form.py` fails if the upload form quotes a
  different length than the manuscript actually has. The length is
  measured on the text that is actually pasted — the abstract rendered
  out of its LaTeX, paragraph breaks included — rather than on a
  normalisation that deleted the mathematics and so counted a string
  nobody would submit.
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
