# Changelog

All notable changes to DeepCatch are recorded here. The project follows
[semantic versioning](https://semver.org/); the Zenodo DOI for each
release is minted from the corresponding GitHub release.

## [2.3.0] — 2026-11-23

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

- **The release is archived and citable: 10.5281/zenodo.24453044.**
  Zenodo's GitHub integration is enabled for this repository, and the
  v2.3.0 release was archived and published automatically: the record
  is <https://zenodo.org/records/24453044> (v2.3.0, MIT, the source
  archive `rollroyces/deepcatch-v2.3.0.zip`), the DOI is registered
  with DataCite, and the GitHub release carries the submission PDF.
  The DOI is declared in `CITATION.cff`, in `MODEL.md`'s BibTeX entry,
  in the methods paper's availability section and in the preprint's
  data-availability section, and
  `test/test_release_doi_is_consistent.py` checks that all four agree,
  that none of them carries a placeholder, and that no DOI appears in
  an availability statement without being declared. The v2.2.0 release
  predates the integration and still carries no DOI.
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

### Added

- **Both manuscripts now carry the AI-use disclosure.** bioRxiv expects
  it, arXiv asks for it, and PLOS's Ethical Publishing Practice requires
  it ("must be clearly reported in a dedicated section of the Methods,
  or in the Acknowledgements section") and treats its absence as
  misrepresentation of methods and contributions. The paragraph names
  the tools the repository's own records name -- `TEAM.md` §1.2: the
  OpenAI / Anthropic / Google providers behind the multi-persona audit
  agents -- plus MiniMax, which ran the current audit-and-fix round.
  It says what the tools did (wrote and reviewed analysis code, ran the
  scripts, performed the audits recorded in `AUDIT_2_FINDINGS.md`,
  drafted and edited parts of the manuscript) and what they did not do
  (generate, alter or interpret any primary data, figure or statistical
  result), and that every number is produced by the committed code from
  the committed artifact. The disclosure took the submission PDF from 8
  pages to 9, and every document that states the page count was
  updated with it.
- **The GSE124600 cohort is now measured, not quoted.** The GEO survey
  described GSE124600 from the series description; this measures it.
  `results/geo_gse124600_cohort.json` records what the three deposited
  supplementary files actually contain -- 205,250 CpG loci per matrix,
  **274 plasma sample columns** across the two plasma files (268 + 6) and
  66 tissue columns (33 colorectal cancer tissue + 33 adjacent mucosa) --
  with each file's byte count and sha256, the stage-label distribution
  counted from the GEO sample records, and the series record's own design
  quote and the quality-control arithmetic behind the cohort size.
  `test/test_geo_gse124600_cohort_is_consistent.py` keeps that artifact and
  the documents that quote the cohort honest: it checks the artifact's own
  arithmetic, pins the cohort size and the stage distribution, and finds
  the documents by scanning the tree for the accession, so a document that
  quotes GSE124600 without the artifact's numbers fails.
- **The affiliation is now checked as a string, not only as "is an
  employer named".** The author guard compared the name string and the
  ORCID across six surfaces, and asked whether one of them named an
  employer -- but it never compared the affiliation itself, and the
  surfaces did not agree: the manuscript's author block said
  `Independent Researcher`, `CITATION.cff` said "Independent Researcher
  (no institutional affiliation; verified email domain umich.edu)", the
  generated abstract said `Independent Researcher`, and the submission
  sheet's form values said "Independent Researcher, Hong Kong" -- a
  locality that appears nowhere in the PDF the form uploads. That matters
  because the form's Author field and the paper's title page are read
  side by side, and `docs/SCHOLAR_VISIBILITY.md` §1 says Google Scholar
  matches a paper to a profile by the affiliation string. The guard now
  compares the declared affiliation across the manuscript, the citation
  metadata, the generated abstract, the JOSS paper and the submission
  sheet's two form rows. It ignores a trailing parenthetical -- the
  citation metadata qualifies its value that way, and that qualifier is
  metadata about the declaration rather than a different declaration --
  and nothing else, so a locality added to one surface only is caught.
  The three new checks move the suite from 1024 to 1027 tests, so every
  document that states the count moves with them and the submission PDF
  and the arXiv bundle are rebuilt.
- **Zenodo was searched for the first time, and it holds a second external
  cohort this repository can score.** The goal asked for processed cfDNA
  features that are public and need no IRB, and the survey had covered
  GEO over 99 series but only two known Zenodo records. Eight queries
  against the Zenodo records API (`scripts/zenodo_cfdna_survey.py`,
  artifact `results/zenodo_cfdna_survey.json`) brought back **142
  distinct records** in the first pages, and one of them is exactly what
  the goal names: **10.5281/zenodo.21269430**, per-sample hg38 fragment
  BED files for 12 plasma samples -- 6 clear-cell renal cell carcinoma
  patients and 6 matched healthy controls -- with a sample metadata table
  carrying the pT stage.

  `scripts/external_cohort_zenodo_ccrcc.py` extracts the repository's
  five benchmark channels from those files and scores them with the model
  trained on FinaleDB publications 6 + 8, importing the protocol from
  `scripts/external_validation_finaledb.py` rather than re-implementing
  it -- the same training-derived threshold (0.976602) comes out of both.
  The result is in `results/external_cohort_zenodo_ccrcc.json`: **AUC
  0.9444**, DeLong 95% CI [0.8156, 1.0], **exact permutation p =
  0.0043** over all 924 label assignments, and the ranking is stable
  across the three harmonisations (0.9444 / 0.9722 / 0.9722).

  This is a second external family -- a different cancer type, a
  different laboratory, a different data source -- and three honest
  limits travel with the number: all six patients are **pT3a (stage
  III)** and two have distant metastasis, so it is a late-stage cohort;
  the assay is about 30x shallower than the training cohort's; and the
  **absolute threshold does not travel** (at the training cohort's
  99%-specificity threshold the cohort's specificity is 0.5, because the
  whole cohort sits high on the score axis -- the controls' median score
  is 0.869 against 0.003 for the training cohort's). An AUC is a rank
  statistic, which is why the rank is the honest number to quote for a
  held-out cohort and a probability threshold taken from another cohort
  is not. `test/test_zenodo_cfdna_survey_and_ccrcc_cohort.py` keeps the
  artifacts, the documents and the pins together.

  A fourth limit was measured afterwards and travels with the number:
  **the cohort's depth is a confound**. The two arms were sequenced to
  different depths (median 1,106,290 fragments per patient against
  931,526 per control), and **depth alone separates them at AUC
  0.8333** -- so the question the headline has to answer is whether the
  model's separation IS the depth. It is not: regressing the score on
  the depth and scoring what is left gives **AUC 0.9722** (partial
  Spearman of score against label given depth 0.6777, against 0.7725
  uncontrolled). The confound is real, is disclosed, and the separation
  survives it. The artifact's `depth_control` block holds the numbers.

  Reading the deposited fragments needed one change to the repository's
  own code: the files are standard BED files and carry a header line
  (`#chrom start end name mapq strand`), which `read_fragments` could
  not parse, so it now skips a leading `#` line. The FinaleDB files do
  not have one, so the training cohort's features are unchanged.
- **The stage-labelled cohort's methylation matrix was opened, and the
  answer is measured rather than planned.** `ULTRA_EARLY_READINESS.md`
  §1.3 called the stage breakdown the one readiness gap and said the
  public cohort with stage labels (GSE124600) "would need a methylation
  feature extractor and its own model". The matrices were downloaded and
  opened instead (`scripts/geo_gse124600_methylation_cohort.py`,
  artifact `results/methylation_cohort_geo_gse124600.json`):
  * the matrix is a methylated-CpG **CALL** matrix with normalised
    counts -- 91.1% of the 205,250 x 274 cells are exactly 0 and the
    maximum is 6,949.25, so a zero means "not called methylated in this
    sample", not "0% methylated";
  * the per-locus summaries are **anti-correlated** with cancer (AUC
    0.27-0.41); the only class signal is the number of CpGs called
    methylated per sample (15,637 in the healthy controls against 20,130
    in the cancer plasma, AUC 0.7688, permutation p < 0.0001) -- a
    burden/coverage measure, the same shape as this repository's own
    FinaleMe null, now reproduced on an independent and larger deposit
    that carries stage labels;
  * the stage breakdown is a **negative** answer: the called-CpG count
    does not rise with stage (Spearman rho 0.0961, p = 0.26 over the
    142 cancer samples; the per-stage AUCs against the controls are
    0.7130 / 0.8223 / 0.7188 / 0.8288). The one public cohort with
    per-sample stage labels does not deliver a Stage I result either, so
    the readiness gap is closed with a measurement rather than a plan --
    and it stays open.

  The checks live in `test/test_geo_gse124600_cohort_is_consistent.py`,
  next to the cohort they belong to, and they pin the finding, the stage
  answer, and the rule that a document quoting the finding has to carry
  the negative answer with it.
- **The Zenodo survey now reads the deposits' cohort tables, not just
  their descriptions.** A survey that stops at the record's description
  reports a labelled cohort as unlabelled. SAUNA's deposit
  (10.5281/zenodo.15097509) carries a `metadata.tsv` with each sample's
  group, coverage and **tumour fraction**: 9 samples (7 CRC, 2 benign)
  at 15-33x coverage, one patient sampled three times with the tumour
  fraction rising 0.2903 -> 0.3271 -> 0.4029. That is the closest
  public thing to this repository's early-detection question -- a
  fragmentomics cohort with a GRADED cancer label at more than the
  training cohort's depth, and a longitudinal series in public data. It is not
  extracted here (13.1 GB, about 148 M fragments per sample -- 296 M
  reads, counted from the files) and no
  number from it is quoted yet; it is recorded with its labels, depths
  and tumour fractions in `results/zenodo_cfdna_survey.json` and
  described in `docs/PUBLIC_DATA_SOURCES.md` §6.2. The survey artifact
  now records which deposits are EXPECTED to carry a cohort table, and
  the guard fails when one is declared and not read -- the earlier note
  in the survey that said the labels were not in the record was wrong,
  and the guard is what makes that kind of claim checkable.
- **A third external cohort, the first stage-graded *fragment* cohort, and
  a negative result with a measured cause.** The Cleveland Clinic PHTS
  screening cohort (Zenodo 10.5281/zenodo.15109495; Liu et al., *Cell
  Reports Medicine* 2024) carries 99 plasma samples with the depositor's
  own stage column, its pre/post split and its diagnosis dates. It has
  been unpacked, converted and scored with the same protocol as the ccRCC
  cohort (`scripts/prepare_zenodo_phts_cohort.py`,
  `scripts/external_cohort_zenodo_phts.py`; artifacts
  `results/zenodo_phts_prep.json`, `results/external_cohort_zenodo_phts.json`):
  * the deposit ships **BEDPE**, not fragment BED. The extractor reads
    columns 1-3, which in a BEDPE are the first mate's coordinates, so
    fed as-is it reads a "fragment" one read long (about 146 bp instead
    of 167 bp) and the fragment-size distribution collapses onto the
    read-length peak -- silently, with every channel wrong by the same
    factor. The preparation script converts each member to the fragment
    form and refuses to report success unless every sample shows the
    nucleosome ladder. All 99 pass; 302,801,515 fragments.
  * the deposit is mapped to **hg19** (chr1 runs to 249,239,706; hg38's
    ends at 248,956,422), while the feature grid is hg38 -- the grid the
    training cohort's features use, so it is the right one. The report
    records the cost: 0.433% of the fragments are dropped by it.
  * the "pre-diagnosis" arm is **not** a lead-time arm. The depositor's
    `date_diff` column puts every one of those nine draws 1 to 7 days
    before the diagnosis, so the longest lead time anywhere in the cohort
    is 7 days and the tumour was already there. An arm called
    "pre-diagnosis" is not an early-detection arm.
  * the tumour-naive model does **not** separate the classes: AUC
    **0.4506** [0.3343, 0.5669], sampled permutation p = 0.80, and the
    null survives the harmonisation choice -- three variants that let the
    cohort normalise itself give 0.4506 / 0.4563 / 0.4563, and the
    control that normalises it by the **training** cohort's statistics
    gives 0.5388 [0.4232, 0.6544], still including 0.5. The stage
    breakdown is flat -- the pre-diagnosis arm (n=9) at chance, stage 0
    0.556, stages 0-1 0.4119, stages 2-3 0.575, all staged 0.4621, no
    stage monotonicity (rho 0.172, p 0.297). This is the first
    stage-graded fragment cohort the repository has scored, so the stage
    question now has an answer rather than a plan.
  * and the null is **not about the transfer**: a model fitted and
    cross-validated inside the cohort reaches AUC 0.5465, CI [0.4319,
    0.6612], and read one channel at a time every interval over the whole
    cancer arm includes 0.5 (0.4367-0.5759). This cohort has no class
    signal in these channels to find, so it is not a valid external test
    in either direction -- it shows where the channels stop working
    (mostly subclinical disease in a screening population), not that the
    model fails.
  * the effect the depositor published does show up, in the arm with the
    highest tumour burden: read one channel at a time, the
    second-neoplasm arm separates (the short/long ratio is higher there
    -- AUC **0.6455** [0.5008, 0.7901], in the direction the depositor's
    own paper reports -- the fragments are shorter, and the
    di-nucleosome proportion is lower). Over the whole cancer arm the
    same channel's interval includes 0.5 again (0.5714 [0.4574, 0.6854]):
    the effect is in the highest-burden arm and not in the cohort.

  The checks live in `test/test_zenodo_cfdna_survey_and_ccrcc_cohort.py`,
  next to the cohorts they belong to, and they pin the finding, the null,
  the lead time, and the rule that a document quoting the numbers has to
  carry the limits with them.
- **A cohort with enough tumour burden does separate -- and the model
  reads the burden.** The SAUNA deposit (Zenodo 10.5281/zenodo.15097509)
  is nine colorectal plasma samples at 15-33x, mapped to hg38, with
  the depositor's own per-sample tumour fraction: 0.0049 to 0.4029, and
  one patient sampled three times as it rises 0.2903 -> 0.3271 -> 0.4029.
  It is the closest public thing to this repository's early-detection
  question -- the training cohort's own cancer type, the build the
  feature grid uses, a graded dose rather than a label -- and it has now
  been unpacked, verified and scored with the same protocol as the other
  held-out cohorts (`scripts/prepare_zenodo_sauna_cohort.py`,
  `scripts/external_cohort_zenodo_sauna.py`; artifacts
  `results/zenodo_sauna_prep.json`, `results/external_cohort_zenodo_sauna.json`):
  * the tumour-naive model's score tracks the tumour fraction: Spearman
    rho **0.6611** over the nine (sampled permutation p = 0.0641) and
    **0.8154** within the seven CRC samples (p = **0.0355**);
  * the raw channel carries more of it than the model does -- the
    short/long ratio's rho against the tumour fraction is 0.7333 overall
    (p = 0.0329) and **0.9643** within the CRC samples (p = **0.00305**);
  * all four harmonisation variants agree on the sign and the size
    (rho 0.66-0.75), so this is not an artefact of how the cohort is
    normalised;
  * the one patient sampled three times has a short/long ratio that
    rises monotonically (0.20709 -> 0.23023 -> 0.26068) and a median
    fragment length that falls (297 -> 295 -> 291 bp). The model's own
    score is saturated at 1.000000 on all three, so it cannot be
    monotone -- the channel moves under it, the score is at its ceiling;
  * the binary contrast is thin and is not the readout: CRC vs benign is
    0.7143 [0.3101, 1.0] on seven against two.

  Against the PHTS cohort's null (above), this is the other half of the
  picture: in a cohort with enough tumour burden the channels read it,
  and in a screening population whose malignant samples are mostly
  subclinical there is nothing to read. Nine samples and one patient
  sampled three times is still a small measurement, and the library's
  fragment-size shape differs from the training cohort's (short/long
  ratio 0.0849 against 0.2182), so the gradient is the honest number
  here and the absolute threshold is not.

  The extractor had to change to make the measurement possible at all:
  these files hold up to 307 million fragments each, and reading one
  whole -- with the chromosome column as Python strings -- needs more
  memory than the machine has. It now reads chunk by chunk, which gives
  exactly the same numbers (re-extracting three ccRCC samples with the
  new code reproduces the committed feature files bit for bit).
- **The end-motif tables are scored, and the confound control is
  what they show.** Zenodo 10.5281/zenodo.17253504 publishes per-sample
  4-mer end-motif frequencies for 223 cfDNA samples -- 87 healthy
  controls from two studies and 136 cancer samples (breast 50, NSCLC 56,
  metastatic colorectal 30) in two library preparations. It is the only
  deposit in the survey carrying **the channel the PHTS paper's own
  finding rests on** -- fragment end motifs -- with cancer/healthy
  labels and per-sample values, and it is 2.3 MB rather than gigabytes
  (`scripts/prepare_zenodo_endmotif_cohort.py`,
  `scripts/external_cohort_zenodo_endmotif.py`; artifacts
  `results/zenodo_endmotif_prep.json`,
  `results/external_cohort_zenodo_endmotif.json`).
  * The tables' motif rows are in the same lexicographic ACGT order as
    the repository's own motif index, which the preparation script
    verifies rather than assumes.
  * The numbers look strong: the published claim's statistic (the
    fraction of 4-mer end-motif counts beginning CC or CA) is higher in
    the cancer arm -- median 0.2637 against 0.2466, AUC 0.8207,
    permutation p = 5e-05 -- and a cross-validated model inside the
    cohort separates cancer from healthy at AUC 0.9958 [0.9884, 1.0].
  * And the confound control rules them out as cancer signals. The
    healthy arm comes from two studies and the cancer arm from three
    others, so a library-level difference is not separable from a class
    difference here -- and the control says how large it is: **the same
    256 features separate BC from NSCLC, two cancer groups from
    different studies, at AUC 0.9996**, higher than the
    cancer-versus-healthy separation, and the two healthy studies differ
    on the CC-/CA-end fraction at AUC 0.7057. The per-group medians say
    the same thing (HI25 0.1586, mCRC 0.1407, HI37 0.2477, BC 0.2645,
    NSCLC 0.2678): the three groups that sit high are the ones whose
    libraries sit high.
  * So none of the class AUCs in that artifact can be read as a cancer
    signal. An AUC without that control is not a result, and the
    artifact says so.
  * The repository's own motif model, trained on the 98 training-cohort
    samples that carry a motif vector (48 healthy, 50 hepatocellular
    carcinoma), lands at **AUC 0.3413** here -- inverted. That is what a
    mismatch of definitions looks like: the deposit says the motifs
    come from the 5' fragment boundary but does not say whether both
    strands' ends were counted per fragment, which is what the
    repository's extraction does. The ambiguity is stated in the
    artifact, and the inversion is the evidence for it.

### Fixed

- **The repository's evidence that bioRxiv's affiliation rule "is not
  enforced in practice" was measuring the wrong thing.** The rule was read
  off a statistic -- 115 of 5,480 January-2026 preprints are
  single-authored -- and taken to show that a submission without an
  institutional affiliation would get through. It does not: a
  single-authored preprint can perfectly well carry an institutional
  affiliation (one PhD student is a single author), so the statistic never
  measured affiliation-less submissions at all. A submission without one
  was rejected on 2026-11-25 with exactly that reason -- *"bioRxiv
  requires authors to add an institutional affiliation ... an organization
  that can provide oversight of research activities"* (manuscript ID
  BIORXIV/2026/725086; nothing was posted, no DOI was minted, and the
  record does not appear in bioRxiv's API). The row that recorded the weak
  evidence now records the rejection instead, and the local submission
  checklist states the consequence: the affiliation is mandatory for the
  bioRxiv path, so the open question is no longer *whether* to declare one
  but *which* organisation, and only if the one declared is true. The
  arXiv fallback, which asks for no institutional affiliation, is no
  longer a formality.
- **The rejected submission was also the wrong manuscript.** Its title was
  the v2.0-era "DeepCatch: Performance-Weighted Multi-Modal Fusion with
  Cumulative Evidence Tracking for Pan-Cancer Detection from Cell-Free
  DNA", not the current title, so the next attempt has to upload the
  committed `paper/biorxiv_submission_v2.3.0.pdf` and use the current
  title in the form. `docs/JOURNAL_SUBMISSION.md` now says so next to
  the file it names.
- **The arXiv fallback was documented as if arXiv had no rules of its
  own.** arXiv is the only venue that asks for no institutional
  affiliation, so after the bioRxiv rejection it is the route that does
  not depend on the affiliation being settled -- but it has three
  requirements that would have stopped the submission, and none of them
  was written down. Its metadata rules say "abstracts longer than 1920
  characters will not be accepted" and the manuscript's abstract is
  2,940, which bioRxiv accepts and arXiv will not; the submission guide
  says "Authors are expected to self-submit" and that a third-party
  submission "without prior authorization will be removed", so on arXiv
  the author presses "Submit Article"; and the same guide warns that a
  new submitter "may be required to find endorsements", which is a real
  risk for a first submission. The authors-field policy adds the sharpest
  one: a claimed affiliation "should be current in the conventional
  sense" and misrepresentation is "possible grounds for immediate and
  permanent suspension" -- so choosing arXiv defers the affiliation
  question rather than removing it. `docs/JOURNAL_SUBMISSION.md` gains
  an arXiv section that states each requirement in arXiv's own words,
  and `docs/SCHOLAR_VISIBILITY.md` section 3.2 records that an arXiv
  preprint's DOI is deposited to **DataCite** (the `10.48550` prefix), not
  Crossref, so the registry to pick when the work is added to the ORCID
  record depends on which server posted it. The DOI checker learned the
  same thing: it sent every identifier to Crossref except a Zenodo one,
  so a real arXiv DOI -- which arXiv deposits to DataCite -- would have
  been reported as unregistered the first time one entered the tree.
  `docs/SCHOLAR_VISIBILITY.md` §3.2 also records the cheaper route for
  the works list: ORCID's **Add DOI** fills the record in from the
  registry, needs no trusted-party grant, and if the ORCID iD is given
  at submission the publisher passes it on and the work arrives by
  itself.
- **The JOSS paper was in the format JOSS stopped accepting.** JOSS
  changed its submission format in late 2025, and its own review bot now
  checks the paper file for a *State of the field*, a *Software design*, a
  *Research impact statement* and an *AI usage disclosure* section: a
  paper in the old format is rejected in pre-review, and this
  repository's guard was still listing the old eight-section format. The
  manuscript now carries all twelve sections. The four new ones are
  written from the repository's own facts -- the open-data landscape and
  what is still missing from it, the design principle that a number
  cannot exist without an artifact behind it, the impact stated honestly
  (developer-led use, the companion-pipeline integration, public
  artifacts with DOIs on both ends, and no documented external use yet),
  and the AI disclosure with the three elements JOSS's policy asks for.
  The guard requires the twelve headings, and it now checks that the
  disclosure is a real disclosure: it must name a model provider and
  confirm the human review, because JOSS treats an incomplete disclosure
  as a research-integrity matter rather than a formatting one. Both
  failure modes are negative-tested.
- **The JOSS pre-review gates were never written down.** JOSS screens a
  submission against four gates before review -- public development
  history of more than six months, demonstrated research impact,
  open-source practices for a single-author project, and iterative
  development over time -- and a paper that fails any one of them is
  desk-rejected with a note saying what is missing.
  `docs/JOURNAL_SUBMISSION.md` now records each gate with the evidence
  for this repository: the repository was created on 2026-04-28, so the
  six-month threshold was cleared on 2026-10-28, and the commit history
  spans six months of 2026; the impact gate is marked partly met and
  names the preprint as the thing that moves it, because a preprint
  server is explicitly not a previous publication and the bioRxiv
  preprint can therefore be posted before, during or after the JOSS
  review; the open-source-practices gate needed a stated support and
  governance expectation, which `CONTRIBUTING.md` now has. The scope
  rule that the paper must not focus on new research results is
  addressed by a sentence in the Summary that says what the measurements
  are and where the findings themselves are discussed.
- **The bioRxiv upload was blocked on an account, not on the ORCID.** The
  ORCID authorization *did* succeed -- the callback returned to
  `submit.biorxiv.org/?orcidcode=0003&orcidid=0009-0008-9113-769X` -- and
  the page then refused the author: *"No account is associated with
  orcid, please login and associate orcid via the edit personal
  information page."* The ORCID side can succeed while the bioRxiv side
  still refuses, because the iD has to be associated with an existing
  bioRxiv account from that account's Personal Info area. The working
  order is email -> "Check for Account" -> create or sign in ->
  associate the iD -> submit, and the registration form's first field is
  the one thing only the author can type. `docs/JOURNAL_SUBMISSION.md`
  and `docs/SCHOLAR_VISIBILITY.md` section 3.1 now carry that sequence,
  so the next attempt starts from the email field instead of from the
  ORCID button.
- **The affiliation decision was documented twice, and the two documents
  lettered the same two options in OPPOSITE order.**
  `docs/JOURNAL_SUBMISSION.md` called "(a)" the stale ORCID entry and
  "(b)" the current one; §7 of the local submission checklist called
  "(a)" the current employment and "(b)" the stale one. A reader who had
  read one and then consulted the other would have applied the wrong
  ending -- and the wrong ending is the difference between declaring an
  institution and not declaring it. Both now use the same letters, each
  option names the exact string it ends with ("University of Michigan,
  Ann Arbor, Michigan, USA" or "Independent Researcher"), and neither is
  presented as a default: which of the two is true is the author's fact,
  not the repository's. The evidence behind the decision was re-checked
  against the ORCID public API on 2026-11-24 -- one open-ended
  University of Michigan employment, no start date, no end date, and zero
  works -- unchanged since 2026-11-17.
- **Two of the surfaces that declare the affiliation were not guarded.**
  `test/test_author_affiliation_is_consistent.py` compared the manuscript,
  `CITATION.cff`, the generated abstract, the JOSS paper, the model card
  and the submission sheet's two form rows -- but `paper/METHODS_PAPER.md`
  and `paper/PAPER.md` declare the affiliation in prose on their author's
  footnote line and were not compared at all. The affiliation decision
  moves them, so a partial edit could have left one of them naming a
  different institution than the paper it accompanies. The guard now
  requires the declared string to appear in both author blocks (the prose
  is searched rather than parsed: a manuscript that names another
  employer, or a truncated one, does not contain it), and both cases are
  negative-tested. `docs/SCHOLAR_VISIBILITY.md` §4 now lists every
  surface that declares the affiliation and every surface that describes
  the author's situation, marking the ones the guard compares, so the
  change is a checklist rather than a search.
- **The GEO survey was too broad, and it had missed the category the
  goal asked about.** `docs/PUBLIC_DATA_SOURCES.md` §4 listed four
  accessions and concluded that GEO hosts "bigWig coverage tracks or
  5hmC/ATAC peaks, not per-sample fragment files". Re-run over **99
  cfDNA-related GEO series** (2026-11-24), that is true for fragments
  and false for methylation: GEO hosts per-sample **methylation
  matrices**, and for plasma cfDNA cohorts it hosts them **with
  per-sample stage labels**. GSE124600 (blood-based colorectal cancer
  detection) carries 142 cancer and 132 control plasma samples (274 in
  total) with `crc grade: Stage_I`, `Stage_II`, `Stage_III` and
  `Stage_IV` on the cancer plasma, plus per-sample
  methylation matrices -- the first public cohort carrying the stage
  labels the readiness document says the open-data cohort lacks. It is
  a *methylation* cohort, so it cannot be scored with this repository's
  fragmentomics model; it would need a methylation feature extractor
  and its own model. GEO still holds no fragment-length matrix for a
  cfDNA cancer cohort (the one recent fragment study deposits coverage
  tracks and uses serum), which is why the held-out cohort comes from
  FinaleDB, whose fragment files are served directly.
- **The GEO survey's cohort size was wrong, and the matrices turned out to
  be complete.** The entry above gave the wrong number of control plasma
  samples for GSE124600. The series record says plasma from 147 CRC
  patients and 134 controls, of which 5 patients and 2 controls failed
  quality control, so the deposited cohort is **142 cancer + 132 control =
  274 plasma samples**; the survey had also left `Stage_IV` (5 samples) out
  of the stage set it listed. Downloading the three deposited files and
  counting their columns shows **205,250 CpG
  loci** in every matrix and **274 plasma sample columns** in total -- 268
  in `sumP_..._dup1.5_...` and 6 in `sumP_..._nodup1.5_...` -- which
  covers every plasma sample in the series, plus 66 tissue columns (33
  colorectal cancer tissue + 33 adjacent mucosa). The measurements, the
  file hashes and the provenance are in
  `results/geo_gse124600_cohort.json`, and
  `test/test_geo_gse124600_cohort_is_consistent.py` keeps the documents and
  the artifact from drifting apart again.
  *(The earlier note that the 1.9 MB supplementary file was a truncated
  download was wrong too: macOS `zcat` cannot read a `.gz` file and fails
  with `can't stat`, which looks exactly like a truncated transfer.
  `gzip -dc` reads the file, and its byte count matches the FTP
  `Content-Length` exactly.)*
- **The submission sheet's form values and the Scholar document described
  an affiliation the paper does not carry.** The bioRxiv and arXiv sheets
  told the author to paste "Yu Ching Lam (Independent Researcher, Hong
  Kong)" into the Author field while the uploaded manuscript declares
  `Independent Researcher` -- so the record and the PDF it describes would
  have carried two different affiliations for the same author. And
  `docs/SCHOLAR_VISIBILITY.md` §4 said the repository says
  "Independent Researcher" "in the manuscripts and in `CITATION.cff`, and
  the Zenodo record archived that same string", which was wrong about the
  citation metadata and about the Zenodo record: both carry the longer
  qualified string. That document's whole purpose is to make every
  surface an index reads say the same affiliation, so describing the
  strings wrongly was the worst kind of error for it to contain. The
  form values now carry the manuscript's string, §4 now lists what each
  surface declares -- including the docs bylines' locality, which is
  contact information rather than a declared affiliation -- and the
  locality stays where it belongs.
- **The count guard had a blind spot, and a stale count was sitting in
  it.** `RESULTS.md` states the standalone test directory's size twice:
  in the Section 8 status row, and in the paragraph under it ("`pytest
  test/` collects **N** (the standalone test files)"). The guard checked
  the row, so the sentence below it stayed at 790 while the row above it
  said 796 -- a reader who took the number from the paragraph got a
  count nine tests behind the one in the table above. The guard now
  checks both places, and the sentence says 799 with the date it was
  measured.
- **The docs site stopped deploying because the Pages action looked for
  the artifact before the artifact list had caught up.** The deploy step
  finds the artifact by NAME through the artifact list API, and that
  list can lag a second or two behind the upload's own success
  message: on 2026-11-23 the upload finalised the artifact and printed
  its ID, and the lookup ran 19 ms later and reported "Found 0
  artifact(s)" for an artifact that existed in the same run. The
  workflow now waits for the artifact to become visible before it
  deploys it, and the deploys are serialised so that two pushes cannot
  deploy at once.
- **Two manuscripts carried no ORCID, and the author identity was only
  checked in three of the six places it appears.** Google Scholar
  matches a paper to a profile by the NAME string, the affiliation and
  the email -- it does not read an ORCID to find an author -- and it
  follows the ORCID to link the work to the profile afterwards. So a
  name that differs by a hyphen or an initial splits the author's works
  across two profiles, invisibly. `paper/PAPER.md` and
  `paper/METHODS_PAPER.md` carried no ORCID at all, and the guard only
  compared three author blocks. Both manuscripts now carry the ORCID,
  and the guard compares the name string and the ORCID across all six
  surfaces an index reads (the LaTeX manuscript, the JOSS paper, both
  markdown manuscripts, `CITATION.cff` and `AUTHORS`).
- **`docs/SCHOLAR_VISIBILITY.md`: what the three indexes need, and
  what only the author can do.** It records the measured state of each
  source (the ORCID name matches every surface, its works list is
  empty, the Zenodo record is complete with the ORCID, the preprint
  does not exist yet) and the ordered steps: post the preprint, add the
  two DOIs to the ORCID record with Search & Link against DataCite,
  and claim exactly one Scholar profile with the same affiliation the
  preprint carries.
- **Three documents described the submission PDF with a figure count the
  PDF does not have.** The PDF embeds three figures -- the panel-LLR
  performance figure, the pooled ROC and the sensitivity-at-specificity
  panel -- and the venue map said "2 figures embedded"; the paper
  README hedged with "2 figures - 3 embedded images"; the next-steps
  note said "2 figures". A count in prose about a binary artefact is
  exactly what drifts: nothing regenerates it. `docs/JOURNAL_SUBMISSION.md`,
  `paper/README.md` and `NEXT_STEPS.md` now state 8 pages and 3 figures,
  and `test/test_submission_description_matches_the_pdf.py` reads the
  page count from the committed PDF and the figure count from the
  manuscript's own `\includegraphics` calls, so a prose count that
  stops matching the artefact fails the build.
- **The external AUCs were quoted without saying what stage the held-out
  cohort was at.** 0.9444, 0.708 and the pooled 0.7385 are the
  project's headline generalisation number, and nothing in the
  repository said what stage those samples were at -- so a reader could
  take them for an early-detection result. FinaleDB's sample record
  carries no stage field, so the stage was read from the records that
  describe each publication: **publication 1's GEO extract protocol
  states that all 52 of its cancer samples come from individuals with a
  clinical diagnosis of Stage IV cancer**, and **publication 9 is a
  metastatic cohort** by its paper's own title ("high concordance with
  metastatic tumors"). Publication 7's per-sample annotation lives in
  its supplementary tables, which PMC serves behind a proof-of-work and
  reCAPTCHA gate; that row is recorded as unresolved rather than
  guessed, and a collaborator with journal access can finish it. The
  consequence is stated wherever the numbers are read -- the
  external-validation document, the preprint (results and
  limitations), the methods paper (§3.2.2 and a new limitation), the
  readiness document and the public-data survey: the external AUCs
  measure **late-stage discrimination** and carry no statement about
  Stage I detection. The characterisation lives in
  `results/external_cohort_stage_annotation.json`, with the record it
  came from, the sentence it says and the date it was read, and
  `test/test_external_cohort_stage_is_stated_where_it_matters.py`
  keeps the documents and the artifact from drifting apart.
- **The model card pointed at a commit 134 commits behind the release.**
  `MODEL.md` said "DeepCatch v2.3.0 (commit `7f9613b` on `main`)", and
  `7f9613b` is the 2026-11-16 paper commit -- the release is the tag
  `v2.3.0` at `170a46a`. A commit hash in a version pointer is a
  promise that cannot be kept: `main` moves, the hash stays, and a
  reader who checks it out gets a tree whose numbers are not the
  numbers the card states. The version is named by its release tag
  (and its DOI) now, which cannot drift, and a new guard fails if a
  current-facing document points at a version with a commit hash. The
  changelog and the dated audit reports are exempt: naming the old
  commit is the point of those files.
- **The DOI check reported the release's own DOI as unregistered.**
  `scripts/check_dois.py` trims the trailing punctuation of a captured
  DOI but not markdown emphasis, so a DOI written in bold
  (**10.5281/zenodo.24453044**) was read as an identifier ending in
  asterisks -- one that resolves to nothing -- and the check that exists
  to catch dead links failed on the DOI the release had just minted.
  The trimming drops the emphasis markers now, the extraction is a
  function of its own so that it can be tested, and seven forms are
  pinned: bare, bold, bold with a full stop, italic, a trailing full
  stop, a markdown link and a LaTeX `\url{}`.
- **The list of documents a quoted figure must appear in was a dict
  literal, and a repeated key silently dropped half of it.** The guard
  reads `EXTRA_DOCS` to learn which documents quote which figure, and
  seven keys in it were written twice -- the second time for a document
  added later. A Python dict keeps the LAST value for a repeated key,
  so the earlier list was gone: **31 document/figure protections were
  lost**, including the documents the superseded-figure guard is scoped
  to, so a superseded number could have come back into them unnoticed.
  The entries are now PAIRS merged by union, so a repeated key adds to
  the protection instead of replacing it, and
  `test/test_every_figure_copy_is_protected.py` fails on any dict
  literal with a duplicate key anywhere in the repository -- the
  duplicate is silently lossy and nothing about it looks like a
  mistake when reading the source.
- **The guard checked that named documents quote a figure, never that
  unnamed ones do.** A document that quotes a figure but was never
  named is a copy that can drift silently: the artifact is regenerated,
  the named documents are updated, and this one keeps the old number.
  That gap is now mechanical: the same test walks every tracked
  document, and a document that quotes a pinned figure must be named
  for it or carry an excuse that says what its number is instead. The
  excuses are checked in both directions -- one that stops being needed,
  stops quoting the figure, or loses its reason fails.
- **Three documents were named for a figure they do not quote.**
  `USAGE.md`, `docs/ALPHAGENOME_PROXY_RESULTS.md` and
  `docs/CADD_PER_SUBGROUP_LLR.md` were named for the ovarian
  per-cancer AUC because they contain the digits `0.989` -- in a stale
  fusion rendering, an AlphaGenome proxy AUC and a subgroup AUC
  respectively. The name was satisfied by the wrong number, which hid
  the first of them: `USAGE.md` carried **0.989** for the cross-repo
  fusion, the older run's value, while the artifact holds 0.9884.
  `USAGE.md` and `METHYLATION_PROJECT.md` now quote the artifact's
  value, and the two documents that quote something else are excused
  instead of named.
- **A document quoted the companion repository's figures without naming
  the commit it read them from.** `METHYLATION_PROJECT.md` states the
  companion pipeline's cross-study AUC (0.974-0.978) as a range; every
  other document that quotes the companion now names its commit, and a
  range is a copy that goes wrong the moment the other repository's
  numbers move. It names `eb1529e` now.
- **The guards read the suite's size with a three-digit pattern, and the
  suite passed a thousand tests.** `(\d{3})` matches "998" and not
  "1002", so the moment the count crossed four digits the README badge
  looked as if it stated no count at all -- the count guard failed with
  "the README no longer states the suite's test count", and the arXiv
  bundle and submission-PDF checks were matching counts that were no
  longer there. The patterns read three or four digits now, in the four
  places that read a count. The suite is **1002 tests** (816 core + 186
  deep-learning), and every document that states it has been updated
  together -- the counts are moved by an anchored replacement, never by
  a raw find-and-replace, which is how a research figure was corrupted
  once already.
- **The documented command for rebuilding the held-out cohort selected
  nothing.** `scripts/fetch_external_cohort.py` caps how many samples each
  publication contributes, and its quota defaulted to empty while
  `--max-per-publication` defaulted to 0 -- so the command the
  reproduction driver tells a reader to run selected no samples and
  stopped with "nothing selected". The `--help` example
  (`9=240,7=31,1=12`) was worse: it ran, and it selected a cohort that
  is not the published one (publication 1 contributed 12 samples where
  the committed cohort has 52). The default is now the quota that
  reproduces the committed cohort, and a guard compares it with the
  per-publication counts in the published artifact -- offline, so it needs
  neither the API nor the S3 bucket.
- **The FinaleDB API is up, and the repository said it was down.**
  `docs/PRETRAINING.md`, `paper/REPRODUCE.md` and `paper/REPRODUCE.sh`
  all stated that the API had been degraded since 2026-09-21. The
  repository's own probe (`scripts/_verify_publication_readiness.py`)
  reports the verdict `ok` as of 2026-11-22 -- every endpoint answers
  200 and `/seqrun` enumerates 2,579 runs -- and a real fetch downloaded
  and extracted samples from the public S3 bucket. The S3 probe reports
  `partial`: the fragment objects it checks are reachable, some other
  keys answer 403. The claims are corrected, with the date they were
  measured.
- **Three documents quoted the headline figures without being pinned to
  them.** `paper/METHODS_PAPER_STATUS.md` restates the methods paper's
  §3.2 rows, `wiki/Jiang-4mer-Validation.md` quotes the Jiang
  publication's per-publication AUC, and `docs/V3_DESIGN.md` quotes the
  panel LLR anchor -- and none of the three was named in the list of
  documents a quoted number must appear in. A regenerated artifact could
  have moved them without anything noticing. They are named now, each pin
  read in its document first: a document that merely ROUNDS to the same
  number (a range's endpoint, or a threshold written as a target) is
  not named, because the quote test matches at three decimals and such a
  document would pass the check without ever quoting the figure.
- **The submission PDF was never checked for the things a preprint
  server rejects.** There is a guard for whether the PDF is the
  manuscript (the figures, the numbers); there was none for whether the
  file is one a server can process. `test/test_submission_pdf_is_submittable.py`
  now requires every font in the PDF to be embedded -- an unembedded one
  is substituted when the PDF is rendered, so the paper comes out in the
  wrong typeface -- the PDF not to be encrypted (bioRxiv's submission
  help lists locked or encrypted PDFs among the formats it does not
  support), its text to be extractable (a scanned page cannot be
  proof-read or indexed), and its size to be inside arXiv's 10 MB limit.
  The PDF this repository builds is 8 pages, 265 KB, with all sixteen
  fonts embedded.
- **The arXiv source bundle was never checked against arXiv's own rules.**
  arXiv does not compile a submission with the author's toolchain: it
  unpacks the archive and compiles it with its own TeX Live, with no
  shell escape and no external programs. A bundle that works locally can
  still come back as "no output produced", and the author learns it from
  an error screen rather than from a check. The bundle this repository
  builds is deliberately plain -- `article`, `hyperref`, `geometry`,
  `amsmath`, `booktabs`, `graphicx`, and the bibliography written inline
  -- and `test/test_arxiv_bundle_is_submittable.py` now keeps it that
  way: no package arXiv cannot run, no shell escape, no absolute path,
  every figure it names a member of the bundle, no build leftovers, no
  file above arXiv's 10 MB limit, and a compiled `.bbl` whenever the
  manuscript reads its references from a bibliography file.
- **A broken arXiv bundle would have been uploaded as if it were fine.**
  The paper workflow verified the bundle with `continue-on-error: true`,
  so the job stayed green and the broken archive was uploaded anyway --
  the one outcome the step's own comment said it was there to prevent.
  The step now fails the job (the submission PDF is uploaded before it,
  so the PDF artefact is unaffected), and the paper job runs the bundle
  guards, which it did not: the bundle exists only in that job, so the
  checks had no job to run in.
- **The three demo scripts were one directory short of the repository
  root.** Each derives the root from its own location -- `docs/demo/terminal`
  is three levels down, so the path needs three `..` -- and each used two.
  `docs/demo/terminal/../..` is `docs/`, so the install demo installed the
  package inside `docs/`, ran the venv's python against `docs/test/...`, and
  printed its success anyway: the install failure was swallowed by `|| true`.
  A script that cannot find the repository root is now a guard failure, and
  the check is mechanical -- it runs the same `cd` the script runs and
  compares the result with the root.
- **The install demo's entry-point step found nothing, every time.** It
  filtered the console scripts by their *value*
  (`src.fragmentomics.train_tumor_naive:main`) while the name is what
  carries `deepcatch`, so the step printed "CLI entry points (0)" on every
  run while the three scripts `pyproject.toml` declares went unmentioned. A
  guard now runs the demo's own filter against the declared entry points and
  fails if it does not find them.
- **A CPU-only install could not collect the suite at all.** `torch` is an
  optional dependency -- the `api` extra provides it -- and four modules
  needed it at import time: `src/foundation/test_integration.py` imported
  it outright, and `src/tissue_deconv/model.py` and
  `src/priming/response_predictor.py` set the name to `None` and then
  subclassed `nn.Module`, so they raised `AttributeError` at import. A
  collection error is not a skip: pytest interrupted the whole run, so a
  reader who installed the runtime and ran the documented test command got
  no test results at all, not even the CPU ones. The modules now skip
  (`pytest.importorskip`) or raise a clear ImportError that the packages'
  `__init__` catches, so the packages import without torch and the CPU
  tests still run. A guard collects the suite with `import torch` blocked
  and fails if anything errors.
- **The JOSS paper told a reviewer to install the runtime and then run
  the suite.** `paper.md` -- the file a JOSS reviewer reads -- said "install
  with `pip install -e .`" and, a few lines later, "the suite is 978
  tests". The base install has no pytest, so the reviewer who follows the
  instruction and runs the suite gets "No module named pytest", on the
  very claim the paper makes about reproducibility. The paragraph now
  names the dev extra for the test suite and the api extra for the
  deep-learning modules, and lists the dependencies the runtime install
  actually brings (it also said "the CPU fragmentomics path needs only
  NumPy, SciPy, scikit-learn, pandas and statsmodels": `statsmodels` is
  in no dependency set at all, and `matplotlib`, `seaborn`, `requests` and
  `tqdm` are). A document that quotes the suite's size must now also say
  how to get a test runner.
- **The repository's stub torch made every "is torch installed?" probe
  lie.** `src/multimodal_fusion/advanced_fusion.py` installs a placeholder
  torch in `sys.modules` when the real one is missing, so that the CPU
  fusion path can import -- but from then on `import torch` succeeded for
  the whole process, in every module. A probe that only asked "is torch
  importable?" therefore believed torch was installed, and the deep
  learning tests ran against a placeholder instead of skipping. The stub
  now carries a marker, the probes honour it, and the stub also provides
  `torch.utils.data` (without it, `from torch.utils.data import ...` failed
  with "No module named 'torch.utils'; 'torch' is not a package"). A CPU
  install now runs the gate as 685 passed + 54 skipped + 2 deselected, and
  the README states that outcome next to the one with torch.
- **The stub torch leaked into the session and broke the real torch's own
  imports.** The stub installs `torch`, `torch.nn`, `torch.nn.functional`
  and `torch.nn.utils` into `sys.modules` when torch is missing.
  `test/test_biomedical_review_fixes.py` simulates a torch-free
  environment by putting a fake torch there and restoring the real one
  afterwards -- but it restored only three of the keys, so the stub's
  `torch.nn.utils` survived the test. When the stub was extended with
  `torch.utils` (so that `from torch.utils.data import ...` can work in
  a CPU-only install), the leak grew into a crash: the real torch's
  `torch.save` does `from torch.utils.serialization import config`, and
  with the stub's `torch.utils` in `sys.modules` -- which is not a
  package -- that raised "No module named 'torch.utils'; 'torch' is not
  a package", and two checkpoint tests failed in the reproduction
  driver's gate. The stub now lists every key it installs, the test
  restores all of them, and the deep-learning CI job runs the
  reproduction driver's gate set in one session: the DL job ran the
  foundation module on its own and the core job ran those tests on their
  own, so the interaction was covered by nothing.
- **Every quickstart installed the package without a test runner, then ran
  the test suite.** `pytest` is a development dependency: `pyproject.toml`
  puts it in the `dev` extra, not in the runtime `dependencies`. The
  README's install block, the Chinese README, `CONTRIBUTING.md`,
  `docs/USER_GUIDE.md`, `docs/REPRODUCE.md`, `USAGE.md`, the docs site's
  quickstart and the install demo script all installed the base
  distribution and then ran `pytest`, which fails on a clean virtualenv
  with `No module named pytest` before a single test runs. The demo was
  worse twice over: it installed with the ambient `pip` rather than the
  venv's own, and then ran the venv's python. `paper/REPRODUCE.sh` and
  `scripts/reproduce_all.sh` now check that pytest is importable and say
  what to install.
  `test/test_documented_install_provides_what_it_runs.py` requires that a
  documented pytest run is preceded by an install that provides pytest
  (or by a check that says how to get it), that `pytest` really is a
  development dependency, that a document's "installed automatically by
  `pip install -e .`" list is the runtime dependency list, that a
  documented command's claimed outcome is what the command selects, and
  that the two READMEs state the same test numbers. Red-green verified:
  six of its tests fail against the documents as they were.
- **The wiki promised packages the base install does not bring.**
  `wiki/Getting-Started.md` listed `torch` and `statsmodels` under "Core
  Python dependencies (installed automatically by `pip install -e .`)".
  Neither is a runtime dependency -- `torch` is in the `api` extra and
  `statsmodels` is in no extra at all -- while the list omitted
  `requests` and `tqdm`, which are. The list is the runtime list now, and
  the optional packages are listed separately with the extra that provides
  each. The READMEs made the same claim as "Minimum dependencies" and
  told the reader to install by hand the packages the runtime install
  already brings, including `statsmodels` -- which is optional (the GC
  correction falls back to a quadratic regression without it, and the
  published cross-study numbers use the committed GC artifact rather than
  a smoothing computed at run time). A list that claims to be what the
  install provides can now only name packages that the install provides.
- **Three test counts in the READMEs had drifted.** The English README
  claimed "15 passed" for the publication-readiness suite, which collects
  16, and "702 passed, 19 skipped, 13 deselected" for the fast-path gate,
  which adds up to 734 while the gate collects 768. The Chinese README
  had drifted further: its badge said 919 where the English one said
  953, its table said the gate collects 734 where the English one said
  768, and its tree said 46 test files where the English one said 57. The
  count guard compared the documents with each other, so it could not see
  a drift inside one of them -- and the Chinese README was not in its
  list at all. It is now, and a documented command's claimed outcome is
  compared with what the command actually selects.
- **The methods paper's per-cancer sensitivity columns were not pinned to
  the artifact that holds them.** Section 3.3 quotes each cancer's
  sens@95/98/99 from `results/per_cancer_sens_at_spec.json`, and the
  artifact carries all three for all five cancers -- but no guard required
  the rows to quote it. A blanket search-and-replace of the suite's test
  count rewrote OV's sens@95 from 0.964 to 0.967 in the table and in the
  paragraph under it, and every guard stayed green. The five rows are
  pinned now, and the paragraph's number is checked at the place it sits
  (the paragraph is one long line, so a per-line presence check cannot
  tell which of its numbers was meant). The rows' AUC column and
  bootstrap CI are still outside the guard: they come from the paper's
  5-seed run, and the committed artifact holds a single-seed AUC and its
  own bootstrap draw, so pinning them there would fail on numbers that
  are not wrong.

- **Five modules documented their maths in docstrings Python mangled.**
  The LaTeX was written in ordinary strings, so `\frac` was read as a form
  feed followed by "rac", and `\mid`, `\sum` and `\epsilon` produced
  invalid-escape warnings on every import; Python has said for years that
  such sequences "will not work in the future". The docstrings are raw
  strings now, so the rendered text is the LaTeX that was written.
  `test/test_no_invalid_escape_sequences.py` compiles every tracked module
  and fails on an invalid escape sequence, so the mangling cannot come
  back. It was found by installing the repository from a fresh clone --
  the first time the warnings were visible, because the working copy had
  been importing them long before.
- **The top-level runner's own artifact was four months stale.**
  `results/final_cross_validated_results.json` was written on 2026-04-30 by
  a quick-mode run of the script as it stood then, so running the
  documented command today produced a different set of numbers from the
  artifact the repository carries. The artifact is now the one the current
  code produces (a full-mode run, 5 seeds), and its own metadata says so.
  Nothing quoted it -- the paper's figures come from the cross-study
  artifacts that `paper/REPRODUCE.sh` regenerates -- but a reader who ran
  `bash RUN_ALL.sh` and compared the output with the committed file had no
  way to tell which was which.
- **The guide's runtime claim was stale and named a flag that does not
  exist.** It said the full pipeline takes "~20-30 minutes" on a MacBook
  Pro and that "the `--demo` mode runs in under 2 minutes" -- `--demo` is
  not a flag the runner accepts (it now errors on an unknown flag), and
  the measured full run takes 7.7 minutes on an Apple M4, with the
  reduced-scale run at 3.2 minutes. The sentence now gives the measured
  number, the machine and the date.
- **Three DOIs in the repository resolved to nothing.** A DOI that does not
  resolve is a dead link in a document a reader follows, and two of them
  looked exactly like working DOIs: the Mathios 2021 lung-cancer
  fragmentome paper had two digits transposed (`...-21394-w` for
  `...-24994-w`), the JOSS manuscript still quoted the PNAS paper whose
  identifier is unregistered (the other manuscripts describe the Cancer
  Discovery paper, whose identifier resolves), and the Finaletoolkit
  preprint's identifier could not be located. The JOSS manuscript also gave
  the Newman 2014 paper another paper's pages.
  `scripts/check_dois.py` now extracts every DOI from the tracked tree and
  asks the registry that issued it -- Crossref for the journal literature,
  DataCite for the Zenodo records, because a Zenodo DOI is not in Crossref
  -- and `test/test_dois_in_the_tree_are_registered.py` runs it, skipping
  when a registry is unreachable rather than failing. All 42 identifiers in
  the tree resolve.
- **The LaTeX manuscript gave a paper the wrong pages.** Its inline
  bibliography gave the Jiang 2020 end-motif paper the pages "692-681" --
  reversed, and not the published range; the methods manuscript's
  bibliography file had the correct "664-673". Two manuscripts, two
  versions of the same paper. Checked against the DOI registry, the other
  eleven entries agree with the registry. The citation guard now compares
  the two bibliographies for every work they share, so a page range (or a
  volume, or an article number) cannot drift apart again.
- **The methods manuscript cited numbered references and listed none.** The
  numbers it used were also inconsistent with each other: FinaleMe was
  `[8]` in one sentence and `[11]` in four others, and the held-out
  publications were `[18]`, `[19]`, `[20]` against a bibliography file
  that ordered them differently. A reader of a candidate PLOS or
  Bioinformatics manuscript saw numbers pointing at nothing.
  `paper/METHODS_PAPER.md` now carries a reference section, numbered by
  first appearance, and every citation points at the right work (the
  Galleri row also said "Klein 2021" for a paper by Liu and colleagues).
  The works the manuscript cited but the bibliography file did not carry
  were added to it.
- **Three DOIs in the bibliography resolved to nothing.** Checked against
  the DOI registry: the Heitzer review's identifier was unassigned (now
  `10.1038/s41576-018-0071-5`), the Jiang 2020 entry named the PNAS paper
  with an identifier that does not resolve (the manuscript describes the
  Cancer Discovery paper -- the 222-patient cohort -- whose identifier
  does), and the Finaletoolkit preprint's identifier could not be located
  at all (its DOI is now dropped rather than kept as a link to nothing).
  `test/test_manuscript_citations_resolve.py` checks that every numbered
  citation resolves to a listed reference, that every listed reference is
  cited, that the LaTeX manuscript's `\\cite{}` keys resolve against its
  own inline bibliography, and that the compiled PDF carries no `[?]`
  marker. `paper/README.md` no longer claims the orphan
  `references_final.bib` is the manuscript's bibliography, and the upload
  checklist no longer tells the author to upload it.
- **A release could cite its DOI in one place and not in another.**
  `ZENODO_SETUP.md` says the minted DOI goes into four places: the
  citation metadata, the model card's BibTeX entry, the manuscript's
  data and code availability section and the methods manuscript's
  availability section. Four files, one hand -- and nothing checked that
  they agreed. A half-done release is worse than no release: a reader
  who follows the citation metadata and a reader who follows the
  manuscript get different answers, and one of them is a link to
  nothing.
  `test/test_release_doi_is_consistent.py` now requires a declared DOI
  to appear in all four, forbids a placeholder DOI in any of them (the
  repository once carried
  `10.5281/zenodo.to-be-assigned-on-first-release`), and flags a DOI
  that appears in an availability statement but in neither the citation
  metadata nor the two public archives this repository cites. Red-green
  verified for all three cases -- and the red-green check is what caught
  a bug in the guard itself: its DOI pattern had a stray bracket, so it
  never matched anything. A guard that cannot fire is worse than no
  guard, because it looks like a check.
- **The quickstart notebook did not run.** The README's second demo is a
  recording of `notebooks/deepcatch_quickstart.ipynb`, and its histogram
  cell wrote to a relative path -- `figures/fsi_hist.png` -- while no
  directory of that name exists in the repository, so the cell failed
  with a FileNotFoundError wherever the notebook was run from. The
  committed outputs came from a machine where the directory happened to
  exist, so nothing could see it: the notebook is not a test file, and no
  test ran it. The cell now writes to a scratch directory (never into the
  repository), and `test/test_quickstart_notebook_runs.py` executes every
  computing cell of the notebook with the working directory the notebook
  assumes, so the demo cannot rot again unnoticed. Red-green verified:
  with the old path restored, the test fails on cell 5 with exactly the
  FileNotFoundError a reader would have seen.
- **The guide's table of API response fields was checked against nothing.**
  The user guide tells a reader what `POST /predict` returns --
  `classical.mds_simpson`, `neural.probability`, `explanation.top_motifs`
  and the three accepted method names -- and the API's schema declares them
  by hand, so renaming a field in `api/schemas.py` changed the API
  without changing the sentence that describes it.
  `test/test_api_contract_matches_the_guide.py` now reads the field list
  out of the guide itself (so it cannot drift into a second copy) and
  requires every component to be a declared field of the schema, and every
  method the guide shows to be a member of `PredictionMethod`. The field
  lookup matches the *declaration*, not the name anywhere in the file: the
  schema's docstrings carry JSON examples with the same names, and a
  looser search passed even after the field had been renamed -- which the
  red-green check caught.
- **The collaborator intake template was outside every guard.** The issue
  template a potential collaborator opens names the files they are asked to
  look at (`src/fragmentomics/finale_features.py`,
  `scripts/external_validation_finaledb.py`,
  `scripts/cross_study_loso.py`), and it was the one user-facing document
  no guard read. It is now in the documented-commands guard's scope.
- **The LaTeX manuscript carried no ORCID, while the checklist said it did.**
  `docs/JOURNAL_SUBMISSION.md` lists the ORCID as present in the manuscript
  and in `CITATION.cff`; the LaTeX author block had neither, so a reviewer
  reading the paper and the citation metadata saw two different author
  records. The author block now carries the same ORCID the repository
  already publishes, and `test/test_author_affiliation_is_consistent.py`
  keeps the three author blocks -- the LaTeX manuscript, the JOSS paper and
  `CITATION.cff` -- from drifting apart: the JOSS manuscript keeps the
  sections JOSS asks for, the three declare the same affiliation (an
  employer named in one has to be named in all three), and the ORCID is the
  same in each. It does not decide whether an institutional affiliation is
  correct -- that is the author's fact, and `docs/JOURNAL_SUBMISSION.md`
  keeps the question open -- it only makes sure the three statements agree
  whichever way it is settled. Red-green verified: a removed JOSS section,
  an employer in one author block and a missing ORCID each fail the test
  that should catch them.
- **Two paper builds at once are not two independent jobs.** The TeX Live
  install in the paper workflow holds the apt lock, so a second build waits
  behind the first for tens of minutes -- one run sat in "Install LaTeX"
  for 45 minutes while the next one finished that step in two seconds,
  because by then the packages were installed. The workflow now runs one
  paper build per branch and the newer commit's build cancels the
  superseded one: the PDF that matters is the one built from the commit
  that is still the head.
- **The public docs site's quickstart could not be run.** The site (GitHub
  Pages — the first thing a reader sees) promised "256 tests, all pass" for
  a suite that collects 934, and its second command was
  `python scripts/real_tcga_validation.py --report`: the script is at the
  repository root, not under `scripts/`, and it has no `--report` flag. Both
  are fixed, and two guards now cover the site: the test-count guard reads
  it as one of the documents that state the suite's size, and the
  documented-commands guard reads its `<pre>` blocks as commands a reader
  runs. Red-green verified: with the old text restored, one guard names
  the missing script and the other names the disagreement between the
  site's 256 and the suite's 934.
- **The public docs site still advertised 256 tests.** The quickstart's
  comment was corrected, but the header's statistic — the number a reader
  sees first, on the site GitHub Pages serves — was still the one from
  several generations ago. It is corrected now, and the test-count guard
  reads the site's markup (the number and the word "tests" are separate
  elements, so a stale number could hide behind its own markup). The link
  guard now covers the site as well: a link there that names a file of
  this repository, whether written relatively or as an absolute
  `github.com/rollroyces/deepcatch/blob/main/...` URL, has to resolve.
  Both are red-green verified — with the old text restored, the guards
  name the site's 256 and the dead link.
- **The citation file still declared a DOI that does not exist.**
  `CITATION.cff` had been cleaned of its placeholder `doi:` key — a DOI that
  resolved to nothing, which rendered in GitHub's "Cite this repository" box
  as a link to a record that does not exist — but the same defect survived
  in another form: an `identifiers:` entry with `type: doi` whose value was
  the string "to be assigned (preprint)". Not a DOI, but sitting in a
  machine-readable field that says `type: doi`, so a citation tool would
  record it as one. The entry is removed, the reason is written into the
  file, and `test/test_citation_metadata_is_real.py` now requires every DOI
  the file declares to look like one: a `doi:` key or an `identifiers:`
  entry of type `doi` must carry a DOI-shaped value, while declaring no
  identifier at all is the honest state until a release carries one.
- **Seven of the model card's links led nowhere.** `MODEL.md` names the
  documents it points at — `docs/PRETRAINING.md`, `docs/PRETRAIN_BUG.md`,
  `docs/CROSS_STUDY_BENCHMARK.md`, `docs/FUSION_ISOTONIC.md`,
  `docs/SENS_AT_SPEC_ABLATION.md`, `docs/SPARSE_AWARE_ABLATION.md`,
  `docs/HARMONIZATION.md` — but the link *targets* were the bare file names,
  which resolve to the repository root, where those files do not live. Every
  one of them was written correctly in the link text and incorrectly in the
  target, which is why the source read fine. The model card is what a
  reviewer opens before the paper, and in a repository whose argument is
  that the evidence can be checked, its links are the argument. `wiki/Home.md`
  had the same defect in the other direction: a licence badge pointing at
  `LICENSE`, which is the repository's file, not a wiki page. Both are fixed,
  and `test/test_internal_doc_links_resolve.py` now checks every relative
  link in every tracked markdown file — resolving a wiki page the way the
  wiki resolves it, and skipping a link that is an example inside a code
  block.
- **The Docker instruction replaced the published artifacts.** The README
  mounted the repository's own `results/` directory into the container,
  whose entry point runs the validation suite and writes its output there:
  following the instruction overwrote the artifacts the repository cites.
  The mount now points at a directory of the reader's own, and a guard
  checks that no documented command bind-mounts a tracked directory of this
  repository. Adding the two guards moved the suite from 926 to 930 tests,
  so the documents that quote the count were updated together and the
  submission PDF and the arXiv bundle were rebuilt.

- **The arXiv source archive was stale, and nothing would have said so.**
  `dist/deepcatch_v2_3_0.tar.gz` is a build product, not a tracked file, so
  nothing kept it honest: the copy on this machine still carried a manuscript
  that said "883 tests", long after the committed manuscript said 919 and then
  924. A stale bundle is not obviously stale — it compiles, it produces the
  same eight pages and it embeds every figure — so it is only visible by
  comparing it with the manuscript it was built from. The paper workflow now
  uploads the bundle as the `deepcatch-arxiv-source` artifact, which means
  the archive that gets submitted is the one built from the committed
  manuscript rather than a local copy that can be stale, and
  `test/test_arxiv_bundle_is_current.py` fails when a local bundle disagrees
  with a fresh build of the same kind. Adding it moved the suite from 924 to
  926 tests, so the documents that quote the count were updated together and
  both the submission PDF and the bundle were rebuilt.

- **The user guide documented a command line that does not exist.**
  `docs/USER_GUIDE.md` told the reader to run `evaluate.py`,
  `agent2-multimodal-fusion/evaluate.py`,
  `agent3-longitudinal/run_final_fixed.py` and
  `validation_framework_tests.py`, with flags, sample outputs and clinical
  verdicts that no file in the repository backed; the README linked to
  `docs/TUMOR_NAIVE.md`, which was never written; the wiki's install step
  ran `pip install -r requirements.txt`, a file the repository does not
  have. A reader who followed the guide got a usage error on the first
  command, and nothing in CI could see it: the unit tests call functions,
  not command lines, and the documentation guards check numbers, not
  paths. Section 2 now describes the three interfaces the release
  actually ships — the REST API (`uvicorn api.main:app`, `POST /predict`),
  the cohort adapter (`scripts/adapter_local_cohort.py`) and the modules
  under `src/` — and the sections around it were corrected against the
  code: where a run writes its output, how the ctDNA levels and the
  sequencing depth are configured, which error the cohort adapter
  raises, and how the test suite is run.
- **`RUN_ALL.sh` could replace the published artifacts and report success
  after a failure.** The runner wrote straight into `results/`, so a
  validation run overwrote the repository's published artifacts; its
  dependency step installed a `requirements.txt` that does not exist; and
  `set -e` swallowed the exit code of the cross-validation, so a failed
  run still ended in "VALIDATION PIPELINE COMPLETE". It also passed
  `--quick` to `real_tcga_validation.py`, which has no such flag, so the
  TCGA step could never run, and it ran the Jiang step without the
  `--input` that script requires, so that step could never run either.
  The runner now writes to a scratch directory (`RUN_ALL_OUTPUT` overrides
  it) and the repository's published artifacts are untouched; it installs
  the package (`pip install -e .`); the cross-validation's exit code is
  captured; the TCGA step receives a path it can write to and a reduced
  scale it understands; the Jiang step runs against a copy the user
  supplies (`JIANG_INPUT`) and says so plainly when there is none; and an
  unknown flag is an error instead of a silent full run.
- **A documented command can no longer name a file that is not there.**
  `test/test_documented_commands_exist.py` checks every script path in the
  user-facing documents against the tree: a path a reader is told to run
  has to resolve exactly, and a name quoted in a comment only has to
  exist somewhere. It is red-green verified — against the previous guide
  it reports the twenty-four fictional paths it found there, against the
  corrected guide it passes. Adding it moved the suite from 919 to 924
  tests, so every
  document that quotes the count was updated together (919 → 924, the CI
  split 733 + 186 → 738 + 186, the fast-path gate 734 → 739) and the
  submission PDF was rebuilt.

- **Three documents cited artifacts that do not hold their figures.**
  A row-level audit — every table row that names an artifact, checked
  against that artifact's own numbers — found three groups of rows
  naming the wrong source. (1) The Fragmentomics row of the
  cross-platform table quoted the HCC_J OvR AUC as **0.7532**; that is
  the value `results/per_cancer_sens_at_spec.json` held on 2026-09-23.
  The artifact was regenerated on 2026-11-16 and the row's own source
  now reads **0.7221** (CI 0.653–0.792). Nothing caught it because
  that row was not in the number guard. (2) The n=5 claim table of the
  sens-at-spec ablation document cited
  `results/sens_at_spec_ablation.json` for its paired deltas, but that
  file was overwritten by the n=20 run — its `n_seeds` is 20 and its
  paired Δ is now −0.0005, CI [−0.0015, +0.0005], not the n=5
  +0.0075 the table quotes. The n=5 figures are real, but they are
  derived from the two n=5 per-run arrays, so the table now names
  those arrays and says which run they are. (3) Four documents (the
  bioRxiv checklist, the methods paper's limitation table,
  `paper/README.md` and the portfolio roadmap) cited
  `results/methylation_channel_diagnosis.json` for the methylation
  null's figures — 0.9444, 103,863 CpGs, p=1.0, median difference 0.0
  — but that artifact does not hold them; they are in
  `results/methylation_cohort_channel.json`. The diagnosis artifact
  holds the MECHANISM (the depth figures, the coverage baselines, the
  decoder degeneracy), so each citation now names the artifact that
  actually holds what it is cited for.
  All three are guarded now: the cross-platform row by a pin that
  requires the current figure and **forbids** the superseded 0.688 CI
  bound; the methylation figures by pins that name the documents that
  quote them; and the n=5 derived figures by a recomputation test
  instead of a pin, because they are derived, no artifact stores them,
  and the n=5 per-run artifacts are the historical leaky-CV record —
  the leaky `StratifiedKFold` they were produced with is gone from the
  code, so the driver cannot regenerate them and must not be asked to.
- **The number guard could not see a count written with a thousands
  separator.** Its parser matched only DECIMAL numbers, so a figure the
  documents write as "103,863" was invisible to it: it could drift
  freely while the guard stayed green — the same blind spot as a figure
  re-rendered at a coarser precision, a rendering the guard cannot
  parse. The parser now matches a separated integer and matches it
  **exactly** (a count like that is written exactly, so an off-by-one
  must not pass). A plain unseparated integer stays invisible: the
  convention is three decimals or an exact separated count, nothing
  looser.
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

- **Neither manuscript had an acknowledgements section, and both target
  journals expect one.** PLOS Computational Biology and Bioinformatics
  Advances ask for it, so its absence was a submission blocker that the
  repository's own journal checklist had carried as open for months.
  Both manuscripts now have one — `paper/paper.tex` and
  `paper/METHODS_PAPER.md` — thanking the public data providers in the
  wording the JOSS manuscript already used (`paper.md`). What stays the
  author's decision is whether to add a disclosure of the AI-assisted
  audits: PLOS's Ethical Publishing Practice requires AI contributions
  to be reported with the tool names, how they were used, how their
  output was validated and which parts of the work they touched, and
  bioRxiv's FAQ says such use "should be detailed in manuscripts". Two
  ready-to-paste variants sit in the local, gitignored checklist.
- **`CITATION.cff` carried a DOI that resolves to nothing.** The
  placeholder `10.5281/zenodo.to-be-assigned-on-first-release` was
  rendered by GitHub's "Cite this repository" box as a link to a
  record that does not exist, and its "on-first-release" wording had
  become false — the v2.2.0 release was a release and got no DOI
  either. The `doi:` key is now absent until a real one is minted, and
  `ZENODO_SETUP.md` says to add it rather than substitute it.
- **The Zenodo integration was never enabled, and the repository did not
  say so.** The Zenodo API returns zero records for this repository, so
  the v2.2.0 release (2026-08-11) was never archived and carries no
  DOI. A release published with the integration off gets a GitHub
  release and no DOI, so `ZENODO_SETUP.md` now states the verified
  state, gives the one-line query that proves it, and puts "enable the
  toggle, verify, then publish" in front of the tag step.
- **bioRxiv's affiliation requirement is a registration gate, not a
  formality.** Its submission guide asks authors to "supply an
  affiliation so there is an entity providing oversight we can
  contact", and openRxiv's 2026-10-20 statement repeats it. The
  manuscript and `CITATION.cff` declare "Independent Researcher", which
  is not an entity — so the ORCID/affiliation inconsistency the
  checklist already flagged sits on the critical path of the upload,
  and is now documented with its primary sources.

- **The author's email address was in the public tree.** The April 2026
  v1.0 draft (`paper/deepcatch_final.tex`, deprecated since and marked
  "do not use") carried it in the author block and the abstract's
  contact line. The current tree no longer carries it — the history
  still does, and a public history cannot be taken back. A new guard,
  `test/test_public_tree_has_no_personal_contact.py`, keeps it out: the
  rule is about the *kind* of address (a consumer mailbox provider is
  a personal mailbox) rather than the address itself, because writing
  the address into the test would put it back into the tree the test
  exists to protect. Institutional contacts are left alone — the FLARE
  provider's address is how a collaborator requests the data, and
  `deepcatch+sponsors@` is a role alias.
- **The suite's size was stated four different ways, and none of them
  was right.** The English README said "685 collected" next to a
  breakdown that summed to 659, the Chinese README said "659
  collected" and "599 passed", and the two tree comments said 43 and
  42 test files holding 646 and 622 tests. Measured at HEAD on 2026-11-17 (commit `2f8f68b`): **873
  tests** (687 core + 186 deep-learning), the fast-path gate collects
  **688** and runs **652 passed + 23 skipped + 13 deselected**, `test/`
  holds 45 test files and 645 tests, `src/` discovery still collects
  228. Every document that states a count now states the measured one.
  `docs/PORTFOLIO_ROADMAP.md` also claimed "4 CI workflows green";
  there are three workflow files and ten jobs.

- **The methods manuscript offered two checkpoint files the repository does
  not contain.** Its Availability section listed
  `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` and
  `checkpoints/foundation_pretrained_SYNTHETIC_v0.pt` under
  "**Pretrained checkpoints:**" as though they shipped with the
  repository. There is no `checkpoints/` directory at all — the only model
  file the repository carries is `model/motif_model_checkpoint.pt` — so a
  reviewer following the Availability section found nothing there, which
  is the same defect as the Zenodo placeholder DOI this repository has
  already removed. The checkpoints are regenerated by
  `scripts/pretrain_production_finaledb.py`, and **no reported number
  depends on them**: the encoder behind the §3.1 fusion comparison is
  trained in-run with `pretrained=False`, so §3 reproduces from a clean
  checkout. The Availability section now says exactly that, and
  `test/test_manuscript_availability_claims.py` keeps it that way.
- **RESULTS.md's test-status table was three generations stale.** It said
  "28 tests in `test/`" while the suite collected 645, "6/6 jobs green"
  when there are ten, and the note beneath it described a README badge
  reading "228/228" that no longer existed. The table now states the
  measured numbers — **650** in `test/`, **43** in
  `src/foundation/test_integration.py`, **878** across `pytest src test`,
  **10** jobs across 3 workflows, and the pipeline repository's **126**
  tests measured in an isolated worktree at its HEAD `eb1529e` — and is
  marked as repo state rather than a research result, because the
  research figures in Sections 1-7 remain the 2026-08-24 snapshot the
  banner describes. The citation block pointed at v2.2.0 under the old
  title and now matches `MODEL.md`'s entry at v2.3.0. Two claims that
  had become false were corrected: the deep-learning modules are no
  longer "smoke-tested only" (the methylation channel runs end to end on
  real public fragments and is reported as a null, the GNN loads real
  loci, and the fusion layer is benchmarked on real data), and the
  methylation channel is no longer an "active extension project" whose
  answer was expected in Q1 2027 — it has been measured, and the answer
  is that the reachable library measures depth rather than methylation.
  Sections 11 and 12 were also out of order.
- **Two more count claims were unguarded and are now pinned.** The
  README's "Full `src/` discovery" row and RESULTS.md's test-status row
  stated numbers that no check covered, so either could drift silently.
  Both are now compared against a real collection, and the CI job count
  in RESULTS.md against the workflow files.

- **Twenty-seven scripts hardcoded one machine's home directory.** Every
  script that needs the FinaleDB feature cache, the multi-class labels, the
  pre-trained checkpoints or the FinaleMe installation pointed at
  `/Users/hermes/...`. Those paths do not exist for a reader, and several
  had stopped existing for the author too: the repository had moved, so
  `ROOT = Path("/Users/hermes/deepcatch")` in six CADD scripts, the two
  `validation/cfdna/` scripts and a handful of defaults were pointing at a
  directory that was not there — broken for everyone, not just for a
  reader. All of them now resolve through one module,
  `scripts/_paths.py`, which looks the sibling pipeline up next to this
  repository, then in the home directory, with the environment variables
  overriding both. Every resolved path is byte-for-byte the path that
  worked before, except the repository root, which now points at the
  repository. The seven test files that carried the same paths — the
  skip-guards that decide whether the sibling-dependent tests run or skip
  — were migrated the same way, so a reader who clones the sibling
  repository into their home directory gets those tests to run rather
  than skip.
- **The reproduction driver and the demo scripts had the same problem.**
  `paper/REPRODUCE.sh` defaulted its feature, labels and sister-pipeline
  paths to `/Users/hermes/...`, so a reader who followed the driver's own
  instruction — clone the sibling repository *next to* this one — and ran
  it hit paths that do not exist on their machine. The defaults are now
  derived from the repository's own location. The same goes for the three
  recorded demo scripts under `docs/demo/terminal/` and for
  `scripts/reproduce_all.sh`.
- **The driver never said where the missing feature cache comes from.**
  The extracted FinaleDB features are in *neither* repository: the
  pipeline repository gitignores `data/features/` (~600 MB for 627
  samples) and `data/raw/` (~100 GB), because both are rebuilt from the
  public FinaleDB API. So "set DEEPCATCH_FEATURES_DIR to the extracted
  FinaleDB features" left the reader hunting for them. The preflight now
  says which repository does not have them, how to build them (`python
  run_cross_study.py --parallel 8 --max-mb 500` in the pipeline
  checkout) and where to point the variable; the methods manuscript's
  Availability section says the same, and `paper/REPRODUCE.md` documents
  the sibling-relative defaults and the one command that produces the
  cache.
- **The paths a reader is told to use are now guarded.**
  `test/test_availability_claims.py` requires that a named checkpoint
  either exists or sits in a paragraph that says it is not distributed,
  that **no tracked script** points at a home directory (a container's
  `/home/node` is not one — it is part of the image), and that the
  missing-cache message names the remedy rather than only the environment
  variable.

- **Twenty documents handed the reader another machine's path.** The
  copy-paste commands in the documentation invoked one machine's
  virtualenv by absolute path — `env -u PYTHONPATH
  /Users/hermes/deepcatch/.venv/bin/python` — including
  `docs/COLLABORATOR_QUICKSTART.md` and
  `docs/COLLABORATOR_DATA_INTERFACE.md`, which are the entry points for
  the wet-lab collaborators this project is recruiting. A reader who
  pasted one got "no such file or directory" and no clue what to do
  instead. The commands are now relative to the repository, the sibling
  pipeline is written as `../cfdna-fragmentomics-pipeline`, and the
  local caches as `~/.cache/…`. Dated records are left alone: the
  changelog, the audit reports and the generated run reports name the
  machine they were produced on, which is the point of a record.
  `test/test_availability_claims.py` checks the documents as well as the
  scripts.

- **The README's quick start could not work, because it never said where
  the data comes from.** It told the reader to run the cross-study
  benchmark — the headline number — and to validate the install with
  `./.venv/bin/python`, but its installation section created neither the
  virtualenv nor the companion repository, and never mentioned that the
  feature cache is in *neither* repository. `USAGE.md` had the correct
  sequence; the README, which is the front door, did not. Both READMEs
  now install the companion pipeline first, build the cache with
  `run_cross_study.py`, create the virtualenv, and say which artifacts
  come from where. A fresh clone of the public repository was checked
  against this: every file the instructions name is present in it.

- **The arXiv route had no source bundle, and a bundle made by hand would
  have compiled with three missing figures.** arXiv takes a LaTeX source
  archive, not a PDF, and the manuscript reaches outside the repository
  for three figures — so copying `paper/paper.tex` alone produces an
  archive whose `\includegraphics` paths do not resolve, and the
  missing figures are invisible until a reader complains.
  `scripts/make_arxiv_bundle.py` flattens the figures next to the
  manuscript, writes `dist/deepcatch_v2_3_0.tar.gz`, and with
  `--verify` compiles the archive standalone and compares it with the
  committed PDF — 8 pages, the same. This is also the fallback that
  needs no institutional affiliation, which is the field bioRxiv's
  registration asks for. A guard now checks that every figure the
  manuscript includes is in the repository.

- **The collaborator path did not work, and nothing had ever run it.**
  The quick-start's optional third step — merge with the open-data
  labels and run the cross-study benchmark — failed with "Found array
  with 0 sample(s)", and three separate defects sat behind that:

  * The adapter passed the *labels'* directory as the benchmark's
    features directory. The pipeline repository keeps its labels at the
    top level and its feature cache under `data/features/`, so the
    loader found no feature files at all.
  * The benchmark reads the FSD channel from `<sample>.fsd.json`, while
    the collaborator interface asks for `<sample>.fsd_histogram.npy` — so
    a cohort that followed the documentation exactly was invisible to the
    benchmark, which drops any sample whose FSD it cannot read. The
    adapter now derives the `.json` from the `.npy`, so the converted
    cohort is directly consumable.
  * The run wrote the benchmark document and the per-cancer artifact to
    their **published** paths, so a collaborator's optional run replaced
    the manuscript's numbers with a single-publication result. Both
    outputs now go to the adapter's own directory.

  The step also had no test at all — the test module's docstring
  claimed it was "smoke-tested with `--quick` synthetic data" and
  nothing exercised it. There is now one, and it runs the whole step on
  a synthetic cohort and asserts that the repository's published files
  come out byte-identical.

  The merged run covers the open-data publications only: the benchmark
  selects publications from the FinaleDB registry, so a study with no
  FinaleDB publication id is dropped from it. The adapter says so in
  its report and the quick-start now says it too, instead of implying
  that the local cohort joins the benchmark.

  `results/local_cohort/` is now gitignored: the adapter's output
  symlinks into a collaborator's own data files and its manifest names
  their lab, contact and IRB number.

  A fourth defect surfaced in CI, where the sibling checkout is absent:
  the GC-corrected path of the cross-study benchmark imported the
  sibling's loader at module level, so the whole path died with a
  `ModuleNotFoundError` naming a module the reader has never heard of.
  It now falls back to the repo-local implementation, which is
  line-for-line identical to the sibling's — verified bit-for-bit on
  the same cohort. The per-study harmoniser still comes from the
  sibling, which is a documented prerequisite, and the adapter now says
  so plainly instead of letting it fail with a traceback.

- **The JOSS manuscript quoted two figures that no guard pinned.**
  `paper.md` is the paper a JOSS reviewer reads, and it quotes 0.708 (the
  second held-out publication) and 0.8731 (the lower bound of the first
  one's DeLong interval). Both are pinned for the methods paper and for
  the READMEs — but not for `paper.md`, and the interval's lower bound
  was pinned nowhere at all, so a document could carry a stale interval
  and nothing would notice. Every document that quotes either figure is
  now named, `paper.md` included.

- **The GC-corrected pooled row's standard deviation was unpinned.** The
  row quotes "0.9684 ± 0.0027" and only the AUC was pinned, so the std
  could drift. The two runs' stds render differently (0.0027 vs
  0.0013), so the row check now requires the GC pair and forbids the
  legacy pair, and a row that mixed the two runs fails. Verified by
  breaking it: a row with the GC AUC and the legacy std turns the
  guard red.

- **Twenty-one documents quoted a pinned figure without being named for
  it.** The quotable-number list names the documents that must quote each
  figure, and the guard checks those documents — so a document that
  quotes a figure but is not named is a copy that can drift silently.
  An audit that scanned every tracked document with the guard's own
  matcher found 21 of them, among them the preprint's own abstract
  (0.7385), the JOSS manuscript, the landing page and the reproduce
  guide. All are named now: 26 numeric entries gained their missing
  documents and 32 label-addressed entries gained rows, so every copy of
  a benchmark figure is checked.

  One collision needed the mechanism that already exists for it:
  `docs/ALPHAMISSENSE_REAL_LLR.md` quotes 0.9785 as the CADD Top-K=200
  AUC while also quoting the with-motifs pooled AUC, so it is allowed to
  carry that scoped value, exactly as `results/README.md` is allowed to
  carry 0.9915.

- **The installation instructions had three small mismatches.** The
  Dockerfile's version label still read `1.0.0` while the repository is
  at 2.3.0; the Dockerfile's header told the reader to build an image
  called `deepcatch-validation:latest` while the README said
  `deepcatch:latest`; and the README presented the image as "Docker:"
  without saying that it runs the bioinformatics validation suite rather
  than the pipeline. It also never said how to import the library: the
  installed distribution is `deepcatch`, its modules are importable as
  `src.fragmentomics`, and the CLI entry points are installed as
  commands. Checked by building the wheel and installing it into a
  clean virtual environment: the wheel carries `src/` and `validation`,
  `import src.fragmentomics` works, and all three console scripts are
  present.

- **The API advertised a version the repository had left behind.** The
  OpenAPI metadata and the service's root response both carried a
  hardcoded `2.1.0` while the repository was at 2.3.0, and the clinical
  report carried the same stale string. A reader who runs the API —
  which the demo recording does — sees that version in `/openapi.json`
  and `/`, so it contradicts the code serving the request. The version is
  now read from the installed distribution's metadata, or from
  `pyproject.toml` in a source checkout (`src/fragmentomics/version.py`),
  and a guard fails the build if it is hardcoded again.

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
