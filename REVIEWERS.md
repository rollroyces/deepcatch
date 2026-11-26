# DeepCatch v2.3.0 — Seeking Expert Review

**Open-source cfDNA benchmark: panel-based MRD detection plus a
tumor-naive fragmentomics channel with held-out external validation.
All data open-access.**
**Repo: github.com/rollroyces/deepcatch | Preprint PDF:
`paper/biorxiv_submission_v2.3.0.pdf` | Venue map:
`docs/JOURNAL_SUBMISSION.md`**

---

## What we built

A panel-based ultra-sensitive cfDNA detection pipeline that aggregates
evidence across a patient's real TCGA tumor mutations — the same architecture
as Signatera, CAPP-Seq, and other MRD assays.  Three scoring methods
benchmarked at five ctDNA fractions (10% → 0.1%) with fixed-specificity
metrics. Context-aware sequencing error model (CpG/homopolymer stranding).
Fragmentomics cross-validated on independent real cfDNA (FLARE/GSE317007).

**Key results** (every number traceable to a committed artifact):

- **Panel LLR** AUC 0.921 at 0.1% ctDNA on 20 real TCGA-LUAD patients.
  Duplex-UMI (1e-4 error) + 50k× depth → sens@95% = 1.000.
- **Tumor-naive fragmentomics** (627 open FinaleDB samples): pooled OOF
  AUC 0.9684 ± 0.0027, and **0.9444 [0.8731–1.0000]** on a held-out
  publication that contributes no training data (0.708 on a second, 0.7385
  pooled over the held-out studies carrying both classes).
- **Transfer, measured**: holding one publication out of *training* costs
  most of the pooled AUC (0.681 in the single-type direction, 0.936 the
  other way) — the boundary that a pooled cross-validation hides.
- **The methylation channel is a null**, reported with its mechanism: its
  AUC equals the library-depth baseline exactly, and at the 103,863 loci
  where both arms have coverage the median cancer-vs-control difference is
  0.0 (each locus tested on the samples that cover it; permutation p = 1.0).
- **Reproducibility is enforced**: the benchmark's PCA is seeded, and a
  regression test asserts that two runs agree exactly.

All code open source. All data from GDC (TCGA-LUAD MAFs), GEO
(GSE317007, GSE185307) and the public FinaleDB API/S3 (no IRB, no data-use
agreement).  1040 tests across the two CI jobs (854 core + 186 deep-learning).

---

## What we're asking for

We are seeking feedback from researchers in computational cfDNA biology,
liquid biopsy assay development, and MRD trial design on the following:

### 1. Simulation design & honesty

- Is the context-aware error model (CpG 10×, homopolymer 5×) realistic for
  modern targeted-panel duplex sequencing?
- Should we add fragment-length-aware scoring (tumor fragments ~134bp vs
  normal ~167bp) to the simulation?  Does the literature support it?
- Are there published per-trinucleotide error rates we should calibrate against?

**Suggested reviewers:** Florent Mouliere (Amsterdam UMC — cfDNA fragmentomics),
Nicholas Turner (Royal Marsden — Signatera MRD trials),
Cristian Tomasetti (City of Hope — DELFI fragmentomics)

### 2. Scoring methods

- The Fisher-method panel score (-log₁₀ Poisson p-value per locus) is the
  CAPP-Seq standard.  At 0.1% ctDNA it gives AUC 0.834 vs LLR's 0.921.
  Is this consistent with published CAPP-Seq performance at similar limits?
- Strand-concordance scoring (biallelic signal validation): are we modeling
  strand asymmetry correctly?  Should the noise model be strand-partitioned
  differently?

**Suggested reviewers:** Ash Alizadeh / Maximilian Diehn (Stanford — CAPP-Seq),
Jillian Phallen (Johns Hopkins — targeted error-correction sequencing)

### 3. Assay sweep → production specification

- We sweep error rates (2e-3 to 1e-5) and depths (5k× to 50k×) and identify
  duplex-UMI consensus (1e-4) + 50k× as the production spec that achieves
  AUC 1.000 at 0.1% ctDNA.  Is this consistent with Signatera/Guardant
  commercial specs?
- What panel size (100-500 loci) is realistic for a per-patient tumor-informed
  design given the median LUAD mutation count (199)?

**Suggested reviewers:** Alexey Aleshin (Natera — Signatera),
Richard Lanman (Guardant Health)

### 4. Fragmentomics validation

- The FLARE dataset (GSE317007, 12 HNSCC samples, ONT sequencing, 256 4-mer
  end-motif frequencies) was cross-validated with our fragmentomics module:
  CG-depletion and AT-enrichment patterns confirmed against the published
  Jiang HCC signature.  Without a larger labeled dataset, how would you
  validate the fragmentomics → cancer detection pipeline further?

**Suggested reviewers:** Loris De Cecco (IRCSS Milan — FLARE authors),
Peiyong Jiang (CUHK — 4-mer end-motif HCC detection)

### 5. External validation and transfer

- The tumor-naive classifier was trained on publications 6 + 8 and tested
  on publications 1, 7 and 9 (never trained on, public API/S3, no IRB):
  0.9444 on the held-out whole-genome study, 0.708 on a second, 0.7385
  pooled over the studies carrying both classes.  The held-out
  **whole-exome** study transfers only partly (breast 0.694, prostate
  0.713) — is the platform boundary
  (whole-genome features vs an exome coverage profile) the right
  explanation, or is it a cohort-composition artefact?
- A case-only held-out study is reported as **sensitivity at a
  training-fixed threshold** rather than an AUC.  Is that the most
  defensible summary for a reviewer, or should such a study be excluded
  from the headline entirely?
- **Four more external cohorts have since been scored** — a second
  external family (ccRCC, 12 samples, AUC 0.9444 [0.8156, 1.0], exact
  permutation p = 0.0043), a colorectal cohort with a **graded** label
  (SAUNA, 9 samples carrying a per-sample tumour fraction), a
  hereditary-cancer **screening** cohort (Cleveland Clinic PHTS, 99
  samples, stage-graded), and the only public cohort carrying
  **fragment end motifs** with cancer/healthy labels (223 samples).
  Three questions come out of them, and we would rather a reviewer asked
  them than the other way round:
  - *Does the model track tumour burden?* In SAUNA — colorectal, the
    training cohort's own cancer type, hg38, at the training cohort's
    depth — the score tracks the tumour fraction (Spearman rho 0.8154
    within the seven CRC samples, permutation p = 0.0355; the
    short/long ratio carries more of it, rho 0.9643), and the one
    patient sampled three times has a short/long ratio that rises with
    the fraction (0.207, 0.230, 0.261).  This is the first external
    result here on a graded cancer label rather than a late-stage
    binary contrast.
  - *Does it work in a screening population?* In PHTS — where the
    malignant samples are mostly subclinical — it does not: AUC 0.4506
    [0.3343, 0.5669], and **not even a model fitted and cross-validated
    inside that cohort separates it** (0.5465 [0.4319, 0.6612]).  So
    that null is about the cohort, not about the model, and we report it
    as such.  Note also that the "pre-diagnosis" arm there is drawn
    1-7 days **before** the diagnosis, not months: the longest lead time
    in the whole cohort is seven days, so it is not a lead-time cohort.
  - *Is the external AUC a sequencing-depth artefact?* Depth alone
    separates neither held-out study's classes above chance (0.5342 on
    publication 7, 0.5271 on publication 1), so the headline is not a
    depth effect — although part of the score is depth-correlated, and
    removing it lowers publication 7 to 0.7436.  The ccRCC cohort does
    carry a depth confound (depth alone reaches 0.8333 there), and the
    model's separation survives it (0.9722 with the depth removed), so
    we publish that baseline next to the number.
  - *What about end motifs — the channel the cfDNA field actually
    argues from?* The end-motif cohort's class AUCs look strong
    (0.8207 on the published CC-/CA-end statistic, permutation p =
    5e-05; 0.9958 within the cohort), and the confound control is what
    decides how they read: **the same 256 features separate two cancer
    groups from different studies at 0.9996**, so the study structure
    in those libraries is larger than any class effect and the class
    AUCs cannot be read as cancer signals.  Our own motif model lands
    at 0.3413 there — inverted — which is what a mismatch of motif
    definitions looks like.
  - Is that the right way to report a cohort whose signals are
    confounded with study: as a null about the cohort, or should such
    cohorts be left out of the external evidence altogether?

**Suggested reviewers:** Y.M. Dennis Lo (CUHK — pioneered cfDNA),
Victor Velculescu (Johns Hopkins — CancerSEEK/DELFI)

### 6. The null results

- The methylation channel measures library depth rather than methylation
  on the reachable public libraries (size-selected at ~300 bp, ~1
  observation per CpG).  We report that with a depth-ablation control on
  a genuinely MspI-digested library.  Is that the right call, or should
  the channel be withheld from the manuscript until an independent deep
  MspI cohort exists?
- The neural foundation encoder loses to a logistic baseline at n=40
  (0.554 vs 0.904) and late fusion wins (0.9825).  We frame that as an
  estimator-class result rather than an architecture result.  Would a
  reviewer accept that framing, or expect the architecture comparison to
  be run at a realistic n?

**Suggested reviewers:** an ML/statistics reviewer for the estimator-class
argument; a methylation reviewer for the null

---

## Reproducibility

All results are generated by a single command:
```bash
python real_tcga_validation.py --n-patients 20 --seeds 5 --cancer-types LUAD
```
Output: `results/real_tcga_validation.json` with all AUCs, sweep, and
per-seed metrics.  Requires only numpy/scipy/scikit-learn — the GDC MAFs
are downloaded automatically on first run (open access).

Full test suite: `pytest src test -q` — 1040 collected on Python 3.11-3.14,
torch 2.13, torch_geometric.

---

## Contact

**Yu Ching Lam** — DeepCatch Project
**Repo:** github.com/rollroyces/deepcatch
**PR for review:** github.com/rollroyces/deepcatch/pull/2

Feedback via GitHub Issues or email.  If you find this work useful, a
GitHub star helps us reach more reviewers.
