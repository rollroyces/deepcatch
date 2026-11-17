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
agreement).  882 tests across the two CI jobs (696 core + 186 deep-learning).

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

Full test suite: `pytest src test -q` — 882 collected on Python 3.11-3.14,
torch 2.13, torch_geometric.

---

## Contact

**Yu Ching Lam** — DeepCatch Project
**Repo:** github.com/rollroyces/deepcatch
**PR for review:** github.com/rollroyces/deepcatch/pull/2

Feedback via GitHub Issues or email.  If you find this work useful, a
GitHub star helps us reach more reviewers.
