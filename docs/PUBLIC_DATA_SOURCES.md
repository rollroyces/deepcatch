# Public data sources for external validation

A survey of where else DeepCatch's features could be validated on public,
IRB-free data — what exists, what it would cost, and why each candidate
is or is not usable today. Written while building the external
validation on FinaleDB publications 1/7/9 (see
[`EXTERNAL_VALIDATION.md`](EXTERNAL_VALIDATION.md)).

The short version: **FinaleDB is the only source that serves the exact
artifact the pipeline needs — per-sample fragment files in a uniform
format — without a data-access committee.** Everything else either
needs a DAC, serves a different measurement, or is a re-packaging of
data already in FinaleDB.

## 1. FinaleDB (used)

| Property | Value |
|---|---|
| Access | Public REST API + public S3 bucket, no registration |
| API | `http://finaledb.research.cchmc.org/api/v1/{publication,seqrun,summary}` |
| Objects | `entries/<seqrun_id>/hg38/<seqrun_id>.hg38.frag.tsv.bgz` |
| Cohort | 2,579 seqruns / 2,505 samples, 8 publications, 1994-2024 |
| Labels | per sample: tissue, disease, publication |
| Publications | 1 Snyder 2016 · 3 Burnham 2018 (urine) · 4 Rabinowitz 2019 (no disease labels) · 5 Viswanathan 2018 · 6 Jiang 2015 · 7 Sun 2019 · 8 Cristiano 2019 · 9 Adalsteinsson 2017 |
| Caveat | the training cohort's control arm is real donor plasma (FinaleDB: `disease: Healthy`, `tissue: blood plasma`), but the per-study z-score is class-mixed and the two training studies have very different case mixes |

**Object availability is partial.** Fragment files (`.frag.tsv.bgz`),
insert-size metrics and the `.bigWig` coverage objects return 200;
`.bedGraph.bgz` objects return 403. The readiness probe used to test a
single key and reported "private" — see
`scripts/_verify_publication_readiness.py`, which now probes several
kinds and reports `partial`.

**Two usable axes of independence** were found:

- *Study*: publications 1, 7 and 9 were never trained on. Publication
  7 is whole-genome and transfers (0.9444); publication 1 is also
  whole-genome and transfers less well (0.708, thirteen tumour types
  against nine controls); publication 9 is whole-exome and transfers
  only partly (breast 0.6938, prostate 0.7127).
- *Assay*: publication 9 is whole-exome, everything else is
  whole-genome, so the assay shift is measurable rather than assumed.

**The held-out cohort is the complete set, not a convenient subset.**
FinaleDB publishes 8 studies. Six of them serve blood-plasma cfDNA;
the other three are excluded for reasons that are properties of the
data, not of our convenience (measured from the public API,
2026-11-16):

| Publication | Study | Plasma cfDNA? | Why it is not in the held-out cohort |
|---|---|---|---|
| 1 | Snyder 2016 (Cell) | yes — 60 samples | **used**: held-out, matched (46 cancer + 12 controls) |
| 3 | Burnham 2018 | **no — urine** (141 seqruns, `tissue: urine`) | a different analyte; urinary cfDNA fragments are much shorter, so this is an assay-transfer study, not a plasma cancer cohort |
| 4 | Rabinowitz 2019 | no disease labels (74 seqruns, prenatal NIPD) | a fetal monogenic-disorder study, not a cancer cohort |
| 5 | Viswanathan 2018 (CRPC) | **no plasma seqruns** (0) | FinaleDB serves no plasma samples for it |
| 6 | Jiang 2015 | yes | **training** cohort |
| 7 | Sun 2019 (Genome Res) | yes — 31 samples | **used**: held-out, matched (18 cancer + 13 controls) |
| 8 | Cristiano 2019 (DELFI) | yes — 540 samples | **training** cohort |
| 9 | Adalsteinsson 2017 | yes — 1,508 samples, **all cancer, no controls** | **used**: held-out, case-only (reported as sensitivity at the training threshold) |

So publications 1, 7 and 9 are *every* plasma cfDNA study FinaleDB serves
that the model never trained on. The cohort cannot be enlarged by
choosing differently: publication 9 has 1,508 samples but no controls
at all, so it can only ever be reported case-only, and publication 1 has
all 58 of its labelled plasma samples in the cohort.

## 2. Zenodo — the CRAG release (surveyed, not used)

| Property | Value |
|---|---|
| Record | [6914806](https://doi.org/10.5281/zenodo.6914806) (concept 3928545, v4) |
| Title | CRAG: de novo characterization of cfDNA fragmentation hotspots in plasma WGS (Zhou et al., *Genome Medicine* 2022) |
| Licence | CC-BY-4.0, open |
| Contents | `liver.tar` 3.3 GB, `breast_1..4.tar` 2.5–3.1 GB (per-sample chr22 fragment files), `BH01.chr22.frag.bed.gz` 90 MB, `Supplemantary_Table 2_Hotspots.zip` 44 MB, `Supplemantary_Materials_1_source_code.zip` 6 MB, plus GC / chromatin-state / mappability references |

**Why it is not a second independent cohort:** its samples are
FinaleDB's. The label table inside the source-code archive
(`CRAG-master/unsupervised_analysis/Sun_data_ID.xlsx`) is exactly
FinaleDB publication 6's composition — 225 samples: 90 HCC, 67 HBV, 36
cirrhosis, 32 healthy — with no stage column. So CRAG re-packages the
same data with a different file layout; running an "independent"
validation on it would be double-counting publication 6.

What it *is* good for: the multi-gigabyte archives are a cross-check on
the feature extractor's behaviour on the same fragments through a
different pipeline, and the hotspot table is a useful comparison for
the fragmentation features.

## 3. Zenodo — FinaleMe (used, for the methylation path)

| Property | Value |
|---|---|
| Record | [14013719](https://doi.org/10.5281/zenodo.14013719) (concept 7647045) |
| Title | Code and de-identified datasets for FinaleMe (Liu et al., *Nat Commun* 2024) |
| Contents | `FinaleMe-main.zip` 43 MB (source), the WGBS prior bigwig 324 MB, the CG-motif bedgraph 155 MB, the CpG index 130 MB, pretrained HMM models, and the reference panels (`meth_wgbs.tar` 1.1 GB, `ref_panel_wgbs.tar` 5.1 GB, `HD_45/46.WGS.b37.tsv.gz` 1.5–1.7 GB) |
| Caveat | the released pretrained HMM models are binary-incompatible with the released source (see `docs/METHYLATION_PIPELINE.md`) |

## 4. GEO (surveyed again, 2026-11-24 -- and the first survey was too broad)

The earlier survey of this section listed four accessions and concluded
that GEO hosts "bigWig coverage tracks or 5hmC/ATAC peaks, not per-sample
fragment files". Re-run over **99 cfDNA-related GEO series**
(`"cell-free DNA"` plus `fragment*` / `fragmentomic*` / `"DNA methylation"`,
`db=gds`, filtered to series), the conclusion is true for fragments and
**false for methylation**: GEO does host per-sample **methylation
matrices**, and for plasma cfDNA cohorts it hosts them **with per-sample
stage labels**.

| Accession | Study | Deposited artifact | Labels | Usable here? |
|---|---|---|---|---|
| **GSE124600** | blood-based colorectal cancer detection (Kang et al.) | per-sample **methylation matrices** (`sumP_..._umepm.txt.gz`): 205,250 CpG loci x **274 plasma** samples (268 + 6 columns across the two deposited plasma files), plus 66 tissue samples | plasma, CRC patient / healthy, and **crc grade: Stage_I / Stage_II / Stage_III / Stage_IV** on the cancer plasma (142 labelled, 132 controls unlabelled) | a *methylation* cohort with stage labels -- not the fragmentomics channel this repository scores |
| GSE97932 | cell-free DNA methylation markers, lung cancer | one combined processed matrix (485k CpGs) | cancer / normal | same caveat |
| GSE97923 | the same study, colorectal | one combined processed matrix | cancer / normal | same caveat |
| GSE118690 | genome-scale methylation, tissue-of-origin of plasma cfDNA | per-sample + combined matrices | tissue-of-origin labels | same caveat |
| GSE328856 | G-quadruplex-forming cfDNA fragment, HCC (2026) | per-sample **bigWig** + raw reads | serum, HCC / healthy | coverage-only, and **serum**, not plasma |
| GSE231296 | plasma 5hmC in EGFR-mutant NSCLC | per-sample bigWig + narrowPeak | — | coverage-only, no fragments |
| GSE282357 | cfDNA chromatin in breast cancer (abemaciclib) | per-sample bigWig (hg19) | — | coverage-only |
| GSE290667 | cfDNA 5hmC for peritoneal metastasis | bigWig + TSV | — | coverage-only |
| GSE246110 | cfDNA 5hmC in gastric cancer | processed TSVs (gene-level rpms) | — | gene-level, not fragments |

Two findings worth keeping:

- **The stage-labelled public plasma cohort exists.** GSE124600's cancer
  plasma samples carry `crc grade: Stage_I`, `Stage_II`, `Stage_III` and
  `Stage_IV` -- **142 cancer plasma** samples are labelled, and the
  remaining **132 plasma** samples are the controls, which carry no stage
  label, for 274 plasma samples in total. That is the first public cohort
  with the stage labels `ULTRA_EARLY_READINESS.md` §1.3 says the open-data
  cohort lacks. It is a *methylation* cohort, so it cannot be scored with
  this repository's fragmentomics model; it would need a methylation feature
  extractor and its own model, which is a project rather than a week of
  work.
- **The methylation matrices are complete, and verified rather than
  assumed.** Downloading the three deposited files (14,341,039 +
  1,897,647 + 7,333,412 bytes, each a valid gzip) and counting their
  columns gives **205,250 CpG loci** in every matrix, **274 plasma sample
  columns** in total -- 268 in `sumP_..._dup1.5_...` and 6 in
  `sumP_..._nodup1.5_...`, which is exactly the 274 plasma samples GEO
  records -- and 66 tissue columns (33 colorectal cancer tissue + 33
  adjacent mucosa). The per-sample stage labels and the per-sample
  methylation values therefore live in the same series, which is what makes
  this cohort usable for the stage-stratified analysis §1.3 asks for. The
  measurements, the file hashes and the provenance are in
  `results/geo_gse124600_cohort.json`.
  *(Re-checking this needs one caveat: macOS `zcat` cannot read a `.gz`
  file -- it fails with `can't stat`, which looks exactly like a truncated
  download. `gzip -dc` reads the file, and every byte count above matches
  the FTP `Content-Length` exactly. An earlier note in this repository
  called that 1.9 MB file truncated; it was not.)*
- **No fragment-length matrix for a cfDNA cancer cohort exists in
  GEO.** The one recent fragment study deposits coverage tracks (and
  uses serum, not plasma). The raw reads are in SRA, but downloading
  and re-aligning a cohort to recover fragments is a different project
  (and the repositories' own "bigWig" products are already the aligned
  summary). For a fragment-level method this is the wrong artifact --
  which is why the held-out cohort in §1 comes from FinaleDB, whose
  fragment files are served directly.

### 4.1 What the GSE124600 methylation matrix actually contains

The survey above established that GSE124600 carries stage-labelled plasma
samples and per-sample methylation matrices; it did not open the matrices.
`scripts/geo_gse124600_methylation_cohort.py` opens them (205,250 loci x
274 plasma samples, plus a sample metadata table read from the GEO sample
records), and the answer is in
`results/methylation_cohort_geo_gse124600.json`:

- **The matrix is a methylated-CpG CALL matrix with normalised counts,
  not a fractional methylome.** 91.1% of the cells are exactly 0 and the
  maximum is 6,949.25, so a zero means "this CpG was not called
  methylated in this sample", not "this CpG is 0% methylated".
- **The per-locus summaries are anti-correlated with cancer** (AUC 0.27
  to 0.41 over the shared-locus means, the mean of the non-zero values
  and the total signal) -- the signal is spread over more loci, which is
  what a normalised count matrix does.
- **The only class signal is the number of CpGs called methylated per
  sample**: 15,637 (median) in the healthy controls against 20,130 in
  the cancer plasma, **AUC 0.7688** (permutation p < 0.0001). That is a
  burden/coverage measure, and it is what this repository's own
  methylation channel found on FinaleMe, where the decode emits hard 0/1
  state assignments and the "methylation" the channel consumed was the
  fraction of covered loci called methylated. An independent, larger,
  stage-labelled deposit reproduces the null.
- **And the stage breakdown is a negative answer.** The called-CpG count
  does *not* rise with stage (Spearman rho = **0.0961**, p = 0.26 over
  the 142 cancer samples; the per-stage AUCs against the healthy
  controls are 0.7130 / 0.8223 / 0.7188 / 0.8288 for Stage_I / II / III /
  IV, which are noisy at these arm sizes). So the one public cohort with
  per-sample stage labels **does not** deliver a stage-stratified
  early-detection result, and `ULTRA_EARLY_READINESS.md` §1.3's readiness
  gap is closed with a measurement rather than a plan.

## 5. EGA (surveyed, blocked by a DAC)

| Study | Dataset | Samples | Access |
|---|---|---|---|
| Mouliere et al. 2021, *Genome-wide cell-free DNA termini* | EGAS00001005738 | 295 resectable oesophageal adenocarcinoma | DAC request + ethics review |
| Mouliere et al. 2019 (NSCLC immunotherapy) | EGAD50000002661 | 244 NSCLC, 488 fastq (2.2 TB) | DAC request |

`data/EGA_DBGAP_DATA_ACCESS_TEMPLATE.md` in this repo is the prepared
request. These are exactly the cohorts that would make the assay-shift
finding moot — the same lab, the same protocol, a different disease
population — but they are not IRB-free.

## 6. Zenodo — the cfDNA search (2026-11-24)

The GEO survey (§4) answered the fragment/methylation question for GEO.
Zenodo held two known records (§2, §3) and **no search**, so it was
searched: eight queries against the Zenodo records API with
`type=dataset`, 25 records per page. The queries, their hit counts on
the day, and the records that were read are in
`results/zenodo_cfdna_survey.json`, written by
`scripts/zenodo_cfdna_survey.py`. Zenodo's search is an OR over the
query terms, so the hit counts are large (3,748–31,903 dataset hits,
except one query at 34) and **142 distinct records** came back across
the first pages; the candidates below were read by hand.

| Record | What it holds | Usable here? |
|---|---|---|
| **[10.5281/zenodo.21269430](https://doi.org/10.5281/zenodo.21269430)** | 12 plasma samples (6 clear-cell renal cell carcinoma, 6 healthy controls) as **per-sample hg38 fragment BED files** + a sample metadata table with pT stage; 105,849,322 bytes | **yes — extracted and scored, see §6.1** |
| [10.5281/zenodo.15097509](https://doi.org/10.5281/zenodo.15097509) | the SAUNA paper's cfDNA fragments as per-sample `*.hg38.frag.tsv.bgz` files (9 samples, 15-33x coverage, 13.1 GB) plus a `metadata.tsv` | the same artifact type **and a richer label**: per-sample group, coverage and tumour fraction -- see §6.2 |
| [10.5281/zenodo.17253504](https://doi.org/10.5281/zenodo.17253504) | per-sample 2/3/4-mer **end-motif frequency tables** (21 files, 2.3 MB) for healthy, breast, metastatic colorectal and NSCLC cohorts | processed features, but motif frequencies only — not the five-channel profile this repository scores |
| [10.5281/zenodo.23195882](https://doi.org/10.5281/zenodo.23195882) | the LEAF-1 pretraining corpus: 211 GB of cfDNA fragments hosted on Hugging Face; this record holds the manifest and checksums | a manifest, not the data — and its fragments come from FinaleDB and the classic plasma studies this repository already has |
| [10.5281/zenodo.23196827](https://doi.org/10.5281/zenodo.23196827) | chromatin landscapes and cfDNA fragmentation in colorectal cancer development (6.7 GB) | the fragments are bundled with chromatin data, not released per sample |
| [10.5281/zenodo.15738112](https://doi.org/10.5281/zenodo.15738112) | head-and-neck immunotherapy end-motif study: a 35 GB `Frag_Files.tar.gz` | per-sample fragments, but 35 GB for the archive and no per-sample label table in the record |
| [10.5281/zenodo.15109495](https://doi.org/10.5281/zenodo.15109495) | cfDNA fragmentomics in PTEN hamartoma tumour syndrome (7.6 GB) | the fragments are inside a BAM archive, not per-sample BED files |

So the Zenodo side of the goal — processed cfDNA features that are public
and re-fetchable without an IRB — is **real and larger than GEO's**, and
one of the records is directly usable by this repository's own
extractor.

### 6.1 A second external family: ccRCC plasma from Zenodo

The held-out cohort's three FinaleDB publications are one cohort shape.
This is another: a **different cancer type** (clear-cell renal cell
carcinoma), a **different laboratory and protocol** (an in-house
clinical study), and a **different data source** (a Zenodo deposit of
per-sample fragment BED files rather than a FinaleDB publication).

`scripts/external_cohort_zenodo_ccrcc.py` extracts the repository's five
benchmark channels from the 12 deposited fragment files and scores the
cohort with the model trained on FinaleDB publications 6 + 8, using the
**same protocol as `scripts/external_validation_finaledb.py`** (imported
from it, not re-implemented: per-study z-score harmonisation on the
training side, per-publication on the held-out side, PCA fitted on the
training cohort only, logistic regression fitted on the training cohort
only, the held-out cohort scored once). The result is in
`results/external_cohort_zenodo_ccrcc.json`.

| | |
|---|---|
| Cohort | 6 ccRCC patients vs 6 healthy controls, matched within one study |
| Stage | **all six patients are pT3a (stage III)**; two have distant metastasis (M1 / Mx) |
| **AUC** | **0.9444**, DeLong 95% CI [0.8156, 1.0] |
| Exact permutation p | **0.0043** (all 924 label assignments; the scores never see the held-out labels) |
| Harmonisation sensitivity | 0.9444 (per study / per publication) and 0.9722 (both global variants) — the ranking does not depend on the choice |
| Assay depth | median 1,106,290 fragments per patient and 931,526 per control, against about 30 M in the training cohort |

Four things are true about this number, and only the first is a result:

1. **The ranking transfers.** A model that never saw this cancer type,
   this laboratory or this data source ranks all six patients above
   four of the six controls, at p = 0.0043 by an exact permutation test.
   That is a second external family with the same qualitative behaviour
   as the FinaleDB held-out cohorts (0.9444 on publication 7, 0.708 on
   publication 1).
2. **The absolute threshold does not transfer.** At the training
   cohort's own 99%-specificity threshold (0.9766), all six patients
   are above it — and so are three of the six controls, so the
   specificity on this cohort is **0.5**. The reason is visible in the
   scores: this cohort's controls have a median score of 0.869 against
   0.003 for the training cohort's controls, i.e. the whole cohort sits
   high on the score axis. That is a cohort-level shift (a different
   library preparation), not a cancer signal, and an AUC is a rank
   statistic — which is why the rank is the honest number to quote for
   a held-out cohort and a probability threshold taken from another
   cohort is not.
3. **It is a late-stage cohort.** All six are pT3a; two have distant
   metastasis. Like the repository's other held-out cohorts, it says
   nothing about early detection.
4. **It is small.** n = 6 vs 6. The point estimate is the finding; the
   interval is wide, and the specificity grid is not resolvable at this
   n (six controls quantise specificity to multiples of 1/6, so every
   value above 83.3% is the same threshold).

**What it changes.** The pooled external number in the manuscript
(0.7385 over publications 1 and 7) is a statement about one cohort
shape; this is a second shape with the same behaviour. Whether it is
folded into that pooled number, or reported alongside it as a separate
cohort, is a decision about the manuscript's headline and is the
author's — the artifact and this section exist so the choice can be
made on the numbers.

**One thing to watch if it goes into the text.** This cohort's AUC
renders as `0.9444`, and so do two other figures in this repository: the
held-out FinaleDB publication 7's AUC (`0.9444`, CI [0.8731, 1.0000],
n=31) and the methylation channel's leave-one-cohort-out AUC
(`0.9444`). Three different measurements, one rendering — the interval
and the cohort size are what tell them apart. The repository's figure
guard cannot separate them either: a rendering shared by more than one
pinned figure is skipped, because a document that writes "0.9444"
cannot be resolved to one of them by its digits. So if this cohort
enters the manuscript, the number has to arrive with its cohort
attached ("0.9444 on the ccRCC cohort, n=6 vs 6") and never as a bare
0.9444 standing next to the other one.

**One change to the repository's code came out of it.** The deposited
fragment files are standard BED files and carry a header line
(`#chrom start end name mapq strand`), which the fragment reader
could not parse; `read_fragments` now skips a leading `#` line, which
the FinaleDB files do not have. Without that, the cohort is not
readable by the repository's own extractor.

### 6.2 A public MRD-style cohort with a graded cancer label: SAUNA

The second fragment deposit is the one that comes closest to the
readiness document's question. Zenodo 15097509 ("SAUNA -- Simulated
Annealing for Unique Nucleosome Arrangements") publishes 9
medium-coverage (15-33x) cfDNA files as per-sample hg38 fragment files,
and -- this is the part the record's description does not tell you -- a
`metadata.tsv` with a per-sample cohort table:

| sample | group | patient | coverage | tumour fraction |
|---|---|---|---|---|
| SeCT-20 | CRC | SeCT20 | 20.3x | 0.01735 |
| SeCT-22 | CRC | SeCT22 | 30.2x | 0.007067 |
| SeCT-26_t1 | CRC | SeCT26 | 29.5x | 0.2903 |
| SeCT-26_t2 | CRC | SeCT26 | 28.5x | 0.3271 |
| SeCT-26_t3 | CRC | SeCT26 | 30.9x | 0.4029 |
| SeCT-47 | CRC | SeCT47 | 32.6x | 0.004879 |
| SeCT-58 | Benign | SeCT58 | 15.7x | 0.006191 |
| SeCT-61 | CRC | SeCT61 | 28.4x | 0.0249 |
| SeCT-82 | Benign | SeCT82 | 28.3x | 0.00724 |

(The table is read from the deposit by `scripts/zenodo_cfdna_survey.py`
and recorded in `results/zenodo_cfdna_survey.json`. An earlier note in
this survey said the labels were not in the record; that was wrong --
the labels are in the deposit's own metadata file, and a survey that
reads only the record's description reports a cohort as unlabelled when
it is not.)

Why this cohort is the closest public thing to the repository's
early-detection claim:

- **It is labelled with a graded cancer signal**, not a stage: every
  sample carries its own ctDNA tumour fraction. Five of the seven CRC
  samples sit between 0.5% and 2.5% tumour fraction -- the regime where
  a classifier is supposed to work, and where this repository's
  tumour-naive channel has never been tested on real post-resection
  plasma.
- **One patient is sampled three times**, with the tumour fraction
  rising 0.2903 -> 0.3271 -> 0.4029 -- a longitudinal series in public
  data, which is the shape a monitoring claim needs.
- **The depth matches the training cohort's** (15-33x, against the
  training cohort's about 30 M fragments per sample), so this cohort does
  not carry the depth confound the ccRCC cohort does.

What is not done here: the deposit is 13.1 GB and its nine samples hold
about 300 M fragments each, so extracting the repository's five channels
from them is a matter of hours, and the binary contrast (7 CRC vs 2
benign) is thin. It is recorded as the cohort that can close the
early-detection question in public data, with the labels, the depths and
the tumour fractions to do it with. No number from it is quoted yet.

### 6.3 A screening cohort with stage labels: the Cleveland Clinic PHTS plasma

Zenodo 15109495 ("Cell-free DNA Fragmentomics and Second Malignant Neoplasm
Risk in Patients with PTEN Hamartoma Tumor Syndrome"; related articles
10.1016/j.xcrm.2023.101384 and 10.1016/j.xpro.2024.103357) is the one
deposit in the survey whose label table the survey did not pick up, and
picking it up is worth doing: the file `phtscfx_meta.csv` (15 KB) carries
**99 plasma samples with a hereditary-cancer-screening design**, and the
labels in it are the ones this repository's fragmentomics channels have
never been tested on.

| Field | What it holds |
|---|---|
| `cancer` | No 50 / Yes 49 |
| `group2` | No Cancer 50 / SMN 26 / PMN 23 |
| `group3` | No Cancer 50 / **Pre-SMN 4, Post-SMN 22, Pre-PMN 5, Post-PMN 18** |
| `stage` | NA 60 / **0: 5, 1: 22, 2: 10, 3: 2** (39 labelled) |
| `ca_1st`, `ca_2nd` | the organ of the first and the second primary (Breast 27 / Thyroid 14 / Endometrial 6 / …) |
| `draw_after_smn` | Yes 23 / No 2 / Unknown 1 / NA 73 |
| `date_diff` | days between the plasma draw and the diagnosis; **−7 to −1 in the "pre" arms, +1 to +37 in the "post" arms** (NA in the 50 no-cancer samples) |

**Why this is the cohort the stage question has been waiting for.** The
stage breakdown in §6.1 of the readiness document has no stage-graded
*fragment* cohort to run on: GSE124600 carries stage labels but ships
methylation matrices, and the methylation channel is a coverage signal
rather than a feature (§3.2 of the readiness document). This deposit
ships **per-sample fragment data** (`02-bamtobed.tar.gz`, 7.6 GB) and the
labels to go with it, and it has now been run (`results/external_cohort_zenodo_phts.json`):

- **Stage.** Cancer-vs-no-cancer AUC within each stage against the 50
  no-cancer samples — the fragmentomics channel's first stage-stratified
  readout on real plasma. **Measured: no separation at any stage** (see
  below).
- **Pre-malignant.** Pre-SMN/Pre-PMN (n=9) against the 50 no-cancer
  samples: plasma drawn *before* a second primary appeared. **Measured:
  at chance** — and the `date_diff` column says why that is not an
  early-detection question (see below).

**The label table, derived and in the shape the extractor reads.** The
depositor's `phtscfx_meta.csv` is a clinical table, not a cohort table,
so it is converted into the tab-separated form the cohort extractor
takes (`sample_id` as the key, `label` 1/0 from `cancer`, and every
clinical column kept as published) by
`scripts/prepare_phts_label_table.py`, and written next to the cohort's
features under `data/external_cohort/` — the same place the ccRCC
cohort's features live, which is local rather than versioned (the
`data/` tree is ignored; the results are what the repository commits).
Keeping `date_diff` and `draw_after_smn` is what makes the arm labels
readable: the labels say *which* arm a sample is in, and the dates say
*when the plasma was drawn relative to the diagnosis*.

**Scored, and three things the run settled before it produced a number.**
The cohort is unpacked, converted and scored
(`scripts/prepare_zenodo_phts_cohort.py` → `results/zenodo_phts_prep.json`;
`scripts/external_cohort_zenodo_phts.py` →
`results/external_cohort_zenodo_phts.json`) with the same protocol as the
ccRCC cohort: train on FinaleDB 6+8, per-study harmonisation on the
training side and per-publication on the held-out side, PCA(200) +
LogisticRegression fit on the training cohort only, the held-out cohort
scored once and never fitted on.

**The deposit ships BEDPE, not fragment BED.** One row per read *pair*,
with the two mates' coordinates in separate columns; the extractor reads
columns 1-3, which in a BEDPE are the first mate's coordinates. Fed
as-is it reads a "fragment" that is one read long (~146 bp instead of
~167 bp) and the fragment-size distribution collapses onto the
read-length peak — silently, with every channel wrong by the same
factor. The preparation script converts each member to the fragment form
(the pair's outer coordinates) and refuses to report success unless
every sample shows the nucleosome ladder: the mono-nucleosome peak as
the mode, and a visible di-nucleosome shoulder. All 99 samples pass;
302,801,515 fragments in total.

**The deposit is hg19; the grid is hg38.** The fragments run to
chr1:249,239,706 — past hg38's chr1 end (248,956,422) and inside hg19's
(249,250,621) — which is what the depositor's own data-availability
statement says too. The extractor's bin grid is hg38, and it has to be:
that is the grid the training cohort's features were computed on. The
preparation report records the cost — the grid drops 1,310,997 of
302,801,515 fragments (0.433%), those starting past an hg38 chromosome
end.

**The "pre-diagnosis" arm is a days-around-diagnosis arm, not a
lead-time arm.** The depositor's `date_diff` column puts every "pre"
draw **1 to 7 days before** the diagnosis (n=9, median 2 days) and every
"post" draw 1 to 37 days after it (n=40, median 9). The longest
lead time anywhere in the cohort is **7 days**. So the tumour was already
there when the "pre" plasma was drawn: this is a cohort of samples drawn
*around* a diagnosis, and a separation in that arm would not be an
early-detection result. The artifact carries the intervals per arm so
that nobody has to take the arm names on trust.

**The result: the model does not transfer, and the effect does.** The
tumour-naive model — trained on FinaleDB 6+8, never fitted on this
cohort — does **not** separate cancer from no cancer here: AUC
**0.4506**, CI95 [0.3343, 0.5669], sampled label-permutation p = 0.80
(an exact enumeration is not possible at n=99; the test is sampled and
seeded, and the artifact says so). The null survives the harmonisation
choice: three variants that let the cohort normalise itself give
0.4506 / 0.4563 / 0.4563, and the control that normalises it by the
**training** cohort's statistics instead — so its absolute position in
feature space is preserved — gives 0.5388 [0.4232, 0.6544], still
including 0.5. The stage breakdown is flat: the pre-diagnosis arm
(n=9) is at chance (AUC 0.5, CI [0.3052, 0.6948]), stage 0 (n=5) 0.556,
stages 0-1 (n=27) 0.4119, stages 2-3 (n=12) 0.575, all staged (n=39)
0.4621, and the stage gradient is not monotone (Spearman rho 0.172,
p 0.297). This is the first stage-graded *fragment* cohort this
repository has scored, so the stage question now has an answer rather
than a plan: no stage-stratified separation either.

**But the effect the depositor reported is in the data.** Read one
channel at a time against the same no-cancer arm, the second-neoplasm
arm separates on its own: the short/long ratio is higher there (0.2135
against 0.1941; AUC **0.6455**, CI [0.5008, 0.7901] — an interval that
excludes 0.5), the fragments are shorter (median 171 bp against 174),
and the di-nucleosome proportion is lower (0.0258 against 0.0280, AUC
0.3964) — every one of them in the direction Liu et al. report. The
model's score does not carry that effect into a separation, so the null
is a **transfer failure, not an absence of signal**, and the mechanism
is visible: the model reads this whole cohort as cancer-like. Its
no-cancer controls score at the top of the model's range (median score
0.9997, against the training cohort's own controls' 0.0027), so the
per-cohort harmonisation pushes every sample into the saturated region
and the AUC is a ranking among near-identical values.

That is the result worth publishing: not "the model works on a second
cohort" — it does not — but "the effect reproduces in the channels, the
model does not transfer, and here is the measurement of why".

**The survey missed this label table until 2026-11-25.** Its reader only
split on tabs, and the deposit's table is a CSV — so the richest label
table in the survey was recorded as a deposit without one. The reader now
sniffs the delimiter and the map knows this deposit; the survey artifact
carries the table with all 99 rows and its clinical columns.

## 7. What this means for the paper

- The external validation is on **public, re-fetchable data with a
  scripted provenance** (`scripts/fetch_external_cohort.py` + the
  manifest), which is the strongest available substitute for a clinical
  cohort.
- The **assay-platform boundary** is now an explicit, measured
  result rather than a caveat: cross-study transfer is strongest when
  the assay is comparable (0.9444 on the matched whole-genome
  publication) and blunted when it is not (0.6938 / 0.7127 on the
  whole-exome one).
- The **stage breakdown** has now been measured on two public
  stage-labelled plasma cohorts instead of being planned, and both
  answers are negative: GSE124600's methylation matrices (the called-CpG
  count does not rise with stage, rho 0.0961, p 0.26) and the Cleveland
  Clinic PHTS fragments (§6.3; AUC 0.4506 pooled, no stage
  stratification, and the "pre-diagnosis" arm drawn 1 to 7 days before
  the diagnosis). The gap itself is sharper than it looked, and it is a
  gap in the *data* rather than in the analysis: FinaleDB carries no
  stage labels, and the only stage-annotated public fragment resources
  are either inside the CRAG re-packaging of publication 6 (no stage
  column) or behind an EGA
  DAC -- and the held-out cohort's own publication 7 is one of them: its
  11 CRC patients and 14 liver-transplant recipients were submitted to the
  EGA under `EGAS00001003117` and `EGAS00001003116`. But the held-out cohort's own publications do describe their
  stage in public records that are reachable without a DAC -- GEO's
  extract protocol for publication 1 (GSE71378) says its 52 cancer
  samples are **Stage IV**, and publication 9 is a **metastatic**
  cohort by its paper's title. So the open cohort carries no Stage I
  case that the public records identify
  (`results/external_cohort_stage_annotation.json`), and the external
  AUCs are a late-stage transfer result rather than an early-detection
  one. Jiang 2015's
  own supplementary table would carry the HCC stage, but the article is
  not open access (Europe PMC refuses the supplementary bundle:
  "Article with id PMC4372002 is not open access one") and PNAS blocks
  automated downloads of `suppl_file` (HTTP 403), so the stage
  annotation is only reachable through a publisher session.
- **The model does not transfer to a third external family, and the
  reason is measured.** On the PHTS cohort the tumour-naive model scores
  every sample at the top of its range -- the cohort's no-cancer
  controls' median score is 0.9997 against the training cohort's own
  controls' 0.0027 -- so the AUC is a ranking among near-identical
  values, while the effect the depositor published is present in the
  channels. A null that comes with a measured mechanism is worth more
  than a caveat, and it is the shape a reviewer will ask for.
