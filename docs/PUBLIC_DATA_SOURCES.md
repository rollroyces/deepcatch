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
| [10.5281/zenodo.15097509](https://doi.org/10.5281/zenodo.15097509) | the SAUNA paper's cfDNA fragments as per-sample `*.hg38.frag.tsv.bgz` files (10 files, 13.1 GB) | the same artifact type, so the same extractor would read it; not extracted here (the labels are not in the record) |
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

**One change to the repository's code came out of it.** The deposited
fragment files are standard BED files and carry a header line
(`#chrom start end name mapq strand`), which the fragment reader
could not parse; `read_fragments` now skips a leading `#` line, which
the FinaleDB files do not have. Without that, the cohort is not
readable by the repository's own extractor.

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
- The **stage breakdown** remains the one readiness gap, and it is
  sharper than it looked: FinaleDB carries no stage labels, and the
  only stage-annotated public fragment resources are either inside the
  CRAG re-packaging of publication 6 (no stage column) or behind an EGA
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
