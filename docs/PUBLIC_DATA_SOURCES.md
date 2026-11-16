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
| Caveat | the *healthy* arm of the training cohort is a modelled TF=0 arm, not real donor plasma |

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

## 4. GEO (surveyed, not usable for this pipeline)

GEO hosts cfDNA studies, but the deposited artifacts are usually
**bigWig coverage tracks or 5hmC/ATAC peaks**, not per-sample fragment
files:

| Accession | Study | Deposited artifact | Usable? |
|---|---|---|---|
| GSE231296 | plasma 5hmC in EGFR-mutant NSCLC | per-sample bigWig + narrowPeak | coverage-only, no fragments |
| GSE282357 | cfDNA chromatin in breast cancer (abemaciclib) | per-sample bigWig (hg19) | coverage-only |
| GSE290667 | cfDNA 5hmC for peritoneal metastasis | bigWig + TSV | coverage-only |
| GSE246110 | cfDNA 5hmC in gastric cancer | processed TSVs (gene-level rpms) | gene-level, not fragments |

The raw reads are in SRA, but downloading and re-aligning a cohort to
recover fragments is a different project (and the repositories' own
"bigWig" products are already the aligned summary). For a
fragment-level method this is the wrong artifact.

## 5. EGA (surveyed, blocked by a DAC)

| Study | Dataset | Samples | Access |
|---|---|---|---|
| Mouliere et al. 2021, *Genome-wide cell-free DNA termini* | EGAS00001005738 | 295 resectable oesophageal adenocarcinoma | DAC request + ethics review |
| Mouliere et al. 2019 (NSCLC immunotherapy) | EGAD50000002661 | 244 NSCLC, 488 fastq (2.2 TB) | DAC request |

`data/EGA_DBGAP_DATA_ACCESS_TEMPLATE.md` in this repo is the prepared
request. These are exactly the cohorts that would make the assay-shift
finding moot — the same lab, the same protocol, a different disease
population — but they are not IRB-free.

## 6. What this means for the paper

- The external validation is on **public, re-fetchable data with a
  scripted provenance** (`scripts/fetch_external_cohort.py` + the
  manifest), which is the strongest available substitute for a clinical
  cohort.
- The **assay-platform boundary** is now an explicit, measured
  result rather than a caveat: cross-study transfer is strongest when
  the assay is comparable (0.9444 on the matched whole-genome
  publication) and blunted when it is not (0.6938 / 0.7127 on the
  whole-exome one).
- The **stage breakdown** remains the one readiness gap: FinaleDB
  carries no stage labels, and the only stage-annotated public
  fragment resources are either inside the CRAG re-packaging of
  publication 6 (no stage column) or behind an EGA DAC. Jiang 2015's
  own supplementary table would carry the HCC stage, but the article is
  not open access (Europe PMC refuses the supplementary bundle:
  "Article with id PMC4372002 is not open access one") and PNAS blocks
  automated downloads of `suppl_file` (HTTP 403), so the stage
  annotation is only reachable through a publisher session.
