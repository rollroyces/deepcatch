# DeepCatch Results — v2.3.0

This directory holds the artifacts every reported number is read from.
Each figure below names the file it comes from; `paper/REPRODUCE.sh`
regenerates them and compares the regenerated numbers against the
committed ones.

---

## Headline results

| Result | Figure | Artifact |
|---|---|---|
| Tumor-naive fragmentomics, cross-study (real public FinaleDB plasma, 627 samples) | **0.9684 ± 0.0027** | `cross_study_finallydb_gc_corrected.json` |
| **Held-out external validation** (publications never trained on) | **0.9444** [0.8731, 1.0] on publication 7 (WGS, n=31); **0.708** [0.5148, 0.9012] on publication 1 (WGS, n=52); **0.7385** pooled over the publications carrying both classes; publication 9 is case-only and reported as sensitivity 0.6569 at the training threshold. Depth alone separates neither publication above chance (0.5342 / 0.5271), so the transfer is not a library-depth artefact | `external_validation_finaledb.json` |
| Leave-one-study-out transfer | 0.6805 (Jiang→Cristiano) / 0.9361 (Cristiano→Jiang) against in-domain 0.9789 / 0.969 | `cross_study_loso.json` |
| Cross-repo fusion (10 seeds) | naive average **0.9884 ± 0.0010**; LR fusion 0.9884 ± 0.0009; isotonic 0.9849 ± 0.0014 | `fusion_ablation_isotonic.json` |
| Per-cancer-type AUC | LUAD 0.9825 · OV 0.9930 · CRC 0.9754 · BRCA 0.9738 · HCC 0.7221 — the spread is real, and that is the interesting part | `per_cancer_sens_at_spec.json` |
| Decision-curve operating points | tumor-naive AUC 0.9773; naive average 0.9915 | `decision_curve.json` |
| Methylation channel on real public fragments | **a null**: 103,863 loci, decoded-β AUC 0.9444 = its coverage baseline (difference 0.0, p = 1.0) | `methylation_cohort_channel.json` |
| Mutation-informed panel (real TCGA-LUAD mutations, 5,738, + CADD v1.7) | LLR AUC 0.921 at 0.1% ctDNA, Sens@95% 0.770 — **the plasma reads are simulated** | `real_tcga_validation.json` |

Three tumor-naive benchmarks exist, each in its own artifact and not
interchangeable: 0.9684 GC-corrected (the headline), 0.9770 legacy
no-GC, and 0.9761 GC-corrected with the 4-mer motif block.

---

## Honest limitations

- ❌ **No clinical validation on a real patient plasma cohort.** The
  external validation above is on *public* plasma cfDNA from FinaleDB —
  a different kind of evidence. A wet-lab cohort is the missing piece,
  and it is what the collaboration notice in the top-level `README.md`
  asks for.
- ❌ **The mutation-informed channel's plasma reads are simulated**
  (Poisson-sampled at a target ctDNA fraction); only the mutations
  themselves are real. The tumor-naive channel runs on real plasma.
- ❌ **No validated tissue-of-origin result.** The module
  (`src/tissue_deconv/`) exists and is wired into the fusion, but its
  29-tissue reference atlas is synthetic because the real cfSort atlas
  is not in this repository, and no tissue-of-origin artifact is
  committed.
- ❌ **No head-to-head against Guardant360 or Grail** on the same blood
  draws.
- ❌ **No ComBat / limma empirical-Bayes harmonization comparison** (the
  `sva` / `pycombat` dependency is not taken).
- ❌ Production and regulatory framing; single-source funding; no second
  author.

What *is* reproducible: every figure in the table above, from the
artifacts named in it.

---

## How to reproduce

```bash
# The paper's reproduction driver: regenerates every quoted artifact and
# compares the regenerated numbers with the committed ones.
cd /path/to/deepcatch
bash paper/REPRODUCE.sh

# Quick variant (fewer seeds, looser tolerances)
bash paper/REPRODUCE.sh --quick

# The older end-to-end pipeline
bash RUN_ALL.sh
```

`paper/REPRODUCE.sh` needs the extracted FinaleDB feature cache
(`DEEPCATCH_FEATURES_DIR`) and the labels file (`DEEPCATCH_LABELS`); see
`paper/REPRODUCE.md`.

---

## File index

| File | Description |
|---|---|
| `cross_study_finallydb_gc_corrected.json` | **Headline** tumor-naive cross-study benchmark (GC-corrected) |
| `cross_study_finallydb.json` | Legacy tumor-naive benchmark, no GC correction |
| `cross_study_finallydb_with_motifs.json` | The same benchmark with 4-mer end motifs |
| `cross_study_finallydb_shuffled_control.json` | Shuffled-label null control |
| `cross_study_loso.json` | Leave-one-study-out transfer |
| `external_validation_finaledb.json` | Held-out external validation on publications 1 / 7 / 9 |
| `external_cohort_zenodo_ccrcc.json` | A **second external family**: ccRCC plasma from Zenodo 10.5281/zenodo.21269430 (AUC 0.9444, exact permutation p = 0.0043; all six patients are pT3a, the absolute threshold does not travel, and the cohort's depth confound is measured: depth alone separates the arms at AUC 0.8333, the model's separation survives it at 0.9722) |
| `zenodo_cfdna_survey.json` | The Zenodo-side public-data survey: the queries, their hit counts and the candidate records that were read |
| `zenodo_endmotif_prep.json` | The end-motif tables unpacked: 223 samples (87 healthy, 136 cancer across breast, NSCLC and metastatic colorectal) with their per-sample 4-mer frequencies in the repository's own motif order |
| `external_cohort_zenodo_endmotif.json` | The end-motif cohort, scored **with its confound control**: the class AUCs look strong (0.8207 on the published CC-/CA-end statistic, 0.9958 within the cohort) and the same 256 features separate two CANCER groups from different studies at 0.9996, so the study structure is larger than any class effect and none of the class AUCs can be read as a cancer signal |
| `zenodo_sauna_prep.json` | The SAUNA cohort's **preparation**: nine colorectal plasma files at 15-33x, hg38, 148-307 M fragments each, every one showing the nucleosome ladder (these libraries are not size-selected, so in four of the nine the di-nucleosome peak is the taller one) |
| `external_cohort_zenodo_sauna.json` | The cohort with a **graded** cancer label: the tumour-naive model's score tracks the tumour fraction (Spearman rho 0.6611 over the nine, p = 0.0641; **0.8154** within the seven CRC samples, p = **0.0355**), the short/long ratio carries more of it still (rho **0.9643** within the CRC samples, p = **0.00305**), and the one patient sampled three times has a short/long ratio that rises monotonically with the tumour fraction. The binary contrast (7 CRC vs 2 benign, 0.7143) is far too thin to quote alone -- the gradient is the readout |
| `zenodo_phts_prep.json` | The Cleveland Clinic PHTS cohort's **preparation**: the deposit's files are BEDPE and had to be converted to fragment form, the deposit is mapped to **hg19** while the feature grid is hg38 (0.433% of fragments dropped), and every sample shows the nucleosome ladder |
| `zenodo_immunotherapy_prep.json` | The immunotherapy cohort's **preparation**: 184 of the deposit's 185 longitudinal plasma samples (68 patients, six academic centres, pembrolizumab), their 256 4-mer end-motif frequencies in the repository's motif order, the study's own rMDS per 500 kb bin, and the DELFI channel -- with the study's filters (50-350 bp, MAPQ >= 30) applied, because the deposited files carry every length and every MAPQ, and the ENCODE blacklist exclusion not applied, because the track is not in the deposit |
| `external_cohort_zenodo_immunotherapy.json` | The treatment-response cohort, scored **with its confound control**: the study's own metric finds **0 of 5,527** bins differential between immunotherapy responders and non-responders at FDR < 0.1 (the paper reports 1,080), and the classifier is at AUC **0.3373** -- 0.4162 with the study's own batch correction -- while the same features separate two collection centres that **share** a response label at **0.9813**, and the centres' responder fractions run from 0.35 to 0.79. The cohort's structure is the centre; the depth is not what does the separating either (depth alone 0.4530, and with the depth regressed out 0.3356) |
| `external_cohort_zenodo_phts.json` | A **third external family**, and the first stage-graded *fragment* cohort: 99 PHTS screening plasmas from Zenodo 10.5281/zenodo.15109495. The tumour-naive model does **not** separate them (AUC 0.4506, CI95 [0.3343, 0.5669], sampled permutation p = 0.80) and a model fitted inside the cohort does not separate either (AUC 0.5465, CI95 [0.4319, 0.6612]) — so this cohort has no class signal to find and is not a valid external test in either direction. The depositor's published effect shows up in the highest-burden arm only (short/long ratio AUC 0.6455, CI95 [0.5008, 0.7901]) The "pre-diagnosis" arm (n=9) was drawn 1-7 days before the diagnosis, so the cohort's longest lead time is 7 days: it is not a lead-time arm |
| `geo_gse124600_cohort.json` | The GEO survey's GSE124600 record, **measured** from the deposited files: 142 cancer + 132 control = 274 plasma samples with `crc grade: Stage_I` / `Stage_II` / `Stage_III` / `Stage_IV`, and the methylation matrices' coverage |
| `methylation_cohort_geo_gse124600.json` | What the GSE124600 methylation matrix **contains**: a methylated-CpG call matrix, the only class signal being the called-CpG count (AUC 0.7688), and the negative answer on the stage breakdown (the count does not rise with stage, rho 0.0961) |
| `per_cancer_sens_at_spec.json` | Per-cancer-type AUC, sens@spec and PPV |
| `fusion_ablation_isotonic.json` | Cross-repo fusion ablation (10 seeds) |
| `decision_curve.json` | Decision-curve operating points |
| `methylation_cohort_channel.json` | Methylation channel measured on real public fragments (a null) |
| `methylation_channel_diagnosis.json` | The depth/coverage diagnostics the null rests on |
| `cross_platform_readiness_probe.json` | The **measured** cross-platform readiness: which of FinaleMe's components are on this machine and which are not, so the state in `docs/CROSS_PLATFORM_READINESS.md` is a number rather than prose |
| `harmonization_check.json` | Harmonisation and true-confound controls |
| `real_tcga_validation.json` | Mutation-informed panel validation |
| `cadd_*.json` | CADD scoring and per-subgroup panel results |
| `publication_readiness.json` | FinaleDB infrastructure and per-publication readiness |
| `*.png`, `ensemble_waterfall.png` | Figures |
| `REPRODUCE_REPORT.md` | The last run of the **older** `scripts/reproduce_all.sh` benchmark pipeline (2026-09-28). The numbers the manuscript quotes are proven by `paper/REPRODUCE.sh`, which re-runs every quoted script and then compares the regenerated numbers against the committed ones — a report from the older pipeline is not evidence that anything was reproduced |

`*.json` artifacts are gitignored by default and the ones this project
quotes are force-added; the raw validation JSONs are regenerated by the
scripts.

---

*Last updated: 2026-11-24 (v2.3.0). Every figure above is quoted from
the artifact named beside it.*
