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
| **Held-out external validation** (publications never trained on) | **0.9444** [0.8731, 1.0] on publication 7 (WGS, n=31); **0.708** [0.5148, 0.9012] on publication 1 (WGS, n=52); **0.7385** pooled over the publications carrying both classes; publication 9 is case-only and reported as sensitivity 0.6569 at the training threshold | `external_validation_finaledb.json` |
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
| `geo_gse124600_cohort.json` | The GEO survey's GSE124600 record, **measured** from the deposited files: 142 cancer + 132 control = 274 plasma samples with `crc grade: Stage_I` / `Stage_II` / `Stage_III` / `Stage_IV`, and the methylation matrices' coverage |
| `per_cancer_sens_at_spec.json` | Per-cancer-type AUC, sens@spec and PPV |
| `fusion_ablation_isotonic.json` | Cross-repo fusion ablation (10 seeds) |
| `decision_curve.json` | Decision-curve operating points |
| `methylation_cohort_channel.json` | Methylation channel measured on real public fragments (a null) |
| `methylation_channel_diagnosis.json` | The depth/coverage diagnostics the null rests on |
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
