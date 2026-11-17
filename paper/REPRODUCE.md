# DeepCatch — Reproduction Guide

> **One-shot reproduction for the methods-paper numbers (`paper/METHODS_PAPER.md`).**

All numerical claims in §3.1, §3.2, §3.3, §3.4, and the *honest framing* language of
§4.2 are reproducible with the commands below. Every script writes a JSON artifact
that lives in `results/`; the paper text is a digest of those JSONs.

---

## One-bash driver

```bash
# from the repo root
bash paper/REPRODUCE.sh            # full reproduction (5 seeds), writes into the repo
bash paper/REPRODUCE.sh --quick    # 3 seeds, writes to a scratch dir, repo untouched
bash paper/REPRODUCE.sh --no-test  # skip the pytest gate
```

The driver does three things, in this order:

1. **Re-runs** every script whose artifact the paper quotes.
2. **Proves** the regenerated numbers match the committed ones, number by number
   (`scripts/reproduce_check.py`). A script that succeeds but writes a different
   number fails the run — that is the failure mode a re-run alone cannot catch.
3. **Gates** on the regression tests.

### The steps, in order

| # | Script | Writes | Paper section it backs |
|---|---|---|---|
| 1 | `scripts/cross_study_finallydb.py --no-gc-correction` | `results/cross_study_finallydb.json`, `docs/CROSS_STUDY_BENCHMARK.md` | §3.2 legacy pooled rows, true-confound control |
| 2 | `scripts/cross_study_finallydb.py --gc-correction` | `results/cross_study_finallydb_gc_corrected.json`, `docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md`, `results/per_cancer_sens_at_spec.json` | §3.2 headline, §3.3 |
| 3 | `scripts/cross_study_finallydb.py --include-motifs` | `results/cross_study_finallydb_with_motifs.json`, `docs/CROSS_STUDY_BENCHMARK_with_motifs.md` | §2.1 the with-motifs benchmark — **skipped without the cached motif features** (see below) |
| 4 | `scripts/cross_study_finallydb_shuffled_control.py` | `results/cross_study_finallydb_shuffled_control.json` | §3.2 shuffled-label null control |
| 5 | `scripts/cross_study_loso.py` | `results/cross_study_loso.json`, `docs/CROSS_STUDY_LOSO.md` | §3.2.1 leave-one-study-out transfer |
| 6 | `scripts/external_validation_finaledb.py` | `results/external_validation_finaledb.json`, `docs/EXTERNAL_VALIDATION.md` | §3.2.2 external validation on held-out public publications |
| 7 | `src/fragmentomics/fusion_ablation.py` | `results/fusion_ablation_isotonic.json` | §3.1 / fusion ablation |
| 8 | `scripts/fusion_regularized_smoke.py` | `results/fusion_regularized_smoke.json` | §2.2 the regularized late fusion / Bayesian model average that replaced the neural encoder |
| 9 | `src/fragmentomics/decision_curve_cli.py` | `results/decision_curve.json` | `RESULTS.md` §4 decision-curve operating points |
| 10 | `scripts/methylation_cohort_channel.py` | `results/methylation_cohort_channel.json` | §3.2.4 the methylation channel (reported as a **null**) |
| 11 | `real_tcga_validation.py` | `results/real_tcga_validation.json` | the panel-LLR AUC (0.921) — the abstract's first figure — **skipped without the GDC masked-MAF cache** |
| 12 | `scripts/plot_figures.py` | `docs/figures/fig1..fig4*.png` | Figures 1–4 |
| 13 | `scripts/ultraearly_readiness.py` | `results/ultraearly_readiness.json`, `docs/ULTRA_EARLY_READINESS_AUTO.md` | the readiness verdicts |
| 14 | `scripts/foundation_real_smoke.py --seeds 20 --loss ce` | `results/sens_at_spec_ce_n20.json` | the paired synthetic-plasma honest baseline behind §3.1 — **skipped without the GDC masked-MAF cache** |

Between the steps and the gate sit the **acceptance floors** (see below) and
the **reproduction check**.

> **Four steps need inputs that are not in the repository.**
> Step 3 (the with-motifs sweep) needs the cached 4-mer end-motif features
> (`<sample>.motifs.npy`), which live in the author's local feature
> directory and are tracked in neither repository: without them the motif
> block is zero-filled for every sample and the run produces a *different*
> benchmark. The step is skipped and says so, and the reproduction check
> then reports that artifact as NOT regenerated and ends with "the
> reproduction is incomplete" — the correct verdict, not a silent pass.
> Step 10 (the methylation channel) needs java, the FinaleMe 0.58 jar and
> the cached per-sample CpG matrices, and is skipped the same way. Step 11
> (the panel-LLR validation) needs the GDC masked-MAF cache in
> `validation/tcga/tcga_cache/`, and step 14 (the paired synthetic-plasma
> baseline behind §3.1) needs the same cache — it is also skipped
> without it, and the check then reports
> `results/sens_at_spec_ce_n20.json` as not regenerated.

> **Why two cross-study runs?** `--gc-correction` defaults to ON, so a driver
> that passes no flag writes the GC-corrected run into the legacy path. That
> silently replaced the legacy baseline and broke the comparison the paper
> quotes. The driver therefore runs the sweep twice, with the flags stated
> explicitly, and the two configurations never overwrite each other.

> **The methylation step can be skipped.** It needs `java`, the FinaleMe 0.58
> jar and the cached per-sample CpG matrices. When any of those is missing the
> step prints `SKIP` and the reproduction check then reports the methylation
> channel as *not regenerated* — which is the correct verdict, not a pass.

### `--quick` never touches the repo

`--quick` writes to a temporary directory and reports the result from there.
The published artifacts are not overwritten with a reduced-seed approximation
of themselves. The full run is the one that reproduces the committed numbers
digit for digit.

---

## Data the driver needs

| What | Default | Override |
|---|---|---|
| Extracted FinaleDB features | `/Users/hermes/cfdna-fragmentomics-pipeline/data/features` | `DEEPCATCH_FEATURES_DIR` |
| Multi-class labels | `/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv` | `DEEPCATCH_LABELS` |
| Held-out external cohort | `data/external_cohort` | `DEEPCATCH_EXTERNAL_DIR` |
| Methylation CpG matrices | `../methylation_cohort` | `DEEPCATCH_METHYLATION_DIR` |

The external cohort is fetched from the public FinaleDB API + S3 by
`scripts/fetch_external_cohort.py` — no IRB, no institutional account.

---

## The screening-grade specificity grid

`cross_study_finallydb.py --include-screening` extends the per-cancer
sens@spec grid with the screening-grade operating points {0.995, 0.999}.
The driver leaves that flag **off**, because the committed
`results/per_cancer_sens_at_spec.json` was produced on the legacy grid
{0.95, 0.98, 0.99}; turning the flag on adds fields, and the regenerated
artifact would then no longer match the committed one.

The screening-grade numbers the readiness doc quotes come from the other
route — `scripts/per_cancer_sens_at_spec.py --specificities 0.95 0.98 0.99
0.995 0.999` (see `ULTRA_EARLY_READINESS.md` §1.2) — and the pooled
screening operating points are already carried by
`results/cross_study_loso.json` and `results/external_validation_finaledb.json`.

> This combination used to crash. `SCREENING_SPEC_GRID` was defined with a
> trailing comma, which made it a one-element tuple *containing* the grid;
> the per-cancer standalone writer then passed that tuple where a float
> specificity was expected and died with `unsupported operand type(s) for
> -: 'float' and 'tuple'`. Fixed in `scripts/cross_study_finallydb.py`;
> pinned by `test/test_screening_spec_grid.py`.

---

## Acceptance floors and the reproduction check

After the steps, the driver asserts three floors on the regenerated numbers:

| Floor | Range | Meaning |
|---|---|---|
| legacy pooled-harmonized AUC | `[0.965, 0.985]` | the §3.2 legacy row is in the right place |
| legacy true-confound AUC | `[0.480, 0.510]` | the study-vs-cancer confound still collapses |
| GC-corrected pooled-harmonized AUC | `[0.960, 0.980]` | the headline row is in the right place |

Then the reproduction check compares the regenerated artifacts against the
committed ones:

```bash
python scripts/reproduce_check.py --baseline <committed>/results --current results
```

It compares **numbers**, not files: the artifacts carry a `generated_at`
timestamp, so byte comparison would fail for a perfectly reproducible pipeline.
The list of numbers compared is the authoritative quotable-number list in
`scripts/benchmark_numbers.py` — the same list the manuscript drift guard
(`test/test_paper_numbers_match_artifacts.py`) enforces against the prose. One
list, so the check and the guard cannot disagree about what is quotable.

A drifted number, a missing artifact, or an empty baseline all fail the run.

**Two tolerance bands.** The check separates *smooth* numbers (AUCs) from
*threshold statistics* (sens@spec). A sensitivity-at-specificity is a step
function of the score distribution: the operating point is fixed by the
healthy tail, so one healthy sample moves the threshold and one positive
moves the sensitivity by 1/n_pos. Between a 3-seed and a 5-seed run the AUCs
move by ~3e-4 while a sens@spec point moves by up to ~3e-2. They are judged
against their own band (`--threshold-tolerance`) so the wide band cannot
become a loophole for the AUCs.

| Mode | Seeds | Smooth tolerance | Threshold tolerance |
|---|---|---:|---:|
| full | 42 13 7 99 1234 (the committed list) | 1e-6 | 1e-6 |
| `--quick` | 42 7 1234 | 5e-3 | 5e-2 |

The full run is the one that proves exact reproducibility: same seed list,
same numbers, to the last digit.

---

## What each artifact is the source of truth for

| Paper section | JSON artifact | Headline value |
|---|---|---|
| §3.1 (paired synthetic-plasma) | `results/sens_at_spec_ce_n20.json` | foundation AUC 0.5635 ± 0.0113 (n=20 paired seeds) |
| §3.1.1 (estimator class) | `results/fusion_regularized_smoke.json` | late fusion 0.9825 ± 0.0094, BMA 0.9650 ± 0.0089 |
| §3.2 per-cohort rows | `results/cross_study_finallydb_gc_corrected.json:per_cohort` (headline) / `results/cross_study_finallydb.json:per_cohort` (legacy) | GC-corrected: jiang 0.9661 ± 0.0031, cristiano 0.9698 ± 0.0016; legacy: jiang 0.9789 ± 0.0027, cristiano 0.9704 ± 0.0023 |
| §3.2 pooled rows | `results/cross_study_finallydb_gc_corrected.json:pooled` (headline) / `results/cross_study_finallydb.json:pooled` (legacy) | GC-corrected harmonized 0.9684 ± 0.0027, sens@99 = 0.700; legacy 0.9770 ± 0.0013, sens@99 = 0.785 |
| §3.2 true-confound control | `results/cross_study_finallydb.json:true_confound_control` | harmonized 0.49 vs no_harmonize 1.00 |
| §3.2.1 LOSO transfer | `results/cross_study_loso.json` | jiang → cristiano 0.6805, cristiano → jiang 0.9361 |
| §3.2.2 external validation | `results/external_validation_finaledb.json` | pub 7 (Sun 2019) AUC 0.9444 [0.8731–1.0000]; pub 1 (Snyder 2016) 0.708 [0.5148, 0.9012]; pooled 0.7385 [0.6281, 0.8488]; case-only pub 9 reported as sensitivity 0.6569 at the training-derived threshold 0.9766 |
| §3.2.4 methylation channel | `results/methylation_cohort_channel.json` | **null**: decoded β AUC 0.9444 = coverage baseline 0.9444 || §3.3 per-cancer (bootstrap CIs) | `results/cross_study_finallydb_gc_corrected.json:per_cancer` | top-5 OvR AUC + sens@spec |
| §3.3 per-cancer (DeLong CIs + PPV) | `results/per_cancer_sens_at_spec.json` | top-5 OvR with DeLong 95% CIs + ppv_at_prevalence |
| §3.4 ablations | `results/sens_at_spec_ablation_n20.json`, `results/sparse_aware_ablation.json`, `results/cross_study_finallydb.json:true_confound_control` | see §3.4 of the paper |

> **Note on §3.3 CIs.** The per-cancer sens@spec table in §3.3 quotes the **DeLong**
> placement-value 95% CIs (interpolated from `per_cancer_sens_at_spec.json`).
> The companion artifact `cross_study_finallydb.json:per_cancer.*.per_specificity`
> uses **bootstrap percentile** 95% CIs over the same OvR scores for the
> cross-study-fit row. The point estimates agree; the two CIs are
> slightly different in their width (DeLong is narrower at n_pos ≥ 10, bootstrap
> is wider but assumption-free) — both are reported in their respective
> JSONs for transparency. The DeLong row is the canonical read.

---

## What this script does NOT do

- It does **not** regenerate the pretrained checkpoints; those are
  reproducible via `scripts/pretrain_production_finaledb.py` separately.
- It does **not** run the §3.4 foundation ablations — those follow the
  per-script drivers documented in `docs/SENS_AT_SPEC_ABLATION.md`
  and `docs/SPARSE_AWARE_ABLATION.md`.
- It does **not** contact the FinaleDB REST API — the cross-study pool runs
  against the pre-extracted feature cache (FinaleDB API has been degraded
  since 2026-09-21; see `docs/PRETRAINING.md`). The held-out external cohort
  is the exception: fetch it once with `scripts/fetch_external_cohort.py`.
- It does **not** rebuild the CpG matrices from raw fragments. That is the
  FinaleMe BAM+matrix stage (≈12 min per sample, whole genome); see
  `docs/METHYLATION_PIPELINE.md` §4. The driver reuses the cached matrices.

---

## Acceptance checklist

After `bash paper/REPRODUCE.sh` completes successfully, the following must hold:

- [x] every artifact in the driver's `ARTIFACTS` list exists and is valid JSON
- [x] `docs/CROSS_STUDY_BENCHMARK.md` and `docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md` refreshed
- [x] Pytest exit code `0` on the 7-file test gate
- [x] the three acceptance floors above hold
- [x] the reproduction check reports every quoted number reproduced

If any of those fail, open an issue with the JSON path and the failing test
name. The numbers in §3 of `paper/METHODS_PAPER.md` are *not* safe to cite
until they hold.

---

*Companion to `paper/METHODS_PAPER.md`.*
