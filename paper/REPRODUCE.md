# DeepCatch — Reproduction Guide

> **One-shot reproduction for the methods-paper numbers (`paper/METHODS_PAPER.md`).**

All numerical claims in §3.1, §3.2, §3.3, §3.4, and the *honest framing* language of
§4.2 are reproducible with the commands below. Every script writes a JSON artifact
that lives in `results/`; the paper text is a digest of those JSONs.

---

## One-bash driver

```bash
# from the repo root (/Users/hermes/deepcatch)
bash paper/REPRODUCE.sh           # full reproduction (default)
bash paper/REPRODUCE.sh --quick   # reduced seed / PCA list (CI sanity check)
bash paper/REPRODUCE.sh --no-test # skip the pytest gate
```

`paper/REPRODUCE.sh` runs, in order:

1. **`scripts/cross_study_finallydb.py`** — pooled OOF AUC (§3.2 per-cohort rows,
   pooled harmonized / no_harmonize rows, true-confound control) + per-cancer OvR
   sens@spec with bootstrap 95% CIs. Writes
   `results/cross_study_finallydb.json` and `docs/CROSS_STUDY_BENCHMARK.md`.
2. **`scripts/cross_study_finallydb.py --out-per-cancer-json
   results/per_cancer_sens_at_spec.json`** — standalone per-cancer sens@spec with
   **DeLong** placement-value 95% CIs + PPV@prev grid for the same OvR scores
   (§3.3). Writes `results/per_cancer_sens_at_spec.json`.
3. **Pytest gate** — the six test files in `test/` and `src/foundation/test_integration.py`
   that exercise the cross-study + per-cancer math, plus the harmony / pretrain /
   projection regression suites. Must remain 136+ green.

---

## Step-by-step manual reproduction (what the script does under the hood)

The driver's full command list (assumes the repo venv at `.venv/` is active):

```bash
cd /Users/hermes/deepcatch
PYTHONPATH=  ./.venv/bin/python scripts/cross_study_finallydb.py \
    --top-cancer-n 5 \
    --out-per-cancer-json results/per_cancer_sens_at_spec.json \
    --out-json results/cross_study_finallydb.json \
    --out-md    docs/CROSS_STUDY_BENCHMARK.md
```

(Equivalent pre-built lines, as documented inline in
`scripts/cross_study_finallydb.py:main`.) Then the test gate:

```bash
PYTHONPATH= ./.venv/bin/python -m pytest \
    test/test_biomedical_review_fixes.py \
    test/test_sparse_aware_projection.py \
    test/test_finaledb_pretrained_loader.py \
    test/test_pretrain_bug_fix.py \
    test/test_harmonization_check.py \
    test/test_per_cancer_sens_at_spec.py \
    src/foundation/test_integration.py \
    --timeout=60 -q
```

---

## What each artifact is the source of truth for

| Paper section | JSON artifact | Headline value |
|---|---|---|
| §3.1 (paired synthetic-plasma) | `results/sens_at_spec_ce_n20.json` | foundation AUC 0.5588 ± 0.0094 (n=20 paired seeds) |
| §3.2 per-cohort rows | `results/cross_study_finallydb.json:per_cohort` | jiang 0.9791 ± 0.0028, cristiano 0.9693 ± 0.0022 |
| §3.2 pooled rows | `results/cross_study_finallydb.json:pooled` | harmonized 0.9738 ± 0.0018, sens@99 = 0.782 |
| §3.2 true-confound control | `results/cross_study_finallydb.json:true_confound_control` | harmonized 0.49 vs no_harmonize 1.00 |
| §3.3 per-cancer (bootstrap CIs) | `results/cross_study_finallydb.json:per_cancer` | top-5 OvR AUC + sens@spec |
| §3.3 per-cancer (DeLong CIs + PPV) | `results/per_cancer_sens_at_spec.json` | top-5 OvR with DeLong 95% CIs + ppv_at_prevalence |
| §3.4 ablations | `results/sens_at_spec_ablation_n20.json`, `results/sparse_aware_ablation.json`, `results/cross_study_finallydb.json:true_confound_control` | see §3.4 of the paper |

> **Note on §3.3 CIs.** The per-cancer sens@spec table in §3.3 quotes the **DeLong**
> placement-value 95% CIs (interpolated from `per_cancer_sens_at_spec.json`).
> The companion artifact `cross_study_finallydb.json:per_cancer.*.per_specificity`
> uses **bootstrap percentile** 95% CIs over the same OvR scores for the
> cross-study-fit row. The point estimates agree; the two CIs are
> slightly different in their width (DeLong is narrower at n_pos ≥ 10, bootstrap
> is wider but assumption-free) — both are reported in their respective
> JSONs for transparency. The DeLong row is the canonical reads.

---

## What this script does NOT do

- It does **not** regenerate the pretrained checkpoints; those are
  reproducible via `scripts/pretrain_production_finaledb.py` separately.
- It does **not** run the foundation-model ablations of §3.4 — those
  follow the per-script drivers documented in `docs/SENS_AT_SPEC_ABLATION.md`
  and `docs/SPARSE_AWARE_ABLATION.md`.
- It does **not** contact the FinaleDB REST API — the cross-study pool runs
  against the pre-extracted feature cache at
  `/Users/hermes/cfdna-fragmentomics-pipeline/data/features/` (FinaleDB API
  has been degraded since 2026-09-21; see `docs/PRETRAINING.md`).
- It does **not** generate plots — `RUN_ALL.sh` covers those.

---

## Acceptance checklist

After `bash paper/REPRODUCE.sh` completes successfully, the following must hold:

- [x] `results/cross_study_finallydb.json` exists, valid JSON, `schema_version == "1.0"`.
- [x] `results/per_cancer_sens_at_spec.json` exists, valid JSON, top-5 cancer rows present.
- [x] `docs/CROSS_STUDY_BENCHMARK.md` refreshed with the run timestamp.
- [x] Pytest exit code `0` on the 7-file test gate above.
- [x] Pooled-harmonized `auc_mean ∈ [0.97, 0.98]`, `auc_std ≤ 0.005` (sanity
  floor — see the JSON directly).
- [x] True-confound `harmonized.auc_mean ∈ [0.48, 0.51]` (anti-batches floor).

If any of those fail, open an issue with the JSON path and the failing test
name. The numbers in §3 of `paper/METHODS_PAPER.md` are *not* safe to cite
until they hold.

---

*Last updated: 2026-09-23. Companion to `paper/METHODS_PAPER.md`.*
