# Reproducing every published number in this repo

> **One-command reproduction**: `cd deepcatch && bash scripts/reproduce_all.sh`
>
> **Honest framing**: this pipeline reproduces the **open-data methods
> benchmark** numbers published in `docs/CROSS_STUDY_BENCHMARK.md`,
> `docs/SENS_AT_SPEC_ABLATION.md`, `docs/SPARSE_AWARE_ABLATION.md`,
> `docs/HARMONIZATION.md`, and the foundation smoke's per-seed JSON.
> It is **NOT** clinical validation, **NOT** external cohort
> validation, and **NOT** a regulatory pipeline.

---

## 0. Prerequisites

- macOS or Linux (Windows works via `reproduce_all.py`)
- Python ≥ 3.10 (project uses 3.14.4)
- `git`, `bash` (4+), `curl` for some downloads
- 17 GB free disk (the FinaleDB + TCGA caches are large but read-only)

The full pipeline (excluding the long cross-study benchmark) takes
~10–15 minutes wall-clock on an M-class MacBook with CPU torch.

---

## 1. Installation

```bash
git clone https://github.com/rollroyces/deepcatch.git
cd deepcatch

# Create venv — required so deps don't conflict with Hermes's global
# Python on macOS (see deepcatch skill pitfall #1).
python3 -m venv .venv
.venv/bin/pip install -e .

# Optional deps for the foundation / GNN steps:
.venv/bin/pip install -e .[api]
```

---

## 2. Pre-flight check

```bash
env -u PYTHONPATH ./.venv/bin/python scripts/verify_environment.py
```

Honors deepcatch skill pitfall **#1 (macOS PYTHONPATH hijack)** by
unstashing the global var before each python invocation.

The script prints a markdown table with a status icon per dependency:

| Status | Check | Found | Required | Note |
|---|---|---|---|---|
| ✅ python≥3.11 | 3.14.4 | 3.11 |  |
| ✅ repo layout | OK | scripts/ + test/ + src/ + docs/ + pyproject.toml |  |
| ✅ deepcatch/.venv | /path/.venv/bin/python | … |  |
| ✅ results/ writable | OK | … |  |
| ✅ FinaleDB features dir | /…/cfdna-fragmentomics-pipeline/data/features (5048 files) | … |  |
| ✅ numpy / scipy / scikit-learn / pandas / torch / pytest | … | … | … |

If any row is ❌, fix it before running the reproduce pipeline. The
script also supports `--json` for machine-readable checks.

For non-Mac users: set `FINALEDB_FEATURES_DIR=/your/path/to/features`
before invoking the reproduce script.

---

## 3. One-command reproduction

```bash
bash scripts/reproduce_all.sh
```

What you get:

- `results/cross_study.json` — pooled cross-study (Jiang + Cristiano) AUC + DeLong CIs
- `results/cross_study_finallydb.json` — alias kept for legacy readers
- `results/per_cancer.json` — per-cancer sens@spec + PPV@prev table
- `results/harmonization.json` — per-study harmonization A/B test
- `results/foundation_smoke.json` — foundation model + LR + naive-avg
  per-seed AUCs (under `--quick` mode, ~30s/run)
- `results/foundation_smoke_ce.json` / `_sens.json` / `_sparse.json` —
  the three smoke variants used by the two ablations
- `results/focal_bce_ablation.json` — CE vs focal-BCE paired t-test
- `results/sparse_aware_ablation.json` — Linear vs SparseAware paired t-test
- `results/pytest_tests.json` — synthesized test summary
- `results/<slug>.log` — per-step stdout/stderr
- `results/REPRODUCE_REPORT.md` — final summary table with wall-clock
- `docs/CROSS_STUDY_BENCHMARK.md` — overwritten by the cross-study step

Honors deepcatch skill pitfalls:

- **#1 (PYTHONPATH hijack)**: every python call is prefixed `env -u PYTHONPATH ./.venv/bin/python`
- **#4 (pipefail)**: the bash script uses `set -euo pipefail`
- **#7 (editable-install paths)**: the repo-root resolver walks up
  until it finds `pyproject.toml + scripts/` together

---

## 4. Wall-clock budget

| Step | Script | Wall (M-class MacBook) |
|---|---|---|
| Pre-flight | `verify_environment.py` | ~3 s |
| Cross-study benchmark | `cross_study_finallydb.py` | ~5 min |
| Per-cancer sens@spec | `per_cancer_sens_at_spec.py --synthetic` | ~30 s |
| Harmonization check | `harmonization_check.py` | ~10 s |
| Foundation smoke (3×) | `foundation_real_smoke.py --quick` ×3 | ~3 min |
| Focal-BCE ablation | `sens_at_spec_ablation.py` | ~30 s |
| Sparse-aware ablation | `sparse_aware_ablation.py` | ~30 s |
| pytest suite | `pytest test/ src/foundation/test_integration.py` | ~60 s |
| **Total** | | **~10–15 min** |

The `--quick` flag cuts the foundation smoke to ~30 s/run (2 seeds ×
3 folds × 2 ensemble). The non-quick defaults take ~5–10 min/run.

---

## 5. Output schema

Each `results/<slug>.json` is the artifact of record — if you change
the schema, you must re-run the pipeline. The scripts never read each
other's outputs at runtime; the only dependency is that the
**ablation** steps read the **foundation smoke** JSONs:

```
focal_bce_ablation.json   reads   foundation_smoke_ce.json + foundation_smoke_sens.json
sparse_aware_ablation.json reads   foundation_smoke.json (CE) + foundation_smoke_sparse.json
```

| JSON | Schema summary | Source doc |
|---|---|---|
| `cross_study.json` | `{"per_cohort":[...], "pooled":{...}, "per_cancer":{...}, "true_confound_control":{...}}` | `docs/CROSS_STUDY_BENCHMARK.md` |
| `per_cancer.json` | `{"per_cancer":[{cancer, n_pos, auc, sens_at_95_ci, sens_at_99_ci, ...}], "pooled":{...}}` | (this file + the per-cancer section in `docs/CROSS_STUDY_BENCHMARK.md`) |
| `harmonization.json` | `{"raw":{pooled_auc_per_seed, ...}, "harmonized":{...}, "delta":{...}, "verdict": "..."}` | `docs/HARMONIZATION.md` |
| `foundation_smoke.json` | `{"foundation_auc_per_seed": [...], "lr_baseline_auc_per_seed": [...], ...}` | the per-seed arrays feed the ablation paired t-tests |
| `focal_bce_ablation.json` | `{"ce_foundation_aucs": [...], "sens_foundation_aucs": [...], "paired_diff": float, "p": float}` | `docs/SENS_AT_SPEC_ABLATION.md` |
| `sparse_aware_ablation.json` | same shape as focal-BCE | `docs/SPARSE_AWARE_ABLATION.md` |
| `pytest_tests.json` | `{"step":"pytest_tests", "tests_passed": N, "deselected": [...]}` | the per-test JSONs under each test module are the canonical contract |

---

## 6. Re-running individual steps

Every script accepts a `--out` flag (or hardcodes a path you can copy
after the run) and can be invoked standalone:

```bash
# Cross-study benchmark alone
env -u PYTHONPATH ./.venv/bin/python scripts/cross_study_finallydb.py \
    --features-dir $FINALEDB_FEATURES_DIR \
    --out-json results/cross_study.json \
    --out-md docs/CROSS_STUDY_BENCHMARK.md

# Per-cancer sens@spec alone (synthetic fixture)
env -u PYTHONPATH ./.venv/bin/python scripts/per_cancer_sens_at_spec.py \
    --synthetic --n 240 --out results/per_cancer.json

# Foundation smoke alone (real TCGA panel-LLR if TCGA_CACHE_DIR set)
env -u PYTHONPATH ./.venv/bin/python scripts/foundation_real_smoke.py \
    --quick --n-patients 20 --seeds 5 \
    --features-dir validation/tcga/tcga_cache \
    --out results/foundation_smoke.json

# Focal-BCE ablation
env -u PYTHONPATH ./.venv/bin/python scripts/sens_at_spec_ablation.py \
    --ce-json results/foundation_smoke_ce.json \
    --sens-json results/foundation_smoke_sens.json \
    --out results/focal_bce_ablation.json
```

---

## 7. Selective re-runs

The bash script accepts subset / skip flags:

```bash
# Just the pre-flight + cross-study benchmark + report
bash scripts/reproduce_all.sh --only verify_environment,cross_study

# Skip the slow pytest step
bash scripts/reproduce_all.sh --skip-tests

# Print the plan without running anything
bash scripts/reproduce_all.sh --dry-run
```

The Python equivalent (`scripts/reproduce_all.py`) accepts the same
flags:

```bash
env -u PYTHONPATH ./.venv/bin/python scripts/reproduce_all.py --only per_cancer,harmonization --dry-run
```

---

## 8. Cross-platform (Windows / Alpine)

If `bash` is not on the PATH, use the Python wrapper:

```bash
env -u PYTHONPATH ./.venv/bin/python scripts/reproduce_all.py
```

It runs each script via `subprocess.run` with explicit
`env={"PYTHONPATH": ""}` per deepcatch skill pitfall #1 and produces
the same `results/REPRODUCE_REPORT.md` artifact.

---

## 9. What this pipeline does NOT do

- **Clinical validation**: this is open-data methods benchmarking on
  FinaleDB publications + TCGA-LUAD. There is no held-out clinical
  cohort.
- **Regulatory-grade QC**: no ISO 13485 / 510(k) controls.
- **Per-sample audit trail**: only the per-seed JSONs are preserved.
- **Network downloads**: the cross-study step reads a pre-existing
  FinaleDB feature cache (set `FINALEDB_FEATURES_DIR` to override).

If you need a clinical-grade validation pipeline, you need an
external held-out cohort + a regulatory sponsor. This repo's
contribution is the **open-data methods benchmark** the published
numbers come from.

---

## 10. Troubleshooting checklist

| Symptom | Fix |
|---|---|
| `ModuleNotFoundError: numpy` despite `pip install` | You forgot `env -u PYTHONPATH`. Hermes's global PYTHONPATH hijacks the venv (deepcatch skill pitfall #1). |
| `bash: declare: -A: invalid option` on macOS | `/bin/bash` is 3.2.57 on macOS. The reproduce script uses bash 3.2-safe parallel arrays — verify you're on the latest version. |
| `feature cache not found` in cross-study step | Set `FINALEDB_FEATURES_DIR=/path/to/cfdna-fragmentomics-pipeline/data/features`. |
| pytest `Timeout (>60.0s)` on `test_foundation_smoke.py` | Those tests spawn `foundation_real_smoke.py` as a subprocess (~80s each). The reproduce script uses `--timeout=180` for them. |
| `cross_study_finallydb.py exceeded 900s` in pytest | One test (`test_end_to_end_script_writes_both_outputs`) runs the cross-study script as a subprocess; the reproduce script deselects it (`--deselect`) and runs it standalone instead. |
| Foundation AUC ≈ 0.69 but expected 0.85 | The `--quick` mode caps seeds/folds/ensemble — this is the trade-off for fitting in the 10–15 min budget. Run with default `--n-folds 5 --seeds 5 --n-ensemble 3` for the published numbers (~5–10 min/run). |
| Finaledb pretrained loader test skips | The `checkpoints/foundation_pretrained_finaledb_PRODUCTION.pt` file is gitignored. Re-run `scripts/pretrain_production_finaledb.py` (~2 s on CPU) to regenerate. |

---

## 11. What "reproduces" means here

The pipeline regenerates:

1. **All headline numbers in `docs/CROSS_STUDY_BENCHMARK.md`** (5-seed
   × 5-fold pooled OOF + per-cancer sens@spec with DeLong CIs).
2. **The harmonization verdict in `docs/HARMONIZATION.md`** (raw vs
   per-study z-score A/B test).
3. **The focal-BCE vs CE verdict in `docs/SENS_AT_SPEC_ABLATION.md`**
   (paired t-test on per-seed arrays).
4. **The sparse-aware vs Linear verdict in
   `docs/SPARSE_AWARE_ABLATION.md`** (paired t-test on per-seed arrays).
5. **The foundation smoke's per-seed arrays** (the JSON is the
   artifact; the docs cite the AUC range).

It does **not** regenerate:

- The CADD weighting papers (`docs/CADD_*.md`) — those use the GDC
  API and a 382-patient WXS cohort; the reproduce pipeline covers
  the cross-study FinaleDB only.
- The methylation channel — that lives in a sibling repo
  (`rollroyces/deepcatch-methylation`).

---

## 12. License & data

- **Code**: MIT (see `LICENSE`).
- **FinaleDB open data**: Jiang 2015 (pub 6) + Cristiano 2019 (pub 8)
  are the two cross-study cohorts. Pre-extracted under
  `/Users/hermes/cfdna-fragmentomics-pipeline/data/features/`; the
  `FINALEDB_FEATURES_DIR` env var overrides the path for non-Mac
  setups.
- **TCGA-LUAD MAFs**: GDC open-access somatic mutations, cached under
  `validation/tcga/tcga_cache/` (~70 MB).

---

## 13. One-command TL;DR

```bash
cd deepcatch
env -u PYTHONPATH ./.venv/bin/python scripts/verify_environment.py
bash scripts/reproduce_all.sh
cat results/REPRODUCE_REPORT.md
```