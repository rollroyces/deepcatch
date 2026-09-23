# Reproduction report

- Generated: 2026-09-23T08:12:54Z
- repo_root: `/Users/hermes/deepcatch`
- python:    `/Users/hermes/deepcatch/.venv/bin/python` (Python 3.14.4)
- total wall: 0.0s

## Step results

| Step | Description | Wall (s) | Exit | Output |
|---|---|---:|---|---|
| 🟡 `cross_study` | cross_study_finallydb.py (FinaleDB Jiang+Cristiano) | 0.0 | DRY | `cross_study.json (MISSING)` |
| 🟡 `per_cancer` | per_cancer_sens_at_spec.py --synthetic | 0.0 | DRY | `per_cancer.json` |
| 🟡 `harmonization` | harmonization_check.py | 0.0 | DRY | `harmonization.json` |
| 🟡 `focal_bce_ablation` | sens_at_spec_ablation.py (CE vs focal-BCE) | 0.0 | DRY | `focal_bce_ablation.json` |
| 🟡 `sparse_aware_ablation` | sparse_aware_ablation.py (Linear vs SparseAware) | 0.0 | DRY | `sparse_aware_ablation.json` |
| 🟡 `foundation_smoke` | foundation_real_smoke.py --quick | 0.0 | DRY | `foundation_smoke.json` |
| 🟡 `pytest_tests` | pytest test/ + src/foundation/test_integration.py | 0.0 | DRY | `pytest_tests.json (MISSING)` |

## Per-step logs

Each step writes `results/<slug>.log` with the full stdout/stderr of the script.

## Honest framing

This script reproduces the **open-data methods benchmark** numbers published in
`docs/CROSS_STUDY_BENCHMARK.md`, `docs/SENS_AT_SPEC_ABLATION.md`,
`docs/SPARSE_AWARE_ABLATION.md`, `docs/HARMONIZATION.md`, and the foundation
smoke's `docs/FOUNDATION_REAL_SMOKE.md`-shaped numbers.

It does **not** perform:
- Clinical validation on an external held-out cohort
- Regulatory-grade quality controls
- Per-sample audit/reproducibility beyond what the per-seed JSONs capture

## Re-running individual steps

```bash
# Cross-study benchmark alone
env -u PYTHONPATH ./.venv/bin/python scripts/cross_study_finallydb.py \
    --features-dir $FINALEDB_FEATURES_DIR \
    --out-json results/cross_study.json \
    --out-md docs/CROSS_STUDY_BENCHMARK.md

# Per-cancer sens@spec alone (synthetic fixture)
env -u PYTHONPATH ./.venv/bin/python scripts/per_cancer_sens_at_spec.py \
    --synthetic --n 240 --out results/per_cancer.json
```

