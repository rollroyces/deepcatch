#!/usr/bin/env bash
# DeepCatch — methods-paper reproduction driver
# =============================================
# Re-runs every script whose JSON artifact is quoted in paper/METHODS_PAPER.md §3.
# Writes JSONs to results/, refreshes docs/CROSS_STUDY_BENCHMARK.md,
# then runs the regression test gate.
#
# Usage:
#   bash paper/REPRODUCE.sh              # full run
#   bash paper/REPRODUCE.sh --quick      # reduced seed / PCA list (CI sanity check)
#   bash paper/REPRODUCE.sh --no-test    # skip the pytest gate
#
# Exit codes:
#   0 = all artifacts produced + tests green
#   nonzero on any failure (set -euo pipefail)
#
# Requirements:
#   - Python 3.x with the repo venv at .venv/
#   - FinaleDB features already pre-extracted to
#     /Users/hermes/cfdna-fragmentomics-pipeline/data/features/
#     (FinaleDB REST API is in degraded state as of 2026-09-21)

set -euo pipefail

# Locate the repo root (parent of paper/, where this script lives).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

# CLI flags
QUICK=0
NO_TEST=0
for arg in "$@"; do
    case "$arg" in
        --quick) QUICK=1 ;;
        --no-test) NO_TEST=1 ;;
        -h|--help)
            head -25 "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *)
            echo "[REPRODUCE] Unknown flag: $arg" >&2
            exit 2
            ;;
    esac
done

# Use the repo venv if it exists; otherwise fall back to whatever python is on PATH.
if [[ -x ".venv/bin/python" ]]; then
    PYTHON="./.venv/bin/python"
else
    PYTHON="$(command -v python3)"
fi
echo "[REPRODUCE] Using python: ${PYTHON}"

# Clear PYTHONPATH to avoid leaking an outer project's modules into DeepCatch.
export PYTHONPATH=

# Seeds: full = [42, 13, 7, 99, 1234]; --quick = [42, 7, 1234]
SEEDS=(42 13 7 99 1234)
if [[ "${QUICK}" == "1" ]]; then
    SEEDS=(42 7 1234)
fi

mkdir -p results docs

echo "================================================================"
echo "[1/3] Cross-study FinaleDB sweep (cross_study_finallydb.py)"
echo "      seeds: ${SEEDS[*]}"
echo "================================================================"

# Build the arg list.
SEED_ARGS=()
for s in "${SEEDS[@]}"; do SEED_ARGS+=(--seeds "$s"); done

# Cross-study benchmark. The `--out-per-cancer-json` flag writes a
# standalone JSON with per-cancer sens@spec (DeLong CIs + PPV@prev).
# The flag's emit path was previously buggy (n_neg summed across cancers
# — 7×264 = 1848 for the harmonized pooled cohort). Fixed upstream:
# each OvR uses the FULL pooled healthy count (~264), the POOLED row
# uses the actual pooled n_pos/n_neg from the cohort section.
"${PYTHON}" scripts/cross_study_finallydb.py \
    --features-dir /Users/hermes/cfdna-fragmentomics-pipeline/data/features \
    --labels-multiclass /Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv \
    --out-json    results/cross_study_finallydb.json \
    --out-md      docs/CROSS_STUDY_BENCHMARK.md \
    --out-per-cancer-json results/per_cancer_sens_at_spec.json \
    --top-cancer-n 5 \
    "${SEED_ARGS[@]}"

echo
echo "[REPRODUCE] Wrote:"
echo "  - results/cross_study_finallydb.json"
echo "  - results/per_cancer_sens_at_spec.json"
echo "  - docs/CROSS_STUDY_BENCHMARK.md"

# Acceptance floors (sanity; both must hold for the §3 numbers to be cite-able).
echo
echo "[SANITY] Acceptance floors:"
HARMONIZED_AUC=$( "${PYTHON}" -c "
import json
with open('results/cross_study_finallydb.json') as f:
    j = json.load(f)
print(j['pooled']['harmonized']['auc_mean'])
" )
TCC_AUC=$( "${PYTHON}" -c "
import json
with open('results/cross_study_finallydb.json') as f:
    j = json.load(f)
print(j['true_confound_control']['cancer_jiang_healthy_cristiano']['harmonized']['auc_mean'])
" )

H_OK=$( "${PYTHON}" -c "print(1 if 0.965 <= float('${HARMONIZED_AUC}') <= 0.985 else 0)" )
T_OK=$( "${PYTHON}" -c "print(1 if 0.48 <= float('${TCC_AUC}') <= 0.51 else 0)" )
echo "  pooled-harmonized auc_mean = ${HARMONIZED_AUC}  (floor [0.965, 0.985]): ${H_OK}"
echo "  true-confound-harmonized = ${TCC_AUC}            (floor [0.480, 0.510]): ${T_OK}"
if [[ "${H_OK}" != "1" || "${T_OK}" != "1" ]]; then
    echo "[REPRODUCE] !! Acceptance floor failed; the numbers in paper §3.2 are not cite-able yet." >&2
    exit 3
fi

if [[ "${NO_TEST}" == "1" ]]; then
    echo
    echo "[REPRODUCE] --no-test set; skipping pytest gate."
    exit 0
fi

echo
echo "================================================================"
echo "[2/3] Pytest gate (7-file cross-study + per-cancer + pretrain + harmony + projection)"
echo "================================================================"

"${PYTHON}" -m pytest \
    test/test_biomedical_review_fixes.py \
    test/test_sparse_aware_projection.py \
    test/test_finaledb_pretrained_loader.py \
    test/test_pretrain_bug_fix.py \
    test/test_harmonization_check.py \
    test/test_per_cancer_sens_at_spec.py \
    src/foundation/test_integration.py \
    --timeout=60 -q

echo
echo "================================================================"
echo "[3/3] Done."
echo "================================================================"
echo "  results/cross_study_finallydb.json     (§3.2 source of truth)"
echo "  results/per_cancer_sens_at_spec.json   (§3.3 source of truth)"
echo "  docs/CROSS_STUDY_BENCHMARK.md          (companion narrative)"
echo "  paper/METHODS_PAPER.md §3              (the paper)"
echo "================================================================"
