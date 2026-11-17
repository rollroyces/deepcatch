#!/usr/bin/env bash
# DeepCatch — methods-paper reproduction driver
# =============================================
# Re-runs every script whose JSON artifact is quoted in paper/METHODS_PAPER.md
# and paper/paper.tex, then PROVES that the regenerated numbers match the
# committed ones, then runs the regression gate.
#
# The proof step matters. Re-running a script is not the same as showing the
# committed numbers are reproducible: a script can succeed and still write a
# different number than the one in the manuscript. `scripts/reproduce_check.py`
# compares the freshly generated artifacts against the committed ones number
# by number, so a silent drift fails the run instead of waiting for a reader
# to notice it.
#
# Usage:
#   bash paper/REPRODUCE.sh              # full run (5 seeds), writes into the repo
#   bash paper/REPRODUCE.sh --quick      # 3 seeds, writes to a scratch dir
#   bash paper/REPRODUCE.sh --no-test    # skip the pytest gate
#
# Exit codes:
#   0 = every artifact regenerated, every quoted number reproduced, tests green
#   nonzero on any failure (set -euo pipefail)
#
# Requirements:
#   - Python 3.x with the repo venv at .venv/ (a system python3 also works)
#   - FinaleDB features already pre-extracted to $DEEPCATCH_FEATURES_DIR
#     (the FinaleDB REST API is in a degraded state as of 2026-09-21)
#   - the held-out external cohort fetched into $DEEPCATCH_EXTERNAL_DIR
#     (scripts/fetch_external_cohort.py — public API + S3, no IRB needed)
#   - for the methylation channel only: java + FinaleMe 0.58 + the cached
#     per-sample CpG matrices. If any of those is missing the step is
#     SKIPPED and says so. It is not a failure: the other channels are
#     independent of it.
#
# Environment overrides:
#   DEEPCATCH_FEATURES_DIR    directory of extracted FinaleDB feature matrices
#   DEEPCATCH_LABELS          path to labels_multiclass.tsv
#   DEEPCATCH_EXTERNAL_DIR    directory holding the held-out external cohort
#   DEEPCATCH_METHYLATION_DIR directory holding matrices/ + cohort_labels.tsv
#   DEEPCATCH_SISTER_PIPELINE checkout of cfdna-fragmentomics-pipeline
#
# The sister-pipeline checkout is a REAL prerequisite, not a convenience:
# the cross-study, LOSO and external-validation steps all import the
# sibling's `honest_benchmark.load5` and `train_classifier._harmonize`, so
# that the training and external features are assembled by exactly the
# same code path. Clone it next to this repository:
#
#   git clone https://github.com/rollroyces/cfdna-fragmentomics-pipeline
#
# and point the driver at it if it is not in the default place:
#
#   export DEEPCATCH_SISTER_PIPELINE=/path/to/cfdna-fragmentomics-pipeline
#
# Runtime (measured on the development machine, full mode): ~20 minutes, of
# which the two cross-study sweeps are ~12 minutes.

set -euo pipefail

# Locate the repo root (parent of paper/, where this script lives).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

# ---------------------------------------------------------------- flags ----
QUICK=0
NO_TEST=0
for arg in "$@"; do
    case "$arg" in
        --quick) QUICK=1 ;;
        --no-test) NO_TEST=1 ;;
        -h|--help)
            head -45 "${BASH_SOURCE[0]}"
            exit 0
            ;;
        *)
            echo "[REPRODUCE] Unknown flag: $arg" >&2
            exit 2
            ;;
    esac
done

# -------------------------------------------------------------- python ----
if [[ -x ".venv/bin/python" ]]; then
    PYTHON="./.venv/bin/python"
else
    PYTHON="$(command -v python3)"
fi
echo "[REPRODUCE] Using python: ${PYTHON}"

# Clear PYTHONPATH to avoid leaking an outer project's modules into DeepCatch.
export PYTHONPATH=

# --------------------------------------------------------------- paths ----
FEATURES_DIR="${DEEPCATCH_FEATURES_DIR:-/Users/hermes/cfdna-fragmentomics-pipeline/data/features}"
LABELS_TSV="${DEEPCATCH_LABELS:-/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv}"
EXTERNAL_DIR="${DEEPCATCH_EXTERNAL_DIR:-${REPO_ROOT}/data/external_cohort}"
METHYLATION_DIR="${DEEPCATCH_METHYLATION_DIR:-${REPO_ROOT}/../methylation_cohort}"

# One temp root for everything this run needs but must not leave behind:
# the baseline snapshot, the methylation work dir, and (in quick mode) the
# scratch output tree.
RUN_TMP="$(mktemp -d -t deepcatch-reproduce-XXXXXX)"
BASELINE_DIR="${RUN_TMP}/baseline"
METH_WORK="${RUN_TMP}/methylation_work"
mkdir -p "${BASELINE_DIR}"

# Full mode regenerates the committed artifacts in place. Quick mode writes
# to a scratch directory so a smoke run can never overwrite the published
# numbers with a reduced-seed approximation of them.
SCRATCH=""
if [[ "${QUICK}" == "1" ]]; then
    SCRATCH="${RUN_TMP}/scratch"
    RESULTS_DIR="${SCRATCH}/results"
    DOCS_DIR="${SCRATCH}/docs"
    FIG_DIR="${SCRATCH}/docs/figures"
    echo "[REPRODUCE] quick mode: writing to ${SCRATCH} (the repo is not touched)"
else
    RESULTS_DIR="${REPO_ROOT}/results"
    DOCS_DIR="${REPO_ROOT}/docs"
    FIG_DIR="${DOCS_DIR}/figures"
fi
mkdir -p "${RESULTS_DIR}" "${DOCS_DIR}" "${FIG_DIR}"

# Seeds: full = [42, 13, 7, 99, 1234] (the list the committed artifacts were
# produced with); --quick = [42, 7, 1234] (a subset, for speed).
SEEDS=(42 13 7 99 1234)
REPRO_TOLERANCE=1e-6
REPRO_THRESHOLD_TOLERANCE=1e-6
if [[ "${QUICK}" == "1" ]]; then
    SEEDS=(42 7 1234)
    # A 3-seed mean is not the 5-seed mean, so the committed numbers are not
    # expected to match to the last digit in quick mode. 5e-3 absorbs seed
    # sampling for the smooth numbers (AUCs move by ~3e-4 between a 3-seed
    # and a 5-seed run) and still catches a real regression.
    REPRO_TOLERANCE=5e-3
    # sens@spec operating points are step functions of the score
    # distribution -- the threshold is fixed by the healthy tail, so one
    # healthy sample moves the operating point and one positive moves the
    # sensitivity by 1/n_pos. They are an order of magnitude more
    # seed-sensitive than the AUC and need their own band. The full run
    # (identical seed list) reproduces them exactly.
    REPRO_THRESHOLD_TOLERANCE=5e-2
fi
# NOTE: the seed list must be ONE argument (`--seeds 42 13 7 99 1234`),
# not one argument per seed. The scripts use argparse `nargs="+"`, and
# argparse keeps the LAST occurrence of a repeated option -- so
# `--seeds 42 --seeds 7 --seeds 1234` silently runs ONE seed and reports
# a ±0.0000 std, which looks like a suspiciously perfect result.
SEED_ARGS=(--seeds "${SEEDS[@]}")

# The artifacts this driver regenerates. The committed copy of each is
# snapshotted before the run and compared against afterwards.
ARTIFACTS=(
    cross_study_finallydb.json
    cross_study_finallydb_gc_corrected.json
    cross_study_finallydb_with_motifs.json
    per_cancer_sens_at_spec.json
    cross_study_finallydb_shuffled_control.json
    cross_study_loso.json
    external_validation_finaledb.json
    fusion_ablation_isotonic.json
    fusion_regularized_smoke.json
    real_tcga_validation.json
    decision_curve.json
    methylation_cohort_channel.json
    sens_at_spec_ce_n20.json
    ultraearly_readiness.json
)

# Cleanup: only ever remove the temp root this script created, and only
# while it still carries the mktemp prefix it was made with.
cleanup() {
    case "${RUN_TMP:-}" in
        */deepcatch-reproduce-*)
            [ -n "${RUN_TMP}" ] && [ -d "${RUN_TMP}" ] && rm -rf "${RUN_TMP}" || true
            ;;
    esac
    return 0
}
trap cleanup EXIT

for a in "${ARTIFACTS[@]}"; do
    if [[ -f "${REPO_ROOT}/results/${a}" ]]; then
        cp "${REPO_ROOT}/results/${a}" "${BASELINE_DIR}/"
    fi
done
echo "[REPRODUCE] Snapshotted $(ls -1 "${BASELINE_DIR}" | wc -l | tr -d ' ') committed artifact(s) as the reproduction baseline"

if [[ ! -d "${FEATURES_DIR}" ]]; then
    echo "[REPRODUCE] !! feature directory not found: ${FEATURES_DIR}" >&2
    echo "[REPRODUCE]    set DEEPCATCH_FEATURES_DIR to the extracted FinaleDB features." >&2
    exit 4
fi

# The sister-pipeline checkout is a hard prerequisite: the cross-study,
# LOSO and external-validation steps import its load5 and _harmonize.
# Without the check below, step 1 dies with "No module named
# 'honest_benchmark'" -- which reads like a broken repo rather than a
# missing prerequisite.
SISTER_PIPELINE="${DEEPCATCH_SISTER_PIPELINE:-/Users/hermes/cfdna-fragmentomics-pipeline}"
if [[ ! -f "${SISTER_PIPELINE}/scripts/honest_benchmark.py" ]]; then
    echo "[REPRODUCE] !! sister-pipeline checkout not found: ${SISTER_PIPELINE}" >&2
    echo "[REPRODUCE]    clone https://github.com/rollroyces/cfdna-fragmentomics-pipeline" >&2
    echo "[REPRODUCE]    and set DEEPCATCH_SISTER_PIPELINE to its path." >&2
    exit 4
fi

echo
echo "================================================================"
echo "[1/14] Cross-study FinaleDB sweep — legacy (no GC correction)"
echo "================================================================"
# `--gc-correction` DEFAULTS TO ON, so a driver that passes no flag writes
# the GC-corrected run into the legacy path. That silently replaced the
# legacy baseline and broke the comparison the paper quotes; hence the two
# explicit runs below.
#
# `--out-per-cancer-json ""` on THIS run is load-bearing. The script's
# default for that flag is the published `results/per_cancer_sens_at_spec.json`,
# which belongs to the GC-corrected run below. Leaving the flag off here
# made the legacy sweep overwrite the committed artifact with legacy
# per-cancer numbers (BRCA 0.9738 -> 0.9794) — a silent drift of the very
# artifact the paper quotes. The legacy run's own per-cancer numbers live in
# its main JSON (`per_cancer`), which is what the legacy narrative quotes.
echo "      seeds: ${SEEDS[*]}"
"${PYTHON}" scripts/cross_study_finallydb.py \
    --features-dir "${FEATURES_DIR}" \
    --labels-multiclass "${LABELS_TSV}" \
    --no-gc-correction \
    --out-json    "${RESULTS_DIR}/cross_study_finallydb.json" \
    --out-md      "${DOCS_DIR}/CROSS_STUDY_BENCHMARK.md" \
    --out-per-cancer-json "" \
    --shuffled-label-control-out "" \
    --top-cancer-n 5 \
    "${SEED_ARGS[@]}"

echo
echo "================================================================"
echo "[2/14] Cross-study FinaleDB sweep — GC-corrected (headline)"
echo "================================================================"
# `--include-screening` is deliberately OFF here. The committed
# `results/per_cancer_sens_at_spec.json` was produced on the legacy
# grid {0.95, 0.98, 0.99}; turning the flag on adds sens_at_995 /
# sens_at_999 fields, so the regenerated artifact would no longer match
# the committed one. The screening-grade grid is a separate route —
# `scripts/per_cancer_sens_at_spec.py --specificities 0.95 0.98 0.99
# 0.995 0.999` (see ULTRA_EARLY_READINESS.md §1.2) — and its pooled
# numbers are already carried by results/cross_study_loso.json and
# results/external_validation_finaledb.json.
"${PYTHON}" scripts/cross_study_finallydb.py \
    --features-dir "${FEATURES_DIR}" \
    --labels-multiclass "${LABELS_TSV}" \
    --gc-correction \
    --out-json    "${RESULTS_DIR}/cross_study_finallydb_gc_corrected.json" \
    --out-md      "${DOCS_DIR}/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md" \
    --out-per-cancer-json "${RESULTS_DIR}/per_cancer_sens_at_spec.json" \
    --shuffled-label-control-out "" \
    --top-cancer-n 5 \
    "${SEED_ARGS[@]}"

echo
echo "================================================================"
echo "[3/14] Cross-study FinaleDB sweep — with 4-mer end motifs"
echo "================================================================"
# The third tumor-naive benchmark: the same pipeline with the 256-dim
# 4-mer end-motif vector added. Eight documents quote its pooled AUC,
# so the driver regenerates it like the other two sweeps.
#
# The motif features (`<sample>.motifs.npy`) are NOT part of the feature
# cache that ships with the sibling checkout: they live in the author's
# local feature directory and are tracked in neither repository. Without
# them the motif block is zero-filled for every sample and the run
# produces a DIFFERENT benchmark (the committed artifact was produced
# with 98 cached motif files), so the step is skipped and says so rather
# than writing a number that cannot match.
MOTIF_COUNT=$(find "${FEATURES_DIR}" -maxdepth 1 -name '*.motifs.npy' 2>/dev/null | wc -l | tr -d ' ')
if [[ "${MOTIF_COUNT}" == "0" ]]; then
    echo "[REPRODUCE] SKIP: no *.motifs.npy in ${FEATURES_DIR}"
    echo "[REPRODUCE]       the with-motifs benchmark cannot be regenerated"
    echo "[REPRODUCE]       without the cached 4-mer end-motif features, and"
    echo "[REPRODUCE]       the reproduction check below will report"
    echo "[REPRODUCE]       results/cross_study_finallydb_with_motifs.json as"
    echo "[REPRODUCE]       NOT regenerated and end with 'the reproduction"
    echo "[REPRODUCE]       is incomplete' -- that is the correct verdict."
else
    echo "[REPRODUCE]       ${MOTIF_COUNT} motif feature file(s) found"
    "${PYTHON}" scripts/cross_study_finallydb.py \
        --features-dir "${FEATURES_DIR}" \
        --labels-multiclass "${LABELS_TSV}" \
        --include-motifs \
        --out-json    "${RESULTS_DIR}/cross_study_finallydb_with_motifs.json" \
        --out-md      "${DOCS_DIR}/CROSS_STUDY_BENCHMARK_with_motifs.md" \
        --out-per-cancer-json "" \
        --shuffled-label-control-out "" \
        --top-cancer-n 5 \
        "${SEED_ARGS[@]}"
fi

echo
echo "================================================================"
echo "[4/14] Shuffled-label null control"
echo "================================================================"
# Single seed (2), matching the committed artifact: the control is a sanity
# check, and a second seed adds minutes without adding evidence.
"${PYTHON}" scripts/cross_study_finallydb_shuffled_control.py \
    --features-dir "${FEATURES_DIR}" \
    --labels-multiclass "${LABELS_TSV}" \
    --seeds 2 \
    --out-json "${RESULTS_DIR}/cross_study_finallydb_shuffled_control.json"

echo
echo "================================================================"
echo "[5/14] Leave-one-study-out transfer (train on one publication)"
echo "================================================================"
"${PYTHON}" scripts/cross_study_loso.py \
    --features-dir "${FEATURES_DIR}" \
    --labels-multiclass "${LABELS_TSV}" \
    --publications 6 8 \
    --pca 200 \
    "${SEED_ARGS[@]}" \
    --out-json "${RESULTS_DIR}/cross_study_loso.json" \
    --out-md   "${DOCS_DIR}/CROSS_STUDY_LOSO.md"

echo
echo "================================================================"
echo "[6/14] External validation on held-out public publications"
echo "================================================================"
if [[ ! -d "${EXTERNAL_DIR}" ]]; then
    echo "[REPRODUCE] !! external cohort not found: ${EXTERNAL_DIR}" >&2
    echo "[REPRODUCE]    run scripts/fetch_external_cohort.py first (public API + S3)." >&2
    exit 4
fi
"${PYTHON}" scripts/external_validation_finaledb.py \
    --features-dir "${FEATURES_DIR}" \
    --labels-multiclass "${LABELS_TSV}" \
    --external-dir "${EXTERNAL_DIR}" \
    --pca 200 \
    --out-json "${RESULTS_DIR}/external_validation_finaledb.json" \
    --out-md   "${DOCS_DIR}/EXTERNAL_VALIDATION.md"

echo
echo "================================================================"
echo "[7/14] Fusion ablation (tumor-naive + mutation-informed)"
echo "================================================================"
PYTHONPATH=src "${PYTHON}" src/fragmentomics/fusion_ablation.py \
    --features-dir "${FEATURES_DIR}" \
    --labels "${LABELS_TSV}" \
    --pca-n 200 \
    --seeds 10 \
    --target-auc 0.92 \
    --out "${RESULTS_DIR}/fusion_ablation_isotonic.json"

echo
echo "================================================================"
echo "[8/14] Regularized late fusion / Bayesian model averaging"
echo "================================================================"
# The estimator that replaced the neural encoder, with the permuted-
# label and arm-swap controls it is reported with. Ten documents quote
# its 0.9825, so the driver regenerates it like the other benchmarks.
"${PYTHON}" scripts/fusion_regularized_smoke.py \
    --out-json "${RESULTS_DIR}/fusion_regularized_smoke.json"

echo
echo "================================================================"
echo "[9/14] Decision curve (clinician-facing operating points)"
echo "================================================================"
# The operating points RESULTS.md section 4 quotes. Same cohort, same
# PCA (200 components) and the same five seeds as the decision curve
# was always run with -- but the PCA is seeded now, so this output is
# reproducible instead of one arbitrary draw.
PYTHONPATH=src "${PYTHON}" src/fragmentomics/decision_curve_cli.py \
    --features-dir "${FEATURES_DIR}" \
    --labels "${LABELS_TSV}" \
    --pca-n 200 \
    --seeds 5 \
    --out "${RESULTS_DIR}/decision_curve.json"

echo
echo "================================================================"
echo "[10/14] Methylation channel (FinaleMe on real public fragments)"
echo "================================================================"
# This step needs java, the FinaleMe 0.58 jar and the cached per-sample CpG
# matrices. When any of those is missing the step is skipped and reported
# as skipped — the other channels do not depend on it, and a missing local
# tool is not a reproducibility failure.
if [[ ! -d "${METHYLATION_DIR}/matrices" ]]; then
    echo "[REPRODUCE] SKIP: no cached CpG matrices at ${METHYLATION_DIR}/matrices"
    echo "[REPRODUCE]       (the methylation channel needs the FinaleMe BAM+matrix stage;"
    echo "[REPRODUCE]        see docs/METHYLATION_PIPELINE.md §4)"
    echo "[REPRODUCE]       the reproduction check below will report the methylation"
    echo "[REPRODUCE]       channel as NOT regenerated — that is the correct verdict."
elif ! command -v java >/dev/null 2>&1; then
    echo "[REPRODUCE] SKIP: java is not on PATH (needed for FinaleMe)"
    echo "[REPRODUCE]       the reproduction check below will report the methylation"
    echo "[REPRODUCE]       channel as NOT regenerated — that is the correct verdict."
else
    # The work dir is a pure output dir (the HMM is trained fresh into it
    # and every sample is decoded into it), so it is pointed at a run-local
    # directory rather than at the shared cache: a reproduction run must not
    # mutate the cache it reads from.
    "${PYTHON}" scripts/methylation_cohort_channel.py \
        --matrices-dir "${METHYLATION_DIR}/matrices" \
        --labels       "${METHYLATION_DIR}/cohort_labels.tsv" \
        --work-dir     "${METH_WORK}" \
        --out-json     "${RESULTS_DIR}/methylation_cohort_channel.json"
fi

echo
echo "================================================================"
echo "[11/14] Mutation-informed panel validation (real TCGA-LUAD mutations)"
echo "================================================================"
# The panel-LLR AUC (0.921 at 0.1% ctDNA) is the first figure in the
# preprint's abstract and the readiness verdict depends on it. It needs
# the GDC masked-MAF cache; without it the step is skipped and says so,
# and the reproduction check then reports that artifact as NOT
# regenerated.
MAF_COUNT=$(find "${REPO_ROOT}/validation/tcga/tcga_cache" -maxdepth 1 \
    -name 'gdc_TCGA-LUAD_*.maf.gz' 2>/dev/null | wc -l | tr -d ' ')
if [[ "${MAF_COUNT}" == "0" ]]; then
    echo "[REPRODUCE] SKIP: no gdc_TCGA-LUAD_*.maf.gz in"
    echo "[REPRODUCE]       ${REPO_ROOT}/validation/tcga/tcga_cache — the GDC"
    echo "[REPRODUCE]       masked-MAF cache is not present, so the panel"
    echo "[REPRODUCE]       validation cannot be regenerated here."
else
    echo "[REPRODUCE]       ${MAF_COUNT} MAF file(s) in the GDC cache"
    "${PYTHON}" real_tcga_validation.py \
        --no-download \
        --cache-dir "${REPO_ROOT}/validation/tcga/tcga_cache" \
        --output "${RESULTS_DIR}/real_tcga_validation.json" \
        --seeds 5 --n-patients 20
fi

echo
echo "================================================================"
echo "[12/14] Figures"
echo "================================================================"
"${PYTHON}" scripts/plot_figures.py \
    --results-dir "${RESULTS_DIR}" \
    --out-dir "${FIG_DIR}"

echo
echo "================================================================"
echo "[13/14] Acceptance floors"
echo "================================================================"
read_harm() {
    "${PYTHON}" -c "
import json, sys
with open(sys.argv[1]) as f:
    j = json.load(f)
print(j['pooled']['harmonized']['auc_mean'])
" "$1"
}
read_confound() {
    "${PYTHON}" -c "
import json, sys
with open(sys.argv[1]) as f:
    j = json.load(f)
print(j['true_confound_control']['cancer_jiang_healthy_cristiano']['harmonized']['auc_mean'])
" "$1"
}

LEGACY_AUC="$(read_harm "${RESULTS_DIR}/cross_study_finallydb.json")"
LEGACY_CONF="$(read_confound "${RESULTS_DIR}/cross_study_finallydb.json")"
GC_AUC="$(read_harm "${RESULTS_DIR}/cross_study_finallydb_gc_corrected.json")"

LEGACY_OK=$("${PYTHON}" -c "print(1 if 0.965 <= float('${LEGACY_AUC}') <= 0.985 else 0)")
CONF_OK=$("${PYTHON}" -c "print(1 if 0.48 <= float('${LEGACY_CONF}') <= 0.51 else 0)")
GC_OK=$("${PYTHON}" -c "print(1 if 0.960 <= float('${GC_AUC}') <= 0.980 else 0)")
echo "  legacy  pooled-harmonized AUC = ${LEGACY_AUC}  (floor [0.965, 0.985]): ${LEGACY_OK}"
echo "  legacy  true-confound AUC     = ${LEGACY_CONF}  (floor [0.480, 0.510]): ${CONF_OK}"
echo "  GC-corr pooled-harmonized AUC = ${GC_AUC}  (floor [0.960, 0.980]): ${GC_OK}"
if [[ "${LEGACY_OK}" != "1" || "${CONF_OK}" != "1" || "${GC_OK}" != "1" ]]; then
    echo "[REPRODUCE] !! acceptance floor failed; the §3 numbers are not cite-able yet." >&2
    exit 3
fi

echo
echo "================================================================"
echo "[14/14] Paired synthetic-plasma honest baseline"
echo "================================================================"
# The methods paper's section 3.1 table quotes the paired-cohort benchmark
# (20 patients x 2 arms, real TCGA-LUAD panel + real mutation-derived
# fragmentomics). Its artifact was committed from an older run of the
# pipeline and no step regenerated it, so the table was a record rather
# than a reproduction. This step re-measures it with the documented
# command, so the table's numbers are compared like the rest.
#
# It needs the GDC masked-MAF cache (validation/tcga/tcga_cache) and
# torch; without either, the step is skipped and says so, and the
# reproduction check then reports the artifact as NOT regenerated --
# the correct verdict rather than a silent pass.
MAF_COUNT=$(find "${REPO_ROOT}/validation/tcga/tcga_cache" -maxdepth 1 \
    -name 'gdc_TCGA-LUAD_*.maf.gz' 2>/dev/null | wc -l | tr -d ' ')
if [[ "${MAF_COUNT}" == "0" ]]; then
    echo "[REPRODUCE] SKIP: no gdc_TCGA-LUAD_*.maf.gz in"
    echo "[REPRODUCE]       ${REPO_ROOT}/validation/tcga/tcga_cache -- the"
    echo "[REPRODUCE]       paired synthetic-plasma baseline cannot be"
    echo "[REPRODUCE]       regenerated without the GDC cache."
elif ! command -v java >/dev/null 2>&1; then
    echo "[REPRODUCE] SKIP: java is not on PATH (the fragmentomics decoder needs it)"
else
    echo "[REPRODUCE]       ${MAF_COUNT} MAF file(s) in the GDC cache"
    "${PYTHON}" scripts/foundation_real_smoke.py \
        --seeds 20 --loss ce \
        --out "${RESULTS_DIR}/sens_at_spec_ce_n20.json"
fi

# ---------------------------------------------- derived readouts ---------
# Re-emit the ultra-early readiness verdict from the freshly refreshed
# JSON artifacts so docs/ULTRA_EARLY_READINESS_AUTO.md is always in
# lock-step with results/*.json. This is offline (no model weights, no
# network) and reads only cached JSON.
#
# It runs BEFORE the proof step on purpose: the readiness verdict is
# derived from the artifacts the run just regenerated, so it belongs in
# the comparison rather than after it.
echo
echo "================================================================"
echo "[READOUT] Ultra-early readiness verdict"
echo "================================================================"
"${PYTHON}" scripts/ultraearly_readiness.py \
    --results-dir "${RESULTS_DIR}" \
    --out-json "${RESULTS_DIR}/ultraearly_readiness.json" \
    --out-md "${DOCS_DIR}/ULTRA_EARLY_READINESS_AUTO.md" \
    || echo "[REPRODUCE] !! ultra-early readiness CLI failed (non-fatal)"

# ---------------------------------------------------------- the proof -----
# Compare the freshly generated artifacts against the committed ones. This
# is the step that makes the driver a reproduction check rather than a
# re-run: a number that moved fails here.
echo
echo "================================================================"
echo "[PROVE] Regenerated numbers vs the committed artifacts"
echo "================================================================"
"${PYTHON}" scripts/reproduce_check.py \
    --baseline "${BASELINE_DIR}" \
    --current "${RESULTS_DIR}" \
    --tolerance "${REPRO_TOLERANCE}" \
    --threshold-tolerance "${REPRO_THRESHOLD_TOLERANCE}"

# ------------------------------------------------------------ gate ------
if [[ "${NO_TEST}" == "1" ]]; then
    echo
    echo "[REPRODUCE] --no-test set; skipping the pytest gate."
else
    echo
    echo "================================================================"
    echo "[GATE] Pytest gate (cross-study + per-cancer + pretrain + harmony + projection)"
    echo "================================================================"
    # `--timeout` needs the pytest-timeout plugin. It is not part of the
    # declared dependencies, so passing it unconditionally made pytest reject
    # the unknown argument and the gate fail on a clean checkout. Add it
    # only when the plugin is importable.
    GATE_TESTS=(
        test/test_biomedical_review_fixes.py
        test/test_sparse_aware_projection.py
        test/test_finaledb_pretrained_loader.py
        test/test_pretrain_bug_fix.py
        test/test_harmonization_check.py
        test/test_per_cancer_sens_at_spec.py
        src/foundation/test_integration.py
    )
    GATE_ARGS=()
    if "${PYTHON}" -c "import pytest_timeout" 2>/dev/null; then
        GATE_ARGS+=(--timeout=60)
    else
        echo "[REPRODUCE] pytest-timeout not installed; running the gate without --timeout"
    fi
    # NOTE: `${GATE_ARGS[@]+"${GATE_ARGS[@]}"}` rather than plain
    # `"${GATE_ARGS[@]}"`. macOS ships bash 3.2, where expanding an EMPTY
    # array under `set -u` raises "unbound variable" and kills the run.
    "${PYTHON}" -m pytest "${GATE_TESTS[@]}" ${GATE_ARGS[@]+"${GATE_ARGS[@]}"} -q
fi

echo
echo "================================================================"
echo "[REPRODUCE] Done."
echo "================================================================"
echo "  ${RESULTS_DIR}/cross_study_finallydb.json              (legacy, no GC)"
echo "  ${RESULTS_DIR}/cross_study_finallydb_gc_corrected.json (headline)"
echo "  ${RESULTS_DIR}/cross_study_finallydb_with_motifs.json   (with 4-mer motifs)"
echo "  ${RESULTS_DIR}/per_cancer_sens_at_spec.json"
echo "  ${RESULTS_DIR}/cross_study_finallydb_shuffled_control.json"
echo "  ${RESULTS_DIR}/cross_study_loso.json"
echo "  ${RESULTS_DIR}/external_validation_finaledb.json"
echo "  ${RESULTS_DIR}/fusion_ablation_isotonic.json"
echo "  ${RESULTS_DIR}/fusion_regularized_smoke.json"
echo "  ${RESULTS_DIR}/real_tcga_validation.json"
echo "  ${RESULTS_DIR}/decision_curve.json"
echo "  ${RESULTS_DIR}/methylation_cohort_channel.json"
echo "  ${RESULTS_DIR}/ultraearly_readiness.json"
echo "  ${RESULTS_DIR}/sens_at_spec_ce_n20.json"
echo "  ${DOCS_DIR}/CROSS_STUDY_BENCHMARK.md"
echo "  ${DOCS_DIR}/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md"
echo "  ${DOCS_DIR}/CROSS_STUDY_LOSO.md"
echo "  ${DOCS_DIR}/EXTERNAL_VALIDATION.md"
echo "  ${FIG_DIR}/fig1..fig4*.png"
echo "  paper/METHODS_PAPER.md §3 — the paper these artifacts back"
echo "================================================================"
if [[ "${QUICK}" == "1" ]]; then
    echo "[REPRODUCE] quick mode: the repo was NOT modified; scratch output discarded."
fi
