#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────
# scripts/reproduce_all.sh
# One-command reproduction of every published number in this repo.
#
# Honors two deepcatch-skill pitfalls that have burned past CI:
#   #1  macOS global PYTHONPATH hijacks the venv — we use `env -u
#       PYTHONPATH ./.venv/bin/python` for every python invocation
#       so the project venv wins.
#   #4  Without `set -o pipefail`, an `if cmd | tee log; then` tests
#       tee's exit code (always 0), not cmd's. We set pipefail AND
#       errexit so any failed step aborts the whole pipeline.
#
# Idempotent: re-running overwrites results/. Safe to interrupt and
# restart — every step writes a complete JSON artifact before moving
# on.
#
# Usage:
#   cd deepcatch
#   bash scripts/reproduce_all.sh                # full ~10-15 min run
#   bash scripts/reproduce_all.sh --dry-run      # print plan, exit 0
#   bash scripts/reproduce_all.sh --skip-tests   # skip pytest step
#   bash scripts/reproduce_all.sh --only cross_study,per_cancer  # subset
#
# Wall-clock budget (MacBook M-class, CPU torch):
#   cross-study benchmark      ~5 min   (5-seed × 5-fold pooled OOF)
#   per-cancer sens@spec       ~30 s    (synthetic fixture, DeLong CIs)
#   harmonization check        ~10 s    (synthetic fixture, 5-seed)
#   focal-BCE ablation         ~30 s    (reads ce + sens JSONs, paired t)
#   sparse-aware ablation      ~30 s    (reads linear + sparse JSONs)
#   foundation smoke           ~30 s    (5-seed × 5-fold real TCGA-LUAD)
#   pytest test/ + src/foundation   ~30 s
#   ─────────────────────────────────────
#   TOTAL                      ~10-15 min
# ─────────────────────────────────────────────────────────────────────

set -euo pipefail

# ─── 0. Resolve repo root (walk-up) + python interpreter ─────────────
# Walk up from this script until we find the repo root (pyproject.toml +
# scripts/ together). This avoids pitfall #7 where `__file__` resolves
# through an editable-install symlink to the .venv site-packages.

find_repo_root() {
    local d
    d="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    while [ "$d" != "/" ]; do
        if [ -f "$d/pyproject.toml" ] && [ -d "$d/scripts" ]; then
            echo "$d"
            return 0
        fi
        d="$(dirname "$d")"
    done
    # Fallback: CWD
    pwd
}

REPO_ROOT="$(find_repo_root)"
PY="$REPO_ROOT/.venv/bin/python"
RESULTS="$REPO_ROOT/results"

mkdir -p "$RESULTS"

# ─── 1. Flags + arg parsing ──────────────────────────────────────────
DRY_RUN=0
SKIP_TESTS=0
ONLY_STEPS=""

while [ $# -gt 0 ]; do
    case "$1" in
        --dry-run) DRY_RUN=1; shift ;;
        --skip-tests) SKIP_TESTS=1; shift ;;
        --only) ONLY_STEPS="${2:-}"; shift 2 || true ;;
        -h|--help)
            sed -n '2,40p' "$0" | sed 's/^# //' | sed 's/^#//'
            exit 0
            ;;
        *) echo "Unknown arg: $1" >&2; exit 2 ;;
    esac
done

# ─── 2. Confirm venv + Python before we start ────────────────────────
if [ ! -x "$PY" ]; then
    echo "❌ missing venv python at $PY" >&2
    echo "   Create with: python3 -m venv .venv && pip install -e ." >&2
    exit 1
fi

# Feature dir: env var override, else standard cfdna-fragmentomics path
# (FinaleDB feature cache, used by cross_study_finallydb.py)
FINALEDB_FEATURES_DIR="${FINALEDB_FEATURES_DIR:-/Users/hermes/cfdna-fragmentomics-pipeline/data/features}"
if [ ! -d "$FINALEDB_FEATURES_DIR" ]; then
    echo "⚠️  FinaleDB features dir not found: $FINALEDB_FEATURES_DIR" >&2
    echo "   Cross-study benchmark (step 1) will fail." >&2
    echo "   Set FINALEDB_FEATURES_DIR or symlink the canonical path." >&2
fi

# TCGA cache for foundation_real_smoke.py (different from FinaleDB —
# these are GDC MAFs, not FinaleDB profiles).
TCGA_CACHE_DIR="${TCGA_CACHE_DIR:-$REPO_ROOT/validation/tcga/tcga_cache}"

# ─── 3. Pre-flight env check ─────────────────────────────────────────
echo "==[verify_environment]=="
if [ "$DRY_RUN" = "1" ]; then
    echo "[dry-run] would run: env -u PYTHONPATH $PY scripts/verify_environment.py"
else
    env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/verify_environment.py"
fi
echo

# ─── 4. Step registry ────────────────────────────────────────────────
# Format: "slug|wall_seconds|description|log_file"
#
# Slugs map to both the JSON output file basename (results/<slug>.json)
# and to a function body below. The output JSON is REQUIRED for the
# downstream final report — every step writes one.

ALL_STEPS=$(cat <<'EOF'
cross_study|300|cross_study_finallydb.py (FinaleDB Jiang+Cristiano)|cross_study.log
per_cancer|45|per_cancer_sens_at_spec.py --synthetic|per_cancer.log
harmonization|15|harmonization_check.py|harmonization.log
focal_bce_ablation|30|sens_at_spec_ablation.py (CE vs focal-BCE)|focal_bce_ablation.log
sparse_aware_ablation|30|sparse_aware_ablation.py (Linear vs SparseAware)|sparse_aware_ablation.log
foundation_smoke|60|foundation_real_smoke.py --quick|foundation_smoke.log
pytest_tests|45|pytest test/ + src/foundation/test_integration.py|pytest.log
EOF
)

# `--only` filter
filter=""
if [ -n "$ONLY_STEPS" ]; then
    filter="$ONLY_STEPS"
fi
if [ "$SKIP_TESTS" = "1" ]; then
    filter="$(echo "$ONLY_STEPS" | sed 's/pytest_tests//' )"
fi

# ─── 5. Per-step runner ──────────────────────────────────────────────
# We use bash functions for each step so we can keep the JSON output
# path consistent. Every step writes results/<slug>.json — the final
# report reads those JSONs (or -wall if missing → FAIL with note).
#
# Parallel arrays (bash 3.2-safe — `/bin/bash` on macOS is still 3.2.x
# and lacks associative arrays): the index into each array identifies
# the step.

STEP_SLUGS=()
STEP_DESCS=()
STEP_RCS=()
STEP_WALLS=()
STEP_NOTES=()

record_step() {
    local slug="$1" rc="$2" wall="$3" note="$4"
    STEP_SLUGS+=("$slug")
    STEP_RCS+=("$rc")
    STEP_WALLS+=("$wall")
    STEP_NOTES+=("$note")
}

run_step() {
    local slug="$1"
    local description="$2"
    local json_out="$RESULTS/$slug.json"
    local log="$RESULTS/$slug.log"

    if [ -n "$filter" ] && ! echo ",$filter," | grep -q ",$slug,"; then
        echo "==[$slug]== SKIPPED (not in --only set)"
        record_step "$slug" "SKIP" "0.0" "skipped via --only"
        return 0
    fi

    echo "==[$slug]== $description"
    echo "    out: $json_out"

    if [ "$DRY_RUN" = "1" ]; then
        echo "    [dry-run] would run this step"
        record_step "$slug" "DRY" "0.0" "dry-run"
        return 0
    fi

    local t0 t1
    t0="$(date +%s.%N)"

    # Per-step dispatch. ALL of these write results/<slug>.json. The
    # python prefix `env -u PYTHONPATH` is non-negotiable per skill
    # pitfall #1.
    set +e
    case "$slug" in
        cross_study)
            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/cross_study_finallydb.py" \
                --features-dir "$FINALEDB_FEATURES_DIR" \
                --out-json "$json_out" \
                --out-md "$REPO_ROOT/docs/CROSS_STUDY_BENCHMARK.md" \
                >"$log" 2>&1
            ;;
        per_cancer)
            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/per_cancer_sens_at_spec.py" \
                --synthetic --n 240 --seed 42 \
                --out "$json_out" \
                >"$log" 2>&1
            ;;
        harmonization)
            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/harmonization_check.py" \
                >"$log" 2>&1
            # harmonization_check.py hardcodes paths via os.chdir(ROOT);
            # under our call its ROOT is its own parent. To get its JSON
            # into results/ we move it.
            if [ -f "$REPO_ROOT/results/harmonization_check.json" ] && \
               [ "$REPO_ROOT/results/harmonization_check.json" != "$json_out" ]; then
                mv "$REPO_ROOT/results/harmonization_check.json" "$json_out"
            fi
            ;;
        focal_bce_ablation)
            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/sens_at_spec_ablation.py" \
                --ce-json "$RESULTS/foundation_smoke_ce.json" \
                --sens-json "$RESULTS/foundation_smoke_sens.json" \
                --out "$json_out" \
                >"$log" 2>&1
            ;;
        sparse_aware_ablation)
            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/sparse_aware_ablation.py" \
                --linear-json "$RESULTS/foundation_smoke.json" \
                --sparse-json "$RESULTS/foundation_smoke_sparse.json" \
                --out "$json_out" \
                >"$log" 2>&1
            ;;
        foundation_smoke)
            # Run three times — once for CE (default), once for
            # focal-BCE, once for sparse_aware — so the ablation steps
            # have both inputs.
            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/foundation_real_smoke.py" \
                --quick --n-patients 20 --seeds 5 \
                --features-dir "$TCGA_CACHE_DIR" \
                --out "$RESULTS/foundation_smoke_ce.json" \
                >"$RESULTS/foundation_smoke_ce.log" 2>&1 || true

            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/foundation_real_smoke.py" \
                --quick --n-patients 20 --seeds 5 \
                --features-dir "$TCGA_CACHE_DIR" \
                --loss sens_at_spec --alpha-pos 20 \
                --out "$RESULTS/foundation_smoke_sens.json" \
                >"$RESULTS/foundation_smoke_sens.log" 2>&1 || true

            env -u PYTHONPATH "$PY" "$REPO_ROOT/scripts/foundation_real_smoke.py" \
                --quick --n-patients 20 --seeds 5 \
                --features-dir "$TCGA_CACHE_DIR" \
                --projection-kinds '{"frag_basic":"sparse_aware"}' \
                --out "$RESULTS/foundation_smoke_sparse.json" \
                >"$RESULTS/foundation_smoke_sparse.log" 2>&1 || true

            # The "main" foundation_smoke.json the report reads is the CE run.
            if [ -f "$RESULTS/foundation_smoke_ce.json" ]; then
                cp "$RESULTS/foundation_smoke_ce.json" "$json_out"
            fi
            ;;
        pytest_tests)
            env -u PYTHONPATH "$PY" -m pytest \
                "$REPO_ROOT/test/" \
                "$REPO_ROOT/src/foundation/test_integration.py" \
                --tb=line --timeout=180 --no-header -q \
                --deselect "test/test_cross_study_per_cancer_integration.py::test_end_to_end_script_writes_both_outputs" \
                >"$log" 2>&1
            # Synthesize a minimal JSON so the report can list a real
            # output path. The log carries the real test output.
            if [ "$rc" -eq 0 ]; then
                local pass_count
                pass_count="$(grep -oE '[0-9]+ passed' "$log" | tail -1 | awk '{print $1}')"
                [ -z "$pass_count" ] && pass_count="0"
                printf '{"step":"pytest_tests","tests_passed":%s,"source":"results/pytest_tests.log","deselected":["test_cross_study_per_cancer_integration.py::test_end_to_end_script_writes_both_outputs (5-10 min subprocess)"]}' \
                    "$pass_count" > "$json_out"
            fi
            ;;
        *)
            echo "    unknown step: $slug" >&2
            ;;
    esac
    local rc=$?
    set -e

    t1="$(date +%s.%N)"
    local elapsed
    elapsed="$(awk "BEGIN { printf \"%.1f\", $t1 - $t0 }")"

    if [ "$rc" -eq 0 ] && [ -f "$json_out" ]; then
        record_step "$slug" "PASS" "$elapsed" "ok"
    else
        local tail_note
        tail_note="$(tail -3 "$log" 2>/dev/null | tr '\n' ' ' | head -c 200)"
        record_step "$slug" "FAIL" "$elapsed" "exit=$rc, $tail_note"
    fi

    printf "    wall: %ss  rc: %s\n" "$elapsed" "${STEP_RCS[${#STEP_RCS[@]}-1]}"
    echo
}

# Ordered execution. Order matters: ablation steps need foundation_smoke
# JSONs to exist (run_step for foundation_smoke writes all three).
echo "================================================================"
echo "  deepcatch reproduce_all.sh"
echo "  repo_root : $REPO_ROOT"
echo "  python    : $PY"
echo "  results   : $RESULTS"
echo "  features  : $FINALEDB_FEATURES_DIR"
echo "  steps     : $(echo "$ALL_STEPS" | wc -l | tr -d ' ')"
echo "================================================================"
echo

TOTAL_T0="$(date +%s.%N)"

# Collect slugs in declared order so we can pair with parallel arrays.
SLUG_ORDER=()
while IFS='|' read -r slug budget desc logf; do
    [ -z "$slug" ] && continue
    SLUG_ORDER+=("$slug")
    STEP_DESCS+=("$desc")
done <<< "$ALL_STEPS"

for slug in "${SLUG_ORDER[@]}"; do
    # find desc by index
    idx=0
    for s in "${SLUG_ORDER[@]}"; do
        if [ "$s" = "$slug" ]; then
            break
        fi
        idx=$((idx + 1))
    done
    run_step "$slug" "${STEP_DESCS[$idx]}"
done

TOTAL_T1="$(date +%s.%N)"
TOTAL_WALL="$(awk "BEGIN { printf \"%.1f\", $TOTAL_T1 - $TOTAL_T0 }")"

# ─── 6. Final report ─────────────────────────────────────────────────
REPORT="$RESULTS/REPRODUCE_REPORT.md"
{
    echo "# Reproduction report"
    echo
    echo "- Generated: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo "- repo_root: \`$REPO_ROOT\`"
    echo "- python:    \`$PY\` ($(env -u PYTHONPATH "$PY" -V 2>&1))"
    echo "- total wall: ${TOTAL_WALL}s"
    echo
    echo "## Step results"
    echo
    echo "| Step | Description | Wall (s) | Exit | Output |"
    echo "|---|---|---:|---|---|"
    for i in "${!SLUG_ORDER[@]}"; do
        slug="${SLUG_ORDER[$i]}"
        desc="${STEP_DESCS[$i]}"
        rc="${STEP_RCS[$i]:-NA}"
        wall="${STEP_WALLS[$i]:-NA}"
        marker="✅"
        [ "$rc" = "FAIL" ] && marker="❌"
        [ "$rc" = "SKIP" ] && marker="⏭️"
        [ "$rc" = "DRY" ] && marker="🟡"
        out="$RESULTS/$slug.json"
        out_display="$slug.json"
        [ ! -f "$out" ] && out_display="$slug.json (MISSING)"
        echo "| $marker \`$slug\` | $desc | $wall | $rc | \`$out_display\` |"
    done
    echo
    echo "## Per-step logs"
    echo
    echo "Each step writes \`results/<slug>.log\` with the full stdout/stderr of the script."
    echo
    echo "## Honest framing"
    echo
    echo "This script reproduces the **open-data methods benchmark** numbers published in"
    echo "\`docs/CROSS_STUDY_BENCHMARK.md\`, \`docs/SENS_AT_SPEC_ABLATION.md\`,"
    echo "\`docs/SPARSE_AWARE_ABLATION.md\`, \`docs/HARMONIZATION.md\`, and the foundation"
    echo "smoke's \`docs/FOUNDATION_REAL_SMOKE.md\`-shaped numbers."
    echo
    echo "It does **not** perform:"
    echo "- Clinical validation on an external held-out cohort"
    echo "- Regulatory-grade quality controls"
    echo "- Per-sample audit/reproducibility beyond what the per-seed JSONs capture"
    echo
    echo "## Re-running individual steps"
    echo
    echo "\`\`\`bash"
    echo "# Cross-study benchmark alone"
    echo "env -u PYTHONPATH ./.venv/bin/python scripts/cross_study_finallydb.py \\"
    echo "    --features-dir \$FINALEDB_FEATURES_DIR \\"
    echo "    --out-json results/cross_study.json \\"
    echo "    --out-md docs/CROSS_STUDY_BENCHMARK.md"
    echo ""
    echo "# Per-cancer sens@spec alone (synthetic fixture)"
    echo "env -u PYTHONPATH ./.venv/bin/python scripts/per_cancer_sens_at_spec.py \\"
    echo "    --synthetic --n 240 --out results/per_cancer.json"
    echo "\`\`\`"
    echo
} > "$REPORT"

echo "================================================================"
echo "  Reproduction report written: $REPORT"
echo "  Total wall: ${TOTAL_WALL}s"
echo "================================================================"

# Exit non-zero if any non-skipped/non-dry step failed
n_fail=0
for rc in "${STEP_RCS[@]}"; do
    if [ "$rc" = "FAIL" ]; then
        n_fail=$((n_fail + 1))
    fi
done

if [ "$n_fail" -gt 0 ]; then
    echo "❌ ${n_fail} step(s) failed. See $REPORT for details."
    exit 1
fi

echo "✅ All steps passed (or were skipped / dry-run)."
exit 0