#!/usr/bin/env bash
#
# DeepCatch Master Validation Runner
# ================================
# Runs the complete validation pipeline with ONE command.
#
# Usage:
#   bash RUN_ALL.sh                 # Full run
#   bash RUN_ALL.sh --quick         # Quick run with reduced data
#   bash RUN_ALL.sh --with-plots    # Also render the comparison plots
#   bash RUN_ALL.sh --python=PATH   # Use a specific interpreter
#
# Environment:
#   RUN_ALL_OUTPUT=<dir>   Write the run's outputs to <dir> instead of a
#                          scratch directory
#   JIANG_INPUT=<file>     Run step 4 on your own copy of the Jiang 4-mer
#                          dataset (not shipped with this repository)
#
# Requirements:
#   - Python 3.10+ with numpy, scipy, scikit-learn
#   - Optional: matplotlib (plots), torch + torch-geometric (GNN channel)
#
# Outputs -- written to the output directory, which is a scratch
# directory by default so a run never replaces the repository's
# published artifacts under results/; set RUN_ALL_OUTPUT to choose one:
#   final_cross_validated_results.json   cross-validated AUCs
#   roc_comparison.png                   (--with-plots)
#   sensitivity_vs_vaf.png               (--with-plots)
#   ensemble_waterfall.png               (--with-plots)
#   real_tcga_validation.json            step 3, when the TCGA cache is present
#   real_tcga_performance.png            step 3, with --with-plots
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Parse flags
QUICK_FLAG=""
SKIP_PLOTS="--skip-plots"
PYTHON_CMD=""

for arg in "$@"; do
    case "$arg" in
        --quick) QUICK_FLAG="--quick" ;;
        --with-plots) SKIP_PLOTS="" ;;
        --skip-plots) SKIP_PLOTS="--skip-plots" ;;   # plots are off by default
        --python=*) PYTHON_CMD="${arg#*=}" ;;
        *)
            echo -e "${RED}[ERROR] Unknown option: $arg${NC}"
            echo "Usage: bash RUN_ALL.sh [--quick] [--with-plots] [--python=PATH]"
            exit 2
            ;;
    esac
done

# ------------------------------------------------------------------
# Environment setup
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}========================================${NC}"
echo -e "${BLUE}${BOLD}  DeepCatch Master Validation Pipeline${NC}"
echo -e "${BLUE}${BOLD}========================================${NC}"
echo ""

# Detect Python
if [ -n "$PYTHON_CMD" ]; then
    PYTHON="$PYTHON_CMD"
else
    # Try common Python commands
    for cmd in python3 python python3.11 python3.10; do
        if command -v "$cmd" &> /dev/null; then
            PYTHON="$cmd"
            break
        fi
    done
fi

if [ -z "${PYTHON:-}" ]; then
    echo -e "${RED}[ERROR] Python not found. Please install Python 3.9+.${NC}"
    echo "  Or specify: bash RUN_ALL.sh --python=/path/to/python"
    exit 1
fi

echo -e "${GREEN}[INFO] Using Python: $PYTHON ($($PYTHON --version 2>&1))${NC}"

# Check numpy
echo -n "[INFO] Checking numpy... "
if $PYTHON -c "import numpy" 2>/dev/null; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${RED}MISSING${NC}"
    echo ""
    echo "Required packages not found. Install them first:"
    echo "  pip install -e ."
    echo ""
    exit 1
fi

# Check scipy
echo -n "[INFO] Checking scipy... "
if $PYTHON -c "import scipy" 2>/dev/null; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${RED}MISSING${NC}"
    exit 1
fi

# Check sklearn
echo -n "[INFO] Checking sklearn... "
if $PYTHON -c "import sklearn" 2>/dev/null; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${RED}MISSING${NC}"
    exit 1
fi

# Check matplotlib (optional)
echo -n "[INFO] Checking matplotlib... "
if $PYTHON -c "import matplotlib" 2>/dev/null; then
    echo -e "${GREEN}OK${NC}"
else
    echo -e "${YELLOW}MISSING (plots disabled)${NC}"
    SKIP_PLOTS="--skip-plots"
fi

# Where this run writes. The artifacts under results/ are the
# repository's published ones -- a validation run must not replace them,
# so the output goes to a scratch directory unless the caller names one.
OUTPUT_DIR="${RUN_ALL_OUTPUT:-$(mktemp -d -t deepcatch-run-all-XXXXXX)}"
mkdir -p "$OUTPUT_DIR"
mkdir -p validation/tcga/results

echo ""
echo -e "${BOLD}Pipeline Configuration:${NC}"
echo "  Quick mode:    ${QUICK_FLAG:-OFF}"
echo "  Plots:         $([ -z "$SKIP_PLOTS" ] && echo 'ENABLED' || echo 'DISABLED')"
echo "  Output:        ${OUTPUT_DIR}"
echo ""

# ------------------------------------------------------------------
# Step 1: Install dependencies (if needed)
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[1/5] Installing dependencies...${NC}"

$PYTHON -m pip install -e . --quiet 2>/dev/null || {
    echo -e "${YELLOW}[WARN] The package could not be installed (no network?). Continuing with the packages already present.${NC}"
}

echo ""

# ------------------------------------------------------------------
# Step 2: Full cross-validation
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[2/5] Running full cross-validation...${NC}"
echo ""

FULL_EXIT=0
$PYTHON run_full_validation.py \
    $QUICK_FLAG \
    $SKIP_PLOTS \
    --output "${OUTPUT_DIR}" \
    --seeds 5 || FULL_EXIT=$?
if [ $FULL_EXIT -ne 0 ]; then
    echo -e "${RED}[ERROR] Full validation failed with exit code $FULL_EXIT${NC}"
    exit $FULL_EXIT
fi

echo ""
echo -e "${GREEN}[OK] Full validation complete${NC}"
echo ""

# ------------------------------------------------------------------
# Step 3: Real data validation (TCGA) — optional
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[3/5] Running real data validation...${NC}"
echo ""

TCGA_SCRIPT="$(dirname $0)/real_tcga_validation.py"
if [ -f "$TCGA_SCRIPT" ]; then
    # real_tcga_validation.py has no --quick flag of its own; the reduced
    # scale is expressed as a smaller cohort instead.
    TCGA_QUICK_ARGS=""
    if [ -n "$QUICK_FLAG" ]; then
        TCGA_QUICK_ARGS="--n-patients 5"
    fi
    $PYTHON real_tcga_validation.py \
        $TCGA_QUICK_ARGS \
        --output "${OUTPUT_DIR}/real_tcga_validation.json" \
        --seeds 5 || echo -e "${YELLOW}[WARN} TCGA validation non-fatal error${NC}"
else
    echo -e "${YELLOW}[WARN] real_tcga_validation.py not found — skipping${NC}"
fi

echo ""
echo -e "${GREEN}[OK] TCGA validation step complete${NC}"
echo ""

# ------------------------------------------------------------------
# Step 4: Jiang 4-mer analysis — optional
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[4/5] Running Jiang 4-mer analysis...${NC}"
echo ""

JIANG_SCRIPT="$(dirname $0)/run_jiang_analysis.py"
if [ -n "${JIANG_INPUT:-}" ] && [ -f "$JIANG_INPUT" ]; then
    # The Jiang 4-mer dataset is not part of this repository (it was
    # requested from the source lab -- see data/JIANG_DATA_REQUEST.md),
    # so this step only runs against a copy the user provides.
    $PYTHON run_jiang_analysis.py \
        --input "$JIANG_INPUT" \
        ${JIANG_ARGS:-} \
        --output "${OUTPUT_DIR}/jiang_reanalysis/" || echo -e "${YELLOW}[WARN] Jiang analysis non-fatal error${NC}"
elif [ -f "$JIANG_SCRIPT" ]; then
    echo -e "${YELLOW}[SKIP] No Jiang 4-mer dataset in this repository (see data/JIANG_DATA_REQUEST.md).${NC}"
    echo -e "${YELLOW}       Set JIANG_INPUT=<path> to run this step on your own copy.${NC}"
else
    echo -e "${YELLOW}[WARN] run_jiang_analysis.py not found — skipping${NC}"
fi

echo ""
echo -e "${GREEN}[OK] Jiang analysis step complete${NC}"
echo ""

# ------------------------------------------------------------------
# Step 5: Supplementary checks
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[5/5] Running supplementary checks...${NC}"
echo ""

# Verify Python syntax across all .py files
BROKEN=0
for f in $(find . -name '*.py' -not -path './.git/*' | sort); do
    $PYTHON -c "import ast; ast.parse(open('$f').read())" 2>/dev/null || {
        echo -e "  ${RED}⚠ Syntax error in $f${NC}"
        BROKEN=$((BROKEN+1))
    }
done
if [ $BROKEN -eq 0 ]; then
    echo -e "  ${GREEN}✓ All Python files pass syntax check${NC}"
else
    echo -e "  ${YELLOW}⚠ $BROKEN file(s) have syntax errors${NC}"
fi

echo ""

# ------------------------------------------------------------------
# Step 6: Summary
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[Summary]${NC}"
echo ""

echo -e "${BLUE}${BOLD}========================================${NC}"
echo -e "${BLUE}${BOLD}  VALIDATION PIPELINE COMPLETE${NC}"
echo -e "${BLUE}${BOLD}========================================${NC}"
echo ""

# Print summary of results if available
if [ -f "${OUTPUT_DIR}/final_cross_validated_results.json" ]; then
    echo -e "${BOLD}Cross-Validated Results:${NC}"
    $PYTHON -c "
import json
with open('${OUTPUT_DIR}/final_cross_validated_results.json') as f:
    data = json.load(f)
for exp_name, exp_data in data.get('experiments', {}).items():
    agg = exp_data.get('aggregated', {})
    auc = agg.get('auc', {})
    print(f\"  {exp_name}: AUC={auc.get('mean',0):.4f} ± {auc.get('std',0):.4f}\")
" 2>/dev/null || true
    echo ""
fi

if [ -f "${OUTPUT_DIR}/real_tcga_validation.json" ]; then
    echo -e "${BOLD}TCGA Validation:${NC}"
    $PYTHON -c "
import json
with open('${OUTPUT_DIR}/real_tcga_validation.json') as f:
    data = json.load(f)
for model in ['variant_caller', 'ml_classifier']:
    entries = data.get('results', {}).get(model) or []
    auc = next((e for e in entries if e.get('metric') == 'auc'), None)
    if auc:
        print(f\"  {model}: AUC={auc['mean']:.4f} ± {auc['std']:.4f}\")
" 2>/dev/null || true
    echo ""
fi

echo -e "${BOLD}Generated Files:${NC}"
for f in "${OUTPUT_DIR}"/*.json "${OUTPUT_DIR}"/*.md "${OUTPUT_DIR}"/*.png; do
    if [ -f "$f" ]; then
        size=$(du -h "$f" | cut -f1)
        echo -e "  ${GREEN}✓${NC} $f ($size)"
    fi
done

echo ""
echo -e "${BOLD}Quick Commands:${NC}"
echo "  Output directory: ${OUTPUT_DIR}"
echo "  View results:    cat ${OUTPUT_DIR}/final_cross_validated_results.json | python -m json.tool"
echo "  Rerun with plots: bash RUN_ALL.sh --with-plots"
echo "  Quick run:       bash RUN_ALL.sh --quick"
echo ""
echo -e "${GREEN}Done! 🦾${NC}"

exit 0
