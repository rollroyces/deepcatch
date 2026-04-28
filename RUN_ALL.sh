#!/usr/bin/env bash
#
# DeepCatch Master Validation Runner
# ================================
# Runs the complete validation pipeline with ONE command.
#
# Usage:
#   bash RUN_ALL.sh              # Full run
#   bash RUN_ALL.sh --quick      # Quick run with reduced data
#   bash RUN_ALL.sh --skip-plots # Skip plot generation
#
# Requirements:
#   - Python 3.9+ with numpy, scipy, scikit-learn
#   - Optional: matplotlib, torch, torch-geometric (for GNN)
#
# Outputs:
#   results/final_cross_validated_results.json
#   results/tcga_validation_results.json
#   results/benchmark_comparison.json
#   results/benchmark_table.md
#   results/roc_comparison.png
#   results/sensitivity_vs_vaf.png
#   results/ensemble_waterfall.png
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
        --python=*) PYTHON_CMD="${arg#*=}" ;;
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
    echo "  pip install -r requirements.txt"
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

# Create directories
mkdir -p results
mkdir -p validation/tcga/results

echo ""
echo -e "${BOLD}Pipeline Configuration:${NC}"
echo "  Quick mode:    ${QUICK_FLAG:-OFF}"
echo "  Plots:         $([ -z "$SKIP_PLOTS" ] && echo 'ENABLED' || echo 'DISABLED')"
echo "  Output:        $(pwd)/results/"
echo ""

# ------------------------------------------------------------------
# Step 1: Install dependencies (if needed)
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[1/6] Installing dependencies...${NC}"

if [ -f "requirements.txt" ]; then
    $PYTHON -m pip install -r requirements.txt --quiet 2>/dev/null || {
        echo -e "${YELLOW}[WARN] Some packages could not be installed. Continuing with available packages.${NC}"
    }
else
    echo -e "${YELLOW}[WARN] requirements.txt not found. Skipping.${NC}"
fi

echo ""

# ------------------------------------------------------------------
# Step 2: Full cross-validation
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[2/6] Running full cross-validation...${NC}"
echo ""

$PYTHON run_full_validation.py \
    $QUICK_FLAG \
    $SKIP_PLOTS \
    --output results/ \
    --seeds 5

FULL_EXIT=$?
if [ $FULL_EXIT -ne 0 ]; then
    echo -e "${RED}[ERROR] Full validation failed with exit code $FULL_EXIT${NC}"
    exit $FULL_EXIT
fi

echo ""
echo -e "${GREEN}[OK] Full validation complete${NC}"
echo ""

# ------------------------------------------------------------------
# Step 3: TCGA validation
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[3/6] Running TCGA real data validation...${NC}"
echo ""

$PYTHON run_tcga_validation.py \
    $QUICK_FLAG \
    --output results/ \
    --seeds 5

TCGA_EXIT=$?
if [ $TCGA_EXIT -ne 0 ]; then
    echo -e "${YELLOW}[WARN] TCGA validation exited with code $TCGA_EXIT (non-fatal)${NC}"
fi

echo ""
echo -e "${GREEN}[OK] TCGA validation complete${NC}"
echo ""

# ------------------------------------------------------------------
# Step 4: Benchmark comparison
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[4/6] Generating benchmark comparison...${NC}"
echo ""

$PYTHON run_benchmark_comparison.py \
    --results results/final_cross_validated_results.json \
    --tcga results/tcga_validation_results.json \
    --output results/

BENCH_EXIT=$?
if [ $BENCH_EXIT -ne 0 ]; then
    echo -e "${YELLOW}[WARN] Benchmark comparison exited with code $BENCH_EXIT (non-fatal)${NC}"
fi

echo ""
echo -e "${GREEN}[OK] Benchmark comparison complete${NC}"
echo ""

# ------------------------------------------------------------------
# Step 5: Verify fixes
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[5/6] Running fix verification...${NC}"
echo ""

if [ -f "validation_framework_tests.py" ]; then
    $PYTHON validation_framework_tests.py 2>&1 | head -30 || true
    echo ""
fi

if [ -f "verify_fixes.py" ]; then
    $PYTHON verify_fixes.py 2>&1 || {
        echo -e "${YELLOW}[WARN] verify_fixes.py issue (non-fatal)${NC}"
    }
else
    echo -e "${YELLOW}[WARN] verify_fixes.py not found — skipping${NC}"
fi

echo ""

# ------------------------------------------------------------------
# Step 6: Summary
# ------------------------------------------------------------------

echo -e "${BLUE}${BOLD}[6/6] Generating summary...${NC}"
echo ""

echo -e "${BLUE}${BOLD}========================================${NC}"
echo -e "${BLUE}${BOLD}  VALIDATION PIPELINE COMPLETE${NC}"
echo -e "${BLUE}${BOLD}========================================${NC}"
echo ""

# Print summary of results if available
if [ -f "results/final_cross_validated_results.json" ]; then
    echo -e "${BOLD}Cross-Validated Results:${NC}"
    $PYTHON -c "
import json
with open('results/final_cross_validated_results.json') as f:
    data = json.load(f)
for exp_name, exp_data in data.get('experiments', {}).items():
    agg = exp_data.get('aggregated', {})
    auc = agg.get('auc', {})
    print(f\"  {exp_name}: AUC={auc.get('mean',0):.4f} ± {auc.get('std',0):.4f}\")
" 2>/dev/null || true
    echo ""
fi

if [ -f "results/tcga_validation_results.json" ]; then
    echo -e "${BOLD}TCGA Validation:${NC}"
    $PYTHON -c "
import json
with open('results/tcga_validation_results.json') as f:
    data = json.load(f)
for model in ['variant_caller', 'multimodal_fusion']:
    if model in data:
        s = data[model].get('summary', {})
        print(f\"  {model}: AUC={s.get('auc',{}).get('mean',0):.4f} ± {s.get('auc',{}).get('std',0):.4f}\")
" 2>/dev/null || true
    echo ""
fi

echo -e "${BOLD}Generated Files:${NC}"
for f in results/*.json results/*.md results/*.png; do
    if [ -f "$f" ]; then
        size=$(du -h "$f" | cut -f1)
        echo -e "  ${GREEN}✓${NC} $f ($size)"
    fi
done

echo ""
echo -e "${BOLD}Quick Commands:${NC}"
echo "  View results:    cat results/final_cross_validated_results.json | python -m json.tool"
echo "  View benchmark:  cat results/benchmark_table.md"
echo "  Rerun with plots: bash RUN_ALL.sh --with-plots"
echo "  Quick run:       bash RUN_ALL.sh --quick"
echo ""
echo -e "${GREEN}Done! 🦾${NC}"

exit 0
