"""
Central configuration for the DeepCatch Python validation pipeline.

All constants, paths, and hyperparameters live here so every module
can import them consistently.  Matches validation/node/realisticDownsample.js
and the FINAL_REAL_DATA_REPORT.md parameter table exactly.
"""

import os
from pathlib import Path

# ── Reproducibility ──────────────────────────────────────────────────────
SEED = 42
SEEDS = [42, 123, 456, 789, 1024]

# ── CV / Bootstrap ───────────────────────────────────────────────────────
N_BOOTSTRAP = 2000
N_FOLDS = 5

# ── Sequencing Parameters (match realisticDownsample.js) ─────────────────
SEQUENCING_DEPTH = 50000         # baseline sequencing depth ×
ERROR_RATE = 0.0001              # baseline per-base error rate
N_BACKGROUND_SITES = 5000        # background sites for specificity estimation
N_LOCI = 50                      # multi-locus monitoring for CET

# ── ctDNA fractions tested ───────────────────────────────────────────────
CTDNA_LEVELS = [0.01, 0.005, 0.0025, 0.001, 0.0005,
                0.00025, 0.0001, 0.00005, 0.00001]

# ── Cancer Types (sourced from COSMIC v99 + TCGA PanCancer) ──────────────
CANCER_TYPES = ['LUAD', 'COADREAD', 'BRCA', 'PRAD',
                'STAD', 'LIHC', 'PAAD', 'OV', 'BLCA', 'HNSC']

# ── Paths ────────────────────────────────────────────────────────────────
ROOT_DIR = Path(__file__).resolve().parent.parent.parent  # cancer-screening/
RESULTS_DIR = ROOT_DIR / 'results'
RESULTS_NODE_DIR = RESULTS_DIR / 'node'
RESULTS_PY_DIR = RESULTS_DIR / 'py'

# Ensure results directory exists
RESULTS_PY_DIR.mkdir(parents=True, exist_ok=True)

# ── Input data paths ─────────────────────────────────────────────────────
TCGA_INPUT_PATH = RESULTS_NODE_DIR / 'real_tcga_data.json'
DOWNSAMPLED_PATH = RESULTS_NODE_DIR / 'real_downsampled.json'
H2H_RESULTS_PATH = RESULTS_NODE_DIR / 'real_headToHead_results.json'
CET_RESULTS_PATH = RESULTS_NODE_DIR / 'real_cet_results.json'

# ── Output paths ─────────────────────────────────────────────────────────
PY_DOWNSAMPLED_PATH = RESULTS_PY_DIR / 'real_downsampled.json'
PY_H2H_PATH = RESULTS_PY_DIR / 'head_to_head_results.json'
PY_CET_PATH = RESULTS_PY_DIR / 'cet_results.json'
PY_TOO_PATH = RESULTS_PY_DIR / 'too_results.json'
PY_COMPARISON_PATH = RESULTS_PY_DIR / 'clinical_comparison.json'


if __name__ == '__main__':
    print("DeepCatch configuration:")
    print(f"  SEED={SEED}, N_BOOTSTRAP={N_BOOTSTRAP}, N_FOLDS={N_FOLDS}")
    print(f"  SEQUENCING_DEPTH={SEQUENCING_DEPTH}×, ERROR_RATE={ERROR_RATE}")
    print(f"  CTDNA_LEVELS: {CTDNA_LEVELS}")
    print(f"  CANCER_TYPES: {CANCER_TYPES}")
    print(f"  Results dir: {RESULTS_PY_DIR}")
