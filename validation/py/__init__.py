"""
DeepCatch Python Validation Pipeline

Core packages for reproducing ALL key claims from FINAL_REAL_DATA_REPORT.md
using Python (numpy/scipy/sklearn) — the data science standard.

Export key classes and functions for easy importing.
"""

from .config import (SEED, SEEDS, N_BOOTSTRAP, N_FOLDS, SEQUENCING_DEPTH,
                     ERROR_RATE, CTDNA_LEVELS, CANCER_TYPES)
from .tcga_loader import load_tcga_data
from .statistical_tests import (bootstrap_ci, delong_test, bonferroni_correct,
                                benjamini_hochberg, compute_auc)
from .realistic_downsample import (downsample_to_cfdna, apply_chip,
                                   apply_variable_shedding, apply_trinucleotide_errors,
                                   apply_variable_input, apply_batch_effects,
                                   apply_inflammatory_elevation)
from .performance_weighted_fusion import (performance_weighted_fusion,
                                          simple_average_fusion, selective_fusion)
from .head_to_head import run_head_to_head
from .cet_validation import run_cet_validation
from .too_validation import run_too_validation
from .compare_published import generate_clinical_comparison

__all__ = [
    # Config
    'SEED', 'SEEDS', 'N_BOOTSTRAP', 'N_FOLDS', 'SEQUENCING_DEPTH',
    'ERROR_RATE', 'CTDNA_LEVELS', 'CANCER_TYPES',
    # Data
    'load_tcga_data',
    # Statistics
    'bootstrap_ci', 'delong_test', 'bonferroni_correct', 'benjamini_hochberg',
    'compute_auc',
    # Downsampling
    'downsample_to_cfdna', 'apply_chip', 'apply_variable_shedding',
    'apply_trinucleotide_errors', 'apply_variable_input', 'apply_batch_effects',
    'apply_inflammatory_elevation',
    # Fusion
    'performance_weighted_fusion', 'simple_average_fusion', 'selective_fusion',
    # Validations
    'run_head_to_head', 'run_cet_validation', 'run_too_validation',
    'generate_clinical_comparison',
]
