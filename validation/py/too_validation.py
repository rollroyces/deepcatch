"""
Tissue-of-Origin (TOO) Prediction Validation

Mirrors validation/node/realCET.js and the FINAL_REAL_DATA_REPORT.

Simulates multi-class logistic regression on methylation + fragmentomic
features to predict tissue of origin.

⚠️ HONEST NOTE: 100% accuracy on simulation is MEANINGLESS.
Grail clinical TOO = 88.7% across 50+ cancer types.
This module is included for completeness but reports honestly
that simulation TOO is not comparable to clinical TOO.

Features modelled: 9 methylation + 3 fragmentomic = 12 features
Cancer types: LUAD, COADREAD, BRCA, PRAD, STAD, LIHC, PAAD, OV (8)

References:
  Jamshidi 2022 Cancer Cell (Grail TOO 88.7%)
  Cohen 2018 Science (CancerSEEK TOO 83%)
  Cristiano 2019 Nature (DELFI fragmentomics)
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import accuracy_score, confusion_matrix

from .config import SEED, CANCER_TYPES, PY_TOO_PATH

logger = logging.getLogger(__name__)

# Feature configurations
N_METHYLATION_FEATURES = 9
N_FRAGMENTOMIC_FEATURES = 3
N_TOTAL_FEATURES = N_METHYLATION_FEATURES + N_FRAGMENTOMIC_FEATURES
DEFAULT_N_CANCER_TYPES = 8
DEFAULT_N_PER_TYPE = 50


# ── Simulated Feature Generation ─────────────────────────────────────────
def _generate_too_features(
    n_cancer_types: int,
    n_per_type: int,
    rng: np.random.RandomState,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic methylation + fragmentomic features.

    Methylation features: tissue-specific CpG island methylation patterns.
    Fragmentomic features: fragment size distribution, end motif, coverage.

    Each cancer type has a unique methylation signature at specific CpG sites.
    """
    cancer_types = CANCER_TYPES[:n_cancer_types]
    n_samples = n_per_type * n_cancer_types

    X = np.zeros((n_samples, N_TOTAL_FEATURES))
    y = np.zeros(n_samples, dtype=int)

    for i, ct in enumerate(cancer_types):
        start = i * n_per_type
        end = start + n_per_type

        # Methylation features: cancer-type-specific patterns
        # Each type has 2-3 "characteristic" methylation marks at specific loci
        meth_bases = np.zeros(N_METHYLATION_FEATURES)
        primary_methyl_idx = i % N_METHYLATION_FEATURES  # primary marker
        secondary_methyl_idx = (i + 3) % N_METHYLATION_FEATURES  # secondary marker
        meth_bases[primary_methyl_idx] = 0.85
        meth_bases[secondary_methyl_idx] = 0.60

        for j in range(n_per_type):
            # Methylation: high signal at characteristic loci + noise
            meth = meth_bases + rng.normal(0, 0.08, N_METHYLATION_FEATURES)
            meth = np.clip(meth, 0, 1)

            # Fragmentomic: type-specific fragmentation patterns
            frag_mean_size = 155 + rng.normal(0, 12)  # bp
            end_motif_ratio = 0.30 + rng.normal(0, 0.05)
            coverage_cv = 0.18 + rng.normal(0, 0.04)

            X[start + j] = np.concatenate([
                meth,
                [frag_mean_size / 200, end_motif_ratio, coverage_cv],
            ])

        y[start:end] = i

    return X, y


# ── Train TOO model with multi-class logistic regression ──────────────────
def _train_too_model(X: np.ndarray, y: np.ndarray,
                     n_splits: int = 5) -> Dict:
    """Train and evaluate TOO classifier."""
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)

    model = LogisticRegression(
        multi_class='multinomial',
        solver='lbfgs',
        max_iter=5000,
        random_state=SEED,
    )

    # Cross-validated accuracy
    cv_scores = cross_val_score(model, X, y, cv=skf, scoring='accuracy')
    cv_accuracy = float(np.mean(cv_scores))

    # Fit on full data for confusion matrix and per-class metrics
    model.fit(X, y)
    y_pred = model.predict(X)
    cm = confusion_matrix(y, y_pred)

    # Per-class accuracy
    per_class = {}
    for i, ct in enumerate(CANCER_TYPES[:len(np.unique(y))]):
        class_mask = y == i
        n_class = class_mask.sum()
        class_correct = (y_pred[class_mask] == i).sum()
        per_class[ct] = {
            'accuracy': float(class_correct / n_class),
            'n_samples': int(n_class),
            'n_correct': int(class_correct),
        }

    return {
        'cv_accuracy': cv_accuracy,
        'cv_accuracy_std': float(np.std(cv_scores)),
        'overall_accuracy': float(accuracy_score(y, y_pred)),
        'per_class': per_class,
        'confusion_matrix': cm.tolist(),
        'model_coefficients': {
            'intercept_shape': list(model.intercept_.shape),
            'coef_shape': list(model.coef_.shape),
        },
    }


# ── Main TOO Validation ───────────────────────────────────────────────────
def run_too_validation(
    n_cancer_types: int = DEFAULT_N_CANCER_TYPES,
    n_per_type: int = DEFAULT_N_PER_TYPE,
    seed: int = SEED,
) -> Dict:
    """
    Multi-class logistic regression on simulated methylation + fragmentomic
    features to predict tissue of origin.

    ⚠️ WARNING: This is SIMULATION ONLY. Features are synthetic with known
    ground truth. Real TOO on heterogeneous cancer types with mixed sample
    quality is dramatically harder.

    Returns: Performance dict with honest caveats.
    """
    rng = np.random.RandomState(seed)

    logger.info(f"TOO: {n_cancer_types} cancer types, {n_per_type} per type "
                f"({n_cancer_types * n_per_type} total)")

    X, y = _generate_too_features(n_cancer_types, n_per_type, rng)
    results = _train_too_model(X, y)

    # Honest comparison note
    grail_too_accuracy = 0.887
    cancerseeek_too_accuracy = 0.83

    output = {
        'metadata': {
            'generated': True,
            'validation_type': 'SIMULATION ONLY',
            'cancer_types': CANCER_TYPES[:n_cancer_types],
            'n_per_type': n_per_type,
            'features': f'{N_METHYLATION_FEATURES} methylation + {N_FRAGMENTOMIC_FEATURES} fragmentomic',
            'seed': seed,
        },
        'performance': results,
        'clinical_comparison': {
            'grail_too_accuracy': grail_too_accuracy,
            'grail_citation': 'Jamshidi 2022 Cancer Cell',
            'grail_note': 'Clinical TOO across 50+ cancer types',
            'cancerseeek_too_accuracy': cancerseeek_too_accuracy,
            'cancerseeek_citation': 'Cohen 2018 Science',
            'cancerseeek_note': 'Clinical TOO across 8 cancer types',
            'deepcatch_note': 'SIMULATION ONLY — not directly comparable',
        },
        'honest_assessment': (
            '❌ DEEPCATCH TOO IS NOT PROVEN. '
            'This module uses synthetic features with known ground truth — '
            'meaningless for Nature-level publication. '
            f'Grail achieved 88.7% TOO accuracy on CLINICAL samples across 50+ types. '
            f'CancerSEEK achieved 83% TOO on CLINICAL samples across 8 types. '
            'DeepCatch TOO accuracy on real heterogeneous samples is unknown.'
        ),
    }

    # Save
    with open(PY_TOO_PATH, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    logger.info(f"Saved TOO results to {PY_TOO_PATH}")

    logger.info(f"  CV Accuracy: {results['cv_accuracy']*100:.1f}% "
                f"(±{results['cv_accuracy_std']*100:.1f}%)")
    logger.info(f"  ⚠️  Grail clinical TOO: {grail_too_accuracy*100:.1f}%")
    logger.info(f"  ⚠️  SIMULATION ONLY — NOT comparable to clinical TOO")

    return output


# ── Demo ──────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    print("=" * 60)
    print("Tissue-of-Origin Validation — Demo")
    print("=" * 60)

    results = run_too_validation(n_cancer_types=5, n_per_type=30, seed=42)

    perf = results['performance']
    print(f"\nCV Accuracy: {perf['cv_accuracy']*100:.1f}%")
    print(f"Per-class accuracy:")
    for ct, stats in perf['per_class'].items():
        print(f"  {ct}: {stats['accuracy']*100:.1f}% ({stats['n_correct']}/{stats['n_samples']})")
    print(f"\n{results['honest_assessment']}")
    print("\n✅ TOO validation complete.")
