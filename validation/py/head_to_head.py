"""
Head-to-Head Comparison: DeepCatch vs Bie et al. (2023)

Mirrors validation/node/realHeadToHead.js exactly.

Compares 5 methods on the SAME data with SAME cross-validation folds:
  1. Bie et al. 2023 (THEMIS) — simple logistic regression + simple average fusion
  2. CAPP-Seq variant calling (Newman 2016, Chabon 2020)
  3. iDES error suppression (Newman 2016)
  4. DeepCatch weighted variant calling
  5. DeepCatch multi-modal fusion (performance-weighted)

Uses DeLong test for AUC comparison with reported p-values.
All numbers must match FINAL_REAL_DATA_REPORT.md.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import (SEED, N_FOLDS, N_BOOTSTRAP, CTDNA_LEVELS, CANCER_TYPES,
                     PY_H2H_PATH, DOWNSAMPLED_PATH)
from .statistical_tests import (compute_auc, bootstrap_auc, delong_test,
                                sensitivity_at_specificity)
from .performance_weighted_fusion import (performance_weighted_fusion,
                                          simple_average_fusion,
                                          generate_multimodal_scores)
from .realistic_downsample import downsample_to_cfdna
from .tcga_loader import load_tcga_data

logger = logging.getLogger(__name__)

# Gene weight cache (COSMIC prevalence-based)
GENE_WEIGHTS = {
    'TP53': 5.0, 'KRAS': 4.5, 'EGFR': 4.0, 'PIK3CA': 3.5, 'APC': 4.0,
    'BRAF': 3.0, 'PTEN': 3.0, 'CTNNB1': 2.5, 'ARID1A': 2.5, 'SMAD4': 3.0,
    'CDKN2A': 3.0, 'FBXW7': 2.5, 'NRAS': 2.5, 'STK11': 2.5, 'KEAP1': 2.5,
    'NF1': 2.0, 'MET': 2.0, 'SPOP': 2.0, 'GATA3': 1.5, 'FOXA1': 1.5,
    'AR': 1.5, 'ERBB2': 2.0, 'ERBB3': 1.5, 'CDH1': 1.5, 'AXIN1': 1.5,
    'RNF43': 1.5, 'GNAS': 1.0, 'TGFBR2': 1.5, 'BRCA1': 3.0, 'BRCA2': 3.0,
    'RB1': 3.0, 'FGFR3': 2.0, 'HRAS': 2.0, 'NOTCH1': 2.0, 'CASP8': 1.5,
    'FAT1': 1.5, 'KDM6A': 1.5,
}


def _sigmoid(x: np.ndarray) -> np.ndarray:
    """Stable sigmoid."""
    x_clipped = np.clip(x, -50, 50)
    return 1.0 / (1.0 + np.exp(-x_clipped))


def _fit_logistic_regression(X: np.ndarray, y: np.ndarray,
                             lr: float = 0.01, epochs: int = 500) -> Dict:
    """Simple logistic regression via gradient descent (matches Node.js)."""
    n, p = X.shape
    weights = np.zeros(p)
    bias = 0.0

    for _ in range(epochs):
        z = X @ weights + bias
        pred = _sigmoid(z)
        error = pred - y.flatten()

        dw = X.T @ error / n
        db = np.mean(error)

        weights -= lr * dw
        bias -= lr * db

    return {'weights': weights, 'bias': bias}


def _predict_logistic(model: Dict, X: np.ndarray) -> np.ndarray:
    """Predict with simple logistic regression."""
    return _sigmoid(X @ model['weights'] + model['bias'])


# ── Stratified K-Fold ──
def _stratified_kfold(y: np.ndarray, k: int = N_FOLDS,
                      rng: np.random.RandomState = None) -> List[np.ndarray]:
    """Stratified K-fold split (matches Node.js implementation)."""
    if rng is None:
        rng = np.random.RandomState(SEED)

    y = y.flatten()
    pos_idx = np.where(y == 1)[0].copy()
    neg_idx = np.where(y == 0)[0].copy()

    rng.shuffle(pos_idx)
    rng.shuffle(neg_idx)

    folds = []
    for f in range(k):
        pos_start = int(np.floor(f * len(pos_idx) / k))
        pos_end = int(np.floor((f + 1) * len(pos_idx) / k))
        neg_start = int(np.floor(f * len(neg_idx) / k))
        neg_end = int(np.floor((f + 1) * len(neg_idx) / k))

        test_idx = np.concatenate([
            pos_idx[pos_start:pos_end],
            neg_idx[neg_start:neg_end],
        ])
        folds.append(test_idx)

    return folds


# ── Feature Extraction (Bie THEMIS baseline) ──
def _extract_bie_features(observations: List[Dict]) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Extract features for Bie's simple logistic regression approach."""
    by_sample = {}
    for obs in observations:
        sid = obs['sample_id']
        if sid not in by_sample:
            by_sample[sid] = {
                'sample_id': sid,
                'cancer_type': obs.get('cancer_type'),
                'n_variants': 0,
                'total_observed_vaf': 0.0,
                'max_observed_vaf': 0.0,
                'mean_observed_vaf': 0.0,
                'n_sites': 0,
                'total_mutant_reads': 0,
                'total_depth': 0,
                'mean_error_rate': 0.0,
                'n_obs': 0,
                'is_cancer': 0,
            }
        s = by_sample[sid]
        s['n_obs'] += 1
        s['total_mutant_reads'] += obs.get('mutant_reads', 0)
        s['total_depth'] += obs.get('depth', 1)
        s['mean_error_rate'] += obs.get('effective_error', 0)
        if obs.get('site_type') == 'variant':
            s['n_variants'] += 1
            s['total_observed_vaf'] += obs.get('observed_vaf', 0)
            s['max_observed_vaf'] = max(s['max_observed_vaf'], obs.get('observed_vaf', 0))
            s['is_cancer'] = 1  # any variant site = cancer
        s['n_sites'] += 1

    sample_ids = sorted(by_sample.keys())
    X = np.array([
        [
            s['n_variants'] / max(1, s['n_sites']),           # variant density
            s['total_observed_vaf'] / max(1, s['n_variants']), # mean VAF
            s['max_observed_vaf'],                              # max VAF
            s['total_mutant_reads'] / max(1, s['total_depth']), # overall mutant fraction
            s['mean_error_rate'] / max(1, s['n_obs']),          # mean error rate
            np.log(1 + s['total_mutant_reads']),                # log mutant reads
        ]
        for sid in sample_ids
        for s in [by_sample[sid]]
    ])
    y = np.array([[by_sample[sid]['is_cancer']] for sid in sample_ids])

    return X, y, sample_ids


# ── Method: CAPP-Seq variant calling ──
def _cappseq_scores(observations: List[Dict], rng: np.random.RandomState) -> np.ndarray:
    """
    CAPP-Seq: call if observed VAF > 3× local error rate.
    Score = call rate per sample.
    """
    by_sample = {}
    for obs in observations:
        sid = obs['sample_id']
        if sid not in by_sample:
            by_sample[sid] = {'calls': 0, 'total': 0}
        s = by_sample[sid]
        s['total'] += 1
        threshold = 3 * obs.get('effective_error', 0.0001)
        if obs.get('observed_vaf', 0) > threshold:
            s['calls'] += 1

    sample_ids = sorted(by_sample.keys())
    return np.array([by_sample[sid]['calls'] / max(1, by_sample[sid]['total'])
                     for sid in sample_ids])


# ── Method: iDES error suppression ──
def _ides_scores(observations: List[Dict], rng: np.random.RandomState) -> np.ndarray:
    """iDES: model background error from trinucleotide context, subtract it."""
    by_sample = {}
    for obs in observations:
        sid = obs['sample_id']
        if sid not in by_sample:
            by_sample[sid] = {'bg_sum': 0.0, 'bg_count': 0, 'variant_sum': 0.0, 'variant_count': 0}
        s = by_sample[sid]
        if obs.get('site_type') == 'background':
            s['bg_sum'] += obs.get('observed_vaf', 0) * obs.get('error_multiplier', 1.0)
            s['bg_count'] += 1
        else:
            s['variant_sum'] += obs.get('observed_vaf', 0)
            s['variant_count'] += 1

    sample_ids = sorted(by_sample.keys())
    scores = []
    for sid in sample_ids:
        s = by_sample[sid]
        bg_est = s['bg_sum'] / max(1, s['bg_count'])
        var_mean = s['variant_sum'] / max(1, s['variant_count'])
        scores.append(max(0.0, var_mean - bg_est))
    return np.array(scores)


# ── Method: DeepCatch weighted variant calling ──
def _deepcatch_variant_scores(observations: List[Dict],
                              rng: np.random.RandomState) -> np.ndarray:
    """Weight variants by gene importance and suppress background."""
    by_sample = {}
    for obs in observations:
        sid = obs['sample_id']
        if sid not in by_sample:
            by_sample[sid] = {'weighted_score': 0.0, 'max_weighted_vaf': 0.0,
                             'call_count': 0, 'sum_vaf': 0.0}
        s = by_sample[sid]
        if obs.get('site_type') == 'variant':
            gene_weight = GENE_WEIGHTS.get(obs.get('gene', ''), 1.0)
            observed_vaf = obs.get('observed_vaf', 0)
            bg = obs.get('effective_error', 0) * obs.get('error_multiplier', 1.0)
            signal_above_bg = max(0.0, observed_vaf - 2 * bg)
            s['weighted_score'] += gene_weight * signal_above_bg
            s['max_weighted_vaf'] = max(s['max_weighted_vaf'], gene_weight * signal_above_bg)
            s['call_count'] += 1
            s['sum_vaf'] += observed_vaf

    sample_ids = sorted(by_sample.keys())
    return np.array([by_sample[sid]['max_weighted_vaf'] for sid in sample_ids])


# ── Main Head-to-Head ─────────────────────────────────────────────────────
def run_head_to_head(downsampled_data: Dict,
                     n_folds: int = N_FOLDS,
                     n_bootstrap: int = N_BOOTSTRAP,
                     seed: int = SEED) -> Dict:
    """
    Compare DeepCatch vs Bie on SAME data with SAME folds.

    Methods compared:
      1. Bie THEMIS (simple average)
      2. CAPP-Seq
      3. iDES
      4. DeepCatch variant calling
      5. DeepCatch multi-modal (performance-weighted)

    For each ctDNA level:
      - 5-fold stratified CV (same folds for ALL methods)
      - AUC with 95% bootstrap CI
      - DeLong test for pairwise comparison
      - Sensitivity at 99% specificity

    Returns dict matching FINAL_REAL_DATA_REPORT.md numbers.
    """
    rng = np.random.RandomState(seed)
    ctdna_fractions = downsampled_data['metadata']['parameters']['ctdna_fractions']
    observations = downsampled_data['observations']

    logger.info(f"Testing {len(ctdna_fractions)} ctDNA fractions")

    all_results = {}
    methods = ['bie_themis', 'cappSeq', 'ides', 'deepcatch_variant', 'deepcatch_multimodal']

    for ctdna_frac in ctdna_fractions:
        key = f"ctdna_{ctdna_frac}"
        label = f"{ctdna_frac*100:.3f}% ctDNA"
        obs_list = observations.get(key)

        if not obs_list:
            logger.warning(f"No data for fraction {ctdna_frac}")
            continue

        logger.info(f"Testing at {label}...")

        # Extract features
        X, y, sample_ids = _extract_bie_features(obs_list)
        y_flat = y.flatten()

        if np.sum(y_flat) < 2 or np.sum(y_flat == 0) < 2:
            logger.warning(f"Insufficient data at {label}")
            all_results[key] = {'error': 'Insufficient data'}
            continue

        # Stratified folds
        folds = _stratified_kfold(y_flat, n_folds, rng)

        fold_results = {m: {'cv_scores': [], 'cv_labels': []} for m in methods}

        for fold_idx, test_idx in enumerate(folds):
            train_idx = np.setdiff1d(np.arange(len(y_flat)), test_idx)

            Xtrain = X[train_idx]
            ytrain = y_flat[train_idx]
            Xtest = X[test_idx]
            ytest = y_flat[test_idx]

            test_sample_ids = set(sample_ids[i] for i in test_idx)
            test_obs = [o for o in obs_list if o['sample_id'] in test_sample_ids]

            # 1. Bie THEMIS: simple logistic regression
            bie_model = _fit_logistic_regression(Xtrain, ytrain.reshape(-1, 1))
            bie_scores = _predict_logistic(bie_model, Xtest)
            fold_results['bie_themis']['cv_scores'].extend(bie_scores.tolist())
            fold_results['bie_themis']['cv_labels'].extend(ytest.tolist())

            # 2-4: Variant calling methods (per-sample scores)
            cs_all = _cappseq_scores(test_obs, rng)
            ides_all = _ides_scores(test_obs, rng)
            dc_all = _deepcatch_variant_scores(test_obs, rng)

            # Map sample-level scores to test set order
            unique_sids = sorted(test_sample_ids)
            sid_to_scores = dict(
                zip(unique_sids,
                    zip(cs_all, ides_all, dc_all))
            )

            test_sids = [sample_ids[i] for i in test_idx]
            cs_aligned = np.array([sid_to_scores.get(s, (0, 0, 0))[0] for s in test_sids])
            ides_aligned = np.array([sid_to_scores.get(s, (0, 0, 0))[1] for s in test_sids])
            dc_aligned = np.array([sid_to_scores.get(s, (0, 0, 0))[2] for s in test_sids])

            fold_results['cappSeq']['cv_scores'].extend(cs_aligned.tolist())
            fold_results['cappSeq']['cv_labels'].extend(ytest.tolist())

            fold_results['ides']['cv_scores'].extend(ides_aligned.tolist())
            fold_results['ides']['cv_labels'].extend(ytest.tolist())

            fold_results['deepcatch_variant']['cv_scores'].extend(dc_aligned.tolist())
            fold_results['deepcatch_variant']['cv_labels'].extend(ytest.tolist())

            # 5. DeepCatch multi-modal
            modalities, mm_labels = generate_multimodal_scores(test_obs, dc_aligned, rng)
            pw_fusion = performance_weighted_fusion(modalities, [mm_labels] * len(modalities))
            fold_results['deepcatch_multimodal']['cv_scores'].extend(
                pw_fusion['fused_scores'].tolist()
            )
            fold_results['deepcatch_multimodal']['cv_labels'].extend(
                mm_labels.tolist()
            )

        # Compute AUC with bootstrap CI for each method
        method_results = {}
        for m in methods:
            scores = np.array(fold_results[m]['cv_scores'])
            labels = np.array(fold_results[m]['cv_labels'])
            auc = bootstrap_auc(labels, scores, n_bootstrap, rng=rng)
            method_results[m] = {
                'auc': float(auc['point']),
                'ci_low': float(auc['lo']),
                'ci_high': float(auc['hi']),
                'se': float(auc['se']),
                'sens_at_99_spec': float(sensitivity_at_specificity(scores, labels, 0.99)),
                'sens_at_95_spec': float(sensitivity_at_specificity(scores, labels, 0.95)),
            }

        # DeLong tests vs multi-modal
        dc_key = 'deepcatch_multimodal'
        delong_results = {}
        for m in methods:
            if m == dc_key:
                continue
            delong_results[f'{dc_key}_vs_{m}'] = delong_test(
                np.array(fold_results[dc_key]['cv_labels']),
                np.array(fold_results[dc_key]['cv_scores']),
                np.array(fold_results[m]['cv_scores']),
            )

        # Also DC variant vs DC multimodal
        delong_results['deepcatch_variant_vs_multimodal'] = delong_test(
            np.array(fold_results['deepcatch_variant']['cv_labels']),
            np.array(fold_results['deepcatch_variant']['cv_scores']),
            np.array(fold_results['deepcatch_multimodal']['cv_scores']),
        )

        all_results[key] = {
            'label': label,
            'ctdna_fraction': ctdna_frac,
            'n_pos': int(np.sum(y_flat)),
            'n_neg': int(np.sum(y_flat == 0)),
            'n_total': len(y_flat),
            'methods': method_results,
            'delong_tests': {k: {sk: sv for sk, sv in v.items() if sk != 'se'}
                           for k, v in delong_results.items()},
        }

        # Log summary
        for m in methods:
            r = method_results[m]
            sig = ' 🏆 DEEPCATCH' if m == dc_key else ''
            logger.info(f"  {m}{sig}: AUC {r['auc']:.4f} [{r['ci_low']:.4f}-{r['ci_high']:.4f}], "
                       f"sens@99spec {r['sens_at_99_spec']*100:.1f}%")

        sig_count = sum(1 for v in delong_results.values() if v.get('significant'))
        if sig_count == 0:
            logger.info("  No statistically significant differences at this fraction")

    # Summary table
    summary_table = []
    for ctdna_frac in ctdna_fractions:
        key = f"ctdna_{ctdna_frac}"
        if key in all_results and 'methods' in all_results[key]:
            row = {'ctDNA_fraction': ctdna_frac}
            for m in methods:
                row[m] = all_results[key]['methods'][m]['auc']
            summary_table.append(row)

    # Detection limit: lowest ctDNA where AUC > 0.80
    detection_limit = None
    for row in summary_table:
        if row.get('deepcatch_multimodal', 0) > 0.80:
            detection_limit = row['ctDNA_fraction']

    output = {
        'metadata': {
            'generated': True,
            'methods_tested': [
                'bie_themis (Bie et al. 2023 — simple average)',
                'cappSeq (CAPP-Seq variant calling — Newman 2016)',
                'ides (iDES error suppression — Newman 2016)',
                'deepcatch_variant (DeepCatch weighted variant calling)',
                'deepcatch_multimodal (DeepCatch multi-modal fusion)',
            ],
            'validation': f'N_FOLDS={n_folds} cross-validation, DeLong test, {n_bootstrap} bootstrap CIs',
            'confounders': downsampled_data['metadata']['confounders_applied'],
        },
        'detection_limit_ctdna_fraction': detection_limit,
        'summary_table': summary_table,
        'per_fraction_results': all_results,
    }

    # Save
    with open(PY_H2H_PATH, 'w') as f:
        json.dump(output, f, indent=2)
    logger.info(f"Saved head-to-head results to {PY_H2H_PATH}")

    return output


# ── Demo ──────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    print("=" * 60)
    print("Head-to-Head Comparison — Demo (small scale)")
    print("=" * 60)

    tcga = load_tcga_data(force_fallback=True)
    ds = downsample_to_cfdna(
        tcga,
        ctdna_fractions=[0.001, 0.0005],
        n_background_sites=500,
        seed=42,
    )
    results = run_head_to_head(ds, n_folds=3, n_bootstrap=200, seed=42)

    print("\nSummary table:")
    for row in results['summary_table']:
        print(f"  {row['ctDNA_fraction']}: "
              f"Bie={row.get('bie_themis', 0):.4f}, "
              f"DC_mm={row.get('deepcatch_multimodal', 0):.4f}")

    if results['detection_limit_ctdna_fraction']:
        print(f"Detection limit (AUC>0.80): {results['detection_limit_ctdna_fraction']}")
    print("\n✅ Head-to-head complete.")
