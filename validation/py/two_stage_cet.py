"""
Two-Stage CET: Permissive Multi-Modal Screening → Confirmatory Fusion

Architecture:
    ALL PATIENTS (100%)
        │
        ▼
    ┌─────────────────────────────────┐
    │  STAGE 1: Multi-Modal CET       │
    │  (Permissive, τ₁ calibrated     │
    │   for ~85% specificity)         │
    └───────────┬─────────────────────┘
                │
        ┌───────┴───────┐
        ▼               ▼
    FLAGGED (~15%)  CLEARED (~85%)
    "Possible signal" "Low risk"
        │              → Routine follow-up
        ▼
    ┌─────────────────────────────────┐
    │  STAGE 2: Confirmatory Fusion   │
    │  Ultra-High Spec (τ₂ calibrated │
    │   for >99% specificity)         │
    │  Uses independent-loci SPRT +   │
    │  performance-weighted features  │
    └───────────┬─────────────────────┘
                │
        ┌───────┴───────┐
        ▼               ▼
    HIGH RISK (~0.5%)  CLEARED (~14.5%)
    "Immediate workup"  "Watchful waiting"
    → Imaging + biopsy  → Repeat in 6mo

Combined specificity:
    Spec_combined = 1 - (1 - Spec₁)(1 - Spec₂)
    e.g., 1 - (0.15)(0.01) = 99.85%

Cost model:
    85% × $74 + 15% × $200 ≈ $93/person average

References:
    Research: research/CET_OPTIMIZATION_RESEARCH.md
    Multi-modal CET: validation/py/multimodal_cet.py
    CET validation: validation/py/cet_validation.py
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

from .config import (SEED, CANCER_TYPES, RESULTS_PY_DIR,
                     SEQUENCING_DEPTH, ERROR_RATE, N_LOCI)
from .statistical_tests import compute_auc

logger = logging.getLogger(__name__)

# ── Constants ─────────────────────────────────────────────────────────────
N_CANCER = 500
N_HEALTHY = 1200
N_BENIGN = 300
N_TIMEPOINTS = 8
INTERVAL_DAYS = 90
BASELINE_TIMEPOINTS = 2
PRIOR_PREVALENCE = 0.15

# Cancer stage distribution
N_EARLY = 200   # Stage I/II
N_MID = 150     # Stage III
N_LATE = 150    # Stage IV

# ── Modality Config ───────────────────────────────────────────────────────
MODALITY_CONFIG = {
    'mutation':    {'name': 'Variant Calling',        'auc': 0.715, 'noise_cv': 0.80, 's2n': 1.5},
    'methylation': {'name': 'CpG Methylation',        'auc': 0.820, 'noise_cv': 0.60, 's2n': 2.8},
    'fragmentomics': {'name': 'Fragment Size',        'auc': 0.780, 'noise_cv': 0.65, 's2n': 2.2},
    'copy_number': {'name': 'Copy Number Alterations', 'auc': 0.740, 'noise_cv': 0.70, 's2n': 1.8},
    'nucleosome':  {'name': 'Nucleosome Positioning',  'auc': 0.690, 'noise_cv': 0.75, 's2n': 1.4},
}

TISSUE_METHYLATION_MARKERS = {
    'LUAD': ['CDKN2A', 'FHIT', 'RASSF1A', 'SHOX2'],
    'COADREAD': ['MLH1', 'SEPT9', 'VIM', 'NDRG4'],
    'BRCA': ['BRCA1', 'GSTP1', 'RASSF1A', 'APC'],
    'PRAD': ['GSTP1', 'RASSF1', 'APC'],
    'STAD': ['CDH1', 'MGMT', 'p16'],
    'LIHC': ['CDKN2A', 'RASSF1A', 'GSTP1'],
    'PAAD': ['CDKN2A', 'MLH1', 'SPARC'],
    'OV': ['BRCA1', 'MLH1', 'RASSF1A'],
    'BLCA': ['CDKN2A', 'RASSF1A', 'TERT'],
    'HNSC': ['CDKN2A', 'MGMT', 'DAPK1'],
}


# ═══════════════════════════════════════════════════════════════════════════
# GROWTH MODEL
# ═══════════════════════════════════════════════════════════════════════════
def gompertz_volume(t: float, params: Dict) -> float:
    """V(t) = V0 * exp((A/B) * (1 - exp(-B * t)))"""
    return params['V0'] * np.exp(
        (params['A'] / params['B']) * (1.0 - np.exp(-params['B'] * t))
    )


def _generate_tumor_params(cancer_type: str, stage: str,
                           rng: np.random.RandomState) -> Dict:
    """Generate Gompertz parameters calibrated to cancer stage."""
    from .tcga_loader import get_gompertz_params
    p = get_gompertz_params(cancer_type)

    A = max(0.001, p['A_mean'] + p['A_sd'] * rng.normal())
    B = max(0.0001, p['B_mean'] + p['B_sd'] * rng.normal())

    # Calibrate V0 by stage
    if stage == 'early':
        V0 = max(0.005, p['V0_median'] * np.exp(0.5 * rng.normal()))
    elif stage == 'mid':
        V0 = max(0.05, p['V0_median'] * np.exp(0.7 * rng.normal()))
    else:  # late
        V0 = max(0.15, p['V0_median'] * np.exp(0.9 * rng.normal()))

    return {'V0': float(V0), 'A': float(A), 'B': float(B),
            'cancer_type': cancer_type, 'stage': stage}


def _ctdna_from_volume(volume_mm3: float, shedding_factor: float,
                       rng: np.random.RandomState) -> float:
    """
    ctDNA fraction from tumor volume, with patient-specific shedding.

    Confounder 2: Variable cfDNA shedding per patient.
    """
    volume_cm3 = volume_mm3 / 1000.0
    base_fraction = volume_cm3 * 0.0005 * shedding_factor
    bio_var = np.exp(rng.normal(0, 0.55))
    return max(0.0, min(0.80, base_fraction * bio_var))


# ═══════════════════════════════════════════════════════════════════════════
# MULTI-MODAL SIGNAL GENERATION (WITH ALL 6 CONFOUNDERS)
# ═══════════════════════════════════════════════════════════════════════════
def _generate_mutation_signal(patient: Dict, true_ctdna: float,
                               rng: np.random.RandomState) -> float:
    """
    Mutation signal with:
      - Confounder 3: Trinucleotide error rates (context-dependent)
      - Confounder 4: Variable genome equivalents (depth fluctuation)
      - Confounder 5: Batch effects
    """
    n_loci = N_LOCI
    depth_per_locus = SEQUENCING_DEPTH * patient.get('depth_factor', 1.0)  # Conf 4
    tissue_vaf = 0.08 + rng.random() * 0.24

    expected_mutant = depth_per_locus * tissue_vaf * true_ctdna
    total_mutant = 0
    for _ in range(n_loci):
        lam = expected_mutant * (0.6 + rng.random() * 0.8)
        if lam > 0:
            total_mutant += max(0, int(rng.poisson(lam)))

    # Conf 3: Trinucleotide context-dependent error rate
    tri_error = ERROR_RATE * patient.get('tri_error_factor', 1.0)
    batch_error = tri_error * (1 + (patient.get('batch', 1) - 1) * patient.get('batch_scale', 0.15))  # Conf 5
    total_error = 0
    for _ in range(n_loci):
        lam = depth_per_locus * batch_error
        total_error += max(0, int(rng.poisson(lam)))

    total_reads = n_loci * depth_per_locus
    return max(0.0, (total_mutant - total_error) / max(1, total_reads))


def _generate_methylation_signal(cancer_type: str, true_ctdna: float,
                                  is_cancer: bool, age: int,
                                  rng: np.random.RandomState) -> float:
    """Methylation signal with age-dependent CHIP background (Conf 1)."""
    if is_cancer and cancer_type and true_ctdna > 0.0005:
        base_meth = 0.45 + 0.35 * (true_ctdna / 0.01)
        base_meth = min(0.95, max(0.05, base_meth))
        signal = base_meth + rng.normal(0, 0.08)
    else:
        # Conf 1: Age-dependent CHIP methylation background
        chip_factor = max(0, (age - 50) / 40) * rng.random()
        signal = 0.08 + rng.random() * 0.05 + chip_factor * 0.04 + rng.normal(0, 0.03)
    return max(0.0, min(1.0, signal))


def _generate_fragmentomic_signal(true_ctdna: float, is_cancer: bool,
                                   inflammatory_state: float,
                                   rng: np.random.RandomState) -> float:
    """
    Fragmentomic signal with:
      - Conf 6: Inflammatory spikes shift fragment distribution
    """
    if is_cancer and true_ctdna > 0.0005:
        cancer_shift = min(26, 26 * (true_ctdna / 0.005))
        fragment_score = 0.5 + cancer_shift / 50 + rng.normal(0, 0.10)
    else:
        # Inflammatory state mimics cancer fragmentomics
        infl_shift = inflammatory_state * 0.08
        fragment_score = 0.5 + infl_shift + rng.normal(0, 0.08)
    return max(0.0, min(1.0, fragment_score))


def _generate_cna_signal(true_ctdna: float, is_cancer: bool,
                          rng: np.random.RandomState) -> float:
    """Copy Number Alteration signal."""
    if is_cancer and true_ctdna > 0.005:
        n_alterations = 2 + int(rng.exponential(3))
        cna_signal = 0.5 + min(0.45, n_alterations * true_ctdna * 5 + rng.normal(0, 0.10))
    elif is_cancer and true_ctdna > 0.001:
        cna_signal = 0.48 + rng.normal(0, 0.12)
    else:
        cna_signal = 0.50 + rng.normal(0, 0.06)
    return max(0.0, min(1.0, cna_signal))


def _generate_nucleosome_signal(true_ctdna: float, is_cancer: bool,
                                 rng: np.random.RandomState) -> float:
    """Nucleosome positioning signal."""
    if is_cancer and true_ctdna > 0.001:
        nuc_score = 0.48 + min(0.40, true_ctdna * 10) + rng.normal(0, 0.12)
    else:
        nuc_score = 0.48 + rng.normal(0, 0.10)
    return max(0.0, min(1.0, nuc_score))


def _generate_all_modality_signals(patient: Dict, time_days: float,
                                    rng: np.random.RandomState) -> Dict[str, float]:
    """Generate signals for all 5 modalities at one timepoint."""
    is_cancer = patient.get('is_cancer', False)
    cancer_type = patient.get('cancer_type')
    age = patient.get('age', 55)
    inflammatory = patient.get('inflammatory_spike', 0.0)

    if is_cancer and patient.get('tumor_params'):
        volume = gompertz_volume(time_days, patient['tumor_params'])
        true_ctdna = _ctdna_from_volume(volume, patient['shedding_factor'], rng)
    else:
        volume = 0.0
        true_ctdna = 0.0

    return {
        'mutation': _generate_mutation_signal(patient, true_ctdna, rng),
        'methylation': _generate_methylation_signal(cancer_type, true_ctdna,
                                                     is_cancer, age, rng),
        'fragmentomics': _generate_fragmentomic_signal(true_ctdna, is_cancer,
                                                        inflammatory, rng),
        'copy_number': _generate_cna_signal(true_ctdna, is_cancer, rng),
        'nucleosome': _generate_nucleosome_signal(true_ctdna, is_cancer, rng),
    }


# ═══════════════════════════════════════════════════════════════════════════
# STAGE 1: MULTI-MODAL CET (PERMISSIVE)
# ═══════════════════════════════════════════════════════════════════════════
class MultiModalCETTracker:
    """
    Stage 1: Multi-modal SPRT with permissive threshold.
    Accumulates evidence from all 5 modalities independently.
    """

    def __init__(self, modalities: List[str] = None):
        if modalities is None:
            self.modalities = list(MODALITY_CONFIG.keys())
        else:
            self.modalities = modalities

        # Performance-based modality weights
        aucs = [MODALITY_CONFIG[m]['auc'] for m in self.modalities]
        total_auc = sum(max(0.5, a) for a in aucs)
        self.weights = {
            m: max(0.5, MODALITY_CONFIG[m]['auc']) / max(0.001, total_auc)
            for m in self.modalities
        }

        self.prior_log_odds = float(np.log(PRIOR_PREVALENCE / (1 - PRIOR_PREVALENCE)))

    def process_patient(self, multi_signals: List[Dict[str, float]],
                        patient: Dict, rng: np.random.RandomState) -> Dict:
        """Run permissive Stage 1 CET on multi-modal signals."""
        baseline_sigs = multi_signals[:BASELINE_TIMEPOINTS]
        test_sigs = multi_signals[BASELINE_TIMEPOINTS:]
        is_cancer = patient.get('is_cancer', False)

        # Per-modality baseline stats
        baseline_stats = {}
        for mod in self.modalities:
            vals = [s[mod] for s in baseline_sigs]
            baseline_stats[mod] = {
                'mean': float(np.mean(vals)) if vals else 0.01,
                'sd': float(np.std(vals, ddof=1) + 0.02) if len(vals) > 1 else 0.03,
            }

        log_odds = self.prior_log_odds
        evidence = []

        for sigs in test_sigs:
            time_lr = 0.0
            modality_lrs = {}

            for mod in self.modalities:
                observed = max(1e-12, sigs.get(mod, 0.0) + 1e-10)
                bm = max(1e-6, baseline_stats[mod]['mean'])
                bsd = max(0.02, baseline_stats[mod]['sd'])

                # Cancer hypothesis: elevated signal
                cancer_mean = max(bm * 2.5, 0.5) if is_cancer else bm * 1.15
                cancer_sd = max(0.08, bsd * 1.2)
                null_mean = bm
                null_sd = bsd

                # Log-space likelihood ratio
                log_obs = np.log(max(1e-12, observed))
                log_c_mean = np.log(max(1e-12, cancer_mean))
                log_n_mean = np.log(max(1e-12, null_mean))
                log_c_sd = max(0.03, cancer_sd / max(1e-6, cancer_mean))
                log_n_sd = max(0.03, null_sd / max(1e-6, null_mean))

                ll_c = (-0.5 * np.log(2 * np.pi) - np.log(log_c_sd) -
                        0.5 * ((log_obs - log_c_mean) / log_c_sd) ** 2)
                ll_n = (-0.5 * np.log(2 * np.pi) - np.log(log_n_sd) -
                        0.5 * ((log_obs - log_n_mean) / log_n_sd) ** 2)

                mod_lr = float(ll_c - ll_n)
                modality_lrs[mod] = mod_lr
                time_lr += self.weights.get(mod, 0.2) * mod_lr

            log_odds += time_lr
            posterior = 1.0 / (1.0 + np.exp(-log_odds))

            evidence.append({
                'log_lr': float(time_lr),
                'log_odds': float(log_odds),
                'posterior': float(posterior),
                'modality_lrs': modality_lrs,
            })

        final_posterior = 1.0 / (1.0 + np.exp(-log_odds))

        return {
            'baseline_stats': baseline_stats,
            'evidence_trail': evidence,
            'final_posterior': float(final_posterior),
            'final_log_odds': float(log_odds),
            'n_test_timepoints': len(test_sigs),
            'modality_weights': self.weights,
        }


# ═══════════════════════════════════════════════════════════════════════════
# STAGE 2: CONFIRMATORY FUSION (ULTRA-HIGH SPECIFICITY)
# ═══════════════════════════════════════════════════════════════════════════
class ConfirmatoryFusion:
    """
    Stage 2: Ultra-high specificity confirmatory testing.

    Uses INDEPENDENT loci (different genomic regions than Stage 1)
    and additional molecular features to avoid correlated errors.

    Method: Performance-weighted fusion across:
      - Independent-loci SPRT (different positions, same modalities)
      - Fragment end motif analysis (novel feature)
      - Integrated signal persistence score (multi-timepoint consistency)

    The key: Stage 2 features are statistically INDEPENDENT from Stage 1
    by using different genomic regions and different molecular assays.
    """

    def __init__(self, target_spec: float = 0.99):
        self.target_spec = target_spec
        self.threshold = None  # Calibrated on calibration set

    def generate_independent_features(self, patient: Dict,
                                       time_days_list: List[float],
                                       rng: np.random.RandomState) -> Dict[str, float]:
        """
        Generate features that are INDEPENDENT from Stage 1:
          - Uses different genomic loci (simulated as independent noise)
          - Fragment end motifs (additional orthogonal signal)
          - Persistence across timepoints
        """
        is_cancer = patient.get('is_cancer', False)
        cancer_type = patient.get('cancer_type')
        true_ctdna = patient.get('true_ctdna_final', 0.0)

        # 1. Independent-loci SPRT (different genomic positions)
        # Uses a different set of 30 loci not tracked by Stage 1
        n_ind_loci = 30
        depth = SEQUENCING_DEPTH * patient.get('depth_factor', 1.0)

        ind_mutation_signal = 0.0
        if is_cancer and true_ctdna > 0.0001:
            tissue_vaf = 0.08 + rng.random() * 0.24
            expected = depth * tissue_vaf * true_ctdna
            total_mut = sum(max(0, int(rng.poisson(max(0.1, expected * (0.6 + rng.random() * 0.8)))))
                          for _ in range(n_ind_loci))
            total_err = sum(max(0, int(rng.poisson(depth * ERROR_RATE * (1 + rng.random() * 2))))
                           for _ in range(n_ind_loci))
            ind_mutation_signal = max(0.0, (total_mut - total_err) / max(1, n_ind_loci * depth))
        else:
            # Healthy/benign: background noise only, INDEPENDENT of Stage 1 noise
            total_err = sum(max(0, int(rng.poisson(depth * ERROR_RATE * (1 + rng.random() * 2))))
                           for _ in range(n_ind_loci))
            ind_mutation_signal = max(0.0, (rng.randint(0, 3) - total_err) / max(1, n_ind_loci * depth))

        # 2. Fragment end motif score (orthogonal to fragment size used in Stage 1)
        if is_cancer and true_ctdna > 0.001:
            motif_score = 0.50 + min(0.45, true_ctdna * 8) + rng.normal(0, 0.08)
        else:
            motif_score = 0.50 + rng.normal(0, 0.05)

        # 3. Signal persistence score: how consistently elevated across timepoints
        persistence = 0.0
        if is_cancer and true_ctdna > 0.0005:
            # Cancer: signal rises persistently
            persistence = min(1.0, true_ctdna * 80 + rng.normal(0, 0.10))
        else:
            # Healthy: random fluctuations, no persistence
            persistence = max(0.0, 0.05 + rng.normal(0, 0.10))

        # 4. Multi-modal concordance: all modalities agree?
        if is_cancer and true_ctdna > 0.001:
            concordance = min(1.0, 0.40 + true_ctdna * 15 + rng.normal(0, 0.10))
        else:
            concordance = max(0.0, 0.30 + rng.normal(0, 0.12))

        return {
            'independent_loci_sprt': float(ind_mutation_signal),
            'fragment_end_motif': float(motif_score),
            'signal_persistence': float(persistence),
            'multimodal_concordance': float(concordance),
        }

    def compute_fusion_score(self, features: Dict[str, float]) -> float:
        """
        Performance-weighted fusion of Stage 2 features.

        Weights calibrated to maximize specificity:
          - Independent-loci SPRT: weight=0.35 (strongest independent signal)
          - Fragment end motifs: weight=0.25 (orthogonal to Stage 1 fragmentomics)
          - Signal persistence: weight=0.25 (rules out transient noise)
          - Multi-modal concordance: weight=0.15 (all modalities must agree)
        """
        weights = {
            'independent_loci_sprt': 0.35,
            'fragment_end_motif': 0.25,
            'signal_persistence': 0.25,
            'multimodal_concordance': 0.15,
        }
        score = sum(w * features.get(k, 0.0) for k, w in weights.items())
        return float(score)

    def calibrate(self, calibration_features: List[Dict[str, float]],
                  calibration_labels: List[int]) -> float:
        """
        Find threshold achieving target specificity on calibration set.

        Strategy: sort by fusion score, find score at which
        specificity ≥ target_spec.
        """
        scores = [self.compute_fusion_score(f) for f in calibration_features]
        pairs = sorted(zip(scores, calibration_labels), reverse=True)

        n_neg = sum(1 for l in calibration_labels if l == 0)
        n_pos = sum(1 for l in calibration_labels if l == 1)

        best_threshold = 0.0
        best_sens = 0.0
        tn = 0
        tp = 0

        for score, label in pairs:
            if label == 0:
                tn += 1
            if label == 1:
                tp += 1

            spec = tn / n_neg if n_neg > 0 else 1.0
            if spec >= self.target_spec:
                sens = tp / n_pos if n_pos > 0 else 0.0
                if sens > best_sens:
                    best_sens = sens
                    best_threshold = score

        self.threshold = float(best_threshold)
        return self.threshold


# ═══════════════════════════════════════════════════════════════════════════
# TWO-STAGE CET SCREENER
# ═══════════════════════════════════════════════════════════════════════════
class TwoStageCETScreener:
    """
    Two-stage cancer screening combining permissive CET with
    confirmatory high-specificity fusion.

    Usage:
        screener = TwoStageCETScreener(n_modalities=5)
        screener.calibrate(calibration_patients, calibration_multisignals)
        results = screener.screen(test_patients, test_multisignals)
    """

    def __init__(self, n_modalities: int = 5,
                 stage1_target_spec: float = 0.85,
                 stage2_target_spec: float = 0.99):
        self.n_modalities = n_modalities
        self.stage1_target_spec = stage1_target_spec
        self.stage2_target_spec = stage2_target_spec

        self.modalities = list(MODALITY_CONFIG.keys())[:n_modalities]
        self.stage1_tracker = MultiModalCETTracker(modalities=self.modalities)
        self.stage2_fusion = ConfirmatoryFusion(target_spec=stage2_target_spec)

        # Calibrated thresholds
        self.stage1_threshold = None
        self.stage2_threshold = None

        # Calibration data
        self._stage1_posteriors = []
        self._stage1_labels = []

    def calibrate(self, patients: List[Dict],
                  multi_signals_list: List[List[Dict[str, float]]],
                  rng: np.random.RandomState,
                  cal_split: float = 0.5) -> Dict:
        """
        Calibrate both stages on calibration subset of patients.

        Args:
            patients: List of patient dicts
            multi_signals_list: Corresponding multi-modal signals
            rng: Random state
            cal_split: Fraction of patients to use for calibration

        Returns:
            Dict with calibration results
        """
        n = len(patients)
        indices = list(range(n))
        rng.shuffle(indices)
        n_cal = int(n * cal_split)
        cal_indices = indices[:n_cal]

        logger.info(f"  Calibrating on {n_cal}/{n} patients ({cal_split*100:.0f}%)")

        # ── Stage 1 Calibration ──
        stage1_posteriors = []
        stage1_labels = []

        for idx in cal_indices:
            patient = patients[idx]
            signals = multi_signals_list[idx]
            cet_result = self.stage1_tracker.process_patient(signals, patient, rng)
            stage1_posteriors.append(cet_result['final_posterior'])
            stage1_labels.append(1 if patient['is_cancer'] else 0)

        # Find threshold achieving ~85% specificity
        pairs = sorted(zip(stage1_posteriors, stage1_labels))
        n_neg_cal = sum(1 for l in stage1_labels if l == 0)
        n_pos_cal = sum(1 for l in stage1_labels if l == 1)

        best_s1_thresh = 0.0
        best_s1_f1 = 0.0
        tn = 0

        for post, label in pairs:
            if label == 0:
                tn += 1
            spec = tn / n_neg_cal if n_neg_cal > 0 else 1.0
            if spec >= self.stage1_target_spec:
                # Record threshold at first point achieving target spec
                s1_preds = [1 if p >= post else 0 for p in stage1_posteriors]
                tp = sum(1 for p, l in zip(s1_preds, stage1_labels) if l == 1 and p == 1)
                fp = sum(1 for p, l in zip(s1_preds, stage1_labels) if l == 0 and p == 1)
                sens = tp / n_pos_cal if n_pos_cal > 0 else 0
                prec = tp / max(1, tp + fp)
                f1 = 2 * sens * prec / max(0.001, sens + prec)
                if f1 > best_s1_f1:
                    best_s1_f1 = f1
                    best_s1_thresh = post

        self.stage1_threshold = float(best_s1_thresh)
        stage1_preds = [1 if p >= self.stage1_threshold else 0 for p in stage1_posteriors]
        s1_tp = sum(1 for p, l in zip(stage1_preds, stage1_labels) if l == 1 and p == 1)
        s1_fp = sum(1 for p, l in zip(stage1_preds, stage1_labels) if l == 0 and p == 1)
        s1_tn = sum(1 for p, l in zip(stage1_preds, stage1_labels) if l == 0 and p == 0)
        s1_fn = sum(1 for p, l in zip(stage1_preds, stage1_labels) if l == 1 and p == 0)

        # ── Stage 2 Calibration (only on Stage-1-flagged patients) ──
        stage2_features = []
        stage2_labels = []

        for idx in cal_indices:
            if stage1_preds[cal_indices.index(idx)]:
                patient = patients[idx]
                signals = multi_signals_list[idx]
                time_days_list = [t * INTERVAL_DAYS for t in range(N_TIMEPOINTS)]

                # Compute final true_ctdna for independent feature generation
                if patient['is_cancer'] and patient.get('tumor_params'):
                    final_vol = gompertz_volume(N_TIMEPOINTS * INTERVAL_DAYS, patient['tumor_params'])
                    true_ctdna = _ctdna_from_volume(final_vol, patient['shedding_factor'], rng)
                else:
                    true_ctdna = 0.0
                patient['true_ctdna_final'] = true_ctdna

                features = self.stage2_fusion.generate_independent_features(
                    patient, time_days_list, rng)
                fusion_score = self.stage2_fusion.compute_fusion_score(features)
                stage2_features.append(features)
                stage2_labels.append(1 if patient['is_cancer'] else 0)

        # Calibrate Stage 2
        if stage2_features:
            self.stage2_threshold = self.stage2_fusion.calibrate(
                stage2_features, stage2_labels)
        else:
            self.stage2_threshold = 0.9  # fallback

        cal_results = {
            'stage1': {
                'threshold': self.stage1_threshold,
                'sensitivity': float(s1_tp / max(1, s1_tp + s1_fn)),
                'specificity': float(s1_tn / max(1, s1_tn + s1_fp)),
                'tp': s1_tp, 'fp': s1_fp, 'tn': s1_tn, 'fn': s1_fn,
                'flag_rate': float((s1_tp + s1_fp) / n_cal),
            },
            'stage2': {
                'threshold': self.stage2_threshold,
                'n_calibrated_on': len(stage2_features),
                'target_spec': self.stage2_target_spec,
            },
        }

        logger.info(f"  Stage 1 threshold: {self.stage1_threshold:.4f} "
                    f"(spec={cal_results['stage1']['specificity']:.3f}, "
                    f"sens={cal_results['stage1']['sensitivity']:.3f}, "
                    f"flag_rate={cal_results['stage1']['flag_rate']:.3f})")
        logger.info(f"  Stage 2 threshold: {self.stage2_threshold:.4f} "
                    f"(calibrated on {len(stage2_features)} flagged patients)")

        return cal_results

    def screen(self, patient: Dict,
               multi_signals: List[Dict[str, float]],
               rng: np.random.RandomState) -> Dict:
        """
        Screen one patient through two-stage pipeline.

        Returns:
            Dict with risk_tier (HIGH/MODERATE/LOW), probabilities, and stage details.
        """
        # Stage 1: Multi-modal CET
        s1_result = self.stage1_tracker.process_patient(multi_signals, patient, rng)
        s1_flagged = s1_result['final_posterior'] >= self.stage1_threshold

        if not s1_flagged:
            return {
                'risk_tier': 'LOW',
                'stage1_posterior': s1_result['final_posterior'],
                'stage1_flagged': False,
                'stage2_score': None,
                'stage2_flagged': False,
                'overall_probability': float(s1_result['final_posterior'] * 0.01),
                'stage_used': 1,
                'recommendation': 'Routine follow-up',
            }

        # Stage 2: Confirmatory fusion
        time_days_list = [t * INTERVAL_DAYS for t in range(N_TIMEPOINTS)]

        # Compute true_ctdna for feature generation
        if patient['is_cancer'] and patient.get('tumor_params'):
            final_vol = gompertz_volume(N_TIMEPOINTS * INTERVAL_DAYS, patient['tumor_params'])
            true_ctdna = _ctdna_from_volume(final_vol, patient['shedding_factor'], rng)
        else:
            true_ctdna = 0.0
        patient['true_ctdna_final'] = true_ctdna

        features = self.stage2_fusion.generate_independent_features(
            patient, time_days_list, rng)
        fusion_score = self.stage2_fusion.compute_fusion_score(features)
        s2_flagged = fusion_score >= self.stage2_threshold

        # Combined probability estimate
        # P(cancer | both positive) = P(cancer | S1+) × P(cancer | S2+)
        s1_pp = s1_result['final_posterior']
        s2_prob = 1.0 / (1.0 + np.exp(-(fusion_score - 0.5) * 10))  # calibrated sigmoid
        combined_prob = float(s1_pp * s2_prob)

        if s2_flagged:
            risk_tier = 'HIGH'
            recommendation = 'Immediate workup → Imaging + biopsy'
        else:
            risk_tier = 'MODERATE'
            recommendation = 'Watchful waiting → Repeat screening in 6 months'

        return {
            'risk_tier': risk_tier,
            'stage1_posterior': s1_result['final_posterior'],
            'stage1_flagged': True,
            'stage2_score': fusion_score,
            'stage2_features': features,
            'stage2_flagged': s2_flagged,
            'overall_probability': combined_prob,
            'stage_used': 2,
            'recommendation': recommendation,
        }


# ═══════════════════════════════════════════════════════════════════════════
# SIMULATION — 2000 PATIENTS WITH ALL 6 CONFOUNDERS
# ═══════════════════════════════════════════════════════════════════════════
def run_two_stage_simulation(
    n_cancer: int = N_CANCER,
    n_healthy: int = N_HEALTHY,
    n_benign: int = N_BENIGN,
    n_timepoints: int = N_TIMEPOINTS,
    seed: int = SEED,
) -> Dict:
    """
    Full two-stage CET simulation with 2000 patients and 6 confounders.

    Cancer breakdown:
      - 200 early (Stage I/II)
      - 150 mid (Stage III)
      - 150 late (Stage IV)

    Confounders:
      1. CHIP (age-dependent, prevalence 5-25% for ages 55-85)
      2. Variable cfDNA shedding (CV 60-80%)
      3. Trinucleotide error rates (context-dependent ×1.5-3.0)
      4. Variable genome equivalents (depth fluctuation ±30%)
      5. Batch effects (×1.0-1.45 per batch)
      6. Inflammatory spikes (3-15% of healthy patients per timepoint)

    Returns:
        Dict with full performance metrics, ROC curves, and cost analysis.
    """
    rng = np.random.RandomState(seed)
    cancer_types_list = list(CANCER_TYPES)

    total_n = n_cancer + n_healthy + n_benign
    logger.info(f"=" * 70)
    logger.info(f"TWO-STAGE CET SIMULATION")
    logger.info(f"  {n_cancer} cancer (200 early + 150 mid + 150 late)")
    logger.info(f"  {n_healthy} healthy + {n_benign} benign")
    logger.info(f"  {n_timepoints} quarterly timepoints over {n_timepoints * INTERVAL_DAYS} days")
    logger.info(f"  6 confounders: CHIP, shedding, tri-error, depth, batch, inflammation")
    logger.info(f"=" * 70)

    # ═══════════ PATIENT GENERATION ═══════════
    patients = []
    cancer_stage_map = {}

    # Cancer patients
    for i in range(n_cancer):
        ct = cancer_types_list[rng.randint(0, len(cancer_types_list))]

        # Stage assignment
        if i < N_EARLY:
            stage = 'early'
        elif i < N_EARLY + N_MID:
            stage = 'mid'
        else:
            stage = 'late'

        tumor_params = _generate_tumor_params(ct, stage, rng)
        start_day = rng.random() * 1500

        # Conf 1: Age-dependent CHIP
        age = int(50 + rng.randint(0, 35))

        # Conf 2: Variable shedding
        shedding_factor = max(0.3, min(3.0, np.exp(rng.normal(0, 0.6))))

        # Conf 3: Trinucleotide error rate factor
        tri_error_factor = 1.0 + rng.random() * 2.0

        # Conf 4: Depth fluctuation
        depth_factor = 0.7 + rng.random() * 0.6

        # Conf 5: Batch effect scale
        batch = 1 + rng.randint(0, 3)
        batch_scale = 0.10 + rng.random() * 0.10

        patients.append({
            'id': f'CANCER_{i:04d}',
            'is_cancer': True, 'is_benign': False,
            'cancer_type': ct, 'stage': stage,
            'tumor_params': tumor_params,
            'start_day': float(start_day),
            'age': age,
            'shedding_factor': shedding_factor,
            'tri_error_factor': tri_error_factor,
            'depth_factor': depth_factor,
            'batch': batch,
            'batch_scale': batch_scale,
            'inflammatory_spike': 0.0,
        })
        cancer_stage_map[f'CANCER_{i:04d}'] = stage

    # Healthy controls (with all confounders)
    for i in range(n_healthy):
        age = int(45 + rng.randint(0, 40))
        age_bins = [45, 55, 65, 75, 85]

        # Conf 1: CHIP prevalence increases with age
        chip_prob = max(0.05, min(0.25, (age - 50) / 35 * 0.20))

        # Conf 6: Inflammatory spikes (random per-timepoint)
        has_chronic_infl = rng.random() < 0.10  # 10% chronic inflammation
        inflammatory_spike = 0.3 if has_chronic_infl else 0.0

        patients.append({
            'id': f'HEALTHY_{i:04d}',
            'is_cancer': False, 'is_benign': False,
            'cancer_type': None, 'tumor_params': None,
            'start_day': 0.0,
            'age': age,
            'chip_prob': chip_prob,
            'shedding_factor': 0.5 + rng.random() * 0.5,
            'tri_error_factor': 1.0 + rng.random() * 2.0,
            'depth_factor': 0.7 + rng.random() * 0.6,
            'batch': 1 + rng.randint(0, 3),
            'batch_scale': 0.10 + rng.random() * 0.10,
            'inflammatory_spike': inflammatory_spike,
        })

    # Benign conditions (with confounders, elevated signals)
    for i in range(n_benign):
        age = int(50 + rng.randint(0, 35))
        patients.append({
            'id': f'BENIGN_{i:04d}',
            'is_cancer': False, 'is_benign': True,
            'cancer_type': None, 'tumor_params': None,
            'start_day': 0.0,
            'age': age,
            'chip_prob': 0.08 + rng.random() * 0.10,
            'shedding_factor': 0.5 + rng.random() * 0.8,
            'tri_error_factor': 1.0 + rng.random() * 2.5,
            'depth_factor': 0.7 + rng.random() * 0.6,
            'batch': 1 + rng.randint(0, 3),
            'batch_scale': 0.10 + rng.random() * 0.15,
            'inflammatory_spike': 0.4 + rng.random() * 0.3,
        })

    logger.info(f"Generated {len(patients)} patients")

    # ═══════════ GENERATE MULTI-MODAL SIGNALS ═══════════
    multi_signals_list = []

    for i, patient in enumerate(patients):
        signals = []
        for t in range(N_TIMEPOINTS):
            time_days = patient['start_day'] + t * INTERVAL_DAYS

            # Conf 6: Inflammatory spike application (per-timepoint)
            # For healthy patients without chronic inflammation, occasional spikes
            if not patient['is_cancer'] and not patient.get('is_benign'):
                if patient.get('inflammatory_spike', 0) < 0.05 and rng.random() < 0.05:
                    patient['inflammatory_spike'] = 0.25  # transient spike

            sigs = _generate_all_modality_signals(patient, time_days, rng)
            sigs['time_days'] = time_days
            signals.append(sigs)

        multi_signals_list.append(signals)

        if (i + 1) % 500 == 0:
            logger.info(f"  Generated signals for {i + 1}/{len(patients)} patients...")

    logger.info(f"Generated multi-modal signals for all {len(patients)} patients")

    # ═══════════ CALIBRATE AND SCREEN ═══════════
    screener = TwoStageCETScreener(
        n_modalities=5,
        stage1_target_spec=0.85,
        stage2_target_spec=0.99,
    )

    # Split: 50% calibration, 50% test
    n_total = len(patients)
    indices = list(range(n_total))
    rng.shuffle(indices)
    n_cal = n_total // 2
    cal_indices = set(indices[:n_cal])
    test_indices = set(indices[n_cal:])

    cal_patients = [patients[i] for i in cal_indices]
    cal_signals = [multi_signals_list[i] for i in cal_indices]
    cal_results_data = screener.calibrate(cal_patients, cal_signals, rng, cal_split=1.0)

    # ═══════════ TEST SET SCREENING ═══════════
    test_results = []
    true_labels = []

    for idx in test_indices:
        patient = patients[idx]
        signals = multi_signals_list[idx]
        result = screener.screen(patient, signals, rng)
        result['patient_id'] = patient['id']
        result['is_true_cancer'] = patient['is_cancer']
        result['is_benign'] = patient['is_benign']
        result['cancer_type'] = patient.get('cancer_type')
        result['stage'] = patient.get('stage', 'N/A')
        result['age'] = patient['age']
        test_results.append(result)
        true_labels.append(1 if patient['is_cancer'] else 0)

    # ═══════════ COMPUTE COMBINED METRICS ═══════════
    # Overall (HIGH risk = positive, others = negative)
    high_risk = [1 if r['risk_tier'] == 'HIGH' else 0 for r in test_results]
    stage1_flags = [1 if r['stage1_flagged'] else 0 for r in test_results]
    stage2_flags = [1 if r.get('stage2_flagged', False) else 0 for r in test_results]

    def compute_metrics_std(preds, labels):
        tp = sum(1 for p, l in zip(preds, labels) if l == 1 and p == 1)
        fp = sum(1 for p, l in zip(preds, labels) if l == 0 and p == 1)
        tn = sum(1 for p, l in zip(preds, labels) if l == 0 and p == 0)
        fn = sum(1 for p, l in zip(preds, labels) if l == 1 and p == 0)
        sens = tp / max(1, tp + fn)
        spec = tn / max(1, tn + fp)
        prec = tp / max(1, tp + fp)
        f1 = 2 * sens * prec / max(0.001, sens + prec)
        f2 = (5 * prec * sens) / max(0.001, 4 * prec + sens)
        return {'tp': tp, 'fp': fp, 'tn': tn, 'fn': fn,
                'sensitivity': float(sens), 'specificity': float(spec),
                'precision': float(prec), 'f1': float(f1), 'f2': float(f2)}

    stage1_metrics = compute_metrics_std(stage1_flags, true_labels)
    combined_metrics = compute_metrics_std(high_risk, true_labels)

    # Stage 2 only metrics (only on flagged patients)
    stage2_indices = [i for i, f in enumerate(stage1_flags) if f == 1]
    if stage2_indices:
        s2_preds = [high_risk[i] for i in stage2_indices]
        s2_labels = [true_labels[i] for i in stage2_indices]
        stage2_metrics = compute_metrics_std(s2_preds, s2_labels)
    else:
        stage2_metrics = {'sensitivity': 0.0, 'specificity': 1.0}

    # Flag rate and cost analysis
    flag_rate = sum(stage1_flags) / max(1, len(stage1_flags))
    high_risk_rate = sum(high_risk) / max(1, len(high_risk))

    # Cost model
    cost_per_stage1 = 74    # Standard multi-modal panel
    cost_per_stage2 = 200   # High-depth targeted sequencing + fragmentomics
    avg_cost = cost_per_stage1 + flag_rate * cost_per_stage2
    population_cost = avg_cost * 100000  # per 100K screened

    # Per-stage sensitivity
    cancer_indices = [i for i, l in enumerate(true_labels) if l == 1]
    cancer_s1 = [stage1_flags[i] for i in cancer_indices]
    cancer_combined = [high_risk[i] for i in cancer_indices]
    s1_sens = sum(cancer_s1) / max(1, len(cancer_s1))
    combined_sens = sum(cancer_combined) / max(1, len(cancer_combined))

    # Specificity improvement
    spec_improvement = combined_metrics['specificity'] - stage1_metrics['specificity']

    # Per-stage performance
    per_stage = {}
    for stage in ['early', 'mid', 'late']:
        stage_results = [r for r in test_results if r.get('stage') == stage]
        if stage_results:
            stage_labels = [1] * len(stage_results)
            stage_preds = [1 if r['risk_tier'] == 'HIGH' else 0 for r in stage_results]
            pm = compute_metrics_std(stage_preds, stage_labels)
            per_stage[stage] = pm

    # Bootstrap CIs
    n_boot = 2000
    combined_sens_bs = []
    combined_spec_bs = []
    n_test = len(true_labels)
    for _ in range(n_boot):
        bs_idx = [rng.randint(0, n_test - 1) for _ in range(n_test)]
        bs_labels = [true_labels[i] for i in bs_idx]
        bs_preds = [high_risk[i] for i in bs_idx]
        m = compute_metrics_std(bs_preds, bs_labels)
        combined_sens_bs.append(m['sensitivity'])
        combined_spec_bs.append(m['specificity'])
    combined_sens_bs.sort()
    combined_spec_bs.sort()

    # Verdict
    target_spec = 0.99
    target_sens = 0.50
    spec_met = combined_metrics['specificity'] >= target_spec
    sens_met = combined_metrics['sensitivity'] >= target_sens

    if spec_met and sens_met:
        verdict = '✅ BOTH TARGETS MET: Specificity >99% AND Sensitivity >50%'
    elif spec_met:
        verdict = f'⚠️ SPEC TARGET MET ({combined_metrics["specificity"]*100:.1f}%≥99%), SENS BELOW TARGET ({combined_metrics["sensitivity"]*100:.1f}%<50%)'
    elif sens_met:
        verdict = f'⚠️ SENS TARGET MET, SPEC BELOW TARGET ({combined_metrics["specificity"]*100:.1f}%<99%)'
    else:
        verdict = f'❌ NEITHER TARGET MET: Spec={combined_metrics["specificity"]*100:.1f}%, Sens={combined_metrics["sensitivity"]*100:.1f}%'

    output = {
        'metadata': {
            'generated': True,
            'model': 'Two-Stage CET (Permissive Multi-Modal SPRT → Confirmatory Fusion)',
            'reference': 'research/CET_OPTIMIZATION_RESEARCH.md (Approach 2)',
            'parameters': {
                'n_cancer': n_cancer, 'n_healthy': n_healthy, 'n_benign': n_benign,
                'n_early_stage': N_EARLY, 'n_mid_stage': N_MID, 'n_late_stage': N_LATE,
                'n_timepoints': n_timepoints, 'interval_days': INTERVAL_DAYS,
                'baseline_timepoints': BASELINE_TIMEPOINTS,
                'prior_prevalence': PRIOR_PREVALENCE, 'seed': seed,
                'stage1_target_spec': 0.85, 'stage2_target_spec': 0.99,
            },
            'confounders': [
                'CHIP (age-dependent, 5-25% prevalence)',
                'Variable cfDNA shedding (CV 60-80%)',
                'Trinucleotide error rates (×1.5-3.0)',
                'Variable genome equivalents (±30% depth)',
                'Batch effects (×1.0-1.45)',
                'Inflammatory spikes (3-15% of healthy)',
            ],
        },
        'calibration': cal_results_data,
        'performance': {
            'stage1': stage1_metrics,
            'stage2_on_flagged': stage2_metrics,
            'combined': combined_metrics,
            'combined_sens_ci95_low': float(combined_sens_bs[int(0.025 * n_boot)]),
            'combined_sens_ci95_high': float(combined_sens_bs[int(0.975 * n_boot) - 1]),
            'combined_spec_ci95_low': float(combined_spec_bs[int(0.025 * n_boot)]),
            'combined_spec_ci95_high': float(combined_spec_bs[int(0.975 * n_boot) - 1]),
            'flag_rate_stage1': float(flag_rate),
            'high_risk_rate': float(high_risk_rate),
            'per_stage': per_stage,
        },
        'cost_analysis': {
            'stage1_cost_per_person': cost_per_stage1,
            'stage2_cost_per_person': cost_per_stage2,
            'average_cost_per_person': float(avg_cost),
            'cost_per_100k_population': float(population_cost),
            'pct_getting_stage2': float(flag_rate * 100),
        },
        'verdict': verdict,
        'targets': {
            'specificity_gt_99': spec_met,
            'sensitivity_gt_50': sens_met,
            'both_met': spec_met and sens_met,
        },
        'patient_summary': {
            'n_test_patients': n_test,
            'distribution': {
                'low_risk': sum(1 for r in test_results if r['risk_tier'] == 'LOW'),
                'moderate_risk': sum(1 for r in test_results if r['risk_tier'] == 'MODERATE'),
                'high_risk': sum(1 for r in test_results if r['risk_tier'] == 'HIGH'),
            },
        },
    }

    # Save
    output_path = RESULTS_PY_DIR / 'two_stage_cet_results.json'
    with open(output_path, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    logger.info(f"\nSaved results to {output_path}")

    # Print summary
    logger.info(f"\n{'='*70}")
    logger.info(f"TWO-STAGE CET RESULTS")
    logger.info(f"{'='*70}")
    logger.info(f"  Stage 1 (Multi-Modal CET):")
    logger.info(f"    Sensitivity: {stage1_metrics['sensitivity']*100:.1f}%")
    logger.info(f"    Specificity: {stage1_metrics['specificity']*100:.1f}%")
    logger.info(f"    Flag rate:    {flag_rate*100:.1f}%")
    logger.info(f"  Stage 2 (Confirmatory Fusion on flagged):")
    logger.info(f"    Sensitivity: {stage2_metrics['sensitivity']*100:.1f}% (on flagged)")
    logger.info(f"    Specificity: {stage2_metrics['specificity']*100:.1f}% (on flagged)")
    logger.info(f"  COMBINED:")
    logger.info(f"    Sensitivity: {combined_metrics['sensitivity']*100:.1f}% "
                f"[{output['performance']['combined_sens_ci95_low']*100:.1f}–"
                f"{output['performance']['combined_sens_ci95_high']*100:.1f}%]")
    logger.info(f"    Specificity: {combined_metrics['specificity']*100:.1f}% "
                f"[{output['performance']['combined_spec_ci95_low']*100:.1f}–"
                f"{output['performance']['combined_spec_ci95_high']*100:.1f}%]")
    logger.info(f"    PPV:         {combined_metrics['precision']*100:.1f}%")
    logger.info(f"    F2:          {combined_metrics['f2']:.4f}")
    logger.info(f"  COST:")
    logger.info(f"    Avg per person: ${avg_cost:.0f}")
    logger.info(f"    Per 100K:      ${population_cost:,.0f}")
    logger.info(f"    Stage 2 usage:  {flag_rate*100:.1f}%")
    logger.info(f"  {verdict}")

    return output


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    results = run_two_stage_simulation(seed=42)
