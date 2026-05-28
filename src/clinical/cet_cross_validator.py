#!/usr/bin/env python3
"""
Nested Cross-Validation CET Validator (P0-A1)

Implements a rigorous nested cross-validation framework for motif-based
cancer early test (CET) feature selection. Addresses the pre-filter leakage
problem in the current pipeline where Mann-Whitney U p-values are computed
on the ENTIRE dataset before cross-validation — allowing the classifier to
indirectly observe test-fold information during feature selection.

Architecture
------------

Outer CV (5-fold stratified)
├── Outer Fold 1
│   ├── Inner CV (3-fold) on trainₒᵤₜₑᵣ
│   │   ├── Feature selection on trainᵢₙₙₑᵣ
│   │   ├── Train LogisticRegression
│   │   └── Validate on testᵢₙₙₑᵣ
│   └── Select top-k motifs (from inner CV consensus)
│       Train LR on trainₒᵤₜₑᵣ → Predict testₒᵤₜₑᵣ
├── ... (repeat for all outer folds)
└── Aggregate: mean AUC, motif selection frequency, stability

This ensures the feature selection NEVER sees test-fold data, providing
an unbiased estimate of generalization performance.

References
----------
- Cawley & Talbot (2010): "On Over-fitting in Model Selection and
  Subsequent Selection Bias in Performance Evaluation"
- Varma & Simon (2006): "Bias in error estimation when using
  cross-validation for model selection"
- Krstajic et al. (2014): "Cross-validation pitfalls when selecting
  and assessing regression and classification models"
"""

from __future__ import annotations

import json
import logging
import sys
from collections import Counter
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

logger = logging.getLogger(__name__)


# ── Constants ─────────────────────────────────────────────────────────────
_DEFAULT_OUTER_FOLDS = 5
_DEFAULT_INNER_FOLDS = 3
_DEFAULT_K_MOTIFS = 50
_DEFAULT_RANDOM_STATE = 42

# ── Try to import config constants; fall back gracefully ──────────────────
try:
    from validation.py.config import SEED, CANCER_TYPES, N_FOLDS, CTDNA_LEVELS
except ImportError:
    SEED = 42
    CANCER_TYPES = ['LUAD', 'COADREAD', 'BRCA', 'PRAD',
                    'STAD', 'LIHC', 'PAAD', 'OV', 'BLCA', 'HNSC']
    N_FOLDS = 5
    CTDNA_LEVELS = [0.01, 0.005, 0.0025, 0.001, 0.0005,
                    0.00025, 0.0001, 0.00005, 0.00001]


# ═══════════════════════════════════════════════════════════════════════════
# DATA STRUCTURES
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class MotifRanking:
    """
    Result of ranking a single motif for differential analysis.

    Attributes
    ----------
    motif_idx : int
        Zero-based motif index in the original feature matrix.
    p_value : float
        Mann-Whitney U two-sided p-value.
    cliffs_delta : float
        Non-parametric effect size in [-1, 1].
    composite_score : float
        Combined score ranking motifs (higher = more discriminative).
        Computed as: -log10(p_value + 1e-12) × |cliffs_delta|.
    adj_p_value : float or None
        Benjamini-Hochberg FDR-corrected p-value (None if FDR not applied).
    """
    motif_idx: int
    p_value: float
    cliffs_delta: float
    composite_score: float
    adj_p_value: Optional[float] = None

    def to_dict(self) -> Dict:
        """Serialize to plain dict for JSON output."""
        return {
            'motif_idx': self.motif_idx,
            'p_value': float(self.p_value),
            'cliffs_delta': float(self.cliffs_delta),
            'composite_score': float(self.composite_score),
            'adj_p_value': float(self.adj_p_value) if self.adj_p_value is not None else None,
        }


@dataclass
class EnrichmentProfile:
    """
    Per-fold diagnostics for a single outer cross-validation fold.

    Attributes
    ----------
    fold : int
        Zero-based outer fold index.
    n_train : int
        Number of training samples.
    n_test : int
        Number of test samples.
    train_auc : float
        AUC on training set (monitor for overfitting).
    test_auc : float
        AUC on the held-out test fold.
    selected_motif_indices : list of int
        Indices of motifs selected in this fold.
    best_k : int or None
        Optimal k from inner CV (if nested); None for non-nested.
    inner_aucs : list of float or None
        Inner CV AUCs across inner folds (if nested).
    """
    fold: int
    n_train: int
    n_test: int
    train_auc: float
    test_auc: float
    selected_motif_indices: List[int] = field(default_factory=list)
    best_k: Optional[int] = None
    inner_aucs: Optional[List[float]] = None

    def to_dict(self) -> Dict:
        """Serialize to plain dict."""
        return {
            'fold': self.fold,
            'n_train': self.n_train,
            'n_test': self.n_test,
            'train_auc': float(self.train_auc),
            'test_auc': float(self.test_auc),
            'n_selected_motifs': len(self.selected_motif_indices),
            'selected_motif_indices': [int(m) for m in self.selected_motif_indices],
            'best_k': self.best_k,
            'inner_aucs': [float(a) for a in self.inner_aucs] if self.inner_aucs else None,
        }


# ═══════════════════════════════════════════════════════════════════════════
# EFFECT SIZE: CLIFF'S DELTA
# ═══════════════════════════════════════════════════════════════════════════

def compute_cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    """
    Compute Cliff's delta — a non-parametric effect size measure.

    Cliff's delta quantifies the degree of overlap between two
    distributions without assuming normality.  It is the probability
    that a randomly chosen observation from *x* exceeds a randomly
    chosen observation from *y*, minus the reverse probability.

    .. math::

        delta = P(x_i > y_j) - P(x_i < y_j)

    Interpretation (Romano et al. 2006):

    - |delta| < 0.147: negligible
    - 0.147 ≤ |delta| < 0.33: small
    - 0.33 ≤ |delta| < 0.474: medium
    - |delta| ≥ 0.474: large

    Parameters
    ----------
    x : np.ndarray (1-D)
        Sample from first distribution (e.g., cancer group).
    y : np.ndarray (1-D)
        Sample from second distribution (e.g., healthy group).

    Returns
    -------
    float
        Cliff's delta value in [-1, 1].

    Notes
    -----
    Implementation uses the rank-biserial formulation for O(n log n)
    complexity by sorting the concatenated arrays.

    References
    ----------
    Cliff, N. (1993). Dominance statistics: Ordinal analyses to answer
    ordinal questions. *Psychological Bulletin*, 114(3), 494–509.
    """
    x = np.asarray(x, dtype=np.float64).ravel()
    y = np.asarray(y, dtype=np.float64).ravel()

    nx, ny = len(x), len(y)
    if nx == 0 or ny == 0:
        return 0.0

    # Count how many y-values each x exceeds, and vice-versa
    # Using sorting for efficiency: O((nx+ny) log(nx+ny))
    combined = np.concatenate([x, y])
    group = np.concatenate([np.ones(nx), np.zeros(ny)])

    order = np.argsort(combined)
    sorted_group = group[order]

    # For each x, count y's with smaller value
    cumulative_y = np.cumsum(1 - sorted_group)

    # Dominance: x > y count minus x < y count
    x_positions = np.where(sorted_group == 1)[0]
    if len(x_positions) == 0:
        return 0.0

    x_gt_y = np.sum(cumulative_y[x_positions])
    x_lt_y = nx * ny - x_gt_y

    delta = (x_gt_y - x_lt_y) / (nx * ny)
    return float(delta)


# ═══════════════════════════════════════════════════════════════════════════
# MOTIF RANKER
# ═══════════════════════════════════════════════════════════════════════════

class MotifRanker:
    """
    Rank motifs by their discriminative power between two groups.

    Supports Mann-Whitney U as the primary statistical test with
    Cliff's delta as a complementary effect-size measure.  A composite
    score combines both into a single ranking criterion.

    .. rubric:: Composite Score

    .. math::

        score_j = -\\log_{10}(p_j + 10^{-12}) \\times |\\delta_j|

    where :math:`p_j` is the Mann-Whitney p-value and :math:`\\delta_j`
    is Cliff's delta for motif *j*.  This balances statistical significance
    with biological relevance (effect size).

    Parameters
    ----------
    method : {'mannwhitney', 'cliffs_delta'}
        Primary ranking method.  ``'mannwhitney'`` (default) uses the
        composite score; ``'cliffs_delta'`` ranks by absolute effect size.
    apply_fdr : bool
        Whether to apply Benjamini-Hochberg FDR correction (default True).

    Attributes
    ----------
    method : str
        Active ranking method.
    apply_fdr : bool
        FDR correction toggle.
    """

    def __init__(self, method: str = 'mannwhitney', apply_fdr: bool = True):
        if method not in ('mannwhitney', 'cliffs_delta'):
            raise ValueError(f"Unknown method '{method}'. "
                             f"Use 'mannwhitney' or 'cliffs_delta'.")
        self.method = method
        self.apply_fdr = apply_fdr

    def rank(self, X: np.ndarray, y: np.ndarray) -> List[MotifRanking]:
        """
        Rank all motifs by their discriminative power.

        For each motif (column of *X*), computes:
        1. Mann-Whitney U two-sided p-value (cancer vs. healthy)
        2. Cliff's delta effect size
        3. Composite score = -log10(p + ε) × |delta|
        4. (Optionally) Benjamini-Hochberg adjusted p-value

        Parameters
        ----------
        X : np.ndarray of shape (n_samples, n_motifs)
            Motif frequency matrix.  Each column is one motif's abundance
            across samples.
        y : np.ndarray of shape (n_samples,)
            Binary labels: 1 = cancer, 0 = healthy/control.

        Returns
        -------
        list of MotifRanking
            Rankings for all motifs, sorted by composite_score descending.
        """
        from scipy.stats import mannwhitneyu

        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64).ravel()

        # Check for NaN values and impute with column mean
        nan_count = int(np.sum(np.isnan(X)))
        if nan_count > 0:
            logger.warning("Input X contains %d NaN values. Replacing with column mean.", nan_count)
            col_means = np.nanmean(X, axis=0)
            nan_mask = np.isnan(X)
            X = np.where(nan_mask, col_means[np.newaxis, :], X)

        n_samples, n_motifs = X.shape
        cancer_idx = np.where(y == 1)[0]
        control_idx = np.where(y == 0)[0]

        if len(cancer_idx) == 0 or len(control_idx) == 0:
            raise ValueError("Both classes must have at least one sample "
                             f"(got {len(cancer_idx)} cancer, {len(control_idx)} control).")

        rankings: List[MotifRanking] = []

        for j in range(n_motifs):
            col = X[:, j]

            # Skip constant columns
            if np.std(col) < 1e-12:
                rankings.append(MotifRanking(
                    motif_idx=j,
                    p_value=1.0,
                    cliffs_delta=0.0,
                    composite_score=0.0,
                ))
                continue

            x_cancer = col[cancer_idx]
            x_control = col[control_idx]

            # Mann-Whitney U (two-sided)
            try:
                u_stat, p_val = mannwhitneyu(x_cancer, x_control, alternative='two-sided')
            except ValueError:
                p_val = 1.0

            # Cliff's delta
            delta = compute_cliffs_delta(x_cancer, x_control)

            # Composite score
            eps = 1e-12
            composite = -np.log10(max(p_val, eps)) * abs(delta)

            rankings.append(MotifRanking(
                motif_idx=j,
                p_value=float(p_val),
                cliffs_delta=float(delta),
                composite_score=float(composite),
            ))

        # Apply FDR correction if requested
        if self.apply_fdr:
            self._apply_benjamini_hochberg(rankings)

        # Sort by composite score descending
        rankings.sort(key=lambda r: r.composite_score, reverse=True)

        return rankings

    @staticmethod
    def _apply_benjamini_hochberg(rankings: List[MotifRanking]) -> None:
        """
        Apply Benjamini-Hochberg FDR correction to p-values in-place.

        Procedure:
        1. Sort p-values ascending.
        2. For each rank *i* (1-indexed), compute critical value:
           :math:`q_i = (i / m) \\times \\alpha` where :math:`\\alpha = 0.05`.
        3. Find the largest *i* where :math:`p_{(i)} \\leq q_i`.
        4. All motifs up to that rank are FDR-significant.

        Parameters
        ----------
        rankings : list of MotifRanking
            Rankings to adjust in-place.  After this call,
            ``adj_p_value`` is set for every entry.
        """
        m = len(rankings)
        if m == 0:
            return

        # Sort by p-value ascending
        sorted_by_p = sorted(rankings, key=lambda r: r.p_value)

        # BH procedure
        alpha = 0.05
        p_vals = np.array([r.p_value for r in sorted_by_p])
        ranks = np.arange(1, m + 1, dtype=np.float64)
        bh_critical = ranks / m * alpha

        # Adjusted p-values: p_raw * m / rank  (capped at 1.0)
        # But enforce monotonicity: adj_p[i] = min(adj_p[i], adj_p[i+1], ...)
        adj_p_raw = np.minimum(p_vals * m / ranks, 1.0)

        # Cumulative min from the right (monotonicity enforcement)
        for i in range(m - 2, -1, -1):
            adj_p_raw[i] = min(adj_p_raw[i], adj_p_raw[i + 1])

        # Assign back
        for r, adj_p in zip(sorted_by_p, adj_p_raw):
            r.adj_p_value = float(adj_p)

    def select_top_k(self, rankings: List[MotifRanking], k: int) -> List[int]:
        """
        Return the indices of the top *k* motifs from a ranking.

        Parameters
        ----------
        rankings : list of MotifRanking
            Ranked motifs (assumed sorted by composite_score descending).
        k : int
            Number of top motifs to select.

        Returns
        -------
        list of int
            Indices of top-k motifs.
        """
        return [r.motif_idx for r in rankings[:k]]


# ═══════════════════════════════════════════════════════════════════════════
# NESTED CET VALIDATOR
# ═══════════════════════════════════════════════════════════════════════════

class NestedCETValidator:
    """
    Nested cross-validation validator for motif-based CET.

    Wraps the feature-selection-and-classification pipeline in a
    doubly-nested CV framework to prevent information leakage:
    feature selection (motif ranking) is performed *only on training
    data* within each outer fold.  The test fold is completely held
    out from ranking, selection, and model fitting.

    .. rubric:: Nested CV Procedure

    1. Outer split: 5-fold stratified CV.
    2. For each outer fold:
       a. **Inner CV** (3-fold) on outer training fold:
          - Split outer-train into inner-train / inner-val.
          - Rank motifs on inner-train only.
          - Train LogisticRegression on top-k inner-train.
          - Evaluate on inner-val.
          - Average inner AUC across inner-folds for each k candidate.
       b. Select optimal *k* (highest mean inner AUC).
       c. Rank motifs on full outer-training fold.
       d. Train LogisticRegression on top-*k* outer-training features.
       e. Predict and evaluate on the held-out outer test fold.
    3. Aggregate outer-fold AUCs.

    Parameters
    ----------
    outer_folds : int
        Number of outer CV folds (default 5).
    inner_folds : int
        Number of inner CV folds for hyperparameter / k selection (default 3).
    k_motifs : int or None
        Number of top motifs to select.  If None, uses the inner CV
        to select the optimal k from ``k_range``.
    k_range : range or list of int
        Candidate k values for inner CV optimization (only used when
        ``k_motifs is None``).
    random_state : int
        Seed for reproducibility.

    Attributes
    ----------
    outer_folds : int
    inner_folds : int
    k_motifs : int or None
    random_state : int
    feature_importance_ranks : list of int or None
        Motif indices ranked by selection frequency across outer folds.
    """

    def __init__(
        self,
        outer_folds: int = _DEFAULT_OUTER_FOLDS,
        inner_folds: int = _DEFAULT_INNER_FOLDS,
        k_motifs: int = _DEFAULT_K_MOTIFS,
        k_range: Optional[Union[range, List[int]]] = None,
        random_state: int = _DEFAULT_RANDOM_STATE,
    ):
        self.outer_folds = outer_folds
        self.inner_folds = inner_folds
        self.k_motifs = k_motifs
        self.k_range = k_range if k_range is not None else range(5, 200, 5)
        if isinstance(self.k_range, range):
            self.k_range = list(self.k_range)
        self.random_state = random_state
        self.feature_importance_ranks: Optional[List[int]] = None

        logger.debug(
            "NestedCETValidator initialized: outer=%d, inner=%d, k_motifs=%s, seed=%d",
            outer_folds, inner_folds, k_motifs, random_state,
        )

    def validate(self, X: np.ndarray, y: np.ndarray) -> Dict:
        """
        Run nested cross-validation and return aggregated results.

        Parameters
        ----------
        X : np.ndarray of shape (n_samples, n_motifs)
            Motif frequency feature matrix.
        y : np.ndarray of shape (n_samples,)
            Binary labels (1 = cancer, 0 = control).

        Returns
        -------
        dict
            Aggregated results with keys:

            - ``outer_aucs``: AUC per outer fold (list of float)
            - ``mean_auc``: mean AUC across outer folds
            - ``std_auc``: standard deviation of outer AUCs
            - ``feature_importance_ranks``: motifs ranked by selection
              frequency (most frequently selected first)
            - ``per_fold_results``: list of ``EnrichmentProfile`` dicts
              with per-fold diagnostics
            - ``n_motifs``: total number of motifs in X
            - ``n_samples``: total number of samples
        """
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64).ravel()

        n_samples, n_motifs = X.shape

        outer_skf = StratifiedKFold(
            n_splits=self.outer_folds,
            shuffle=True,
            random_state=self.random_state,
        )

        outer_aucs: List[float] = []
        per_fold_results: List[EnrichmentProfile] = []
        selection_counter: Counter = Counter()

        use_inner_cv = self.k_motifs is None

        for fold_idx, (train_idx, test_idx) in enumerate(outer_skf.split(X, y)):
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # ── Determine k (fixed or via inner CV) ──
            if use_inner_cv:
                best_k = self._inner_cv_select_k(X_train, y_train)
            else:
                best_k = self.k_motifs

            # ── Rank motifs on outer training data only ──
            ranker = MotifRanker(method='mannwhitney', apply_fdr=True)
            rankings = ranker.rank(X_train, y_train)
            selected_indices = ranker.select_top_k(rankings, best_k)

            # Track selection frequency
            for mi in selected_indices:
                selection_counter[mi] += 1

            # ── Train classifier on selected features ──
            X_train_sel = X_train[:, selected_indices]
            X_test_sel = X_test[:, selected_indices]

            clf = LogisticRegression(
                C=1.0,
                solver='liblinear',
                random_state=self.random_state + fold_idx,
                max_iter=5000,
            )
            clf.fit(X_train_sel, y_train)

            # ── Predict and evaluate ──
            train_probs = clf.predict_proba(X_train_sel)[:, 1]
            test_probs = clf.predict_proba(X_test_sel)[:, 1]

            train_auc = roc_auc_score(y_train, train_probs)
            test_auc = roc_auc_score(y_test, test_probs)

            outer_aucs.append(float(test_auc))

            pf = EnrichmentProfile(
                fold=fold_idx,
                n_train=len(train_idx),
                n_test=len(test_idx),
                train_auc=float(train_auc),
                test_auc=float(test_auc),
                selected_motif_indices=selected_indices,
                best_k=best_k,
                inner_aucs=None,  # inner CV details stored if needed
            )
            per_fold_results.append(pf)

            logger.info(
                "Fold %d/%d | train AUC=%.4f | test AUC=%.4f | n_selected=%d | k=%d",
                fold_idx + 1, self.outer_folds, train_auc, test_auc,
                len(selected_indices), best_k,
            )

        # ── Feature importance: motifs ranked by selection frequency ──
        ranked_motifs = [
            motif_idx for motif_idx, _ in
            selection_counter.most_common()
        ]
        self.feature_importance_ranks = ranked_motifs

        mean_auc = float(np.mean(outer_aucs))
        std_auc = float(np.std(outer_aucs, ddof=1))

        logger.info(
            "Nested CV complete: mean AUC=%.4f ± %.4f (across %d folds)",
            mean_auc, std_auc, self.outer_folds,
        )

        return {
            'outer_aucs': [float(a) for a in outer_aucs],
            'mean_auc': mean_auc,
            'std_auc': std_auc,
            'feature_importance_ranks': ranked_motifs,
            'per_fold_results': [pf.to_dict() for pf in per_fold_results],
            'n_motifs': n_motifs,
            'n_samples': n_samples,
        }

    def _inner_cv_select_k(self, X_train: np.ndarray, y_train: np.ndarray) -> int:
        """
        Select optimal k via inner cross-validation.

        For each candidate k in ``self.k_range``:
        1. Split ``X_train`` into inner folds.
        2. For each inner fold: rank on inner-train, train on top-k, eval.
        3. Average inner AUCs across folds.
        4. Pick k with highest mean inner AUC.

        Parameters
        ----------
        X_train : np.ndarray
            Outer-fold training data.
        y_train : np.ndarray
            Outer-fold training labels.

        Returns
        -------
        int
            Optimal k (number of motifs to select).
        """
        inner_skf = StratifiedKFold(
            n_splits=self.inner_folds,
            shuffle=True,
            random_state=self.random_state,
        )

        best_k = self.k_range[0]
        best_mean_auc = 0.0
        k_scores: Dict[int, float] = {}

        for k in self.k_range:
            inner_aucs: List[float] = []

            for inner_train_idx, inner_val_idx in inner_skf.split(X_train, y_train):
                X_it = X_train[inner_train_idx]
                X_iv = X_train[inner_val_idx]
                y_it = y_train[inner_train_idx]
                y_iv = y_train[inner_val_idx]

                # Rank on inner train only
                ranker = MotifRanker(method='mannwhitney', apply_fdr=True)
                rankings = ranker.rank(X_it, y_it)
                selected = ranker.select_top_k(rankings, min(k, X_it.shape[1]))

                if not selected:
                    inner_aucs.append(0.5)
                    continue

                clf = LogisticRegression(
                    C=1.0, solver='liblinear',
                    max_iter=5000, random_state=self.random_state,
                )
                clf.fit(X_it[:, selected], y_it)
                probs = clf.predict_proba(X_iv[:, selected])[:, 1]
                inner_aucs.append(roc_auc_score(y_iv, probs))

            mean_auc = float(np.mean(inner_aucs))
            k_scores[k] = mean_auc

            if mean_auc > best_mean_auc:
                best_mean_auc = mean_auc
                best_k = k

        logger.debug("Inner CV: best_k=%d (AUC=%.4f) from %d candidates",
                     best_k, best_mean_auc, len(self.k_range))
        return best_k

    def compare_leakage(self, X: np.ndarray, y: np.ndarray) -> Dict:
        """
        Compare nested (unbiased) vs. non-nested (leaky) CV performance.

        **Nested**: feature selection inside each outer fold (correct,
        no leakage).  **Non-nested**: feature selection on the full
        dataset, then standard CV (pre-filter leakage, optimistic bias).

        The difference:

        .. math::

            \\Delta_{\\text{leak}} = \\text{AUC}_{\\text{non-nested}}
            - \\text{AUC}_{\\text{nested}}

        estimates the optimism from pre-filter leakage.  A positive delta
        means the current non-nested pipeline overestimates performance.

        Parameters
        ----------
        X : np.ndarray of shape (n_samples, n_motifs)
            Motif frequency feature matrix.
        y : np.ndarray of shape (n_samples,)
            Binary labels.

        Returns
        -------
        dict
            Comparison results with keys:

            - ``nested_mean_auc``: unbiased AUC from nested CV
            - ``nested_std_auc``: std of nested AUCs
            - ``non_nested_mean_auc``: biased AUC from non-nested CV
            - ``non_nested_std_auc``: std of non-nested AUCs
            - ``leakage_delta``: difference (non_nested − nested), positive = optimistic
            - ``leakage_pct``: leakage as percentage of nested AUC
            - ``nested_fold_aucs``: per-fold nested AUCs
            - ``non_nested_fold_aucs``: per-fold non-nested AUCs
        """
        logger.info("Running nested CV...")
        nested_results = self.validate(X, y)

        logger.info("Running non-nested (leaky) CV for comparison...")
        non_nested_results = self._run_non_nested_cv(X, y)

        nested_mean = nested_results['mean_auc']
        non_nested_mean = non_nested_results['mean_auc']
        leakage_delta = non_nested_mean - nested_mean
        leakage_pct = (leakage_delta / max(nested_mean, 0.001)) * 100.0

        logger.info(
            "Leakage comparison: nested=%.4f, non_nested=%.4f, delta=%.4f (%.1f%%)",
            nested_mean, non_nested_mean, leakage_delta, leakage_pct,
        )

        return {
            'nested_mean_auc': nested_mean,
            'nested_std_auc': nested_results['std_auc'],
            'non_nested_mean_auc': non_nested_mean,
            'non_nested_std_auc': non_nested_results['std_auc'],
            'leakage_delta': float(leakage_delta),
            'leakage_pct': float(leakage_pct),
            'nested_fold_aucs': nested_results['outer_aucs'],
            'non_nested_fold_aucs': non_nested_results['outer_aucs'],
        }

    def _run_non_nested_cv(self, X: np.ndarray, y: np.ndarray) -> Dict:
        """
        Run standard (non-nested) CV with pre-filtering — the "leaky" approach.

        Feature selection is performed once on the ENTIRE dataset, then
        a standard stratified K-fold CV is run using those pre-selected
        features.  This exposes the classifier to test-fold information
        during feature selection, producing an optimistic bias.

        Parameters
        ----------
        X : np.ndarray
        y : np.ndarray

        Returns
        -------
        dict
            With same structure as :meth:`validate` return.
        """
        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64).ravel()

        # ── Pre-filter: rank on entire dataset (THE LEAK) ──
        ranker = MotifRanker(method='mannwhitney', apply_fdr=True)
        rankings = ranker.rank(X, y)
        selected_indices = ranker.select_top_k(rankings, self.k_motifs)

        # ── Standard CV on pre-selected features ──
        skf = StratifiedKFold(
            n_splits=self.outer_folds,
            shuffle=True,
            random_state=self.random_state,
        )

        outer_aucs: List[float] = []
        per_fold_results: List[EnrichmentProfile] = []

        for fold_idx, (train_idx, test_idx) in enumerate(skf.split(X, y)):
            X_train = X[train_idx][:, selected_indices]
            X_test = X[test_idx][:, selected_indices]
            y_train, y_test = y[train_idx], y[test_idx]

            clf = LogisticRegression(
                C=1.0, solver='liblinear',
                max_iter=5000,
                random_state=self.random_state + fold_idx,
            )
            clf.fit(X_train, y_train)

            train_probs = clf.predict_proba(X_train)[:, 1]
            test_probs = clf.predict_proba(X_test)[:, 1]

            train_auc = roc_auc_score(y_train, train_probs)
            test_auc = roc_auc_score(y_test, test_probs)

            outer_aucs.append(float(test_auc))

            per_fold_results.append(EnrichmentProfile(
                fold=fold_idx,
                n_train=len(train_idx),
                n_test=len(test_idx),
                train_auc=float(train_auc),
                test_auc=float(test_auc),
                selected_motif_indices=selected_indices,
                best_k=self.k_motifs,
            ))

        mean_auc = float(np.mean(outer_aucs))
        std_auc = float(np.std(outer_aucs, ddof=1))

        return {
            'outer_aucs': [float(a) for a in outer_aucs],
            'mean_auc': mean_auc,
            'std_auc': std_auc,
            'per_fold_results': [pf.to_dict() for pf in per_fold_results],
        }

    def find_optimal_k(
        self,
        X: np.ndarray,
        y: np.ndarray,
        k_range: Optional[Union[range, List[int]]] = None,
    ) -> Dict:
        """
        Knee-elbow analysis to find the optimal number of motifs (k).

        Performs nested CV for each candidate k in *k_range*, then applies
        a knee-point detection algorithm on the mean-AUC vs. k curve.
        The optimal k balances model performance with parsimony.

        Knee detection method:
        1. Fit a line between the first and last point of the (k, AUC) curve.
        2. For each k, compute the orthogonal distance from that point to
           the connecting line.
        3. The k maximizing this distance is the knee/elbow point.

        Parameters
        ----------
        X : np.ndarray of shape (n_samples, n_motifs)
        y : np.ndarray of shape (n_samples,)
        k_range : range or list of int, optional
            Candidate k values.  Defaults to ``range(5, 200, 5)``.

        Returns
        -------
        dict
            Analysis results with keys:

            - ``k_values``: list of k candidates tested
            - ``auc_means``: mean AUC for each k
            - ``auc_stds``: standard deviation of AUC for each k
            - ``optimal_k``: recommended k (knee point)
            - ``optimal_auc``: AUC at the optimal k
            - ``knee_distances``: orthogonal distance from line for each k
        """
        if k_range is None:
            k_range = list(range(5, 200, 5))
        if isinstance(k_range, range):
            k_range = list(k_range)

        X = np.asarray(X, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64).ravel()

        k_values: List[int] = []
        auc_means: List[float] = []
        auc_stds: List[float] = []

        logger.info(
            "find_optimal_k: testing %d candidates from k=%d to k=%d...",
            len(k_range), k_range[0], k_range[-1],
        )

        for k in k_range:
            validator = NestedCETValidator(
                outer_folds=self.outer_folds,
                inner_folds=self.inner_folds,
                k_motifs=k,
                random_state=self.random_state,
            )
            result = validator.validate(X, y)
            k_values.append(k)
            auc_means.append(result['mean_auc'])
            auc_stds.append(result['std_auc'])

            logger.info(
                "  k=%d: AUC=%.4f ± %.4f",
                k, result['mean_auc'], result['std_auc'],
            )

        # ── Knee-point detection ──
        # Normalize k and AUC to [0, 1] for distance computation
        k_arr = np.array(k_values, dtype=np.float64)
        auc_arr = np.array(auc_means, dtype=np.float64)

        k_min, k_max = k_arr[0], k_arr[-1]
        k_norm = (k_arr - k_min) / max(k_max - k_min, 1)
        auc_min, auc_max = auc_arr[0], auc_arr[-1]
        auc_norm = (auc_arr - auc_min) / max(auc_max - auc_min, 1e-12)

        # Line from (k_norm[0], auc_norm[0]) to (k_norm[-1], auc_norm[-1])
        p1 = np.array([k_norm[0], auc_norm[0]])
        p2 = np.array([k_norm[-1], auc_norm[-1]])
        line_vec = p2 - p1
        line_len_sq = np.dot(line_vec, line_vec)
        line_len_sq = max(line_len_sq, 1e-12)

        distances = []
        for i in range(len(k_arr)):
            p = np.array([k_norm[i], auc_norm[i]])
            # Project p onto line segment, clamp
            t = max(0.0, min(1.0, np.dot(p - p1, line_vec) / line_len_sq))
            proj = p1 + t * line_vec
            dist = np.linalg.norm(p - proj)
            distances.append(float(dist))

        knee_idx = int(np.argmax(distances))
        optimal_k = k_values[knee_idx]
        optimal_auc = auc_means[knee_idx]

        logger.info(
            "Knee-elbow optimal: k=%d (AUC=%.4f)",
            optimal_k, optimal_auc,
        )

        return {
            'k_values': k_values,
            'auc_means': auc_means,
            'auc_stds': auc_stds,
            'optimal_k': optimal_k,
            'optimal_auc': optimal_auc,
            'knee_distances': distances,
        }


# ═══════════════════════════════════════════════════════════════════════════
# SYNTHETIC DATA GENERATION (FOR DEMO)
# ═══════════════════════════════════════════════════════════════════════════

def _generate_synthetic_motif_data(
    n_motifs: int = 256,
    n_cancer: int = 50,
    n_control: int = 50,
    n_informative: int = 30,
    effect_strength: float = 0.8,
    random_state: int = SEED,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic motif frequency data for testing.

    Produces a motif×sample matrix with a known set of informative
    motifs that discriminate cancer from control.  Non-informative
    motifs are pure noise.

    Parameters
    ----------
    n_motifs : int
        Total number of motifs (columns).
    n_cancer : int
        Number of cancer samples.
    n_control : int
        Number of control (healthy) samples.
    n_informative : int
        Number of motifs with true discriminative signal.
    effect_strength : float
        Multiplier for cancer vs. control separation in informative motifs.
    random_state : int
        Seed for reproducibility.

    Returns
    -------
    X : np.ndarray of shape (n_cancer + n_control, n_motifs)
        Motif frequency matrix.
    y : np.ndarray of shape (n_cancer + n_control,)
        Binary labels (1 = cancer).
    """
    rng = np.random.RandomState(random_state)

    n_total = n_cancer + n_control
    X = np.zeros((n_total, n_motifs), dtype=np.float64)
    y = np.zeros(n_total, dtype=np.float64)
    y[:n_cancer] = 1.0

    # ── Background noise for all motifs ──
    # Motifs have different baseline frequencies (log-normal)
    baseline_mu = rng.lognormal(mean=0.0, sigma=0.6, size=n_motifs)
    baseline_cv = 0.15  # coefficient of variation

    for j in range(n_motifs):
        # All samples get baseline with poisson noise
        baseline = baseline_mu[j]
        X[:, j] = rng.poisson(lam=max(1, baseline), size=n_total).astype(np.float64)

        # Add biological variability
        X[:, j] *= np.exp(rng.normal(0, baseline_cv, n_total))

    # ── Informative motifs: elevated in cancer ──
    informative_indices = rng.choice(n_motifs, size=n_informative, replace=False)

    for j in informative_indices:
        # Cancer samples get boosted signal
        boost = rng.uniform(1.5, 1.5 + effect_strength * 3.0, n_cancer)
        X[:n_cancer, j] *= boost

        # Some motifs are also depressed in cancer (reverse effect)
        if rng.random() < 0.2:
            X[:n_cancer, j] = X[:n_cancer, j] / max(boost.mean(), 1.0)
            X[n_cancer:, j] *= rng.uniform(1.3, 3.0, n_control)

    # Add a few constant/noisy columns to test robustness
    # (handled by ranker which skips constant columns)

    return X, y


# ═══════════════════════════════════════════════════════════════════════════
# DEMO / MAIN
# ═══════════════════════════════════════════════════════════════════════════

def run_demo(random_state: int = SEED) -> Dict:
    """
    Run a standalone demonstration of the nested CV CET validator.

    Generates synthetic motif data (256 motifs, 100 samples with 30
    truly informative), then:

    1. Runs nested CV → prints mean ± std AUC
    2. Compares leakage vs. non-nested → prints delta
    3. Prints top 10 motifs by selection frequency
    4. Runs find_optimal_k and shows recommended k

    Parameters
    ----------
    random_state : int
        Random seed.

    Returns
    -------
    dict
        Demo results dictionary.
    """
    logger.info("=" * 70)
    logger.info("CET NESTED CROSS-VALIDATION DEMO")
    logger.info("=" * 70)

    # ── Generate synthetic data ──
    n_motifs = 256
    n_cancer = 50
    n_control = 50
    n_informative = 30

    logger.info(
        "Generating synthetic data: %d motifs, %d samples (%d cancer + %d control), "
        "%d truly informative motifs",
        n_motifs, n_cancer + n_control, n_cancer, n_control, n_informative,
    )

    X, y = _generate_synthetic_motif_data(
        n_motifs=n_motifs,
        n_cancer=n_cancer,
        n_control=n_control,
        n_informative=n_informative,
        random_state=random_state,
    )

    logger.info("  X shape: %s, y shape: %s, class balance: %.2f",
                X.shape, y.shape, y.mean())

    # ── 1. Nested CV ──
    logger.info("\n--- 1. Nested Cross-Validation ---")
    validator = NestedCETValidator(
        outer_folds=5,
        inner_folds=3,
        k_motifs=50,
        random_state=random_state,
    )
    nested_results = validator.validate(X, y)

    logger.info(
        "  Mean AUC: %.4f ± %.4f",
        nested_results['mean_auc'], nested_results['std_auc'],
    )
    logger.info("  Per-fold AUCs: %s",
                [f"{a:.4f}" for a in nested_results['outer_aucs']])

    # ── 2. Leakage comparison ──
    logger.info("\n--- 2. Leakage Comparison ---")
    leakage_results = validator.compare_leakage(X, y)

    logger.info(
        "  Nested AUC:     %.4f ± %.4f",
        leakage_results['nested_mean_auc'], leakage_results['nested_std_auc'],
    )
    logger.info(
        "  Non-nested AUC: %.4f ± %.4f",
        leakage_results['non_nested_mean_auc'], leakage_results['non_nested_std_auc'],
    )
    logger.info(
        "  Leakage delta:  %.4f (%.1f%% of nested AUC)",
        leakage_results['leakage_delta'], leakage_results['leakage_pct'],
    )

    if leakage_results['leakage_delta'] > 0.01:
        logger.info(
            "  ⚠️  Positive leakage detected — non-nested CV overestimates performance!"
        )
    else:
        logger.info("  ✅ Negligible leakage — nested and non-nested agree.")

    # ── 3. Top motifs ──
    logger.info("\n--- 3. Top Motifs by Selection Frequency ---")
    top_10 = nested_results['feature_importance_ranks'][:10]
    logger.info("  Rank | Motif ID | Selection Count (out of %d folds)",
                validator.outer_folds)

    selection_counts = Counter()
    for pf in nested_results['per_fold_results']:
        for mi in pf['selected_motif_indices']:
            selection_counts[mi] += 1

    for rank, motif_idx in enumerate(top_10, 1):
        count = selection_counts.get(motif_idx, 0)
        logger.info("  %4d | %8d | %d", rank, motif_idx, count)

    # ── 4. Optimal k ──
    logger.info("\n--- 4. Knee-Elbow Optimal k Analysis ---")
    k_results = validator.find_optimal_k(X, y, k_range=range(5, 200, 20))

    logger.info(
        "  Recommended k: %d (AUC=%.4f)",
        k_results['optimal_k'], k_results['optimal_auc'],
    )

    # Print k vs AUC table (abbreviated)
    logger.info("  k    | AUC")
    for k_val, auc_val in zip(k_results['k_values'], k_results['auc_means']):
        marker = " ← KNEE" if k_val == k_results['optimal_k'] else ""
        logger.info("  %3d  | %.4f%s", k_val, auc_val, marker)

    logger.info("\n" + "=" * 70)
    logger.info("DEMO COMPLETE")
    logger.info("=" * 70)

    return {
        'nested_cv': nested_results,
        'leakage_comparison': leakage_results,
        'top_10_motifs': top_10,
        'optimal_k_analysis': k_results,
    }


# ═══════════════════════════════════════════════════════════════════════════
# EXPORTS
# ═══════════════════════════════════════════════════════════════════════════

__all__ = [
    # Data structures
    "MotifRanking",
    "EnrichmentProfile",
    # Effect size
    "compute_cliffs_delta",
    # Feature ranking
    "MotifRanker",
    # Cross-validation
    "NestedCETValidator",
    # Demo
    "run_demo",
    "_generate_synthetic_motif_data",
]


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%H:%M:%S',
    )
    results = run_demo(random_state=SEED)
    # Also dump full JSON to stdout for inspection
    print("\n--- Full Results (JSON) ---")
    print(json.dumps({
        'nested_mean_auc': results['nested_cv']['mean_auc'],
        'nested_std_auc': results['nested_cv']['std_auc'],
        'leakage_delta': results['leakage_comparison']['leakage_delta'],
        'top_10_motifs': results['top_10_motifs'],
        'optimal_k': results['optimal_k_analysis']['optimal_k'],
        'optimal_k_auc': results['optimal_k_analysis']['optimal_auc'],
    }, indent=2))
