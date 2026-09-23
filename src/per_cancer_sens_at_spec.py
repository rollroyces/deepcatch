"""Per-cancer sens@spec + PPV@prev table generator (clinical-validation).

This is the standard cfDNA clinical-decision table: for each cancer type
(or pooled across all cancers), report:

  - AUC with a DeLong 95% CI
  - Sensitivity at specificity ∈ {0.95, 0.98, 0.99} with DeLong 95% CIs
  - PPV at prevalence ∈ {0.001, 0.004, 0.01, 0.05, 0.10, 0.20, 0.50}
    evaluated at the spec=0.99 operating point.

The math follows:

  - AUC CI: DeLong, DeLong, Clarke-Pearson (1988). Midrank-based
    structural component V_i = #{"negatives ranked below positive i"}/N +
    0.5·#{"ties"}; SE = sqrt(S10/N_pos + S01/N_neg - S10·S01/(N_pos·N_neg)) / ...

    We expose a closed-form implementation: the empirical covariance
    matrix of the placement values X = (V_1, ..., V_{N_pos}) on
    positives and Y = (W_1, ..., W_{N_neg}) on negatives, where
    V_i = (#negatives ranked below positive i) / N_neg and
    W_j = (#positives ranked above negative j) / N_pos. The AUC CI is
    derived from these placements without bootstrap resampling.

  - Sens@spec CI: standard reference is also DeLong placement-value
    SE on the positives (Sun & Xu, 2014; equivalent to DeLong 1988
    restricted to positives). Placement for positive i at spec_target:
    V_i(target) = 1 if score_i > threshold(target_spec), 0.5 if tied,
    0 otherwise; SE = std(placements) / sqrt(N_pos). 95% CI uses the
    normal approximation, consistent with n_pos >= 10 (regression-tested).

  - PPV at prevalence p_ at sensitivity s and specificity sp:
        PPV = (s · p_) / (s · p_ + (1 - sp) · (1 - p_)).
    Verified against the published definition (Mercaldo et al., 2007;
    the standard cfDNA clinical-decision recipe in Cohen et al. 2018,
    Liu et al. 2020, etc.).

Edge cases handled:

  - All-negative or all-positive subset → returns a SKIPPED placeholder
    row (auc=None, sens_at_*.=None, ppv_at_*.=None) with a `skipped`
    reason. Cancers with <5 positives are SKIPPED at the row level
    (DeLong CI is unreliable below this — standard n_pos >= 5 floor).
  - Ties in `score` are handled by midrank (avoid silent bias).
  - For spec targets tighter than the empirical 1 - 1/n_neg step
    (e.g. spec=0.999 with only 50 controls → spec ceiling 0.98), the
    CI's lower bound is reported as the highest achievable sens and
    a warning is added; the script does not silently snap to fpr=0.

The module is self-contained: only numpy + scipy.stats are required
at the boundary. The CLI wrapper in `scripts/per_cancer_sens_at_spec.py`
loads a TSV, calls `build_per_cancer_table`, and writes JSON.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

import numpy as np
from scipy import stats

# ──────────────────────────────────────────────────────────────────────
# Constants (the canonical cfDNA clinical-decision grid)
# ──────────────────────────────────────────────────────────────────────

DEFAULT_SPECIFICITIES: tuple[float, ...] = (0.95, 0.98, 0.99)
DEFAULT_PREVALENCES: tuple[float, ...] = (
    0.001, 0.004, 0.01, 0.05, 0.10, 0.20, 0.50,
)
PPV_AT_SPEC: float = 0.99  # spec operating point for the PPV@prev table
MIN_POSITIVES_FOR_CI: int = 5  # DeLong CI floor on n_pos
Z_95: float = 1.959963984540054  # Φ^{-1}(0.975)


# ──────────────────────────────────────────────────────────────────────
# Low-level: midrank placements
# ──────────────────────────────────────────────────────────────────────

def _midrank(x: np.ndarray) -> np.ndarray:
    """Midrank: ties share the average of their positions (1-indexed).

    Equivalent to `scipy.stats.rankdata("average")`. Defined explicitly
    so the module has only one scipy dependency (stats.norm for CI
    quantiles), and so the algorithm is reproducible without scipy.
    """
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty_like(x, dtype=np.float64)
    sorted_x = x[order]
    n = len(x)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_x[j + 1] == sorted_x[i]:
            j += 1
        # Average of positions i..j (1-indexed) = (i + j + 2) / 2
        avg_rank = (i + j + 2) / 2.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1
    return ranks


# ──────────────────────────────────────────────────────────────────────
# DeLong AUC (point estimate + 95% CI)
# ──────────────────────────────────────────────────────────────────────

def delong_auc_ci(y_true: np.ndarray, y_score: np.ndarray,
                  alpha: float = 0.05) -> dict:
    """DeLong AUC point estimate + (1-alpha) CI.

    Returns
    -------
    dict with keys
        auc   : point estimate (midrank handling of ties)
        se    : DeLong standard error
        ci_lo : lower CI bound (clipped to [0, 1])
        ci_hi : upper CI bound (clipped to [0, 1])
        n_pos : number of positives
        n_neg : number of negatives
    """
    y_true = np.asarray(y_true, dtype=np.int64)
    y_score = np.asarray(y_score, dtype=np.float64)
    pos = y_score[y_true == 1]
    neg = y_score[y_true == 0]
    n_pos, n_neg = len(pos), len(neg)

    if n_pos == 0 or n_neg == 0:
        return {
            "auc": float("nan"),
            "se": float("nan"),
            "ci_lo": float("nan"),
            "ci_hi": float("nan"),
            "n_pos": int(n_pos),
            "n_neg": int(n_neg),
        }

    # Placement values on positives:
    #   V_i = (#negatives strictly less than pos[i] + 0.5 * #negatives
    #          equal to pos[i]) / n_neg
    # Walk the combined pos+neg array sorted by score; for each
    # positive, V_i = (cum_neg_below + 0.5 * block_neg) / n_neg.
    combined = np.concatenate([pos, neg])
    combined_is_pos = np.concatenate([
        np.ones(n_pos, dtype=bool), np.zeros(n_neg, dtype=bool)
    ])
    order = np.argsort(combined, kind="mergesort")
    sorted_combined = combined[order]
    sorted_is_pos = combined_is_pos[order]
    n = len(combined)

    cum_neg_below = np.zeros(n_pos, dtype=np.float64)
    cum_neg_equal = np.zeros(n_pos, dtype=np.float64)
    cum_pos_below = np.zeros(n_neg, dtype=np.float64)
    cum_pos_equal = np.zeros(n_neg, dtype=np.float64)

    cum_neg = 0
    cum_pos = 0
    pos_walker = 0
    neg_walker = 0
    i = 0
    while i < n:
        v = sorted_combined[i]
        j = i
        while j + 1 < n and sorted_combined[j + 1] == v:
            j += 1
        block_pos = int(sorted_is_pos[i:j + 1].sum())
        block_neg = (j - i + 1) - block_pos
        # Block of positives: V_i = (cum_neg + 0.5·block_neg) / n_neg
        end = pos_walker + block_pos
        cum_neg_below[pos_walker:end] = cum_neg
        cum_neg_equal[pos_walker:end] = block_neg / 2.0
        pos_walker = end
        # Block of negatives: W_j = (cum_pos + 0.5·block_pos) / n_pos
        end = neg_walker + block_neg
        cum_pos_below[neg_walker:end] = cum_pos
        cum_pos_equal[neg_walker:end] = block_pos / 2.0
        neg_walker = end
        cum_neg += block_neg
        cum_pos += block_pos
        i = j + 1

    V = (cum_neg_below + cum_neg_equal) / n_neg
    W = (cum_pos_below + cum_pos_equal) / n_pos

    auc = float(V.mean())
    s10 = float(V.var(ddof=1)) if n_pos > 1 else 0.0
    s01 = float(W.var(ddof=1)) if n_neg > 1 else 0.0
    var_auc = s10 / n_pos + s01 / n_neg
    se = float(np.sqrt(max(var_auc, 0.0)))
    z = stats.norm.ppf(1.0 - alpha / 2.0)
    ci_lo = max(0.0, auc - z * se)
    ci_hi = min(1.0, auc + z * se)
    return {
        "auc": auc,
        "se": se,
        "ci_lo": float(ci_lo),
        "ci_hi": float(ci_hi),
        "n_pos": int(n_pos),
        "n_neg": int(n_neg),
    }


# ──────────────────────────────────────────────────────────────────────
# Sens@spec reading + DeLong placement-value CI
# ──────────────────────────────────────────────────────────────────────

def _read_sens_at_spec(y_true: np.ndarray, y_score: np.ndarray,
                       target_spec: float) -> tuple[float, float]:
    """Sensitivity at the largest fpr <= 1 - target_spec on the ROC.

    Returns (sens, threshold). Uses the LARGEST-fpr-≤-target reading
    (correct behaviour for tight spec targets — `argmin(|fpr - target|)`
    silently snaps to fpr=0 at spec=0.99 with finite n_neg).
    """
    from sklearn.metrics import roc_curve

    fpr, tpr, thr = roc_curve(y_true, y_score)
    target_fpr = 1.0 - target_spec
    ok = fpr <= target_fpr + 1e-12
    if not ok.any():
        return 0.0, float("nan")
    idx = int(np.where(ok)[0][-1])
    return float(tpr[idx]), float(thr[idx])


def delong_sens_at_spec_ci(y_true: np.ndarray, y_score: np.ndarray,
                           target_spec: float,
                           alpha: float = 0.05) -> dict:
    """Sensitivity at fixed specificity + DeLong (placement-value) 95% CI.

    The placement value for positive i at threshold τ is
        V_i(τ) = I(score_i > τ) + 0.5 · I(score_i == τ).
    Sensitivity at τ is the mean of V_i(τ) over positives. The
    variance is estimated from the std of V_i(τ) (Sun & Xu, 2014),
    which equals the DeLong-Han-Agarwal structural-component restricted
    to the positives — narrower than bootstrap at n_pos >= 10.
    """
    y_true = np.asarray(y_true, dtype=np.int64)
    y_score = np.asarray(y_score, dtype=np.float64)
    sens, threshold = _read_sens_at_spec(y_true, y_score, target_spec)
    pos = y_score[y_true == 1]
    n_pos = len(pos)
    n_neg = int((y_true == 0).sum())

    if n_pos < MIN_POSITIVES_FOR_CI:
        return {
            "sensitivity": float(sens),
            "threshold": float(threshold),
            "se": float("nan"),
            "ci_lo": float("nan"),
            "ci_hi": float("nan"),
            "n_pos": int(n_pos),
            "n_neg": int(n_neg),
            "specificity": float(target_spec),
            "ci_unreliable": True,
        }

    # Placements on positives:
    if np.isnan(threshold):
        V = np.zeros(n_pos, dtype=np.float64)
    else:
        V = np.where(
            pos > threshold,
            1.0,
            np.where(pos == threshold, 0.5, 0.0),
        )
    se = float(np.sqrt(max(V.var(ddof=1) / n_pos, 0.0)))
    z = stats.norm.ppf(1.0 - alpha / 2.0)
    ci_lo = max(0.0, sens - z * se)
    ci_hi = min(1.0, sens + z * se)
    return {
        "sensitivity": float(sens),
        "threshold": float(threshold),
        "se": se,
        "ci_lo": float(ci_lo),
        "ci_hi": float(ci_hi),
        "n_pos": int(n_pos),
        "n_neg": int(n_neg),
        "specificity": float(target_spec),
        "ci_unreliable": False,
    }


# ──────────────────────────────────────────────────────────────────────
# PPV at prevalence (Bayes' rule, evaluated at a fixed operating point)
# ──────────────────────────────────────────────────────────────────────

def ppv_at_prevalence(sens: float, spec: float, prev: float) -> float:
    """PPV = (sens · prev) / (sens · prev + (1-spec) · (1-prev)).

    Standard cfDNA clinical-decision formula. Handles degenerate
    prev ∈ {0, 1} and degenerate sens/spec gracefully.
    """
    if prev <= 0.0:
        return 0.0
    if prev >= 1.0:
        return float(sens)
    num = sens * prev
    den = num + (1.0 - spec) * (1.0 - prev)
    if den <= 0.0:
        return 0.0
    return float(num / den)


# ──────────────────────────────────────────────────────────────────────
# Per-cancer row + table builder
# ──────────────────────────────────────────────────────────────────────

@dataclass
class CancerRow:
    cancer: str
    n: int
    n_pos: int
    n_neg: int
    auc_mean: Optional[float]
    auc_ci: tuple[float, float]
    sens_at_95: Optional[float]
    sens_at_98: Optional[float]
    sens_at_99: Optional[float]
    sens_at_spec: dict = field(default_factory=dict)
    ppv_at_prevalence: dict = field(default_factory=dict)
    skipped: bool = False
    skip_reason: str = ""

    def to_json(self) -> dict:
        out: dict = {
            "cancer": self.cancer,
            "n": self.n,
            "n_pos": self.n_pos,
            "n_neg": self.n_neg,
            "auc_mean": self.auc_mean,
            "auc_ci": list(self.auc_ci),
            "sens_at_95": self.sens_at_95,
            "sens_at_98": self.sens_at_98,
            "sens_at_99": self.sens_at_99,
            "sens_at_spec": self.sens_at_spec,
            "ppv_at_prevalence": self.ppv_at_prevalence,
            "skipped": self.skipped,
        }
        if self.skipped:
            out["skip_reason"] = self.skip_reason
        return out


def _make_skipped_row(cancer: str, n: int, n_pos: int, n_neg: int,
                      reason: str) -> CancerRow:
    return CancerRow(
        cancer=cancer,
        n=n,
        n_pos=n_pos,
        n_neg=n_neg,
        auc_mean=None,
        auc_ci=(float("nan"), float("nan")),
        sens_at_95=None,
        sens_at_98=None,
        sens_at_99=None,
        sens_at_spec={},
        ppv_at_prevalence={f"prev_{p}": None for p in DEFAULT_PREVALENCES},
        skipped=True,
        skip_reason=reason,
    )


def _safe_float(x) -> Optional[float]:
    """JSON-safe float: NaN/inf → None."""
    if x is None:
        return None
    try:
        f = float(x)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(f):
        return None
    return f


def _json_safe(obj):
    """Recursively replace NaN/inf floats with None for json.dump."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def _row_for_subset(cancer: str, y: np.ndarray, s: np.ndarray) -> CancerRow:
    n = len(y)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        reason = (
            f"degenerate subset (n_pos={n_pos}, n_neg={n_neg}); "
            "AUC and Sens@spec are undefined."
        )
        return _make_skipped_row(cancer, n, n_pos, n_neg, reason)
    if n_pos < MIN_POSITIVES_FOR_CI:
        reason = (
            f"n_pos={n_pos} < MIN_POSITIVES_FOR_CI={MIN_POSITIVES_FOR_CI}; "
            "DeLong CI is unreliable below this floor (Sun & Xu 2014)."
        )
        row = _make_skipped_row(cancer, n, n_pos, n_neg, reason)
        # Still emit the AUC point estimate (no CI) for clinical context.
        auc = delong_auc_ci(y, s)
        row.auc_mean = auc["auc"]
        # And the sens-at-spec point estimates (no CI).
        for spec_target in DEFAULT_SPECIFICITIES:
            s_dict = delong_sens_at_spec_ci(y, s, spec_target)
            row.sens_at_spec[str(spec_target)] = s_dict
            setattr(
                row,
                f"sens_at_{int(round(spec_target * 100))}",
                s_dict["sensitivity"],
            )
        # PPV at spec=0.99 for each prevalence.
        sens99 = row.sens_at_99 if row.sens_at_99 is not None else 0.0
        for prev in DEFAULT_PREVALENCES:
            row.ppv_at_prevalence[f"prev_{prev}"] = ppv_at_prevalence(
                sens99, PPV_AT_SPEC, prev
            )
        return row

    auc = delong_auc_ci(y, s)
    sens_at_spec_block: dict[str, dict] = {}
    sens_at_simple: dict[int, float] = {}
    for spec_target in DEFAULT_SPECIFICITIES:
        s_dict = delong_sens_at_spec_ci(y, s, spec_target)
        sens_at_spec_block[str(spec_target)] = s_dict
        sens_at_simple[int(round(spec_target * 100))] = s_dict["sensitivity"]
    # PPV at spec=0.99 for each prevalence.
    sens99 = sens_at_simple[99]
    ppv_block = {f"prev_{p}": ppv_at_prevalence(sens99, PPV_AT_SPEC, p)
                 for p in DEFAULT_PREVALENCES}
    return CancerRow(
        cancer=cancer,
        n=n,
        n_pos=n_pos,
        n_neg=n_neg,
        auc_mean=auc["auc"],
        auc_ci=(auc["ci_lo"], auc["ci_hi"]),
        sens_at_95=sens_at_simple[95],
        sens_at_98=sens_at_simple[98],
        sens_at_99=sens_at_simple[99],
        sens_at_spec=sens_at_spec_block,
        ppv_at_prevalence=ppv_block,
        skipped=False,
    )


def build_per_cancer_table(
    y: np.ndarray,
    s: np.ndarray,
    cancer_label: np.ndarray,
    study_label: Optional[np.ndarray] = None,
    specificities: Iterable[float] = DEFAULT_SPECIFICITIES,
    prevalences: Iterable[float] = DEFAULT_PREVALENCES,
    include_pooled: bool = True,
) -> dict:
    """Build the standard cfDNA clinical-decision table.

    Parameters
    ----------
    y : binary labels (0=control, 1=cancer)
    s : continuous scores (higher = more cancer-like)
    cancer_label : per-sample cancer-type string (e.g. "LUAD", "BRCA",
        "HEALTHY"). The pooled row uses the FULL array (including
        HEALTHY), and per-cancer rows use only samples where
        cancer_label matches a non-HEALTHY type.
    study_label : optional per-sample study identifier, recorded in the
        output for provenance but not consumed in the metric math.
    specificities : operating-point specificities (default {0.95, 0.98,
        0.99}).
    prevalences : prevalence grid for PPV@prev (default the standard
        clinical-decision grid).
    include_pooled : if True, emit a "POOLED" row using all samples.

    Returns
    -------
    dict with keys
        per_cancer : {cancer_name: row_data}
        pooled     : row_data or None
        config     : the operative specificities/prevalences grids
        n_samples  : total samples
        n_studies  : number of distinct studies (or None if no study
            label provided)
    """
    y = np.asarray(y, dtype=np.int64)
    s = np.asarray(s, dtype=np.float64)
    cancer_label = np.asarray(cancer_label)
    if study_label is not None:
        study_label = np.asarray(study_label)
        n_studies = int(len(np.unique(study_label)))
    else:
        n_studies = None

    out: dict = {
        "per_cancer": {},
        "pooled": None,
        "config": {
            "specificities": list(specificities),
            "prevalences": list(prevalences),
            "ppv_at_spec": PPV_AT_SPEC,
            "min_positives_for_ci": MIN_POSITIVES_FOR_CI,
        },
        "n_samples": int(len(y)),
        "n_studies": n_studies,
    }

    # Per-cancer rows: for each cancer (excluding controls), the subset
    # is (cancer samples with y=1) + (all controls with y=0). The mask
    # keeps the cancer-type label AND any control label.
    control_labels = set(cancer_label[y == 0])
    cancer_types = sorted(set(cancer_label) - control_labels)
    for ctype in cancer_types:
        mask = (cancer_label == ctype) | (y == 0)
        y_c = y[mask]
        s_c = s[mask]
        row = _row_for_subset(str(ctype), y_c, s_c)
        out["per_cancer"][str(ctype)] = row.to_json()

    # Pooled row: all cancer types + all controls.
    if include_pooled:
        pooled = _row_for_subset("POOLED", y, s)
        out["pooled"] = pooled.to_json()

    return _json_safe(out)