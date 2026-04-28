# New Issues Discovered in Second Review

**Review Date:** 2026-04-28
**Source:** Independent verification of fixes against validation_framework.py API

---

## CRITICAL: Integration Layer Gap Between validation_framework.py and Fixed Files

### N1: Missing Functions in validation_framework.py — 4 Fixed Files Would Crash at Import

**Severity:** CRITICAL (blocks all evaluation pipelines)
**Affected files:** 
- `agent1-variant-calling/evaluate_fixed.py`
- `agent2-multimodal-fusion/evaluate_fixed.py`
- `agent3-longitudinal/run_final_fixed.py`
- `agent3-longitudinal/improved_methods_fixed.py` (partially)

**Root cause:** The validation framework was designed with classes (`BootstrapCI`, `run_with_seeds`), but the fixed files call non-existent wrapper functions (`bootstrap_metrics`, `multi_seed_summary`, `interpolate_sensitivity_at_specificity`, `bootstrap_confidence_interval`). This is a systematic naming mismatch — the framework and the files were developed without a shared API contract.

**Functions called but not existing:**

| Function Called | Called In | Validation Framework Equivalent |
|----------------|-----------|-------------------------------|
| `bootstrap_metrics(y_true, y_score, n_bootstrap=500, seed=42)` | evaluate_fixed.py (agent1, agent2), improved_methods_fixed.py, run_final_fixed.py | `BootstrapCI(seed=seed).compute(y_true, y_pred, y_score)` |
| `interpolate_sensitivity_at_specificity(fpr, tpr, target_specificity)` | evaluate_fixed.py (agent1, agent2) | **NOTHING** — must be newly added |
| `multi_seed_summary(...)` | evaluate_fixed.py (agent1, agent2), run_final_fixed.py | `run_with_seeds(...)` |
| `bootstrap_confidence_interval(...)` | run_final_fixed.py | `BootstrapCI().compute(...)` |

**Impact:** ImportError or AttributeError at runtime. None of the evaluation pipelines can execute.

**Fix options:**
1. **Recommended:** Add `interpolate_sensitivity_at_specificity()` to validation_framework.py. For `bootstrap_metrics` and `multi_seed_summary`, either add thin wrappers or fix callers to use the existing `BootstrapCI` and `run_with_seeds` APIs directly.
2. **Alternative:** Rewrite all callers to use the existing validation_framework.py API.

---

## MODERATE: Wrong API Usage in Fixed Files

### N2: `DataSplitter` Constructor and Return Value Mismatches

**Severity:** MODERATE (prevents execution in 2 files)
**Affected files:** `agent1-variant-calling/evaluate_fixed.py`, `agent3-longitudinal/run_final_fixed.py`

**Problem 1 — Wrong constructor:**
```python
# Called (WRONG):
DataSplitter(train_ratio=0.6, val_ratio=0.2, test_ratio=0.2, seed=seed)

# Actual signature:
DataSplitter(seed: int = 42)
```

The `train_ratio`, `val_ratio`, `test_ratio` are arguments to `split()`, not the constructor.

**Problem 2 — Wrong return value handling:**
```python
# Called (WRONG):
split = splitter.split(np.arange(n), labels)
cal_trajectories = [all_trajs[i] for i in split.val_idx]
test_trajectories = [all_trajs[i] for i in split.test_idx]

# Actual return:
# splitter.split() returns (X_train, y_train, X_val, y_val, X_test, y_test)
# — a 6-tuple of numpy arrays, NOT an object with .val_idx/.test_idx
```

**Fix:** Unpack the tuple and use indices directly or repackage.

---

### N3: `ThresholdCalibrator` Constructor and Return Value Mismatches

**Severity:** MODERATE (prevents execution in run_final_fixed.py)
**Affected file:** `agent3-longitudinal/run_final_fixed.py`

**Problem 1 — Wrong constructor:**
```python
# Called (WRONG):
ThresholdCalibrator(target_metric='specificity', target_value=0.999)

# Actual signature:
ThresholdCalibrator(criterion: str = "youden")
```

**Problem 2 — Wrong return value treatment:**
```python
# Called (WRONG):
cal_result = calibrator.calibrate(cet_scores_cal, cet_labels_cal)
cet_threshold = float(cal_result['threshold'])  # cal_result treated as dict

# Actual return:
# calibrate() returns float (the optimal threshold value)
```

**Fix:**
```python
calibrator = ThresholdCalibrator(criterion='specificity')
cet_threshold = calibrator.calibrate(cet_scores_cal, cet_labels_cal, target_specificity=0.999)
cal_metrics = calibrator.calibration_metrics_  # dict with sensitivity, specificity
```

---

## MODERATE: Paper-Code Inconsistency

### N4: Paper Describes CET Algorithm that Was Removed from Code

**Severity:** MODERATE (paper describes non-existent functionality)
**Affected file:** `paper/main_fixed.tex`, Section 3.3 and 4.4

**Problem:** The paper's Results section describes CET as using "a streak bonus for consecutive measurements above baseline, and a trend bonus for positive log-linear slopes." The Methods section shows the formula:

```
S_t = Σ [log P(m|λ_grow)/P(m|λ_stable) + β_s·𝕀(streak≥3) + β_t·β̂_log-linear]
```

The Limitations section (minor #13) says "CET bonus weights are arbitrary."

**BUT** the fixed code (`improved_methods_fixed.py`) REMOVED all bonuses as part of fix C9. CET now uses pure SPRT. There are no β_s or β_t parameters.

**Status:** Paper and code describe different algorithms. Needs reconciliation.

---

## MINOR: Dead Import in improved_methods_fixed.py

### N5: Non-functional Import Causes ImportError

**Severity:** MINOR (but blocks file loading)
**Affected file:** `agent3-longitudinal/improved_methods_fixed.py`, line ~25

**Problem:**
```python
from validation_framework import bootstrap_metrics
```

`bootstrap_metrics` doesn't exist in validation_framework.py. This import would raise ImportError, preventing the entire file from loading — even though the CET/BOCD/Kalman classes don't use any validation_framework functions (they implement bootstrap inline).

**Fix:** Remove the import line. The file doesn't actually need anything from validation_framework.

---

## MINOR: Missing Function Implementation

### N6: `interpolate_sensitivity_at_specificity()` Not Implemented Anywhere

**Severity:** MINOR (but blocks 2 files)
**Affected files:** `agent1-variant-calling/evaluate_fixed.py`, `agent2-multimodal-fusion/evaluate_fixed.py`

**Problem:** Both evaluate_fixed.py files call `interpolate_sensitivity_at_specificity(fpr, tpr, target_specificity)` but this function exists nowhere in the codebase. It's not in validation_framework.py, not in the fixed files themselves, and not in any other import path.

**Required implementation** (to add to validation_framework.py):
```python
def interpolate_sensitivity_at_specificity(fpr: np.ndarray, tpr: np.ndarray, target_specificity: float) -> float:
    """Interpolate sensitivity at a target specificity using linear interpolation
    between bracketing ROC curve points."""
    specificity = 1.0 - fpr
    # Find points bracketing target specificity
    idx = np.searchsorted(specificity, target_specificity)
    if idx == 0:
        return float(tpr[0])
    if idx >= len(specificity):
        return float(tpr[-1])
    # Linear interpolation
    w = (target_specificity - specificity[idx-1]) / (specificity[idx] - specificity[idx-1])
    return float(tpr[idx-1] + w * (tpr[idx] - tpr[idx-1]))
```

---

## Summary

| ID | Severity | Description | Files Affected | Fix Time |
|----|----------|-------------|----------------|----------|
| N1 | CRITICAL | Missing function exports from validation_framework | 4 files | 1-2 hours |
| N2 | MODERATE | Wrong DataSplitter constructor and return usage | 2 files | 30 min |
| N3 | MODERATE | Wrong ThresholdCalibrator constructor and return usage | 1 file | 20 min |
| N4 | MODERATE | Paper-code CET algorithm mismatch | 1 file | 30 min |
| N5 | MINOR | Dead import blocks improved_methods_fixed.py | 1 file | 5 min |
| N6 | MINOR | interpolate_sensitivity_at_specificity unimplemented | 2 files | 30 min |

**Total estimated fix time:** 3-4 hours for all issues.

**Note:** If the integration is fixed properly, all 6 of these issues become moot and the original 15 fixes (C1-C17) would be verified as correct. These are mechanical API mismatch issues, not conceptual problems.

---

*New Issues Report — Second Review Agent — April 28, 2026*
