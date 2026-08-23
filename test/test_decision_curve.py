"""Tests for decision_curve + per_specificity_table."""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from src.fragmentomics.decision_curve import (  # noqa: E402
    net_benefit, decision_curve, per_specificity_table,
)


def test_net_benefit_at_high_sensitivity():
    """A perfect classifier should give net_benefit = prevalence at low thresholds."""
    y = np.array([0, 0, 0, 1, 1, 1])
    s = np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])  # perfect
    prev = 0.5
    # At threshold 0.5: predicts all positive correctly, NB = TP/N - FP/N*1
    # = 0.5 - 0 = 0.5 = prevalence
    nb, nb_all, nb_none = net_benefit(y, s, 0.5)
    assert abs(nb - prev) < 1e-9
    # Treat-none baseline is exactly 0
    assert nb_none == 0.0


def test_net_benefit_treat_all_baseline():
    """Treat-all baseline = prevalence - (1-prevalence) * w."""
    y = np.array([1] * 5 + [0] * 5)  # prevalence 0.5
    s = np.array([0.5] * 10)  # model predicts 1 at threshold ≤ 0.5
    # At threshold 0.5, w = 1; model = treat-all in this case
    nb_m, nb_a, _ = net_benefit(y, s, 0.5)
    # NB_all = 0.5 - 0.5 * 1 = 0
    assert abs(nb_a - 0.0) < 1e-9
    # Model predicts pos for everyone = TP=5, FP=5
    # NB_m = 5/10 - 5/10 * 1 = 0
    assert abs(nb_m) < 1e-9


def test_decision_curve_returns_arrays():
    y = np.array([1] * 30 + [0] * 70)
    rng = np.random.default_rng(0)
    s = np.where(y == 1, rng.normal(0.7, 0.2, 100),
                 rng.normal(0.3, 0.2, 100))
    out = decision_curve(y, s, thresholds=np.linspace(0.05, 0.45, 9))
    assert len(out["thresholds"]) == 9
    assert len(out["nb_model"]) == 9
    assert all(0.0 <= nb <= 1.0 for nb in out["nb_model"])
    assert out["clinical_value_range"] is not None


def test_per_specificity_table_at_least_one_row():
    y = np.array([1] * 5 + [0] * 5)
    s = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    rows = per_specificity_table(y, s, specificities=[0.5, 0.9])
    assert len(rows) == 2
    for r in rows:
        assert 0.0 <= r["sensitivity"] <= 1.0
        assert 0.0 <= r["specificity"] <= 1.0
        # Threshold is in [0, 1] OR is +inf (no threshold achieves this spec)
        thr = r["operating_threshold"]
        assert (0.0 <= thr <= 1.0) or (thr == float("inf")) or np.isnan(thr)


def test_per_specificity_monotone_in_spec():
    """Higher specificity requirement → sensitivity should not increase."""
    y = np.array([1] * 30 + [0] * 70)
    rng = np.random.default_rng(0)
    s = np.where(y == 1, rng.normal(0.6, 0.2, 100),
                 rng.normal(0.4, 0.2, 100))
    rows = per_specificity_table(y, s,
                                 specificities=[0.5, 0.7, 0.9, 0.95, 0.99])
    sens = [r["sensitivity"] for r in rows]
    # Sensitivity should be non-increasing with stricter specificity
    for i in range(len(sens) - 1):
        assert sens[i + 1] <= sens[i] + 1e-9, (
            f"non-monotone: {sens} at specificities {rows}")


def test_net_benefit_degenerate_cohort():
    """No-positive cohort → all NBs are 0 (undefined, returned as 0)."""
    y = np.zeros(10)
    s = np.full(10, 0.5)
    nb_m, nb_a, nb_n = net_benefit(y, s, 0.3)
    assert nb_m == 0.0 and nb_a == 0.0 and nb_n == 0.0


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-v"]))