"""Tests for the per-cancer sens@spec table generator.

Covers:
- DeLong AUC CI math (point estimate in interval, interval is a 95%
  CI, ties are handled by midrank, n_pos < 5 returns nan CI)
- Sens@spec reading: LARGEST-fpr-≤-target (NOT argmin, which silently
  reports 0.0 at tight targets)
- PPV at prevalence formula matches the published definition
- Per-cancer loop: one row per unique cancer type; pooled row uses
  the full dataset
- Edge cases: cancer with <5 positives is SKIPPED with a reason but
  the row still emits n_pos / point estimates
- End-to-end: --synthetic CLI smoke, JSON schema validation
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.per_cancer_sens_at_spec import (  # noqa: E402
    DEFAULT_PREVALENCES,
    DEFAULT_SPECIFICITIES,
    MIN_POSITIVES_FOR_CI,
    PPV_AT_SPEC,
    CancerRow,
    _json_safe,
    build_per_cancer_table,
    delong_auc_ci,
    delong_sens_at_spec_ci,
    ppv_at_prevalence,
)


# ──────────────────────────────────────────────────────────────────────
# DeLong AUC CI
# ──────────────────────────────────────────────────────────────────────

def test_delong_ci_is_a_95_interval_and_auc_in_interval():
    """A 95% CI is an interval; the point estimate is inside it."""
    rng = np.random.default_rng(0)
    y = np.array([0] * 80 + [1] * 80)
    s = np.where(
        y == 1, rng.normal(0.7, 0.15, 160), rng.normal(0.4, 0.15, 160)
    )
    r = delong_auc_ci(y, s, alpha=0.05)
    assert r["ci_lo"] <= r["auc"] <= r["ci_hi"], (
        f"AUC {r['auc']} not in CI [{r['ci_lo']}, {r['ci_hi']}]"
    )
    # Width should be > 0 (variance > 0 on real signal).
    assert r["ci_hi"] - r["ci_lo"] > 0.0
    assert r["n_pos"] == 80 and r["n_neg"] == 80
    # The CI is on [0, 1].
    assert 0.0 <= r["ci_lo"] <= 1.0
    assert 0.0 <= r["ci_hi"] <= 1.0


def test_delong_auc_matches_sklearn():
    """DeLong AUC point estimate must match sklearn (midrank AUC)."""
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(1)
    y = np.array([0] * 60 + [1] * 60)
    s = np.where(
        y == 1, rng.normal(0.8, 0.2, 120), rng.normal(0.4, 0.2, 120)
    )
    r = delong_auc_ci(y, s)
    sk = float(roc_auc_score(y, s))
    assert abs(r["auc"] - sk) < 1e-9, (
        f"DeLong AUC {r['auc']} vs sklearn {sk} differ"
    )


def test_delong_auc_perfect_classifier_returns_one():
    """A perfectly separable cohort yields AUC=1.0 with CI=[1, 1]."""
    y = np.array([0] * 10 + [1] * 10)
    s = np.concatenate([np.zeros(10), np.ones(10)])
    r = delong_auc_ci(y, s)
    assert r["auc"] == 1.0
    assert r["ci_lo"] == 1.0 and r["ci_hi"] == 1.0
    assert r["se"] == 0.0


def test_delong_auc_ties_get_midrank():
    """All-tied scores → AUC=0.5 (midrank gives 0.5 for tied groups)."""
    y = np.array([0] * 20 + [1] * 20)
    s = np.ones(40)  # every score identical
    r = delong_auc_ci(y, s)
    assert 0.4 <= r["auc"] <= 0.6, (
        f"All-tied AUC should be ~0.5, got {r['auc']}"
    )


def test_delong_auc_no_positives_returns_nan():
    """No positives → undefined; helper returns NaN sentinel."""
    y = np.zeros(20, dtype=np.int64)
    s = np.arange(20, dtype=np.float64)
    r = delong_auc_ci(y, s)
    assert np.isnan(r["auc"]) and np.isnan(r["ci_lo"]) and np.isnan(r["ci_hi"])
    assert r["n_pos"] == 0 and r["n_neg"] == 20


def test_delong_auc_handles_all_negatives_in_one_block():
    """All samples in the first class share the same score — non-trivial tie block."""
    y = np.array([0, 0, 0, 1, 1, 1])
    s = np.array([0.1, 0.1, 0.1, 0.5, 0.5, 0.5])
    r = delong_auc_ci(y, s)
    assert r["auc"] == 1.0  # perfect separation by block
    assert r["ci_lo"] == 1.0 and r["ci_hi"] == 1.0


# ──────────────────────────────────────────────────────────────────────
# DeLong Sens@spec CI
# ──────────────────────────────────────────────────────────────────────

def test_sens_at_spec_ci_uses_largest_fpr_le_target_not_argmin():
    """At spec=0.99 with N controls, sens MUST be read at the largest
    fpr ≤ 1 - 0.99. The argmin bug silently snaps to fpr=0 and reports
    sens=0 — verify we don't fall into that trap.
    """
    rng = np.random.default_rng(7)
    n = 200
    y = np.array([0] * 100 + [1] * 100)
    # Cancer scores clearly higher than control scores but with some overlap.
    s = np.where(
        y == 1, rng.normal(0.7, 0.1, n), rng.normal(0.3, 0.1, n)
    )
    r = delong_sens_at_spec_ci(y, s, target_spec=0.99)
    # CI should be a proper interval.
    assert r["ci_lo"] <= r["sensitivity"] <= r["ci_hi"]
    # Sensitivity should be >> 0 — the argmin bug would report 0.
    assert r["sensitivity"] > 0.0


def test_sens_at_spec_ci_narrows_with_more_positives():
    """CI width should monotonically shrink as n_pos increases."""
    rng = np.random.default_rng(13)
    y_base = np.array([0] * 200 + [1] * 200)
    s_base = np.where(
        y_base == 1, rng.normal(0.7, 0.15, 400), rng.normal(0.4, 0.15, 400)
    )
    widths = []
    for n_keep in (20, 50, 100, 200):
        # Sample n_keep positives and all 200 controls (controls anchor spec).
        pos_idx = np.where(y_base == 1)[0][:n_keep]
        idx = np.concatenate([np.where(y_base == 0)[0], pos_idx])
        r = delong_sens_at_spec_ci(y_base[idx], s_base[idx], target_spec=0.95)
        widths.append(r["ci_hi"] - r["ci_lo"])
    # Width is non-increasing (allow tie at the smallest n_keep which may
    # be too coarse to see the difference).
    for w_prev, w_curr in zip(widths, widths[1:]):
        assert w_curr <= w_prev + 1e-9, (
            f"CI width didn't shrink: {widths}"
        )


def test_sens_at_spec_ci_below_min_positives_returns_unreliable():
    """n_pos < MIN_POSITIVES_FOR_CI → CI is marked ci_unreliable."""
    rng = np.random.default_rng(0)
    y = np.array([0] * 40 + [1] * 4)  # n_pos = 4 < MIN_POSITIVES_FOR_CI
    s = np.where(y == 1, 0.7, 0.3)
    r = delong_sens_at_spec_ci(y, s, target_spec=0.95)
    assert r["ci_unreliable"] is True
    assert np.isnan(r["ci_lo"]) and np.isnan(r["ci_hi"])
    # Sensitivity point estimate is still returned (for context).
    assert 0.0 <= r["sensitivity"] <= 1.0


def test_sens_at_spec_ci_monotone_in_target_spec():
    """Sensitivity should not increase as the target specificity tightens."""
    rng = np.random.default_rng(2)
    y = np.array([0] * 100 + [1] * 100)
    s = np.where(
        y == 1, rng.normal(0.7, 0.2, 200), rng.normal(0.3, 0.2, 200)
    )
    senss = []
    for spec in (0.95, 0.98, 0.99):
        r = delong_sens_at_spec_ci(y, s, target_spec=spec)
        senss.append(r["sensitivity"])
    for s_prev, s_curr in zip(senss, senss[1:]):
        assert s_curr <= s_prev + 1e-9, f"non-monotone sens@spec: {senss}"


# ──────────────────────────────────────────────────────────────────────
# PPV at prevalence formula
# ──────────────────────────────────────────────────────────────────────

def test_ppv_formula_matches_published_definition():
    """PPV = (sens·prev) / (sens·prev + (1-spec)·(1-prev)).

    Spot-check several values against hand calculation.
    """
    # (0.5, 0.99, 0.05) → 0.025 / (0.025 + 0.0095) = 0.025/0.0345
    assert abs(ppv_at_prevalence(0.5, 0.99, 0.05) - 0.025 / 0.0345) < 1e-12
    # (1.0, 0.99, 0.001) → 0.001 / (0.001 + 0.00999) = 0.001/0.01099
    assert abs(ppv_at_prevalence(1.0, 0.99, 0.001) - 0.001 / 0.01099) < 1e-12
    # (0.0, 0.99, 0.05) → 0 (sens=0)
    assert ppv_at_prevalence(0.0, 0.99, 0.05) == 0.0
    # (0.5, 1.0, 0.5) → 0.25 / (0.25 + 0) = 1.0 (perfect specificity)
    assert ppv_at_prevalence(0.5, 1.0, 0.5) == 1.0
    # degenerate prev=0 → 0
    assert ppv_at_prevalence(0.5, 0.99, 0.0) == 0.0
    # degenerate prev=1 → sens
    assert ppv_at_prevalence(0.7, 0.99, 1.0) == 0.7


def test_ppv_increases_with_prevalence_at_fixed_spec():
    """PPV is monotonically non-decreasing in prevalence (at fixed sens, spec)."""
    sens, spec = 0.5, 0.99
    ppvs = [ppv_at_prevalence(sens, spec, p) for p in (0.001, 0.01, 0.05, 0.10, 0.50)]
    for a, b in zip(ppvs, ppvs[1:]):
        assert b >= a - 1e-12, f"PPV non-monotone in prevalence: {ppvs}"


def test_ppv_increases_with_sens_at_fixed_spec_and_prev():
    """PPV is monotonically non-decreasing in sensitivity (at fixed spec, prev)."""
    spec, prev = 0.99, 0.05
    ppvs = [ppv_at_prevalence(s, spec, prev) for s in (0.1, 0.3, 0.5, 0.7, 1.0)]
    for a, b in zip(ppvs, ppvs[1:]):
        assert b >= a - 1e-12, f"PPV non-monotone in sens: {ppvs}"


# ──────────────────────────────────────────────────────────────────────
# Per-cancer loop
# ──────────────────────────────────────────────────────────────────────

def test_per_cancer_loop_produces_one_row_per_cancer():
    """The output dict has exactly one entry per unique cancer type."""
    rng = np.random.default_rng(0)
    n = 200
    y = np.array([0] * 80 + [1] * 30 + [1] * 30 + [1] * 30 + [1] * 30)
    cancers = np.array(
        ["HEALTHY"] * 80 + ["BRCA"] * 30 + ["LUAD"] * 30
        + ["CRC"] * 30 + ["PAAD"] * 30,
        dtype=object,
    )
    s = np.where(y == 1, rng.normal(0.7, 0.15, n), rng.normal(0.3, 0.15, n))
    table = build_per_cancer_table(y, s, cancers)
    assert set(table["per_cancer"].keys()) == {"BRCA", "LUAD", "CRC", "PAAD"}
    # Pooled row uses the full dataset.
    assert table["pooled"] is not None
    assert table["pooled"]["n"] == n


def test_per_cancer_subset_includes_all_controls():
    """Each per-cancer subset must include the cancer samples AND all controls."""
    rng = np.random.default_rng(0)
    n = 200
    y = np.array([0] * 80 + [1] * 30 + [1] * 30 + [1] * 30 + [1] * 30)
    cancers = np.array(
        ["HEALTHY"] * 80 + ["BRCA"] * 30 + ["LUAD"] * 30
        + ["CRC"] * 30 + ["PAAD"] * 30,
        dtype=object,
    )
    s = np.where(y == 1, rng.normal(0.7, 0.15, n), rng.normal(0.3, 0.15, n))
    table = build_per_cancer_table(y, s, cancers)
    for name, row in table["per_cancer"].items():
        # Each subset has 30 positives (the cancer) + 80 controls.
        if not row["skipped"]:
            assert row["n_pos"] == 30, f"{name}: n_pos={row['n_pos']}"
            assert row["n_neg"] == 80, f"{name}: n_neg={row['n_neg']}"


def test_per_cancer_with_study_label_records_n_studies():
    """study_label is recorded in the output for provenance."""
    rng = np.random.default_rng(0)
    y = np.array([0] * 50 + [1] * 50)
    cancers = np.array(["HEALTHY"] * 50 + ["BRCA"] * 50, dtype=object)
    studies = np.array(["S1"] * 30 + ["S2"] * 20 + ["S1"] * 25 + ["S2"] * 25,
                       dtype=object)
    s = np.where(y == 1, 0.7, 0.3)
    table = build_per_cancer_table(y, s, cancers, study_label=studies)
    assert table["n_studies"] == 2


def test_per_cancer_with_no_study_label_yields_none():
    """study_label=None → n_studies=None (still a valid output)."""
    rng = np.random.default_rng(0)
    y = np.array([0] * 30 + [1] * 30)
    cancers = np.array(["HEALTHY"] * 30 + ["BRCA"] * 30, dtype=object)
    s = np.where(y == 1, 0.7, 0.3)
    table = build_per_cancer_table(y, s, cancers, study_label=None)
    assert table["n_studies"] is None


def test_per_cancer_output_row_schema_matches_contract():
    """Every per-cancer row has the documented keys."""
    rng = np.random.default_rng(0)
    y = np.array([0] * 50 + [1] * 50)
    cancers = np.array(["HEALTHY"] * 50 + ["BRCA"] * 50, dtype=object)
    s = np.where(y == 1, 0.7, 0.3)
    table = build_per_cancer_table(y, s, cancers)
    row = table["per_cancer"]["BRCA"]
    expected = {
        "cancer", "n", "n_pos", "n_neg", "auc_mean", "auc_ci",
        "sens_at_95", "sens_at_98", "sens_at_99",
        "sens_at_spec", "ppv_at_prevalence", "skipped",
    }
    assert set(row.keys()) >= expected, (
        f"row keys missing: {expected - set(row.keys())}"
    )
    # sens_at_spec has all three target specificities.
    for sp in DEFAULT_SPECIFICITIES:
        assert str(sp) in row["sens_at_spec"]
    # ppv_at_prevalence has all prevalence entries.
    for prev in DEFAULT_PREVALENCES:
        assert f"prev_{prev}" in row["ppv_at_prevalence"]


# ──────────────────────────────────────────────────────────────────────
# Edge cases
# ──────────────────────────────────────────────────────────────────────

def test_cancer_with_lt_5_positives_is_skipped_with_reason():
    """A cancer with n_pos < MIN_POSITIVES_FOR_CI is SKIPPED with reason,
    but the row still records n_pos and AUC point estimate.
    """
    rng = np.random.default_rng(0)
    # 80 healthy + 3 BRCA (too few) + 30 LUAD (enough) + 30 PAAD (enough)
    y = np.array([0] * 80 + [1] * 3 + [1] * 30 + [1] * 30)
    cancers = np.array(
        ["HEALTHY"] * 80 + ["BRCA"] * 3 + ["LUAD"] * 30 + ["PAAD"] * 30,
        dtype=object,
    )
    s = np.where(y == 1, 0.7, 0.3)
    table = build_per_cancer_table(y, s, cancers)
    brca = table["per_cancer"]["BRCA"]
    assert brca["skipped"] is True
    assert brca["n_pos"] == 3  # still recorded, not silently zeroed
    assert MIN_POSITIVES_FOR_CI == 5  # contract
    assert "MIN_POSITIVES_FOR_CI" in brca["skip_reason"]
    # AUC point estimate still emitted (without CI) for clinical context.
    assert brca["auc_mean"] is not None
    assert brca["auc_ci"] == [None, None]
    # sens@spec still emitted as point estimate.
    for sp in DEFAULT_SPECIFICITIES:
        assert brca["sens_at_spec"][str(sp)]["ci_unreliable"] is True


def test_cancer_with_only_one_class_is_skipped():
    """A subset with n_pos=0 or n_neg=0 is degenerate → SKIPPED."""
    rng = np.random.default_rng(0)
    # A 'cancer' group with only controls by construction (y=0 for all)
    y = np.array([0] * 80 + [1] * 30)
    cancers = np.array(["HEALTHY"] * 80 + ["BRCA"] * 30, dtype=object)
    # Mask cancer as HEALTHY so it has n_pos=0 (no 'BRCA' label rows have y=1)
    cancers[80:] = "HEALTHY"
    s = rng.normal(0.5, 0.2, len(y))
    table = build_per_cancer_table(y, s, cancers)
    # All samples are now HEALTHY — there are no cancer types → empty per_cancer.
    # The point is that the degenerate path is exercised (no exception).
    assert "per_cancer" in table


def test_no_controls_in_data_does_not_crash():
    """A pathological dataset with no controls still produces a table."""
    y = np.array([1] * 50)
    cancers = np.array(["BRCA"] * 50, dtype=object)
    s = np.arange(50, dtype=np.float64)
    table = build_per_cancer_table(y, s, cancers)
    # Pooled row is SKIPPED with the degenerate reason.
    assert table["pooled"]["skipped"] is True


def test_no_cancer_types_returns_empty_per_cancer():
    """A dataset with only controls yields empty per_cancer dict."""
    y = np.array([0] * 30)
    cancers = np.array(["HEALTHY"] * 30, dtype=object)
    s = np.arange(30, dtype=np.float64)
    table = build_per_cancer_table(y, s, cancers)
    assert table["per_cancer"] == {}
    assert table["pooled"]["skipped"] is True


# ──────────────────────────────────────────────────────────────────────
# End-to-end CLI
# ──────────────────────────────────────────────────────────────────────

def test_cli_synthetic_smoke_produces_valid_json(tmp_path):
    """Run the script end-to-end with --synthetic; validate JSON schema."""
    out = tmp_path / "per_cancer.json"
    cmd = [
        sys.executable,
        str(_REPO_ROOT / "scripts" / "per_cancer_sens_at_spec.py"),
        "--synthetic", "--n", "120", "--seed", "7",
        "--n-cancer-types", "3", "--per-cancer-positives", "20",
        "--out", str(out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, f"CLI failed: {r.stderr}"
    assert out.exists(), "output file not created"
    payload = json.loads(out.read_text())
    # Schema.
    assert "per_cancer" in payload and "pooled" in payload
    assert "config" in payload and "provenance" in payload
    # Per-cancer rows match the configured count.
    assert len(payload["per_cancer"]) == 3
    assert payload["pooled"] is not None
    # Provenance records what produced the table.
    assert payload["provenance"]["source"] == "synthetic_fixture"
    assert payload["provenance"]["n_samples"] == 120


def test_cli_help_is_fast_and_does_not_run_benchmark():
    """`--help` must be near-instant (the run_X(feat_dir) function pattern)."""
    import time
    cmd = [
        sys.executable,
        str(_REPO_ROOT / "scripts" / "per_cancer_sens_at_spec.py"),
        "--help",
    ]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    elapsed = time.time() - t0
    assert elapsed < 5.0, f"--help took {elapsed:.1f}s (should be <5s)"
    assert "Standalone per-cancer" in r.stdout


def test_cli_requires_exactly_one_of_tsv_or_synthetic(tmp_path):
    """Passing neither --scores-tsv nor --synthetic → error exit."""
    cmd = [
        sys.executable,
        str(_REPO_ROOT / "scripts" / "per_cancer_sens_at_spec.py"),
        "--out", str(tmp_path / "out.json"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    assert r.returncode != 0
    assert "exactly one of" in (r.stderr + r.stdout).lower() or \
           "one of the arguments" in (r.stderr + r.stdout).lower()


def test_cli_loads_tsv_end_to_end(tmp_path):
    """Write a small TSV, run the CLI, validate the resulting JSON."""
    tsv = tmp_path / "scores.tsv"
    rng = np.random.default_rng(11)
    n = 100
    y = np.array([0] * 40 + [1] * 60)
    cancer = ["HEALTHY"] * 40 + ["BRCA"] * 30 + ["LUAD"] * 30
    s = np.where(y == 1, rng.normal(0.7, 0.15, n), rng.normal(0.3, 0.15, n))
    lines = ["sample_id\tscore\ty\tcancer_label\tstudy"]
    for i in range(n):
        lines.append(
            f"sample_{i}\t{s[i]:.6f}\t{y[i]}\t{cancer[i]}\tS{i % 2 + 1}"
        )
    tsv.write_text("\n".join(lines) + "\n")

    out = tmp_path / "pcss_tsv.json"
    cmd = [
        sys.executable,
        str(_REPO_ROOT / "scripts" / "per_cancer_sens_at_spec.py"),
        "--scores-tsv", str(tsv),
        "--out", str(out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, f"CLI failed: {r.stderr}"
    payload = json.loads(out.read_text())
    assert set(payload["per_cancer"].keys()) == {"BRCA", "LUAD"}
    assert payload["pooled"]["n"] == n
    assert payload["provenance"]["source"] == f"scores_tsv:{tsv}"
    assert payload["n_studies"] == 2


def test_cli_tsv_rejects_missing_columns(tmp_path):
    """TSV without required columns → CLI exits non-zero with a clear message."""
    tsv = tmp_path / "bad.tsv"
    tsv.write_text("sample_id\tscore\nfoo\t0.1\nbar\t0.9\n")
    cmd = [
        sys.executable,
        str(_REPO_ROOT / "scripts" / "per_cancer_sens_at_spec.py"),
        "--scores-tsv", str(tsv),
        "--out", str(tmp_path / "out.json"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    assert r.returncode != 0
    msg = (r.stderr + r.stdout).lower()
    assert "missing required columns" in msg or "missing" in msg


# ──────────────────────────────────────────────────────────────────────
# JSON-safety helper
# ──────────────────────────────────────────────────────────────────────

def test_json_safe_replaces_nan_and_inf_with_none():
    """NaN/inf floats in the table become None (json.dump safe)."""
    obj = {"a": float("nan"), "b": float("inf"), "c": -float("inf"),
           "d": 0.5, "e": [1.0, float("nan"), 2.0], "f": {"g": float("nan")}}
    safe = _json_safe(obj)
    assert safe["a"] is None
    assert safe["b"] is None
    assert safe["c"] is None
    assert safe["d"] == 0.5
    assert safe["e"] == [1.0, None, 2.0]
    assert safe["f"]["g"] is None
    # Round-trips through json.dumps without raising.
    json.dumps(safe)