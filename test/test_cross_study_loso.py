"""Tests for the leave-one-study-out (LOSO) cross-publication benchmark.

`scripts/cross_study_loso.py` measures the transfer that the pooled
cross-study benchmark cannot: train on one publication, test on the
other.

Two layers of test:
1. Unit tests of the helpers (no external data, fast).
2. Artifact tests over `results/cross_study_loso.json`, skipped when the
   artifact is absent (the LOSO run needs the sister pipeline's
   pre-extracted FinaleDB features, which are not committed here).

The artifact tests assert *structural* properties and the direction
ordering, not a single magic AUC, so a legitimate re-run on different
features cannot be failed by a hard-coded constant.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from scripts.cross_study_loso import (  # noqa: E402
    auc_with_ci,
    component_sweep,
    fit_pca_lr,
    harmonize,
    in_domain_cv,
    per_cancer_ovr,
    sens_at_spec,
    spec_grid,
)

ARTIFACT = _REPO_ROOT / "results" / "cross_study_loso.json"
SCRIPT = _REPO_ROOT / "scripts" / "cross_study_loso.py"


# --------------------------------------------------------------------------- #
# Import safety (guards the CI collection failure)
# --------------------------------------------------------------------------- #
def test_sister_pipeline_is_not_imported_at_module_level():
    """The sibling checkout is absent on CI, so it must not be imported
    at module scope -- otherwise collection of this file fails there."""
    for lineno, line in enumerate(SCRIPT.read_text().splitlines(), 1):
        if "from honest_benchmark import" in line:
            assert line[:1].isspace(), (
                f"line {lineno}: module-level import of the sibling "
                "pipeline breaks collection where the sibling repo is absent"
            )
        if "from cross_study_finallydb import" in line:
            raise AssertionError(
                f"line {lineno}: cross_study_finallydb imports the sibling "
                "pipeline at its own module level and must be imported "
                "lazily via importlib"
            )


def test_module_imports_without_the_sibling_pipeline():
    """Importing the script must not require the sibling checkout."""
    import subprocess
    import sys
    code = (
        "import sys; sys.path.insert(0, %r); sys.path.insert(0, %r); "
        "import scripts.cross_study_loso as m; print(m.in_domain_cv.__name__)"
        % (str(_REPO_ROOT), str(_REPO_ROOT / "scripts"))
    )
    env = {"PATH": "/usr/bin:/bin", "DEEPCATCH_SISTER_PIPELINE": "/nonexistent"}
    out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                         text=True, env=env)
    assert out.returncode == 0, out.stderr
    assert "in_domain_cv" in out.stdout


# --------------------------------------------------------------------------- #
# Determinism
# --------------------------------------------------------------------------- #
def test_in_domain_cv_is_deterministic():
    """The in-domain reference must reproduce run-to-run.

    The headline benchmark's PCA is unseeded and its per-publication AUC
    drifts ~0.003 between runs; this reference is seeded, so it must not.
    """
    rng = np.random.default_rng(7)
    n = 80
    y = np.r_[np.zeros(40, int), np.ones(40, int)]
    X = np.vstack([
        rng.normal(0, 1, (40, 30)),
        rng.normal(1.2, 1, (40, 30)),
    ])
    a = in_domain_cv(X, y, [42, 13], 10)
    b = in_domain_cv(X, y, [42, 13], 10)
    assert a == b


def test_fit_pca_lr_is_deterministic():
    rng = np.random.default_rng(8)
    X = rng.normal(size=(50, 15))
    y = np.r_[np.zeros(25, int), np.ones(25, int)]
    p1, m1 = fit_pca_lr(X, y, 8)
    p2, m2 = fit_pca_lr(X, y, 8)
    assert np.allclose(p1.transform(X), p2.transform(X))
    assert np.allclose(m1.predict_proba(p1.transform(X)),
                       m2.predict_proba(p2.transform(X)))


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def test_harmonize_self_zscore_centres_each_publication():
    rng = np.random.default_rng(0)
    Xtr = rng.normal(loc=10.0, scale=1.0, size=(40, 8))
    Xte = rng.normal(loc=-5.0, scale=3.0, size=(30, 8))
    A, B = harmonize(Xtr.copy(), Xte.copy(), "self_zscore")
    # Each publication is centred on its own statistics.
    assert abs(float(A.mean())) < 1e-9
    assert abs(float(B.mean())) < 1e-9
    # The shift between publications is therefore removed.
    assert abs(float(A.mean() - B.mean())) < 1e-9


def test_harmonize_train_global_uses_no_test_statistics():
    rng = np.random.default_rng(1)
    Xtr = rng.normal(loc=0.0, scale=1.0, size=(40, 8))
    Xte = rng.normal(loc=100.0, scale=1.0, size=(30, 8))
    A, B = harmonize(Xtr.copy(), Xte.copy(), "train_global")
    assert abs(float(A.mean())) < 1e-9
    # The test publication keeps its own offset: the train scaler is applied.
    assert float(B.mean()) > 50.0


def test_harmonize_rejects_unknown_mode():
    with pytest.raises(ValueError):
        harmonize(np.zeros((5, 2)), np.zeros((5, 2)), "nope")


def test_sens_at_spec_is_monotone_non_increasing():
    rng = np.random.default_rng(2)
    y = np.array([0, 0, 0, 1, 1, 1])
    score = np.array([0.1, 0.2, 0.3, 0.7, 0.8, 0.9])
    assert sens_at_spec(y, score, 0.95) == pytest.approx(1.0)
    assert sens_at_spec(y, score, 0.99) <= sens_at_spec(y, score, 0.95)


def test_spec_grid_covers_screening_operating_points():
    rng = np.random.default_rng(3)
    y = np.r_[np.zeros(200, int), np.ones(100, int)]
    score = np.r_[rng.normal(0, 1, 200), rng.normal(3, 1, 100)]
    grid = spec_grid(y, score)
    assert {"0.95", "0.98", "0.99", "0.995", "0.999"} <= set(grid)
    assert grid["0.999"] <= grid["0.995"] <= grid["0.95"]


def test_auc_with_ci_brackets_the_point_estimate():
    rng = np.random.default_rng(4)
    y = np.r_[np.zeros(150, int), np.ones(80, int)]
    score = np.r_[rng.normal(0, 1, 150), rng.normal(2.5, 1, 80)]
    res = auc_with_ci(y, score)
    assert res["ci_lo"] <= res["auc"] <= res["ci_hi"]
    assert res["n_pos"] == 80 and res["n_neg"] == 150


def test_per_cancer_ovr_scores_each_type_against_healthy():
    # per_cancer_ovr skips types with n < 5.
    y = np.array([0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1])
    score = np.array([0.1, 0.2, 0.15, 0.05, 0.12, 0.18,
                      0.8, 0.9, 0.85, 0.95, 0.88, 0.92])
    disease = np.array(["HEALTHY"] * 6 + ["LUAD"] * 6)
    rows = per_cancer_ovr(y, score, disease)
    assert "LUAD" in rows
    assert not rows["LUAD"].get("skipped")
    assert rows["LUAD"]["n_pos"] == 6
    assert rows["LUAD"]["auc"] > 0.9


def test_component_sweep_refits_and_is_monotone_in_ks():
    """The sweep must refit the classifier per k, not reuse one fit."""
    rng = np.random.default_rng(5)
    n, p = 60, 20
    y = np.r_[np.zeros(30, int), np.ones(30, int)]
    X = np.vstack([
        rng.normal(0, 1, (30, p)),
        rng.normal(1.5, 1, (30, p)),
    ])
    from sklearn.decomposition import PCA
    pca = PCA(n_components=10, random_state=0).fit(X)
    sweep = component_sweep(pca, X, y, X, y, (2, 5, 10))
    assert set(sweep) == {"2", "5", "10"}
    assert all(0.0 <= v <= 1.0 for v in sweep.values())
    # A trivially separable cohort stays high at every k.
    assert min(sweep.values()) > 0.9


# --------------------------------------------------------------------------- #
# Artifact
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def loso():
    if not ARTIFACT.exists():
        pytest.skip("results/cross_study_loso.json not present "
                    "(requires the sister pipeline's FinaleDB features)")
    with open(ARTIFACT) as f:
        return json.load(f)


def test_artifact_schema(loso):
    for key in ("config", "cohort", "loso", "in_domain",
                "within_study_cross_type_control", "interpretation"):
        assert key in loso, f"missing top-level key: {key}"
    assert loso["cohort"]["n_samples_with_features"] == 627
    assert loso["cohort"]["n_features"] > 1000


def test_both_directions_present_in_both_harmonization_modes(loso):
    combos = {(d["train_study"], d["test_study"], d["harmonization"])
              for d in loso["loso"] if not d.get("skipped")}
    assert ("jiang", "cristiano", "self_zscore") in combos
    assert ("cristiano", "jiang", "self_zscore") in combos
    assert ("jiang", "cristiano", "train_global") in combos
    assert ("cristiano", "jiang", "train_global") in combos


def test_transfer_is_worse_than_in_domain(loso):
    """The structural claim LOSO exists to make.

    A model trained on one publication must score worse on the held-out
    publication than the same pipeline scores within a publication.
    """
    in_domain = {s["study"]: s["cv_auc_mean"] for s in loso["in_domain"]}
    loso_auc = {(d["train_study"], d["test_study"]): d["pooled_test_auc"]["auc"]
                for d in loso["loso"] if not d.get("skipped")
                and d["harmonization"] == "self_zscore"}
    for (tr, te), auc in loso_auc.items():
        assert auc < in_domain[tr], (
            f"{tr}->{te} LOSO AUC {auc} is not below the in-domain "
            f"{tr} AUC {in_domain[tr]}"
        )
        assert auc < in_domain[te]


def test_direction_ordering_holds(loso):
    """Pan-cancer training must beat single-type training on transfer."""
    loso_auc = {(d["train_study"], d["test_study"]): d["pooled_test_auc"]["auc"]
                for d in loso["loso"] if not d.get("skipped")
                and d["harmonization"] == "self_zscore"}
    jc = loso_auc[("jiang", "cristiano")]
    cj = loso_auc[("cristiano", "jiang")]
    assert cj > jc, (
        f"expected train-Cristiano/test-Jiang ({cj}) to beat "
        f"train-Jiang/test-Cristiano ({jc})"
    )


def test_within_study_cross_type_control_transfers(loso):
    """Unseen cancer type within one study is NOT the blocker.

    If this control were also at chance, the Jiang->Cristiano drop would
    be uninterpretable (study shift vs cancer-type specificity).
    """
    ctrl = loso["within_study_cross_type_control"]
    rows = ctrl.get("cristiano", {})
    measured = [r["auc"] for r in rows.values()
                if isinstance(r, dict) and not r.get("skipped")]
    assert len(measured) >= 3
    assert min(measured) > 0.80, (
        "within-study unseen-cancer-type transfer collapsed to "
        f"{min(measured)} -- the LOSO drop cannot then be attributed "
        "to the study holdout"
    )
