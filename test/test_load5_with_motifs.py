"""Tests for scripts/load5_with_motifs.py and the --include-motifs flag.

These tests pin four contracts:

1. ``load5_with_motifs`` adds exactly +256 dims over ``load5_only``
   (63_246 + 256 = 63_502-dim feature vector per sample).
2. ``load5_with_motifs_optional`` keeps all 5-DELFI samples, zero-filling
   the motif block when motifs.npy is absent (graceful skip, head-to-head
   comparison vs the OFF-default benchmark).
3. The motif block per sample sums to ~1.0 (the .motifs.npy files are
   stored as a 4-mer frequency distribution).
4. ``scripts/cross_study_finallydb.py --include-motifs`` exits cleanly on
   --help and the CLI surface exposes the flag.

The end-to-end "5-seed x 5-fold with --include-motifs produces a different
pooled AUC" check is a separate test marked slow and skipped by default
in tight CI; run it manually with:
    pytest test/test_load5_with_motifs.py -m slow -v
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
_SCRIPTS = _REPO_ROOT / "scripts"
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_SCRIPTS))

from load5_with_motifs import (  # noqa: E402
    DEFAULT_FEATURES_DIR,
    MOTIF_DIM,
    load5_only,
    load5_with_motifs,
    load5_with_motifs_optional,
)

# Cross-platform skip: tests need /Users/hermes/cfdna-fragmentomics-pipeline/data/features
# (gitignored). CI runners don't have it.
_FEAT_DIR = Path(os.environ.get(
    "CFDNA_FEATURES_DIR",
    "/Users/hermes/cfdna-fragmentomics-pipeline/data/features",
))
_HAS_FEATURES = (
    _FEAT_DIR.is_dir()
    and any(_FEAT_DIR.glob("*.delfi_5mb_ratio.npy"))
    and any(_FEAT_DIR.glob("*.motifs.npy"))
)

pytestmark = pytest.mark.skipif(
    not _HAS_FEATURES,
    reason=(
        f"Local FinaleDB features cache not at {_FEAT_DIR} "
        "(data/features/ is gitignored; CI lacks cohort)"
    ),
)


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #

def _build_labels_from_dir(feat_dir: Path, max_n: int = 50):
    """Build a synthetic labels dict from the files in `feat_dir`.

    Uses any sample with a motifs.npy as label=1 (pseudo-cancer) and the
    next ten non-motif samples as label=0 (pseudo-healthy). This keeps
    the fixture fast (≤50 samples) and exercises the loader's actual
    file paths without mocking.
    """
    motif_samples = sorted(
        p.name[:-len(".motifs.npy")]
        for p in feat_dir.glob("*.motifs.npy")
    )[:max_n]
    # Pick a few non-motif samples for the optional-loader test
    non_motif_samples = sorted(
        p.name[:-len(".delfi_5mb_ratio.npy")]
        for p in feat_dir.glob("*.delfi_5mb_ratio.npy")
        if not (feat_dir / f"{p.name[:-len('.delfi_5mb_ratio.npy')]}.motifs.npy").exists()
    )[:10]
    labels = {s: 1 for s in motif_samples}
    labels.update({s: 0 for s in non_motif_samples})
    return labels


@pytest.fixture
def feat_dir() -> Path:
    return _FEAT_DIR


@pytest.fixture
def small_labels(feat_dir) -> dict:
    return _build_labels_from_dir(feat_dir, max_n=50)


# --------------------------------------------------------------------------- #
# Unit tests for the loader
# --------------------------------------------------------------------------- #

def test_load5_only_dim_is_63246(small_labels, feat_dir):
    """The 5-channel DELFI + FSD vector is 63,246-dim per sample."""
    X, y, st = load5_only(small_labels, {}, str(feat_dir))
    if X.shape[0] == 0:
        pytest.skip("no samples with all 5 DELFI channels in this cache")
    assert X.shape[1] == 63_246, (
        f"expected 63,246 dims for load5_only, got {X.shape[1]}"
    )
    assert y.shape[0] == X.shape[0]
    assert st.shape[0] == X.shape[0]


def test_load5_with_motifs_includes_256_dim_motif_vector(small_labels, feat_dir):
    """load5_with_motifs adds exactly +256 dims over load5_only.

    Per the loader docstring, X.shape = (n_kept, 63_246 + 256) = (n_kept, 63_502).
    """
    X5, _, _ = load5_only(small_labels, {}, str(feat_dir))
    Xm, _, _ = load5_with_motifs(small_labels, {}, str(feat_dir))
    if Xm.shape[0] == 0:
        pytest.skip("no samples with motifs.npy in this cache")
    assert Xm.shape[1] == X5.shape[1] + MOTIF_DIM == 63_502, (
        f"expected 63,502 dims (63,246 + {MOTIF_DIM}), "
        f"load5_only={X5.shape[1]} load5_with_motifs={Xm.shape[1]}"
    )


def test_load5_with_motifs_validates_sum_to_one(small_labels, feat_dir):
    """The motif block (last 256 dims) sums to ~1.0 per sample.

    .motifs.npy files are stored as 4-mer frequency distributions
    (sum to 1.0). This test pins the convention; if the upstream
    feature extractor ever switches to raw counts, this test fails.
    """
    X, _, _ = load5_with_motifs(small_labels, {}, str(feat_dir))
    if X.shape[0] == 0:
        pytest.skip("no samples with motifs.npy in this cache")
    motif_block = X[:, -MOTIF_DIM:]
    sums = motif_block.sum(axis=1)
    # All motif blocks should sum to ~1.0 (frequency distribution).
    assert np.allclose(sums, 1.0, atol=1e-3), (
        f"motif block sums not ~1.0; min={sums.min():.4f} "
        f"max={sums.max():.4f} mean={sums.mean():.4f}"
    )


def test_load5_with_motifs_handles_missing_motif_file(small_labels, feat_dir):
    """load5_with_motifs_optional keeps all 5-DELFI samples by zero-filling
    the motif block when motifs.npy is absent. This is the head-to-head
    comparison mode for the --include-motifs benchmark.
    """
    # Build a labels dict where ALL samples have the 5 DELFI channels
    # but some are missing motifs.npy
    all_5chan_samples = sorted(
        p.name[:-len(".delfi_5mb_ratio.npy")]
        for p in feat_dir.glob("*.delfi_5mb_ratio.npy")
        if (feat_dir /
            f"{p.name[:-len('.delfi_5mb_ratio.npy')]}.fsd.json").exists()
    )
    labels = {s: 1 for s in all_5chan_samples[:50]}
    X_optional, y_optional, _ = load5_with_motifs_optional(
        labels, {}, str(feat_dir),
    )
    X_strict, y_strict, _ = load5_with_motifs(
        labels, {}, str(feat_dir),
    )
    if X_optional.shape[0] == 0:
        pytest.skip("no 5-channel samples in this cache")
    # The optional loader must keep MORE samples than strict
    # (strict drops motifs-missing samples; optional zero-fills).
    assert X_optional.shape[0] >= X_strict.shape[0]
    # And must produce 63,502-dim vectors.
    assert X_optional.shape[1] == 63_502, (
        f"expected 63,502 dims, got {X_optional.shape[1]}"
    )
    # Verify zero-fill: count how many motif blocks sum to 0 (zero-filled)
    # vs ~1 (real motifs).
    motif_block = X_optional[:, -MOTIF_DIM:]
    sums = motif_block.sum(axis=1)
    n_zero_filled = int((sums == 0).sum())
    n_real = int((np.abs(sums - 1.0) < 1e-3).sum())
    # At least one of each OR all-zero (single-mode cohort)
    assert n_zero_filled + n_real == X_optional.shape[0]
    if X_strict.shape[0] > 0:
        # Some samples were zero-filled (i.e., not in strict subset)
        assert n_zero_filled > 0, (
            "optional loader returned same sample count as strict; "
            "expected at least one zero-filled motif block"
        )


def test_load5_with_motifs_strict_drops_missing(small_labels, feat_dir):
    """load5_with_motifs (strict) drops any sample missing motifs.npy.

    This is the contract that reduces the cohort to the ~98 samples
    with motifs cached locally.
    """
    X_strict, y_strict, _ = load5_with_motifs(small_labels, {}, str(feat_dir))
    X5, _, _ = load5_only(small_labels, {}, str(feat_dir))
    if X5.shape[0] == 0:
        pytest.skip("no 5-DELFI samples in this cache")
    # Strict must return FEWER samples than load5_only (it drops no-motif)
    assert X_strict.shape[0] <= X5.shape[0]
    # And every strict sample must have a real motif block
    if X_strict.shape[0] > 0:
        sums = X_strict[:, -MOTIF_DIM:].sum(axis=1)
        assert np.allclose(sums, 1.0, atol=1e-3), (
            f"strict loader returned a sample whose motif block "
            f"does not sum to 1.0; min sum={sums.min():.4f}"
        )


# --------------------------------------------------------------------------- #
# CLI surface tests for cross_study_finallydb.py --include-motifs
# --------------------------------------------------------------------------- #

def test_cross_study_pipeline_accepts_include_motifs_flag():
    """`cross_study_finallydb.py --help` must:
    - exit 0 fast (≤ 30s — guards against the --help-triggers-work bug)
    - show --include-motifs in the help output
    - show --motif-cohort-mode in the help output
    """
    py = sys.executable  # works in any venv
    t0 = time.time()
    r = subprocess.run(
        [py, str(_SCRIPTS / "cross_study_finallydb.py"), "--help"],
        capture_output=True, text=True, timeout=30,
        cwd=str(_REPO_ROOT),
        env={**os.environ, "PYTHONPATH": ""},
    )
    elapsed = time.time() - t0
    assert r.returncode == 0, f"--help failed: {r.stderr}"
    assert elapsed < 30, f"--help took {elapsed:.1f}s — too slow"
    assert "--include-motifs" in r.stdout, (
        "--include-motifs not in --help output; flag not wired into argparse"
    )
    assert "--motif-cohort-mode" in r.stdout, (
        "--motif-cohort-mode not in --help output"
    )


# --------------------------------------------------------------------------- #
# Slow / end-to-end test (manually run; default skipped in CI)
# --------------------------------------------------------------------------- #

@pytest.mark.slow
def test_with_motifs_pooled_auc_differs_from_baseline(feat_dir):
    """End-to-end: --include-motifs produces a different pooled AUC
    than the OFF-default baseline.

    This test runs the full 5-seed x 5-fold benchmark TWICE (once with
    --include-motifs, once without) and compares pooled AUC. Total cost
    ~15 min on M-series. Marked @pytest.mark.slow so tight CI skips it.

    Pass criterion: with-motifs and without-motifs pooled AUCs are both
    computed; the test reports the delta honestly. A delta of 0 is a
    real outcome ("motifs add zero value") that we record but do not
    treat as a failure.
    """
    labels_tsv = (
        feat_dir.parent / "labels_multiclass.tsv"
    )
    if not labels_tsv.exists():
        labels_tsv = Path(
            "/Users/hermes/cfdna-fragmentomics-pipeline/labels_multiclass.tsv"
        )
    if not labels_tsv.exists():
        pytest.skip(f"labels_multiclass.tsv not found at {labels_tsv}")

    py = sys.executable
    base_args = [
        "--features-dir", str(feat_dir),
        "--labels-multiclass", str(labels_tsv),
        "--seeds", "42", "13", "7", "99", "1234",
        "--out-per-cancer-json", "",
    ]

    # Run WITHOUT --include-motifs (existing baseline)
    base_json = _REPO_ROOT / "results" / "_test_load5_with_motifs_baseline.json"
    r1 = subprocess.run(
        [py, str(_SCRIPTS / "cross_study_finallydb.py"),
         *base_args,
         "--out-json", str(base_json),
         "--out-md", str(_REPO_ROOT / "docs" / "_test_load5_with_motifs_baseline.md")],
        capture_output=True, text=True, timeout=1200,
        cwd=str(_REPO_ROOT),
        env={**os.environ, "PYTHONPATH": ""},
    )
    assert r1.returncode == 0, f"baseline run failed: {r1.stderr[-500:]}"

    # Run WITH --include-motifs (new flag)
    motifs_json = _REPO_ROOT / "results" / "_test_load5_with_motifs_with_motifs.json"
    r2 = subprocess.run(
        [py, str(_SCRIPTS / "cross_study_finallydb.py"),
         *base_args,
         "--include-motifs",
         "--motif-cohort-mode", "zero_fill",
         "--out-json", str(motifs_json),
         "--out-md", str(_REPO_ROOT / "docs" / "_test_load5_with_motifs_with_motifs.md")],
        capture_output=True, text=True, timeout=1200,
        cwd=str(_REPO_ROOT),
        env={**os.environ, "PYTHONPATH": ""},
    )
    assert r2.returncode == 0, f"with-motifs run failed: {r2.stderr[-500:]}"

    with open(base_json) as f:
        base = json.load(f)
    with open(motifs_json) as f:
        motifs = json.load(f)

    base_auc = base["pooled"]["harmonized"]["auc_mean"]
    motifs_auc = motifs["pooled"]["harmonized"]["auc_mean"]
    delta = motifs_auc - base_auc

    # The test PASSES when both runs complete; the delta is reported
    # via a print for the human reviewer to inspect. A delta of zero
    # is acceptable (motifs add zero value on this cohort) — not a test
    # failure.
    print(
        f"\n[test_with_motifs_pooled_auc_differs_from_baseline] "
        f"baseline pooled AUC (load5, 63,246-dim): {base_auc:.4f} "
        f"± {base['pooled']['harmonized']['auc_std']:.4f}\n"
        f"with-motifs pooled AUC (load5+256 motifs, 63,502-dim): "
        f"{motifs_auc:.4f} ± {motifs['pooled']['harmonized']['auc_std']:.4f}\n"
        f"DELTA: {delta:+.4f} "
        f"({'POSITIVE: motifs help' if delta > 0 else 'NEGATIVE: motifs add noise' if delta < 0 else 'ZERO: motifs add no value'})"
    )

    # Cleanup the intermediate test artifacts
    for p in (base_json, motifs_json,
              _REPO_ROOT / "docs" / "_test_load5_with_motifs_baseline.md",
              _REPO_ROOT / "docs" / "_test_load5_with_motifs_with_motifs.md"):
        try:
            p.unlink()
        except FileNotFoundError:
            pass

    assert base["config"]["motif_features"]["included"] is False
    assert motifs["config"]["motif_features"]["included"] is True
    assert motifs["config"]["motif_features"]["dim"] == 256
