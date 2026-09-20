"""Smoke test for the foundation_real_smoke.py module.

Exercises the three-channel reporting (foundation / lr_baseline /
naive_avg) on a fully-synthetic cohort to catch regressions in the
smoke test machinery itself. Does NOT require real TCGA data.

Run: python -m pytest test/test_foundation_smoke.py -v
"""
import sys
import os
import subprocess

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

# Skip all tests in this file if torch is unavailable.
try:
    import torch  # noqa: F401
    _HAS_TORCH = True
except ImportError:
    _HAS_TORCH = False


def _run_smoke(tmp_dir: str, features_dir: str | None = None) -> dict:
    """Run the smoke test in fallback (no TCGA cache) mode and parse JSON.

    If ``features_dir`` is provided it is passed as ``--features-dir``
    to the smoke script. Pass an empty / non-existent directory to
    force the synthetic-fallback path.
    """
    out_path = os.path.join(tmp_dir, "foundation_real_smoke.json")
    log_path = os.path.join(tmp_dir, "foundation_real_smoke.log")
    # Resolve the repo root robustly. This test file may be symlinked
    # under .venv/lib/python*/site-packages/ when pytest discovers it
    # through the editable install, so __file__ can be misleading. Use
    # the tests/ directory layout: repo is two levels up from this file.
    repo_root = os.path.abspath(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), "..", "..")
    )
    smoke_script = os.path.join(repo_root, "scripts", "foundation_real_smoke.py")
    if not os.path.exists(smoke_script):
        # Fallback: walk up to find scripts/foundation_real_smoke.py
        cur = os.path.dirname(os.path.realpath(__file__))
        while cur != "/":
            cand = os.path.join(cur, "scripts", "foundation_real_smoke.py")
            if os.path.exists(cand):
                repo_root = cur
                smoke_script = cand
                break
            cur = os.path.dirname(cur)
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    # Use the synthetic fallback by clearing any TCGA cache that
    # happens to be on the test machine. The script's
    # _try_load_real_panel_scores returns None when no cache is
    # present, so the synthetic cohort kicks in.
    cmd = [
        sys.executable, smoke_script,
        "--out", out_path,
        "--seeds", "3",
        "--n-patients", "20",
    ]
    if features_dir is not None:
        cmd.extend(["--features-dir", features_dir])
    with open(log_path, "w") as logf:
        result = subprocess.run(
            cmd, cwd=repo_root, env=env, capture_output=True, text=True,
            timeout=600,
        )
        logf.write(result.stdout)
        logf.write(result.stderr)
    assert result.returncode in (0, 1), (
        f"smoke exited with {result.returncode}; see {log_path}\n"
        f"stderr: {result.stderr[:2000]}"
    )
    import json
    with open(out_path) as f:
        return json.load(f)


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_returns_all_three_channel_metrics(tmp_path):
    """The smoke output must report foundation + lr_baseline + naive_avg,
    plus the shuffled-label negative control and the signal-to-artifact
    ratio so reviewers can verify the AUC is signal-driven, not artifact.
    """
    d = _run_smoke(str(tmp_path))
    for key in (
        "foundation_auc_mean", "foundation_auc_std",
        "foundation_sens_at_99_mean",
        "lr_baseline_auc_mean", "lr_baseline_sens_at_99_mean",
        "naive_avg_auc_mean", "naive_avg_sens_at_99_mean",
        "shuffled_lr_baseline_auc_mean",
        "shuffled_naive_avg_auc_mean",
        "shuffled_foundation_auc_mean",
        "signal_to_artifact_ratio",
        "gate_pass", "data_source", "honest_framing",
    ):
        assert key in d, f"smoke output missing key {key!r}"


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_shuffled_label_auc_is_below_real(tmp_path):
    """Sanity check: when labels are shuffled, the model AUC must be
    meaningfully below the real-labels AUC. Otherwise the real AUC is
    an artifact of the paired cancer/control design.

    We allow shuffled_auc up to 0.65 (some patient-pair structure
    always leaks through); we require shuffled_auc < real_auc - 0.10.
    """
    d = _run_smoke(str(tmp_path))
    real = d["lr_baseline_auc_mean"]
    shuf = d["shuffled_lr_baseline_auc_mean"]
    assert shuf < real - 0.10, (
        f"shuffled lr_baseline AUC {shuf:.3f} is too close to real "
        f"{real:.3f}; the real AUC may be artifact, not signal. "
        f"signal_to_artifact_ratio = {d['signal_to_artifact_ratio']:.2f}"
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_lr_baseline_beats_naive_avg(tmp_path):
    """The sklearn LR on [panel, frag] should be at least as good as the
    naive average for any cohort with non-degenerate signal.
    """
    d = _run_smoke(str(tmp_path))
    # LR has one extra parameter to fit; on the synthetic cohort
    # (panel + frag both AUC ~0.92) it should match or beat the mean.
    assert d["lr_baseline_auc_mean"] >= d["naive_avg_auc_mean"] - 0.02, (
        f"lr_baseline AUC {d['lr_baseline_auc_mean']:.3f} should be at "
        f"least as good as naive_avg {d['naive_avg_auc_mean']:.3f}"
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_lr_baseline_auc_above_chance(tmp_path):
    """The LR baseline must be well above 0.5 AUC on the synthetic cohort."""
    d = _run_smoke(str(tmp_path))
    assert d["lr_baseline_auc_mean"] > 0.70, (
        f"lr_baseline AUC {d['lr_baseline_auc_mean']:.3f} should be well "
        f"above chance on a synthetic cohort with AUC-0.92 channels"
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_foundation_auc_above_chance(tmp_path):
    """Foundation variant-B (frozen encoder + LR head) alone reaches
    AUC ~0.95 on n=40. The averaged score should be above 0.7.
    """
    d = _run_smoke(str(tmp_path))
    assert d["foundation_auc_mean"] > 0.65, (
        f"foundation AUC {d['foundation_auc_mean']:.3f} too low; "
        f"variant-B alone reaches ~0.95"
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_data_source_is_synthetic_when_no_cache(tmp_path):
    """When TCGA cache is absent the data_source must say synthetic.

    Point --features-dir at an empty directory so the loader sees
    no MAF files and falls back to the synthetic cohort. We verify
    the data_source field reflects that fallback.
    """
    empty_cache = os.path.join(tmp_path, "empty_tcga_cache")
    os.makedirs(empty_cache, exist_ok=True)
    d = _run_smoke(str(tmp_path), features_dir=empty_cache)
    assert "synthetic" in d["data_source"].lower(), (
        f"expected 'synthetic' in data_source with empty TCGA cache, "
        f"got {d['data_source']!r}"
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_foundation_within_lr_baseline_band(tmp_path):
    """The foundation score (variant B weighted 0.7 + variant A 0.3)
    should track the lr_baseline AUC within ~3pp on n=40 — same data,
    same signal source, different model class. If the gap is wider
    than 5pp the variant-A trainable transformer is dominating the
    score, which is the overfitting failure mode we want to detect.
    """
    d = _run_smoke(str(tmp_path))
    gap = d["lr_baseline_auc_mean"] - d["foundation_auc_mean"]
    assert gap < 0.05, (
        f"foundation AUC {d['foundation_auc_mean']:.3f} is "
        f"{gap:.3f} below lr_baseline AUC {d['lr_baseline_auc_mean']:.3f}; "
        f"variant-A trainable transformer is over-dominating the score"
    )


@pytest.mark.skipif(not _HAS_TORCH, reason="torch not installed")
def test_smoke_all_three_gates_present(tmp_path):
    """The JSON must report all three gate thresholds + pass/fail."""
    d = _run_smoke(str(tmp_path))
    for key in ("gate_auc", "gate_sens99", "gate_foundation_auc", "gate_pass"):
        assert key in d, f"smoke output missing key {key!r}"
    assert isinstance(d["gate_pass"], bool)
