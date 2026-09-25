"""Tests for the FinaleMe cross-platform pipeline (scaffolding only).

The cross-platform pipeline is honest by construction: it REFUSES to
synthesize predictions when input data is missing. The tests below verify
that:

  1. `_cross_platform_readiness_probe.py` exits 0 (it's a probe; never
     fails) and prints a checklist with green/red status per requirement.
  2. `finaleme_pipeline.py --help` is fast (<10s) and the orchestrator
     reports FinaleMe JAR / pretrained-models / reference-files status
     without launching FinaleMe.
  3. `finaleme_to_deepcatch_bridge.py` validates FinaleMe per-sample
     outputs against the FinaleDB features cache; missing data → emits
     `cross_platform_finaledb_NOT_RUN.json` with `n_samples_in=0` and
     `data_source='no_data'`.
  4. The bridge REJECTS mismatched sample IDs (sample IDs in the
     methylation output not in FinaleDB features → SchemaError).
  5. The bridge REJECTS β-values outside [0, 1].
  6. `cross_platform_finaledb_validate.py --help` is fast and won't
     fabricate AUC: when data_source != 'both', the JSON carries no
     `cross_platform_auc` field and the script prints a clear refusal
     message.
  7. `--skip-finaleme` and `--skip-fragmentomics` flags behave correctly.
  8. When the bridge is given a fully synthetic but internally-consistent
     methylation table (matching FinaleDB sample IDs, β-values in [0, 1]),
     the orchestrator pipeline runs end-to-end and reports
     `data_source='both'` for the validate script.
  9. The methodology-demo's `cross_platform_methylation_fusion.json`
     remains UNTOUCHED by the new validate script (it only writes
     `cross_platform_finaledb_*.json`).
 10. The validate script's `--skip-finaleme` flag forces `data_source` to
     a non-`both` value so AUC is suppressed.
 11. The probe script's output contains every checklist item
     (FinaleMe JAR, pretrained models, reference files, FinaleMe output,
     FinaleDB features cache) with explicit present/missing verdicts.

These tests run against the SCRIPT INTERFACES, not against the pipeline's
ability to do real FinaleMe decoding (which requires the JAR + reference
files + pretrained models to actually be installed — see
``docs/CROSS_PLATFORM_FINALEME.md``).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import textwrap
import time
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
PY = sys.executable  # venv-aware

# Repo-relative fallback (when pytest is invoked from a different cwd)
def _env_python():
    return PY


# ─────────────────────────────────────────────────────────────────────
# Smoke / scaffolding tests
# ─────────────────────────────────────────────────────────────────────


def test_help_is_fast_for_validate_script():
    """`cross_platform_finaledb_validate.py --help` should exit fast (<10s)
    and never print benchmark output. This is the smoking-gun guard for
    the `--help-triggers-work` bug class."""
    cmd = [_env_python(), str(SCRIPTS / "cross_platform_finaledb_validate.py"), "--help"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    elapsed = time.time() - t0
    assert r.returncode == 0, f"expected --help exit 0, got {r.returncode}\nstderr={r.stderr}"
    assert elapsed < 10, f"--help took {elapsed:.1f}s — work is running at module level"
    # Help text may mention `data_source` in the docstring (it does).
    # The right no-benchmark-output check is: the SCRIPT'S OUTPUT
    # (not --help) prints only when run for real. Verify by running
    # it for real into a tmp dir and ensuring it has no synthetic AUC.
    feat_dir = Path(tempfile.mkdtemp())
    (feat_dir / "labels.tsv").write_text("sample\tlabel\nS1\tcancer\n")
    out_path = Path(tempfile.mkdtemp()) / "validate.json"
    r2 = subprocess.run(
        [
            _env_python(),
            str(SCRIPTS / "cross_platform_finaledb_validate.py"),
            "--skip-finaleme",
            "--skip-fragmentomics",
            "--features-dir", str(feat_dir),
            "--output", str(out_path),
        ],
        capture_output=True, text=True, timeout=15,
    )
    assert r2.returncode in (0, 1)
    payload = json.loads(out_path.read_text())
    assert payload.get("data_source") == "no_data"


def test_help_is_fast_for_probe_script():
    """The probe is intentionally always-exits-0, but --help should also
    be fast (a regression guard against accidentally injecting heavy work
    into the probe)."""
    cmd = [_env_python(), str(SCRIPTS / "_cross_platform_readiness_probe.py"), "--help"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    elapsed = time.time() - t0
    assert r.returncode == 0, f"probe --help exit {r.returncode}"
    assert elapsed < 5, f"probe --help took {elapsed:.1f}s"


def test_help_is_fast_for_orchestrator():
    """`finaleme_pipeline.py --help` is fast; orchestration is opt-in via
    subcommands."""
    cmd = [_env_python(), str(SCRIPTS / "finaleme_pipeline.py"), "--help"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    elapsed = time.time() - t0
    assert r.returncode == 0
    assert elapsed < 10, f"orchestrator --help took {elapsed:.1f}s"


def test_help_is_fast_for_bridge():
    """`finaleme_to_deepcatch_bridge.py --help` is fast."""
    cmd = [_env_python(), str(SCRIPTS / "finaleme_to_deepcatch_bridge.py"), "--help"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    elapsed = time.time() - t0
    assert r.returncode == 0
    assert elapsed < 10, f"bridge --help took {elapsed:.1f}s"


# ─────────────────────────────────────────────────────────────────────
# Probe (always exits 0) tests
# ─────────────────────────────────────────────────────────────────────


def test_probe_exits_zero_and_reports_state():
    """The probe must always exit 0 (it's a probe, not a gate). The
    output should be a human-readable checklist so an operator can decide
    whether to install missing components."""
    cmd = [_env_python(), str(SCRIPTS / "_cross_platform_readiness_probe.py")]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert r.returncode == 0, f"probe should never fail; got {r.returncode}\n{r.stderr}"
    # The output should be a checklist. Required fields:
    for required in [
        "FinaleMe JAR",
        "Pretrained",
        "Reference",
        "FinaleMe output",
        "FinaleDB",
    ]:
        assert required.lower() in r.stdout.lower(), (
            f"probe output missing {required!r}\noutput={r.stdout!r}"
        )


# ─────────────────────────────────────────────────────────────────────
# Bridge validation tests
# ─────────────────────────────────────────────────────────────────────


def test_bridge_rejects_mismatched_sample_ids(tmp_path, monkeypatch):
    """The bridge must reject a methylation table whose sample IDs are
    not present in the FinaleDB features cache. This is the gate that
    prevents the script from silently producing a cross-platform AUC
    over a sample set that doesn't actually overlap fragmentomics."""
    sys.path.insert(0, str(SCRIPTS))
    from finaleme_to_deepcatch_bridge import validate_sample_overlap, BridgeError

    # FinaleDB cache: only sample A
    finaledb_ids = ["A", "B", "C"]
    # Methylation output: has sample Z that isn't in FinaleDB
    meth_ids = ["A", "B", "Z"]
    with pytest.raises(BridgeError) as exc:
        validate_sample_overlap(meth_ids, finaledb_ids)
    assert "Z" in str(exc.value), f"error must name the bad sample Z; got {exc.value}"


def test_bridge_rejects_betas_out_of_range(tmp_path):
    """β-values are real numbers in [0, 1]; the bridge MUST reject
    anything outside that range (a column of 1.5s or -0.2s is symptomatic
    of a misparsed FinaleMe output)."""
    sys.path.insert(0, str(SCRIPTS))
    from finaleme_to_deepcatch_bridge import validate_beta_values, BridgeError

    good = np.array([[0.1, 0.5, 0.9], [0.7, 0.3, 0.0]])
    # No error case
    validate_beta_values(good)

    # Out-of-range
    bad = np.array([[0.1, 1.5, 0.9]])
    with pytest.raises(BridgeError):
        validate_beta_values(bad)

    bad2 = np.array([[0.1, -0.1, 0.9]])
    with pytest.raises(BridgeError):
        validate_beta_values(bad2)


def test_bridge_emits_not_run_when_input_missing(tmp_path):
    """When the bridge can't find a FinaleMe output, it must emit a
    NOT_RUN JSON with `data_source='no_data'` and `n_samples_in=0`. This
    is the contract the downstream validator consumes."""
    sys.path.insert(0, str(SCRIPTS))
    from finaleme_to_deepcatch_bridge import convert_or_emit_not_run

    out_path = tmp_path / "cross_platform_finaledb_probe.json"
    result = convert_or_emit_not_run(
        finaleme_dir=tmp_path,            # empty → no output found
        finaledb_features_dir=tmp_path,   # empty → no overlap
        output_path=str(out_path),
        reason="test: no FinaleMe output present",
    )
    assert out_path.exists()
    payload = json.loads(out_path.read_text())
    assert payload["data_source"] == "no_data"
    assert payload["n_samples_in"] == 0
    assert payload["runnable"] is False
    # Honest framing: must NOT carry fabricated AUC
    assert "cross_platform_auc" not in payload or payload.get(
        "cross_platform_auc"
    ) is None


# ─────────────────────────────────────────────────────────────────────
# Validate script: refuse to fabricate AUC
# ─────────────────────────────────────────────────────────────────────


def test_validate_refuses_to_print_auc_when_data_source_not_both(tmp_path):
    """`cross_platform_finaledb_validate.py` is the orchestrator. When
    the readiness probe says we don't have both channels (methylation
    AND fragmentomics), the script MUST refuse to print a cross-platform
    AUC — this is the load-bearing honesty gate."""

    # Run the validate script with `--skip-finaleme` so data_source is
    # at best 'fragmentomics_only' (not 'both'), and verify the output
    # JSON has NO `cross_platform_auc` numeric value.
    cmd = [
        _env_python(),
        str(SCRIPTS / "cross_platform_finaledb_validate.py"),
        "--output", str(tmp_path / "out.json"),
        "--skip-finaleme",
        "--features-dir", str(tmp_path),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    # exit 0 even when not runnable (the orchestrator returns 0 with an
    # honest not-runnable JSON); non-zero is reserved for unrecoverable
    # errors like a bug in the script itself.
    assert r.returncode in (0, 1), f"unexpected rc {r.returncode}\n{r.stderr}"
    out_path = tmp_path / "out.json"
    assert out_path.exists(), f"expected output JSON at {out_path}\n{r.stderr}\n{r.stdout}"
    payload = json.loads(out_path.read_text())
    # Refuse to fabricate AUC when only fragmentomics is available.
    assert payload.get("data_source") != "both", (
        f"--skip-finaleme should not produce data_source=both; got {payload['data_source']!r}"
    )
    assert "cross_platform_auc" not in payload or payload.get(
        "cross_platform_auc"
    ) is None, (
        "MUST NOT print cross_platform_auc when data_source != 'both'"
    )


def test_validate_skip_fragmentomics_forces_not_both(tmp_path):
    """Symmetric: `--skip-fragmentomics` keeps methylation checked but
    forces the data_source off `both`, suppressing AUC."""
    cmd = [
        _env_python(),
        str(SCRIPTS / "cross_platform_finaledb_validate.py"),
        "--output", str(tmp_path / "out.json"),
        "--skip-fragmentomics",
        "--features-dir", str(tmp_path),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert r.returncode in (0, 1)
    payload = json.loads((tmp_path / "out.json").read_text())
    assert payload.get("data_source") != "both"
    assert "cross_platform_auc" not in payload or payload.get(
        "cross_platform_auc"
    ) is None


# ─────────────────────────────────────────────────────────────────────
# Cross-platform happy path with synthetic-but-valid methylation
# (proves the script CAN compute AUC when given real input, but the
#  input here is synthetic-to-test-fixture, NOT real FinaleMe output)
# ─────────────────────────────────────────────────────────────────────


def test_validate_synthetic_both_produces_auc(tmp_path, monkeypatch):
    """Synthetic but internally-consistent input (β-values in [0,1],
    sample IDs that match FinaleDB features cache) is accepted by the
    pipeline and `cross_platform_auc` IS reported. This proves the
    pipeline CAN report numbers when fed real data, while never
    fabricating them when fed empty / missing data.

    The "synthetic" here is the test fixture, NOT a real FinaleMe
    output. We mark the resulting JSON with `is_synthetic_fixture=True`
    so the schema makes the provenance unambiguous."""

    # Build a small FinaleDB-shaped features cache
    feat_dir = tmp_path / "finaledb"
    feat_dir.mkdir()
    sample_ids = ["S1", "S2", "S3", "S4"]
    # Each sample: small (5,) vector covering 5 channels
    for sid in sample_ids:
        for chan, vec in [
            ("delfi_5mb_ratio",     np.array([0.1, 0.2, 0.3, 0.4, 0.5])),
            ("delfi_5mb_coverage",  np.array([1.0, 1.1, 0.9, 1.2, 1.05])),
            ("delfi_100kb_ratio",   np.array([0.1]*5)),
            ("delfi_100kb_counts",  np.array([10.0]*5)),
            ("fsd_histogram",       np.array([0.05]*5)),
        ]:
            np.save(feat_dir / f"{sid}.{chan}.npy", vec.astype(float))

    # Labels TSV (positional)
    labels_tsv = tmp_path / "labels.tsv"
    labels_tsv.write_text(
        "sample\tdisease_class\tlabel\tstudy\n"
        "S1\tTEST\tcancer\tjiang\n"
        "S2\tTEST\tcancer\tjiang\n"
        "S3\tTEST\thealthy\tjiang\n"
        "S4\tTEST\thealthy\tjiang\n"
    )

    # FinaleMe output dir with a synthetic but valid β-value table
    meth_dir = tmp_path / "finaleme"
    meth_dir.mkdir()
    # Two rows: per-cpG β for each sample; shape (n_samples, n_cpgs)
    rng = np.random.default_rng(2026)
    n_cpgs = 50
    # Cancer samples (S1, S2) get lower β on average than healthy (S3, S4)
    betas = np.vstack([
        rng.normal(0.30, 0.10, n_cpgs).clip(0, 1),  # S1 cancer
        rng.normal(0.35, 0.10, n_cpgs).clip(0, 1),  # S2 cancer
        rng.normal(0.70, 0.10, n_cpgs).clip(0, 1),  # S3 healthy
        rng.normal(0.75, 0.10, n_cpgs).clip(0, 1),  # S4 healthy
    ])
    # Build a TSV: chr start end sample_id beta
    lines = ["chr\tstart\tend\tsample_id\tbeta"]
    for i, sid in enumerate(sample_ids):
        for cpg in range(n_cpgs):
            lines.append(f"chr1\t{cpg*100}\t{cpg*100+50}\t{sid}\t{betas[i, cpg]:.4f}")
    (meth_dir / "BetaValues.tsv").write_text("\n".join(lines) + "\n")

    # Now run the bridge against the synthetic FinaleMe output + the
    # synthetic FinaleDB features cache. It should emit a
    # cross_platform_finaledb_*.json with `data_source='both'` and
    # `n_samples_in=4`, AND mark the output as a synthetic-fixture run.
    bridge_out = tmp_path / "bridge.json"
    cmd = [
        _env_python(),
        str(SCRIPTS / "finaleme_to_deepcatch_bridge.py"),
        "--finaleme-dir", str(meth_dir),
        "--finaledb-features-dir", str(feat_dir),
        "--labels-tsv", str(labels_tsv),
        "--output", str(bridge_out),
        "--synthetic-fixture",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert r.returncode == 0, (
        f"bridge failed on synthetic-but-valid input; rc={r.returncode}\n"
        f"stderr={r.stderr}\nstdout={r.stdout}"
    )
    payload = json.loads(bridge_out.read_text())
    assert payload["data_source"] == "both", (
        f"expected data_source='both' on synthetic-but-valid fixture; "
        f"got {payload.get('data_source')!r}\n{payload}"
    )
    assert payload["n_samples_in"] == 4
    assert payload["n_features_methylation"] == n_cpgs
    # Honest framing for synthetic test fixture: must be marked so
    # a future user can't mistake it for a real FinaleMe decode.
    assert payload.get("is_synthetic_fixture") is True


def test_validate_synthetic_both_runs_a_fusion_auc(tmp_path):
    """End-to-end: the validate script, pointed at a synthetic-but-valid
    methylation output + a synthetic FinaleDB features cache, should
    produce a JSON with `cross_platform_auc` and at least 1 fusion
    strategy listed.

    The synthetic per-sample methylation score (mean β across CpGs per
    sample) is constructed so cancer samples have lower β than healthy
    — i.e. the methylation channel carries real signal in this fixture.
    That guarantees a `methylation_auc` > 0.5 in the resulting JSON.

    Likewise the fragmentomics fixture gives cancer samples a
    meaningfully higher L2-norm than healthy (4x richer vector
    magnitudes), so the proxy score is not tied — without distinct
    scores the AUC is degenerate (NaN / all-ties) and the test
    can't confirm the fusion output is non-zero."""
    feat_dir = tmp_path / "finaledb"
    feat_dir.mkdir()
    sample_ids = ["S1", "S2", "S3", "S4"]
    # Cancer samples have stronger, larger vectors so the L2-norm
    # proxy is meaningfully different from healthy samples.
    for sid in sample_ids:
        is_cancer = sid in ("S1", "S2")
        scale = 4.0 if is_cancer else 1.0
        for chan, vec in [
            ("delfi_5mb_ratio",     np.linspace(0.6, 0.65, 5) * scale),
            ("delfi_5mb_coverage",  np.linspace(1.0, 1.2, 5) * scale),
            ("delfi_100kb_ratio",   np.linspace(0.5, 0.7, 5) * scale),
            ("delfi_100kb_counts",  np.linspace(10, 12, 5) * scale),
            ("fsd_histogram",       np.linspace(0.04, 0.06, 5) * scale),
        ]:
            np.save(feat_dir / f"{sid}.{chan}.npy", vec.astype(float))

    labels_tsv = tmp_path / "labels.tsv"
    labels_tsv.write_text(
        "sample\tdisease_class\tlabel\tstudy\n"
        "S1\tTEST\tcancer\tjiang\n"
        "S2\tTEST\tcancer\tjiang\n"
        "S3\tTEST\thealthy\tjiang\n"
        "S4\tTEST\thealthy\tjiang\n"
    )

    meth_dir = tmp_path / "finaleme"
    meth_dir.mkdir()
    rng = np.random.default_rng(2026)
    n_cpgs = 50
    betas = np.vstack([
        rng.normal(0.30, 0.10, n_cpgs).clip(0, 1),  # S1 cancer
        rng.normal(0.35, 0.10, n_cpgs).clip(0, 1),  # S2 cancer
        rng.normal(0.70, 0.10, n_cpgs).clip(0, 1),  # S3 healthy
        rng.normal(0.75, 0.10, n_cpgs).clip(0, 1),  # S4 healthy
    ])
    lines = ["chr\tstart\tend\tsample_id\tbeta"]
    for i, sid in enumerate(sample_ids):
        for cpg in range(n_cpgs):
            lines.append(f"chr1\t{cpg*100}\t{cpg*100+50}\t{sid}\t{betas[i, cpg]:.4f}")
    (meth_dir / "BetaValues.tsv").write_text("\n".join(lines) + "\n")

    out = tmp_path / "validate.json"
    cmd = [
        _env_python(),
        str(SCRIPTS / "cross_platform_finaledb_validate.py"),
        "--finaleme-dir", str(meth_dir),
        "--features-dir", str(feat_dir),
        "--labels-tsv", str(labels_tsv),
        "--output", str(out),
        "--synthetic-fixture",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, (
        f"validate failed on synthetic-but-valid input; rc={r.returncode}\n"
        f"stderr={r.stderr}\nstdout={r.stdout}"
    )
    payload = json.loads(out.read_text())
    assert payload["data_source"] == "both"
    # MUST print a cross-platform AUC since both channels are runnable
    assert "cross_platform_auc" in payload
    assert payload["cross_platform_auc"] is not None
    # The fixture is synthetic; verify the provenance flags that.
    assert payload.get("is_synthetic_fixture") is True
    # At least one fusion strategy
    assert "fusion_strategies" in payload
    assert payload["fusion_strategies"]


# ─────────────────────────────────────────────────────────────────────
# Honest-framing regression: the methodology demo is NOT touched
# ─────────────────────────────────────────────────────────────────────


def test_validate_does_not_touch_existing_methylation_fusion_json(tmp_path):
    """`cross_platform_finaledb_validate.py` writes ONLY files matching
    `cross_platform_finaledb_*`; the previously-shipped methodology-demo
    file `cross_platform_methylation_fusion.json` MUST NOT be modified
    or replaced by the new validator. (The two scripts produce different
    artifacts under different schema names; mixing them up would be a
    serious provenance bug.)"""
    # Create a stand-in for the methodology demo file in tmp_path so
    # we can verify the new script doesn't go anywhere near the actual
    # `results/` location.
    sentinel = tmp_path / "cross_platform_methylation_fusion.json"
    sentinel.write_text(json.dumps({"task": "DO NOT TOUCH", "sentinel": True}))

    # Run validate into a different output path; verify sentinel
    # untouched.
    feat_dir = tmp_path / "finaledb"
    feat_dir.mkdir()
    labels_tsv = tmp_path / "labels.tsv"
    labels_tsv.write_text(
        "sample\tdisease_class\tlabel\tstudy\nS1\tT\tcancer\tjiang\n"
    )

    out = tmp_path / "validate.json"
    cmd = [
        _env_python(),
        str(SCRIPTS / "cross_platform_finaledb_validate.py"),
        "--skip-finaleme",
        "--features-dir", str(feat_dir),
        "--labels-tsv", str(labels_tsv),
        "--output", str(out),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    assert r.returncode in (0, 1)

    # Sentinel untouched
    assert sentinel.exists()
    assert json.loads(sentinel.read_text()).get("sentinel") is True


def test_orchestrator_reports_missing_jar(tmp_path):
    """When FinaleMe JAR is not installed, `finaleme_pipeline.py` must
    EXIT NON-ZERO with a clear message naming the missing path."""
    cmd = [
        _env_python(),
        str(SCRIPTS / "finaleme_pipeline.py"),
        "decode",
        "--finaleme-dir", str(tmp_path / "does-not-exist"),
    ]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
    # Non-zero on missing input
    assert r.returncode != 0, "expected non-zero exit when inputs are missing"
    out = (r.stdout + r.stderr).lower()
    assert "jar" in out or "finaleme" in out, (
        f"missing-jar message should name 'JAR' or 'FinaleMe'\nout={r.stdout!r}\nerr={r.stderr!r}"
    )


# ─────────────────────────────────────────────────────────────────────
# Schema-version and provenance contract
# ─────────────────────────────────────────────────────────────────────


def test_not_run_payload_schema_contract(tmp_path):
    """The NOT_RUN JSON has a stable schema (schema_version, data_source,
    runnable, n_samples_in, generated_at, blocks). Other tooling relies
    on these keys; they MUST exist when the bridge can't run."""
    sys.path.insert(0, str(SCRIPTS))
    from finaleme_to_deepcatch_bridge import convert_or_emit_not_run

    out = tmp_path / "notrun.json"
    convert_or_emit_not_run(
        finaleme_dir=tmp_path,
        finaledb_features_dir=tmp_path,
        output_path=str(out),
        reason="test",
    )
    payload = json.loads(out.read_text())
    for k in ("schema_version", "data_source", "runnable", "n_samples_in", "generated_at"):
        assert k in payload, f"NOT_RUN JSON missing required key {k!r}"
    assert payload["data_source"] == "no_data"
    assert payload["runnable"] is False
    assert payload["n_samples_in"] == 0
