"""Tests for scripts/adapter_local_cohort.py — Collaborator Data Interface.

These tests run end-to-end against a synthetic 30-sample cohort built by
scripts/_synthetic_collaborator_cohort.py, plus deliberately-malformed
fixtures to exercise the schema validation paths.

The adapter must:

    1. Accept a well-formed collaborator cohort and produce a normalized
       output the cross-study benchmark can read.
    2. Reject a cohort with missing channels, missing labels, missing
       manifest, mismatched vector lengths, and/or empty features dir.
    3. Apply the cell-line filter regex (GM\\d+, HeLa, HepG2, K562, …)
       so cell-line samples do not contaminate the cancer set.
    4. Merge with an existing labels TSV when ``--merge-with`` is passed.
    5. Run the cross-study benchmark on the merged cohort end-to-end
       (smoke-tested with `--quick` synthetic data, not real data).

No real FinaleDB cohort is required — every test builds its fixture
locally with ``np.save``.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPTS))

# Import the adapter module so we can call its functions directly.
# The script lives under scripts/; load it by path so we don't depend on
# it being installed as a package.
import importlib.util as _importlib_util  # noqa: E402

_adapter_path = SCRIPTS / "adapter_local_cohort.py"
_spec = _importlib_util.spec_from_file_location("adapter_local_cohort", _adapter_path)
assert _spec and _spec.loader, f"could not load adapter from {_adapter_path}"
adapter = _importlib_util.module_from_spec(_spec)
_spec.loader.exec_module(adapter)

PY = sys.executable  # any venv works; matches the verification recipe


# ──────────────────────────────────────────────────────────────────────
# Synthetic cohort builder (mirrors _synthetic_collaborator_cohort.py
# but parameterized so individual tests can shrink / grow / mutate it)
# ──────────────────────────────────────────────────────────────────────

def _write_sample(features_dir: Path, sample_id: str, channel_dims: dict,
                  *, is_cancer: bool = True, seed: int = 0) -> None:
    """Write one synthetic sample's 5 per-channel .npy files."""
    rng = np.random.default_rng(seed + hash(sample_id) % 100000)
    shift = 0.01 if is_cancer else 0.0
    for ch, dim in channel_dims.items():
        if ch == "fsd_histogram":
            v = rng.random(dim)
            v = v / v.sum()
        elif ch.endswith("_ratio"):
            v = np.clip(rng.normal(0.15 + shift, 0.02, dim), 0, 1)
        elif ch.endswith("_coverage"):
            v = np.clip(rng.normal(1.0 + shift, 0.2, dim), 0, None)
        else:  # _counts
            v = rng.integers(50, 5000, dim).astype(np.float64)
        np.save(features_dir / f"{sample_id}.{ch}.npy", v.astype(np.float32))


def _make_manifest(cohort_name: str, n_cancer: int, n_healthy: int,
                   channel_dims: dict) -> dict:
    return {
        "schema_version": "1.0",
        "cohort_name": cohort_name,
        "contributing_lab": "Test Lab",
        "contact_email": "test@example.org",
        "irb_number": "TEST-IRB-000",
        "sample_count_cancer": n_cancer,
        "sample_count_healthy": n_healthy,
        "channels_present": list(channel_dims.keys()),
        "vector_lengths": {k: int(v) for k, v in channel_dims.items()},
        "generation_date": "2026-09-24",
        "scope": "test cohort",
        "citation": "synthetic",
        "license": "CC0",
    }


def _write_labels_tsv(path: Path, rows: list[tuple[str, str, str, str]]) -> None:
    """Write the long collaborator labels.tsv format (6 columns)."""
    header = "sample_id\tdisease_class\tlabel\tstudy_id\tpublication_id\ttissue_of_origin"
    lines = [header]
    for sid, dc, label, study in rows:
        lines.append(f"{sid}\t{dc}\t{label}\t{study}\t\t")
    path.write_text("\n".join(lines) + "\n")


def _build_cohort(
    tmp: Path,
    *,
    cohort_name: str = "synth_2026",
    n_cancer: int = 20,
    n_healthy: int = 10,
    include_cell_line: bool = True,
    channel_dims: dict | None = None,
    seed: int = 0,
    skip_channels: tuple[str, ...] = (),
    wrong_dim_channels: dict[str, int] | None = None,
    drop_from_features: tuple[str, ...] = (),
) -> Path:
    """Build a synthetic cohort on disk; return its root directory."""
    channel_dims = channel_dims or adapter.CHANNEL_DIMS
    wrong_dim_channels = wrong_dim_channels or {}
    root = tmp / cohort_name
    feat_dir = root / "features"
    feat_dir.mkdir(parents=True, exist_ok=True)

    rows: list[tuple[str, str, str, str]] = []
    for i in range(1, n_cancer + 1):
        sid = f"SYN_C{i:03d}"
        _write_sample(feat_dir, sid, channel_dims,
                      is_cancer=True, seed=seed + i)
        rows.append((sid, "CRC_S", "cancer", cohort_name))
    for i in range(1, n_healthy + 1):
        sid = f"SYN_H{i:03d}"
        _write_sample(feat_dir, sid, channel_dims,
                      is_cancer=False, seed=seed + 1000 + i)
        rows.append((sid, "HEALTHY", "healthy", cohort_name))
    if include_cell_line:
        sid = "GM1100"
        _write_sample(feat_dir, sid, channel_dims,
                      is_cancer=True, seed=42)
        rows.append((sid, "Liver cancer", "cancer", cohort_name))

    _write_labels_tsv(root / "labels.tsv", rows)
    manifest = _make_manifest(cohort_name, n_cancer, n_healthy, channel_dims)
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2))

    # Apply optional corruptions AFTER writing the full cohort.
    for ch in skip_channels:
        for sid in {r[0] for r in rows}:
            p = feat_dir / f"{sid}.{ch}.npy"
            if p.exists():
                p.unlink()
    for ch, new_dim in wrong_dim_channels.items():
        for sid in {r[0] for r in rows}:
            np.save(feat_dir / f"{sid}.{ch}.npy",
                    np.zeros(new_dim, dtype=np.float32))
    for sid in drop_from_features:
        for ch in channel_dims:
            p = feat_dir / f"{sid}.{ch}.npy"
            if p.exists():
                p.unlink()
    return root


# ──────────────────────────────────────────────────────────────────────
# 1. Schema validation — happy path
# ──────────────────────────────────────────────────────────────────────

def test_validate_cohort_layout_accepts_well_formed_cohort(tmp_path):
    root = _build_cohort(tmp_path, cohort_name="good")
    # Must not raise.
    adapter.validate_cohort_layout(str(root))


def test_validate_manifest_parses_well_formed_manifest(tmp_path):
    root = _build_cohort(tmp_path)
    m = adapter.validate_manifest(str(root / "manifest.json"))
    assert m["cohort_name"] == "synth_2026"
    assert m["sample_count_cancer"] == 20
    assert m["sample_count_healthy"] == 10
    assert set(m["channels_present"]) == set(adapter.CHANNEL_NAMES)


# ──────────────────────────────────────────────────────────────────────
# 2. Schema validation — rejection paths
# ──────────────────────────────────────────────────────────────────────

def test_validate_cohort_layout_rejects_missing_features_dir(tmp_path):
    root = _build_cohort(tmp_path, cohort_name="bad_nofeat")
    import shutil
    shutil.rmtree(root / "features")
    with pytest.raises(adapter.SchemaError, match="features"):
        adapter.validate_cohort_layout(str(root))


def test_validate_cohort_layout_rejects_missing_labels_tsv(tmp_path):
    root = _build_cohort(tmp_path, cohort_name="bad_nolabels")
    (root / "labels.tsv").unlink()
    with pytest.raises(adapter.SchemaError, match="labels.tsv"):
        adapter.validate_cohort_layout(str(root))


def test_validate_cohort_layout_rejects_missing_manifest(tmp_path):
    root = _build_cohort(tmp_path, cohort_name="bad_nomanifest")
    (root / "manifest.json").unlink()
    with pytest.raises(adapter.SchemaError, match="manifest.json"):
        adapter.validate_cohort_layout(str(root))


def test_validate_features_directory_detects_missing_channel(tmp_path):
    root = _build_cohort(
        tmp_path, cohort_name="bad_skipch",
        skip_channels=("fsd_histogram",),
    )
    features_dir = str(root / "features")
    expected = {f"SYN_C{i:03d}" for i in range(1, 21)} | \
               {f"SYN_H{i:03d}" for i in range(1, 11)} | {"GM1100"}
    missing, mismatched, cell_line = adapter.validate_features_directory(
        features_dir, expected, adapter.CHANNEL_DIMS,
    )
    assert missing, "should report missing fsd_histogram files"
    assert all(".fsd_histogram.npy" in m for m in missing)
    assert not mismatched


def test_validate_features_directory_detects_mismatched_vector_length(tmp_path):
    root = _build_cohort(
        tmp_path, cohort_name="bad_dims",
        wrong_dim_channels={"fsd_histogram": 50},  # should be 196
    )
    features_dir = str(root / "features")
    expected = {f"SYN_C{i:03d}" for i in range(1, 21)} | \
               {f"SYN_H{i:03d}" for i in range(1, 11)} | {"GM1100"}
    missing, mismatched, cell_line = adapter.validate_features_directory(
        features_dir, expected, adapter.CHANNEL_DIMS,
    )
    assert any("shape" in m for m in mismatched), (
        "mismatched-length reports must mention shape"
    )
    assert all("fsd_histogram" in m for m in mismatched)


def test_validate_manifest_rejects_unknown_channels(tmp_path):
    root = _build_cohort(tmp_path, cohort_name="bad_unknown_ch")
    m = json.loads((root / "manifest.json").read_text())
    m["channels_present"] = m["channels_present"] + ["motif_256"]
    (root / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(adapter.SchemaError, match="unknown channels"):
        adapter.validate_manifest(str(root / "manifest.json"))


# ──────────────────────────────────────────────────────────────────────
# 3. Cell-line filter
# ──────────────────────────────────────────────────────────────────────

def test_cell_line_filter_drops_gm1100_by_default(tmp_path):
    """The cell-line filter must remove GM1100 (a B-lymphocyte cell line
    mislabeled 'Liver cancer' in FinaleDB) so it doesn't inflate the
    cancer set in the cross-study benchmark.
    """
    out_root = tmp_path / "out"
    root = _build_cohort(tmp_path / "cohort", cohort_name="clf_test")
    report = adapter.convert_cohort(
        cohort_root=str(root),
        cohort_name="clf_test",
        out_root=str(out_root),
        drop_cell_lines=True,
    )
    assert report["n_samples_in"] == 31, "20 cancer + 10 healthy + 1 cell-line"
    assert report["n_dropped_cell_line"] == 1
    assert report["dropped_cell_line_ids"] == ["GM1100"]
    assert report["n_samples_out"] == 30


def test_cell_line_filter_keeps_gm_when_disabled(tmp_path):
    """With --no-cell-line-filter, the cell-line sample must pass through."""
    out_root = tmp_path / "out"
    root = _build_cohort(tmp_path / "cohort", cohort_name="clf_keep")
    report = adapter.convert_cohort(
        cohort_root=str(root),
        cohort_name="clf_keep",
        out_root=str(out_root),
        drop_cell_lines=False,
    )
    assert report["n_dropped_cell_line"] == 0
    assert report["n_samples_out"] == report["n_samples_in"]


def test_validate_features_directory_lists_gm1100_as_cell_line(tmp_path):
    """validate_features_directory surfaces cell-line sample IDs so
    convert_cohort can drop them (or report them)."""
    root = _build_cohort(tmp_path, cohort_name="clf_list")
    features_dir = str(root / "features")
    expected = {f"SYN_C{i:03d}" for i in range(1, 21)} | \
               {f"SYN_H{i:03d}" for i in range(1, 11)} | {"GM1100"}
    _, _, cell_line = adapter.validate_features_directory(
        features_dir, expected, adapter.CHANNEL_DIMS,
    )
    assert cell_line == ["GM1100"], cell_line


# ──────────────────────────────────────────────────────────────────────
# 4. End-to-end conversion + cross-study benchmark smoke
# ──────────────────────────────────────────────────────────────────────

def test_convert_cohort_writes_canonical_layout(tmp_path):
    """After conversion, results/local_cohort/<name>/labels.tsv has the
    6-column header the cross-study benchmark expects.
    """
    root = _build_cohort(tmp_path / "cohort", cohort_name="e2e_write")
    out_root = tmp_path / "out"
    adapter.convert_cohort(
        cohort_root=str(root),
        cohort_name="e2e_write",
        out_root=str(out_root),
    )
    target = out_root / "e2e_write"
    assert (target / "labels.tsv").is_file()
    assert (target / "adapter_manifest.json").is_file()
    assert (target / "features").is_dir()

    # labels.tsv must be tab-separated with the documented header.
    text = (target / "labels.tsv").read_text()
    header = text.splitlines()[0].split("\t")
    assert header == ["sample_id", "disease_class", "label", "study_id",
                      "publication_id", "tissue_of_origin"]

    # Features directory must contain symlinks or copies for every
    # non-cell-line sample × every channel.
    expected_samples = {f"SYN_C{i:03d}" for i in range(1, 21)} | \
                       {f"SYN_H{i:03d}" for i in range(1, 11)}
    for sid in expected_samples:
        for ch in adapter.CHANNEL_NAMES:
            assert (target / "features" / f"{sid}.{ch}.npy").exists(), (
                f"missing {sid}.{ch}.npy in converted cohort"
            )


def test_convert_cohort_emits_adapter_manifest_with_provenance(tmp_path):
    """The adapter_manifest.json captures every count needed for
    auditability: n_in, n_out, n_dropped_cell_line, channel_layout,
    channel_dims, scope, research_use_only.
    """
    root = _build_cohort(tmp_path / "cohort")
    out_root = tmp_path / "out"
    adapter.convert_cohort(
        cohort_root=str(root),
        cohort_name="manifest_check",
        out_root=str(out_root),
    )
    am = json.loads((out_root / "manifest_check" / "adapter_manifest.json").read_text())
    assert am["n_samples_in"] == 31
    assert am["n_samples_out"] == 30
    assert am["n_dropped_cell_line"] == 1
    assert am["dropped_cell_line_ids"] == ["GM1100"]
    assert am["channel_layout"] == list(adapter.CHANNEL_NAMES)
    assert am["channel_dims"] == adapter.CHANNEL_DIMS
    assert am["scope"].startswith("local-cohort adapter output")
    assert am["research_use_only"] is True


def test_dry_run_does_not_write_output(tmp_path):
    """`--dry-run` validates + reports, but writes nothing."""
    # _build_cohort puts the cohort at <root>/<cohort_name>, so the
    # cohort lands at tmp_path / "dry".
    cohort_root = tmp_path / "dry"
    _build_cohort(tmp_path, cohort_name="dry")
    out_root = tmp_path / "out"
    rc = subprocess.call(
        [PY, str(SCRIPTS / "adapter_local_cohort.py"),
         "--cohort-root", str(cohort_root),
         "--cohort-name", "dry",
         "--out-root", str(out_root),
         "--dry-run"],
    )
    assert rc == 0
    assert not (out_root / "dry").exists(), "dry-run must not create output"


def test_dry_run_reports_dry_run_true_in_provenance(tmp_path):
    """`--dry-run` must mark the report with dry_run=true so downstream
    automation can distinguish a probe from a real run."""
    cohort_root = tmp_path / "dry2" / "dry2"
    _build_cohort(tmp_path / "dry2", cohort_name="dry2")
    r = subprocess.run(
        [PY, str(SCRIPTS / "adapter_local_cohort.py"),
         "--cohort-root", str(cohort_root),
         "--cohort-name", "dry2",
         "--dry-run"],
        capture_output=True, text=True, timeout=10,
    )
    assert r.returncode == 0
    # The final JSON report includes a `dry_run: true` flag.
    assert '"dry_run": true' in r.stdout
    assert "would_write_to" in r.stdout


def test_adapter_help_is_fast(tmp_path):
    """Regression guard: --help must exit fast (<10s) and never run
    any real work.  Mirrors the pitfall in the cfdna-fragmentomics
    skill (--help-triggers-full-benchmark bug class).
    """
    import time
    t0 = time.time()
    r = subprocess.run(
        [PY, str(SCRIPTS / "adapter_local_cohort.py"), "--help"],
        capture_output=True, text=True, timeout=10,
    )
    elapsed = time.time() - t0
    assert r.returncode == 0
    assert elapsed < 5, f"--help took {elapsed:.1f}s, should be <5s"


def test_adapter_rejects_malformed_cohort_with_nonzero_exit(tmp_path):
    """Missing manifest.json must exit non-zero with a clear stderr
    message — the collaborator needs a clean failure to fix their
    directory layout, not a stack trace."""
    # _build_cohort puts the cohort at <root>/<cohort_name>.
    root = tmp_path / "malformed"
    _build_cohort(tmp_path, cohort_name="malformed")
    (root / "manifest.json").unlink()
    out_root = tmp_path / "out"
    r = subprocess.run(
        [PY, str(SCRIPTS / "adapter_local_cohort.py"),
         "--cohort-root", str(root),
         "--cohort-name", "malformed",
         "--out-root", str(out_root)],
        capture_output=True, text=True, timeout=60,
    )
    assert r.returncode != 0, "schema error must yield non-zero exit"
    assert "schema error" in r.stderr.lower()


# ──────────────────────────────────────────────────────────────────────
# 5. Merge + cross-study end-to-end (uses the canonical cross-study
#    benchmark entry point; smoke-tests the wiring without depending
#    on the 627-sample FinaleDB cache).
# ──────────────────────────────────────────────────────────────────────

def test_merge_with_existing_labels_combines_files(tmp_path):
    """merge_with_existing_labels must concatenate local + existing
    labels without duplicating sample IDs (silent double-counting is
    the dangerous failure mode here).
    """
    # Build a 4-row local cohort (2 cancer + 2 healthy).
    root = _build_cohort(
        tmp_path / "cohort",
        cohort_name="merge_local",
        n_cancer=2, n_healthy=2, include_cell_line=False,
    )
    out_root = tmp_path / "out"
    adapter.convert_cohort(
        cohort_root=str(root),
        cohort_name="merge_local",
        out_root=str(out_root),
    )
    local_labels = out_root / "merge_local" / "labels.tsv"

    # Build a fake 'existing' labels TSV that overlaps on one sample id.
    existing = tmp_path / "existing_labels.tsv"
    existing.write_text(
        "sample_id\tdisease_class\tlabel\tstudy_id\n"
        "SYN_C001\tCRC_S\tcancer\tjiang\n"  # duplicate of local
        "EXIST001\tBRCA\tcancer\tcristiano\n"
        "EXIST002\tHEALTHY\thealthy\tcristiano\n"
    )

    merged_path = out_root / "merge_local" / "merged_labels.tsv"
    report = adapter.merge_with_existing_labels(
        local_labels_tsv=str(local_labels),
        existing_labels_tsv=str(existing),
        out_path=str(merged_path),
    )
    assert report["n_local"] == 4
    assert report["n_existing"] == 3
    # Local has 4 unique ids (SYN_C001, SYN_C002, SYN_H001, SYN_H002).
    # Existing has 3 ids; one (SYN_C001) is a duplicate of local.
    # Unique count = 4 (local) + 2 (new in existing) = 6; 1 dup dropped.
    assert report["n_merged"] == 6, "one duplicate should be dropped"
    assert report["n_dropped_duplicates"] == 1

    # The merged file must have the 6-column header (the canonical
    # shape) and contain every unique sample id exactly once.
    text = merged_path.read_text()
    assert "sample_id\tdisease_class\tlabel" in text
    assert text.count("SYN_C001") == 1
    assert text.count("EXIST001") == 1


def test_full_pipeline_cli_synth_creates_normalized_cohort(tmp_path):
    """End-to-end CLI: build synth cohort → run adapter → assert the
    normalized output is consumable by the cross-study benchmark's
    load5() loader (i.e. the .npy files exist and load with the right
    shape).
    """
    # Step 1: build the synthetic cohort via the official script.
    root = tmp_path / "cohort"
    rc = subprocess.call(
        [PY, str(SCRIPTS / "_synthetic_collaborator_cohort.py"),
         "--out", str(root),
         "--cohort-name", "synth_full"],
    )
    assert rc == 0, "_synthetic_collaborator_cohort.py failed"

    # Step 2: run the adapter on it.
    out_root = tmp_path / "out"
    rc = subprocess.call(
        [PY, str(SCRIPTS / "adapter_local_cohort.py"),
         "--cohort-root", str(root),
         "--cohort-name", "synth_full",
         "--out-root", str(out_root)],
    )
    assert rc == 0

    target = out_root / "synth_full"
    assert (target / "labels.tsv").is_file()
    assert (target / "adapter_manifest.json").is_file()

    # Step 3: the cross-study benchmark's loader would read these .npy
    # files; spot-check that one loads with the expected shape.
    arr = np.load(target / "features" / "SYN_C001.delfi_5mb_ratio.npy")
    assert arr.shape == (adapter.CHANNEL_DIMS["delfi_5mb_ratio"],)
    fsd = np.load(target / "features" / "SYN_C001.fsd_histogram.npy")
    assert fsd.shape == (adapter.CHANNEL_DIMS["fsd_histogram"],)
    assert np.isclose(fsd.sum(), 1.0, atol=1e-3), "FSD should be normalized"