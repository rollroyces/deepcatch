#!/usr/bin/env python3
"""Adapter for local collaborator cohorts (Collaborator Data Interface).

External clinical collaborators (or any third party who wants to
contribute a plasma cohort to the benchmark) provide their data in the
standardized layout documented at
``docs/COLLABORATOR_DATA_INTERFACE.md``:

    <cohort_root>/features/<sample_id>.<channel>.npy      (5 channels)
    <cohort_root>/labels.tsv                              (TSV)
    <cohort_root>/manifest.json                           (cohort metadata)

This adapter:

    1. **Validates** the directory layout, the .npy vector lengths,
       and the label/feature counts (mandatory step — silent bugs in
       cohort conversion have destroyed real benchmark results).
    2. **Filters out cell-line technical controls** (e.g. GM1100 is
       mislabeled "Liver cancer" in FinaleDB — its fragmentation is
       nothing like in-vivo plasma cfDNA and inflates the cancer
       set; the regex is from the cfdna-fragmentomics skill).
    3. **Writes a normalized cohort** to
       ``results/local_cohort/<cohort_name>/`` in the SAME shape the
       cross-study benchmark (``scripts/cross_study_finallydb.py``)
       expects: a features directory + a labels TSV with the canonical
       ``sample / disease_class / label / study`` columns.  When
       ``--merge-with`` is passed the adapter extends an existing
       labels file (deepcatch or cfdna-fragmentomics-pipeline) with
       the local samples and runs the cross-study benchmark on the
       combined cohort.
    4. **Emits a JSON provenance manifest** describing what was
       processed: input counts, dropped counts, kept counts, channel
       layout, and the resolved ``--merge-with`` target.

The adapter is intentionally **read-only on the input** — it never
mutates the collaborator's directory.  All writes go under
``results/local_cohort/<cohort_name>/`` so re-running produces a
deterministic, gitignored output.

Why a directory adapter, not a Python API?
    The whole point of the Collaborator Data Interface is that external
    collaborators don't have to know anything about the deepcatch /
    cfdna-fragmentomics-pipeline internals.  They produce a directory
    on disk that satisfies the schema; we do the conversion.  This
    script is the contract surface.

Usage:
    # 1) Sanity-check a cohort layout
    env -u PYTHONPATH ./.venv/bin/python \\
        scripts/adapter_local_cohort.py \\
        --cohort-root /path/to/cohort \\
        --dry-run

    # 2) Convert + validate + write normalized output
    env -u PYTHONPATH ./.venv/bin/python \\
        scripts/adapter_local_cohort.py \\
        --cohort-root /path/to/cohort \\
        --cohort-name my_lab_2026

    # 3) Merge with existing labels + run cross-study benchmark
    env -u PYTHONPATH ./.venv/bin/python \\
        scripts/adapter_local_cohort.py \\
        --cohort-root /path/to/cohort \\
        --cohort-name my_lab_2026 \\
        --merge-with ../cfdna-fragmentomics-pipeline/labels_multiclass.tsv \\
        --cross-study
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from typing import Dict, List, Optional, Tuple

import numpy as np


# Cell-line filter regex — keep in sync with
#   cfdna-fragmentomics skill: GM\d+, HeLa, HepG2, K562, etc.
#   scripts/cross_study_finallydb.py
CELL_LINE_RE = re.compile(
    r"^(GM\d+|HeLa|HepG2|K562|HL60|Jurkat|Raji|MCF7|U937|THP1|HEK293|"
    r"HCT116|SW480|A549|GM12878)",
    re.I,
)

# Channel contract — must match src/fragmentomics/tumor_naive_adapter.py.
# If this changes, both adapters and the cross-study benchmark break in
# lockstep; update all three.
CHANNEL_NAMES = [
    "delfi_5mb_ratio",
    "delfi_5mb_coverage",
    "delfi_100kb_ratio",
    "delfi_100kb_counts",
    "fsd_histogram",
]
CHANNEL_DIMS: Dict[str, int] = {
    "delfi_5mb_ratio": 631,
    "delfi_5mb_coverage": 631,
    "delfi_100kb_ratio": 30894,
    "delfi_100kb_counts": 30894,
    "fsd_histogram": 196,
}


# ──────────────────────────────────────────────────────────────────────
# Errors
# ──────────────────────────────────────────────────────────────────────

class SchemaError(ValueError):
    """Raised when a collaborator's cohort violates the schema.

    The error message is intentionally short and human-readable — the
    collaborator is NOT a deepcatch developer and shouldn't need to
    read Python tracebacks to fix their directory layout.
    """


# The canonical label vocabulary lives in the package, so every reader of
# a labels file in this repository resolves labels the same way. A label
# that is not a case or a control word is a schema error, not a control
# to be invented.
_SRC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from fragmentomics.labels import (  # noqa: E402
    LabelVocabularyError, looks_like_header, resolve_columns,
    resolve_label,
)


# ──────────────────────────────────────────────────────────────────────
# Validation
# ──────────────────────────────────────────────────────────────────────

def _required_paths(cohort_root: str) -> Tuple[str, str, str]:
    """Return (features_dir, labels_tsv, manifest_json) paths."""
    features_dir = os.path.join(cohort_root, "features")
    labels_tsv = os.path.join(cohort_root, "labels.tsv")
    manifest_json = os.path.join(cohort_root, "manifest.json")
    return features_dir, labels_tsv, manifest_json


def validate_cohort_layout(cohort_root: str) -> None:
    """Validate that the cohort root exists and has the expected files.

    Raises SchemaError with a human-readable message on any violation.
    """
    if not os.path.isdir(cohort_root):
        raise SchemaError(
            f"cohort_root does not exist or is not a directory: {cohort_root!r}. "
            f"Expected <cohort_root>/{{features/,labels.tsv,manifest.json}}."
        )

    features_dir, labels_tsv, manifest_json = _required_paths(cohort_root)
    if not os.path.isdir(features_dir):
        raise SchemaError(
            f"missing features/ subdirectory under {cohort_root!r}. "
            f"Each sample must have 5 per-channel .npy files in features/."
        )
    if not os.path.isfile(labels_tsv):
        raise SchemaError(
            f"missing labels.tsv at {labels_tsv!r}. "
            f"See docs/COLLABORATOR_DATA_INTERFACE.md for the column schema."
        )
    if not os.path.isfile(manifest_json):
        raise SchemaError(
            f"missing manifest.json at {manifest_json!r}. "
            f"The manifest records cohort metadata (lab, IRB, sample counts)."
        )


def validate_manifest(manifest_path: str) -> dict:
    """Validate manifest.json against the schema.

    Returns the parsed dict on success; raises SchemaError on missing
    fields.  Only fields that affect downstream behavior are required;
    optional fields (citation, license) are NOT enforced.
    """
    with open(manifest_path) as f:
        m = json.load(f)
    required_top = {"cohort_name", "channels_present", "vector_lengths"}
    missing = sorted(required_top - set(m))
    if missing:
        raise SchemaError(
            f"manifest.json missing required fields: {missing}. "
            f"See docs/COLLABORATOR_DATA_INTERFACE.md §manifest.json."
        )

    channels = m["channels_present"]
    if not isinstance(channels, list) or not channels:
        raise SchemaError(
            "manifest.json channels_present must be a non-empty list of "
            "channel names. Got: " + repr(channels)
        )
    vec_lengths = m["vector_lengths"]
    if not isinstance(vec_lengths, dict):
        raise SchemaError(
            "manifest.json vector_lengths must be a dict of "
            "channel_name -> int. Got: " + repr(vec_lengths)
        )
    # Warn (but don't reject) if collaborator included channels we don't
    # know how to load.  The cross-study benchmark requires all 5.
    unknown = sorted(set(channels) - set(CHANNEL_NAMES))
    if unknown:
        raise SchemaError(
            f"manifest.json channels_present includes unknown channels: "
            f"{unknown}. Expected subset of {CHANNEL_NAMES}."
        )
    return m


def validate_features_directory(
    features_dir: str,
    expected_samples: set[str],
    vec_lengths: dict,
) -> Tuple[List[str], List[str], List[str]]:
    """Walk features/, check that every expected sample has all channels
    present and the right vector lengths.

    Returns (missing_files, mismatched_lengths, sample_with_cell_line_id).
    Each entry is a human-readable string ready for the collaborator.
    """
    missing: List[str] = []
    mismatched: List[str] = []
    cell_line_samples: List[str] = []

    for sid in sorted(expected_samples):
        for ch in CHANNEL_NAMES:
            path = os.path.join(features_dir, f"{sid}.{ch}.npy")
            if not os.path.isfile(path):
                missing.append(f"{sid}.{ch}.npy")
                continue
            try:
                arr = np.load(path)
            except Exception as e:  # corrupted .npy
                mismatched.append(f"{sid}.{ch}.npy (load error: {e})")
                continue
            expected_dim = vec_lengths.get(ch, CHANNEL_DIMS[ch])
            if arr.shape != (expected_dim,):
                mismatched.append(
                    f"{sid}.{ch}.npy (shape {arr.shape}, expected ({expected_dim},))"
                )
        if CELL_LINE_RE.match(sid):
            cell_line_samples.append(sid)
    return missing, mismatched, cell_line_samples


# ──────────────────────────────────────────────────────────────────────
# Labels TSV
# ──────────────────────────────────────────────────────────────────────

def load_labels_tsv(path: str) -> List[dict]:
    """Load a collaborator labels.tsv into a list of dict rows.

    Required columns: sample_id, disease_class, label.
    Optional columns: study_id, publication_id, tissue_of_origin.

    Returns a list of dicts with the canonical keys.  Empty rows are
    skipped.  The header is detected by name match on the first
    required column — older files without headers are tolerated as long
    as they have >=3 tab-separated columns.
    """
    rows: List[dict] = []
    required = ("sample_id", "disease_class", "label")
    optional = ("study_id", "publication_id", "tissue_of_origin")
    with open(path) as f:
        headered = None
        label_col = study_col = None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            if headered is None:
                # The header is decided ONCE, from the first row, and the
                # column indices come from it. Deciding per row dropped
                # any sample whose id happened to be "sample"; asking a
                # DATA row for its column names picked up its disease
                # class as the label.
                headered = looks_like_header(parts)
                label_col, study_col = resolve_columns(parts,
                                                       headered=headered)
                if headered:
                    continue
            if label_col is None or label_col >= len(parts):
                continue
            row = {
                "sample_id": parts[0],
                "disease_class": parts[1] if len(parts) > 1 else "",
            }
            try:
                row["label"] = "cancer" if resolve_label(parts[label_col],
                                                         parts[0]) else "healthy"
            except LabelVocabularyError as exc:
                raise SchemaError(str(exc)) from exc
            if study_col is not None and study_col < len(parts):
                row["study_id"] = parts[study_col]
            if len(parts) > 4:
                row["publication_id"] = parts[4]
            if len(parts) > 5:
                row["tissue_of_origin"] = parts[5]
            rows.append(row)
    return rows


def write_labels_tsv(path: str, rows: List[dict]) -> None:
    """Write the canonical 4-column labels TSV the cross-study benchmark
    reads (``sample / disease_class / label / study``).

    The collaborator interface allows optional columns; the cross-study
    benchmark only needs the first 4.  Extra columns are preserved if
    they were present in the input, but the first 4 are guaranteed.
    """
    cols = ["sample_id", "disease_class", "label", "study_id",
            "publication_id", "tissue_of_origin"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            line = "\t".join(str(r.get(c, "")) for c in cols)
            f.write(line + "\n")


# ──────────────────────────────────────────────────────────────────────
# Conversion
# ──────────────────────────────────────────────────────────────────────

# The fragment-size histogram: 196 bins of 5 bp, from 20 bp to 1000 bp.
_FSD_BIN_LOW = 20
_FSD_BIN_WIDTH = 5
_FSD_N_BINS = 196


def _write_fsd_json(features_dir: str, sample_id: str,
                    histogram_path: str) -> None:
    """Write `<sample>.fsd.json` from the collaborator's histogram .npy.

    The benchmark's feature loaders -- the sibling pipeline's
    `honest_benchmark.load5` and the repo-local fallback in
    `scripts/load5_with_motifs.py` -- read the FSD channel from
    `<sample>.fsd.json`, a JSON file with a `size_bins` map. The
    collaborator interface documented here asks for
    `<sample>.fsd_histogram.npy` instead, because a laboratory can
    export a normalised histogram far more easily than it can produce
    the pipeline's JSON.

    So the adapter writes both: the .npy the interface specifies, and
    the .json the benchmark can read. Without it a converted cohort is
    invisible to the benchmark, because a sample with no readable FSD
    is dropped by its loader -- silently, as an empty cohort.
    """
    import numpy as np

    hist = np.load(histogram_path).astype(float).ravel()
    if hist.size != _FSD_N_BINS:
        raise ValueError(
            f"{sample_id}: the FSD histogram has {hist.size} bins, expected "
            f"{_FSD_N_BINS} (5 bp bins from {_FSD_BIN_LOW} bp to "
            f"{_FSD_BIN_LOW + _FSD_BIN_WIDTH * _FSD_N_BINS} bp)"
        )
    total = float(hist.sum())
    if total <= 0:
        raise ValueError(f"{sample_id}: the FSD histogram sums to {total}")

    edges = [_FSD_BIN_LOW + _FSD_BIN_WIDTH * i for i in range(_FSD_N_BINS)]
    mass = hist / total
    size_bins = {f"{lo}-{lo + _FSD_BIN_WIDTH}": float(m)
                 for lo, m in zip(edges, mass)}

    def quantile(q: float) -> float:
        """The bin's upper edge where the cumulative mass crosses q."""
        cumulative = 0.0
        for (lo, hi), m in zip(zip(edges, edges[1:] + [edges[-1] + _FSD_BIN_WIDTH]), mass):
            cumulative += m
            if cumulative >= q:
                return float(hi)
        return float(edges[-1] + _FSD_BIN_WIDTH)

    def fraction(lo_bp: float, hi_bp: float) -> float:
        total_in_band = 0.0
        for (lo, hi), m in zip(zip(edges, edges[1:] + [edges[-1] + _FSD_BIN_WIDTH]), mass):
            overlap = min(hi, hi_bp) - max(lo, lo_bp)
            if overlap > 0:
                total_in_band += m * (overlap / (hi - lo))
        return float(total_in_band)

    mode_index = int(np.argmax(mass))
    short_fraction = fraction(100.0, 150.0)
    long_fraction = fraction(150.0, 220.0)
    payload = {
        "sample": sample_id,
        "mode": "frag",
        # The collaborator's histogram is normalised, so the fragment
        # count behind it is not recoverable from the file.
        "fragment_count": None,
        "median_length": quantile(0.5),
        "mean_length": float(sum(m * (lo + _FSD_BIN_WIDTH / 2)
                                 for lo, m in zip(edges, mass))),
        "mode_length": float(edges[mode_index] + _FSD_BIN_WIDTH / 2),
        "p10": quantile(0.10),
        "p25": quantile(0.25),
        "p75": quantile(0.75),
        "p90": quantile(0.90),
        "short_fraction_100_150": short_fraction,
        "long_fraction_150_220": long_fraction,
        "short_long_ratio": (short_fraction / long_fraction
                              if long_fraction > 0 else None),
        "size_bins": size_bins,
        "source": "collaborator fsd_histogram.npy, normalised to sum=1",
    }
    with open(os.path.join(features_dir, f"{sample_id}.fsd.json"), "w") as f:
        json.dump(payload, f, indent=2)


def convert_cohort(
    cohort_root: str,
    cohort_name: str,
    out_root: str,
    drop_cell_lines: bool = True,
    dry_run: bool = False,
) -> dict:
    """Validate + convert one cohort directory.

    Returns a dict with provenance: n_in, n_out, n_filtered_cell_line,
    n_missing_features, channel_layout.  When ``dry_run`` is True the
    function still validates the cohort and computes the report, but
    does NOT create the normalized cohort directory on disk.  This is
    the contract used by the CLI's ``--dry-run`` flag.
    """
    validate_cohort_layout(cohort_root)
    features_dir, labels_tsv, manifest_json = _required_paths(cohort_root)
    manifest = validate_manifest(manifest_json)

    # Load TSV rows; preserve original order.
    rows = load_labels_tsv(labels_tsv)
    if not rows:
        raise SchemaError(
            f"labels.tsv at {labels_tsv!r} has no data rows."
        )

    expected_samples = {r["sample_id"] for r in rows}
    vec_lengths = manifest["vector_lengths"]
    missing, mismatched, cell_line_samples = validate_features_directory(
        features_dir, expected_samples, vec_lengths
    )

    if missing or mismatched:
        problems = missing + mismatched
        # Show at most 10 problems in the error to keep it readable.
        show = problems[:10]
        more = "" if len(problems) <= 10 else (
            f"\n  ... and {len(problems) - 10} more"
        )
        raise SchemaError(
            f"schema validation failed for {cohort_root!r}: "
            f"{len(problems)} problem(s).\n  - " + "\n  - ".join(show) + more
        )

    # Filter cell lines (and warn about which ones were dropped).
    dropped_ids: List[str] = []
    if drop_cell_lines:
        keep_rows = []
        for r in rows:
            if CELL_LINE_RE.match(r["sample_id"]):
                dropped_ids.append(r["sample_id"])
                continue
            keep_rows.append(r)
    else:
        keep_rows = list(rows)

    # Map manifest fields onto the canonical output labels shape.
    out_rows: List[dict] = []
    for r in keep_rows:
        # The cross-study benchmark's labels loader (in
        # scripts/cross_study_finallydb.py) expects the third column
        # to be 'label' (cancer/healthy).  It uses 'study' as the
        # 4th column for per-study harmonization — the local cohort's
        # study_id maps there.
        out_rows.append({
            "sample_id": r["sample_id"],
            "disease_class": r["disease_class"],
            "label": r["label"],
            "study_id": r.get("study_id", "") or cohort_name,
            "publication_id": r.get("publication_id", ""),
            "tissue_of_origin": r.get("tissue_of_origin", ""),
        })

    # Write the normalized cohort (skipped when dry_run).
    target_dir = os.path.join(out_root, cohort_name)
    target_features = os.path.join(target_dir, "features")

    if not dry_run:
        os.makedirs(target_features, exist_ok=True)

        # Symlink the per-sample .npy files (read-only — never copy/mutate
        # the collaborator's originals).  If symlinks are unsupported on
        # the platform, fall back to copy.
        for sid in sorted({r["sample_id"] for r in out_rows}):
            for ch in CHANNEL_NAMES:
                src = os.path.join(features_dir, f"{sid}.{ch}.npy")
                dst = os.path.join(target_features, f"{sid}.{ch}.npy")
                if os.path.lexists(dst):
                    continue
                try:
                    os.symlink(os.path.abspath(src), dst)
                except (OSError, NotImplementedError):
                    import shutil
                    shutil.copyfile(src, dst)
            # The benchmark reads the FSD channel from a .json, so the
            # converted cohort carries that too (see _write_fsd_json).
            _write_fsd_json(
                target_features, sid,
                os.path.join(features_dir, f"{sid}.fsd_histogram.npy"),
            )

        # Canonical 4-col labels.tsv for the cross-study benchmark.
        write_labels_tsv(
            os.path.join(target_dir, "labels.tsv"),
            out_rows,
        )

    # Adapter output manifest — separate from the collaborator's
    # input manifest.  The output manifest records what we did
    # (provenance), not what they told us about their cohort.
    out_manifest = {
        "schema_version": "1.0",
        "cohort_name": cohort_name,
        "source_cohort_root": os.path.abspath(cohort_root),
        "source_manifest": manifest,
        "n_samples_in": len(rows),
        "n_samples_out": len(out_rows),
        "n_dropped_cell_line": len(dropped_ids),
        "dropped_cell_line_ids": dropped_ids,
        "n_missing_features": len(missing),
        "n_mismatched_lengths": len(mismatched),
        "channel_layout": list(CHANNEL_NAMES),
        "channel_dims": CHANNEL_DIMS,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scope": "local-cohort adapter output. NOT a clinical-validation "
                 "cohort. See research_use_only: true.",
        "research_use_only": True,
        "dry_run": dry_run,
    }
    if not dry_run:
        with open(os.path.join(target_dir, "adapter_manifest.json"), "w") as f:
            json.dump(out_manifest, f, indent=2)

    return out_manifest


# ──────────────────────────────────────────────────────────────────────
# Merge + cross-study entry points
# ──────────────────────────────────────────────────────────────────────

def merge_with_existing_labels(
    local_labels_tsv: str,
    existing_labels_tsv: str,
    out_path: str,
) -> dict:
    """Combine local + existing labels into one TSV.

    The cross-study benchmark reads the standard
    ``sample / disease_class / label / study`` columns.  Both input
    files must already be in that shape — use the conversion step
    above first if the local file is in the longer collaborator
    schema.

    Returns a dict with merge provenance (counts, source files).
    """
    local_rows = load_labels_tsv(local_labels_tsv)
    existing_rows = load_labels_tsv(existing_labels_tsv)

    # Preserve unique sample IDs (the cross-study benchmark loads by
    # sample id, so duplicates would silently double-count).
    seen: set[str] = set()
    merged: List[dict] = []
    dropped_dups = 0
    for r in local_rows + existing_rows:
        sid = r["sample_id"]
        if sid in seen:
            dropped_dups += 1
            continue
        seen.add(sid)
        merged.append(r)

    write_labels_tsv(out_path, merged)
    return {
        "n_local": len(local_rows),
        "n_existing": len(existing_rows),
        "n_merged": len(merged),
        "n_dropped_duplicates": dropped_dups,
        "out_path": out_path,
    }


def _features_dir_for_labels(labels_tsv: str) -> str:
    """Where the feature files of an existing labels TSV live.

    The pipeline repository keeps its labels at the top level
    (`labels_multiclass.tsv`) and its feature cache under
    `data/features/`. The adapter used to pass the *labels'* directory,
    which is one level too high: the benchmark's loader looks for
    `<features_dir>/<sample_id>.<channel>.npy`, found nothing, and
    failed with "Found array with 0 sample(s)" — an error that says
    nothing about the cause. Both layouts are tried.
    """
    root = os.path.dirname(labels_tsv)
    for candidate in (os.path.join(root, "data", "features"), root):
        if os.path.isdir(candidate):
            return candidate
    return os.path.join(root, "data", "features")


def _run_cross_study(merged_labels_tsv: str, merged_features_dir: str,
                     out_json: str) -> int:
    """Invoke scripts/cross_study_finallydb.py on the merged cohort.

    Returns the cross-study script's exit code (0 on success).
    Silently passes --publications "8" (Cristiano 2019) since the
    local cohort does not have a FinaleDB publication id — the
    benchmark will pick the publication column from the merged TSV.
    """
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    script = os.path.join(repo_root, "scripts", "cross_study_finallydb.py")
    import subprocess
    out_dir = os.path.dirname(os.path.abspath(out_json)) or "."
    cmd = [
        sys.executable, script,
        "--labels-multiclass", merged_labels_tsv,
        "--features-dir", merged_features_dir,
        "--out-json", out_json,
        "--publications", "8",
        # The benchmark writes a benchmark document and a per-cancer
        # artifact by default, and those two are the repository's
        # PUBLISHED files — the ones the manuscript quotes. A
        # collaborator's optional run must not overwrite them with a
        # single-publication result, so both go to the adapter's own
        # output directory (the per-cancer artifact is skipped: it
        # belongs to the published run, not to this one).
        "--out-md", os.path.join(out_dir, "CROSS_STUDY_BENCHMARK.md"),
        "--out-per-cancer-json", "",
    ]
    print(f"[adapter_local_cohort] running: {' '.join(cmd)}")
    return subprocess.call(cmd)


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────

def _build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--cohort-root", required=True,
        help="Path to the collaborator cohort root (features/, labels.tsv, "
             "manifest.json).",
    )
    ap.add_argument(
        "--cohort-name", default="local_cohort",
        help="Name used for the output subdirectory under results/local_cohort/. "
             "Default: 'local_cohort'.",
    )
    ap.add_argument(
        "--out-root", default="results/local_cohort",
        help="Root directory for the normalized cohort output. Default: "
             "results/local_cohort (relative to the deepcatch repo root).",
    )
    ap.add_argument(
        "--no-cell-line-filter", action="store_true",
        help="By default, samples whose IDs match the cell-line regex "
             "(GM\\d+, HeLa, HepG2, K562, ...) are dropped. Pass this flag "
             "to KEEP them — useful when validating the filter itself.",
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="Validate the cohort and print a JSON report to stdout, but "
             "don't write any output files.",
    )
    ap.add_argument(
        "--merge-with", default=None,
        help="Optional path to an existing labels TSV (e.g. cfdna-fragmentomics-pipeline/"
             "labels_multiclass.tsv). When passed, the adapter writes a combined "
             "labels file under results/local_cohort/<cohort_name>/merged_labels.tsv "
             "and (if --cross-study) invokes the cross-study benchmark on it.",
    )
    ap.add_argument(
        "--cross-study", action="store_true",
        help="Run scripts/cross_study_finallydb.py on the merged cohort. "
             "Implies --merge-with (the merged labels TSV is the input).",
    )
    ap.add_argument(
        "--cross-study-out",
        default="results/cross_study_local_cohort.json",
        help="Output JSON for the cross-study benchmark when --cross-study is set.",
    )
    return ap


def main() -> int:
    ap = _build_argparser()
    args = ap.parse_args()

    if args.cross_study and not args.merge_with:
        ap.error("--cross-study requires --merge-with")

    cohort_root = os.path.abspath(args.cohort_root)
    out_root = os.path.abspath(args.out_root)
    target_dir = os.path.join(out_root, args.cohort_name)

    try:
        report = convert_cohort(
            cohort_root=cohort_root,
            cohort_name=args.cohort_name,
            out_root=out_root,
            drop_cell_lines=not args.no_cell_line_filter,
            dry_run=args.dry_run,
        )
    except SchemaError as e:
        # Print to stderr and exit non-zero — the collaborator gets a
        # human-readable error and a clear exit code.
        sys.stderr.write(f"schema error: {e}\n")
        return 2

    if args.dry_run:
        report["would_write_to"] = target_dir
        print(json.dumps(report, indent=2))
        return 0

    print(f"[adapter_local_cohort] wrote normalized cohort to {target_dir}/")
    print(f"[adapter_local_cohort] n_in={report['n_samples_in']} "
          f"n_out={report['n_samples_out']} "
          f"n_dropped_cell_line={report['n_dropped_cell_line']}")

    # Optional merge + cross-study.
    if args.merge_with:
        local_labels = os.path.join(target_dir, "labels.tsv")
        existing_labels = os.path.abspath(args.merge_with)
        merged_path = os.path.join(target_dir, "merged_labels.tsv")
        merge_report = merge_with_existing_labels(
            local_labels, existing_labels, merged_path,
        )
        print(f"[adapter_local_cohort] merged with {existing_labels}: "
              f"n_local={merge_report['n_local']} "
              f"n_existing={merge_report['n_existing']} "
              f"n_merged={merge_report['n_merged']}")
        report["merge"] = merge_report

        if args.cross_study:
            # The merged TSV lives in the deepcatch results dir, but
            # the features still live under the cfdna-fragmentomics-pipeline
            # features dir — the cross-study benchmark's load5() reads by
            # sample id and looks for <features_dir>/<sample>.delfi_*.npy.
            # For the merged cohort, point the benchmark at the
            # original features dir of the existing labels, since the
            # local cohort's features are not visible to load5's glob.
            rc = _run_cross_study(
                merged_labels_tsv=merged_path,
                merged_features_dir=_features_dir_for_labels(existing_labels),
                out_json=os.path.abspath(args.cross_study_out),
            )
            report["cross_study_exit_code"] = rc
            # Say plainly what that run did and did not include: the
            # benchmark selects publications from the FinaleDB registry
            # (--publications 8), so a study with no FinaleDB
            # publication id is dropped from it. The normalised cohort
            # above is the deliverable; this run is the open-data
            # benchmark over the merged label file.
            report["cross_study_scope"] = (
                "open-data publications only (--publications 8); the local "
                "study's samples are dropped by the benchmark's publication "
                "filter, so they are not part of this run"
            )
            print("[adapter_local_cohort] NOTE: that run covers the "
                  "open-data publications only -- the benchmark selects "
                  "publications from the FinaleDB registry, so the local "
                  "study's samples are not part of it. The normalised "
                  "cohort above is the deliverable.")

    # Always emit the final report to stdout as JSON (one line so it is
    # easy to grep; the JSON contains all the provenance).
    print("---ADAPTER_REPORT---")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())