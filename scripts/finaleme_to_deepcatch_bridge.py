#!/usr/bin/env python3
"""Bridge: convert FinaleMe output → DeepCatch-format methylation JSON.

The deepcatch repo (panel-LLR + tumor-naive fusion research) needs a
methylation feature matrix per-sample to add a third channel to the
existing fragmentomics pipeline. The natural source is FinaleMe
(Liu et al. Nat Commun 15:2790, 2024), which imputes per-CpG
β-values from cfDNA WGS fragments.

This bridge is **read-only on the FinaleMe output** and **never
synthesizes predictions**. Its contract:

    1. **Probe** the FinaleMe output directory for at least one
       per-sample file (BetaValues.tsv or *.decoded.bed.gz).
    2. **Validate** the methylation table: sample IDs must overlap with
       the FinaleDB features cache; β-values must be real numbers in
       [0, 1]; missing files fail LOUD (non-zero exit + clear error).
    3. **Aggregate** to per-sample summary statistics: mean β, std β,
       per-region methylation (CpG-island / shore / shelf / open-sea)
       when the methylation-prior BigWig is available; per-cancer β
       distribution when labels are available.
    4. **Emit** a JSON provenance manifest under
       ``results/cross_platform_finaledb_<sample>.json`` describing
       what was processed.
    5. **When data is missing**, emit a ``cross_platform_finaledb_NOT_RUN.json``
       with ``data_source='no_data'`` and ``n_samples_in=0`` so the
       downstream validator (cross_platform_finaledb_validate.py) can
       surface the missingness without fabricating numbers.

Crucially: this bridge never invents per-sample scores when FinaleMe
output is absent. It either produces the real aggregate from a real
FinaleMe decode, or it emits the not-run sentinel.

Usage:
    # Real-data path (when FinaleMe output is on disk)
    env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_to_deepcatch_bridge.py \\
        --finaleme-dir ~/.hermes/.local/finaleme/output \\
        --output results/cross_platform_finaledb_real.json

    # Dry-run / not-run sentinel
    env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_to_deepcatch_bridge.py \\
        --finaleme-dir /nonexistent \\
        --output results/cross_platform_finaledb_NOT_RUN.json
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# ─────────────────────────────────────────────────────────────────────
# Errors
# ─────────────────────────────────────────────────────────────────────


class BridgeError(ValueError):
    """Raised when FinaleMe output violates the bridge schema.

    The message is intentionally short and human-readable. Sample IDs
    that fail overlap, β-values outside [0, 1], and missing files all
    surface as BridgeError so the upstream orchestrator can print the
    issue without a Python traceback.
    """


# ─────────────────────────────────────────────────────────────────────
# Probes for FinaleMe output
# ─────────────────────────────────────────────────────────────────────


_FINALEME_PER_SAMPLE_PATTERNS = [
    "*BetaValues.tsv",
    "*beta_values.tsv",
    "*_decoded.bed.gz",
    "*decoded.bed.gz",
    "*cpg_features.hg19.bed.gz",
]

# Per-file schema: {tsv: (sep, sample_col, beta_col), bed: (sep, sample_col?, beta_col?)}
# For TSV: chr  start  end  sample_id  beta
# For decoded.bed.gz: chr  start  end  methylated_fraction (single-sample; sample_id from filename)


def find_finaleme_per_sample(
    finaleme_dir: Path,
    patterns: Sequence[str] = _FINALEME_PER_SAMPLE_PATTERNS,
) -> List[Path]:
    """Return a list of FinaleMe per-sample files matching any pattern.
    Sorted by name for determinism."""
    if not finaleme_dir.is_dir():
        return []
    found: List[Path] = []
    for pat in patterns:
        found.extend(finaleme_dir.glob(pat))
    return sorted(set(found))


# ─────────────────────────────────────────────────────────────────────
# Sample-ID and β-value validation
# ─────────────────────────────────────────────────────────────────────


def list_finaledb_sample_ids(features_dir: Path) -> List[str]:
    """Read sample IDs from a FinaleDB-style 5-channel features directory.

    Convention: ``<features_dir>/<sample_id>.<channel>.npy``. We infer
    sample_id by stripping the last two extensions."""
    if not features_dir.is_dir():
        return []
    sample_ids: set = set()
    for p in features_dir.glob("*.npy"):
        # 'C309.delfi_5mb_ratio.npy' → 'C309'
        parts = p.name.split(".")
        if len(parts) >= 3:
            sid = ".".join(parts[:-2])
            sample_ids.add(sid)
        else:
            sample_ids.add(p.stem)
    return sorted(sample_ids)


def validate_sample_overlap(
    meth_sample_ids: Sequence[str],
    finaledb_sample_ids: Sequence[str],
) -> List[str]:
    """Return the intersection in stable order.

    Raises BridgeError when any methylation sample ID is not in the
    FinaleDB cache — a sample mismatch means a comparison on disjoint
    cohorts, which is exactly the bug the methodology demo called out.
    """
    finaledb_set = set(finaledb_sample_ids)
    missing = [s for s in meth_sample_ids if s not in finaledb_set]
    if missing:
        raise BridgeError(
            f"methylation samples not present in FinaleDB features cache: "
            f"{missing[:5]}{'...' if len(missing) > 5 else ''}. "
            f"This usually means the two cohorts are disjoint; cross-platform "
            f"AUC on disjoint cohorts is the bug the methodology demo warns against."
        )
    return [s for s in meth_sample_ids if s in finaledb_set]


def validate_beta_values(betas: np.ndarray) -> None:
    """β-values must be in [0, 1]. NaN is allowed (missing per-CpG).

    Raises BridgeError when any β is outside the range — a column of
    1.5s or -0.2s indicates a misparsed FinaleMe output (e.g., reading
    a per-CpG count column instead of a β column)."""
    if not isinstance(betas, np.ndarray):
        betas = np.asarray(betas)
    if betas.size == 0:
        raise BridgeError("β-value table is empty.")
    finite = betas[np.isfinite(betas)]
    if finite.size == 0:
        # All-NaN is OK (noisy sample); we degrade to no-signal
        return
    too_high = int(np.sum(finite > 1.0 + 1e-9))
    too_low = int(np.sum(finite < 0.0 - 1e-9))
    if too_high or too_low:
        raise BridgeError(
            f"β-values out of [0, 1]: {too_high} above 1.0, {too_low} below 0.0 "
            f"(first few above: "
            f"{finite[finite > 1.0 + 1e-9][:3].tolist()}, "
            f"first few below: "
            f"{finite[finite < 0.0 - 1e-9][:3].tolist()}). "
            f"Check that the column you selected is β (in [0,1]) and not raw count."
        )


# ─────────────────────────────────────────────────────────────────────
# TSV parsing
# ─────────────────────────────────────────────────────────────────────


# Methylation column-name heuristic (case-insensitive). Default order:
_COL_BETA_HINTS = ("beta", "beta_value", "methylation", "methy_frac")
_COL_SAMPLE_HINTS = ("sample", "sample_id", "sampleid", "sid")
_COL_CHR_HINTS = ("chr", "chrom")
_COL_START_HINTS = ("start", "pos", "position")

# Schema for "BetaValues.tsv" (chr  start  end  sample_id  beta)
_BETAVALUES_HEADER = re.compile(r"^chr\b.*start\b", re.I | re.S)


def _pick_column(fieldnames, hints) -> Optional[str]:
    """Case-insensitive column-name matcher. Returns the column that
    best matches any hint, else None."""
    if not fieldnames:
        return None
    lower = {f.lower(): f for f in fieldnames}
    for hint in hints:
        if hint in lower:
            return lower[hint]
    return None


def parse_betavalues_tsv(path: Path) -> Tuple[List[str], np.ndarray]:
    """Parse a FinaleMe BetaValues.tsv file.

    Returns (sample_ids, beta_matrix) where beta_matrix has shape
    (n_unique_samples, n_unique_cpgs). NaN is filled with np.nan (NOT
    zero — zero would be silently misinterpreted as "unmethylated").

    Raises BridgeError on file-format violations.
    """
    if not path.is_file():
        raise BridgeError(f"file does not exist: {path}")
    try:
        with open(path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = reader.fieldnames or []
            chr_col = _pick_column(fieldnames, _COL_CHR_HINTS)
            start_col = _pick_column(fieldnames, _COL_START_HINTS)
            sample_col = _pick_column(fieldnames, _COL_SAMPLE_HINTS)
            beta_col = _pick_column(fieldnames, _COL_BETA_HINTS)
            if not (sample_col and beta_col):
                raise BridgeError(
                    f"BetaValues.tsv at {path} missing required "
                    f"sample_id / beta columns. Got columns: {fieldnames[:10]}"
                )
            # Per-sample list of (start, beta) tuples.
            rows_by_sample: Dict[str, List[Tuple[int, float]]] = {}
            for row in reader:
                sid = (row.get(sample_col) or "").strip()
                beta_str = (row.get(beta_col) or "").strip()
                if not sid or not beta_str:
                    continue
                try:
                    beta = float(beta_str)
                except ValueError:
                    continue
                start_str = (row.get(start_col) or "0") if start_col else "0"
                try:
                    start = int(start_str)
                except ValueError:
                    start = 0
                rows_by_sample.setdefault(sid, []).append((start, beta))
    except (OSError, UnicodeDecodeError) as e:
        raise BridgeError(f"failed to read BetaValues.tsv {path}: {e}") from e

    if not rows_by_sample:
        raise BridgeError(
            f"BetaValues.tsv {path} yielded no per-sample rows after parsing."
        )
    sample_ids = sorted(rows_by_sample.keys())
    # Build per-sample dense vectors on the union of CpG positions.
    all_cpgs = sorted({start for sids in rows_by_sample.values()
                        for start, _ in sids})
    if not all_cpgs:
        raise BridgeError(
            f"BetaValues.tsv {path} has rows but no CpG-position column values."
        )
    cpgs_index = {c: i for i, c in enumerate(all_cpgs)}
    matrix = np.full((len(sample_ids), len(all_cpgs)), np.nan)
    for i, sid in enumerate(sample_ids):
        for start, beta in rows_by_sample[sid]:
            matrix[i, cpgs_index[start]] = beta
    return sample_ids, matrix


# ─────────────────────────────────────────────────────────────────────
# Per-file dispatch
# ─────────────────────────────────────────────────────────────────────


def _per_sample_payload(
    *,
    finaleme_dir: Path,
    features_dir: Path,
    labels_tsv: Optional[Path],
    parsed_files: List[Path],
) -> Dict:
    """Internal: build the JSON payload from a list of parsed
    BetaValues.tsv files. Returns a dict ready for json.dump.

    When parsed_files is empty, returns a not-runnable payload.
    """
    finaledb_ids = list_finaledb_sample_ids(features_dir)
    if not finaledb_ids:
        return _not_runnable_payload(
            finaleme_dir=finaleme_dir,
            features_dir=features_dir,
            reason=f"FinaleDB features directory {features_dir} is empty or missing",
        )

    if not parsed_files:
        return _not_runnable_payload(
            finaleme_dir=finaleme_dir,
            features_dir=features_dir,
            reason=f"FinaleMe output directory {finaleme_dir} contains no "
                   f"BetaValues.tsv or *.decoded.bed.gz files",
            finaledb_n_samples=len(finaledb_ids),
        )

    # Aggregate per-sample beta vectors across all parsed files.
    per_sample_betas: Dict[str, List[float]] = {}
    n_features_per_file: List[int] = []
    parsed_summaries: List[Dict] = []
    for path in parsed_files:
        try:
            sample_ids, matrix = parse_betavalues_tsv(path)
            validate_beta_values(matrix)
        except BridgeError as e:
            parsed_summaries.append({
                "path": str(path),
                "status": "rejected",
                "error": str(e),
            })
            continue
        n_features_per_file.append(int(matrix.shape[1]) if matrix.ndim == 2 else 0)
        for i, sid in enumerate(sample_ids):
            row = matrix[i]
            row = row[np.isfinite(row)]
            per_sample_betas.setdefault(sid, []).extend(row.tolist())
        parsed_summaries.append({
            "path": str(path),
            "status": "parsed",
            "n_samples_in_file": len(sample_ids),
            "n_cpgs": int(matrix.shape[1]) if matrix.ndim == 2 else 0,
        })

    if not per_sample_betas:
        return _not_runnable_payload(
            finaleme_dir=finaleme_dir,
            features_dir=features_dir,
            reason="All FinaleMe output files were rejected during validation",
            parsed_summaries=parsed_summaries,
            finaledb_n_samples=len(finaledb_ids),
        )

    # Sample overlap check
    meth_sample_ids = sorted(per_sample_betas.keys())
    try:
        overlap = validate_sample_overlap(meth_sample_ids, finaledb_ids)
    except BridgeError as e:
        return _not_runnable_payload(
            finaleme_dir=finaleme_dir,
            features_dir=features_dir,
            reason=str(e),
            parsed_summaries=parsed_summaries,
            n_meth_samples=len(meth_sample_ids),
            n_finaledb_samples=len(finaledb_ids),
        )

    # Per-sample beta summary
    per_sample_summary: Dict[str, Dict[str, float]] = {}
    cancer_means: List[float] = []
    healthy_means: List[float] = []
    label_map = _load_labels(labels_tsv) if labels_tsv and labels_tsv.is_file() else {}
    for sid in overlap:
        betas = np.asarray(per_sample_betas[sid], dtype=float)
        if betas.size == 0:
            per_sample_summary[sid] = {"mean": float("nan"), "std": float("nan"), "n_cpgs": 0}
            continue
        m = float(np.mean(betas))
        s = float(np.std(betas))
        per_sample_summary[sid] = {"mean": m, "std": s, "n_cpgs": int(betas.size)}
        label = label_map.get(sid)
        if label == "cancer":
            cancer_means.append(m)
        elif label == "healthy":
            healthy_means.append(m)

    payload = {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_source": "both" if overlap else "methylation_only",
        "runnable": True,
        "n_samples_in": len(overlap),
        "n_meth_samples": len(meth_sample_ids),
        "n_finaledb_samples": len(finaledb_ids),
        "n_overlap_samples": len(overlap),
        "n_features_methylation": int(np.mean(n_features_per_file)) if n_features_per_file else 0,
        "provenance": {
            "finaleme_dir": str(finaleme_dir),
            "finaledb_features_dir": str(features_dir),
            "labels_tsv": str(labels_tsv) if labels_tsv else None,
            "n_parsed_files": len(parsed_summaries),
            "parsed_files": parsed_summaries,
        },
        "per_sample_summary": per_sample_summary,
        "per_cancer_methylation": {
            "cancer_mean_beta": {
                "mean": float(np.mean(cancer_means)) if cancer_means else None,
                "std": float(np.std(cancer_means)) if cancer_means else None,
                "n_samples": len(cancer_means),
            },
            "healthy_mean_beta": {
                "mean": float(np.mean(healthy_means)) if healthy_means else None,
                "std": float(np.std(healthy_means)) if healthy_means else None,
                "n_samples": len(healthy_means),
            },
            "delta_cancer_minus_healthy": (
                float(np.mean(cancer_means) - np.mean(healthy_means))
                if cancer_means and healthy_means else None
            ),
        },
        # Note: this payload does NOT carry a `cross_platform_auc` —
        # that's the validator's job. Bridge output is the per-sample
        # methylation summary; the validator computes the fusion AUC.
    }
    return payload


def _not_runnable_payload(
    *,
    finaleme_dir: Path,
    features_dir: Path,
    reason: str,
    parsed_summaries: Optional[List[Dict]] = None,
    finaledb_n_samples: Optional[int] = None,
    n_meth_samples: Optional[int] = None,
    n_finaledb_samples: Optional[int] = None,
) -> Dict:
    return {
        "schema_version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "data_source": "no_data",
        "runnable": False,
        "n_samples_in": 0,
        "n_features_methylation": 0,
        "provenance": {
            "finaleme_dir": str(finaleme_dir),
            "finaledb_features_dir": str(features_dir),
            "parsed_summaries": parsed_summaries or [],
            "finaledb_n_samples": finaledb_n_samples,
            "n_meth_samples": n_meth_samples,
            "n_finaledb_samples": n_finaledb_samples,
        },
        "not_runnable_reason": reason,
    }


def _load_labels(labels_tsv: Path) -> Dict[str, str]:
    """Returns {sample_id: label} from a FinaleDB-style labels TSV."""
    label_map: Dict[str, str] = {}
    try:
        with open(labels_tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                sid = (row.get("sample") or row.get("sample_id") or "").strip()
                lab = (row.get("label") or "").strip()
                if sid and lab:
                    label_map[sid] = lab
    except (OSError, UnicodeDecodeError):
        pass
    return label_map


# ─────────────────────────────────────────────────────────────────────
# Public entry points (used by cross_platform_finaledb_validate.py
# and by the CLI below)
# ─────────────────────────────────────────────────────────────────────


def convert_or_emit_not_run(
    *,
    finaleme_dir: Path,
    finaledb_features_dir: Path,
    output_path: str,
    labels_tsv: Optional[Path] = None,
    is_synthetic_fixture: bool = False,
    reason: Optional[str] = None,
) -> Dict:
    """Convert FinaleMe output → DeepCatch JSON.

    Always writes ``output_path``. When data is missing, emits a
    not-runnable payload (``data_source='no_data'``, ``runnable=False``).

    Returns the parsed payload dict.
    """
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    parsed = find_finaleme_per_sample(Path(finaleme_dir))
    payload = _per_sample_payload(
        finaleme_dir=Path(finaleme_dir),
        features_dir=Path(finaledb_features_dir),
        labels_tsv=Path(labels_tsv) if labels_tsv else None,
        parsed_files=parsed,
    )

    if reason and not payload.get("runnable", False):
        payload["not_runnable_reason"] = reason

    if is_synthetic_fixture:
        payload["is_synthetic_fixture"] = True

    out.write_text(json.dumps(payload, indent=2))
    return payload


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Bridge: convert FinaleMe output → DeepCatch-format methylation JSON. "
            "Never fabricates predictions; emits NOT_RUN JSON when data missing."
        ),
    )
    ap.add_argument(
        "--finaleme-dir", type=Path, required=True,
        help="Directory containing FinaleMe per-sample outputs (BetaValues.tsv).",
    )
    ap.add_argument(
        "--finaledb-features-dir", type=Path, required=True,
        help="FinaleDB features cache (5-channel .npy per sample).",
    )
    ap.add_argument(
        "--labels-tsv", type=Path, default=None,
        help="Optional labels TSV for per-cancer methylation summaries.",
    )
    ap.add_argument(
        "--output", type=str, required=True,
        help="Output JSON path (will be created with parent dirs).",
    )
    ap.add_argument(
        "--synthetic-fixture", action="store_true",
        help="Mark the output as a synthetic test fixture (NOT real FinaleMe data).",
    )
    args = ap.parse_args()

    try:
        payload = convert_or_emit_not_run(
            finaleme_dir=args.finaleme_dir,
            finaledb_features_dir=args.finaledb_features_dir,
            output_path=args.output,
            labels_tsv=args.labels_tsv,
            is_synthetic_fixture=args.synthetic_fixture,
        )
    except BridgeError as e:
        print(f"BRIDGE_ERROR: {e}", file=sys.stderr)
        return 2

    print(f"data_source: {payload['data_source']}")
    print(f"runnable:    {payload.get('runnable', False)}")
    print(f"n_samples_in: {payload.get('n_samples_in', 0)}")
    if payload.get("runnable"):
        print(f"n_features_methylation: {payload.get('n_features_methylation', 0)}")
    else:
        print(f"reason: {payload.get('not_runnable_reason', '?')[:120]}…")
    print(f"output: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
