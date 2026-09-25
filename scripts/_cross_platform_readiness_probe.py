#!/usr/bin/env python3
"""Cross-platform readiness probe.

Prints the current state of the cross-platform pipeline inputs:

* FinaleMe JAR (finaledb-finaleme / Zenodo hybrid recipe in
  ``~/.hermes/.local/finaleme/``)
* FinaleMe pretrained models (``~/.hermes/.local/finaleme/models/``)
* FinaleMe reference files (1.4 GB total: hg19.2bit, methylation
  prior, CpG motif, mappability BED)
* FinaleMe output directory (per-sample BetaValues.tsv files)
* FinaleDB features cache (the local on-disk 5-channel features)

This script is a PROBE — it ALWAYS exits 0. It is intentionally
non-destructive. Operators read its output to decide whether to
install the missing components. See
``docs/CROSS_PLATFORM_FINALEME.md`` for the runnable state recipe.

Usage:
    env -u PYTHONPATH ./.venv/bin/python scripts/_cross_platform_readiness_probe.py
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Tuple

# ─────────────────────────────────────────────────────────────────────
# Standard paths (deepcatch convention: ~/.hermes/.local/finaleme/*)
# ─────────────────────────────────────────────────────────────────────

FINALEME_HOME = Path("~/.hermes/.local/finaleme").expanduser()

# Candidate JAR locations: try the canonical production path first,
# then the hybrid JAR that the cfdna-fragmentomics skill's recipe builds.
FINALEME_JAR_CANDIDATES: List[Path] = [
    FINALEME_HOME / "FinaleMe-0.61-jar-with-dependencies.jar",
    FINALEME_HOME / "FinaleMe-hybrid-jar-with-dependencies.jar",
    FINALEME_HOME / "FinaleMe-0.58.1-jar-with-dependencies.jar",
]

# Pretrained models
FINALEME_MODELS_DIR = FINALEME_HOME / "models"
FINALEME_PRETRAINED_FILES: List[str] = [
    "healthy_WGS.mincg7.example.hmm_model",
    "cancer_WGS.mincg7.example.hmm_model",
]

# Reference files (1.4 GB total)
FINALEME_REF_FILES: List[Tuple[str, Path]] = [
    ("hg19.2bit",                       FINALEME_HOME / "hg19.2bit"),
    ("methylation prior bw",            FINALEME_HOME / "wgbs_buffyCoat_jensen2015GB.methy.hg19.bw"),
    ("CpG motif bedgraph",              FINALEME_HOME / "CG_motif.hg19.common_chr.pos_only.bedgraph.gz"),
    ("mappability BED",                 FINALEME_HOME / "mappability.bed"),
    ("chrom sizes",                     FINALEME_HOME / "hg19.chrom.sizes"),
]

# FinaleMe output directory (per-sample BetaValues.tsv / bed.gz)
DEFAULT_FINALEME_OUTPUT_DIR = FINALEME_HOME / "output"

# FinaleDB features cache (deepcatch's cfdna-fragmentomics sibling)
DEFAULT_FINALEDB_FEATURES_DIR = Path(
    "/Users/hermes/cfdna-fragmentomics-pipeline/data/features"
)

# ─────────────────────────────────────────────────────────────────────
# Color helpers (with safe fallback if NO_COLOR is set or no TTY)
# ─────────────────────────────────────────────────────────────────────


def _green(s: str) -> str:
    if not sys.stdout.isatty() or os.environ.get("NO_COLOR"):
        return s
    return f"\033[32m{s}\033[0m"


def _red(s: str) -> str:
    if not sys.stdout.isatty() or os.environ.get("NO_COLOR"):
        return s
    return f"\033[31m{s}\033[0m"


def _yellow(s: str) -> str:
    if not sys.stdout.isatty() or os.environ.get("NO_COLOR"):
        return s
    return f"\033[33m{s}\033[0m"


def _bold(s: str) -> str:
    if not sys.stdout.isatty() or os.environ.get("NO_COLOR"):
        return s
    return f"\033[1m{s}\033[0m"


# ─────────────────────────────────────────────────────────────────────
# Probe helpers
# ─────────────────────────────────────────────────────────────────────


def _check_jar() -> Tuple[bool, str]:
    """Return (found, detail) for any FinaleMe JAR candidate."""
    for cand in FINALEME_JAR_CANDIDATES:
        if cand.is_file():
            size_mb = cand.stat().st_size / 1024 / 1024
            return True, f"{cand} ({size_mb:.1f} MB)"
    return False, (
        f"not found. Tried {len(FINALEME_JAR_CANDIDATES)} candidate path(s) "
        f"under {FINALEME_HOME}. See docs/CROSS_PLATFORM_FINALEME.md."
    )


def _check_pretrained() -> Tuple[bool, str]:
    if not FINALEME_MODELS_DIR.is_dir():
        return False, f"{FINALEME_MODELS_DIR} does not exist."
    missing = [f for f in FINALEME_PRETRAINED_FILES
               if not (FINALEME_MODELS_DIR / f).is_file()]
    if missing:
        return False, (
            f"missing {missing} in {FINALEME_MODELS_DIR}. "
            f"Zenodo record 14013719 is the source. See "
            f"docs/CROSS_PLATFORM_FINALEME.md."
        )
    return True, str(FINALEME_MODELS_DIR)


def _check_reference_files() -> Tuple[bool, str]:
    missing = [(name, path) for name, path in FINALEME_REF_FILES if not path.is_file()]
    if missing:
        names = ", ".join(name for name, _ in missing)
        return False, f"missing {names} (1.4 GB total)."
    return True, "all 5 reference files present"


def _check_finaleme_output(finaleme_dir: Path) -> Tuple[bool, str]:
    if not finaleme_dir.is_dir():
        return False, f"{finaleme_dir} does not exist (or is not a directory)."
    # Look for any BetaValues.tsv, *.decoded.bed.gz, *.cpg_features.hg19.bed.gz
    found = list(finaleme_dir.glob("*BetaValues.tsv")) + \
            list(finaleme_dir.glob("*decoded.bed.gz")) + \
            list(finaleme_dir.glob("*cpg_features.hg19.bed.gz"))
    if not found:
        return False, (
            f"{finaleme_dir} exists but contains no FinaleMe output "
            f"(BetaValues.tsv / *.decoded.bed.gz / *.cpg_features.hg19.bed.gz)."
        )
    return True, f"{len(found)} output file(s) in {finaleme_dir}"


def _check_finaledb_features(features_dir: Path) -> Tuple[bool, str]:
    if not features_dir.is_dir():
        return False, (
            f"{features_dir} does not exist. The FinaleDB features cache "
            f"is at /Users/hermes/cfdna-fragmentomics-pipeline/data/features/."
        )
    # Count .npy files; there should be 5 per sample
    npy = list(features_dir.glob("*.npy"))
    if not npy:
        return False, f"{features_dir} exists but contains no .npy files."
    sample_ids = {p.name.rsplit(".", 2)[0] for p in npy}
    return True, (
        f"{len(sample_ids)} samples × {len(npy) // max(len(sample_ids), 1)} channels "
        f"({len(npy)} total .npy files) in {features_dir}"
    )


def _check_java() -> Tuple[bool, str]:
    """Check for any Java 21+ runtime available on PATH."""
    if shutil.which("java") is None:
        return False, "no `java` binary on $PATH"
    try:
        out = subprocess.run(
            ["java", "-version"], capture_output=True, text=True, timeout=10
        )
        if out.returncode == 0:
            version = (
                (out.stderr or out.stdout).splitlines()[0]
                if out.stderr or out.stdout
                else "unknown"
            )
            return True, version
        return False, f"`java -version` exited {out.returncode}"
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return False, "java -version did not respond"


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(
        description=(
            "Probe cross-platform pipeline readiness. "
            "ALWAYS exits 0; this is a probe, not a gate."
        ),
    )
    ap.add_argument(
        "--finaleme-output-dir", type=Path, default=DEFAULT_FINALEME_OUTPUT_DIR,
        help="Path to the FinaleMe output directory (BetaValues.tsv / .bed.gz).",
    )
    ap.add_argument(
        "--features-dir", type=Path, default=DEFAULT_FINALEDB_FEATURES_DIR,
        help="Path to the FinaleDB features cache (.npy files).",
    )
    args = ap.parse_args()

    now = datetime.now(timezone.utc).isoformat()

    print(_bold("FinaleMe cross-platform readiness probe"))
    print(f"generated_at: {now}")
    print(f"finaleme_home: {FINALEME_HOME}")
    print(f"finaleme_output_dir: {args.finaleme_output_dir}")
    print(f"features_dir: {args.features_dir}")
    print()

    # Checklist rows: (label, found, detail)
    rows: List[Tuple[str, bool, str]] = []
    rows.append(("Java runtime",       *_check_java()))
    rows.append(("FinaleMe JAR",        *_check_jar()))
    rows.append(("Pretrained HMM models", *_check_pretrained()))
    rows.append(("Reference files",     *_check_reference_files()))
    rows.append(("FinaleMe output",     *_check_finaleme_output(args.finaleme_output_dir)))
    rows.append(("FinaleDB features",   *_check_finaledb_features(args.features_dir)))

    print(_bold("Checklist"))
    for label, ok, detail in rows:
        marker = _green("[OK]") if ok else _red("[MISSING]")
        print(f"  {marker:<10}  {label:<22}  {detail}")
    print()

    # Verdict: data_source
    finaleme_output_ok, meth_detail = _check_finaleme_output(args.finaleme_output_dir)
    finaledb_ok, frag_detail = _check_finaledb_features(args.features_dir)
    if finaleme_output_ok and finaledb_ok:
        data_source = "both"
        verdict = _green(
            "data_source = both. Run scripts/cross_platform_finaledb_validate.py "
            "for a real cross-platform AUC."
        )
    elif finaleme_output_ok:
        data_source = "methylation_only"
        verdict = _yellow(
            "data_source = methylation_only. Cross-platform AUC will NOT be printed "
            "(fragmentomics channel missing — install FinaleDB features or move "
            "to a machine that has them)."
        )
    elif finaledb_ok:
        data_source = "fragmentomics_only"
        verdict = _yellow(
            "data_source = fragmentomics_only. Cross-platform AUC will NOT be printed "
            "(methylation channel missing — install FinaleMe JAR + pretrained models "
            "+ reference files; run scripts/finaleme_pipeline.py decode to produce "
            "per-CpG β-values)."
        )
    else:
        data_source = "no_data"
        verdict = _red(
            "data_source = no_data. Both channels missing — no metrics possible."
        )

    print(_bold("Verdict"))
    print(f"  data_source: {data_source}")
    print(f"  {verdict}")
    print()

    # Helpful next-step links
    if data_source != "both":
        print(_bold("Next steps"))
        if not _check_pretrained()[0]:
            print(
                "  • Download pretrained HMM models from Zenodo record 14013719 "
                f"into {FINALEME_MODELS_DIR}/"
            )
        if not _check_jar()[0]:
            print(
                f"  • Install FinaleMe JAR at {FINALEME_HOME}/ "
                "(see docs/CROSS_PLATFORM_FINALEME.md §Prerequisites)."
            )
        if not _check_reference_files()[0]:
            print(
                "  • Download reference files (~1.4 GB) — "
                "see cfdna-fragmentomics skill §Reference files."
            )

    # Always exit 0 — probe never fails
    return 0


if __name__ == "__main__":
    sys.exit(main())
