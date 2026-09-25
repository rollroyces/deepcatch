#!/usr/bin/env python3
"""FinaleMe pipeline orchestrator (Java 21 + Maven + pretrained models).

FinaleMe (Liu et al. Nat Commun 15:2790, 2024) is a Java tool that
imputes per-CpG methylation β-values from cfDNA WGS fragments. It runs
in three steps:

    Step 1: extract CpG-fragment features from a BAM/tabix input
            (``CpgFeatureMatrixBuilder``, ~14 min per sample on M4)
    Step 2: train per-CpG HMM models (skipped when pretrained models
            are available from Zenodo record 14013719)
    Step 3: decode per-CpG β-values using the trained (or pretrained)
            HMM model (``FinaleMe``, ~32 s per sample for chr22)

This script orchestrates the three steps. It deliberately
**never fabricates predictions**: when the FinaleMe JAR or pretrained
models or reference files are missing, it emits a JSON
provenance report and exits non-zero with a clear message naming
which input is missing. The downstream bridge
(``scripts/finaleme_to_deepcatch_bridge.py``) consumes the final
``*.decoded.bed.gz`` outputs.

Subcommands:
    probe        Probe JAR / pretrained models / reference files (exit 0).
    features     Run Step 1 (CpgFeatureMatrixBuilder) on a BAM/tabix input.
    decode       Run Step 3 (FinaleMe decode with pretrained models); skip
                 Step 2 training when models are present.

The cross-platform validator (``scripts/cross_platform_finaledb_validate.py``)
calls ``probe`` to decide whether the methylation channel is runnable.

Prerequisites (full recipe in ``docs/CROSS_PLATFORM_FINALEME.md``):
    * Java 21 (`brew install --cask zulu@21`; or the tarball install at
      ``~/.hermes/.local/jdk/`` per the deepcatch skill)
    * FinaleMe JAR (built from source with `mvn clean package`, or
      the hybrid JAR per the cfdna-fragmentomics skill — Zenodo
      pretrained models deserialize against a hybrid JAR that combines
      v0.61's streaming decoder with v0.58.1's TreeMap-based HMM
      classes, with two legacy-package remaps).
    * Reference files (1.4 GB total: hg19.2bit, methylation prior bw,
      CpG motif, mappability BED, chrom sizes).
    * Pretrained HMM models from Zenodo record 14013719
      (~150 kB, two files).

Usage:
    env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_pipeline.py probe
    env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_pipeline.py \\
        features --input sample.bam --output-dir /tmp/finaleme/
    env -u PYTHONPATH ./.venv/bin/python scripts/finaleme_pipeline.py \\
        decode --features-file features.bed.gz \\
        --output-dir /tmp/finaleme/ \\
        --model both
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ─────────────────────────────────────────────────────────────────────
# Standard paths
# ─────────────────────────────────────────────────────────────────────

FINALEME_HOME = Path("~/.hermes/.local/finaleme").expanduser()
FINALEME_JAR_CANDIDATES: List[Path] = [
    FINALEME_HOME / "FinaleMe-0.61-jar-with-dependencies.jar",
    FINALEME_HOME / "FinaleMe-hybrid-jar-with-dependencies.jar",
    FINALEME_HOME / "FinaleMe-0.58.1-jar-with-dependencies.jar",
]
FINALEME_MODELS_DIR = FINALEME_HOME / "models"
HEALTHY_MODEL_FILE = "healthy_WGS.mincg7.example.hmm_model"
CANCER_MODEL_FILE = "cancer_WGS.mincg7.example.hmm_model"

# JVM heap ceiling for FinaleMe on M4 / 16 GB
DEFAULT_JVM_HEAP = "12G"
DEFAULT_PARALLEL_GC_THREADS = 10

# Reference-file locations
REFERENCE_FILES: Dict[str, Path] = {
    "hg19.2bit":                  FINALEME_HOME / "hg19.2bit",
    "methy_prior_bw":             FINALEME_HOME / "wgbs_buffyCoat_jensen2015GB.methy.hg19.bw",
    "cpg_motif_bedgraph":         FINALEME_HOME / "CG_motif.hg19.common_chr.pos_only.bedgraph.gz",
    "mappability_bed":            FINALEME_HOME / "mappability.bed",
    "chrom_sizes":                FINALEME_HOME / "hg19.chrom.sizes",
}

# FinaleMe Java classes
FINALEME_STEP1_CLASS = "edu.northwestern.epifluidlab.finaleme.utils.CpgFeatureMatrixBuilder"
FINALEME_STEP3_CLASS = "edu.northwestern.epifluidlab.finaleme.hmm.FinaleMe"

REQUIRED_INPUTS: List[str] = [
    "java_runtime",            # `java` on PATH
    "finaleme_jar",            # any of the JAR candidates present
    "pretrained_models",       # both healthy + cancer HMM models
    "reference_files",         # all 5 reference files
]


# ─────────────────────────────────────────────────────────────────────
# Errors
# ─────────────────────────────────────────────────────────────────────


class MissingInputError(RuntimeError):
    """Raised when a required FinaleMe input is missing.

    The orchestrator catches this and prints a clear ``MISSING``
    message naming the absent input, then exits 1.
    """


class JavaNotFoundError(MissingInputError):
    pass


class JarNotFoundError(MissingInputError):
    pass


class PretrainedModelsMissing(MissingInputError):
    pass


class ReferenceFilesMissing(MissingInputError):
    pass


# ─────────────────────────────────────────────────────────────────────
# Probe
# ─────────────────────────────────────────────────────────────────────


@dataclass
class FinaleMeReadiness:
    """Per-input readiness verdict (matches the probe JSON schema)."""

    java_runtime: Dict = field(default_factory=dict)
    finaleme_jar: Dict = field(default_factory=dict)
    pretrained_models: Dict = field(default_factory=dict)
    reference_files: Dict = field(default_factory=dict)
    generated_at: str = ""

    def is_runnable(self) -> bool:
        return all([
            self.java_runtime.get("present"),
            self.finaleme_jar.get("present"),
            self.pretrained_models.get("present"),
            self.reference_files.get("present"),
        ])

    def to_dict(self) -> Dict:
        d = asdict(self)
        d["runnable"] = self.is_runnable()
        return d


def probe_java() -> Dict:
    """Detect `java` on $PATH. Returns `{present, version, path}`."""
    path = shutil.which("java")
    if path is None:
        return {"present": False, "error": "no `java` on $PATH"}
    try:
        out = subprocess.run(
            ["java", "-version"], capture_output=True, text=True, timeout=10
        )
        if out.returncode != 0:
            return {"present": False, "error": f"`java -version` exit {out.returncode}"}
        version_line = (out.stderr or out.stdout).splitlines()[0] if out.stderr or out.stdout else "?"
        return {"present": True, "version": version_line, "path": path}
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        return {"present": False, "error": f"java -version: {e}"}


def probe_jar() -> Dict:
    """Find a FinaleMe JAR among the candidate paths."""
    for cand in FINALEME_JAR_CANDIDATES:
        if cand.is_file():
            return {
                "present": True,
                "path": str(cand),
                "size_mb": round(cand.stat().st_size / 1024 / 1024, 1),
            }
    return {
        "present": False,
        "tried": [str(p) for p in FINALEME_JAR_CANDIDATES],
        "error": (
            "no FinaleMe JAR found. See docs/CROSS_PLATFORM_FINALEME.md "
            "§Prerequisites for the v0.61 / hybrid JAR recipe."
        ),
    }


def probe_pretrained_models() -> Dict:
    """Both pretrained HMM files must be present."""
    if not FINALEME_MODELS_DIR.is_dir():
        return {
            "present": False,
            "path": str(FINALEME_MODELS_DIR),
            "error": f"models directory does not exist: {FINALEME_MODELS_DIR}",
        }
    h = FINALEME_MODELS_DIR / HEALTHY_MODEL_FILE
    c = FINALEME_MODELS_DIR / CANCER_MODEL_FILE
    if h.is_file() and c.is_file():
        return {
            "present": True,
            "healthy_model": str(h),
            "cancer_model": str(c),
        }
    missing = []
    if not h.is_file():
        missing.append(HEALTHY_MODEL_FILE)
    if not c.is_file():
        missing.append(CANCER_MODEL_FILE)
    return {
        "present": False,
        "path": str(FINALEME_MODELS_DIR),
        "missing": missing,
        "error": (
            f"missing pretrained models: {missing}. "
            f"Source: Zenodo record 14013719."
        ),
    }


def probe_reference_files() -> Dict:
    """All five reference files must be present."""
    missing = [(name, str(path)) for name, path in REFERENCE_FILES.items()
               if not path.is_file()]
    if missing:
        return {
            "present": False,
            "missing": missing,
            "error": (
                f"missing {len(missing)} reference file(s) "
                f"(~1.4 GB total). See deepcatch skill §FinaleMe."
            ),
        }
    return {
        "present": True,
        "files": {n: str(p) for n, p in REFERENCE_FILES.items()},
    }


def probe_all() -> FinaleMeReadiness:
    r = FinaleMeReadiness()
    r.java_runtime = probe_java()
    r.finaleme_jar = probe_jar()
    r.pretrained_models = probe_pretrained_models()
    r.reference_files = probe_reference_files()
    r.generated_at = datetime.now(timezone.utc).isoformat()
    return r


# ─────────────────────────────────────────────────────────────────────
# Subcommand handlers
# ─────────────────────────────────────────────────────────────────────


def cmd_probe(args: argparse.Namespace) -> int:
    r = probe_all()
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r.to_dict(), indent=2))
    print(json.dumps(r.to_dict(), indent=2))
    print(f"probe_report: {out}")
    return 0


def cmd_features(args: argparse.Namespace) -> int:
    """Run FinaleMe Step 1 (CpgFeatureMatrixBuilder) — produces
    <prefix>.cpg_features.hg19.bed.gz.

    Uses tabix mode (-fragmentInputTabix -fragStrandColumn N) by
    default — BAM mode hangs on M4 / 16 GB.
    """
    r = probe_all()
    if not r.java_runtime.get("present"):
        raise JavaNotFoundError(r.java_runtime.get("error", "java missing"))
    if not r.finaleme_jar.get("present"):
        raise JarNotFoundError(r.finaleme_jar.get("error", "FinaleMe JAR missing"))
    if not r.reference_files.get("present"):
        raise ReferenceFilesMissing(r.reference_files.get("error"))

    input_path = Path(args.input)
    if not input_path.is_file():
        raise MissingInputError(f"input not found: {input_path}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (input_path.stem.replace(".frag.bed.gz", "") + ".cpg_features.hg19.bed.gz")

    jar_path = r.finaleme_jar["path"]
    cmd = [
        "java",
        f"-Xmx{args.jvm_heap}",
        "-XX:+UseParallelGC",
        f"-XX:ParallelGCThreads={args.parallel_gc_threads}",
        "-cp", jar_path,
        FINALEME_STEP1_CLASS,
        str(REFERENCE_FILES["hg19.2bit"]),
        str(REFERENCE_FILES["cpg_motif_bedgraph"]),
        str(REFERENCE_FILES["cpg_motif_bedgraph"]),
        str(input_path),
        str(output_path),
        "-fragmentInputTabix",
        "-fragStrandColumn", str(args.frag_strand_column),
        "-valueWigs", f"methyPrior:0:{REFERENCE_FILES['methy_prior_bw']}",
        "-inferMethyFromValueWig",
        "-useNoChrPrefixBam",
        "-t", str(args.parallel_gc_threads),
    ]
    print(" ".join(cmd))
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        return rc
    print(f"features_output: {output_path}")
    return 0


def cmd_decode(args: argparse.Namespace) -> int:
    """Run FinaleMe Step 3 (decode per-CpG β using pretrained models).

    Always uses -decodeModeOnly to skip Step 2 training. Decodes with
    BOTH healthy + cancer HMM models so the LR baseline has a
    two-channel methylation feature (cancer mean β − healthy mean β)
    plus the direct cancer model probability.
    """
    r = probe_all()
    if not r.java_runtime.get("present"):
        raise JavaNotFoundError(r.java_runtime.get("error", "java missing"))
    if not r.finaleme_jar.get("present"):
        raise JarNotFoundError(r.finaleme_jar.get("error", "FinaleMe JAR missing"))
    if not r.pretrained_models.get("present"):
        raise PretrainedModelsMissing(r.pretrained_models.get("error"))
    if not r.reference_files.get("present"):
        raise ReferenceFilesMissing(r.reference_files.get("error"))

    features_path = Path(args.features_file)
    if not features_path.is_file():
        raise MissingInputError(f"features file not found: {features_path}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    jar_path = r.finaleme_jar["path"]
    models_to_run: List[Tuple[str, Path, Path]] = []
    if args.model in ("healthy", "both"):
        out = output_dir / f"{features_path.stem}_healthy_decoded.bed.gz"
        models_to_run.append(("healthy", Path(r.pretrained_models["healthy_model"]), out))
    if args.model in ("cancer", "both"):
        out = output_dir / f"{features_path.stem}_cancer_decoded.bed.gz"
        models_to_run.append(("cancer", Path(r.pretrained_models["cancer_model"]), out))

    n_ok = 0
    for label, model_path, out_path in models_to_run:
        cmd = [
            "java",
            f"-Xmx{args.jvm_heap}",
            "-cp", jar_path,
            FINALEME_STEP3_CLASS,
            str(model_path),
            str(features_path),
            str(out_path),
            "-decodeModeOnly",
            "-bwOutput",
            "-chromSizeFile", str(REFERENCE_FILES["chrom_sizes"]),
        ]
        print(f"[{label}] " + " ".join(cmd))
        rc = subprocess.run(cmd).returncode
        if rc != 0:
            print(f"[{label}] decode failed rc={rc}", file=sys.stderr)
            return rc
        n_ok += 1

    print(f"decoded_with_{n_ok}_model(s); outputs in {output_dir}")
    return 0


# ─────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "FinaleMe pipeline orchestrator. Probes for required inputs "
            "(FinaleMe JAR, pretrained models, reference files); never "
            "synthesizes predictions when missing."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=False)

    p_probe = sub.add_parser("probe", help="Probe inputs (always exits 0).")
    p_probe.add_argument(
        "--output", type=str, default="results/finaleme_readiness.json",
        help="Path to write the probe JSON.",
    )
    p_probe.set_defaults(func=cmd_probe)

    p_features = sub.add_parser("features", help="Run FinaleMe Step 1.")
    p_features.add_argument("--input", type=str, required=True,
                            help="Input BAM/tabix fragment file.")
    p_features.add_argument("--output-dir", type=str, required=True,
                            help="Directory for output .cpg_features.hg19.bed.gz.")
    p_features.add_argument("--jvm-heap", type=str, default=DEFAULT_JVM_HEAP,
                            help=f"JVM -Xmx (default {DEFAULT_JVM_HEAP}).")
    p_features.add_argument("--parallel-gc-threads", type=int,
                            default=DEFAULT_PARALLEL_GC_THREADS,
                            help=f"Parallel-GC threads (default {DEFAULT_PARALLEL_GC_THREADS}).")
    p_features.add_argument("--frag-strand-column", type=int, default=6,
                            help="Column index for fragment strand (1-based).")
    p_features.set_defaults(func=cmd_features)

    p_decode = sub.add_parser("decode", help="Run FinaleMe Step 3 (decode).")
    p_decode.add_argument("--features-file", type=str, required=True,
                          help="Per-sample features file (.cpg_features.hg19.bed.gz).")
    p_decode.add_argument("--output-dir", type=str, required=True,
                          help="Directory for output *.decoded.bed.gz.")
    p_decode.add_argument("--model", type=str, default="both",
                          choices=("healthy", "cancer", "both"),
                          help="Which pretrained model(s) to decode with.")
    p_decode.add_argument("--jvm-heap", type=str, default=DEFAULT_JVM_HEAP,
                          help=f"JVM -Xmx (default {DEFAULT_JVM_HEAP}).")
    p_decode.set_defaults(func=cmd_decode)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if not getattr(args, "command", None):
        parser.print_help()
        return 0
    try:
        return int(args.func(args))
    except MissingInputError as e:
        print(f"MISSING: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
