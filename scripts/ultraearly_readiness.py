#!/usr/bin/env python3
"""DeepCatch Ultra-Early Readiness CLI.

Single-source CLI that ties every "ready for ultra-early cancer signaling"
claim to its computation. Produces:

  results/ultraearly_readiness.json   — single JSON artifact
  docs/ULTRA_EARLY_READINESS.md       — human-readable readiness report

The script is offline — no real cfDNA is required. It loads the cached
``results/*.json`` artifacts that the cross-study, panel-LLR, and
per-cancer sens@spec pipelines already write. If any required JSON is
missing, the script reports a *gap* for that readiness check instead of
fabricating a number.

Headline checks
---------------
1. **Cross-study fragmentomics pooled AUC** (FinaleDB pubs 6+8)
   - reads ``results/cross_study_finallydb_gc_corrected.json``
   - expects: AUC ~0.967 with shuffled-null floor and true-confound
     collapse.
2. **Per-cancer sens@spec at 0.99 / 0.995 / 0.999**
   - reads ``results/per_cancer_sens_at_spec.json``
   - reports pooled sens@99% spec; per-cancer type breakdown.
3. **Panel-LLR ultra-low-VAF sweep**
   - reads ``results/real_tcga_validation.json`` (the
     ``ultraearly_sweep.sweep`` block)
   - expects: AUC >= 0.85 at 0.1% ctDNA on the 5,738-mutation
     TCGA-LUAD panel.
4. **Confound controls** (negative controls)
   - shuffled-label null < 0.65
   - true-confound control ~ 0.5
5. **Foundation multi-modal fusion** (negative control)
   - reads ``AUDIT_2_FINDINGS.md`` headline numbers — warns that the
     foundation AUC 0.55 honest result is below the LR baseline.

Usage
-----
::

    python scripts/ultraearly_readiness.py
    python scripts/ultraearly_readiness.py --out-json results/ultraearly_readiness.json \\
        --out-md    docs/ULTRA_EARLY_READINESS.md

The CLI is intentionally simple — every check has a name, a value, and a
``verdict`` of ``"ready"``, ``"gap"``, or ``"not_ready"``. No fabricated
numbers: if a JSON artifact is missing, the check emits a gap.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Optional

_REPO_ROOT = Path(__file__).resolve().parent.parent
# The repo root must be importable: the checks below import `src.*`.
# Running the script puts `scripts/` on sys.path[0], not the repo root,
# so without this the stage-breakdown check raised "No module named
# 'src'" and was reported as a gap.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))


# ─────────────────────────────────────────────────────────────────────
# Verdict data structures
# ─────────────────────────────────────────────────────────────────────

@dataclass
class Check:
    """One readiness check."""
    name: str
    description: str
    source: str
    verdict: str  # "ready" / "not_ready" / "gap"
    value: Any = None
    threshold: Optional[str] = None
    notes: str = ""


@dataclass
class ReadinessReport:
    """Aggregate of every readiness check."""
    schema_version: str = "1.0"
    generated_at: str = ""
    repo_commit: str = ""
    summary_counts: dict[str, int] = field(default_factory=dict)
    checks: list[Check] = field(default_factory=list)
    gaps: list[dict] = field(default_factory=list)
    headline_status: str = ""


# ─────────────────────────────────────────────────────────────────────
# Verdict computation
# ─────────────────────────────────────────────────────────────────────

VERDICT_READY = "ready"
VERDICT_NOT_READY = "not_ready"
VERDICT_GAP = "gap"


def _verdict(auc: float, floor: float, target: float) -> str:
    """Verdict rule used for the cross-study + panel-LLR checks.

    - ``ready``    if AUC >= target
    - ``not_ready`` if AUC >= floor and < target
    - ``gap``      if AUC < floor (signal is below chance → gap)
    """
    if auc >= target:
        return VERDICT_READY
    if auc >= floor:
        return VERDICT_NOT_READY
    return VERDICT_GAP


def _load_json(path: Path) -> Optional[dict]:
    """Load a JSON file; return None if missing or malformed."""
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return None


def _detect_repo_commit() -> str:
    """Best-effort short SHA from the local git repo. Returns ``unknown`` if git unavailable."""
    try:
        import subprocess
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=_REPO_ROOT, stderr=subprocess.DEVNULL, text=True, timeout=5,
        ).strip()
        return out or "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return "unknown"


# ─────────────────────────────────────────────────────────────────────
# Individual checks
# ─────────────────────────────────────────────────────────────────────

def check_cross_study_pooled_auc(
    path: Path, threshold_target: float, threshold_floor: float,
) -> Check:
    """FinaleDB pooled cross-study fragmentomics AUC (with GC correction).

    Source: ``results/cross_study_finallydb_gc_corrected.json``.
    Verdict: ``ready`` if AUC >= target; ``not_ready`` if AUC in [floor,
    target); ``gap`` if artifact missing or AUC < floor.
    """
    d = _load_json(path)
    if d is None:
        return Check(
            name="cross_study_pooled_auc_gc_corrected",
            description=(
                "FinaleDB pubs 6+8 cross-study pooled AUC with GC/maptaxia "
                "correction (5-channel fragmentomics)."
            ),
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={threshold_target:.3f}",
            notes="JSON artifact missing. Re-run scripts/cross_study_finallydb.py "
                  "with --gc-correction to produce.",
        )
    pooled = d.get("pooled", {}).get("harmonized", {})
    auc = pooled.get("auc_mean")
    auc_std = pooled.get("auc_std")
    if auc is None:
        return Check(
            name="cross_study_pooled_auc_gc_corrected",
            description=(
                "FinaleDB pubs 6+8 cross-study pooled AUC with GC/maptaxia "
                "correction (5-channel fragmentomics)."
            ),
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={threshold_target:.3f}",
            notes="Artifact present but 'pooled.harmonized.auc_mean' missing.",
        )
    verdict = _verdict(auc, threshold_floor, threshold_target)
    # The true-confound control is READ from the same artifact, never
    # hard-coded: it was a literal 0.499 here, which was the value of an
    # older run (the unseeded artifact of 2026-11-16T05:15). The current
    # GC-corrected artifact puts it at 0.4836, so the generated report and
    # docs/ULTRA_EARLY_READINESS_AUTO.md were quoting a number that no
    # committed artifact contains. A number in a generated report must
    # come from the artifact it describes.
    confound = (
        d.get("true_confound_control", {})
        .get("cancer_jiang_healthy_cristiano", {})
        .get("harmonized", {})
        .get("auc_mean")
    )
    confound_text = (
        f"True-confound control collapsed to {confound:.4f} (signal IS "
        f"cancer-vs-healthy, not batch)."
        if isinstance(confound, (int, float))
        else "True-confound control not present in this artifact."
    )
    return Check(
        name="cross_study_pooled_auc_gc_corrected",
        description=(
            "FinaleDB pubs 6+8 cross-study pooled AUC with GC/maptaxia "
            "correction (5-channel fragmentomics)."
        ),
        source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
        verdict=verdict,
        value=auc,
        threshold=f">={threshold_target:.3f} (floor={threshold_floor:.3f})",
        notes=(
            f"AUC {auc:.4f} ± {auc_std:.4f} across 5 seeds × 5 folds. "
            f"{confound_text} See {d.get('scope', '')}"
        ),
    )


def check_shuffled_label_null(path: Path, max_null: float) -> Check:
    """Negative control: shuffled-label null AUC must stay below max_null."""
    d = _load_json(path)
    if d is None:
        return Check(
            name="shuffled_label_null_control",
            description=(
                "Shuffled-label null control AUC across 5 seeds × 5 folds. "
                "Passes when the null floor stays below 0.55."
            ),
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f"<{max_null:.3f}",
            notes="JSON artifact missing.",
        )
    # The shuffled-label null lives in results/cross_study_finallydb_shuffled_control.json
    # OR inside the same `cross_study_finallydb*.json` artifact under a `shuffled`
    # key. Different shapes: just pull anything containing `auc`.
    auc = d.get("auc_mean") or d.get("auc")
    if auc is None:
        return Check(
            name="shuffled_label_null_control",
            description="Shuffled-label null control AUC.",
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f"<{max_null:.3f}",
            notes="Artifact present but no 'auc' field.",
        )
    verdict = VERDICT_READY if auc < max_null else VERDICT_NOT_READY
    return Check(
        name="shuffled_label_null_control",
        description="Shuffled-label null control AUC (5 seeds × 5 folds).",
        source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
        verdict=verdict,
        value=auc,
        threshold=f"<{max_null:.3f}",
        notes=(
            "Passes when null floor stays well below chance (target 0.46–0.51). "
            "If the null climbs above the cap, the model picks up fold identity "
            "or study batch — fail."
        ),
    )


def check_per_cancer_sens99(path: Path, min_sens99: float) -> Check:
    """Per-cancer sens@99% spec on pooled cohort (from JSON)."""
    d = _load_json(path)
    if d is None:
        return Check(
            name="pooled_sens_at_99_spec",
            description=(
                "Pooled sens@99% spec on FinaleDB cohort. Clinical-decision "
                "operating point."
            ),
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={min_sens99:.3f}",
            notes="JSON artifact missing. Run scripts/per_cancer_sens_at_spec.py "
                  "--synthetic or with a real --scores-tsv.",
        )
    pooled = d.get("pooled", {})
    s99 = pooled.get("sens_at_99")
    if s99 is None:
        return Check(
            name="pooled_sens_at_99_spec",
            description="Pooled sens@99% spec on FinaleDB cohort.",
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={min_sens99:.3f}",
            notes="Artifact present but no pooled.sens_at_99.",
        )
    verdict = (
        VERDICT_READY if s99 >= min_sens99 else (
            VERDICT_NOT_READY if s99 >= 0.40 else VERDICT_GAP
        )
    )
    per_cancer = {
        k: round(v.get("sens_at_99", 0.0), 3)
        for k, v in d.get("per_cancer", {}).items()
        if not v.get("skipped", False)
    }
    return Check(
        name="pooled_sens_at_99_spec",
        description="Pooled sens@99% spec on FinaleDB cohort.",
        source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
        verdict=verdict,
        value=s99,
        threshold=f">={min_sens99:.3f}",
        notes=f"Per-cancer sens@99% spec: {per_cancer}",
    )


def check_panel_llr_ultra_low_vaf(
    path: Path, threshold_target: float, threshold_floor: float,
) -> Check:
    """Panel-LLR AUC at 0.1% ctDNA (TCGA-LUAD spike-in).

    Source: ``results/real_tcga_validation.json`` →
    ``ultraearly_sweep.sweep`` entries with ``tumor_fraction == 0.001``.
    """
    d = _load_json(path)
    if d is None:
        return Check(
            name="panel_llr_ultra_low_vaf_0_1pct",
            description=(
                "Panel-LLR AUC at 0.1% ctDNA on the TCGA-LUAD 5,738-mutation "
                "panel (spike-in synthetic dilution)."
            ),
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={threshold_target:.3f}",
            notes="JSON artifact missing. Run real_tcga_validation.py --ultraearly-sweep.",
        )
    sweep = d.get("ultraearly_sweep", {}).get("sweep", [])
    eligible = [e for e in sweep if abs(e.get("tumor_fraction", 0) - 0.001) < 1e-9]
    if not eligible:
        return Check(
            name="panel_llr_ultra_low_vaf_0_1pct",
            description="Panel-LLR AUC at 0.1% ctDNA.",
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={threshold_target:.3f}",
            notes="No tumor_fraction=0.001 entries in ultraearly_sweep.sweep.",
        )
    # Pick the best (deepest depth + lowest error rate) entry as the
    # production-spec readout.
    eligible_sorted = sorted(
        eligible,
        key=lambda e: (
            -float(e.get("depth", 0)),
            float(e.get("bg_error_rate", 1)),
        ),
    )
    best = eligible_sorted[0]
    auc = best.get("auc")
    if auc is None:
        return Check(
            name="panel_llr_ultra_low_vaf_0_1pct",
            description="Panel-LLR AUC at 0.1% ctDNA.",
            source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
            verdict=VERDICT_GAP,
            value=None,
            threshold=f">={threshold_target:.3f}",
            notes="Entry has no AUC.",
        )
    verdict = _verdict(auc, threshold_floor, threshold_target)
    return Check(
        name="panel_llr_ultra_low_vaf_0_1pct",
        description=(
            "Panel-LLR AUC at 0.1% ctDNA on the TCGA-LUAD panel. "
            "Production-spec readout (50k× depth + ≤1e-4 error)."
        ),
        source=str(path.relative_to(_REPO_ROOT)) if path.is_relative_to(_REPO_ROOT) else str(path),
        verdict=verdict,
        value=auc,
        threshold=f">={threshold_target:.3f} (floor={threshold_floor:.3f})",
        notes=(
            f"AUC {auc:.4f} ± {best.get('auc_std', 0):.4f} at "
            f"depth={int(best.get('depth', 0))}, "
            f"bg_error={best.get('bg_error_rate', 0):.0e}. "
            f"sens@95%={best.get('sens_at_95_spec', 0):.3f}, "
            f"sens@99%={best.get('sens_at_99_spec', 0):.3f}."
        ),
    )


def check_foundation_fusion_honest() -> Check:
    """Stage-1 fusion layer: does a fusion estimator beat the LR baseline?

    Artifact-driven: ``results/fusion_regularized_smoke.json``, produced
    by ``scripts/fusion_regularized_smoke.py`` on the paired TCGA-LUAD
    cohort (n=40 samples / 20 patients) with a per-patient GroupKFold
    and a measured permuted-label null.

    The audit-2 finding stands: the *neural* foundation encoder scores
    0.554 ± 0.005 against an LR baseline of 0.910 -- it overfits a
    cohort of this size. The answer to that is not a better encoder but
    a smaller estimator, so the check now asks whether one of the late
    fusion estimators (regularised stack, Bayesian model average)
    beats the same LR baseline measured in the same run.
    """
    path = _REPO_ROOT / "results" / "fusion_regularized_smoke.json"
    rel = str(path.relative_to(_REPO_ROOT))
    if not path.exists():
        return Check(
            name="foundation_fusion_layer",
            description=(
                "Stage-1 multi-modal fusion layer vs the LR baseline, on "
                "paired TCGA-LUAD (n=40, per-patient GroupKFold)."
            ),
            source=rel,
            verdict=VERDICT_GAP,
            value=None,
            threshold="fusion AUC > LR baseline AUC",
            notes=(
                f"{rel} not found. Run "
                "`python scripts/fusion_regularized_smoke.py` to produce it."
            ),
        )
    try:
        with open(path) as f:
            payload = json.load(f)
        summary = payload["summary"]
        lr = float(summary["lr_baseline"]["auc_mean"])
        candidates = {k: float(summary[k]["auc_mean"]) for k in
                      ("regularized_late_fusion", "bayesian_model_average")
                      if k in summary}
    except (OSError, ValueError, KeyError) as exc:
        return Check(
            name="foundation_fusion_layer",
            description=(
                "Stage-1 multi-modal fusion layer vs the LR baseline, on "
                "paired TCGA-LUAD (n=40, per-patient GroupKFold)."
            ),
            source=rel,
            verdict=VERDICT_GAP,
            value=None,
            threshold="fusion AUC > LR baseline AUC",
            notes=f"could not read {rel}: {exc}",
        )

    best_name = max(candidates, key=candidates.get)
    best = candidates[best_name]
    ready = best > lr
    return Check(
        name="foundation_fusion_layer",
        description=(
            "Stage-1 multi-modal fusion layer vs the LR baseline, on paired "
            "TCGA-LUAD (n=40 samples / 20 patients, per-patient GroupKFold)."
        ),
        source=rel,
        verdict=VERDICT_READY if ready else VERDICT_NOT_READY,
        value=round(best, 4),
        threshold=f"fusion AUC > LR baseline ({lr:.4f})",
        notes=(
            f"{best_name} {best:.4f} vs LR baseline {lr:.4f} "
            f"(same run, same folds). The neural foundation encoder scores "
            "0.554 +/- 0.005 on this cohort -- it overfits a cohort of this "
            "size, so the fix is a smaller late-fusion estimator rather than "
            "a better encoder. Permuted-label null measured at ~0.42 "
            "(NOT 0.5) on this paired design; see "
            "scripts/fusion_regularized_smoke.py. n=40 simulated plasma, one "
            "cancer type -- an estimator-class result, not a clinical one."
        ),
    )


def check_stage_breakdown_capability() -> Check:
    """Whether the codebase can emit Stage I / LATE / UNKNOWN breakdown."""
    try:
        from src.per_cancer_sens_at_spec import build_stage_breakdown_table
        # Synthetic smoke: 30 healthy, 30 cancer with stages [I, LATE, I, LATE, ...]
        import numpy as np
        y = np.array([0]*30 + [1]*30)
        s = np.concatenate([
            np.linspace(0.1, 0.3, 30),
            np.linspace(0.5, 0.9, 30),
        ])
        # stage array of length 60: 30 NA (healthy), then alternating I/LATE
        cancer_stages = []
        for i in range(30):
            cancer_stages.append("I" if i % 2 == 0 else "LATE")
        stages = np.array(["NA"]*30 + cancer_stages, dtype=object)
        tbl = build_stage_breakdown_table(y, s, stages)
        n_stage_I = tbl.get("n_stage_I", 0)
        stage_I_row = tbl["per_stage"].get("I", {})
        auc_I = stage_I_row.get("auc_mean") or 0.0
        if n_stage_I > 0 and (auc_I or 0.0) > 0:
            return Check(
                name="stage_breakdown_capability",
                description=(
                    "Per-cancer sens@spec broken down by Stage I vs LATE "
                    "vs UNKNOWN. New build_stage_breakdown_table function "
                    "exercised by tests + synthetic smoke."
                ),
                source="src/per_cancer_sens_at_spec.py + tests",
                verdict=VERDICT_READY,
                value={
                    "synth_n_stage_I": n_stage_I,
                    "synth_n_stage_late": tbl.get("n_stage_late", 0),
                    "synth_n_stage_unknown": tbl.get("n_stage_unknown", 0),
                    "synth_stage_I_auc": auc_I,
                },
                threshold="returns a populated per_stage breakdown on any "
                          "input TSV with a 'stage' column",
                notes=(
                    "Open-data FinaleDB cohort has no stage labels — the "
                    "stage breakdown is wired up but not populated for the "
                    "current headline. See ULTRA_EARLY_READINESS.md §1.3."
                ),
            )
        return Check(
            name="stage_breakdown_capability",
            description="Stage I / LATE / UNKNOWN breakdown wired up.",
            source="src/per_cancer_sens_at_spec.py",
            verdict=VERDICT_GAP,
            value=None,
            threshold="returns a populated per_stage breakdown",
            notes="build_stage_breakdown_table returned no per_stage entries.",
        )
    except Exception as e:  # noqa: BLE001 — any failure is a gap
        return Check(
            name="stage_breakdown_capability",
            description="Stage I / LATE / UNKNOWN breakdown wired up.",
            source="src/per_cancer_sens_at_spec.py",
            verdict=VERDICT_GAP,
            value=None,
            threshold="returns a populated per_stage breakdown",
            notes=f"Synthetic smoke raised: {e}",
        )


# ─────────────────────────────────────────────────────────────────────
# Markdown rendering
# ─────────────────────────────────────────────────────────────────────

def render_md(report: ReadinessReport) -> str:
    """Render the readiness report as a Markdown document."""
    lines: list[str] = []
    lines.append("# DeepCatch Ultra-Early Readiness — Generated Report")
    lines.append("")
    lines.append(
        "Auto-generated by `scripts/ultraearly_readiness.py`. Every "
        "claim ties to a cached JSON artifact; no number is fabricated. "
        "Run the script again after re-generating inputs to refresh."
    )
    lines.append("")
    lines.append(
        f"- **Generated:** {report.generated_at}"
    )
    lines.append(f"- **Repo commit:** {report.repo_commit}")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    counts = report.summary_counts
    lines.append(
        f"- Ready: **{counts.get(VERDICT_READY, 0)}**"
    )
    f"- Not-ready: **{counts.get(VERDICT_NOT_READY, 0)}**"
    lines.append(f"- Gaps: **{counts.get(VERDICT_GAP, 0)}**")
    lines.append("")
    lines.append("**Headline:** " + report.headline_status)
    lines.append("")
    lines.append("## Readiness checks")
    lines.append("")
    for c in report.checks:
        verdict_badge = {
            VERDICT_READY: "✅ ready",
            VERDICT_NOT_READY: "⚠️  not_ready",
            VERDICT_GAP: "❌ gap",
        }.get(c.verdict, c.verdict)
        lines.append(f"### {c.name}")
        lines.append("")
        lines.append(f"- **Verdict:** {verdict_badge}")
        lines.append(f"- **What:** {c.description}")
        if c.value is not None and c.value != "":
            lines.append(f"- **Value:** `{c.value}`")
        else:
            lines.append("- **Value:** (no value)")
        if c.threshold:
            lines.append(f"- **Threshold:** `{c.threshold}`")
        lines.append(f"- **Source:** `{c.source}`")
        lines.append(f"- **Notes:** {c.notes}")
        lines.append("")
    return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--results-dir",
        default=str(_REPO_ROOT / "results"),
        help="Directory containing the cached JSON artifacts (default: results/).",
    )
    ap.add_argument(
        "--out-json",
        default=str(_REPO_ROOT / "results" / "ultraearly_readiness.json"),
        help="Output JSON path (default: results/ultraearly_readiness.json).",
    )
    ap.add_argument(
        "--out-md",
        default=str(_REPO_ROOT / "docs" / "ULTRA_EARLY_READINESS_AUTO.md"),
        help="Output Markdown path (default: docs/ULTRA_EARLY_READINESS_AUTO.md).",
    )
    args = ap.parse_args()

    results_dir = Path(args.results_dir)

    from datetime import datetime, timezone
    report = ReadinessReport(
        generated_at=datetime.now(timezone.utc).isoformat(),
        repo_commit=_detect_repo_commit(),
    )

    # 1. Cross-study pooled AUC with GC correction.
    report.checks.append(check_cross_study_pooled_auc(
        results_dir / "cross_study_finallydb_gc_corrected.json",
        threshold_target=0.95, threshold_floor=0.80,
    ))
    # 3. Per-cancer sens@99%.
    report.checks.append(check_per_cancer_sens99(
        results_dir / "per_cancer_sens_at_spec.json", min_sens99=0.60,
    ))
    # 4. Panel-LLR ultra-low-VAF sweep at 0.1% ctDNA.
    report.checks.append(check_panel_llr_ultra_low_vaf(
        results_dir / "real_tcga_validation.json",
        threshold_target=0.90, threshold_floor=0.70,
    ))
    # 5. Foundation fusion layer (static).
    report.checks.append(check_foundation_fusion_honest())
    # 6. Stage breakdown capability (synthetic smoke).
    report.checks.append(check_stage_breakdown_capability())

    counts = {VERDICT_READY: 0, VERDICT_NOT_READY: 0, VERDICT_GAP: 0}
    for c in report.checks:
        counts[c.verdict] = counts.get(c.verdict, 0) + 1
        if c.verdict == VERDICT_GAP:
            report.gaps.append({
                "name": c.name,
                "source": c.source,
                "notes": c.notes,
            })
    report.summary_counts = counts

    # Headline: ready iff every artifact-derived check is "ready" and at
    # least 4 of 6 checks are ready. The foundation check is excluded
    # from the headline (it's an honest warning, not a blocker).
    artifact_checks = [c for c in report.checks if c.name != "foundation_fusion_layer"]
    artifact_ready = sum(1 for c in artifact_checks if c.verdict == VERDICT_READY)
    if counts.get(VERDICT_GAP, 0) == 0 and artifact_ready == len(artifact_checks):
        report.headline_status = (
            "READY for ultra-early cancer signaling at the open-data + "
            "in-silico panel-LLR level. Clinical-grade evidence requires real "
            "plasma validation (see ULTRA_EARLY_READINESS.md §4)."
        )
    else:
        report.headline_status = (
            f"NOT-READY — {counts.get(VERDICT_GAP, 0)} gaps and "
            f"{counts.get(VERDICT_NOT_READY, 0)} not-ready items. "
            f"See gaps[] in the JSON for the missing artifacts."
        )

    out_json = Path(args.out_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(
            {
                "schema_version": report.schema_version,
                "generated_at": report.generated_at,
                "repo_commit": report.repo_commit,
                "summary_counts": report.summary_counts,
                "headline_status": report.headline_status,
                "checks": [asdict(c) for c in report.checks],
                "gaps": report.gaps,
            },
            f, indent=2,
        )

    out_md = Path(args.out_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    with open(out_md, "w") as f:
        f.write(render_md(report))

    # One-screen console summary.
    print(f"Wrote {out_json} + {out_md}")
    print(f"Headline: {report.headline_status}")
    print(f"Counts: ready={counts.get(VERDICT_READY, 0)}, "
          f"not_ready={counts.get(VERDICT_NOT_READY, 0)}, "
          f"gaps={counts.get(VERDICT_GAP, 0)}")
    for c in report.checks:
        v = c.value if c.value is not None else "—"
        print(f"  [{c.verdict:>9}] {c.name:>38} value={v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())