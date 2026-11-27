"""The decision-curve CLI must read the study column from the header.

What went wrong
---------------
`src/fragmentomics/decision_curve_cli.py` resolved the LABEL with the
header-aware resolver (`load_labels_tsv`) but read the STUDY
positionally -- column 2 of each line. In the canonical FinaleDB layout
(`sample  disease_class  label  study`) column 2 is the *label*, so the
study became "cancer"/"healthy". The per-study harmonisation then
z-scored each class separately, which removes the signal the classifier
is supposed to find, and the run silently reported a chance-level AUC
(0.4841 against 0.9773) instead of failing.

The adapter's own docstring warns about exactly this ("reading the
study positionally breaks on the second (it picks up the label)"), and
the bug was invisible because the only labels file the CLI was ever run
against was the headerless `labels_cross_study.tsv`, where the
positional read happens to be right. Pointed at the canonical
`labels_multiclass.tsv` -- the file every other script in the
repository reads -- it produced a nonsense result.

This test builds a small synthetic cohort whose features carry a
class signal, writes the labels in the canonical headered layout, and
requires the CLI to recover that signal. It then requires the headered
and headerless layouts to give the same answer, because the tool must
not depend on which file it is handed.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

DIM_5MB = 40
DIM_100KB = 40
DIM_FSD = 30


def _write_sample_features(features_dir: Path, sample_id: str, label: int,
                           rng: np.random.Generator) -> None:
    """Write one sample's five feature files (the adapter's layout)."""
    centre = 0.6 if label == 1 else 0.4
    np.save(features_dir / f"{sample_id}.delfi_5mb_ratio.npy",
            rng.normal(centre, 0.2, DIM_5MB))
    np.save(features_dir / f"{sample_id}.delfi_5mb_coverage.npy",
            rng.normal(centre, 0.2, DIM_5MB))
    np.save(features_dir / f"{sample_id}.delfi_100kb_ratio.npy",
            rng.normal(centre, 0.2, DIM_100KB))
    np.save(features_dir / f"{sample_id}.delfi_100kb_counts.npy",
            rng.normal(centre, 0.2, DIM_100KB))
    size_bins = {
        f"{100 + 4 * i}-{104 + 4 * i}": float(rng.normal(0.5, 0.1))
        for i in range(DIM_FSD)
    }
    with open(features_dir / f"{sample_id}.fsd.json", "w") as f:
        json.dump({"size_bins": size_bins}, f)


def _build_cohort(tmp_path: Path):
    """Two studies, 20 samples each (10 cancer, 10 healthy)."""
    feat_dir = tmp_path / "features"
    feat_dir.mkdir()
    rng = np.random.default_rng(7)
    rows = []  # (sample_id, disease_class, label, study)
    for study in ("studyA", "studyB"):
        for i in range(10):
            rows.append((f"{study}C{i:02d}", "CANCER", "cancer", study))
            rows.append((f"{study}H{i:02d}", "HEALTHY", "healthy", study))
    for sample_id, _disease, label, _study in rows:
        _write_sample_features(feat_dir, sample_id, 1 if label == "cancer" else 0,
                               rng)
    return feat_dir, rows


def _write_labels(path: Path, rows, headered: bool) -> None:
    lines = ["sample\tdisease_class\tlabel\tstudy"] if headered else []
    for sample_id, disease, label, study in rows:
        if headered:
            lines.append(f"{sample_id}\t{disease}\t{label}\t{study}")
        else:
            lines.append(f"{sample_id}\t{label}\t{study}")
    path.write_text("\n".join(lines) + "\n")


def _run_cli(feat_dir: Path, labels: Path, out: Path) -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(_REPO_ROOT / "src")
    proc = subprocess.run(
        [sys.executable, "-u",
         str(_REPO_ROOT / "src" / "fragmentomics" / "decision_curve_cli.py"),
         "--features-dir", str(feat_dir),
         "--labels", str(labels),
         "--pca-n", "10",
         "--seeds", "2",
         "--out", str(out)],
        capture_output=True, text=True, cwd=str(_REPO_ROOT), env=env,
        timeout=600,
    )
    assert proc.returncode == 0, (
        f"the decision-curve CLI failed (rc={proc.returncode}):\n"
        f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
    )
    return json.loads(out.read_text())


def _run_module_cli(module: str, feat_dir: Path, labels: Path,
                    out: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, "-m", module,
         "--features-dir", str(feat_dir),
         "--labels", str(labels),
         "--pca-n", "10",
         "--seeds", "2",
         "--out", str(out)],
        capture_output=True, text=True, cwd=str(_REPO_ROOT),
        timeout=600,
    )
    assert proc.returncode == 0, (
        f"{module} CLI failed (rc={proc.returncode}):\n"
        f"stdout:\n{proc.stdout[-2000:]}\nstderr:\n{proc.stderr[-2000:]}"
    )
    return json.loads(out.read_text())


def _fusion_auc(payload: dict) -> float:
    return next(
        m["mean"] for m in payload["per_strategy"]["tumor_naive"]
        if m["metric"] == "auc"
    )


def test_the_headered_labels_layout_keeps_the_signal(tmp_path):
    """The canonical `sample disease_class label study` layout must work.

    With the study column read positionally, the harmonisation groups
    by class instead of by study, the class signal is z-scored away and
    the AUC collapses to chance.
    """
    feat_dir, rows = _build_cohort(tmp_path)
    labels = tmp_path / "labels_multiclass.tsv"
    _write_labels(labels, rows, headered=True)

    payload = _run_cli(feat_dir, labels, tmp_path / "out.json")
    auc = payload["per_strategy"]["tumor_naive"]["auc"]
    assert auc > 0.7, (
        f"the tumor-naive AUC is {auc:.4f} on a cohort with a clear class "
        f"signal -- the study column was probably read positionally, so "
        f"the harmonisation grouped by class and removed the signal"
    )

    tn = _run_module_cli("src.fragmentomics.train_tumor_naive",
                         feat_dir, labels, tmp_path / "tn.json")
    assert tn["studies"] == ["studyA", "studyB"], (
        f"train_tumor_naive parsed studies={tn['studies']}; expected real "
        f"study names rather than label values"
    )
    assert tn["auc_mean"] > 0.7, (
        f"train_tumor_naive AUC is {tn['auc_mean']:.4f} on a cohort with a "
        f"clear class signal"
    )

    fusion = _run_module_cli("src.fragmentomics.fusion_ablation",
                             feat_dir, labels, tmp_path / "fusion.json")
    assert _fusion_auc(fusion) > 0.7, (
        f"fusion_ablation tumor_naive AUC is {_fusion_auc(fusion):.4f} on a "
        f"cohort with a clear class signal"
    )


def test_the_two_label_layouts_agree(tmp_path):
    """The same cohort must give the same answer in either layout."""
    feat_dir, rows = _build_cohort(tmp_path)
    headered = tmp_path / "labels_headered.tsv"
    headerless = tmp_path / "labels_headerless.tsv"
    _write_labels(headered, rows, headered=True)
    _write_labels(headerless, rows, headered=False)

    a = _run_cli(feat_dir, headered, tmp_path / "a.json")
    b = _run_cli(feat_dir, headerless, tmp_path / "b.json")
    assert a["per_strategy"]["tumor_naive"]["auc"] == pytest.approx(
        b["per_strategy"]["tumor_naive"]["auc"], abs=1e-9), (
        "the headered and headerless labels files describe the same cohort "
        "but the decision curve gave different answers"
    )

    tn_a = _run_module_cli("src.fragmentomics.train_tumor_naive",
                           feat_dir, headered, tmp_path / "tn_a.json")
    tn_b = _run_module_cli("src.fragmentomics.train_tumor_naive",
                           feat_dir, headerless, tmp_path / "tn_b.json")
    assert tn_a["studies"] == tn_b["studies"] == ["studyA", "studyB"]
    assert tn_a["auc_mean"] == pytest.approx(tn_b["auc_mean"], abs=1e-9), (
        "the headered and headerless labels files describe the same cohort "
        "but train_tumor_naive gave different answers"
    )

    fus_a = _run_module_cli("src.fragmentomics.fusion_ablation",
                            feat_dir, headered, tmp_path / "fus_a.json")
    fus_b = _run_module_cli("src.fragmentomics.fusion_ablation",
                            feat_dir, headerless, tmp_path / "fus_b.json")
    assert _fusion_auc(fus_a) == pytest.approx(_fusion_auc(fus_b), abs=1e-9), (
        "the headered and headerless labels files describe the same cohort "
        "but fusion_ablation gave different answers"
    )
