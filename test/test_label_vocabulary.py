"""Every reader of a labels file must resolve labels the same way.

What went wrong
---------------
Two readers of the same labels file resolved the label cell by exact
string equality -- ``1 if cell == "cancer" else 0`` -- so every other
spelling ("Cancer", "tumor", "control") became a **healthy control**,
silently. One of them writes the headline cross-study AUC; the other
writes the external validation. Reading the same file through the
package's own resolver gave a *different* label vector, which is what
exposed it.

The collaborator interface promises a schema error for a label that is
neither a case nor a control; the adapter did not deliver one, so a
collaborator's "tumor" rows would have been filed as healthy donors.

The vocabulary now lives in one place (`src/fragmentomics/labels.py`)
and every reader goes through it: an unknown label is an error to be
fixed, never a control to be invented.
"""
from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from fragmentomics.labels import (  # noqa: E402
    LabelVocabularyError, resolve_label,
)


def _load(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("raw,expected", [
    ("cancer", 1), ("Cancer", 1), ("CANCER", 1), (" tumor ", 1),
    ("positive", 1), ("case", 1), ("1", 1), ("true", 1),
    ("healthy", 0), ("Healthy", 0), ("control", 0), ("normal", 0),
    ("0", 0), ("false", 0),
])
def test_the_shared_resolver_reads_the_case_and_control_words(raw, expected):
    assert resolve_label(raw) == expected


@pytest.mark.parametrize("raw", ["CRC", "unknown", "0.5", "", "yes"])
def test_the_shared_resolver_refuses_anything_else(raw):
    with pytest.raises(LabelVocabularyError):
        resolve_label(raw)


def _write_labels(path: Path, rows, headered=True):
    lines = []
    if headered:
        lines.append("sample\tdisease_class\tlabel\tstudy")
    for row in rows:
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")


ROWS = [
    ("C001", "CRC", "cancer", "jiang"),
    ("C002", "CRC", "Cancer", "jiang"),
    ("C003", "CRC", "tumor", "jiang"),
    ("H001", "HEALTHY", "healthy", "cristiano"),
    ("H002", "HEALTHY", "control", "cristiano"),
]


def test_the_cross_study_reader_resolves_labels_like_the_package(tmp_path):
    """The headline reader and the package's resolver must agree."""
    csf = _load("csf_label_test", str(_REPO_ROOT / "scripts" /
                                     "cross_study_finallydb.py"))
    path = tmp_path / "labels.tsv"
    _write_labels(path, ROWS)

    labels, studies, _disease, _pub = csf.load_labels_multiclass(
        str(path), requested_publications={"6", "8"})
    assert labels == {"C001": 1, "C002": 1, "C003": 1,
                      "H001": 0, "H002": 0}


def test_the_cross_study_reader_refuses_an_unknown_label(tmp_path):
    csf = _load("csf_label_test2", str(_REPO_ROOT / "scripts" /
                                      "cross_study_finallydb.py"))
    path = tmp_path / "labels.tsv"
    _write_labels(path, [("C001", "CRC", "CRC", "jiang")])

    with pytest.raises(LabelVocabularyError):
        csf.load_labels_multiclass(str(path),
                                   requested_publications={"6", "8"})


def test_the_external_validation_reader_refuses_an_unknown_label(tmp_path):
    ext = _load("ext_label_test", str(_REPO_ROOT / "scripts" /
                                      "external_validation_finaledb.py"))
    path = tmp_path / "labels.tsv"
    _write_labels(path, [("C001", "CRC", "CRC", "jiang")])

    # The label reading is exercised on its own: the cohort loader
    # imports the sibling pipeline, which is absent in CI.
    with pytest.raises(LabelVocabularyError):
        ext._read_labels(str(path), {"6", "8"})


def test_the_collaborator_adapter_reports_an_unknown_label_as_a_schema_error(
        tmp_path):
    """The collaborator interface promises a schema error; deliver one."""
    adapter = _load("adapter_label_test", str(_REPO_ROOT / "scripts" /
                                             "adapter_local_cohort.py"))
    path = tmp_path / "labels.tsv"
    _write_labels(path, [("C001", "CRC", "CRC", "jiang")])

    with pytest.raises(adapter.SchemaError) as excinfo:
        adapter.load_labels_tsv(str(path))
    assert "CRC" in str(excinfo.value)


def test_the_collaborator_adapter_normalises_the_case_and_control_words(
        tmp_path):
    adapter = _load("adapter_label_test2", str(_REPO_ROOT / "scripts" /
                                              "adapter_local_cohort.py"))
    path = tmp_path / "labels.tsv"
    _write_labels(path, ROWS)

    rows = adapter.load_labels_tsv(str(path))
    assert {r["sample_id"]: r["label"] for r in rows} == {
        "C001": "cancer", "C002": "cancer", "C003": "cancer",
        "H001": "healthy", "H002": "healthy",
    }


# The three documented layouts. Every reader must answer identically on
# each of them.
LAYOUTS = {
    "canonical headered": (
        "sample\tdisease_class\tlabel\tstudy\n"
        "C001\tCRC\tcancer\tjiang\n"
        "C002\tCRC\ttumor\tjiang\n"
        "H001\tHEALTHY\thealthy\tcristiano\n"
        "H002\tHEALTHY\tcontrol\tcristiano\n"
    ),
    "headerless three-column": (
        "C001\tcancer\tjiang\n"
        "C002\ttumor\tjiang\n"
        "H001\thealthy\tcristiano\n"
        "H002\tcontrol\tcristiano\n"
    ),
    "headerless four-column": (
        "C001\tCRC\tcancer\tjiang\n"
        "C002\tCRC\ttumor\tjiang\n"
        "H001\tHEALTHY\thealthy\tcristiano\n"
        "H002\tHEALTHY\tcontrol\tcristiano\n"
    ),
}

EXPECTED = {"C001": 1, "C002": 1, "H001": 0, "H002": 0}


def _read_with_every_reader(path):
    """The label vector each reader produces for the same file."""
    csf = _load("csf_agree", str(_REPO_ROOT / "scripts" /
                                "cross_study_finallydb.py"))
    ext = _load("ext_agree", str(_REPO_ROOT / "scripts" /
                                 "external_validation_finaledb.py"))
    adapter = _load("adapter_agree", str(_REPO_ROOT / "scripts" /
                                        "adapter_local_cohort.py"))
    from fragmentomics.tumor_naive_adapter import load_labels_and_studies

    cross, _studies, _disease, _pub = csf.load_labels_multiclass(
        str(path), requested_publications={"6", "8"})
    external, _s, _d, _p = ext._read_labels(str(path), {"6", "8"})
    collaborator = {r["sample_id"]: (1 if r["label"] == "cancer" else 0)
                    for r in adapter.load_labels_tsv(str(path))}
    package, _s = load_labels_and_studies(str(path))
    return {
        "cross_study_finallydb": cross,
        "external_validation": external,
        "adapter_local_cohort": collaborator,
        "tumor_naive_adapter": package,
    }


@pytest.mark.parametrize("layout", sorted(LAYOUTS))
def test_every_reader_resolves_the_same_file_identically(tmp_path, layout):
    """Four readers, one file, one answer.

    They used to disagree: the package's reader took the DISEASE CLASS
    for a headerless four-column file, and the two script readers
    skipped a three-column file entirely -- silently, returning an
    empty cohort. Reading the same file through every reader and
    comparing the OUTPUT is what exposes this class.
    """
    path = tmp_path / "labels.tsv"
    path.write_text(LAYOUTS[layout])

    vectors = _read_with_every_reader(path)
    for reader, vector in vectors.items():
        assert vector == EXPECTED, (
            f"{reader} read the {layout} layout as {vector}, not {EXPECTED}")
