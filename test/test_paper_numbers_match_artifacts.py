"""The manuscript must quote the numbers the artifacts actually contain.

Why this file exists
---------------------
The benchmark artifacts were regenerated from a seeded pipeline in
v2.3.0, and the prose was updated to match. Nothing enforced that, so
the next regeneration would have silently desynchronised the manuscript
from the evidence again -- which is exactly what had happened before
(the committed JSON came from a PCA-200 run, the document from a PCA-20
run, and the manuscript quoted both).

This test reads the headline numbers **out of the artifacts** and
requires them to appear in the documents that cite them. Regenerate an
artifact without updating the prose and this goes red.

The list of quotable numbers is not written here: it lives in
``scripts/benchmark_numbers.py`` (``QUOTED_NUMBERS`` /
``QUOTED_SPECIAL``), which is also the script that prints every quoted
figure with its artifact. One list, so the test and the printer can
never disagree about what is quotable.

It is deliberately a text-presence check rather than a parser: the point
is not to re-derive the manuscript, it is to make drift impossible to
miss.
"""
from __future__ import annotations

import importlib.util
import io
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))


def _load_numbers_module():
    spec = importlib.util.spec_from_file_location(
        "benchmark_numbers", _REPO_ROOT / "scripts" / "benchmark_numbers.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_numbers = _load_numbers_module()


def _read(rel: str) -> str:
    return io.open(_REPO_ROOT / rel, encoding="utf-8").read()


def _cases():
    """(label, value, documents) for every quoted number."""
    for label, value, _artifact, documents in _numbers.quoted_values():
        yield label, value, documents


@pytest.mark.parametrize("label,value,documents", list(_cases()))
def test_quoted_number_is_present_wherever_it_is_cited(label, value,
                                                       documents):
    rendered = f"{value:.4f}"
    for doc in documents:
        text = _read(doc)
        assert rendered in text, (
            f"{doc} does not quote {rendered} (from {label}). Update the "
            f"prose or the artifact.")


def test_methylation_null_is_quoted_with_its_permutation_p_value():
    """The null must be quoted together with the test that establishes it."""
    tex = _read("paper/paper.tex")
    assert "0.9444" in tex
    assert "103,863" in tex, (
        "the matched-locus CpG count that the null rests on is missing")


def test_paired_fusion_delta_matches_the_artifact():
    """The preprint quotes a paired delta; it must be the artifact's."""
    artifact = _numbers.load("fusion_ablation_isotonic.json")
    if artifact is None:
        pytest.skip("fusion artifact not present")
    per_seed = artifact["per_seed_per_strategy"]
    tn = [s["tumor_naive"]["auc"] for s in per_seed]
    na = [s["naive_average"]["auc"] for s in per_seed]
    delta = sum(b - a for a, b in zip(tn, na)) / len(tn)
    rendered = f"{delta:.4f}".lstrip("-")
    assert rendered in _read("paper/paper.tex"), (
        f"the paired fusion delta quoted in the preprint is not the one in the "
        f"artifact ({rendered})")


# Figures that the v2.3.0 regeneration superseded. They must not reappear
# in a current-facing document; the CHANGELOG and the dated audit records
# are allowed to name them, because describing the old value is the
# point of those files.
SUPERSEDED = [
    "0.9670",   # the old GC-corrected pooled AUC
    "0.9747",   # the old legacy pooled AUC
    "0.9753",   # the companion pipeline's unseeded tumor-naive figure
    # 0.9886 is deliberately absent: it was the old naive-average fusion
    # AUC, but it is also the CURRENT per-cancer OV AUC, so a string check
    # cannot distinguish them. The presence check covers the fusion
    # number instead (the documents must quote 0.9884).
    "0.4657",   # the old shuffled-label null
    "0.9225",   # the old PCA-20 pooled AUC
    "0.9458",   # the old PCA-20 per-publication AUC
]
CURRENT_DOCS = [
    "README.md",
    "MODEL.md",
    "REVIEWERS.md",
    "TEAM.md",
    "USAGE.md",
    "paper/paper.tex",
    "paper/PAPER.md",
    "paper/METHODS_PAPER.md",
    "paper/README.md",
    "paper.md",
    "docs/CROSS_STUDY_BENCHMARK.md",
    "docs/CROSS_STUDY_BENCHMARK_GC_CORRECTED.md",
    "docs/CROSS_STUDY_BENCHMARK_with_motifs.md",
    "docs/CROSS_STUDY_LOSO.md",
    "docs/EXTERNAL_VALIDATION.md",
    "docs/FUSION_ISOTONIC.md",
    "docs/HARMONIZATION.md",
    "docs/METHYLATION_PIPELINE.md",
    "docs/site/index.html",
]


@pytest.mark.parametrize("value", SUPERSEDED)
def test_superseded_figures_do_not_reappear(value):
    """A regenerated number must not leave its predecessor behind."""
    for doc in CURRENT_DOCS:
        path = _REPO_ROOT / doc
        if not path.exists():
            continue
        text = io.open(path, encoding="utf-8").read()
        assert value not in text, (
            f"{doc} still quotes the superseded figure {value}")
