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
import os
import re
import subprocess
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
    "ULTRA_EARLY_READINESS.md",
    "docs/EXTERNAL_VALIDATION.md",
    "docs/CROSS_STUDY_LOSO.md",
    "docs/PORTFOLIO_ROADMAP.md",
    "docs/METHYLATION_PIPELINE.md",
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
    "ULTRA_EARLY_READINESS.md",
    "docs/PORTFOLIO_ROADMAP.md",
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


def test_figures_are_rendered_from_the_headline_artifact():
    """The figures must plot the benchmark the manuscript quotes.

    `scripts/plot_figures.py` once read the legacy no-GC artifact while
    the text quoted the GC-corrected run, so Figure 1 in the submission
    PDF showed a 0.9770 curve next to a 0.9684 sentence. This asserts the
    script's pooled-artifact constant is the headline one.
    """
    spec = importlib.util.spec_from_file_location(
        "plot_figures", _REPO_ROOT / "scripts" / "plot_figures.py")
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:  # matplotlib absent in a minimal env
        pytest.skip(f"plot_figures not importable here: {exc}")
    assert Path(module.POOLED_JSON).name == (
        "cross_study_finallydb_gc_corrected.json"), (
        "the figures must be rendered from the headline (GC-corrected) "
        f"benchmark, not {Path(module.POOLED_JSON).name}")
    assert Path(module.PER_CANCER_JSON).name == "per_cancer_sens_at_spec.json"


def test_abstract_fits_the_biorxiv_character_limit():
    """bioRxiv rejects an abstract over 3000 characters.

    The abstract grew with the v2.3.0 results and passed 3000 without
    anyone checking; the upload would have been refused. This keeps the
    margin visible.
    """
    tex = _read("paper/paper.tex")
    match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    assert match, "paper.tex has no abstract"
    plain = re.sub(r"\\[a-zA-Z]+\{([^}]*)\}", r"\1", match.group(1))
    plain = re.sub(r"\\[a-zA-Z]+", "", plain)
    plain = re.sub(r"\$[^$]*\$", "0", plain)
    length = len(" ".join(plain.split()))
    assert length <= 3000, (
        f"the abstract is {length} characters; bioRxiv rejects anything "
        f"over 3000")
    assert length <= 2950, (
        f"the abstract is {length} characters -- only {3000 - length} of "
        "margin left under bioRxiv's limit")


def _abstract_length() -> int:
    """The abstract's length, measured the way the test above measures it."""
    tex = _read("paper/paper.tex")
    match = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S)
    assert match, "paper.tex has no abstract"
    plain = re.sub(r"\\[a-zA-Z]+\{([^}]*)\}", r"\1", match.group(1))
    plain = re.sub(r"\\[a-zA-Z]+", "", plain)
    plain = re.sub(r"\$[^$]*\$", "0", plain)
    return len(" ".join(plain.split()))


def test_changelog_states_the_current_abstract_length():
    """The release note tells the reader how long the shipped abstract is.

    The abstract is edited after the note is written often enough -- it
    went from 2894 to 2897 when the dataset names were made precise -- and
    a release note that states the wrong length sends the author to
    bioRxiv with the wrong number in hand.
    """
    length = _abstract_length()
    changelog = _read("CHANGELOG.md")
    match = re.search(r"it is now (\d+) with every\s+number intact", changelog)
    assert match, (
        "the CHANGELOG no longer states the abstract's length; update it "
        "and this test")
    stated = int(match.group(1))
    assert stated == length, (
        f"the CHANGELOG says the abstract is {stated} characters; "
        f"paper/paper.tex makes it {length}")


def test_biorxiv_form_agrees_with_the_manuscript():
    """The upload form must quote the manuscript's own length and numbers.

    The form is gitignored (it carries the corresponding author's
    details), so this only runs on the author's machine -- which is where
    it matters, because that is where the form is filled in.
    """
    form = _REPO_ROOT / "paper" / "BIORXIV_SUBMISSION.md"
    if not form.is_file():
        pytest.skip(
            "paper/BIORXIV_SUBMISSION.md is gitignored -- the check runs on "
            "the author's machine")
    checker = _REPO_ROOT / "scripts" / "check_biorxiv_form.py"
    proc = subprocess.run(
        [sys.executable, str(checker)],
        capture_output=True, text=True, timeout=60,
        env={k: v for k, v in os.environ.items() if k != "PYTHONPATH"},
    )
    assert proc.returncode == 0, (
        f"the bioRxiv submission form has drifted from the manuscript:\n"
        f"{proc.stdout}\n{proc.stderr}")
