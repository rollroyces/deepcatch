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
    # A document may quote a number at the precision its table uses, so
    # the check compares each written number against the artifact value
    # at the precision it was written at. A bare substring search would
    # not do: "0.970" cannot tell the GC-corrected publication-8 AUC
    # (0.9698) from the legacy one (0.9704), and both round to it.
    for doc in documents:
        text = _read(doc)
        assert _numbers.quotes(text, value), (
            f"{doc} does not quote {value:.6f} (from {label}). Update the "
            f"prose or the artifact.")


@pytest.mark.parametrize("doc,anchor,required,forbidden,source",
                         list(_numbers.table_row_values()))
def test_table_row_quotes_one_run_only(doc, anchor, required, forbidden,
                                       source):
    """A table row must not mix numbers from two different runs.

    Every number in a row can exist somewhere in the repository while the
    row itself mixes two artifacts -- the document-level checks above
    stay green because each number is present somewhere. This finds the
    row that contains the anchor and requires it to quote the right
    numbers, to not quote the same quantity from another run, and to
    name that run's artifact in its source column.
    """
    text = _read(doc)
    lines = [line for line in text.splitlines() if anchor in line]
    assert lines, f"{doc} has no table row containing {anchor!r}"
    assert len(lines) == 1, (
        f"{doc} has {len(lines)} rows containing {anchor!r}; the guard "
        f"cannot tell which one to check")
    row = lines[0]
    for artifact_name, value in required:
        assert _numbers.quotes(row, value), (
            f"{doc}: the row {anchor!r} does not quote {value:.6f} from "
            f"{artifact_name}. It quotes: {row.strip()}")
    for artifact_name, value in forbidden:
        assert not _numbers.quotes(row, value), (
            f"{doc}: the row {anchor!r} quotes {value:.6f} from "
            f"{artifact_name}, a different run than the rest of the row. "
            f"It quotes: {row.strip()}")
    if source:
        assert source in row, (
            f"{doc}: the row {anchor!r} must name {source} as its source. "
            f"It quotes: {row.strip()}")


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
#
# The check is a literal one, so a value that merely shares a prefix
# with a current figure does not trip it: the old GC-corrected pooled
# AUC "0.9670" is banned, but the per-cancer PAAD AUC "0.967" and the
# clean-panel LLR AUC "0.922" are different numbers and stay.
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
    # Superseded by the v2.3.0 correction of the external validation:
    # the first fetch took 9 of publication 1's 60 samples, and the
    # harmonisation was a single pooled z-score rather than the
    # per-publication one the script documents.
    "0.9915",   # the old publication-7 external AUC
    "0.8819",   # the old pooled external AUC
    "0.4454",   # the old publication-9 external sensitivity
    # Superseded by the seeded LOSO run.
    "0.9701",   # the old in-domain Cristiano AUC
    # Superseded by the GC-corrected true-confound control. It was still
    # hard-coded in the readiness generator until this release.
    "0.499",
    # The old mutation-only fusion AUC, from the pre-regeneration
    # fusion artifact. It is the only one of the stale fusion figures
    # that can be banned literally: the others ("0.9887", "0.937",
    # "0.845") are prefixes of legitimate current numbers, so they are
    # pinned by the table-row checks below instead.
    "0.9242",
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
    # CITATION.cff is what GitHub's "Cite this repository" button hands a
    # reader, and it was still quoting the superseded pooled AUC (0.9670)
    # because no check covered a .cff file at all.
    "CITATION.cff",
    # The Simplified-Chinese README carried the same stale generation as
    # the English one (the fusion table from a pre-regeneration run) with
    # no check covering it at all.
    "README.zh-CN.md",
    # RESULTS.md is the repository's headline results document and it
    # was covered by no check at all: its TL;DR quoted the legacy
    # tumor-naive benchmark, a fusion AUC the fusion artifact no longer
    # holds, and a "0.499 -> 0.500" harmonization pair that no artifact
    # backs, while its Section 4 decision-curve table predated the
    # seeded PCA and could not be reproduced by the committed code.
    "RESULTS.md",
]


@pytest.mark.parametrize("value", SUPERSEDED)
def test_superseded_figures_do_not_reappear(value):
    """A regenerated number must not leave its predecessor behind.

    The check is a literal one. It has to be: the superseded values are
    recorded at four decimals, and a current figure that merely shares
    their first three ("0.967" for the per-cancer PAAD AUC, "0.922" for
    the clean-panel LLR AUC) is a different number, not the old one.
    """
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


def test_the_chinese_readme_carries_every_figure_the_english_one_does():
    """The translation must not omit a figure the original states.

    README.zh-CN.md is a translation of README.md. It drifted a whole
    generation once, quoting a pre-regeneration fusion table and a test
    badge from two releases earlier, and nothing checked it. This
    requires the translation to carry every number the English original
    states, so a figure added to the original cannot be forgotten in
    the translation.

    Numbers are compared as written: a figure quoted at three decimals
    in one file and four in the other is still the same figure, so the
    check is on the value, not the string.
    """
    import re

    def figures(text):
        out = set()
        for match in re.finditer(r"(?<![\d.])\d+\.\d{3,}(?![\d])", text):
            out.add(round(float(match.group(0)), 4))
        return out

    english = _read("README.md")
    chinese = _read("README.zh-CN.md")
    missing = figures(english) - figures(chinese)
    assert not missing, (
        "README.zh-CN.md does not carry these figures from README.md: "
        f"{sorted(missing)}")


def test_citation_metadata_agrees_with_the_manuscript():
    """CITATION.cff must describe the paper that is actually being published.

    GitHub's "Cite this repository" button hands this file to a reader,
    and nothing checked it: it quoted a title the manuscript does not
    use, and a pooled AUC (0.9670) that the v2.3.0 regeneration had
    already superseded. The title is now the manuscript's own, and the
    figures it quotes are checked against the artifacts by the
    superseded list above.
    """
    tex = _read("paper/paper.tex")
    match = re.search(r"\\title\{(.*?)\}\s*\n", tex, re.S)
    assert match, "paper.tex has no title"
    title = re.sub(r"\\textbf\{(.*)\}", r"\1", match.group(1)).strip()
    title = " ".join(title.split())
    cff = _read("CITATION.cff")
    stated = re.search(r'^title:\s*"(.*?)"', cff, re.M)
    assert stated, "CITATION.cff has no title"
    assert stated.group(1) == title, (
        "CITATION.cff quotes a different title from the manuscript:\n"
        f"  CITATION.cff: {stated.group(1)}\n"
        f"  paper.tex:     {title}")

    # The affiliation the manuscript declares is the one the citation
    # metadata must carry -- a reviewer clicks through from one to the
    # other.
    assert "Independent Researcher" in cff, (
        "CITATION.cff no longer declares the affiliation the manuscript "
        "declares")
    assert "Independent Researcher" in tex, (
        "paper.tex no longer declares an affiliation; CITATION.cff must "
        "be updated to match")


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
