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
SUPERSEDED = _numbers.SUPERSEDED


def _load_abstract_module():
    """The renderer that turns the manuscript's abstract into paste-ready text."""
    spec = importlib.util.spec_from_file_location(
        "generate_preprint_abstract",
        _REPO_ROOT / "scripts" / "generate_preprint_abstract.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_abstract = _load_abstract_module()


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


def test_the_guard_can_see_a_thousands_separated_integer():
    """A count written "103,863" must be visible to the number guard.

    The guard's parser matched only DECIMAL numbers, so a figure the
    documents quote with a thousands separator was invisible to it: it
    could drift freely while the guard stayed green. That is the same
    blind spot as a figure re-rendered at a coarser precision -- a
    rendering the guard cannot parse.

    A count like this is written EXACTLY, so it is matched exactly:
    an off-by-one must not pass.
    """
    assert _numbers.quotes("over 103,863 loci", 103863.0), (
        "the guard cannot see a thousands-separated integer")
    assert not _numbers.quotes("over 103,863 loci", 103864.0), (
        "a thousands-separated integer was matched with a tolerance")
    assert not _numbers.quotes("over 103,863 loci", 103.0), (
        "the comma-separated count was parsed as its leading part")
    assert not _numbers.quotes("over 103,863 loci", 863.0), (
        "the comma-separated count was parsed as its trailing part")
    # A plain integer without a separator is NOT a benchmark rendering:
    # the guard's convention is three decimals or an exact separated
    # count, nothing looser. "n = 103863" must not be read as a
    # figure that pins a value.
    assert not _numbers.quotes("n = 103863", 103863.0), (
        "an unseparated integer was treated as an exact benchmark figure")


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


@pytest.mark.parametrize("value,documents",
                         sorted(_numbers.scoped_superseded_documents().items()))
def test_scoped_superseded_figures_do_not_reappear(value, documents):
    """A superseded value is banned only where it would be a regression.

    Several superseded values are ALSO the rendering of a legitimate
    current figure -- 0.9915 is the decision curve's naive-average AUC,
    0.9670 is the OV per-cancer DeLong CI upper bound, 0.9701 one of the
    pooled per-seed AUCs. A global literal ban would forbid a correct
    document, so each value is banned only in the documents that quote
    the quantity it replaced.
    """
    for doc in documents:
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

    The length is measured on the string that is actually pasted: the
    abstract rendered out of its LaTeX with the maths intact. The
    previous measure deleted the maths (`$...$` became `0`), so it
    under-counted and could not see an abstract that was long because
    of its formulae.
    """
    length = _abstract_length()
    assert length <= 3000, (
        f"the abstract is {length} characters; bioRxiv rejects anything "
        f"over 3000")
    assert length <= 2950, (
        f"the abstract is {length} characters -- only {3000 - length} of "
        "margin left under bioRxiv's limit")


def _abstract_length() -> int:
    """The abstract's length, measured the way the test above measures it.

    This is the length of the text that goes into the bioRxiv form: the
    LaTeX rendered to plain text -- not the LaTeX source, and not a lossy
    normalisation of it. `scripts/generate_preprint_abstract.py` is the
    same renderer the paste-ready document uses, so the form, that
    document and this test cannot disagree about the number.
    """
    tex = _read("paper/paper.tex")
    return _abstract.abstract_length(tex)


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


def test_the_methylation_cohort_split_is_the_artifacts_split():
    """The documents must state the cohort's case/control split it has.

    Six documents said the 31 real samples were "17 colorectal cancer
    and 14 liver-transplant controls"; the artifact holds 18 and 13.
    The composition is a pair of integers, which the number guard cannot
    see (it matches decimals), so it is checked here.
    """
    artifact = _numbers.load("methylation_cohort_channel.json")
    if artifact is None:
        pytest.skip("methylation artifact not present")
    per_sample = artifact["cohort"]["per_sample"]
    cases = sum(1 for row in per_sample.values() if row.get("label") == 1)
    controls = sum(1 for row in per_sample.values() if row.get("label") == 0)
    assert cases and controls, "the artifact carries no labelled samples"

    for doc in ("paper/paper.tex", "paper/METHODS_PAPER.md", "README.md",
                "ULTRA_EARLY_READINESS.md", "docs/METHYLATION_PIPELINE.md"):
        text = _read(doc)
        assert f"{cases} colorectal" in text, (
            f"{doc} does not state the artifact's case count ({cases} "
            f"colorectal cancer)")
        assert f"{controls} liver-transplant" in text, (
            f"{doc} does not state the artifact's control count ({controls} "
            f"liver-transplant controls)")

    chinese = _read("README.zh-CN.md")
    assert f"{cases} 例结直肠癌" in chinese, (
        f"README.zh-CN.md does not state the artifact's case count ({cases})")
    assert f"{controls} 例肝移植" in chinese, (
        f"README.zh-CN.md does not state the artifact's control count "
        f"({controls})")

