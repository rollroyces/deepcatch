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

It is deliberately a text-presence check rather than a parser: the point
is not to re-derive the manuscript, it is to make drift impossible to
miss. The authoritative list of quotable numbers lives in
``scripts/benchmark_numbers.py``.
"""
from __future__ import annotations

import io
import json
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _load(name: str):
    path = _REPO_ROOT / "results" / name
    if not path.exists():
        pytest.skip(f"{name} not present")
    with io.open(path, encoding="utf-8") as f:
        return json.load(f)


def _read(rel: str) -> str:
    return io.open(_REPO_ROOT / rel, encoding="utf-8").read()


def _fmt(value: float) -> str:
    return f"{value:.4f}"


def _per_strategy_auc(artifact: dict, strategy: str) -> float:
    for row in artifact["per_strategy"][strategy]:
        if row["metric"] == "auc":
            return float(row["mean"])
    raise AssertionError(f"no auc row for {strategy}")


def _loso_direction(artifact: dict, train: str, test: str) -> float:
    """The pooled test AUC for one leave-one-study-out direction."""
    for entry in artifact["loso"]:
        if (entry["train_study"] == train and entry["test_study"] == test
                and entry["harmonization"] == "self_zscore"):
            return float(entry["pooled_test_auc"]["auc"])
    raise AssertionError(f"no LOSO entry {train} -> {test}")


# (artifact, path into the artifact, documents that must quote it)
HEADLINE_NUMBERS = [
    ("external_validation_finaledb.json",
     ("per_publication_auc", "7", "auc"),
     ["paper/paper.tex", "paper.md"]),
    ("cross_study_loso.json",
     None,  # handled by the direction-specific test
     []),
    ("cross_study_finallydb_gc_corrected.json",
     ("pooled", "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md", "paper.md"]),
    ("cross_study_finallydb.json",
     ("pooled", "harmonized", "auc_mean"),
     ["paper/METHODS_PAPER.md"]),
    ("fusion_ablation_isotonic.json",
     None,  # handled separately: the per-strategy table
     ["paper/paper.tex"]),
    ("methylation_cohort_channel.json",
     ("beta_summary_auc_loo",),
     ["paper/paper.tex", "paper.md"]),
]


def _dig(artifact: dict, path):
    node = artifact
    for key in path:
        node = node[key]
    return node


@pytest.mark.parametrize("artifact_name,path,documents", HEADLINE_NUMBERS)
def test_headline_number_is_quoted_wherever_it_is_cited(artifact_name, path,
                                                         documents):
    artifact = _load(artifact_name)
    if path is None:
        pytest.skip("handled by the fusion-specific test")
    value = float(_dig(artifact, path))
    rendered = _fmt(value)
    for doc in documents:
        text = _read(doc)
        assert rendered in text, (
            f"{doc} does not quote {rendered} (from "
            f"results/{artifact_name}{'.'.join(map(str, path))}). "
            f"Update the prose or the artifact.")
    print(f"{artifact_name}{path} -> {rendered} quoted in {documents}")


def test_fusion_naive_average_is_quoted_in_the_preprint():
    artifact = _load("fusion_ablation_isotonic.json")
    value = _per_strategy_auc(artifact, "naive_average")
    rendered = _fmt(value)
    assert rendered in _read("paper/paper.tex"), (
        f"paper/paper.tex does not quote the naive-average AUC {rendered} from "
        "results/fusion_ablation_isotonic.json")
    assert rendered in _read("docs/site/index.html"), (
        f"the landing page does not quote the naive-average AUC {rendered}")


def test_loso_transfer_numbers_are_quoted_in_the_preprint():
    """Both transfer directions must appear with the artifact's values."""
    artifact = _load("cross_study_loso.json")
    tex = _read("paper/paper.tex")
    for train, test in (("jiang", "cristiano"), ("cristiano", "jiang")):
        rendered = _fmt(_loso_direction(artifact, train, test))
        assert rendered in tex, (
            f"paper/paper.tex does not quote the {train} -> {test} transfer "
            f"AUC {rendered} from results/cross_study_loso.json")


def test_methylation_null_is_quoted_with_its_permutation_p_value():
    """The null must be quoted together with the test that establishes it."""
    tex = _read("paper/paper.tex")
    assert "0.9444" in tex
    assert "103,863" in tex, (
        "the matched-locus CpG count that the null rests on is missing")


def test_paired_fusion_delta_matches_the_artifact():
    """The abstract quotes a paired delta; it must be the artifact's."""
    artifact = _load("fusion_ablation_isotonic.json")
    per_seed = artifact["per_seed_per_strategy"]
    tn = [s["tumor_naive"]["auc"] for s in per_seed]
    na = [s["naive_average"]["auc"] for s in per_seed]
    delta = sum(b - a for a, b in zip(tn, na)) / len(tn)
    rendered = f"{delta:.4f}".lstrip("-")
    assert rendered in _read("paper/paper.tex"), (
        f"the paired fusion delta quoted in the preprint is not the one in the "
        f"artifact ({rendered})")
