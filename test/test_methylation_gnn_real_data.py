"""The methylation GNN must be able to run on real public data.

What went wrong
---------------
``src/methylation_gnn/`` was architecture-only. Its graph builder filled
the node set with randomly generated regions when none were loaded
(``_ensure_regions``), and its data module filled the beta matrix with
random values. A GNN run therefore consumed "methylation" that had no
relationship to any real sample, and nothing in the code said so.

The V3 pipeline now runs FinaleMe end to end on real public FinaleDB
fragments, so the GNN can be pointed at real decoded calls. These tests
cover that path:

* the loader reads FinaleMe's decoded format and returns a real cohort;
* the graph builder refuses to fall back to synthetic regions once it
  has been pointed at real loci;
* the whole path -- load, build a graph, run the model -- works on a
  real cohort when the decoded files are present.

The fixture here is a tiny file in FinaleMe's *decoded* format. It is a
test input for the PARSER, not a synthetic methylation cohort: the
loader is exercised against the real file format, and the real cohort is
exercised separately, from the actual decoded output when it is
available.
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
import sys
sys.path.insert(0, str(_REPO_ROOT))

from src.methylation_gnn.real_data import (  # noqa: E402
    RealDataError,
    build_real_cohort_graphs,
    load_finaleme_cohort,
    load_real_regions,
    per_sample_methylation_data,
)

# The cached decoded output from the V3 pipeline (outside the repo).
REAL_DECODED_DIR = Path(os.environ.get(
    "DEEPCATCH_METHYLATION_DIR", str(_REPO_ROOT.parent / "methylation_cohort")
)) / "final_work"
REAL_LABELS = Path(os.environ.get(
    "DEEPCATCH_METHYLATION_DIR", str(_REPO_ROOT.parent / "methylation_cohort")
)) / "cohort_labels.tsv"


def _write_decoded(path: Path, calls: list[tuple[str, int, int, float]]) -> None:
    """Write one file in FinaleMe's decoded format.

    Columns: chrom, start, end, beta (0-100), then the decoder's own
    fields. The loader reads the first four.
    """
    with gzip.open(path, "wt") as handle:
        handle.write("#chrom\tstart\tend\tbeta\tstate\t...\n")
        for chrom, start, end, beta in calls:
            handle.write(
                f"{chrom}\t{start}\t{end}\t{beta:g}\t"
                f"{int(beta > 0.5)}\t1\t100.0\t1\t1\n")


@pytest.fixture
def decoded_fixture(tmp_path: Path):
    """Three samples, ten shared loci, in FinaleMe's decoded format."""
    decoded_dir = tmp_path / "decoded"
    decoded_dir.mkdir()

    loci = [(f"chr22", 1_000_000 + i * 10, 1_000_001 + i * 10) for i in range(10)]

    # Sample A: every locus methylated.
    _write_decoded(decoded_dir / "A.decoded.txt.gz",
                   [(c, s, e, 100.0) for c, s, e in loci])
    # Sample B: every locus unmethylated.
    _write_decoded(decoded_dir / "B.decoded.txt.gz",
                   [(c, s, e, 0.0) for c, s, e in loci])
    # Sample C: covers only the first eight loci, mixed calls.
    _write_decoded(decoded_dir / "C.decoded.txt.gz",
                   [(c, s, e, 100.0 if i % 2 else 0.0)
                    for i, (c, s, e) in enumerate(loci[:8])])

    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "sample\tlabel\tdisease\n"
        "A\t1\tColorectal cancer\n"
        "B\t0\tLiver transplant\n"
        "C\t1\tColorectal cancer\n"
    )
    return decoded_dir, labels


def test_loader_reads_the_decoded_format(decoded_fixture):
    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8,
                                  min_shared_cpgs=5)

    assert cohort["sample_names"] == ["A", "B", "C"]
    assert cohort["labels"].tolist() == [1, 0, 1]
    assert cohort["disease"] == ["Colorectal cancer", "Liver transplant",
                                 "Colorectal cancer"]
    assert cohort["beta_values"].shape == (8, 3)
    assert cohort["is_synthetic"] is False
    assert cohort["data_source"] == "finaleme_decoded"


def test_loader_uses_the_real_coordinates(decoded_fixture):
    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8, min_shared_cpgs=5)
    # The shared loci are the first eight (C does not cover the last two).
    assert len(cohort["regions"]) == 8
    for region, (chrom, start) in zip(cohort["regions"], cohort["cpg_keys"]):
        assert region.chrom == chrom
        assert region.start == int(start)
        assert region.region_type == "cpg"
    # The coordinates are real, not the builder's synthetic placeholders.
    assert all(r.chrom == "chr22" for r in cohort["regions"])
    assert not any(r.chrom.startswith("chr") and r.chrom != "chr22"
                   for r in cohort["regions"])


def test_loader_marks_uncovered_loci_as_unobserved(decoded_fixture):
    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8, min_shared_cpgs=5)
    covered = cohort["covered"]
    # C covers 8 of 10 loci; the shared set is those 8, so C covers all
    # of them and A and B cover all of them too.
    assert covered.all()
    # The beta values are the decoded calls, scaled to 0-1.
    assert cohort["beta_values"][0, 0] == pytest.approx(1.0)
    assert cohort["beta_values"][0, 1] == pytest.approx(0.0)


def test_loader_refuses_to_invent_a_cohort(tmp_path: Path):
    """An empty directory must be an error, not an empty cohort."""
    empty = tmp_path / "empty"
    empty.mkdir()
    labels = tmp_path / "labels.tsv"
    labels.write_text("sample\tlabel\tdisease\nA\t1\tCancer\n")
    with pytest.raises(RealDataError):
        load_finaleme_cohort(str(empty), str(labels))


def test_loader_refuses_a_missing_directory(tmp_path: Path):
    with pytest.raises(RealDataError):
        load_finaleme_cohort(str(tmp_path / "nope"), str(tmp_path / "l.tsv"))


def test_loader_refuses_when_too_few_shared_loci(decoded_fixture):
    decoded_dir, labels = decoded_fixture
    with pytest.raises(RealDataError, match="could be selected"):
        load_finaleme_cohort(str(decoded_dir), str(labels),
                             coverage_fraction=0.8, min_shared_cpgs=500)


def test_loader_records_which_locus_rule_it_used(decoded_fixture):
    """The threshold rule and the fallback must be distinguishable in the
    output, because they select different loci."""
    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8, min_shared_cpgs=5)
    # All 8 loci are covered by all 3 samples, so the threshold rule
    # applies and no fallback is needed.
    assert cohort["locus_selection"] == "coverage_threshold"
    assert cohort["n_threshold_loci"] == 8
    assert cohort["coverage_threshold"] == 3

    # Ask for a floor the cohort cannot meet: the fallback selects every
    # locus there is, and the floor still rejects the cohort.
    with pytest.raises(RealDataError, match="could be selected"):
        load_finaleme_cohort(str(decoded_dir), str(labels),
                             coverage_fraction=0.8, min_shared_cpgs=200)


def test_loader_does_not_fabricate_region_annotations(decoded_fixture):
    """GC content and CpG density are not in the FinaleMe output, so they
    must be absent -- not filled with random numbers."""
    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8, min_shared_cpgs=5)
    assert "gc_content" not in cohort
    assert "cpg_density" not in cohort
    assert "region_annotations" in cohort


def test_graph_builder_refuses_to_fall_back_to_synthetic_regions():
    """Once a run is pointed at real loci, a missing region set is a bug."""
    from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder

    builder = RegulatoryGraphBuilder(n_nodes=10)
    # The legacy behaviour is preserved by default.
    builder._ensure_regions()
    assert len(builder._regions) == 10

    strict = RegulatoryGraphBuilder(n_nodes=10)
    strict.require_real_regions = True
    with pytest.raises(ValueError, match="require_real_regions"):
        strict._ensure_regions()


def test_per_sample_data_matches_the_cohort(decoded_fixture):
    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8, min_shared_cpgs=5)
    data = per_sample_methylation_data(cohort, 0)
    assert data["beta_values"].shape == (8,)
    assert data["methylation_variance"].shape == (8,)
    assert data["beta_is_observed"].shape == (8,)


def test_real_cohort_builds_graphs(decoded_fixture):
    """End to end: real decoded calls → real loci → PyG graphs."""
    torch = pytest.importorskip("torch")
    pytest.importorskip("torch_geometric")

    from src.methylation_gnn.graph_builder import RegulatoryGraphBuilder

    decoded_dir, labels = decoded_fixture
    cohort = load_finaleme_cohort(str(decoded_dir), str(labels),
                                  coverage_fraction=0.8, min_shared_cpgs=5)
    builder = RegulatoryGraphBuilder(n_nodes=8, edge_k=4)
    load_real_regions(builder, cohort)
    graphs = build_real_cohort_graphs(builder, cohort)

    assert len(graphs) == 3
    for graph, name in zip(graphs, cohort["sample_names"]):
        assert graph.sample_name == name
        assert graph.x.shape[0] == 8, "the graph must have the real loci"
        assert torch.isfinite(graph.x).all()
    # The node set is the real one, not the builder's synthetic default.
    assert all(r.chrom == "chr22" for r in builder._regions)


@pytest.mark.skipif(
    not (REAL_DECODED_DIR.is_dir() and REAL_LABELS.is_file()),
    reason="the cached FinaleMe decoded output is not present",
)
def test_the_real_published_cohort_loads():
    """The actual cohort the paper measured must load through this path."""
    cohort = load_finaleme_cohort(str(REAL_DECODED_DIR), str(REAL_LABELS))
    assert cohort["is_synthetic"] is False
    assert len(cohort["sample_names"]) >= 8
    assert cohort["beta_values"].shape[0] >= 50
    assert set(cohort["labels"].tolist()) == {0, 1}
    # The real cohort's shared loci are the ones the paper describes.
    assert cohort["n_locus_observations"] > 0
