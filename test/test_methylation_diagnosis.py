"""Tests for the methylation channel diagnosis and library check.

These two scripts exist to decide whether the methylation channel carries
signal or is a restatement of library complexity, so their arithmetic is
pinned down here:

- the boundary estimator turns "this observation sat at a fragment end"
  into an unmethylated probability, and must ignore CpGs that were seen
  too rarely to have a state;
- the separation helper must not report a feature that separates in the
  *opposite* direction as a feature that does not separate;
- the matched-locus test must find a real difference when one is planted
  and must not invent one when there is none;
- the library check must recognise a size-selected library and must
  tolerate the ``chr`` prefix difference between the fragment files and
  the reference CpG bedgraphs.

No java, no network and no real data are needed.
"""
from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        name, _REPO_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


diag = _load("methylation_channel_diagnosis")
lib = _load("methylation_library_check")


MATRIX_HEADER = ("chr\tstart\tend\treadName\tFragLen\tFrag_strand\t"
                 "methy_stat\tNorm_Frag_cov\tbaseQ\tOffset_frag\t"
                 "Dist_frag_end\tkmer")


def _write_matrix(path: Path, rows) -> None:
    """rows: list of (start, dist_to_frag_end)."""
    with gzip.open(path, "wt") as f:
        f.write(MATRIX_HEADER + "\n")
        for i, (start, dist) in enumerate(rows):
            f.write(f"chr22\t{start}\t{start + 1}\tfrag{i}\t150\t+\tm\t"
                    f"0.518072\t40\t50\t{dist}\tACGT\n")


def _write_decoded(path: Path, betas) -> None:
    with gzip.open(path, "wt") as f:
        f.write("#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\t"
                "total_count_predict\tmethy_perc_obs\tmethy_count_obs\t"
                "total_count_obs\n")
        for i, b in enumerate(betas):
            f.write(f"chr22\t{1000 + i * 10}\t{1001 + i * 10}\t"
                    f"{b}\t1\t1\t{b}\t1\t1\n")


# ---------------------------------------------------------------- diagnosis
def test_read_matrix_counts_observations_and_cpgs(tmp_path):
    path = tmp_path / "s.chr22_matrix.txt.gz"
    _write_matrix(path, [(100, 10), (100, 20), (200, 0)])
    info = diag.read_matrix(str(path))
    assert info["n_observations"] == 3
    assert info["n_cpgs"] == 2
    assert info["mean_obs_per_cpg"] == pytest.approx(1.5)
    # one of the three observations sat at a fragment end
    assert info["n_observations_at_fragment_end"] == 1
    assert info["boundary_observation_rate"] == pytest.approx(1 / 3, abs=1e-5)


def test_read_matrix_reports_a_truncated_file_as_empty(tmp_path):
    path = tmp_path / "broken.chr22_matrix.txt.gz"
    path.write_bytes(b"\x1f\x8b\x08\x00 truncated")
    info = diag.read_matrix(str(path))
    assert info["n_observations"] == 0
    assert info["n_cpgs"] == 0


def test_boundary_estimator_ignores_cpgs_seen_only_once():
    """A single observation cannot tell "interior" from "at the end"."""
    per_cpg = {("chr22", 1): [0], ("chr22", 2): [0, 5], ("chr22", 3): [7, 7]}
    est = diag.boundary_estimates(per_cpg, min_obs=2)
    assert est["n_cpgs_used"] == 2
    # the CpG seen twice with one end-observation is half unmethylated
    assert est["mean_unmethylated_fraction"] == pytest.approx(0.25)


def test_boundary_estimator_reports_nothing_without_identified_cpgs():
    est = diag.boundary_estimates({("chr22", 1): [3]}, min_obs=2)
    assert est["n_cpgs_used"] == 0
    assert est["mean_unmethylated_fraction"] is None


def test_separation_is_direction_agnostic():
    """A feature that is LOWER in the positive class still separates."""
    y = np.array([0, 0, 0, 1, 1, 1])
    x = np.array([9.0, 8.0, 7.0, 3.0, 2.0, 1.0])
    sep = diag.separation(x, y)
    assert sep["auc"] < 0.5
    assert sep["auc_directionless"] == pytest.approx(1.0)
    assert sep["direction"] == "lower_feature_means_positive_class"


def test_separation_reports_none_for_a_single_class():
    y = np.array([0, 0, 0])
    x = np.array([1.0, 2.0, 3.0])
    assert diag.separation(x, y)["auc"] is None


def test_read_decoded_reports_decoder_degeneracy(tmp_path):
    path = tmp_path / "s.decoded.txt.gz"
    _write_decoded(path, [100.0, 100.0, 0.0, 55.0])
    dec = diag.read_decoded(str(path))
    assert dec["n_cpgs"] == 4
    assert dec["fraction_exactly_1"] == pytest.approx(0.5)
    assert dec["fraction_exactly_0"] == pytest.approx(0.25)
    assert dec["fraction_intermediate"] == pytest.approx(0.25)


def test_matched_locus_test_finds_a_planted_difference():
    # Eight samples per arm: with only two, the permutation null is
    # symmetric and every draw is as extreme as the observation.
    per_cpg = {}
    labels = {}
    for arm, value in (("case", 0.9), ("ctrl", 0.1)):
        for s in range(8):
            name = f"{arm}{s}"
            labels[name] = 1 if arm == "case" else 0
            per_cpg[name] = {}
            for c in range(100):
                per_cpg[name][("chr22", c)] = value
    out = diag.matched_locus_test(per_cpg, labels, n_permutations=200)
    assert out["n_cpgs_tested"] == 100
    assert out["median_difference_cancer_minus_control"] == pytest.approx(0.8)
    assert out["permutation_p_value"] < 0.01


def test_matched_locus_test_finds_nothing_when_there_is_nothing():
    per_cpg = {}
    labels = {}
    for arm in ("case", "ctrl"):
        for s in range(8):
            name = f"{arm}{s}"
            labels[name] = 1 if arm == "case" else 0
            per_cpg[name] = {}
            for c in range(100):
                per_cpg[name][("chr22", c)] = 0.5
    out = diag.matched_locus_test(per_cpg, labels, n_permutations=200)
    assert out["median_difference_cancer_minus_control"] == pytest.approx(0.0)
    assert out["permutation_p_value"] == pytest.approx(1.0)


def test_matched_locus_test_refuses_to_score_too_few_cpgs():
    per_cpg = {"case": {("chr22", 1): 0.9}, "ctrl": {("chr22", 1): 0.1}}
    out = diag.matched_locus_test(per_cpg, {"case": 1, "ctrl": 0})
    assert "note" in out


# ------------------------------------------------------------ library check
def test_normalise_chrom_ignores_the_chr_prefix():
    assert lib.normalise_chrom("chr22") == "22"
    assert lib.normalise_chrom("22") == "22"
    assert lib.normalise_chrom("chrX") == "X"


def test_load_cpg_positions_matches_across_prefixes(tmp_path):
    path = tmp_path / "cpg.bedgraph.gz"
    with gzip.open(path, "wt") as f:
        f.write("chr22\t100\t101\t1\n")
        f.write("chr22\t200\t201\t1\n")
        f.write("chr1\t300\t301\t1\n")
    assert lib.load_cpg_positions(str(path), "22") == [100, 200]
    assert lib.load_cpg_positions(str(path), "chr22") == [100, 200]


def _write_fragments(path: Path, rows) -> None:
    with gzip.open(path, "wt") as f:
        for i, (chrom, start, end) in enumerate(rows):
            f.write(f"{chrom}\t{start}\t{end}\tfrag{i}\t60\t+\n")


def test_check_file_recognises_a_size_selected_library(tmp_path):
    """Every fragment under the ceiling means the ends are truncation points."""
    path = tmp_path / "s.frag.tsv.bgz"
    rows = [("22", 1000 + i * 10, 1000 + i * 10 + 150) for i in range(100)]
    _write_fragments(path, rows)
    out = lib.check_file(str(path), "22", [])
    assert out["n_fragments"] == 100
    assert out["fraction_at_or_above_300bp"] == 0.0
    assert out["length_profile"]["150-200"] == pytest.approx(1.0)


def test_check_file_reports_a_long_tail(tmp_path):
    path = tmp_path / "s.frag.tsv.bgz"
    rows = [("22", 1000, 1200), ("22", 2000, 2600)]
    _write_fragments(path, rows)
    out = lib.check_file(str(path), "22", [])
    assert out["fraction_at_or_above_300bp"] == pytest.approx(0.5)


def test_check_file_measures_the_digestion_signature(tmp_path):
    """Fragment ends that land on a CpG beat the chance rate."""
    cpg = tmp_path / "cpg.bedgraph.gz"
    with gzip.open(cpg, "wt") as f:
        for p in range(1000, 2000, 10):
            f.write(f"chr22\t{p}\t{p + 1}\t1\n")
    path = tmp_path / "s.frag.tsv.bgz"
    rows = [("22", 1000 + i * 10, 1010 + i * 10) for i in range(50)]
    _write_fragments(path, rows)
    out = lib.check_file(str(path), "22",
                         lib.load_cpg_positions(str(cpg), "22"))
    sig = out["fragment_ends_on_a_cpg"]
    assert sig["start_rate"] == pytest.approx(1.0)
    assert sig["start_enrichment_over_chance"] > 1.0


def test_check_file_handles_an_empty_file(tmp_path):
    path = tmp_path / "empty.frag.tsv.bgz"
    with gzip.open(path, "wt") as f:
        f.write("")
    out = lib.check_file(str(path), "22", [])
    assert out["n_fragments"] == 0
