"""Tests for the FinaleMe pipeline helpers.

Two things are tested here, both without any external dependency:

- the fragment-file conversion FinaleMe needs (5-column FinaleDB schema
  with bare chromosome names -> 6-column bedpe with a "chr" prefix);
- the decoder-input builder's handling of a missing prior ("NaN", upper
  case -- Java's ``Double.parseDouble`` rejects "nan").

The bigwig-backed part of the decoder-input builder is skipped where
pyBigWig is not installed; the CI job that runs this file does not
install it.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_REPO_ROOT / "scripts"))


def test_frag_to_bedpe_conversion(tmp_path):
    """FinaleDB's 5-column rows become the 6-column bedpe schema."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "fm_bam", _REPO_ROOT / "scripts" / "finaleme_bam_from_frag.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    src = tmp_path / "S1.frag.tsv.bgz"
    with gzip.open(src, "wt") as f:
        f.write("1\t100\t250\t30\t+\n")
        f.write("2\t300\t400\t0\t-\n")
        f.write("chr3\t10\t20\t60\t+\n")

    class _FakeFa:
        index = {"chr1": (0, 1000, 50), "chr2": (0, 1000, 50),
                 "chr3": (0, 1000, 50)}

        def fetch(self, chrom, start, end):
            return "ACGT" * ((end - start) // 4)

    sam = tmp_path / "out.sam"
    n = mod.write_sam(src, _FakeFa(), sam)
    assert n == 3
    lines = [l for l in sam.read_text().splitlines() if not l.startswith("@")]
    assert lines, "no alignment lines written"
    for line in lines:
        f = line.split("\t")
        assert len(f) == 11, f"expected 11 SAM fields, got {len(f)}"
        assert f[2].startswith("chr"), "chromosome names must carry the chr prefix"
        # SAM order: QNAME FLAG RNAME POS MAPQ CIGAR ...
        assert f[5].endswith("M"), "CIGAR must be a length"
        assert f[5] != f[4], (
            "CIGAR must be the fragment length, not the MAPQ (a BAM built "
            "that way is rejected by htsjdk")


def test_decoder_input_writes_nan_for_missing_priors(tmp_path):
    """A CpG with no prior must be "NaN" (upper case), not "nan"."""
    pyBigWig = pytest.importorskip("pyBigWig")

    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "fm_dec", _REPO_ROOT / "scripts" / "finaleme_decoder_input.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    matrix = tmp_path / "m.txt.gz"
    header = ["chr", "start", "end", "readName", "FragLen", "Frag_strand",
              "methy_stat", "Norm_Frag_cov", "baseQ", "Offset_frag",
              "Dist_frag_end", "kmer"]
    with gzip.open(matrix, "wt") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(["chr1", "100", "101", "r1", "150", "+", "m",
                            "0.65", "40", "10", "10", "AGTC"]) + "\n")
        f.write("\t".join(["chr2", "200", "201", "r2", "150", "+", "m",
                            "0.65", "40", "10", "10", "AGTC"]) + "\n")

    class _FakeBw:
        def __init__(self, *a, **k):
            pass

        def values(self, chrom, start, end):
            if chrom == "chr1":
                return [12.5]
            return [float("nan")]

        def close(self):
            pass

    mod.pyBigWig = None
    import sys as _sys
    _sys.modules["pyBigWig"] = type("m", (), {"open": staticmethod(_FakeBw)})

    out = tmp_path / "d.txt.gz"
    stats = mod.build(matrix, Path("/does/not/matter.bw"), out)
    assert stats["rows"] == 2
    assert stats["with_prior"] == 1
    rows = gzip.open(out, "rt").read().splitlines()
    assert rows[0].endswith("methyPrior")
    assert rows[1].endswith("12.5000")
    assert rows[2].endswith("NaN")
    assert not rows[2].endswith("nan")
