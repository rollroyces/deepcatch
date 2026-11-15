"""Tests for the methylation cohort-channel scoring.

The functions under test are the ones that decide whether the channel
is real or an artefact, so they are worth pinning down:

- ``summarise`` computes the per-sample beta summary and can cap the
  number of CpGs, which is what removes the coverage confound;
- the scoring reports the beta summary *and* the coverage baseline, so a
  reader can see whether the channel adds anything over library
  complexity.

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

_spec = importlib.util.spec_from_file_location(
    "methylation_cohort_channel",
    _REPO_ROOT / "scripts" / "methylation_cohort_channel.py")
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)

summarise = _mod.summarise


def _write_decoded(path: Path, betas, chromosome: str = "chr22") -> None:
    with gzip.open(path, "wt") as f:
        f.write("#chr\tstart\tend\tmethy_perc_predict\tmethy_count_predict\t"
                "total_count_predict\tmethy_perc_obs\tmethy_count_obs\t"
                "total_count_obs\n")
        for i, b in enumerate(betas):
            f.write(f"{chromosome}\t{1000 + i * 10}\t{1001 + i * 10}\t"
                    f"{b}\t1\t1\t{b}\t1\t1\n")


def test_summarise_reads_beta_percentages(tmp_path):
    path = tmp_path / "s.decoded.txt.gz"
    _write_decoded(path, [100.0, 50.0, 0.0, 100.0])
    s = summarise(str(path))
    assert s["n_cpgs"] == 4
    assert s["beta_mean"] == pytest.approx(0.625)
    # median of [0.0, 0.5, 1.0, 1.0]
    assert s["beta_median"] == pytest.approx(0.75)
    assert s["fraction_methylated"] == pytest.approx(0.5)
    assert s["fraction_fully_methylated"] == pytest.approx(0.5)


def test_summarise_caps_the_cpg_count_uniformly(tmp_path):
    """Depth matching keeps every k-th CpG, not a prefix."""
    path = tmp_path / "s.decoded.txt.gz"
    _write_decoded(path, [0.0] * 10 + [100.0] * 10)
    s = summarise(str(path), max_cpgs=5)
    assert s["n_cpgs"] == 5
    assert s["n_cpgs_before_match"] == 20
    # A uniform stride over a file that is half zero and half one must
    # keep a mix; a prefix would keep only zeros.
    assert 0.0 < s["beta_mean"] < 1.0


def test_summarise_reports_a_truncated_file_as_empty(tmp_path):
    """A file still being written must not raise or be counted."""
    path = tmp_path / "broken.decoded.txt.gz"
    path.write_bytes(b"\x1f\x8b\x08\x00 truncated")
    s = summarise(str(path))
    assert s["n_cpgs"] == 0


def test_summarise_handles_an_empty_file(tmp_path):
    path = tmp_path / "empty.decoded.txt.gz"
    with gzip.open(path, "wt") as f:
        f.write("#chr\tstart\tend\n")
    s = summarise(str(path))
    assert s["n_cpgs"] == 0
