"""The methylation channel's numbers must not depend on the hash seed.

What went wrong
---------------
The shared-locus ("common CpG set") analysis collected the shared loci in
a ``set`` of ``"chr:pos"`` strings and then iterated that set directly::

    vals = [m[k] for k in common if k in m]

Python randomises string hashing per process (``PYTHONHASHSEED``), so
the set iterates in a different order on every run. The order reaches
the numbers: the float sum behind ``np.mean`` is not associative, so the
last bits of ``beta_mean`` move, and the leave-one-out AUC computed from
those summaries moves with it. Measured on the same 31-sample cohort,
the common-locus AUC came back as **0.8120**, **0.8162** and **0.8291**
under three hash seeds, while every headline number (decoded-beta AUC
0.9444, the coverage baseline 0.9444, the per-CpG median difference 0.0
over 103,863 CpGs) stayed put.

This test pins the ordering rule directly, so the artifact cannot depend
on the hash seed again. It does not need java, FinaleMe or the matrices.
"""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "methylation_cohort_channel",
        _REPO_ROOT / "scripts" / "methylation_cohort_channel.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_mod = _load_module()


def _loci(n=500):
    return {f"chr1:{i * 10 + 1}": float(i % 7) / 7.0 for i in range(n)}


def test_common_values_are_returned_in_sorted_key_order():
    loci = _loci()
    values = _mod._ordered_common_values(set(loci), loci)
    expected = [loci[k] for k in sorted(loci)]
    assert values == expected


def test_common_values_do_not_depend_on_the_set_iteration_order():
    """The same loci, inserted in different orders, must give the same
    values in the same order."""
    loci = _loci()
    keys = list(loci)
    a = _mod._ordered_common_values(set(keys), loci)
    b = _mod._ordered_common_values(set(reversed(keys)), loci)
    assert a == b


def test_common_values_ignore_keys_the_map_does_not_cover():
    loci = {"chr1:1": 0.1, "chr1:2": 0.2}
    common = {"chr1:1", "chr1:2", "chr1:3"}
    assert _mod._ordered_common_values(common, loci) == [0.1, 0.2]


def test_the_script_does_not_iterate_a_set_directly_for_the_summary():
    """A source guard: the shared-locus values must come from the
    ordering helper, not from a bare `for k in common`."""
    source = (_REPO_ROOT / "scripts" / "methylation_cohort_channel.py").read_text()
    assert "for k in common" not in source, (
        "the shared-locus summary iterates the locus set directly again; "
        "use _ordered_common_values() so the numbers do not depend on "
        "PYTHONHASHSEED"
    )


@pytest.mark.parametrize("seed", ["1", "987654321"])
def test_helper_output_is_identical_across_hash_seeds(seed):
    """Run the helper in a fresh interpreter with a fixed hash seed and
    confirm the output is byte-identical to the default-seed run."""
    program = (
        "import importlib.util, json, sys\n"
        "spec = importlib.util.spec_from_file_location('m', sys.argv[1])\n"
        "m = importlib.util.module_from_spec(spec)\n"
        "spec.loader.exec_module(m)\n"
        "loci = {f'chr1:{i*10+1}': float(i % 7)/7.0 for i in range(500)}\n"
        "vals = m._ordered_common_values(set(loci), loci)\n"
        "print(json.dumps([repr(v) for v in vals]))\n"
    )
    env = {**os.environ, "PYTHONHASHSEED": seed}
    env.pop("PYTHONPATH", None)
    proc = subprocess.run(
        [sys.executable, "-c", program,
         str(_REPO_ROOT / "scripts" / "methylation_cohort_channel.py")],
        capture_output=True, text=True, env=env, timeout=120)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip(), "the helper produced no output"

    env_default = {**os.environ}
    env_default.pop("PYTHONHASHSEED", None)
    env_default.pop("PYTHONPATH", None)
    base = subprocess.run(
        [sys.executable, "-c", program,
         str(_REPO_ROOT / "scripts" / "methylation_cohort_channel.py")],
        capture_output=True, text=True, env=env_default, timeout=120)
    assert base.returncode == 0, base.stderr
    assert proc.stdout == base.stdout, (
        "the shared-locus values changed with PYTHONHASHSEED"
    )
