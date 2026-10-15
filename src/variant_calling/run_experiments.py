"""Real implementation for run_experiments.py.

Original location: src/agent1-variant-calling/run_experiments.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent1-variant-calling/` directory that was never committed to git.
"""

from __future__ import annotations

"""Orchestrate a small variant-calling benchmark on synthetic reads."""
import numpy as np

from .bayesian_caller import call_variant
from .baseline_callers import simple_threshold_caller
from .evaluate import evaluate_calls
from .simulate import simulate_reads


def run(n_variants: int = 200, n_sites: int = 1000, seed: int = 0) -> dict:
    reads, true_pos = simulate_reads(n_variants=n_variants,
                                      n_sites=n_sites, seed=seed)
    called = np.zeros(n_sites, dtype=bool)
    for site_idx, (ref_n, alt_n) in enumerate(reads):
        call = call_variant("chr1", site_idx, "A", "T", ref_n, alt_n)
        if call.p_variant > 0.5:
            called[site_idx] = True
    eval_ = evaluate_calls(true_pos, called)
    return {
        "n_sites": n_sites,
        "n_variants_truth": int(true_pos.sum()),
        "n_called": int(called.sum()),
        "precision": eval_.precision,
        "recall": eval_.recall,
        "f1": eval_.f1,
    }


if __name__ == "__main__":  # pragma: no cover
    import json
    print(json.dumps(run(), indent=2))


__all__ = ["run"]
