"""Real implementation for run_all.py.

Original location: src/agent6-ensemble/run_all.py.py (committed as broken symlink,
now restored as a real module).

This is a minimal, stdlib + numpy + scipy + sklearn implementation that
provides the surface advertised by the README. Domain-specific knobs (real
trained weights, GPU kernels, hyperparameter search) are intentionally
omitted — see docs/ and METHYLATION_PROJECT.md for the full design.

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent6-ensemble/` directory that was never committed to git.
"""

from __future__ import annotations

"""End-to-end driver for the ensemble pipeline.

Runs all available ensembles on a synthetic cohort and reports operating
points. Intended as a smoke test — see docs/V3_DESIGN.md for the full
methylation β-value extension plan.
"""
from .ensemble_core import EnsembleMember, StackedEnsemble, TwoStageScreener
from .decision_threshold import threshold_at_specificity, sweep_specificities
from .simulate_performance import simulate_cohort, measure_at_specificity


def _dummy_member(seed: int) -> EnsembleMember:
    import numpy as np
    rng = np.random.default_rng(seed)

    def _predict(X):
        return X @ rng.normal(0, 1, size=(X.shape[1], 1)) + rng.normal(0, 0.1)
    return EnsembleMember(name=f"member_{seed}", predict_fn=_predict)


def run_all(n_pos: int = 200, n_neg: int = 800,
            target_specificity: float = 0.99,
            n_iter: int = 50) -> dict:
    X, y = simulate_cohort(n_pos=n_pos, n_neg=n_neg)
    members = [_dummy_member(s) for s in range(3)]
    ens = StackedEnsemble(members)
    ens.fit_global_weights(X[: n_pos + n_neg // 2], y[: n_pos + n_neg // 2], n_iter=n_iter)
    screener = TwoStageScreener(ens, target_specificity=target_specificity)
    screener.fit_threshold(X[n_pos + n_neg // 2:], y_neg := X[y == 0])
    ops = sweep_specificities(ens.predict(X[y == 1]), ens.predict(X[y == 0]))
    return {
        "operating_points": ops,
        "threshold": screener._threshold,
        "sensitivity_at_target": float(
            (ens.predict(X[y == 1]).flatten() >= screener._threshold).mean()),
    }


if __name__ == "__main__":  # pragma: no cover
    import json
    out = run_all()
    print(json.dumps({"threshold": out["threshold"],
                      "n_op_points": len(out["operating_points"])},
                     indent=2))


__all__ = ["run_all"]
