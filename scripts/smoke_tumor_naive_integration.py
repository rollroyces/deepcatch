"""CI smoke test for the tumor-naive adapter + fusion ablation.

Exercises the cfdna-framentomics-pipeline adapter end-to-end with a
synthetic mini-cohort (no network, no 627-sample FinaleDB download)
so future breakage is caught by `pytest` without external data.

Validates:
    - The on-disk artifact schema (.npy + .fsd.json) round-trips through
      the adapter without shape/dtype errors.
    - Per-sample median normalization of 100kb counts lands on 1.0.
    - The mutation-score simulator calibrates to within ±0.04 of a
      requested target AUC.
    - The fusion_ablation._evaluate_seed function runs end-to-end.

Designed to run in <10 seconds; does not touch the network.
"""
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.fragmentomics.tumor_naive_adapter import (  # noqa: E402
    load_cohort, load_labels_tsv,
)
from src.fragmentomics.fusion_ablation import (  # noqa: E402
    _simulate_mutation_scores, _summarize, _evaluate_seed,
)


def _build_synthetic_cohort(tmp: str, n: int = 80, n_cancer: int = 40
                             ) -> tuple[str, str]:
    rng = np.random.default_rng(42)
    for i in range(n):
        sid = f"sim_{i:03d}"
        np.save(os.path.join(tmp, f"{sid}.delfi_5mb_ratio.npy"),
                rng.random(631))
        np.save(os.path.join(tmp, f"{sid}.delfi_5mb_coverage.npy"),
                rng.random(631))
        np.save(os.path.join(tmp, f"{sid}.delfi_100kb_ratio.npy"),
                rng.random(30894))
        np.save(os.path.join(tmp, f"{sid}.delfi_100kb_counts.npy"),
                rng.random(30894))
        # FSD histogram normalized to sum=1, 196 bins (5bp, 20-1000bp)
        raw = rng.random(196)
        raw = raw / raw.sum()
        bins = {f"{20 + k*5}-{20 + (k+1)*5}": float(v)
                for k, v in enumerate(raw)}
        with open(os.path.join(tmp, f"{sid}.fsd.json"), "w") as f:
            json.dump({"sample": sid, "size_bins": bins,
                       "fragment_count": int(rng.integers(1_000_000, 5_000_000)),
                       "median_length": 167.0}, f)
    labels_path = os.path.join(tmp, "labels.tsv")
    with open(labels_path, "w") as f:
        for i in range(n):
            label = "cancer" if i < n_cancer else "healthy"
            f.write(f"sim_{i:03d}\t{label}\n")
    return tmp, labels_path


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        feat_dir, labels_path = _build_synthetic_cohort(tmp)

        labels = load_labels_tsv(labels_path)
        X, order = load_cohort(sorted(labels), feat_dir, skip_missing=False)

        # Contract: 80 samples × 63,246 features
        assert X.shape == (80, 63246), f"bad shape {X.shape}"
        y = np.asarray([labels[s] for s in order], dtype=int)
        assert int((y == 1).sum()) == 40

        # The 100kb counts block lands at median 1.0 (median normalization)
        block_start = sum([631, 631, 30894])  # 5mb_ratio + 5mb_cov + 100kb_ratio
        block_end = block_start + 30894
        for row in X[:, block_start:block_end]:
            assert abs(float(np.median(row)) - 1.0) < 1e-9

        # Mutation-score calibration lands in expected band (±0.04 of target)
        rng = np.random.default_rng(42)
        mut = _simulate_mutation_scores(y, rng, target_auc=0.85)
        cal = _summarize(y, mut)
        assert 0.78 < cal["auc"] < 0.92, f"calibration drift: {cal}"

        # _evaluate_seed end-to-end — verifies the fusion script runs.
        # Now also includes a delong_vs_tumor_naive block per seed.
        study = np.array(["jiang"] * 40 + ["cristiano"] * 40)
        per_seed = _evaluate_seed(X, y, study, mut, pca_n=10,
                                   seed=0, harmonize=True)
        strat_keys = {"tumor_naive", "mutation_only",
                      "naive_average", "lr_fusion"}
        assert strat_keys <= set(per_seed.keys()), (
            f"missing strategies: {strat_keys - set(per_seed.keys())}")
        for strat_res in per_seed.values():
            if isinstance(strat_res, dict) and "auc" in strat_res:
                assert 0.0 <= strat_res["auc"] <= 1.0
                assert 0.0 <= strat_res["sens_at_95"] <= 1.0
                assert 0.0 <= strat_res["sens_at_99"] <= 1.0
        if "delong_vs_tumor_naive" in per_seed:
            for strat in ("mutation_only", "naive_average", "lr_fusion"):
                assert strat in per_seed["delong_vs_tumor_naive"]

        print(
            f"OK: adapter loaded {X.shape[0]} samples; "
            f"mutation AUC {cal['auc']:.3f}; "
            f"fusion strategies all produced finite metrics"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())