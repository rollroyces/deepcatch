"""End-to-end AUC gate for the DeepCatch tumor-naive adapter.

Mirrors `scripts/auc_reproducibility_gate.py` in the pipeline repo,
but on the *DeepCatch* side: it builds the same synthetic cohort and
runs the DeepCatch adapter (load_cohort + 5-seed LR-on-PCA) end-to-end.

If a regression in the adapter breaks the signal-pickup path (e.g.
median-normalization flipped, scaling bug, the wrong channel index
in the load loop), this gate catches it.

Runs in <30 seconds, no network required.
"""
import json
import os
import sys
import tempfile

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, PROJECT_ROOT)

from src.fragmentomics.tumor_naive_adapter import load_cohort, load_labels_tsv  # noqa

RNG = np.random.default_rng(20240820)
FLOOR_AUC = 0.80


def _make_synthetic_cohort(features_dir: str, n_cancer: int = 40,
                             n_healthy: int = 40) -> list[str]:
    """Same synthetic cohort as the pipeline's AUC gate."""
    samples = []
    for i in range(n_cancer):
        s = f"canc_{i:03d}"
        np.save(os.path.join(features_dir, f"{s}.delfi_5mb_ratio.npy"),
                RNG.random(631))
        np.save(os.path.join(features_dir, f"{s}.delfi_5mb_coverage.npy"),
                RNG.random(631) + 0.5)
        samples.append(s)
    for i in range(n_healthy):
        s = f"hlth_{i:03d}"
        np.save(os.path.join(features_dir, f"{s}.delfi_5mb_ratio.npy"),
                RNG.random(631))
        np.save(os.path.join(features_dir, f"{s}.delfi_5mb_coverage.npy"),
                RNG.random(631) + 0.5)
        samples.append(s)
    noise_100 = RNG.random((30894,))
    signal_shift = np.zeros(30894)
    signal_shift[:50] = 0.20
    for i in range(n_cancer):
        s = f"canc_{i:03d}"
        np.save(os.path.join(features_dir, f"{s}.delfi_100kb_ratio.npy"),
                noise_100 + signal_shift)
        np.save(os.path.join(features_dir, f"{s}.delfi_100kb_counts.npy"),
                RNG.random(30894) * 1000 + 100)
    for i in range(n_healthy):
        s = f"hlth_{i:03d}"
        np.save(os.path.join(features_dir, f"{s}.delfi_100kb_ratio.npy"),
                noise_100)
        np.save(os.path.join(features_dir, f"{s}.delfi_100kb_counts.npy"),
                RNG.random(30894) * 1000 + 100)
    for s in samples:
        raw = RNG.random(196)
        raw = raw / raw.sum()
        bins = {f"{20+k*5}-{20+(k+1)*5}": float(v) for k, v in enumerate(raw)}
        with open(os.path.join(features_dir, f"{s}.fsd.json"), "w") as f:
            json.dump({"sample": s, "size_bins": bins,
                       "fragment_count": int(RNG.integers(1_000_000, 5_000_000)),
                       "median_length": 167.0}, f)
    labels_path = os.path.join(features_dir, "labels.tsv")
    with open(labels_path, "w") as f:
        for s in samples:
            label = "cancer" if s.startswith("canc") else "healthy"
            f.write(f"{s}\t{label}\n")
    return samples


def main() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        _make_synthetic_cohort(tmp)
        labels = load_labels_tsv(os.path.join(tmp, "labels.tsv"))
        X, order = load_cohort(sorted(labels), tmp, skip_missing=False)
        y = np.asarray([labels[s] for s in order], dtype=int)
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
        yt, ys = [], []
        for tr, te in cv.split(X, y):
            sc = StandardScaler().fit(X[tr])
            Xtr = sc.transform(X[tr]); Xte = sc.transform(X[te])
            max_pca = min(Xtr.shape[0], Xtr.shape[1])
            pca = PCA(n_components=min(30, max_pca), random_state=0).fit(Xtr)
            m = LogisticRegression(max_iter=20000, tol=1e-8,
                                       random_state=0).fit(
                pca.transform(Xtr), y[tr])
            ys.extend(m.predict_proba(pca.transform(Xte))[:, 1].tolist())
            yt.extend(y[te].tolist())
        auc = float(roc_auc_score(yt, ys))
        print(f"[adapter_auc_gate] synthetic cohort: {X.shape[0]} samples, "
              f"AUC={auc:.4f}")
        if auc < FLOOR_AUC:
            print(f"[adapter_auc_gate] FAIL: AUC {auc:.4f} < floor {FLOOR_AUC:.4f}.")
            return 1
        print(f"[adapter_auc_gate] PASS: AUC {auc:.4f} >= floor {FLOOR_AUC:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())