"""Real implementation for train.py.

Original location: src/agent2-multimodal-fusion/train.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent2-multimodal-fusion/` directory that was never committed to git.
"""

from __future__ import annotations

"""Minimal training loop for a multi-modal fusion classifier.

Uses scikit-learn's ``LogisticRegression`` as a stand-in for the production
Transformer fusion model. Designed to be drop-in replaceable once the
trained checkpoint is committed.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression

from ..synthetic_data.pipeline import SyntheticCohortPipeline
from ..synthetic_data.data_generator import MultiModalDataGenerator


def collect_features(n_cancer: int = 100, n_healthy: int = 100,
                     seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    pipe = SyntheticCohortPipeline(
        n_cancer=n_cancer, n_healthy=n_healthy,
        generator=MultiModalDataGenerator(seed=seed),
    )
    rows, labels = [], []
    for x, y in pipe:
        rows.append(x); labels.append(y)
    return np.stack(rows), np.array(labels)


def train_classifier(n_cancer: int = 100, n_healthy: int = 100,
                     seed: int = 0, **lr_kwargs) -> LogisticRegression:
    X, y = collect_features(n_cancer=n_cancer, n_healthy=n_healthy, seed=seed)
    lr_kwargs.setdefault("max_iter", 200)
    clf = LogisticRegression(**lr_kwargs)
    clf.fit(X, y)
    return clf


__all__ = ["collect_features", "train_classifier"]
