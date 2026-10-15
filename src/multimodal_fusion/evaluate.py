"""Real implementation for evaluate.py.

Original location: src/agent2-multimodal-fusion/evaluate.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent2-multimodal-fusion/` directory that was never committed to git.
"""

from __future__ import annotations

"""Compute AUC + sens@fixed-spec for the multimodal fusion classifier."""
from dataclasses import dataclass

import numpy as np
from sklearn.metrics import roc_auc_score

from .train import collect_features
from ..ensemble.decision_threshold import threshold_at_specificity


@dataclass
class FusionEval:
    auc: float
    threshold: float
    sensitivity_at_target_spec: float


def evaluate(target_specificity: float = 0.99, seed: int = 1) -> FusionEval:
    from .train import train_classifier
    clf = train_classifier(seed=seed)
    X_val, y_val = collect_features(seed=seed + 100)
    scores = clf.predict_proba(X_val)[:, 1]
    auc = float(roc_auc_score(y_val, scores))
    t = threshold_at_specificity(scores[y_val == 0], target_specificity)
    sens = float((scores[y_val == 1] >= t).mean())
    return FusionEval(auc=auc, threshold=t, sensitivity_at_target_spec=sens)


__all__ = ["FusionEval", "evaluate"]
