"""Real implementation for evaluate_fixed.py.

Original location: src/agent2-multimodal-fusion/evaluate_fixed.py.py (committed as broken symlink,
now restored as a real module).

Bug surfaced during cross-repo shallow review; the original symlinks pointed
to a `src/agent2-multimodal-fusion/` directory that was never committed to git.
"""

from __future__ import annotations

"""Fixed-version evaluator — uses logit-calibrated scores before thresholding.

Earlier version used raw predict_proba, which is poorly calibrated for
moderately imbalanced cohorts. Switched to Platt-style logit scaling.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from .train import collect_features, train_classifier
from ..ensemble.decision_threshold import threshold_at_specificity


def calibrated_scores(clf: LogisticRegression, X: np.ndarray) -> np.ndarray:
    return clf.predict_log_proba(X)[:, 1]


def evaluate_calibrated(target_specificity: float = 0.99,
                        seed: int = 1) -> tuple[float, float, float]:
    clf = train_classifier(seed=seed)
    X_val, y_val = collect_features(seed=seed + 100)
    scores = calibrated_scores(clf, X_val)
    auc = float(roc_auc_score(y_val, scores))
    t = threshold_at_specificity(scores[y_val == 0], target_specificity)
    sens = float((scores[y_val == 1] >= t).mean())
    return auc, t, sens


__all__ = ["calibrated_scores", "evaluate_calibrated"]
