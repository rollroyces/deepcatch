"""
XGBoost classifier for cancer screening -- first-line triage model.

Ported from ``cancer-screening/src/ml_models/xgboost_classifier.py`` (SHA f4b13827).

Fast, interpretable baseline for cancer detection from multi-omics features.
SHAP-based explainability for every prediction.

NOTE: Requires ``pip install 'deepcatch[ultra]'`` (xgboost, shap).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import os
import pickle
import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

try:
    import xgboost as xgb
    import shap
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
    from sklearn.isotonic import IsotonicRegression
    HAS_DEPS = True
except ImportError:
    HAS_DEPS = False


class XGBoostCancerClassifier:
    """
    XGBoost-based cancer classifier for liquid biopsy screening.

    Training:
        model = XGBoostCancerClassifier()
        model.train(features_df, labels_series, config)

    Prediction:
        results = model.predict(test_features)
    """

    DEFAULT_PARAMS = {
        'objective': 'binary:logistic',
        'eval_metric': 'aucpr',
        'max_depth': 6,
        'learning_rate': 0.05,
        'subsample': 0.8,
        'colsample_bytree': 0.8,
        'min_child_weight': 3,
        'reg_alpha': 0.1,
        'reg_lambda': 1.0,
        'scale_pos_weight': 50,
        'n_estimators': 1000,
        'early_stopping_rounds': 50,
        'verbosity': 0,
        'tree_method': 'hist',
    }

    def __init__(
        self,
        output_dir: str = "models",
        seed: int = 42,
    ):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.seed = seed
        self.model = None
        self.calibrator = None
        self.shap_explainer = None
        self.feature_names = []
        self.config = {}

        logger.info("XGBoostCancerClassifier initialized")

    def train(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        config: Optional[dict] = None,
        eval_fraction: float = 0.2,
        n_folds: int = 5,
    ) -> dict:
        """Train XGBoost with cross-validation and early stopping."""
        if not HAS_DEPS:
            raise ImportError(
                "Install dependencies: pip install 'deepcatch[ultra]' "
                "(xgboost, shap, scikit-learn)"
            )

        params = {**self.DEFAULT_PARAMS, **(config or {})}
        self.config = params
        self.feature_names = list(X.columns)

        logger.info(
            f"Training XGBoost: {X.shape[0]} samples, {X.shape[1]} features, "
            f"positive rate={y.mean():.3f}"
        )

        from sklearn.model_selection import train_test_split
        X_train, X_val, y_train, y_val = train_test_split(
            X, y, test_size=eval_fraction, stratify=y, random_state=self.seed
        )

        dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=self.feature_names)
        dval = xgb.DMatrix(X_val, label=y_val, feature_names=self.feature_names)

        evals = [(dtrain, 'train'), (dval, 'eval')]
        booster = xgb.train(
            params,
            dtrain,
            num_boost_round=params['n_estimators'],
            evals=evals,
            early_stopping_rounds=params.get('early_stopping_rounds', 50),
            verbose_eval=100,
        )

        cv_results = xgb.cv(
            {k: v for k, v in params.items() if k not in ['n_estimators', 'early_stopping_rounds']},
            dtrain,
            num_boost_round=params['n_estimators'],
            nfold=n_folds,
            stratified=True,
            metrics=['auc', 'aucpr'],
            seed=self.seed,
        )

        self.model = xgb.XGBClassifier(
            n_estimators=booster.best_iteration,
            max_depth=params['max_depth'],
            learning_rate=params['learning_rate'],
            subsample=params['subsample'],
            colsample_bytree=params['colsample_bytree'],
            min_child_weight=params['min_child_weight'],
            reg_alpha=params['reg_alpha'],
            reg_lambda=params['reg_lambda'],
            scale_pos_weight=params['scale_pos_weight'],
            use_label_encoder=False,
            eval_metric='aucpr',
            random_state=self.seed,
        )
        self.model.fit(X_train, y_train)

        raw_probs = self.model.predict_proba(X_val)[:, 1]
        self.calibrator = IsotonicRegression(out_of_bounds='clip')
        self.calibrator.fit(raw_probs, y_val)

        val_probs = self.calibrator.predict(raw_probs)
        metrics = {
            'best_iteration': booster.best_iteration,
            'best_score': booster.best_score,
            'cv_auc_mean': float(cv_results['test-aucpr-mean'].iloc[-1]),
            'cv_auc_std': float(cv_results['test-aucpr-stdv'].iloc[-1]),
            'val_auc_roc': float(roc_auc_score(y_val, val_probs)),
            'val_auc_pr': float(average_precision_score(y_val, val_probs)),
            'val_n_samples': len(y_val),
            'val_n_cancer': int(y_val.sum()),
        }

        self.shap_explainer = shap.TreeExplainer(self.model)

        self.save()
        logger.info(f"Training complete. Metrics: {metrics}")

        return metrics

    def predict(
        self,
        X: pd.DataFrame,
        return_shap: bool = True,
        return_explanation: bool = True,
    ) -> dict:
        """Predict cancer probability with SHAP explanations."""
        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")

        if len(X.shape) == 1:
            X = X.to_frame().T

        raw_probs = self.model.predict_proba(X)[:, 1]

        if self.calibrator:
            probs = self.calibrator.predict(raw_probs)
        else:
            probs = raw_probs

        risk_levels = []
        for p in probs:
            if p >= 0.5:
                risk_levels.append('high')
            elif p >= 0.1:
                risk_levels.append('medium')
            else:
                risk_levels.append('low')

        results = {
            'cancer_probability': probs,
            'raw_score': raw_probs,
            'risk_level': np.array(risk_levels),
        }

        if return_shap and self.shap_explainer:
            shap_values = self.shap_explainer.shap_values(X)
            results['shap_values'] = shap_values

            if return_explanation:
                explanations = self._generate_explanations(X, shap_values)
                results['explanation'] = explanations

        return results

    def _generate_explanations(
        self,
        X: pd.DataFrame,
        shap_values: np.ndarray,
    ) -> list:
        """Generate human-readable explanations from SHAP values."""
        explanations = []

        for i in range(len(X)):
            row_shap = shap_values[i]
            row_features = X.iloc[i]
            feature_names = self.feature_names

            sorted_idx = np.argsort(np.abs(row_shap))[::-1]

            top_features = []
            for j in sorted_idx[:5]:
                fname = feature_names[j]
                shap_val = row_shap[j]
                feat_val = float(row_features.iloc[j])

                direction = "increased" if shap_val > 0 else "decreased"
                top_features.append(
                    f"  {fname}: {feat_val:.4f} (SHAP: {shap_val:+.4f}, {direction} cancer risk)"
                )

            explanations.append(
                f"Sample {i}: Cancer probability driven by:\n" + "\n".join(top_features)
            )

        return explanations

    def feature_importance(self, top_k: int = 30) -> pd.DataFrame:
        """Get top-k most important features."""
        if self.model is None:
            raise ValueError("Model not trained.")

        importance = self.model.feature_importances_
        df = pd.DataFrame({
            'feature': self.feature_names,
            'importance': importance,
        }).sort_values('importance', ascending=False)

        return df.head(top_k)

    def save(self, path: Optional[str] = None) -> str:
        """Save model, calibrator, and config to disk."""
        path = path or str(self.output_dir / "xgb_cancer_classifier.pkl")

        artifacts = {
            'model': self.model,
            'calibrator': self.calibrator,
            'config': self.config,
            'feature_names': self.feature_names,
            'seed': self.seed,
        }

        with open(path, 'wb') as f:
            pickle.dump(artifacts, f)

        logger.info(f"Model saved to {path}")
        return path

    def load(self, path: str) -> None:
        """Load model artifacts from disk."""
        with open(path, 'rb') as f:
            artifacts = pickle.load(f)

        self.model = artifacts['model']
        self.calibrator = artifacts['calibrator']
        self.config = artifacts['config']
        self.feature_names = artifacts['feature_names']
        self.seed = artifacts['seed']

        self.shap_explainer = shap.TreeExplainer(self.model)
        logger.info(f"Model loaded from {path}")


def simulate_training_data(
    n_samples: int = 5000,
    seed: int = 42,
) -> tuple:
    """Generate simulated multi-omics features for testing."""
    np.random.seed(seed)

    n_cancer = int(n_samples * 0.02)
    n_healthy = n_samples - n_cancer

    mutation_genes = [
        'TP53', 'KRAS', 'EGFR', 'PIK3CA', 'BRAF', 'PTEN',
        'APC', 'CDKN2A', 'RB1', 'ATM'
    ]

    features = {}

    for gene in mutation_genes:
        healthy_vaf = np.random.exponential(1e-6, n_healthy)
        cancer_vaf = np.random.exponential(5e-5, n_cancer)
        features[f'vaf_{gene}'] = np.concatenate([healthy_vaf, cancer_vaf])

    methylation_genes = ['RASSF1A', 'SEPT9', 'SHOX2', 'BRCA1', 'GSTP1']
    for gene in methylation_genes:
        healthy_beta = np.random.beta(2, 8, n_healthy)
        cancer_beta = np.random.beta(8, 2, n_cancer)
        features[f'beta_{gene}'] = np.concatenate([healthy_beta, cancer_beta])

    for gene in ['MYC', 'ERBB2', 'EGFR', 'PTEN']:
        healthy_cnv = np.random.normal(0, 0.1, n_healthy)
        cancer_cnv = np.random.normal(0.5 if gene in ['MYC', 'ERBB2'] else -0.5, 0.3, n_cancer)
        features[f'cnv_log2_{gene}'] = np.concatenate([healthy_cnv, cancer_cnv])

    proteins = ['CA125', 'CEA', 'CA19-9', 'PSA', 'AFP']
    for protein in proteins:
        healthy_prot = np.random.lognormal(0, 0.5, n_healthy)
        cancer_prot = np.random.lognormal(1.5, 0.8, n_cancer)
        features[f'protein_{protein}'] = np.concatenate([healthy_prot, cancer_prot])

    features['age'] = np.concatenate([
        np.random.randint(30, 75, n_healthy),
        np.random.randint(40, 80, n_cancer),
    ])
    features['bmi'] = np.concatenate([
        np.random.normal(25, 4, n_healthy),
        np.random.normal(27, 5, n_cancer),
    ])
    features['pack_years'] = np.concatenate([
        np.random.exponential(2, n_healthy),
        np.random.exponential(15, n_cancer),
    ])

    df = pd.DataFrame(features)
    labels = np.concatenate([np.zeros(n_healthy), np.ones(n_cancer)])

    shuffle_idx = np.random.permutation(n_samples)
    return (
        df.iloc[shuffle_idx].reset_index(drop=True),
        pd.Series(labels[shuffle_idx]),
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    print("Generating simulated training data...")
    X, y = simulate_training_data(n_samples=5000, seed=42)
    print(f"Data shape: {X.shape}, Positive rate: {y.mean():.3f}")

    print("Training XGBoost classifier...")
    clf = XGBoostCancerClassifier(output_dir="models")
    metrics = clf.train(X, y, n_folds=3, eval_fraction=0.2)

    print(f"\nTraining complete!")
    print(f"  CV AUC-PR: {metrics['cv_auc_mean']:.4f} +/- {metrics['cv_auc_std']:.4f}")
    print(f"  Validation AUC-ROC: {metrics['val_auc_roc']:.4f}")
    print(f"  Validation AUC-PR: {metrics['val_auc_pr']:.4f}")

    print("\nTop 10 important features:")
    importance = clf.feature_importance(top_k=10)
    print(importance.to_string(index=False))
