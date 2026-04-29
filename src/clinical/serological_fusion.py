#!/usr/bin/env python3
"""
Clinical Serological Fusion Module

Integrates traditional serum biomarkers (pepsinogen, gastrin-17, H. pylori)
with cfDNA-based DeepCatch predictions via a learnable shallow fusion branch.

Initially designed for gastric cancer (STAD), extensible to other cancer types
with their respective serum markers.

Reference: THEMIS gastric cancer framework (Bie et al. 2023)
"""
import numpy as np
from typing import Dict, List, Optional, Tuple
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


class SerologicalFusion:
    """
    Shallow fusion branch combining cfDNA predictions with serum biomarkers.
    
    Uses logistic regression to learn optimal combination weights between
    cfDNA-only scores and clinical serum markers.
    
    Markers:
    - PG I/II ratio (pepsinogen): gastric atrophy indicator
    - Gastrin-17 (G-17): gastric acid secretion marker
    - H. pylori IgG serology: infection status
    """
    
    def __init__(self, random_state: int = 42):
        self.scaler = StandardScaler()
        self.model = LogisticRegression(
            C=0.1, class_weight='balanced', random_state=random_state
        )
        self._fitted = False
        self.feature_names = ['cfdna_score', 'pg_ratio', 'gastrin_17', 'hpylori_igg']
    
    def fit(
        self,
        cfdna_scores: np.ndarray,
        serum_markers: Dict[str, np.ndarray],
        labels: np.ndarray
    ) -> 'SerologicalFusion':
        """
        Fit the serological fusion model.
        
        Parameters
        ----------
        cfdna_scores : (n,) array
            DeepCatch cfDNA-only prediction scores.
        serum_markers : dict
            Keys: 'pg_ratio', 'gastrin_17', 'hpylori_igg'.
            Values: (n,) arrays of serum marker values.
        labels : (n,) array
            Binary labels (1=cancer, 0=healthy).
        """
        features = np.column_stack([
            cfdna_scores,
            serum_markers.get('pg_ratio', np.zeros_like(cfdna_scores)),
            serum_markers.get('gastrin_17', np.zeros_like(cfdna_scores)),
            serum_markers.get('hpylori_igg', np.zeros_like(cfdna_scores)),
        ])
        
        features = self.scaler.fit_transform(features)
        self.model.fit(features, labels)
        self._fitted = True
        return self
    
    def predict_proba(
        self,
        cfdna_scores: np.ndarray,
        serum_markers: Dict[str, np.ndarray]
    ) -> np.ndarray:
        """Return cancer probability after serological fusion."""
        if not self._fitted:
            raise RuntimeError("Model not fitted. Call fit() first.")
        
        features = np.column_stack([
            cfdna_scores,
            serum_markers.get('pg_ratio', np.zeros_like(cfdna_scores)),
            serum_markers.get('gastrin_17', np.zeros_like(cfdna_scores)),
            serum_markers.get('hpylori_igg', np.zeros_like(cfdna_scores)),
        ])
        
        features = self.scaler.transform(features)
        return self.model.predict_proba(features)[:, 1]
    
    def get_weights(self) -> Dict[str, float]:
        """Return learned feature weights."""
        if not self._fitted:
            return {}
        return dict(zip(self.feature_names, self.model.coef_[0]))


class IntegrativeScoringSystem:
    """
    Integrative Scoring System combining cfDNA + serological markers.
    
    Toggleable: when serum data is available, activates serological fusion.
    When not available, falls back to cfDNA-only scoring.
    """
    
    def __init__(self, enable_serological: bool = False):
        self.enable_serological = enable_serological
        self.sero_fusion = SerologicalFusion() if enable_serological else None
    
    def score(
        self,
        cfdna_scores: np.ndarray,
        serum_markers: Optional[Dict[str, np.ndarray]] = None
    ) -> np.ndarray:
        """
        Compute integrated cancer risk scores.
        
        If serological fusion is enabled AND serum data provided, combines
        both. Otherwise returns cfDNA-only scores.
        """
        if self.enable_serological and self.sero_fusion is not None and serum_markers is not None:
            if not self.sero_fusion._fitted:
                return cfdna_scores  # Fall back to cfDNA-only
            return self.sero_fusion.predict_proba(cfdna_scores, serum_markers)
        return cfdna_scores


__all__ = [
    "SerologicalFusion",
    "IntegrativeScoringSystem",
]
