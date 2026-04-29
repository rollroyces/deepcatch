"""
Clinical Module for DeepCatch.

Serological fusion integrating traditional serum biomarkers
(pepsinogen, gastrin-17, H. pylori) with cfDNA predictions.
"""

from .serological_fusion import SerologicalFusion, IntegrativeScoringSystem

__all__ = [
    "SerologicalFusion",
    "IntegrativeScoringSystem",
]
