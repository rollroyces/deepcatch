"""
Clinical Module for DeepCatch.

- Serological fusion integrating traditional serum biomarkers
  (pepsinogen, gastrin-17, H. pylori) with cfDNA predictions.
- Frequency input: loads pre-computed motif frequency vectors
  (Jiang 4-mer data) for CET analysis.
- Plot generation: publication-quality visualisations (volcano,
  heatmap, ROC, feature importance).
- Nested cross-validation CET validator for unbiased motif-based
  feature selection and performance estimation.
"""

from .frequency_input import FrequencyDataset, PlotGenerator
from .serological_fusion import SerologicalFusion, IntegrativeScoringSystem
from .cet_cross_validator import (
    MotifRanking,
    EnrichmentProfile,
    MotifRanker,
    NestedCETValidator,
    compute_cliffs_delta,
    run_demo,
)
from .clinical_interpretation import ClinicalReportGenerator

__all__ = [
    # Frequency input
    "FrequencyDataset",
    "PlotGenerator",
    # Serological fusion
    "SerologicalFusion",
    "IntegrativeScoringSystem",
    # CET cross-validation
    "MotifRanking",
    "EnrichmentProfile",
    "MotifRanker",
    "NestedCETValidator",
    "compute_cliffs_delta",
    "run_demo",
    # Clinical interpretation
    "ClinicalReportGenerator",
]
