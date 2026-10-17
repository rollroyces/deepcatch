"""cfRNA read-side companion module for the methylation pathway.

This module is the **read-side** counterpart to ``mrnavax_rare`` (mRNA
therapeutic design side). It lifts plasma cfRNA persistence predictions
into a real-time chronic-disease-management (CDM) feature schema that
downstream clinical-decision tooling can act on.

The coupling keeps no runtime dependency on ``mrnavax_rare`` -- the
contract is a strict :class:`CfRNAPersistenceInput` TypedDict shape.

Ported from ``rollroyces/deepcatch-methylation`` (SHA 69f3ba32) -- see
``docs/ultrasensitive/PHASE0_PLAN.md`` and ``METHYLATION_PROJECT.md``
for context.
"""

from .cfrna_coupling import (
    ANOMALY_RATIO_HIGH,
    ANOMALY_RATIO_LOW,
    CfRNAPersistenceInput,
    DosingRecommendation,
    MIN_HALFLIFE_HOURS,
    RealTimeCDMObservation,
    couple_cfrna_to_cdm,
)
from .methylation_reader import (
    MethylationArrayInput,
    TISSUE_OF_ORIGIN_MARKERS,
    build_cfrna_persistence_input,
    extract_cfrna_signature,
    make_synthetic_liver_beta_values,
    make_synthetic_non_liver_beta_values,
)

__all__ = [
    # cfrna coupling (read-side)
    "ANOMALY_RATIO_HIGH",
    "ANOMALY_RATIO_LOW",
    "CfRNAPersistenceInput",
    "DosingRecommendation",
    "MIN_HALFLIFE_HOURS",
    "RealTimeCDMObservation",
    "couple_cfrna_to_cdm",
    # methylation array reader
    "MethylationArrayInput",
    "TISSUE_OF_ORIGIN_MARKERS",
    "build_cfrna_persistence_input",
    "extract_cfrna_signature",
    "make_synthetic_liver_beta_values",
    "make_synthetic_non_liver_beta_values",
]

__version__ = "0.11.0"
