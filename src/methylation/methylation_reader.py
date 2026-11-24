"""Methylation-array -> cfRNA persistence input read pipeline (v0.11).

This module is the **read-side** complement to the cfRNA coupling in
``deepcatch.methylation.cfrna_coupling``. It defines a standardized read
pipeline that ingests methylation array data (Illumina EPIC 850K or
processed beta values) and produces a
:class:`CfRNAPersistenceInput` that the coupling can consume.

Architecture
------------

::

    methylation array (IDAT / processed beta values)
        -> ``MethylationArrayInput`` (TypedDict)
            -> ``extract_cfrna_signature()`` (tissue-of-origin identification)
                -> ``build_cfrna_persistence_input()`` (normalizer)
                    -> ``CfRNAPersistenceInput`` (TypedDict)
                        -> ``couple_cfrna_to_cdm()`` (anomaly + dosing)
                            -> ``RealTimeCDMObservation`` (CDM feature schema)

Honest disclosure
-----------------

The v0.11 read pipeline is a **structural scaffold** that defines the
input contract for cfRNA persistence coupling. It does NOT call a
real methylation-array analysis pipeline (e.g. minfi, methylKit,
SeSAMe) and does NOT extract tissue-of-origin markers from raw IDAT
files. Real production deployment would need minfi / methylKit / SeSAMe,
a tissue-of-origin marker panel, m5C detection, and ddPCR-quantified
cfRNA concentration.

Ported from ``rollroyces/deepcatch-methylation`` v0.11 (SHA 69f3ba32).
"""

from __future__ import annotations

import math
from typing import Mapping, TypedDict


TISSUE_OF_ORIGIN_MARKERS: dict = {
    "liver": ("liver_marker_001", "liver_marker_002", "liver_marker_003"),
    "lung": ("lung_marker_001", "lung_marker_002", "lung_marker_003"),
    "cns": ("cns_marker_001", "cns_marker_002", "cns_marker_003"),
    "muscle": ("muscle_marker_001", "muscle_marker_002"),
    "cardiac": ("cardiac_marker_001", "cardiac_marker_002"),
    "default": ("generic_marker_001",),
}


class MethylationArrayInput(TypedDict):
    """The standard input shape for the v0.11 methylation read pipeline."""

    sample_id: str
    tissue: str
    beta_values: Mapping[str, float]
    cfrna_concentration_pg_ml: float
    sampling_time_post_dose_hours: float
    dose_cfrna_at_t0_pg_ml: float


_MIN_MARKERS_FOR_TISSUE_CALL: int = 2
_TISSUE_BETA_THRESHOLD: float = 0.5


def extract_cfrna_signature(
    beta_values: Mapping[str, float],
    tissue: str,
) -> tuple:
    """Identify the tissue-of-origin from methylation beta values."""
    markers = TISSUE_OF_ORIGIN_MARKERS.get(tissue, TISSUE_OF_ORIGIN_MARKERS["default"])
    hit_count = 0
    beta_sum = 0.0
    for marker in markers:
        if marker in beta_values:
            beta = float(beta_values[marker])
            if beta >= _TISSUE_BETA_THRESHOLD:
                hit_count += 1
                beta_sum += beta
    if hit_count >= _MIN_MARKERS_FOR_TISSUE_CALL:
        return tissue, beta_sum / hit_count
    return "unknown", 0.0


def build_cfrna_persistence_input(
    methylation_input: MethylationArrayInput,
) -> "CfRNAPersistenceInput":
    """Build a :class:`CfRNAPersistenceInput` from methylation data."""
    from .cfrna_coupling import CfRNAPersistenceInput  # noqa: PLC0415

    tissue = methylation_input["tissue"].strip().lower()
    beta_values = methylation_input.get("beta_values", {}) or {}
    cfrna_conc = float(methylation_input.get("cfrna_concentration_pg_ml", 0.0))

    called_tissue, mean_beta = extract_cfrna_signature(beta_values, tissue)

    if cfrna_conc <= 0.0:
        cfrna_conc = mean_beta * float(
            methylation_input.get("dose_cfrna_at_t0_pg_ml", 0.0)
        )

    if called_tissue == tissue and cfrna_conc > 0.0:
        persistence_category = "medium"
    elif cfrna_conc > 0.0:
        persistence_category = "low"
    else:
        persistence_category = "low"

    half_life_h = 8.0

    return CfRNAPersistenceInput(
        predicted_half_life_hours=half_life_h,
        persistence_category=persistence_category,
        target_tissue=tissue,
        codon_optimization="balanced",
        feature_breakdown={
            "mean_beta": mean_beta,
            "called_tissue": 1.0 if called_tissue == tissue else 0.0,
            "concentration_pg_ml": cfrna_conc,
        },
    )


def make_synthetic_liver_beta_values(tissue: str = "liver") -> dict:
    """Generate synthetic beta values for a tissue-derived cfRNA sample."""
    markers = TISSUE_OF_ORIGIN_MARKERS.get(tissue, TISSUE_OF_ORIGIN_MARKERS["default"])
    return {m: 0.85 for m in markers}


def make_synthetic_non_liver_beta_values(tissue: str = "lung") -> dict:
    """Generate synthetic beta values for a non-liver cfRNA sample.

    The liver markers are below the tissue-call threshold (<= 0.5)
    but the tissue-specific markers are above (>= 0.5) -- exercises
    the "tissue mismatch" path in ``extract_cfrna_signature``.
    """
    liver_markers = TISSUE_OF_ORIGIN_MARKERS["liver"]
    tissue_markers = TISSUE_OF_ORIGIN_MARKERS.get(
        tissue, TISSUE_OF_ORIGIN_MARKERS["default"]
    )
    return {
        **{m: 0.30 for m in liver_markers},
        **{m: 0.85 for m in tissue_markers},
    }


__all__ = [
    "MethylationArrayInput",
    "TISSUE_OF_ORIGIN_MARKERS",
    "build_cfrna_persistence_input",
    "extract_cfrna_signature",
    "make_synthetic_liver_beta_values",
    "make_synthetic_non_liver_beta_values",
]
