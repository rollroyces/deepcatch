"""cfRNA -> real-time CDM feature coupling (v0.11.0).

This module is the **read-side** counterpart to
``mrnavax_rare.safety.cfrna.predict_cfrna_persistence``. It consumes the
prediction shape emitted by the design side and lifts it into a
real-time chronic-disease-management (CDM) feature schema that
downstream clinical-decision tooling can act on.

Why a separate module?
----------------------

``mrnavax-rare`` is the mRNA-therapeutic design side -- it predicts how
long a designed CDS should persist in plasma. ``deepcatch.methylation``
is the cfRNA read side -- it observes how much cfRNA is actually
circulating in the patient's plasma sample and flags anomalies against
the prediction. The two packages share NO runtime code (no shared
imports); they couple over a strict :class:`CfRNAPersistenceInput`
shape. This keeps the read side deployable independently of the design
side (e.g. for sites that already have mrnavax-rare outputs in a
research database but don't need the full design pipeline).

Honest disclosure
-----------------

The v0.10.0 anomaly thresholds (``ANOMALY_RATIO_LOW = 0.5``,
``ANOMALY_RATIO_HIGH = 2.0``) are **research-use defaults**. They
will need clinical calibration before any production use. The
dosing-recommendation rules are heuristic; they are NOT a substitute
for clinician judgment.

Ported from ``rollroyces/deepcatch-methylation`` v0.10.0 (SHA 69f3ba32).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping, TypedDict


class CfRNAPersistenceInput(TypedDict):
    """The minimal input shape expected from mrnavax-rare.

    Mirrors ``mrnavax_rare.safety.cfrna.CfRNAPersistencePrediction.as_dict()``.
    Defined as a :class:`TypedDict` so this package does not need to
    import mrnavax-rare at runtime -- the contract is structural, not
    nominal.
    """

    predicted_half_life_hours: float
    persistence_category: str
    target_tissue: str
    codon_optimization: str
    feature_breakdown: Mapping[str, float]


ANOMALY_RATIO_LOW: float = 0.5
ANOMALY_RATIO_HIGH: float = 2.0
MIN_HALFLIFE_HOURS: float = 1.0

DosingRecommendation = str


@dataclass(frozen=True)
class RealTimeCDMObservation:
    """One cfRNA plasma observation lifted into the real-time CDM feature schema."""

    observed_cfRNA_concentration_pg_ml: float
    sampling_time_post_dose_hours: float
    dose_cfRNA_at_t0_pg_ml: float
    predicted_half_life_hours: float
    predicted_persistence_category: str
    expected_residual_fraction: float
    observed_residual_fraction: float
    observed_vs_expected_ratio: float
    anomaly_flag: bool
    dosing_recommendation: DosingRecommendation
    notes: tuple = field(default_factory=tuple)


__all__ = [
    "ANOMALY_RATIO_HIGH",
    "ANOMALY_RATIO_LOW",
    "CfRNAPersistenceInput",
    "DosingRecommendation",
    "MIN_HALFLIFE_HOURS",
    "RealTimeCDMObservation",
    "couple_cfrna_to_cdm",
]


def _validate_input(prediction: CfRNAPersistenceInput) -> tuple:
    """Validate the input shape; return ``(half_life, category)``."""
    half_life = float(prediction["predicted_half_life_hours"])
    category = str(prediction["persistence_category"]).strip().lower()
    if category not in {"high", "medium", "low"}:
        raise ValueError(
            f"Unknown persistence_category: {category!r}; expected one of "
            f"'high', 'medium', 'low'."
        )
    if half_life <= 0.0:
        raise ValueError(
            f"predicted_half_life_hours must be positive; got {half_life!r}."
        )
    if not math.isfinite(half_life):
        raise ValueError(
            f"predicted_half_life_hours must be finite; got {half_life!r}."
        )
    return half_life, category


def _compute_residuals(
    observed_conc: float,
    dose_at_t0: float,
    sampling_time_h: float,
    half_life_h: float,
) -> tuple:
    """Compute ``(expected, observed, ratio)``."""
    if dose_at_t0 <= 0.0:
        raise ValueError(
            f"dose_cfRNA_at_t0_pg_ml must be positive; got {dose_at_t0!r}."
        )
    if sampling_time_h < 0.0:
        raise ValueError(
            f"sampling_time_post_dose_hours must be non-negative; "
            f"got {sampling_time_h!r}."
        )
    effective_half_life = max(half_life_h, MIN_HALFLIFE_HOURS)
    expected_residual = math.exp(-math.log(2) * sampling_time_h / effective_half_life)
    observed_residual = observed_conc / dose_at_t0
    if expected_residual <= 0.0:
        ratio = float("inf")
    else:
        ratio = observed_residual / expected_residual
    return expected_residual, observed_residual, ratio


def _recommend(
    ratio: float,
    anomaly: bool,
    half_life_h: float,
    category: str,
) -> tuple:
    """Map ``(ratio, anomaly, half_life, category)`` to ``(rec, notes)``."""
    notes: list = []

    if half_life_h < MIN_HALFLIFE_HOURS:
        notes.append(
            f"predicted half-life {half_life_h} h is below MIN_HALFLIFE_HOURS "
            f"= {MIN_HALFLIFE_HOURS} h -- dosing held pending re-prediction."
        )
        return ("hold", tuple(notes))

    if anomaly:
        if ratio < ANOMALY_RATIO_LOW:
            notes.append(
                f"observed_vs_expected_ratio = {ratio:.3f} below "
                f"ANOMALY_RATIO_LOW = {ANOMALY_RATIO_LOW} -- faster clearance "
                f"than predicted; potential under-dosing."
            )
            if category == "low":
                notes.append(
                    "compound signal: predicted persistence is also 'low' -- "
                    "recommend escalating to clinical review."
                )
                return ("escalate", tuple(notes))
            return ("increase", tuple(notes))
        if ratio > ANOMALY_RATIO_HIGH:
            notes.append(
                f"observed_vs_expected_ratio = {ratio:.3f} above "
                f"ANOMALY_RATIO_HIGH = {ANOMALY_RATIO_HIGH} -- slower clearance "
                f"than predicted; potential accumulation / saturation."
            )
            if category == "low":
                notes.append(
                    "compound signal: predicted persistence is 'low' but "
                    "actual is much higher -- review for compliance / dosing "
                    "interval errors."
                )
                return ("escalate", tuple(notes))
            return ("decrease", tuple(notes))

    return ("maintain", tuple(notes))


def couple_cfrna_to_cdm(
    prediction: CfRNAPersistenceInput,
    observed_cfRNA_concentration_pg_ml: float,
    sampling_time_post_dose_hours: float,
    dose_cfRNA_at_t0_pg_ml: float,
) -> RealTimeCDMObservation:
    """Couple a mrnavax-rare cfRNA prediction to a real-time CDM observation."""
    half_life_h, category = _validate_input(prediction)
    expected_residual, observed_residual, ratio = _compute_residuals(
        observed_conc=observed_cfRNA_concentration_pg_ml,
        dose_at_t0=dose_cfRNA_at_t0_pg_ml,
        sampling_time_h=sampling_time_post_dose_hours,
        half_life_h=half_life_h,
    )
    anomaly = not (ANOMALY_RATIO_LOW <= ratio <= ANOMALY_RATIO_HIGH)
    rec, notes = _recommend(ratio, anomaly, half_life_h, category)
    return RealTimeCDMObservation(
        observed_cfRNA_concentration_pg_ml=float(observed_cfRNA_concentration_pg_ml),
        sampling_time_post_dose_hours=float(sampling_time_post_dose_hours),
        dose_cfRNA_at_t0_pg_ml=float(dose_cfRNA_at_t0_pg_ml),
        predicted_half_life_hours=half_life_h,
        predicted_persistence_category=category,
        expected_residual_fraction=expected_residual,
        observed_residual_fraction=observed_residual,
        observed_vs_expected_ratio=ratio,
        anomaly_flag=anomaly,
        dosing_recommendation=rec,
        notes=notes,
    )
