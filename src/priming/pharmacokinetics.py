#!/usr/bin/env python3
"""
Priming Agent Pharmacokinetics / Pharmacodynamics
===================================================

1-compartment PK model with first-order elimination for 5 priming agent types.
Models cfDNA clearance suppression and ctDNA concentration boost.

Adapted from:
- Martin-Alonso et al. (2024) Science: priming agents transiently suppress cfDNA clearance
- Gabrielsson & Weiner (2016) "Pharmacokinetic and Pharmacodynamic Data Analysis"

All half-lives and PK parameters sourced from published literature.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import AGENT_SPECS, AGENT_NAMES

logger = logging.getLogger(__name__)


class PKModel:
    """1-compartment pharmacokinetic model for priming agents.

    Simulates the time course of priming agent concentration and its
    effect on cfDNA clearance rate.

    Parameters
    ----------
    agent : str
        Agent type name (scFv, liposome, nanoparticle, polymeric_micelle, dendrimer).

    Notes
    -----
    Uses first-order elimination: dC/dt = -k_el * C
    where k_el = CL / Vd (clearance / volume of distribution)
    """

    def __init__(self):
        pass

    def simulate(
        self,
        agent: str,
        dose_mg: float,
        patient_weight_kg: float = 70.0,
        liver_function: float = 1.0,
        n_timepoints: int = 100,
        duration_hours: float = 48.0,
    ) -> dict:
        """Simulate PK time course for a single dose.

        Parameters
        ----------
        agent : str
            Agent type name.
        dose_mg : float
            Administered dose in mg.
        patient_weight_kg : float
            Patient weight in kg.
        liver_function : float
            Liver function score (0=severely impaired, 1=normal).
        n_timepoints : int
            Number of time points in output.
        duration_hours : float
            Simulation duration in hours.

        Returns
        -------
        dict with keys:
            time_hours : ndarray
            concentration_mg_L : ndarray
            clearance_rate : ndarray (cfDNA clearance rate multiplier)
            ctDNA_boost_factor : ndarray
            peak_concentration : float
            time_to_peak_h : float
            auc_mg_h_L : float
            effective_half_life_h : float
            cfdna_suppression_duration_h : float
        """
        if agent not in AGENT_NAMES:
            raise ValueError(f"Unknown agent: {agent}. Choose from {AGENT_NAMES}")

        if patient_weight_kg <= 0:
            logger.warning(f"patient_weight_kg={patient_weight_kg} ≤ 0; clamping to 1")
            patient_weight_kg = 1.0

        if dose_mg <= 0:
            logger.warning(f"dose_mg={dose_mg} ≤ 0; returning zero-effect result")
            return self._zero_result(n_timepoints, duration_hours)

        liver_function = np.clip(liver_function, 0.0, 1.0)

        spec = AGENT_SPECS[agent]
        k_el_nominal = spec["clearance_rate_L_h"] / spec["volume_of_distribution_L"]

        # Adjust elimination rate for liver function (impaired liver → slower clearance)
        liver_factor = 0.3 + 0.7 * liver_function  # range [0.3, 1.0]
        k_el = k_el_nominal * liver_factor

        # Weight-based volume adjustment
        vd = spec["volume_of_distribution_L"] * (patient_weight_kg / 70.0)
        f = spec["bioavailability"]

        # Initial concentration (IV bolus approximation)
        c0 = (f * dose_mg) / vd

        time = np.linspace(0, duration_hours, n_timepoints)
        concentration = c0 * np.exp(-k_el * time)

        # cfDNA clearance rate: 1 - normalized agent concentration
        # When agent is high → clearance suppressed → ctDNA accumulates
        clearance_rate = 1.0 - np.clip(concentration / (c0 + 1e-8), 0, 1)
        clearance_rate *= liver_factor  # baseline clearance also liver-dependent

        # ctDNA boost factor: inversely related to clearance
        # Maximum boost when clearance is maximally suppressed
        max_suppression = 1.0 - np.min(clearance_rate)
        ctDNA_boost = 1.0 + (max_suppression * 10.0) * (1.0 - clearance_rate)

        # Metrics
        peak_conc = float(np.max(concentration))
        time_to_peak = float(time[np.argmax(concentration)])
        auc = float(np.trapezoid(concentration, time))

        # Effective half-life (time to drop below 50% of peak)
        half_mask = concentration >= (peak_conc / 2.0)
        effective_half_life = float(time[half_mask][-1]) if np.any(half_mask) else 0.0

        # Duration of cfDNA suppression (>20% reduction in clearance)
        suppression_mask = clearance_rate < 0.8
        suppression_duration = (
            float(time[suppression_mask][-1] - time[suppression_mask][0])
            if np.any(suppression_mask) and np.sum(suppression_mask) > 1
            else 0.0
        )

        return {
            "time_hours": time,
            "concentration_mg_L": concentration,
            "clearance_rate": clearance_rate,
            "ctDNA_boost_factor": ctDNA_boost,
            "peak_concentration": peak_conc,
            "time_to_peak_h": time_to_peak,
            "auc_mg_h_L": auc,
            "effective_half_life_h": effective_half_life,
            "cfdna_suppression_duration_h": suppression_duration,
        }

    def _zero_result(self, n_timepoints: int, duration_hours: float) -> dict:
        """Return zero-effect simulation result."""
        time = np.linspace(0, duration_hours, n_timepoints)
        return {
            "time_hours": time,
            "concentration_mg_L": np.zeros(n_timepoints),
            "clearance_rate": np.ones(n_timepoints),
            "ctDNA_boost_factor": np.ones(n_timepoints),
            "peak_concentration": 0.0,
            "time_to_peak_h": 0.0,
            "auc_mg_h_L": 0.0,
            "effective_half_life_h": 0.0,
            "cfdna_suppression_duration_h": 0.0,
        }


class OptimalDosingSchedule:
    """Compute optimal dosing schedule for priming agents.

    Maximizes ctDNA boost while respecting toxicity constraints.
    Uses grid search over dose and timing parameters.

    Parameters
    ----------
    config : PrimingConfig or None
        Configuration object.
    """

    def __init__(self, config=None):
        from .config import PrimingConfig

        self.config = config or PrimingConfig()
        self.pk_model = PKModel()

    def compute(self, agent: str, patient_data: dict) -> dict:
        """Compute optimal dosing schedule for a patient.

        Parameters
        ----------
        agent : str
            Agent type name.
        patient_data : dict
            Patient features including weight, liver function, etc.
            Required keys: weight_kg, liver_function, tumor_type, tumor_stage.

        Returns
        -------
        dict with:
            optimal_dose_mg : float
            optimal_timing_h : float (time before blood draw)
            predicted_ctDNA_boost : float
            predicted_toxicity_risk : float
            schedule : list[dict] (dose + timing steps)
            effect_curve : ndarray (time points)
            is_feasible : bool
            recommendation : str
        """
        if agent not in AGENT_NAMES:
            raise ValueError(f"Unknown agent: {agent}. Choose from {AGENT_NAMES}")

        # Extract patient data with fallbacks
        weight = patient_data.get("weight_kg", 70.0)
        liver_fn = patient_data.get("liver_function", 1.0)
        tumor_type = patient_data.get("tumor_type", 0)
        tumor_stage = patient_data.get("tumor_stage", 1)
        baseline_cfdna = patient_data.get("baseline_cfdna_ng_ml", 20.0)

        if weight <= 0:
            weight = 70.0
            logger.warning("weight_kg ≤ 0, using default 70 kg")

        # Dose search: 5-50 mg/kg range in 10 steps
        spec = AGENT_SPECS[agent]
        dose_per_kg_candidates = np.linspace(0.5, 50.0, 20)
        best_score = -np.inf
        best_result = None

        ref_half_life = spec["half_life_hours"]

        for dose_per_kg in dose_per_kg_candidates:
            dose_mg = dose_per_kg * weight
            sim = self.pk_model.simulate(
                agent=agent,
                dose_mg=dose_mg,
                patient_weight_kg=weight,
                liver_function=liver_fn,
            )

            # Score = boost × feasibility penalty
            boost = float(np.max(sim["ctDNA_boost_factor"]))
            suppression_dur = sim["cfdna_suppression_duration_h"]

            # Toxicity proxy: dose above 5 mg/kg starts increasing risk
            toxicity_risk = self._estimate_toxicity(dose_per_kg, agent, liver_fn)
            if toxicity_risk > 0.8:
                continue  # Skip unsafe doses

            # Feasibility: suppression should last long enough for blood draw
            # but not too long (>24h impractical)
            timing_feasibility = np.clip(suppression_dur / 6.0, 0, 1) if 1 <= suppression_dur <= 24 else 0.3

            score = boost * (1.0 - 0.3 * toxicity_risk) * timing_feasibility

            if score > best_score:
                best_score = score
                best_result = {
                    "dose_mg": dose_mg,
                    "dose_per_kg": dose_per_kg,
                    "simulation": sim,
                    "toxicity_risk": toxicity_risk,
                }

        if best_result is None:
            # Fallback: minimal dose
            dose_mg = 0.5 * weight
            sim = self.pk_model.simulate(
                agent=agent, dose_mg=dose_mg,
                patient_weight_kg=weight, liver_function=liver_fn,
            )
            tox = self._estimate_toxicity(0.5, agent, liver_fn)
            best_result = {"dose_mg": dose_mg, "dose_per_kg": 0.5, "simulation": sim, "toxicity_risk": tox}

        sim = best_result["simulation"]
        peak_time = sim["time_to_peak_h"]

        # Optimal blood draw: at peak ctDNA boost
        boost_idx = np.argmax(sim["ctDNA_boost_factor"])
        optimal_timing = float(sim["time_hours"][boost_idx])
        predicted_boost = float(sim["ctDNA_boost_factor"][boost_idx])

        # Calculate schedule steps
        schedule_steps = []
        # Pre-dose: baseline
        schedule_steps.append({
            "time_h": -2.0,
            "action": "baseline_blood_draw",
            "cfDNA_clearance_pct": 100.0,
        })
        # Administer
        schedule_steps.append({
            "time_h": 0.0,
            "action": "administer_priming_agent",
            "dose_mg": round(best_result["dose_mg"], 1),
        })
        # Peak blood draw
        schedule_steps.append({
            "time_h": round(optimal_timing, 2),
            "action": "peak_blood_draw",
            "predicted_ctDNA_boost": round(predicted_boost, 2),
        })

        return {
            "optimal_dose_mg": round(best_result["dose_mg"], 2),
            "optimal_dose_per_kg": round(best_result["dose_per_kg"], 2),
            "optimal_timing_h": round(optimal_timing, 2),
            "predicted_ctDNA_boost": round(predicted_boost, 2),
            "predicted_toxicity_risk": round(best_result["toxicity_risk"], 3),
            "schedule": schedule_steps,
            "effect_curve": sim["ctDNA_boost_factor"],
            "is_feasible": predicted_boost >= 2.0,
            "recommendation": self._generate_recommendation(
                agent, predicted_boost, best_result["toxicity_risk"]
            ),
        }

    def _estimate_toxicity(
        self, dose_per_kg: float, agent: str, liver_function: float
    ) -> float:
        """Estimate toxicity risk (0-1) based on dose and agent type."""
        spec = AGENT_SPECS[agent]

        # Base toxicity from protein binding (higher PB → more off-target)
        base_tox = spec["protein_binding"] * 0.4

        # Dose-dependent toxicity (sigmoid)
        dose_tox = 1.0 / (1.0 + np.exp(-(dose_per_kg - 10.0) / 5.0))

        # Liver function modifier
        liver_tox = (1.0 - liver_function) * 0.5

        return float(np.clip(base_tox + 0.6 * dose_tox + liver_tox, 0.0, 1.0))

    def _generate_recommendation(
        self, agent: str, predicted_boost: float, toxicity_risk: float
    ) -> str:
        """Generate a human-readable recommendation."""
        if predicted_boost >= 5.0 and toxicity_risk < 0.3:
            return (
                f"Strong candidate for {agent} priming. "
                f"Predicted {predicted_boost:.1f}x ctDNA boost with low toxicity risk ({toxicity_risk:.2f}). "
                f"Proceed with standard protocol."
            )
        elif predicted_boost >= 2.0 and toxicity_risk < 0.6:
            return (
                f"Moderate candidate for {agent} priming. "
                f"Predicted {predicted_boost:.1f}x ctDNA boost with acceptable toxicity ({toxicity_risk:.2f}). "
                f"Consider monitoring liver enzymes post-administration."
            )
        else:
            return (
                f"Limited benefit expected from {agent} priming. "
                f"Predicted {predicted_boost:.1f}x boost with toxicity risk {toxicity_risk:.2f}. "
                f"Consider alternative agent or higher-sensitivity detection approach."
            )
