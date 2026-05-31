#!/usr/bin/env python3
"""
Priming Agents Data Generation
================================

Synthetic data generators for priming scenarios:
- Patient profile generation
- Priming effect simulation
- Clinical trial dataset generation

All PK parameters sourced from published literature.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from .config import AGENT_SPECS, AGENT_NAMES

logger = logging.getLogger(__name__)

# ── Tumor Type Constants ─────────────────────────────────────────

TUMOR_TYPES = {
    0: "lung",
    1: "colorectal",
    2: "breast",
    3: "pancreatic",
    4: "liver",
    5: "gastric",
    6: "ovarian",
    7: "esophageal",
    8: "head_neck",
    9: "other",
}

# Baseline ctDNA levels by stage (ng/mL) — from Cohen et al. 2018, Phallen et al. 2017
BASELINE_CFDNA_BY_STAGE = {
    0: (5, 20),
    1: (10, 40),
    2: (15, 80),
    3: (30, 200),
    4: (50, 500),
}


def generate_patient_profiles(
    n_patients: int, seed: int = 42
) -> List[dict]:
    """Generate synthetic patient profiles for priming simulations.

    Parameters
    ----------
    n_patients : int
        Number of patient profiles to generate.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    list[dict] : List of patient data dictionaries.
    """
    rng = np.random.RandomState(seed)
    profiles = []

    for i in range(n_patients):
        # Demographics
        age = int(rng.uniform(35, 85))
        weight_kg = rng.normal(72, 15)
        weight_kg = float(np.clip(weight_kg, 35, 150))
        bmi = rng.normal(26, 5)
        bmi = float(np.clip(bmi, 16, 50))

        # Organ function (correlated with age)
        liver_fn = float(np.clip(rng.normal(0.9 - (age - 40) * 0.003, 0.15), 0.1, 1.0))
        renal_fn = float(np.clip(rng.normal(0.95 - (age - 40) * 0.002, 0.12), 0.1, 1.0))

        # Tumor characteristics
        tumor_type = int(rng.choice(list(TUMOR_TYPES.keys())))
        # 60% stage I-II, 40% stage III-IV (screening population bias)
        tumor_stage = int(rng.choice([1, 2, 3, 4], p=[0.35, 0.25, 0.20, 0.20]) if rng.random() < 0.9 else 0)

        # Tumor fraction (low for early stage, higher for late)
        if tumor_stage == 0:
            tumor_fraction = 0.0
        elif tumor_stage == 1:
            tumor_fraction = float(10 ** rng.uniform(-4, -2))  # 0.01% - 1%
        elif tumor_stage == 2:
            tumor_fraction = float(10 ** rng.uniform(-3, -1))  # 0.1% - 10%
        elif tumor_stage == 3:
            tumor_fraction = float(10 ** rng.uniform(-2, -0.5))  # 1% - 31%
        else:
            tumor_fraction = float(10 ** rng.uniform(-1.5, 0))  # 3% - 100%

        # cfDNA baseline
        cfdna_range = BASELINE_CFDNA_BY_STAGE.get(tumor_stage, (10, 100))
        baseline_cfdna = float(rng.uniform(cfdna_range[0], cfdna_range[1]))

        # Lab values
        albumin = float(rng.normal(42, 5))
        albumin = float(np.clip(albumin, 22, 55))
        bilirubin = float(np.clip(rng.exponential(8), 2, 45))
        creatinine = float(np.clip(rng.normal(80, 25), 35, 450))
        alt = float(np.clip(rng.exponential(20), 5, 180))
        ast = float(np.clip(rng.exponential(22), 5, 190))
        platelets = float(rng.normal(250, 70))
        platelets = float(np.clip(platelets, 55, 580))
        neutrophils = float(rng.normal(4.5, 2))
        neutrophils = float(np.clip(neutrophils, 1.0, 18))
        hemoglobin = float(rng.normal(135, 15))
        hemoglobin = float(np.clip(hemoglobin, 75, 175))
        crp = float(np.clip(rng.exponential(5), 0, 180))

        prior_treatment = int(rng.choice([0, 1, 2, 3], p=[0.55, 0.25, 0.10, 0.10]))
        performance_status = int(np.clip(rng.poisson(1), 0, 4))

        profile = {
            "patient_id": f"PT-{i + 1:04d}",
            "age": age,
            "weight_kg": weight_kg,
            "bmi": bmi,
            "liver_function": liver_fn,
            "renal_function": renal_fn,
            "tumor_type": tumor_type,
            "tumor_stage": tumor_stage,
            "tumor_fraction": tumor_fraction,
            "baseline_cfdna_ng_ml": baseline_cfdna,
            "albumin_g_L": albumin,
            "bilirubin_umol_L": bilirubin,
            "creatinine_umol_L": creatinine,
            "alt_IU_L": alt,
            "ast_IU_L": ast,
            "platelet_count": platelets,
            "neutrophil_count": neutrophils,
            "hemoglobin_g_L": hemoglobin,
            "crp_mg_L": crp,
            "prior_treatment": prior_treatment,
            "performance_status": performance_status,
        }
        profiles.append(profile)

    return profiles


def apply_priming_effect(
    features: np.ndarray,
    agent: str,
    dose_mg: float,
    patient_weight_kg: float = 70.0,
    liver_function: float = 1.0,
    noise_level: float = 0.05,
    seed: Optional[int] = None,
) -> dict:
    """Simulate the effect of a priming agent on cfDNA features.

    Parameters
    ----------
    features : np.ndarray
        Original cfDNA features (pre-priming), shape (n_features,).
    agent : str
        Agent type name.
    dose_mg : float
        Administered dose in mg.
    patient_weight_kg : float
        Patient weight.
    liver_function : float
        Liver function score.
    noise_level : float
        Relative noise level to add.
    seed : int, optional
        Random seed.

    Returns
    -------
    dict with:
        pre_priming_features : ndarray
        post_priming_features : ndarray
        boost_factor : float
        peak_time_h : float
        signal_to_noise : float
        agent : str
        dose_mg : float
    """
    rng = np.random.RandomState(seed)
    features = np.asarray(features, dtype=np.float64)

    if agent not in AGENT_NAMES:
        raise ValueError(f"Unknown agent: {agent}. Choose from {AGENT_NAMES}")

    if dose_mg <= 0:
        return {
            "pre_priming_features": features.copy(),
            "post_priming_features": features.copy(),
            "boost_factor": 1.0,
            "peak_time_h": 0.0,
            "signal_to_noise": 1.0,
            "agent": agent,
            "dose_mg": 0.0,
        }

    spec = AGENT_SPECS[agent]

    # Compute expected boost from PK model
    from .pharmacokinetics import PKModel

    pk = PKModel()
    sim = pk.simulate(
        agent=agent,
        dose_mg=dose_mg,
        patient_weight_kg=patient_weight_kg,
        liver_function=liver_function,
    )

    boost_factor = float(np.max(sim["ctDNA_boost_factor"]))
    peak_time = sim["time_to_peak_h"]

    # boost_factor affects ctDNA features proportionally to their magnitude
    # Larger features (true ctDNA signal) get boosted more
    boosted = features.copy()

    # ctDNA signal amplification: features above 50th percentile get boosted
    threshold = np.percentile(np.abs(features), 50) if len(features) > 1 else 0
    signal_mask = np.abs(features) >= threshold

    # True ctDNA features are boosted by ~boost_factor
    # Background features get some boost too (noise amplification)
    boosted[signal_mask] *= boost_factor
    boosted[~signal_mask] *= (1.0 + 0.1 * (boost_factor - 1.0))  # noise gets 10% of boost

    # Add biological/technical noise
    noise = rng.normal(0, noise_level * np.std(features) + 1e-8, size=features.shape)
    boosted += noise

    # Compute effective SNR
    signal_power = np.var(boosted[signal_mask]) if np.any(signal_mask) else 0
    noise_power = np.var(boosted[~signal_mask]) if np.any(~signal_mask) else 1e-10
    effective_snr = float(signal_power / max(noise_power, 1e-10))

    return {
        "pre_priming_features": features.copy(),
        "post_priming_features": boosted,
        "boost_factor": boost_factor,
        "peak_time_h": peak_time,
        "signal_to_noise": effective_snr,
        "agent": agent,
        "dose_mg": dose_mg,
    }


def simulate_clinical_trial(
    n_patients: int = 100,
    n_agents: int = 5,
    n_features: int = 50,
    seed: int = 42,
) -> tuple:
    """Generate a complete simulated clinical trial dataset.

    Creates paired pre/post priming measurements with ground truth labels
    for ctDNA detection.

    Parameters
    ----------
    n_patients : int
        Number of patients.
    n_agents : int
        Number of agent types to test per patient.
    n_features : int
        Number of cfDNA features per sample.
    seed : int
        Random seed.

    Returns
    -------
    tuple:
        patients : list[dict]  — patient profiles
        trial_data : dict  — structured trial dataset
        summary : dict  — trial summary statistics
    """
    rng = np.random.RandomState(seed)
    patients = generate_patient_profiles(n_patients, seed=seed)

    agents = list(AGENT_NAMES)[: min(n_agents, len(AGENT_NAMES))]

    pre_samples = []   # shape (n_patients, n_features)
    post_samples = []  # shape (n_patients * n_agents, n_features)
    labels = []        # cancer (1) or healthy (0)
    boost_factors = [] # per patient-agent
    patient_ids = []
    agent_ids = []
    stages = []

    for i, patient in enumerate(patients):
        # Generate baseline cfDNA features
        tumor_fraction = patient["tumor_fraction"]
        stage = patient["tumor_stage"]

        # Base features: random with signal proportional to tumor fraction
        base_features = rng.normal(0, 1, n_features)
        if tumor_fraction > 0:
            # Inject ctDNA signal proportional to tumor fraction
            signal_features = rng.normal(0, 2, n_features)
            base_features = (1 - tumor_fraction) * base_features + tumor_fraction * signal_features

        pre_samples.append(base_features)

        # Cancer label: 1 if tumor present (stage > 0 and fraction > 0.0001)
        has_cancer = stage > 0 and tumor_fraction > 1e-4
        labels.append(1.0 if has_cancer else 0.0)

        # Apply priming agents
        for agent in agents:
            # Dose: 1-20 mg/kg
            dose_per_kg = rng.uniform(1, 20)
            dose_mg = dose_per_kg * patient["weight_kg"]

            result = apply_priming_effect(
                features=base_features,
                agent=agent,
                dose_mg=dose_mg,
                patient_weight_kg=patient["weight_kg"],
                liver_function=patient["liver_function"],
                seed=seed + i * 10 + agents.index(agent),
            )

            post_samples.append(result["post_priming_features"])
            boost_factors.append(result["boost_factor"])
            patient_ids.append(patient["patient_id"])
            agent_ids.append(agent)
            stages.append(stage)

    pre_samples = np.array(pre_samples)
    post_samples = np.array(post_samples)
    labels = np.array(labels)
    boost_factors = np.array(boost_factors)

    # Summary statistics
    cancer_mask = labels == 1.0
    summary = {
        "n_patients": n_patients,
        "n_agents": len(agents),
        "n_features": n_features,
        "agents_tested": agents,
        "n_with_cancer": int(np.sum(cancer_mask)),
        "n_without_cancer": int(np.sum(~cancer_mask)),
        "mean_boost_factor": float(np.mean(boost_factors)),
        "median_boost_factor": float(np.median(boost_factors)),
        "max_boost_factor": float(np.max(boost_factors)),
        "boost_by_stage": {},
    }

    # Boost by stage
    stage_arr = np.array(stages)
    for s in range(5):
        mask = stage_arr == s
        if np.any(mask):
            summary["boost_by_stage"][int(s)] = {
                "mean": float(np.mean(boost_factors[mask])),
                "median": float(np.median(boost_factors[mask])),
                "n": int(np.sum(mask)),
            }

    trial_data = {
        "pre_samples": pre_samples,
        "post_samples": post_samples,
        "labels": labels,
        "boost_factors": boost_factors,
        "patient_ids": patient_ids,
        "agent_ids": agent_ids,
        "stages": stages,
    }

    return patients, trial_data, summary
