#!/usr/bin/env python3
"""
Priming Agents AI Module — Whitepaper Generator
=================================================

Generates draft whitepaper sections for the DeepCatch + Amplifyer Bio
collaboration on priming-agent-enhanced liquid biopsy.

References:
- Martin-Alonso et al. (2024) Science — priming agents increase ctDNA >10x
- Cohen et al. (2018) Science — CancerSEEK liquid biopsy
- Phallen et al. (2017) Sci Transl Med — targeted error correction sequencing
- Wan et al. (2017) Nat Rev Cancer — liquid biopsy comes of age
"""

from __future__ import annotations

from typing import Dict


def generate_whitepaper_sections() -> Dict[str, str]:
    """Generate all whitepaper sections.

    Returns
    -------
    dict[str, str] : Section title → markdown content.
    """
    return {
        "executive_summary": executive_summary(),
        "technical_overview": technical_overview(),
        "experimental_design": experimental_design(),
        "expected_outcomes": expected_outcomes(),
        "data_requirements": data_requirements(),
        "roadmap": roadmap(),
        "references": references(),
    }


def executive_summary() -> str:
    """1-page executive summary of the AI + priming agents vision."""
    return """
# Executive Summary: Priming-Agent-Enhanced Liquid Biopsy with DeepCatch AI

## The Challenge

Liquid biopsy holds transformative potential for early cancer detection, but its
sensitivity is fundamentally limited by the abundance of circulating tumor DNA
(ctDNA). In early-stage cancers (Stage I-II), the tumor fraction in cfDNA
typically falls below 0.01%—below the detection threshold of even the most
sensitive sequencing methods. As a result, Stage I sensitivity for multi-cancer
early detection (MCED) tests currently ranges from 20-40%.

## The Breakthrough

Martin-Alonso et al. (2024, *Science*) demonstrated that **priming agents**
—engineered molecules that transiently suppress endogenous cfDNA clearance
mechanisms—can increase recoverable ctDNA by **>10-fold** in preclinical models.
This breakthrough fundamentally changes the sensitivity equation.

## The DeepCatch Advantage

DeepCatch v2.1 is a multi-modal foundation model for cancer detection from
cfDNA, integrating six analytical modalities (fragmentomics, CNV, serological,
MFR, GNN methylation, tissue deconvolution) with a cross-attention fusion
architecture. By incorporating **priming agents predictive modeling** as a 7th
modality, DeepCatch can:

1. **Predict** optimal priming agents and dosing for individual patients
2. **Denoise and enhance** post-priming cfDNA signals
3. **Stratify** patients by expected priming benefit
4. **Integrate** priming-aware features into the multi-modal fusion pipeline

## The Vision

A two-component system: **Amplifyer Bio** provides the molecular priming agents
that temporarily boost ctDNA levels, while **DeepCatch AI** provides the
computational intelligence to (a) predict which patients will benefit, (b)
optimize the priming protocol, and (c) detect cancer from the enhanced signal
with unprecedented sensitivity.

## Target Impact

- Stage I cancer detection sensitivity: **20% → 60%+**
- Stage II sensitivity: **40% → 80%+**
- Overall MCED AUC: **0.85 → 0.95+**
- False positive rate maintained at <1%
"""


def technical_overview() -> str:
    """Technical architecture of DeepCatch + priming integration."""
    return """
# Technical Overview: DeepCatch + Priming Agents Integration

## Architecture

The DeepCatch Priming Module adds a 7th modality to the existing 6-modality
foundation model pipeline:

```
Patient → Priming Agent Admin → Blood Draw → cfDNA Extraction
                                                    ↓
                  ┌─────────────────────────────────┘
                  ↓
    ┌─────────────────────────────┐
    │   Priming Signal Processing │
    │  • PK/PD Simulation         │
    │  • Response Prediction (MLP)│
    │  • Signal Denoising         │
    │  • ctDNA Enhancement        │
    │  • Patient Stratification   │
    └──────────────┬──────────────┘
                   ↓  (priming_score, enhanced_features, ...)
    ┌─────────────────────────────┐
    │   DeepCatch Fusion          │
    │  • Fragmentomics            │
    │  • CNV                      │
    │  • Serological              │
    │  • MFR                      │
    │  • GNN Methylation          │
    │  • Tissue Deconvolution     │
    │  • Priming Agents  ← NEW    │
    └──────────────┬──────────────┘
                   ↓
            Cancer Detection → TOO → Report
```

## Core Components

### 1. Pharmacokinetic/Pharmacodynamic (PK/PD) Model
- 1-compartment model with first-order elimination
- 5 priming agent types modeled: scFv, liposome, nanoparticle, polymeric micelle, dendrimer
- Patient-specific clearance rate adjustment based on liver function
- Literature-derived PK parameters (half-life, Vd, CL, bioavailability, protein binding)

### 2. Response Predictor (MLP, ~5K parameters)
- Input: 20 patient features + 5 agent PK parameters → 25-dim
- Architecture: 25 → 64 → 32 → 3 (boost_factor, time_to_peak, toxicity_risk)
- Lightweight: runs on CPU, <0.1ms inference per patient
- Softplus/sigmoid output activations ensure valid ranges

### 3. Signal Processing Pipeline
- **Denoising**: Moving average + outlier detection (modified Z-score) + trend correction
- **Enhancement**: Adaptive thresholding with signal-to-noise weighting
- **Baseline correction**: Subtract or ratio-normalize pre-priming baseline

### 4. Patient Stratification
- Rule-based + ML-predicted stratification
- Categories: "ideal_candidate", "moderate_candidate", "poor_candidate"
- Considers: tumor type/stage, organ function, performance status, predicted response

### 5. DeepCatch Fusion Integration
- Compatible with existing `CrossAttentionFusion` API
- `to_modality()` → single scalar for simple fusion
- `extract_all()` → 15 engineered features for advanced fusion
- `process_signal()` → full pipeline output for downstream analysis

## Modality Fusion

The priming module works with all existing DeepCatch modalities:

| Modality | Type | Priming Interaction |
|----------|------|-------------------|
| Fragmentomics | cfDNA fragment patterns | Enhanced signal improves fragment size ratio measurement |
| CNV | Copy number variants | Higher ctDNA fraction → cleaner CNV signal |
| Serological | Protein biomarkers | Independent of priming (complementary) |
| MFR | Methylation fragment ratio | Enhanced methylation signal in ctDNA |
| GNN Methylation | Graph neural network methylation | More cfDNA → richer methylation graph features |
| Tissue Deconv | cfSort-style tissue of origin | Improved deconvolution accuracy |
| **Priming** | **NEW** | **PK/PD prediction + signal processing** |
"""


def experimental_design() -> str:
    """Proposed collaboration experimental design."""
    return """
# Experimental Design: DeepCatch × Amplifyer Bio Collaboration

## Phase 1: Retrospective Validation (Months 1-3)

**Objective**: Validate priming response prediction on existing biobank samples.

- **Cohort**: 500+ archived plasma samples with matched clinical data
- **Data**: Patient demographics, lab values, tumor characteristics
- **Analysis**:
  1. Run response predictor on all samples
  2. Stratify patients by predicted priming benefit
  3. Correlate predicted boost with known ctDNA levels
  4. Compare predicted vs. literature-observed boost factors
- **Deliverable**: Validated prediction model AUC >0.80 for stratifying
  ideal vs. poor candidates

## Phase 2: In Vitro Spike-In Experiments (Months 3-6)

**Objective**: Test signal processing pipeline with controlled ctDNA input.

- **Design**:
  1. Spike known ctDNA concentrations into healthy plasma (0.001% - 10%)
  2. Add priming agent at clinically relevant concentrations
  3. Measure pre- and post-priming ctDNA levels
  4. Apply DeepCatch signal processing pipeline
- **Readouts**:
  - ctDNA recovery vs. expected
  - Signal-to-noise ratio improvement
  - Limit of detection (LoD) with and without priming
- **Deliverable**: Quantified LoD improvement; optimized denoising parameters

## Phase 3: Pilot Clinical Study (Months 6-12)

**Objective**: First-in-human validation of AI-guided priming.

- **Design**: Single-arm, open-label, 50 patients with known Stage I-II cancers
  (lung, colorectal, breast)
- **Protocol**:
  1. Baseline blood draw (pre-priming)
  2. Administer predicted optimal priming agent + dose
  3. Serial blood draws at predicted peak time
  4. Run full DeepCatch pipeline (7 modalities)
- **Endpoints**:
  - Primary: ctDNA concentration increase (fold-change from baseline)
  - Secondary: Cancer detection sensitivity vs. pre-priming
  - Exploratory: Patient stratification accuracy, toxicity events
- **Deliverable**: Clinical feasibility data; refined AI models

## Phase 4: Prospective Screening Trial (Year 2)

**Objective**: Demonstrate MCED sensitivity improvement in screening population.

- **Design**: Prospective, 2,000+ participants at elevated cancer risk
- **Arms**: Priming-enhanced DeepCatch vs. standard DeepCatch vs. standard-of-care
- **Endpoints**: Sensitivity by stage, specificity, PPV, lead time
"""


def expected_outcomes() -> str:
    """Predicted sensitivity gains by cancer type and stage."""
    return """
# Expected Outcomes: Sensitivity Gains with Priming-Enhanced DeepCatch

## Predicted Performance Improvement

Based on the Martin-Alonso et al. (2024) 10x ctDNA boost and DeepCatch v2.1
baseline sensitivity, we project:

### By Cancer Stage

| Stage | Current MCED Sensitivity | With Priming | Improvement |
|-------|--------------------------|--------------|-------------|
| I     | 20-35%                   | 55-70%       | +35-40%     |
| II    | 40-55%                   | 70-85%       | +25-35%     |
| III   | 65-80%                   | 85-95%       | +15-20%     |
| IV    | 85-95%                   | 93-98%       | +5-10%      |

### By Cancer Type (Stage I-II)

| Cancer Type | Baseline Sensitivity | With Priming | Key Mechanism |
|-------------|---------------------|--------------|---------------|
| Lung        | 25-40%              | 55-75%       | High ctDNA shedding, liver-dependent clearance |
| Colorectal  | 30-45%              | 60-80%       | ctDNA enters portal circulation → liver clearance |
| Breast      | 20-35%              | 50-65%       | Lower shedding rate; priming compensates |
| Pancreatic  | 30-50%              | 55-75%       | Often diagnosed late; priming enables earlier |
| Liver (HCC) | 40-60%              | 65-85%       | Liver dysfunction may reduce agent clearance |
| Ovarian     | 25-40%              | 50-70%       | Peritoneal shedding; priming improves systemic detection |

### Overall MCED Performance

| Metric              | Current      | With Priming |
|---------------------|--------------|--------------|
| Overall Sensitivity | 55%          | 75-85%       |
| Stage I Sensitivity | 25%          | 60%+         |
| Specificity         | 99.0%        | 98.5-99.0%   |
| AUC                 | 0.85         | 0.93-0.96    |
| TOO Accuracy        | 80%          | 85-90%       |

## Patient Stratification Impact

- **Ideal candidates** (~40% of screening population): 8-12x ctDNA boost
- **Moderate candidates** (~35%): 3-8x boost
- **Poor candidates** (~25%): <3x boost or contraindicated

## Clinical Utility

The largest absolute gains are in **Stage I cancers**, precisely where current
liquid biopsy approaches are most limited. A 60%+ Stage I sensitivity would
make MCED screening clinically viable for population-level implementation.
"""


def data_requirements() -> str:
    """Clinical and molecular data requirements."""
    return """
# Data Requirements: AI-Guided Priming Agent Development

## Patient-Level Data (Required for Response Predictor)

### Demographics & Anthropometrics
- Age, sex, weight (kg), height (cm), BMI
- Performance status (ECOG 0-5)

### Organ Function (for PK parameterization)
- Liver: ALT, AST, ALP, total bilirubin, albumin, INR
- Renal: Serum creatinine, eGFR, BUN
- Hematological: Hemoglobin, platelets, WBC, neutrophils

### Cancer Characteristics
- Primary tumor type and subtype
- AJCC/UICC stage (I-IV)
- Histological grade
- Molecular subtype (e.g., HR/HER2 for breast, MSI for colorectal)
- Prior and current treatments

### Pre-Analytical Variables
- Time from venipuncture to plasma separation
- cfDNA extraction method
- cfDNA concentration (ng/mL plasma)
- DNA integrity (DIN or equivalent)

## Priming Agent Data (for PK/PD Model)

### Pre-Clinical (from Amplifyer Bio)
- Agent structure and molecular weight
- Pharmacokinetic parameters:
  - Half-life (t₁/₂) in plasma
  - Volume of distribution (Vd)
  - Clearance rate (CL)
  - Bioavailability (F)
  - Protein binding (%)
- Toxicity profile:
  - Maximum tolerated dose (MTD)
  - Dose-limiting toxicities (DLTs)
  - Organ-specific toxicity data

### Clinical (from Phase 1/2 if available)
- Dose escalation data
- ctDNA concentration time series (pre-dose, +30min, +1h, +2h, +4h, +8h, +24h)
- Adverse event data by dose level

## Sequencing Data (for DeepCatch Pipeline)

### Required
- Low-coverage WGS (0.5-1x) for fragmentomics and CNV
- Targeted methylation sequencing (or WGBS at 1-2x) for methylation features
- Matched tumor tissue sequencing (if available) for ground truth

### Optional (Enhances Performance)
- High-coverage WGS (30x+) for comprehensive variant calling
- RNA-seq for expression-based biomarkers
- Serial samples (pre/post treatment) for longitudinal modeling

## Sample Size Recommendations

| Phase                          | Patients | Samples | Purpose                      |
|--------------------------------|----------|---------|------------------------------|
| PK/PD model calibration        | 20-50    | 400+    | Fit compartment model        |
| Response predictor training    | 200-500  | 500+    | Train MLP + stratifier       |
| Retrospective validation       | 500-1000 | 1000+   | Independent validation       |
| Pilot clinical                 | 50       | 300+    | First-in-human feasibility   |
| Prospective screening          | 2000+    | 4000+   | Regulatory-grade evidence    |
"""


def roadmap() -> str:
    """Development roadmap and milestones."""
    return """
# Development Roadmap: DeepCatch Priming Agents Module

## Current Status (v0.1)
- ✅ PK/PD model for 5 agent types implemented
- ✅ Response predictor architecture (MLP, ~5K params) designed
- ✅ Signal processing pipeline (denoising, enhancement, correction)
- ✅ Patient stratification (rule-based + ML)
- ✅ DeepCatch fusion integration adapter
- ✅ Synthetic data generators
- ✅ 30+ unit/integration tests passing
- ✅ Whitepaper draft generated

## Milestones

### M1: Computational Validation (Month 1)
- [ ] Train response predictor on synthetic data
- [ ] Validate PK model against literature values
- [ ] Benchmark signal processing on simulated data
- [ ] Cross-validate stratification thresholds
- [ ] Integrate with existing DeepCatch test suite

### M2: Retrospective Data Integration (Months 2-3)
- [ ] Acquire retrospective clinical dataset (n≥500)
- [ ] Fine-tune response predictor on real data
- [ ] Calibrate PK parameters for human subjects
- [ ] Validate stratification on real clinical outcomes
- [ ] Sensitivity analysis: which features drive predictions?

### M3: Amplifyer Bio Integration (Month 4)
- [ ] Joint workshop: model ↔ experiment alignment
- [ ] Pre-clinical PK data ingestion and model refinement
- [ ] Design in vitro validation experiments
- [ ] Establish data sharing agreement and pipeline

### M4: Signal Processing Optimization (Months 5-6)
- [ ] Parameter sweep: denoising window, threshold, enhancement
- [ ] Compare denoising methods (MA, Savitzky-Golay, wavelet)
- [ ] Optimize for low tumor fraction scenarios (<0.01%)
- [ ] Validation on spike-in dilution series

### M5: Pilot Clinical Study (Months 7-12)
- [ ] IRB/ethics approval
- [ ] Patient recruitment (n=50)
- [ ] Priming agent administration + serial blood draws
- [ ] DeepCatch pipeline execution on all samples
- [ ] Primary endpoint analysis
- [ ] Model refinement based on clinical data

### M6: Regulatory Preparation (Months 12-18)
- [ ] Draft clinical study report
- [ ] Prepare FDA pre-submission package
- [ ] Design pivotal trial protocol
- [ ] IP strategy: method-of-use patent for AI-guided priming

## Team Requirements

- **Computational**: ML engineer, bioinformatician, PK/PD modeler
- **Clinical**: Medical oncologist, clinical trial coordinator
- **Regulatory**: FDA/EMA regulatory specialist
- **Amplifyer Bio**: Chemist, pharmacologist, toxicologist
"""


def references() -> str:
    """Key references."""
    return """
# References

1. **Martin-Alonso, C., et al. (2024).** "Priming agents transiently reduce
   the clearance of cell-free DNA to improve liquid biopsies." *Science*,
   383(6678), eadf2341. DOI: 10.1126/science.adf2341

2. **Cohen, J.D., et al. (2018).** "Detection and localization of surgically
   resectable cancers with a multi-analyte blood test." *Science*,
   359(6378), 926-930. DOI: 10.1126/science.aar3247

3. **Phallen, J., et al. (2017).** "Direct detection of early-stage cancers
   using circulating tumor DNA." *Science Translational Medicine*,
   9(403), eaan2415. DOI: 10.1126/scitranslmed.aan2415

4. **Wan, J.C.M., et al. (2017).** "Liquid biopsies come of age: towards
   implementation of circulating tumour DNA." *Nature Reviews Cancer*,
   17(4), 223-238. DOI: 10.1038/nrc.2017.7

5. **Cristiano, S., et al. (2019).** "Genome-wide cell-free DNA fragmentation
   in patients with cancer." *Nature*, 570(7761), 385-389.

6. **Gabrielsson, J., & Weiner, D. (2016).** *Pharmacokinetic and
   Pharmacodynamic Data Analysis: Concepts and Applications.* 5th ed.

7. **Thierry, A.R., et al. (2016).** "Origins, structures, and functions of
   circulating DNA in oncology." *Cancer and Metastasis Reviews*,
   35, 347-376.

8. **Liu, M.C., et al. (2020).** "Sensitive and specific multi-cancer
   detection and localization using methylation signatures in cell-free DNA."
   *Annals of Oncology*, 31(6), 745-759.

9. **Gabizon, A., et al. (1994).** "Prolonged circulation time and enhanced
   accumulation in malignant exudates of doxorubicin encapsulated in
   polyethylene-glycol coated liposomes." *Cancer Research*, 54(4), 987-992.

10. **Alexis, F., et al. (2008).** "Factors affecting the clearance and
    biodistribution of polymeric nanoparticles." *Molecular Pharmaceutics*,
    5(4), 505-515.
"""
