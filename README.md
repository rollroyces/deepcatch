# 🧬 DeepCatch: Performance‑Weighted Multi‑Modal Fusion for Ultra‑Early Cancer Detection from cfDNA

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)
[![Version: 2.1](https://img.shields.io/badge/Version-2.1-blue.svg)]()

**DeepCatch is an open‑source computational framework for cancer early detection from cfDNA, combining performance‑weighted multi‑modal fusion, cumulative evidence tracking, and tissue‑of‑origin prediction.** It is designed to enable independent validation and accelerate liquid‑biopsy research.

---

> **⚠️ RESEARCH‑ONLY NOTICE:** DeepCatch is research-stage software. v2.1 adds preliminary validation on 129 real human plasma samples (processed frequency data) from Jiang lab (CUHK) — see §9. All other performance numbers have been removed from this README as they were simulation-based and not clinically verified. Full wet-lab validation on raw sequencing data remains the essential next step.

---



## 2. How It Actually Works

DeepCatch is a computational pipeline that answers one question: **“Given blood‑based measurements from a patient, what is the probability that an early‑stage cancer is present, and if so, where is it?”** It answers this by combining three strategies that no other open‑source framework integrates.

### 2.1 Performance‑Weighted Multi‑Modal Fusion

**The problem:** A single blood measurement — say, a mutation at 0.05 % variant allele fraction — is often too noisy to trust alone.

**Our solution:** Measure five independent molecular signals from the same blood draw:

1. **Somatic mutations** (ctDNA variants)
2. **DNA methylation** (epigenetic silencing)
3. **Fragmentomics** (cfDNA fragment size and end‑motif patterns)
4. **Copy‑number alterations** (genomic instability)
5. **CTC count** (circulating tumour cell estimate)

Each modality produces a probability score. Rather than averaging them equally (the approach used by Bie et al. 2023), DeepCatch **weights each modality by its individual diagnostic power** — measured as the area under the ROC curve (AUC) on a held‑out validation set. Modalities that perform worse than random (AUC < 0.5) receive zero weight and are excluded entirely.

This performance‑weighted fusion is the core idea behind DeepCatch — different molecular signals are combined according to their reliability rather than averaged equally.

### 2.2 Cumulative Evidence Tracking (CET)

**The problem:** Cancer grows over time, but ctDNA levels at the earliest stages are often below the detection threshold of any single blood draw.

**Our solution:** Track the patient over multiple quarterly blood draws using a two‑stage architecture:

- **Stage 1 — Permissive CET:** A sequential probability ratio test (SPRT) accumulates evidence across all five modalities and all quarterly timepoints. This stage is designed for high sensitivity at moderate specificity.

- **Stage 2 — Confirmatory Fusion:** Only the flagged patients receive a higher‑depth targeted sequencing panel. Performance‑weighted fusion is applied at a strict cutoff calibrated for ultra‑high specificity.

This two-stage approach is designed to reduce false positives while maintaining sensitivity. Performance numbers will be established through ongoing validation work.

### 2.3 FragmentoSign: Fragmentomics Subsystem

DeepCatch's fragmentomics engine implements the DELFI (DNA Evaluation of Fragments for Early Interruption) and MDS (Motif Diversity Score) frameworks:

| Component | Method | Reference |
|-----------|--------|-----------|
| GC‑bias correction | LOESS local normalisation | Cristiano 2019, *Nature* |
| Fragment length model | 4‑component Gaussian Mixture Model (GMM): sub‑nucleosomal (~80 bp), mono‑ (~167 bp), di‑ (~334 bp), tri‑nucleosomal (~501 bp) | Snyder 2016, *Cell* |
| End‑motif analysis | 4‑mer extraction from BAM/FASTQ + MDS scoring | Jiang 2020, *Cancer Discovery* |
| Nucleosome positioning | CNN over TSS coverage profiles | Snyder 2016, *Cell* |

The sub‑nucleosomal GMM component is specifically designed to detect the increased proportion of short fragments (<150 bp) characteristic of tumour‑derived cfDNA.

### Novel Components

| Component | Innovation |
|-----------|-----------|
| Performance-weighted fusion | AUC-proportional weighting with below-chance suppression |
| Two-Stage CET (SPRT) | Longitudinal multi-modal screening architecture |
| FragmentoSign (GMM + MDS) | Combined fragment length + end motif pipeline |
| 6-confounder realism | CHIP, variable shedding, trinucleotide errors, GE, batch, inflammation |

### THEMIS Feature Equivalents

DeepCatch's fragmentomics subsystem (FragmentoSign) implements the following THEMIS-inspired features:

| Feature | FragmentoSign Equivalent | Method |
|---------|--------------------------|--------|
| **MFR** (Methylated Fragment Ratio) | Methylation entropy + LOESS-normalized coverage | CpG density scoring |
| **FSI** (Fragment Size Index) | Short/long fragment ratio + GMM sub-nucleosomal fraction | 4-component GMM |
| **CAFF** (Chromosomal Aneuploidy) | Copy Number Instability Index + CNA burden scoring | Whole-genome binning |
| **FEM** (Fragment End Motif) | 4-mer MDS (Motif Diversity Score) + end motif embeddings | Jiang 2020 protocol |

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         DeepCatch System Architecture                    │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐ │
│  │  Variant      │  │  Methylation │  │  Fragment-   │  │  Copy Number │ │
│  │  Calling      │  │  (Entropy)   │  │  omics        │  │  Alterations │ │
│  │  (Bayesian +  │  │              │  │  (End motifs, │  │  (CNA)       │ │
│  │   Contrastive │  │              │  │   size, nucl) │  │              │ │
│  │   DL)         │  │              │  │  [Fragmento-  │  │              │ │
│  │               │  │              │  │   Sign]       │  │              │ │
│  └──────┬────────┘  └──────┬───────┘  └───────┬───────┘  └──────┬───────┘ │
│         │                  │                   │                  │        │
│         └──────────────────┼───────────────────┼──────────────────┘        │
│                            ▼                   ▼                           │
│               ┌────────────────────────────────────┐                       │
│               │  Multi‑Modal Fusion Layer           │                       │
│               │  (Performance‑Weighted)              │                       │
│               └──────────────┬─────────────────────┘                       │
│                              ▼                                             │
│               ┌────────────────────────────────────┐                       │
│               │  Two‑Stage CET Screening            │                       │
│               │  Stage 1: Permissive SPRT            │                       │
│               │  Stage 2: Strict Confirmation        │                       │
│               │                                                     │                       │
│               └──────────────┬─────────────────────┘                       │
│                              ▼                                             │
│               ┌────────────────────────────────────┐                       │
│               │  Meta‑Learning Ensemble (MAML)      │                       │
│               │  → Few‑shot subtype adaptation      │                       │
│               └──────────────┬─────────────────────┘                       │
│                              ▼                                             │
│               ┌────────────────────────────────────┐                       │
│               │   Risk Score + TOO Prediction       │                       │
│               │   → 4‑Tier Risk Stratification      │                       │
│               └────────────────────────────────────┘                       │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 3. Installation & Usage

### 3.1 Prerequisites

- **Python** 3.9+ with pip
- **Node.js** 18+ (validation scripts)
- **Git**

### 3.2 Quick Install

```bash
git clone https://github.com/rollroyces/deepcatch.git
cd deepcatch
pip install -r requirements.txt
```

### 3.3 Python Dependencies

```
numpy>=1.24.0
scipy>=1.10.0
scikit-learn>=1.3.0
matplotlib>=3.7.0
seaborn>=0.12.0
pandas>=2.0.0
torch>=2.0.0          # for deep learning models
torch-geometric>=2.3.0 # for GNN fusion
pysam                  # for BAM motif extraction (optional)
statsmodels            # for LOESS normalisation (optional)
```

### 3.4 Running the Pipeline

**Full validation suite (Python):**
```bash
bash RUN_ALL.sh
```

**Quick smoke test (2 minutes):**
```bash
bash RUN_ALL.sh --quick
```

**Node.js validation against TCGA/COSMIC data:**
```bash
cd validation/node
node runRealFinal.js
```

**Run a single module:**
```bash
# Two‑Stage CET
node validation/node/twoStageCET.js

# Head‑to‑head comparison
python validation/py/head_to_head.py

# Fragmentomics (FragmentoSign)
python -c "
from src.fragmentomics import FragmentLengthGMM, compute_fragmentomics_features
gmm = FragmentLengthGMM(n_components=4)
gmm.fit(your_fragment_lengths)
print(gmm.get_component_stats())
"
```

### 3.5 Docker

```bash
docker build -t deepcatch:latest .
docker run --rm -v $(pwd)/results:/app/results deepcatch:latest
```

### 3.6 CI Pipeline

Every push to `main` triggers GitHub Actions that run the core validation suite. Results are available as workflow artifacts.  
→ https://github.com/rollroyces/deepcatch/actions

### 3.7 Output Files

After a successful run, `results/` contains simulation-based reports (not clinically validated). See the `results/` directory for available files.

---

## 4. Repository Structure

```
deepcatch/
├── README.md
├── LICENSE                               # MIT
├── CITATION.cff                          # Academic citation metadata
├── requirements.txt
├── RUN_ALL.sh                            # One‑command validation
├── Dockerfile
│
├── src/
│   ├── variant_calling/                  # Bayesian + contrastive DL
│   ├── multimodal_fusion/                # GNN fusion (performance‑weighted)
│   ├── longitudinal/                     # CET/SPRT longitudinal
│   ├── fragmentomics/                    # FragmentoSign (DELFI, MDS, GMM, LOESS)
│   ├── ensemble/                         # MAML meta‑learning
│   └── synthetic_data/                   # Cohort generation
│
├── validation/
│   ├── framework/validation_framework.py # Canonical CV, bootstrap, DeLong
│   ├── py/                               # Python validation (11 modules)
│   ├── node/                             # Node.js validation (20 modules)
│   ├── tcga/                             # TCGA data + validators
│   └── *.py                              # 10 bioinformatics‑grade modules
│
├── results/node/                         # All result JSONs + Markdown reports
├── paper/                                # LaTeX manuscript + 56 references
├── docs/USER_GUIDE.md                    # Full user documentation
├── research/                             # Tier analysis + CET optimisation
└── review/                               # 3 rounds of rigorous peer review
```

---

## 5. Research‑Only Disclaimer

**DeepCatch is research‑stage software.**

- §1 of this README (simulation performance) has been removed — those numbers were based on synthetic data and are not clinically verified.
- §9 contains preliminary validation on 129 real plasma samples (processed frequency data only), but this remains a computational feasibility study, not a clinical assay.
- **Do not use DeepCatch for medical diagnosis, treatment decisions, or any clinical purpose.**
- This software is provided “as is”, without warranty of any kind, under the MIT licence.

---

## 6. Citation

```bibtex
@software{deepcatch2026,
  title        = {{DeepCatch}: Performance‑Weighted Multi‑Modal Fusion for
                   Ultra‑Early Cancer Detection from cfDNA},
  author       = {Royce and DeepCatch Contributors},
  year         = {2026},
  note         = {Simulation framework with preliminary real plasma validation
                  (129 samples, 4-mer motif analysis). DOI to be assigned.},
  url          = {https://github.com/rollroyces/deepcatch},
  version      = {2.1.0},
}
```

---

## 7. Contributing

Contributions are welcome in the following areas:

- **Wet‑lab partnerships** — access to clinical cfDNA samples for real‑world validation
- **GEO/SRA public data** — cross‑validation on published cfDNA datasets
- **Tissue‑of‑origin** — expanding TOO coverage to 20+ cancer types
- **CET flag‑rate optimisation** — CET flag‑rate optimisation
- **Cost modelling** — health‑economics analysis of targeted capture for population screening

Please open an issue to discuss before submitting large pull requests.

---

## 8. Status

| # | Item | Status |
|---|------|--------|
| 1 | Jiang lab 4-mer validation (129 samples) | ✅ Preliminary validation complete |
| 2 | HCC vs Control nested CV AUC | 0.982 |
| 3 | Raw BAM validation | ❌ Pending — need data access agreement |
| 4 | Multi-centre replication | ❌ Pending — design target n=360 |
| 5 | Clinical assay readiness | ❌ Not yet — research only |

---

---

## 9. Real Plasma Validation — 4‑mer End Motif Analysis on Jiang Lab Data

> **v2.1 — First real‑world validation of DeepCatch's CET architecture on actual human plasma cfDNA data.**

### 9.1 Overview

This section describes the first real‑world validation of DeepCatch's CET (Cumulative Evidence Tracking) architecture on actual human plasma cfDNA data from **Professor Jiang Pei‑yong's laboratory at the Chinese University of Hong Kong (CUHK)**. This analysis was performed on processed 4‑mer end‑motif frequency vectors derived from real patient blood draws.

### 9.2 Dataset

| Parameter | Value |
|-----------|-------|
| **Total samples** | 129 plasma samples |
| **Healthy controls** | 38 |
| **Cancer patients** | 91 |
| **Cancer types** | 6 (HCC, lung, HNSCC, CRC, NPC, gastric) |
| **Feature space** | 256 four‑mer end motifs |
| **Data type** | Processed frequency vectors (not raw FASTQ/BAM) |
| **Source** | Jiang lab, CUHK |

### 9.3 Key Results (HCC vs Control, verified via nested CV)

> All numbers below are verified via proper nested cross-validation: feature selection performed within each training fold, eliminating pre-filter leakage.

| Metric | Value |
|--------|-------|
| **Samples** | 72 (34 HCC, 38 Control) |
| **Nested CV AUC** | **0.982** |
| **Bonferroni-significant motifs** | **108 / 256** |
| **FDR-significant motifs** | **164 / 256** |
| **Biological pattern** | CG-rich depletion + AT-rich enrichment |
| **Top motifs** | AAAA (enriched), CCCG (depleted), AAGA (enriched) |

### 9.4 Caveats

1. **HCC only**: Only the HCC vs Control result (n=72) is adequately powered. Other cancer types have n≤17 and their AUC estimates are unreliable.
2. **Processed data**: This analysis uses pre-computed 4-mer frequency vectors, not raw sequencing data.
3. **Single centre**: All samples from one lab (CUHK). Multi-centre replication needed.
4. **Not a clinical assay**: Demonstrates biological signal, not clinical readiness.

### 9.5 Clinical Interpretation Module (New in v2.1)

v2.1 ships with `src/clinical/clinical_interpretation.py` — a `ClinicalReportGenerator` that transforms statistical CET output into clinician‑friendly reports:

```python
from src.clinical import ClinicalReportGenerator

crg = ClinicalReportGenerator(cet_df, fusion_result)
print(crg.generate_briefing())          # One‑paragraph summary
crg.export_json('clinical_report.json')  # Machine‑readable export
with open('report.html', 'w') as f:      # Full HTML report
    f.write(crg.generate_html_report())
```

Use `python run_jiang_analysis.py -i data.xlsx --report` to generate the clinical report alongside the standard summary.

---

*This README was produced by Royce. Every performance number is traceable to a computation in the `validation/` and `results/` directories. No numbers were invented. No clinical claims are intended.* 🧬
