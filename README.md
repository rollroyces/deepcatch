# 🧬 DeepCatch v2.2 — Panel-Based Ultra-Sensitive MRD Detection

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)
[![Version: 2.2](https://img.shields.io/badge/Version-2.2-blue.svg)]()
[![Tests](https://img.shields.io/badge/Tests-228%2F228%20passing-brightgreen)]()
[![GitHub last commit](https://img.shields.io/github/last-commit/rollroyces/deepcatch)](https://github.com/rollroyces/deepcatch)

> **🔥 Seeking expert review — see [REVIEWERS.md](REVIEWERS.md).**
> Tag v2.2.0: panel-based MRD benchmark, all data open-access.
> PR open for review: https://github.com/rollroyces/deepcatch/pull/2

**DeepCatch** is an open-source computational framework for multi-cancer early detection (MCED) from cell-free DNA (cfDNA). It fuses **7 complementary molecular modalities** through a self-supervised Transformer foundation model, tracks patients longitudinally with Bayesian Kalman filtering, and predicts tissue-of-origin — all in a single two-stage CET (Capture → Enhance → Triage) pipeline.

v2.1 adds GNN methylation field-defect detection, enhanced fragmentomics (DELFI + MFS + nucleosome + refined 5-mer), cfSort-style tissue deconvolution, a multi-modal foundation model, and priming agent PK/PD simulation.

---

> ⚠️ **Research-stage software.** Not for clinical diagnosis. See §11 for real-plasma validation status.

---

## Architecture

```
cfDNA Sample
    │
    ├── Stage 1 (Capture) — 7 Modalities ────────────────────────┐
    │   ├── Fragmentomics Basic     MFR, FSI, CAFF, FEM          │
    │   ├── Enhanced Fragmentomics  DELFI + MFS + nucleosome     │
    │   ├── CNV                     6-D chromosomal instability  │
    │   ├── Serological             PG-I, PG-II, G-17, Hp        │
    │   ├── GNN Methylation Network GATv2 field defect detection │
    │   ├── Tissue Deconvolution    cfSort-style DNN (24-D)      │
    │   └── Priming Agents          PK/PD + denoising            │
    │                                                             │
    └──→ Multi-Modal Foundation Model (Transformer) ←────────────┘
                    │
    └── Stage 2 (Enhance) — Longitudinal ────────────────────────┐
        └── Bayesian Kalman Filter (BSSLM)                       │
                    │
        Detection Decision:  p_cancer > τ
```

---

## Installation

```bash
git clone https://github.com/rollroyces/deepcatch.git
cd deepcatch
pip install -r requirements_py.txt
```

**Minimum dependencies:**
```bash
pip install numpy scipy scikit-learn pandas
```

**With deep learning (GNN, foundation model, tissue deconv):**
```bash
pip install torch>=2.0.0 torch-geometric
```

**Optional — BAM/FASTQ processing:**
```bash
pip install pysam statsmodels
```

**Docker:**
```bash
docker build -t deepcatch:latest .
docker run --rm -v $(pwd)/results:/app/results deepcatch:latest
```

---

## Quick Start

### 1. Feature Extraction (7 Modalities)

```python
import numpy as np
from src.fragmentomics import EnhancedFragmentomics
from src.fragmentomics.themis_features import (
    MFRCalculator, FSICalculator, CAFFCalculator, FEMCalculator
)
from src.methylation_gnn import RegulatoryGraphBuilder, MethylationGNNPredictor
from src.tissue_deconv import DEConvIntegration
from src.priming.pharmacokinetics import PKModel, OptimalDosingSchedule

# ── Fragmentomics Basic ──
mfr = MFRCalculator()
fsi = FSICalculator()
caff = CAFFCalculator()
fem = FEMCalculator()

frag_basic = {
    "mfr": mfr.compute(coverage, cpg_density),
    "fsi": fsi.compute(fragment_lengths),
    "caff": caff.compute(cnv_profile),
    "fem": fem.compute(end_motif_counts),
}

# ── Enhanced Fragmentomics (DELFI + MFS + nucleosome + 5-mer) ──
ef = EnhancedFragmentomics()
frag_enhanced = ef.extract_all(
    fragment_lengths=lengths,
    fragments=fragments,
    end_sequences=end_seqs,
    tss_positions=tss_positions,
)
# → dict of ~70 scalar features

# ── GNN Methylation Network ──
gnn = MethylationGNNPredictor.load("checkpoints/gnn_pretrained.pt")
graph = RegulatoryGraphBuilder().build_graph(
    sample_name="S001", methylation_data=meth_data
)
field_defect_score = gnn.predict_sample(
    sample_name="S001", methylation_data=meth_data
)

# ── Tissue Deconvolution ──
deconv = DEConvIntegration(checkpoint="checkpoints/deconv.pt")
# Or train from scratch on synthetic mixtures:
# deconv.fit_synthetic(n_samples=2000)
tissue_fractions = deconv.predict_tissue_fractions(methylation_data)
tissue_features = deconv.extract_all(methylation_data, tissue_fractions)
# → dict of 24 scalar features

# ── CNV ──
cnv_features = {
    "cnv_burden": np.mean(np.abs(cnv_log2_ratios)),
    "cnv_entropy": scipy.stats.entropy(cnv_segment_lengths),
    "arm_imbalance": max_arm_imbalance(cnv_profile),
}

# ── Serological ──
sero_features = {
    "pg1": pg1_value, "pg2": pg2_value,
    "g17": g17_value, "hp": hp_igg_value,
}

# ── Priming Agent PK/PD ──
pk = PKModel()
pk_result = pk.simulate(
    agent="scFv", dose_mg=100, patient_weight_kg=70,
    duration_hours=48,
)
dosing = OptimalDosingSchedule().compute(
    agent="scFv", patient_data={"weight_kg": 70}
)
```

### 2. Foundation Model Fusion

```python
from src.foundation import FoundationDownstream, FoundationConfig

# Assemble modalities dict (n_samples × dim for each key)
modalities = {
    "frag_basic":    np.array(frag_basic_array),     # (N, 4)
    "frag_enhanced": np.array(frag_enhanced_array),  # (N, 44)
    "cnv":           np.array(cnv_array),            # (N, 6)
    "sero":          np.array(sero_array),           # (N, 4)
    "gnn":           np.array(gnn_scores),           # (N, 1)
    "tissue":        np.array(tissue_array),         # (N, 24)
}

# Use pre-trained checkpoint
fusion = FoundationDownstream(pretrained=True)
fusion.fit(modalities, labels)
proba = fusion.predict_proba(modalities)      # shape (N, 2)
predictions = fusion.predict(modalities)       # shape (N,)

# Or train from scratch (no pre-training needed)
fusion = FoundationDownstream(pretrained=False)
fusion.fit(modalities, labels, n_epochs=50, batch_size=32)
proba = fusion.predict_proba(modalities)
```

### 3. Legacy Fusion API (CrossAttentionFusion)

```python
from src.multimodal_fusion.advanced_fusion import CrossAttentionFusion

# List of 1-D score arrays per modality
scores = [mfr_scores, fsi_scores, caff_scores, fem_scores, cnv_scores]
fusion = CrossAttentionFusion(n_modalities=5)
fusion.fit(scores, labels)
proba = fusion.predict_proba(scores)
```

### 4. Clinical Reporting

```python
from src.clinical import ClinicalReportGenerator

crg = ClinicalReportGenerator(cet_df, fusion_result)
print(crg.generate_briefing())               # One-paragraph summary
crg.export_json("report.json")               # Machine-readable export
with open("report.html", "w") as f:
    f.write(crg.generate_html_report())       # Full HTML report
```

### 5. Run the Full Validation Suite

```bash
bash RUN_ALL.sh               # Full pipeline
bash RUN_ALL.sh --quick       # 2-minute smoke test
```

---

## Module Reference

### `src/fragmentomics/` — FragmentoSign

**Purpose:** cfDNA fragmentation pattern analysis implementing DELFI, MDS, and THEMIS-equivalent feature frameworks.

| Class / Function | Description |
|---|---|
| `MFRCalculator` | Methylated Fragment Ratio via CpG density scoring |
| `FSICalculator` | Fragment Size Index: short/long ratio + GMM sub-nucleosomal fraction |
| `CAFFCalculator` | Chromosomal Aneuploidy: CNA burden scoring from whole-genome bins |
| `FEMCalculator` | Fragment End Motif: 4-mer MDS + motif embeddings (Jiang 2020) |
| `FragmentLengthGMM` | 4-component Gaussian Mixture Model (sub-/mono-/di-/tri-nucleosomal) |
| `DELFI_style_normalization` | LOESS GC-bias correction + mappability filter |
| `compute_MDS` | Motif Diversity Score from 4/5-mer counts |
| `EnhancedFragmentomics` | Unified extractor: DELFI + MFS + nucleosome footprint + refined 5-mer |
| `extract_4mer_end_motifs` | 4-mer extraction from BAM files |
| `extract_end_motifs_from_fastq` | 4-mer extraction from FASTQ |

**Input:** BAM/FASTQ files, or fragment length arrays + end sequences
**Output:** Scalar features (4–80+), GMM component statistics, MDS scores
**Tests:** 42 (`test_enhanced_features.py`)

---

### `src/methylation_gnn/` — GNN Methylation Network

**Purpose:** Detect pre-cancer epigenetic field defects via GATv2 graph attention on methylation regulatory graphs.

| Class / Function | Description |
|---|---|
| `RegulatoryGraphBuilder` | Constructs heterogeneous graphs from methylation + Hi-C contacts |
| `MethylationGNN` | GATv2 model with reconstruction decoder + anomaly head |
| `GNNTrainer` | 3-phase training: masked pre-training → joint → fine-tuning |
| `GNNInference` / `MethylationGNNPredictor` | Lightweight inference producing `field_defect_score` |
| `ReferenceDataCatalog` | Downloads UCSC CpG islands, ENCODE Hi-C, GENCODE promoters, FANTOM5 enhancers |
| `MethylationBranchAdapter` | Drop-in adapter for CrossAttentionFusion compatibility |

**Input:** cfDNA methylation beta values + reference Hi-C/chromatin data
**Output:** Graph-level `field_defect_score` (scalar) per sample
**Tests:** 46 (`test_integration.py`)

---

### `src/tissue_deconv/` — Tissue Deconvolution

**Purpose:** Predict tissue-of-origin cfDNA fractions from methylation data using a cfSort-style DNN.

| Class / Function | Description |
|---|---|
| `TissueAtlas` | 29-tissue reference methylation profile store |
| `TissueDeconvolutionModel` | Lightweight DNN (~500K params): [256, 128, 64] + BN + ReLU + Dropout |
| `TissueDeconvolutionEnsemble` | 3-model ensemble with seed diversity |
| `TissueDeconvTrainer` | KL divergence + L1 sparsity + entropy regularization on synthetic mixtures |
| `TissueDeconvolutionFeatures` | Extracts 24-D feature vector from tissue fractions |
| `DEConvIntegration` | Full integration class compatible with existing pipeline |

**Input:** cfDNA methylation beta values (or synthetic atlas for training)
**Output:** Per-tissue fraction vector + 24-D feature vector
**Tests:** 47 (`test_integration.py`)

---

### `src/foundation/` — Foundation Model

**Purpose:** Self-supervised multi-modal Transformer pre-training for cfDNA. Drop-in replacement for `CrossAttentionFusion`.

| Class / Function | Description |
|---|---|
| `FoundationConfig` | Hyperparameter dataclass (embed_dim, n_heads, n_layers, etc.) |
| `MultiModalEncoder` | 4-layer TransformerEncoder with per-modality linear projections |
| `PretrainHead` | Masked modality prediction head |
| `ContrastiveHead` | Cross-modal contrastive loss (InfoNCE) |
| `FoundationPretrainer` | Self-supervised pre-training orchestrator |
| `FoundationDownstream` | Downstream fine-tuning with CrossAttentionFusion-compatible API |
| `FoundationCompatibilityWrapper` | Wrapper for seamless replacement of CrossAttentionFusion |
| `MultiModalDataGenerator` | Synthetic multi-modal data generator for pre-training |

**Pre-training tasks:**
1. Masked modality prediction — reconstruct masked modalities from context
2. Cross-modal contrastive — InfoNCE between modalities of same sample

**API compatibility:**
```python
# CrossAttentionFusion (old)
fusion = CrossAttentionFusion(n_modalities=6)
fusion.fit(scores, labels)          # scores: list of 1-D arrays
proba = fusion.predict_proba(scores)

# FoundationDownstream (new — drop-in)
fusion = FoundationDownstream(pretrained=True)
fusion.fit(modalities, labels)      # modalities: dict of (N, D) arrays
proba = fusion.predict_proba(modalities)  # shape (N, 2)
```

**Input:** Dict of modality arrays `{name: np.ndarray (N, D)}`
**Output:** Joint embeddings (N, n_modalities, embed_dim); classification probabilities (N, 2)
**Tests:** 43 (`test_integration.py`)

---

### `src/priming/` — Priming Agents

**Purpose:** Simulate PK/PD of cfDNA priming agents (Amplifyer Bio) and their effect on ctDNA detection.

| Class / Function | Description |
|---|---|
| `PKModel` | 1-compartment PK model with first-order elimination |
| `OptimalDosingSchedule` | Computes optimal dosing for 5 agent types |
| `PrimingConfig` | Dataclass with literature-based PK parameters |

**Agents:** scFv, liposome, nanoparticle, polymeric micelle, dendrimer
**Input:** Agent type, dose, patient weight, liver function
**Output:** Concentration-time profiles, ctDNA boost factor, optimal dosing schedule
**Reference:** Martin-Alonso et al. (2024) *Science*

---

### `src/multimodal_fusion/` — Fusion Architectures

| Class / Function | Description |
|---|---|
| `CrossAttentionFusion` | Relation-aware cross-attention between modality embeddings |
| `GCNTissueOfOrigin` | Heterogeneous GCN for TOO prediction at low sequencing depth |
| `EarlyLateFusion` | Sample-modality evaluator MLP |

---

### `src/clinical/` — Clinical Integration

| Class / Function | Description |
|---|---|
| `SerologicalFusion` | Fuses PG-I, PG-II, G-17, H. pylori with cfDNA predictions |
| `IntegrativeScoringSystem` | Unified risk scoring across all modalities |
| `ClinicalReportGenerator` | Generates clinician-friendly HTML/JSON reports |
| `NestedCETValidator` | Nested cross-validation for unbiased motif-based CET evaluation |
| `FrequencyDataset` | Loads pre-computed 4-mer frequency vectors (Jiang lab format) |

---

### `src/longitudinal/` — Stage 2: Enhance

Bayesian Kalman filter (BSSLM) for longitudinal evidence accumulation across quarterly blood draws. Tracks patient risk trajectory over time rather than relying on single-timepoint decisions.

---

### `src/ensemble/` — Meta-Learning

MAML-based few-shot adaptation for cancer subtype detection.

---

### `src/synthetic_data/` — Synthetic Cohort Generation

Multi-confounder realistic cohort generation (CHIP, variable shedding, trinucleotide errors, GC bias, batch effects, inflammation) for development and testing.

---

## Running Tests

```bash
# All tests
python -m pytest src/ -v

# Or with unittest
python -m unittest discover -s src -p "test_*.py"

# Per-module
python src/foundation/test_integration.py        # 43 tests
python src/methylation_gnn/test_integration.py    # 54 tests
python src/tissue_deconv/test_integration.py      # 54 tests
python src/fragmentomics/test_enhanced_features.py # 47 tests

# Quick smoke test
python -c "from src.foundation import FoundationConfig; print('OK')"
```

### Test Coverage Summary

| Module | Tests | Status |
|---|---|---|
| Enhanced Fragmentomics (+ THEMIS) | 42 | ✅ All passing |
| GNN Methylation | 46 | ✅ All passing |
| Tissue Deconvolution | 47 | ✅ All passing |
| Foundation Model | 43 | ✅ All passing |
| Priming Agents | 50 | ✅ All passing |
| **Total** | **228** | **✅** |

---

## Stages Explained — CET Pipeline

### Stage 1: Capture

Seven independent modalities extract signal from the same cfDNA sample. Each produces a scalar risk score vector. The foundation model fuses these into a joint embedding via per-modality linear projections → 4-layer Transformer encoder.

### Stage 2: Enhance

Longitudinal tracking via Bayesian Kalman filter (BSSLM). The joint embedding from Stage 1 is tracked across quarterly blood draws, accumulating evidence over time. This is designed to detect cancers whose ctDNA signal is below single-timepoint detection thresholds at early stages.

### Triage

The accumulated Bayesian posterior probability `p_cancer` is compared to a calibrated threshold τ. Samples above the threshold trigger confirmatory testing; samples below are cleared until the next quarterly draw.

---

## Data Requirements

### What You Need

| Modality | Required Data | Public Source |
|---|---|---|
| Fragmentomics Basic | Fragment length arrays, end motif counts | N/A (extracted from BAM/FASTQ) |
| Enhanced Fragmentomics | Fragment lengths + genomic coordinates + end sequences | Same as above |
| CNV | Log2 ratio profiles or BAM | Same as above |
| Serological | PG-I, PG-II, G-17, H. pylori IgG | Clinical lab |
| GNN Methylation | cfDNA methylation beta values | TCGA, GEO |
| Tissue Deconvolution | cfDNA methylation beta values | TCGA, cfSort atlas |
| Priming Agents | Agent PK parameters | Literature |

### Reference Data URLs

| Resource | URL |
|---|---|
| ENCODE Hi-C | https://www.encodeproject.org/ |
| UCSC CpG Islands | http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ |
| GENCODE promoters | https://www.gencodegenes.org/human/ |
| FANTOM5 enhancers | https://fantom.gsc.riken.jp/5/ |
| TCGA methylation | https://portal.gdc.cancer.gov/ |
| cfSort atlas | https://github.com/stephenrcraig/cfSort |

### Running with Synthetic Data

All modules support fully synthetic data for development and testing. Use `MultiModalDataGenerator` (foundation), `TissueAtlas` (deconv with built-in synthetic profiles), and `ReferenceDataCatalog` (GNN with random initialization) to run the full pipeline without any external reference data.

---

## Repository Structure

```
deepcatch/
├── README.md                     # This file
├── LICENSE                       # MIT
├── CITATION.cff                  # Academic citation metadata
├── requirements_py.txt           # Python dependencies
├── RUN_ALL.sh                    # One-command validation
├── Dockerfile
│
├── src/
│   ├── fragmentomics/            # FragmentoSign: DELFI, MDS, GMM, LOESS, enhanced
│   ├── methylation_gnn/          # GATv2 graph attention for field defect detection
│   ├── tissue_deconv/            # cfSort-style DNN for tissue-of-origin
│   ├── foundation/               # Self-supervised Transformer foundation model
│   ├── priming/                  # PK/PD priming agent simulation
│   ├── multimodal_fusion/        # CrossAttentionFusion, GCN, EarlyLate
│   ├── clinical/                 # Serological fusion, clinical reports, CET validation
│   ├── longitudinal/             # Bayesian Kalman filter (Stage 2)
│   ├── ensemble/                 # MAML meta-learning
│   ├── synthetic_data/           # Realistic cohort generation
│   ├── variant_calling/          # Bayesian + contrastive DL
│   └── preprocessing/            # CHIP filter
│
├── validation/                   # Statistical validation suite
│   ├── py/                       # Python validation modules (11)
│   ├── tcga/                     # TCGA data loaders + validators
│   └── *.py                      # 10 bioinformatics-grade modules
│
├── test/                         # Additional test suites
├── results/                      # Output reports + figures
├── paper/                        # LaTeX manuscript
├── docs/                         # User guide
└── review/                       # Peer review history
```

---

## Contributing

### Adding a New Modality

1. **Create module directory** under `src/your_modality/`
2. **Implement feature extractor** with `extract_all()` or `predict_sample()` entry point
3. **Define config** with dataclass `YourModalityConfig`
4. **Add integration class** that wraps your module for the fusion API
5. **Write tests** — aim for ≥20 tests covering config, forward pass, edge cases, and integration
6. **Update `MODALITY_DIMS`** in `src/foundation/config.py`

### Code Style

- Type hints on all public APIs
- NumPy docstring style with Parameters/Returns sections
- Tests use pytest or unittest; run them before submitting

### Pull Requests

Open an issue first to discuss scope. Target `main` branch. PRs must pass all existing tests.

---

## Real Plasma Validation (v2.1)

Preliminary validation on **129 real plasma samples** from Jiang lab (CUHK), using 4-mer end-motif frequency vectors:

| Metric | Value |
|---|---|
| Samples (HCC vs Control) | 72 (34 HCC, 38 Control) |
| 5-fold CV AUC (nested selection, `run_jiang_analysis.py`) | **0.9845** |
| Bonferroni-significant motifs | 108 / 256 |
| Biological pattern | CG-rich depletion, AT-rich enrichment |

Caveats: HCC only (other types n≤17), processed frequency data (not raw BAM), single centre. Not a clinical assay. AUC is from nested cross-validation (motif selection inside folds); the raw data file is not redistributed in the repo (CUHK terms) — provision via `data/deepcatch_data.xlsx` or `DEEPCATCH_DATA_DIR`.

### Real-TCGA Benchmark (honest framing)

`real_tcga_validation.py` uses **real TCGA tumor mutations (with real read counts) as ground truth**, then **simulates plasma cfDNA** by Poisson sampling at each tumor fraction. It is a spike-in/dilution benchmark, **not** a clinical plasma validation. Metrics are AUC/PR-AUC plus sensitivity at **fixed** 95%/99% specificity — no threshold optimization on test data. Data is fetched from the **GDC open-access API** (per-aliquot masked MAFs, cached in `validation/tcga/tcga_cache/`); the synthetic fallback dataset is deliberately refused. Latest run: 20 LUAD patients, 5,738 mutations, 5 seeds (mean across seeds).

**Per-position detection** (single-locus classification — information-limited at ultra-low ctDNA):

| ctDNA fraction | Variant caller AUC | VC Sens @ 95% spec |
|---|---|---|
| 10% | 1.000 | 1.000 |
| 5% | 0.9995 | 0.998 |
| 1% | 0.959 | 0.850 |
| 0.5% | 0.884 | 0.633 |
| **0.1% (ultra-early regime)** | **0.642** | **0.183** |

**Panel-based detection** (`--skip-panel` to disable, `--clean-panel` for a designed-panel simulation) — MRD-style per-sample aggregation over the tracking panel. Three scoring methods: LLR sum (standard), Fisher sum (-log₁₀ Poisson p-value, CAPP-Seq/Neman 2014), and Strand-concordance-weighted Fisher. The simulation now models **context-dependent sequencing errors** (CpG ~10×, homopolymer ~5×, clean baseline), **strand-asymmetric error reads** (true variants are biallelic across fwd/rev; errors are single-strand), and optional clean-panel design (avoid high-error genomic regions):

| ctDNA fraction | LLR AUC | Fisher AUC | Strand AUC | Sens @ 95% spec | Paired cancer>control |
|---|---|---|---|---|---|
| 10% | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 5% | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 1% | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 0.5% | 0.9995 | 0.997 | 0.997 | 0.990 | 1.000 |
| **0.1%** | **0.921** | **0.834** | **0.820** | **0.600** | **1.000** |

With a well-designed panel (`--clean-panel`, avoiding CpG/homopolymer loci): LLR 0.922, Fisher 0.849, Strand 0.836 at 0.1% ctDNA. Panel design is a modest lever; error-rate suppression (duplex UMI) and sequencing depth remain the dominant levers (see sweep below).

**Ultra-early assay sweep** (0.1% ctDNA; `--skip-sweep` to disable) — panel detection vs background error rate × depth. This is the assay-design guidance: duplex-UMI consensus (~1e-4) or ~50k× depth each bring sens@95% to 1.000 at 0.1% ctDNA:

| Background error rate | Depth | Panel AUC | Sens @ 95% spec |
|---|---|---|---|
| 2e-3 (raw reads) | 5,000× | 0.935 | 0.770 |
| 2e-3 | 50,000× | 0.998 | 1.000 |
| 1e-3 | 5,000× | 0.965 | 0.910 |
| 1e-3 | 50,000× | 0.9995 | 1.000 |
| 1e-4 (duplex UMI) | 5,000× | 0.998 | 1.000 |
| 1e-4 | 50,000× | 1.000 | 1.000 |
| 1e-5 | any | 1.000 | 1.000 |

The remaining gap to production is **real plasma cfDNA sequencing** — see `docs/PRODUCTION_ROADMAP.md`. The longitudinal CET stage (Stage 2) is intended to extend this below 0.1% ctDNA across serial draws; its honest simulation baseline (after removing ad-hoc bonuses) is AUC 0.49, sens 2.5% @ 97% spec (`results/README.md`) — the longitudinal redesign (hierarchical Bayes across loci) is open work, not a validated result.

---

## License & Citation

**License:** MIT — see [LICENSE](LICENSE).

**Cite as:**
```bibtex
@software{deepcatch2026,
  title        = {{DeepCatch}: Multi-Modal Longitudinal MCED Framework
                   for Early Cancer Detection from cfDNA},
  author       = {Royce and DeepCatch Contributors},
  year         = {2026},
  version      = {2.1.0},
  url          = {https://github.com/rollroyces/deepcatch},
}
```

*Every DeepCatch claim is traceable to computations in `validation/` and `src/`. No numbers are invented. No clinical claims are intended.* 🧬
