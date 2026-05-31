# DeepCatch v2.1 — GNN Methylation Network Extension

## Overview

Added a Graph Neural Network (GNN) branch to DeepCatch Stage 1 (Capture) pipeline for detecting pre-cancer epigenetic field defects through methylation network reconstruction.

**Date:** 2026-05-31  
**Author:** Research Sub-agent (implementation) / Royce (design)  
**Module:** `src/methylation_gnn/`

## Motivation

Current DeepCatch v2.0 detects cancer through cfDNA fragmentation patterns, but has a fundamental detection floor at tumor fractions below 0.01%. Pre-cancer epigenetic field defects — hypermethylation at tumor suppressor CpG islands, global hypomethylation, enhancer dysregulation — manifest before any cancer cell exists and are detectable through methylation network analysis.

A GNN is the natural architecture because:
1. Methylation is inherently relational (3D chromatin topology)
2. Hi-C data provides ground-truth physical interaction edges
3. Field defects manifest at the subnetwork level
4. Graph anomaly detection identifies coordinated dysregulation

## New Files

```
src/methylation_gnn/
├── __init__.py              # Module exports (17 public names)
├── config.py                # Hyperparameters: GNNConfig, NODE_FEATURE_SPEC
├── graph_builder.py         # Graph construction: RegulatoryGraphBuilder
├── gnn_model.py             # GATv2Conv + dual head: MethylationGNN
├── gnn_trainer.py           # 3-phase training: GNNTrainer
├── gnn_inference.py         # Inference pipeline: GNNInference, MethylationGNNPredictor
├── data.py                  # Reference data URLs + preprocessing
├── integration.py           # Fusion adapter: ModularArmsBuilder
└── test_integration.py      # Comprehensive test suite
```

## Architecture

```
Methylation GNN Architecture:
  Node Features (N × 20)
    → Feature Projection (20→64)
    → GATv2Conv × 3 layers (64→128→256)
    → Dual Head:
        ├── Reconstruction Decoder (256→128→64→20)
        └── Anomaly Head (256→128→64→1) → field_defect_score

Training:
  Phase 1: Self-supervised masked node prediction (30% nodes)
  Phase 2: Joint reconstruction + anomaly (supervised)
```

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| **GATv2Conv** over GCN/GAT | Dynamic attention computes attention after linear transform; more expressive for biological graphs |
| **Dual head** (reconstruction + anomaly) | Reconstruction learns normal network; anomaly learns diagnostic significance |
| **3 layers** | Sufficient receptive field for 50K-node regulatory graphs |
| **Heterogeneous edges** (5 types) | Hi-C, co-methylation, co-fragmentation, proximity, regulatory domain |
| **50K nodes** | Covers CpG islands + enhancers + promoters; tractable on most GPUs |
| **~2M parameters** | Fits on Apple Silicon MPS and 16GB GPU |

## Integration

The GNN branch produces a `field_defect_score ∈ [0,1]` that feeds into the existing `CrossAttentionFusion` as a **5th modality**:

```python
# v2.0: 4 modalities
fusion = CrossAttentionFusion(n_modalities=4)
fusion.fit([frag, cnv, sero, mfr], labels)

# v2.1: 5 modalities with GNN
from src.methylation_gnn import extend_fusion_with_gnn
gnn_scores = adapter.predict_batch(samples)
all_scores = extend_fusion_with_gnn([frag, cnv, sero, mfr], gnn_scores)
fusion = CrossAttentionFusion(n_modalities=5)
fusion.fit(all_scores, labels)
```

## Files Modified

- `src/methylation_gnn/__init__.py` — New module (created)
- `src/methylation_gnn/config.py` — New (GNNConfig dataclass with validation)
- `src/methylation_gnn/graph_builder.py` — New (RegulatoryGraphBuilder)
- `src/methylation_gnn/gnn_model.py` — New (MethylationGNN, FieldDefectLoss)
- `src/methylation_gnn/gnn_trainer.py` — New (GNNTrainer with 3-phase training)
- `src/methylation_gnn/gnn_inference.py` — New (GNNInference, MethylationGNNPredictor)
- `src/methylation_gnn/data.py` — New (ReferenceDataCatalog, preprocessing)
- `src/methylation_gnn/integration.py` — New (ModularArmsBuilder, MethylationBranchAdapter)
- `src/methylation_gnn/test_integration.py` — New (comprehensive test suite)

## Dependencies Added

```txt
torch >= 2.0.0          # PyTorch with MPS/CUDA support
torch_geometric >= 2.5.0 # PyG for graph neural networks
```

## Test Coverage

`test_integration.py` covers:
- **Config:** defaults, validation, serialization, node feature spec
- **Graph Builder:** synthetic region creation, BED parsing, node features (all sources, missing data), edge fallback, PyG graph construction
- **GNN Model:** creation, forward pass, pretrain forward, predict, homogeneous fallback, loss function, smallest model
- **Trainer:** creation, pretrain epoch, full training cycle, checkpoint save/load, predict
- **Inference:** load, predict, predict_details, high-level predictor
- **Integration:** extend_fusion_with_gnn, ModularArmsBuilder (with/without GNN), fusion compatibility
- **Data:** reference catalog, download URLs, synthetic data generation, methylation preprocessing, co-methylation matrix, BED parsing, CpG region generation

Run: `python src/methylation_gnn/test_integration.py`

## Remaining Work (Data Acquisition)

The following are NOT part of this PR — they require network access and large data downloads:

1. Download UCSC CpG island track for hg38
2. Download GENCODE v44 promoter annotations (GTF)
3. Download FANTOM5 enhancer atlas (BED)
4. Download GM12878 Hi-C from Rao et al. 2014 (GSE63525)
5. Download ENCODE chromatin marks (bigWig)
6. Download TCGA methylation β-value matrices (Illumina 450K/850K)

All URLs and instructions are catalogued in `data.py::ReferenceDataCatalog`.

## Rollback

To revert to v2.0 behavior: simply don't import from `src.methylation_gnn`, and use `n_modalities=4` in fusion layers. The new module has no side effects on existing code.
