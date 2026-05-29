# 🏗️ Pipeline Architecture

DeepCatch is built around a **two‑stage CET (Cumulative Evidence Tracking)** architecture with **performance‑weighted multi‑modal fusion**. This page explains each component in detail and provides a complete CLI reference for `run_jiang_analysis.py`.

---

## Table of Contents

1. [Two‑Stage CET Architecture](#two-stage-cet-architecture)
2. [Performance‑Weighted Multi‑Modal Fusion](#performance-weighted-multi-modal-fusion)
3. [FragmentoSign Subsystem](#fragmentosign-subsystem)
4. [Feature Engineering](#feature-engineering)
5. [Nested Cross‑Validation](#nested-cross-validation)
6. [CLI Reference](#cli-reference)

---

## Two‑Stage CET Architecture

The core innovation of DeepCatch is a **longitudinal screening architecture** designed to catch cancer signals that are too weak for a single blood draw.

### Stage 1 — Permissive CET (Enrichment)

Each modality is tested independently per timepoint using non‑parametric tests:

| Test | When Used | Purpose |
|------|-----------|---------|
| **Mann‑Whitney U** | Per‑motif frequency comparison | Rank‑based enrichment — robust to non‑normal distributions |
| **Cliff's delta** | Effect size | Direction‑aware magnitude (range: −1 to +1) |
| **Benjamini‑Hochberg FDR** | Multiple testing correction | Controls false discovery rate at α = 0.05 |

The output is a ranked list of features by **composite score**:

```text
composite_score = −log₁₀(p_value) × |Cliff's delta|
```

Features with p ≥ α receive score = 0 and are excluded from fusion.

> **Design rationale:** Stage 1 is intentionally permissive — we want high sensitivity. False positives are managed in Stage 2.

### Stage 2 — Confirmatory Fusion

Only the top‑k ranked features are fed into **L2‑regularised logistic regression** (C=10.0, liblinear solver, max 5000 iterations):

```text
P(cancer | features) = σ(β₀ + β₁·motif₁ + β₂·motif₂ + ... + βₖ·motifₖ)
```

- **Why C=10.0?** Weaker regularisation preserves more discriminative information from correlated motif features — critical for high‑signal data like 4‑mer frequencies.
- **Why logistic regression?** Interpretable coefficients, calibrated probabilities, and well‑behaved with small‑to‑moderate sample sizes.

### Optimal‑k Selection (Elbow Method)

When `--optimal-k` is enabled, the pipeline finds the knee point of the cumulative composite score curve. This avoids both under‑fitting (too few motifs) and over‑fitting (noisy, low‑rank motifs diluting the signal).

---

## Performance‑Weighted Multi‑Modal Fusion

*(Conceptual — the 4‑mer pipeline uses a single modality. Multi‑modal fusion is the full DeepCatch framework.)*

In the full multi‑modal pipeline, five molecular signals are measured from the same blood draw:

1. **Somatic mutations** — ctDNA variant allele fractions
2. **DNA methylation** — entropy‑based CpG density scoring
3. **Fragmentomics** — fragment size + end‑motif patterns (FragmentoSign)
4. **Copy‑number alterations** — CNA burden scoring
5. **CTC count** — circulating tumour cell estimate

Each modality produces a probability score. Instead of averaging equally (Bie et al. 2023 approach), DeepCatch **weights each modality by its individual AUC** measured on a held‑out validation set. Modalities with AUC < 0.5 receive zero weight and are excluded.

```text
Fusion Score = Σᵢ (wᵢ × pᵢ),   where wᵢ = max(0, AUCᵢ − 0.5)
```

---

## FragmentoSign Subsystem

DeepCatch's fragmentomics engine — **FragmentoSign** — implements three key frameworks:

| Component | Method | Reference |
|-----------|--------|-----------|
| **GC‑bias correction** | LOESS local normalisation | Cristiano 2019, *Nature* |
| **Fragment length model** | 4‑component Gaussian Mixture Model (GMM) | Snyder 2016, *Cell* |
| **End‑motif analysis** | 4‑mer extraction + MDS (Motif Diversity Score) | Jiang 2020, *Cancer Discovery* |

### GMM Components

| Component | Mean Size | Biological Origin |
|-----------|-----------|-------------------|
| Sub‑nucleosomal | ~80 bp | Tumour‑derived short fragments |
| Mono‑nucleosomal | ~167 bp | Normal nucleosome‑wrapped DNA |
| Di‑nucleosomal | ~334 bp | Two‑nucleosome fragments |
| Tri‑nucleosomal | ~501 bp | Three‑nucleosome fragments |

The sub‑nucleosomal fraction is the key cancer signal — tumour‑derived cfDNA is enriched for fragments shorter than 150 bp.

### MDS (Motif Diversity Score)

The MDS quantifies end‑motif diversity across 256 possible 4‑mers. Cancer plasma typically shows **reduced MDS** — a restricted set of cleavage patterns — compared to healthy controls. This is computed as:

```text
MDS = −Σᵢ fᵢ × log₂(fᵢ)    (Shannon entropy of 4‑mer frequencies)
```

---

## Feature Engineering

The `run_jiang_analysis.py` pipeline applies two transformations by default:

### 1. Rank Transformation (`--rank-features`)

Each sample's 256 motif frequencies are converted to within‑sample ranks:

```python
# For each sample: replace frequency with rank (1 = lowest, 256 = highest)
ranked[j] = argsort(argsort(frequencies[j]))
```

**Why?** Ranks are robust to outliers, scale‑invariant, and capture relative motif abundance patterns rather than absolute frequencies — which can vary due to technical factors (library prep, sequencing depth).

### 2. CG/AT Composition Ratio Features (`--ratio-features`)

Additional features derived from nucleotide composition:

| Feature | Formula | Biological Meaning |
|---------|---------|-------------------|
| `cg_ratio` | Σ(CG motifs) / Σ(all motifs) | Global CpG content — lower in cancer |
| `at_ratio` | Σ(AT motifs) / Σ(all motifs) | AT‑rich motif abundance — higher in cancer |
| `cg_at_ratio` | cg_ratio / at_ratio | Combined index — most discriminative single feature |

Both transformations can be disabled:
```bash
python run_jiang_analysis.py -i data.xlsx --no-rank-features --no-ratio-features
```

---

## Nested Cross‑Validation

Standard CV (cross‑validation) can produce **optimistically biased** AUC estimates when feature selection and model fitting share the same data. Nested CV solves this by splitting into **outer** and **inner** loops:

```
┌────────────────────────────────────────────────────────┐
│                   Nested Cross-Validation               │
│                                                         │
│  Outer Loop (5 folds):                                  │
│  ┌───────────────────────────────────────────────────┐ │
│  │ Fold 1: [Train] [Train] [Train] [Train] │ Test   │ │
│  │   Inner Loop (3 folds):                            │ │
│  │   ┌─────────────────────────────────────┐         │ │
│  │   │ [Inner Train] │ [Inner Val]         │         │ │
│  │   │ → Select top-k + tune C in inner CV │         │ │
│  │   │ → Fit final model on all inner data  │         │ │
│  │   │ → Evaluate on held-out outer Test    │         │ │
│  │   └─────────────────────────────────────┘         │ │
│  └───────────────────────────────────────────────────┘ │
│  ... repeat for all 5 outer folds ...                   │
│                                                         │
│  Final AUC = mean(outer fold AUCs) ± std                 │
└────────────────────────────────────────────────────────┘
```

**Key guarantee:** Feature selection (k motifs) and hyperparameter tuning (C) happen entirely within the inner loop for each outer fold — the outer test set is never seen during training.

Enable nested CV:
```bash
python run_jiang_analysis.py -i data.xlsx --nested-cv
```

---

## CLI Reference — `run_jiang_analysis.py`

### Required

| Flag | Short | Type | Description |
|------|-------|------|-------------|
| `--input` | `-i` | path | Frequency data file (`.xlsx`, `.csv`, `.npy`, `.npz`) |

### Cancer Type Selection

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--cancer-type` | | string | `None` | Cancer type column value (e.g., `HCC`, `LC`) |
| `--control-label` | | string | `None` | Label for control group (auto‑detected if binary) |
| `--labels` | `-l` | path | `None` | Separate label file |

### Analysis Parameters

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--top-k` | | int | `50` | Number of top motifs for logistic regression fusion |
| `--alpha` | | float | `0.05` | Significance threshold for FDR correction |
| `--lr-C` | | float | `10.0` | LogisticRegression C (inverse regularisation strength) |
| `--select-by` | | choice | `p_value` | Feature selection method: `p_value` or `variance` |
| `--seed` | | int | `42` | Random seed for reproducibility |

### Feature Engineering

| Flag | Description |
|------|-------------|
| `--rank-features` | Convert motif frequencies to per‑sample ranks (default: on) |
| `--no-rank-features` | Use raw frequencies instead |
| `--ratio-features` | Add CG/AT composition ratio features (default: on) |
| `--no-ratio-features` | Disable composition ratio features |

### Validation & Output

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--nested-cv` | | flag | `False` | Enable nested cross‑validation (5 outer × 3 inner folds) |
| `--optimal-k` | | flag | `False` | Find optimal k via elbow method |
| `--plot` | `-p` | flag | `False` | Generate 4 plots (volcano, heatmap, ROC, importance) |
| `--report` | `-r` | flag | `False` | Generate clinical interpretation report (HTML + JSON) |
| `--output` | `-o` | path | `results/jiang_reanalysis/` | Output directory |
| `--verbose` | `-v` | flag | `False` | Enable debug logging |

### Example Commands

```bash
# Minimal: HCC vs Control, default parameters
python run_jiang_analysis.py -i deepcatch_data.xlsx --cancer-type HCC

# Full analysis: nested CV + optimal-k + plots + report
python run_jiang_analysis.py -i deepcatch_data.xlsx \
  --cancer-type HCC --control-label Control \
  --top-k 50 --nested-cv --optimal-k --plot --report

# Run all pairwise comparisons (no --cancer-type)
python run_jiang_analysis.py -i deepcatch_data.xlsx --top-k 30

# Custom output directory
python run_jiang_analysis.py -i data.csv -o my_results/ --plot

# Verbose mode for debugging
python run_jiang_analysis.py -i data.xlsx -v --top-k 20
```

---

## Related Pages

- 🏠 **[Home →](Home)**
- 🚀 **[Getting Started →](Getting-Started)**
- 🧪 **[Jiang 4‑mer Validation →](Jiang-4mer-Validation)**
