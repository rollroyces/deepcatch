# 🚀 Getting Started with DeepCatch

This guide walks you through installing DeepCatch and running your first 4‑mer CET analysis in under 5 minutes.

---

## 📋 Prerequisites

| Requirement | Minimum Version | Notes |
|-------------|----------------|-------|
| **Python** | 3.9+ | Recommended: 3.10 or 3.11 |
| **pip** | 21.0+ | Comes with Python ≥ 3.9 |
| **Git** | 2.0+ | For cloning the repo |
| **Disk space** | ~ 1 GB | For dependencies + sample data |

Core Python dependencies (installed automatically via `requirements.txt`):

```text
numpy >= 1.24.0
scipy >= 1.10.0
scikit-learn >= 1.3.0
matplotlib >= 3.7.0
seaborn >= 0.12.0
pandas >= 2.0.0
torch >= 2.0.0          # Deep learning models (optional for basic CET)
statsmodels             # LOESS normalisation (optional)
```

---

## 💿 Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/rollroyces/deepcatch.git
cd deepcatch
```

### Step 2: Install Dependencies

```bash
pip install -r requirements.txt
```

> **Note:** If `requirements.txt` is not present (e.g., early release), install the core packages manually:
> ```bash
> pip install numpy scipy pandas scikit-learn matplotlib seaborn
> ```

### Step 3: Verify Installation

```bash
python -c "from src.clinical.frequency_input import FrequencyDataset; print('✅ DeepCatch ready.')"
```

Expected output:

```
✅ DeepCatch ready.
```

---

## ⚡ Quick Demo — Jiang 4‑mer HCC Analysis

The fastest way to see DeepCatch in action is the Jiang 4‑mer CET pipeline. This runs the full analysis on 129 real plasma samples:

```bash
python run_jiang_analysis.py \
  -i results/prof_jiang_4mer_analysis/deepcatch_data.xlsx \
  --cancer-type HCC \
  --control-label Control \
  --top-k 50 \
  --nested-cv \
  --optimal-k \
  --plot \
  --report
```

**What this does (in order):**

1. Loads 256 4‑mer motif frequencies from the Excel file
2. Filters to HCC vs Control (72 samples: 34 HCC, 38 Control)
3. Runs Mann‑Whitney U per motif (CET Step 2)
4. Computes Cliff's delta effect sizes
5. Applies Benjamini‑Hochberg FDR correction
6. Finds the optimal number of top motifs (elbow method)
7. Logistic regression fusion on top‑k motifs
8. Nested cross‑validation (5×3 folds) for unbiased AUC
9. Generates 4 plots (volcano, heatmap, ROC, feature importance)
10. Produces a clinical interpretation report (HTML + JSON)

**Runtime:** ~30 seconds on a modern laptop.

---

## 📂 Expected Output Files

After a successful run, `results/jiang_reanalysis/` (default) will contain:

```
results/jiang_reanalysis/
├── cet_motif_results.csv         # Per‑motif p‑values, effect sizes, scores
├── fusion_coefficients.csv       # Logistic regression coefficients
├── nested_cv_results.csv         # Nested CV AUC with 95% CI
├── summary_report.md             # Full Markdown summary report
├── clinical_report.html          # Clinician‑friendly HTML report (--report)
├── clinical_report.json          # Machine‑readable export (--report)
├── clinical_briefing.txt         # One‑paragraph text summary (--report)
└── plots/
    ├── volcano.png               # Volcano plot (‑log₁₀(p) vs effect size)
    ├── heatmap.png               # Top‑30 motif heatmap by sample
    ├── roc_curve.png             # ROC curve with AUC
    └── feature_importance.png    # Top‑20 motif coefficients
```

---

## 🔧 Troubleshooting

| Problem | Likely Cause | Solution |
|---------|-------------|----------|
| `ModuleNotFoundError: No module named 'sklearn'` | Missing dependency | `pip install scikit-learn` |
| `ImportError: src.clinical.cet_cross_validator` | Module not on path | Run from repo root: `cd deepcatch` |
| `ValueError: Need ≥ 3 samples per class` | Too few samples | Use a cancer type with ≥ 3 cases |
| `FileNotFoundError: deepcatch_data.xlsx` | Wrong path | Use full path: `results/prof_jiang_4mer_analysis/deepcatch_data.xlsx` |
| Nested CV prints warning | `cet_cross_validator` not installed | Falls back to sklearn GridSearchCV — still valid |
| Plots look empty or missing | Matplotlib backend issue | Set `export MPLBACKEND=Agg` before running |
| `pysam` import error | Optional BAM dependency | Only needed for raw FASTQ/BAM input; ignore for frequency data |

### Still stuck?

- Check the [Pipeline Architecture](Pipeline-Architecture) page for CLI flag details
- Open an issue on [GitHub](https://github.com/rollroyces/deepcatch/issues)
- Run with `--verbose` / `-v` flag for detailed debug output:
  ```bash
  python run_jiang_analysis.py -i data.xlsx -v
  ```

---

## 🎯 Next Steps

- 📖 Understand the pipeline: **[Pipeline Architecture →](Pipeline-Architecture)**
- 🧪 Dive into the validation results: **[Jiang 4‑mer Validation →](Jiang-4mer-Validation)**
- 🏠 Back to **[Home →](Home)**
