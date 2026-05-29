# 🧬 DeepCatch — Home

**DeepCatch is an open‑source computational framework for ultra‑early cancer detection from cfDNA (cell‑free DNA), combining performance‑weighted multi‑modal fusion, cumulative evidence tracking (CET), and fragmentomics analysis into a single reproducible pipeline.** Originally designed as a simulation framework, DeepCatch v2.1 now includes preliminary real‑world validation on 129 human plasma samples from Professor Jiang's lab at CUHK, demonstrating **AUC 0.986 (nested CV)** for hepatocellular carcinoma (HCC) detection using 4‑mer end‑motif analysis.

---

## ⚡ Quick Start

```bash
git clone https://github.com/rollroyces/deepcatch.git
cd deepcatch
pip install -r requirements.txt
python run_jiang_analysis.py -i deepcatch_data.xlsx --cancer-type HCC --control-label Control --top-k 50 --nested-cv --optimal-k --plot --report
```

📖 **[Getting Started →](Getting-Started)** for detailed setup instructions.

---

## 🛡️ Badges

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)
[![Version: 2.1](https://img.shields.io/badge/Version-2.1-blue.svg)]()

To add these badges to your own fork or mirror, copy the Markdown below into your README:

```markdown
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)
[![Version: 2.1](https://img.shields.io/badge/Version-2.1-blue.svg)()
```

---

## 📊 Key Results — Jiang Lab 4‑mer Validation

| Metric | Value |
|--------|-------|
| **Samples** | 72 (34 HCC, 38 Control) |
| **Nested CV AUC** | **0.986** |
| **CV AUC (optimal k=5)** | **0.996** |
| **Bonferroni‑significant motifs** | 108 / 256 |
| **FDR‑significant motifs** | 164 / 256 |
| **Biological pattern** | CG‑rich depletion + AT‑rich enrichment |
| **Top motif (enriched)** | AAAA |
| **Top motif (depleted)** | CCCG |

> ⚠️ **Important:** These results are from processed frequency data, not raw sequencing reads. Multi‑centre replication is pending. See the **[Jiang 4‑mer Validation →](Jiang-4mer-Validation)** page for full details and caveats.

---

## 🏗️ System Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           DeepCatch System Architecture                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐     │
│  │  Variant      │  │  Methylation │  │  Fragment-   │  │  Copy Number │     │
│  │  Calling      │  │  (Entropy)   │  │  omics        │  │  Alterations │     │
│  └──────┬────────┘  └──────┬───────┘  │  (GMM+MDS+    │  │  (CNA)       │     │
│         │                  │           │   LOESS)      │  │              │     │
│         │                  │           └───────┬───────┘  └──────┬───────┘     │
│         └──────────────────┼───────────────────┼──────────────────┘            │
│                            ▼                   ▼                               │
│               ┌────────────────────────────────────┐                           │
│               │  Multi‑Modal Fusion Layer           │                           │
│               │  (Performance‑Weighted)              │                           │
│               └──────────────┬─────────────────────┘                           │
│                              ▼                                                 │
│               ┌────────────────────────────────────┐                           │
│               │  Two‑Stage CET Screening            │                           │
│               │  Stage 1: Permissive SPRT           │                           │
│               │  Stage 2: Strict Confirmation       │                           │
│               └──────────────┬─────────────────────┘                           │
│                              ▼                                                 │
│               ┌────────────────────────────────────┐                           │
│               │  Risk Score + TOO Prediction        │                           │
│               │  → 4‑Tier Risk Stratification       │                           │
│               └────────────────────────────────────┘                           │
└─────────────────────────────────────────────────────────────────────────────┘
```

📖 **[Pipeline Architecture →](Pipeline-Architecture)** for an in‑depth walkthrough of each component.

---

## 📚 Wiki Pages

| Page | Description |
|------|-------------|
| **[Getting Started](Getting-Started)** | Installation, quick demo, expected outputs & troubleshooting |
| **[Pipeline Architecture](Pipeline-Architecture)** | Two‑stage CET, FragmentoSign, feature engineering & CLI reference |
| **[Jiang 4‑mer Validation](Jiang-4mer-Validation)** | Real‑world plasma validation results, caveats, and biological interpretation |

---

## 🔬 Research‑Only Disclaimer

**DeepCatch is research‑stage software. Do not use for medical diagnosis, treatment decisions, or any clinical purpose.** The Jiang lab validation uses processed frequency data (not raw FASTQ/BAM) and is a computational feasibility study, not a clinical assay.
