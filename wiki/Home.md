# 🧬 DeepCatch — Home

**DeepCatch is an open‑source computational framework for ultra‑early cancer detection from cfDNA (cell‑free DNA), combining performance‑weighted multi‑modal fusion, cumulative evidence tracking (CET), and fragmentomics analysis into a single reproducible pipeline.** Originally designed as a simulation framework, DeepCatch v2.1 also carries preliminary real‑world work on 129 plasma samples from Professor Jiang's lab at CUHK. The hepatocellular carcinoma (HCC) vs control comparison inside that cohort is 34 patients against 38 controls: the committed 5‑fold stratified run reports **AUC 0.985** at 94.1% sensitivity (95% specificity), and a nested cross‑validation of the same comparison reports **0.986** — that run's input and output are not committed here, so the tree supports the first figure and not the second. See the [Jiang 4‑mer Validation →](Jiang-4mer-Validation) page.

---

## ⚡ Quick Start

```bash
git clone https://github.com/rollroyces/deepcatch.git
cd deepcatch
pip install -e .
python run_jiang_analysis.py -i deepcatch_data.xlsx --cancer-type HCC --control-label Control --top-k 50 --nested-cv --optimal-k --plot --report
```

📖 **[Getting Started →](Getting-Started)** for detailed setup instructions.

---

## 🛡️ Badges

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/rollroyces/deepcatch/blob/main/LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-green.svg)](https://www.python.org/)
[![Version: 2.3.0](https://img.shields.io/badge/Version-2.3.0-blue.svg)]()

To add these badges to your own fork or mirror, copy the Markdown below into your README:

```markdown
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/rollroyces/deepcatch/blob/main/LICENSE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-green.svg)](https://www.python.org/)
[![Version: 2.3.0](https://img.shields.io/badge/Version-2.3.0-blue.svg)()
```

---

## 📊 Key Results — Jiang Lab 4‑mer Validation

| Metric | Value |
|--------|-------|
| **Samples** | 72 (34 HCC, 38 Control), out of 129 plasma samples in the pan‑cancer cohort |
| **CV AUC (k=50)** | **0.985** — the committed run (`results/prof_jiang_4mer_analysis/part2_cet_results.txt`) |
| **Nested CV AUC** | 0.986 — from a run whose input and output are not committed here |
| **CV AUC (optimal k=5)** | 0.996 — same run, likewise not committed |
| **Bonferroni‑significant motifs** | 108 / 256 |
| **FDR‑significant motifs** | 164 / 256 |
| **Biological pattern** | CG‑rich depletion + AT‑rich enrichment |
| **Top motif (enriched)** | AAAA |
| **Top motif (depleted)** | CCCG |

> ⚠️ **Important:** These results are from processed frequency data, not raw sequencing reads. Multi‑centre replication is pending. The committed artifact supports **0.985 / 94.1%** (n=72); the nested‑CV and optimal‑k figures above come from a run whose input and output are not in this repository. This page describes the v2.1 Jiang‑lab work: the current manuscript ([`paper/paper.tex`](../paper/paper.tex), v2.3.0) makes **no claim** about these 129 samples, and whether the collaborating laboratory is credited in it is the author's decision. See the **[Jiang 4‑mer Validation →](Jiang-4mer-Validation)** page for full details and caveats.

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
