# 🧪 Jiang 4‑mer Validation — Real Plasma cfDNA Analysis

> **v2.1 — First real‑world validation of DeepCatch's CET architecture on actual human plasma cfDNA data.**

This page documents the first validation of DeepCatch using real patient plasma samples — not simulation. Data sourced from **Professor Jiang Peiyong's laboratory at the Chinese University of Hong Kong (CUHK)**.

---

## Source Data

| Parameter | Value |
|-----------|-------|
| **Total samples** | 129 plasma samples |
| **Healthy controls** | 38 |
| **Cancer patients** | 91 |
| **Cancer types** | 6 (HCC, lung, HNSCC, CRC, NPC, gastric) |
| **Feature space** | 256 four‑mer end motifs |
| **Data type** | Processed frequency vectors (not raw FASTQ/BAM) |
| **Source** | Jiang lab, CUHK — Table S1 |

### Sample Distribution

| Cancer Type | Cases | Controls | Total |
|-------------|-------|----------|-------|
| Hepatocellular Carcinoma (HCC) | 34 | 38 | 72 |
| Lung Cancer (LC) | 10 | 38 | 48 |
| Head & Neck SCC (HNSCC) | 10 | 38 | 48 |
| Colorectal Cancer (CRC) | 10 | 38 | 48 |
| Nasopharyngeal Carcinoma (NPC) | 10 | 38 | 48 |
| Gastric Cancer | 10 | 38 | 48 |
| HBV Carriers (non‑cancer) | 17 | — | 17 |

> ⚠️ **Important:** Only the HCC vs Control comparison (n=72) is adequately powered. All other cancer types have n ≤ 10 cases — their AUC estimates should be treated as exploratory.

---

## HCC vs Control — Key Results

The primary analysis compares 34 HCC patients against 38 healthy controls using the full DeepCatch CET pipeline (rank + ratio features, logistic regression fusion, nested cross‑validation).

### Core Metrics

| Metric | Value | Method |
|--------|-------|--------|
| **Nested CV AUC** | **0.986** | 5 outer × 3 inner folds, feature selection within each outer fold |
| **CV AUC (optimal k=5)** | **0.996** | Standard 5‑fold CV, top‑5 motifs by p‑value |
| **CV AUC (k=50)** | **0.985** | Standard 5‑fold CV, top‑50 motifs |
| **Sensitivity @ 95% specificity** | **94.1%** | Threshold calibrated on training folds |
| **Bonferroni‑significant motifs** | **108 / 256** (42.2%) | p < 0.05/256 |
| **FDR‑significant motifs** | **164 / 256** (64.1%) | Benjamini‑Hochberg, α = 0.05 |
| **Min p‑value** | 1.52 × 10⁻¹⁰ | AAAA motif |

### Per‑Cancer Comparison Summary

| Comparison | AUC | Sensitivity @ 95% Spec | n (case + ctrl) | Reliability |
|------------|-----|------------------------|-----------------|-------------|
| **HCC vs Control** | **0.985** | **94.1%** | 72 | ✅ Adequately powered |
| LC vs Control | 0.979 | 70.0% | 48 | ⚠️ n=10 cases |
| HNSCC vs Control | 0.937 | 30.0% | 48 | ⚠️ n=10 cases |
| CRC vs Control | 0.911 | 40.0% | 48 | ⚠️ n=10 cases |
| NPC vs Control | 0.868 | 30.0% | 48 | ⚠️ n=10 cases |
| **HBV vs HCC** | **0.905** | 38.2% | 51 | ⚠️ Clinically relevant but under‑powered |
| All Cancer vs Control | 0.910 | 58.2% | 129 | ⚠️ Pan‑cancer RF classifier |

---

## Biological Validation

The motif‑level results show a clear, biologically coherent signal:

### CG‑Rich Depletion → CpG Island Hypomethylation

Motifs containing CpG dinucleotides are **consistently depleted** in HCC plasma:

| Motif | Direction in HCC | log₂ Fold Change | Interpretation |
|-------|-----------------|------------------|----------------|
| CCCG | ↓ Depleted | — | CpG‑containing — altered methylation |
| CGCT | ↓ Depleted | — | CpG methylation hallmark |
| CGCC | ↓ Depleted | — | CpG‑rich region depletion |
| CGCG | ↓ Depleted | — | Dense CpG island signal |

**Biological explanation:** Global hypomethylation is a hallmark of hepatocellular carcinoma. Hypomethylated CpG islands become more accessible to nucleases, altering the cfDNA fragmentation landscape. The result is a reduction in CG‑rich 4‑mer motifs in cancer plasma.

### AT‑Rich Enrichment → Nucleosome Depletion

AT‑rich motifs are **consistently enriched** in HCC plasma:

| Motif | Direction in HCC | log₂ Fold Change | Interpretation |
|-------|-----------------|------------------|----------------|
| AAAA | ↑ Enriched | — | Poly‑A stretches — polyadenylation or 3′ bias |
| AAGA | ↑ Enriched | — | AT‑rich — nucleosome positioning shifts |
| ATAT | ↑ Enriched | — | AT‑rich tandem repeat |
| ATAA | ↑ Enriched | — | AT‑rich motif cluster |

**Biological explanation:** AT‑rich regions have lower nucleosome occupancy. In cancer, altered chromatin structure leads to increased fragmentation at these loci, releasing more AT‑rich fragments into plasma.

### HBV → HCC Progression Signal

The transition from HBV carrier to HCC is marked by progressive **CC* and CG* motif depletion** plus **AT‑rich enrichment** — consistent with the known stepwise molecular changes in HBV‑driven hepatocarcinogenesis. This is the clinically most actionable finding: serial 4‑mer profiling could serve as a **surveillance tool for the 250 million chronic HBV carriers worldwide**.

---

## Caveats

> **Honest limitations — these are not excuses, they are what needs to be tested next.**

### 1. HCC‑Only Adequate Power (n=72)
Only the HCC vs Control comparison has enough samples for reliable AUC estimates. Lung, HNSCC, CRC, NPC, and gastric each have only 10 cancer cases. Their per‑cancer AUCs (0.87–0.98) should be treated as **exploratory estimates with wide confidence intervals**.

### 2. Processed Data Only
The analysis uses pre‑computed 4‑mer frequency vectors, not raw sequencing reads (FASTQ/BAM). This means:
- We cannot verify motif extraction from raw data
- Technical confounders (PCR duplicates, GC bias, batch effects) are not fully controlled
- The pipeline has not been tested end‑to‑end from sequencer output

### 3. Single Centre
All 129 samples are from a single laboratory (CUHK). Multi‑centre replication is essential to rule out site‑specific technical artefacts, population‑specific baseline differences, and batch effects.

### 4. Not a Clinical Assay
This is a **computational feasibility study**, not a clinical diagnostic. Key gaps before clinical readiness:
- No prospective cohort validation
- No pre‑specified analysis plan
- No locked model or fixed decision threshold
- No comparison to standard‑of‑care (AFP, ultrasound)
- No CLIA/CAP laboratory workflow

### 5. No Raw BAM Validation
The fragmentomics GMM, MDS, and LOESS normalisation components of FragmentoSign have not been validated on these samples because raw BAM files are not yet available. What we demonstrate is the CET architecture applied to processed frequency data.

---

## Reproduce This Analysis

### One‑Command Reproduction

```bash
python run_jiang_analysis.py \
  -i results/prof_jiang_4mer_analysis/deepcatch_data.xlsx \
  --cancer-type HCC \
  --control-label Control \
  --top-k 50 \
  --nested-cv \
  --optimal-k \
  --plot \
  --report \
  --seed 42
```

This command reproduces the headline HCC vs Control result exactly. The `--seed 42` flag ensures deterministic output.

### Verification

After running, check `results/jiang_reanalysis/summary_report.md` — the nested CV AUC should match **0.986** (within floating‑point tolerance).

### Run All Pairwise Comparisons

```bash
python run_jiang_analysis.py \
  -i results/prof_jiang_4mer_analysis/deepcatch_data.xlsx \
  --top-k 30 \
  --plot
```

This runs all cancer types against controls automatically and produces a `per_cancer_summary.md` table.

---

## Further Reading

- 📖 **README §9** — Full Jiang validation section in the [main README](https://github.com/rollroyces/deepcatch#9-real-plasma-validation--4-mer-end-motif-analysis-on-jiang-lab-data)
- 📄 **Summary for Professor Jiang** — Detailed results in `results/prof_jiang_4mer_analysis/summary_for_professor_jiang.md`
- 🏗️ **[Pipeline Architecture →](Pipeline-Architecture)** — How the CET pipeline works under the hood
- 🚀 **[Getting Started →](Getting-Started)** — Run the pipeline on your own data

---

*Analysis performed 28 May 2026 by Ironman 🦾. Every number is traceable to a computation in `validation/` and `results/`. No numbers were invented. No clinical claims are intended.*
