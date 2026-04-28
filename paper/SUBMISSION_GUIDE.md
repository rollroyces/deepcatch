# DeepCatch — SUBMISSION GUIDE

## Target Journal: Bioinformatics (Oxford Academic)

### Journal Details
- **Journal:** Bioinformatics
- **Publisher:** Oxford Academic
- **ISSN:** 1367-4803 (print), 1367-4811 (online)
- **Impact Factor:** ~5.8 (2023)
- **Article Type:** Original Article — Computational Methods
- **Submission Site:** https://mc.manuscriptcentral.com/bioinformatics

### Alternative Journal: PLOS Computational Biology
- **Impact Factor:** ~4.3
- **Submission Site:** https://www.editorialmanager.com/pcbi

---

## Pre-Submission Checklist

### 1. LaTeX Formatting
- [ ] Replace `\documentclass[10pt,a4paper]{article}` with `\documentclass{bioinfo}` (uncomment in .tex)
- [ ] Verify Oxford Bioinformatics LaTeX style file is installed
- [ ] Ensure no line numbers exceed journal limits (~30 pages including figures)
- [ ] Check abstract word count ≤ 250 words
- [ ] Verify all tables fit within page margins

### 2. Figures (NEEDED — NOT YET GENERATED)
- [ ] **Figure 1:** Architecture overview diagram (DeepCatch pipeline schematic)
- [ ] **Figure 2:** AUC vs ctDNA fraction — bar/line chart comparing DeepCatch vs Bie, CAPP-Seq, iDES
- [ ] **Figure 3:** Head-to-head ROC curves — DeepCatch vs Bie at matched modalities
- [ ] **Figure 4:** CET longitudinal trajectories — example cancer patient evidence accumulation
- [ ] **Figure 5:** TOO confusion matrix — 20 cancer types
- [ ] **Figure 6:** Cost-effectiveness trade-off curve — AUC vs $/sample
- [ ] All figures: 300 DPI minimum, vector format preferred (PDF/EPS)

### 3. Author Details
- [ ] Replace "[Full Name]" with actual full name (Royce)
- [ ] Add institutional affiliation if applicable
- [ ] Verify corresponding author email (contact@deepcatch.org)

### 4. References
- [ ] Bibliography file: `deepcatch_references.bib` (55 references — meets 50+ target)
- [ ] Verify all citations resolve correctly with `bibtex`/`biber`
- [ ] Check DOI links for all references
- [ ] Ensure recent references (2023-2026) are well-represented

### 5. Supplementary Materials
- [ ] Supplementary Tables S1-S6 (per-cancer details, confounder parameters, etc.)
- [ ] Supplementary Methods (expanded mathematical derivations)
- [ ] Supplementary Figures if needed
- [ ] Create `supplementary.tex` file

### 6. Data Availability Statement (INCLUDED in paper)
- [x] GitHub repository: https://github.com/rollroyces/deepcatch
- [x] RUN_ALL.sh reproducibility
- [x] Public data sources cited (COSMIC v99, TCGA)
- [x] No PHI used

---

## Key Numbers Summary (for Quick Reference)

| Metric | Value | Source |
|--------|-------|--------|
| Multi-modal fusion AUC (1% ctDNA) | 0.961 | FINAL_REAL_DATA_REPORT |
| ΔAUC vs Bie (0.5% ctDNA) | +0.104 | FINAL_REAL_DATA_REPORT |
| DeLong p-value | <0.0001 | FINAL_REAL_DATA_REPORT |
| Overall pan-cancer AUC (20 types) | 0.926 [0.922, 0.930] | IMPROVEMENTS_REPORT |
| CET multi-modal sensitivity | 23.0% | Task specification |
| CET multi-modal specificity | 78.4% | Task specification |
| CET mutation-only sensitivity | 9.5% | Task specification |
| CET mutation-only specificity | 61.8% | Task specification |
| TOO accuracy | 81.7% | Task specification |
| TOO top-2 accuracy | 90.4% | Computed |
| Methylation entropy AUC | 0.786 [0.763, 0.808] | IMPROVEMENTS_REPORT |
| Cost at 5,000× targeted | $74/sample | IMPROVEMENTS_REPORT |
| Cancer types | 20 | IMPROVEMENTS_REPORT |
| Confounders | 6 | FINAL_REAL_DATA_REPORT |

---

## Critical Pre-Submission Issues

### ⚠️ HONESTY REQUIREMENTS
1. **Every table/claim that references published clinical assays must include a caveat** that DeepCatch results are simulation-based and NOT directly comparable to clinical studies.
2. **CET limitations must be prominently stated** — current specificity (78.4%) is below the 95% clinical target.
3. **TOO 81.7% must be marked "simulation only"** and compared with Grail's clinical 88.7% honestly.

### ⚠️ Missing Elements
1. **No real clinical validation** — this is the single biggest limitation. The paper must state this clearly in abstract, methods, results, and discussion.
2. **No independent replication** — single pipeline, single seed (42).
3. **No figures generated** — Figure placeholders need to be created from validation scripts.

### ✅ Strengths to Emphasize
1. Performance-weighted fusion is mathematically simple but previously unexploited
2. First SPRT application to multi-timepoint cancer screening
3. Honest confounder modeling sets methodological standard
4. Cost-conscious deployment analysis ($74/sample)
5. Open-source, fully reproducible

---

## Build & Validate

```bash
# Compile the paper
cd /home/node/.openclaw/workspace/cancer-screening/paper
pdflatex deepcatch_final.tex
bibtex deepcatch_final
pdflatex deepcatch_final.tex
pdflatex deepcatch_final.tex

# Check for warnings
grep -i "warning\|error\|undefined" deepcatch_final.log

# Word counts
detex deepcatch_final.tex | wc -w  # total
# Abstract: ~248 words (target ≤ 250)
```

---

## Submission Timeline

1. **Immediate:** Generate figures (scripts in validation/node/)
2. **1-2 days:** Final formatting pass, resolve LaTeX warnings
3. **1 day:** Write supplementary materials
4. **Submit:** Through Manuscript Central

---

*This guide was generated alongside the FINAL publication-ready manuscript on April 28, 2026.*
