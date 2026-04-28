# DeepCatch Submission Checklist
## Target: Bioinformatics (Oxford Academic) — Computational Methods Track

---

## Manuscript Files

- [x] **main_final.tex** — Complete LaTeX manuscript (~45 KB)
- [x] **references_final.bib** — Bibliography with 40+ citations (~20 KB)
- [x] **abstract_final.txt** — Plain-text abstract (240 words, within 250 limit)
- [x] **cover_letter_final.txt** — Cover letter for editors
- [ ] **supplementary.tex** — Supplementary information (use existing `supplementary.tex` from `paper/`)
- [ ] **figures/** — Generate from `results/node/*.json` using matplotlib/seaborn

---

## Honesty Rules Verification

| # | Rule | Status | Evidence |
|---|------|--------|----------|
| 1 | Every number matches FINAL_REAL_DATA_REPORT | ✅ VERIFIED | AUC 0.9610, ΔAUC +0.1434, p<0.0001 all from Section 3 |
| 2 | Every limitation acknowledged | ✅ VERIFIED | Discussion Section 4.1 lists 8 limitations |
| 3 | "Simulation" appears prominently | ✅ VERIFIED | Abstract, Introduction, Results, Discussion all flag simulation |
| 4 | No clinical claims without clinical data | ✅ VERIFIED | Comparison table has "NOT Directly Comparable" header |
| 5 | Depth/cost caveats included | ✅ VERIFIED | 50,000× vs 5,000× clinical standard noted |
| 6 | CET failure discussed openly | ✅ VERIFIED | 2.5% sensitivity reported; 3 reasons for failure analyzed |

---

## Key Numbers Audit

### From FINAL_REAL_DATA_REPORT → Paper

| Data Point | Report Value | Paper Location | Match? |
|-----------|-------------|----------------|--------|
| Multi-Modal AUC @ 1% ctDNA | 0.9610 | Table 2 | ✅ |
| ΔAUC vs Bie @ 1% ctDNA | +0.1434 | Table 5 | ✅ |
| p-value @ 1% ctDNA | <0.0001 | Table 5 | ✅ |
| CET Sensitivity | 2.5% | Table 6 | ✅ |
| CET Specificity | 97.0% | Table 6 | ✅ |
| CET AUC | 0.4926 | Not in main paper (in supplement) | ✅ |
| Variant sens @ 99% spec (1% ctDNA) | 72.8% | Table 4 | ✅ |
| Cancer types | 8 | Table 1 | ✅ |
| Confounders | 6 | Table 2 | ✅ |

### From headToHead_results.json → Paper

| Data Point | JSON Value | Paper Location | Match? |
|-----------|-----------|----------------|--------|
| DC 5-mod AUC | 0.9784 | Table 3 | ✅ |
| DC 5-mod CI | [0.9721, 0.9837] | Table 3 | ✅ |
| Bie 4-mod AUC | 0.9668 | Table 3 | ✅ |
| Bie 4-mod CI | [0.9592, 0.9740] | Table 3 | ✅ |
| DC5 vs Bie4 ΔAUC | +0.0116 | Table 3 | ✅ |
| DC5 vs Bie4 p-value | <0.0001 | Table 4 | ✅ |

---

## Before Submission

### Critical
- [ ] **Generate figures** from `results/node/*.json` data files
- [ ] **Multi-seed validation** — Run with seeds 43, 44, 45, 46, 47; report mean ± SD
- [ ] **Supplementary information** — Update `supplementary.tex` with real-data results
- [ ] **Author list** — Finalize author names and affiliations
- [ ] **Corresponding author email** — Verify `royce@deepcatch.org` is active

### Recommended
- [ ] **Language editing** — Professional scientific editing for non-native English
- [ ] **Reference formatting** — Verify all DOIs and journal styles
- [ ] **Figure quality** — Ensure 300+ DPI for all figures
- [ ] **Data deposition** — Upload simulation data to Zenodo/Figshare with DOI
- [ ] **Code release** — Create GitHub repository with README and RUN_ALL.sh
- [ ] **ORCID** — Add ORCID IDs for all authors

### Pre-submission Checklist
- [ ] All author approvals obtained
- [ ] Conflict of interest statement finalized
- [ ] Funding acknowledgments added (if applicable)
- [ ] Data availability statement verified
- [ ] Code availability statement verified
- [ ] Manuscript formatted per Bioinformatics author guidelines
- [ ] Page limits checked (Bioinformatics: ~8-10 pages for Original Article)
- [ ] Abstract ≤ 250 words (current: ~240)
- [ ] Keywords added: liquid biopsy, circulating tumor DNA, early cancer detection, multi-modal fusion, longitudinal screening, variant calling, simulation study

---

## Journal-Specific Requirements

### Bioinformatics (Oxford Academic)
- [ ] Use `\documentclass{bioinfo}` (currently `article` — change before submission)
- [ ] Author names: Firstname Lastname^affiliation^ format
- [ ] Graphical abstract (optional but recommended)
- [ ] Supplementary data as separate PDF
- [ ] Data availability section required

### OR: PLOS Computational Biology
- [ ] Use PLOS LaTeX template
- [ ] Author summary (non-technical, 150-200 words)
- [ ] Supporting information captions in main text
- [ ] Data availability statement (PLOS requires deposition)

---

## Figures to Generate

| Figure | Data Source | Description |
|--------|-----------|-------------|
| Fig 1 | Architecture schematic | DeepCatch pipeline overview (4 panels) |
| Fig 2 | `headToHead_results.json` | AUC vs ctDNA fraction (all methods) |
| Fig 3 | `headToHead_results.json` + `fusion_results.json` | ΔAUC bar chart with significance stars |
| Fig 4 | `real_cet_results.json` | CET score trajectories + ROC |
| Fig 5 | `published_comparison.json` | Clinical comparison forest plot |
| Fig S1 | `real_tcga_data.json` | Per-cancer-type shedding rate distributions |
| Fig S2 | `variant_calling_results.json` | VAF sensitivity curves |

---

## Submission Checklist Summary

```
manuscript:   ✅ main_final.tex
references:   ✅ references_final.bib
abstract:     ✅ abstract_final.txt
cover letter: ✅ cover_letter_final.txt
checklist:    ✅ SUBMISSION_CHECKLIST.md (this file)
supplement:   ⚠️ needs update with real-data results
figures:      ⚠️ needs generation from results JSON
multi-seed:   ⚠️ not yet performed
```

---

*Last updated: 2026-04-28 — Ironman 🦾*
