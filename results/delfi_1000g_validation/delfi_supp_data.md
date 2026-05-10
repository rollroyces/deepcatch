# DELFI Published Statistics — Compiled from Literature

## Source: Cristiano et al. 2019, Nature 570:385-389
### "Genome-wide cell-free DNA fragmentation in patients with cancer"

---

## Study Overview
- **Sample size**: 236 cancer patients + 245 healthy controls
- **Cancer types**: 7 (breast, colorectal, lung, ovarian, pancreatic, gastric, bile duct)
- **Sequencing**: Low-coverage (~1-2×) whole-genome sequencing of cfDNA
- **Key innovation**: Genome-wide fragmentation profiles (DELFI) for cancer detection
- **ML method**: Gradient tree boosting
- **Performance**: AUC 0.94; 57–99% sensitivity across 7 cancers at 98% specificity

---

## Fragment-Level Statistics (from paper main text, figures, and supplementary)

### Fragment Size Distribution
- **Healthy cfDNA**: Peak at ~167 bp (nucleosome protection), with 10-bp periodicity
- **Cancer cfDNA**: 
  - More variable fragment sizes
  - Increased proportion of short fragments (<150 bp)
  - Decreased proportion at nucleosome peak (~167 bp)
  - Altered 10-bp periodicity patterns (oscillations)
- **Fragment size ratio**: Cancer patients showed significantly different short-to-long fragment ratios compared to healthy controls across multiple genomic windows

### Genome-Wide Fragmentation Profiles
- DELFI analyzed fragmentation patterns in **5-Mb windows** across the genome
- Each window had ~50,000 fragments
- **504 non-overlapping windows** across autosomes
- Features extracted per window:
  - Coverage (normalized fragment counts)
  - Fragment size distribution (short/long ratios)
  - Fragment end motifs
  - Copy number alterations

### MDS (Motif Diversity Score) Reference Values
*From Jiang et al. 2020, Cancer Discovery 10:664-673 (end-motif profiling paper)*:

| Category | MDS Range | Description |
|---|---|---|
| **Healthy cfDNA** | 0.92–0.96 | Normal nucleosome-mediated fragmentation |
| **Cancer patient cfDNA** | 0.88–0.92 | Altered nuclease activity, more selective fragmentation |
| **Advanced cancer** | 0.85–0.88 | Strongly aberrant fragmentation patterns |
| **Synthetic/random** | >0.97 | No biological fragmentation signal |

### Tissue-Specific Fragmentation
- Different cancer types exhibit distinct fragmentation patterns
- Lung cancer: most divergent fragmentation profiles
- Colorectal cancer: strong fragmentation signal in specific genomic regions
- Breast cancer: moderate fragmentation alterations

### Fragment End Motif Patterns
- **CCCA** motif: enriched in cfDNA (associated with DNASE1L3 cleavage preference)
- **CCAG/CCTG/CTGG**: DNASE1 family cleavage signatures
- Cancer-associated shifts in preferred cleavage motifs
- 4-mer end motifs provide information complementary to fragment size

### Performance by Cancer Stage
| Stage | Sensitivity (at 98% specificity) |
|---|---|
| Stage I | 57% (average across cancers) |
| Stage II | 67% |
| Stage III | 78% |
| Stage IV | 90% |
| All stages | 73% |

### Combination with Mutation Analysis
- DELFI alone: 73% sensitivity at 98% specificity
- DELFI + mutations: 91% sensitivity at 98% specificity
- Demonstrates value of multi-modal analysis (fragmentomics + genomics)

### Technical Validation
- **Reproducibility**: High correlation between technical replicates (R² > 0.95)
- **Library preparation**: Compatible with standard WGS library prep
- **Input requirement**: As low as 1 ng cfDNA
- **Sequencing depth**: ~1-2× coverage sufficient

---

## Comparison Studies

### Mathios et al. 2021 (Nat Commun) — DELFI for Lung Cancer
- **Sample**: 365 screening + 431 validation
- **Performance**: AUC 0.98 (SCLC vs NSCLC classification); 94% detection at 80% spec
- **Key finding**: DELFI distinguish SCLC from NSCLC with high accuracy
- **Clinical integration**: Combined with clinical risk factors improved performance

### Mazzone et al. 2024 (Cancer Discov) — DELFI Clinical Validation
- **Sample**: 958 for clinical validation (lung cancer)
- **Performance**: 51% sens @ 80% spec (all stages); 40% stage I
- **Real-world implementation**: Demonstrated clinical feasibility

---

## Key Biological Insights for DeepCatch Validation

1. **cfDNA fragmentation is non-random**: Nucleosome positioning, nuclease preferences (DNASE1, DNASE1L3, DFFB), and chromatin structure all influence where DNA breaks
2. **Cancer alters fragmentation**: Both the nucleosome landscape and nuclease activity change in cancer, producing measurable shifts in fragment patterns
3. **MDS is a robust single-number summary**: Motif Diversity Score captures the selectivity of fragmentation — lower MDS indicates more selective (typically cancer-associated) cleavage
4. **Whole-blood WGS vs cfDNA**: Whole-blood WGS (like 1000 Genomes) fragments differently than cfDNA due to different nucleosome structures and different DNase environments
5. **Synthetic data**: Random fragmentation produces MDS close to 1.0 — clearly distinguishable from biological samples

---

## References
1. Cristiano, S. et al. (2019). Nature 570:385-389. PMID: 31142840
2. Jiang, P. et al. (2020). Cancer Discovery 10(5):664-673. PMID: 32111602
3. Mathios, D. et al. (2021). Nat Commun 12:5060. PMID: 34417454
4. Mazzone, P.J. et al. (2024). Cancer Discov. PMID: 38270515
