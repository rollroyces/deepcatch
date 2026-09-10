# DeepCatch Methylation Project — Scope & Plan

**Status:** Active planning (started 2026-09-10)
**Author:** Yu Ching Lam (Independent Researcher)
**ORCID:** [0009-0008-9113-769X](https://orcid.org/0009-0008-9113-769X)
**Scope:** CpG methylation-based cfDNA cancer detection as a complementary channel to the existing fragmentomics pipeline.

---

## 1. Why this project exists

The fragmentomics pipeline (`rollroyces/cfdna-fragmentomics-pipeline`) achieves **AUC 0.974-0.978** on a 627-sample cross-study cohort using only **fragment-length features** (FSD + DELFI 5Mb/100kb ratios). The DeepCatch framework (`rollroyces/deepcatch`) adds a **mutation-channel fusion** that lifts headline to **AUC 0.989** on a synthetic cohort.

**What fragmentomics cannot do:**
- Identify **tissue of origin** (FOI) beyond what the cancer-vs-healthy signal already encodes.
- Achieve the >50% sensitivity at 99.5% specificity GRAIL's methylation-based Galleri reports across 50+ cancer types (Liu 2020, CCGA-1).
- Detect **pre-cancer epigenetic field defects** — methylation changes that precede any fragmentomic shift.

**What methylation adds (per round-4 journal reviewer audit):**
> "Methylation is the largest single signal in any MCED context."

The framework already has **2,200 lines of methylation scaffolding** (`src/methylation_gnn/`, GATv2-based pre-cancer field-defect detector) and **46 unit tests**, but **zero production wiring**: the headline fusion pipeline (`fusion_ablation.py`, `multimodal_fusion/advanced_fusion.py`) never imports it. There is no real-data methylation feature in either repo.

This project closes that gap.

---

## 2. Headline research question

> **Does CpG methylation add signal beyond fragmentomics on a 627-sample cross-study cohort?**

The honest answer is currently **unknown** because no head-to-head has ever been run. Three outcomes are scientifically interesting:

| Outcome | Interpretation |
|---|---|
| Methylation AUC ≪ Fragmentomics AUC | Fragmentomics dominates; methylation provides little marginal value at this sample size |
| Methylation AUC ≈ Fragmentomics AUC | Both signals extract overlapping biology; combination gives modest lift (ΔAUC ~0.01-0.02) |
| Methylation AUC ≫ Fragmentomics AUC | Methylation carries tissue-specific signal fragmentomics cannot access; combination gives large lift (ΔAUC ≥0.05) |

GRAIL's CCGA-1 numbers (sensitivity 0.50, specificity 0.99 across 50+ cancer types; Liu 2020) suggest methylation wins on broad pan-cancer, but their cohort and assay are materially different.

---

## 3. Strategy: WGS → methylation imputation (no new assay needed)

The critical practical decision: **do not require new methylation data.** FinaleMe (Liu et al., *Nature Communications* 15:2790, 2024; MIT-licensed; https://github.com/epifluidlab/FinaleMe) **predicts single-CpG methylation status from plain cfDNA WGS** without bisulfite conversion.

| Property | Value |
|---|---|
| Input | FinaleDB-style WGS fragments (the data we already ingest) |
| Output | Per-CpG methylation probability + tissue-of-origin score |
| Validated auROC | 0.91 on fragments with ≥5 CpGs in CpG-rich regions |
| Training data | Sun 2015 / Sun 2018 WGBS (dbGaP phs000846 / phs001417) |
| Cost | Zero new raw bytes; the WGS we already have is sufficient |

This is the **only realistic path to methylation features on a 627-sample cohort within 6 months**. All alternative paths (targeted bisulfite sequencing, EM-seq, WGBS) require either (a) DAC approval (months), (b) purchasing arrays ($100K+), or (c) waiting for a public WGBS cfDNA cohort at the scale of CCGA — which does not exist.

**Caveat:** FinaleMe auROC drops in CpG-poor regions and has documented cancer-WGS false positives. It is **not** a substitute for true bisulfite sequencing, but it is the best public option for adding methylation to the existing fragmentomics cohort.

---

## 4. Phased plan

### Phase 0 — Validation: does FinaleMe work on FinaleDB fragments? (2-3 weeks)

**Goal:** Verify that FinaleMe runs cleanly on FinaleDB WGS BED/bigWig inputs and produces methylation features with a meaningful cancer-vs-healthy signal.

**Tasks:**
1. Install FinaleMe + FinaleToolkit (`pip install FinaleToolkit`).
2. Run FinaleMe on a small subset of FinaleDB samples (e.g., 50 samples: 25 HCC + 25 healthy from Jiang 2018).
3. Compute per-CpG methylation features:
   - Mean β-value per chromosome arm
   - CpG-island-level methylation
   - Tissue-of-origin score (FinaleMe's built-in TOO classifier)
4. Train an L2-LR baseline on these features (5-fold CV, single-study) and report AUC.
5. **Honest target:** AUC > 0.6 (i.e. above chance). If FinaleMe doesn't even produce above-chance signal on the existing data, stop and revisit the data strategy.

**Deliverable:** `results/finaleme_smoke_validation.json` + a short report documenting the subset, pipeline, and AUC.

**Exit criterion:** continue to Phase 1 only if AUC > 0.65 on a 50-sample subset.

### Phase 1 — Methylation feature pipeline on the full 627 cohort (4-6 weeks)

**Goal:** Produce a methylation-feature matrix that aligns with the existing fragmentomics cohort.

**Tasks:**
1. Run FinaleMe on all 627 FinaleDB samples (Cristiano 2019 + Jiang 2018 + others).
2. Aggregate to per-sample feature matrix:
   - Per-chromosome mean β-value (24 features)
   - Per-CpG-island coverage statistics (~28K CpG islands → 50K features after variance filtering)
   - TOO score vector (5 features: WBC / hepatocyte / colon / lung / breast contributions)
3. Same 5-fold pooled OOF CV protocol as fragmentomics (Cristiano 2019 train/val split, Jiang 2018 held-out).
4. **Headline metric:** cross-study AUC for methylation-only baseline (LR no-PCA, C=1).

**Deliverable:**
- `scripts/finaleme_extract.py` — WGS → methylation features
- `scripts/methylation_baseline.py` — LR baseline + cross-study AUC
- `results/methylation_baseline.json`

**Target:** Cross-study AUC > 0.75. (Lower than fragmentomics AUC 0.974 is expected because FinaleMe is an imputation, not a direct assay. The question is whether the signal is large enough to be useful.)

### Phase 2 — Head-to-head comparison + fusion (3-4 weeks)

**Goal:** Quantify the methylation contribution above fragmentomics.

**Tasks:**
1. Combine fragmentomics features (5 channels: FSD + DELFI 5Mb + DELFI 100kb) with methylation features (Phase 1 output) into a single feature matrix.
2. Run 3 model configurations:
   - **Fragmentomics-only** (re-run existing pipeline for fresh baseline)
   - **Methylation-only** (Phase 1 output)
   - **Fragmentomics + Methylation** (concatenated feature matrix, LR no-PCA, C=1)
3. 5-fold pooled OOF CV; report AUC, ΔAUC vs fragmentomics-only, 95% DeLong CI, paired t-test, multiple-testing correction.
4. Run the same experiment with `MethylationBranchAdapter` (from existing scaffolding) plugged into the multimodal fusion layer.

**Deliverable:**
- `scripts/methylation_vs_fragmentomics.py` — head-to-head comparison
- `scripts/methylation_fusion.py` — combined feature + fusion ablation
- `results/methylation_head_to_head.json` + `results/methylation_fusion.json`
- An honest **"did methylation help?"** section in the project report, regardless of direction.

**Target:** Document the ΔAUC honestly. If positive: ≥0.005 AUC lift with p<0.05. If negative or null: still publish — negative results matter.

### Phase 3 — Optional: tissue-of-origin (TOO) ablation (2-3 weeks, only if Phase 1 shows TOO AUC > 0.6)

**Goal:** Test whether the methylation TOO score alone can identify cancer type (HCC vs CRC vs LUAD vs BRCA, etc.) within the cross-study cohort.

**Tasks:**
1. Train per-cancer-type methylation baseline using only FinaleMe's TOO output.
2. Report per-cancer-type AUC + confusion matrix.
3. Compare against fragmentomics per-cancer-type AUC (which is currently undocumentable due to lack of `cancer_type` metadata in labels).

**Deliverable:** `results/methylation_too.json`

**Decision rule:** Skip Phase 3 if Phase 1 shows TOO AUC < 0.6, or if sample size per cancer type is < 30.

### Phase 4 — Validation of the existing GNN scaffolding (2-4 weeks, deferred until a methylation expert collaborator is found)

**Goal:** Train the existing `src/methylation_gnn/` GATv2 model on a real methylation cohort.

**Tasks:**
1. Source a public methylation dataset with at least 100 cancer + 100 healthy samples (TCGA-LIHC is the most likely candidate — 380 tumors + 50 matched normal).
2. Wire the dataset through `RegulatoryGraphBuilder` with real β-values and real region BED files (download UCSC CpG islands, GENCODE v44 promoters as documented in `methylation_gnn/data.py`).
3. Train `MethylationGNN` (3× GATv2Conv, dual head) for full 2-phase pretrain + finetune.
4. Compare against the Phase 1 LR baseline.

**Deliverable:**
- `scripts/methylation_gnn_validate.py` — train + evaluate on real data
- `checkpoints/methylation_gnn/finetune_best.pt` — first real checkpoint
- `results/methylation_gnn_validation.json` — AUC vs LR baseline

**Exit criterion:** This phase requires a methylation-bioinformatician collaborator (per TEAM.md §6.2). Do not block on it; the Phase 0-2 work does not require this.

---

## 5. What's already done (no work needed)

| Component | Location | Status |
|---|---|---|
| Methylation GNN architecture | `src/methylation_gnn/gnn_model.py` (640 lines) | Done; smoke-tested only |
| Graph builder | `src/methylation_gnn/graph_builder.py` (762 lines) | Done; needs real BED files for production |
| Trainer (2-phase pretrain + finetune) | `src/methylation_gnn/gnn_trainer.py` (686 lines) | Done; smoke-tested only |
| Inference wrapper | `src/methylation_gnn/gnn_inference.py` (366 lines) | Done; needs checkpoint |
| Fusion adapter | `src/methylation_gnn/integration.py` (408 lines) | Done; not wired into headline fusion |
| Reference-data catalog | `src/methylation_gnn/data.py` (619 lines) | URL list only; no downloader |
| 46 smoke tests | `src/methylation_gnn/test_integration.py` | All synthetic; need real-data thresholds |
| Fragmentomics pipeline (cross-study baseline) | `cfdna-fragmentomics-pipeline/scripts/honest_benchmark.py` | Production; AUC 0.974-0.978 |
| 5-fold pooled OOF protocol | `cfdna-fragmentomics-pipeline/scripts/evaluate_cv.py` | Production |
| Adapter to fragmentomics output | `deepcatch/src/fragmentomics/tumor_naive_adapter.py` | Production; extends to methylation trivially |
| FinaleDB ingestion | `cfdna-fragmentomics-pipeline/scripts/fetch_finaledb.py` | Production; the WGS we already have |

The scaffolding exists. The data exists (FinaleDB WGS). What's missing is **wiring** + **FinaleMe feature extraction** + **validation**.

---

## 6. Open questions

1. **FinaleMe limitations on cancer samples:** The 0.91 auROC was reported on healthy + breast/prostate. Performance on HCC, CRC, LUAD may differ. Phase 0 will quantify this.
2. **Does FinaleMe's TOO score generalize across cancer types?** It was trained on a limited tissue panel. Pan-cancer TOO accuracy at the 627-cohort scale is unknown.
3. **Per-cancer-type methylation baseline:** Is the signal uniform across cancer types, or does methylation carry strong signal for some (HCC, CRC) and weak for others (breast, prostate)?
4. **Should we attempt the GNN scaffold (Phase 4) without a methylation expert?** The risk is getting stuck on debug cycles that a domain expert would resolve quickly.

---

## 7. Timeline & resource budget

| Phase | Wall-clock | Compute | Disk | New code |
|---|---|---|---|---|
| Phase 0 (smoke validation) | 2-3 weeks | 1 CPU-day | 5 GB | ~300 lines |
| Phase 1 (full 627 methylation baseline) | 4-6 weeks | 10 CPU-days | 20 GB | ~600 lines |
| Phase 2 (head-to-head + fusion) | 3-4 weeks | 5 CPU-days | 5 GB | ~400 lines |
| Phase 3 (TOO ablation, conditional) | 2-3 weeks | 2 CPU-days | 5 GB | ~300 lines |
| Phase 4 (GNN training, deferred) | 2-4 weeks | 50 GPU-hours (if collaborator has GPU) | 100 GB | ~200 lines |

**Total minimum path (Phases 0-2):** 9-13 weeks wall-clock, no new raw-data downloads.

**Total with Phase 4:** 11-17 weeks + a methylation-expert collaborator.

**Disk budget:** ~130 GB additional to the existing 500 GB FinaleDB cache (no raw methylation download needed; FinaleMe runs on existing WGS fragments).

**Compute budget:** ~17 CPU-days + optional 50 GPU-hours for Phase 4.

---

## 8. Success criteria

This project is successful if, at the end of Phase 2:

1. ✅ Methylation features are extracted from the full 627 cohort using FinaleMe.
2. ✅ Cross-study AUC for methylation-only baseline is reported.
3. ✅ Head-to-head comparison vs fragmentomics is published (positive OR negative).
4. ✅ Combined feature AUC is reported with proper DeLong CI and multiple-testing correction.
5. ✅ All results are reproducible from a fresh clone (following the existing `cfdna-fetch` → `cfdna-fsd` → `cfdna-delfi` → `cfdna-finaleme` workflow).
6. ✅ The honest finding is documented, regardless of direction.

This project is NOT successful if:

- ❌ Headline methylation AUC is reported without a confidence interval.
- ❌ Methylation features are tuned on the test set.
- ❌ The fusion result is reported without showing the methylation-only baseline.
- ❌ Negative results are suppressed.

---

## 9. Repo structure

This project will live in `rollroyces/deepcatch` as a new top-level module: `src/methylation/`.

Final structure:

```
src/methylation/
├── __init__.py
├── finaleme_extract.py       # WGS → methylation features (Phase 1)
├── methylation_baseline.py   # LR baseline + 5-fold CV (Phase 1)
├── methylation_vs_fragmentomics.py  # head-to-head (Phase 2)
├── methylation_fusion.py     # combined feature + fusion ablation (Phase 2)
├── too_ablation.py           # tissue-of-origin (Phase 3)
├── gnn_validate.py           # GNN on real data (Phase 4)
├── test/
│   ├── test_finaleme_extract.py
│   ├── test_methylation_baseline.py
│   ├── test_methylation_fusion.py
│   └── test_too_ablation.py
├── results/                   # JSON outputs (mirroring cfdna-fragmentomics-pipeline/results/)
└── METHYLATION_PROJECT.md     # this file (also at repo root)

scripts/                       # NEW top-level scripts (cross-repo entry points)
├── cfdna-finaleme             # CLI wrapper around finaleme_extract.py
└── cfdna-methylation-baseline # CLI wrapper around methylation_baseline.py
```

A separate `rollroyces/deepcatch-methylation` repo is **not** warranted at this stage — the project is an extension of DeepCatch, not a separate framework. Splitting repos adds CI/maintenance burden without research benefit. If the project grows past Phase 4 (the GNN training), we can revisit repo splitting.

---

## 10. Honest constraints

- **No institutional affiliation** (per user, 2026-09-10). I am working solo.
- **No methylation-bioinformatician collaborator** (yet). This blocks Phase 4 but not Phases 0-2.
- **No GPU access on this hardware.** Phase 4 will require a collaborator with GPU resources, or borrowing cloud compute.
- **FinaleMe is an imputation, not an assay.** The signal ceiling is lower than true bisulfite methylation. This is a feature-engineering project, not a "validate GRAIL's claims" project.
- **Headline numbers will not beat Galleri.** Galleri uses targeted bisulfite on >100K CpGs with thousands of patients. We have 627 WGS samples and imputed methylation. The contribution is methodological: showing what signal is extractable from a 627-sample cohort with imputation, not competing with industry.

---

## 11. References

- **Liu et al. 2020** (CCGA-1): Liu MC, et al. "Sensitive and specific multi-cancer detection with the multi-cancer early detection (MCED) test." *Annals of Oncology* 31:745-759. DOI 10.1016/j.annonc.2020.02.011
- **FinaleMe:** Liu Z, et al. (2024). "FinaleMe: predicting methylation status from cfDNA WGS." *Nature Communications* 15:2790. DOI 10.1038/s41467-024-47196-6. Code: https://github.com/epifluidlab/FinaleMe
- **FinaleDB:** Zheng H, Zhu Z, Liu Z. (2021). "FinaleDB: a browser and database of cell-free DNA fragmentation." *Bioinformatics* 37(16):2502-2503. DOI 10.1093/bioinformatics/btaa999
- **Loyfer 2023:** Loyfer N, et al. (2023). "A DNA methylation atlas of normal human cell types." *Nature* 613:585-594. DOI 10.1038/s41586-022-05580-6
- **Moss 2018:** Moss J, et al. (2018). "Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease." *Nature Communications* 9:5068. DOI 10.1038/s41467-018-07466-6
- **TCGA Methylation:** GDC Methylation Array Pipeline. https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_Pipeline/
- **cfMethDB:** Sun et al. (2025). *Genomics, Proteomics & Bioinformatics*. DOI 10.1093/gpbjnl/qzaf092

---

## 12. Status

| Item | Status |
|---|---|
| Inventory of existing methylation scaffolding | ✅ Done (see /tmp/methylation_inventory.md, 2026-09-10) |
| Inventory of public methylation data sources | ✅ Done (see /tmp/methylation_data_sources.md, 2026-09-10) |
| This project plan | ✅ Done |
| Phase 0 (FinaleMe smoke validation) | 🔲 Not started |
| Phase 1 (627-cohort methylation baseline) | 🔲 Not started |
| Phase 2 (head-to-head + fusion) | 🔲 Not started |
| Phase 3 (TOO ablation) | 🔲 Not started (conditional) |
| Phase 4 (GNN validation) | 🔲 Not started (deferred until collaborator) |

Last updated: 2026-09-10
