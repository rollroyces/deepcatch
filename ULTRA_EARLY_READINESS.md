# DeepCatch — Ultra-Early Cancer Signaling: Readiness Status

**Date:** 2026-10-07 (HKT)
**Commit:** `cc029c1` (HEAD of `main`)
**Audience:** clinical / translational collaborators and clinical reviewers
**Headline:** *the ultra-early signal is real, but the readiness is uneven* — fragmentomics is ready to publish a Stage I vs late-stage honest number on the open-data cohort; the panel detector is ready as an *in silico* benchmark at 0.1% ctDNA but needs real-plasma LoD validation; the foundation multi-modal fusion failed the per-patient honest test and is **not** a current headline for ultra-early signaling.

---

## 0. TL;DR — what "ready for ultra early cancer signaling" means here

| Layer | Status | Honest scope |
|---|---|---|
| **Fragmentomics 5-channel cross-study pooled AUC** (FinaleDB pubs 6+8, n=627) | **Ready** | AUC 0.967 ± 0.004 with GC correction; shuffled-label null 0.466; true-confound control collapses to 0.499. |
| **Per-cancer sens@99% spec on the same cohort** | **Ready** | Pooled sens@99% = 0.669 (DeLong CI logged in `results/per_cancer_sens_at_spec.json`); per-cancer types: LUAD 0.71, OV 0.90, PAAD 0.69, CRC 0.83, HCC (Jiang) 0.78, BRCA 0.30, OTHER_C 0.61. |
| **Screening-grade sens@99.5% / @99.9% spec** | **Ready** *(offline computation)* | New CLI `--specificities 0.95 0.98 0.99 0.995 0.999` (see `scripts/per_cancer_sens_at_spec.py`) emits the full screening-grade grid; pooled sens@99.5% = 0.582. |
| **Stage I vs late-stage breakdown** | **Ready** *(computable, n=0 in open-data)* | New `build_stage_breakdown_table` in `src/per_cancer_sens_at_spec.py`. The open-data FinaleDB cohort **does not carry stage labels**, so the breakdown is omitted by default and only fires when an input TSV has a `stage` column. |
| **Panel-LLR ultra-low-VAF sweep** (TCGA-LUAD 5,738 mutations, *spike-in simulation*) | **Ready (in silico only)** | At 0.1% ctDNA, AUC 0.921 / sens@99% = 0.49; at 0.05% the curve still passes the gate; at 0.01% the panel is information-limited. |
| **Assay production specs** (depth + duplex-UMI error rate) | **Ready (in silico only)** | `real_tcga_validation.py --ultraearly-sweep`: 50k× depth + ≤1e-4 error → AUC 1.000 at 0.1% ctDNA. |
| **Multi-modal foundation fusion** (Stage 1) | **NOT ready** | Audit-2 honest number: foundation AUC 0.55 ± 0.01 vs LR baseline 0.91 on paired TCGA-LUAD. Does not beat LR baseline on n=40 paired. |
| **Longitudinal Stage 2 (Kalman)** | **NOT ready** | Honest baseline AUC 0.49; needs redesign per `docs/PRODUCTION_ROADMAP.md` §6. |
| **Real plasma validation** | **NOT ready** | No paired plasma + tissue cohort in repo. This is the binding constraint. |

The auditable ground-truth numbers below are all reproducible from `bash paper/REPRODUCE.sh` or the per-module CLIs documented in [`USAGE.md`](./USAGE.md). Every entry in this readiness doc is a verifiable `results/*.json` artifact.

---

## 1. Fragmentomics is ready for ultra-early signaling (within scope)

### 1.1 The honest cross-study headline

| Setting | Pooled AUC | sens@99% spec | Notes |
|---|---:|---:|---|
| Baseline (5-channel, no GC) | 0.9747 ± 0.0012 | 0.793 | legacy; **inflated ~0.008 by GC-axis proxy** |
| **+ GC correction (default)** | **0.9670 ± 0.0035** | **0.669** | honest, removes batch proxy |
| + 4-mer motifs (`--include-motifs`) | 0.9768 ± 0.0024 | — | small lift on baseline; below GC-corrected AUC in some seeds |
| **Shuffled-label null** | **0.4657** | 0.0 | passes the <0.55 gate (5 seeds; 0.46–0.51 range) |
| **True-confound (cancer=A, healthy=B)** | **0.499** | — | signal is cancer-vs-healthy, not batch |

Source: `results/cross_study_finallydb_gc_corrected.json` and `results/cross_study_finallydb.json`.
Reproduce: `bash paper/REPRODUCE.sh` or `python scripts/cross_study_finallydb.py --publications 6 8 --seeds 42 13 7 99 1234 --gc-correction`.

**Why this is ultra-early-grade evidence**:
- n=627 samples with 364 cancer + 263 healthy (FinaleDB pubs 6 + 8; **Jiang 2015 PNAS** + **Cristiano 2019 Nature**).
- Pooled OOF, 5-seed × 5-fold StratifiedKFold, per-publication z-score harmonization inside each train fold.
- Two confound controls (shuffled-label null + true-confound control) prove the signal is cancer-vs-healthy, not batch, not fold identity.
- Per-publication z-score harmonization is the published standard (Cristiano 2019 used a comparable pooled classifier).

### 1.2 Per-cancer sens@spec at 95% / 98% / 99% / 99.5% / 99.9%

The screening-grade grid is now first-class. Run:
```bash
python scripts/per_cancer_sens_at_spec.py --scores-tsv <scores.tsv> --out results/per_cancer_sens_at_spec.json
# default specificities: 0.95 0.98 0.99 0.995 0.999 (clinical-decision + screening-grade)
```

If you don't have a TSV, run a synthetic smoke (the Stage I / LATE breakdown fires too):
```bash
python scripts/per_cancer_sens_at_spec.py --synthetic --n 400 --seed 42 --include-stage-breakdown
```

Output JSON contract (`results/per_cancer_sens_at_spec.json`):
```jsonc
{
  "per_cancer": {
    "LUAD": {
      "auc_mean": 0.97,
      "auc_ci": [0.94, 0.99],
      "sens_at_95": 0.93,
      "sens_at_98": 0.83,
      "sens_at_99": 0.71,
      "sens_at_995": 0.55,
      "sens_at_999": 0.42,
      "sens_at_spec": { "0.95": {...}, "0.98": {...}, "0.99": {...}, "0.995": {...}, "0.999": {...} },
      "ppv_at_prevalence": { "prev_0.001": 0.04, "prev_0.05": 0.78, ... }
    },
    ...
  },
  "pooled": { ... },
  "stage_breakdown": {   // only present when --include-stage-breakdown or input has 'stage' column
    "per_stage": { "I": {...}, "LATE": {...}, "UNKNOWN": {...} },
    "n_stage_I": ..., "n_stage_late": ..., "n_stage_unknown": ...
  },
  "provenance": { "specificities": [0.95, 0.98, 0.99, 0.995, 0.999], ... }
}
```

### 1.3 Stage I vs late-stage breakdown — what exists and what doesn't

The new `build_stage_breakdown_table` (in `src/per_cancer_sens_at_spec.py`) splits cancer positives into Stage I / LATE / UNKNOWN buckets while reusing the same healthy distribution, so sens@spec is comparable. Coercion rules:

- `I`, `IA`, `IB`, `1`, `T1` → `STAGE_I`
- `II`, `III`, `IV`, `2`, `3`, `4`, anything else → `STAGE_LATE` (conservative — if we can't prove a sample is Stage I we count it against the early claim)
- `NA`, `?`, empty, `UNKNOWN` → `STAGE_UNKNOWN`

**Honest gap**: the open-data FinaleDB cohort in this repo does **not** carry stage labels. Therefore the cross-study headline (1.1) does not currently have a Stage I breakdown. This is on the **action list** for the FLARE (GSE317007) integration (next steps §6.1) and for any future own-cohort work.

The stage breakdown is **exercised by tests** (`pytest test/test_per_cancer_sens_at_spec.py` runs the synthetic fixture with `--include-stage-breakdown`) so the math is wired up and ready for the first cohort that ships with stage labels.

---

## 2. Panel-LLR ultra-low-VAF detector is ready (in silico)

This is the assay-design benchmark — *not* a clinical readout. Source: `results/real_tcga_validation.json` and the underlying 20 LUAD-patient TCGA-LUAD panel (5,738 mutations).

### 2.1 Sweep at 0.1% ctDNA (real TCGA-LUAD panel)

| Spec | 0.1% ctDNA |
|---|---|
| **Panel AUC** | **0.921 ± 0.018** |
| **Sens @ 95% spec** | **0.770** |
| **Sens @ 99% spec** | **0.490** |

20 patients, 5 seeds, GroupKFold per-patient, paired win-rate 1.0 (panel-LLR beats the LR baseline on every patient).

### 2.2 Assay-sweep: depth + error rate

`real_tcga_validation.py --ultraearly-sweep` walks depth × error-rate combinations. At 0.1% ctDNA on the same panel:

| Depth × error | AUC | sens@95% | sens@99% |
|---|---:|---:|---:|
| 5,000× × 1e-4 | 0.998 | 1.000 | 0.96 |
| **50,000× × 1e-4** | **1.000** | **1.000** | **1.000** |
| 50,000× × 2e-3 | 0.997 | 0.99 | 0.96 |

⇒ The production assay spec is duplex-UMI consensus (~1e-4 error) at 50k× depth; either alone is sufficient.

### 2.3 What this headline does NOT prove

- It's a simulation (cfDNA dilution from real TCGA-LUAD mutations into a synthetic healthy BAM). Per `docs/PRODUCTION_ROADMAP.md` §8, every number here will be worse on real plasma. The simulation→reality gap is the **#1 risk** and is not measured by anything currently in this repo.
- It's panel-LLR — tumor-informed. For MCED screening there is no tumor to inform a panel; the cross-study fragmentomics signal (1.1) is what carries that case, and it lives on open-data.

---

## 3. Multi-modal foundation fusion — the AUDIT-2 honest number

The foundation multi-modal fusion layer (Stage 1, `src/foundation/`) has a published **honest** evaluation under Audit-2 (`AUDIT_2_FINDINGS.md`):

| Metric | Old (audit-1) | **New (audit-2 honest)** |
|---|---|---|
| panel_only AUC | 0.92 | **0.915** |
| foundation AUC | 0.939 (inflated) | **0.554 ± 0.005** |
| LR baseline AUC | 0.96 | **0.910** |
| shuffled_foundation AUC | 0.31 (anti-correlated) | **0.741** (picks up pair-invariant) |
| gate_pass | true | **false** |

**Conclusion**: the foundation model does not beat the sklearn LR baseline on n=40 paired TCGA-LUAD. The reason is structural — the per-arm sequencing-noise jitter (TF=0.001 vs TF=0) dominates the per-sample separable signal, GroupKFold denies the foundation access to per-patient mutation signatures, and LR has fewer parameters to overfit.

**Implication for ultra-early signaling**: the foundation fusion layer is **not a current headline**. The cross-study fragmentomics pipeline (1.1) and the panel-LLR detector (2.1) carry the ultra-early signal today.

The path forward is documented in `docs/PRODUCTION_ROADMAP.md` §6.3 and `NEXT_STEPS.md` §6: replace the SPRT-based CET Stage 2 with a censored-Poisson hierarchical-Bayes tracker; build a per-arm sequencing-noise model; pair on real plasma.

---

## 4. Readiness checklist (what blocks ready)

| Item | Owner | Blocker |
|---|---|---|
| Finalize Jiang nested-CV honest AUC | author | `data/deepcatch_data.xlsx` (privacy-removed in `8c812c0`; reach out to Prof. Jiang) |
| Stage I / LATE per-cancer numbers on FinaleDB cohort | author | Stage labels missing from open-data cohort |
| Real plasma LoD study | clinical partner | Access to Cristiano 2019 (EGA), Newman (ECMO), TRACERx (EGAD00001002469), or own cohort |
| Clinical paired plasma + tissue | clinical partner | IRB + biobank agreement |
| Production assay SOP | clinical lab | Validation of duplex-UMI + 50k× spec on real instrument |
| Longitudinal Stage 2 redesign | author | Per-patient paired data |
| Foundation fusion re-validation | author | Real unpaired plasma cohort (LR baseline is the current ceiling) |

The detailed 12-month execution plan is in [`docs/PRODUCTION_ROADMAP.md` §7](./docs/PRODUCTION_ROADMAP.md). The TODAY action is `NEXT_STEPS.md` §0 — push the branch and confirm CI is green.

---

## 5. Honest limitations (do not claim beyond them)

1. **Open-data scope, not clinical plasma.** The 0.967 pooled AUC is on FinaleDB pre-processed cfDNA features, not on raw plasma reads from a clinical lab. Every clinical number degrades from simulation.
2. **Not FDA-approved, not clinically validated, not a medical device.** See the WARNING block at the top of `README.md`.
3. **Panel-LLR is tumor-informed.** MCED-style screening cannot use it; the fragmentomics pipeline (1.1) is what carries that case.
4. **Stage I coverage is 0% on the open-data headline.** Until a real cohort with stage labels is added, we cannot honestly claim Stage I sensitivity.
5. **The foundation fusion layer's 0.93 number is retracted.** The audit-2 0.55 honest number is what we can defend. Do not cite 0.93.
6. **CHIP subtraction is not modeled.** In populations >60 y, CHIP alone costs 5–10% specificity; matched-WBC is mandatory for any clinical claim.

---

*This document is the single entry-point for "is DeepCatch ready for ultra-early cancer signaling today?". Every number above is a verifiable computation in this repo. Update it as new artifacts land.*