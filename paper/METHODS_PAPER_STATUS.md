# Methods Paper Skeleton — Status Meta-Doc

**File:** `paper/METHODS_PAPER.md`
**Date:** 2026-09-23 (HKT, UTC+08:00)
**Commit:** `fb1e523` (HEAD of `main`)
**Status:** SKETCH — methods-paper skeleton, NOT a finished manuscript

---

## What's in this skeleton

| Section | Status | Word count (approx) | Notes |
|---|---|---:|---|
| Title + authors + venue target + abstract | ✅ drafted | ~360 | venue = PLOS Comp Bio / Bioinformatics Advances |
| 1. Introduction (motivation, related work, gap) | ✅ drafted | ~390 | cites DELFI, CAPP-Seq, Galleri, FinaleMe, FinaleDB, Jiang 2020 |
| 2.1 Fragmentomics features | ✅ drafted | ~140 | six-modality schema |
| 2.2 Multi-modal fusion architecture | ✅ drafted | ~210 | foundation encoder + cross-attention + LR head |
| 2.3 Panel-based detection | ✅ drafted | ~80 | pointer to companion paper |
| 2.4 Validation protocol | ✅ drafted | ~270 | 5-seed × 5-fold GroupKFold + shuffled-label + true-confound + DeLong CIs + bit-identical sanity |
| 3.1 Honest baseline AUC on paired synthetic cohort | ✅ filled | ~220 | n=20, AUC ≈ 0.55 honest under GroupKFold; reads sensibly |
| 3.2 Cross-study AUC on FinaleDB pub 6+8 | 🟡 TODO (pointer only) | ~80 | waits on `docs/CROSS_STUDY_BENCHMARK.md` |
| 3.3 Per-cancer sens@spec | 🟡 TODO (pointer only) | ~50 | waits on `results/per_cancer_sens_at_spec/` |
| 3.4 Honest null/negative ablations | ✅ filled | ~150 | focal-BCE null, sparse-aware negative/marginal, harmonization negative |
| 4.1 What this work does | ✅ drafted | ~80 | engineering substrate, not clinical assay |
| 4.2 What this work does NOT do | ✅ drafted | ~150 | clinical validation, FDA pathway, prospective screening |
| 4.3 Limitations | ✅ drafted | ~260 | small n, simulated healthy, no real plasma, FinaleDB API down, encoder size |
| 4.4 Future work | ✅ drafted | ~150 | real-plasma IRB cohort, cross-platform, methylation, larger pretrain |
| 5. Availability | ✅ drafted | ~200 | github URL, license, data pointers, reproduce command |

**Total: ~2,800 words** (well under the 5,000-word ceiling).

---

## What's a stub (placeholder, not yet content)

- **§3.2 Cross-study AUC table** — the cross-study sweep hasn't landed yet; section is a pointer to the forthcoming `docs/CROSS_STUDY_BENCHMARK.md`. Once that sibling subagent's work lands, paste the per-cancer AUC table here with DeLong CIs.
- **§3.3 Per-cancer sens@spec** — depends on `results/per_cancer_sens_at_spec/*.json`. Once those JSONs exist, fill in the per-cancer sensitivity at 95% / 99% specificity with DeLong CIs.
- **§2.3 Panel-based detection** — pointer to `paper/PAPER.md` (companion benchmark paper) is intentional; this methods paper is about fragmentomics, not panel scoring.

## What's TODO (post-skeleton)

- [ ] Fill in §3.2 once `docs/CROSS_STUDY_BENCHMARK.md` is written.
- [ ] Fill in §3.3 once `results/per_cancer_sens_at_spec/*.json` are produced.
- [ ] Author list (`[co-authors: TBD]` placeholder). Yu Ching Lam is the lead; collaborators TBD.
- [ ] Decide final venue between PLOS Computational Biology and Bioinformatics Advances (both are plausible; PLOS Comp Bio is broader scope, Bioinformatics Advances is methods-focused).
- [ ] Add a Figure 1 — the multi-modal encoder architecture diagram. Source: `paper/PAPER.md` Figure 1 (panel-LLR scoring) is already there; this paper needs a different figure (encoder + cross-attention + pre-training phases).
- [ ] Add a Figure 2 — the validation protocol schematic (5-seed × 5-fold GroupKFold + shuffled-label control).
- [ ] Add a Table 1 — the six-modality schema with concrete per-modality dim and example feature.
- [ ] Run the LaTeX port (`paper/METHODS_PAPER.tex`) for the venue submission if the venue requires LaTeX (Bioinformatics Advances does; PLOS Comp Bio accepts either).

---

## Cross-references and provenance

This skeleton was built by reading the following documents and synthesizing their honest framing:

- `docs/PRETRAINING.md` — pre-training pipeline (PRODUCTION_CONFIG, 200-sample real FinaleDB, Phase 1 + Phase 2 losses)
- `docs/PRETRAIN_BUG.md` (referenced from PRETRAINING.md) — synthetic-bypass bug + regression test
- `docs/SPARSE_AWARE_ABLATION.md` — null/negative ablation on sparse-aware projection
- `docs/SENS_AT_SPEC_ABLATION.md` — null ablation on focal-BCE loss (supersedes n=5 with n=20 honest paired test)
- `AUDIT_2_FINDINGS.md` (repo root) — consolidated audit-2 with P0-A (CV leak), P0-B (synthetic frag replacement), P0-C (shuffled-label control bug)
- `paper/PAPER.md` (companion bioRxiv benchmark paper) — different paper, different venue, different scope (panel-LLR MRD vs fragmentomics methods)
- `paper/paper.tex` — LaTeX source of the companion paper; NOT modified by this skeleton
- `paper/BIORXIV_SUBMISSION.md` — submission package for the companion paper; NOT modified
- `results/sens_at_spec_ce_n20.json`, `results/sparse_aware_*.json`, `results/harmonization_check.json`, `results/pretrain_production_finaledb.json` — JSON artifacts cited in §3

The methods paper is a **different paper target** from the bioRxiv submission. The bioRxiv submission (paper/PAPER.md, paper/paper.tex) is the **benchmark paper** — open-source panel-LLR MRD pipeline, AUC 0.921 at 0.1% ctDNA on real TCGA-LUAD mutations. The methods paper (paper/METHODS_PAPER.md, this file) is the **fragmentomics framework paper** — engineering substrate for cfDNA fragmentomics research, honest negative ablations, no clinical claims.

---

## Verification commands

The parent task asks that the following pytest command remain green:

```bash
cd /Users/hermes/deepcatch
env -u PYTHONPATH ./.venv/bin/python -m pytest \
    test/test_biomedical_review_fixes.py \
    test/test_sparse_aware_projection.py \
    test/test_finaledb_pretrained_loader.py \
    test/test_pretrain_bug_fix.py \
    test/test_harmonization_check.py \
    test/test_per_cancer_sens_at_spec.py \
    src/foundation/test_integration.py \
    --timeout=60 -q
```

This skeleton does not modify any Python source — only adds three new markdown / bib files under `paper/`. The test suite should remain at its current passing count.

---

*Last updated: 2026-09-23.*
