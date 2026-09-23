# Methods Paper Skeleton — Status Meta-Doc

**File:** `paper/METHODS_PAPER.md`
**Date:** 2026-09-23 (HKT, UTC+08:00)
**Commit:** `9f6662f` (HEAD of `main`)
**Status:** DRAFT — methods-paper draft with grounded §3 numbers.

---

## What's in this draft

| Section | Status | Word count (approx) | Notes |
|---|---|---:|---|
| Title + authors + venue target + abstract | ✅ drafted | ~440 | abstract now reports the §3.2 pooled AUC + §3.3 per-cancer OV + true-confound numerics |
| 1. Introduction (motivation, related work, gap) | ✅ drafted | ~390 | cites DELFI, CAPP-Seq, Galleri, FinaleMe, FinaleDB, Jiang 2020 |
| 2.1 Fragmentomics features | ✅ drafted | ~140 | six-modality schema |
| 2.2 Multi-modal fusion architecture | ✅ drafted | ~210 | foundation encoder + cross-attention + LR head |
| 2.3 Panel-based detection | ✅ drafted | ~80 | pointer to companion paper |
| 2.4 Validation protocol | ✅ drafted | ~270 | 5-seed × 5-fold GroupKFold + shuffled-label + true-confound + DeLong CIs + bit-identical sanity |
| 3.1 Honest baseline AUC on paired synthetic cohort | ✅ filled | ~220 | n=20, AUC ≈ 0.55 honest under GroupKFold |
| 3.2 Cross-study AUC on FinaleDB pub 6+8 | ✅ filled | ~430 | pooled 0.974 ± 0.002, sens@99=0.782; Jiang alone 0.979, Cristiano alone 0.969; true-confound table (0.49 vs 1.00) |
| 3.3 Per-cancer sens@spec (DeLong CIs + PPV@prev) | ✅ filled | ~470 | top-5 cancers, DeLong CIs at 95/98/99% specificity; prevalence-floor PPV summary table |
| 3.4 Honest null/negative ablations | ✅ filled | ~190 | focal-BCE null, sparse-aware negative/marginal, harmonization re-classified MIXED/PROTECTIVE on the cross-study confound control |
| 3.5 Comparison to published cfDNA benchmarks | ✅ filled (NEW) | ~430 | DELFI, CancerSEEK, Galleri, CAPP-Seq, FinaleMe comparator table + honest-positioning paragraph |
| 4.1 What this work does | ✅ drafted | ~80 | engineering substrate, not clinical assay |
| 4.2 What this work does NOT do | ✅ drafted (touched) | ~230 | added HCC_J outlier reading, per-cancer denominator limits, cross-study framing |
| 4.3 Limitations | ✅ drafted | ~260 | small n, simulated healthy, no real plasma, FinaleDB API down, encoder size |
| 4.4 Future work | ✅ drafted | ~150 | real-plasma IRB cohort, cross-platform, methylation, larger pretrain |
| 5. Availability | ✅ drafted | ~210 | github URL, license, data pointers, one-bash reproduction pointer |
| `paper/REPRODUCE.md` | ✅ NEW | ~190 | one-bash driver + step-by-step manual commands + acceptance checklist |
| `paper/REPRODUCE.sh` | ✅ NEW (executable) | (script) | runs cross_study_finallydb.py + pytest gate; accepts `--quick`, `--no-test` |

**Total: ~4,470 words** (well under the 5,000-word ceiling; measured by `wc -w paper/METHODS_PAPER.md` = 4,428 — `wc -w` under-counts markdown table cells).

---

## What's filled (vs the previous skeleton)

- **§3.2** — now has the per-cohort rows (Jiang 0.9791 ± 0.0028, Cristiano 0.9693 ± 0.0022), pooled rows (harmonized 0.9738 ± 0.0018, no_harmonize 0.9650 ± 0.0019), sens@spec at 95/98/99%, and the four-row **true-confound control** table. Source: `results/cross_study_finallydb.json`. Companion narrative: `docs/CROSS_STUDY_BENCHMARK.md`.
- **§3.3** — top-5 cancer OvR table (HCC_J, LUAD, BRCA, PAAD, OV) with AUC ± std, DeLong 95% CI on AUC, sens@95/98/99% with DeLong CI on sens@99, n_cancer + n_total. Followed by a **prevalence-floor PPV summary table** at spec=99% across the standard cfDNA grid {0.001, 0.004, 0.01, 0.05, 0.10, 0.20, 0.50}. Source: `results/per_cancer_sens_at_spec.json` (DeLong CIs; cross-references `results/cross_study_finallydb.json` for bootstrap CIs and for the cohort-size accounting).
- **§3.5** (NEW) — comparator table for DELFI, CancerSEEK, Galleri, CAPP-Seq, FinaleMe, and DeepCatch — plus a three-point honest-positioning paragraph. **Not a clinical claim.**
- **§4.2** — reclassified "harmonization" limitation from "NEGATIVE: signal is partly study-confounded" to "MIXED, but the true-confound control proves harmonization is **PROTECTIVE** for the cross-study claim".
- **§5** — added the `paper/REPRODUCE.sh` pointer in addition to `RUN_ALL.sh`.

## What's still a stub (or TODO)

- **Author list** (`[co-authors: TBD]` placeholder). Yu Ching Lam is the lead; collaborators TBD.
- **Final venue** between PLOS Computational Biology and Bioinformatics Advances (both plausible).
- **Figure 1** — the multi-modal encoder architecture diagram. Source: `paper/PAPER.md` Figure 1 (panel-LLR scoring) is already there; this paper needs a different figure (encoder + cross-attention + pre-training phases).
- **Figure 2** — the validation protocol schematic (5-seed × 5-fold GroupKFold + shuffled-label control).
- **Table 1** — the six-modality schema with concrete per-modality dim and example feature (currently inline in §2.1).
- **LaTeX port** (`paper/METHODS_PAPER.tex`) for the venue submission (Bioinformatics Advances requires LaTeX; PLOS Comp Bio accepts either).

## What's intentionally NOT in this draft

- **Clinical operating points.** sens@99 = 0.782 is quoted honestly as a pooled-internal-CV number; it is not extrapolated to a held-out clinical cohort, MCED-population prevalence, or FDA pathway. The §3.5 paragraph explicitly rules out a head-to-head reading.
- **A pretraining-cohort AUC.** The §3.2 number is on the cross-study pool, not on the 200-sample pretraining cohort. The pretraining cohort's purpose is the encoder (MODALITY-token + masked-prediction), not the headline AUC — see `docs/PRETRAINING.md`.
- **Paired-design ablations of §3.4.** Those run under their own drivers (`docs/SENS_AT_SPEC_ABLATION.md`, `docs/SPARSE_AWARE_ABLATION.md`); the paper quotes the verdict, not the full setup.

## Cross-references and provenance

This draft synthesized the following sources:

- `docs/CROSS_STUDY_BENCHMARK.md` — the running cross-study narrative (already shipped by sibling subagent on commit `9f6662f`).
- `results/cross_study_finallydb.json` — the JSON artifact behind §3.2 (pooled, per_cohort, true_confound_control, per_cancer).
- `results/per_cancer_sens_at_spec.json` — the JSON artifact behind §3.3 (DeLong CIs + PPV@prev).
- `paper/METHODS_PAPER_CITATIONS.bib` — references for [5], [6], [7], [9], [10], [11], [16], DeLong 1988.
- `docs/HARMONIZATION.md` — backs the §3.4 harmonization re-classification.
- `docs/SENS_AT_SPEC_ABLATION.md`, `docs/SPARSE_AWARE_ABLATION.md` — back the §3.4 row.
- `docs/PRETRAINING.md`, `docs/PRETRAIN_BUG.md` — back §2.2.
- `AUDIT_2_FINDINGS.md` (repo root) — backs §3.1's "honest framing" paragraph.
- `MODEL.md` — backs §2.2 PRODUCTION_CONFIG.

---

## Verification commands

The parent task asks that the following pytest command remain green after any
edit:

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

This draft does not modify any Python source — only adds `paper/REPRODUCE.md`,
`paper/REPRODUCE.sh`, and updates the markdown docs in `paper/`. The test
suite remains at its current passing count.

The paper-text reproduction is now the single-command driver
`bash paper/REPRODUCE.sh`; the per-section JSON mapping is documented in
`paper/REPRODUCE.md` §"What each artifact is the source of truth for".

---

## Remaining TODOs (post-draft)

- [ ] Confirm co-authors (or finalize single-author manuscript).
- [ ] Decide final venue between PLOS Computational Biology and Bioinformatics Advances.
- [ ] Add Figure 1 (encoder architecture).
- [ ] Add Figure 2 (validation-protocol schematic).
- [ ] Add Table 1 (six-modality schema).
- [ ] Run the LaTeX port (`paper/METHODS_PAPER.tex`) once Figure 1/2 + Table 1 land.
- [ ] Cross-link `paper/REPRODUCE.sh` from `RUN_ALL.sh` so a single top-level
      driver regenerates everything.

---

*Last updated: 2026-09-23.*
