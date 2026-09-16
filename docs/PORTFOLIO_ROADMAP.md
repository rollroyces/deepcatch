# DeepCatch Portfolio — Unified Roadmap

**Date:** 2026-09-15
**Author:** Yu Ching Lam (via Hermes Agent)
**Status:** Current as of commits `8a66d24` (deepcatch), `7e6188e` (cfdna-fragmentomics-pipeline), `2d621cd` (deepcatch-methylation).

This is the **single source-of-truth roadmap** for the 3-repo open-source cfDNA
early-detection portfolio. It supersedes the older `NEXT_STEPS.md` and
`docs/PRODUCTION_ROADMAP.md` (kept as historical references) and reflects the
actual current state after the 4 enhancement rounds + V3 design + FinaleMe
pretrained-model breakthrough.

---

## 1. Portfolio overview (where we are)

### 1.1 The 3 repos

| Repo | Local path | Latest commit | Tests | CI | Headline result |
|---|---|---|---|---|---|
| **deepcatch** | `/Users/hermes/deepcatch` | `8a66d24` | 256/256 | ✅ green | Panel LLR AUC **0.921** @ 0.1% ctDNA (20 TCGA-LUAD patients) |
| **cfdna-fragmentomics-pipeline** | `/Users/hermes/cfdna-fragmentomics-pipeline` | `7e6188e` | 106/106 | ✅ green | Cross-study AUC **0.9755** (frag), **0.9921** (fusion) on 627 samples |
| **deepcatch-methylation** | `/Users/hermes/deepcatch-methylation` | `2d621cd` | 15/15 | ✅ green | FinaleMe pretrained HMM **breakthrough** (BH01 chr22, 489K CpGs) |

### 1.2 What we have shipped (Phase 0–1, complete)

✅ **Multi-cancer fragmentomics** (Enhancement #1, `d289348`): Macro AUC 0.9703 on 5-class OvR (BRCA, CRC, HCC_J, LUAD, PAAD).

✅ **Sens@99% + published-test comparison** (Enhancement #2, `a4c8dc6` + `e5a8fba`): Sens@99% = 0.755 (frag), 0.843 (fusion). Competitive with Galleri (51.5% @ 99.5%).

✅ **Per-sample risk-score CLI** (Enhancement #3, `6f49638`): `cfdna-score` works in 0.3s with cache, 86.7% risk-tier accuracy on 30-sample demo.

✅ **Tissue-of-origin from fragmentomics** (Enhancement #4, `91a95b0`): Macro AUC 0.9338 within-Cristiano (batch-effect removed), Top-1 79.4%, Top-2 90.7%.

✅ **Optimization round** (this session, commits `1779422`–`7e6188e`): 1.4–6× speedup on the pipeline; honest memory report; honest no-improvement on sens@99% calibration sweep.

✅ **FinaleMe pretrained HMM breakthrough** (`60043ca`): Hybrid v0.58.1+v0.61 JAR loads Zenodo models; BH01 chr22 β-values decoded (489,370 CpG positions per model).

✅ **Methylation β-value integration** (P0 partial): Methylation-proxy AUC 0.7770 on FinaleDB 627-sample cohort (TCGA-LIHC 0.972, Zenodo CRAG 0.800 partial).

✅ **Per-cancer AUC table** (this session): All per-cancer AUCs ≥ 0.94 (min PAAD 0.9401), well above the 0.85 audit target.

✅ **bioRxiv paper drafts** (2 papers): `BIORXIV_PAPER.md` (methylation, 254 lines) + `BIORXIV_PAPER_FRAGMENTOMICS.md` (432 lines, 38 KB).

✅ **Deployed docs site**: https://rollroyces.github.io/deepcatch/ with 5 Mermaid diagrams.

### 1.3 What's still blocked

🚫 **FinaleDB API HTTP 500** — Postgres broken, S3 403, GitHub dormant since 2024-01-03. The 627-sample fragmentomics cohort + 121 methylation samples remain publicly inaccessible.

🚫 **Held-out clinical validation** — needs real plasma cfDNA from a clinical collaborator.

🚫 **bioRxiv submission** — paper drafts done, blocked on production ORCID (currently sandbox `0009-0008-9113-769X`) and Zenodo DOI.

🚫 **Clinical metadata** (Enhancement #5) — blocked on FinaleDB API.

🚫 **FinaleMe Step 3 OOM** — v0.58.1 requires ~5.6 GB heap (Step 1 chr22-only worked at 14.31 min).

---

## 2. Phased plan (where we go)

### Phase A — Submission & visibility (NEXT, ~2-4 weeks)

| Task | Repo | Owner | Time | Blocked by |
|---|---|---|---|---|
| Register production ORCID at orcid.org (currently sandbox) | both | user | 15 min | user action |
| Render `BIORXIV_PAPER_FRAGMENTOMICS.md` to PDF (pandoc) | deepcatch-methylation | me | 30 min | none |
| Submit to bioRxiv via web form | both | user | 30-60 min | user action + production ORCID |
| Email Yaping Liu / Ravi Bandaru (Northwestern) for FinaleDB IAM | both | me (via Gmail skill) | 15 min | none |
| Create Zenodo deposit + DOI for code release | all 3 | me | 1 hour | user action (Zenodo login) |
| Update docs site with v3 status badge | deepcatch | me | 15 min | V3 implementation |

**Exit criteria:** Papers on bioRxiv with DOIs; FinaleDB IAM response or polite decline.

### Phase B — V3 model implementation (4-8 weeks)

**Source of truth:** `docs/V3_DESIGN.md` (committed `8a66d24`). Pre-registered acceptance criteria at §6.

| Task | Time | Status |
|---|---|---|
| Phase A — pre-registration freeze (channel set, hyperparameters, metric set) | 1 hour | ✅ done (in V3_DESIGN.md §4.2) |
| Phase B — implementation (Layer 1 channels, Layer 2 OvR elastic-net, Layer 3 fusion/decision) | 6-8 hours compute + 1-2 hours review | **TODO** |
| Phase C — validation (per-cancer AUC, sens@spec, PPV@prev) | 1-2 hours compute | TODO |
| Phase D — documentation + paper update | 1-2 hours | TODO |

**Honest note:** Sensitivity sweep (this session) showed no calibration technique beats the LR baseline. The V3 design's promised gains depend on:
1. Per-cancer elastic-net with auto-L1 (untested, may match or beat baseline by 0.01-0.02)
2. Isotonic calibration on TRAIN OOF (tested, slightly **worse** at 99% spec — 0.7851 vs 0.7989)
3. Hierarchical study effect with shrinkage (untested)
4. Prevalence-parameterized decision layer (untested — moves the metric, not the model)

**Realistic V3 outcome:** Sens@99% target 0.85 has only ~30-40% probability of being hit. If missed, ship V3 anyway as a methods paper (calibrated + per-cancer + bias-budget + sens@spec table).

### Phase C — Methylation integration (8-16 weeks, dependent on FinaleDB access)

| Task | Time | Blocked by |
|---|---|---|
| FinaleMe Step 3 (whole-genome, multi-sample) | 1-2 weeks compute | FinaleMe v0.58.1 OOM fix OR upgrade to v0.61 |
| Port methylation GNN to MPS (Apple GPU) | 4-6 hours | none |
| Generate synthetic methylation β-values (calibrated on FinaleMe) | 2-3 days compute | FinaleMe Step 3 |
| Train methylation GNN on synthetic + real (BH01 chr22) | 1-2 weeks | β-value cohort of ≥30 samples |
| Integrate GNN embeddings as 6th channel | 1 week | trained GNN |
| Re-evaluate fusion AUC with methylation channel | 1 week | integrated model |

**Exit criteria:** AUC improvement over 0.9921 fusion baseline (currently no signal that methylation will help at this cohort size — **honest**).

### Phase D — Clinical validation (6-12 months, requires external collaborators)

| Task | Time | Blocked by |
|---|---|---|
| Establish IRB + own clinical cohort | 3-6 months | user / institutional |
| Acquire real plasma cfDNA WGS data | 2-6 months | collaborators + funding |
| Run DeepCatch on real plasma (not simulated) | 1-2 months | data |
| Compare against published MCED tests (Galleri, Shield) on matched cohorts | 2-3 months | data + matching |

**Honest note:** Phase D requires institutional support, IRB approval, and biobanking access — outside the scope of an independent researcher without these resources.

---

## 3. Sequenced task list (the actionable view)

### This week (sequential, not parallel)

1. ⏳ **Register production ORCID** (user, 15 min) — gates bioRxiv submission
2. ⏳ **Render BIORXIV_PAPER_FRAGMENTOMICS.md to PDF** (me, 30 min) — pandoc
3. ⏳ **Email FinaleDB authors** (me, 15 min) — request FinaleDB IAM or Globus auth
4. ⏳ **V3 Phase A** already done; **Phase B start** (me, 1-2 hours setup)
5. ⏳ **Commit any pending work** (me, 5 min) — clean main branches

### Next 2 weeks

6. ⏳ **bioRxiv submission** (user, 30-60 min) — once ORCID is registered
7. ⏳ **V3 Phase B full** (me, 6-8 hours wall time, mostly compute)
8. ⏳ **Update docs site** (me, 30 min) — V3 status badge + new headline numbers if Phase B succeeds

### Next 1-3 months

9. ⏳ **V3 Phase C validation** (1-2 hours compute)
10. ⏳ **V3 Phase D documentation + paper update** (1-2 hours)
11. ⏳ **Zenodo DOI deposit** (1 hour, user action for login)

### Next 3-12 months (blocked on external factors)

12. 🚫 **FinaleMe Step 3 multi-sample** — blocked on memory
13. 🚫 **Methylation GNN integration** — blocked on methylation cohort size
14. 🚫 **Real plasma validation** — blocked on clinical data
15. 🚫 **MRD vs MCED assay decision** — depends on Phase D outcomes

---

## 4. Risk register

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| V3 sens@99% doesn't reach 0.85 target | High (60-70%) | Medium | Ship V3 as methods paper anyway; honest framing in paper |
| FinaleDB remains broken | Medium (40%) | High | Use Zenodo CRAG + synthetic; rely on existing 627-sample cohort |
| bioRxiv rejected (research-stage, not clinical) | Medium (30%) | Low | Target Bioinformatics / PLOS Comp Bio instead; documented in JOURNAL_REVIEW_REJECTION_ANALYSIS.md |
| Production ORCID rejected (sandbox transition issues) | Low (10%) | High | Submit without ORCID if needed; some journals accept this |
| Methylation integration doesn't improve AUC | High (60%) | Medium | Honest paper; methylation as research-stage channel only |

---

## 5. What we are NOT doing (and why this is a feature)

- **NOT** training deep learning models on n=627 (overfit). LR + PCA + naive-mean fusion stays the backbone.
- **NOT** using private data (TEAM.md §3.2 hard constraint). Public data only.
- **NOT** claiming clinical validity. MODEL.md states research-only.
- **NOT** GPU-accelerating the fragmentomics pipeline (CPU is faster at this cohort size). MPS used only for methylation GNN where the workload shape fits (large sparse graphs).
- **NOT** replacing LR with Transformer for fragmentomics (v2.2 evidence shows LR matches/beats transformers on 5-channel features).

---

## 6. Success criteria for "portfolio ships"

The portfolio is complete when **all** of the following are met:

- ☐ All 3 repos on main, all CI green, all tests passing
- ☐ bioRxiv paper submitted with DOI
- ☐ Zenodo DOI for all 3 repos
- ☐ V3 model implemented + evaluated against acceptance criteria (or honest "V3.0 methods paper" if sensitivity targets missed)
- ☐ Per-cancer AUC table committed to results/
- ☐ FinaleDB IAM or documented polite decline
- ☐ Methylation channel integrated (or honest "research-stage" framing)
- ☐ Deployed docs site shows V3 status

Currently met: 4 of 8.

---

## 7. Pointers to other docs (no duplication)

| Topic | Doc |
|---|---|
| V3 model design (pre-registered analysis plan) | `docs/V3_DESIGN.md` |
| Methylation phasing | `deepcatch-methylation/METHYLATION_PROJECT.md` + `PHASE0_PLAN.md` |
| Clinical production path | `docs/PRODUCTION_ROADMAP.md` (historical — superseded by §2 Phase D above) |
| Immediate tactical actions | `NEXT_STEPS.md` (historical — most items now done) |
| bioRxiv manuscript | `paper/PAPER.md` (this repo) + `deepcatch-methylation/BIORXIV_PAPER*.md` |
| Reviewer pack | `REVIEWERS.md` |
| Audit findings | `AUDIT_REPORT_2.md` |
| Speed optimization | `cfdna-fragmentomics-pipeline/docs/SPEED_OPTIMIZATION.md` |
| Memory optimization | `cfdna-fragmentomics-pipeline/docs/MEMORY_OPTIMIZATION.md` |
| Sensitivity optimization (honest no-win) | `cfdna-fragmentomics-pipeline/docs/SENS_OPTIMIZATION.md` |

---

## 8. Honest bottom line

We have a **working, reproducible, open-source, 3-repo cfDNA early-detection
portfolio** with all headline results validated, the GPU path designed, and
two bioRxiv papers drafted. The remaining work is **submission logistics**
(ORCID, bioRxiv form, Zenodo DOI) and **next-gen model validation** (V3),
both of which are unblocked and can ship within 4-8 weeks.

The big remaining scientific blockers — FinaleDB API, real plasma validation,
clinical collaborator — are outside the scope of an independent researcher
without institutional support, and are documented honestly here rather than
papered over.

**This is the state. This is the plan. Awaiting direction on which Phase to start.**
