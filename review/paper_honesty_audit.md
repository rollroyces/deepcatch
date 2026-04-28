# Paper Honesty Audit — Claim-by-Claim Verification

**Audit Date:** 2026-04-28
**Source:** `paper/main_fixed.tex` vs `paper/CLAIMS_LOG.md` vs Fixed Code
**Method:** For each OLD→NEW claim change in CLAIMS_LOG.md, verify presence and accuracy in main_fixed.tex

---

## Overall Assessment

The paper has been substantially corrected from the original. Every inflated claim has been either removed or properly qualified. The framing shifted from "established detection capability" to "proof-of-concept simulation study." The Limitations section is now the FIRST subsection of Discussion, organized by severity — a strong signal of scientific honesty.

**Overall Verdict: ✅ HONEST — Paper properly represents the evidence with one paper-code mismatch found.**

---

## Claim-by-Claim Verification

| # | CLAIMS_LOG OLD→NEW | In Paper? | Location in main_fixed.tex | Correct? | Notes |
|---|-------------------|-----------|---------------------------|----------|-------|
| 1 | Title: "Ultra-Early Pan-Cancer Detection at 0.001% VAF" → "A Proof-of-Concept Simulation Study" | ✅ | \title{} | ✅ | Title line directly matches CLAIMS_LOG |
| 2 | Abstract: no simulation flag → "This is a proof-of-concept simulation study; no real patient data were used" | ✅ | Abstract, first paragraph (bold) | ✅ | Bold emphasis in source |
| 3 | Abstract: "17% sensitivity" → "17% (5/30; 95% CI: 6-35%) of simulated variants" | ✅ | Abstract, Results section | ✅ | CI added, "simulated" qualifier |
| 4 | Abstract: "AUC 0.692, 11.9% improvement" → "AUC 0.692 (95% CI: 0.58-0.80; n=90); statistical significance not established" | ✅ | Abstract, Results section | ✅ | CI, n, and significance caveat all present |
| 5 | Abstract: "100% sensitivity at 99.95% specificity" → "1.000 sensitivity (95% CI: 0.996-1.000) at 0.9995 specificity in a single simulation run" | ✅ | Abstract, Results section | ✅ | "1.000" not "100%", CI, "single run" |
| 6 | Abstract: "1.55-fold improvement" → "representing improved detection... achieved 64.3% sensitivity (95% CI: 63.2-65.4%)" | ✅ | Abstract, Results section | ✅ | "1.55×" removed; side-by-side reporting |
| 7 | Abstract: "enabling tumor detection at approximately 3mm³" → "under model assumptions to estimated tumor volumes of ~3mm³" | ✅ | Abstract, Results section | ✅ | "under model assumptions" qualifier |
| 8 | Abstract: "reduces screening costs ~40%" → "may reduce per-cancer-found costs by ~40% under simplified cost model assumptions" | ✅ | Abstract, Results section | ✅ | "may", "simplified", "assumptions" |
| 9 | Abstract: "6-18 month lead-time" → "model-based projections suggest advance... actual lead times depend on biology" | ✅ | Abstract, Results section | ✅ | Specific numbers removed |
| 10 | Abstract closing: no validation requirement → "prospective clinical validation with real patient cohorts is required" | ✅ | Abstract, Conclusions | ✅ | Ends with validation requirement |
| 11 | Results: contrastive data leakage → added caveats about small n, single seed | ✅ | Results §3.1, Table 1 footnote | ✅ | "further validation with independent test sets" |
| 12 | Table 1: "17% sensitivity" → added "95% CI: 6-35%", footnote about single seed | ✅ | Table 1 | ✅ | CI column added; footnote present |
| 13 | Results: GNN "11.9% improvement" → "statistical significance is not established", CI on ΔAUC | ✅ | Results §3.2 | ✅ | "CI spans [-0.03, +0.18]" |
| 14 | Results: CET "100% sensitivity" → "1.000 sensitivity [0.996-1.000], single simulation run" | ✅ | Results §3.3 | ✅ | Both CI and single-run qualifier |
| 15 | Table 3: Temporal Transformer "100%/100%" → "trajectory classifier, not clinical cancer detector" | ✅ | Table 2 footnote | ✅ | Explicit caveat about 0% STABLE |
| 16 | MAML section entirely: "99% 1-shot accuracy" → REMOVED | ✅ | Results §3.4.3 | ✅ | "Preliminary exploration only; same data for both — uninterpretable" |
| 17 | Benchmark: "5.3× Grail" → REMOVED; separated "Clinical" from "Simulation" with "not directly comparable" | ✅ | Table 4 | ✅ | Two sections clearly separated |
| 18 | Results: "17-34× earlier detection" → REMOVED | ✅ | Not found in paper | ✅ | Completely removed |
| 19 | Results: "projected 92%/98.5%" → "under stated degradation assumptions... may differ substantially" | ✅ | Discussion, Limitations | ✅ | Couched as assumptions |
| 20 | Discussion: limitations buried → Limitations FIRST, severity-graded | ✅ | Discussion §4.1 | ✅ | 5 CRITICAL, 6 MAJOR, 4 MINOR |
| 21 | Methods: no validation standards → NEW "Validation Standards" section | ✅ | Methods §4.1 | ✅ | Documents splits, seeds, CIs, missing tests |
| 22 | Supplementary: sensitivity analysis → NEW sensitivity analysis results | ✅ | Appendix, Table S7 | ✅ | 5 break scenarios documented |
| 23 | CET bonus weights: arbitrariness noted | ⚠️ SEE BELOW | Results §3.3 and Methods §4.4.2 | ⚠️ | Paper still describes bonuses; code removed them |
| 24 | Single-timepoint: pooled measurements warning added | ✅ | Methods §4.1 | ✅ | "may inflate effective sample size" |

---

## ⚠️ ISSUE: Claim 23 — CET Algorithm Mismatch

### Paper says (Results §3.3):
> "CET computes a running evidence score incorporating a sequential probability ratio, a streak bonus for consecutive measurements above baseline, and a trend bonus for positive log-linear slopes."

### Paper says (Methods §4.4.2):
> S_t = Σ [log P(m|λ_grow)/P(m|λ_stable) + β_s·𝕀(streak≥3) + β_t·β̂_log-linear]
> where β_s = 1.5 and β_t = 2.0

### Paper says (Limitations, minor #13):
> "CET bonus weights are arbitrary. The streak bonus weight and trend bonus weight are fixed values without formal calibration."

### Fixed code says (improved_methods_fixed.py, C9 fix):
```python
# FIXED C9: Pure Sequential Probability Ratio Test (SPRT)
# REMOVED arbitrary bonuses:
#   - streak_bonus = 0.5 * min(n, 5)  (REMOVED)
#   - trend_bonus = max(0, slope) * 3.0  (REMOVED)
# Now uses ONLY the log-likelihood ratio
```

### Assessment:

The fixed code REMOVED the bonuses entirely as part of fix C9 (arbitrary bonus weights). CET now uses pure SPRT — mathematically cleaner and more defensible. However, the **paper was never updated to reflect this change**. The paper still describes the OLD algorithm with bonuses.

**This is a genuine paper-code inconsistency.** Two options:
1. Update the paper to describe CET as "pure SPRT" and remove all bonus mentions — matches the fixed code, simpler, more scientifically defensible
2. Update the code to keep the bonuses but document them honestly — re-introduces C9 concern about arbitrary hyperparameters

**Recommendation:** Option 1 (match paper to code). Pure SPRT is scientifically cleaner. The paper should state: "CET uses a pure Sequential Probability Ratio Test without heuristic bonuses, accumulating Poisson log-likelihood ratios across quarterly measurements."

---

## Additional Paper Audit Checks

### Target Journal Appropriateness

**Paper claims:** PLOS Computational Biology / Bioinformatics

**Assessment:** ✅ APPROPRIATE. This is a methods-focused, simulation-only study. PLOS Computational Biology accepts well-executed simulation studies with proper validation. The honest framing as "proof-of-concept" aligns with the journal's scope. The original target of Nature Medicine/Lancet was wildly inappropriate for simulation-only work.

### Limitation Prominence

**Paper:** Limitations are the FIRST subsection of Discussion (§4.1), before Clinical Implications.

**Assessment:** ✅ EXCELLENT. Most papers bury limitations at the end. Putting them first is a strong signal of scientific integrity. The severity-grading (CRITICAL/MAJOR/MINOR) is clear and actionable. Every limitation includes an explicit statement about what it means for result interpretation.

### CI Coverage

**Assessment:** ✅ All headline claims now include 95% CIs. The supplementary Table S1 (CI Summary) provides a single reference point. The paper explicitly notes that some CIs are based on small samples and should be interpreted cautiously.

### Bayesian Results Qualification

**Assessment:** ✅ Fair. The paper communicates that Bayesian-informed results remain simulation-based and that prior choices influence posteriors in small-sample regimes.

### Honesty Patterns

**What the paper does well:**
- Bold "no real patient data" declaration in abstract
- Every quantitative claim has a CI
- Limitations come before implications
- "5.3× Grail", "17-34× earlier", "99% 1-shot" claims fully removed
- Test-data calibration openly disclosed as a current limitation
- Single-seed numbers clearly flagged
- Absence of statistical tests (DeLong, McNemar) honestly noted

**What the paper could improve:**
- Update CET algorithm description to match the fixed code (remove bonus terms)
- Consider estimating ±0.02-0.05 AUC variance from different GNN graph seeds
- Add a statement that "Results should not be cited without noting the simulation-only nature of this work"

---

## Verdict

**Paper Honesty: ✅ HONEST (one paper-code mismatch found)**

The paper has been transformed from an overclaimed manuscript targeting Nature/Lancet into an honest proof-of-concept simulation study targeting PLOS Computational Biology. The single remaining issue (CET algorithm mismatch) is a mechanical synchronization problem — the paper was updated for the old code, then the code removed bonuses, but the paper wasn't updated to match. This is a 30-minute fix.

**The paper is ready for submission after:**
1. Fixing the CET algorithm description (remove bonus terms, describe pure SPRT)
2. Fixing the code integration issues (so code and paper describe the same system)

---

*Paper Honesty Audit — Second Review Agent — April 28, 2026*
