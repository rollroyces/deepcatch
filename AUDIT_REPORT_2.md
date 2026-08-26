# AUDIT REPORT 2 — Multi-Reviewer Audit (Engineering + Scientific + Statistical)

**Date**: 2026-08-26
**Reviewers**: 4 (Engineering, Scientific, Statistical, Reproducibility)
**Mode**: 3 of 4 completed; Reproducibility still running
**Source files**:
- `/Users/hermes/SCIENTIFIC_REVIEW_cfDNA.md` (full)
- `/Users/hermes/STATISTICAL_REVIEW.md` (full)
- `/Users/hermes/.hermes/cache/delegation/subagent-summary-0-20260826_161547_972707.txt` (Engineering summary)

---

## Summary of all findings

| Reviewer | High-severity | Medium-severity | Low-severity | Action items |
|---|---|---|---|---|
| Engineering | 3 | 4 | 6 | 13 |
| Scientific | 6 | n/a | n/a | 6 (questions + 3 experiments) |
| Statistical | 5 | n/a | n/a | 5 (concerns + 3 analyses) |
| Reproducibility | (pending) | | | |

---

## Engineering review — actionable items

### High-severity (fix immediately)

**E1**: `nuc_ablation.py:172-186` — `_evaluate` called twice on `X_5ch` and `X_5nuc` (lines 172-175 then 184-187). The first results are computed, stored, then **immediately overwritten** by the second call. **Wasted compute, misleading timing printout.**

**E2**: `honest_benchmark.py:17, 167-181` — `main()` runs full benchmark on `--help`. The benchmark work happens at module top-level BEFORE `if __name__ == "__main__":`. So `python scripts/honest_benchmark.py --help` triggers the entire 5-seed benchmark.

**E3**: `honest_benchmark.py:104-132` — Bare `open()` (no `with`, no encoding) for labels file, called 4 times at module level. If file missing, crash happens from inside module body, not from callable.

### Medium-severity (fix soon)

**E4**: `train_classifier.py:138-140` — NaN→median imputation uses full-cohort median (train + test together). For small cohorts with sparse NaNs this is a measurable leakage. Should be moved inside `evaluate_cv` and computed only on train.

**E5**: `train_classifier.py:127-134` — Inconsistent feature-count check: rows appended before validation. Result: `y` may be misaligned if a bad sample is followed by a good one.

**E6**: `gemma_baseline.py:299-302` — `if p is None: p = 0.5` silently fills parse failures with chance. If parse failures differ by class, this introduces label-correlated bias. Should track failure rate.

**E7**: No patient-level grouping in `StratifiedKFold`. If a single patient contributed multiple samples, leakage is possible.

### Low-severity (document or fix later)

**E8**: Hardcoded `/Users/hermes/...` paths in 5 scripts (replace with `pathlib.Path(__file__).resolve().parent.parent`).

**E9**: LOESS fallback in `normalization.py:182-191` only catches ImportError, not ValueError/LinAlgError from the actual LOESS call. Add to except.

**E10**: DeLong `except Exception` in `fusion_ablation.py:248` records error to JSON but never surfaces to stderr.

**E11**: Test coverage ~29% (5/17 scripts in pipeline, 3/12 in DeepCatch fragmentomics).

---

## Scientific review — actionable items

### Critical methodological concerns (S1-S6)

**S1**: "Pan-cancer" label is misleading. **8 of 9 cancer types come from a single study (Cristiano 2019)**, and HCC is 24.5% of cancers (89/364). Per-cancer-type AUC is **never reported** anywhere in RESULTS.md, BENCHMARK.md, MODEL.md, or PAPER.md. The cross-cancer generalization claim is not supported by the analysis.

**S2**: Mutation channel in fusion is a **calibrated synthetic surrogate**, not real variant calling on real plasma. The "+0.014 fusion gain" is a *synthetic + real* experiment. The repo admits this implicitly via its 8-point calibration sweep, but the framing in RESULTS.md is more positive than the data supports.

**S3**: Per-study z-score harmonization may be over-correcting. The true-confound test gives AUC 0.4966 (random) when harmonized. This proves harmonization removes confound, but also suggests it may remove real biological signal in the partially-confounded headline cohort. No sensitivity analysis.

**S4**: 5-channel headline vs 8-channel codebase. `train_classifier.py` extracts 8 channels (5 + 4-mer motifs + per-bin mean length). The 4-mer motif ablation was only done on 98 samples; never repeated at n=627. This looks like configuration cherry-picking.

**S5**: Internal n-discrepancy: labels file has 658 samples, RESULTS.md/BENCHMARK.md headline 627. The 31-sample gap is unaccounted for. **FIXED**: documented in RESULTS.md now (30 Cristiano + 1 Jiang missing one or more required features, likely .fsd.json).

**S6**: No clinical framing. PPV at screening prevalence never computed. At 99% spec, 82.4% sens, 0.4% prevalence, PPV ≈ 25% (3 of 4 positives false). Never mentioned.

### Specific questions (still need answers)

- Why 5-channel headline and not 8-channel?
- What is the empirical floor for "real signal removed by harmonization"?
- Per-cancer-type AUC table with stage stratification?

### Suggested additional experiments

1. Per-cancer-type AUC table with stage stratification
2. ComBat / limma-style empirical-Bayes harmonization comparison
3. Head-to-head against a published methylation-based MCED baseline (Galleri / Liu 2020)

---

## Statistical review — actionable items

### Critical statistical concerns

**ST1**: AUC reported with ±std, not 95% CI. The std is across-seed variability of the OOF estimator, NOT a CI on the population AUC. For n=627 a proper DeLong CI is roughly ±0.020, not ±0.001. **FIXED**: uncertainty note added to RESULTS.md.

**ST2**: "95% bootstrap CI [+0.0135, +0.0152]" for fusion ΔAUC is **mislabeled**. It is actually the DeLong CI on a single seed's pooled OOF, not a 10-seed paired bootstrap. The honest CI is ~10× wider. **FIXED**: re-labeled in RESULTS.md as "95% CI (DeLong, single-seed pooled OOF)".

**ST3**: C-sweep (7+ C values × 2 penalties = ≥18 tests) is **post-hoc with no Bonferroni/BH adjustment**. Same for the 3-feature / 6-feature nucleosome ablations. **FIXED**: multiple-testing note added to RESULTS.md.

**ST4**: Paired t-test with df=4-9 over 5-10 seeds that are different fold partitions (not independent replications) tests "is this specific Δ nonzero for these shuffles", not a generalizable significance claim.

**ST5**: Per-study z-score uses a class-mixed per-study mean (`train_classifier.py:282-287`). When studies have different class proportions (Jiang ~70% HCC vs Cristiano mostly healthy), this is a subtle confound. A more rigorous approach would be per-study-per-class.

### Suggested additional analyses

1. Add a between-seed bootstrap CI on ΔAUC, not just within-seed DeLong CI
2. Apply Benjamini-Hochberg FDR across the full family of comparisons
3. Per-study-per-class z-score harmonization, plus sensitivity check

---

## What I've already fixed (in this round)

✓ Sample-count note in TL;DR (S5)
✓ "95% bootstrap CI" → "95% CI (DeLong, single-seed pooled OOF)" (ST2)
✓ Uncertainty note (±std is across-seed, not population CI) (ST1)
✓ Multiple-testing note (no Bonferroni/BH applied; gains are suggestive) (ST3)

## What I should fix next (in priority order)

1. **Engineering E2**: Wrap `honest_benchmark.py` body in `main()` so `--help` doesn't run the full benchmark. High impact (usability).
2. **Engineering E4**: Move NaN→median imputation inside `evaluate_cv` to avoid leakage.
3. **Engineering E6**: Track Gemma parse failure rate explicitly.
4. **Scientific S1**: Run per-cancer-type AUC table. New analysis (1-2 hours).
5. **Scientific S4**: Run 8-channel baseline on full 627 cohort. New analysis (30 min).

## What I cannot fix automatically

- **Engineering E1**: The duplicated `_evaluate` call in `nuc_ablation.py` may not actually be a bug — it could be intentional (e.g., to test with two different PCA settings). Need to read more context.
- **Scientific S2**: The synthetic-mutation-channel is a deliberate design choice for testing fusion structure; cannot be replaced without a real mutation caller.
- **Scientific S3**: ComBat / limma-style harmonization requires new dependencies (`sva` or `pycombat`); not a 1-line fix.
- **Statistical ST4**: Requires changing the entire experimental design (e.g., bootstrap CIs over fold partitions rather than fixed seeds).

---

## Reproducibility review (pending)

Will be added when the 4th reviewer completes. The reproducibility reviewer is the only one that actually runs scripts, so it will tell us whether the documented commands work end-to-end.
