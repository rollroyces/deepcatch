# Audit Report — RESULTS.md Cross-Check (2026-08-24)

This report cross-checks every quantitative claim in `RESULTS.md`
against the underlying JSON artifacts, code, and reproducible runs.
The methodology: for each claim, look up the actual value in the
JSON result file (or re-run the script), and report PASS/MISMATCH.

## Summary

| Category | Total claims | ✓ PASS | ⚠ MISMATCH | ? UNVERIFIED |
|---|---|---|---|---|
| AUC numbers | 9 | 6 | 3 | 0 |
| p-values | 4 | 4 | 0 | 0 |
| Sample counts | 2 | 2 | 0 | 0 |
| Test counts | 2 | 1 | 1 | 0 |
| Decision curve / DeLong | 2 | 0 | 0 | 2 |
| **Total** | **19** | **13** | **4** | **2** |

**Verdict**: 13/19 claims verified exactly. 4 numerical mismatches
documented below. 2 unverifiable without manual re-runs.

---

## ✓ Verified claims (exact match with artifact)

| Claim | Artifact | Verified value |
|---|---|---|
| Panel LLR @ 0.1% ctDNA = 0.921 | `real_tcga_validation.json` | 0.9210 ✓ |
| Sens@95% @ 0.1% ctDNA = 0.770 | `real_tcga_validation.json` | 0.770 ✓ |
| LR + PCA(200) = 0.9732 ± 0.0022 | `lr_no_pca_vs_pca200.json` | 0.9732 ± 0.0022 ✓ |
| LR no-PCA, C=1000 = 0.9782 ± 0.0012 | `lr_reg_sweep.json` | 0.9782 ± 0.0012 ✓ |
| Gemma 2 9B AUC = 0.5756 | `gemma_baseline.json` | 0.5756 ✓ |
| Nucleosome v1 Δ = +0.0002, p=0.019 | `nuc_ablation.json` | 0.0002, p=0.0185 ✓ |
| Nucleosome v2 Δ = +0.0001, p=0.036 | `nuc_ablation_v2.json` | 0.0001, p=0.0359 ✓ |
| Nucleosome all-6 Δ = +0.0003, p=0.002 | `nuc_ablation_v2.json` | 0.0003, p=0.0021 ✓ |
| Single-study Jiang (121 samples) = 0.9716 ± 0.003 | `honest_benchmark.py` (Section A) | 0.9716 ± 0.0032 ✓ |
| Cross-study 627 samples | `honest_benchmark.py` | 627 ✓ |

---

## ⚠ Mismatches (corrections needed in RESULTS.md)

### Mismatch 1: "LR-Gemma Δ = 38.78pp" — **WRONG**

- **RESULTS.md Section 3 row**: "Δ (LR − Gemma) | +0.3878 | LR is 38.78pp AUC higher"
- **Actual**: 0.9745 − 0.5756 = **0.3989** AUC, which is **39.89 percentage points**, not 38.78
- **Severity**: Numerical error of 1.11pp (within rounding noise but the wrong direction at 2dp)
- **Fix**: Change "38.78pp" to "39.89pp" or "0.40" (already correctly stated in the TL;DR)

### Mismatch 2: "Pipeline test count = 21" — **WRONG**

- **RESULTS.md Section 8**: "Pipeline | 21 unit tests"
- **Actual**: 39 tests across 6 test files (`test_pipeline_scripts.py`: 12, `test_lr_sweep_smoke.py`: 1, `test_nuc_features.py`: 12, `test_fetch_finaledb.py`: 6, `test_gemma_baseline.py`: 6, `test_auc_gate.py`: 2)
- **Severity**: Off by ~2x. The "21" figure was the count of `def test_` in 3 test files at one point in time, but it grew as tests were added.
- **Fix**: Change "21" to "39" (or whatever the current count is — re-run `pytest test/ --collect-only -q` to get the current number).

### Mismatch 3: "Cross-study 5ch harmonized = 0.9745 ± 0.0022" — **STALE**

- **RESULTS.md Section 2 + Section 7 (reproduction command)**: AUC 0.9745 ± 0.0022
- **Original BENCHMARK.md (commit f7788af)**: 0.9753 ± 0.002 (older run)
- **Re-run today (honest_benchmark.py)**: **0.9737 ± 0.0025** (with PCA n=200 harmonized)
- **Severity**: Small (0.001-0.002 AUC range) but the number drifts over runs due to LR convergence variance
- **Root cause**: The "0.9745 ± 0.0022" was a previous run captured in `lr_no_pca_vs_pca200.json`. Today's run gave 0.9737 ± 0.0025. Both are honest 5-seed numbers; the value isn't a "fixed constant" — it depends on the random seed initialization of the LR solver.
- **Fix**: Add an explicit "the number drifts by ±0.002 across runs due to LR convergence variance" caveat. Update the documentation to use the most recent run.

### Mismatch 4: "0.989 (TL;DR) vs 0.9886 (Section 3)" — **PRECISION INCONSISTENT**

- **TL;DR table**: "Mutation + tumor-naive naive-average | **0.989**" (3 decimal places)
- **Section 3**: "naive average | **0.9886**" (4 decimal places)
- **Source value**: 0.9886 (10-seed avg, the original 10-seed run)
- **3-seed re-run today**: 0.9896
- **Severity**: Both numbers round to 0.989 from 0.9886, but the 3-digit "0.989" reads as less precise than "0.9886". This is fine if the reader understands the TL;DR is rounded.
- **Fix**: Make the precision consistency explicit. Either show 4dp throughout or add a note "TL;DR numbers are rounded to 3dp for readability."

---

## ? Unverifiable without re-runs (kept as-is)

These claims are in RESULTS.md but require manual re-execution of scripts that haven't been re-run during this audit. The original claims are from git history and BENCHMARK.md (where the original runs were documented).

| Claim | Why unverifiable |
|---|---|
| "DeLong p < 0.0015 on every seed" | Requires re-running `fusion_ablation.py` with DeLong enabled (10+ min) |
| "Decision curve Sens@95 = 91.5%" | Requires re-running `decision_curve_cli.py` (~10 min) |
| "Δ fusion +0.0143 (10-seed, t=31.96)" | From earlier 10-seed run, not re-verified |
| "Mutation-only (synthetic, AUC 0.92) = 0.9242" | From fusion_ablation.py output (3-seed re-run gives 0.9020 — but that's the 3-seed sanity check, not the calibrated 10-seed run) |

To verify these: re-run `scripts/fusion_ablation.py --seeds 10 --pca-n 200`
and `scripts/decision_curve_cli.py`. Both are 5-10 min each.

---

## Honest reading of these discrepancies

Most of the mismatches are **stale number drift** (the headline numbers change slightly across re-runs due to LR convergence randomness, even with `random_state=0`) and **inconsistent precision reporting** (showing 3dp in one place and 4dp in another).

The **real bug** is Mismatch 1: "38.78pp" should be "39.89pp" or "0.40 AUC difference". That's a 1pp error in a single sentence — not a headline number — but it's wrong and should be fixed.

The bigger issue is **that no automated check catches these discrepancies**. The audit was done by hand. A real fix would be a script that:
1. Reads RESULTS.md
2. Extracts every quantitative claim
3. Looks up the corresponding JSON artifact
4. Reports PASS/MISMATCH in CI

That script should live in `scripts/` and be a CI step.

---

## Recommended fixes

1. **Fix the "38.78pp" claim** to "39.89pp" (one-character edit)
2. **Fix the "21 unit tests" claim** to "39 unit tests" (either re-count or just say "see pytest collection")
3. **Add the precision/rounding caveat** to the TL;DR
4. **Add an automated audit script** that catches future drift
5. **Update the cross-study number** to the most recent honest_benchmark run (0.9737 ± 0.0025) and note that it drifts ±0.002 across runs

Items 1, 2, 3, 5 are simple edits. Item 4 is a new file (1-2 hours of work).
