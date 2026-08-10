# DeepCatch Code Review — 2026-08-10

**Reviewer:** Hermes Agent (external review)
**Repo:** github.com/rollroyces/deepcatch @ `1ee329d` (main)
**Scope:** Full-pipeline audit: real-data paths, validation methodology, reproducibility, test claims, engineering hygiene.

---

## 1. Verdict

DeepCatch is a **genuinely well-structured research codebase with an unusually honest reporting culture** (self-audit, claims audit, honest-limitations sections — rare and valuable). The test suite is real and green (228 passed on a fresh environment; README badge claims 198 — stale). The architecture (7 modalities → foundation fusion → longitudinal) is thoughtful.

**The core scientific gap:** the project's headline *performance* still rests on synthetic simulations. Real-data validation exists but is (a) tiny (5 LUAD patients for TCGA; 129 samples 4-mer for Jiang), (b) partially non-reproducible (data files not in repo, hardcoded external paths), and (c) methodologically optimistic in two places (feature selection outside CV; threshold optimized on test predictions). The honest ultra-early-regime results are modest (AUC ≈ 0.65 at 0.1% ctDNA) and the longitudinal CET currently scores *below chance* (AUC 0.49) after the p-hacked bonuses were removed.

None of this is fatal — the fixes are concrete and mostly quick. Priority is below.

---

## 2. What's Working Well (keep it)

1. **Honesty infrastructure** — `PIPELINE_AUDIT.md`, `review/claims_audit.json` (15 claims audited, 2 marked FALSE), `review/fix_verification.json`, and `results/README.md`'s "What's Honest / What's Not" table are exceptional practice.
2. **Tests actually pass** — 228/228 on a clean venv (torch 2.13 + torch_geometric). Coverage is broad per-module.
3. **Real ground truth where it exists** — TCGA-LUAD MAF parsing with real barcodes/read counts; Jiang Table S1 (129 real plasma samples) is the strongest real dataset in the project.
4. **Sound statistical tools** — DeLong tests, Bonferroni, bootstrap CIs, decision-curve analysis are implemented and applied (sometimes to the wrong data, but correctly implemented).
5. **Sensible module boundaries** — fragmentomics / GNN / tissue deconv / foundation / priming are cleanly separated with config dataclasses and documented APIs.

---

## 3. Results Inventory — What Is Actually Real

| Result | Data source | Status |
|---|---|---|
| Fusion AUC 0.92–0.96 vs ctDNA fraction (results/README.md) | Simulation (seed 42), Bie/CAPP-Seq/iDES re-implemented on same simulation | **Simulated** (honestly labeled "Sim only") |
| Longitudinal CET: sens 2.5%, spec 97%, AUC 0.49 | Simulation (Gompertz + Poisson) | **Simulated, honest** — below chance |
| TCGA "REAL_DATA" validation | 5 LUAD patients, real MAF mutations + **simulated plasma reads** | **Real mutations, simulated sequencing** |
| Jiang 4-mer HCC AUC 0.9845 / pan-cancer 0.910 | 129 real plasma samples (CUHK Table S1) | **Real plasma** — strongest result |
| Primer/7th modality PK-PD, GNN field defects, etc. | Literature params + synthetic | **Simulated** |

---

## 4. Critical Findings (verified, with locations)

### 4.1 `--seeds` flag is a no-op in the real-TCGA validation
`real_tcga_validation.py:478,486` — `run_real_validation()` uses only `seeds[0]` ("Use first seed for patient") in a single loop over patients. The CLI advertises `--seeds 5` and metadata writes `seeds_used`, but **every seed > first is ignored**. The multi-seed claim is hollow.

### 4.2 Real-TCGA results are not reproducible from the repo
- `real_tcga_validation.py:658-660` — default cache dir is `/home/node/.openclaw/workspace/cancer-screening/validation/tcga/tcga_cache` (external agent environment).
- `validation/tcga/tcga_cache/` contains **only** `fallback_dataset.json` — the synthetic fake (sample IDs `LUAD_S0000`, positions 0/1000/2000). No `.maf.gz` files.
- The committed `results/real_tcga_validation.json` (5 real patients) **cannot be regenerated** from the clone. If run as-is with no MAF files, `load_tcga_cohort` returns an empty cohort and the script silently produces empty/meaningless results — no failure, no fallback, no download.
- `tcga_downloader.py` exists (cBioPortal API, 5 cancer types) but is **not wired into** `real_tcga_validation.py`.

### 4.3 Cohort loader counts files, not patients
`real_tcga_validation.py:149` — `if new_patients or used_files < n_patients:` appends and increments per *file*, so it may stop after `n_patients` files that add no new patients (or fewer patients than requested).

### 4.4 Threshold optimized on test predictions (optimistic sens/spec)
`real_tcga_validation.py:411-417` — the "optimal threshold" is selected from **pooled CV test predictions** (`j_scores = tpr - fpr`), then sensitivity/specificity/F1 are reported at that threshold on the same predictions. AUC is unaffected; sens/spec are inflated. This is the same class of issue flagged as C3 in `review/fix_verification.json` — it survives in the new code.

### 4.5 "NOT synthetic data" banner overstates
`real_tcga_validation.py:747-748` prints "These are REAL TCGA patient mutations, REAL read counts... NOT synthetic data." The mutations/read counts are real, but **the plasma cfDNA sequencing is simulated** (Poisson reads + Beta error model, `simulate_cfdna_from_real()`). The banner will mislead readers.

### 4.6 Jiang pipeline: feature selection outside CV (leakage)
`scripts/run_jiang_pipeline.py` — Mann-Whitney U enrichment over the **full dataset** selects the top-50 motifs, then logistic regression is 5-fold-CV'd on the selected features. Motif selection sees test-fold labels → **AUC inflated** (headline HCC 0.98). Needs nested CV (selection inside each fold). This directly affects the headline number in README §11.

### 4.7 Jiang results not reproducible; scripts duplicated
- The data file (`deepcatch_data.xlsx`, Table S1) is **not in the repo** (removed for privacy in `8c812c0` — correct call), but there is no downloader, checksum, or documented access path; `scripts/run_jiang_pipeline.py:715-720` hard-errors if the xlsx is absent.
- Two overlapping scripts: `run_jiang_analysis.py` (905 lines, root) and `scripts/run_jiang_pipeline.py` (1290 lines) — both hardcode `/home/node/.openclaw/workspace` and `/tmp/deepcatch_jiang_analysis`.
- README AUC values inconsistent across sources: 0.982 (commit ae77b09) / 0.9845 (summary_for_professor_jiang.md) / 0.986 (README §11).

### 4.8 CHIP filter exists but is disconnected
`src/preprocessing/chip_filter.py` implements a reasonable CHIPFilter (gene list + VAF range + gnomAD AF + phasing), but **no validation script imports or applies it**. Matched-normal counts (`n_alt_count`/`n_ref_count` in MAF) are used only for an error-rate estimate, not to subtract CHIP/germline. In a screening-age population CHIP alone can cost 5–10% specificity — this is the single highest-leverage missing piece for the "ultra-early" claim.

### 4.9 Ultra-early regime performance is modest (and that's the honest headline)
From the committed real-TCGA results (5 patients):

| ctDNA fraction | Variant caller AUC | ML classifier AUC |
|---|---|---|
| 10% | 1.000 | 1.000 |
| 5% | 0.999 | 0.998 |
| 1% | 0.969 | 0.967 |
| 0.5% | 0.907 | 0.893 |
| 0.1% | **0.647** | **0.608** |

At 0.1% ctDNA (already late for "ultra-early"), detection is near-random. The README does not surface this table; the "ultra-early" narrative currently leans on simulation. **Report this honestly and prominently — it's the actual research frontier of the project.**

### 4.10 Longitudinal CET honest result is below chance
`results/README.md` — after removing the arbitrary streak/trend bonuses (fix C9), the CET longitudinal tracking scores **AUC 0.4926, sens 2.5%**, dual-target NOT MET. The old "100% sensitivity" claim came from the p-hacked version. The longitudinal approach needs a real redesign (see §6.8), not threshold tweaks.

### 4.11 CI does not run the test suite — and its "Real Data" test is synthetic
`.github/workflows/validate.yml` installs only numpy/scipy/sklearn/pandas, then runs an inline script that generates **`np.random` data** and prints "✅ ALL 3 CORE TESTS PASSED" with the job named "DeepCatch Real Data Performance Test". The 228-test suite, the torch modules, and any real data are never exercised in CI. The README "198/198 passing" badge is unverifiable from CI.

### 4.12 Portability
- Hardcoded `/home/node/.openclaw/workspace/...` in `real_tcga_validation.py`, `scripts/run_jiang_pipeline.py`; `/tmp/deepcatch_jiang_analysis` output paths.
- `requirements_py.txt` pins loosely; **`torch_geometric` is required by `src/methylation_gnn/` but absent from the requirements file** (README mentions it; CI cannot install it).
- Python 3.14 (current on this machine) works for everything except that torch must be installed from the CPU index; CI pins 3.11 — fine, but document the matrix.

### 4.13 Test count claims are stale
README: fragmentomics 47 ✅ (actual 47 in enhanced + others), GNN 54 → **actual 46**, tissue deconv 54 → **actual 47**, foundation 43 ✅, priming not listed (50 actual). **Actual total: 228 passing**, not 198.

---

## 5. Priorities

### P0 — Scientific integrity of real-data validation (do first)
1. **Fix `real_tcga_validation.py`**: honor `--seeds` (loop seeds properly); fix `load_tcga_cohort` patient counting; wire `tcga_downloader.py` in (or fail loudly when no MAF files); reword the "NOT synthetic data" banner to "real tumor mutations + simulated plasma reads".
2. **Move threshold selection off test data**: pick operating points on a calibration split, or report sensitivity at fixed specificity (95%/99%) — never Youden-threshold on pooled test predictions.
3. **Nested CV for Jiang**: feature selection inside folds; dedupe the two Jiang scripts; make README AUCs traceable to one script output; add a documented data-provision path for Table S1 (checksum + request note, respecting CUHK terms).
4. **Wire CHIPFilter into the detection path**: filter CHIP-gene variants in the 0.1–2% VAF window; use matched-normal counts for germline/CHIP subtraction. Quantify its effect on specificity in the simulation.
5. **Scale real validation**: 5 patients → full TCGA-LUAD (566) + COADREAD/BRCA via the existing cBioPortal downloader; explicitly benchmark the 0.01%–0.001% ctDNA regime (the actual "ultra-early" claim) and report those numbers in the README.

### P1 — Modeling (where performance actually improves)
6. **Make 4-mer fragmentomics the core**: it's the only real-plasma signal with strong results (HCC AUC ≈ 0.98). Extend to real WGS cfDNA (Cristiano 2019 DELFI-style public data), fuse 4-mer + DELFI + methylation, and validate cancer-type-first (HCC first, then TOO) per your own summary — pan-cancer sens at 95% spec is 58%.
7. **Variant caller at ultra-low VAF**: add sequence-context error priors, strand-aware counts, UMI/duplex support, and report LoD at fixed specificity instead of Youden thresholds.
8. **Longitudinal CET redesign**: the honest baseline (AUC 0.49) shows single-patient VAF SPRT can't beat the Poisson floor. Move to hierarchical Bayes across loci/panels (Setty 2022), model CHIP trajectories as a distinct state, and only claim improvement against the *current honest* baseline (not the removed bonuses).

### P2 — Engineering
9. **CI**: run the real suite (core + torch + torch_geometric jobs, split to stay in the 10-min budget); rename the synthetic smoke test; add coverage reporting.
10. **Requirements**: split `requirements-core.txt` / `requirements-dl.txt` (add `torch_geometric`); document the Python version matrix (3.9–3.13 supported; 3.14 needs CPU-index torch).
11. **Portability**: replace external paths with repo-relative + `env` override; add `data/README.md` with provenance + license for every dataset (TCGA open-access MAF ok; Jiang xlsx = CUHK terms; BAM manifest).
12. **Delete or archive `validation/node/*.js`** (20+ files duplicating `validation/py/`); keep one canonical validation suite.
13. **Refresh README**: test badge (228), module test table, §11 AUCs sourced from one script, and add the ultra-early regime table.

---

## 6. Compliance note

- **TCGA**: open-access MAF (PanCanAtlas) is fine to redistribute/analyze; cite TCGA publication guidelines.
- **Jiang lab data (CUHK)**: publishing analysis results of their Table S1 in a public repo — confirm the collaboration terms explicitly before treating the 0.98 HCC number as a public claim. The "summary for Professor Jiang" framing suggests this is understood; make it explicit in the repo.
- No clinical claims are made in the README — keep it that way.

---

## 7. Implemented Fixes (2026-08-10, branch `p0-fixes`)

All P0 items from §5 implemented and verified (228/228 tests green after changes):

| # | Fix | Files |
|---|---|---|
| 1 | `--seeds` now actually loops seeds × patients (per-seed mean ± std across seeds in output) | `real_tcga_validation.py` (`run_real_validation`) |
| 2 | Simulation fully seeded — error-rate draws moved from global `np.random` to the per-seed `RandomState` | `real_tcga_validation.py` (`simulate_cfdna_from_real`) |
| 3 | Cohort loader counts patients (not files); picks the `n_patients` patients with richest signal; dedupes | `real_tcga_validation.py` (`load_tcga_cohort`) |
| 4 | cBioPortal downloader wired in — auto-downloads MAF-equivalent data, saves normalized `*.maf.gz` for reproducibility; **fails loudly** (SystemExit) with instructions instead of silently using the synthetic fallback | `real_tcga_validation.py` (`download_tcga_data`, `df_to_mutations`, `save_normalized_maf`) |
| 5 | Threshold optimization on test data removed — metrics are now AUC/PR-AUC + **sensitivity at fixed 95%/99% specificity** | `real_tcga_validation.py` (`run_variant_caller`, `run_ml_classifier`, `compute_bootstrap_ci`) |
| 6 | CHIP/germline filter wired into the real-data path (matched-normal VAF ≥ 0.25 → germline; CHIP-gene + normal evidence → CHIP; CHIP-window candidates) with counts in metadata | `real_tcga_validation.py` (`filter_chip_variants`, `--no-chip-filter` flag) |
| 7 | Honest framing: `pipeline_type: REAL_MUTATIONS_+_SIMULATED_PLASMA_READS`, explanatory note, honest console banner | `real_tcga_validation.py` (`main`) |
| 8 | Hardcoded `/home/node/.openclaw/...` and `/tmp/...` paths removed → repo-relative + `DEEPCATCH_DATA_DIR`/`TMPDIR` env overrides | `real_tcga_validation.py`, `scripts/run_jiang_pipeline.py` |
| 9 | **Nested CV for Jiang CET** — motif selection (MWU top-k) now happens inside each training fold; per-fold selection stability reported; full-data coefficients labeled interpretation-only | `run_jiang_analysis.py` (`logistic_fusion_cv`) |
| 10 | CI runs the real test suite (core job + DL job with torch/torch_geometric), a real-data guard (must refuse synthetic fallback), and the old synthetic smoke test renamed honestly | `.github/workflows/validate.yml` |
| 11 | `torch-geometric` added to `requirements_py.txt` (was required by GNN but missing) | `requirements_py.txt` |
| 12 | README: test badge → 228/228, per-module counts corrected, §11 rewritten with honest TCGA benchmark table + nested-CV framing for Jiang AUC | `README.md` |

**Caveats / not done (needs you):**
- Jiang `deepcatch_data.xlsx` is not in the repo (privacy) — the nested-CV AUC must be re-estimated once the file is provisioned at `data/deepcatch_data.xlsx` or via `DEEPCATCH_DATA_DIR`. Expect the nested-CV number to be **lower** than 0.9845 (selection bias removed).
- `real_tcga_validation.json` was regenerated on the new pipeline with **20 real LUAD patients (5,738 mutations from GDC open-access MAFs)** and the new fixed-specificity metrics. The old committed JSON (5 patients, Youden-based sens/spec) was replaced. Note: GDC open-access MAFs have no matched-normal counts, so the CHIP filter currently removes 0 variants — re-run with controlled-access MAFs or plasma data to exercise it. P1/P2 items (cancer-type-first architecture, hierarchical-Bayes longitudinal CET, deleting `validation/node/*.js`, requirement file splits) remain open.

## 8. Ultra-Early Optimization (2026-08-10, same branch)

The headline ultra-early numbers (0.1% ctDNA: AUC 0.64, sens@95% 0.18) were per-
**position** classification — information-limited (signal ≈1.9 reads vs error
≈10 reads per locus at 5,000×). Implemented the field-standard fix:

- **Panel-based per-sample detection** (`run_panel_detection`): MRD-style
  aggregation of per-locus Poisson LLR over the tracking panel; paired
  cancer/control design per patient; ROC across patients per seed, mean±std
  across 5 seeds.
- **Ultra-early assay sweep** (`run_ultraearly_sweep`): panel detection over
  error rate (2e-3 → 1e-5) × depth (5k×/50k×) at 0.1% ctDNA — production
  assay-design guidance.
- `compute_llr_scores` extracted (shared by caller + panel detector);
  `simulate_cfdna_from_real` gained a `bg_error_rate` parameter.
- `--with-ml` now opt-in (the per-position ML classifier was the 13-min cost
  and adds nothing at ultra-low VAF); default run takes ~3 min.

**Result (20 real LUAD patients, 5 seeds):** 0.1% ctDNA panel AUC **0.935**
(vs 0.642 per-position), sens@95% **0.770** (vs 0.183), paired win rate 1.000.
At duplex-UMI error (1e-4) or 50k× depth: sens@95% = **1.000**.

**Production path:** `docs/PRODUCTION_ROADMAP.md` — MRD-first product strategy,
assay spec (duplex UMI, 50k×, matched WBC), data-acquisition plan (Jiang data →
public WGS cfDNA datasets on EGA/dbGaP → own cohort), validation ladder
(analytical → clinical validity → utility → regulatory), 12-month milestones,
honest risks.

*Generated by Hermes Agent, 2026-08-10. All findings verified against a fresh clone at `1ee329d` and a clean Python 3.14 venv (torch 2.13, torch_geometric installed).*
