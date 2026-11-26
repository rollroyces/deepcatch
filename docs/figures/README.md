# DeepCatch validation figures

Four standard cfDNA validation figures, generated from existing JSON
artifacts. Pure plotting (no model re-runs, no recomputation of metrics).
Matplotlib only — no seaborn, no heavy deps.

Two further figures are **schematics** rather than plots: the
multi-modal encoder architecture and the validation protocol. Each is a
self-contained HTML file with an inline SVG, plus the PNG beside it for
the manuscript. They are drawn, not computed — every constant on them
is read out of the code they describe, and the HTML is the source of
truth (the PNG is its export).

## Files

| File | Description | Source JSON |
|---|---|---|
| `fig1_pooled_roc.png` | Pooled cross-study per-seed AUC strip plot + pooled sens@spec operating points (harmonized vs no_harmonize). No fabricated smooth ROC curve — the JSON stores per-seed AUC + sens@spec only, not the underlying `(y_true, y_score)` pair. | `results/cross_study_finallydb.json` (`pooled.harmonized`, `pooled.no_harmonize`, `cohort`) |
| `fig2_per_cancer_roc.png` | Per-cancer per-seed AUC distribution for the top-5 cancers by sample count (LUAD, BRCA, OV, PAAD, HCC_J). Mean ± std annotated; per-seed points jittered. | `results/cross_study_finallydb.json:per_cancer.{LUAD,BRCA,OV,PAAD,HCC_J}` |
| `fig3_calibration.png` | Operating-point curve (sens vs 1−spec) at the three stored spec points (0.95, 0.98, 0.99), with the chance diagonal for reference. **Not** a probability calibration curve — the LR scores are not calibrated probabilities; this is the same sens@spec grid as fig 1b plotted in ROC space. | `results/cross_study_finallydb.json:pooled.{harmonized,no_harmonize}.sens_at_spec_*` |
| `fig4_sens_operating_points.png` | Grouped bar chart of per-cancer sensitivity at spec ∈ {0.95, 0.98, 0.99} for the top-5 cancers. | `results/per_cancer_sens_at_spec.json:per_cancer.{LUAD,BRCA,OV,PAAD,HCC_J}.sens_at_{95,98,99}` |
| `fig5_encoder_architecture.html` / `.png` | **Schematic**: the multi-modal encoder — six modality slots → per-modality `LinearProjection` → `[CLS]` token + 4-layer transformer → the two pre-training phases → the PRODUCTION checkpoint → the 0.7 / 0.3 blend → the loss → the cancer score. Drawn from `src/foundation/config.py`, `src/foundation/model.py` and `scripts/pretrain_production_finaledb.py`. | — (drawn from the code) |
| `fig6_validation_protocol.html` / `.png` | **Schematic**: the validation protocol — paired cohort → per-patient 5-seed × 5-fold GroupKFold → out-of-fold scores → reported metrics, with the three controls (true-confound, shuffled-label, bit-identical sanity) running alongside. Drawn from §2.4 of `paper/METHODS_PAPER.md`. | — (drawn from the paper) |

## Reproducibility

Regenerate every figure with one command:

```bash
env -u PYTHONPATH ./.venv/bin/python \
    ./scripts/plot_figures.py
```

This writes the four PNGs into `docs/figures/`. Each of those four is
≤85 KB (verified — total <300 KB for all 4), well below the 500 KB
disk-pressure cap on this Mac mini. The two schematics are exported
separately (2x device scale, 2432 px wide, 150 KB and 137 KB) and their
HTML is the source of truth; re-export with:

```bash
python3 scripts/export_diagram_png.py \
    docs/figures/fig5_encoder_architecture.html \
    docs/figures/fig5_encoder_architecture.png 2
python3 scripts/export_diagram_png.py \
    docs/figures/fig6_validation_protocol.html \
    docs/figures/fig6_validation_protocol.png 2
```

(The helper renders the HTML and screenshots the `<svg>` element only,
so the editorial header and footer are dropped by design. It needs
Playwright with Chromium.)

Validate-only (no files written, used by `test/test_plot_figures.py`):

```bash
env -u PYTHONPATH ./.venv/bin/python \
    ./scripts/plot_figures.py --validate-only
```

## Honesty contract

These figures do **not** invent metrics. Per the
`cfdna-early-detection-validation` skill (per-seed-mean, never single-
seed AUC, and the standard "no synthetic ROC curve" rule), the pooled
OOF predictions are stored only as summary statistics in the JSON.
Re-running the pooled OOF to recover `(y_true, y_score)` would take
~5 minutes; we deliberately avoid that cost and instead plot the
honest summary:

- **Fig 1**: per-seed AUC distribution + the published sens@spec grid.
- **Fig 2**: per-cancer per-seed AUC distribution (top-5 cancers).
- **Fig 3**: sens vs 1−spec at the three stored spec points (NOT a
  probability calibration; LR scores are uncalibrated).
- **Fig 4**: per-cancer sens@spec operating points, straight from the
  JSON sens_at_spec block.

Every caption makes the data source explicit so the reader can trace
each number back to its JSON file.

## Disk budget

| Figure | Size |
|---|---:|
| fig1_pooled_roc.png | ~68 KB |
| fig2_per_cancer_roc.png | ~72 KB |
| fig3_calibration.png | ~79 KB |
| fig4_sens_operating_points.png | ~54 KB |
| **Total** | **~273 KB** |

Each PNG stays well under the 500 KB cap from the brief.
