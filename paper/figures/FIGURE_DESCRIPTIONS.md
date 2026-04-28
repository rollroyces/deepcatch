================================================================================
DeepCatch Paper: Figure Descriptions
================================================================================

All figures are placeholders. Actual figures should be generated from the
agent results using matplotlib/seaborn and placed in this directory.

FILE NAMING: figures/fig{N}_{description}.png

================================================================================

Figure 1: DeepCatch System Architecture Overview
------------------------------------------------
figures/fig1_architecture.png

DESCRIPTION: A schematic overview of the DeepCatch pipeline with four panels:
  Panel A: Variant Calling — Bayesian + contrastive learning pipeline with
    input (sequencing reads) → position-specific error modeling → contrastive
    embedding → variant probability score
  Panel B: Multi-Modal Fusion — 6 modality encoders (ctDNA, methylation,
    fragmentomics, copy number, CTC, miRNA) feeding into the heterogeneous
    GNN with 9 edge types → global attention pooling → cancer score
  Panel C: Longitudinal Tracking — Timeline showing quarterly blood draws →
    CET score accumulation → detection threshold crossing at ~306 days
  Panel D: Ensemble & Risk Stratification — All detection signals → stacked
    ensemble with isotonic calibration → 4-tier risk output → clinical action

STYLE: Clean schematic with color-coded components. Blue = single-modality
encoders, green = fusion/integration, red = output/decision layer.

================================================================================

Figure 2: Ultra-Low VAF Variant Detection Performance
------------------------------------------------------
figures/fig2_variant_calling.png

DESCRIPTION: Three-panel figure showing variant calling results:
  Panel A: ROC curves for all 7 callers (SimpleVAF, MuTect2-like, VarDict-like,
    UMI Simple, Fisher Exact, Bayesian, Contrastive Ensemble) with AUC annotations
  Panel B: Bar chart of per-VAF sensitivity for the contrastive ensemble.
    X-axis: VAF levels (5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 0.01). Y-axis: sensitivity.
    Bars show fraction of 30 variants detected at each VAF
  Panel C: Precision-recall curves (more informative than ROC given the extreme
    class imbalance with 180 positives vs 99,820 negatives)

DATA SOURCE: agent1-variant-calling/evaluation_report.md
             agent1-variant-calling/caller_comparison.csv
             agent1-variant-calling/per_vaf_sensitivity.csv

================================================================================

Figure 3: Multi-Modal Fusion Performance
------------------------------------------
figures/fig3_multimodal_fusion.png

DESCRIPTION: Three-panel figure:
  Panel A: Bar chart comparing AUC for all models. Two groups: single-modality
    baselines (6 bars) and fusion methods (3 bars). GNN heterogeneous bar
    highlighted in bold color. Values: CTC 0.619, Variants 0.609,
    Fragmentomics 0.542, CN 0.524, miRNA 0.474, Methylation 0.464,
    Late Fusion 0.614, Cross-Attention 0.595, GNN Hetero 0.692
  Panel B: Overlaid ROC curves for single modalities (thin, light) and GNN
    (thick, dark). AUC annotations. Inset shows zoom of 0.95–1.00 specificity
  Panel C: Schematic of the heterogeneous graph with all 6 node types and
    9 edge types labeled. Node sizes proportional to per-node attention weights
    learned by the global pooling layer

DATA SOURCE: agent2-multimodal-fusion/results/results.json

================================================================================

Figure 4: Longitudinal Detection via Cumulative Evidence Tracking
------------------------------------------------------------------
figures/fig4_longitudinal.png

DESCRIPTION: Four-panel figure:
  Panel A: Representative VAF trajectories for 3 patient types (healthy,
    cancer, benign) over 730 days. Points = quarterly measurements with
    error bars showing Poisson 95% CI. Cancer trajectory shows exponential
    growth; benign shows transient spike at day 270
  Panel B: Corresponding CET score trajectories for the same 3 patients.
    Horizontal dashed line = detection threshold (τ=7.6). Red star marks
    detection at day 306 for cancer patient
  Panel C: ROC curves comparing CET, BOCD-v2, Kalman-Adaptive, Transformer,
    and single-timepoint baseline. CET dominates (top-left corner)
  Panel D: Bar chart of sensitivity at fixed 99.95% specificity for all
    methods. CET = 100%, BOCD-v2 = 65.2%, Kalman = 98.5%, Transformer = 100%,
    Single-TP = 64.3%

DATA SOURCE: agent3-longitudinal/results/final_results.json

================================================================================

Figure 5: Ensemble Strategy and Risk Stratification
-----------------------------------------------------
figures/fig5_ensemble.png

DESCRIPTION: Four-panel figure:
  Panel A: Heatmap showing ensemble sensitivity gain as a function of
    inter-detector correlation (ρ, x-axis) and number of detectors
    (y-axis). Color = sensitivity gain over best individual detector.
    Gain is positive (blue) for ρ < 0.5 and negative (red) for ρ ≥ 0.7
  Panel B: Learned ensemble weights as pie chart or horizontal bar chart:
    Variant Calling 39.8%, Fragmentomics 28.1%, Methylation 22.8%,
    Multi-Modal Fusion 8.1%, Longitudinal 1.2%
  Panel C: Risk score distribution histogram (log scale) with 4-tier
    boundaries marked. Red ticks above = cancer cases. Tier annotations
    show population %, action, and cancer capture rate
  Panel D: MAML few-shot adaptation curve. X-axis = number of examples
    (1, 3, 5, 10). Y-axis = balanced accuracy. Three lines: MAML (blue,
    >99% from 1-shot), transfer learning (green, intermediate), random
    init (red, poor with <50 examples)

DATA SOURCE: agent6-ensemble/simulation_results.json
             agent6-ensemble/SUMMARY.md

================================================================================

Figure 6: Benchmark Comparison and Clinical Impact
---------------------------------------------------
figures/fig6_benchmark.png

DESCRIPTION: Two panels:
  Panel A: Comparison table/chart of DeepCatch vs published methods
    (CancerSEEK, GRAIL/Galleri, DELFI, PanSeer) with annotations for
    sensitivity, specificity, number of cancer types, and detection stage.
    Note caveat that DeepCatch results are from synthetic data
  Panel B: Timeline diagram showing the estimated diagnostic window shift.
    Tumor volume on log scale (x-axis), progression timeline with key
    milestones: tumor initiation → DeepCatch detection (~3 mm³) → current
    ctDNA detection (~100 mm³) → imaging detection (~10,000 mm³) → clinical
    symptoms → late-stage diagnosis. Annotated with lead-time estimates
    (6–18 months)

================================================================================
