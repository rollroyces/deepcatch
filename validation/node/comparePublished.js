#!/usr/bin/env node
/**
 * comparePublished.js — PHASE 5: Honest Comparison vs Real Published Results
 * 
 * Extracts DeepCatch's performance and compares HONESTLY against published clinical assays.
 * 
 * CRITICAL RULES:
 * - If DeepCatch has no clinical validation → state this clearly
 * - If DeepCatch outperforms but at higher sequencing depth → note this
 * - If DeepCatch's TOO is simulation-only → compare with Grail's CLINICAL 88.7%
 * - Every comparison must note important caveats
 * 
 * Published Reference Data (from peer-reviewed literature):
 */
const fs = require('fs');
const path = require('path');

const HEADTOHEAD_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'real_headToHead_results.json');
const CET_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'real_cet_results.json');
const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'published_comparison.json');

// ── Published Clinical Results (from peer-reviewed literature) ──
// Every number has a citation
const PUBLISHED_ASSAYS = [
  {
    name: 'Guardant360 (Guardant Health)',
    type: 'ctDNA NGS (74 genes, tumor-informed available)',
    citation: 'Odegaard 2018 Clin Cancer Res; Lanman 2015 PLoS ONE; multiple publications',
    clinical_validation: true,
    sample_size: '>200,000 clinical tests',
    sensitivity_overall: 0.853,
    specificity: 0.996,
    lod_ctdna: 0.0001, // 0.01%
    cancer_types: 50,
    too_accuracy: null, // Tumor-informed, N/A
    sequencing_depth: '5,000× (clinical standard)',
    notes: 'Commercially available. FDA-approved companion diagnostic for osimertinib. Performance in advanced/metastatic setting primarily.',
    stage_I_sensitivity: null,
    source: 'Guardant360 technical specifications and peer-reviewed publications',
  },
  {
    name: 'FoundationOne Liquid CDx (Foundation Medicine)',
    type: 'Hybrid capture NGS (324 genes)',
    citation: 'Woodhouse 2020 PLoS ONE; multiple FDA submissions',
    clinical_validation: true,
    sample_size: '>30,000 clinical tests',
    sensitivity_overall: 0.837,
    specificity: 0.995,
    lod_ctdna: 0.001, // 0.1%
    cancer_types: 50,
    too_accuracy: null, // Tumor-informed
    sequencing_depth: '~5,000×',
    notes: 'FDA-approved. Higher LOD than Guardant360. Better for tumor mutational burden assessment.',
    stage_I_sensitivity: null,
    source: 'FoundationOne CDx technical specifications',
  },
  {
    name: 'Grail Galleri (MCED)',
    type: 'Targeted methylation sequencing (>100,000 CpG regions)',
    citation: 'Jamshidi 2022 Cancer Cell; Liu 2020 Ann Oncol; Klein 2021 Ann Oncol',
    clinical_validation: true,
    sample_size: '15,254 (CCGA substudy), >140,000 (NHS-Galleri trial)',
    sensitivity_overall: 0.515, // 51.5% overall at 99.5% specificity
    specificity: 0.995,
    lod_ctdna: null, // Methylation-based, not VAF-dependent
    cancer_types: 50,
    too_accuracy: 0.887, // 88.7% tissue-of-origin accuracy
    sequencing_depth: '~30× WGBS equivalent (targeted)',
    notes: 'Commercially available ($949). NHS trial with 140,000 participants. Sensitivity increases with stage: 16.7% (I) → 90.5% (IV). Best-in-class for multi-cancer early detection breadth.',
    stage_I_sensitivity: 0.167,
    source: 'Klein 2021 Ann Oncol; Jamshidi 2022 Cancer Cell',
  },
  {
    name: 'CancerSEEK (Thrive/Exact Sciences)',
    type: 'Multi-analyte: ctDNA mutations (61 genes) + protein biomarkers (8 proteins)',
    citation: 'Cohen 2018 Science; Lennon 2020 Science',
    clinical_validation: true,
    sample_size: '1,005 cancer + 812 healthy (discovery); 10,006 women (DETECT-A)',
    sensitivity_overall: 0.70,
    specificity: 0.99,
    lod_ctdna: null, // Not reported as VAF
    cancer_types: 8,
    too_accuracy: 0.83, // 83% TOO accuracy (among detected)
    sequencing_depth: '~30,000× (targeted amplicon)',
    notes: 'DETECT-A study: 10,006 women screened, doubled the number of cancers detected by standard screening. Combined with PET-CT for confirmation.',
    stage_I_sensitivity: 0.43,
    source: 'Cohen 2018 Science; Lennon 2020 Science',
  },
  {
    name: 'DELFI (Delfi Diagnostics)',
    type: 'Genome-wide fragmentomics (low-coverage WGS)',
    citation: 'Cristiano 2019 Nature; Mathios 2021 Nat Commun; Mazzone 2024 Cancer Discov',
    clinical_validation: true,
    sample_size: '958 (lung cancer clinical validation)',
    sensitivity_overall: 0.73,
    specificity: 0.98,
    lod_ctdna: null, // Fragmentomics-based, not VAF
    cancer_types: 7,
    too_accuracy: 0.75, // Estimated from Cristiano 2019
    sequencing_depth: '1-2× WGS (low coverage)',
    notes: 'Low-cost ($100-200). No targeted enrichment needed. FDA breakthrough device for lung cancer screening.',
    stage_I_sensitivity: 0.57,
    source: 'Cristiano 2019 Nature; Mathios 2021 Nat Commun',
  },
  {
    name: 'PanSeer (Singlera Genomics)',
    type: 'Targeted methylation sequencing (595 regions)',
    citation: 'Chen 2020 Nat Commun',
    clinical_validation: true,
    sample_size: '605 asymptomatic (191 later diagnosed) + 223 cancer',
    sensitivity_overall: 0.88,
    specificity: 0.96,
    lod_ctdna: 0.00001, // 0.001% reported
    cancer_types: 5,
    too_accuracy: null,
    sequencing_depth: 'Targeted bisulfite PCR',
    notes: 'ONLY assay to demonstrate pre-symptomatic detection (up to 4 years before diagnosis) using longitudinal archived samples. GOLD STANDARD for longitudinal validation.',
    stage_I_sensitivity: null,
    pre_diagnosis_sensitivity: 0.95, // 95% up to 4 years pre-diagnosis
    source: 'Chen 2020 Nat Commun',
  },
  {
    name: 'Bie et al. 2023 (THEMIS)',
    type: 'Multi-modal: methylation + fragmentomics + CNA from single enzymatic assay',
    citation: 'Bie 2023 Nat Commun',
    clinical_validation: false,
    sample_size: '780 cancer + 497 healthy',
    sensitivity_overall: null, // Not reported as overall; 73% early-stage at 99% spec
    specificity: 0.99,
    lod_ctdna: 0.001,
    cancer_types: 7,
    too_accuracy: null,
    sequencing_depth: 'WMS (whole methylome)',
    notes: 'Academic validation only (no clinical trial). Used for head-to-head in our study.',
    stage_I_sensitivity: 0.73, // At 99% specificity
    source: 'Bie 2023 Nat Commun',
  },
];

// ── Extract DeepCatch performance from head-to-head results ──
function extractDeepCatchPerformance(headToHead) {
  const perf = {
    variant_calling: {},
    multimodal: {},
    detection_limit: headToHead.detection_limit_ctdna_fraction,
    summary: [],
  };

  if (headToHead.per_fraction_results) {
    for (const [key, result] of Object.entries(headToHead.per_fraction_results)) {
      if (result.error || !result.methods) continue;

      const dcVar = result.methods.deepcatch_variant;
      const dcMulti = result.methods.deepcatch_multimodal;

      const frac = result.ctdna_fraction;
      perf.summary.push({
        ctdna_fraction: frac,
        variant_auc: dcVar?.auc,
        variant_sens_95spec: dcVar?.sens_at_95_spec,
        variant_sens_99spec: dcVar?.sens_at_99_spec,
        multimodal_auc: dcMulti?.auc,
        multimodal_sens_95spec: dcMulti?.sens_at_95_spec,
        multimodal_sens_99spec: dcMulti?.sens_at_99_spec,
      });

      if (frac === 0.01) {
        perf.variant_calling.at_1pct = dcVar;
        perf.multimodal.at_1pct = dcMulti;
      }
      if (frac === 0.001) {
        perf.variant_calling.at_0_1pct = dcVar;
        perf.multimodal.at_0_1pct = dcMulti;
      }
      if (frac === 0.0001) {
        perf.variant_calling.at_0_01pct = dcVar;
        perf.multimodal.at_0_01pct = dcMulti;
      }
    }
  }

  return perf;
}

// ── Compare DeepCatch to Published ──
function generateComparisons(dcPerf, cetResults) {
  const comparisons = [];

  // 1. Guardant360 comparison
  comparisons.push({
    comparison: 'DeepCatch vs Guardant360',
    metric: 'Limit of Detection (ctDNA fraction)',
    guardant360: { value: 0.0001, units: 'ctDNA fraction (0.01%)', citation: 'Lanman 2015 PLoS ONE' },
    deepcatch: { value: dcPerf.detection_limit, units: 'ctDNA fraction (lowest with AUC>0.80)',
      interpretation: dcPerf.detection_limit <= 0.0001 ? 'Comparable or better than Guardant360' : 'Higher LOD (worse) than Guardant360' },
    caveats: [
      'DeepCatch LOD is simulation-based; Guardant360 LOD is clinical',
      'Guardant360 uses molecular barcoding (UMIs) with error correction; DeepCatch uses in silico error suppression',
      'Guardant360 has >200,000 clinical samples; DeepCatch has 0',
    ],
    honest_assessment: dcPerf.detection_limit <= 0.0001 ?
      'DeepCatch SIMULATION shows comparable LOD. This must be validated in real patient samples.' :
      `DeepCatch LOD (${(dcPerf.detection_limit*100).toFixed(2)}%) is ${dcPerf.detection_limit/0.0001}x higher than Guardant360 clinical LOD (0.01%). This is a significant gap.`,
  });

  // 2. Grail Galleri comparison
  const dcSensAt99Spec = dcPerf.summary.find(s => s.multimodal_sens_99spec > 0)?.multimodal_sens_99spec || 0;
  comparisons.push({
    comparison: 'DeepCatch vs Grail Galleri (MCED)',
    metric: 'Sensitivity at 99.5% Specificity',
    grail: { value: 0.515, units: 'Overall sensitivity', citation: 'Klein 2021 Ann Oncol; Jamshidi 2022 Cancer Cell' },
    deepcatch_multimodal: { value: dcSensAt99Spec, units: 'Estimated sensitivity at 99% spec at optimal ctDNA fraction', note: 'From head-to-head simulation' },
    caveats: [
      'Grail Galleri is methylation-based (proprietary); DeepCatch is mutation + multi-modal',
      'Grail has clinical data from 15,254-subject CCGA study; DeepCatch has simulation only',
      'Grail Galleri is FDA breakthrough device and commercially available; DeepCatch is a research concept',
      'Grail detected 51.5% at 99.5% specificity across 50+ cancer types; DeepCatch simulation covers 8 types',
      'Grail TOO accuracy: 88.7% (CLINICAL); DeepCatch TOO: SIMULATION ONLY',
    ],
    honest_assessment: 'DIRECT COMPARISON NOT POSSIBLE: Grail Galleri is a validated clinical test with >15,000 patients. DeepCatch has zero clinical patients. Any comparison would be apples-to-oranges. DeepCatch\'s simulation performance cannot be compared to clinical reality without wet-lab validation.',
  });

  // 3. TOO comparison
  comparisons.push({
    comparison: 'Tissue-of-Origin (TOO) Accuracy',
    metric: 'TOO Accuracy',
    grail: { value: 0.887, units: '88.7% (CLINICAL)', citation: 'Jamshidi 2022 Cancer Cell' },
    cancerseeek: { value: 0.83, units: '83% (CLINICAL)', citation: 'Cohen 2018 Science' },
    deepcatch: { value: null, units: 'NOT MEASURED ON REAL DATA',
      note: 'DeepCatch TOO was only validated on simulated data. Cannot report TOO accuracy without real multi-class labeled data.' },
    caveats: [
      'Grail: Clinical TOO across 50+ cancer types with 88.7% accuracy',
      'CancerSEEK: Clinical TOO across 8 cancer types with 83% accuracy',
      'DeepCatch: TOO simulation used synthetic data with known ground truth — meaningless for Nature-level publication',
    ],
    honest_assessment: '❌ DEEPCATCH TOO IS NOT PROVEN: Previous TOO validation used simulation-only data. Real TOO accuracy on heterogeneous cancer types with mixed sample quality is unknown. This is a critical gap for publication.',
  });

  // 4. CET vs PanSeer longitudinal
  if (cetResults) {
    const cSens = cetResults.performance?.sensitivity || 0;
    const cSpec = cetResults.performance?.specificity_overall || 0;

    comparisons.push({
      comparison: 'DeepCatch CET vs PanSeer (Longitudinal)',
      metric: 'Pre-diagnosis Detection',
      panseer: { value: 0.95, units: 'Sensitivity 1-4 years pre-diagnosis', citation: 'Chen 2020 Nat Commun' },
      deepcatch_cet: {
        sensitivity: cSens,
        specificity: cSpec,
        median_detection_days: cetResults.performance?.median_detection_days,
        target_met: cetResults.targets?.both_met || false,
        note: 'Gompertz growth model + 8 quarterly timepoints, 700 simulated patients',
      },
      caveats: [
        'PanSeer: REAL patient data from Taizhou Longitudinal Study (123,115 subjects, up to 4 years pre-diagnosis)',
        'DeepCatch CET: SIMULATION only with Gompertz growth model',
        'PanSeer used archived blood samples collected years before cancer diagnosis',
        'DeepCatch has zero longitudinal patient data',
      ],
      honest_assessment: cetResults.targets?.both_met ?
        `DeepCatch CET simulation MET the dual target (sens≥70% spec≥95%). However, this is simulated — PanSeer demonstrated 95% pre-diagnosis sensitivity in REAL patients. Without longitudinal clinical validation, DeepCatch CET remains a theoretical concept only.` :
        `DeepCatch CET simulation did NOT meet targets (sens=${(cSens*100).toFixed(0)}%, spec=${(cSpec*100).toFixed(0)}%). PanSeer achieved 95% pre-diagnosis sensitivity in REAL patients. DeepCatch CET is not competitive.`,
    });
  }

  // 5. Sequencing cost/depth comparison
  comparisons.push({
    comparison: 'Sequencing Requirements',
    metric: 'Cost and Depth',
    assays: {
      'Guardant360': { depth: '5,000×', cost: '$5,800 (clinical)', coverage: '74 genes' },
      'FoundationOne Liquid': { depth: '5,000×', cost: '$5,800 (clinical)', coverage: '324 genes' },
      'Grail Galleri': { depth: '~30× WGBS', cost: '$949', coverage: '>100,000 CpG regions' },
      'DELFI': { depth: '1-2× WGS', cost: '$100-200', coverage: 'Genome-wide fragments' },
      'DeepCatch (projected)': { depth: '50,000×', cost: 'Unknown (academic)', coverage: '>50 genes + multi-modal' },
    },
    caveats: [
      'DeepCatch requires 10× higher sequencing depth than Guardant360 (50,000× vs 5,000×)',
      'This 10× depth requirement may not be clinically or economically viable',
      'DeepCatch\'s cost advantage (if any) comes from using cheaper in silico error suppression vs molecular barcoding',
      'Commercial assays benefit from economies of scale (millions of samples); DeepCatch does not',
    ],
    honest_assessment: '⚠️ DEEPCATCH REQUIRES 10× MORE SEQUENCING: At 50,000× depth vs Guardant360\'s clinical 5,000×, DeepCatch\'s LOD advantage may be partly or entirely attributable to increased sequencing depth, not better algorithms. A fair comparison must hold sequencing depth constant.',
  });

  return comparisons;
}

// ── Generate Final Table ──
function generateSummaryTable(dcPerf, cetResults) {
  const rows = [];

  PUBLISHED_ASSAYS.forEach(assay => {
    rows.push({
      assay: assay.name,
      type: assay.type,
      clinical_validation: assay.clinical_validation,
      sample_size: assay.sample_size,
      sensitivity: assay.sensitivity_overall,
      specificity: assay.specificity,
      lod_ctdna: assay.lod_ctdna,
      too_accuracy: assay.too_accuracy,
      cancer_types: assay.cancer_types,
      depth: assay.sequencing_depth,
      notes: assay.notes,
    });
  });

  // Add DeepCatch rows
  const bestDcFrac = dcPerf.summary.find(s => s.ctdna_fraction === 0.001) || 
                     dcPerf.summary[0] || {};

  rows.push({
    assay: 'DeepCatch (variant calling)',
    type: 'Weighted multi-gene variant calling with trinucleotide error suppression',
    clinical_validation: false,
    sample_size: '0 clinical, simulation only',
    sensitivity: bestDcFrac.variant_sens_95spec || null,
    specificity: bestDcFrac.variant_sens_95spec !== undefined ? 'Simulated (at target specificity)' : 'Not validated on real data',
    lod_ctdna: dcPerf.detection_limit,
    too_accuracy: null,
    cancer_types: 8,
    depth: '50,000× (simulation)',
    notes: '❌ SIMULATION ONLY. No clinical validation. Cannot be directly compared to clinical assays.',
  });

  rows.push({
    assay: 'DeepCatch (multi-modal fusion)',
    type: 'Performance-weighted fusion: variant calling + fragmentomics + methylation',
    clinical_validation: false,
    sample_size: '0 clinical, simulation only',
    sensitivity: bestDcFrac.multimodal_sens_95spec || null,
    specificity: bestDcFrac.multimodal_sens_95spec !== undefined ? 'Simulated (at target specificity)' : 'Not validated',
    lod_ctdna: dcPerf.detection_limit,
    too_accuracy: 'SIMULATION ONLY — not validated',
    cancer_types: 8,
    depth: '50,000× (simulation)',
    notes: '❌ SIMULATION ONLY. AUC from performance-weighted fusion of simulated modalities. Not clinically validated.',
  });

  if (cetResults) {
    rows.push({
      assay: 'DeepCatch CET (longitudinal)',
      type: 'Hierarchical Bayes Cumulative Evidence Tracking with Gompertz growth',
      clinical_validation: false,
      sample_size: '700 simulated patients (0 clinical)',
      sensitivity: cetResults.performance?.sensitivity,
      specificity: cetResults.performance?.specificity_overall,
      lod_ctdna: null,
      too_accuracy: null,
      cancer_types: 8,
      depth: 'N/A (longitudinal)',
      notes: `❌ SIMULATION ONLY. ${cetResults.targets?.both_met ? 'Meets dual target in simulation' : 'Does not meet targets in simulation'}. No longitudinal patient data.`,
    });
  }

  return rows;
}

// ── MAIN ──
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH REAL-DATA VALIDATION — PHASE 5: Published Comparison');
  console.log('='.repeat(70));
  console.log();

  let dcPerf = null;
  let cetResults = null;

  // Load head-to-head results
  try {
    const h2h = JSON.parse(fs.readFileSync(HEADTOHEAD_PATH, 'utf8'));
    dcPerf = extractDeepCatchPerformance(h2h);
    console.log('✅ Loaded head-to-head results for DeepCatch performance extraction');
  } catch (err) {
    console.log(`⚠️  Could not load head-to-head results: ${err.message}`);
    dcPerf = { detection_limit: null, summary: [], variant_calling: {}, multimodal: {} };
  }

  // Load CET results
  try {
    cetResults = JSON.parse(fs.readFileSync(CET_PATH, 'utf8'));
    console.log('✅ Loaded CET results');
    console.log(`   CET: sens=${(cetResults.performance?.sensitivity * 100).toFixed(1)}%, spec=${(cetResults.performance?.specificity_overall * 100).toFixed(1)}%, AUC=${cetResults.performance?.auc?.toFixed(4)}`);
  } catch (err) {
    console.log(`⚠️  Could not load CET results: ${err.message}`);
  }

  console.log();

  // Generate comparisons
  const comparisons = generateComparisons(dcPerf, cetResults);
  const summaryTable = generateSummaryTable(dcPerf, cetResults);

  // ── HONEST VERDICT ──
  let honestAssessment;
  if (!dcPerf.detection_limit) {
    honestAssessment = '❌ CANNOT ASSESS: DeepCatch head-to-head results not available. No comparison possible.';
  } else if (dcPerf.detection_limit > 0.001) {
    honestAssessment = `❌ DEEPCATCH LOD IS TOO HIGH: Detection limit (${(dcPerf.detection_limit*100).toFixed(2)}% ctDNA) is ${(dcPerf.detection_limit/0.0001).toFixed(0)}× worse than Guardant360 clinical LOD (0.01%). Without wet-lab validation showing comparable or better LOD, DeepCatch cannot claim clinical utility. The multi-modal fusion advantage at higher ctDNA fractions is not sufficient to overcome poor sensitivity at clinically relevant ctDNA fractions.`;
  } else {
    honestAssessment = `⚠️ PARTIALLY PROMISING: DeepCatch SIMULATION shows competitive LOD (${(dcPerf.detection_limit*100).toFixed(2)}% ctDNA) but this is simulation-only. The gap between simulation and clinical reality is large — factors like sample degradation, PCR artifacts, GC bias, and inter-patient variability are not fully captured. DeepCatch requires: (1) wet-lab validation on real patient samples, (2) head-to-head comparison against Guardant360 or CAPP-Seq on the same patient samples, (3) demonstration that performance advantage persists at matched sequencing depth.`;
  }

  const output = {
    metadata: {
      generated: new Date().toISOString(),
      published_assay_references: PUBLISHED_ASSAYS.map(a => ({ name: a.name, citation: a.citation })),
      comparison_caveats: 'ALL DeepCatch results are SIMULATION-BASED. No clinical validation has been performed. Direct comparison to clinical assays is NOT scientifically valid — this comparison is provided for context only.',
    },
    published_assays: PUBLISHED_ASSAYS,
    deepcatch_extracted_performance: dcPerf,
    deepcatch_cet_performance: cetResults?.performance || null,
    head_to_head_comparisons: comparisons,
    summary_table: summaryTable,
    honest_assessment: honestAssessment,
    requirements_for_validation: [
      '1. Test DeepCatch on real patient plasma samples (n ≥ 200 cancer, n ≥ 200 healthy)',
      '2. Head-to-head on same samples against Guardant360 or CAPP-Seq',
      '3. Match sequencing depth to clinical standard (5,000×) for fair comparison',
      '4. Validate TOO on multi-class real data with known primary tumors',
      '5. Longitudinal cohort: ≥500 patients with serial blood draws over ≥2 years',
      '6. Independent validation at a separate institution',
      '7. Pre-register analysis plan to prevent p-hacking',
    ],
    publication_readiness: {
      can_publish_as_commentary: true,
      can_publish_as_methods_paper: 'Only if wet-lab validation is added',
      can_publish_as_clinical_validation: false,
      simulation_only: true,
      next_step: 'Partner with clinical collaborators for real sample validation',
    },
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));

  console.log('='.repeat(70));
  console.log('HONEST ASSESSMENT');
  console.log('='.repeat(70));
  console.log();
  console.log(honestAssessment);
  console.log();
  console.log('📋 Requirements for clinical validation:');
  output.requirements_for_validation.forEach(r => console.log(`   ${r}`));
  console.log();
  console.log(`💾 Saved published comparison to ${path.basename(OUTPUT_PATH)}`);
  console.log('\n✅ Phase 5 complete.');
  console.log('='.repeat(70));
})();
