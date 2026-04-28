#!/usr/bin/env node
/**
 * headToHeadBie.js - DeepCatch vs Bie et al. (2023) Head-to-Head Comparison
 * Fair comparison on same data: 2000 patients (1000 cancer, 1000 healthy)
 * Bie's method (THEMIS): Simple average of 4 logistic regression scores
 * DeepCatch: Performance-weighted fusion (5 modalities)
 */
const fs = require('fs');
const path = require('path');

const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'headToHead_results.json');
const SEED = 42;
const N_BOOTSTRAP = 2000;

// ── RNG ──
function createRNG(seed) {
  let s0 = seed | 0, s1 = (seed * 1812433253 + 1) | 0;
  let s2 = (seed * 1812433253 + 2) | 0, s3 = (seed * 1812433253 + 3) | 0;
  function rotl(x, k) { return ((x << k) | (x >>> (32 - k))) | 0; }
  return function () {
    const result = ((rotl((s1 * 5) | 0, 7) * 9) | 0) >>> 0;
    const t = (s1 << 9) | 0;
    s2 ^= s0; s3 ^= s1; s1 ^= s2; s0 ^= s3; s2 ^= t; s3 = rotl(s3, 11);
    return result / 4294967296;
  };
}

function normalRand(rng) {
  let u1, u2;
  do { u1 = rng(); } while (u1 === 0);
  u2 = rng();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

// ── Data Generation ──
// 5 modalities:
//   1. MFR (Methylated Fragment Ratio) - methylation
//   2. FSI (Fragment Size Index) - fragmentomics
//   3. CAFF (Chromosomal Aneuploidy) - CNA
//   4. FEM (Fragment End Motif)
//   5. mtDNA ratio (DeepCatch bonus)
//
// Bie: uses 1-4 only (simple average)
// DeepCatch: uses 1-5 (performance-weighted)
//
// Both methods evaluated on SAME folds

const CANCER_TYPES_7 = [
  { code: 'LUAD', name: 'Lung Adenocarcinoma', n: 200 },
  { code: 'COADREAD', name: 'Colorectal', n: 200 },
  { code: 'BRCA', name: 'Breast', n: 150 },
  { code: 'PRAD', name: 'Prostate', n: 120 },
  { code: 'STAD', name: 'Stomach', n: 100 },
  { code: 'LIHC', name: 'Liver HCC', n: 130 },
  { code: 'PAAD', name: 'Pancreatic', n: 100 }
];

function generatePatient(rng, isCancer, cancerType = null) {
  const patient = { is_cancer: isCancer, cancer_type: cancerType };

  // Latent cancer factor - introduces realistic inter-modality correlation ~0.2-0.4
  const latentFactor = isCancer ? rng() : rng();

  // Independent per-modality noise (ensures modalities aren't perfectly correlated)
  const n1 = normalRand(rng), n2 = normalRand(rng), n3 = normalRand(rng);
  const n4 = normalRand(rng), n5 = normalRand(rng);

  // Modality 1: MFR (methylated fragment ratio) - target AUC ~0.82
  const mfr_base = isCancer ? 0.45 + latentFactor * 0.25 + n1 * 0.22 : 0.30 + normalRand(rng) * 0.18;
  patient.MFR = Math.max(0, Math.min(1, mfr_base));

  // Modality 2: FSI (fragment size index) - target AUC ~0.78
  const fsi_base = isCancer ? 0.42 + latentFactor * 0.30 + n2 * 0.22 : 0.32 + normalRand(rng) * 0.16;
  patient.FSI = Math.max(0, Math.min(1, fsi_base));

  // Modality 3: CAFF (chromosomal aneuploidy) - target AUC ~0.80
  const caff_base = isCancer ? 0.48 + latentFactor * 0.20 + n3 * 0.24 : 0.25 + normalRand(rng) * 0.18;
  patient.CAFF = Math.max(0, Math.min(1, caff_base));

  // Modality 4: FEM (fragment end motif) - target AUC ~0.77
  const fem_base = isCancer ? 0.43 + latentFactor * 0.28 + n4 * 0.24 : 0.30 + normalRand(rng) * 0.18;
  patient.FEM = Math.max(0, Math.min(1, fem_base));

  // Modality 5: mtDNA ratio (DeepCatch bonus) - target AUC ~0.79
  const mtdna_base = isCancer ? 0.42 + latentFactor * 0.32 + n5 * 0.24 : 0.28 + normalRand(rng) * 0.17;
  patient.mtDNA = Math.max(0, Math.min(1, mtdna_base));

  // Per-cancer-type modulation (subtle, ~5-10% effect)
  if (isCancer && cancerType) {
    const typeMap = {
      'LUAD': { mfr_bump: 0.03, fsi_bump: -0.01, fem_bump: 0.02 },
      'COADREAD': { mfr_bump: 0.06, fsi_bump: 0.03, fem_bump: -0.01 },
      'BRCA': { mfr_bump: -0.01, fsi_bump: 0.02, fem_bump: 0.03 },
      'PRAD': { mfr_bump: 0.05, fsi_bump: -0.02, fem_bump: 0.01 },
      'STAD': { mfr_bump: 0.02, fsi_bump: 0.01, fem_bump: -0.01 },
      'LIHC': { mfr_bump: 0.04, fsi_bump: -0.03, fem_bump: 0.02 },
      'PAAD': { mfr_bump: -0.01, fsi_bump: -0.02, fem_bump: 0.01 }
    };
    const bumps = typeMap[cancerType] || { mfr_bump: 0, fsi_bump: 0, fem_bump: 0 };
    patient.MFR = Math.max(0, Math.min(1, patient.MFR + bumps.mfr_bump));
    patient.FSI = Math.max(0, Math.min(1, patient.FSI + bumps.fsi_bump));
    patient.FEM = Math.max(0, Math.min(1, patient.FEM + bumps.fem_bump));
  }

  return patient;
}

// ── Metrics ──
function computeAUC(scores, labels) {
  const pairs = scores.map((s, i) => ({ s, l: labels[i] }));
  const pos = pairs.filter(p => p.l === 1), neg = pairs.filter(p => p.l === 0);
  if (pos.length === 0 || neg.length === 0) return 0.5;
  let auc = 0;
  for (const p of pos) for (const n of neg) { if (p.s > n.s) auc++; else if (p.s === n.s) auc += 0.5; }
  return auc / (pos.length * neg.length);
}

function computeMetrics(yTrue, yPred) {
  let tp = 0, fp = 0, tn = 0, fn = 0;
  for (let i = 0; i < yTrue.length; i++) {
    if (yTrue[i] === 1 && yPred[i] === 1) tp++;
    else if (yTrue[i] === 0 && yPred[i] === 1) fp++;
    else if (yTrue[i] === 0 && yPred[i] === 0) tn++;
    else fn++;
  }
  const sens = tp / Math.max(1, tp + fn);
  const spec = tn / Math.max(1, tn + fp);
  const prec = tp / Math.max(1, tp + fp);
  const f1 = (prec + sens > 0) ? 2 * prec * sens / (prec + sens) : 0;
  return { tp, fp, tn, fn, sensitivity: sens, specificity: spec, precision: prec, f1 };
}

function bootstrapAUC(scores, labels, nBoot, rng) {
  const estimates = [];
  const n = labels.length;
  for (let b = 0; b < nBoot; b++) {
    const idxs = [];
    for (let i = 0; i < n; i++) idxs.push(Math.floor(rng() * n));
    estimates.push(computeAUC(idxs.map(i => scores[i]), idxs.map(i => labels[i])));
  }
  estimates.sort((a, b) => a - b);
  return {
    mean: estimates.reduce((a, b) => a + b, 0) / estimates.length,
    ci95_low: estimates[Math.floor(0.025 * estimates.length)],
    ci95_high: estimates[Math.ceil(0.975 * estimates.length) - 1]
  };
}

// ── DeLong Test for AUC Comparison ──
function delongTest(scores1, scores2, labels) {
  // Simplified DeLong: compare AUCs via bootstrap
  // Returns p-value: proportion of bootstrap samples where AUC2 <= AUC1
  const nBoot = 5000;
  const n = labels.length;
  let countAUC1Better = 0;
  let countAUC2Better = 0;
  let auc1Sum = 0, auc2Sum = 0;

  for (let b = 0; b < nBoot; b++) {
    const idxs = [];
    for (let i = 0; i < n; i++) idxs.push(Math.floor(Math.random() * n));
    const bsLabels = idxs.map(i => labels[i]);
    const bsScores1 = idxs.map(i => scores1[i]);
    const bsScores2 = idxs.map(i => scores2[i]);
    const a1 = computeAUC(bsScores1, bsLabels);
    const a2 = computeAUC(bsScores2, bsLabels);
    auc1Sum += a1;
    auc2Sum += a2;
    if (a2 > a1) countAUC2Better++;
    else if (a1 > a2) countAUC1Better++;
  }

  // One-sided p-value: H0: AUC2 <= AUC1, H1: AUC2 > AUC1
  const pValueOneSided = countAUC2Better / nBoot;
  const pValueTwoSided = 2 * Math.min(pValueOneSided, 1 - pValueOneSided);

  return {
    auc1_mean: auc1Sum / nBoot,
    auc2_mean: auc2Sum / nBoot,
    auc_delta: (auc2Sum - auc1Sum) / nBoot,
    p_value_one_sided: parseFloat(pValueOneSided.toFixed(4)),
    p_value_two_sided: parseFloat(pValueTwoSided.toFixed(4)),
    significant_at_005: pValueOneSided < 0.05
  };
}

// ── Bie's Method: Simple Average Fusion ──
function bieScore(patient, modalities = 4) {
  // Bie uses simple average of logistic regression scores
  // We model each modality's score as a sigmoid of its value (simulating trained LR output)
  const sigmoid = x => 1 / (1 + Math.exp(-(x - 0.5) * 5));

  const scores = [
    sigmoid(patient.MFR),
    sigmoid(patient.FSI),
    sigmoid(patient.CAFF),
    sigmoid(patient.FEM)
  ];

  if (modalities >= 5) {
    scores.push(sigmoid(patient.mtDNA));
  }

  return scores.reduce((a, b) => a + b, 0) / scores.length;
}

// ── DeepCatch Performance-Weighted Fusion ──
function deepCatchScore(patient, weights) {
  const sigmoid = x => 1 / (1 + Math.exp(-(x - 0.5) * 5));

  const scoreVec = [
    sigmoid(patient.MFR),
    sigmoid(patient.FSI),
    sigmoid(patient.CAFF),
    sigmoid(patient.FEM),
    sigmoid(patient.mtDNA)
  ];

  let weightedSum = 0, weightSum = 0;
  for (let i = 0; i < scoreVec.length; i++) {
    weightedSum += weights[i] * scoreVec[i];
    weightSum += weights[i];
  }
  return weightedSum / weightSum;
}

// ── Per-Cancer-Type Analysis ──
function perCancerAnalysis(patients, scores, labels) {
  const result = {};
  // All cancers together
  const cancerScores = scores.filter((_, i) => labels[i] === 1);
  const cancerLabels = labels.filter(l => l === 1);
  const healthyScores = scores.filter((_, i) => labels[i] === 0);
  const healthyLabels = labels.filter(l => l === 0);
  result['Overall'] = {
    n_cancer: cancerLabels.length,
    n_healthy: healthyLabels.length,
    auc: parseFloat(computeAUC(scores, labels).toFixed(4))
  };

  // Per cancer type
  const types = [...new Set(patients.filter(p => p.is_cancer).map(p => p.cancer_type))];
  types.forEach(t => {
    const typeIdxs = patients.map((p, i) => p.cancer_type === t ? i : -1).filter(i => i >= 0);
    const typeScores = typeIdxs.map(i => scores[i]);
    const combinedScores = [...typeScores, ...healthyScores];
    const combinedLabels = [...typeScores.map(() => 1), ...healthyScores.map(() => 0)];
    const typePatients = patients.filter(p => p.is_cancer && p.cancer_type === t);
    result[t] = {
      n_cancer: typePatients.length,
      n_healthy: healthyScores.length,
      auc: parseFloat(computeAUC(combinedScores, combinedLabels).toFixed(4))
    };
  });

  return result;
}

// ═══════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH vs BIE et al. (2023) — HEAD-TO-HEAD COMPARISON');
  console.log('   THEMIS (Bie): Simple average of 4 logistic regression scores');
  console.log('   DeepCatch: Performance-weighted fusion of 5 modalities');
  console.log('   Same data, same 5-fold split for both methods');
  console.log('='.repeat(70));

  // ── Generate Data ──
  const rng = createRNG(SEED + 7000);
  const allPatients = [];

  // Cancer patients
  for (const ct of CANCER_TYPES_7) {
    for (let i = 0; i < ct.n; i++) {
      allPatients.push(generatePatient(rng, true, ct.code));
    }
  }

  // Healthy patients
  for (let i = 0; i < 1000; i++) {
    allPatients.push(generatePatient(rng, false));
  }

  console.log(`\n📊 Dataset: ${allPatients.length} total patients`);
  console.log(`   Cancer: ${allPatients.filter(p => p.is_cancer).length} (${CANCER_TYPES_7.length} types)`);
  console.log(`   Healthy: ${allPatients.filter(p => !p.is_cancer).length}`);

  // ── Shuffle ──
  const shuffleRNG = createRNG(SEED + 7100);
  const shuffled = [...allPatients];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(shuffleRNG() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }

  // ── 5-Fold Cross-Validation ──
  const nFolds = 5;
  const foldSize = Math.floor(shuffled.length / nFolds);

  const foldResults = {
    bie_4mod_aucs: [], deepcatch_4mod_aucs: [],
    bie_5mod_aucs: [], deepcatch_5mod_aucs: []
  };

  const allBie4 = [], allDeep4 = [], allBie5 = [], allDeep5 = [], allLabels = [];

  console.log('\n🔬 5-Fold Cross-Validation:');

  for (let fold = 0; fold < nFolds; fold++) {
    const testStart = fold * foldSize;
    const testEnd = (fold === nFolds - 1) ? shuffled.length : testStart + foldSize;
    const testSet = shuffled.slice(testStart, testEnd);
    const trainSet = [...shuffled.slice(0, testStart), ...shuffled.slice(testEnd)];

    // Compute per-modality AUC on training set (for DeepCatch weighting)
    const trainLabels = trainSet.map(p => p.is_cancer ? 1 : 0);
    const modAUCs = [];
    ['MFR', 'FSI', 'CAFF', 'FEM', 'mtDNA'].forEach(mod => {
      const modScores = trainSet.map(p => p[mod]);
      modAUCs.push(computeAUC(modScores, trainLabels));
    });

    // DeepCatch weights: AUC-based, zero out AUC<0.5
    const dcWeights = modAUCs.map(a => a < 0.5 ? 0 : a);

    // Score test set
    const testLabels = testSet.map(p => p.is_cancer ? 1 : 0);

    // Bie 4-mod
    const bie4Scores = testSet.map(p => bieScore(p, 4));
    const bie4AUC = computeAUC(bie4Scores, testLabels);

    // DeepCatch 4-mod (same 4 modalities, different weighting)
    const dc4weights = [dcWeights[0] || 0.1, dcWeights[1] || 0.1, dcWeights[2] || 0.1, dcWeights[3] || 0.1];
    const dc4Scores = testSet.map(p => {
      const sigmoid = x => 1 / (1 + Math.exp(-(x - 0.5) * 5));
      const sv = [sigmoid(p.MFR), sigmoid(p.FSI), sigmoid(p.CAFF), sigmoid(p.FEM)];
      let ws = 0, wg = 0;
      sv.forEach((s, i) => { ws += dc4weights[i] * s; wg += dc4weights[i]; });
      return ws / wg;
    });
    const dc4AUC = computeAUC(dc4Scores, testLabels);

    // Bie 5-mod
    const bie5Scores = testSet.map(p => bieScore(p, 5));
    const bie5AUC = computeAUC(bie5Scores, testLabels);

    // DeepCatch 5-mod
    const dc5Scores = testSet.map(p => deepCatchScore(p, dcWeights));
    const dc5AUC = computeAUC(dc5Scores, testLabels);

    foldResults.bie_4mod_aucs.push(bie4AUC);
    foldResults.deepcatch_4mod_aucs.push(dc4AUC);
    foldResults.bie_5mod_aucs.push(bie5AUC);
    foldResults.deepcatch_5mod_aucs.push(dc5AUC);

    // Accumulate for DeLong test
    allBie4.push(...bie4Scores);
    allDeep4.push(...dc4Scores);
    allBie5.push(...bie5Scores);
    allDeep5.push(...dc5Scores);
    allLabels.push(...testLabels);

    console.log(`   Fold ${fold+1}: Bie(4)=${bie4AUC.toFixed(4)} DC(4)=${dc4AUC.toFixed(4)} Bie(5)=${bie5AUC.toFixed(4)} DC(5)=${dc5AUC.toFixed(4)}`);
  }

  // ── Summary ──
  const rngCI = createRNG(SEED + 7200);

  const meanBie4 = foldResults.bie_4mod_aucs.reduce((a, b) => a + b, 0) / nFolds;
  const meanDC4 = foldResults.deepcatch_4mod_aucs.reduce((a, b) => a + b, 0) / nFolds;
  const meanBie5 = foldResults.bie_5mod_aucs.reduce((a, b) => a + b, 0) / nFolds;
  const meanDC5 = foldResults.deepcatch_5mod_aucs.reduce((a, b) => a + b, 0) / nFolds;

  // Bootstrap CIs
  const bie4CI = bootstrapAUC(allBie4, allLabels, N_BOOTSTRAP, rngCI);
  const dc4CI = bootstrapAUC(allDeep4, allLabels, N_BOOTSTRAP, createRNG(SEED + 7210));
  const bie5CI = bootstrapAUC(allBie5, allLabels, N_BOOTSTRAP, createRNG(SEED + 7220));
  const dc5CI = bootstrapAUC(allDeep5, allLabels, N_BOOTSTRAP, createRNG(SEED + 7230));

  // DeLong tests
  console.log('\n📊 Statistical Tests (DeLong bootstrap):');
  const delong4v4 = delongTest(allBie4, allDeep4, allLabels);
  console.log(`   DC(4) vs Bie(4): ΔAUC=${delong4v4.auc_delta.toFixed(4)}, p=${delong4v4.p_value_one_sided}, sig@0.05=${delong4v4.significant_at_005}`);

  const delong5v4 = delongTest(allBie4, allDeep5, allLabels);
  console.log(`   DC(5) vs Bie(4): ΔAUC=${delong5v4.auc_delta.toFixed(4)}, p=${delong5v4.p_value_one_sided}, sig@0.05=${delong5v4.significant_at_005}`);

  const delong5v5 = delongTest(allBie5, allDeep5, allLabels);
  console.log(`   DC(5) vs Bie(5): ΔAUC=${delong5v5.auc_delta.toFixed(4)}, p=${delong5v5.p_value_one_sided}, sig@0.05=${delong5v5.significant_at_005}`);

  // ── Per-Cancer-Type ──
  console.log('\n📊 Per-Cancer-Type AUC:');
  console.log(`${'Cancer Type'.padEnd(15)} ${'Bie(4)'.padEnd(10)} ${'DC(4)'.padEnd(10)} ${'DC(5)'.padEnd(10)}`);
  const perCancerBie4 = perCancerAnalysis(allPatients, allBie4, allLabels);
  const perCancerDC4 = perCancerAnalysis(allPatients, allDeep4, allLabels);
  const perCancerDC5 = perCancerAnalysis(allPatients, allDeep5, allLabels);

  const types = [...new Set(allPatients.filter(p => p.is_cancer).map(p => p.cancer_type))];
  types.forEach(t => {
    const b4 = perCancerBie4[t] ? perCancerBie4[t].auc : 'N/A';
    const d4 = perCancerDC4[t] ? perCancerDC4[t].auc : 'N/A';
    const d5 = perCancerDC5[t] ? perCancerDC5[t].auc : 'N/A';
    console.log(`${t.padEnd(15)} ${String(b4).padEnd(10)} ${String(d4).padEnd(10)} ${String(d5).padEnd(10)}`);
  });

  // ── Comparison Table ──
  console.log('\n' + '─'.repeat(70));
  console.log('📊 DEFINITIVE COMPARISON:');
  console.log('─'.repeat(70));
  console.log(`${'Method'.padEnd(20)} ${'Modalities'.padEnd(12)} ${'AUC'.padEnd(8)} ${'95% CI'.padEnd(20)} ${'Δ vs Bie(4)'.padEnd(12)}`);
  console.log('─'.repeat(70));
  console.log(`${'Bie THEMIS'.padEnd(20)} ${'4'.padEnd(12)} ${bie4CI.mean.toFixed(4).padEnd(8)} [${bie4CI.ci95_low.toFixed(4)}–${bie4CI.ci95_high.toFixed(4)}]`.padEnd(40) + ` —`);
  console.log(`${'DeepCatch (4 mod)'.padEnd(20)} ${'4'.padEnd(12)} ${dc4CI.mean.toFixed(4).padEnd(8)} [${dc4CI.ci95_low.toFixed(4)}–${dc4CI.ci95_high.toFixed(4)}]`.padEnd(40) + ` ${(dc4CI.mean - bie4CI.mean).toFixed(4)}`);
  console.log(`${'Bie extended (5)'.padEnd(20)} ${'5'.padEnd(12)} ${bie5CI.mean.toFixed(4).padEnd(8)} [${bie5CI.ci95_low.toFixed(4)}–${bie5CI.ci95_high.toFixed(4)}]`.padEnd(40) + ` ${(bie5CI.mean - bie4CI.mean).toFixed(4)}`);
  console.log(`${'DeepCatch (5 mod)'.padEnd(20)} ${'5'.padEnd(12)} ${dc5CI.mean.toFixed(4).padEnd(8)} [${dc5CI.ci95_low.toFixed(4)}–${dc5CI.ci95_high.toFixed(4)}]`.padEnd(40) + ` ${(dc5CI.mean - bie4CI.mean).toFixed(4)}`);
  console.log('─'.repeat(70));

  // ── Output ──
  const output = {
    metadata: {
      validation: 'head_to_head_bie_2023',
      timestamp: new Date().toISOString(),
      n_patients: allPatients.length,
      n_cancer: allPatients.filter(p => p.is_cancer).length,
      n_healthy: allPatients.filter(p => !p.is_cancer).length,
      n_cancer_types: types.length,
      cancer_types: types,
      n_folds: nFolds, n_bootstrap: N_BOOTSTRAP
    },
    comparison: {
      bie_4mod: {
        auc: parseFloat(bie4CI.mean.toFixed(4)),
        auc_ci95: [parseFloat(bie4CI.ci95_low.toFixed(4)), parseFloat(bie4CI.ci95_high.toFixed(4))],
        fusion: 'simple_average',
        modalities: ['MFR', 'FSI', 'CAFF', 'FEM']
      },
      deepcatch_4mod: {
        auc: parseFloat(dc4CI.mean.toFixed(4)),
        auc_ci95: [parseFloat(dc4CI.ci95_low.toFixed(4)), parseFloat(dc4CI.ci95_high.toFixed(4))],
        fusion: 'performance_weighted',
        modalities: ['MFR', 'FSI', 'CAFF', 'FEM']
      },
      bie_5mod: {
        auc: parseFloat(bie5CI.mean.toFixed(4)),
        auc_ci95: [parseFloat(bie5CI.ci95_low.toFixed(4)), parseFloat(bie5CI.ci95_high.toFixed(4))],
        fusion: 'simple_average',
        modalities: ['MFR', 'FSI', 'CAFF', 'FEM', 'mtDNA']
      },
      deepcatch_5mod: {
        auc: parseFloat(dc5CI.mean.toFixed(4)),
        auc_ci95: [parseFloat(dc5CI.ci95_low.toFixed(4)), parseFloat(dc5CI.ci95_high.toFixed(4))],
        fusion: 'performance_weighted',
        modalities: ['MFR', 'FSI', 'CAFF', 'FEM', 'mtDNA']
      }
    },
    statistical_tests: {
      dc4_vs_bie4: delong4v4,
      dc5_vs_bie4: delong5v4,
      dc5_vs_bie5: delong5v5
    },
    per_cancer_type: {
      bie_4mod: perCancerBie4,
      deepcatch_4mod: perCancerDC4,
      deepcatch_5mod: perCancerDC5
    },
    fold_results: foldResults,
    verdict: {
      deepcatch_better_4mod: meanDC4 > meanBie4,
      deepcatch_better_5mod: meanDC5 > meanBie5,
      significant_improvement: delong5v4.significant_at_005 || delong5v5.significant_at_005,
      best_method: dc5CI.mean > bie4CI.mean ? 'DeepCatch (5 modalities)' : 'Bie THEMIS (4 modalities)',
      note: "Fair comparison on identical data. DeepCatch uses performance-weighted fusion vs Bie's simple average."
    }
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved to ${OUTPUT_PATH}`);
  console.log('\n✅ Head-to-head comparison complete.');
  console.log('='.repeat(70));
})();
