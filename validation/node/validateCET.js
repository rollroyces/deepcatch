#!/usr/bin/env node
/**
 * validateCET.js - Longitudinal CET with Proper Threshold Calibration
 * IMPROVEMENT: 4 threshold methods + lambda parameter tuning
 */
const fs = require('fs');
const path = require('path');

const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'cet_results.json');
const SEED = 42;
const N_BOOTSTRAP = 2000;

// Seeded RNG
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

function poisson(lambda, rng) {
  if (lambda < 30) {
    const L = Math.exp(-lambda);
    let k = 0, p = 1;
    do { k++; p *= rng(); } while (p > L);
    return k - 1;
  }
  const x = normalRand(rng) * Math.sqrt(lambda) + lambda;
  return Math.max(0, Math.round(x));
}

// Parameters
const N_CANCER = 200;
const N_HEALTHY = 400;
const N_BENIGN = 100;
const N_QUARTERS = 8;
const SEQUENCING_DEPTH = 50000;
const ERROR_RATE = 0.0001;
const BASELINE_VAF = 0.000003;
const HEALTHY_VAF = 0.000001;
const DOUBLING_TIME_MIN = 150;
const DOUBLING_TIME_MAX = 300;

// Lambda values to test
const LAMBDA_VALUES = [0.001, 0.005, 0.01, 0.05, 0.1];

function simulatePatient(type, rng, patientIdx) {
  const patient = { id: `${type}_P${patientIdx}`, type, quarters: [] };
  let baseVAF;
  if (type === 'cancer') baseVAF = BASELINE_VAF * (0.5 + rng() * 1.5);
  else if (type === 'healthy') baseVAF = HEALTHY_VAF * (0.3 + rng() * 1.4);
  else baseVAF = HEALTHY_VAF * 2;
  const doublingTime = DOUBLING_TIME_MIN + rng() * (DOUBLING_TIME_MAX - DOUBLING_TIME_MIN);

  for (let q = 0; q < N_QUARTERS; q++) {
    const tDays = q * 90;
    let trueVAF;
    if (type === 'cancer') {
      const growthRate = Math.log(2) / doublingTime;
      trueVAF = baseVAF * Math.exp(growthRate * tDays);
      trueVAF *= (0.8 + rng() * 0.4);
    } else if (type === 'healthy') {
      trueVAF = baseVAF + normalRand(rng) * baseVAF * 0.3;
      trueVAF = Math.max(0, trueVAF);
    } else {
      const spikeFactor = Math.exp(-((q - 2.5) ** 2) / 3) * (5 + rng() * 15);
      trueVAF = baseVAF + spikeFactor * HEALTHY_VAF * 10;
      trueVAF = Math.max(0, trueVAF);
    }
    const trueReads = SEQUENCING_DEPTH * trueVAF;
    const mutantReads = poisson(Math.max(0.1, trueReads), rng);
    patient.quarters.push({
      quarter: q, days: tDays, true_vaf: trueVAF,
      mutant_reads: mutantReads, depth: SEQUENCING_DEPTH,
      observed_vaf: mutantReads / SEQUENCING_DEPTH
    });
  }
  return patient;
}

function computeCETScore(patient, lambdaCancerPerRead = 0.00001) {
  let cumulativeScore = 0;
  const timepoints = [];
  patient.quarters.forEach(q => {
    const lambda0 = ERROR_RATE;
    const lambda1 = ERROR_RATE + lambdaCancerPerRead;
    const logLR = q.mutant_reads * Math.log(lambda1 / lambda0) - q.depth * (lambda1 - lambda0);
    cumulativeScore += logLR;
    timepoints.push({
      quarter: q.quarter, log_lr: logLR, cumulative_score: cumulativeScore,
      mutant_reads: q.mutant_reads, observed_vaf: q.observed_vaf
    });
  });
  return { cumulative_score: cumulativeScore, timepoints };
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
  const npv = tn / Math.max(1, tn + fn);
  const f1 = (prec + sens > 0) ? 2 * prec * sens / (prec + sens) : 0;
  const f2 = (4 * prec + sens > 0) ? (5 * prec * sens) / (4 * prec + sens) : 0;
  return { tp, fp, tn, fn, sensitivity: sens, specificity: spec, precision: prec, npv, f1, f2 };
}

function computeAUC(scores, labels) {
  const pairs = scores.map((s, i) => ({ s, l: labels[i] }));
  const pos = pairs.filter(p => p.l === 1), neg = pairs.filter(p => p.l === 0);
  if (pos.length === 0 || neg.length === 0) return 0.5;
  let auc = 0, nPos = pos.length, nNeg = neg.length;
  for (const p of pos) for (const n of neg) { if (p.s > n.s) auc++; else if (p.s === n.s) auc += 0.5; }
  return auc / (nPos * nNeg);
}

function bootstrapCI(values, nBoot, rng) {
  const estimates = [];
  const n = values.length;
  for (let b = 0; b < nBoot; b++) {
    const sample = [];
    for (let i = 0; i < n; i++) sample.push(values[Math.floor(rng() * n)]);
    estimates.push(sample.reduce((a, b) => a + b, 0) / sample.length);
  }
  estimates.sort((a, b) => a - b);
  const mean = estimates.reduce((a, b) => a + b, 0) / estimates.length;
  const lo = estimates[Math.floor(0.025 * estimates.length)];
  const hi = estimates[Math.ceil(0.975 * estimates.length) - 1];
  return { mean, ci95_low: lo, ci95_high: hi };
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
  const mean = estimates.reduce((a, b) => a + b, 0) / estimates.length;
  const lo = estimates[Math.floor(0.025 * estimates.length)];
  const hi = estimates[Math.ceil(0.975 * estimates.length) - 1];
  return { mean, ci95_low: lo, ci95_high: hi };
}

// ---- Threshold Calibration Methods ----

// METHOD A: Youden's J Index
function youdenThreshold(scores, labels) {
  const sorted = scores.map((s, i) => ({ s, l: labels[i] })).sort((a, b) => a.s - b.s);
  const totalPos = labels.filter(l => l === 1).length;
  const totalNeg = labels.filter(l => l === 0).length;
  let bestJ = -1, bestThresh = 0;
  for (let i = 0; i < sorted.length; i++) {
    const thresh = sorted[i].s;
    const tp = sorted.filter(p => p.s >= thresh && p.l === 1).length;
    const tn = sorted.filter(p => p.s < thresh && p.l === 0).length;
    const sens = tp / Math.max(1, totalPos);
    const spec = tn / Math.max(1, totalNeg);
    const J = sens + spec - 1;
    if (J > bestJ) { bestJ = J; bestThresh = thresh; }
  }
  // Evaluate metrics at this threshold
  const preds = scores.map(s => s >= bestThresh ? 1 : 0);
  const metrics = computeMetrics(labels, preds);
  return { threshold: parseFloat(bestThresh.toFixed(6)), youden_j: parseFloat(bestJ.toFixed(4)), ...metrics };
}

// METHOD B: F2 Score (weights recall 2x over precision)
function f2Threshold(scores, labels) {
  const sorted = scores.map((s, i) => ({ s, l: labels[i] })).sort((a, b) => a.s - b.s);
  const totalPos = labels.filter(l => l === 1).length;
  const totalNeg = labels.filter(l => l === 0).length;
  let bestF2 = -1, bestThresh = 0;
  for (let i = 0; i < sorted.length; i++) {
    const thresh = sorted[i].s;
    const tp = sorted.filter(p => p.s >= thresh && p.l === 1).length;
    const fp = sorted.filter(p => p.s >= thresh && p.l === 0).length;
    const fn = totalPos - tp;
    const prec = tp / Math.max(1, tp + fp);
    const rec = tp / Math.max(1, tp + fn);
    const f2 = (4 * prec + rec > 0) ? (5 * prec * rec) / (4 * prec + rec) : 0;
    if (f2 > bestF2) { bestF2 = f2; bestThresh = thresh; }
  }
  const preds = scores.map(s => s >= bestThresh ? 1 : 0);
  const metrics = computeMetrics(labels, preds);
  return { threshold: parseFloat(bestThresh.toFixed(6)), f2_score: parseFloat(bestF2.toFixed(4)), ...metrics };
}

// METHOD C: Cost-Sensitive Threshold
function costSensitiveThreshold(scores, labels, costFP = 5000, costFN = 200000) {
  const sorted = scores.map((s, i) => ({ s, l: labels[i] })).sort((a, b) => a.s - b.s);
  const totalPos = labels.filter(l => l === 1).length;
  const totalNeg = labels.filter(l => l === 0).length;
  let bestCost = Infinity, bestThresh = 0;
  for (let i = 0; i < sorted.length; i++) {
    const thresh = sorted[i].s;
    const tp = sorted.filter(p => p.s >= thresh && p.l === 1).length;
    const fp = sorted.filter(p => p.s >= thresh && p.l === 0).length;
    const fn = totalPos - tp;
    const cost = fp * costFP + fn * costFN;
    if (cost < bestCost) { bestCost = cost; bestThresh = thresh; }
  }
  const preds = scores.map(s => s >= bestThresh ? 1 : 0);
  const metrics = computeMetrics(labels, preds);
  const expectedCostPer100K = (bestCost / labels.length) * 100000;
  return {
    threshold: parseFloat(bestThresh.toFixed(6)),
    cost_per_sample: parseFloat((bestCost / labels.length).toFixed(2)),
    expected_cost_per_100K: parseFloat(expectedCostPer100K.toFixed(2)),
    cost_FP: costFP,
    cost_FN: costFN,
    ...metrics
  };
}

// METHOD D: Multi-Tier Thresholds (percentile-based)
function multiTierThresholds(scores, labels) {
  const sorted = scores.map((s, i) => ({ s, l: labels[i] })).sort((a, b) => a.s - b.s);
  const percentiles = {
    P90: sorted[Math.floor(0.90 * sorted.length)],
    P97_5: sorted[Math.floor(0.975 * sorted.length)],
    P99_7: sorted[Math.floor(0.997 * sorted.length)],
    P99_9: sorted[Math.floor(0.999 * sorted.length)]
  };

  const tiers = {};
  for (const [name, p] of Object.entries(percentiles)) {
    const thresh = p.s;
    const preds = scores.map(s => s >= thresh ? 1 : 0);
    const metrics = computeMetrics(labels, preds);
    const totalPos = labels.filter(l => l === 1).length;
    const totalNeg = labels.filter(l => l === 0).length;
    const cost = metrics.fp * 5000 + metrics.fn * 200000;
    tiers[name] = {
      threshold: parseFloat(thresh.toFixed(6)),
      percentile: name === 'P90' ? 90 : name === 'P97_5' ? 97.5 : name === 'P99_7' ? 99.7 : 99.9,
      sensitivity: parseFloat(metrics.sensitivity.toFixed(4)),
      specificity: parseFloat(metrics.specificity.toFixed(4)),
      ppv: parseFloat(metrics.precision.toFixed(4)),
      npv: parseFloat(metrics.npv.toFixed(4)),
      f1: parseFloat(metrics.f1.toFixed(4)),
      f2: parseFloat(metrics.f2.toFixed(4)),
      expected_cost_per_100K: parseFloat(((cost / labels.length) * 100000).toFixed(2)),
      risk_level: name === 'P90' ? 'Low Risk (screening)' :
                  name === 'P97_5' ? 'Borderline (referral)' :
                  name === 'P99_7' ? 'Elevated (urgent workup)' : 'High Risk (immediate action)',
      recommended_action: name === 'P90' ? 'Repeat in 6 months' :
                          name === 'P97_5' ? 'Confirmatory testing' :
                          name === 'P99_7' ? 'Imaging + biopsy' : 'Emergency oncology referral'
    };
  }
  return tiers;
}

// ---- Lambda Tuning ----
function tuneLambda(patients, calibSet, testSet, lambdaValues, calibrationMethod, rng) {
  const results = [];
  for (const lambda of lambdaValues) {
    // Compute CET scores at final quarter for all patients with this lambda
    const calibScores = [];
    const calibLabels = [];
    const testScores = [];
    const testLabels = [];

    calibSet.forEach(p => {
      const cet = computeCETScore(p, lambda);
      calibScores.push(cet.timepoints[N_QUARTERS - 1].cumulative_score);
      calibLabels.push(p.type === 'cancer' ? 1 : 0);
    });

    testSet.forEach(p => {
      const cet = computeCETScore(p, lambda);
      testScores.push(cet.timepoints[N_QUARTERS - 1].cumulative_score);
      testLabels.push(p.type === 'cancer' ? 1 : 0);
    });

    // Find threshold using the specified method on calibration set
    let threshold;
    if (calibrationMethod === 'youden') {
      threshold = youdenThreshold(calibScores, calibLabels).threshold;
    } else if (calibrationMethod === 'f2') {
      threshold = f2Threshold(calibScores, calibLabels).threshold;
    } else if (calibrationMethod === 'cost') {
      threshold = costSensitiveThreshold(calibScores, calibLabels).threshold;
    } else {
      threshold = youdenThreshold(calibScores, calibLabels).threshold;
    }

    // Evaluate on test set
    const testPreds = testScores.map(s => s >= threshold ? 1 : 0);
    const metrics = computeMetrics(testLabels, testPreds);
    const auc = computeAUC(testScores, testLabels);
    const aucCI = bootstrapAUC(testScores, testLabels, N_BOOTSTRAP, rng);

    results.push({
      lambda: parseFloat(lambda.toFixed(4)),
      threshold: parseFloat(threshold.toFixed(6)),
      auc: parseFloat(aucCI.mean.toFixed(4)),
      auc_ci95_low: parseFloat(aucCI.ci95_low.toFixed(4)),
      auc_ci95_high: parseFloat(aucCI.ci95_high.toFixed(4)),
      sensitivity: parseFloat(metrics.sensitivity.toFixed(4)),
      specificity: parseFloat(metrics.specificity.toFixed(4)),
      ppv: parseFloat(metrics.precision.toFixed(4)),
      npv: parseFloat(metrics.npv.toFixed(4)),
      f1: parseFloat(metrics.f1.toFixed(4)),
      f2: parseFloat(metrics.f2.toFixed(4)),
      tp: metrics.tp, fp: metrics.fp, tn: metrics.tn, fn: metrics.fn
    });
  }

  // Best lambda by AUC
  const bestByAUC = results.reduce((b, r) => r.auc > b.auc ? r : b, results[0]);
  const bestByF2 = results.reduce((b, r) => r.f2 > b.f2 ? r : b, results[0]);

  return { all_results: results, best_by_auc: bestByAUC, best_by_f2: bestByF2 };
}

// ---- Main ----
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH NODE.JS - IMPROVED CET CALIBRATION');
  console.log('   4 Threshold Methods + λ Parameter Tuning');
  console.log('='.repeat(70));

  const rng = createRNG(SEED + 4000);
  console.log(`\n⚙️  Simulating ${N_CANCER + N_HEALTHY + N_BENIGN} patients...`);

  const allPatients = [];
  for (let i = 0; i < N_CANCER; i++) allPatients.push(simulatePatient('cancer', rng, i));
  for (let i = 0; i < N_HEALTHY; i++) allPatients.push(simulatePatient('healthy', rng, i));
  for (let i = 0; i < N_BENIGN; i++) allPatients.push(simulatePatient('benign', rng, i));

  // Shuffle
  const rngShuffle = createRNG(SEED + 4500);
  const shuffled = [...allPatients];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(rngShuffle() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }

  // Split: 30% calibration, 70% test
  const nCalib = Math.floor(shuffled.length * 0.3);
  const calibPatients = shuffled.slice(0, nCalib);
  const testPatients = shuffled.slice(nCalib);

  // Default lambda for calibration comparison
  const DEFAULT_LAMBDA = 0.01;

  // Collect CET scores at final quarter
  const calibScores = [], calibLabels = [];
  calibPatients.forEach(p => {
    const cet = computeCETScore(p, DEFAULT_LAMBDA);
    calibScores.push(cet.timepoints[N_QUARTERS - 1].cumulative_score);
    calibLabels.push(p.type === 'cancer' ? 1 : 0);
  });

  const testScores = [], testLabels = [];
  testPatients.forEach(p => {
    const cet = computeCETScore(p, DEFAULT_LAMBDA);
    testScores.push(cet.timepoints[N_QUARTERS - 1].cumulative_score);
    testLabels.push(p.type === 'cancer' ? 1 : 0);
  });

  const rngCI = createRNG(SEED + 5000);
  const testAUC = bootstrapAUC(testScores, testLabels, N_BOOTSTRAP, rngCI);
  console.log(`\n📊 Dataset: ${calibPatients.length} calibration, ${testPatients.length} test`);
  console.log(`   CET AUC (final quarter, λ=${DEFAULT_LAMBDA}): ${testAUC.mean.toFixed(4)} [${testAUC.ci95_low.toFixed(4)}–${testAUC.ci95_high.toFixed(4)}]`);

  // === METHOD A: Youden's J ===
  console.log('\n🔬 METHOD A: Youden\'s J Index');
  const youden = youdenThreshold(calibScores, calibLabels);
  const testPredsYouden = testScores.map(s => s >= youden.threshold ? 1 : 0);
  const youdenEval = computeMetrics(testLabels, testPredsYouden);
  console.log(`   Calibration threshold: ${youden.threshold.toFixed(6)}`);
  console.log(`   Test: Sens=${(youdenEval.sensitivity*100).toFixed(1)}% Spec=${(youdenEval.specificity*100).toFixed(1)}% PPV=${(youdenEval.precision*100).toFixed(1)}% F2=${youdenEval.f2.toFixed(4)}`);

  // === METHOD B: F2 Score ===
  console.log('\n🔬 METHOD B: F2 Score (recall-weighted)');
  const f2 = f2Threshold(calibScores, calibLabels);
  const testPredsF2 = testScores.map(s => s >= f2.threshold ? 1 : 0);
  const f2Eval = computeMetrics(testLabels, testPredsF2);
  console.log(`   Calibration threshold: ${f2.threshold.toFixed(6)}`);
  console.log(`   Test: Sens=${(f2Eval.sensitivity*100).toFixed(1)}% Spec=${(f2Eval.specificity*100).toFixed(1)}% PPV=${(f2Eval.precision*100).toFixed(1)}% F2=${f2Eval.f2.toFixed(4)}`);

  // === METHOD C: Cost-Sensitive ===
  console.log('\n🔬 METHOD C: Cost-Sensitive ($5K FP, $200K FN)');
  const cost = costSensitiveThreshold(calibScores, calibLabels);
  const testPredsCost = testScores.map(s => s >= cost.threshold ? 1 : 0);
  const costEval = computeMetrics(testLabels, testPredsCost);
  const testCost = (costEval.fp * 5000 + costEval.fn * 200000) / testLabels.length * 100000;
  console.log(`   Calibration threshold: ${cost.threshold.toFixed(6)}`);
  console.log(`   Calibration cost/sample: $${cost.cost_per_sample}`);
  console.log(`   Test: Sens=${(costEval.sensitivity*100).toFixed(1)}% Spec=${(costEval.specificity*100).toFixed(1)}% PPV=${(costEval.precision*100).toFixed(1)}% F2=${costEval.f2.toFixed(4)}`);
  console.log(`   Expected cost per 100K: $${testCost.toFixed(0)}`);

  // === METHOD D: Multi-Tier ===
  console.log('\n🔬 METHOD D: Multi-Tier Risk Stratification');
  const tiers = multiTierThresholds(calibScores, calibLabels);
  const tierResults = {};
  for (const [name, tier] of Object.entries(tiers)) {
    const testPreds = testScores.map(s => s >= tier.threshold ? 1 : 0);
    const testEval = computeMetrics(testLabels, testPreds);
    tierResults[name] = {
      ...tier,
      test_sensitivity: parseFloat(testEval.sensitivity.toFixed(4)),
      test_specificity: parseFloat(testEval.specificity.toFixed(4)),
      test_ppv: parseFloat(testEval.precision.toFixed(4))
    };
    console.log(`   ${name.padEnd(6)} (${tier.risk_level.padEnd(25)}): Thresh=${tier.threshold.toFixed(6)}, Test Sens=${(testEval.sensitivity*100).toFixed(1)}%, Spec=${(testEval.specificity*100).toFixed(1)}%, PPV=${(testEval.precision*100).toFixed(1)}%`);
  }

  // === Lambda Tuning ===
  console.log('\n🔧 LAMBDA PARAMETER TUNING:');
  console.log(`   λ values: ${LAMBDA_VALUES.join(', ')}`);

  const lambdaResults = {};
  for (const method of ['youden', 'f2', 'cost']) {
    console.log(`\n   --- Using ${method} calibration ---`);
    const tuning = tuneLambda(allPatients, calibPatients, testPatients, LAMBDA_VALUES, method, createRNG(SEED + 5100 + (method === 'youden' ? 0 : method === 'f2' ? 100 : 200)));
    lambdaResults[method] = tuning;

    console.log(`   ${'λ'.padEnd(10)} ${'AUC'.padEnd(10)} ${'Sens%'.padEnd(10)} ${'Spec%'.padEnd(10)} ${'F2'.padEnd(10)} ${'PPV%'.padEnd(10)}`);
    tuning.all_results.forEach(r => {
      console.log(`   ${String(r.lambda).padEnd(10)} ${r.auc.toFixed(4).padEnd(10)} ${(r.sensitivity*100).toFixed(1).padEnd(10)} ${(r.specificity*100).toFixed(1).padEnd(10)} ${r.f2.toFixed(4).padEnd(10)} ${(r.ppv*100).toFixed(1)}`);
    });
    console.log(`   🏆 Best by AUC: λ=${tuning.best_by_auc.lambda} (AUC=${tuning.best_by_auc.auc.toFixed(4)}, F2=${tuning.best_by_auc.f2.toFixed(4)})`);
    console.log(`   🏆 Best by F2:  λ=${tuning.best_by_f2.lambda} (F2=${tuning.best_by_f2.f2.toFixed(4)}, AUC=${tuning.best_by_f2.auc.toFixed(4)})`);
  }

  // === Quarterly tracking for the best lambda ===
  const bestLambda = lambdaResults['cost'].best_by_f2.lambda;
  console.log(`\n📈 Longitudinal tracking with λ=${bestLambda} (best by F2 with cost calibration):`);

  const quarterlyResults = {};
  for (let q = 0; q < N_QUARTERS; q++) {
    const qScores = [], qLabels = [];
    testPatients.forEach(p => {
      const cet = computeCETScore(p, bestLambda);
      qScores.push(cet.timepoints[q].cumulative_score);
      qLabels.push(p.type === 'cancer' ? 1 : 0);
    });
    const opt = lambdaResults['cost'];
    // Find threshold from calibration set
    const qCalibScores = [], qCalibLabels = [];
    calibPatients.forEach(p => {
      const cet = computeCETScore(p, bestLambda);
      qCalibScores.push(cet.timepoints[q].cumulative_score);
      qCalibLabels.push(p.type === 'cancer' ? 1 : 0);
    });
    const qThresh = costSensitiveThreshold(qCalibScores, qCalibLabels).threshold;
    const qPreds = qScores.map(s => s >= qThresh ? 1 : 0);
    const qMetrics = computeMetrics(qLabels, qPreds);
    const qAUC = computeAUC(qScores, qLabels);

    // Single timepoint comparison
    const singleScores = testPatients.map(p => p.quarters[q].observed_vaf);
    const singleAUC = computeAUC(singleScores, qLabels);

    quarterlyResults[q] = {
      quarter: q, months: q * 3,
      cet_auc: parseFloat(qAUC.toFixed(4)),
      single_auc: parseFloat(singleAUC.toFixed(4)),
      cet_sensitivity: parseFloat(qMetrics.sensitivity.toFixed(4)),
      cet_specificity: parseFloat(qMetrics.specificity.toFixed(4)),
      cet_ppv: parseFloat(qMetrics.precision.toFixed(4)),
      cet_f1: parseFloat(qMetrics.f1.toFixed(4)),
      cet_f2: parseFloat(qMetrics.f2.toFixed(4)),
      threshold: parseFloat(qThresh.toFixed(6))
    };
    console.log(`   Q${q+1} (${q*3}mo): CET AUC=${qAUC.toFixed(4)}, Single AUC=${singleAUC.toFixed(4)}, Sens=${(qMetrics.sensitivity*100).toFixed(1)}%, Spec=${(qMetrics.specificity*100).toFixed(1)}%, F2=${qMetrics.f2.toFixed(4)}`);
  }

  // === Output ===
  const output = {
    metadata: {
      validation_type: 'improved_longitudinal_cet',
      version: '2.0.0',
      n_cancer: N_CANCER, n_healthy: N_HEALTHY, n_benign: N_BENIGN,
      n_quarters: N_QUARTERS, total_duration_months: N_QUARTERS * 3,
      calibration_fraction: 0.3, test_fraction: 0.7,
      sequencing_depth: SEQUENCING_DEPTH, error_rate: ERROR_RATE,
      default_lambda: DEFAULT_LAMBDA, lambda_values_tested: LAMBDA_VALUES,
      n_calibration: calibPatients.length, n_test: testPatients.length,
      n_bootstrap: N_BOOTSTRAP,
      timestamp: new Date().toISOString()
    },
    baseline_cet_auc: {
      auc: parseFloat(testAUC.mean.toFixed(4)),
      ci95_low: parseFloat(testAUC.ci95_low.toFixed(4)),
      ci95_high: parseFloat(testAUC.ci95_high.toFixed(4))
    },
    calibration_methods: {
      youden_j: {
        calibration: { threshold: youden.threshold, youden_j: youden.youden_j, sensitivity: youden.sensitivity, specificity: youden.specificity },
        test: { sensitivity: parseFloat(youdenEval.sensitivity.toFixed(4)), specificity: parseFloat(youdenEval.specificity.toFixed(4)),
                ppv: parseFloat(youdenEval.precision.toFixed(4)), npv: parseFloat(youdenEval.npv.toFixed(4)),
                f1: parseFloat(youdenEval.f1.toFixed(4)), f2: parseFloat(youdenEval.f2.toFixed(4)) }
      },
      f2_score: {
        calibration: { threshold: f2.threshold, f2_score: f2.f2_score, sensitivity: f2.sensitivity, specificity: f2.specificity },
        test: { sensitivity: parseFloat(f2Eval.sensitivity.toFixed(4)), specificity: parseFloat(f2Eval.specificity.toFixed(4)),
                ppv: parseFloat(f2Eval.precision.toFixed(4)), npv: parseFloat(f2Eval.npv.toFixed(4)),
                f1: parseFloat(f2Eval.f1.toFixed(4)), f2: parseFloat(f2Eval.f2.toFixed(4)) }
      },
      cost_sensitive: {
        calibration: { threshold: cost.threshold, cost_per_sample: cost.cost_per_sample, 
                       sensitivity: cost.sensitivity, specificity: cost.specificity },
        test: { sensitivity: parseFloat(costEval.sensitivity.toFixed(4)), specificity: parseFloat(costEval.specificity.toFixed(4)),
                ppv: parseFloat(costEval.precision.toFixed(4)), npv: parseFloat(costEval.npv.toFixed(4)),
                f1: parseFloat(costEval.f1.toFixed(4)), f2: parseFloat(costEval.f2.toFixed(4)),
                expected_cost_per_100K: parseFloat(testCost.toFixed(0)) }
      },
      multi_tier: tierResults
    },
    lambda_tuning: {
      youden: lambdaResults['youden'],
      f2: lambdaResults['f2'],
      cost: lambdaResults['cost'],
      overall_best_lambda: parseFloat(bestLambda.toFixed(4)),
      overall_best_auc: parseFloat(lambdaResults['cost'].best_by_auc.auc.toFixed(4)),
      overall_best_f2: parseFloat(lambdaResults['cost'].best_by_f2.f2.toFixed(4))
    },
    quarterly_results_best_lambda: quarterlyResults,
    comparison_summary: {
      note: 'Comparison of calibration methods on test set performance',
      methods: [
        {
          method: 'Youden\'s J',
          sensitivity: parseFloat(youdenEval.sensitivity.toFixed(4)),
          specificity: parseFloat(youdenEval.specificity.toFixed(4)),
          f2: parseFloat(youdenEval.f2.toFixed(4))
        },
        {
          method: 'F2 Score (recall-weighted)',
          sensitivity: parseFloat(f2Eval.sensitivity.toFixed(4)),
          specificity: parseFloat(f2Eval.specificity.toFixed(4)),
          f2: parseFloat(f2Eval.f2.toFixed(4))
        },
        {
          method: 'Cost-Sensitive',
          sensitivity: parseFloat(costEval.sensitivity.toFixed(4)),
          specificity: parseFloat(costEval.specificity.toFixed(4)),
          f2: parseFloat(costEval.f2.toFixed(4))
        }
      ]
    }
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved to ${OUTPUT_PATH}`);
  console.log('\n✅ Improved CET calibration complete.');
  console.log('='.repeat(70));
})();
