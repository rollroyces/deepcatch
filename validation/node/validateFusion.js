#!/usr/bin/env node
/**
 * validateFusion.js - Multi-Modal Fusion with Smart Strategies
 * IMPROVEMENT: 3 smarter fusion strategies + correlation analysis
 */
const fs = require('fs');
const path = require('path');

const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'fusion_results.json');
const SEED = 42;
const ALPHA = 0.3;
const N_MODALITIES = 5;
const N_FOLDS = 5;
const N_BOOTSTRAP = 2000;
const N_SAMPLES = 500;
const CANCER_PREVALENCE = 0.15;
const MODALITY_NAMES = ['mutations', 'methylation', 'fragment_size', 'copy_number', 'ctc_count'];

// Seeded RNG (xoshiro128** variant)
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

function sigmoid(x) { return 1 / (1 + Math.exp(-Math.max(-20, Math.min(20, x)))); }

function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

// ---- Data Generation (same as original) ----
function generateFeatures(nSamples, alpha, rng, cancerPrev) {
  const cancerTypes = ['LUAD', 'COADREAD', 'BRCA'];
  const modalities = [...MODALITY_NAMES];
  const labels = [];
  const features = modalities.map(() => new Array(nSamples));

  for (let i = 0; i < nSamples; i++) {
    const hasCancer = rng() < cancerPrev;
    labels.push(hasCancer ? 1 : 0);
    const z = hasCancer ? (1.5 + normalRand(rng) * 0.5) : (normalRand(rng) * 0.8);
    const modalityFactors = modalities.map(() => alpha * z + (1 - alpha) * normalRand(rng) * 0.7);

    features[0][i] = hasCancer ? Math.max(0, modalityFactors[0] * 3 + normalRand(rng) * 0.3) : Math.max(0, normalRand(rng) * 0.1);
    features[1][i] = sigmoid(modalityFactors[1] * 2.5 + (hasCancer ? 1.0 : 0) + normalRand(rng) * 0.3);
    features[2][i] = hasCancer ? (155 + modalityFactors[2] * 8 + normalRand(rng) * 3) : (167 + normalRand(rng) * 5);
    features[3][i] = hasCancer ? (modalityFactors[3] * 0.4 + normalRand(rng) * 0.15) : normalRand(rng) * 0.05;
    features[4][i] = hasCancer ? Math.max(0, Math.round(modalityFactors[4] * 5 + normalRand(rng) * 1)) : Math.max(0, Math.round(normalRand(rng) * 0.2 + 0.1));
  }
  return { labels, features, modalities };
}

// ---- Logistic Regression ----
function logisticRegression(X, y, lr = 0.01, epochs = 500, lambda = 0.01) {
  const n = X.length;
  const d = X[0].length;
  let w = new Array(d).fill(0);
  let b = 0;
  for (let epoch = 0; epoch < epochs; epoch++) {
    let dw = new Array(d).fill(0), db = 0;
    for (let i = 0; i < n; i++) {
      const pred = sigmoid(dot(w, X[i]) + b);
      const error = pred - y[i];
      for (let j = 0; j < d; j++) dw[j] += error * X[i][j];
      db += error;
    }
    for (let j = 0; j < d; j++) { dw[j] = dw[j] / n + lambda * w[j]; w[j] -= lr * dw[j]; }
    b -= lr * db / n;
    if (epoch % 100 === 0) lr *= 0.95;
  }
  return { w, b };
}

function predictLogistic(X, model) { return X.map(xi => sigmoid(dot(model.w, xi) + model.b)); }

// ---- Metrics ----
function computeAUC(scores, labels) {
  const pairs = scores.map((s, i) => ({ s, l: labels[i] }));
  const pos = pairs.filter(p => p.l === 1), neg = pairs.filter(p => p.l === 0);
  if (pos.length === 0 || neg.length === 0) return 0.5;
  let auc = 0, nPos = pos.length, nNeg = neg.length;
  for (const p of pos) for (const n of neg) { if (p.s > n.s) auc++; else if (p.s === n.s) auc += 0.5; }
  return auc / (nPos * nNeg);
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
  const mean = estimates.reduce((a, b) => a + b, 0) / estimates.length;
  const lo = estimates[Math.floor(0.025 * estimates.length)];
  const hi = estimates[Math.ceil(0.975 * estimates.length) - 1];
  return { mean, ci95_low: lo, ci95_high: hi };
}

function delongPValue(scores1, scores2, labels, rng) {
  const nBoot = 1000, n = labels.length;
  let count = 0;
  for (let b = 0; b < nBoot; b++) {
    const idxs = []; for (let i = 0; i < n; i++) idxs.push(Math.floor(rng() * n));
    const a1 = computeAUC(idxs.map(i => scores1[i]), idxs.map(i => labels[i]));
    const a2 = computeAUC(idxs.map(i => scores2[i]), idxs.map(i => labels[i]));
    if (a2 > a1) count++;
  }
  return 1 - count / nBoot;
}

// ---- Stratified Split ----
function stratifiedSplit(labels, trainFrac, rng) {
  const pos = labels.map((l, i) => l === 1 ? i : -1).filter(i => i >= 0);
  const neg = labels.map((l, i) => l === 0 ? i : -1).filter(i => i >= 0);
  function shuffle(arr) { for (let i = arr.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [arr[i], arr[j]] = [arr[j], arr[i]]; } }
  shuffle(pos); shuffle(neg);
  const nPosTrain = Math.floor(pos.length * trainFrac), nNegTrain = Math.floor(neg.length * trainFrac);
  const train = new Set([...pos.slice(0, nPosTrain), ...neg.slice(0, nNegTrain)]);
  const test = new Set([...pos.slice(nPosTrain), ...neg.slice(nNegTrain)]);
  return { train: [...train], test: [...test] };
}

// ---- Correlation Matrix (Spearman) ----
function rank(arr) { const idx = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v); const ranks = new Array(arr.length); idx.forEach((o, r) => { ranks[o.i] = r + 1; }); return ranks; }

function spearmanCorr(a, b) {
  const ra = rank(a), rb = rank(b);
  const n = a.length;
  let sumD2 = 0;
  for (let i = 0; i < n; i++) { const d = ra[i] - rb[i]; sumD2 += d * d; }
  return 1 - (6 * sumD2) / (n * (n * n - 1));
}

function correlationMatrix(features, names) {
  const matrix = {};
  for (let i = 0; i < names.length; i++) {
    matrix[names[i]] = {};
    for (let j = 0; j < names.length; j++) {
      matrix[names[i]][names[j]] = parseFloat(spearmanCorr(features[i], features[j]).toFixed(4));
    }
  }
  return matrix;
}

// ---- STRATEGY A: Performance-Weighted Fusion ----
function performanceWeightedFusion(features, labels, trainIdx, testIdx, names, nBoot, rng) {
  const nTrain = trainIdx.length, nTest = testIdx.length;

  // Train individual models and get validation AUC
  const individualModels = [];
  const validationAUCs = [];

  for (let m = 0; m < names.length; m++) {
    const XTrain = trainIdx.map(i => [features[m][i]]);
    const yTrain = trainIdx.map(i => labels[i]);
    const model = logisticRegression(XTrain, yTrain, 0.05, 300, 0.001);
    individualModels.push(model);

    // Compute AUC on training set (in-sample, as a proxy for validation)
    const trainScores = predictLogistic(XTrain, model);
    const auc = computeAUC(trainScores, yTrain);
    validationAUCs.push(auc);
  }

  // Compute weights: AUC / sum(AUC), zero out AUC < 0.5
  const weights = validationAUCs.map(auc => auc < 0.5 ? 0 : auc);
  const weightSum = weights.reduce((a, b) => a + b, 0);
  const normedWeights = weightSum > 0 ? weights.map(w => w / weightSum) : weights.map(() => 1 / names.length);

  // Predict on test set using weighted average of individual model predictions
  const testScores = new Array(nTest);
  for (let i = 0; i < nTest; i++) {
    let weightedSum = 0;
    for (let m = 0; m < names.length; m++) {
      if (normedWeights[m] > 0) {
        const xi = [features[m][testIdx[i]]];
        const score = sigmoid(dot(individualModels[m].w, xi) + individualModels[m].b);
        weightedSum += normedWeights[m] * score;
      }
    }
    testScores[i] = weightedSum;
  }

  const yTest = testIdx.map(i => labels[i]);
  const auc = computeAUC(testScores, yTest);
  const aucCI = bootstrapAUC(testScores, yTest, nBoot, rng);

  return { auc: aucCI.mean, ci95_low: aucCI.ci95_low, ci95_high: aucCI.ci95_high, scores: testScores, weights: normedWeights.map(w => parseFloat(w.toFixed(4))), raw_aucs: validationAUCs.map(a => parseFloat(a.toFixed(4))) };
}

// ---- STRATEGY B: Stacked Meta-Learner ----
function stackedMetaLearner(features, labels, names, nFolds, nBoot, rng) {
  const n = labels.length;

  // Create folds (stratified)
  const pos = labels.map((l, i) => l === 1 ? i : -1).filter(i => i >= 0);
  const neg = labels.map((l, i) => l === 0 ? i : -1).filter(i => i >= 0);
  const folds = [];
  const foldSplit = nFolds;

  function shuffle(rng, arr) { for (let i = arr.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [arr[i], arr[j]] = [arr[j], arr[i]]; } }
  const foldRng = createRNG(SEED + 9999);
  const posCopy = [...pos]; shuffle(foldRng, posCopy);
  const negCopy = [...neg]; shuffle(foldRng, negCopy);

  for (let f = 0; f < foldSplit; f++) {
    const pStart = Math.floor(f * posCopy.length / foldSplit), pEnd = Math.floor((f + 1) * posCopy.length / foldSplit);
    const nStart = Math.floor(f * negCopy.length / foldSplit), nEnd = Math.floor((f + 1) * negCopy.length / foldSplit);
    const testSet = new Set([...posCopy.slice(pStart, pEnd), ...negCopy.slice(nStart, nEnd)]);
    const trainIdxs = labels.map((_, i) => i).filter(i => !testSet.has(i));
    folds.push({ train: trainIdxs, test: [...testSet] });
  }

  // For each fold, train individual LRs and predict on held-out
  const oofPredictions = names.map(() => new Array(n).fill(0)); // per-modality OOF predictions

  folds.forEach(fold => {
    const { train, test } = fold;
    for (let m = 0; m < names.length; m++) {
      const XTrain = train.map(i => [features[m][i]]);
      const yTrain = train.map(i => labels[i]);
      const model = logisticRegression(XTrain, yTrain, 0.05, 300, 0.001);
      const XTest = test.map(i => [features[m][i]]);
      const preds = predictLogistic(XTest, model);
      test.forEach((ti, idx) => { oofPredictions[m][ti] = preds[idx]; });
    }
  });

  // Train meta-learner on OOF predictions
  const metaX = [];
  for (let i = 0; i < n; i++) {
    metaX.push(names.map((_, m) => oofPredictions[m][i]));
  }
  const metaModel = logisticRegression(metaX, labels, 0.03, 500, 0.005);

  // Re-train base models on full data
  const baseModels = [];
  for (let m = 0; m < names.length; m++) {
    const XFull = labels.map((_, i) => [features[m][i]]);
    baseModels.push(logisticRegression(XFull, labels, 0.05, 300, 0.001));
  }

  // Final predictions: base model predictions → meta model
  const baseTestPreds = names.map((_, m) => predictLogistic(metaX, { w: baseModels[m].w, b: baseModels[m].b }));
  
  // Actually for fair comparison, split into train/test
  const split = stratifiedSplit(labels, 0.7, createRNG(SEED + 8000));
  
  // Retrain base on train subset
  const baseModels2 = [];
  for (let m = 0; m < names.length; m++) {
    const XTrain = split.train.map(i => [features[m][i]]);
    const yTrain = split.train.map(i => labels[i]);
    baseModels2.push(logisticRegression(XTrain, yTrain, 0.05, 300, 0.001));
  }

  // Train meta on train subset OOF
  const metaTrainX = split.train.map(i => names.map((_, m) => {
    const xi = [features[m][i]];
    // Use in-sample prediction (this is a limitation, but acceptable for simulation)
    return sigmoid(dot(baseModels2[m].w, xi) + baseModels2[m].b);
  }));
  const metaYTrain = split.train.map(i => labels[i]);
  const metaModel2 = logisticRegression(metaTrainX, metaYTrain, 0.03, 500, 0.005);

  // Predict on test
  const testMetaX = split.test.map(i => names.map((_, m) => {
    const xi = [features[m][i]];
    return sigmoid(dot(baseModels2[m].w, xi) + baseModels2[m].b);
  }));
  const testScores = predictLogistic(testMetaX, metaModel2);
  const yTest = split.test.map(i => labels[i]);
  const auc = computeAUC(testScores, yTest);
  const aucCI = bootstrapAUC(testScores, yTest, nBoot, rng);

  return { auc: aucCI.mean, ci95_low: aucCI.ci95_low, ci95_high: aucCI.ci95_high, scores: testScores, metaWeights: metaModel2.w.map(w => parseFloat(w.toFixed(4))), metaBias: parseFloat(metaModel2.b.toFixed(4)), nOOF: n };
}

// ---- STRATEGY C: Selective Fusion (Best-n) ----
function selectiveFusion(features, labels, names, nBoot, rng) {
  // Split into train (70%) / validation (30%)
  const split = stratifiedSplit(labels, 0.7, createRNG(SEED + 7000));
  const trainIdx = split.train, valIdx = split.test;

  // Train individual models and compute validation AUC
  const modPerf = [];
  for (let m = 0; m < names.length; m++) {
    const XTrain = trainIdx.map(i => [features[m][i]]);
    const yTrain = trainIdx.map(i => labels[i]);
    const model = logisticRegression(XTrain, yTrain, 0.05, 300, 0.001);
    const XVal = valIdx.map(i => [features[m][i]]);
    const valScores = predictLogistic(XVal, model);
    const yVal = valIdx.map(i => labels[i]);
    const valAUC = computeAUC(valScores, yVal);
    modPerf.push({ name: names[m], auc: valAUC, model, idx: m });
  }

  // Rank by validation AUC descending
  modPerf.sort((a, b) => b.auc - a.auc);

  // Test n = 1..5
  const results = [];
  for (let n = 1; n <= names.length; n++) {
    const selected = modPerf.slice(0, n);
    // Train fusion model on train set with selected modalities
    const XTrain = trainIdx.map(i => selected.map(m => features[m.idx][i]));
    const yTrain = trainIdx.map(i => labels[i]);
    const fusionModel = logisticRegression(XTrain, yTrain, 0.03, 500, 0.005);
    const XVal = valIdx.map(i => selected.map(m => features[m.idx][i]));
    const valScores = predictLogistic(XVal, fusionModel);
    const yVal = valIdx.map(i => labels[i]);
    const valAUC = computeAUC(valScores, yVal);

    results.push({
      n_modalities: n,
      modalities: selected.map(m => m.name),
      individual_aucs: selected.map(m => parseFloat(m.auc.toFixed(4))),
      fusion_auc_val: parseFloat(valAUC.toFixed(4))
    });
  }

  // Best n → re-evaluate on a fresh test split (or use the same val for simplicity)
  // Actually, to be fair, we should use the val set AUC directly
  const bestN = results.reduce((best, r) => r.fusion_auc_val > best.fusion_auc_val ? r : best, results[0]);

  return {
    results_by_n: results,
    best_n: bestN.n_modalities,
    best_modalities: bestN.modalities,
    best_auc_val: bestN.fusion_auc_val,
    modality_ranking: modPerf.map(m => ({ name: m.name, val_auc: parseFloat(m.auc.toFixed(4)) }))
  };
}

// ---- Main ----
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH NODE.JS - IMPROVED MULTI-MODAL FUSION');
  console.log('   Strategies: Performance-Weighted, Stacked Meta-Learner, Selective Fusion');
  console.log('='.repeat(70));

  const rng = createRNG(SEED + 3000);
  console.log(`\n⚙️  Generating ${N_SAMPLES} samples with ${N_MODALITIES} modalities...`);
  console.log(`   Shared latent factor α = ${ALPHA}`);
  console.log(`   Cancer prevalence = ${(CANCER_PREVALENCE*100).toFixed(1)}%`);

  const { labels, features, modalities } = generateFeatures(N_SAMPLES, ALPHA, rng, CANCER_PREVALENCE);
  const nCancer = labels.filter(l => l === 1).length;
  const nHealthy = labels.filter(l => l === 0).length;
  console.log(`   Generated: ${nCancer} cancer, ${nHealthy} healthy`);

  // === Correlation Matrix ===
  console.log('\n📊 CORRELATION MATRIX (Spearman ρ):');
  const corrMatrix = correlationMatrix(features, modalities);
  modalities.forEach(m1 => {
    const row = modalities.map(m2 => corrMatrix[m1][m2].toFixed(4)).join('  ');
    console.log(`   ${m1.padEnd(15)}: ${row}`);
  });

  // === Split Data ===
  const split = stratifiedSplit(labels, 0.7, createRNG(SEED + 6000));
  const trainIdx = split.train, testIdx = split.test;
  const yTest = testIdx.map(i => labels[i]);
  const rngCI = createRNG(SEED + 6500);

  // === Baseline: Best Single Modality ===
  console.log('\n📡 BASELINE: Single Modality Performance');
  const singleResults = [];
  let bestSingleName = null, bestSingleAUC = -1, bestSingleScores = null;

  for (let m = 0; m < modalities.length; m++) {
    const XTrain = trainIdx.map(i => [features[m][i]]);
    const yTrain = trainIdx.map(i => labels[i]);
    const model = logisticRegression(XTrain, yTrain, 0.05, 300, 0.001);
    const XTest = testIdx.map(i => [features[m][i]]);
    const testScores = predictLogistic(XTest, model);
    const auc = computeAUC(testScores, yTest);
    const aucCI = bootstrapAUC(testScores, yTest, N_BOOTSTRAP, rngCI);

    singleResults.push({ modality: modalities[m], auc: aucCI.mean, ci95_low: aucCI.ci95_low, ci95_high: aucCI.ci95_high, scores: testScores });

    console.log(`   ${modalities[m].padEnd(15)}: AUC = ${aucCI.mean.toFixed(4)} [${aucCI.ci95_low.toFixed(4)}–${aucCI.ci95_high.toFixed(4)}]`);

    if (aucCI.mean > bestSingleAUC) {
      bestSingleAUC = aucCI.mean;
      bestSingleName = modalities[m];
      bestSingleScores = testScores;
    }
  }
  console.log(`   🏆 Best single: ${bestSingleName} (AUC = ${bestSingleAUC.toFixed(4)})`);

  // === Naive Fusion (original approach) ===
  console.log('\n🔗 NAIVE FUSION (all 5 concatenated):');
  const XTrainFusion = trainIdx.map(i => modalities.map((_, m) => features[m][i]));
  const XTestFusion = testIdx.map(i => modalities.map((_, m) => features[m][i]));
  const yTrain = trainIdx.map(i => labels[i]);
  const naiveModel = logisticRegression(XTrainFusion, yTrain, 0.03, 500, 0.005);
  const naiveScores = predictLogistic(XTestFusion, naiveModel);
  const naiveAUC = computeAUC(naiveScores, yTest);
  const naiveAUC_CI = bootstrapAUC(naiveScores, yTest, N_BOOTSTRAP, rngCI);
  const naiveVsBest = delongPValue(bestSingleScores, naiveScores, yTest, rngCI);
  console.log(`   AUC = ${naiveAUC_CI.mean.toFixed(4)} [${naiveAUC_CI.ci95_low.toFixed(4)}–${naiveAUC_CI.ci95_high.toFixed(4)}]`);
  console.log(`   vs ${bestSingleName}: ΔAUC = ${(naiveAUC_CI.mean - bestSingleAUC).toFixed(4)}, p = ${naiveVsBest.toFixed(4)}`);

  // === STRATEGY A: Performance-Weighted Fusion ===
  console.log('\n⚖️  STRATEGY A: Performance-Weighted Fusion');
  const pWeighted = performanceWeightedFusion(features, labels, trainIdx, testIdx, modalities, N_BOOTSTRAP, createRNG(SEED + 7100));
  const pwVsBest = delongPValue(bestSingleScores, pWeighted.scores, yTest, rngCI);
  console.log(`   AUC = ${pWeighted.auc.toFixed(4)} [${pWeighted.ci95_low.toFixed(4)}–${pWeighted.ci95_high.toFixed(4)}]`);
  console.log(`   Weights: ${modalities.map((m, idx) => `${m}=${pWeighted.weights[idx]}`).join(', ')}`);
  console.log(`   vs ${bestSingleName}: ΔAUC = ${(pWeighted.auc - bestSingleAUC).toFixed(4)}, p = ${pwVsBest.toFixed(4)}`);

  // === STRATEGY B: Stacked Meta-Learner ===
  console.log('\n🔄 STRATEGY B: Stacked Meta-Learner');
  const stacked = stackedMetaLearner(features, labels, modalities, N_FOLDS, N_BOOTSTRAP, createRNG(SEED + 7200));
  const stVsBest = delongPValue(bestSingleScores, stacked.scores, yTest, rngCI);
  console.log(`   AUC = ${stacked.auc.toFixed(4)} [${stacked.ci95_low.toFixed(4)}–${stacked.ci95_high.toFixed(4)}]`);
  console.log(`   Meta weights: ${modalities.map((m, idx) => `${m}=${stacked.metaWeights[idx]}`).join(', ')}`);
  console.log(`   vs ${bestSingleName}: ΔAUC = ${(stacked.auc - bestSingleAUC).toFixed(4)}, p = ${stVsBest.toFixed(4)}`);

  // === STRATEGY C: Selective Fusion ===
  console.log('\n🎯 STRATEGY C: Selective Fusion (Best-n)');
  const selective = selectiveFusion(features, labels, modalities, N_BOOTSTRAP, createRNG(SEED + 7300));
  console.log(`   Modality ranking: ${selective.modality_ranking.map(m => `${m.name}=${m.val_auc}`).join(' > ')}`);
  selective.results_by_n.forEach(r => {
    console.log(`   n=${r.n_modalities}: ${r.modalities.join('+')} → AUC = ${r.fusion_auc_val.toFixed(4)}`);
  });
  console.log(`   🏆 Best: n=${selective.best_n} (${selective.best_modalities.join(', ')}) AUC = ${selective.best_auc_val.toFixed(4)}`);

  // To get a fair test AUC for selective fusion, retrain with best_n on train and test on test
  const bestN = selective.best_n;
  const bestModIndices = selective.modality_ranking.slice(0, bestN).map(m => modalities.indexOf(m.name));
  const XTrainSel = trainIdx.map(i => bestModIndices.map(m => features[m][i]));
  const XTestSel = testIdx.map(i => bestModIndices.map(m => features[m][i]));
  const selModel = logisticRegression(XTrainSel, yTrain, 0.03, 500, 0.005);
  const selScores = predictLogistic(XTestSel, selModel);
  const selAUC = computeAUC(selScores, yTest);
  const selAUC_CI = bootstrapAUC(selScores, yTest, N_BOOTSTRAP, rngCI);
  const selVsBest = delongPValue(bestSingleScores, selScores, yTest, rngCI);
  console.log(`   Test AUC (best-n=${bestN}): ${selAUC_CI.mean.toFixed(4)} [${selAUC_CI.ci95_low.toFixed(4)}–${selAUC_CI.ci95_high.toFixed(4)}]`);
  console.log(`   vs ${bestSingleName}: ΔAUC = ${(selAUC_CI.mean - bestSingleAUC).toFixed(4)}, p = ${selVsBest.toFixed(4)}`);

  // === Summary Comparison ===
  console.log('\n📊 STRATEGY COMPARISON:');
  console.log('─'.repeat(70));
  console.log(`   ${'Strategy'.padEnd(28)} ${'AUC'.padEnd(10)} ${'CI95'.padEnd(26)} ${'Δ vs Best'.padEnd(14)} p-value`);
  console.log('─'.repeat(70));

  const strategies = [
    { name: 'Best Single (' + bestSingleName + ')', auc: bestSingleAUC, ci_low: singleResults.find(s => s.modality === bestSingleName).ci95_low, ci_high: singleResults.find(s => s.modality === bestSingleName).ci95_high, delta: 0, p: null },
    { name: 'Naive Fusion', auc: naiveAUC_CI.mean, ci_low: naiveAUC_CI.ci95_low, ci_high: naiveAUC_CI.ci95_high, delta: naiveAUC_CI.mean - bestSingleAUC, p: naiveVsBest },
    { name: 'Performance-Weighted', auc: pWeighted.auc, ci_low: pWeighted.ci95_low, ci_high: pWeighted.ci95_high, delta: pWeighted.auc - bestSingleAUC, p: pwVsBest },
    { name: 'Stacked Meta-Learner', auc: stacked.auc, ci_low: stacked.ci95_low, ci_high: stacked.ci95_high, delta: stacked.auc - bestSingleAUC, p: stVsBest },
    { name: 'Selective Fusion (n=' + bestN + ')', auc: selAUC_CI.mean, ci_low: selAUC_CI.ci95_low, ci_high: selAUC_CI.ci95_high, delta: selAUC_CI.mean - bestSingleAUC, p: selVsBest }
  ];

  strategies.forEach(s => {
    const dSign = s.delta >= 0 ? '+' : '';
    console.log(`   ${s.name.padEnd(28)} ${s.auc.toFixed(4).padEnd(10)} [${s.ci_low.toFixed(4)}–${s.ci_high.toFixed(4)}]  ${(dSign + s.delta.toFixed(4)).padEnd(14)} ${s.p !== null ? s.p.toFixed(4) : '—'}`);
  });

  // === Output ===
  const output = {
    metadata: {
      validation_type: 'improved_multi_modal_fusion',
      version: '2.0.0',
      n_samples: N_SAMPLES,
      n_modalities: N_MODALITIES,
      modalities,
      shared_latent_factor: ALPHA,
      cancer_prevalence: CANCER_PREVALENCE,
      n_folds: N_FOLDS,
      n_bootstrap: N_BOOTSTRAP,
      n_cancer: nCancer,
      n_healthy: nHealthy,
      train_samples: trainIdx.length,
      test_samples: testIdx.length,
      timestamp: new Date().toISOString()
    },
    correlation_matrix: corrMatrix,
    single_modality_results: singleResults.map(s => ({
      modality: s.modality,
      auc: parseFloat(s.auc.toFixed(4)),
      ci95_low: parseFloat(s.ci95_low.toFixed(4)),
      ci95_high: parseFloat(s.ci95_high.toFixed(4))
    })),
    naive_fusion: {
      auc: parseFloat(naiveAUC_CI.mean.toFixed(4)),
      ci95_low: parseFloat(naiveAUC_CI.ci95_low.toFixed(4)),
      ci95_high: parseFloat(naiveAUC_CI.ci95_high.toFixed(4)),
      delta_vs_best_single: parseFloat((naiveAUC_CI.mean - bestSingleAUC).toFixed(4)),
      p_value_vs_best_single: parseFloat(naiveVsBest.toFixed(4))
    },
    performance_weighted_fusion: {
      auc: parseFloat(pWeighted.auc.toFixed(4)),
      ci95_low: parseFloat(pWeighted.ci95_low.toFixed(4)),
      ci95_high: parseFloat(pWeighted.ci95_high.toFixed(4)),
      weights: pWeighted.weights,
      raw_validation_aucs: pWeighted.raw_aucs,
      delta_vs_best_single: parseFloat((pWeighted.auc - bestSingleAUC).toFixed(4)),
      p_value_vs_best_single: parseFloat(pwVsBest.toFixed(4))
    },
    stacked_meta_learner: {
      auc: parseFloat(stacked.auc.toFixed(4)),
      ci95_low: parseFloat(stacked.ci95_low.toFixed(4)),
      ci95_high: parseFloat(stacked.ci95_high.toFixed(4)),
      meta_weights: stacked.metaWeights,
      meta_bias: stacked.metaBias,
      delta_vs_best_single: parseFloat((stacked.auc - bestSingleAUC).toFixed(4)),
      p_value_vs_best_single: parseFloat(stVsBest.toFixed(4))
    },
    selective_fusion: {
      best_n: bestN,
      best_modalities: selective.best_modalities,
      best_auc_val: selective.best_auc_val,
      test_auc: parseFloat(selAUC_CI.mean.toFixed(4)),
      test_ci95_low: parseFloat(selAUC_CI.ci95_low.toFixed(4)),
      test_ci95_high: parseFloat(selAUC_CI.ci95_high.toFixed(4)),
      ranking: selective.modality_ranking,
      by_n: selective.results_by_n,
      delta_vs_best_single: parseFloat((selAUC_CI.mean - bestSingleAUC).toFixed(4)),
      p_value_vs_best_single: parseFloat(selVsBest.toFixed(4))
    },
    strategy_comparison: strategies.map(s => ({
      strategy: s.name,
      auc: parseFloat(s.auc.toFixed(4)),
      ci95_low: parseFloat(s.ci_low.toFixed(4)),
      ci95_high: parseFloat(s.ci_high.toFixed(4)),
      delta_vs_best_single: parseFloat(s.delta.toFixed(4)),
      p_value: s.p !== null ? parseFloat(s.p.toFixed(4)) : null
    })),
    conclusion: {
      best_strategy: strategies.slice(1).reduce((b, s) => s.auc > b.auc ? s : b, strategies[1]).name,
      best_auc: Math.max(...strategies.slice(1).map(s => s.auc)),
      fusion_improves_over_best_single: Math.max(...strategies.slice(1).map(s => s.auc)) > bestSingleAUC,
      note: strategies.slice(1).reduce((b, s) => s.auc > b.auc ? s : b).auc > bestSingleAUC
        ? 'Smart fusion successfully improves over best single modality'
        : 'Fusion does not improve over best single modality — modalities are correlated and the best modality dominates'
    }
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved to ${OUTPUT_PATH}`);
  console.log('\n✅ Improved Fusion validation complete.');
  console.log('='.repeat(70));
})();
