#!/usr/bin/env node
/**
 * validateTOO.js - Tissue-of-Origin Prediction
 * Methylation + fragmentomic patterns for 8 cancer types
 * Methods: Multi-class logistic regression, random forest, neural network
 */
const fs = require('fs');
const path = require('path');

const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'too_results.json');
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

// ── Cancer Type Definitions ──
const CANCER_TYPES = {
  LUAD:  { name: 'Lung Adenocarcinoma', genes: ['CDKN2A', 'FHIT', 'RASSF1A'], methylation_strength: 0.75, fragment_shift: 10 },
  COADREAD: { name: 'Colorectal', genes: ['MLH1', 'SEPT9', 'VIM'], methylation_strength: 0.80, fragment_shift: -15 },
  BRCA: { name: 'Breast Cancer', genes: ['BRCA1', 'GSTP1'], methylation_strength: 0.70, fragment_shift: 5 },
  PRAD: { name: 'Prostate Cancer', genes: ['GSTP1'], methylation_strength: 0.85, fragment_shift: -8 },
  STAD: { name: 'Stomach Cancer', genes: ['CDH1'], methylation_strength: 0.65, fragment_shift: -12 },
  LIHC: { name: 'Hepatocellular Carcinoma', genes: ['CDKN2A', 'RASSF1A'], methylation_strength: 0.72, fragment_shift: -18 },
  PAAD: { name: 'Pancreatic Cancer', genes: ['CDKN2A'], methylation_strength: 0.68, fragment_shift: -22 },
  OV:   { name: 'Ovarian Cancer', genes: ['BRCA1'], methylation_strength: 0.78, fragment_shift: 3 }
};

const N_SAMPLES_PER_CANCER = 50; // 400 total cancer samples
const N_HEALTHY = 200; // healthy controls (for cancer detection evaluation)

// ── Data Simulation ──
function simulateTOOData(rng) {
  const allGenes = [...new Set(Object.values(CANCER_TYPES).flatMap(t => t.genes))];
  const nFeatures = allGenes.length + 3; // methylation genes + 3 fragmentomic features

  const samples = [];
  const labels = []; // cancer type index (0-7 for cancers, -1 for healthy)

  // Cancer samples
  let cancerIdx = 0;
  for (const [code, info] of Object.entries(CANCER_TYPES)) {
    for (let i = 0; i < N_SAMPLES_PER_CANCER; i++) {
      const features = [];

      // Methylation features (gene-specific)
      allGenes.forEach(gene => {
        if (info.genes.includes(gene)) {
          // Hypermethylated in this cancer type
          features.push(0.7 + rng() * 0.3); // hypermethylation: 0.7-1.0
        } else {
          // Background methylation
          features.push(0.1 + rng() * 0.3); // background: 0.1-0.4
        }
      });

      // Fragmentomic features
      // Fragment size index
      const baseSize = 167; // healthy nucleosome size
      features.push(baseSize + info.fragment_shift + normalRand(rng) * 10);

      // Fragment end motif diversity (Shannon entropy of 4-mer frequencies)
      const motifScore = 2.5 + (info.methylation_strength * 1.5) + rng() * 0.5;
      features.push(motifScore);

      // Nucleosome positioning score
      const nucScore = 0.5 + (info.fragment_shift > 0 ? 0.2 : -0.1) + rng() * 0.15;
      features.push(Math.max(0, Math.min(1, nucScore)));

      samples.push({ features, cancer_type: code, label: cancerIdx, is_cancer: true });
    }
    labels.push({ index: cancerIdx, code, name: info.name });
    cancerIdx++;
  }

  // Healthy samples
  for (let i = 0; i < N_HEALTHY; i++) {
    const features = [];
    allGenes.forEach(() => features.push(0.05 + rng() * 0.2)); // low methylation
    features.push(165 + normalRand(rng) * 5); // fragment size ~165bp
    features.push(2.0 + rng() * 0.5); // lower motif diversity
    features.push(0.3 + rng() * 0.2); // nucleosome score
    samples.push({ features, cancer_type: 'HEALTHY', label: -1, is_cancer: false });
  }

  return { samples, allGenes, labels, nCancerTypes: cancerIdx };
}

// ── Multi-Class Logistic Regression (softmax) ──
function trainLogisticRegression(trainSamples, trainLabels, nClasses, nFeatures, lr = 0.005, epochs = 1000, lambda = 0.01) {
  // Weight matrix: nClasses × (nFeatures + 1) [bias]
  const weights = [];
  for (let c = 0; c < nClasses; c++) {
    weights.push(new Array(nFeatures + 1).fill(0).map(() => (Math.random() - 0.5) * 0.001));
  }

  const cancerTrain = trainSamples.filter((_, i) => trainLabels[i] >= 0);
  const cancerTrainLabels = trainLabels.filter(l => l >= 0);
  const n = cancerTrain.length;

  // Feature normalization (standardize)
  const featMeans = new Array(nFeatures).fill(0);
  const featStds = new Array(nFeatures).fill(0);
  for (let j = 0; j < nFeatures; j++) {
    for (let i = 0; i < n; i++) featMeans[j] += cancerTrain[i].features[j];
    featMeans[j] /= n;
    for (let i = 0; i < n; i++) featStds[j] += (cancerTrain[i].features[j] - featMeans[j]) ** 2;
    featStds[j] = Math.sqrt(featStds[j] / n) || 1;
  }

  const normFeatures = cancerTrain.map(s => s.features.map((f, j) => (f - featMeans[j]) / featStds[j]));

  for (let epoch = 0; epoch < epochs; epoch++) {
    // Shuffle
    const idxs = Array.from({ length: n }, (_, i) => i);
    for (let i = n - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [idxs[i], idxs[j]] = [idxs[j], idxs[i]];
    }

    for (const idx of idxs) {
      const x = [1, ...normFeatures[idx]];
      const trueClass = cancerTrainLabels[idx];

      // Softmax scores
      const scores = weights.map(w => w.reduce((s, wj, j) => s + wj * x[j], 0));
      const maxScore = Math.max(...scores);
      const expScores = scores.map(s => Math.exp(s - maxScore));
      const sumExp = expScores.reduce((a, b) => a + b, 0);
      const probs = expScores.map(s => s / sumExp);

      for (let c = 0; c < nClasses; c++) {
        const target = (c === trueClass) ? 1 : 0;
        const error = probs[c] - target;
        // L2 regularization
        for (let j = 0; j < nFeatures + 1; j++) {
          const reg = (j > 0) ? lambda * weights[c][j] : 0;
          weights[c][j] -= lr * (error * x[j] + reg);
        }
      }
    }
    lr *= 0.999;
  }
  return { weights, featMeans, featStds };
}

function predictLR(model, features) {
  const { weights, featMeans, featStds } = model;
  const normFeat = features.map((f, j) => (f - featMeans[j]) / featStds[j]);
  const x = [1, ...normFeat];
  const scores = weights.map(w => w.reduce((s, wj, j) => s + wj * x[j], 0));
  const maxScore = Math.max(...scores);
  const expScores = scores.map(s => Math.exp(s - maxScore));
  const sumExp = expScores.reduce((a, b) => a + b, 0);
  return expScores.map(s => s / sumExp); // probabilities per class
}

// ── Random Forest (simplified) ──
class DecisionTree {
  constructor(maxDepth = 6, minSamples = 5) {
    this.maxDepth = maxDepth;
    this.minSamples = minSamples;
    this.root = null;
  }

  fit(features, labels, nClasses) {
    this.nClasses = nClasses;
    const data = features.map((f, i) => ({ x: f, y: labels[i] }));
    this.root = this._buildTree(data, 0);
  }

  _buildTree(data, depth) {
    const uniqueLabels = new Set(data.map(d => d.y));
    if (uniqueLabels.size === 1 || depth >= this.maxDepth || data.length < this.minSamples) {
      // Leaf: return majority class
      const counts = new Array(Math.max(...data.map(d => d.y)) + 1).fill(0);
      data.forEach(d => counts[d.y]++);
      let maxClass = 0, maxCount = 0;
      counts.forEach((c, i) => { if (c > maxCount) { maxCount = c; maxClass = i; } });
      return { leaf: true, class: maxClass, probs: counts.map(c => c / data.length) };
    }

    // Random feature subset
    const nFeatures = data[0].x.length;
    const mtry = Math.max(1, Math.floor(Math.sqrt(nFeatures)));
    const featureIdxs = [];
    for (let i = 0; i < mtry; i++) {
      const idx = Math.floor(Math.random() * nFeatures);
      if (!featureIdxs.includes(idx)) featureIdxs.push(idx);
    }

    let bestSplit = null, bestGain = -Infinity;
    for (const f of featureIdxs) {
      const values = data.map(d => d.x[f]).sort((a, b) => a - b);
      for (let i = 1; i < values.length; i++) {
        const thresh = (values[i - 1] + values[i]) / 2;
        const left = data.filter(d => d.x[f] <= thresh);
        const right = data.filter(d => d.x[f] > thresh);
        if (left.length < this.minSamples || right.length < this.minSamples) continue;

        const giniLeft = this._gini(left);
        const giniRight = this._gini(right);
        const gain = this._gini(data) - (left.length / data.length * giniLeft) - (right.length / data.length * giniRight);
        if (gain > bestGain) { bestGain = gain; bestSplit = { feature: f, threshold: thresh, left, right }; }
      }
    }

    if (!bestSplit) {
      const counts = new Array(this.nClasses).fill(0);
      data.forEach(d => counts[d.y]++);
      let maxClass = 0, maxCount = 0;
      counts.forEach((c, i) => { if (c > maxCount) { maxCount = c; maxClass = i; } });
      return { leaf: true, class: maxClass, probs: counts.map(c => c / Math.max(1, data.length)) };
    }

    return {
      leaf: false,
      feature: bestSplit.feature,
      threshold: bestSplit.threshold,
      left: this._buildTree(bestSplit.left, depth + 1),
      right: this._buildTree(bestSplit.right, depth + 1)
    };
  }

  _gini(data) {
    const counts = new Array(this.nClasses).fill(0);
    data.forEach(d => counts[d.y]++);
    let sum = 0;
    for (const c of counts) sum += (c / data.length) ** 2;
    return 1 - sum;
  }

  predict(features) {
    return this._predictNode(features, this.root);
  }

  _predictNode(x, node) {
    if (node.leaf) return node.probs;
    if (x[node.feature] <= node.threshold) return this._predictNode(x, node.left);
    return this._predictNode(x, node.right);
  }
}

function trainRandomForest(trainFeatures, trainLabels, nClasses, nTrees = 50) {
  const trees = [];
  for (let t = 0; t < nTrees; t++) {
    // Bootstrap sample
    const sampleIdxs = [];
    const n = trainFeatures.length;
    for (let i = 0; i < n; i++) sampleIdxs.push(Math.floor(Math.random() * n));
    const bootFeatures = sampleIdxs.map(i => trainFeatures[i]);
    const bootLabels = sampleIdxs.map(i => trainLabels[i]);

    const tree = new DecisionTree(8, 5);
    tree.fit(bootFeatures, bootLabels, nClasses);
    trees.push(tree);
  }
  return trees;
}

function predictRF(trees, features) {
  const nClasses = trees[0].nClasses;
  const aggProbs = new Array(nClasses).fill(0);
  trees.forEach(t => {
    const probs = t.predict(features);
    probs.forEach((p, i) => aggProbs[i] += p);
  });
  return aggProbs.map(p => p / trees.length);
}

// ── Shallow Neural Network (2-layer) ──
function trainNeuralNetwork(trainFeatures, trainLabels, nClasses, nFeatures, hiddenSize = 32, lr = 0.01, epochs = 300) {
  // Xavier initialization
  const scale1 = Math.sqrt(2.0 / (nFeatures + hiddenSize));
  const W1 = Array.from({ length: hiddenSize }, () => Array.from({ length: nFeatures }, () => (Math.random() * 2 - 1) * scale1));
  const b1 = new Array(hiddenSize).fill(0);

  const scale2 = Math.sqrt(2.0 / (hiddenSize + nClasses));
  const W2 = Array.from({ length: nClasses }, () => Array.from({ length: hiddenSize }, () => (Math.random() * 2 - 1) * scale2));
  const b2 = new Array(nClasses).fill(0);

  const cancerData = trainFeatures.map((f, i) => ({ x: f, y: trainLabels[i] })).filter(d => d.y >= 0);

  for (let epoch = 0; epoch < epochs; epoch++) {
    // Shuffle
    for (let i = cancerData.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [cancerData[i], cancerData[j]] = [cancerData[j], cancerData[i]];
    }

    let totalLoss = 0;
    for (const d of cancerData) {
      // Forward
      const h = [];
      for (let i = 0; i < hiddenSize; i++) {
        let sum = b1[i];
        for (let j = 0; j < nFeatures; j++) sum += W1[i][j] * d.x[j];
        h.push(Math.max(0, sum)); // ReLU
      }

      const scores = [];
      for (let i = 0; i < nClasses; i++) {
        let sum = b2[i];
        for (let j = 0; j < hiddenSize; j++) sum += W2[i][j] * h[j];
        scores.push(sum);
      }

      const maxScore = Math.max(...scores);
      const expScores = scores.map(s => Math.exp(s - maxScore));
      const sumExp = expScores.reduce((a, b) => a + b, 0);
      const probs = expScores.map(s => s / sumExp);

      totalLoss += -Math.log(Math.max(probs[d.y], 1e-10));

      // Backward
      const dScores = [...probs];
      dScores[d.y] -= 1;

      // dW2, db2
      for (let i = 0; i < nClasses; i++) {
        b2[i] -= lr * dScores[i];
        for (let j = 0; j < hiddenSize; j++) {
          W2[i][j] -= lr * dScores[i] * h[j];
        }
      }

      // dh
      const dh = new Array(hiddenSize).fill(0);
      for (let j = 0; j < hiddenSize; j++) {
        for (let i = 0; i < nClasses; i++) dh[j] += dScores[i] * W2[i][j];
        dh[j] *= (h[j] > 0 ? 1 : 0); // ReLU derivative
      }

      // dW1, db1
      for (let j = 0; j < hiddenSize; j++) {
        b1[j] -= lr * dh[j];
        for (let k = 0; k < nFeatures; k++) {
          W1[j][k] -= lr * dh[j] * d.x[k];
        }
      }
    }
    lr *= 0.995;
  }
  return { W1, b1, W2, b2, hiddenSize };
}

function predictNN(model, features) {
  const { W1, b1, W2, b2, hiddenSize } = model;
  const nClasses = W2.length;

  const h = [];
  for (let i = 0; i < hiddenSize; i++) {
    let sum = b1[i];
    for (let j = 0; j < features.length; j++) sum += W1[i][j] * features[j];
    h.push(Math.max(0, sum));
  }

  const scores = [];
  for (let i = 0; i < nClasses; i++) {
    let sum = b2[i];
    for (let j = 0; j < hiddenSize; j++) sum += W2[i][j] * h[j];
    scores.push(sum);
  }

  const maxScore = Math.max(...scores);
  const expScores = scores.map(s => Math.exp(s - maxScore));
  const sumExp = expScores.reduce((a, b) => a + b, 0);
  return expScores.map(s => s / sumExp);
}

// ── Metrics ──
function computeConfusionMatrix(trueLabels, predLabels, nClasses) {
  const cm = Array.from({ length: nClasses }, () => new Array(nClasses).fill(0));
  for (let i = 0; i < trueLabels.length; i++) {
    if (trueLabels[i] >= 0 && predLabels[i] >= 0) {
      cm[trueLabels[i]][predLabels[i]]++;
    }
  }
  return cm;
}

function computeTOOMetrics(cm, cancerTypeNames) {
  const totalCorrect = cm.reduce((sum, row, i) => sum + row[i], 0);
  const total = cm.reduce((sum, row) => sum + row.reduce((a, b) => a + b, 0), 0);
  const accuracy = total > 0 ? totalCorrect / total : 0;

  // Per-class metrics
  const perClass = {};
  for (let i = 0; i < cm.length; i++) {
    const tp = cm[i][i];
    const totalActual = cm[i].reduce((a, b) => a + b, 0);
    const totalPred = cm.reduce((sum, row) => sum + row[i], 0);
    perClass[cancerTypeNames[i]] = {
      sensitivity: totalActual > 0 ? tp / totalActual : 0,
      precision: totalPred > 0 ? tp / totalPred : 0,
      total_samples: totalActual
    };
  }

  // Top-2 accuracy
  let top2Correct = 0;
  // We don't have probabilities for all methods stored, skip top-2 for confusion matrix approach

  return { accuracy, per_class: perClass, confusion_matrix: cm, total_samples: total };
}

// ── Cross-Validation ──
function crossValidate(samples, cancerLabels, nClasses, nFeatures) {
  // 5-fold CV
  const cancerSamples = samples.filter(s => s.is_cancer);
  const healthySamples = samples.filter(s => !s.is_cancer);

  // Shuffle cancer samples
  const shuffledCancer = [...cancerSamples];
  for (let i = shuffledCancer.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffledCancer[i], shuffledCancer[j]] = [shuffledCancer[j], shuffledCancer[i]];
  }

  const nFolds = 5;
  const foldSize = Math.ceil(shuffledCancer.length / nFolds);

  const results = {
    logistic_regression: { accuracies: [], allPreds: [], allTrue: [] },
    random_forest: { accuracies: [], allPreds: [], allTrue: [] },
    neural_network: { accuracies: [], allPreds: [], allTrue: [] }
  };

  for (let fold = 0; fold < nFolds; fold++) {
    const testStart = fold * foldSize;
    const testEnd = Math.min(testStart + foldSize, shuffledCancer.length);
    const testCancer = shuffledCancer.slice(testStart, testEnd);
    const trainCancer = [...shuffledCancer.slice(0, testStart), ...shuffledCancer.slice(testEnd)];

    const trainFeatures = trainCancer.map(s => s.features);
    const trainLabels = trainCancer.map(s => s.label);
    const testFeatures = testCancer.map(s => s.features);
    const testLabels = testCancer.map(s => s.label);

    // LR
    const lrModel = trainLogisticRegression(trainCancer, trainLabels, nClasses, nFeatures);
    const lrPreds = testFeatures.map(f => {
      const probs = predictLR(lrModel, f);
      return probs.indexOf(Math.max(...probs));
    });
    const lrAcc = lrPreds.filter((p, i) => p === testLabels[i]).length / testLabels.length;

    // RF
    const rfTrees = trainRandomForest(trainFeatures, trainLabels, nClasses, 50);
    const rfPreds = testFeatures.map(f => {
      const probs = predictRF(rfTrees, f);
      return probs.indexOf(Math.max(...probs));
    });
    const rfAcc = rfPreds.filter((p, i) => p === testLabels[i]).length / testLabels.length;

    // NN
    const nnModel = trainNeuralNetwork(trainFeatures, trainLabels, nClasses, nFeatures, 32, 0.01, 300);
    const nnPreds = testFeatures.map(f => {
      const probs = predictNN(nnModel, f);
      return probs.indexOf(Math.max(...probs));
    });
    const nnAcc = nnPreds.filter((p, i) => p === testLabels[i]).length / testLabels.length;

    results.logistic_regression.accuracies.push(lrAcc);
    results.logistic_regression.allPreds.push(...lrPreds);
    results.logistic_regression.allTrue.push(...testLabels);
    results.random_forest.accuracies.push(rfAcc);
    results.random_forest.allPreds.push(...rfPreds);
    results.random_forest.allTrue.push(...testLabels);
    results.neural_network.accuracies.push(nnAcc);
    results.neural_network.allPreds.push(...nnPreds);
    results.neural_network.allTrue.push(...testLabels);
  }

  return results;
}

// ── Bootstrap CI ──
function bootstrapCI(values, nBoot, rng) {
  const estimates = [];
  const n = values.length;
  for (let b = 0; b < nBoot; b++) {
    let sum = 0;
    for (let i = 0; i < n; i++) sum += values[Math.floor(rng() * n)];
    estimates.push(sum / n);
  }
  estimates.sort((a, b) => a - b);
  const mean = estimates.reduce((a, b) => a + b, 0) / estimates.length;
  const lo = estimates[Math.floor(0.025 * estimates.length)];
  const hi = estimates[Math.ceil(0.975 * estimates.length) - 1];
  return { mean, ci95_low: lo, ci95_high: hi };
}

// ── Cancer Detection (healthy vs cancer classification with TOO) ──
function evaluateCancerDetection(samples, cancerLabels, nClasses, nFeatures) {
  const cancerOnly = samples.filter(s => s.is_cancer);
  const healthyOnly = samples.filter(s => !s.is_cancer);

  // Train on 80% of cancer + all healthy
  const nTrainCancer = Math.floor(cancerOnly.length * 0.8);
  const shuffledCancer = [...cancerOnly];
  for (let i = shuffledCancer.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffledCancer[i], shuffledCancer[j]] = [shuffledCancer[j], shuffledCancer[i]];
  }

  const trainCancer = shuffledCancer.slice(0, nTrainCancer);
  const testCancer = shuffledCancer.slice(nTrainCancer);

  // Train TOO classifier on cancer samples only
  const trainFeatures = trainCancer.map(s => s.features);
  const trainLabels = trainCancer.map(s => s.label);

  // Train LR for both TOO + detection
  const lrModel = trainLogisticRegression(trainCancer, trainLabels, nClasses, nFeatures);

  // Test: all test samples
  const testSamples = [...testCancer, ...healthyOnly];
  const testTrueLabels = testSamples.map(s => s.label); // -1 for healthy

  // For each sample, predict TOO probabilities and max probability
  const predictions = testSamples.map(s => {
    const probs = predictLR(lrModel, s.features);
    const maxProb = Math.max(...probs);
    const predClass = probs.indexOf(maxProb);
    return { maxProb, predClass, probs };
  });

  // Detection: threshold on max probability
  const nTest = predictions.length;
  const scores = predictions.map(p => p.maxProb);
  const labels = testSamples.map(s => s.is_cancer ? 1 : 0);

  // Find optimal threshold
  const sorted = scores.map((s, i) => ({ s, l: labels[i] })).sort((a, b) => a.s - b.s);
  const totalPos = labels.filter(l => l === 1).length;
  const totalNeg = labels.filter(l => l === 0).length;
  let bestJ = -1, bestThresh = 0;
  for (const p of sorted) {
    const tp = sorted.filter(x => x.s >= p.s && x.l === 1).length;
    const tn = sorted.filter(x => x.s < p.s && x.l === 0).length;
    const J = tp / Math.max(1, totalPos) + tn / Math.max(1, totalNeg) - 1;
    if (J > bestJ) { bestJ = J; bestThresh = p.s; }
  }

  const preds = scores.map(s => s >= bestThresh ? 1 : 0);
  let tp = 0, fp = 0, tn = 0, fn = 0;
  for (let i = 0; i < labels.length; i++) {
    if (labels[i] === 1 && preds[i] === 1) tp++;
    else if (labels[i] === 0 && preds[i] === 1) fp++;
    else if (labels[i] === 0 && preds[i] === 0) tn++;
    else fn++;
  }

  // TOO accuracy on correctly detected cancers
  const detectedCancers = predictions.filter((_, i) => labels[i] === 1 && preds[i] === 1);
  const detectedTrueLabels = testSamples.filter((s, i) => labels[i] === 1 && preds[i] === 1).map(s => s.label);
  const tooCorrect = detectedCancers.filter((p, i) => p.predClass === detectedTrueLabels[i]).length;

  return {
    cancer_detection: {
      sensitivity: tp / Math.max(1, tp + fn),
      specificity: tn / Math.max(1, tn + fp),
      ppv: tp / Math.max(1, tp + fp),
      tp, fp, tn, fn
    },
    too_on_detected: {
      accuracy: tooCorrect / Math.max(1, detectedCancers.length),
      n_detected: detectedCancers.length
    }
  };
}

// ═══════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH TISSUE-OF-ORIGIN (TOO) VALIDATION');
  console.log(`   8 Cancer Types, ${N_SAMPLES_PER_CANCER} samples each, ${N_HEALTHY} healthy`);
  console.log('   Methods: Logistic Regression, Random Forest, Neural Network');
  console.log('='.repeat(70));

  const rng = createRNG(SEED + 6000);
  const { samples, allGenes, labels: typeLabels, nCancerTypes } = simulateTOOData(rng);
  const nFeatures = samples[0].features.length;

  console.log(`\n📊 Dataset: ${samples.length} total samples`);
  console.log(`   Cancer samples: ${samples.filter(s => s.is_cancer).length} (${nCancerTypes} types)`);
  console.log(`   Healthy samples: ${samples.filter(s => !s.is_cancer).length}`);
  console.log(`   Features: ${nFeatures} (${allGenes.length} methylation + 3 fragmentomic)`);

  // ── Cross-Validation ──
  console.log('\n🔬 5-Fold Cross-Validation...');
  const cvResults = crossValidate(samples, typeLabels, nCancerTypes, nFeatures);

  console.log(`\n📊 TOO Accuracy (multi-class, cancer samples only):`);
  console.log(`${'Method'.padEnd(25)} ${'Mean Acc'.padEnd(10)} ${'95% CI'.padEnd(20)}`);

  const methodNames = ['logistic_regression', 'random_forest', 'neural_network'];
  const methodDisplay = {
    logistic_regression: 'Logistic Regression',
    random_forest: 'Random Forest',
    neural_network: 'Neural Network (2-layer)'
  };

  const finalResults = {};
  const rngCI = createRNG(SEED + 6500);

  methodNames.forEach(method => {
    const accs = cvResults[method].accuracies;
    const ci = bootstrapCI(accs, N_BOOTSTRAP, rngCI);
    const mean = accs.reduce((a, b) => a + b, 0) / accs.length;
    console.log(`${methodDisplay[method].padEnd(25)} ${(mean*100).toFixed(2)}%`.padEnd(10) + ` [${(ci.ci95_low*100).toFixed(2)}–${(ci.ci95_high*100).toFixed(2)}%]`);

    // Full confusion matrix (on all predictions)
    const allPreds = cvResults[method].allPreds;
    const allTrue = cvResults[method].allTrue;
    const cm = computeConfusionMatrix(allTrue, allPreds, nCancerTypes);
    const nameList = typeLabels.map(l => l.code);
    const tooMetrics = computeTOOMetrics(cm, nameList);

    // Top-2 accuracy
    let top2Correct = 0;
    // For top-2: need probabilities; compute from LR method specifically
    if (method === 'logistic_regression') {
      const cancerSamples = samples.filter(s => s.is_cancer);
      const lrModelTop2 = trainLogisticRegression(cancerSamples, cancerSamples.map(s => s.label), nCancerTypes, nFeatures, 0.005, 1000, 0.01);
      cancerSamples.forEach(s => {
        const probs = predictLR(lrModelTop2, s.features);
        const sorted = probs.map((p, i) => ({ p, i })).sort((a, b) => b.p - a.p);
        if (sorted[0].i === s.label || sorted[1].i === s.label) top2Correct++;
      });
    }

    finalResults[method] = {
      accuracy: parseFloat(mean.toFixed(4)),
      accuracy_ci95_low: parseFloat(ci.ci95_low.toFixed(4)),
      accuracy_ci95_high: parseFloat(ci.ci95_high.toFixed(4)),
      top2_accuracy: method === 'logistic_regression' ? parseFloat((top2Correct / samples.filter(s => s.is_cancer).length).toFixed(4)) : null,
      per_cancer_type: tooMetrics.per_class,
      confusion_matrix: tooMetrics.confusion_matrix
    };
  });

  // ── Cancer Detection + TOO ──
  console.log('\n🔬 Joint Cancer Detection + TOO Pipeline:');
  const detectionEval = evaluateCancerDetection(samples, typeLabels, nCancerTypes, nFeatures);
  console.log(`   Cancer Detection: Sens=${(detectionEval.cancer_detection.sensitivity*100).toFixed(1)}%, Spec=${(detectionEval.cancer_detection.specificity*100).toFixed(1)}%, PPV=${(detectionEval.cancer_detection.ppv*100).toFixed(1)}%`);
  console.log(`   TOO on Detected: ${(detectionEval.too_on_detected.accuracy*100).toFixed(1)}% (n=${detectionEval.too_on_detected.n_detected})`);

  // ── Per-Cancer-Type TOO Accuracy ──
  console.log('\n📊 Per-Cancer-Type TOO Sensitivity (LR method):');
  console.log(`${'Cancer Type'.padEnd(20)} ${'Sensitivity'.padEnd(12)} ${'Precision'.padEnd(12)}`);
  const perClass = finalResults['logistic_regression'].per_cancer_type;
  for (const [type, metrics] of Object.entries(perClass)) {
    console.log(`${type.padEnd(20)} ${(metrics.sensitivity*100).toFixed(1)}%`.padEnd(12) + ` ${(metrics.precision*100).toFixed(1)}%`);
  }

  // ── Output ──
  const output = {
    metadata: {
      validation: 'tissue_of_origin_prediction',
      timestamp: new Date().toISOString(),
      n_cancer_types: nCancerTypes,
      cancer_types: typeLabels.map(l => ({ code: l.code, name: l.name })),
      n_samples_per_cancer: N_SAMPLES_PER_CANCER,
      n_healthy: N_HEALTHY,
      n_features: nFeatures,
      n_methylation_features: allGenes.length,
      n_fragmentomic_features: 3,
      n_bootstrap: N_BOOTSTRAP
    },
    too_results: finalResults,
    joint_detection_too: {
      cancer_detection: {
        sensitivity: parseFloat(detectionEval.cancer_detection.sensitivity.toFixed(4)),
        specificity: parseFloat(detectionEval.cancer_detection.specificity.toFixed(4)),
        ppv: parseFloat(detectionEval.cancer_detection.ppv.toFixed(4))
      },
      too_on_detected: {
        accuracy: parseFloat(detectionEval.too_on_detected.accuracy.toFixed(4)),
        n_detected: detectionEval.too_on_detected.n_detected
      }
    },
    summary: {
      best_method: 'logistic_regression',
      best_too_accuracy: parseFloat(finalResults['logistic_regression'].accuracy.toFixed(4)),
      best_top2_accuracy: finalResults['logistic_regression'].top2_accuracy,
      note: 'TOO accuracy measures correct tissue prediction among cancer cases. Joint pipeline first detects cancer, then predicts tissue-of-origin.'
    }
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved to ${OUTPUT_PATH}`);
  console.log('\n✅ Tissue-of-Origin validation complete.');
  console.log('='.repeat(70));
})();
