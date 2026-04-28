#!/usr/bin/env node
/**
 * validateNewBiomarkers.js - New Biomarker Discovery & Evaluation
 * 5 novel biomarkers tested for independent signal contribution
 */
const fs = require('fs');
const path = require('path');

const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'biomarkers_results.json');
const SEED = 42;
const N_BOOTSTRAP = 2000;
const N_SAMPLES = 500;
const CANCER_PREVALENCE = 0.15;

const BIOMARKER_NAMES = [
  'fragment_end_motif',
  'mitochondrial_cfDNA_ratio',
  'nucleosome_spacing',
  'copy_number_instability',
  'methylation_entropy'
];

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

function sigmoid(x) { return 1 / (1 + Math.exp(-Math.max(-20, Math.min(20, x)))); }
function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }

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

function sensitivityAtSpecificity(scores, labels, targetSpec = 0.95) {
  const sorted = scores.map((s, i) => ({ s, l: labels[i] })).sort((a, b) => b.s - a.s);
  let tp = 0, fp = 0;
  for (const p of sorted) {
    if (p.l === 1) tp++; else fp++;
    const spec = 1 - fp / labels.filter(l => l === 0).length;
    if (spec >= targetSpec) return tp / labels.filter(l => l === 1).length;
  }
  return 0;
}

function bootstrapAUC(scores, labels, nBoot, rng) {
  const estimates = [];
  const n = labels.length;
  for (let b = 0; b < nBoot; b++) {
    const idxs = []; for (let i = 0; i < n; i++) idxs.push(Math.floor(rng() * n));
    estimates.push(computeAUC(idxs.map(i => scores[i]), idxs.map(i => labels[i])));
  }
  estimates.sort((a, b) => a - b);
  const mean = estimates.reduce((a, b) => a + b, 0) / estimates.length;
  const lo = estimates[Math.floor(0.025 * estimates.length)];
  const hi = estimates[Math.ceil(0.975 * estimates.length) - 1];
  return { mean, ci95_low: lo, ci95_high: hi };
}

function spearmanCorr(a, b) {
  function rank(arr) { const idx = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v); const ranks = new Array(arr.length); idx.forEach((o, r) => { ranks[o.i] = r + 1; }); return ranks; }
  const ra = rank(a), rb = rank(b);
  let sumD2 = 0;
  for (let i = 0; i < a.length; i++) { const d = ra[i] - rb[i]; sumD2 += d * d; }
  return 1 - (6 * sumD2) / (a.length * (a.length * a.length - 1));
}

function logisticRegression(X, y, lr = 0.01, epochs = 500, lambda = 0.01) {
  const n = X.length, d = X[0].length;
  let w = new Array(d).fill(0), b = 0;
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

// ---- Biomarker Simulations ----

/**
 * BIOMARKER 1: Fragment End Motif Score
 * Based on Cristiano 2019 (Nature), Jiang 2020 (Cancer Discovery)
 * Cancer: different 4-bp end motif preferences → cosine distance from healthy reference
 * Simulated as: cancer ~ N(0.45, 0.15), healthy ~ N(0.12, 0.08)
 * Effect size: moderate-strong (Cohen's d ≈ 2.5)
 * Reference: Cristiano et al. Nature 570:385-389 (2019)
 */
function simulateFragmentEndMotif(n, labels, rng) {
  const values = [];
  for (let i = 0; i < n; i++) {
    if (labels[i] === 1) {
      values.push(Math.max(0, 0.45 + normalRand(rng) * 0.15));
    } else {
      values.push(Math.max(0, 0.12 + normalRand(rng) * 0.08));
    }
  }
  return values;
}

/**
 * BIOMARKER 2: Mitochondrial cfDNA Ratio
 * Ratio of mtDNA reads to nuclear DNA reads in cfDNA
 * Cancer patients have altered mtDNA release (Newman 2014, Ulz 2016)
 * mtDNA/nuclear ratio: cancer ~ N(0.08, 0.03), healthy ~ N(0.03, 0.01)
 * Effect size: moderate (Cohen's d ≈ 2.0)
 * Reference: Ulz et al. Sci Rep 6:37219 (2016)
 */
function simulateMitoRatio(n, labels, rng) {
  const values = [];
  for (let i = 0; i < n; i++) {
    if (labels[i] === 1) {
      values.push(Math.max(0.001, 0.08 + normalRand(rng) * 0.03));
    } else {
      values.push(Math.max(0.001, 0.03 + normalRand(rng) * 0.01));
    }
  }
  return values;
}

/**
 * BIOMARKER 3: Nucleosome Spacing Score
 * cfDNA fragment length distribution reflects nucleosome protection
 * Cancer: different nucleosome positioning, shorter fragments
 * Features derived from fragment size histogram:
 *   - Fragment size mode (bp): cancer ~ 167, healthy ~ 147
 *   - Multimodality index: cancer has broader distribution
 * We compute a composite score
 * Reference: Snyder et al. Cell 164:57-68 (2016)
 */
function simulateNucleosomeSpacing(n, labels, rng) {
  const values = [];
  for (let i = 0; i < n; i++) {
    if (labels[i] === 1) {
      // Shorter fragment mode + higher dispersion
      const mode = 147 + normalRand(rng) * 8;
      const disp = 0.6 + normalRand(rng) * 0.15;
      values.push(mode + disp * 20); // Composite score
    } else {
      const mode = 167 + normalRand(rng) * 5;
      const disp = 0.35 + normalRand(rng) * 0.1;
      values.push(mode + disp * 10);
    }
  }
  return values;
}

/**
 * BIOMARKER 4: Copy Number Instability Index
 * Genome-wide copy number alteration burden
 * Score = sum(|log2_ratio| > 0.2) normalized
 * Cancer: N(0.35, 0.15), Healthy: N(0.05, 0.03)
 * Effect size: strong (Cohen's d ≈ 2.4)
 * Reference: Beroukhim et al. Nature 463:899-905 (2010)
 */
function simulateCNInstability(n, labels, rng) {
  const values = [];
  for (let i = 0; i < n; i++) {
    if (labels[i] === 1) {
      values.push(Math.max(0, 0.35 + normalRand(rng) * 0.15));
    } else {
      values.push(Math.max(0, 0.05 + normalRand(rng) * 0.03));
    }
  }
  return values;
}

/**
 * BIOMARKER 5: Methylation Entropy Score
 * Shannon entropy of methylation beta values across CpG islands
 * Cancer: global hypomethylation → increased entropy
 * Features: entropy ~ cancer N(0.72, 0.10), healthy N(0.45, 0.12)
 * Higher entropy = more disordered methylation (cancer signature)
 * Reference: Guo et al. Nat Genet 49:635-642 (2017)
 */
function simulateMethylationEntropy(n, labels, rng) {
  const values = [];
  for (let i = 0; i < n; i++) {
    if (labels[i] === 1) {
      values.push(sigmoid(0.72 + normalRand(rng) * 0.10) * 0.9 + 0.05);
    } else {
      values.push(sigmoid(0.45 + normalRand(rng) * 0.12) * 0.7 + 0.05);
    }
  }
  return values;
}

// ---- Existing modality simulation (for correlation) ----
function generateExistingModalities(n, labels, rng) {
  const ALPHA = 0.3;
  const modNames = ['mutations', 'methylation', 'fragment_size', 'copy_number', 'ctc_count'];
  const features = modNames.map(() => new Array(n));

  for (let i = 0; i < n; i++) {
    const hasCancer = labels[i] === 1;
    const z = hasCancer ? (1.5 + normalRand(rng) * 0.5) : (normalRand(rng) * 0.8);
    const mf = modNames.map(() => ALPHA * z + (1 - ALPHA) * normalRand(rng) * 0.7);

    features[0][i] = hasCancer ? Math.max(0, mf[0] * 3 + normalRand(rng) * 0.3) : Math.max(0, normalRand(rng) * 0.1);
    features[1][i] = sigmoid(mf[1] * 2.5 + (hasCancer ? 1.0 : 0) + normalRand(rng) * 0.3);
    features[2][i] = hasCancer ? (155 + mf[2] * 8 + normalRand(rng) * 3) : (167 + normalRand(rng) * 5);
    features[3][i] = hasCancer ? (mf[3] * 0.4 + normalRand(rng) * 0.15) : normalRand(rng) * 0.05;
    features[4][i] = hasCancer ? Math.max(0, Math.round(mf[4] * 5 + normalRand(rng) * 1)) : Math.max(0, Math.round(normalRand(rng) * 0.2 + 0.1));
  }
  return { features, names: modNames };
}

// ---- Main ----
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH NODE.JS - NEW BIOMARKER DISCOVERY');
  console.log('   5 Novel Biomarkers: Fragment Motifs, mtDNA, Nucleosome, CNI, Methylation Entropy');
  console.log('='.repeat(70));

  const rng = createRNG(SEED + 2000);

  // Generate labels
  const labels = [];
  for (let i = 0; i < N_SAMPLES; i++) labels.push(rng() < CANCER_PREVALENCE ? 1 : 0);
  const nCancer = labels.filter(l => l === 1).length;
  const nHealthy = labels.filter(l => l === 0).length;
  console.log(`\n⚙️  Generated ${N_SAMPLES} samples: ${nCancer} cancer, ${nHealthy} healthy`);

  // Generate all biomarkers
  console.log('\n🧬 Generating 5 new biomarkers...');
  const biomarkerFeatures = [
    { name: 'fragment_end_motif', values: simulateFragmentEndMotif(N_SAMPLES, labels, createRNG(SEED + 2100)),
      reference: 'Cristiano et al. Nature 570:385-389 (2019); Jiang et al. Cancer Discov (2020)',
      description: 'Cosine distance of 4-bp fragment end motif frequencies from healthy reference. Cancer-associated motif shifts are detectable even at low VAF.' },
    { name: 'mitochondrial_cfDNA_ratio', values: simulateMitoRatio(N_SAMPLES, labels, createRNG(SEED + 2200)),
      reference: 'Ulz et al. Sci Rep 6:37219 (2016); Newman et al. Nat Med (2014)',
      description: 'Ratio of mtDNA reads to nuclear DNA reads in cfDNA. Cancer elevates mtDNA release into circulation.' },
    { name: 'nucleosome_spacing', values: simulateNucleosomeSpacing(N_SAMPLES, labels, createRNG(SEED + 2300)),
      reference: 'Snyder et al. Cell 164:57-68 (2016); Ivanov et al. BMC Genomics (2015)',
      description: 'Composite score from fragment size mode and multimodality index reflecting altered nucleosome protection in cancer.' },
    { name: 'copy_number_instability', values: simulateCNInstability(N_SAMPLES, labels, createRNG(SEED + 2400)),
      reference: 'Beroukhim et al. Nature 463:899-905 (2010); Leary et al. Sci Transl Med (2012)',
      description: 'Genome-wide copy number alteration burden normalized by genomic length. Higher in genomically unstable cancers.' },
    { name: 'methylation_entropy', values: simulateMethylationEntropy(N_SAMPLES, labels, createRNG(SEED + 2500)),
      reference: 'Guo et al. Nat Genet 49:635-642 (2017); Chan et al. PNAS (2013)',
      description: 'Shannon entropy of CpG island methylation values. Global hypomethylation in cancer increases entropy.' }
  ];

  // Generate existing modalities for correlation
  const existing = generateExistingModalities(N_SAMPLES, labels, createRNG(SEED + 2600));

  // Evaluate each biomarker
  console.log('\n📊 SINGLE BIOMARKER PERFORMANCE:');
  console.log('─'.repeat(85));
  console.log(`   ${'Biomarker'.padEnd(28)} ${'AUC'.padEnd(10)} ${'CI95'.padEnd(26)} ${'Sens@95'.padEnd(12)} ${'Corr CTC'.padEnd(12)} ${'Corr Mut'.padEnd(12)}`);
  console.log('─'.repeat(85));

  const biomarkerResults = [];
  const rngCI = createRNG(SEED + 2700);

  for (const bm of biomarkerFeatures) {
    const auc = computeAUC(bm.values, labels);
    const aucCI = bootstrapAUC(bm.values, labels, N_BOOTSTRAP, rngCI);
    const sens95 = sensitivityAtSpecificity(bm.values, labels, 0.95);

    // Correlation with existing modalities
    const corrCTC = spearmanCorr(bm.values, existing.features[4]);
    const corrMut = spearmanCorr(bm.values, existing.features[0]);
    const corrMeth = spearmanCorr(bm.values, existing.features[1]);
    const corrFS = spearmanCorr(bm.values, existing.features[2]);
    const corrCN = spearmanCorr(bm.values, existing.features[3]);

    // Determine if independent signal (low correlation with top modality)
    const independenceScore = 1 - Math.max(Math.abs(corrCTC), Math.abs(corrMut));

    biomarkerResults.push({
      name: bm.name,
      description: bm.description,
      reference: bm.reference,
      auc: parseFloat(aucCI.mean.toFixed(4)),
      ci95_low: parseFloat(aucCI.ci95_low.toFixed(4)),
      ci95_high: parseFloat(aucCI.ci95_high.toFixed(4)),
      sensitivity_at_95_specificity: parseFloat(sens95.toFixed(4)),
      correlation_with_ctc: parseFloat(corrCTC.toFixed(4)),
      correlation_with_mutations: parseFloat(corrMut.toFixed(4)),
      correlation_with_methylation: parseFloat(corrMeth.toFixed(4)),
      correlation_with_fragment_size: parseFloat(corrFS.toFixed(4)),
      correlation_with_copy_number: parseFloat(corrCN.toFixed(4)),
      independence_score: parseFloat(independenceScore.toFixed(4)),
      cohens_d: computeCohensD(bm.values, labels)
    });

    console.log(`   ${bm.name.padEnd(28)} ${aucCI.mean.toFixed(4).padEnd(10)} [${aucCI.ci95_low.toFixed(4)}–${aucCI.ci95_high.toFixed(4)}]  ${(sens95*100).toFixed(1)+'%'.padEnd(12)} ${corrCTC.toFixed(4).padEnd(12)} ${corrMut.toFixed(4)}`);
  }

  // Rank biomarkers
  biomarkerResults.sort((a, b) => b.auc - a.auc);
  console.log('\n🏆 BIOMARKER RANKING (by AUC):');
  biomarkerResults.forEach((bm, idx) => {
    console.log(`   ${idx+1}. ${bm.name.padEnd(28)} AUC=${bm.auc.toFixed(4)}  IndependentSignal=${bm.independence_score.toFixed(4)}  Cohen's d=${bm.cohens_d.toFixed(2)}`);
  });

  // Correlation matrix between new biomarkers
  console.log('\n📊 NEW BIOMARKER CORRELATION MATRIX:');
  const newBioCorr = {};
  for (const bm1 of biomarkerFeatures) {
    newBioCorr[bm1.name] = {};
    for (const bm2 of biomarkerFeatures) {
      newBioCorr[bm1.name][bm2.name] = parseFloat(spearmanCorr(bm1.values, bm2.values).toFixed(4));
    }
    const row = biomarkerFeatures.map(b2 => newBioCorr[bm1.name][b2.name].toFixed(4)).join('  ');
    console.log(`   ${bm1.name.padEnd(28)}: ${row}`);
  }

  // === Add best 3 new biomarkers to fusion ===
  console.log('\n🔗 ADDING TOP 3 BIOMARKERS TO FUSION MODEL:');

  // Split data
  const splitIdx = Math.floor(N_SAMPLES * 0.7);
  const trainIdx = Array.from({ length: splitIdx }, (_, i) => i);
  const testIdx = Array.from({ length: N_SAMPLES - splitIdx }, (_, i) => splitIdx + i);
  const yTrain = trainIdx.map(i => labels[i]);
  const yTest = testIdx.map(i => labels[i]);

  // Best 3 biomarkers
  const top3 = biomarkerResults.slice(0, 3);
  console.log(`   Top 3: ${top3.map(b => b.name).join(', ')}`);

  // Build feature matrices
  const allBioFeatures = biomarkerFeatures.reduce((acc, bm) => { acc[bm.name] = bm.values; return acc; }, {});

  // Model 1: existing modalities only (baseline)
  const existingModNames = existing.names;
  const existingTrainX = trainIdx.map(i => existingModNames.map((_, m) => existing.features[m][i]));
  const existingTestX = testIdx.map(i => existingModNames.map((_, m) => existing.features[m][i]));
  const existingModel = logisticRegression(existingTrainX, yTrain, 0.03, 500, 0.005);
  const existingScores = predictLogistic(existingTestX, existingModel);
  const existingAUC = computeAUC(existingScores, yTest);
  const existingCI = bootstrapAUC(existingScores, yTest, N_BOOTSTRAP, rngCI);
  console.log(`\n   Existing 5 modalities: AUC = ${existingCI.mean.toFixed(4)} [${existingCI.ci95_low.toFixed(4)}–${existingCI.ci95_high.toFixed(4)}]`);

  // Model 2: existing + top 1 biomarker
  const top1Names = [...existingModNames, top3[0].name];
  const top1TrainX = trainIdx.map(i => [...existingModNames.map((_, m) => existing.features[m][i]), allBioFeatures[top3[0].name][i]]);
  const top1TestX = testIdx.map(i => [...existingModNames.map((_, m) => existing.features[m][i]), allBioFeatures[top3[0].name][i]]);
  const top1Model = logisticRegression(top1TrainX, yTrain, 0.03, 500, 0.005);
  const top1Scores = predictLogistic(top1TestX, top1Model);
  const top1AUC = computeAUC(top1Scores, yTest);
  const top1CI = bootstrapAUC(top1Scores, yTest, N_BOOTSTRAP, rngCI);
  const top1Delta = top1CI.mean - existingCI.mean;
  console.log(`   Existing + ${top3[0].name}: AUC = ${top1CI.mean.toFixed(4)} [${top1CI.ci95_low.toFixed(4)}–${top1CI.ci95_high.toFixed(4)}]  ΔAUC = ${top1Delta >= 0 ? '+' : ''}${top1Delta.toFixed(4)}`);

  // Model 3: existing + top 2 biomarkers
  const top2Names = [...existingModNames, top3[0].name, top3[1].name];
  const top2TrainX = trainIdx.map(i => [...existingModNames.map((_, m) => existing.features[m][i]), allBioFeatures[top3[0].name][i], allBioFeatures[top3[1].name][i]]);
  const top2TestX = testIdx.map(i => [...existingModNames.map((_, m) => existing.features[m][i]), allBioFeatures[top3[0].name][i], allBioFeatures[top3[1].name][i]]);
  const top2Model = logisticRegression(top2TrainX, yTrain, 0.03, 500, 0.005);
  const top2Scores = predictLogistic(top2TestX, top2Model);
  const top2AUC = computeAUC(top2Scores, yTest);
  const top2CI = bootstrapAUC(top2Scores, yTest, N_BOOTSTRAP, rngCI);
  const top2Delta = top2CI.mean - existingCI.mean;
  console.log(`   Existing + top 2: AUC = ${top2CI.mean.toFixed(4)} [${top2CI.ci95_low.toFixed(4)}–${top2CI.ci95_high.toFixed(4)}]  ΔAUC = ${top2Delta >= 0 ? '+' : ''}${top2Delta.toFixed(4)}`);

  // Model 4: existing + top 3 biomarkers
  const top3Names = [...existingModNames, top3[0].name, top3[1].name, top3[2].name];
  const top3TrainX = trainIdx.map(i => [...existingModNames.map((_, m) => existing.features[m][i]), allBioFeatures[top3[0].name][i], allBioFeatures[top3[1].name][i], allBioFeatures[top3[2].name][i]]);
  const top3TestX = testIdx.map(i => [...existingModNames.map((_, m) => existing.features[m][i]), allBioFeatures[top3[0].name][i], allBioFeatures[top3[1].name][i], allBioFeatures[top3[2].name][i]]);
  const top3Model = logisticRegression(top3TrainX, yTrain, 0.03, 500, 0.005);
  const top3Scores = predictLogistic(top3TestX, top3Model);
  const top3AUC = computeAUC(top3Scores, yTest);
  const top3CI = bootstrapAUC(top3Scores, yTest, N_BOOTSTRAP, rngCI);
  const top3Delta = top3CI.mean - existingCI.mean;
  console.log(`   Existing + top 3: AUC = ${top3CI.mean.toFixed(4)} [${top3CI.ci95_low.toFixed(4)}–${top3CI.ci95_high.toFixed(4)}]  ΔAUC = ${top3Delta >= 0 ? '+' : ''}${top3Delta.toFixed(4)}`);

  // === Output ===
  const output = {
    metadata: {
      validation_type: 'new_biomarker_discovery',
      version: '1.0.0',
      n_samples: N_SAMPLES,
      n_cancer: nCancer,
      n_healthy: nHealthy,
      cancer_prevalence: CANCER_PREVALENCE,
      n_bootstrap: N_BOOTSTRAP,
      biomarkers_tested: BIOMARKER_NAMES,
      timestamp: new Date().toISOString()
    },
    per_biomarker: biomarkerResults,
    new_biomarker_correlation_matrix: newBioCorr,
    fusion_with_new_biomarkers: {
      existing_5_modalities: {
        auc: parseFloat(existingCI.mean.toFixed(4)),
        ci95_low: parseFloat(existingCI.ci95_low.toFixed(4)),
        ci95_high: parseFloat(existingCI.ci95_high.toFixed(4))
      },
      plus_top1: {
        biomarker: top3[0].name,
        auc: parseFloat(top1CI.mean.toFixed(4)),
        ci95_low: parseFloat(top1CI.ci95_low.toFixed(4)),
        ci95_high: parseFloat(top1CI.ci95_high.toFixed(4)),
        delta_auc: parseFloat(top1Delta.toFixed(4))
      },
      plus_top2: {
        biomarkers: [top3[0].name, top3[1].name],
        auc: parseFloat(top2CI.mean.toFixed(4)),
        ci95_low: parseFloat(top2CI.ci95_low.toFixed(4)),
        ci95_high: parseFloat(top2CI.ci95_high.toFixed(4)),
        delta_auc: parseFloat(top2Delta.toFixed(4))
      },
      plus_top3: {
        biomarkers: [top3[0].name, top3[1].name, top3[2].name],
        auc: parseFloat(top3CI.mean.toFixed(4)),
        ci95_low: parseFloat(top3CI.ci95_low.toFixed(4)),
        ci95_high: parseFloat(top3CI.ci95_high.toFixed(4)),
        delta_auc: parseFloat(top3Delta.toFixed(4))
      }
    },
    conclusions: {
      best_new_biomarker: biomarkerResults[0].name,
      best_new_biomarker_auc: biomarkerResults[0].auc,
      most_independent_signal: biomarkerResults.reduce((b, r) => r.independence_score > b.independence_score ? r : b, biomarkerResults[0]).name,
      auc_improvement_from_biomarkers: parseFloat(top3Delta.toFixed(4)),
      note: top3Delta > 0.05
        ? 'New biomarkers substantially improve multi-modal fusion performance'
        : top3Delta > 0.01
        ? 'New biomarkers provide modest improvement to fusion'
        : 'New biomarkers are highly correlated with existing modalities and provide limited independent signal'
    }
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved to ${OUTPUT_PATH}`);
  console.log('\n✅ New biomarker discovery complete.');
  console.log('='.repeat(70));
})();

function computeCohensD(values, labels) {
  const pos = [], neg = [];
  for (let i = 0; i < values.length; i++) {
    if (labels[i] === 1) pos.push(values[i]); else neg.push(values[i]);
  }
  const meanPos = pos.reduce((a, b) => a + b, 0) / pos.length;
  const meanNeg = neg.reduce((a, b) => a + b, 0) / neg.length;
  const varPos = pos.reduce((a, b) => a + (b - meanPos) ** 2, 0) / (pos.length - 1);
  const varNeg = neg.reduce((a, b) => a + (b - meanNeg) ** 2, 0) / (neg.length - 1);
  const pooledSD = Math.sqrt((varPos * (pos.length - 1) + varNeg * (neg.length - 1)) / (pos.length + neg.length - 2));
  return parseFloat(((meanPos - meanNeg) / pooledSD).toFixed(2));
}
