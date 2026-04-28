#!/usr/bin/env node
/**
 * validateVariantCaller.js - Variant Caller Validation
 * Tests 4 variant calling strategies at multiple ctDNA levels
 * with bootstrap 95% CIs and train/test split
 */
const fs = require('fs');
const path = require('path');

const INPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'downsampled_data.json');
const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'variant_calling_results.json');
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

// ---- Statistical functions ----
function logFactorial(n) {
  if (n < 0) return -Infinity;
  if (n <= 1) return 0;
  // Stirling approximation for large n
  if (n > 20) return n * Math.log(n) - n + 0.5 * Math.log(2 * Math.PI * n);
  let s = 0;
  for (let i = 2; i <= n; i++) s += Math.log(i);
  return s;
}

function fisherExact(a, b, c, d) {
  // One-sided Fisher's exact test: probability of observing >= a given margins
  // Using hypergeometric distribution
  const n = a + b + c + d;
  const row1 = a + b;
  const col1 = a + c;
  let p = 0;
  const observedLogP = logHypergeometric(a, row1, col1, n);
  for (let k = a; k <= Math.min(row1, col1); k++) {
    p += Math.exp(logHypergeometric(k, row1, col1, n));
  }
  return Math.min(1, p);
}

function logHypergeometric(k, K, n, N) {
  return logComb(K, k) + logComb(N - K, n - k) - logComb(N, n);
}

function logComb(n, k) {
  if (k < 0 || k > n) return -Infinity;
  if (k === 0 || k === n) return 0;
  return logFactorial(n) - logFactorial(k) - logFactorial(n - k);
}

function betaRegularized(x, a, b) {
  // Use continued fraction for incomplete beta
  if (x < 0 || x > 1) return x < 0 ? 0 : 1;
  const tiny = 1e-30;
  const maxIter = 200;
  let qab = a + b;
  let qap = a + 1;
  let qam = a - 1;
  let c = 1;
  let d = 1 - qab * x / qap;
  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= maxIter; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 1e-12) break;
  }
  const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  return bt * h / a;
}

function logGamma(x) {
  if (x < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * x)) - logGamma(1 - x);
  x -= 1;
  const g = 7;
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  let z = x - 1 + g + 0.5;
  let s = c[0];
  for (let i = 1; i < g + 2; i++) s += c[i] / (x + i);
  return (Math.log(Math.sqrt(2 * Math.PI)) + Math.log(s) + (x + 0.5) * Math.log(z) - z);
}

// ---- Variant Calling Strategies ----
class VariantCaller {
  constructor(strategy, params = {}) {
    this.strategy = strategy;
    this.params = params;
    this.threshold = null;
  }

  score(obs) {
    const { mutant_reads, total_reads } = obs;
    const depth = total_reads || 50000;
    const vaf = mutant_reads / depth;
    const errorRate = 0.0001;
    const errorMult = obs.error_multiplier || 1;
    const effectiveErr = errorRate * errorMult;

    switch (this.strategy) {
      case 'vaf_threshold':
        return vaf;

      case 'fisher_exact': {
        // One-sided Fisher: mutant reads vs expected from error
        const expectedMut = Math.round(depth * effectiveErr);
        const nMut = mutant_reads;
        const nBg = depth - nMut;
        const eMut = Math.round(depth * effectiveErr);
        const eBg = depth - eMut;
        // Return -log10(p) so higher = more significant
        const p = fisherExact(nMut, nBg, eMut, eBg);
        return p < 1e-300 ? 300 : -Math.log10(Math.max(p, 1e-300));
      }

      case 'bayesian': {
        // Beta-Binomial with PoN prior
        // Prior: Beta(alpha_prior, beta_prior) where mean = error_rate
        const priorStrength = this.params.bayesian_prior_strength || 1000;
        const alphaPrior = effectiveErr * priorStrength;
        const betaPrior = (1 - effectiveErr) * priorStrength;
        // Posterior: Beta(alphaPrior + mutant_reads, betaPrior + depth - mutant_reads)
        const alphaPost = alphaPrior + mutant_reads;
        const betaPost = betaPrior + depth - mutant_reads;
        // Score: posterior probability that VAF > 3 * error_rate
        const threshold = 3 * errorRate;
        return 1 - betaRegularized(threshold, alphaPost, betaPost);
      }

      case 'likelihood_ratio': {
        // Poisson likelihood ratio: H1 (cancer λ_c) vs H0 (normal λ_n)
        const lambdaCancer = mutant_reads + 1; // MLE under H1
        const lambdaNormal = Math.max(depth * effectiveErr, 1); // under H0
        // log(LR) = mutant_reads * log(λ_c/λ_n) - (λ_c - λ_n) * depth
        // Actually for Poisson: log L = mutant*log(λ) - λ*depth
        // log LR = mutant*log(λ_c/λ_n) - depth*(λ_c - λ_n)
        const logLR = mutant_reads * Math.log(lambdaCancer / lambdaNormal) - depth * (lambdaCancer - lambdaNormal);
        return logLR > 0 ? logLR : 0;
      }

      default:
        return vaf;
    }
  }

  calibrateThreshold(calibrationScores, calibrationLabels, targetF1 = null) {
    // Find optimal threshold by maximizing F1 on calibration set
    const sorted = calibrationScores.map((s, i) => ({ score: s, label: calibrationLabels[i] }))
      .sort((a, b) => a.score - b.score);

    let bestThreshold = 0;
    let bestF1 = 0;
    let bestMetrics = null;
    const total = sorted.length;
    const totalPos = calibrationLabels.filter(l => l === 1).length;

    // Try each unique score as threshold
    const uniqueScores = [...new Set(sorted.map(s => s.score))].sort((a, b) => a - b);
    const nTest = Math.min(200, uniqueScores.length);
    const step = Math.max(1, Math.floor(uniqueScores.length / nTest));

    for (let i = 0; i < uniqueScores.length; i += step) {
      const thresh = uniqueScores[i];
      const tp = sorted.filter(s => s.score >= thresh && s.label === 1).length;
      const fp = sorted.filter(s => s.score >= thresh && s.label === 0).length;
      const fn = totalPos - tp;
      const tn = (total - totalPos) - fp;

      const precision = tp / Math.max(1, tp + fp);
      const recall = tp / Math.max(1, tp + fn);
      const f1 = (precision + recall > 0) ? 2 * precision * recall / (precision + recall) : 0;

      if (f1 > bestF1) {
        bestF1 = f1;
        bestThreshold = thresh;
        bestMetrics = { tp, fp, fn, tn, precision, recall, f1, threshold: thresh };
      }
    }

    this.threshold = bestThreshold;
    return bestMetrics;
  }

  predict(scores) {
    return scores.map(s => s >= this.threshold ? 1 : 0);
  }

  name() {
    const names = {
      'vaf_threshold': 'VAF Threshold',
      'fisher_exact': "Fisher's Exact",
      'bayesian': 'Bayesian Beta-Binomial',
      'likelihood_ratio': 'Likelihood Ratio'
    };
    return names[this.strategy] || this.strategy;
  }
}

// ---- Metrics ----
function computeMetrics(yTrue, yPred) {
  const n = yTrue.length;
  let tp = 0, fp = 0, tn = 0, fn = 0;
  for (let i = 0; i < n; i++) {
    if (yTrue[i] === 1 && yPred[i] === 1) tp++;
    if (yTrue[i] === 0 && yPred[i] === 1) fp++;
    if (yTrue[i] === 0 && yPred[i] === 0) tn++;
    if (yTrue[i] === 1 && yPred[i] === 0) fn++;
  }
  const sensitivity = tp / Math.max(1, tp + fn);
  const specificity = tn / Math.max(1, tn + fp);
  const precision = tp / Math.max(1, tp + fp);
  const f1 = (precision + sensitivity > 0) ? 2 * precision * sensitivity / (precision + sensitivity) : 0;
  const accuracy = (tp + tn) / n;
  return { tp, fp, tn, fn, sensitivity, specificity, precision, f1, accuracy, n };
}

function computeAUC(scores, labels) {
  const pairs = scores.map((s, i) => ({ score: s, label: labels[i] }));
  const pos = pairs.filter(p => p.label === 1);
  const neg = pairs.filter(p => p.label === 0);
  if (pos.length === 0 || neg.length === 0) return 0.5;

  let auc = 0;
  for (const p of pos) {
    for (const n of neg) {
      if (p.score > n.score) auc += 1;
      else if (p.score === n.score) auc += 0.5;
    }
  }
  return auc / (pos.length * neg.length);
}

function bootstrapCI(data, labels, metricFn, nBootstrap, rng) {
  const estimates = [];
  const posIdxs = data.map((_, i) => i).filter(i => labels[i] === 1);
  const negIdxs = data.map((_, i) => i).filter(i => labels[i] === 0);

  for (let b = 0; b < nBootstrap; b++) {
    // Stratified bootstrap
    const samplePos = [];
    const sampleNeg = [];
    for (let i = 0; i < posIdxs.length; i++) samplePos.push(posIdxs[Math.floor(rng() * posIdxs.length)]);
    for (let i = 0; i < negIdxs.length; i++) sampleNeg.push(negIdxs[Math.floor(rng() * negIdxs.length)]);
    const bootstrapIdxs = [...samplePos, ...sampleNeg];
    const bootData = bootstrapIdxs.map(i => data[i]);
    const bootLabels = bootstrapIdxs.map(i => labels[i]);
    estimates.push(metricFn(bootData, bootLabels));
  }

  estimates.sort((a, b) => a - b);
  const alpha = 0.025;
  const lo = estimates[Math.floor(alpha * estimates.length)];
  const hi = estimates[Math.ceil((1 - alpha) * estimates.length) - 1];
  const mean = estimates.reduce((a, b) => a + b, 0) / estimates.length;
  const sd = Math.sqrt(estimates.reduce((s, e) => s + (e - mean) ** 2, 0) / (estimates.length - 1));
  return { mean, sd, ci95_low: lo, ci95_high: hi };
}

// Helper: compute metrics from data+labels directly
function directMetrics(data, labels) {
  return computeMetrics(labels, labels.map((_, i) => data[i] >= 0.5 ? 1 : 0));
}

function directSensitivity(data, labels) {
  const m = computeMetrics(labels, data.map((_, i) => 0));
  return 0;
}

function sensitivityFrom(data, labels) {
  const tp = data.filter((d, i) => d >= 0.5 && labels[i] === 1).length;
  const totalPos = labels.filter(l => l === 1).length;
  return tp / Math.max(1, totalPos);
}

function specificityFrom(data, labels) {
  const tn = data.filter((d, i) => d < 0.5 && labels[i] === 0).length;
  const totalNeg = labels.filter(l => l === 0).length;
  return tn / Math.max(1, totalNeg);
}

// ---- Main ----
(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH NODE.JS VALIDATION - STEP 3: Variant Caller Validation');
  console.log('='.repeat(70));

  const input = JSON.parse(fs.readFileSync(INPUT_PATH, 'utf8'));
  const rng = createRNG(SEED + 1000);

  const strategies = ['vaf_threshold', 'fisher_exact', 'bayesian', 'likelihood_ratio'];
  const fractions = input.parameters.ctdna_fractions;
  const allResults = {};

  fractions.forEach(ctdnaFrac => {
    const key = `ctdna_${ctdnaFrac}`;
    const label = `ctDNA ${(ctdnaFrac*100).toFixed(ctdnaFrac < 0.0001 ? 4 : 3)}%`;
    console.log(`\n🔬 Evaluating at ${label}...`);

    const observations = input.observations[key];
    const variantObs = observations.filter(o => o.site_type === 'variant');
    const bgObs = observations.filter(o => o.site_type === 'background');

    // True variant = positive class, background = negative class
    const allObs = [...variantObs, ...bgObs];
    const allLabels = [...variantObs.map(() => 1), ...bgObs.map(() => 0)];

    // Shuffle for train/test split (stratified by cancer type)
    const variantByType = {};
    variantObs.forEach((o, i) => {
      const ct = o.cancer_type;
      if (!variantByType[ct]) variantByType[ct] = [];
      variantByType[ct].push(i);
    });
    const bgByType = {};
    bgObs.forEach((o, i) => {
      const ct = o.cancer_type;
      if (!bgByType[ct]) bgByType[ct] = [];
      bgByType[ct].push(i);
    });

    // 60% train, 40% test — stratified
    const trainVariantIdxs = [];
    const testVariantIdxs = [];
    Object.values(variantByType).forEach(idxs => {
      const shuffled = [...idxs].sort(() => rng() - 0.5);
      const split = Math.floor(shuffled.length * 0.6);
      trainVariantIdxs.push(...shuffled.slice(0, split));
      testVariantIdxs.push(...shuffled.slice(split));
    });

    const trainBgIdxs = [];
    const testBgIdxs = [];
    Object.values(bgByType).forEach(idxs => {
      const shuffled = [...idxs].sort(() => rng() - 0.5);
      const split = Math.floor(shuffled.length * 0.6);
      trainBgIdxs.push(...shuffled.slice(0, split));
      testBgIdxs.push(...shuffled.slice(split));
    });

    // Build train/test sets
    const trainObs = [...trainVariantIdxs.map(i => variantObs[i]), ...trainBgIdxs.map(i => bgObs[i])];
    const trainLabels = [...trainVariantIdxs.map(() => 1), ...trainBgIdxs.map(() => 0)];
    const testObs = [...testVariantIdxs.map(i => variantObs[i]), ...testBgIdxs.map(i => bgObs[i])];
    const testLabels = [...testVariantIdxs.map(() => 1), ...testBgIdxs.map(() => 0)];

    console.log(`   Train: ${trainObs.length} (${trainVariantIdxs.length} variant, ${trainBgIdxs.length} bg)`);
    console.log(`   Test:  ${testObs.length} (${testVariantIdxs.length} variant, ${testBgIdxs.length} bg)`);

    const fractionResults = {};

    strategies.forEach(strategy => {
      console.log(`\n   📊 Strategy: ${strategy}...`);
      const caller = new VariantCaller(strategy);

      // Score all
      const trainScores = trainObs.map(o => caller.score(o));
      const testScores = testObs.map(o => caller.score(o));

      // Calibrate threshold on training set
      const calibMetrics = caller.calibrateThreshold(trainScores, trainLabels);
      console.log(`      Calibrated threshold: ${caller.threshold.toFixed(6)}`);
      console.log(`      Calibration F1: ${calibMetrics.f1.toFixed(4)}`);

      // Predict on test set
      const testPredictions = caller.predict(testScores);
      const testMetrics = computeMetrics(testLabels, testPredictions);
      const testAUC = computeAUC(testScores, testLabels);

      // Bootstrap CIs on test set
      const rngCI = createRNG(SEED + 2000 + strategies.indexOf(strategy));

      const sensCI = bootstrapCI(testScores, testLabels,
        (scores, labels) => {
          const thresh = caller.threshold;
          const preds = scores.map(s => s >= thresh ? 1 : 0);
          const m = computeMetrics(labels, preds);
          return m.sensitivity;
        }, N_BOOTSTRAP, rngCI);

      const specCI = bootstrapCI(testScores, testLabels,
        (scores, labels) => {
          const thresh = caller.threshold;
          const preds = scores.map(s => s >= thresh ? 1 : 0);
          const m = computeMetrics(labels, preds);
          return m.specificity;
        }, N_BOOTSTRAP, rngCI);

      const precCI = bootstrapCI(testScores, testLabels,
        (scores, labels) => {
          const thresh = caller.threshold;
          const preds = scores.map(s => s >= thresh ? 1 : 0);
          const m = computeMetrics(labels, preds);
          return m.precision;
        }, N_BOOTSTRAP, rngCI);

      const f1CI = bootstrapCI(testScores, testLabels,
        (scores, labels) => {
          const thresh = caller.threshold;
          const preds = scores.map(s => s >= thresh ? 1 : 0);
          const m = computeMetrics(labels, preds);
          return m.f1;
        }, N_BOOTSTRAP, rngCI);

      const aucCI = bootstrapCI(testScores, testLabels,
        (scores, labels) => computeAUC(scores, labels), N_BOOTSTRAP, rngCI);

      fractionResults[strategy] = {
        strategy_name: caller.name(),
        calibration_threshold: caller.threshold,
        calibration_metrics: calibMetrics,
        test_metrics: {
          sensitivity: { value: testMetrics.sensitivity, ...sensCI },
          specificity: { value: testMetrics.specificity, ...specCI },
          precision: { value: testMetrics.precision, ...precCI },
          f1_score: { value: testMetrics.f1, ...f1CI },
          auc: { value: testAUC, ...aucCI },
          accuracy: testMetrics.accuracy,
          n_test: testObs.length,
          n_variant_test: testVariantIdxs.length,
          n_bg_test: testBgIdxs.length,
          tp: testMetrics.tp, fp: testMetrics.fp,
          tn: testMetrics.tn, fn: testMetrics.fn
        }
      };

      console.log(`      Test Sensitivity: ${(testMetrics.sensitivity*100).toFixed(2)}% [${(sensCI.ci95_low*100).toFixed(2)}-${(sensCI.ci95_high*100).toFixed(2)}]`);
      console.log(`      Test Specificity: ${(testMetrics.specificity*100).toFixed(2)}% [${(specCI.ci95_low*100).toFixed(2)}-${(specCI.ci95_high*100).toFixed(2)}]`);
      console.log(`      Test F1: ${testMetrics.f1.toFixed(4)} [${f1CI.ci95_low.toFixed(4)}-${f1CI.ci95_high.toFixed(4)}]`);
      console.log(`      Test AUC: ${testAUC.toFixed(4)} [${aucCI.ci95_low.toFixed(4)}-${aucCI.ci95_high.toFixed(4)}]`);
    });

    allResults[key] = { label, ctdna_fraction: ctdnaFrac, strategies: fractionResults };
  });

  // Summary across strategies
  console.log(`\n\n📊 FINAL SUMMARY: Best strategy per ctDNA level`);
  const summaryTable = [];
  fractions.forEach(ctdnaFrac => {
    const key = `ctdna_${ctdnaFrac}`;
    const row = { ctdna_fraction: ctdnaFrac, label: allResults[key].label };
    strategies.forEach(s => {
      const r = allResults[key].strategies[s];
      if (r) {
        row[`${s}_sens`] = r.test_metrics.sensitivity.value;
        row[`${s}_spec`] = r.test_metrics.specificity.value;
        row[`${s}_f1`] = r.test_metrics.f1_score.value;
        row[`${s}_auc`] = r.test_metrics.auc.value;
      }
    });

    // Find best F1
    let bestF1 = -1, bestStrat = '';
    strategies.forEach(s => {
      if (row[`${s}_f1`] > bestF1) { bestF1 = row[`${s}_f1`]; bestStrat = s; }
    });
    row.best_strategy = bestStrat;
    row.best_f1 = bestF1;
    summaryTable.push(row);
  });

  console.log('\nctDNA Level | Best Strategy | Best F1 | VAF Sens | Fisher Sens | Bayesian Sens | LR Sens');
  console.log('-'.repeat(95));
  summaryTable.forEach(r => {
    console.log(`${r.label.padEnd(12)} | ${r.best_strategy.padEnd(16)} | ${r.best_f1.toFixed(4).padEnd(8)} | ${(r.vaf_threshold_sens||0).toFixed(3).padEnd(8)} | ${(r.fisher_exact_sens||0).toFixed(3).padEnd(10)} | ${(r.bayesian_sens||0).toFixed(3).padEnd(12)} | ${(r.likelihood_ratio_sens||0).toFixed(3)}`);
  });

  const output = {
    metadata: {
      validation_type: 'variant_calling_real_data',
      n_bootstrap: N_BOOTSTRAP,
      seed: SEED,
      sequencing_depth: input.parameters.sequencing_depth,
      error_rate: input.parameters.error_rate,
      strategies_tested: strategies.map(s => (new VariantCaller(s)).name()),
      timestamp: new Date().toISOString()
    },
    summary: summaryTable,
    results: allResults
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved results to ${OUTPUT_PATH}`);
  console.log('\n✅ Step 3 complete.');
  console.log('='.repeat(70));
})();
