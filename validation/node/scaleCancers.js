#!/usr/bin/env node
/**
 * scaleCancers.js — FIX 2: Scale Cancer Types (8 → 20)
 * 
 * Problem: Only 8 cancer types. GRAIL covers 50+. Moldovan covers >10.
 * Adds 12 more cancer types with real TCGA mutation frequencies.
 * Runs performance-weighted multi-modal fusion on all 20 types.
 */

const fs = require('fs');
const path = require('path');

const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'scale_results.json');

const SEED = 42;
const N_BOOTSTRAP = 2000;
const N_SAMPLES_PER_TYPE = 100;  // 100 cancer + 100 healthy per type

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

function sigmoid(x) { return 1 / (1 + Math.exp(-Math.max(-50, Math.min(50, x)))); }

// ── 20 CANCER TYPES WITH REAL TCGA/COSMIC MUTATION FREQUENCIES ──
const ALL_CANCER_TYPES = [
  // Original 8
  {
    code: 'LUAD', name: 'Lung Adenocarcinoma',
    genes: [
      { gene: 'TP53', prevalence: 0.46 }, { gene: 'KRAS', prevalence: 0.33 },
      { gene: 'EGFR', prevalence: 0.14 }, { gene: 'STK11', prevalence: 0.17 },
      { gene: 'KEAP1', prevalence: 0.17 }, { gene: 'BRAF', prevalence: 0.07 },
      { gene: 'NF1', prevalence: 0.08 }, { gene: 'MET', prevalence: 0.07 },
    ],
    tmb: 8.7, shedding: { mean: 0.0032, cv: 1.1 },
    tcga_samples: 566,
  },
  {
    code: 'COADREAD', name: 'Colorectal Adenocarcinoma',
    genes: [
      { gene: 'APC', prevalence: 0.81 }, { gene: 'TP53', prevalence: 0.59 },
      { gene: 'KRAS', prevalence: 0.44 }, { gene: 'PIK3CA', prevalence: 0.20 },
      { gene: 'FBXW7', prevalence: 0.12 }, { gene: 'SMAD4', prevalence: 0.10 },
      { gene: 'BRAF', prevalence: 0.12 }, { gene: 'NRAS', prevalence: 0.07 },
    ],
    tmb: 4.5, shedding: { mean: 0.008, cv: 0.9 },
    tcga_samples: 594,
  },
  {
    code: 'BRCA', name: 'Breast Invasive Carcinoma',
    genes: [
      { gene: 'TP53', prevalence: 0.37 }, { gene: 'PIK3CA', prevalence: 0.36 },
      { gene: 'GATA3', prevalence: 0.12 }, { gene: 'CDH1', prevalence: 0.11 },
      { gene: 'MAP3K1', prevalence: 0.09 }, { gene: 'PTEN', prevalence: 0.06 },
      { gene: 'AKT1', prevalence: 0.04 }, { gene: 'ERBB2', prevalence: 0.06 },
    ],
    tmb: 1.8, shedding: { mean: 0.0012, cv: 1.3 },
    tcga_samples: 1084,
  },
  {
    code: 'PRAD', name: 'Prostate Adenocarcinoma',
    genes: [
      { gene: 'SPOP', prevalence: 0.11 }, { gene: 'TP53', prevalence: 0.09 },
      { gene: 'FOXA1', prevalence: 0.07 }, { gene: 'PTEN', prevalence: 0.06 },
      { gene: 'AR', prevalence: 0.03 }, { gene: 'PIK3CA', prevalence: 0.04 },
      { gene: 'BRCA2', prevalence: 0.04 }, { gene: 'ATM', prevalence: 0.04 },
    ],
    tmb: 0.9, shedding: { mean: 0.0004, cv: 1.4 },
    tcga_samples: 494,
  },
  {
    code: 'STAD', name: 'Stomach Adenocarcinoma',
    genes: [
      { gene: 'TP53', prevalence: 0.49 }, { gene: 'ARID1A', prevalence: 0.22 },
      { gene: 'CDH1', prevalence: 0.16 }, { gene: 'PIK3CA', prevalence: 0.16 },
      { gene: 'KRAS', prevalence: 0.08 }, { gene: 'RHOA', prevalence: 0.11 },
      { gene: 'ERBB2', prevalence: 0.10 }, { gene: 'CTNNB1', prevalence: 0.07 },
    ],
    tmb: 3.3, shedding: { mean: 0.005, cv: 1.0 },
    tcga_samples: 441,
  },
  {
    code: 'LIHC', name: 'Liver Hepatocellular Carcinoma',
    genes: [
      { gene: 'CTNNB1', prevalence: 0.26 }, { gene: 'TP53', prevalence: 0.31 },
      { gene: 'ARID1A', prevalence: 0.10 }, { gene: 'ARID2', prevalence: 0.07 },
      { gene: 'AXIN1', prevalence: 0.08 }, { gene: 'TERT', prevalence: 0.44 },
      { gene: 'KEAP1', prevalence: 0.05 }, { gene: 'RB1', prevalence: 0.04 },
    ],
    tmb: 2.6, shedding: { mean: 0.006, cv: 1.0 },
    tcga_samples: 377,
  },
  {
    code: 'PAAD', name: 'Pancreatic Adenocarcinoma',
    genes: [
      { gene: 'KRAS', prevalence: 0.93 }, { gene: 'TP53', prevalence: 0.72 },
      { gene: 'SMAD4', prevalence: 0.32 }, { gene: 'CDKN2A', prevalence: 0.28 },
      { gene: 'ARID1A', prevalence: 0.08 }, { gene: 'RNF43', prevalence: 0.06 },
      { gene: 'GNAS', prevalence: 0.05 }, { gene: 'BRAF', prevalence: 0.03 },
    ],
    tmb: 2.5, shedding: { mean: 0.007, cv: 0.9 },
    tcga_samples: 185,
  },
  {
    code: 'OV', name: 'Ovarian Serous Cystadenocarcinoma',
    genes: [
      { gene: 'TP53', prevalence: 0.96 }, { gene: 'BRCA1', prevalence: 0.16 },
      { gene: 'BRCA2', prevalence: 0.12 }, { gene: 'NF1', prevalence: 0.10 },
      { gene: 'RB1', prevalence: 0.08 }, { gene: 'CDK12', prevalence: 0.05 },
      { gene: 'CSMD3', prevalence: 0.08 }, { gene: 'FAT3', prevalence: 0.07 },
    ],
    tmb: 2.5, shedding: { mean: 0.010, cv: 0.8 },
    tcga_samples: 489,
  },
  // ── NEW 12 ──
  {
    code: 'CESC', name: 'Cervical Squamous Cell Carcinoma',
    genes: [
      { gene: 'PIK3CA', prevalence: 0.31 }, { gene: 'KRAS', prevalence: 0.08 },
      { gene: 'TP53', prevalence: 0.14 }, { gene: 'EP300', prevalence: 0.15 },
      { gene: 'FBXW7', prevalence: 0.15 }, { gene: 'PTEN', prevalence: 0.09 },
      { gene: 'ARID1A', prevalence: 0.10 }, { gene: 'ERBB2', prevalence: 0.07 },
    ],
    tmb: 4.0, shedding: { mean: 0.004, cv: 1.1 },
    tcga_samples: 304, source: 'TCGA PanCancer (Burk 2017 Nature)',
  },
  {
    code: 'ESCA', name: 'Esophageal Carcinoma',
    genes: [
      { gene: 'TP53', prevalence: 0.83 }, { gene: 'CDKN2A', prevalence: 0.48 },
      { gene: 'NFE2L2', prevalence: 0.18 }, { gene: 'PIK3CA', prevalence: 0.15 },
      { gene: 'NOTCH1', prevalence: 0.12 }, { gene: 'KMT2D', prevalence: 0.14 },
      { gene: 'FAT1', prevalence: 0.10 }, { gene: 'ARID1A', prevalence: 0.08 },
    ],
    tmb: 5.5, shedding: { mean: 0.006, cv: 1.0 },
    tcga_samples: 185, source: 'TCGA PanCancer (Cancer Genome Atlas Network 2017 Nature)',
  },
  {
    code: 'KIRC', name: 'Kidney Renal Clear Cell Carcinoma',
    genes: [
      { gene: 'VHL', prevalence: 0.72 }, { gene: 'PBRM1', prevalence: 0.41 },
      { gene: 'SETD2', prevalence: 0.15 }, { gene: 'BAP1', prevalence: 0.14 },
      { gene: 'MTOR', prevalence: 0.06 }, { gene: 'PTEN', prevalence: 0.04 },
      { gene: 'KDM5C', prevalence: 0.07 }, { gene: 'TP53', prevalence: 0.05 },
    ],
    tmb: 1.3, shedding: { mean: 0.003, cv: 1.2 },
    tcga_samples: 534, source: 'TCGA KIRC (Creighton 2013 Nature)',
  },
  {
    code: 'LGG', name: 'Lower Grade Glioma',
    genes: [
      { gene: 'IDH1', prevalence: 0.75 }, { gene: 'TP53', prevalence: 0.49 },
      { gene: 'ATRX', prevalence: 0.33 }, { gene: 'CIC', prevalence: 0.19 },
      { gene: 'FUBP1', prevalence: 0.13 }, { gene: 'EGFR', prevalence: 0.10 },
      { gene: 'PTEN', prevalence: 0.06 }, { gene: 'PIK3CA', prevalence: 0.06 },
    ],
    tmb: 1.0, shedding: { mean: 0.0005, cv: 1.5 },
    tcga_samples: 515, source: 'TCGA LGG (Brat 2015 NEJM)',
  },
  {
    code: 'SKCM', name: 'Skin Cutaneous Melanoma',
    genes: [
      { gene: 'BRAF', prevalence: 0.52 }, { gene: 'NRAS', prevalence: 0.28 },
      { gene: 'NF1', prevalence: 0.14 }, { gene: 'TP53', prevalence: 0.16 },
      { gene: 'PTEN', prevalence: 0.12 }, { gene: 'KIT', prevalence: 0.06 },
      { gene: 'ARID2', prevalence: 0.07 }, { gene: 'CDKN2A', prevalence: 0.12 },
    ],
    tmb: 11.5, shedding: { mean: 0.008, cv: 0.9 },
    tcga_samples: 470, source: 'TCGA SKCM (Akbani 2015 Cell)',
  },
  {
    code: 'THCA', name: 'Thyroid Carcinoma',
    genes: [
      { gene: 'BRAF', prevalence: 0.60 }, { gene: 'NRAS', prevalence: 0.09 },
      { gene: 'HRAS', prevalence: 0.04 }, { gene: 'KRAS', prevalence: 0.03 },
      { gene: 'RET', prevalence: 0.07 }, { gene: 'TERT', prevalence: 0.15 },
      { gene: 'EIF1AX', prevalence: 0.02 }, { gene: 'PIK3CA', prevalence: 0.03 },
    ],
    tmb: 0.4, shedding: { mean: 0.001, cv: 1.3 },
    tcga_samples: 500, source: 'TCGA THCA (Agrawal 2014 Cell)',
  },
  {
    code: 'UCEC', name: 'Uterine Corpus Endometrial Carcinoma',
    genes: [
      { gene: 'PTEN', prevalence: 0.57 }, { gene: 'PIK3CA', prevalence: 0.42 },
      { gene: 'ARID1A', prevalence: 0.37 }, { gene: 'TP53', prevalence: 0.25 },
      { gene: 'KRAS', prevalence: 0.20 }, { gene: 'CTNNB1', prevalence: 0.22 },
      { gene: 'FBXW7', prevalence: 0.12 }, { gene: 'PPP2R1A', prevalence: 0.12 },
    ],
    tmb: 5.0, shedding: { mean: 0.005, cv: 1.0 },
    tcga_samples: 547, source: 'TCGA UCEC (Levine 2013 Nature)',
  },
  {
    code: 'GBM', name: 'Glioblastoma Multiforme',
    genes: [
      { gene: 'EGFR', prevalence: 0.57 }, { gene: 'PTEN', prevalence: 0.33 },
      { gene: 'TP53', prevalence: 0.28 }, { gene: 'NF1', prevalence: 0.18 },
      { gene: 'PIK3CA', prevalence: 0.12 }, { gene: 'RB1', prevalence: 0.11 },
      { gene: 'PDGFRA', prevalence: 0.10 }, { gene: 'IDH1', prevalence: 0.06 },
    ],
    tmb: 3.0, shedding: { mean: 0.0003, cv: 1.6 },
    tcga_samples: 396, source: 'TCGA GBM (Brennan 2013 Cell)',
  },
  {
    code: 'AML', name: 'Acute Myeloid Leukemia',
    genes: [
      { gene: 'NPM1', prevalence: 0.27 }, { gene: 'FLT3', prevalence: 0.28 },
      { gene: 'DNMT3A', prevalence: 0.23 }, { gene: 'IDH2', prevalence: 0.10 },
      { gene: 'IDH1', prevalence: 0.08 }, { gene: 'TET2', prevalence: 0.10 },
      { gene: 'RUNX1', prevalence: 0.09 }, { gene: 'CEBPA', prevalence: 0.07 },
    ],
    tmb: 1.0, shedding: { mean: 0.040, cv: 0.5 },  // Liquid tumor — high shedding
    tcga_samples: 200, source: 'TCGA LAML (Ley 2013 NEJM)',
  },
  {
    code: 'DLBC', name: 'Diffuse Large B-Cell Lymphoma',
    genes: [
      { gene: 'MYD88', prevalence: 0.29 }, { gene: 'CD79B', prevalence: 0.21 },
      { gene: 'EZH2', prevalence: 0.22 }, { gene: 'KMT2D', prevalence: 0.32 },
      { gene: 'CREBBP', prevalence: 0.28 }, { gene: 'BCL2', prevalence: 0.24 },
      { gene: 'TP53', prevalence: 0.20 }, { gene: 'CARD11', prevalence: 0.10 },
    ],
    tmb: 2.5, shedding: { mean: 0.035, cv: 0.6 },
    tcga_samples: 48, source: 'TCGA DLBC (Schmitz 2018 NEJM)',
  },
  {
    code: 'SARC', name: 'Sarcoma',
    genes: [
      { gene: 'TP53', prevalence: 0.40 }, { gene: 'RB1', prevalence: 0.12 },
      { gene: 'ATRX', prevalence: 0.18 }, { gene: 'CDKN2A', prevalence: 0.10 },
      { gene: 'NF1', prevalence: 0.06 }, { gene: 'PIK3CA', prevalence: 0.05 },
      { gene: 'KMT2C', prevalence: 0.07 }, { gene: 'PTEN', prevalence: 0.04 },
    ],
    tmb: 1.5, shedding: { mean: 0.002, cv: 1.2 },
    tcga_samples: 261, source: 'TCGA SARC (Abeshouse 2017 Cell)',
  },
  {
    code: 'MESO', name: 'Mesothelioma',
    genes: [
      { gene: 'BAP1', prevalence: 0.57 }, { gene: 'NF2', prevalence: 0.39 },
      { gene: 'CDKN2A', prevalence: 0.38 }, { gene: 'TP53', prevalence: 0.15 },
      { gene: 'SETD2', prevalence: 0.09 }, { gene: 'LATS2', prevalence: 0.10 },
      { gene: 'TP53BP1', prevalence: 0.06 }, { gene: 'SETDB1', prevalence: 0.06 },
    ],
    tmb: 1.3, shedding: { mean: 0.003, cv: 1.1 },
    tcga_samples: 87, source: 'TCGA MESO (Bueno 2016 Nat Genet)',
  },
];

// ── Multi-modal Feature Generation ──
// Generates realistic features for each cancer type: variant, methylation, fragment, expression

function generateFeaturesForType(cancerType, nCancer, nHealthy, ctDNAFraction, rng) {
  const samples = [];
  
  // ── Modality 1: Variant Calling Score ──
  function generateVariantScore(isCancer, ctIdx) {
    const depth = 50000;
    const errorRate = 0.0001;
    
    if (!isCancer) {
      // Healthy: background noise + CHIP
      const chipProb = 0.10; // 10% have CHIP
      const hasChip = rng() < chipProb;
      const noiseVAF = hasChip ? 0.005 + rng() * 0.02 : rng() * 0.002;
      const noiseReads = depth * noiseVAF * (0.5 + rng());
      const detection = noiseReads > 2 ? 1 : 0;
      return {
        nVariantsDetected: detection + (rng() < 0.02 ? 1 : 0),
        maxVAF: noiseVAF * (1 + normalRand(rng) * 0.3),
        meanVAF: noiseVAF * 0.8,
        variantScore: sigmoid(-5 + normalRand(rng) * 2 + detection * 3),
      };
    }
    
    // Cancer: real variants + ctDNA fraction
    const nTotalVariants = cancerType.genes.length;
    let nDetected = 0;
    let sumVAF = 0, maxVAF = 0;
    
    cancerType.genes.forEach(g => {
      const hasMutation = rng() < g.prevalence * 0.8; // 80% of prevalence detectable
      if (hasMutation) {
        const trueVAF = ctDNAFraction * (0.05 + rng() * 0.3);
        const observedVAF = trueVAF * (1 + normalRand(rng) * 0.2);
        const expectedReads = depth * observedVAF;
        
        // Poisson detection
        let poissonProb = 0;
        for (let k = 0; k < 3; k++) {
          poissonProb += Math.exp(-expectedReads) * Math.pow(expectedReads, k) / factorial(k);
        }
        const detected = rng() < (1 - poissonProb);
        
        if (detected || observedVAF > 0.001) {
          nDetected++;
          sumVAF += observedVAF;
          maxVAF = Math.max(maxVAF, observedVAF);
        }
      }
    });
    
    const variantScore = sigmoid(
      -3.0 + nDetected * 0.5 + Math.min(maxVAF / ctDNAFraction, 10) * 2.0 + normalRand(rng) * 1.5
    );
    
    return {
      nVariantsDetected: nDetected,
      maxVAF,
      meanVAF: nDetected > 0 ? sumVAF / nDetected : 0,
      variantScore,
    };
  }
  
  // ── Modality 2: Methylation Entropy Score ──
  function generateMethylationScore(isCancer, ctIdx) {
    if (!isCancer) {
      const score = sigmoid(-2 + normalRand(rng) * 1.5);
      return { entropy: 0.35 + rng() * 0.3, methylationScore: score };
    }
    const tmbFactor = cancerType.tmb / 10;
    const entropy = 0.55 + rng() * 0.35 + tmbFactor * 0.1;
    const score = sigmoid(
      -1.5 + entropy * 3.0 + Math.log(1 + ctDNAFraction * 100) * 0.5 + normalRand(rng) * 1.0
    );
    return { entropy, methylationScore: score };
  }
  
  // ── Modality 3: Fragment Size Score ──
  function generateFragmentScore(isCancer, ctIdx) {
    if (!isCancer) {
      const score = sigmoid(-1 + normalRand(rng) * 2);
      return { fragRatio: 0.15 + rng() * 0.15, fragmentScore: score };
    }
    const fragRatio = 0.25 + rng() * 0.35;
    const score = sigmoid(
      -2.0 + fragRatio * 4.0 + Math.log(1 + ctDNAFraction * 50) * 0.8 + normalRand(rng) * 1.2
    );
    return { fragRatio, fragmentScore: score };
  }
  
  // ── Modality 4: CNV Score ──
  function generateCNVScore(isCancer, ctIdx) {
    if (!isCancer) {
      const score = sigmoid(-0.5 + normalRand(rng) * 1.8);
      return { cnvEvents: Math.floor(rng() * 3), cnvScore: score };
    }
    const tmbFactor = cancerType.tmb / 10;
    const cnvEvents = Math.floor(rng() * 5 + tmbFactor * 3);
    const score = sigmoid(
      -2.5 + cnvEvents * 0.6 + Math.log(1 + ctDNAFraction * 30) * 0.6 + normalRand(rng) * 1.3
    );
    return { cnvEvents, cnvScore: score };
  }
  
  // Generate samples
  for (let i = 0; i < nCancer; i++) {
    const v = generateVariantScore(true, ctDNAFraction);
    const m = generateMethylationScore(true, ctDNAFraction);
    const f = generateFragmentScore(true, ctDNAFraction);
    const c = generateCNVScore(true, ctDNAFraction);
    
    samples.push({
      sample_id: `${cancerType.code}_cancer_${i}`,
      cancer_type: cancerType.code,
      label: 1,
      ctDNA_fraction: ctDNAFraction,
      features: {
        variant: v,
        methylation: m,
        fragment: f,
        cnv: c,
      },
    });
  }
  
  for (let i = 0; i < nHealthy; i++) {
    const v = generateVariantScore(false, 0);
    const m = generateMethylationScore(false, 0);
    const f = generateFragmentScore(false, 0);
    const c = generateCNVScore(false, 0);
    
    samples.push({
      sample_id: `healthy_${i}`,
      cancer_type: null,
      label: 0,
      ctDNA_fraction: 0,
      features: {
        variant: v,
        methylation: m,
        fragment: f,
        cnv: c,
      },
    });
  }
  
  return samples;
}

function factorial(n) {
  if (n <= 1) return 1;
  let r = 1;
  for (let i = 2; i <= n; i++) r *= i;
  return r;
}

// ── AUC ──
function computeAUC(scores, labels) {
  const pairs = scores.map((s, i) => ({ s, l: labels[i] }));
  pairs.sort((a, b) => b.s - a.s);
  
  let auc = 0, prevFpr = 0, prevTpr = 0;
  let totalPos = 0, totalNeg = 0;
  labels.forEach(l => { if (l) totalPos++; else totalNeg++; });
  if (totalPos === 0 || totalNeg === 0) return 0.5;
  
  let tp = 0, fp = 0;
  for (let i = 0; i < pairs.length; i++) {
    if (pairs[i].l) tp++; else fp++;
    if (i === pairs.length - 1 || pairs[i].s !== pairs[i + 1]?.s) {
      const tpr = tp / totalPos;
      const fpr = fp / totalNeg;
      auc += (fpr - prevFpr) * (tpr + prevTpr) / 2;
      prevFpr = fpr;
      prevTpr = tpr;
    }
  }
  return auc;
}

// ── Performance-Weighted Fusion ──
// Multimodal scores are combined with weights proportional to each modality's standalone AUC
function performanceWeightedFusion(samples) {
  const rng = createRNG(SEED);
  
  // Split: 70% train, 30% test
  const indices = samples.map((_, i) => i);
  for (let i = indices.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    [indices[i], indices[j]] = [indices[j], indices[i]];
  }
  
  const nTrain = Math.floor(samples.length * 0.7);
  const trainIdx = new Set(indices.slice(0, nTrain));
  
  const trainSet = samples.filter((_, i) => trainIdx.has(i));
  const testSet = samples.filter((_, i) => !trainIdx.has(i));
  
  // Compute per-modality AUC on training set
  const modalities = ['variant', 'methylation', 'fragment', 'cnv'];
  const modalityAUCs = {};
  
  modalities.forEach(mod => {
    const scores = trainSet.map(s => s.features[mod][`${mod}Score`]);
    const labels = trainSet.map(s => s.label);
    modalityAUCs[mod] = computeAUC(scores, labels);
  });
  
  // Weight = (AUC - 0.5)^2 — emphasizes strong modalities, de-emphasizes weak
  const weights = {};
  let totalWeight = 0;
  modalities.forEach(mod => {
    weights[mod] = Math.pow(Math.max(0, modalityAUCs[mod] - 0.5), 2);
    totalWeight += weights[mod];
  });
  
  // Normalize
  modalities.forEach(mod => {
    weights[mod] = totalWeight > 0 ? weights[mod] / totalWeight : 0.25;
  });
  
  // Fuse on test set
  const fusedScores = testSet.map(s => {
    let score = 0;
    modalities.forEach(mod => {
      score += weights[mod] * s.features[mod][`${mod}Score`];
    });
    return score;
  });
  const testLabels = testSet.map(s => s.label);
  
  const fusionAUC = computeAUC(fusedScores, testLabels);
  
  // Bootstrap CI
  const bootAUCs = [];
  const bootRng = createRNG(SEED + 99);
  for (let b = 0; b < N_BOOTSTRAP; b++) {
    const idx = new Array(testSet.length).fill(0).map(() => Math.floor(bootRng() * testSet.length));
    const bs = idx.map(i => fusedScores[i]);
    const bl = idx.map(i => testLabels[i]);
    bootAUCs.push(computeAUC(bs, bl));
  }
  bootAUCs.sort((a, b) => a - b);
  
  return {
    fusionAUC,
    ci95_lo: bootAUCs[Math.floor(N_BOOTSTRAP * 0.025)],
    ci95_hi: bootAUCs[Math.floor(N_BOOTSTRAP * 0.975)],
    modalityAUCs,
    weights,
    nTest: testSet.length,
  };
}

// ── MAIN ──
console.log('='.repeat(70));
console.log('FIX 2: SCALING CANCER TYPES (8 → 20)');
console.log('='.repeat(70));
console.log(`\nTotal cancer types: ${ALL_CANCER_TYPES.length}`);
console.log(`Samples per type: ${N_SAMPLES_PER_TYPE} cancer + ${N_SAMPLES_PER_TYPE} healthy`);
console.log(`Total cancer samples: ${ALL_CANCER_TYPES.length * N_SAMPLES_PER_TYPE}`);
console.log(`Total healthy samples: ${ALL_CANCER_TYPES.length * N_SAMPLES_PER_TYPE}`);

const ctDNAFractions = [0.001, 0.0025, 0.005, 0.01];
const perTypeResults = {};
const overallFusionResults = {};

// Run per-cancer-type fusion
console.log('\n📊 Per-Cancer-Type Multi-Modal Fusion AUC');
console.log('-'.repeat(70));

let overallFusedScores = [];
let overallLabels = [];

ALL_CANCER_TYPES.forEach(ct => {
  console.log(`\n  ${ct.code} (${ct.name}):`);
  console.log(`    Genotype: ${ct.genes.map(g => g.gene).join(', ')}`);
  console.log(`    TMB: ${ct.tmb}/Mb | ctDNA Shedding: ${(ct.shedding.mean*100).toFixed(2)}% (CV: ${ct.shedding.cv})`);
  console.log(`    TCGA samples: ${ct.tcga_samples} | Source: ${ct.source || 'COSMIC v99 + TCGA PanCancer'}`);
  
  perTypeResults[ct.code] = {};
  
  ctDNAFractions.forEach(ctdna => {
    // Re-seed per type+ctDNA combo
    const typeRng = createRNG(SEED + ALL_CANCER_TYPES.indexOf(ct) * 100 + Math.round(ctdna * 10000));
    const samples = generateFeaturesForType(ct, N_SAMPLES_PER_TYPE, N_SAMPLES_PER_TYPE, ctdna, typeRng);
    const result = performanceWeightedFusion(samples);
    
    perTypeResults[ct.code][ctdna] = result;
    
    console.log(`    ctDNA ${(ctdna*100).toFixed(1)}%: AUC = ${result.fusionAUC.toFixed(4)} [${result.ci95_lo.toFixed(4)}, ${result.ci95_hi.toFixed(4)}]`);
    console.log(`      Modality weights: V=${result.weights.variant.toFixed(2)} M=${result.weights.methylation.toFixed(2)} F=${result.weights.fragment.toFixed(2)} C=${result.weights.cnv.toFixed(2)}`);
    
    // Collect for overall fusion
    const testSet = samples;
    const modalities = ['variant', 'methylation', 'fragment', 'cnv'];
    const fusedScore = testSet.map(s => {
      let score = 0;
      modalities.forEach(mod => {
        const w = result.weights[mod] || 0.25;
        score += w * s.features[mod][`${mod}Score`];
      });
      return { score, label: s.label };
    });
    fusedScore.forEach(fs => {
      overallFusedScores.push(fs.score);
      overallLabels.push(fs.label);
    });
  });
});

// Overall performance across all 20 types
console.log('\n📊 OVERALL PERFORMANCE (All 20 Cancer Types Combined)');
console.log('-'.repeat(70));

const overallAUC = computeAUC(overallFusedScores, overallLabels);
const bootRng = createRNG(SEED + 777);
const overallBoot = [];
for (let b = 0; b < N_BOOTSTRAP; b++) {
  const idx = new Array(overallFusedScores.length).fill(0).map(() => Math.floor(bootRng() * overallFusedScores.length));
  const bs = idx.map(i => overallFusedScores[i]);
  const bl = idx.map(i => overallLabels[i]);
  overallBoot.push(computeAUC(bs, bl));
}
overallBoot.sort((a, b) => a - b);

console.log(`  Overall AUC: ${overallAUC.toFixed(4)} [${overallBoot[Math.floor(N_BOOTSTRAP*0.025)].toFixed(4)}, ${overallBoot[Math.floor(N_BOOTSTRAP*0.975)].toFixed(4)}]`);
console.log(`  Samples evaluated: ${overallFusedScores.length}`);
console.log(`  Cancer types: ${ALL_CANCER_TYPES.length}`);

// Per-type summary table
console.log('\n📊 Per-Cancer-Type AUC Summary (at 0.5% ctDNA)');
console.log('-'.repeat(70));
console.log('  Type       Name                              AUC      CI95              TMB');
console.log('  ----       ----                              ---      ----              ---');

ALL_CANCER_TYPES.forEach(ct => {
  const r = perTypeResults[ct.code]?.[0.005];
  if (r) {
    console.log(`  ${ct.code.padEnd(10)} ${ct.name.padEnd(34)} ${r.fusionAUC.toFixed(4)}  [${r.ci95_lo.toFixed(3)},${r.ci95_hi.toFixed(3)}]  ${ct.tmb.toFixed(1)}`);
  }
});

// Best/worst performing types
console.log('\n📊 Best & Worst Performing Cancers (at 0.5% ctDNA)');
console.log('-'.repeat(70));

const ranked = ALL_CANCER_TYPES
  .map(ct => ({
    code: ct.code,
    name: ct.name,
    auc: perTypeResults[ct.code]?.[0.005]?.fusionAUC || 0.5,
    shedding: ct.shedding.mean,
    tmb: ct.tmb,
  }))
  .sort((a, b) => b.auc - a.auc);

console.log('\n  TOP 5:');
ranked.slice(0, 5).forEach((r, i) => {
  console.log(`    ${i+1}. ${r.code} — AUC=${r.auc.toFixed(4)}, shedding=${(r.shedding*100).toFixed(2)}%, TMB=${r.tmb}`);
});

console.log('\n  BOTTOM 5:');
ranked.slice(-5).reverse().forEach((r, i) => {
  console.log(`    ${20-4+i}. ${r.code} — AUC=${r.auc.toFixed(4)}, shedding=${(r.shedding*100).toFixed(2)}%, TMB=${r.tmb}`);
});

// Compare to Grail Galleri
console.log('\n📊 Comparison to Published MCED Assays');
console.log('-'.repeat(70));
console.log(`  DeepCatch (20 types, simulation): AUC=${overallAUC.toFixed(4)}`);
console.log(`  Grail Galleri (50+ types, CLINICAL): Sensitivity=51.5% @ 99.5% Specificity`);
console.log(`  CancerSEEK (8 types, CLINICAL): Sensitivity=70.0% @ 99.0% Specificity`);
console.log(`  ⚠️ Note: DeepCatch is simulation-only. Clinical comparison not meaningful.`);

// ── Output ──
const output = {
  generated: new Date().toISOString(),
  total_cancer_types: ALL_CANCER_TYPES.length,
  cancer_types: ALL_CANCER_TYPES.map(ct => ({
    code: ct.code,
    name: ct.name,
    tcga_samples: ct.tcga_samples,
    tmb: ct.tmb,
    shedding_mean: ct.shedding.mean,
    shedding_cv: ct.shedding.cv,
    top_genes: ct.genes.slice(0, 5).map(g => `${g.gene}(${(g.prevalence*100).toFixed(0)}%)`).join(', '),
    source: ct.source || 'COSMIC v99 + TCGA PanCancer',
  })),
  samples_per_type: N_SAMPLES_PER_TYPE,
  ctDNA_fractions_tested: ctDNAFractions,
  per_type_results: perTypeResults,
  overall: {
    auc: overallAUC,
    ci95_lo: overallBoot[Math.floor(N_BOOTSTRAP * 0.025)],
    ci95_hi: overallBoot[Math.floor(N_BOOTSTRAP * 0.975)],
    n_samples: overallFusedScores.length,
    best_types: ranked.slice(0, 5).map(r => ({ code: r.code, name: r.name, auc: r.auc })),
    worst_types: ranked.slice(-5).reverse().map(r => ({ code: r.code, name: r.name, auc: r.auc })),
  },
  comparison: {
    deepcatch_20_types: { auc: overallAUC, note: 'SIMULATION ONLY' },
    grail_galleri: { sensitivity: 51.5, specificity: 99.5, cancer_types: 50, note: 'CLINICAL VALIDATION' },
    cancerSEEK: { sensitivity: 70.0, specificity: 99.0, cancer_types: 8, note: 'CLINICAL VALIDATION' },
  },
};

fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
console.log(`\n✅ Scale results written to ${OUTPUT_PATH}`);
