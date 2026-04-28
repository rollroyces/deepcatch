#!/usr/bin/env node
/**
 * downsample.js - cfDNA Downsampling Engine
 * Simulates ctDNA dilution to ultra-low fractions with Poisson sampling noise
 */
const fs = require('fs');
const path = require('path');

const DATASET_PATH = path.join(__dirname, '..', 'tcga', 'tcga_cache', 'fallback_dataset.json');
const OUTPUT_PATH = path.join(__dirname, '..', '..', 'results', 'node', 'downsampled_data.json');

const SEQUENCING_DEPTH = 50000;
const ERROR_RATE = 0.0001;
const CTDNA_FRACTIONS = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]; // 1 down to 0.001%
const SEED = 42;

// Seeded RNG (xoshiro128**)
function createRNG(seed) {
  let s0 = seed | 0;
  let s1 = (seed * 1812433253 + 1) | 0;
  let s2 = (seed * 1812433253 + 2) | 0;
  let s3 = (seed * 1812433253 + 3) | 0;
  function rotl(x, k) { return ((x << k) | (x >>> (32 - k))) | 0; }
  return function () {
    const result = ((rotl((s1 * 5) | 0, 7) * 9) | 0) >>> 0;
    const t = (s1 << 9) | 0;
    s2 ^= s0;
    s3 ^= s1;
    s1 ^= s2;
    s0 ^= s3;
    s2 ^= t;
    s3 = rotl(s3, 11);
    return result / 4294967296; // [0, 1)
  };
}

// Poisson random variate (Knuth algorithm)
function poisson(lambda, rng) {
  if (lambda < 30) {
    // Knuth's method for small lambda
    const L = Math.exp(-lambda);
    let k = 0;
    let p = 1;
    do {
      k++;
      p *= rng();
    } while (p > L);
    return k - 1;
  } else {
    // Normal approximation for large lambda
    const x = normalRand(rng) * Math.sqrt(lambda) + lambda;
    return Math.max(0, Math.round(x));
  }
}

function normalRand(rng) {
  // Box-Muller
  let u1, u2;
  do { u1 = rng(); } while (u1 === 0);
  u2 = rng();
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

// Trinucleotide context error multipliers (10x range)
const TRINUC_ERROR_MULTIPLIERS = {
  'C_G': 5.0,    // CpG sites, highest error
  'T_G': 3.0,
  'C_A': 2.5,
  'G_C': 2.0,
  'A_T': 1.5,
  'default': 1.0,
  'T_A': 0.8,
  'G_T': 0.7,
  'A_C': 0.6,
  'C_T': 0.5   // lowest error
};

function getTrinucContext(pos) {
  // Deterministic mapping based on position
  const contexts = Object.keys(TRINUC_ERROR_MULTIPLIERS).filter(k => k !== 'default');
  return contexts[pos % contexts.length];
}

(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH NODE.JS VALIDATION - STEP 2: cfDNA Downsampling');
  console.log('='.repeat(70));

  const raw = JSON.parse(fs.readFileSync(DATASET_PATH, 'utf8'));
  const variants = raw.ground_truth_variants;
  const samples = raw.sample_metadata;

  const rng = createRNG(SEED);

  console.log(`\n⚙️  Parameters:`);
  console.log(`   Sequencing depth: ${SEQUENCING_DEPTH}×`);
  console.log(`   Error rate: ${ERROR_RATE} (${(ERROR_RATE*100).toFixed(3)}%)`);
  console.log(`   ctDNA levels: ${CTDNA_FRACTIONS.map(f => (f*100).toFixed(f < 0.001 ? 4 : f < 0.01 ? 3 : 1)+'%').join(', ')}`);
  console.log(`   Trinuclotide context error range: ${Object.values(TRINUC_ERROR_MULTIPLIERS).filter(v=>typeof v==='number').reduce((a,b)=>Math.min(a,b))}× to ${Object.values(TRINUC_ERROR_MULTIPLIERS).filter(v=>typeof v==='number').reduce((a,b)=>Math.max(a,b))}×`);

  // Also build background site arrays for specificity
  const variantPositions = new Set(variants.map(v => `${v.sample_id}:${v.chrom}:${v.pos}`));
  const backgroundSites = [];
  for (let i = 0; i < 1000; i++) {
    const sample = samples[i % samples.length];
    const pos = 500000 + i * 500;
    backgroundSites.push({ sample_id: sample.sample_id, cancer_type: sample.cancer_type, chrom: 'chr1', pos, is_true_variant: false });
  }

  const allSites = [...variants.map(v => ({ ...v, is_true_variant: true })), ...backgroundSites];
  console.log(`   True variants: ${variants.length}`);
  console.log(`   Background sites: ${backgroundSites.length}`);
  console.log(`   Total sites to downsample: ${allSites.length}`);

  // Downsample each site at each ctDNA level
  const results = {};
  const perFractionStats = {};

  CTDNA_FRACTIONS.forEach(ctdnaFrac => {
    const key = `ctdna_${ctdnaFrac}`;
    const label = `ctDNA ${ctdnaFrac >= 0.01 ? (ctdnaFrac*100).toFixed(1) : (ctdnaFrac*100).toFixed(ctdnaFrac < 0.001 ? 4 : 3)}%`;
    console.log(`\n   🔬 Downsampling at ${label}...`);

    const observations = [];
    let totalTrueMutantReads = 0, totalTrueBgReads = 0;
    let totalBgMutantReads = 0, totalBgBgReads = 0;

    allSites.forEach((site, idx) => {
      const context = getTrinucContext(site.pos);
      const errorMult = TRINUC_ERROR_MULTIPLIERS[context] || TRINUC_ERROR_MULTIPLIERS.default;
      const effectiveErrorRate = ERROR_RATE * errorMult;

      if (site.is_true_variant) {
        // True variant: observed reads = poisson(depth * VAF * ctDNA_fraction)
        const trueLambda = SEQUENCING_DEPTH * site.true_vaf * ctdnaFrac;
        const bgLambda = SEQUENCING_DEPTH * effectiveErrorRate;
        const mutantReads = poisson(trueLambda, rng);
        const bgReads = poisson(bgLambda, rng);
        totalTrueMutantReads += mutantReads;
        totalTrueBgReads += bgReads;
        observations.push({
          site_type: 'variant',
          sample_id: site.sample_id,
          cancer_type: site.cancer_type,
          gene: site.gene,
          chrom: site.chrom,
          pos: site.pos,
          true_vaf: site.true_vaf,
          ctdna_fraction: ctdnaFrac,
          mutant_reads: mutantReads,
          background_reads: bgReads,
          total_reads: SEQUENCING_DEPTH,
          observed_vaf: mutantReads / SEQUENCING_DEPTH,
          expected_vaf: site.true_vaf * ctdnaFrac,
          error_multiplier: errorMult,
          trinuc_context: context
        });
      } else {
        // Background site: only background reads (no true variant)
        const bgLambda = SEQUENCING_DEPTH * effectiveErrorRate;
        const mutantReads = poisson(bgLambda, rng);
        const bgReads = poisson(bgLambda, rng);
        totalBgMutantReads += mutantReads;
        totalBgBgReads += bgReads;
        observations.push({
          site_type: 'background',
          sample_id: site.sample_id,
          cancer_type: site.cancer_type,
          chrom: site.chrom,
          pos: site.pos,
          true_vaf: 0,
          ctdna_fraction: ctdnaFrac,
          mutant_reads: mutantReads,
          background_reads: bgReads,
          total_reads: SEQUENCING_DEPTH,
          observed_vaf: mutantReads / SEQUENCING_DEPTH,
          expected_vaf: 0,
          error_multiplier: errorMult,
          trinuc_context: context
        });
      }
    });

    results[key] = observations;

    const nVariant = variants.length;
    const nBg = backgroundSites.length;
    const meanTrueVaf = variants.reduce((s,v)=>s+v.true_vaf,0)/variants.length;
    perFractionStats[key] = {
      label,
      ctdna_fraction: ctdnaFrac,
      n_variant_sites: nVariant,
      n_background_sites: nBg,
      mean_true_vaf: meanTrueVaf,
      expected_mean_obs_vaf: meanTrueVaf * ctdnaFrac,
      mean_true_mutant_reads: totalTrueMutantReads / nVariant,
      mean_true_bg_reads: totalTrueBgReads / nVariant,
      mean_bg_mutant_reads: totalBgMutantReads / nBg,
      mean_bg_bg_reads: totalBgBgReads / nBg,
      snr_estimate: (totalTrueMutantReads/nVariant) / Math.max(1, (totalTrueBgReads/nVariant))
    };

    console.log(`      Mean mutant reads (true sites): ${perFractionStats[key].mean_true_mutant_reads.toFixed(2)}`);
    console.log(`      Mean mutant reads (bg sites):   ${perFractionStats[key].mean_bg_mutant_reads.toFixed(2)}`);
    console.log(`      SNR estimate: ${perFractionStats[key].snr_estimate.toFixed(4)}`);
  });

  // Output
  const output = {
    parameters: {
      sequencing_depth: SEQUENCING_DEPTH,
      error_rate: ERROR_RATE,
      ctdna_fractions: CTDNA_FRACTIONS,
      trinuc_error_multipliers: TRINUC_ERROR_MULTIPLIERS,
      seed: SEED,
      n_variants: variants.length,
      n_background_sites: backgroundSites.length,
      n_samples: samples.length
    },
    per_fraction_stats: perFractionStats,
    observations: results
  };

  fs.writeFileSync(OUTPUT_PATH, JSON.stringify(output, null, 2));
  console.log(`\n💾 Saved downsampled data to ${OUTPUT_PATH}`);
  console.log(`   File size: ${(fs.statSync(OUTPUT_PATH).size / 1024).toFixed(1)} KB`);
  console.log('\n✅ Step 2 complete.');
  console.log('='.repeat(70));
})();
