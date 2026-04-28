#!/usr/bin/env node
/**
 * loadData.js - Load and parse the real TCGA fallback dataset
 * Produces: parsed summary with variant/background statistics
 */
const fs = require('fs');
const path = require('path');

const DATASET_PATH = path.join(__dirname, '..', 'tcga', 'tcga_cache', 'fallback_dataset.json');

(function main() {
  console.log('='.repeat(70));
  console.log('DEEPCATCH NODE.JS VALIDATION - STEP 1: Load Dataset');
  console.log('='.repeat(70));

  const raw = JSON.parse(fs.readFileSync(DATASET_PATH, 'utf8'));
  const variants = raw.ground_truth_variants;
  const samples = raw.sample_metadata;

  // Build variant index by sample
  const variantMap = {};
  variants.forEach(v => {
    if (!variantMap[v.sample_id]) variantMap[v.sample_id] = [];
    variantMap[v.sample_id].push(v);
  });

  // Cancer type stats
  const cancerTypes = {};
  const variantsPerType = {};
  const samplesPerType = {};
  samples.forEach(s => {
    cancerTypes[s.cancer_type] = true;
    samplesPerType[s.cancer_type] = (samplesPerType[s.cancer_type] || 0) + 1;
  });
  variants.forEach(v => {
    variantsPerType[v.cancer_type] = (variantsPerType[v.cancer_type] || 0) + 1;
  });

  // VAF distribution
  const vafs = variants.map(v => v.true_vaf).sort((a, b) => a - b);
  const vafStats = {
    min: vafs[0],
    q5: percentile(vafs, 5),
    q25: percentile(vafs, 25),
    median: percentile(vafs, 50),
    q75: percentile(vafs, 75),
    q95: percentile(vafs, 95),
    max: vafs[vafs.length - 1],
    mean: vafs.reduce((a, b) => a + b, 0) / vafs.length,
    sd: stddev(vafs)
  };

  // Gene distribution
  const geneCounts = {};
  const geneVAFs = {};
  variants.forEach(v => {
    geneCounts[v.gene] = (geneCounts[v.gene] || 0) + 1;
    if (!geneVAFs[v.gene]) geneVAFs[v.gene] = [];
    geneVAFs[v.gene].push(v.true_vaf);
  });

  // Tumor purity distribution
  const purities = samples.map(s => s.tumor_purity).sort((a, b) => a - b);
  const nMutDist = {};
  samples.forEach(s => { nMutDist[s.n_mutations] = (nMutDist[s.n_mutations] || 0) + 1; });

  // Background positions: all sample/chrom/pos that don't have a variant
  const variantPositions = new Set(variants.map(v => `${v.sample_id}:${v.chrom}:${v.pos}`));
  // Generate some background positions from samples that have room
  const backgroundPositions = [];
  const variantSamples = new Set(variants.map(v => v.sample_id));
  const nonVariantSamples = samples.map(s => s.sample_id).filter(s => !variantSamples.has(s));

  // Use samples without variants + a few random positions from variant samples
  // Generate background positions at predictable offsets
  const maxPos = 200000;
  for (let i = 0; i < 500; i++) {
    const sample = samples[i % samples.length];
    const pos = 1000000 + i * 1000; // offset beyond variant range
    const key = `${sample.sample_id}:chr1:${pos}`;
    if (!variantPositions.has(key)) {
      backgroundPositions.push({
        sample_id: sample.sample_id,
        cancer_type: sample.cancer_type,
        chrom: 'chr1',
        pos: pos,
        is_true_variant: false
      });
    }
  }

  // Summary
  const summary = {
    dataset: 'TCGA fallback dataset',
    total_samples: samples.length,
    total_variants: variants.length,
    cancer_types: Object.keys(cancerTypes),
    cancer_genes: Object.keys(geneCounts),
    samples_per_type: samplesPerType,
    variants_per_type: variantsPerType,
    samples_with_variants: variantSamples.size,
    samples_without_variants: nonVariantSamples.length,
    n_mutations_per_sample: nMutDist,
    variant_vaf_stats: vafStats,
    gene_variant_counts: geneCounts,
    tumor_purity_range: { min: purities[0], max: purities[purities.length - 1] },
    background_positions_generated: backgroundPositions.length
  };

  console.log('\n📊 Dataset Summary:');
  console.log(`   Total samples: ${summary.total_samples}`);
  console.log(`   Total variants: ${summary.total_variants}`);
  console.log(`   Cancer types: ${summary.cancer_types.join(', ')}`);
  console.log(`   Samples per type: ${JSON.stringify(summary.samples_per_type)}`);
  console.log(`   Variants per type: ${JSON.stringify(summary.variants_per_type)}`);
  console.log(`   Samples WITH variants: ${summary.samples_with_variants}`);
  console.log(`   Samples WITHOUT variants: ${summary.samples_without_variants}`);
  console.log(`   n_mutations distribution: ${JSON.stringify(summary.n_mutations_per_sample)}`);
  console.log(`\n📈 VAF Distribution:`);
  console.log(`   Min: ${(vafStats.min*100).toFixed(2)}%`);
  console.log(`   Q5: ${(vafStats.q5*100).toFixed(2)}%`);
  console.log(`   Q25: ${(vafStats.q25*100).toFixed(2)}%`);
  console.log(`   Median: ${(vafStats.median*100).toFixed(2)}%`);
  console.log(`   Q75: ${(vafStats.q75*100).toFixed(2)}%`);
  console.log(`   Q95: ${(vafStats.q95*100).toFixed(2)}%`);
  console.log(`   Max: ${(vafStats.max*100).toFixed(2)}%`);
  console.log(`   Mean ± SD: ${(vafStats.mean*100).toFixed(2)}% ± ${(vafStats.sd*100).toFixed(2)}%`);
  console.log(`\n🧬 Gene Counts:`);
  Object.entries(geneCounts).sort((a, b) => b[1] - a[1]).forEach(([gene, count]) => {
    console.log(`   ${gene}: ${count}`);
  });
  console.log(`\n🩸 Tumor Purity: ${(purities[0]*100).toFixed(1)}% – ${(purities[purities.length-1]*100).toFixed(1)}%`);
  console.log(`\n✅ Step 1 complete.`);
  console.log('='.repeat(70));

  // Save parsed data
  const parsed = { variants, samples, backgroundPositions, summary };
  fs.writeFileSync(
    path.join(__dirname, '..', '..', 'results', 'node', 'parsed_data.json'),
    JSON.stringify(parsed, null, 2)
  );
})();

function percentile(sorted, p) {
  const idx = (p / 100) * (sorted.length - 1);
  const lo = Math.floor(idx);
  const hi = Math.ceil(idx);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo);
}

function stddev(arr) {
  const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
  const variance = arr.reduce((s, x) => s + (x - mean) ** 2, 0) / (arr.length - 1);
  return Math.sqrt(variance);
}
