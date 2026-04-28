#!/usr/bin/env node
/**
 * runAll.js - Master runner for DeepCatch Node.js validation pipeline
 * Executes all validation steps in sequence
 */
const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

const STEPS = [
  { file: 'loadData.js', name: 'Dataset Loading' },
  { file: 'downsample.js', name: 'cfDNA Downsampling' },
  { file: 'validateVariantCaller.js', name: 'Variant Caller Validation' },
  { file: 'validateFusion.js', name: 'Multi-Modal Fusion' },
  { file: 'validateCET.js', name: 'Longitudinal CET' },
  { file: 'compileResults.js', name: 'Results Compilation' }
];

const cwd = __dirname;

console.log('='.repeat(70));
console.log('🦾 DEEPCATCH NODE.JS VALIDATION PIPELINE');
console.log('='.repeat(70));
console.log(`Node.js ${process.version}`);
console.log(`Working directory: ${cwd}`);
console.log(`Start time: ${new Date().toISOString()}`);
console.log('='.repeat(70));

let passed = 0;
let failed = 0;
const startTime = Date.now();

STEPS.forEach((step, i) => {
  const stepStart = Date.now();
  console.log(`\n${'='.repeat(70)}`);
  console.log(`STEP ${i + 1}/${STEPS.length}: ${step.name}`);
  console.log(`Command: node ${step.file}`);
  console.log('='.repeat(70));

  try {
    const result = execSync(`node ${step.file}`, {
      cwd,
      stdio: 'inherit',
      timeout: 300000, // 5 min timeout
      env: { ...process.env, NODE_OPTIONS: '--max-old-space-size=4096' }
    });
    const elapsed = ((Date.now() - stepStart) / 1000).toFixed(1);
    console.log(`\n✅ STEP ${i + 1} PASSED (${elapsed}s)`);
    passed++;
  } catch (err) {
    const elapsed = ((Date.now() - stepStart) / 1000).toFixed(1);
    console.log(`\n❌ STEP ${i + 1} FAILED (${elapsed}s)`);
    console.log(`   Error: ${err.message}`);
    if (err.stderr) console.log(`   Stderr: ${err.stderr.toString().slice(0, 500)}`);
    failed++;
  }
});

const totalElapsed = ((Date.now() - startTime) / 1000).toFixed(1);
console.log(`\n${'='.repeat(70)}`);
console.log(`PIPELINE COMPLETE`);
console.log(`Passed: ${passed}/${STEPS.length} | Failed: ${failed}/${STEPS.length}`);
console.log(`Total time: ${totalElapsed}s`);
console.log('='.repeat(70));

if (failed > 0) process.exit(1);
