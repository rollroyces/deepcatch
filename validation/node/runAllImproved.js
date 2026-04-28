#!/usr/bin/env node
/**
 * runAllImproved.js - Unified Re-Run of All Improved Validations
 */
const { execSync } = require('child_process');
const path = require('path');

const SCRIPT_DIR = __dirname;

const scripts = [
  { name: 'Smart Fusion', file: 'validateFusion.js', desc: '3 strategies + correlation analysis' },
  { name: 'CET Calibration', file: 'validateCET.js', desc: '4 methods + lambda tuning' },
  { name: 'New Biomarkers', file: 'validateNewBiomarkers.js', desc: '5 features + selection' }
];

const HR = '='.repeat(70);

console.log(HR);
console.log('DEEPCATCH IMPROVED VALIDATION - UNIFIED RE-RUN');
console.log(HR);
console.log(`Time: ${new Date().toISOString()}`);
console.log(`Node: ${process.version}`);
console.log(`Scripts: ${scripts.length}`);
console.log(HR);

const results = {};

for (let i = 0; i < scripts.length; i++) {
  const script = scripts[i];
  console.log(`\n\n${'▂'.repeat(70)}`);
  console.log(`[${i + 1}/${scripts.length}] ${script.name}: ${script.desc}`);
  console.log(`${'▂'.repeat(70)}`);

  const scriptPath = path.join(SCRIPT_DIR, script.file);

  try {
    const start = Date.now();
    const output = execSync(`node "${scriptPath}"`, {
      encoding: 'utf-8',
      stdio: ['pipe', 'pipe', 'pipe'],
      maxBuffer: 10 * 1024 * 1024,
      timeout: 120000
    });
    const elapsed = ((Date.now() - start) / 1000).toFixed(1);
    console.log(output);
    results[script.name] = { success: true, elapsedSeconds: parseFloat(elapsed) };
    console.log(`\n✅ ${script.name} completed in ${elapsed}s`);
  } catch (err) {
    results[script.name] = { success: false, error: err.message };
    console.log(`\n❌ ${script.name} FAILED: ${err.message}`);
    if (err.stderr) console.log(`stderr: ${err.stderr}`);
  }
}

console.log(`\n\n${HR}`);
console.log('SUMMARY');
console.log(HR);

for (const script of scripts) {
  const r = results[script.name];
  const status = r.success ? '✅' : '❌';
  const time = r.success ? ` (${r.elapsedSeconds}s)` : '';
  console.log(`   ${status} ${script.name}: ${script.desc}${time}`);
}

const allSuccess = Object.values(results).every(r => r.success);
console.log(`\n${allSuccess ? '🎉 All validations passed!' : '⚠️  Some validations failed. Check above for details.'}`);

// Write summary
const fs = require('fs');
const summaryPath = path.join(__dirname, '..', '..', 'results', 'node', 'improved_run_summary.json');
fs.writeFileSync(summaryPath, JSON.stringify({
  timestamp: new Date().toISOString(),
  node_version: process.version,
  results,
  all_success: allSuccess
}, null, 2));

console.log(`\n📝 Summary saved to ${summaryPath}`);
console.log(HR);
