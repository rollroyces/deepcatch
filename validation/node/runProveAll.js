#!/usr/bin/env node
/**
 * runProveAll.js - Master Runner for DeepCatch Proof Pipeline
 * Executes all 4 missions in sequence, then compiles the proven report.
 */
const { execSync } = require('child_process');
const path = require('path');
const fs = require('fs');

const SCRIPT_DIR = __dirname;
const RESULTS_DIR = path.join(__dirname, '..', '..', 'results', 'node');

// Ensure results directory exists
fs.mkdirSync(RESULTS_DIR, { recursive: true });

const SCRIPTS = [
  { name: 'CET v2 (Specificity Fix)', file: 'validateCET_v2.js', output: 'cet_v2_results.json' },
  { name: 'Tissue-of-Origin', file: 'validateTOO.js', output: 'too_results.json' },
  { name: 'Head-to-Head vs Bie', file: 'headToHeadBie.js', output: 'headToHead_results.json' },
  { name: 'Multi-Cancer (10 types)', file: 'validateMultiCancer.js', output: 'multicancer_results.json' },
];

const startTime = Date.now();

console.log('╔══════════════════════════════════════════════════════════════════════╗');
console.log('║     DEEPCATCH PROVE PIPELINE — Master Runner                       ║');
console.log('║     4 Missions → REAL numbers → Proven Report                      ║');
console.log('╚══════════════════════════════════════════════════════════════════════╝');
console.log(`\n📅 Started: ${new Date().toISOString()}`);
console.log(`🖥️  Node.js: ${process.version}`);
console.log(`📂 Working Dir: ${SCRIPT_DIR}\n`);

let allPassed = true;
const results = [];

for (let i = 0; i < SCRIPTS.length; i++) {
  const script = SCRIPTS[i];
  const scriptPath = path.join(SCRIPT_DIR, script.file);
  const outputPath = path.join(RESULTS_DIR, script.output);

  console.log('─'.repeat(70));
  console.log(`\n🚀 MISSION ${i + 1}/${SCRIPTS.length}: ${script.name}`);
  console.log(`   Script: ${script.file}`);
  console.log(`   Output: ${script.output}\n`);

  try {
    const missionStart = Date.now();
    const stdout = execSync(`/usr/local/bin/node "${scriptPath}"`, {
      cwd: SCRIPT_DIR,
      encoding: 'utf8',
      timeout: 300000, // 5 minutes max per mission
      maxBuffer: 50 * 1024 * 1024 // 50MB
    });
    const elapsed = ((Date.now() - missionStart) / 1000).toFixed(1);
    console.log(stdout);

    // Verify output file exists
    if (fs.existsSync(outputPath)) {
      const stats = fs.statSync(outputPath);
      console.log(`\n✅ MISSION ${i + 1} COMPLETE (${elapsed}s)`);
      console.log(`   Output: ${outputPath} (${(stats.size / 1024).toFixed(1)} KB)`);
      results.push({ mission: script.name, status: 'PASS', file: script.output, size_kb: (stats.size / 1024).toFixed(1) });

      // Quick validation: check JSON is parseable
      try {
        JSON.parse(fs.readFileSync(outputPath, 'utf8'));
      } catch (e) {
        console.log(`   ⚠️  Output file is not valid JSON: ${e.message}`);
        results[results.length - 1].status = 'WARN';
      }
    } else {
      console.log(`\n❌ MISSION ${i + 1} FAILED: Output file not found`);
      results.push({ mission: script.name, status: 'FAIL', file: script.output });
      allPassed = false;
    }
  } catch (err) {
    console.error(`\n❌ MISSION ${i + 1} ERROR:`);
    console.error(`   ${err.message}`);
    if (err.stdout) console.error(`   STDOUT: ${err.stdout.slice(-500)}`);
    if (err.stderr) console.error(`   STDERR: ${err.stderr.slice(-500)}`);
    results.push({ mission: script.name, status: 'ERROR', error: err.message });
    allPassed = false;
  }
}

// ── Compile Report ──
console.log('\n' + '═'.repeat(70));
console.log('\n📊 COMPILING PROVEN VALIDATION REPORT...\n');

try {
  const reportStart = Date.now();
  const stdout = execSync(`/usr/local/bin/node "${path.join(SCRIPT_DIR, 'compileProvenReport.js')}"`, {
    cwd: SCRIPT_DIR,
    encoding: 'utf8',
    timeout: 60000
  });
  const elapsed = ((Date.now() - reportStart) / 1000).toFixed(1);
  console.log(stdout);
  console.log(`\n✅ Report compilation complete (${elapsed}s)`);
} catch (err) {
  console.error(`\n❌ Report compilation error: ${err.message}`);
}

// ── Summary ──
const totalElapsed = ((Date.now() - startTime) / 1000).toFixed(1);

console.log('\n' + '═'.repeat(70));
console.log('╔══════════════════════════════════════════════════════════════════════╗');
console.log('║     PIPELINE COMPLETE                                              ║');
console.log('╚══════════════════════════════════════════════════════════════════════╝');
console.log(`\n📊 Mission Results:`);
results.forEach(r => {
  const icon = r.status === 'PASS' ? '✅' : r.status === 'WARN' ? '⚠️' : '❌';
  console.log(`   ${icon} ${r.mission}: ${r.status} (${r.size_kb || 'N/A'} KB)`);
});

const reportPath = path.join(RESULTS_DIR, 'PROVEN_VALIDATION_REPORT.md');
const reportExists = fs.existsSync(reportPath);
console.log(`\n📄 Proven Report: ${reportExists ? '✅ Generated' : '❌ Missing'}`);
if (reportExists) {
  const reportSize = fs.statSync(reportPath).size;
  console.log(`   Path: ${reportPath} (${(reportSize / 1024).toFixed(1)} KB)`);
}

console.log(`\n⏱️  Total time: ${totalElapsed}s`);
console.log(`📂 Results directory: ${RESULTS_DIR}`);
console.log(`\n📋 Generated files:`);
fs.readdirSync(RESULTS_DIR).forEach(f => {
  if (f.endsWith('.json') || f.endsWith('.md')) {
    const stat = fs.statSync(path.join(RESULTS_DIR, f));
    console.log(`   ${f} (${(stat.size / 1024).toFixed(1)} KB)`);
  }
});

console.log('\n🦾 DeepCatch Prove Pipeline — Complete.\n');
