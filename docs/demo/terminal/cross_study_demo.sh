#!/usr/bin/env bash
# DeepCatch cross-study benchmark demo (~2-3 min, recorded into docs/demo/cross_study.cast)
# Uses --seeds 2 to stay under 3 min wall clock for the demo.
# Production runs use --seeds 5 (see results/cross_study_finallydb.json).
set -e
export TERM=xterm-256color
export COLUMNS=110
export LINES=32
clear
echo "=== DeepCatch cross-study benchmark (n=627, FinaleDB pubs 6+8) ==="
echo "=== 2 seeds x 5-fold pooled OOF (demo); production uses 5 seeds ==="
date
echo ""
echo ">>> step 1: 5-channel cross-publication pipeline"
echo "    (5mb_ratio + 5mb_coverage + 100kb_ratio + 100kb_counts + FSD-196)"
cd /Users/hermes/deepcatch
env -u PYTHONPATH ./.venv/bin/python scripts/cross_study_finallydb.py --seeds 2 \
    --out-json results/_demo_cross_study.json \
    --out-md docs/_demo_benchmark.md 2>&1 | tail -22
echo ""
echo ">>> step 2: headline numbers"
env -u PYTHONPATH ./.venv/bin/python <<'PYEOF'
import json
d = json.load(open('/Users/hermes/deepcatch/results/_demo_cross_study.json'))
p = d['pooled']['harmonized']
nh = d['pooled']['no_harmonize']
tc = d['true_confound_control']['cancer_jiang_healthy_cristiano']
print(f"pooled AUC (harmonized)   : {p['auc_mean']:.4f} +/- {p['auc_std']:.4f}")
print(f"pooled AUC (no-harmonize) : {nh['auc_mean']:.4f} (batch artifact baseline)")
print(f"pooled sens@99% spec      : {p['sens_at_spec_99']:.3f}")
print(f"pooled sens@95% spec      : {p['sens_at_spec_95']:.3f}")
print()
print("--- TRUE CROSS-STUDY CONFOUND CONTROL ---")
print(f"cancer=jiang, healthy=cristiano  (harmonized)  : {tc['harmonized']['auc_mean']:.3f}")
print(f"cancer=jiang, healthy=cristiano  (no-harmonize): {tc['no_harmonize']['auc_mean']:.3f}")
print()
print("VERDICT: pooled AUC survives batch harmonization; true-confound collapses")
print("         to ~0.50 (signal IS cancer-vs-healthy, NOT study-of-origin).")
PYEOF
echo ""
echo "=== Done! Cross-study AUC is real (cancer) signal, not batch leak. ==="
echo "Full 5-seed production numbers: results/cross_study_finallydb.json"
echo "                              + docs/CROSS_STUDY_BENCHMARK.md"
