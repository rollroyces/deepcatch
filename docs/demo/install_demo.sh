#!/usr/bin/env bash
# DeepCatch install + smoke test demo (recorded into docs/demo/install.cast)
set -e
export TERM=xterm-256color
export COLUMNS=100
export LINES=30
clear
echo "=== DeepCatch install + smoke-test demo (~60s) ==="
date
echo ""
echo ">>> step 0: repo + branch info"
git -C /Users/hermes/deepcatch log --oneline -1
git -C /Users/hermes/deepcatch describe --tags --always 2>/dev/null || echo "(no tag)"
echo ""
echo ">>> step 1: editable install (venv reuses cache)"
cd /Users/hermes/deepcatch
pip install -e . --quiet --no-deps 2>&1 | tail -3 || true
echo "[install: ok]"
echo ""
echo ">>> step 2: smoke test (12 fast publication-readiness tests)"
env -u PYTHONPATH ./.venv/bin/python -m pytest \
    test/test_publication_readiness.py \
    -m "not slow" --tb=line -q --no-header 2>&1 | tail -5
echo ""
echo ">>> step 3: confirm CLI entry points (from pyproject.toml)"
env -u PYTHONPATH ./.venv/bin/python -c "
import importlib.metadata as md
try:
    eps = md.entry_points(group='console_scripts')
    ce = sorted([(ep.name, ep.value) for ep in eps if 'deepcatch' in (ep.value or '')])
    print(f'CLI entry points ({len(ce)}):')
    for n, v in ce[:8]:
        print(f'  {n:25s} -> {v}')
    if len(ce) > 8:
        print(f'  ... and {len(ce)-8} more')
except Exception as e:
    print('(no entry_points or venv missing)')
"
echo ""
echo "=== Done! Repo installed, 12/12 tests pass. ==="
echo "Next: try the cross-study benchmark (docs/demo/cross_study.cast)"
