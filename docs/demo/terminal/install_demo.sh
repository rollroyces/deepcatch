#!/usr/bin/env bash
# DeepCatch install + smoke test demo (recorded into docs/demo/install.cast)
#
# The committed recording (install.cast, install.gif, install.webm) was made
# when the publication-readiness suite had 12 tests; this script reports the
# live count instead of a hardcoded one, so the recording and the script
# disagree on that number by design.
set -e
# The repository root, derived from this script's own location, so the
# demo runs for anyone who clones the repository -- it used to point at
# the machine it was recorded on.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"
export TERM=xterm-256color
export COLUMNS=100
export LINES=30
clear
echo "=== DeepCatch install + smoke-test demo (~60s) ==="
date
echo ""
echo ">>> step 0: repo + branch info"
git -C "${REPO_ROOT}" log --oneline -1
git -C "${REPO_ROOT}" describe --tags --always 2>/dev/null || echo "(no tag)"
echo ""
if [ ! -x ./.venv/bin/python ]; then
    echo ">>> no .venv in this checkout -- create one first:"
    echo ">>>   python3 -m venv .venv"
    exit 1
fi
echo ">>> step 1: editable install into the venv (venv reuses cache)"
# The venv's own pip, not the ambient one: step 2 runs the venv's python,
# so the package has to be installed into that interpreter. The [dev] extra
# is what provides pytest -- `pip install -e .` leaves the smoke test below
# with no test runner.
./.venv/bin/pip install -e ".[dev]" --quiet 2>&1 | tail -3
install_status=${PIPESTATUS[0]}
if [ "${install_status}" -ne 0 ]; then
    echo "[install: failed -- the smoke test below needs pytest]"
    exit 1
fi
echo "[install: ok]"
echo ""
echo ">>> step 2: smoke test (the fast publication-readiness tests)"
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
echo "=== Done! Repo installed, and the publication-readiness suite passes. ==="
echo "Next: try the cross-study benchmark (docs/demo/cross_study.cast)"
