#!/usr/bin/env python3
"""Ultra-early readiness CI gate.

Reads ``results/ultraearly_readiness.json`` (emitted by
``scripts/ultraearly_readiness.py``) and exits non-zero if any
artifact-driven readiness check has regressed to "gap" (i.e. the
required results/*.json is missing or invalid). This is a real-data
regression guard: existing artifact-driven checks must not silently
lose their evidence.

Environment variables
---------------------
``FAIL_ON_GAP``  default "true"
    When true, any ``gap`` verdict fails the gate. Set to "false" to
    allow new gaps during a transitional state (e.g. the
    ``cfdna-fragmentomics-pipeline`` sibling repo not being
    available).

Usage
-----
    python3 scripts/ultraearly_readiness_gate.py
    FAIL_ON_GAP=false python3 scripts/ultraearly_readiness_gate.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_JSON = _REPO_ROOT / "results" / "ultraearly_readiness.json"


def main() -> int:
    json_path = Path(os.environ.get("READINESS_JSON", str(DEFAULT_JSON)))
    fail_on_gap = os.environ.get("FAIL_ON_GAP", "true").strip().lower() in {
        "1", "true", "yes", "on",
    }

    if not json_path.exists():
        print(f"❌ FAIL: readiness JSON missing at {json_path}", file=sys.stderr)
        print("   Re-run scripts/ultraearly_readiness.py first.", file=sys.stderr)
        return 1
    try:
        with open(json_path) as f:
            r = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"❌ FAIL: could not load readiness JSON: {e}", file=sys.stderr)
        return 1

    counts = r.get("summary_counts", {})
    headline = r.get("headline_status", "(missing)")
    print(json.dumps(counts, indent=2))
    print(f"Headline: {headline}")

    gaps = r.get("gaps", [])
    if gaps:
        print(f"\n❌ {len(gaps)} readiness gap(s) detected:")
        for g in gaps:
            print(f"   - {g.get('name', '?')}: {g.get('notes', '')[:200]}")

    not_ready = [
        c for c in r.get("checks", [])
        if c.get("verdict") == "not_ready"
    ]
    if not_ready:
        print(f"\n⚠️  {len(not_ready)} check(s) currently 'not_ready':")
        for c in not_ready:
            print(f"   - {c.get('name', '?')}: value={c.get('value')}")

    if fail_on_gap and counts.get("gap", 0) > 0:
        print(
            f"\nFAIL_ON_GAP=true and counts.gap={counts.get('gap', 0)} > 0 → "
            f"failing the build.",
            file=sys.stderr,
        )
        return 1

    print("\n✅ Ultra-early readiness gate passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())