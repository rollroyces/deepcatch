#!/usr/bin/env python3
"""
DeepCatch Complete Python Validation Pipeline
=============================================

Produces ALL key results from FINAL_REAL_DATA_REPORT.md using Python
(numpy/scipy/sklearn). Each phase reproduces the corresponding Node.js
script's output with the same algorithms.

Usage:
    python run_all.py              # Full run (~10-30 min)
    python run_all.py --demo       # Demo run (~2 min, reduced scale)
    python run_all.py --phase 3    # Run specific phase (1-6)
    python run_all.py --skip-downsample  # Skip Phase 2, load from disk

Phases:
    1. Load TCGA/COSMIC data
    2. Downsample with 6 realistic confounders
    3. Head-to-head comparison vs Bie et al. (2023)
    4. CET longitudinal validation
    5. Tissue-of-Origin validation
    6. Clinical assay comparison + final report
"""

import json
import sys
import argparse
import time
import logging
from pathlib import Path

# Add parent directory for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from validation.py.config import (RESULTS_PY_DIR, PY_DOWNSAMPLED_PATH,
                                   PY_H2H_PATH, PY_CET_PATH, PY_TOO_PATH,
                                   PY_COMPARISON_PATH, CTDNA_LEVELS)
from validation.py.tcga_loader import load_tcga_data
from validation.py.statistical_tests import (compute_auc, bootstrap_auc,
                                             delong_test, sensitivity_at_specificity)
from validation.py.realistic_downsample import downsample_to_cfdna
from validation.py.performance_weighted_fusion import performance_weighted_fusion
from validation.py.head_to_head import run_head_to_head
from validation.py.cet_validation import run_cet_validation
from validation.py.too_validation import run_too_validation
from validation.py.compare_published import generate_clinical_comparison

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(name)s] %(levelname)s: %(message)s',
    datefmt='%H:%M:%S',
)
logger = logging.getLogger('run_all')

BANNER = """
╔══════════════════════════════════════════════════════════╗
║          DeepCatch Python Validation Pipeline           ║
║      Reproducing ALL claims from FINAL_REAL_DATA_       ║
║      REPORT.md using numpy/scipy/sklearn                ║
╚══════════════════════════════════════════════════════════╝
"""


def main():
    parser = argparse.ArgumentParser(description='DeepCatch Validation Pipeline')
    parser.add_argument('--demo', action='store_true', help='Run at reduced scale')
    parser.add_argument('--skip-downsample', action='store_true',
                        help='Skip Phase 2, load from disk')
    parser.add_argument('--phase', type=int, default=0,
                        help='Run only specific phase (1-6, 0=all)')
    args = parser.parse_args()

    print(BANNER)

    is_demo = args.demo
    demo_params = {
        'n_background_sites': 500,
        'ctdna_fractions': [0.001, 0.0005, 0.0001],
        'n_folds': 3,
        'n_bootstrap': 200,
        'cet_n_cancer': 50,
        'cet_n_healthy': 100,
        'cet_n_benign': 25,
        'too_n_types': 5,
        'too_n_per_type': 30,
    } if is_demo else {}

    all_start = time.time()

    # ── Phase 1: Load Data ────────────────────────────────────────────
    phase1_start = time.time()
    print("\n" + "=" * 60)
    print("[1/6] Loading TCGA/COSMIC data...")
    print("=" * 60)

    tcga = load_tcga_data(force_fallback=is_demo)
    n_samples = len(tcga['dataset']['samples'])
    n_variants = len(tcga['dataset']['variants'])
    n_cancer = sum(1 for s in tcga['dataset']['samples'] if s['is_cancer'])
    print(f"  ✓ Loaded {n_samples} samples ({n_cancer} cancer, {n_samples - n_cancer} healthy)")
    print(f"  ✓ {n_variants} somatic variants across "
          f"{len(set(v['cancer_type'] for v in tcga['dataset']['variants']))} cancer types")
    print(f"  Phase 1: {time.time() - phase1_start:.1f}s\n")

    if args.phase and args.phase == 1:
        return

    # ── Phase 2: Downsample ───────────────────────────────────────────
    phase2_start = time.time()
    print("\n" + "=" * 60)
    print("[2/6] Downsampling with realistic confounders...")
    print("=" * 60)

    if args.skip_downsample:
        try:
            import json
            with open(PY_DOWNSAMPLED_PATH) as f:
                downsampled = json.load(f)
            print(f"  ✓ Loaded from {PY_DOWNSAMPLED_PATH}")
        except FileNotFoundError:
            print("  ⚠️  No saved data found. Running downsampling anyway...")
            downsampled = downsample_to_cfdna(tcga, **demo_params)
    else:
        downsample_kwargs = {'tcga_data': tcga, **demo_params}
        downsampled = downsample_to_cfdna(**downsample_kwargs)

    logger.info("Phase 2 complete: %.1fs", time.time() - phase2_start)
    if args.phase and args.phase == 2:
        return

    # ── Phase 3: Head-to-Head ─────────────────────────────────────────
    phase3_start = time.time()
    print("\n" + "=" * 60)
    print("[3/6] Head-to-head comparison...")
    print("=" * 60)

    h2h_kwargs = {'downsampled_data': downsampled, **demo_params}
    h2h = run_head_to_head(**h2h_kwargs)

    # Print summary table
    print("\n  AUC vs ctDNA Fraction:")
    print(f"  {'ctDNA':>10s} {'Bie THEMIS':>12s} {'CAPP-Seq':>10s} "
          f"{'iDES':>8s} {'DC Variant':>12s} {'DC Multi-Modal':>15s}")
    print("  " + "-" * 68)
    for row in h2h.get('summary_table', []):
        frac_str = f"{row['ctDNA_fraction']*100:.3f}%"
        print(f"  {frac_str:>10s} "
              f"{row.get('bie_themis', 0):>12.4f} "
              f"{row.get('cappSeq', 0):>10.4f} "
              f"{row.get('ides', 0):>8.4f} "
              f"{row.get('deepcatch_variant', 0):>12.4f} "
              f"{row.get('deepcatch_multimodal', 0):>15.4f}")

    det_limit = h2h.get('detection_limit_ctdna_fraction')
    print(f"\n  Detection limit (AUC > 0.80): ctDNA fraction "
          f"{f'{det_limit*100:.2f}%' if det_limit else 'NOT REACHED'}")

    # DeLong test summary
    sig_count = 0
    total_comps = 0
    for frac_key, frac_result in h2h.get('per_fraction_results', {}).items():
        if 'delong_tests' in frac_result:
            for comp_name, comp in frac_result['delong_tests'].items():
                total_comps += 1
                if comp.get('significant'):
                    sig_count += 1
    print(f"  DeLong tests: {sig_count}/{total_comps} significant (p < 0.05)")

    logger.info("Phase 3 complete: %.1fs", time.time() - phase3_start)
    if args.phase and args.phase == 3:
        return

    # ── Phase 4: CET ──────────────────────────────────────────────────
    phase4_start = time.time()
    print("\n" + "=" * 60)
    print("[4/6] CET longitudinal validation...")
    print("=" * 60)

    cet_kwargs = demo_params if is_demo else {}
    cet_kwargs = {k: v for k, v in cet_kwargs.items() if k.startswith('cet_')}
    # Rename to match function signature
    if 'cet_n_cancer' in cet_kwargs:
        cet_kwargs['n_cancer'] = cet_kwargs.pop('cet_n_cancer')
        cet_kwargs['n_healthy'] = cet_kwargs.pop('cet_n_healthy')
        cet_kwargs['n_benign'] = cet_kwargs.pop('cet_n_benign')
    cet = run_cet_validation(**cet_kwargs)

    logger.info("Phase 4 complete: %.1fs", time.time() - phase4_start)
    if args.phase and args.phase == 4:
        return

    # ── Phase 5: TOO ──────────────────────────────────────────────────
    phase5_start = time.time()
    print("\n" + "=" * 60)
    print("[5/6] Tissue-of-Origin validation...")
    print("=" * 60)

    too_kwargs = {'n_cancer_types': demo_params.get('too_n_types', 8),
                  'n_per_type': demo_params.get('too_n_per_type', 50)} if is_demo else {}
    too = run_too_validation(**too_kwargs)

    logger.info("Phase 5 complete: %.1fs", time.time() - phase5_start)
    if args.phase and args.phase == 5:
        return

    # ── Phase 6: Comparison ───────────────────────────────────────────
    phase6_start = time.time()
    print("\n" + "=" * 60)
    print("[6/6] Clinical assay comparison...")
    print("=" * 60)

    comparison = generate_clinical_comparison(h2h_results=h2h,
                                              cet_results=cet,
                                              too_results=too)

    # ── FINAL REPORT ──────────────────────────────────────────────────
    total_time = time.time() - all_start

    print("\n" + "=" * 60)
    print("FINAL VERDICT")
    print("=" * 60)
    print()

    # Detection limit
    print(f"  1. Detection Limit: "
          f"{'✅' if det_limit else '❌'} "
          f"{f'{det_limit*100:.2f}% ctDNA' if det_limit else 'NOT DETERMINED'}")

    # Best multi-modal AUC
    best_auc = 0.0
    for row in h2h.get('summary_table', []):
        best_auc = max(best_auc, row.get('deepcatch_multimodal', 0))
    print(f"  2. Best Multi-Modal AUC: {best_auc:.4f}")

    # Multi-modal advantage
    sig_h2h = sig_count > 0
    print(f"  3. Multi-Modal Advantage: "
          f"{'✅ Statistically significant' if sig_h2h else '⚠️ Not significant'} "
          f"(DeLong test, {sig_count}/{total_comps} fractions p<0.05)")

    # CET targets
    cet_met = cet.get('targets', {}).get('both_met', False)
    c_sens = cet.get('performance', {}).get('sensitivity', 0)
    c_spec = cet.get('performance', {}).get('specificity_overall', 0)
    print(f"  4. CET Dual Target: {'✅' if cet_met else '❌'} "
          f"(sens={c_sens*100:.1f}%, spec={c_spec*100:.1f}%)")

    # Clinical validation
    print(f"  5. Clinical Validation: ❌ ZERO clinical samples")

    # TOO
    too_acc = too.get('performance', {}).get('cv_accuracy', 0)
    print(f"  6. TOO Accuracy: ❌ Simulation ({too_acc*100:.1f}%) — "
          f"not comparable to Grail 88.7%")

    # Summary
    print(f"\n  {"=" * 58}")
    print(f"  FINAL VERDICT: {comparison['honest_assessment'][:80]}...")
    print(f"  {"=" * 58}")

    print(f"\n  Total runtime: {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"\n  All results saved to: {RESULTS_PY_DIR}/")
    print(f"    - real_downsampled.json")
    print(f"    - head_to_head_results.json")
    print(f"    - cet_results.json")
    print(f"    - too_results.json")
    print(f"    - clinical_comparison.json")

    print(f"\n✅ All Python validations complete. 🦾")

    # ── Platform-specific output check ──
    _check_disk_results()


def _check_disk_results():
    """Verify all result files were written."""
    expected = [PY_DOWNSAMPLED_PATH, PY_H2H_PATH, PY_CET_PATH,
                PY_TOO_PATH, PY_COMPARISON_PATH]
    all_ok = True
    for path in expected:
        if path.exists():
            size_kb = path.stat().st_size / 1024
            print(f"  ✓ {path.name} ({size_kb:.0f} KB)")
        else:
            print(f"  ✗ {path.name} — MISSING!")
            all_ok = False

    if not all_ok:
        logger.warning("Some result files were not written. Check logs above.")


if __name__ == '__main__':
    main()
