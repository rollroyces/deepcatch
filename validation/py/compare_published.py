#!/usr/bin/env python3
"""
Honest Comparison vs Published Clinical Assays

Mirrors validation/node/comparePublished.js — every comparison comes with
caveats about simulation vs clinical reality.

Critical rules:
  - If DeepCatch has no clinical validation → state this clearly
  - If DeepCatch outperforms but at higher depth → note this
  - If TOO is simulation-only → compare honestly with Grail's clinical 88.7%
  - Every comparison must note important caveats

Published assays referenced:
  - Guardant360 (Odegaard 2018 Clin Cancer Res)
  - FoundationOne Liquid CDx (Woodhouse 2020 PLoS ONE)
  - Grail Galleri (Klein 2021 Ann Oncol; Jamshidi 2022 Cancer Cell)
  - CancerSEEK (Cohen 2018 Science; Lennon 2020 Science)
  - DELFI (Cristiano 2019 Nature; Mathios 2021 Nat Commun)
  - PanSeer (Chen 2020 Nat Commun)
  - Bie et al. 2023 THEMIS (Bie 2023 Nat Commun)
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any

import numpy as np

from .config import PY_COMPARISON_PATH, PY_H2H_PATH, PY_CET_PATH, PY_TOO_PATH

logger = logging.getLogger(__name__)

# ── Published Clinical Assay Reference Data ────────────────────────────────
PUBLISHED_ASSAYS = [
    {
        'name': 'Guardant360 (Guardant Health)',
        'type': 'ctDNA NGS (74 genes, tumor-informed available)',
        'citation': 'Odegaard 2018 Clin Cancer Res; Lanman 2015 PLoS ONE',
        'clinical_validation': True,
        'sample_size': '>200,000 clinical tests',
        'sensitivity_overall': 0.853,
        'specificity': 0.996,
        'lod_ctdna': 0.0001,  # 0.01%
        'cancer_types': 50,
        'too_accuracy': None,
        'sequencing_depth': '5,000× (clinical standard)',
        'notes': 'Commercially available. FDA-approved companion diagnostic. '
                 'Performance in advanced/metastatic setting primarily.',
        'stage_I_sensitivity': None,
    },
    {
        'name': 'FoundationOne Liquid CDx (Foundation Medicine)',
        'type': 'Hybrid capture NGS (324 genes)',
        'citation': 'Woodhouse 2020 PLoS ONE',
        'clinical_validation': True,
        'sample_size': '>30,000 clinical tests',
        'sensitivity_overall': 0.837,
        'specificity': 0.995,
        'lod_ctdna': 0.001,  # 0.1%
        'cancer_types': 50,
        'too_accuracy': None,
        'sequencing_depth': '~5,000×',
        'notes': 'FDA-approved. Higher LOD than Guardant360. '
                 'Better for tumor mutational burden assessment.',
        'stage_I_sensitivity': None,
    },
    {
        'name': 'Grail Galleri (MCED)',
        'type': 'Targeted methylation sequencing (>100,000 CpG regions)',
        'citation': 'Klein 2021 Ann Oncol; Jamshidi 2022 Cancer Cell; Liu 2020 Ann Oncol',
        'clinical_validation': True,
        'sample_size': '15,254 (CCGA substudy), >140,000 (NHS-Galleri trial)',
        'sensitivity_overall': 0.515,
        'specificity': 0.995,
        'lod_ctdna': None,  # Methylation-based
        'cancer_types': 50,
        'too_accuracy': 0.887,
        'sequencing_depth': '~30× WGBS equivalent (targeted)',
        'notes': 'Commercially available ($949). Sensitivity: 16.7% (I) → 90.5% (IV). '
                 'Best-in-class MCED breadth. FDA breakthrough device.',
        'stage_I_sensitivity': 0.167,
    },
    {
        'name': 'CancerSEEK (Thrive/Exact Sciences)',
        'type': 'Multi-analyte: ctDNA mutations (61 genes) + protein biomarkers (8 proteins)',
        'citation': 'Cohen 2018 Science; Lennon 2020 Science',
        'clinical_validation': True,
        'sample_size': '1,005 cancer + 812 healthy; 10,006 women (DETECT-A)',
        'sensitivity_overall': 0.70,
        'specificity': 0.99,
        'lod_ctdna': None,
        'cancer_types': 8,
        'too_accuracy': 0.83,
        'sequencing_depth': '~30,000× (targeted amplicon)',
        'notes': 'DETECT-A study doubled cancers detected by standard screening. '
                 'Combined with PET-CT for confirmation.',
        'stage_I_sensitivity': 0.43,
    },
    {
        'name': 'DELFI (Delfi Diagnostics)',
        'type': 'Genome-wide fragmentomics (low-coverage WGS)',
        'citation': 'Cristiano 2019 Nature; Mathios 2021 Nat Commun; Mazzone 2024 Cancer Discov',
        'clinical_validation': True,
        'sample_size': '958 (lung cancer clinical validation)',
        'sensitivity_overall': 0.73,
        'specificity': 0.98,
        'lod_ctdna': None,
        'cancer_types': 7,
        'too_accuracy': 0.75,
        'sequencing_depth': '1-2× WGS (low coverage)',
        'notes': 'Low-cost ($100-200). No targeted enrichment needed. '
                 'FDA breakthrough device for lung cancer screening.',
        'stage_I_sensitivity': 0.57,
    },
    {
        'name': 'PanSeer (Singlera Genomics)',
        'type': 'Targeted methylation sequencing (595 regions)',
        'citation': 'Chen 2020 Nat Commun',
        'clinical_validation': True,
        'sample_size': '605 asymptomatic (191 later diagnosed) + 223 cancer',
        'sensitivity_overall': 0.88,
        'specificity': 0.96,
        'lod_ctdna': 0.00001,
        'cancer_types': 5,
        'too_accuracy': None,
        'sequencing_depth': 'Targeted bisulfite PCR',
        'notes': 'ONLY assay demonstrating pre-symptomatic detection '
                 '(up to 4 years before diagnosis). GOLD STANDARD for longitudinal.',
        'stage_I_sensitivity': None,
        'pre_diagnosis_sensitivity': 0.95,
    },
    {
        'name': 'Bie et al. 2023 (THEMIS)',
        'type': 'Multi-modal: methylation + fragmentomics + CNA from single enzymatic assay',
        'citation': 'Bie 2023 Nat Commun',
        'clinical_validation': False,
        'sample_size': '780 cancer + 497 healthy',
        'sensitivity_overall': None,
        'specificity': 0.99,
        'lod_ctdna': 0.001,
        'cancer_types': 7,
        'too_accuracy': None,
        'sequencing_depth': 'WMS (whole methylome)',
        'notes': 'Academic validation only (no clinical trial). '
                 'Used for head-to-head in our study.',
        'stage_I_sensitivity': 0.73,
    },
]


# ── Extract DeepCatch performance ──
def _extract_dc_performance(h2h: Dict) -> Dict:
    """Extract DeepCatch performance from head-to-head results."""
    perf = {
        'variant_calling': {},
        'multimodal': {},
        'detection_limit': h2h.get('detection_limit_ctdna_fraction'),
        'summary': [],
    }

    if 'per_fraction_results' in h2h:
        for key, result in h2h['per_fraction_results'].items():
            if result.get('error'):
                continue
            dc_var = result.get('methods', {}).get('deepcatch_variant', {})
            dc_multi = result.get('methods', {}).get('deepcatch_multimodal', {})

            frac = result.get('ctdna_fraction', 0)
            perf['summary'].append({
                'ctdna_fraction': frac,
                'variant_auc': dc_var.get('auc'),
                'variant_sens_95spec': dc_var.get('sens_at_95_spec'),
                'variant_sens_99spec': dc_var.get('sens_at_99_spec'),
                'multimodal_auc': dc_multi.get('auc'),
                'multimodal_sens_95spec': dc_multi.get('sens_at_95_spec'),
                'multimodal_sens_99spec': dc_multi.get('sens_at_99_spec'),
            })

    return perf


def generate_clinical_comparison(h2h_results: Optional[Dict] = None,
                                 cet_results: Optional[Dict] = None,
                                 too_results: Optional[Dict] = None) -> Dict:
    """
    Generate honest comparison table vs published clinical assays.

    Args:
        h2h_results: Output from run_head_to_head().
        cet_results: Output from run_cet_validation().
        too_results: Output from run_too_validation().

    Returns:
        Full comparison dict with summary table and honest assessment.
    """
    # Try loading from disk if not provided
    if h2h_results is None:
        try:
            with open(PY_H2H_PATH) as f:
                h2h_results = json.load(f)
            logger.info("Loaded head-to-head results from disk")
        except (FileNotFoundError, json.JSONDecodeError):
            logger.warning("No head-to-head results available for comparison")
            h2h_results = {}

    if cet_results is None:
        try:
            with open(PY_CET_PATH) as f:
                cet_results = json.load(f)
            logger.info("Loaded CET results from disk")
        except (FileNotFoundError, json.JSONDecodeError):
            logger.warning("No CET results available")
            cet_results = None

    dc_perf = _extract_dc_performance(h2h_results)

    # Build summary table
    summary_table = []

    # Published assays
    for assay in PUBLISHED_ASSAYS:
        summary_table.append({
            'assay': assay['name'],
            'type': assay['type'],
            'clinical_validation': assay['clinical_validation'],
            'sample_size': assay['sample_size'],
            'sensitivity': assay['sensitivity_overall'],
            'specificity': assay['specificity'],
            'lod_ctdna': assay['lod_ctdna'],
            'too_accuracy': assay['too_accuracy'],
            'cancer_types': assay['cancer_types'],
            'depth': assay['sequencing_depth'],
            'notes': assay['notes'],
        })

    # DeepCatch rows
    best_dc = next((s for s in dc_perf.get('summary', [])
                     if s.get('ctdna_fraction') == 0.001), None)
    if not best_dc and dc_perf.get('summary'):
        best_dc = dc_perf['summary'][0]

    if best_dc:
        summary_table.append({
            'assay': 'DeepCatch (variant calling)',
            'type': 'Weighted multi-gene variant calling with trinucleotide error suppression',
            'clinical_validation': False,
            'sample_size': '0 clinical, simulation only',
            'sensitivity': best_dc.get('variant_sens_95spec'),
            'specificity': ('Simulated (at target specificity)'
                            if best_dc.get('variant_sens_95spec') is not None
                            else 'Not validated'),
            'lod_ctdna': dc_perf.get('detection_limit'),
            'too_accuracy': None,
            'cancer_types': 8,
            'depth': '50,000× (simulation)',
            'notes': '❌ SIMULATION ONLY. No clinical validation.',
        })

        summary_table.append({
            'assay': 'DeepCatch (multi-modal fusion)',
            'type': 'Performance-weighted fusion: variant calling + fragmentomics + methylation',
            'clinical_validation': False,
            'sample_size': '0 clinical, simulation only',
            'sensitivity': best_dc.get('multimodal_sens_95spec'),
            'specificity': ('Simulated (at target specificity)'
                            if best_dc.get('multimodal_sens_95spec') is not None
                            else 'Not validated'),
            'lod_ctdna': dc_perf.get('detection_limit'),
            'too_accuracy': 'SIMULATION ONLY — not validated',
            'cancer_types': 8,
            'depth': '50,000× (simulation)',
            'notes': '❌ SIMULATION ONLY. AUC from performance-weighted fusion.',
        })

    if cet_results and cet_results.get('performance'):
        p = cet_results['performance']
        summary_table.append({
            'assay': 'DeepCatch CET (longitudinal)',
            'type': 'Hierarchical Bayes Cumulative Evidence Tracking with Gompertz growth',
            'clinical_validation': False,
            'sample_size': '700 simulated patients (0 clinical)',
            'sensitivity': p.get('sensitivity'),
            'specificity': p.get('specificity_overall'),
            'lod_ctdna': None,
            'too_accuracy': None,
            'cancer_types': 8,
            'depth': 'N/A (longitudinal)',
            'notes': ('❌ SIMULATION ONLY. '
                      f"{'Meets dual target in simulation' if cet_results.get('targets', {}).get('both_met') else 'Does not meet targets in simulation'}. "
                      'No longitudinal patient data.'),
        })

    # ── Honest comparisons ──
    comparisons = _generate_detail_comparisons(dc_perf, cet_results, too_results)

    # ── Honest assessment ──
    det_limit = dc_perf.get('detection_limit')
    if det_limit is None:
        honest = ('❌ CANNOT ASSESS: DeepCatch head-to-head results not available. '
                  'No comparison possible.')
    elif det_limit > 0.001:
        honest = (
            f'❌ DEEPCATCH LOD IS TOO HIGH: Detection limit '
            f'({det_limit*100:.2f}% ctDNA) is {det_limit/0.0001:.0f}× worse '
            f'than Guardant360 clinical LOD (0.01%). Without wet-lab validation '
            f'showing comparable or better LOD, DeepCatch cannot claim clinical utility.'
        )
    else:
        honest = (
            f'⚠️ PARTIALLY PROMISING: DeepCatch SIMULATION shows competitive LOD '
            f'({det_limit*100:.2f}% ctDNA) but this is simulation-only. The gap '
            f'between simulation and clinical reality is large. DeepCatch requires: '
            f'(1) wet-lab validation on real patient samples, (2) head-to-head '
            f'comparison against Guardant360 on same samples, (3) demonstration '
            f'that performance advantage persists at matched sequencing depth.'
        )

    output = {
        'metadata': {
            'generated': True,
            'comparison_caveats': (
                'ALL DeepCatch results are SIMULATION-BASED. No clinical validation '
                'has been performed. Direct comparison to clinical assays is NOT '
                'scientifically valid — this comparison is provided for context only.'
            ),
        },
        'published_assays': PUBLISHED_ASSAYS,
        'deepcatch_extracted_performance': dc_perf,
        'deepcatch_cet_performance': cet_results.get('performance') if cet_results else None,
        'head_to_head_comparisons': comparisons,
        'summary_table': summary_table,
        'honest_assessment': honest,
        'requirements_for_validation': [
            '1. Test on real patient plasma samples (n ≥ 200 cancer, n ≥ 200 healthy)',
            '2. Head-to-head on same samples against Guardant360 or CAPP-Seq',
            '3. Match sequencing depth to clinical standard (5,000×) for fair comparison',
            '4. Validate TOO on multi-class real data with known primary tumors',
            '5. Longitudinal cohort: ≥500 patients, serial blood draws over ≥2 years',
            '6. Independent validation at a separate institution',
            '7. Pre-register analysis plan to prevent p-hacking',
        ],
        'publication_readiness': {
            'can_publish_as_commentary': True,
            'can_publish_as_methods_paper': 'Only if wet-lab validation is added',
            'can_publish_as_clinical_validation': False,
            'simulation_only': True,
            'next_step': 'Partner with clinical collaborators for real sample validation',
        },
    }

    # Save
    with open(PY_COMPARISON_PATH, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    logger.info(f"Saved clinical comparison to {PY_COMPARISON_PATH}")

    return output


def _generate_detail_comparisons(dc_perf: Dict,
                                 cet_results: Optional[Dict],
                                 too_results: Optional[Dict]) -> List[Dict]:
    """Generate detailed pairwise comparisons with honest caveats."""
    comparisons = []

    det_limit = dc_perf.get('detection_limit')

    # 1. Guardant360 LOD comparison
    comparisons.append({
        'comparison': 'DeepCatch vs Guardant360',
        'metric': 'Limit of Detection (ctDNA fraction)',
        'guardant360': {
            'value': 0.0001,
            'units': 'ctDNA fraction (0.01%)',
            'citation': 'Lanman 2015 PLoS ONE',
        },
        'deepcatch': {
            'value': det_limit,
            'units': 'ctDNA fraction (lowest with AUC>0.80)',
        },
        'caveats': [
            'DeepCatch LOD is simulation-based; Guardant360 LOD is clinical',
            'Guardant360 uses molecular barcoding (UMIs) with error correction',
            'Guardant360 has >200,000 clinical samples; DeepCatch has 0',
        ],
        'honest_assessment': (
            f'DeepCatch SIMULATION shows {"comparable" if det_limit and det_limit <= 0.0001 else "worse"} LOD. '
            'Must be validated in real patient samples.'
        ),
    })

    # 2. Grail Galleri comparison
    dc_sens_at_99 = None
    for s in dc_perf.get('summary', []):
        if s.get('multimodal_sens_99spec') and s['multimodal_sens_99spec'] > 0:
            dc_sens_at_99 = s['multimodal_sens_99spec']
            break

    comparisons.append({
        'comparison': 'DeepCatch vs Grail Galleri (MCED)',
        'metric': 'Sensitivity at 99.5% Specificity',
        'grail': {
            'value': 0.515,
            'units': 'Overall sensitivity',
            'citation': 'Klein 2021 Ann Oncol',
        },
        'deepcatch': {
            'value': dc_sens_at_99,
            'units': 'Estimated sensitivity at 99% spec (simulation)',
        },
        'caveats': [
            'Grail is methylation-based; DeepCatch is mutation + multi-modal',
            'Grail has clinical data from 15,254-subject CCGA study',
            'Grail FDA breakthrough device and commercially available',
            'Grail TOO accuracy: 88.7% CLINICAL; DeepCatch TOO: SIMULATION ONLY',
        ],
        'honest_assessment': (
            'DIRECT COMPARISON NOT POSSIBLE: Grail is a validated clinical test '
            'with >15,000 patients. DeepCatch has zero clinical patients.'
        ),
    })

    # 3. TOO comparison
    comparisons.append({
        'comparison': 'Tissue-of-Origin (TOO) Accuracy',
        'metric': 'TOO Accuracy',
        'grail': {'value': 0.887, 'units': '88.7% (CLINICAL)', 'citation': 'Jamshidi 2022'},
        'cancerseeek': {'value': 0.83, 'units': '83% (CLINICAL)', 'citation': 'Cohen 2018'},
        'deepcatch': {
            'value': too_results.get('performance', {}).get('cv_accuracy') if too_results else None,
            'units': 'SIMULATION ONLY — not clinically validated',
        },
        'caveats': [
            'Grail: Clinical TOO across 50+ cancer types (88.7% accuracy)',
            'CancerSEEK: Clinical TOO across 8 types (83% accuracy)',
            'DeepCatch: SIMULATION ONLY — meaningless for publication',
        ],
        'honest_assessment': (
            '❌ DEEPCATCH TOO IS NOT PROVEN. Previous TOO validation used '
            'simulation-only data. Real TOO accuracy is unknown.'
        ),
    })

    # 4. CET vs PanSeer
    if cet_results:
        c_perf = cet_results.get('performance', {})
        comparisons.append({
            'comparison': 'DeepCatch CET vs PanSeer (Longitudinal)',
            'metric': 'Pre-diagnosis Detection',
            'panseer': {'value': 0.95, 'units': 'Sensitivity 1-4 years pre-dx',
                        'citation': 'Chen 2020 Nat Commun'},
            'deepcatch_cet': {
                'sensitivity': c_perf.get('sensitivity'),
                'specificity': c_perf.get('specificity_overall'),
                'target_met': cet_results.get('targets', {}).get('both_met', False),
            },
            'caveats': [
                'PanSeer: REAL patient data from Taizhou Longitudinal Study (123,115 subjects)',
                'DeepCatch CET: SIMULATION only with Gompertz growth model',
                'PanSeer used archived blood samples collected years before diagnosis',
            ],
            'honest_assessment': (
                'PanSeer achieved 95% pre-diagnosis sensitivity in REAL patients. '
                'DeepCatch CET is simulation-only — not competitive.'
            ),
        })

    # 5. Cost/depth comparison
    comparisons.append({
        'comparison': 'Sequencing Requirements',
        'metric': 'Cost and Depth',
        'assays': {
            'Guardant360': '5,000× depth, $5,800',
            'Grail Galleri': '~30× WGBS, $949',
            'DELFI': '1-2× WGS, $100-200',
            'DeepCatch': '50,000× depth, Unknown cost',
        },
        'caveats': [
            'DeepCatch requires 10× more depth than Guardant360 (50,000× vs 5,000×)',
            'Commercial assays benefit from economies of scale',
        ],
        'honest_assessment': (
            '⚠️ DEEPCATCH REQUIRES 10× MORE SEQUENCING. At 50,000× depth vs 5,000×, '
            "DeepCatch's LOD advantage may be attributable to increased depth, not better algorithms."
        ),
    })

    return comparisons


# ── Demo ──────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    print("=" * 60)
    print("Clinical Comparison — Demo (with sample data)")
    print("=" * 60)

    # Mock data for demo
    mock_h2h = {
        'detection_limit_ctdna_fraction': 0.00001,
        'summary_table': [
            {'ctDNA_fraction': 0.001,
             'bie_themis': 0.82, 'deepcatch_multimodal': 0.92},
        ],
        'per_fraction_results': {
            'ctdna_0.001': {
                'ctdna_fraction': 0.001,
                'methods': {
                    'deepcatch_variant': {
                        'auc': 0.78, 'sens_at_95_spec': 0.12, 'sens_at_99_spec': 0.08,
                    },
                    'deepcatch_multimodal': {
                        'auc': 0.92, 'sens_at_95_spec': 0.65, 'sens_at_99_spec': 0.55,
                    },
                },
            },
        },
    }

    mock_cet = {
        'performance': {
            'sensitivity': 0.025,
            'specificity_overall': 0.97,
        },
        'targets': {'both_met': False},
    }

    mock_too = {
        'performance': {'cv_accuracy': 0.95},
    }

    result = generate_clinical_comparison(mock_h2h, mock_cet, mock_too)

    print(f"\n{result['honest_assessment']}")
    print(f"\nSummary table ({len(result['summary_table'])} rows):")
    for row in result['summary_table']:
        clinical = '✅' if row.get('clinical_validation') else '❌'
        sens = f"{row.get('sensitivity', 0)*100:.1f}%" if row.get('sensitivity') is not None else 'N/A'
        print(f"  {row['assay']}: sens={sens}, clinical={clinical}")

    print("\n✅ Clinical comparison complete.")
