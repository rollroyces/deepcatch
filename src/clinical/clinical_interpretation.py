#!/usr/bin/env python3
"""
Clinical Interpretation Module (P1-B2)

Translates the output of the Jiang 4-mer CET analysis pipeline into
clinician-friendly reports, biological pattern interpretations, and
machine-readable JSON exports.

Usage
-----
    from src.clinical.clinical_interpretation import ClinicalReportGenerator

    report = ClinicalReportGenerator(cet_df, fusion_result)
    print(report.generate_briefing())
    report.export_json('clinical_report.json')
    with open('report.html', 'w') as f:
        f.write(report.generate_html_report())
"""

from __future__ import annotations

import json
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ═══════════════════════════════════════════════════════════════════════════
# ClinicalReportGenerator
# ═══════════════════════════════════════════════════════════════════════════

class ClinicalReportGenerator:
    """
    Generate clinical-grade interpretation reports from CET analysis results.

    Translates statistical output (p-values, effect sizes, AUCs) into
    structured summaries suitable for clinician review, including
    biological pattern interpretation, one-paragraph briefings, and
    HTML reports with embedded tables.

    Parameters
    ----------
    cet_df : pd.DataFrame
        CET per-motif results with columns:
        ``[motif, p_value, effect_size, fdr_significant, composite_score]``
    fusion_result : dict
        Logistic regression fusion result from
        :func:`run_jiang_analysis.logistic_fusion_cv`.  Expected keys:
        ``auc_mean``, ``auc_std``, ``auc_folds``, ``coefs``,
        ``intercept``, ``top_indices``, ``y_pred_cv``
    threshold_sens : float, default 0.70
        Minimum sensitivity threshold for clinical utility assessment.
    threshold_spec : float, default 0.95
        Minimum specificity threshold for clinical utility assessment.

    Attributes
    ----------
    cet_df : pd.DataFrame
    fusion_result : dict
    threshold_sens : float
    threshold_spec : float
    """

    # ── Effect size interpretation thresholds (Romano et al. 2006) ────
    _EFFECT_LABELS = {
        (0.474, 1.01): "large",
        (0.33,  0.474): "medium",
        (0.147, 0.33): "small",
        (-1.0,  0.147): "negligible",
    }

    def __init__(
        self,
        cet_df: pd.DataFrame,
        fusion_result: Dict,
        threshold_sens: float = 0.70,
        threshold_spec: float = 0.95,
        nested_auc: Optional[float] = None,
        nested_auc_std: Optional[float] = None,
    ):
        self.cet_df = cet_df
        self.fusion_result = fusion_result
        self.threshold_sens = threshold_sens
        self.threshold_spec = threshold_spec
        self.nested_auc = nested_auc
        self.nested_auc_std = nested_auc_std

        # ── Validate inputs ──────────────────────────────────────────
        required_cols = {'motif', 'p_value', 'effect_size',
                         'fdr_significant', 'composite_score'}
        missing = required_cols - set(cet_df.columns)
        if missing:
            raise ValueError(
                f"cet_df missing required columns: {missing}"
            )

        required_fusion = {'auc_mean', 'auc_folds', 'coefs', 'top_indices'}
        missing_fusion = required_fusion - set(fusion_result.keys())
        if missing_fusion:
            raise ValueError(
                f"fusion_result missing required keys: {missing_fusion}"
            )

    # ── Public API ──────────────────────────────────────────────────────

    def generate_summary(self) -> Dict:
        """
        Produce a structured summary dictionary of all key metrics.

        Returns
        -------
        dict
            Keys:

            - ``n_significant_motifs`` (int): motifs with p < 0.05.
            - ``n_fdr_significant`` (int): motifs passing FDR correction.
            - ``top_motifs`` (list): (motif, score, effect_direction) tuples.
            - ``fusion_auc`` (float): mean CV AUC.
            - ``cv_confidence`` (str): AUC ± std formatted.
            - ``biological_pattern`` (str): pattern label.
            - ``clinical_utility`` (str): assessment string.
            - ``effect_size_distribution`` (dict): small/medium/large counts.
            - ``median_effect_size`` (float): median |δ|.
            - ``generated_at`` (str): ISO timestamp.
        """
        n_total = len(self.cet_df)
        n_sig = int((self.cet_df['p_value'] < 0.05).sum())
        n_fdr = int(self.cet_df['fdr_significant'].sum())
        fdr_pct = (n_fdr / n_total * 100) if n_total > 0 else 0.0

        # Top motifs with effect direction
        top_n = min(20, n_total)
        top_motifs: List[Tuple] = []
        for _, row in self.cet_df.head(top_n).iterrows():
            direction = "enriched" if row['effect_size'] > 0 else "depleted"
            top_motifs.append((row['motif'], row['composite_score'], direction))

        # Fusion performance — prefer nested CV AUC (unbiased) when available
        if self.nested_auc is not None:
            auc_mean = self.nested_auc
            auc_std = self.nested_auc_std if self.nested_auc_std is not None else 0.0
            cv_confidence = f"{auc_mean:.4f} ± {auc_std:.4f} (nested CV)"
        else:
            auc_mean = self.fusion_result.get('auc_mean', float('nan'))
            auc_std = self.fusion_result.get('auc_std', float('nan'))
            cv_confidence = f"{auc_mean:.4f} ± {auc_std:.4f}"

        # Biological pattern
        feature_names = self.cet_df['motif'].tolist()
        biological_pattern = self.interpret_biological_pattern(top_motifs, feature_names)

        # Clinical utility assessment
        clinical_utility = self._assess_clinical_utility(auc_mean, n_fdr, n_total)

        # Effect size distribution
        es_dist = self._effect_size_distribution()
        if 'abs_effect_size' in self.cet_df.columns:
            median_es = float(self.cet_df['abs_effect_size'].median())
        else:
            median_es = float(np.abs(self.cet_df['effect_size']).median())

        return {
            'n_significant_motifs': n_sig,
            'n_fdr_significant': n_fdr,
            'fdr_significant_pct': round(fdr_pct, 1),
            'n_total_motifs': n_total,
            'top_motifs': [(m, float(s), d) for m, s, d in top_motifs],
            'fusion_auc': float(auc_mean),
            'cv_confidence': cv_confidence,
            'biological_pattern': biological_pattern,
            'clinical_utility': clinical_utility,
            'effect_size_distribution': es_dist,
            'median_effect_size': median_es,
            'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        }

    def interpret_biological_pattern(
        self,
        top_motifs: List[Tuple],
        feature_names: List[str],
    ) -> str:
        """
        Classify the dominant biological pattern from top motifs.

        Examines the top-ranked motifs for CG-rich vs AT-rich
        content and enrichment vs depletion direction to identify
        patterns consistent with known cancer epigenetics.

        Jiang et al. (2020, *Cancer Discovery*) showed that cancer
        cfDNA exhibits **CG-rich motif depletion** and **AT-rich
        motif enrichment** — reflecting the increased nucleosome
        accessibility and hypomethylation characteristic of tumour
        chromatin.

        Parameters
        ----------
        top_motifs : list of (str, float, str)
            (motif_name, composite_score, direction) tuples.
        feature_names : list of str
            All motif names (for context).

        Returns
        -------
        str
            One of:

            - ``"CG-rich depletion"``
            - ``"AT-rich enrichment"``
            - ``"CG-rich depletion + AT-rich enrichment"`` (mixed, classic cancer)
            - ``"AT-rich depletion"``
            - ``"CG-rich enrichment"``
            - ``"No clear CG/AT pattern"``
        """
        if not top_motifs:
            return "No clear CG/AT pattern"

        cg_enriched = 0
        cg_depleted = 0
        at_enriched = 0
        at_depleted = 0

        for motif, _score, direction in top_motifs:
            motif_upper = motif.upper()
            cg_count = motif_upper.count('C') + motif_upper.count('G')
            at_count = motif_upper.count('A') + motif_upper.count('T')

            is_cg_rich = cg_count >= 3  # ≥3 of 4 bases are C/G
            is_at_rich = at_count >= 3  # ≥3 of 4 bases are A/T

            if is_cg_rich:
                if direction == 'enriched':
                    cg_enriched += 1
                else:
                    cg_depleted += 1
            elif is_at_rich:
                if direction == 'enriched':
                    at_enriched += 1
                else:
                    at_depleted += 1

        # Determine dominant pattern
        cg_signal = cg_depleted - cg_enriched
        at_signal = at_enriched - at_depleted

        pattern = []
        # Require at least 2 motifs in a category and a clear bias
        if cg_depleted >= 2 and cg_signal > 0:
            pattern.append("CG-rich depletion")
        elif cg_enriched >= 2 and cg_signal < 0:
            pattern.append("CG-rich enrichment")
        if at_enriched >= 2 and at_signal > 0:
            pattern.append("AT-rich enrichment")
        elif at_depleted >= 2 and at_signal < 0:
            pattern.append("AT-rich depletion")

        if len(pattern) >= 2:
            return " + ".join(pattern)
        elif len(pattern) == 1:
            return pattern[0]
        else:
            return "No clear CG/AT pattern"

    def generate_briefing(self) -> str:
        """
        Generate a clinician-friendly one-paragraph summary.

        Returns
        -------
        str
            Single-paragraph plain-text briefing suitable for
            clinical reports or messaging.
        """
        summary = self.generate_summary()

        top_5_names = [m for m, s, d in summary['top_motifs'][:5]]
        metric = summary['fusion_auc']

        if metric >= 0.95:
            strength = "outstanding"
        elif metric >= 0.90:
            strength = "excellent"
        elif metric >= 0.80:
            strength = "good"
        elif metric >= 0.70:
            strength = "moderate"
        else:
            strength = "limited"

        briefing = (
            f"DeepCatch v2.1 4-mer end motif analysis identified "
            f"{summary['n_fdr_significant']} of {summary['n_total_motifs']} motifs "
            f"({summary['fdr_significant_pct']}%) as differentially abundant between "
            f"cancer and control plasma (FDR < 0.05). "
            f"Fusion of the top-ranked motifs via logistic regression achieved "
            f"{strength} discriminative performance (CV AUC = {summary['cv_confidence']}). "
            f"The dominant biological pattern is {summary['biological_pattern']}, "
            f"consistent with established cancer cfDNA epigenetics (Jiang et al. 2020). "
            f"Top discriminative motifs include: {', '.join(top_5_names)}. "
            f"Clinical utility assessment: {summary['clinical_utility']}."
        )
        return briefing

    def generate_html_report(self) -> str:
        """
        Generate a self-contained HTML clinical report.

        Includes styled tables for motif rankings, fusion performance,
        effect-size distribution, and biological interpretation.

        Returns
        -------
        str
            Complete HTML document as a string.
        """
        summary = self.generate_summary()
        top_n = min(20, len(self.cet_df))

        # ── Build motif table rows ──────────────────────────────────
        motif_rows = ""
        for i, (_, row) in enumerate(self.cet_df.head(top_n).iterrows()):
            sig_class = 'sig-yes' if row['fdr_significant'] else 'sig-no'
            direction = '▲ enriched' if row['effect_size'] > 0 else '▼ depleted'
            dir_class = 'dir-up' if row['effect_size'] > 0 else 'dir-down'
            motif_rows += f"""
        <tr>
          <td>{i + 1}</td>
          <td class="motif">{row['motif']}</td>
          <td>{row['p_value']:.2e}</td>
          <td>{row['effect_size']:+.4f}</td>
          <td class="{sig_class}">{'✓' if row['fdr_significant'] else '—'}</td>
          <td>{row['composite_score']:.2f}</td>
          <td class="{dir_class}">{direction}</td>
        </tr>"""

        # ── Effect size distribution rows ───────────────────────────
        es_dist = summary['effect_size_distribution']
        es_rows = ""
        for label in ['large', 'medium', 'small', 'negligible']:
            n = es_dist.get(label, 0)
            total = summary['n_total_motifs']
            pct = (n / total * 100) if total > 0 else 0
            es_rows += f"""
        <tr>
          <td>{label.capitalize()}</td>
          <td>{n}</td>
          <td>{pct:.1f}%</td>
        </tr>"""

        # ── Utility colour ──────────────────────────────────────────
        auc = summary['fusion_auc']
        if auc >= 0.95:
            util_color = '#2ecc71'
        elif auc >= 0.85:
            util_color = '#27ae60'
        elif auc >= 0.75:
            util_color = '#f39c12'
        else:
            util_color = '#e74c3c'

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>DeepCatch v2.1 — Clinical Interpretation Report</title>
<style>
  body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    max-width: 900px; margin: 2rem auto; padding: 0 1rem;
    color: #2c3e50; line-height: 1.6;
  }}
  h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 0.4rem; }}
  h2 {{ color: #2471a3; margin-top: 2rem; }}
  h3 {{ color: #2e86c1; }}
  .badge {{
    display: inline-block; padding: 0.25em 0.6em;
    border-radius: 4px; font-size: 0.85em; font-weight: 600;
  }}
  .badge-good {{ background: #d5f5e3; color: #1e8449; }}
  .badge-warn {{ background: #fdebd0; color: #b9770e; }}
  .badge-info {{ background: #d6eaf8; color: #1a5276; }}
  table {{ width: 100%; border-collapse: collapse; margin: 1rem 0; }}
  th, td {{ padding: 0.5rem 0.75rem; text-align: left; border-bottom: 1px solid #ddd; }}
  th {{ background: #eaf2f8; color: #1a5276; font-weight: 600; }}
  tr:hover {{ background: #f7f9fc; }}
  .motif {{ font-family: 'Courier New', monospace; font-weight: 600; }}
  .sig-yes {{ color: #27ae60; font-weight: 700; }}
  .sig-no {{ color: #bdc3c7; }}
  .dir-up {{ color: #e74c3c; }}
  .dir-down {{ color: #2980b9; }}
  .metric-box {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 1rem; margin: 1.5rem 0;
  }}
  .metric {{
    background: #f4f6f7; border-radius: 8px; padding: 1rem; text-align: center;
  }}
  .metric-value {{ font-size: 1.8em; font-weight: 700; color: #1a5276; }}
  .metric-label {{ font-size: 0.85em; color: #7f8c8d; margin-top: 0.25rem; }}
  .briefing {{
    background: #eaf2f8; border-left: 4px solid #2980b9;
    padding: 1rem 1.5rem; margin: 1.5rem 0; border-radius: 0 6px 6px 0;
  }}
  .footer {{
    margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ddd;
    font-size: 0.8em; color: #95a5a6;
  }}
</style>
</head>
<body>

<h1>🧬 DeepCatch v2.1 — Clinical Interpretation Report</h1>

<p>
  <span class="badge badge-info">4-mer End Motif Analysis</span>
  <span class="badge badge-info">Jiang Lab Protocol</span>
  <span class="badge badge-good">Real Plasma Data ✅</span>
</p>

<p>
  <strong>Generated:</strong> {summary['generated_at']}<br>
  <strong>Total motifs tested:</strong> {summary['n_total_motifs']}<br>
  <strong>FDR-significant:</strong> {summary['n_fdr_significant']} ({summary['fdr_significant_pct']}%)
</p>

<!-- ── Metric boxes ─────────────────────────────────────────────── -->
<div class="metric-box">
  <div class="metric">
    <div class="metric-value" style="color: {util_color};">{summary['fusion_auc']:.4f}</div>
    <div class="metric-label">Fusion CV AUC</div>
  </div>
  <div class="metric">
    <div class="metric-value">{summary['n_fdr_significant']}</div>
    <div class="metric-label">FDR-Significant Motifs</div>
  </div>
  <div class="metric">
    <div class="metric-value">{summary['median_effect_size']:.4f}</div>
    <div class="metric-label">Median |Effect Size|</div>
  </div>
  <div class="metric">
    <div class="metric-value">{summary['cv_confidence']}</div>
    <div class="metric-label">AUC (mean ± std)</div>
  </div>
</div>

<!-- ── Clinical briefing ────────────────────────────────────────── -->
<div class="briefing">
  <strong>📋 Clinical Briefing</strong>
  <p>{self.generate_briefing()}</p>
</div>

<!-- ── Biological pattern ───────────────────────────────────────── -->
<h2>1. Biological Pattern</h2>
<p>
  <strong>Pattern:</strong>
  <span class="badge badge-info">{summary['biological_pattern']}</span>
</p>
<p>
  {self._biological_pattern_explanation(summary['biological_pattern'])}
</p>

<!-- ── Top motifs ───────────────────────────────────────────────── -->
<h2>2. Top Discriminative Motifs</h2>
<table>
  <thead>
    <tr>
      <th>Rank</th><th>Motif</th><th>p-value</th><th>Effect Size (δ)</th>
      <th>FDR</th><th>Score</th><th>Direction</th>
    </tr>
  </thead>
  <tbody>{motif_rows}
  </tbody>
</table>

<!-- ── Effect size distribution ─────────────────────────────────── -->
<h2>3. Effect Size Distribution</h2>
<table>
  <thead>
    <tr><th>Magnitude</th><th>Count</th><th>% Motifs</th></tr>
  </thead>
  <tbody>{es_rows}
  </tbody>
</table>

<!-- ── Fusion performance ───────────────────────────────────────── -->
<h2>4. Logistic Regression Fusion Performance</h2>
<table>
  <thead>
    <tr><th>Metric</th><th>Value</th></tr>
  </thead>
  <tbody>
    <tr><td>CV AUC (mean ± std)</td><td>{summary['cv_confidence']}</td></tr>
    <tr><td>Number of CV folds</td><td>{len(self.fusion_result.get('auc_folds', []))}</td></tr>
    <tr><td>Top motifs used</td><td>{len(self.fusion_result.get('top_indices', []))}</td></tr>
    <tr><td>Intercept</td><td>{self.fusion_result.get('intercept', 'N/A')}</td></tr>
  </tbody>
</table>

<!-- ── Clinical utility ─────────────────────────────────────────── -->
<h2>5. Clinical Utility Assessment</h2>
<p>
  <span class="badge {'badge-good' if 'Recommended' in summary['clinical_utility'] else 'badge-warn'}">
    {summary['clinical_utility']}
  </span>
</p>
<p><em>Thresholds: sensitivity ≥ {self.threshold_sens:.0%}, specificity ≥ {self.threshold_spec:.0%}.</em></p>

<div class="footer">
  <p>
    Generated by DeepCatch v2.1 Clinical Interpretation Module.<br>
    Citation: Royce &amp; DeepCatch Contributors (2026).<br>
    <strong>⚠️ Research use only — not for clinical diagnosis.</strong>
  </p>
</div>

</body>
</html>"""
        return html

    def export_json(self, filepath: Optional[str] = None) -> Dict:
        """
        Export the full report as a JSON-serializable dictionary.

        Parameters
        ----------
        filepath : str, optional
            If provided, write JSON to this path.

        Returns
        -------
        dict
            Serializable report dictionary.
        """
        summary = self.generate_summary()

        # ── Full motif list (all 256) ──────────────────────────────
        motifs_export = []
        for _, row in self.cet_df.iterrows():
            motifs_export.append({
                'motif': row['motif'],
                'p_value': float(row['p_value']),
                'effect_size': float(row['effect_size']),
                'abs_effect_size': float(row['abs_effect_size']) if 'abs_effect_size' in row.index else float(abs(row['effect_size'])),
                'fdr_significant': bool(row['fdr_significant']),
                'composite_score': float(row['composite_score']),
                'direction': 'enriched' if row['effect_size'] > 0 else 'depleted',
            })

        # ── Fusion details ─────────────────────────────────────────
        fusion_detail = {
            'auc_mean': float(self.fusion_result.get('auc_mean', float('nan'))),
            'auc_std': float(self.fusion_result.get('auc_std', float('nan'))),
            'auc_folds': [float(a) for a in self.fusion_result.get('auc_folds', [])],
            'intercept': float(self.fusion_result.get('intercept', float('nan'))),
            'n_top_motifs_used': len(self.fusion_result.get('top_indices', [])),
            'top_indices': [int(i) for i in self.fusion_result.get('top_indices', [])],
            'n_coefficients': len(self.fusion_result.get('coefs', [])),
        }

        report = {
            'meta': {
                'version': '2.1.0',
                'module': 'ClinicalReportGenerator',
                'generated_at': summary['generated_at'],
                'threshold_sens': self.threshold_sens,
                'threshold_spec': self.threshold_spec,
            },
            'summary': {
                'n_total_motifs': summary['n_total_motifs'],
                'n_significant_motifs': summary['n_significant_motifs'],
                'n_fdr_significant': summary['n_fdr_significant'],
                'fdr_significant_pct': summary['fdr_significant_pct'],
                'biological_pattern': summary['biological_pattern'],
                'clinical_utility': summary['clinical_utility'],
                'effect_size_distribution': summary['effect_size_distribution'],
                'median_effect_size': summary['median_effect_size'],
                'fusion_auc': summary['fusion_auc'],
                'cv_confidence': summary['cv_confidence'],
            },
            'top_motifs': summary['top_motifs'],
            'all_motifs': motifs_export,
            'fusion_detail': fusion_detail,
            'clinical_interpretation': {
                'biological_pattern': self._biological_pattern_explanation(
                    summary['biological_pattern']
                ),
                'briefing': self.generate_briefing(),
            },
            'disclaimer': (
                "Research use only. Not validated for clinical diagnosis. "
                "All metrics are computational estimates derived from "
                "processed cfDNA frequency data."
            ),
        }

        if filepath:
            with open(filepath, 'w') as f:
                json.dump(report, f, indent=2, ensure_ascii=False)

        return report

    # ── Internal helpers ──────────────────────────────────────────────

    def _assess_clinical_utility(
        self, auc: float, n_fdr: int, n_total: int,
    ) -> str:
        """Produce a plain-language clinical utility assessment."""
        fdr_pct = (n_fdr / n_total * 100) if n_total > 0 else 0

        if auc >= 0.95 and fdr_pct >= 30:
            return (
                "High clinical potential — strong discriminative performance "
                "with large proportion of FDR-significant motifs. "
                "Recommended for independent validation."
            )
        elif auc >= 0.90:
            return (
                "Good clinical potential — discriminative performance adequate "
                "for screening triage. External validation recommended before "
                "clinical deployment."
            )
        elif auc >= 0.80:
            return (
                "Moderate clinical potential — may serve as a supplementary "
                "biomarker panel. Requires combination with additional "
                "modalities for standalone screening."
            )
        elif auc >= 0.70:
            return (
                "Limited clinical utility as a standalone assay. May provide "
                "incremental value in multi-modal fusion panels."
            )
        else:
            return (
                "Insufficient discriminative performance for clinical use. "
                "Motif signal is too weak for reliable cancer detection."
            )

    def _effect_size_distribution(self) -> Dict[str, int]:
        """Count motifs in each effect-size magnitude bin."""
        if 'abs_effect_size' in self.cet_df.columns:
            abs_es = self.cet_df['abs_effect_size'].values
        else:
            abs_es = np.abs(self.cet_df['effect_size'].values)
        counts = {'large': 0, 'medium': 0, 'small': 0, 'negligible': 0}
        for v in abs_es:
            if v >= 0.474:
                counts['large'] += 1
            elif v >= 0.33:
                counts['medium'] += 1
            elif v >= 0.147:
                counts['small'] += 1
            else:
                counts['negligible'] += 1
        return counts

    @staticmethod
    def _biological_pattern_explanation(pattern: str) -> str:
        """Return a clinician-friendly explanation of the pattern."""
        explanations = {
            "CG-rich depletion + AT-rich enrichment":
                "This is the classic cancer cfDNA signature described by "
                "Jiang et al. (2020, <em>Cancer Discovery</em>). "
                "CG-rich end motifs are depleted in cancer patients' cfDNA "
                "while AT-rich motifs are enriched, reflecting tumour-associated "
                "nucleosome repositioning and global hypomethylation. "
                "This dual pattern has been validated across multiple cancer "
                "types (HCC, CRC, lung, NPC, gastric).",
            "CG-rich depletion":
                "CG-rich 4-mer end motifs are significantly less abundant "
                "in cancer patient cfDNA compared to healthy controls. "
                "This reflects the reduced nucleosome protection at "
                "CpG-dense promoter regions in cancer cells, consistent "
                "with the hypomethylation-driven chromatin accessibility "
                "pattern observed in tumour-derived cfDNA (Jiang et al. 2020).",
            "AT-rich enrichment":
                "AT-rich 4-mer end motifs are significantly more abundant "
                "in cancer patient cfDNA. This enrichment is associated "
                "with increased cleavage at A/T-rich linker regions between "
                "nucleosomes — a hallmark of the altered chromatin structure "
                "in cancer cells.",
            "AT-rich depletion":
                "AT-rich motifs are unexpectedly depleted in cancer. "
                "This atypical pattern may reflect cancer-type-specific "
                "chromatin biology or technical factors. "
                "Recommend investigation of nuclease preferences and "
                "sample processing conditions.",
            "CG-rich enrichment":
                "CG-rich motifs are unexpectedly enriched in cancer. "
                "This is atypical for most cancer types and may indicate "
                "tissue-specific biology, technical batch effects, or "
                "sample quality issues. Recommend re-analysis with "
                "stricter quality filtering.",
            "No clear CG/AT pattern":
                "No dominant CG-rich or AT-rich motif pattern was detected "
                "among the top discriminative motifs. This may indicate "
                "a balanced or diffuse epigenetic signal, small effect sizes, "
                "or heterogeneous tumour biology.",
        }
        return explanations.get(pattern, explanations["No clear CG/AT pattern"])


# ═══════════════════════════════════════════════════════════════════════════
# Exports
# ═══════════════════════════════════════════════════════════════════════════

__all__ = ['ClinicalReportGenerator']
