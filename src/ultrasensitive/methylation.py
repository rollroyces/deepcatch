"""
Methylation analysis for cancer detection and tissue-of-origin prediction.

Ported from ``cancer-screening/src/core/methylation.py`` (SHA 4b880edb).

Methylation patterns are among the earliest molecular changes in carcinogenesis.
Cancer-specific methylation signatures can be detected from cfDNA, enabling:
  1. Cancer detection (aberrant methylation vs. normal)
  2. Tissue of origin prediction (which organ the tumor came from)
  3. Tumor fraction estimation (from methylation entropy)

Methods:
  - Bisulfite sequencing (WGBS) -- gold standard
  - TAPS (TET-assisted pyridine borane sequencing) -- higher recovery
  - EPIC/450K arrays -- cost-effective, widely available
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class MethylationProcessor:
    """
    Process bisulfite sequencing data for cancer detection.

    Key concepts:
      - Beta value: fraction of methylated reads at a CpG site
        beta = M / (M + U + alpha)
      - M-value: logit transformation of beta (more normally distributed)
      - Differentially methylated regions (DMRs)
      - Tumor fraction estimated from cancer-specific methylation patterns
    """

    TISSUE_DMR_PANELS = {
        'lung_adenocarcinoma': {
            'genes': ['RASSF1A', 'SHOX2', 'RARB', 'CDKN2A', 'MGMT', 'APC'],
            'hypomethylated': ['SATalpha', 'LINE1'],
            'hyper_beta_shift': 0.2,
        },
        'colorectal': {
            'genes': ['SEPT9', 'VIM', 'CST6', 'GDNF', 'RASSF1A', 'SFRP2'],
            'hypomethylated': ['SATalpha', 'LINE1'],
            'hyper_beta_shift': 0.15,
        },
        'breast': {
            'genes': ['RASSF1A', 'BRCA1_promoter', 'RARB', 'SFRP1', 'MAL'],
            'hypomethylated': [],
            'hyper_beta_shift': 0.18,
        },
        'pancreas': {
            'genes': ['RASSF1A', 'NPTX2', 'UCHL1', 'CLDN5', 'BNC1'],
            'hypomethylated': ['SATalpha'],
            'hyper_beta_shift': 0.22,
        },
        'liver': {
            'genes': ['RASSF1A', 'RASSF2', 'SEPT9', 'NDRG2', 'CDKN2A'],
            'hypomethylated': [],
            'hyper_beta_shift': 0.20,
        },
        'prostate': {
            'genes': ['GSTP1', 'RASSF1A', 'APC', 'RARB', 'PITX2'],
            'hypomethylated': [],
            'hyper_beta_shift': 0.17,
        },
    }

    def __init__(
        self,
        smoothing_constant: float = 100.0,
        min_coverage: int = 5,
        min_cpg_per_gene: int = 3,
        tumor_beta_threshold: float = 0.3,
        normal_beta_threshold: float = 0.3,
    ):
        self.smoothing_constant = smoothing_constant
        self.min_coverage = min_coverage
        self.min_cpg_per_gene = min_cpg_per_gene
        self.tumor_beta_threshold = tumor_beta_threshold
        self.normal_beta_threshold = normal_beta_threshold

        logger.info(
            f"MethylationProcessor initialized: "
            f"alpha={smoothing_constant}, min_cov={min_coverage}"
        )

    def compute_beta_values(
        self,
        methylated_counts: np.ndarray,
        unmethylated_counts: np.ndarray,
    ) -> np.ndarray:
        """Compute beta values from counts with smoothing."""
        m = np.array(methylated_counts, dtype=float)
        u = np.array(unmethylated_counts, dtype=float)
        alpha = self.smoothing_constant

        beta = (m + alpha) / (m + u + 2 * alpha)
        return np.clip(beta, 0.0, 1.0)

    def compute_m_values(self, beta: np.ndarray) -> np.ndarray:
        """Logit-transform beta values."""
        beta = np.clip(beta, 1e-6, 1 - 1e-6)
        m = np.log2(beta / (1 - beta))
        return m

    def load_from_bsseq(self, bsseq_path: str) -> pd.DataFrame:
        """Load bisulfite sequencing data from a Bismark coverage file."""
        df = pd.read_csv(
            bsseq_path,
            sep='\t',
            names=['chrom', 'pos', 'methylated', 'unmethylated'],
            dtype={'chrom': str, 'pos': int, 'methylated': int, 'unmethylated': int},
        )
        df = df[df['methylated'] + df['unmethylated'] >= self.min_coverage]
        df['beta'] = self.compute_beta_values(df['methylated'], df['unmethylated'])
        df['m_value'] = self.compute_m_values(df['beta'])
        logger.info(f"Loaded {len(df)} CpG sites from {bsseq_path}")
        return df

    def compute_tumor_fraction(
        self,
        sample_df: pd.DataFrame,
        reference_normal_df: pd.DataFrame,
    ) -> float:
        """Estimate tumor fraction from methylation."""
        if sample_df.empty or reference_normal_df.empty:
            return 0.0

        merged = sample_df.merge(
            reference_normal_df,
            on=['chrom', 'pos'],
            suffixes=('_sample', '_normal'),
            how='inner',
        )

        if len(merged) < 100:
            logger.warning("Too few overlapping CpG sites for TF estimation")
            return 0.0

        cancer_signal = merged[
            (merged['beta_normal'] < 0.2) & (merged['beta_sample'] > 0.3)
        ]

        if len(cancer_signal) == 0:
            return 0.0

        estimated_tfs = (cancer_signal['beta_sample'] - cancer_signal['beta_normal']) / 0.8
        estimated_tfs = np.clip(estimated_tfs, 0.0, 1.0)

        tumor_fraction = float(np.median(estimated_tfs))
        logger.info(
            f"Estimated tumor fraction from methylation: {tumor_fraction:.2%} "
            f"(from {len(cancer_signal)} informative CpG sites)"
        )
        return tumor_fraction

    def predict_tissue_of_origin(
        self,
        sample_df: pd.DataFrame,
    ) -> dict:
        """Predict which tissue the tumor originated from based on methylation."""
        if sample_df.empty:
            return {}

        tissue_scores = {}

        for tissue, panel in self.TISSUE_DMR_PANELS.items():
            score = 0
            informative_sites = 0

            for gene in panel['genes']:
                gene_sites = sample_df[
                    sample_df.get('gene', pd.Series(dtype=str)) == gene
                ]

                if len(gene_sites) >= self.min_cpg_per_gene:
                    mean_beta = gene_sites['beta'].mean()
                    if mean_beta > self.tumor_beta_threshold:
                        score += 1
                        informative_sites += len(gene_sites)

            if len(panel['genes']) > 0:
                tissue_scores[tissue] = score / len(panel['genes'])
            else:
                tissue_scores[tissue] = 0.0

        ranked = dict(sorted(tissue_scores.items(), key=lambda x: x[1], reverse=True))

        top_tissue = max(ranked, key=ranked.get) if ranked else 'unknown'
        confidence = ranked.get(top_tissue, 0.0)

        logger.info(
            f"Tissue of origin prediction: {top_tissue} (confidence={confidence:.2f})"
        )
        logger.info(f"Full ranking: {ranked}")

        return ranked

    def compute_methylation_entropy(self, sample_df: pd.DataFrame) -> float:
        """Compute methylation entropy as a proxy for epigenetic dysregulation."""
        if sample_df.empty:
            return 0.0

        beta = sample_df['beta'].values
        bins = np.digitize(beta, np.linspace(0, 1, 11))
        counts = np.bincount(bins, minlength=11).astype(float)
        probs = counts / counts.sum()
        probs = probs[probs > 0]

        entropy = -np.sum(probs * np.log2(probs))
        max_entropy = np.log2(11)

        normalized_entropy = entropy / max_entropy

        return float(normalized_entropy)

    def detect_dmrs(
        self,
        sample_df: pd.DataFrame,
        reference_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Detect Differentially Methylated Regions (DMRs)."""
        if reference_df is None:
            logger.info("No reference provided -- using built-in DMR panels")
            return pd.DataFrame()

        merged = sample_df.merge(
            reference_df,
            on=['chrom', 'pos'],
            suffixes=('_sample', '_reference'),
            how='inner',
        )

        if len(merged) < 50:
            return pd.DataFrame()

        merged['delta_beta'] = merged['beta_sample'] - merged['beta_reference']
        merged['window'] = (merged['pos'] // 500) * 500

        dmr_candidates = []
        for (chrom, window), group in merged.groupby(['chrom', 'window']):
            mean_delta = group['delta_beta'].mean()
            n_cpg = len(group)

            if n_cpg >= 3 and abs(mean_delta) > 0.1:
                dmr_candidates.append({
                    'chrom': chrom,
                    'start': int(window),
                    'end': int(window + 500),
                    'mean_delta_beta': round(float(mean_delta), 4),
                    'n_cpg_sites': n_cpg,
                    'direction': 'hyper' if mean_delta > 0 else 'hypo',
                    'p_value': float(stats.ttest_1samp(group['delta_beta'], 0)[1]),
                })

        dmrs = pd.DataFrame(dmr_candidates)
        if not dmrs.empty:
            dmrs['q_value'] = stats.false_discovery_control(
                dmrs['p_value'], method='bh'
            )
            dmrs = dmrs[dmrs['q_value'] < 0.05]

        logger.info(f"Found {len(dmrs)} significant DMRs")
        return dmrs
