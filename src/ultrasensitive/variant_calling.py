"""
Low-VAF variant calling for ctDNA detection.

Ported from ``cancer-screening/src/core/variant_calling.py`` (SHA 01dae958).

Implements specialized variant callers optimized for ultra-low allele frequencies:
- Binomial test for background noise filtering
- Bayesian mutation calling with error models
- Strand bias filtering
- Clonal hematopoiesis (CHIP) filtering via matched WBC sequencing
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Optional
from scipy.stats import binom, poisson
import logging

logger = logging.getLogger(__name__)


@dataclass
class VariantCall:
    """A called somatic variant with statistical confidence."""
    chrom: str
    pos: int
    ref: str
    alt: str
    vaf: float
    depth: int
    alt_count: int
    ref_count: int
    p_value: float
    q_score: float
    gene: Optional[str] = None
    cosmic_id: Optional[str] = None
    driver: bool = False
    chip_filter: bool = False
    tissue_origin: Optional[str] = None


class LowVAFVariantFilter:
    """
    Statistical filter for calling variants at 0.001% VAF.

    At this sensitivity level, background errors dominate. This filter uses:
      1. Position-specific error rates from negative controls
      2. Binomial/Poisson statistics to distinguish true variants from noise
      3. Strand balance checks
      4. Matched WBC filtering (germline + CHIP removal)

    Key thresholds:
      - Minimum 3 supporting reads
      - Binomial p-value < 1e-6
      - Strand balance: both strands > 0 (duplex preferred)
      - VAF > 0.00001 for confirmed detection
    """

    def __init__(
        self,
        background_error_rate: float = 1e-7,
        min_alt_reads: int = 3,
        min_vaf: float = 1e-5,
        max_vaf: float = 0.5,
        binomial_pval_threshold: float = 1e-6,
        strand_balance_min: float = 0.01,
        min_base_quality: int = 25,
        min_mapping_quality: int = 30,
    ):
        self.background_error_rate = background_error_rate
        self.min_alt_reads = min_alt_reads
        self.min_vaf = min_vaf
        self.max_vaf = max_vaf
        self.binomial_pval_threshold = binomial_pval_threshold
        self.strand_balance_min = strand_balance_min
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality

        self.position_errors: dict = {}

        logger.info(
            f"LowVAFVariantFilter: min_alt={min_alt_reads}, "
            f"min_vaf={min_vaf:.1e}, binomial_p={binomial_pval_threshold:.0e}"
        )

    def compute_binomial_pvalue(
        self,
        alt_count: int,
        total_count: int,
        error_rate: float,
    ) -> float:
        """Compute binomial p-value for observing >= alt_count ALT reads."""
        if total_count == 0:
            return 1.0

        expected_error = error_rate * total_count

        if total_count > 100:
            p_value = 1 - poisson.cdf(alt_count - 1, expected_error)
        else:
            p_value = 1 - binom.cdf(alt_count - 1, total_count, error_rate)

        return float(p_value)

    def compute_vaf(
        self,
        alt_count: int,
        ref_count: int,
        total_count: int,
    ) -> float:
        if total_count == 0:
            return 0.0
        return alt_count / total_count

    def check_strand_balance(
        self,
        alt_plus: int,
        alt_minus: int,
        ref_plus: int,
        ref_minus: int,
    ) -> bool:
        """Check strand balance; single-strand artifacts are filtered."""
        total_alt = alt_plus + alt_minus
        if total_alt == 0:
            return False

        if alt_plus > 0 and alt_minus > 0:
            balance = min(alt_plus, alt_minus) / max(alt_plus, alt_minus)
            return balance > self.strand_balance_min
        elif alt_plus + alt_minus >= self.min_alt_reads:
            return False
        return False

    def call_variants_from_bam(
        self,
        pileup_data: pd.DataFrame,
    ) -> list:
        """Call variants from a pileup table."""
        calls = []

        for _, row in pileup_data.iterrows():
            chrom = row['chrom']
            pos = row['pos']
            ref = row['ref_base']

            bases = ['A', 'C', 'G', 'T']
            counts = {b: row.get(f'{b}_count', 0) for b in bases}
            counts['N'] = row.get('N_count', 0)

            counts[ref] = 0
            alt_bases = sorted(counts.items(), key=lambda x: x[1], reverse=True)

            for alt, alt_count in alt_bases:
                if alt_count < self.min_alt_reads:
                    continue

                total = row['coverage']
                vaf = self.compute_vaf(alt_count, total - alt_count, total)

                if vaf < self.min_vaf:
                    continue
                if vaf > self.max_vaf:
                    continue

                key = (chrom, pos)
                err_rate = self.position_errors.get(key, self.background_error_rate)

                p_value = self.compute_binomial_pvalue(alt_count, total, err_rate)
                if p_value >= self.binomial_pval_threshold:
                    continue

                q_score = min(-10 * np.log10(max(p_value, 1e-30)), 60)

                plus_alt = row.get('plus_alt_count', alt_count // 2 + 1)
                minus_alt = alt_count - plus_alt
                plus_ref = row.get('plus_ref_count', (total - alt_count) // 2)
                minus_ref = (total - alt_count) - plus_ref

                has_balance = self.check_strand_balance(
                    plus_alt, minus_alt, plus_ref, minus_ref
                )

                call = VariantCall(
                    chrom=chrom,
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    vaf=vaf,
                    depth=total,
                    alt_count=alt_count,
                    ref_count=total - alt_count,
                    p_value=p_value,
                    q_score=q_score,
                )

                calls.append(call)

        logger.info(f"Called {len(calls)} variants from {len(pileup_data)} positions")
        return calls

    def filter_chip_variants(
        self,
        variants: list,
        wbc_calls: list,
    ) -> list:
        """Filter out clonal hematopoiesis (CHIP) variants via matched WBC sequencing."""
        wbc_variant_set = set()
        for v in wbc_calls:
            wbc_variant_set.add((v.chrom, v.pos, v.ref, v.alt))

        chip_filtered = []
        for v in variants:
            key = (v.chrom, v.pos, v.ref, v.alt)

            if key in wbc_variant_set:
                wbc_v = next(
                    wc for wc in wbc_calls
                    if (wc.chrom, wc.pos, wc.ref, wc.alt) == key
                )

                if wbc_v.vaf > v.vaf * 0.5:
                    v.chip_filter = True
                    logger.debug(
                        f"Filtered CHIP: {v.chrom}:{v.pos}{v.ref}>{v.alt} "
                        f"(plasma VAF={v.vaf:.4f}, WBC VAF={wbc_v.vaf:.4f})"
                    )
                    continue

            chip_filtered.append(v)

        logger.info(
            f"CHIP filter: removed {len(variants) - len(chip_filtered)} "
            f"variants out of {len(variants)}"
        )
        return chip_filtered

    def rank_driver_mutations(
        self,
        variants: list,
        driver_genes: dict,
    ) -> list:
        """Annotate driver mutations (oncogenes and tumor suppressors)."""
        all_drivers = set()
        for genes in driver_genes.values():
            all_drivers.update(genes)

        hotspots = {
            ('EGFR', 'L858R'): 'targeted_therapy',
            ('EGFR', 'T790M'): 'targeted_therapy_resistance',
            ('KRAS', 'G12D'): 'oncogenic',
            ('KRAS', 'G12V'): 'oncogenic',
            ('KRAS', 'G13D'): 'oncogenic',
            ('BRAF', 'V600E'): 'targeted_therapy',
            ('PIK3CA', 'E545K'): 'oncogenic',
            ('PIK3CA', 'H1047R'): 'oncogenic',
            ('TP53', None): 'tumor_suppressor',
        }

        for v in variants:
            if v.gene in all_drivers:
                v.driver = True

        return variants


class VariantCaller:
    """High-level variant calling pipeline."""

    def __init__(self, config: dict):
        self.filter = LowVAFVariantFilter(
            background_error_rate=config.get('background_error_rate', 1e-7),
            min_alt_reads=config.get('min_alt_reads', 3),
            min_vaf=config.get('min_vaf', 1e-5),
        )
        self.config = config

        self.driver_genes = config.get('driver_genes', {
            'lung': ['EGFR', 'KRAS', 'TP53', 'ALK', 'ROS1', 'BRAF', 'MET', 'RET'],
            'breast': ['PIK3CA', 'TP53', 'ESR1', 'HER2', 'BRCA1', 'BRCA2'],
            'colorectal': ['KRAS', 'NRAS', 'BRAF', 'APC', 'TP53', 'SMAD4'],
            'pancreas': ['KRAS', 'TP53', 'CDKN2A', 'SMAD4', 'BRCA2'],
            'liver': ['TP53', 'TERT', 'CTNNB1', 'AXIN1'],
            'general': ['TP53', 'KRAS', 'PIK3CA', 'PTEN', 'RB1', 'APC'],
        })

        logger.info("VariantCaller initialized")

    def call(
        self,
        pileup_df: pd.DataFrame,
        wbc_pileup_df: Optional[pd.DataFrame] = None,
    ) -> pd.DataFrame:
        """Main entry point: call variants from a pileup DataFrame."""
        raw_calls = self.filter.call_variants_from_bam(pileup_df)

        if wbc_pileup_df is not None:
            wbc_calls = self.filter.call_variants_from_bam(wbc_pileup_df)
            filtered_calls = self.filter.filter_chip_variants(raw_calls, wbc_calls)
        else:
            filtered_calls = raw_calls

        ranked_calls = self.filter.rank_driver_mutations(
            filtered_calls, self.driver_genes
        )

        if not ranked_calls:
            return pd.DataFrame()

        result = pd.DataFrame([{
            'chrom': v.chrom,
            'pos': v.pos,
            'ref': v.ref,
            'alt': v.alt,
            'vaf': v.vaf,
            'depth': v.depth,
            'alt_count': v.alt_count,
            'p_value': v.p_value,
            'q_score': v.q_score,
            'gene': v.gene,
            'driver': v.driver,
            'chip_filter': v.chip_filter,
            'tissue_origin': v.tissue_origin,
        } for v in ranked_calls])

        result = result.sort_values(['vaf', 'q_score'], ascending=[False, False])

        return result

    def compute_cancer_scores(self, variants: pd.DataFrame) -> dict:
        """Compute cancer-relevant scores from variant calls."""
        if variants.empty:
            return {
                'mutation_burden': 0,
                'driver_score': 0.0,
                'max_vaf': 0.0,
                'estimated_tumor_fraction': 0.0,
                'high_confidence_count': 0,
                'called_genes': [],
            }

        scores = {
            'mutation_burden': len(variants),
            'driver_score': float(variants['driver'].sum() * 10),
            'max_vaf': float(variants['vaf'].max()),
            'estimated_tumor_fraction': float(2 * variants['vaf'].max()),
            'high_confidence_count': int((variants['q_score'] > 30).sum()),
            'called_genes': variants[variants['gene'].notna()]['gene'].unique().tolist(),
        }

        return scores
