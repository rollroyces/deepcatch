"""
Copy Number Variation (CNV) analysis from shallow whole-genome sequencing.

Ported from ``cancer-screening/src/core/cnv_analyzer.py`` (SHA 7dd7d48f).

CNVs (amplifications and deletions) are detectable in ctDNA at tumor fractions
as low as 1-5%. At 0.001% VAF, CNV is supporting evidence rather than primary signal.

This module provides:
  - Normalized coverage computation
  - GC-bias correction
  - Copy number segmentation (CBS algorithm)
  - Tumor fraction estimation from CNV magnitude

Bug fixes applied during port:
  - ``logging.getLogger(__Logger(__name__))`` -> ``logging.getLogger(__name__)``
  - Replaced legacy dict ordering for hotspot tuples (Python 3.7+ keeps insertion
    order, no change needed).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import signal, stats
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class CNVAnalyzer:
    """
    CNV detection from low-pass WGS (sWGS) data.

    Pipeline:
      1. Bin genome into fixed windows (10kb, 50kb, 100kb)
      2. Count reads per bin, normalize by total coverage
      3. GC-content correction (avoids coverage bias)
      4. Copy number segmentation (Circular Binary Segmentation)
      5. Call amplifications (log2 > 0.3) and deletions (log2 < -0.3)

    Sensitivity: ~5% tumor fraction for 100kb bins
    Complementary to: UMI-based SNV calling (which works at 0.001% VAF)
    """

    AMPLIFICATION_HOTSPOTS = {
        'ERBB2': ('chr17', 39_687_000, 39_737_000),
        'MYC':   ('chr8',  127_735_000, 127_788_000),
        'EGFR':  ('chr7',  55_019_000, 55_093_000),
        'CCND1': ('chr11', 69_165_000, 69_192_000),
        'FGFR1': ('chr8',  38_411_000, 38_469_000),
        'MDM2':  ('chr12', 67_536_000, 67_579_000),
        'KRAS':  ('chr12', 25_248_000, 25_355_000),
        'BRAF':  ('chr7',  140_713_000, 140_773_000),
    }

    DELETION_HOTSPOTS = {
        'PTEN':   ('chr10', 89_693_000, 89_763_000),
        'CDKN2A': ('chr9',  21_967_000, 22_029_000),
        'SMAD4':  ('chr18', 46_075_000, 46_163_000),
        'RB1':    ('chr13', 48_378_000, 48_597_000),
        'APC':    ('chr5',  112_700_000, 112_847_000),
        'BRCA1':  ('chr17', 43_044_000, 43_240_000),
        'BRCA2':  ('chr13', 32_315_000, 32_985_000),
    }

    def __init__(
        self,
        bin_size: int = 50_000,
        genome_size_mb: int = 3000,
        min_reads_per_bin: int = 50,
        gc_window: int = 100,
        segmentation_alpha: float = 1e-4,
        nperm: int = 1000,
    ):
        self.bin_size = bin_size
        self.genome_size_mb = genome_size_mb
        self.min_reads_per_bin = min_reads_per_bin
        self.gc_window = gc_window
        self.segmentation_alpha = segmentation_alpha
        self.nperm = nperm

        logger.info(
            f"CNVAnalyzer: bin_size={bin_size//1000}kb, "
            f"segments={genome_size_mb * 1e6 / bin_size:.0f} bins"
        )

    def bin_coverage(
        self,
        bam_path: str,
        chromosomes: Optional[list] = None,
    ) -> pd.DataFrame:
        """Compute read depth per genomic bin (placeholder; uses pysam in production)."""
        if chromosomes is None:
            chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

        bins = []
        for chrom in chromosomes:
            chrom_length = {'chr1': 248_956_422, 'chr2': 242_193_529}.get(
                chrom, 200_000_000
            )
            for start in range(0, chrom_length, self.bin_size):
                bins.append({
                    'chrom': chrom,
                    'start': start,
                    'end': min(start + self.bin_size, chrom_length),
                    'read_count': 0,
                    'gc_content': 0.5,
                })

        return pd.DataFrame(bins)

    def gc_correction(self, coverage: pd.DataFrame) -> pd.Series:
        """Correct read counts for GC-content bias using sliding-window LOESS."""
        gc = coverage['gc_content'].values
        counts = coverage['read_count'].values

        mask = counts > 0
        gc_valid = gc[mask]
        counts_valid = counts[mask]

        if len(gc_valid) < 100:
            return pd.Series(np.ones(len(coverage)), index=coverage.index)

        from scipy.ndimage import uniform_filter1d
        sorted_idx = np.argsort(gc_valid)
        gc_sorted = gc_valid[sorted_idx]
        counts_sorted = counts_valid[sorted_idx]

        window = max(50, len(gc_sorted) // 20)
        expected = uniform_filter1d(counts_sorted, window, mode='nearest')

        expected_full = np.zeros(len(coverage))
        for i, idx in enumerate(sorted_idx):
            expected_full[idx] = expected[i]

        normalized = np.where(expected_full > 0, counts / expected_full, 1.0)
        normalized = np.clip(normalized, 0.1, 10.0)

        return pd.Series(normalized, index=coverage.index)

    def compute_log2_ratio(
        self,
        normalized_coverage: pd.Series,
        ploidy: int = 2,
    ) -> pd.Series:
        """Convert normalized coverage to log2 copy number ratio."""
        expected = normalized_coverage.median()
        log2_ratio = np.log2(normalized_coverage / expected + 1e-10)
        return log2_ratio

    def circular_binary_segmentation(
        self,
        segments: pd.DataFrame,
    ) -> pd.DataFrame:
        """Circular Binary Segmentation (simplified Python approximation)."""
        log2 = segments['log2_ratio'].values
        n = len(log2)

        if n < 10:
            segments['segment_id'] = 0
            return segments

        z = stats.zscore(log2)
        threshold = stats.norm.ppf(1 - self.segmentation_alpha / 2)

        segment_ids = np.zeros(n, dtype=int)
        current_segment = 0

        for i in range(1, n):
            if abs(log2[i] - log2[i-1]) < 0.2:
                segment_ids[i] = current_segment
            else:
                current_segment += 1
                segment_ids[i] = current_segment

        segments['segment_id'] = segment_ids
        segment_means = segments.groupby('segment_id')['log2_ratio'].transform('mean')
        segments['segment_mean'] = segment_means

        return segments

    def call_cnvs(
        self,
        segmented: pd.DataFrame,
        tumor_fraction: float = 0.1,
    ) -> pd.DataFrame:
        """Call copy number alterations from segmented data."""
        threshold = 0.3 * tumor_fraction
        cnv_threshold = max(threshold, 0.2)

        calls = []
        for seg_id, seg_data in segmented.groupby('segment_id'):
            seg_mean = seg_data['segment_mean'].iloc[0]
            seg_len = len(seg_data)

            if seg_mean > cnv_threshold:
                cnv_type = 'amplification'
                cnv_range = seg_mean
            elif seg_mean < -cnv_threshold:
                cnv_type = 'deletion'
                cnv_range = abs(seg_mean)
            else:
                continue

            chrom = seg_data['chrom'].iloc[0]
            start = seg_data['start'].iloc[0]
            end = seg_data['end'].iloc[-1]

            estimated_cn = 2 * (2 ** seg_mean)
            annotation = self._annotate_hotspot(chrom, start, end)

            calls.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'cnv_type': cnv_type,
                'log2_ratio': round(float(seg_mean), 3),
                'estimated_copy_number': round(float(estimated_cn), 1),
                'length_bp': end - start,
                'n_bins': seg_len,
                'hotspot_annotation': annotation,
            })

        df = pd.DataFrame(calls)
        logger.info(f"Called {len(df)} CNVs from {segmented['segment_id'].nunique()} segments")
        return df

    def _annotate_hotspot(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Optional[str]:
        """Check if CNV overlaps a known cancer-related gene."""
        for gene, (ref_chrom, ref_start, ref_end) in self.AMPLIFICATION_HOTSPOTS.items():
            if chrom == ref_chrom and start <= ref_end and end >= ref_start:
                return f"{gene} amplification"

        for gene, (ref_chrom, ref_start, ref_end) in self.DELETION_HOTSPOTS.items():
            if chrom == ref_chrom and start <= ref_end and end >= ref_start:
                return f"{gene} deletion"

        return None

    def estimate_tumor_fraction_from_cnv(
        self,
        cnv_calls: pd.DataFrame,
    ) -> float:
        """Estimate tumor fraction from CNV magnitude."""
        if cnv_calls.empty:
            return 0.0

        deletions = cnv_calls[cnv_calls['cnv_type'] == 'deletion']

        if deletions.empty:
            amplifications = cnv_calls[cnv_calls['cnv_type'] == 'amplification']
            if amplifications.empty:
                return 0.0

            max_log2 = amplifications['log2_ratio'].max()
            tf = (2 ** max_log2) - 1
        else:
            min_log2 = deletions['log2_ratio'].min()
            tf = 2 * (1 - (2 ** min_log2))

        tf = max(0.0, min(1.0, tf))
        logger.info(f"CNV-estimated tumor fraction: {tf:.2%}")
        return tf

    def full_pipeline(self, bam_path: str) -> dict:
        """Run the complete CNV analysis pipeline."""
        logger.info(f"Running CNV pipeline on: {bam_path}")

        coverage = self.bin_coverage(bam_path)
        gc_corrected = self.gc_correction(coverage)
        coverage['gc_corrected'] = gc_corrected

        log2_ratios = self.compute_log2_ratio(gc_corrected)
        coverage['log2_ratio'] = log2_ratios

        segmented = self.circular_binary_segmentation(coverage)
        cnvs = self.call_cnvs(segmented, tumor_fraction=0.1)
        refined_tf = self.estimate_tumor_fraction_from_cnv(cnvs)

        return {
            'coverage': coverage,
            'segments': segmented,
            'cnv_calls': cnvs,
            'estimated_tumor_fraction': refined_tf,
        }
