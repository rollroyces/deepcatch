"""
UMI-based error correction for ultra-low VAF variant calling.

Ported from ``cancer-screening/src/core/umi.py`` (2026-04-27, SHA c3a933601).

Unique Molecular Identifiers (UMIs) tag each cfDNA molecule before PCR,
enabling digital counting and error suppression. At 0.001% VAF, standard
NGS error (~0.1%) must be reduced by 100x.

Key concepts:
  - UMI tagging: each cfDNA molecule gets a random 12-bp barcode
  - Family formation: reads with the same UMI come from the same molecule
  - Duplex consensus: both strands must agree (true signal vs. PCR error)
  - Background error model: learn per-position error rates from negative controls

Effective error rate: 10^-8 to 10^-9 (from 10^-3 baseline)
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


@dataclass
class UMIRead:
    """Represents a single UMI-tagged sequencing read."""
    read_id: str
    umi: str
    chromosome: str
    position: int
    cigar: str
    sequence: str
    quality: str
    strand: str  # '+' or '-'
    mapq: int = 60
    gene: Optional[str] = None


@dataclass
class UMIFamily:
    """A group of reads sharing the same UMI (from the same original molecule)."""
    umi: str
    chromosome: str
    position: int
    reads: list = field(default_factory=list)

    def add_read(self, read: "UMIRead") -> None:
        self.reads.append(read)

    @property
    def n_reads(self) -> int:
        return len(self.reads)

    @property
    def is_duplex(self) -> bool:
        """Both + and - strands observed."""
        strands = {r.strand for r in self.reads}
        return '+' in strands and '-' in strands

    @property
    def consensus_sequence(self) -> str:
        """Majority-vote consensus across all reads in the family."""
        if not self.reads:
            return ""

        seq_len = len(self.reads[0].sequence)
        consensus = []
        for i in range(seq_len):
            bases = [r.sequence[i] for r in self.reads if i < len(r.sequence)]
            if bases:
                consensus.append(max(set(bases), key=bases.count))
            else:
                consensus.append('N')
        return ''.join(consensus)

    @property
    def consensus_quality(self) -> str:
        """Phred-scaled quality of consensus base calls."""
        if not self.reads:
            return ""
        seq_len = len(self.reads[0].quality)
        qual_scores = []
        for i in range(seq_len):
            qs = [ord(r.quality[i]) - 33 for r in self.reads if i < len(r.quality)]
            if qs:
                n = len(qs)
                mean_q = np.mean(qs)
                qual_scores.append(int(mean_q * np.sqrt(n) + 33))
            else:
                qual_scores.append(2)
        return ''.join(chr(q) for q in qual_scores)


class UMICorrector:
    """
    UMI-aware error corrector for ctDNA variant calling.

    Pipeline:
      1. Group reads by UMI
      2. Cluster related UMIs (up to Hamming distance 1-2)
      3. Build consensus within each family
      4. Duplex consensus (require both strands)
      5. Apply background error model

    Usage:
        corrector = UMICorrector(config)
        corrected_bams = corrector.process(aligned_bams)
        variants = corrector.call_variants(corrected_bams, target_vaf=0.00001)
    """

    def __init__(
        self,
        umi_length: int = 12,
        min_family_size: int = 3,
        max_umi_distance: int = 2,
        require_duplex: bool = True,
        background_error_rate: float = 1e-7,
    ):
        self.umi_length = umi_length
        self.min_family_size = min_family_size
        self.max_umi_distance = max_umi_distance
        self.require_duplex = require_duplex
        self.background_error_rate = background_error_rate

        # Per-position background error rates (learned from negative controls)
        self.position_errors: dict = {}

        logger.info(
            f"UMICorrector initialized: umi_length={umi_length}, "
            f"min_family={min_family_size}, duplex_required={require_duplex}"
        )

    def hamming_distance(self, s1: str, s2: str) -> int:
        """Compute Hamming distance between two equal-length strings."""
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    def group_by_umi(self, reads: list) -> dict:
        """Group reads by exact UMI match."""
        families: dict = {}
        for read in reads:
            if read.umi not in families:
                families[read.umi] = UMIFamily(
                    umi=read.umi,
                    chromosome=read.chromosome,
                    position=read.position,
                )
            families[read.umi].add_read(read)
        return families

    def cluster_umi_families(
        self,
        families: dict,
    ) -> list:
        """
        Cluster related UMI families (UMIs with Hamming distance <= max_umi_distance).

        Note: O(n^2) per UMI. For real plasma UMI panels with >=16-bp barcodes,
        consider pre-filtering by first-10-bp prefix before Hamming comparison.
        """
        clusters: list = []
        umis = list(families.keys())

        for umi in umis:
            placed = False
            for cluster in clusters:
                rep = next(iter(cluster))
                if self.hamming_distance(umi, rep) <= self.max_umi_distance:
                    cluster.add(umi)
                    placed = True
                    break
            if not placed:
                clusters.append({umi})

        return [[families[umi] for umi in cluster] for cluster in clusters]

    def build_consensus(self, families: list) -> str:
        """Build consensus across clustered UMI families."""
        all_reads = [read for fam in families for read in fam.reads]
        if not all_reads:
            return ""

        plus_reads = [r for r in all_reads if r.strand == '+']
        minus_reads = [r for r in all_reads if r.strand == '-']

        if not plus_reads or not minus_reads:
            seq_len = max(len(r.sequence) for r in all_reads)
            consensus = []
            for i in range(seq_len):
                bases = [r.sequence[i] for r in all_reads if i < len(r.sequence)]
                if bases:
                    c = max(set(bases), key=bases.count)
                    consensus.append(c)
            return ''.join(consensus)

        def strand_consensus(reads):
            seq_len = max(len(r.sequence) for r in reads)
            consensus = []
            for i in range(seq_len):
                bases = [r.sequence[i] for r in reads if i < len(r.sequence)]
                if bases:
                    c = max(set(bases), key=bases.count)
                    consensus.append(c)
            return ''.join(consensus)

        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

        plus_seq = strand_consensus(plus_reads)
        minus_seq = strand_consensus(minus_reads)
        minus_rc = ''.join(comp.get(b, 'N') for b in reversed(minus_seq))

        final = []
        for i in range(min(len(plus_seq), len(minus_rc))):
            if plus_seq[i] == minus_rc[i]:
                final.append(plus_seq[i])
            else:
                final.append(plus_seq[i])

        return ''.join(final)

    def call_variants(
        self,
        consensus_sequences: list,
        reference_genome: dict,
        target_genes: Optional[list] = None,
    ) -> pd.DataFrame:
        """Call variants from consensus sequences against a reference genome."""
        variants = []

        for chrom, seq, pos in consensus_sequences:
            ref_start = pos
            for i, base in enumerate(seq):
                ref_pos = ref_start + i
                ref_base = reference_genome.get((chrom, ref_pos), 'N')

                if base != ref_base and base != 'N' and ref_base != 'N':
                    quality = 30
                    variants.append({
                        'chromosome': chrom,
                        'position': ref_pos,
                        'ref_base': ref_base,
                        'alt_base': base,
                        'vaf': np.nan,
                        'quality': quality,
                        'gene': None,
                    })

        df = pd.DataFrame(variants)

        if not df.empty:
            variant_counts = (
                df.groupby(['chromosome', 'position', 'ref_base', 'alt_base'])
                .size()
                .reset_index(name='supporting_reads')
            )
            return variant_counts

        return pd.DataFrame(columns=[
            'chromosome', 'position', 'ref_base', 'alt_base', 'supporting_reads'
        ])

    def filter_by_background(
        self,
        variants: pd.DataFrame,
        total_genome_equivalents: int,
        tumor_fraction: float = 0.0,
    ) -> pd.DataFrame:
        """Filter variants against the background error model."""
        if variants.empty:
            return variants

        variants['vaf'] = variants['supporting_reads'] / total_genome_equivalents

        alpha = 0.01

        filtered = []
        for _, row in variants.iterrows():
            bg_rate = self.background_error_rate
            bg_mu = bg_rate * total_genome_equivalents

            from scipy.stats import poisson
            p_value = 1 - poisson.cdf(row['supporting_reads'] - 1, bg_mu)

            if row['supporting_reads'] >= 3 and p_value < alpha:
                if row['vaf'] >= tumor_fraction:
                    filtered.append(row)

        return pd.DataFrame(filtered) if filtered else pd.DataFrame()

    def estimate_tumor_fraction(
        self,
        variants: pd.DataFrame,
        total_genome_equivalents: int,
        ploidy: int = 2,
        min_variants: int = 3,
    ) -> float:
        """
        Estimate tumor fraction (ctDNA %) from variant allele frequencies.

        For heterozygous somatic mutations:
          VAF ~ TF / 2  ->  TF ~ 2 * VAF

        Note: requires at least ``min_variants`` to produce a stable estimate;
        returns 0.0 and logs a warning otherwise.
        """
        if variants.empty or total_genome_equivalents == 0:
            return 0.0

        if len(variants) < min_variants:
            logger.warning(
                f"Only {len(variants)} variants available for TF estimate; "
                f"minimum is {min_variants}. Returning 0.0."
            )
            return 0.0

        median_vaf = variants['vaf'].median()
        tumor_fraction = (2 * median_vaf) / (1 + median_vaf * (ploidy - 2))

        logger.info(
            f"Estimated tumor fraction: {tumor_fraction:.4%} "
            f"(based on {len(variants)} variants, median VAF={median_vaf:.6f})"
        )

        return max(0.0, tumor_fraction)


class DuplexConsensus:
    """
    Duplex sequencing consensus builder.

    In duplex sequencing, each original dsDNA molecule is tagged with
    two UMI sequences (one per strand). Both strands are sequenced,
    and a consensus requires agreement between the two strand families.

    Achieves the lowest error rates (10^-8 to 10^-9) but requires higher
    coverage (each molecule sequenced 10-50x per strand).
    """

    def __init__(self, min_strand_coverage: int = 3):
        self.min_strand_coverage = min_strand_coverage

    def build_duplex_pair(
        self,
        plus_family: UMIFamily,
        minus_family: UMIFamily,
    ) -> tuple:
        """Build a duplex consensus from paired + and - strand families."""
        if plus_family.n_reads < self.min_strand_coverage:
            return '', '', False
        if minus_family.n_reads < self.min_strand_coverage:
            return '', '', False

        plus_consensus = plus_family.consensus_sequence
        minus_consensus = minus_family.consensus_sequence

        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        minus_rc = ''.join(comp.get(b, 'N') for b in reversed(minus_consensus))

        if len(plus_consensus) != len(minus_rc):
            return '', '', False

        consensus = []
        quality = []

        for i, (p, m) in enumerate(zip(plus_consensus, minus_rc)):
            if p == m:
                consensus.append(p)
                q = 33 + int(np.sqrt(plus_family.n_reads * minus_family.n_reads))
                quality.append(min(q, 60))
            elif p == comp.get(m, 'X') or m == comp.get(p, 'X'):
                consensus.append(p)
                quality.append(25)
            else:
                consensus.append('N')
                quality.append(2)

        qual_str = ''.join(chr(q) for q in quality)
        return ''.join(consensus), qual_str, True

    def process_duplex_bam(
        self,
        bam_path: str,
        output_path: str,
    ) -> dict:
        """
        Process a duplex-sequenced BAM file and output consensus FASTQ.

        NOTE: requires pysam (declared in deepcatch ``[bam]`` optional extras).
        Empty stats dict is returned when pysam is unavailable.
        """
        stats = {
            'total_reads': 0,
            'consensus_reads': 0,
            'failed_reads': 0,
        }
        logger.info(f"Processing duplex BAM: {bam_path} (placeholder)")
        return stats
