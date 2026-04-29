#!/usr/bin/env python3
"""
THEMIS-Inspired cfDNA Fragment Feature Extractors

Implements the four core features from the THEMIS gastric cancer detection
framework (Bie et al. 2023, Nature Communications):

1. MFR — Methylated Fragment Ratio
2. FSI — Fragment Size Index  
3. CAFF — Cancer-associated Fragmental Aneuploidy
4. FEM — Fragment End Motif analysis
"""
import numpy as np
from typing import Dict, Tuple, List, Optional
from collections import Counter


class MFRCalculator:
    """
    Methylated Fragment Ratio (MFR) calculator.
    
    For 1-Mb genome bins, computes the fraction of fully methylated fragments.
    Tumor-derived cfDNA shows aberrant methylation patterns detectable via
    fragment-level methylation density.
    
    Parameters
    ----------
    bin_size : int
        Genome bin size in base pairs (default 1_000_000 for 1-Mb bins).
    """
    
    def __init__(self, bin_size: int = 1_000_000):
        self.bin_size = bin_size
    
    def compute(
        self,
        fragments: List[Dict],
        genome_length: int = 3_000_000_000
    ) -> np.ndarray:
        """
        Compute MFR per 1-Mb genome bin.
        
        Returns (n_bins,) array of methylation ratios.
        """
        n_bins = genome_length // self.bin_size
        methylated = np.zeros(n_bins)
        total = np.zeros(n_bins)
        
        for frag in fragments:
            bin_idx = frag['start'] // self.bin_size
            if bin_idx < n_bins:
                total[bin_idx] += 1
                if frag.get('methylated', False):
                    methylated[bin_idx] += 1
        
        with np.errstate(divide='ignore', invalid='ignore'):
            mfr = np.where(total > 0, methylated / total, 0.0)
        return mfr
    
    def aggregate_score(self, mfr: np.ndarray) -> float:
        """Aggregate bin-level MFR to a single diagnostic score."""
        return float(np.mean(mfr))


class FSICalculator:
    """
    Fragment Size Index (FSI) calculator.
    
    Ratio of short (100-166 bp) to long (169-240 bp) fragments.
    Tumor-derived cfDNA tends to be shorter due to altered nucleosome
    positioning and increased nuclease accessibility.
    
    Reference: Cristiano et al. 2019, Nature 570:385-389
    """
    
    SHORT_MIN = 100
    SHORT_MAX = 166
    LONG_MIN = 169
    LONG_MAX = 240
    
    @classmethod
    def compute(cls, fragment_lengths: np.ndarray) -> Dict[str, float]:
        """
        Compute FSI and derived statistics.
        
        Returns dict with fsi, short_count, long_count, short_frac, long_frac.
        """
        short = np.sum((fragment_lengths >= cls.SHORT_MIN) & 
                       (fragment_lengths <= cls.SHORT_MAX))
        long = np.sum((fragment_lengths >= cls.LONG_MIN) & 
                      (fragment_lengths <= cls.LONG_MAX))
        
        fsi = long / (short + 1e-6)
        total = len(fragment_lengths)
        
        return {
            'fsi': float(fsi),
            'short_count': int(short),
            'long_count': int(long),
            'short_frac': float(short / total if total > 0 else 0),
            'long_frac': float(long / total if total > 0 else 0),
        }


class CAFFCalculator:
    """
    Cancer-associated Fragmental Aneuploidy (CAFF) calculator.
    
    Identifies the 5 most aberrant chromosome arms and computes a
    copy-number alteration score. Cancer genomes exhibit chromosomal
    instability detectable through fragment coverage deviations.
    
    Reference: Bie et al. 2023, Nature Communications
    """
    
    CHROM_ARM_BOUNDARIES = {
        '1p': (0, 121_000_000), '1q': (121_000_000, 248_956_422),
        '2p': (0, 93_000_000), '2q': (93_000_000, 242_193_529),
        '3p': (0, 91_000_000), '3q': (91_000_000, 198_295_559),
        '4p': (0, 50_000_000), '4q': (50_000_000, 190_214_555),
        '5p': (0, 48_000_000), '5q': (48_000_000, 181_538_259),
        '6p': (0, 61_000_000), '6q': (61_000_000, 170_805_979),
        '7p': (0, 59_000_000), '7q': (59_000_000, 159_345_973),
        '8p': (0, 45_000_000), '8q': (45_000_000, 145_138_636),
        '9p': (0, 49_000_000), '9q': (49_000_000, 138_394_717),
        '10p': (0, 40_000_000), '10q': (40_000_000, 133_797_422),
        '11p': (0, 53_000_000), '11q': (53_000_000, 135_086_622),
        '12p': (0, 35_000_000), '12q': (35_000_000, 133_275_309),
        '13q': (0, 114_364_328),
        '14q': (0, 107_043_718),
        '15q': (0, 101_991_189),
        '16p': (0, 36_000_000), '16q': (36_000_000, 90_338_345),
        '17p': (0, 22_000_000), '17q': (22_000_000, 83_257_441),
        '18p': (0, 15_000_000), '18q': (15_000_000, 80_373_285),
        '19p': (0, 26_000_000), '19q': (26_000_000, 58_617_616),
        '20p': (0, 27_000_000), '20q': (27_000_000, 64_444_167),
        '21q': (0, 46_709_983),
        '22q': (0, 50_818_468),
    }
    
    def __init__(self, n_top_arms: int = 5):
        self.n_top_arms = n_top_arms
    
    def compute(
        self,
        per_arm_coverage: Dict[str, float],
        expected_coverage: float = 1.0
    ) -> Dict[str, float]:
        """
        Compute CAFF score from per-chromosome-arm coverage.
        
        Returns dict with caff_score, aberrant_arms, per_arm_z_scores.
        """
        z_scores = {}
        for arm, cov in per_arm_coverage.items():
            z_scores[arm] = abs(cov - expected_coverage) / (expected_coverage + 1e-6)
        
        sorted_arms = sorted(z_scores.items(), key=lambda x: x[1], reverse=True)
        top_arms = sorted_arms[:self.n_top_arms]
        
        caff_score = np.mean([z for _, z in top_arms])
        
        return {
            'caff_score': float(caff_score),
            'aberrant_arms': [arm for arm, _ in top_arms],
            'per_arm_z_scores': {arm: float(z) for arm, z in top_arms},
        }


class FEMCalculator:
    """
    Fragment End Motif (FEM) calculator.
    
    Analyzes 4-mer frequencies at fragment ends, particularly for fragments
    <171 bp. Aberrant nuclease cleavage patterns in cancer produce distinct
    end-motif signatures.
    
    Reference: Jiang et al. 2020, Nature Genetics 52:712-719
    """
    
    BASES = ['A', 'C', 'G', 'T']
    
    @classmethod
    def _generate_all_4mers(cls) -> List[str]:
        return [a+b+c+d for a in cls.BASES for b in cls.BASES 
                for c in cls.BASES for d in cls.BASES]
    
    def __init__(self):
        self.all_4mers = self._generate_all_4mers()
        self.motif_to_idx = {m: i for i, m in enumerate(self.all_4mers)}
    
    def compute(
        self,
        end_sequences: List[str],
        fragment_lengths: Optional[np.ndarray] = None,
        short_threshold: int = 171
    ) -> Dict[str, np.ndarray]:
        """
        Compute FEM frequencies with optional short-fragment filtering.
        
        Returns dict with motif_counts (256,), frequencies (256,), mds, short_mds.
        """
        n_motifs = len(self.all_4mers)
        motif_counts = np.zeros(n_motifs, dtype=np.int64)
        short_counts = np.zeros(n_motifs, dtype=np.int64)
        
        for i, seq in enumerate(end_sequences):
            if len(seq) >= 4:
                motif = seq[:4].upper()
                idx = self.motif_to_idx.get(motif)
                if idx is not None:
                    motif_counts[idx] += 1
                    if fragment_lengths is not None and i < len(fragment_lengths):
                        if fragment_lengths[i] < short_threshold:
                            short_counts[idx] += 1
        
        total = motif_counts.sum()
        short_total = short_counts.sum()
        
        frequencies = motif_counts / total if total > 0 else np.zeros(n_motifs)
        mds = self._compute_mds(motif_counts)
        short_mds = self._compute_mds(short_counts) if short_total > 0 else 0.0
        
        return {
            'motif_counts': motif_counts,
            'frequencies': frequencies,
            'mds': float(mds),
            'short_mds': float(short_mds),
        }
    
    @staticmethod
    def _compute_mds(counts: np.ndarray) -> float:
        """Motif Diversity Score — normalized Simpson diversity."""
        n = len(counts)
        total = counts.sum()
        if total == 0:
            return 0.0
        p = counts / total
        simpson = np.sum(p ** 2)
        return (1 - simpson) / (1 - 1.0 / n)


__all__ = [
    "MFRCalculator",
    "FSICalculator",
    "CAFFCalculator",
    "FEMCalculator",
]
