#!/usr/bin/env python3
"""
CHIP (Clonal Hematopoiesis of Indeterminate Potential) Filter

Detects and removes somatic variants likely originating from white blood
cells using population-frequency databases and read-backed phasing.
"""
import numpy as np
from typing import Dict, List, Set, Tuple

# Known CHIP-associated genes (Genovese 2014 NEJM, Jaiswal 2014 NEJM)
CHIP_GENES = {'DNMT3A', 'TET2', 'ASXL1', 'TP53', 'JAK2', 'SF3B1', 
              'SRSF2', 'PPM1D', 'GNB1', 'CBL', 'IDH2', 'U2AF1',
              'ZRSR2', 'EZH2', 'ETV6', 'RUNX1', 'GNAS', 'CUX1'}

# Population prevalence thresholds (gnomAD v4)
CHIP_MAX_POP_AF = 0.01  # Variants above this in gnomAD are likely germline


class CHIPFilter:
    """
    Filter for clonal hematopoiesis variants.
    
    Removes variants that are:
    1. In known CHIP-associated genes
    2. At low VAF (0.1-2% typical of CHIP)
    3. Present in population databases above threshold
    4. Without read-backed phasing evidence for tumor origin
    """
    
    def __init__(
        self,
        chip_genes: Set[str] = CHIP_GENES,
        max_pop_af: float = CHIP_MAX_POP_AF,
        chip_vaf_range: Tuple[float, float] = (0.001, 0.02)
    ):
        self.chip_genes = chip_genes
        self.max_pop_af = max_pop_af
        self.chip_vaf_min, self.chip_vaf_max = chip_vaf_range
    
    def is_chip_variant(self, variant: Dict) -> bool:
        """
        Check if a variant is likely from CHIP.
        
        Parameters
        ----------
        variant : dict
            Keys: gene, vaf, population_af, phased_with_germline
        
        Returns
        -------
        bool
            True if variant is likely CHIP-derived.
        """
        gene = variant.get('gene', '')
        vaf = variant.get('vaf', 0)
        pop_af = variant.get('population_af', 0)
        phased_germline = variant.get('phased_with_germline', False)
        
        # Check 1: Known CHIP gene
        if gene not in self.chip_genes:
            return False
        
        # Check 2: CHIP-typical VAF range
        if vaf < self.chip_vaf_min or vaf > self.chip_vaf_max:
            return False
        
        # Check 3: Population frequency
        if pop_af > self.max_pop_af:
            return True
        
        # Check 4: Phased with germline → likely CHIP
        if phased_germline:
            return True
        
        return False
    
    def filter(self, variants: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
        """
        Split variants into tumor-derived and CHIP-derived.
        
        Returns (tumor_variants, chip_variants).
        """
        tumor = []
        chip = []
        for v in variants:
            if self.is_chip_variant(v):
                chip.append(v)
            else:
                tumor.append(v)
        return tumor, chip


class NanoparticleEnrichmentSimulator:
    """
    Simulate nanoparticle-based ctDNA enrichment.
    
    Models the effect of selective ctDNA enrichment using:
    - Gold nanoparticles (AuNP): size-selective cfDNA capture
    - Magnetic nano-electrodes: charge-based tumor cfDNA separation
    
    Tumor cfDNA is shorter and has different methylation → can be
    selectively enriched via physical/chemical nanoparticle properties.
    """
    
    def __init__(
        self,
        enrichment_factor: float = 5.0,
        specificity: float = 0.95
    ):
        """
        Parameters
        ----------
        enrichment_factor : float
            Fold enrichment of ctDNA over background cfDNA.
        specificity : float
            Fraction of enriched material that is true ctDNA.
        """
        self.enrichment_factor = enrichment_factor
        self.specificity = specificity
    
    def simulate(
        self,
        true_ctdna_fraction: float,
        n_genome_equivalents: int = 30000
    ) -> Dict[str, float]:
        """
        Simulate post-enrichment ctDNA fraction.
        
        Parameters
        ----------
        true_ctdna_fraction : float
            Native ctDNA fraction before enrichment.
        n_genome_equivalents : int
            Number of genome equivalents in sample.
        
        Returns
        -------
        dict with enriched_ctdna_fraction, sensitivity_gain, effective_ge
        """
        # Pre-enrichment: true_ctdna_fraction
        # Post-enrichment: enriched by nanoparticle capture
        enriched = true_ctdna_fraction * self.enrichment_factor * self.specificity
        enriched = min(enriched, 1.0)  # Cap at 100%
        
        sensitivity_gain = enriched / (true_ctdna_fraction + 1e-10)
        
        return {
            'native_ctdna_fraction': true_ctdna_fraction,
            'enriched_ctdna_fraction': enriched,
            'sensitivity_gain': sensitivity_gain,
            'effective_genome_equivalents': int(n_genome_equivalents * enriched),
        }


__all__ = [
    "CHIPFilter",
    "CHIP_GENES",
    "CHIP_MAX_POP_AF",
    "NanoparticleEnrichmentSimulator",
]
