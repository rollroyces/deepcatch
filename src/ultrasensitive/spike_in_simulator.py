"""
Spike-in simulator for validating 0.001% VAF detection sensitivity.

Ported from ``cancer-screening/src/utils/spike_in_simulator.py`` (SHA 72f878c3).

Generates synthetic cfDNA samples by mixing:
  1. Healthy donor cfDNA (background)
  2. Tumor DNA at known, precise dilutions
  3. Sequencer noise (modeled from real error profiles)

Gold-standard experiment for validating Limit of Detection (LOD).

Usage:
  python -m src.ultrasensitive.spike_in_simulator --vaf 0.001 --n_samples 100
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import argparse
import logging
from pathlib import Path
from dataclasses import dataclass

logger = logging.getLogger(__name__)

# Cancer hotspot mutations for spike-in (hg38 coordinates)
HOTSPOT_MUTATIONS = {
    'EGFR': [
        ('L858R', 'CTG>TTG', 0.15),
        ('T790M', 'ACC>ATC', 0.10),
        ('Exon19 deletion', 'TTAAG>---', 0.20),
        ('G719X', 'GGC>GTC', 0.05),
    ],
    'KRAS': [
        ('G12D', 'GGT>GAT', 0.25),
        ('G12V', 'GGT>GTT', 0.20),
        ('G13D', 'GGC>GAC', 0.08),
        ('Q61H', 'CAA>CAT', 0.05),
    ],
    'TP53': [
        ('R248Q', 'CGC>CAC', 0.06),
        ('R273H', 'CGT>CAT', 0.05),
        ('R175H', 'CGC>CAC', 0.04),
    ],
    'PIK3CA': [
        ('E545K', 'GAG>AAG', 0.04),
        ('H1047R', 'CAT>CGT', 0.06),
    ],
    'BRAF': [
        ('V600E', 'GTG>GAG', 0.08),
    ],
    'PTEN': [
        ('R130*', 'CGA>TGA', 0.02),
    ],
}


@dataclass
class SpikeInConfig:
    """Configuration for a spike-in experiment."""
    tumor_fraction: float
    genome_equivalents: int
    n_samples: int
    noise_profile: str = 'illumina'
    chip_contamination: float = 0.0


class SpikeInSimulator:
    """
    Simulate spike-in experiments for LOD validation.

    Generates synthetic cfDNA samples by spiking tumor mutations
    into healthy cfDNA background at precise dilutions.

    Example: 0.001% VAF
      - 50,000 genome equivalents (10 mL blood)
      - 0.001% VAF = 0.5 mutant molecules
      - True VAF = 1/50,000 = 0.002% (after rounding to 1 molecule)

    Statistical framework:
      - Poisson model for molecule counting
      - Binomial model for variant detection
      - Power analysis for required n_samples
    """

    GENOME_EQUIVALENTS_PER_ML_PLASMA = 5000
    BASES_PER_GENOME = 3_200_000_000
    REGION_SIZE_FOR_VAF = 100

    def __init__(self, config: SpikeInConfig):
        self.config = config
        np.random.seed(int(config.tumor_fraction * 1e6 + config.n_samples))

        logger.info(
            f"SpikeInSimulator: VAF={config.tumor_fraction:.1e}, "
            f"GE={config.genome_equivalents:,}, n={config.n_samples}"
        )

    def compute_molecule_counts(self) -> dict:
        """Compute expected molecule counts from VAF."""
        ge = self.config.genome_equivalents

        region_molecules = ge * self.REGION_SIZE_FOR_VAF / self.BASES_PER_GENOME
        expected_mutant = region_molecules * self.config.tumor_fraction
        expected_ref = region_molecules - expected_mutant

        return {
            'genome_equivalents': ge,
            'region_molecules': region_molecules,
            'expected_mutant': expected_mutant,
            'expected_ref': expected_ref,
            'poisson_detected': int(np.random.poisson(expected_mutant)) if expected_mutant > 0 else 0,
        }

    def generate_sample(self) -> dict:
        """Generate one simulated sample with spike-in mutations."""
        molecules = self.compute_molecule_counts()
        reads_per_molecule = 3
        total_reads = int(molecules['region_molecules'] * reads_per_molecule)

        noise_rates = {
            'illumina': 1e-4,
            'nanopore': 1e-3,
            'pacbio': 1e-2,
        }
        bg_error = noise_rates.get(self.config.noise_profile, 1e-4)

        results = {
            'sample_id': None,
            'total_reads': total_reads,
            'genome_equivalents': self.config.genome_equivalents,
            'target_vaf': self.config.tumor_fraction,
            'mutations': {},
        }

        for gene, mutations in HOTSPOT_MUTATIONS.items():
            for mutation_name, _, allele_freq_in_cancer in mutations:
                this_mutation_fraction = self.config.tumor_fraction * allele_freq_in_cancer
                n_molecules = self.compute_molecule_counts()

                mutant_reads = np.random.poisson(
                    n_molecules['region_molecules'] * this_mutation_fraction
                )

                ref_reads = int(n_molecules['region_molecules']) - mutant_reads

                n_bases = self.REGION_SIZE_FOR_VAF
                false_alt_reads = np.random.binomial(ref_reads, bg_error)

                observed_alt = mutant_reads + false_alt_reads
                observed_ref = ref_reads - false_alt_reads
                total = observed_alt + observed_ref

                observed_vaf = float(observed_alt) / float(total) if total > 0 else 0.0

                results['mutations'][f'{gene}_{mutation_name}'] = {
                    'mutant_reads': int(mutant_reads),
                    'ref_reads': int(observed_ref),
                    'false_alt_reads': int(false_alt_reads),
                    'observed_alt_reads': int(observed_alt),
                    'observed_vaf': observed_vaf,
                    'target_vaf': float(this_mutation_fraction),
                    'called': mutant_reads >= 3 and observed_vaf > self.config.tumor_fraction * 0.5,
                    'chrom': self._get_chrom(gene),
                    'position': self._get_position(gene),
                }

        return results

    def run_spike_in_experiment(self) -> tuple:
        """Run the full spike-in experiment with n_samples replicates."""
        all_results = []

        for i in range(self.config.n_samples):
            sample = self.generate_sample()
            sample['sample_id'] = f"S{i+1:03d}"
            sample['replicate'] = i + 1

            n_mutations = len(sample['mutations'])
            n_called = sum(1 for m in sample['mutations'].values() if m['called'])
            mean_vaf = float(np.mean([m['observed_vaf'] for m in sample['mutations'].values()]))
            max_vaf = float(np.max([m['observed_vaf'] for m in sample['mutations'].values()]))

            sample['n_mutations'] = n_mutations
            sample['n_detected'] = n_called
            sample['detection_rate'] = n_called / max(n_mutations, 1)
            sample['mean_observed_vaf'] = mean_vaf
            sample['max_observed_vaf'] = max_vaf

            all_results.append(sample)

        summary = pd.DataFrame([{
            'sample_id': r['sample_id'],
            'target_vaf': r['target_vaf'],
            'n_mutations': r['n_mutations'],
            'n_detected': r['n_detected'],
            'detection_rate': r['detection_rate'],
            'mean_observed_vaf': r['mean_observed_vaf'],
            'max_observed_vaf': r['max_observed_vaf'],
        } for r in all_results])

        return summary, all_results

    def power_analysis(self) -> dict:
        """Power analysis: how many samples needed to detect at target VAF?"""
        vaf = self.config.tumor_fraction
        ge = self.config.genome_equivalents
        bg_rate = 1e-7  # After UMI correction

        expected_mutant = ge * vaf

        from scipy.stats import poisson

        def detection_prob(lam: float) -> float:
            if lam <= 0:
                return 0.0
            return 1 - sum(poisson.pmf(k, lam) for k in range(3))

        power = detection_prob(expected_mutant)

        if power > 0 and power < 1:
            n_for_80pct = float(np.log(0.2) / np.log(1 - power))
        else:
            n_for_80pct = float('inf')

        return {
            'target_vaf': vaf,
            'expected_mutant_molecules': expected_mutant,
            'single_sample_power': float(power),
            'samples_for_80pct_power': int(np.ceil(n_for_80pct)) if n_for_80pct != float('inf') else None,
            'total_reads_needed': int(
                self.config.genome_equivalents * 30 / self.BASES_PER_GENOME * self.REGION_SIZE_FOR_VAF
            ),
        }

    @staticmethod
    def _get_chrom(gene: str) -> str:
        chromosomes = {
            'EGFR': 'chr7', 'KRAS': 'chr12', 'TP53': 'chr17',
            'PIK3CA': 'chr3', 'BRAF': 'chr7', 'PTEN': 'chr10',
        }
        return chromosomes.get(gene, 'chr17')

    @staticmethod
    def _get_position(gene: str) -> int:
        positions = {
            'EGFR': 55_245_000,
            'KRAS': 25_218_000,
            'TP53': 7_677_000,
            'PIK3CA': 179_240_000,
            'BRAF': 140_719_000,
            'PTEN': 89_692_000,
        }
        return positions.get(gene, 7_677_000)


def run_lod_validation(vaf_range: list, n_samples: int = 20, out_csv: str = None) -> pd.DataFrame:
    """
    Run Limit of Detection validation across a range of VAF levels.

    Args:
        vaf_range: List of target VAFs to test (e.g. ``[1e-2, 1e-3, 1e-4]``).
        n_samples: Number of technical replicates per VAF level.
        out_csv: Optional output CSV path for committed ``results/`` artifacts.

    Returns:
        DataFrame summarizing detection rate, statistical power, and LOD per VAF.
    """
    results = []

    for vaf in vaf_range:
        config = SpikeInConfig(
            tumor_fraction=vaf,
            genome_equivalents=50_000,
            n_samples=n_samples,
        )
        sim = SpikeInSimulator(config)
        summary, _ = sim.run_spike_in_experiment()
        power = sim.power_analysis()

        detection_rate = float(summary['detection_rate'].mean())
        power_estimate = float(power['single_sample_power'])

        results.append({
            'target_vaf': vaf,
            'expected_mutant_molecules': power['expected_mutant_molecules'],
            'detection_rate_mean': detection_rate,
            'detection_rate_std': float(summary['detection_rate'].std()),
            'single_sample_power': power_estimate,
            'samples_for_80pct': power['samples_for_80pct_power'],
            'lod_achieved': detection_rate >= 0.95,
        })

        logger.info(
            f"VAF={vaf:.1e}: detection={detection_rate:.1%}, "
            f"power={power_estimate:.1%}, "
            f"n_for_80pct={power['samples_for_80pct_power']}"
        )

    df = pd.DataFrame(results)
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        logger.info(f"Wrote LOD validation results to {out_csv}")
    return df


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    vaf_range = [1e-2, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6]

    print("=" * 70)
    print("Limit of Detection (LOD) Validation — CANCER-SEEK Spike-In Experiment")
    print("=" * 70)

    results = run_lod_validation(vaf_range, n_samples=20)

    print("\n\nLOD VALIDATION SUMMARY")
    print("=" * 70)
    print(results.to_string(index=False))

    achieved_lod = float(results[results['lod_achieved']]['target_vaf'].max())
    print(f"\nConfirmed LOD (>=95% detection): {achieved_lod:.1e}")
