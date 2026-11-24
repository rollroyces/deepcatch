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
    """Configuration for a spike-in experiment.

    ``sequencing_depth`` is the number of reads covering a target site. It is
    the knob that decides whether a VAF is detectable at all: the simulator
    samples ``depth`` reads at a site and asks how many carry the mutation.
    """
    tumor_fraction: float
    genome_equivalents: int
    n_samples: int
    noise_profile: str = 'illumina'
    chip_contamination: float = 0.0
    sequencing_depth: int = 100_000


class SpikeInSimulator:
    """
    Simulate spike-in experiments for LOD validation.

    Generates synthetic cfDNA samples by spiking tumour mutations into healthy
    cfDNA background at precise dilutions, and asks whether each hotspot would
    be *called* at a given sequencing depth.

    The honest physics this models
    -----------------------------
    At 50,000 genome equivalents (10 mL of plasma) and a VAF of 0.001%, there
    are **0.5 mutant molecules in the entire genome**. No assay can call a
    fraction of a molecule. What a sequencer actually sees at a target site is
    ``sequencing_depth`` reads, of which the expected mutant count is
    ``depth * VAF * allele_frequency_in_tumour``.

    So detection is governed by depth, not by replicate count:

    * 0.001% VAF at 100,000x -> 1 mutant read per site -> below the
      3-supporting-read threshold -> **not called**.
    * 0.001% VAF needs ~300,000x to clear the threshold, and ~450,000x for
      80% power.

    This is why the CANCER-SEEK architecture document specifies 100,000x
    coverage for 0.001% sensitivity and why the UMI error-correction stage
    exists: at this concentration the limiting factor is the number of
    *independent* molecules, not the number of reads.

    Statistical framework:
      - Poisson model for mutant-read counting at a site.
      - Binomial model for platform background error producing false ALT reads.
      - Depth-based power analysis consistent with the calling rule.
    """

    GENOME_EQUIVALENTS_PER_ML_PLASMA = 5000
    BASES_PER_GENOME = 3_200_000_000
    REGION_SIZE_FOR_VAF = 100

    #: Minimum supporting reads required to call a variant.
    MIN_SUPPORTING_READS = 3

    #: Per-base background error rate by platform.
    NOISE_RATES = {
        'illumina': 1e-4,
        'nanopore': 1e-3,
        'pacbio': 1e-2,
    }

    def __init__(self, config: SpikeInConfig):
        self.config = config
        np.random.seed(int(config.tumor_fraction * 1e6 + config.n_samples))

        logger.info(
            f"SpikeInSimulator: VAF={config.tumor_fraction:.1e}, "
            f"GE={config.genome_equivalents:,}, depth={config.sequencing_depth:,}x, "
            f"n={config.n_samples}"
        )

    @property
    def background_error_rate(self) -> float:
        """Per-base background error rate for the configured platform."""
        return self.NOISE_RATES.get(self.config.noise_profile, 1e-4)

    def compute_molecule_counts(self) -> dict:
        """Compute the molecule and read budgets implied by this configuration.

        Two different quantities are reported and they must not be conflated:

        * ``expected_mutant_molecules_genome_wide`` -- ``GE * VAF``. This is how
          many mutant *molecules* exist in the plasma. At 0.001% VAF with
          50,000 GE this is 0.5: less than one molecule.
        * ``expected_mutant_reads`` -- ``depth * VAF``. How many reads carrying
          the mutation a sequencer would observe at a target site.
        """
        ge = self.config.genome_equivalents
        vaf = self.config.tumor_fraction
        depth = self.config.sequencing_depth

        molecules_in_region = ge * self.REGION_SIZE_FOR_VAF / self.BASES_PER_GENOME

        return {
            'genome_equivalents': ge,
            'molecules_in_target_region': molecules_in_region,
            'expected_mutant_molecules_genome_wide': ge * vaf,
            'expected_ref_molecules_genome_wide': ge * (1.0 - vaf),
            'sequencing_depth': depth,
            'expected_mutant_reads': depth * vaf,
            'expected_ref_reads': depth * (1.0 - vaf),
        }

    def generate_sample(self) -> dict:
        """Generate one simulated sample with spike-in mutations."""
        counts = self.compute_molecule_counts()
        depth = self.config.sequencing_depth
        bg_error = self.background_error_rate

        results = {
            'sample_id': None,
            'total_reads': depth,
            'genome_equivalents': self.config.genome_equivalents,
            'sequencing_depth': depth,
            'target_vaf': self.config.tumor_fraction,
            'mutations': {},
        }

        for gene, mutations in HOTSPOT_MUTATIONS.items():
            for mutation_name, _, allele_freq_in_cancer in mutations:
                this_mutation_fraction = self.config.tumor_fraction * allele_freq_in_cancer

                mutant_reads = int(np.random.poisson(depth * this_mutation_fraction))
                ref_reads = depth - mutant_reads

                false_alt_reads = int(np.random.binomial(ref_reads, bg_error))

                observed_alt = mutant_reads + false_alt_reads
                observed_ref = ref_reads - false_alt_reads
                total = observed_alt + observed_ref

                observed_vaf = float(observed_alt) / float(total) if total > 0 else 0.0

                results['mutations'][f'{gene}_{mutation_name}'] = {
                    'mutant_reads': mutant_reads,
                    'ref_reads': observed_ref,
                    'false_alt_reads': false_alt_reads,
                    'observed_alt_reads': observed_alt,
                    'observed_vaf': observed_vaf,
                    'target_vaf': float(this_mutation_fraction),
                    'allele_freq_in_cancer': float(allele_freq_in_cancer),
                    'called': bool(
                        mutant_reads >= self.MIN_SUPPORTING_READS
                        and observed_vaf > this_mutation_fraction * 0.5
                    ),
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
            'sequencing_depth': r['sequencing_depth'],
            'n_mutations': r['n_mutations'],
            'n_detected': r['n_detected'],
            'detection_rate': r['detection_rate'],
            'mean_observed_vaf': r['mean_observed_vaf'],
            'max_observed_vaf': r['max_observed_vaf'],
        } for r in all_results])

        return summary, all_results

    def power_analysis(self, allele_frequency: float = 1.0) -> dict:
        """Detection power at this configuration.

        ``allele_frequency`` is the fraction of tumour cells carrying the
        hotspot. The default of 1.0 is the most generous case (the mutation is
        clonal in the tumour); real hotspots are reported in
        ``HOTSPOT_MUTATIONS`` and are all below 0.25.

        Power is the probability that a single sample yields at least
        ``MIN_SUPPORTING_READS`` mutant reads -- the same threshold
        ``generate_sample`` applies when it decides whether to call a variant.
        Power therefore depends on sequencing depth, not on replicate count;
        the replicate count only reduces the variance of the estimate.
        """
        from scipy.stats import poisson

        vaf = self.config.tumor_fraction
        ge = self.config.genome_equivalents
        depth = self.config.sequencing_depth
        af = float(allele_frequency)

        effective_vaf = vaf * af
        expected_mutant_reads = depth * effective_vaf

        def detection_prob(lam: float) -> float:
            if lam <= 0:
                return 0.0
            return float(1 - sum(poisson.pmf(k, lam) for k in range(self.MIN_SUPPORTING_READS)))

        power = detection_prob(expected_mutant_reads)

        return {
            'target_vaf': vaf,
            'allele_frequency': af,
            'effective_vaf': effective_vaf,
            'expected_mutant_molecules_genome_wide': ge * vaf,
            'expected_mutant_reads': expected_mutant_reads,
            'sequencing_depth': depth,
            'single_sample_power': power,
            'min_supporting_reads': self.MIN_SUPPORTING_READS,
            'depth_for_80pct_power': self.depth_for_power(0.80, allele_frequency=af),
            'depth_for_95pct_power': self.depth_for_power(0.95, allele_frequency=af),
        }

    def depth_for_power(self, target_power: float, allele_frequency: float = 1.0) -> int:
        """Smallest sequencing depth reaching ``target_power`` detection power.

        Solves ``P(Poisson(depth * VAF * af) >= MIN_SUPPORTING_READS) >= target_power``
        for ``depth`` by bisection over the Poisson tail.
        """
        from scipy.stats import poisson

        vaf = self.config.tumor_fraction * float(allele_frequency)
        if vaf <= 0:
            return 0

        min_reads = self.MIN_SUPPORTING_READS

        def power_at(depth: int) -> float:
            lam = depth * vaf
            return float(1 - sum(poisson.pmf(k, lam) for k in range(min_reads)))

        lo, hi = 0, 1
        while power_at(hi) < target_power:
            hi *= 2
            if hi > 50_000_000:
                return -1
        while lo + 1 < hi:
            mid = (lo + hi) // 2
            if power_at(mid) >= target_power:
                hi = mid
            else:
                lo = mid
        return int(hi)

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


def run_lod_validation(
    vaf_range: list,
    n_samples: int = 20,
    out_csv: str = None,
    sequencing_depth: int = 100_000,
    genome_equivalents: int = 50_000,
    detection_threshold: float = 0.95,
) -> pd.DataFrame:
    """
    Run Limit of Detection validation across a range of VAF levels.

    Args:
        vaf_range: List of target VAFs to test (e.g. ``[1e-2, 1e-3, 1e-4]``).
        n_samples: Number of technical replicates per VAF level.
        out_csv: Optional output CSV path for committed ``results/`` artifacts.
        sequencing_depth: Reads covering a target site.
        genome_equivalents: Genome equivalents in the sample (50,000 == 10 mL plasma).
        detection_threshold: Fraction of hotspots that must be called for a
            VAF level to count as detected.

    Returns:
        DataFrame summarizing detection rate, statistical power, the depth
        needed for 80% power, and whether the LOD was reached, per VAF.
    """
    results = []

    for vaf in vaf_range:
        config = SpikeInConfig(
            tumor_fraction=vaf,
            genome_equivalents=genome_equivalents,
            n_samples=n_samples,
            sequencing_depth=sequencing_depth,
        )
        sim = SpikeInSimulator(config)
        summary, _ = sim.run_spike_in_experiment()
        power = sim.power_analysis()

        detection_rate = float(summary['detection_rate'].mean())
        power_estimate = float(power['single_sample_power'])

        results.append({
            'target_vaf': vaf,
            'sequencing_depth': sequencing_depth,
            'expected_mutant_molecules_genome_wide': power['expected_mutant_molecules_genome_wide'],
            'expected_mutant_reads': power['expected_mutant_reads'],
            'detection_rate_mean': detection_rate,
            'detection_rate_std': float(summary['detection_rate'].std()),
            'single_sample_power': power_estimate,
            'depth_for_80pct_power': power['depth_for_80pct_power'],
            'depth_for_95pct_power': power['depth_for_95pct_power'],
            'lod_achieved': detection_rate >= detection_threshold,
        })

        logger.info(
            f"VAF={vaf:.1e}: detection={detection_rate:.1%}, "
            f"power={power_estimate:.1%}, "
            f"depth_for_80pct={power['depth_for_80pct_power']:,}x"
        )

    df = pd.DataFrame(results)
    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        logger.info(f"Wrote LOD validation results to {out_csv}")
    return df


#: Default VAF ladder for the LOD sweep: 1% down to 1e-6 (1e-6 == 1 part
#: per million, below the 0.001% = 1e-5 design floor of CANCER-SEEK).
DEFAULT_VAF_LADDER = (1e-2, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6)


def _parse_vaf_range(raw: str) -> list:
    """Parse a comma-separated VAF ladder from the command line.

    Accepts ``1e-2,1e-3`` and ``1E-2, 1E-3``. Raises ``ValueError`` on an
    unparsable or non-positive entry so a typo fails loudly rather than
    silently running a different experiment.
    """
    values = []
    for token in raw.split(','):
        token = token.strip()
        if not token:
            continue
        value = float(token)
        if value <= 0.0:
            raise ValueError(f"VAF must be positive; got {token!r}.")
        if value > 1.0:
            raise ValueError(f"VAF must be <= 1.0; got {token!r}.")
        values.append(value)
    if not values:
        raise ValueError("No VAF values parsed from --vaf-range.")
    return values


def main(argv=None) -> int:
    """Console-script entry point for ``deepcatch-lod``.

    Returns a process exit code (0 on success) so the function is usable both
    as ``deepcatch-lod`` and from ``python -m src.ultrasensitive.spike_in_simulator``.
    """
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog="deepcatch-lod",
        description=(
            "Limit-of-Detection validation for the CANCER-SEEK ultra-sensitive "
            "pipeline. Simulates tumor DNA spiked into healthy cfDNA at known "
            "VAF levels and reports the detection rate and statistical power "
            "at each level. SYNTHETIC FIXTURE -- not a clinical measurement."
        ),
    )
    parser.add_argument(
        "--vaf-range",
        "--vaf_range",
        dest="vaf_range",
        default=",".join(f"{v:g}" for v in DEFAULT_VAF_LADDER),
        help=(
            "Comma-separated VAF levels to test "
            f"(default: {','.join(f'{v:g}' for v in DEFAULT_VAF_LADDER)})."
        ),
    )
    parser.add_argument(
        "--n-samples",
        "--n_samples",
        dest="n_samples",
        type=int,
        default=20,
        help="Technical replicates per VAF level (default: 20).",
    )
    parser.add_argument(
        "--sequencing-depth",
        "--sequencing_depth",
        dest="sequencing_depth",
        type=int,
        default=100_000,
        help=(
            "Reads covering a target site (default: 100000). Detection at a "
            "given VAF is governed by depth, not by replicate count."
        ),
    )
    parser.add_argument(
        "--genome-equivalents",
        "--genome_equivalents",
        dest="genome_equivalents",
        type=int,
        default=50_000,
        help="Genome equivalents in the sample (default: 50000, about 10 mL plasma).",
    )
    parser.add_argument(
        "--out-csv",
        default=None,
        help="Optional path for the results CSV (e.g. results/lod_validation.csv).",
    )
    args = parser.parse_args(argv)

    if args.n_samples < 1:
        parser.error("--n-samples must be >= 1")
    if args.sequencing_depth < 1:
        parser.error("--sequencing-depth must be >= 1")
    if args.genome_equivalents < 1:
        parser.error("--genome-equivalents must be >= 1")

    try:
        vaf_range = _parse_vaf_range(args.vaf_range)
    except ValueError as exc:
        parser.error(str(exc))

    print("=" * 70)
    print("Limit of Detection (LOD) Validation — CANCER-SEEK Spike-In Experiment")
    print("  (synthetic fixture — NOT a clinical measurement)")
    print(f"  depth={args.sequencing_depth:,}x  "
          f"genome_equivalents={args.genome_equivalents:,}  "
          f"replicates={args.n_samples}")
    print("=" * 70)

    results = run_lod_validation(
        vaf_range,
        n_samples=args.n_samples,
        out_csv=args.out_csv,
        sequencing_depth=args.sequencing_depth,
        genome_equivalents=args.genome_equivalents,
    )

    print("\n\nLOD VALIDATION SUMMARY")
    print("=" * 70)
    print(results.to_string(index=False))

    achieved = results[results['lod_achieved']]['target_vaf']
    if len(achieved):
        print(f"\nLOD reached at this depth (>=95% of hotspots called): "
              f"{float(achieved.max()):.1e}")
    else:
        print("\nNo VAF level reached >=95% hotspot detection at this depth.")
        cheapest = results.sort_values('target_vaf', ascending=False).iloc[0]
        print(f"  Lowest VAF tested: {cheapest['target_vaf']:.1e} "
              f"(detection {cheapest['detection_rate_mean']:.1%}); "
              f"depth for 80% power at that VAF: "
              f"{int(cheapest['depth_for_80pct_power']):,}x")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
