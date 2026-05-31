#!/usr/bin/env python3
"""
Reference Data Utilities for GNN Methylation Network
=====================================================

Provides download URLs and preprocessing functions for the external
reference datasets required by the GNN branch:

- **CpG islands**: UCSC Genome Browser track
- **Promoters**: GENCODE / FANTOM5 annotations
- **Enhancers**: FANTOM5 enhancer atlas / GeneHancer
- **Hi-C**: 4D Nucleome / ENCODE chromatin contact maps

All download functions return the public URL and suggest an output path.
Actual downloading is left to the caller (to avoid network dependencies
in the core module).

Preprocessing functions convert raw methylation β-value matrices to the
node feature format expected by ``RegulatoryGraphBuilder``.

Datasets Summary
----------------

.. list-table::
    :header-rows: 1

    * - Dataset
      - Source
      - Content
      - Use

    * - UCSC CpG Islands
      - `genome.ucsc.edu <http://genome.ucsc.edu/>`_
      - ~28,000 CpG island regions
      - Node definitions

    * - GENCODE v44 (Promoters)
      - `gencodegenes.org <https://www.gencodegenes.org/>`_
      - TSS ± 2kb for all protein-coding genes
      - Node definitions

    * - FANTOM5 Enhancers
      - `fantom.gsc.riken.jp/5/ <https://fantom.gsc.riken.jp/5/>`_
      - Permissive enhancer atlas (65K+ regions)
      - Node definitions + chromatin state

    * - 4DN Hi-C (GM12878)
      - `4dnucleome.org <https://data.4dnucleome.org/>`_
      - 1kb-10kb resolution Hi-C
      - Reference graph edges

    * - ENCODE Chromatin
      - `encodeproject.org <https://www.encodeproject.org/>`_
      - DNase-seq, ChIP-seq (H3K4me3, etc.)
      - Chromatin node features

    * - TCGA Methylation
      - `portal.gdc.cancer.gov <https://portal.gdc.cancer.gov/>`_
      - Illumina 450K/850K β-value matrices
      - Training data + co-methylation edges

    * - Roadmap Epigenomics
      - `egg2.wustl.edu/roadmap <https://egg2.wustl.edu/roadmap/>`_
      - 127 reference epigenomes
      - Chromatin features

    * - GeneHancer
      - `genecards.org <https://www.genecards.org/>`_
      - Enhancer-promoter interactions
      - Regulatory edge types

Public Data URLs
----------------
All data is openly accessible. Some require an API key for programmatic
download (ENCODE, 4DN). TCGA requires GDC authentication.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# ── Public Dataset URLs ─────────────────────────────────────────

@dataclass
class DatasetEntry:
    """Metadata for a reference dataset."""
    name: str
    description: str
    url: str
    format: str
    size_estimate: str
    use_in_gnn: str
    requires_auth: bool = False


REFERENCE_DATASETS: Dict[str, DatasetEntry] = {
    "ucsc_cpg_islands": DatasetEntry(
        name="UCSC CpG Islands (hg38)",
        description="CpG island track from UCSC Genome Browser",
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz",
        format="Tab-separated (hg38 BED-style)",
        size_estimate="~5 MB compressed",
        use_in_gnn="Node definitions: CpG island regions",
    ),
    "ucsc_cpg_islands_hg19": DatasetEntry(
        name="UCSC CpG Islands (hg19)",
        description="CpG island track from UCSC Genome Browser",
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz",
        format="Tab-separated (hg19 BED-style)",
        size_estimate="~5 MB compressed",
        use_in_gnn="Node definitions: CpG island regions",
    ),
    "gencode_promoters": DatasetEntry(
        name="GENCODE v44 Promoters",
        description="TSS ± 2000bp for protein-coding genes",
        url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz",
        format="GTF (gene transfer format)",
        size_estimate="~50 MB compressed",
        use_in_gnn="Node definitions: promoter regions",
    ),
    "fantom5_enhancers": DatasetEntry(
        name="FANTOM5 Enhancer Atlas",
        description="Permissive enhancers from CAGE-seq data",
        url="https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2.bed.gz",
        format="BED",
        size_estimate="~3 MB compressed",
        use_in_gnn="Node definitions: enhancer regions",
    ),
    "encode_ctcf": DatasetEntry(
        name="ENCODE CTCF Binding Sites",
        description="CTCF ChIP-seq peaks across cell types",
        url="https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&target.label=CTCF&biosample_ontology.term_name=GM12878",
        format="BED (narrowPeak/broadPeak)",
        size_estimate="~50 MB across experiments",
        use_in_gnn="Node definitions: CTCF binding sites",
        requires_auth=False,
    ),
    "encode_dnase": DatasetEntry(
        name="ENCODE DNase-seq (GM12878)",
        description="DNase I hypersensitivity in GM12878 lymphoblastoid cells",
        url="https://www.encodeproject.org/experiments/ENCSR000EMT/",
        format="bigWig",
        size_estimate="~200 MB per replicate",
        use_in_gnn="Chromatin node features: dnase_signal",
    ),
    "encode_histone_marks": DatasetEntry(
        name="ENCODE Histone Marks (GM12878)",
        description="H3K4me3, H3K27ac, H3K27me3, H3K9me3 ChIP-seq",
        url="https://www.encodeproject.org/search/?type=Experiment&status=released&assay_title=Histone+ChIP-seq&biosample_ontology.term_name=GM12878",
        format="bigWig",
        size_estimate="~1 GB total",
        use_in_gnn="Chromatin node features: histone mark signals",
    ),
    "hic_rao2014": DatasetEntry(
        name="Rao et al. 2014 Hi-C (GM12878)",
        description="High-resolution Hi-C chromatin contact maps",
        url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525",
        format=".hic / .mcool / contact matrices",
        size_estimate="~10 GB per cell type (1kb resolution)",
        use_in_gnn="Reference graph edges: physical_interaction",
    ),
    "hic_4dn": DatasetEntry(
        name="4D Nucleome Hi-C",
        description="Multi-tissue Hi-C at 1kb resolution",
        url="https://data.4dnucleome.org/",
        format=".hic / .mcool",
        size_estimate="~50 GB total",
        use_in_gnn="Reference graph edges: tissue-specific interactions",
    ),
    "roadmap_epigenomics": DatasetEntry(
        name="Roadmap Epigenomics",
        description="127 reference epigenomes (methylation + histone marks)",
        url="https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/pval/",
        format="bigWig",
        size_estimate="~5 GB per mark",
        use_in_gnn="Chromatin features across tissues",
    ),
    "tcga_methylation": DatasetEntry(
        name="TCGA Methylation (Illumina 450K)",
        description="β-value matrices for 33 cancer types",
        url="https://portal.gdc.cancer.gov/",
        format="IDAT → processed β-value matrix",
        size_estimate="~2 TB total (processable in batches)",
        use_in_gnn="Training labels + methylation features",
        requires_auth=True,
    ),
}


class ReferenceDataCatalog:
    """
    Catalog of reference datasets with download instructions.

    Provides programmatic access to dataset metadata and download URLs.
    Does NOT perform downloads (to avoid network/disk dependencies in
    the core module). Use the standalone utility scripts for actual
    data acquisition.
    """

    @staticmethod
    def list_datasets() -> List[str]:
        """Return names of all registered reference datasets."""
        return list(REFERENCE_DATASETS.keys())

    @staticmethod
    def get_entry(name: str) -> Optional[DatasetEntry]:
        """Get metadata for a named dataset."""
        return REFERENCE_DATASETS.get(name)

    @staticmethod
    def print_catalog() -> str:
        """Human-readable catalog of all datasets."""
        lines = ["Reference Datasets for GNN Methylation Network", "=" * 60, ""]
        for key, entry in REFERENCE_DATASETS.items():
            lines.append(f"  {key}")
            lines.append(f"    Description: {entry.description}")
            lines.append(f"    URL: {entry.url}")
            lines.append(f"    Format: {entry.format}")
            lines.append(f"    Size: {entry.size_estimate}")
            lines.append(f"    Use: {entry.use_in_gnn}")
            if entry.requires_auth:
                lines.append(f"    ⚠️  Requires authentication")
            lines.append("")
        return "\n".join(lines)


# ── Convenience download functions ──────────────────────────────

def download_ucsc_cpg_islands(
    genome: str = "hg38", output_dir: str = "data/reference"
) -> str:
    """Return URL and suggested output path for UCSC CpG island track."""
    url = (
        f"https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/cpgIslandExt.txt.gz"
    )
    return url


def download_gencode_promoters(
    release: str = "44", output_dir: str = "data/reference"
) -> str:
    """Return URL for GENCODE annotation GTF."""
    url = (
        f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
        f"release_{release}/gencode.v{release}.annotation.gtf.gz"
    )
    return url


def download_fantom5_enhancers(output_dir: str = "data/reference") -> str:
    """Return URL for FANTOM5 permissive enhancer BED."""
    return (
        "https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/"
        "human_permissive_enhancers_phase_1_and_2.bed.gz"
    )


def download_encode_hic(cell_type: str = "GM12878") -> str:
    """Return URL hint for ENCODE Hi-C data (search page)."""
    return (
        f"https://www.encodeproject.org/search/"
        f"?type=Experiment&assay_title=Hi-C"
        f"&biosample_ontology.term_name={cell_type}"
    )


# ── Methylation data preprocessing ──────────────────────────────

def preprocess_methylation_betas(
    beta_matrix: np.ndarray,
    cpg_positions: Optional[np.ndarray] = None,
    regions: Optional[List[Tuple[str, int, int]]] = None,
    fill_missing: float = 0.5,
) -> Dict[str, np.ndarray]:
    """
    Preprocess a raw methylation β-value matrix to node features.

    Converts per-CpG β-values to per-region aggregate features:
    - mean β-value
    - methylation variance (epiallelic heterogeneity)
    - Shannon entropy of methylation pattern

    Parameters
    ----------
    beta_matrix : (n_cpgs, n_samples) array
        Raw methylation β-values in [0, 1] with possible NaN.
    cpg_positions : (n_cpgs, 2) array or None
        CpG genomic positions (chrom_idx, position). If None,
        regions must be provided for aggregation.
    regions : list of (chrom, start, end) or None
        Genomic region definitions to aggregate CpGs into.
        Required if cpg_positions is None.
    fill_missing : float
        Value to fill NaN/0 entries (default 0.5 = hemimethylated).

    Returns
    -------
    features : dict
        Keys:
        - 'beta_values': (n_regions, n_samples) mean β per region
        - 'methylation_variance': (n_regions, n_samples) variance
        - 'methylation_entropy': (n_regions, n_samples) Shannon entropy
        - 'cpg_density': (n_regions,) CpGs per 100bp
    """
    n_cpgs, n_samples = beta_matrix.shape

    # Fill missing values
    beta = beta_matrix.copy()
    beta[np.isnan(beta)] = fill_missing
    beta = np.clip(beta, 0.0, 1.0)

    if regions is not None and cpg_positions is not None:
        # Per-region aggregation
        n_regions = len(regions)
        agg_mean = np.zeros((n_regions, n_samples), dtype=np.float32)
        agg_var = np.zeros((n_regions, n_samples), dtype=np.float32)
        agg_entropy = np.zeros((n_regions, n_samples), dtype=np.float32)
        cpg_counts = np.zeros(n_regions, dtype=np.int32)

        for r_idx, (chrom, start, end) in enumerate(regions):
            # Find CpGs in this region
            mask = (
                (cpg_positions[:, 0] == chrom)
                & (cpg_positions[:, 1] >= start)
                & (cpg_positions[:, 1] < end)
            )
            region_cpgs = beta[mask, :]
            n_cpgs_region = region_cpgs.shape[0]

            if n_cpgs_region > 0:
                agg_mean[r_idx] = region_cpgs.mean(axis=0)
                agg_var[r_idx] = region_cpgs.var(axis=0)
                # Shannon entropy per sample
                if n_cpgs_region > 1:
                    # Discretize β into bins for entropy computation
                    for s in range(n_samples):
                        vals = region_cpgs[:, s]
                        hist, _ = np.histogram(vals, bins=10, range=(0, 1))
                        probs = hist / hist.sum()
                        valid = probs > 0
                        agg_entropy[r_idx, s] = -np.sum(
                            probs[valid] * np.log2(probs[valid])
                        ) / np.log2(10)  # normalized to [0, 1]
                else:
                    agg_entropy[r_idx, :] = 0.0

            cpg_counts[r_idx] = n_cpgs_region

        # CpG density (CpGs per 100bp)
        region_lengths = np.array([
            max(1, (end - start) / 100.0) for _, start, end in regions
        ], dtype=np.float32)
        cpg_density = cpg_counts.astype(np.float32) / region_lengths

        logger.info(
            "Preprocessed %d CpGs → %d regions (%d samples)",
            n_cpgs, n_regions, n_samples,
        )

        return {
            "beta_values": agg_mean.astype(np.float32),
            "methylation_variance": agg_var.astype(np.float32),
            "methylation_entropy": agg_entropy.astype(np.float32),
            "cpg_density": cpg_density.astype(np.float32),
        }

    # No regions: return per-CpG features (useful for small-scale tests)
    logger.info(
        "No region aggregation: returning per-CpG features (%d CpGs)", n_cpgs
    )
    return {
        "beta_values": beta.mean(axis=1).astype(np.float32),
        "methylation_variance": beta.var(axis=1).astype(np.float32),
        "methylation_entropy": np.zeros(n_cpgs, dtype=np.float32),
        "cpg_density": np.ones(n_cpgs, dtype=np.float32),  # 1 CpG per region
    }


def build_comethylation_matrix(
    beta_matrix: np.ndarray,
    min_correlation: float = 0.5,
) -> np.ndarray:
    """
    Build co-methylation correlation matrix from β-values.

    Computes pairwise Spearman correlation between CpG sites/regions
    across samples. Values with |ρ| < min_correlation are zeroed out
    (sparsification).

    Parameters
    ----------
    beta_matrix : (n_cpgs, n_samples) array
        Methylation β-values.
    min_correlation : float
        Minimum absolute correlation to keep an edge.

    Returns
    -------
    corr_matrix : (n_cpgs, n_cpgs) float32 array
        Sparse correlation matrix (zero for weak correlations).
    """
    n_cpgs = beta_matrix.shape[0]

    # For large matrices, use chunked computation
    if n_cpgs > 10_000:
        logger.warning(
            "Large co-methylation matrix (%d × %d): may be slow. "
            "Consider pre-filtering CpGs.",
            n_cpgs, n_cpgs,
        )

    beta = beta_matrix.copy()
    beta[np.isnan(beta)] = 0.5
    beta = np.clip(beta, 0.0, 1.0)

    # Rank-transform for Spearman correlation
    from scipy import stats

    ranked = np.zeros_like(beta)
    for i in range(n_cpgs):
        ranked[i] = stats.rankdata(beta[i])

    # Pearson on ranks = Spearman
    corr = np.corrcoef(ranked)
    corr = np.nan_to_num(corr, nan=0.0)

    # Sparsify
    corr[np.abs(corr) < min_correlation] = 0.0
    np.fill_diagonal(corr, 0.0)

    logger.info(
        "Co-methylation matrix: %d × %d, %.2f%% non-zero",
        n_cpgs, n_cpgs,
        100.0 * np.count_nonzero(corr) / (n_cpgs * n_cpgs),
    )

    return corr.astype(np.float32)


# ── BED file utilities ──────────────────────────────────────────

def parse_regions_from_bed(
    bed_path: str, max_regions: int = 100_000
) -> List[Tuple[str, int, int, str, Optional[str]]]:
    """
    Parse a BED file into a list of (chrom, start, end, type, gene) tuples.

    Expected format: chrom start end type [gene] [strand]
    Ignores BED headers (track, browser lines).

    Parameters
    ----------
    bed_path : str
        Path to BED file.
    max_regions : int
        Maximum regions to load.

    Returns
    -------
    regions : list of tuples
    """
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            region_type = parts[3].lower().replace(" ", "_")
            gene = parts[4] if len(parts) > 4 and parts[4] != "." else None
            regions.append((chrom, start, end, region_type, gene))
            if len(regions) >= max_regions:
                break
    return regions


def generate_cpg_regions(
    chrom_sizes: Dict[str, int],
    cpg_island_locations: Optional[List[Tuple[str, int, int]]] = None,
    n_regions: int = 50_000,
) -> List[Tuple[str, int, int, str, None]]:
    """
    Generate synthetic CpG-focused genomic regions for testing.

    Uses chromosome sizes to distribute regions proportionally across
    chromosomes. If CpG island locations are provided, prioritizes those.

    Parameters
    ----------
    chrom_sizes : dict
        Chromosome name → length (e.g., {'chr1': 248956422, ...}).
    cpg_island_locations : list or None
        Known CpG island positions.
    n_regions : int
        Number of regions to generate.

    Returns
    -------
    list of (chrom, start, end, type, None)
    """
    total_genome = sum(chrom_sizes.values())
    regions = []

    # Priority: known CpG islands
    if cpg_island_locations:
        for chrom, start, end in cpg_island_locations[:n_regions]:
            regions.append((chrom, start, end, "cpg_island", None))
        n_remaining = n_regions - len(regions)
    else:
        n_remaining = n_regions

    # Fill remaining with evenly-spaced windows per chromosome
    if n_remaining > 0:
        for chrom, length in sorted(chrom_sizes.items()):
            frac = length / total_genome
            n_chrom_regions = max(1, int(frac * n_remaining))
            step = max(1000, length // (n_chrom_regions + 1))
            for i in range(n_chrom_regions):
                start = i * step
                end = min(start + step // 2, length)
                region_type = "cpg_island" if i % 3 == 0 else "other"
                regions.append((chrom, start, end, region_type, None))
                if len(regions) >= n_regions:
                    break
            if len(regions) >= n_regions:
                break

    return regions[:n_regions]


# ── Synthetic data generation for testing ───────────────────────

def generate_synthetic_methylation_data(
    n_regions: int = 5_000,
    n_samples: int = 100,
    n_cancer: int = 50,
    seed: int = 42,
) -> Dict[str, np.ndarray]:
    """
    Generate synthetic methylation data for development/testing.

    Creates realistic-looking β-value matrices with:
    - Healthy samples: concentrated around 0.3 (unmethylated) or 0.8 (methylated)
    - Cancer samples: global hypomethylation (shift toward 0.2) with
      some CpG islands showing hypermethylation (shift toward 0.9)

    Parameters
    ----------
    n_regions : int
        Number of simulated genomic regions.
    n_samples : int
        Total number of samples.
    n_cancer : int
        Number of cancer samples (rest are healthy).
    seed : int
        Random seed.

    Returns
    -------
    data : dict
        Keys:
        - 'beta_values': (n_regions, n_samples) β-value matrix
        - 'labels': (n_samples,) binary labels (0=healthy, 1=cancer)
        - 'cpg_density': (n_regions,) CpG density
        - 'gc_content': (n_regions,) GC content
    """
    rng = np.random.RandomState(seed)
    n_healthy = n_samples - n_cancer

    # Generate base methylation states (bimodal: unmethylated vs methylated)
    base_state = rng.choice([0.2, 0.8], size=n_regions, p=[0.3, 0.7])

    # Add biological noise
    healthy_beta = np.zeros((n_regions, n_healthy), dtype=np.float32)
    for i in range(n_healthy):
        noise = rng.normal(0, 0.05, size=n_regions)
        healthy_beta[:, i] = np.clip(base_state + noise, 0.0, 1.0)

    # Cancer: global hypomethylation + focal hypermethylation at some regions
    hypermethylation_mask = rng.rand(n_regions) < 0.3
    cancer_beta = np.zeros((n_regions, n_cancer), dtype=np.float32)
    for i in range(n_cancer):
        noise = rng.normal(0, 0.08, size=n_regions)
        cancer_state = base_state.copy()
        cancer_state[~hypermethylation_mask] -= 0.15  # global hypomethylation
        cancer_state[hypermethylation_mask] += 0.1  # focal hypermethylation
        cancer_beta[:, i] = np.clip(cancer_state + noise, 0.0, 1.0)

    # Combine
    beta_values = np.hstack([healthy_beta, cancer_beta]).astype(np.float32)
    labels = np.array([0] * n_healthy + [1] * n_cancer, dtype=np.int64)

    # Region-level features
    cpg_density = rng.uniform(0, 30, size=n_regions).astype(np.float32)
    gc_content = rng.uniform(0.3, 0.7, size=n_regions).astype(np.float32)

    logger.info(
        "Generated synthetic methylation: %d regions, %d samples (%d cancer)",
        n_regions, n_samples, n_cancer,
    )

    return {
        "beta_values": beta_values,
        "labels": labels,
        "cpg_density": cpg_density,
        "gc_content": gc_content,
    }
