#!/usr/bin/env python3
"""
Regulatory Graph Builder — Methylation Network Construction
============================================================

Constructs heterogeneous graphs from methylation data for GNN-based
pre-cancer field defect detection.

Each node is a genomic regulatory region (CpG island, enhancer, promoter,
CTCF site, DHS). Edges come from two sources:

1. **Reference Hi-C** — 3D chromatin contact maps from public data
   (4DN/ENCODE) connecting physically interacting regions.
2. **Co-methylation** — regions whose methylation states are correlated
   across samples (from TCGA reference).
3. **Co-fragmentation** — regions whose cfDNA coverage patterns are
   correlated within the same sample.

Interactive demo
----------------

.. code-block:: python

    from src.methylation_gnn import RegulatoryGraphBuilder, GenomicRegion
    import numpy as np

    # Simulate methylation data for 2000 regions
    n_regions = 2000
    builder = RegulatoryGraphBuilder(n_nodes=n_regions, edge_k=15)

    methylation = {
        "beta_values": np.random.beta(2, 3, size=n_regions),
        "cpg_density": np.random.uniform(0, 30, size=n_regions),
        "gc_content": np.random.uniform(0.3, 0.7, size=n_regions),
    }

    coverage = np.random.poisson(50, size=n_regions).astype(float)

    graph = builder.build_graph(
        sample_name="demo_sample",
        methylation_data=methylation,
        cfDNA_coverage=coverage,
        label=0,
    )
    print(f"Nodes: {graph.x.shape[0]}, Edges: {graph.edge_index.shape[1]}")
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Union

import numpy as np

try:
    import torch
    from torch_geometric.data import Data

    _HAS_PYG = True
except ImportError:
    _HAS_PYG = False
    torch = None
    Data = None

    class Data:  # type: ignore
        """Placeholder for environments without PyG."""

        pass


logger = logging.getLogger(__name__)

# ── Constants ───────────────────────────────────────────────────

REGION_TYPES = {
    "cpg_island": 0,
    "enhancer": 1,
    "promoter": 2,
    "ctcf": 3,
    "dhs": 4,
    "other": 4,  # fallback
}

EDGE_TYPES = {
    "physical_interaction": 0,
    "co_methylation": 1,
    "co_fragmentation": 2,
    "genomic_proximity": 3,
    "regulatory_domain": 4,
}

N_REGION_TYPES = len(set(REGION_TYPES.values()))
N_EDGE_TYPES = len(EDGE_TYPES)


@dataclass
class GenomicRegion:
    """
    A node in the methylation graph representing a regulatory region.

    Attributes
    ----------
    chrom : str
        Chromosome label (e.g., 'chr1', 'chrX').
    start : int
        0-based start coordinate.
    end : int
        0-based end coordinate (exclusive).
    region_type : str
        One of 'cpg_island', 'enhancer', 'promoter', 'ctcf', 'dhs'.
    gene_name : str or None
        Associated gene symbol (if any).
    strand : str
        '+' or '-'.
    """

    chrom: str
    start: int
    end: int
    region_type: str = "other"
    gene_name: Optional[str] = None
    strand: str = "+"

    @property
    def length(self) -> int:
        """Region length in base pairs."""
        return self.end - self.start

    @property
    def midpoint(self) -> int:
        """Region midpoint coordinate."""
        return (self.start + self.end) // 2

    @property
    def type_id(self) -> int:
        """Numeric region type ID (0-4)."""
        return REGION_TYPES.get(self.region_type, 4)

    def distance_to(self, other: "GenomicRegion") -> int:
        """Genomic distance (bp) between midpoints (same chromosome) or inf."""
        if self.chrom != other.chrom:
            return int(1_000_000_000)  # effectively infinite
        return abs(self.midpoint - other.midpoint)


class RegulatoryGraphBuilder:
    """
    Construct heterogeneous regulatory graphs for GNN methylation analysis.

    Builds PyG ``Data`` objects where:
    - **Nodes** = genomic regulatory regions (CpG islands, enhancers, etc.)
    - **Edges** = Hi-C contacts, co-methylation, co-fragmentation
    - **Node features** = methylation β, CpG metrics, chromatin marks,
      fragmentomics signals, region type one-hot

    Gracefully handles missing data: if Hi-C or methylation correlation
    matrices are not available, falls back to co-fragmentation + proximity
    edges. This ensures the builder works with minimal input (cfDNA only).

    Parameters
    ----------
    n_nodes : int
        Maximum number of genomic regions (trim to top-N by relevance).
    edge_k : int
        Maximum edges per node (sparsification for scalability).
    hi_c_file : str or None
        Path to preprocessed Hi-C contact matrix (.npz format).
    reference_genome : str
        Genome build ('hg38' or 'hg19').
    seed : int
        Random seed for reproducibility.
    """

    def __init__(
        self,
        n_nodes: int = 50_000,
        edge_k: int = 20,
        hi_c_file: Optional[str] = None,
        reference_genome: str = "hg38",
        seed: int = 42,
    ):
        self.n_nodes = n_nodes
        self.edge_k = edge_k
        self.hi_c_file = hi_c_file
        self.reference_genome = reference_genome
        self.seed = seed

        self._regions: List[GenomicRegion] = []
        self._hi_c_matrix: Optional[np.ndarray] = None
        self._chromosome_map: Dict[str, List[int]] = {}  # chrom → region indices

        self._rng = np.random.RandomState(seed)

    # ── Region loading ──────────────────────────────────────────

    def load_reference_regions(self, bed_path: str) -> None:
        """
        Load genomic reference regions from a BED file.

        Expected format (tab-separated, no header):
            chrom   start   end   type   gene_name   strand

        Where ``type`` is one of: cpg_island, enhancer, promoter,
        ctcf, dhs.

        Parameters
        ----------
        bed_path : str
            Path to BED file with genomic region annotations.
        """
        self._regions = []
        self._chromosome_map = {}

        with open(bed_path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                region_type = parts[3].lower().replace(" ", "_")
                gene_name = parts[4] if len(parts) > 4 and parts[4] != "." else None
                strand = parts[5] if len(parts) > 5 else "+"

                region = GenomicRegion(
                    chrom=chrom,
                    start=start,
                    end=end,
                    region_type=region_type,
                    gene_name=gene_name,
                    strand=strand,
                )
                self._regions.append(region)

                idx = len(self._regions) - 1
                self._chromosome_map.setdefault(chrom, []).append(idx)

                if len(self._regions) >= self.n_nodes:
                    break

        logger.info(
            "Loaded %d genomic regions from %s", len(self._regions), bed_path
        )
        self.n_nodes = min(self.n_nodes, len(self._regions))

    def load_regions_from_list(
        self, regions: List[GenomicRegion]
    ) -> None:
        """Load regions directly from a list (no BED parsing needed)."""
        self._regions = regions[: self.n_nodes]
        self._chromosome_map = {}
        for i, r in enumerate(self._regions):
            self._chromosome_map.setdefault(r.chrom, []).append(i)
        logger.info("Loaded %d regions from list", len(self._regions))
        self.n_nodes = min(self.n_nodes, len(self._regions))

    def _ensure_regions(self) -> None:
        """Create synthetic regions if none were loaded."""
        if not self._regions:
            logger.warning(
                "No reference regions loaded; generating %d synthetic regions", self.n_nodes
            )
            self._regions = [
                GenomicRegion(
                    chrom=f"chr{(i % 22) + 1}",
                    start=i * 10_000,
                    end=i * 10_000 + 5_000,
                    region_type="other",
                )
                for i in range(self.n_nodes)
            ]
            self._chromosome_map = {}
            for i, r in enumerate(self._regions):
                self._chromosome_map.setdefault(r.chrom, []).append(i)

    # ── Hi-C loading ────────────────────────────────────────────

    def load_hi_c(self, hi_c_matrix: Optional[np.ndarray] = None) -> None:
        """
        Load preprocessed Hi-C contact matrix.

        Parameters
        ----------
        hi_c_matrix : (n_nodes, n_nodes) np.ndarray or None
            Square Hi-C contact matrix at region resolution.
            If None and hi_c_file was specified at init, tries to load from file.
        """
        if hi_c_matrix is not None:
            self._hi_c_matrix = hi_c_matrix
            logger.info("Loaded Hi-C matrix: %s", hi_c_matrix.shape)
            return

        if self.hi_c_file and self.hi_c_file.endswith(".npz"):
            try:
                data = np.load(self.hi_c_file)
                if "contacts" in data:
                    self._hi_c_matrix = data["contacts"]
                else:
                    self._hi_c_matrix = data[list(data.keys())[0]]
                logger.info(
                    "Loaded Hi-C matrix from %s: %s",
                    self.hi_c_file,
                    self._hi_c_matrix.shape,
                )
            except Exception as exc:
                logger.warning("Failed to load Hi-C file %s: %s", self.hi_c_file, exc)
                self._hi_c_matrix = None

    # ── Edge construction ───────────────────────────────────────

    def _build_hic_edges(self) -> Tuple[List[int], List[int], List[float]]:
        """Build edges from Hi-C contact matrix (reference 3D topology)."""
        src, dst, wgt = [], [], []
        if self._hi_c_matrix is None:
            return src, dst, wgt

        n = min(self.n_nodes, self._hi_c_matrix.shape[0])
        for i in range(n):
            row = self._hi_c_matrix[i, :n]
            # top-k Hi-C contacts per region
            if n <= self.edge_k:
                top_indices = np.argsort(row)[::-1][:self.edge_k]
            else:
                top_indices = np.argpartition(row, -self.edge_k)[-self.edge_k:]
            for j in top_indices:
                if j != i and row[j] > 0:
                    src.append(i)
                    dst.append(int(j))
                    wgt.append(float(np.log(row[j] + 1)))

        logger.info("Hi-C edges: %d", len(src))
        return src, dst, wgt

    def _build_comethylation_edges(
        self, methylation_corr: Optional[np.ndarray]
    ) -> Tuple[List[int], List[int], List[float]]:
        """Build edges from co-methylation correlation matrix."""
        src, dst, wgt = [], [], []
        if methylation_corr is None:
            return src, dst, wgt

        n = min(self.n_nodes, methylation_corr.shape[0])
        for i in range(n):
            row = np.abs(methylation_corr[i, :n])
            row[i] = -1  # exclude self
            if n <= self.edge_k:
                top_indices = np.argsort(row)[::-1][:self.edge_k]
            else:
                top_indices = np.argpartition(row, -self.edge_k)[-self.edge_k:]
            for j in top_indices:
                val = float(methylation_corr[i, j])
                if abs(val) > 0.1:  # minimum correlation threshold
                    src.append(i)
                    dst.append(int(j))
                    wgt.append(val)

        logger.info("Co-methylation edges: %d", len(src))
        return src, dst, wgt

    def _build_cofragmentation_edges(
        self, cfDNA_coverage: Optional[np.ndarray]
    ) -> Tuple[List[int], List[int], List[float]]:
        """
        Build edges from cfDNA coverage correlation.

        Regions with similar coverage patterns are likely co-fragmenting,
        which reflects 3D proximity and regulatory co-accessibility.
        """
        src, dst, wgt = [], [], []
        if cfDNA_coverage is None:
            return src, dst, wgt

        n = min(self.n_nodes, len(cfDNA_coverage))
        # Use chromosome grouping for proximity-based co-fragmentation analysis
        window = max(5, min(100, n // 10))  # adaptive window size

        for chrom, indices in self._chromosome_map.items():
            chrom_regions = [idx for idx in indices if idx < n]
            if len(chrom_regions) < 2:
                continue

            for pos, idx_i in enumerate(chrom_regions):
                start = max(0, pos - window)
                end = min(len(chrom_regions), pos + window + 1)
                neighbors = [
                    chrom_regions[k] for k in range(start, end) if k != pos
                ]
                # Limit to edge_k neighbors
                if len(neighbors) > self.edge_k:
                    neighbors = list(
                        self._rng.choice(neighbors, self.edge_k, replace=False)
                    )
                for idx_j in neighbors:
                    cov_i = cfDNA_coverage[idx_i]
                    cov_j = cfDNA_coverage[idx_j]
                    # Inverse coverage difference as similarity proxy
                    sim = 1.0 / (abs(cov_i - cov_j) + 1e-6)
                    w = float(np.log(sim + 1))
                    src.append(idx_i)
                    dst.append(idx_j)
                    wgt.append(w)

        logger.info("Co-fragmentation edges: %d", len(src))
        return src, dst, wgt

    def _build_proximity_edges(self) -> Tuple[List[int], List[int], List[float]]:
        """
        Build edges from linear genomic proximity.

        Regions within 100 kb on the same chromosome are connected,
        reflecting shared regulatory neighborhood effects.
        """
        src, dst, wgt = [], [], []
        n = self.n_nodes
        max_dist = 100_000  # 100 kb

        for chrom, indices in self._chromosome_map.items():
            chrom_regions = sorted(
                [idx for idx in indices if idx < n],
                key=lambda i: self._regions[i].start,
            )
            for pos_i in range(len(chrom_regions)):
                idx_i = chrom_regions[pos_i]
                reg_i = self._regions[idx_i]
                for pos_j in range(pos_i + 1, len(chrom_regions)):
                    idx_j = chrom_regions[pos_j]
                    reg_j = self._regions[idx_j]
                    dist = abs(reg_i.start - reg_j.start)
                    if dist > max_dist:
                        break  # sorted, so all subsequent are farther
                    src.append(idx_i)
                    dst.append(idx_j)
                    wgt.append(float(np.exp(-dist / (max_dist / 3))))

        logger.info("Proximity edges: %d", len(src))
        return src, dst, wgt

    def build_edges(
        self,
        cfDNA_coverage: Optional[np.ndarray] = None,
        methylation_corr: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Build multi-source edge set for the heterogeneous graph.

        Parameters
        ----------
        cfDNA_coverage : (n_nodes,) array or None
            Per-node cfDNA read coverage for co-fragmentation edges.
        methylation_corr : (n_nodes, n_nodes) array or None
            Precomputed methylation correlation matrix.

        Returns
        -------
        edge_index : (2, n_edges) int64 array
            Sparse adjacency list in COO format.
        edge_attr : (n_edges, 1) float32 array
            Edge weights.
        edge_type : (n_edges,) int64 array
            Edge type IDs (0=Hi-C, 1=co-methylation, 2=co-fragmentation,
            3=proximity, 4=regulatory_domain).
        """
        self._ensure_regions()
        n = self.n_nodes

        sources: Dict[int, Tuple[List, List, List]] = {}

        # 1. Hi-C physical interaction edges
        sources[0] = self._build_hic_edges()

        # 2. Co-methylation edges
        sources[1] = self._build_comethylation_edges(methylation_corr)

        # 3. Co-fragmentation edges (from cfDNA)
        sources[2] = self._build_cofragmentation_edges(cfDNA_coverage)

        # 4. Genomic proximity edges
        sources[3] = self._build_proximity_edges()

        # 5. Regulatory domain edges (TAD boundary sharing)
        #    Placeholder: share edges with proximity for now
        sources[4] = ([], [], [])

        # Merge all edge sources
        all_src, all_dst, all_wgt, all_type = [], [], [], []
        for type_id, (s, d, w) in sources.items():
            all_src.extend(s)
            all_dst.extend(d)
            all_wgt.extend(w)
            all_type.extend([type_id] * len(s))

        # If no edges at all, fall back to random edges (prevent empty graph)
        if not all_src:
            logger.warning("No edges produced; generating fallback random edges")
            for i in range(n):
                neighbors = list(self._rng.choice(n, min(self.edge_k, n - 1), replace=False))
                for j in neighbors:
                    if j != i:
                        all_src.append(i)
                        all_dst.append(int(j))
                        all_wgt.append(0.1)
                        all_type.append(3)  # proximity type

        edge_index = np.array([all_src, all_dst], dtype=np.int64)
        edge_attr = np.array(all_wgt, dtype=np.float32).reshape(-1, 1)
        edge_type = np.array(all_type, dtype=np.int64)

        logger.info(
            "Total edges: %d (avg %.1f per node)",
            edge_index.shape[1],
            edge_index.shape[1] / max(1, n),
        )
        return edge_index, edge_attr, edge_type

    # ── Node feature construction ───────────────────────────────

    def build_node_features(
        self,
        methylation_data: Dict[str, np.ndarray],
        chromatin_data: Optional[Dict[str, np.ndarray]] = None,
        fragmentomics_data: Optional[Dict[str, np.ndarray]] = None,
    ) -> np.ndarray:
        """
        Construct the node feature matrix from multi-modal data.

        Feature order (d=20):
        0:  mean_methylation       — β-value per region
        1:  methylation_variance   — epiallelic heterogeneity
        2:  cpg_density            — CpGs per 100bp
        3:  gc_content             — GC fraction
        4:  cpg_obs_exp            — Observed/expected CpG
        5:  coverage_depth         — Normalized cfDNA coverage
        6:  fragment_size_mean     — Mean fragment length
        7:  fragment_short_frac    — Fraction of short fragments
        8:  end_motif_diversity    — 4-mer motif diversity
        9:  methylation_entropy    — Shannon entropy of β-values
        10: dnase_signal           — Chromatin accessibility
        11: h3k4me3_signal         — Active promoter mark
        12: h3k27ac_signal         — Active enhancer mark
        13: h3k27me3_signal        — Polycomb repression
        14: h3k9me3_signal         — Constitutive heterochromatin
        15: region_type_cpg        — Is CpG island (one-hot)
        16: region_type_enhancer   — Is enhancer (one-hot)
        17: region_type_promoter   — Is promoter (one-hot)
        18: region_type_ctcf       — Is CTCF site (one-hot)
        19: region_type_dhs        — Is DHS site (one-hot)

        Missing data is filled with zeros and logged.

        Parameters
        ----------
        methylation_data : dict
            Keys: 'beta_values', optionally 'methylation_entropy',
            'methylation_variance', 'cpg_density', 'cpg_obs_exp', 'gc_content'.
        chromatin_data : dict or None
            Keys: 'dnase', 'h3k4me3', 'h3k27ac', 'h3k27me3', 'h3k9me3'.
        fragmentomics_data : dict or None
            Keys: 'coverage_depth', 'fragment_size_mean',
            'fragment_short_frac', 'end_motif_diversity'.

        Returns
        -------
        x : (n_nodes, 20) float32 array
            Node feature matrix.
        """
        self._ensure_regions()
        n = self.n_nodes

        # Helper: extract feature with fill value
        def _get(arr_dict, key, size):
            if arr_dict and key in arr_dict:
                val = arr_dict[key]
                if isinstance(val, np.ndarray):
                    padded = np.zeros(size, dtype=np.float32)
                    copy_len = min(len(val), size)
                    padded[:copy_len] = val[:copy_len]
                    return padded
                return np.full(size, float(val), dtype=np.float32)
            return np.zeros(size, dtype=np.float32)

        features: List[np.ndarray] = []

        # Methylation features (0-4, 9)
        features.append(_get(methylation_data, "beta_values", n))
        features.append(_get(methylation_data, "methylation_variance", n))
        features.append(_get(methylation_data, "cpg_density", n))
        features.append(_get(methylation_data, "gc_content", n))
        features.append(_get(methylation_data, "cpg_obs_exp", n))

        # Fragmentomics features (5-8)
        frag = fragmentomics_data or {}
        features.append(_get(frag, "coverage_depth", n))
        features.append(_get(frag, "fragment_size_mean", n))
        features.append(_get(frag, "fragment_short_frac", n))
        features.append(_get(frag, "end_motif_diversity", n))

        # Methylation entropy (9)
        features.append(_get(methylation_data, "methylation_entropy", n))

        # Chromatin features (10-14)
        chrom = chromatin_data or {}
        features.append(_get(chrom, "dnase", n))
        features.append(_get(chrom, "h3k4me3", n))
        features.append(_get(chrom, "h3k27ac", n))
        features.append(_get(chrom, "h3k27me3", n))
        features.append(_get(chrom, "h3k9me3", n))

        # Region type one-hot (15-19)
        type_ids = np.array(
            [REGION_TYPES.get(r.region_type, 4) for r in self._regions[:n]],
            dtype=np.int64,
        )
        onehot = np.eye(5, dtype=np.float32)[type_ids]
        for col in range(5):
            features.append(onehot[:, col])

        x = np.column_stack(features).astype(np.float32)

        # Log missing data stats
        n_missing = int(np.sum(np.all(x == 0, axis=1)))
        if n_missing > 0:
            logger.info(
                "Node features: %d/%d nodes have all-zero features (missing data)",
                n_missing, n,
            )

        logger.info("Node feature matrix: %s (%d features)", x.shape, x.shape[1])
        return x

    # ── Full graph construction ─────────────────────────────────

    def build_graph(
        self,
        sample_name: str,
        methylation_data: Dict[str, np.ndarray],
        cfDNA_coverage: Optional[np.ndarray] = None,
        methylation_corr: Optional[np.ndarray] = None,
        chromatin_data: Optional[Dict[str, np.ndarray]] = None,
        fragmentomics_data: Optional[Dict[str, np.ndarray]] = None,
        label: int = 0,
        cancer_type: str = "unknown",
        tumor_fraction: float = 0.0,
    ) -> "Data":
        """
        Build a complete PyG Data object for one cfDNA sample.

        This is the main entry point. It constructs nodes, edges, and
        packs everything into a ``torch_geometric.data.Data``.

        Parameters
        ----------
        sample_name : str
            Sample identifier (for logging and metadata).
        methylation_data : dict
            Methylation features per region. Must contain at minimum
            ``beta_values``.
        cfDNA_coverage : (n_nodes,) array or None
            Per-node cfDNA coverage for co-fragmentation edges.
        methylation_corr : (n_nodes, n_nodes) array or None
            Co-methylation correlation matrix.
        chromatin_data : dict or None
            Chromatin signal features (DNase, histone marks).
        fragmentomics_data : dict or None
            cfDNA fragment-level features (coverage, length, motifs).
        label : int
            0 = healthy, 1 = cancer.
        cancer_type : str
            TCGA abbreviation (e.g., 'LUAD', 'COAD').
        tumor_fraction : float
            Estimated tumor fraction in cfDNA (0 to 1).

        Returns
        -------
        torch_geometric.data.Data
            Graph with ``x``, ``edge_index``, ``edge_attr``, ``edge_type``,
            ``y``, and metadata attributes.

        Raises
        ------
        ImportError
            If torch_geometric is not installed.
        """
        if not _HAS_PYG:
            raise ImportError(
                "torch_geometric is required for graph construction. "
                "Install with: pip install torch_geometric"
            )

        self._ensure_regions()

        # Build features and edges
        x = self.build_node_features(
            methylation_data,
            chromatin_data=chromatin_data,
            fragmentomics_data=fragmentomics_data,
        )
        edge_index, edge_attr, edge_type = self.build_edges(
            cfDNA_coverage=cfDNA_coverage,
            methylation_corr=methylation_corr,
        )

        # Construct graph
        graph = Data(
            x=torch.from_numpy(x),
            edge_index=torch.from_numpy(edge_index),
            edge_attr=torch.from_numpy(edge_attr),
            edge_type=torch.from_numpy(edge_type),
            y=torch.tensor([label], dtype=torch.long),
            sample_name=sample_name,
            cancer_type=cancer_type,
            tumor_fraction=torch.tensor([tumor_fraction], dtype=torch.float32),
        )

        logger.info(
            "Built graph for %s: %d nodes, %d edges, %d features",
            sample_name,
            graph.x.shape[0],
            graph.edge_index.shape[1],
            graph.x.shape[1],
        )
        return graph

    def build_graph_batch(
        self,
        samples: List[Dict],
    ) -> List["Data"]:
        """
        Build graphs for multiple samples at once.

        Parameters
        ----------
        samples : list of dict
            Each dict must have keys matching ``build_graph`` parameters
            (at minimum: ``sample_name``, ``methylation_data``).

        Returns
        -------
        list of Data
            One PyG graph per sample.
        """
        graphs = []
        for sample in samples:
            g = self.build_graph(**sample)
            graphs.append(g)
        return graphs


# ── Aliases for backward compatibility ──────────────────────────

MethylationGraphBuilder = RegulatoryGraphBuilder
