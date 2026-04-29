"""
FragmentoSign: Fragmentomics Subsystem for DeepCatch.

Implements DELFI and MDS frameworks for cfDNA fragment analysis:
- LOESS GC-bias normalization
- 4-mer end motif extraction from BAM/FASTQ
- Gaussian Mixture Model for fragment length distributions
- Nucleosome positioning analysis

References:
- Cristiano et al. 2019, Nature 570:385-389
- Jiang et al. 2020, Nature Genetics 52:712-719
- Snyder et al. 2016, Cell 164:57-68
"""

from .normalization import (
    compute_gc_content,
    loess_normalize,
    DELFI_style_normalization,
    compute_MDS,
)
from .bam_motif_extractor import (
    extract_4mer_end_motifs,
    extract_end_motifs_from_fastq,
    compute_MDS as compute_MDS_from_counts,  # re-exported
)
from .fragment_gmm import (
    FragmentLengthGMM,
    compute_fragmentomics_features,
)

__all__ = [
    # Normalization
    "compute_gc_content",
    "loess_normalize",
    "DELFI_style_normalization",
    "compute_MDS",
    # BAM motif extraction
    "extract_4mer_end_motifs",
    "extract_end_motifs_from_fastq",
    "compute_MDS_from_counts",
    # Fragment GMM
    "FragmentLengthGMM",
    "compute_fragmentomics_features",
]
