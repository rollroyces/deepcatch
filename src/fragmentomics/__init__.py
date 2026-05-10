"""
FragmentoSign: Fragmentomics Subsystem for DeepCatch.

Provides a complete pipeline for cfDNA fragment analysis based on
fragmentation pattern signatures (FragmentoSigns). Implements the
DELFI [1]_ and MDS [2]_ frameworks.

.. rubric:: Public API

**Normalization** (``fragmentomics.normalization``)
    - :func:`compute_gc_content` — GC content of a DNA sequence
    - :func:`loess_normalize` — LOESS-based GC-bias correction
    - :func:`DELFI_style_normalization` — Full DELFI normalization pipeline
      (mappability filter → LOESS GC correction → median-centering)
    - :func:`compute_MDS` — Motif Diversity Score from motif index arrays

**BAM Motif Extraction** (``fragmentomics.bam_motif_extractor``)
    - :func:`extract_4mer_end_motifs` — Extract 4-mer end motifs from BAM
    - :func:`extract_end_motifs_from_fastq` — Extract end motifs from FASTQ
    - :func:`compute_MDS_from_counts` — MDS from raw count arrays
      (re-export of :func:`bam_motif_extractor.compute_MDS`)

**Fragment GMM** (``fragmentomics.fragment_gmm``)
    - :class:`FragmentLengthGMM` — Gaussian Mixture Model for fragment
      length distributions with 4 nucleosomal components
    - :func:`compute_fragmentomics_features` — Full feature extraction
      pipeline: basic statistics + DELFI + GMM features

.. rubric:: Submodules

- ``normalization`` — GC-bias correction and MDS computation
- ``bam_motif_extractor`` — 4-mer end motif extraction from BAM and FASTQ
- ``fragment_gmm`` — GMM-based fragment length analysis

.. rubric:: References

.. [1] Cristiano, S. et al. (2019). Nature 570:385-389. PMID: 31142840
.. [2] Jiang, P. et al. (2020). Cancer Discovery 10(5):664-673. PMID: 32111602
.. [3] Snyder, M.W. et al. (2016). Cell 164:57-68. PMID: 26771485
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
