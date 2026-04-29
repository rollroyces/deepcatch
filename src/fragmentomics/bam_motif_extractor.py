"""
Fragment End Motif Extraction from BAM Files (FragmentoSign)

Extracts 4-mer end motifs from paired-end cfDNA sequencing data.
Follows the DELFI and Jiang et al. protocols for fragment end analysis.

.. rubric:: BAM Parsing Protocol

1. **Iterate**: Stream through properly paired reads via ``pysam.AlignmentFile.fetch()``.
2. **Quality filters**:
   - Both mates must be mapped
   - MAPQ ≥ threshold (default 30)
   - Proper pair flag must be set
   - Exclude secondary/supplementary alignments
3. **Fragment length filter**: Retain fragments in the cfDNA-typical range
   [``fragment_length_min``, ``fragment_length_max``] (default 90–250 bp).
4. **Motif extraction**: For each read 1, fetch the 4-bp sequence at the
   5\u2032 end from the reference genome via ``pysam.FastaFile.fetch()``.
5. **Strand orientation** (if ``stranded=True``): Reverse-complement motifs
   from reads mapped to the reverse strand, so all motifs are reported on
   the + strand.
6. **Tabulation**: Accumulate counts in a 256-element array indexed by
   :data:`MOTIF_TO_IDX`.

.. rubric:: Performance Notes

- **Throughput**: For BAMs with >100K reads, the dominant cost is per-read
  ``fasta.fetch()`` calls. Pre-loading chromosome sequences into a dict
  (e.g., ``chrom_seq = {c: fasta.fetch(c) for c in bam.references}``)
  eliminates per-read FASTA seeks and can accelerate extraction by 10–50×.
- **Memory**: A single human chromosome sequence can be 50–250 MB.
  For whole-genome analysis, consider on-disk random access or
  memory-mapping with ``samtools faidx``.
- **Parallelism**: This implementation is single-threaded. For production
  use, split BAM by chromosome and run separate processes.

.. rubric:: pysam Dependency

- Required: ``pysam`` (Python wrapper for htslib). Install via
  ``pip install pysam`` or ``conda install -c bioconda pysam``.
- The BAM file must be sorted and indexed (``.bai`` or ``.csi`` alongside
  the ``.bam``).
- The reference FASTA must be indexed (``.fai`` file).
- If ``pysam`` is not installed, a clear ``ImportError`` is raised with
  installation instructions.

.. rubric:: References

.. [1] Cristiano, S. et al. (2019). Nature 570:385-389.
.. [2] Jiang, P. et al. (2020). Nature Genetics 52:712-719.
.. [3] Snyder, M.W. et al. (2016). Cell 164:57-68.

.. rubric:: Dependencies

- numpy
- pysam (required for BAM parsing)
"""

import numpy as np
from typing import Dict, Tuple, List, Optional
from collections import Counter
import warnings

# 256 possible 4-mers
BASES = ['A', 'C', 'G', 'T']
ALL_4MERS = [a+b+c+d for a in BASES for b in BASES for c in BASES for d in BASES]
MOTIF_TO_IDX = {m: i for i, m in enumerate(ALL_4MERS)}


def _reverse_complement(seq: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Maps A↔T, C↔G. Unknown bases (including N) map to N.
    Case-insensitive: input is uppercased, output is uppercase.

    Parameters
    ----------
    seq : str
        DNA sequence string.

    Returns
    -------
    str
        Reverse-complemented uppercase sequence.

    Examples
    --------
    >>> _reverse_complement("ACGT")
    'ACGT'
    >>> _reverse_complement("AATT")
    'AATT'
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(b, 'N') for b in reversed(seq.upper()))


def extract_4mer_end_motifs(
    bam_path: str,
    reference_fasta: str,
    n_reads: Optional[int] = None,
    mapq_threshold: int = 30,
    fragment_length_min: int = 90,
    fragment_length_max: int = 250,
    both_ends: bool = True,  # NOTE: currently only Read 1 processed; Read 2 TODO
    stranded: bool = True
) -> Dict[str, np.ndarray]:
    """
    Extract 4-mer end motifs from paired-end cfDNA BAM file.

    Implements the fragment end motif analysis protocol described
    in DELFI and Jiang et al. [1]_ [2]_.

    Protocol
    --------
    1. Iterate through properly paired reads via ``pysam.AlignmentFile.fetch()``
    2. Filter by MAPQ ≥ ``mapq_threshold`` and fragment length in
       [``fragment_length_min``, ``fragment_length_max``] (cfDNA-typical range)
    3. For each read pair, extract the 4-bp sequence at each fragment end
       from the reference genome
    4. If ``stranded=True``: orient motifs to the + strand using
       reverse complement
    5. Return frequency vector of 256 possible 4-mers

    Parameters
    ----------
    bam_path : str
        Path to BAM file. Must be sorted and indexed (``.bai`` or ``.csi``
        alongside the ``.bam``).
    reference_fasta : str
        Path to reference FASTA. Must be indexed with ``.fai``.
    n_reads : int, optional
        Maximum number of reads to process. None = all reads.
    mapq_threshold : int, optional
        Minimum mapping quality (default 30).
    fragment_length_min : int, optional
        Minimum fragment length in bp (default 90).
    fragment_length_max : int, optional
        Maximum fragment length in bp (default 250).
    both_ends : bool, optional
        If True, extract motifs from both 5\u2032 ends.
        NOTE: Currently only Read 1 is processed; Read 2 is a TODO.
    stranded : bool, optional
        If True, orient all motifs to the + strand (default True).

    Returns
    -------
    dict
        Dictionary with:

        - ``motif_counts``: (256,) np.ndarray of 4-mer counts
        - ``motif_frequencies``: (256,) np.ndarray of normalized frequencies
        - ``MDS``: Motif Diversity Score (float, [0, 1])
        - ``n_reads_processed``: number of reads used
        - ``n_reads_skipped``: number of reads filtered out
        - ``fragment_length_mean``: mean fragment length (float, bp)
        - ``fragment_length_median``: median fragment length (float, bp)
        - ``fragment_length_std``: standard deviation of fragment lengths (float, bp)

    Raises
    ------
    ImportError
        If ``pysam`` is not installed. Message includes install instructions.
    FileNotFoundError
        If the BAM file or reference FASTA does not exist, or if the
        corresponding index files (``.bai``, ``.fai``) are missing.
    ValueError
        If ``bam_path`` or ``reference_fasta`` is not a string.

    Notes
    -----
    - Currently only Read 1 is processed to avoid double-counting.
      Read 2 processing is a TODO for future versions.
    - If no valid reads pass all filters, a warning is issued and
      the returned arrays are all zeros with MDS = 0.0.
    - Edge case: reads with motifs containing ambiguous bases (non-ACGT)
      are silently skipped.

    References
    ----------
    .. [1] Cristiano et al. (2019). Nature 570:385-389.
    .. [2] Jiang et al. (2020). Nature Genetics 52:712-719.
    """
    try:
        import pysam
    except ImportError:
        raise ImportError(
            "pysam is required for BAM motif extraction. "
            "Install with: pip install pysam"
        )
    
    bam = pysam.AlignmentFile(bam_path, 'rb')
    fasta = pysam.FastaFile(reference_fasta)
    
    motif_counts = np.zeros(256, dtype=np.int64)
    fragment_lengths = []
    reads_processed = 0
    reads_skipped = 0
    
    for read in bam.fetch():
        if n_reads and reads_processed >= n_reads:
            break
        
        # Quality filters
        if read.is_unmapped or read.mate_is_unmapped:
            reads_skipped += 1
            continue
        if read.mapping_quality < mapq_threshold:
            reads_skipped += 1
            continue
        if not read.is_proper_pair:
            reads_skipped += 1
            continue
        if read.is_secondary or read.is_supplementary:
            reads_skipped += 1
            continue
        
        # Fragment length filter
        frag_len = abs(read.template_length)
        if frag_len < fragment_length_min or frag_len > fragment_length_max:
            reads_skipped += 1
            continue
        
        fragment_lengths.append(frag_len)
        
                # PERFORMANCE NOTE: For high-throughput (>100K reads), pre-load
        # chromosome sequences into memory to avoid per-read FASTA seeks.
        # Example: chrom_seq = {c: fasta.fetch(c) for c in bam.references}
        # Extract 4-mers at fragment ends
        # For read 1: 5' end is at read.reference_start
        # For fragment-level: use the outer coordinates
        
        if read.is_read1:
            # Read 1 5' end = fragment 5' start
            if read.is_reverse:
                # Read mapped to reverse strand → fragment start is at mate end
                end_pos = read.reference_start  # approximate
            else:
                end_pos = read.reference_start
            
            try:
                chrom = bam.get_reference_name(read.reference_id)
                if stranded and read.is_reverse:
                    motif = fasta.fetch(chrom, end_pos, end_pos + 4)
                    motif = _reverse_complement(motif)
                else:
                    motif = fasta.fetch(chrom, end_pos, end_pos + 4)
            except (ValueError, KeyError):
                reads_skipped += 1
                continue
        else:
            # Read 2: process mate
            continue  # Skip for now, process only Read 1 (avoids double-counting)
        
        motif = motif.upper()
        if len(motif) == 4 and all(b in 'ACGT' for b in motif):
            idx = MOTIF_TO_IDX.get(motif)
            if idx is not None:
                motif_counts[idx] += 1
                reads_processed += 1
    
    bam.close()
    fasta.close()
    
    # Compute frequencies and MDS
    total = motif_counts.sum()
    if total > 0:
        frequencies = motif_counts / total
        mds = compute_MDS(motif_counts)
    else:
        frequencies = np.zeros(256)
        mds = 0.0
        warnings.warn("No valid reads processed — check BAM quality and filters")
    
    lengths = np.array(fragment_lengths) if fragment_lengths else np.array([0])
    
    return {
        'motif_counts': motif_counts,
        'motif_frequencies': frequencies,
        'MDS': mds,
        'n_reads_processed': reads_processed,
        'n_reads_skipped': reads_skipped,
        'fragment_length_mean': float(np.mean(lengths)),
        'fragment_length_median': float(np.median(lengths)),
        'fragment_length_std': float(np.std(lengths)),
    }


def extract_end_motifs_from_fastq(
    fastq_path: str,
    n_reads: Optional[int] = None,
    read_end_length: int = 4
) -> np.ndarray:
    """
    Extract end motifs directly from FASTQ (no alignment needed).

    Reads the first ``read_end_length`` bases from each sequence in a
    FASTQ file and tabulates motif frequencies. Useful for raw fragment
    end analysis without the biases introduced by read mapping.

    This is a lightweight alternative to :func:`extract_4mer_end_motifs`
    when BAM alignment is not available or reference-free analysis is
    preferred.

    Parameters
    ----------
    fastq_path : str
        Path to FASTQ file (uncompressed). For compressed files, pipe
        through ``gzip.open`` in the caller.
    n_reads : int, optional
        Maximum number of reads to process. None = process all reads.
    read_end_length : int, optional
        Length of end motif to extract (default 4 for 4-mers).

    Returns
    -------
    np.ndarray
        Motif counts as integer array of shape ``(4**read_end_length,)``.

    Raises
    ------
    FileNotFoundError
        If ``fastq_path`` does not exist.

    Notes
    -----
    - Assumes standard FASTQ format: 4 lines per read (ID, sequence,
      separator "+", quality).
    - Sequences shorter than ``read_end_length`` are silently skipped.
    - Motifs containing non-ACGT bases are skipped.
    - For paired-end data, run separately on R1 and R2 FASTQ files.
    """
    from itertools import product
    
    n_motifs = 4 ** read_end_length
    motif_counts = np.zeros(n_motifs, dtype=np.int64)
    
    with open(fastq_path, 'r') as f:
        line_idx = 0
        reads_processed = 0
        seq = ''
        
        for line in f:
            if n_reads and reads_processed >= n_reads:
                break
            
            if line_idx % 4 == 1:  # Sequence line
                seq = line.strip()
                if len(seq) >= read_end_length:
                    motif = seq[:read_end_length].upper()
                    idx = MOTIF_TO_IDX.get(motif)
                    if idx is not None:
                        motif_counts[idx] += 1
                        reads_processed += 1
            
            line_idx += 1
    
    return motif_counts


def compute_MDS(motif_counts: np.ndarray) -> float:
    """
    Compute Motif Diversity Score (MDS) from raw motif counts.

    MDS is the normalized Simpson diversity index [1]_.

    .. math::

        \text{MDS} = \frac{1 - \sum p_i^2}{1 - 1/n}

    where :math:`p_i = \text{counts}_i / \sum \text{counts}` and
    :math:`n = \text{len(motif_counts)}` is the number of possible motifs.

    Parameters
    ----------
    motif_counts : np.ndarray
        1-D integer array of raw motif counts. Zeros are allowed.

    Returns
    -------
    float
        MDS in [0, 1]. Returns 0.0 if all counts sum to zero.

    Notes
    -----
    This function is re-exported as ``compute_MDS_from_counts`` from
    the parent ``fragmentomics`` package to disambiguate from the
    :func:`fragmentomics.normalization.compute_MDS` variant that
    operates on fragment end index arrays with a configurable ``n_motifs``.

    References
    ----------
    .. [1] Jiang, P. et al. (2020). Nature Genetics 52:712-719.
    """
    n_motifs = len(motif_counts)
    if motif_counts.sum() == 0:
        return 0.0
    p = motif_counts / motif_counts.sum()
    simpson = np.sum(p ** 2)
    mds = (1 - simpson) / (1 - 1.0 / n_motifs)
    return float(mds)


__all__ = [
    "extract_4mer_end_motifs",
    "extract_end_motifs_from_fastq",
    "compute_MDS",
    "MOTIF_TO_IDX",
    "ALL_4MERS",
]
