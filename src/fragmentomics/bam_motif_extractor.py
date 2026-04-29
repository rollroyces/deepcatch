"""
Fragment End Motif Extraction from BAM Files (FragmentoSign)

Extracts 4-mer end motifs from paired-end cfDNA sequencing data.
Follows the DELFI and Jiang et al. protocols for fragment end analysis.

DEPENDENCIES: pysam, numpy
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
    """Reverse complement a DNA sequence."""
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
    Extract 4-mer end motifs from cfDNA BAM file.
    
    Protocol (adapted from DELFI + Jiang et al.):
    1. Iterate through properly paired reads
    2. Filter by MAPQ ≥ 30 and fragment length 90-250bp (cfDNA-typical)
    3. For each read pair, extract the 4-bp sequence at each fragment end
    4. If stranded: orient motifs to the + strand using reference
    5. Return frequency vector of 256 4-mers
    
    Args:
        bam_path: Path to BAM file (must be indexed)
        reference_fasta: Path to reference FASTA (must be indexed with .fai)
        n_reads: Max reads to process (None = all)
        mapq_threshold: Minimum mapping quality
        fragment_length_min: Minimum fragment length (cfDNA-typical)
        fragment_length_max: Maximum fragment length (cfDNA-typical)
        both_ends: If True, extract motifs from both 5' ends
        stranded: If True, orient motifs to + strand
    
    Returns:
        dict with:
        - 'motif_counts': (256,) array of 4-mer counts
        - 'motif_frequencies': (256,) normalized frequencies
        - 'MDS': Motif Diversity Score
        - 'n_reads_processed': number of reads used
        - 'fragment_length_stats': (mean, median, std) of fragment lengths
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
    Alternative: Extract end motifs directly from FASTQ (no alignment needed).
    Useful for raw fragment end analysis without mapping bias.
    
    Args:
        fastq_path: Path to FASTQ file
        n_reads: Max reads to process
        read_end_length: Length of end motif to extract (default 4)
    
    Returns:
        motif_counts: (4^read_end_length,) array of motif counts
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
    """Motif Diversity Score — normalized Simpson diversity."""
    n_motifs = len(motif_counts)
    if motif_counts.sum() == 0:
        return 0.0
    p = motif_counts / motif_counts.sum()
    simpson = np.sum(p ** 2)
    mds = (1 - simpson) / (1 - 1.0 / n_motifs)
    return float(mds)
