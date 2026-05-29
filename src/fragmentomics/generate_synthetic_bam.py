"""
Generate synthetic paired-end BAM file with cfDNA-like fragment properties.

Used for testing bam_motif_extractor.py when real cfDNA BAM data
is not available.

Creates:
- Reads with cfDNA-like fragment length distribution (mode ~167bp)
- Properly paired reads with realistic CIGAR strings
- Sorted and indexed BAM ready for motif extraction
"""

import numpy as np
import pysam
import os
import sys


def create_fasta_index(fasta_path: str) -> str:
    """
    Create a .fai index for a FASTA file.
    Reads the FASTA to determine line lengths and byte offsets.
    """
    fai_path = fasta_path + '.fai'
    if os.path.exists(fai_path):
        print(f"FAI already exists: {fai_path}")
        return fai_path
    
    with open(fasta_path, 'r') as f:
        offsets = {}
        lengths = {}
        current_chrom = None
        current_len = 0
        current_offset = None
        line_bases = None
        line_width = None
        
        pos = 0
        for line in f:
            line_stripped = line.rstrip('\n')
            if line.startswith('>'):
                if current_chrom is not None:
                    lengths[current_chrom] = current_len
                    offsets[current_chrom] = {
                        'offset': current_offset,
                        'line_bases': line_bases,
                        'line_width': line_width
                    }
                current_chrom = line_stripped[1:].split()[0]
                current_len = 0
                current_offset = None
                line_bases = None
                line_width = None
            else:
                seq = line_stripped.upper()
                if current_offset is None:
                    current_offset = pos
                    line_bases = len(seq)
                    line_width = len(line.rstrip('\r\n'))
                current_len += len(seq)
            
            pos += len(line)
        
        if current_chrom is not None:
            lengths[current_chrom] = current_len
            offsets[current_chrom] = {
                'offset': current_offset,
                'line_bases': line_bases,
                'line_width': line_width
            }
    
    with open(fai_path, 'w') as f:
        for chrom in lengths:
            o = offsets[chrom]
            f.write(f"{chrom}\t{lengths[chrom]}\t{o['offset']}\t{o['line_bases']}\t{o['line_width']}\n")
    
    print(f"Created FAI index: {fai_path}")
    return fai_path


def generate_cfDNA_fragments(
    reference_fasta: str,
    chrom: str,
    n_fragments: int,
    read_length: int = 150,
    frag_mean: int = 167,
    frag_std: int = 25,
    seed: int = 42
) -> list:
    """
    Generate synthetic cfDNA-like paired-end reads.
    
    Parameters
    ----------
    reference_fasta : str
        Path to indexed reference FASTA
    chrom : str
        Chromosome to generate reads from
    n_fragments : int
        Number of fragment pairs to generate
    read_length : int
        Length of each read (default 150bp, typical for Illumina)
    frag_mean : int
        Mean fragment length (cfDNA mode ~167bp)
    frag_std : int
        Standard deviation of fragment length
    seed : int
        Random seed for reproducibility
        
    Returns
    -------
    list of dict with read1_seq, read1_pos, read2_seq, read2_pos, is_reverse
    """
    rng = np.random.default_rng(seed)
    
    fasta = pysam.FastaFile(reference_fasta)
    chrom_seq = fasta.fetch(chrom)
    chrom_len = len(chrom_seq)
    fasta.close()
    
    fragments = []
    
    for _ in range(n_fragments):
        # Generate fragment length from (truncated) normal distribution
        frag_len = int(rng.normal(frag_mean, frag_std))
        frag_len = max(read_length + 10, min(1000, frag_len))
        
        # Pick a random start position
        max_start = chrom_len - frag_len - 1
        if max_start <= 0:
            continue
        
        start = int(rng.integers(0, max_start))
        end = start + frag_len
        
        # Randomly designate strand (50% forward, 50% reverse)
        is_reverse = rng.random() > 0.5
        
        # Read 1: 5' end of the fragment
        r1_seq = chrom_seq[start:start + read_length]
        if len(r1_seq) < read_length:
            continue
            
        # Read 2: 5' end of the reverse complement at the other end
        r2_start = end - read_length
        r2_seq = chrom_seq[r2_start:end]
        if len(r2_seq) < read_length:
            continue
        
        # Reverse complement for Read 2 (it reads from the opposite strand)
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        r2_seq_rc = ''.join(complement.get(b, 'N') for b in reversed(r2_seq.upper()))
        
        fragments.append({
            'read1_seq': r1_seq,
            'read1_pos': start,
            'r1_is_reverse': is_reverse,
            'read2_seq': r2_seq_rc,
            'read2_pos': r2_start,
            'r2_is_reverse': not is_reverse,
            'frag_len': frag_len,
        })
    
    return fragments


def write_bam(
    fragments: list,
    reference_fasta: str,
    chrom: str,
    output_bam: str,
    read_group: str = "simulated_cfdna"
):
    """
    Write synthetic fragments to a sorted, indexed BAM file.
    """
    # Open reference to get chrom length
    fasta = pysam.FastaFile(reference_fasta)
    chrom_len = len(fasta.fetch(chrom))
    fasta.close()
    
    # Create header
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': chrom, 'LN': chrom_len}],
        'RG': [{'ID': read_group, 'SM': 'synthetic_cfdna', 'LB': 'lib1', 'PL': 'ILLUMINA'}],
    }
    
    # Collect all reads
    all_reads = []
    for i, f in enumerate(fragments):
        qname = f"frag_{i}"
        flag_r1 = 99 if not f['r1_is_reverse'] else 83  # 99: first in pair, mapped, proper pair; 83: first, mapped, proper pair, reverse
        flag_r2 = 147 if not f['r2_is_reverse'] else 163  # 147: second in pair, mapped, proper pair, reverse
        
        # Read 1
        a1 = pysam.AlignedSegment()
        a1.query_name = qname
        a1.query_sequence = f['read1_seq']
        a1.flag = flag_r1
        a1.reference_id = 0
        a1.reference_start = f['read1_pos']
        a1.mapping_quality = 60
        a1.cigarstring = f"{len(f['read1_seq'])}M"
        a1.next_reference_id = 0
        a1.next_reference_start = f['read2_pos']
        a1.template_length = f['frag_len'] if not f['r1_is_reverse'] else -f['frag_len']
        a1.query_qualities = pysam.qualities_to_qualitystring([40] * len(f['read1_seq']))
        a1.set_tag('RG', read_group)
        
        # Read 2
        a2 = pysam.AlignedSegment()
        a2.query_name = qname
        a2.query_sequence = f['read2_seq']
        a2.flag = flag_r2
        a2.reference_id = 0
        a2.reference_start = f['read2_pos']
        a2.mapping_quality = 60
        a2.cigarstring = f"{len(f['read2_seq'])}M"
        a2.next_reference_id = 0
        a2.next_reference_start = f['read1_pos']
        a2.template_length = -f['frag_len'] if not f['r2_is_reverse'] else f['frag_len']
        a2.query_qualities = pysam.qualities_to_qualitystring([40] * len(f['read2_seq']))
        a2.set_tag('RG', read_group)
        
        all_reads.append(a1)
        all_reads.append(a2)
    
    # Sort by position
    all_reads.sort(key=lambda r: (r.reference_id, r.reference_start))
    
    # Write BAM
    with pysam.AlignmentFile(output_bam, 'wb', header=header) as out:
        for r in all_reads:
            out.write(r)
    
    # Index
    pysam.index(output_bam)
    
    print(f"Wrote {len(all_reads)} reads ({len(fragments)} fragments) to {output_bam}")
    print(f"Created index: {output_bam}.bai")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate synthetic cfDNA BAM for testing")
    parser.add_argument('--fasta', required=True, help='Reference FASTA (must have .fai)')
    parser.add_argument('--chrom', default='chr22', help='Chromosome to use')
    parser.add_argument('--n-fragments', type=int, default=100000, help='Number of fragments')
    parser.add_argument('--read-length', type=int, default=150, help='Read length')
    parser.add_argument('--frag-mean', type=int, default=167, help='Mean fragment length')
    parser.add_argument('--frag-std', type=int, default=25, help='Fragment length std')
    parser.add_argument('--output', required=True, help='Output BAM path')
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    # Ensure FAI exists
    fai_path = args.fasta + '.fai'
    if not os.path.exists(fai_path):
        print(f"Creating FAI index for {args.fasta}...")
        create_fasta_index(args.fasta)
    
    print(f"Generating {args.n_fragments} cfDNA-like fragments from {args.chrom}...")
    fragments = generate_cfDNA_fragments(
        args.fasta,
        args.chrom,
        args.n_fragments,
        read_length=args.read_length,
        frag_mean=args.frag_mean,
        frag_std=args.frag_std,
        seed=args.seed
    )
    
    print(f"Writing BAM to {args.output}...")
    write_bam(fragments, args.fasta, args.chrom, args.output)
    print("Done!")


if __name__ == '__main__':
    main()
