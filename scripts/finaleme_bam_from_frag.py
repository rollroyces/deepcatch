#!/usr/bin/env python3
"""Build a FinaleMe-ready BAM from a FinaleDB fragment file.

Why
---
FinaleMe's per-CpG metrics class (``CpgMultiMetricsStatsV2``) reads a
BAM, and it needs each read to carry its **bases and qualities** -- it
inspects ``read.getReadBases()`` and ``read.getBaseQualities()`` per CpG.
A FinaleDB ``*.frag.tsv.bgz`` has no sequence, so a BAM built straight
from it fails with "Index 10 out of bounds for length 0".

This script bridges the two: it writes a SAM in which every fragment
becomes a perfectly-aligned read pair whose sequence is the reference
bases over the fragment span (the fragment records are exact matches
by construction), then sorts and indexes it with samtools.

Input schema (FinaleDB, 5 columns, bare chromosome names):
    chrom  start  end  mapq  strand

Output: coordinate-sorted, indexed BAM with "chr"-prefixed names, which
is what the hg19 reference files used by FinaleMe (CpG index, CG motif
bedgraph) expect.

Usage:
    python scripts/finaleme_bam_from_frag.py \
        --input  /path/to/EE12345.hg19.frag.tsv.bgz \
        --fa     /path/to/hg19.fa \
        --out-bam /path/to/EE12345.bam
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
from pathlib import Path

QUAL = "I"  # Phred 40, i.e. "this read matches the reference exactly"


class FastaIndex:
    """Minimal random-access reader over a faidx'd FASTA."""

    def __init__(self, fa_path: Path):
        self.fa = Path(fa_path)
        fai = Path(str(self.fa) + ".fai")
        if not fai.exists():
            subprocess.run(["samtools", "faidx", str(self.fa)], check=True)
        self.index: dict[str, tuple[int, int, int]] = {}
        with open(fai) as f:
            for line in f:
                name, length, offset, line_bases, line_width = line.split("\t")[:5]
                self.index[name] = (int(offset), int(length),
                                    int(line_bases))

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open [start, end) on `chrom`."""
        if chrom not in self.index:
            raise KeyError(f"chromosome {chrom} not in the reference index")
        offset, length, line_bases = self.index[chrom]
        if start >= end or end > length:
            return ""
        # Walk the FASTA's fixed-width lines to the first base of `start`.
        line = start // line_bases
        col = start % line_bases
        byte_pos = offset + line * (line_bases + 1) + col
        with open(self.fa, "rb") as f:
            f.seek(byte_pos)
            raw = f.read(end - start)
        return raw.decode("ascii", errors="replace").replace("\n", "")


def write_sam(frag_path: Path, fa: FastaIndex, sam_path: Path) -> int:
    n = 0
    with gzip.open(frag_path, "rt") as fi, open(sam_path, "w") as fo:
        fo.write("@HD\tVN:1.6\tSO:coordinate\n")
        for chrom, (offset, length, line_bases) in sorted(
                fa.index.items(), key=lambda kv: kv[1][1], reverse=True):
            fo.write(f"@SQ\tSN:{chrom}\tLN:{length}\n")
        for i, line in enumerate(fi):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0] if parts[0].startswith("chr") else "chr" + parts[0]
            start, end = int(parts[1]), int(parts[2])
            if end <= start:
                continue
            try:
                seq = fa.fetch(chrom, start, end)
            except KeyError:
                continue
            if not seq:
                continue
            L = len(seq)
            q = f"frag{i}"
            qual = QUAL * L
            fo.write(f"{q}\t99\t{chrom}\t{start + 1}\t60\t{L}M\t=\t{start + 1}"
                     f"\t{L}\t{seq}\t{qual}\n")
            fo.write(f"{q}\t147\t{chrom}\t{start + 1}\t60\t{L}M\t=\t{start + 1}"
                     f"\t-{L}\t{seq}\t{qual}\n")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="FinaleDB .frag.tsv.bgz")
    ap.add_argument("--fa", required=True, help="hg19 reference FASTA")
    ap.add_argument("--out-bam", required=True)
    ap.add_argument("--keep-sam", action="store_true")
    args = ap.parse_args()

    frag = Path(args.input)
    if not frag.is_file():
        print(f"input not found: {frag}", file=sys.stderr)
        return 1
    out_bam = Path(args.out_bam)
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    sam = out_bam.with_suffix(".sam")

    fa = FastaIndex(Path(args.fa))
    n = write_sam(frag, fa, sam)
    print(f"wrote {n} fragment pairs to {sam}", flush=True)

    subprocess.run(["samtools", "sort", "-O", "BAM", "-o", str(out_bam),
                    str(sam)], check=True)
    subprocess.run(["samtools", "index", str(out_bam)], check=True)
    if not args.keep_sam:
        sam.unlink()
    print(f"wrote {out_bam} (+ .bai)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
