#!/usr/bin/env python3
"""Turn a FinaleMe per-CpG matrix into the input the HMM decoder expects.

The decoder (``org.cchmc.epifluidlab.finaleme.hmm.FinaleMe``) reads a
matrix whose last column is the methylation prior at each CpG, as a
percentage:

    chr  start  end  readName  FragLen  Frag_strand  methy_stat
    Norm_Frag_cov  baseQ  Offset_frag  Dist_frag_end  methyPrior

``CpgMultiMetricsStatsV2`` writes 11 columns and, with ``-kmerLen > 0``
(the default), appends a k-mer string as a 12th column -- so the prior
has to be supplied here. The prior comes from the WGBS buffy-coat
methylation bigwig that ships with FinaleMe
(``wgbs_buffyCoat_jensen2015GB.methy.hg19.bw``), the same file the
``-valueWigs methyPrior:0:<bw>`` option uses.

Missing priors are written as ``NaN`` (upper case: Java's
``Double.parseDouble`` accepts ``NaN`` but not ``nan``), and the decoder
skips those CpGs.

Usage:
    python scripts/finaleme_decoder_input.py \
        --matrix  MBC.cpg_matrix.txt.gz \
        --prior   /path/to/wgbs_buffyCoat_jensen2015GB.methy.hg19.bw \
        --out     MBC.decoder_input.txt.gz
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
from pathlib import Path


def build(matrix_path: Path, prior_path: Path, out_path: Path) -> dict:
    import pyBigWig

    bw = pyBigWig.open(str(prior_path))
    n_rows = 0
    n_with_prior = 0
    n_nan = 0
    with gzip.open(matrix_path, "rt") as fi, gzip.open(out_path, "wt",
                                                      compresslevel=6) as fo:
        header = fi.readline().rstrip("\n").split("\t")
        fo.write("\t".join(header[:11] + ["methyPrior"]) + "\n")
        for line in fi:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 11:
                continue
            prior = "NaN"
            try:
                vals = bw.values(parts[0], int(parts[1]), int(parts[2]))
                if vals and vals[0] is not None:
                    val = float(vals[0])
                    if math.isnan(val):
                        n_nan += 1
                    else:
                        prior = f"{val:.4f}"
                        n_with_prior += 1
            except Exception:
                n_nan += 1
            fo.write("\t".join(parts[:11] + [prior]) + "\n")
            n_rows += 1
    bw.close()
    return {"rows": n_rows, "with_prior": n_with_prior, "nan_prior": n_nan}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True,
                    help="per-CpG matrix from CpgMultiMetricsStatsV2 (.txt.gz)")
    ap.add_argument("--prior", required=True,
                    help="WGBS methylation prior bigwig (hg19)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    matrix = Path(args.matrix)
    if not matrix.is_file():
        print(f"matrix not found: {matrix}", file=sys.stderr)
        return 1
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    stats = build(matrix, Path(args.prior), out)
    print(f"wrote {out}: {stats['rows']} rows, {stats['with_prior']} with a "
          f"prior, {stats['nan_prior']} NaN")
    return 0


if __name__ == "__main__":
    sys.exit(main())
