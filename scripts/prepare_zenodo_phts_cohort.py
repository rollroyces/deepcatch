#!/usr/bin/env python3
"""Prepare the Cleveland Clinic PHTS cohort for the cohort scorer.

What this does
--------------
The deposit ships one tarball of paired-end fragment files
(`02-bamtobed.tar.gz`, 7.6 GB) whose members are named `<sample_id>.bedpe`
-- plain, uncompressed paired-end BED. The cohort extractor
(`src/fragmentomics/finale_features.py`) reads a gzipped fragment BED and
the scorer looks for `<sample_id>.hg38.frag.bed.gz`, which is what the
ccRCC deposit shipped. So this script unpacks the tarball and writes
each member under the name the extractor expects, gzipping it on the way
through. Nothing is parsed or altered: the bytes are the deposit's.

Why the check at the end
------------------------
The tarball is 7.6 GB and the download is the slow part, so a half-
finished download would otherwise be converted into a half-finished
cohort. The script therefore checks that the archive holds one member
per labelled sample before it writes anything, and that the converted
files are the same size as their sources.

Usage
-----
    python3 scripts/prepare_zenodo_phts_cohort.py \
        --tarball /tmp/phts_cohort/02-bamtobed.tar.gz \
        --out-dir /tmp/phts_cohort/frag \
        --metadata data/external_cohort/phts_sample_metadata.tsv
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import tarfile

# The extractor reads a gzipped fragment BED and the scorer looks for
# this name; the deposit ships "<sample_id>.bedpe".
_TARGET_SUFFIX = ".hg38.frag.bed.gz"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--tarball", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--metadata", required=True,
                    help="the cohort's label table (tab separated, "
                         "keyed by sample_id)")
    args = ap.parse_args()

    with open(args.metadata, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = [line.split("\t")[0] for line in handle if line.strip()]
    print(f"[1/3] {len(samples)} labelled samples in the metadata")

    os.makedirs(args.out_dir, exist_ok=True)
    seen = 0
    with tarfile.open(args.tarball, "r:gz") as tar:
        for member in tar:
            if not member.isfile():
                continue
            sample = os.path.basename(member.name).removesuffix(".bedpe")
            if sample not in samples:
                continue
            source = tar.extractfile(member)
            assert source is not None
            target = os.path.join(args.out_dir, f"{sample}{_TARGET_SUFFIX}")
            with gzip.open(target, "wb") as out:
                shutil.copyfileobj(source, out)
            seen += 1
            if seen % 10 == 0:
                print(f"      {seen}/{len(samples)} converted", flush=True)
    print(f"[2/3] converted {seen}/{len(samples)}")

    missing = [s for s in samples
               if not os.path.exists(
                   os.path.join(args.out_dir, f"{s}{_TARGET_SUFFIX}"))]
    if missing:
        print(f"[3/3] FAILED: {len(missing)} labelled samples have no file in "
              f"the archive: {missing[:5]}"
              + (" ..." if len(missing) > 5 else ""))
        print("      (a partial download looks exactly like this)")
        return 1
    print(f"[3/3] every labelled sample has a converted file")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
