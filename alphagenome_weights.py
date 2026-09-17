"""
AlphaGenome Atlas AVI weighting for DeepCatch v2.2 panel-LLR.

Provides deterministic per-variant AVI (AlphaGenome Variant Impact) weights
that are used as fixed multipliers in the panel log-likelihood aggregation.

TERMS-OF-SERVICE NOTES (CRITICAL):
--------------------------------------
Per the AlphaGenome Terms of Use, outputs "should not be used for the
training of other machine learning models." We AVOID that violation in two
ways:

  1. AVI scores are consumed ONLY as fixed scalar weights, applied by a
     deterministic hand-coded aggregator (a weighted sum / max of weighted
     per-locus scores). No learned weights, no gradient updates, no model
     parameters are fit using AVI scores as inputs or labels.

  2. The AVI-weighted aggregator is structurally analogous to weighting by
     any published conservation or deleteriousness score (PhyloP, CADD,
     AlphaMissense, SIFT, PolyPhen). The use of such scores as weights in
     deterministic MRD-style aggregation is established practice in the
     ctDNA literature (e.g. CAPP-Seq panel prioritisation).

AVI is also a NON-COMMERCIAL artifact per the Terms.

Data sources, in priority order:
  - Live Atlas API      (preferred when ALPHAGENOME_API_KEY is set)
  - Local Tabix bulk    (when a pre-downloaded avi_snvs_*.tsv.gz cache exists)
  - Deterministic proxy (fallback; clearly flagged as PROXY_USED in results)

The deterministic proxy is NOT a substitute for real AVI scores. It is
derived only from TCGA-VCF-derived annotations we already have on hand
(variant_class, position in the panel). It exists so that the weighted
aggregator can be exercised end-to-end and unit-tested without external
dependencies. Any run that uses the proxy is marked in its output JSON.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

LOG = logging.getLogger("deepcatch.alphagenome_weights")

# ---------------------------------------------------------------------------
# Variant key + weight container
# ---------------------------------------------------------------------------

def variant_key(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Stable GRCh38-style variant key used for AVI lookup."""
    # Strip leading 'chr' for canonical key, but keep chromosome in the key
    chrom_canon = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
    return f"chr{chrom_canon}:{pos}:{ref}>{alt}"


@dataclass
class AVIWeight:
    """Container for one variant's AVI weight and provenance."""
    variant_key: str
    avi_score: float           # AVI on its native scale (Phred-scaled per Atlas paper)
    avi_norm: float            # AVI normalised to [0, 1] for weighting
    source: str                # 'atlas_api' | 'tabix_local' | 'proxy' | 'missing'
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Deterministic proxy (fallback only)
# ---------------------------------------------------------------------------

# Order-of-magnitude mapping mirroring the biological ranking that AVI is
# *trained to* recover (per the Atlas preprint: discriminates pathogenic
# non-coding from benign, integrating AlphaMissense + VEP termination +
# conservation). Higher = more deleterious, as AVI is Phred-scaled.
#
# Values below are NOT published AVI numbers; they are a stable per-variant
# anchor used purely so the weighted-aggregator code path runs.
_PROXY_BY_CLASS: Mapping[str, float] = {
    "Nonsense_Mutation":      30.0,
    "Frame_Shift_Del":        29.5,
    "Frame_Shift_Ins":        29.0,
    "Splice_Site":            27.0,
    "Nonstop_Mutation":       26.0,
    "Translation_Start_Site": 25.0,
    "Missense_Mutation":      18.0,
    "In_Frame_Del":           16.0,
    "In_Frame_Ins":           15.5,
    "Splice_Region":          14.0,
    "5'UTR":                  10.0,
    "3'UTR":                   8.0,
    "Silent":                  3.0,
    "Intron":                  2.0,
    "IGR":                     1.0,
    "RNA":                     1.0,
}
_PROXY_DEFAULT = 5.0


def proxy_avi_for_variant(variant_class: str) -> float:
    """Deterministic per-variant AVI proxy from a MAF Variant_Classification.

    Used ONLY when neither the Atlas API nor a Tabix cache is available.
    Returned on AVI's Phred-like scale. Stable for any given variant_class.
    """
    return _PROXY_BY_CLASS.get(variant_class, _PROXY_DEFAULT)


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def normalize_avi(scores: Sequence[float]) -> np.ndarray:
    """Map AVI scores to weights in (0, 1] with a small floor.

    Floor of 0.05 prevents a mutation with AVI=0 from being entirely
    silenced when its signal is non-zero. Uses min-max scaling on the
    observed cohort; this matches published CADD-style normalisation.
    """
    arr = np.asarray(scores, dtype=float)
    if arr.size == 0:
        return arr
    lo, hi = float(np.min(arr)), float(np.max(arr))
    if hi - lo < 1e-9:
        return np.ones_like(arr)
    norm = (arr - lo) / (hi - lo)
    # floor so every mutation contributes at least 5% of max weight
    return np.clip(norm, 0.05, 1.0)


# ---------------------------------------------------------------------------
# Atlas API client (requires ALPHAGENOME_API_KEY)
# ---------------------------------------------------------------------------

def fetch_avi_via_atlas(variants: Sequence[Tuple[str, int, str, str]],
                        api_key: Optional[str] = None,
                        batch_size: int = 200,
                        max_workers: int = 8,
                        cache_path: Optional[Path] = None) -> Dict[str, AVIWeight]:
    """Fetch AVI scores from the AlphaGenome Atlas gRPC API.

    Args:
        variants: list of (chrom, pos, ref, alt) tuples (chrom may include 'chr').
        api_key: AlphaGenome API key. Falls back to $ALPHAGENOME_API_KEY.
        batch_size: variants per request.
        max_workers: parallel workers for batch queries.
        cache_path: optional JSON path to read/write a persistent cache.

    Returns: mapping variant_key -> AVIWeight (source='atlas_api').
    """
    key = api_key or os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        raise RuntimeError(
            "ALPHAGENOME_API_KEY not set. Obtain a free non-commercial key at "
            "https://deepmind.google.com/science/alphagenome and either export "
            "the env var or pass api_key=..."
        )
    try:
        from alphagenome.atlas import atlas as atlas_mod
        from alphagenome.data import genome
    except ImportError as e:
        raise RuntimeError(
            "alphagenome>=0.9.0 is required for Atlas API access. "
            "Install with: pip install alphagenome"
        ) from e

    client = atlas_mod.create(api_key=key)
    out: Dict[str, AVIWeight] = {}
    t0 = time.time()

    def to_variant(chrom: str, pos: int, ref: str, alt: str):
        chrom_canon = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
        return genome.Variant(
            chromosome=chrom_canon,
            position=int(pos),
            reference_bases=ref,
            alternate_bases=alt,
        )

    # Atlas supports query_variants in one round-trip; we batch anyway for
    # memory bounds.
    for i in range(0, len(variants), batch_size):
        chunk = variants[i : i + batch_size]
        var_objs = [to_variant(*v) for v in chunk]
        try:
            res = client.query_variants(
                variants=var_objs,
                requested_scorers=["AVI"],
                progress_bar=False,
                max_workers=max_workers,
            )
        except Exception as e:
            for chrom, pos, ref, alt in chunk:
                vk = variant_key(chrom, pos, ref, alt)
                out[vk] = AVIWeight(variant_key=vk, avi_score=0.0, avi_norm=0.0,
                                    source="missing", error=str(e))
            continue
        adata = res.get("AVI") if isinstance(res, Mapping) else None
        if adata is None:
            for chrom, pos, ref, alt in chunk:
                vk = variant_key(chrom, pos, ref, alt)
                out[vk] = AVIWeight(variant_key=vk, avi_score=0.0, avi_norm=0.0,
                                    source="missing", error="no_avi_dataset")
            continue
        # The AnnData object exposes per-variant scores; index by variant name.
        # adata.X may be a numpy ndarray, a scipy sparse matrix, or None — be defensive.
        x = adata.X
        try:
            if x is None:
                scores = np.zeros(len(var_objs))
            else:
                # scipy sparse matrices expose .toarray(); numpy arrays do not
                toarray = getattr(x, "toarray", None)
                if callable(toarray):
                    scores = np.asarray(toarray()).flatten()
                else:
                    scores = np.asarray(x).flatten()
        except Exception:
            scores = np.asarray(x).flatten()
        names = list(adata.obs_names) if adata.obs_names is not None else [
            variant_key(*c) for c in chunk
        ]
        for name, score in zip(names, scores):
            out[str(name)] = AVIWeight(
                variant_key=str(name),
                avi_score=float(score),
                avi_norm=0.0,  # set in normalize step
                source="atlas_api",
            )

    elapsed = time.time() - t0
    LOG.info("Atlas API: %d variants in %.1fs (%.1f var/s)",
             len(out), elapsed, len(out) / max(elapsed, 1e-6))

    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump({k: asdict(v) for k, v in out.items()}, f, indent=2)

    return out


# ---------------------------------------------------------------------------
# Tabix local cache
# ---------------------------------------------------------------------------

def fetch_avi_via_tabix(variants: Sequence[Tuple[str, int, str, str]],
                        tabix_path: Path) -> Dict[str, AVIWeight]:
    """Look up AVI scores from a local Tabix-indexed avi_snvs_*.tsv.gz.

    File format expected (Tabix-indexed, one row per SNV):
        #chrom  pos  ref  alt  AVI
        chr1     100  A    T    23.4
        ...

    This is the on-disk format the AlphaGenome downloads page exposes once
    extracted from the 88.5 GB zip. We do not download the zip here — the
    user is expected to have extracted it to tabix_path.
    """
    import pysam  # type: ignore

    tbx = pysam.TabixFile(str(tabix_path))
    out: Dict[str, AVIWeight] = {}
    for chrom, pos, ref, alt in variants:
        chrom_canon = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
        vk = f"chr{chrom_canon}:{pos}:{ref}>{alt}"
        try:
            rows = tbx.fetch(chrom_canon, int(pos) - 1, int(pos))
        except Exception as e:
            out[vk] = AVIWeight(variant_key=vk, avi_score=0.0, avi_norm=0.0,
                                source="missing", error=str(e))
            continue
        match_score: Optional[float] = None
        for row in rows:
            fields = row.split("\t")
            if len(fields) < 5:
                continue
            if fields[2] == ref and fields[3] == alt:
                try:
                    match_score = float(fields[4])
                except ValueError:
                    continue
                break
        if match_score is None:
            out[vk] = AVIWeight(variant_key=vk, avi_score=0.0, avi_norm=0.0,
                                source="missing", error="not_in_tabix")
        else:
            out[vk] = AVIWeight(variant_key=vk, avi_score=match_score,
                                avi_norm=0.0, source="tabix_local")
    return out


# ---------------------------------------------------------------------------
# High-level: weight a cohort of TCGA mutations
# ---------------------------------------------------------------------------

def load_tcga_mutations_with_ref_alt(cache_dir: str) -> List[Dict]:
    """Load TCGA-LUAD mutations and ensure each carries ref/alt alleles.

    The default loader in real_tcga_validation.py drops ref/alt to save
    space; for AVI lookup we need both alleles. This re-reads the gzipped
    MAF files.
    """
    cache_path = Path(cache_dir)
    out: List[Dict] = []
    for maf in sorted(cache_path.glob("*.maf.gz")):
        with gzip.open(maf, "rt", errors="replace") as f:
            lines = f.read().splitlines()
        # Find header
        hdr_idx = None
        for i, line in enumerate(lines):
            if line.startswith("Hugo_Symbol"):
                hdr_idx = i
                break
        if hdr_idx is None:
            continue
        header = lines[hdr_idx].split("\t")
        col = {n: header.index(n) for n in (
            "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome",
            "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
            "Variant_Classification", "t_alt_count", "t_ref_count",
        ) if n in header}
        for line in lines[hdr_idx + 1 :]:
            if not line.strip():
                continue
            f_ = line.split("\t")
            try:
                t_alt = int(f_[col["t_alt_count"]])
                t_ref = int(f_[col["t_ref_count"]])
            except (KeyError, ValueError, IndexError):
                continue
            if t_alt + t_ref < 10:
                continue
            chrom = f_[col["Chromosome"]]
            ref = f_[col["Reference_Allele"]]
            alt = f_[col["Tumor_Seq_Allele2"]]
            if not ref or not alt or len(ref) != 1 or len(alt) != 1:
                # AVI Atlas SNV table covers SNVs only; skip indels here
                continue
            vc = f_[col.get("Variant_Classification", -1)] if "Variant_Classification" in col else ""
            out.append({
                "gene": f_[col["Hugo_Symbol"]],
                "sample": f_[col["Tumor_Sample_Barcode"]][:12],
                "chrom": chrom,
                "pos": int(f_[col["Start_Position"]]),
                "ref": ref,
                "alt": alt,
                "t_alt": t_alt,
                "t_depth": t_alt + t_ref,
                "tumor_vaf": t_alt / (t_alt + t_ref),
                "variant_class": vc,
            })
    return out


def weight_cohort(
    cohort_mutations: Sequence[Dict],
    *,
    api_key: Optional[str] = None,
    tabix_path: Optional[Path] = None,
    cache_path: Optional[Path] = None,
) -> Tuple[Dict[str, AVIWeight], str]:
    """Attach AVI weights to every TCGA mutation in the cohort.

    Returns: (mapping variant_key -> AVIWeight, primary_source)
        primary_source is 'atlas_api', 'tabix_local', 'proxy', or
        'mixed'. The result mapping covers every input mutation (missing
        entries are marked source='missing' with a non-zero error string).
    """
    variants: List[Tuple[str, int, str, str]] = [
        (m["chrom"], m["pos"], m["ref"], m["alt"]) for m in cohort_mutations
    ]

    # Try live Atlas API first if key is available
    if api_key or os.environ.get("ALPHAGENOME_API_KEY"):
        try:
            weights = fetch_avi_via_atlas(variants, api_key=api_key,
                                          cache_path=cache_path)
            return weights, "atlas_api"
        except Exception as e:
            LOG.warning("Atlas API path failed: %s — falling back", e)

    # Then local Tabix
    if tabix_path is not None and tabix_path.exists():
        try:
            weights = fetch_avi_via_tabix(variants, tabix_path)
            n_real = sum(1 for w in weights.values() if w.source == "tabix_local")
            if n_real > 0:
                return weights, "tabix_local"
        except Exception as e:
            LOG.warning("Tabix path failed: %s — falling back", e)

    # Fallback: deterministic proxy (clearly flagged)
    weights = {}
    for m in cohort_mutations:
        vk = variant_key(m["chrom"], m["pos"], m["ref"], m["alt"])
        weights[vk] = AVIWeight(
            variant_key=vk,
            avi_score=proxy_avi_for_variant(m.get("variant_class", "")),
            avi_norm=0.0,
            source="proxy",
        )
    return weights, "proxy"


# ---------------------------------------------------------------------------
# CLI sanity-check
# ---------------------------------------------------------------------------

if __name__ == "__main__":  # pragma: no cover
    import argparse, sys
    p = argparse.ArgumentParser(description="Smoke-test AlphaGenome AVI weighting.")
    p.add_argument("--cache-dir", default="validation/tcga/tcga_cache")
    p.add_argument("--api-key", default=None)
    p.add_argument("--tabix", default=None)
    p.add_argument("--limit", type=int, default=10)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    muts = load_tcga_mutations_with_ref_alt(args.cache_dir)
    print(f"Loaded {len(muts)} TCGA mutations with ref/alt")
    sample = muts[: args.limit]
    weights, source = weight_cohort(sample, api_key=args.api_key,
                                    tabix_path=Path(args.tabix) if args.tabix else None)
    print(f"Primary source: {source}")
    for m, vk in zip(sample, [variant_key(m['chrom'], m['pos'], m['ref'], m['alt']) for m in sample]):
        w = weights.get(vk)
        if w is None:
            print(f"  {vk}  MISSING")
        else:
            print(f"  {vk}  AVI={w.avi_score:6.2f}  source={w.source}  "
                  f"gene={m['gene']:>10}  class={m.get('variant_class','')}")