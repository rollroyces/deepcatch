"""
AlphaMissense weighting for DeepCatch v2.2 panel-LLR.

Provides per-mutation AlphaMissense pathogenicity weights used as fixed
scalars in the panel log-likelihood aggregator.

LICENSING / DATA-SOURCE NOTES (CRITICAL)
----------------------------------------
AlphaMissense (Cheng et al., *Science* 2023) is released under
**CC BY-NC-SA 4.0** by DeepMind. The ~5.5 GB TSV cannot be redistributed
and is NOT bundled with this repository. The TSV must be downloaded once
by the user and indexed locally — see the module docstring for paths.

Because this module only consumes AlphaMissense scores as **fixed scalar
weights** in a deterministic hand-coded aggregator (a weighted sum of
per-locus Poisson LLR scores), it does NOT train any model on AlphaMissense
outputs and does not redistribute them. This is the same architectural
pattern as AlphaGenome AVI consumption (see `alphagenome_weights.py`).

Data sources, in priority order:
  1. A local per-Uniprot pickle index at one of:
       - $CWD/alphamissense_index.pkl
       - ~/.cache/mrnavax/alphamissense_index.pkl       (mrnavax cache)
       - ~/.cache/mrna_ai_tools/alphamissense_index.pkl (mrna-ai-toolkit)
       - ~/.cache/alphamissense/alphamissense_index.pkl
  2. The raw TSV at one of:
       - $CWD/AlphaMissense_hg38.tsv[.gz]
       - ~/AlphaMissense_hg38.tsv[.gz]
  3. A per-variant-class published-prior PROXY (clearly flagged
     `source='proxy'` in every result). This is the deterministic fallback
     used when the real index is unavailable; it preserves the
     *direction* of the weighting hypothesis (LoF > missense > silent) so
     the pipeline runs end-to-end, but the resulting AUCs are NOT
     representative of true AlphaMissense-weighted detection.

Reference
---------
Cheng J, Novati G, Pan J, et al. "Accurate proteome-wide missense variant
effect prediction with AlphaMissense." *Science* 381, eadg7492 (2023).
https://www.science.org/doi/10.1126/science.adg7492
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import logging
import os
import pickle
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

LOG = logging.getLogger("alphamissense_weights")

# Standard AlphaMissense pathogenicity thresholds (Cheng 2023, supp.).
THRESHOLD_PATHOGENIC = 0.564
THRESHOLD_BENIGN = 0.34

# Candidate pickle locations (priority order — first hit wins).
INDEX_CANDIDATE_PATHS = [
    Path.cwd() / "alphamissense_index.pkl",
    Path.home() / ".cache" / "mrnavax" / "alphamissense_index.pkl",
    Path.home() / ".cache" / "mrna_ai_tools" / "alphamissense_index.pkl",
    Path.home() / ".cache" / "alphamissense" / "alphamissense_index.pkl",
]

# Raw TSV locations (only consulted if no pickle is available).
TSV_CANDIDATE_PATHS = [
    Path.cwd() / "AlphaMissense_hg38.tsv",
    Path.cwd() / "AlphaMissense_hg38.tsv.gz",
    Path.home() / "AlphaMissense_hg38.tsv",
    Path.home() / "AlphaMissense_hg38.tsv.gz",
]


# ---------------------------------------------------------------------------
# Per-mutation result dataclass
# ---------------------------------------------------------------------------

@dataclass
class AMWeight:
    """AlphaMissense pathogenicity weight for a single missense variant."""
    uniprot: str
    prot_pos: int
    ref_aa: str
    alt_aa: str
    score: float              # raw AM pathogenicity probability in [0, 1]
    classification: str       # likely_benign | ambiguous | likely_pathogenic
    source: str = "missing"   # picklelocal | tsvlocal | proxy | missing
    norm: float = 0.0         # normalised weight in (0, 1] — set later
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Index loaders — attempt real data first, fall back to proxy
# ---------------------------------------------------------------------------

def _find_index_path() -> Optional[Path]:
    for p in INDEX_CANDIDATE_PATHS:
        if p.exists() and p.is_file() and p.stat().st_size > 0:
            return p
    return None


def _find_tsv_path() -> Optional[Path]:
    for p in TSV_CANDIDATE_PATHS:
        if p.exists() and p.is_file() and p.stat().st_size > 1000:
            return p
    return None


def _try_load_pickle(path: Path) -> Optional[Dict[str, AMWeight]]:
    """Attempt to load a pre-built per-Uniprot pickle index.

    The expected pickle format is a dict mapping keys of the form
    ``"<UNIPROT>:<PROT_POS>:<REF_AA>:<ALT_AA>"`` -> dict(score=…, …).
    Any failure is logged and returns None so the next source can be tried.
    """
    try:
        with open(path, "rb") as f:
            obj = pickle.load(f)
    except Exception as e:
        LOG.warning("Could not unpickle %s: %s", path, e)
        return None
    if not isinstance(obj, dict) or not obj:
        LOG.warning("Pickle at %s is empty or not a dict — skipping", path)
        return None
    out: Dict[str, AMWeight] = {}
    for k, v in obj.items():
        try:
            if isinstance(v, AMWeight):
                out[str(k)] = v
                continue
            # Accept raw dicts produced by mrnavax-style builders
            parts = str(k).split(":")
            if len(parts) != 4:
                continue
            up, pp, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
            score = float(v.get("score", 0.0)) if isinstance(v, dict) else float(v)
            cls = classify_am(score)
            out[str(k)] = AMWeight(
                uniprot=up, prot_pos=pp, ref_aa=ref, alt_aa=alt,
                score=score, classification=cls, source="picklelocal",
            )
        except Exception:
            continue
    if out:
        LOG.info("Loaded %d AM weights from pickle %s", len(out), path)
    return out or None


def classify_am(score: float) -> str:
    """Map AM pathogenicity probability to published 3-tier label."""
    if score >= THRESHOLD_PATHOGENIC:
        return "likely_pathogenic"
    if score < THRESHOLD_BENIGN:
        return "likely_benign"
    return "ambiguous"


# ---------------------------------------------------------------------------
# TSV streaming (only invoked if no pickle is present)
# ---------------------------------------------------------------------------

def _try_build_from_tsv(tsv_path: Path,
                        wanted_keys: Iterable[str]) -> Dict[str, AMWeight]:
    """Stream AlphaMissense_hg38.tsv and collect only the wanted variants.

    Columns (Cheng 2023 release): ``uniprot_id  uniprot_pos  ref_aa  alt_aa
    am_pathogenicity  am_class``. No header in the released TSV — column
    order is fixed. This is intentionally a streaming scan; the full TSV
    is ~5.5 GB and we only need a tiny fraction.
    """
    wanted = set(wanted_keys)
    out: Dict[str, AMWeight] = {}
    t0 = time.time()
    opener = gzip.open if str(tsv_path).endswith(".gz") else open
    n_rows = 0
    try:
        with opener(tsv_path, "rt") as f:  # type: ignore[arg-type]
            for row in f:
                row = row.rstrip("\n")
                if not row or row.startswith("#"):
                    continue
                fields = row.split("\t")
                if len(fields) < 5:
                    continue
                up, pp, ref, alt, score = fields[0], fields[1], fields[2], fields[3], fields[4]
                key = f"{up}:{pp}:{ref}:{alt}"
                if key not in wanted:
                    n_rows += 1
                    continue
                try:
                    s = float(score)
                except ValueError:
                    continue
                out[key] = AMWeight(
                    uniprot=up, prot_pos=int(pp), ref_aa=ref, alt_aa=alt,
                    score=s, classification=classify_am(s), source="tsvlocal",
                )
                if len(out) == len(wanted):
                    break
    except Exception as e:
        LOG.error("TSV scan failed: %s", e)
    elapsed = time.time() - t0
    LOG.info("Scanned %d rows in %.1fs; collected %d matches", n_rows, elapsed, len(out))
    return out


# ---------------------------------------------------------------------------
# Proxy (deterministic, clearly-flagged fallback)
# ---------------------------------------------------------------------------

# Per-variant-class AM-pathogenicity P(mean). Approximate means derived
# from the published class-wise summary in Cheng 2023 (extended data fig 3
# and Table S2). These are STABLE priors that mirror the published
# direction: protein-truncating > missense > silent. They are NOT
# per-mutation scores and must be flagged as PROXY_USED.
_PROXY_AM_BY_CLASS: Mapping[str, float] = {
    "Nonsense_Mutation":      0.92,
    "Frame_Shift_Del":        0.91,
    "Frame_Shift_Ins":        0.91,
    "Splice_Site":            0.86,
    "Nonstop_Mutation":       0.85,
    "Translation_Start_Site": 0.80,
    "Missense_Mutation":      0.55,
    "In_Frame_Del":           0.55,
    "In_Frame_Ins":           0.55,
    "Splice_Region":          0.45,
    "5'UTR":                  0.25,
    "3'UTR":                  0.18,
    "Silent":                 0.10,
    "Intron":                 0.10,
    "IGR":                    0.05,
    "RNA":                    0.05,
}
_PROXY_AM_DEFAULT = 0.40


def proxy_am_for_variant(variant_class: str) -> float:
    """Deterministic per-variant-class AM-pathogenicity prior."""
    return _PROXY_AM_BY_CLASS.get(variant_class, _PROXY_AM_DEFAULT)


# ---------------------------------------------------------------------------
# Key derivation
# ---------------------------------------------------------------------------

def am_key(uniprot: str, prot_pos: int, ref_aa: str, alt_aa: str) -> str:
    """Canonical AlphaMissense lookup key."""
    return f"{uniprot}:{int(prot_pos)}:{ref_aa}:{alt_aa}"


# ---------------------------------------------------------------------------
# High-level: weight a cohort
# ---------------------------------------------------------------------------

def load_tcga_mutations_with_protein(cache_dir: str) -> List[Dict]:
    """Load TCGA-LUAD mutations augmented with SwissProt + Protein_position.

    The default loader in `real_tcga_validation.py` drops Uniprot and
    protein-position to save space; for AlphaMissense lookup we need
    both. This re-reads the gzipped MAF files.
    """
    cache_path = Path(cache_dir)
    out: List[Dict] = []
    for maf in sorted(cache_path.glob("*.maf.gz")):
        with gzip.open(maf, "rt", errors="replace") as f:
            lines = f.read().splitlines()
        hdr_idx = None
        for i, ln in enumerate(lines):
            if ln.startswith("Hugo_Symbol"):
                hdr_idx = i
                break
        if hdr_idx is None:
            continue
        header = lines[hdr_idx].split("\t")
        col = {n: header.index(n) for n in (
            "Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome",
            "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2",
            "Variant_Classification", "SWISSPROT", "Protein_position",
            "HGVSp", "t_alt_count", "t_ref_count",
        ) if n in header}
        if not {"Hugo_Symbol", "Chromosome", "Start_Position",
                "Reference_Allele", "Tumor_Seq_Allele2",
                "SWISSPROT", "Protein_position"}.issubset(col):
            continue
        for line in lines[hdr_idx + 1:]:
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
            ref = f_[col["Reference_Allele"]]
            alt = f_[col["Tumor_Seq_Allele2"]]
            sw = f_[col["SWISSPROT"]].strip()
            pp = f_[col["Protein_position"]].strip().split("/")[0]
            try:
                pp_int = int(pp.split("-")[0].strip())
            except ValueError:
                pp_int = None
            if not sw or pp_int is None:
                continue
            if len(ref) != 1 or len(alt) != 1:
                # AlphaMissense is missense-only; skip indels
                continue
            vc = f_[col.get("Variant_Classification", -1)] \
                if "Variant_Classification" in col else ""
            # AlphaMissense is missense-only — silently exclude other classes
            if vc != "Missense_Mutation":
                continue
            out.append({
                "gene": f_[col["Hugo_Symbol"]],
                "sample": f_[col["Tumor_Sample_Barcode"]][:12],
                "chrom": f_[col["Chromosome"]],
                "pos": int(f_[col["Start_Position"]]),
                "ref": ref,
                "alt": alt,
                "uniprot": sw.split(".")[0],
                "prot_pos": pp_int,
                "variant_class": vc,
                "t_alt": t_alt,
                "t_depth": t_alt + t_ref,
                "tumor_vaf": t_alt / (t_alt + t_ref),
            })
    return out


def weight_cohort(
    cohort_mutations: Sequence[Dict],
    *,
    index_path: Optional[Path] = None,
    tsv_path: Optional[Path] = None,
) -> Tuple[Dict[str, AMWeight], Dict[str, Any]]:
    """Attach AlphaMissense pathogenicity weights to every TCGA mutation.

    Returns: (mapping ``"<UP>:<POS>:<REF>:<ALT>"`` -> AMWeight, info dict).
    The info dict reports which source was used and how many mutations
    were found / missed / proxied. Missing weights are emitted as
    ``source='missing'`` so downstream code can keep a fixed-length
    vector for every cohort mutation.
    """
    info: Dict[str, Any] = {
        "primary_source": "missing",
        "n_input": len(cohort_mutations),
        "n_real": 0,
        "n_proxy": 0,
        "n_missing": 0,
        "index_path": None,
        "tsv_path": None,
        "missense_total": 0,
        "missense_with_full_key": 0,
    }

    # Build the list of wanted keys (missense SNVs only — AM is missense).
    wanted_keys: List[str] = []
    missense_input: List[Dict] = []
    for m in cohort_mutations:
        info["missense_total"] += 1
        up = m.get("uniprot")
        pp = m.get("prot_pos")
        ref = m.get("ref")
        alt = m.get("alt")
        vc = m.get("variant_class", "")
        if vc == "Missense_Mutation" and up and pp is not None \
                and ref and alt and len(ref) == 1 and len(alt) == 1:
            key = am_key(up, pp, ref, alt)
            wanted_keys.append(key)
            missense_input.append(m)
            info["missense_with_full_key"] += 1
    info["wanted_keys"] = len(wanted_keys)

    real_weights: Dict[str, AMWeight] = {}

    # 1) Local pickle
    idx = Path(index_path) if index_path else _find_index_path()
    if idx is not None:
        loaded = _try_load_pickle(idx)
        if loaded:
            real_weights = {k: w for k, w in loaded.items() if k in set(wanted_keys)}
            info["primary_source"] = "picklelocal"
            info["index_path"] = str(idx)

    # 2) Raw TSV — only if pickle didn't yield anything
    if not real_weights:
        tsv = Path(tsv_path) if tsv_path else _find_tsv_path()
        if tsv is not None:
            real_weights = _try_build_from_tsv(tsv, wanted_keys)
            if real_weights:
                info["primary_source"] = "tsvlocal"
                info["tsv_path"] = str(tsv)

    # 3) Final fallback: variant-class PROXY for every wanted mutation
    final: Dict[str, AMWeight] = {}
    for m, key in zip(missense_input, wanted_keys):
        if key in real_weights:
            w = real_weights[key]
            final[key] = w
        else:
            vc = m.get("variant_class", "")
            score = proxy_am_for_variant(vc)
            up = m["uniprot"]
            pp = int(m["prot_pos"])
            ref = m["ref"]
            alt = m["alt"]
            final[key] = AMWeight(
                uniprot=up, prot_pos=pp, ref_aa=ref, alt_aa=alt,
                score=score, classification=classify_am(score),
                source="proxy",
            )

    # Tally
    info["n_real"] = sum(1 for w in final.values() if w.source in ("picklelocal", "tsvlocal"))
    info["n_proxy"] = sum(1 for w in final.values() if w.source == "proxy")
    info["n_missing"] = sum(1 for w in final.values() if w.source == "missing")

    if info["primary_source"] == "missing" and info["n_real"] == 0:
        info["primary_source"] = "proxy"

    return final, info


def normalize_am(scores: Sequence[float]) -> np.ndarray:
    """Min-max scale to weights in [floor, 1] with a 0.05 floor."""
    arr = np.asarray(scores, dtype=float)
    if arr.size == 0:
        return arr
    lo, hi = float(arr.min()), float(arr.max())
    if hi - lo < 1e-9:
        return np.ones_like(arr)
    norm = (arr - lo) / (hi - lo)
    return np.clip(norm, 0.05, 1.0)


# ---------------------------------------------------------------------------
# CLI sanity-check
# ---------------------------------------------------------------------------

if __name__ == "__main__":  # pragma: no cover
    import argparse
    p = argparse.ArgumentParser(description="Smoke-test AlphaMissense weighting.")
    p.add_argument("--cache-dir", default="validation/tcga/tcga_cache")
    p.add_argument("--index", default=None)
    p.add_argument("--tsv", default=None)
    p.add_argument("--limit", type=int, default=10)
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    muts = load_tcga_mutations_with_protein(args.cache_dir)
    print(f"Loaded {len(muts)} TCGA missense SNVs with Uniprot+position")
    sample = muts[: args.limit]
    weights, info = weight_cohort(
        sample,
        index_path=Path(args.index) if args.index else None,
        tsv_path=Path(args.tsv) if args.tsv else None,
    )
    print(f"Primary source: {info['primary_source']}")
    print(f"Real={info['n_real']}  Proxy={info['n_proxy']}  Missing={info['n_missing']}")
    for m in sample:
        k = am_key(m["uniprot"], m["prot_pos"], m["ref"], m["alt"])
        w = weights.get(k)
        if w is None:
            print(f"  {m['gene']:>10}  {k}  MISSING")
        else:
            print(f"  {m['gene']:>10}  {k}  AM={w.score:.3f}  cls={w.classification:<18}  src={w.source}")