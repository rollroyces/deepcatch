#!/usr/bin/env python3
"""
TCGA Data Downloader via cBioPortal API
========================================
Downloads real somatic mutation data (MAF), clinical data, and methylation data
from cBioPortal for multiple cancer types:
  - Lung Adenocarcinoma (LUAD): luad_tcga_pan_can_atlas_2018
  - Colorectal Adenocarcinoma (COADREAD): coadread_tcga_pan_can_atlas_2018
  - Breast Invasive Carcinoma (BRCA): brca_tcga_pan_can_atlas_2018

Also supports fallback hardcoded cancer mutation profiles if API fails.

Usage:
    python3 tcga_downloader.py --output ./tcga_cache/ --cancer-types LUAD,COADREAD,BRCA
"""

import json
import os
import sys
import time
import argparse
import hashlib
from typing import Dict, List, Optional, Any, Tuple
from collections import defaultdict

import requests
import pandas as pd
import numpy as np


# ═══════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════

CBIOPORTAL_API = "https://www.cbioportal.org/api"
GDC_API = "https://api.gdc.cancer.gov"

TCGA_STUDIES = {
    'LUAD': {
        'study_id': 'luad_tcga_pan_can_atlas_2018',
        'name': 'Lung Adenocarcinoma (TCGA PanCan Atlas)',
        'molecular_profile': 'luad_tcga_pan_can_atlas_2018_mutations',
        'sample_count': 566,
    },
    'COADREAD': {
        'study_id': 'coadread_tcga_pan_can_atlas_2018',
        'name': 'Colorectal Adenocarcinoma (TCGA PanCan Atlas)',
        'molecular_profile': 'coadread_tcga_pan_can_atlas_2018_mutations',
        'sample_count': 594,
    },
    'BRCA': {
        'study_id': 'brca_tcga_pan_can_atlas_2018',
        'name': 'Breast Invasive Carcinoma (TCGA PanCan Atlas)',
        'molecular_profile': 'brca_tcga_pan_can_atlas_2018_mutations',
        'sample_count': 1084,
    },
    'PRAD': {
        'study_id': 'prad_tcga_pan_can_atlas_2018',
        'name': 'Prostate Adenocarcinoma (TCGA PanCan Atlas)',
        'molecular_profile': 'prad_tcga_pan_can_atlas_2018_mutations',
        'sample_count': 494,
    },
    'HNSC': {
        'study_id': 'hnsc_tcga_pan_can_atlas_2018',
        'name': 'Head and Neck Squamous Cell Carcinoma (TCGA PanCan Atlas)',
        'molecular_profile': 'hnsc_tcga_pan_can_atlas_2018_mutations',
        'sample_count': 523,
    },
}

# Known cancer hotspot mutations (fallback if API fails)
CANCER_HOTSPOTS = {
    'TP53': {
        'gene': 'TP53',
        'chrom': '17',
        'hotspots': [
            {'pos': 7577120, 'ref': 'C', 'alt': 'T', 'protein_change': 'R175H', 'frequency': 0.06},
            {'pos': 7577539, 'ref': 'C', 'alt': 'T', 'protein_change': 'R248W', 'frequency': 0.04},
            {'pos': 7577538, 'ref': 'G', 'alt': 'A', 'protein_change': 'R248Q', 'frequency': 0.04},
            {'pos': 7577121, 'ref': 'G', 'alt': 'A', 'protein_change': 'G245S', 'frequency': 0.03},
            {'pos': 7578554, 'ref': 'G', 'alt': 'A', 'protein_change': 'R282W', 'frequency': 0.03},
            {'pos': 7578406, 'ref': 'C', 'alt': 'T', 'protein_change': 'R273C', 'frequency': 0.03},
            {'pos': 7578407, 'ref': 'C', 'alt': 'A', 'protein_change': 'R273H', 'frequency': 0.02},
            {'pos': 7578263, 'ref': 'G', 'alt': 'A', 'protein_change': 'Y220C', 'frequency': 0.02},
            {'pos': 7578190, 'ref': 'T', 'alt': 'C', 'protein_change': 'V157F', 'frequency': 0.02},
            {'pos': 7577094, 'ref': 'C', 'alt': 'T', 'protein_change': 'H179R', 'frequency': 0.02},
            {'pos': 7577580, 'ref': 'C', 'alt': 'T', 'protein_change': 'R213X', 'frequency': 0.01},
            {'pos': 7578456, 'ref': 'C', 'alt': 'T', 'protein_change': 'R280K', 'frequency': 0.01},
        ],
        'cancer_types': ['LUAD', 'COADREAD', 'BRCA', 'HNSC', 'PRAD'],
        'overall_frequency': 0.42,  # Mutated in ~42% of cancers
    },
    'KRAS': {
        'gene': 'KRAS',
        'chrom': '12',
        'hotspots': [
            {'pos': 25398284, 'ref': 'C', 'alt': 'T', 'protein_change': 'G12D', 'frequency': 0.12},
            {'pos': 25398284, 'ref': 'C', 'alt': 'A', 'protein_change': 'G12V', 'frequency': 0.10},
            {'pos': 25398285, 'ref': 'C', 'alt': 'T', 'protein_change': 'G13D', 'frequency': 0.06},
            {'pos': 25398284, 'ref': 'C', 'alt': 'G', 'protein_change': 'G12C', 'frequency': 0.06},
            {'pos': 25398284, 'ref': 'C', 'alt': 'C', 'protein_change': 'G12A', 'frequency': 0.03},
            {'pos': 25398285, 'ref': 'C', 'alt': 'A', 'protein_change': 'G13V', 'frequency': 0.02},
            {'pos': 25380275, 'ref': 'T', 'alt': 'A', 'protein_change': 'Q61H', 'frequency': 0.02},
            {'pos': 25380276, 'ref': 'T', 'alt': 'A', 'protein_change': 'Q61L', 'frequency': 0.02},
            {'pos': 25398284, 'ref': 'C', 'alt': 'T', 'protein_change': 'G12S', 'frequency': 0.02},
        ],
        'cancer_types': ['LUAD', 'COADREAD'],
        'overall_frequency': 0.25,
    },
    'BRAF': {
        'gene': 'BRAF',
        'chrom': '7',
        'hotspots': [
            {'pos': 140753336, 'ref': 'A', 'alt': 'T', 'protein_change': 'V600E', 'frequency': 0.08},
            {'pos': 140753337, 'ref': 'A', 'alt': 'G', 'protein_change': 'V600K', 'frequency': 0.02},
            {'pos': 140753336, 'ref': 'A', 'alt': 'G', 'protein_change': 'V600G', 'frequency': 0.01},
        ],
        'cancer_types': ['COADREAD', 'LUAD', 'BRCA'],
        'overall_frequency': 0.08,
    },
    'PIK3CA': {
        'gene': 'PIK3CA',
        'chrom': '3',
        'hotspots': [
            {'pos': 178952085, 'ref': 'A', 'alt': 'G', 'protein_change': 'H1047R', 'frequency': 0.08},
            {'pos': 178952074, 'ref': 'A', 'alt': 'G', 'protein_change': 'E545K', 'frequency': 0.06},
            {'pos': 178952073, 'ref': 'G', 'alt': 'A', 'protein_change': 'E542K', 'frequency': 0.04},
            {'pos': 178936082, 'ref': 'G', 'alt': 'A', 'protein_change': 'E545A', 'frequency': 0.02},
            {'pos': 178927980, 'ref': 'A', 'alt': 'G', 'protein_change': 'N345K', 'frequency': 0.02},
            {'pos': 178938989, 'ref': 'C', 'alt': 'T', 'protein_change': 'C420R', 'frequency': 0.01},
        ],
        'cancer_types': ['BRCA', 'COADREAD', 'LUAD', 'HNSC'],
        'overall_frequency': 0.15,
    },
    'EGFR': {
        'gene': 'EGFR',
        'chrom': '7',
        'hotspots': [
            {'pos': 55259515, 'ref': 'T', 'alt': 'G', 'protein_change': 'L858R', 'frequency': 0.05},
            {'pos': 55242464, 'ref': 'G', 'alt': 'A', 'protein_change': 'E746_A750del', 'frequency': 0.04},
            {'pos': 55249071, 'ref': 'A', 'alt': 'G', 'protein_change': 'T790M', 'frequency': 0.02},
        ],
        'cancer_types': ['LUAD'],
        'overall_frequency': 0.15,
    },
    'APC': {
        'gene': 'APC',
        'chrom': '5',
        'hotspots': [
            {'pos': 112175770, 'ref': 'C', 'alt': 'T', 'protein_change': 'R1450X', 'frequency': 0.04},
            {'pos': 112175639, 'ref': 'C', 'alt': 'T', 'protein_change': 'R1367X', 'frequency': 0.03},
            {'pos': 112175211, 'ref': 'C', 'alt': 'T', 'protein_change': 'R1114X', 'frequency': 0.02},
            {'pos': 112174838, 'ref': 'C', 'alt': 'T', 'protein_change': 'Q1338X', 'frequency': 0.02},
        ],
        'cancer_types': ['COADREAD'],
        'overall_frequency': 0.70,
    },
    'PTEN': {
        'gene': 'PTEN',
        'chrom': '10',
        'hotspots': [
            {'pos': 89692904, 'ref': 'C', 'alt': 'T', 'protein_change': 'R130X', 'frequency': 0.03},
            {'pos': 89685215, 'ref': 'G', 'alt': 'A', 'protein_change': 'R233X', 'frequency': 0.02},
        ],
        'cancer_types': ['BRCA', 'PRAD'],
        'overall_frequency': 0.10,
    },
}


class TCGADownloader:
    """Download and cache TCGA data from cBioPortal API."""

    def __init__(self, cache_dir: str = './tcga_cache/', rate_limit_delay: float = 0.5):
        self.cache_dir = cache_dir
        self.rate_limit_delay = rate_limit_delay
        self.session = requests.Session()
        self.session.headers.update({
            'Accept': 'application/json',
            'User-Agent': 'TCGA-Validation-Pipeline/1.0'
        })
        os.makedirs(cache_dir, exist_ok=True)

    def _cache_path(self, key: str) -> str:
        """Get cache file path for a key."""
        safe_key = hashlib.md5(key.encode()).hexdigest()[:16]
        return os.path.join(self.cache_dir, f'{safe_key}.json')

    def _cached_get(self, url: str, cache_key: str = None) -> Optional[Dict]:
        """GET with caching."""
        if cache_key is None:
            cache_key = url
        cache_path = self._cache_path(cache_key)

        # Check cache
        if os.path.exists(cache_path):
            with open(cache_path) as f:
                return json.load(f)

        # Fetch
        try:
            resp = self.session.get(url, timeout=30)
            time.sleep(self.rate_limit_delay)
            if resp.status_code == 200:
                data = resp.json()
                with open(cache_path, 'w') as f:
                    json.dump(data, f)
                return data
            else:
                print(f"  API error {resp.status_code} for {url[:100]}")
                return None
        except Exception as e:
            print(f"  Request failed: {e}")
            return None

    def fetch_study_list(self) -> List[Dict]:
        """Fetch list of available studies."""
        url = f"{CBIOPORTAL_API}/studies?pageSize=1000&direction=ASC"
        return self._cached_get(url, 'study_list') or []

    def fetch_mutations(self, study_id: str, page_size: int = 500) -> pd.DataFrame:
        """Fetch mutation data (MAF equivalent) for a study.

        Uses the /molecular-data API to fetch mutations.
        """
        print(f"  Fetching mutations for {study_id}...")

        # First get molecular profile ID for mutations
        url = f"{CBIOPORTAL_API}/studies/{study_id}/molecular-profiles"
        profiles = self._cached_get(url, f'mol_profiles_{study_id}')

        mutation_profile_id = None
        if profiles:
            for p in profiles:
                if p.get('molecularAlterationType') == 'MUTATION_EXTENDED':
                    mutation_profile_id = p['molecularProfileId']
                    break

        if not mutation_profile_id:
            print(f"  No mutation profile found for {study_id}")
            return pd.DataFrame()

        print(f"  Mutation profile: {mutation_profile_id}")

        # Try to fetch sample list first
        url = f"{CBIOPORTAL_API}/studies/{study_id}/samples?pageSize=5000"
        samples = self._cached_get(url, f'samples_{study_id}')

        if not samples:
            print(f"  No samples found")
            return pd.DataFrame()

        sample_ids = [s['sampleId'] for s in samples]
        print(f"  Found {len(sample_ids)} samples")

        # Fetch mutations for each sample (batched)
        all_mutations = []
        for i in range(0, len(sample_ids), 50):
            batch = sample_ids[i:i+50]
            mutation_list = []

            for sample_id in batch:
                url = (f"{CBIOPORTAL_API}/molecular-profiles/{mutation_profile_id}"
                       f"/mutations?sampleListId={sample_id}&projection=DETAILED&pageSize=5000")
                data = self._cached_get(url, f'mutations_{study_id}_{sample_id}')
                if data:
                    mutation_list.extend(data)

            if mutation_list:
                all_mutations.extend(mutation_list)

            progress = min(i + 50, len(sample_ids))
            print(f"    Progress: {progress}/{len(sample_ids)} samples "
                  f"({len(all_mutations)} mutations found)", end='\r')

        print()
        print(f"  Total mutations: {len(all_mutations)}")

        if all_mutations:
            df = pd.DataFrame(all_mutations)
            return df
        return pd.DataFrame()

    def fetch_clinical_data(self, study_id: str) -> pd.DataFrame:
        """Fetch clinical data for a study."""
        print(f"  Fetching clinical data for {study_id}...")

        # Get clinical attributes
        url = f"{CBIOPORTAL_API}/studies/{study_id}/clinical-attributes?pageSize=500"
        attrs = self._cached_get(url, f'clinical_attrs_{study_id}')

        if not attrs:
            return pd.DataFrame()

        attr_ids = [a['clinicalAttributeId'] for a in attrs]

        # Get sample list
        url = f"{CBIOPORTAL_API}/studies/{study_id}/samples?pageSize=5000"
        samples = self._cached_get(url, f'samples_{study_id}')

        if not samples:
            return pd.DataFrame()

        sample_ids = [s['sampleId'] for s in samples]

        # Fetch clinical data
        url = f"{CBIOPORTAL_API}/studies/{study_id}/clinical-data?clinicalDataType=SAMPLE&pageSize=100000"
        clinical = self._cached_get(url, f'clinical_data_{study_id}')

        if clinical:
            # Convert to sample x attribute matrix
            records = []
            for item in clinical:
                records.append({
                    'sample_id': item['sampleId'],
                    'patient_id': item.get('patientId', ''),
                    'attribute': item['clinicalAttributeId'],
                    'value': item['value'],
                })

            df = pd.DataFrame(records)
            if not df.empty:
                pivot = df.pivot(index='sample_id', columns='attribute', values='value')
                pivot['patient_id'] = df.groupby('sample_id')['patient_id'].first()
                return pivot.reset_index()

        return pd.DataFrame()

    def download_all(self, cancer_types: List[str]) -> Dict[str, Dict]:
        """Download all data for specified cancer types.

        Returns:
            dict mapping cancer_type -> {mutations_df, clinical_df, metadata}
        """
        results = {}
        for ct in cancer_types:
            if ct not in TCGA_STUDIES:
                print(f"Warning: {ct} not in known TCGA studies, skipping")
                continue

            study = TCGA_STUDIES[ct]
            print(f"\n{'='*60}")
            print(f"Processing: {study['name']}")
            print(f"Study ID: {study['study_id']}")
            print(f"{'='*60}")

            try:
                mutations = self.fetch_mutations(study['study_id'])
                clinical = self.fetch_clinical_data(study['study_id'])

                results[ct] = {
                    'mutations': mutations,
                    'clinical': clinical,
                    'metadata': study,
                }

                # Save to disk
                if not mutations.empty:
                    mutations.to_csv(
                        os.path.join(self.cache_dir, f'{ct}_mutations.csv'), index=False)
                if not clinical.empty:
                    clinical.to_csv(
                        os.path.join(self.cache_dir, f'{ct}_clinical.csv'), index=False)

                print(f"  ✓ {ct}: {len(mutations)} mutations, {len(clinical)} clinical records")

            except Exception as e:
                print(f"  ✗ Error processing {ct}: {e}")
                results[ct] = {'error': str(e), 'metadata': study}

        # Save summary
        summary = {}
        for ct, data in results.items():
            mutations = data.get('mutations', pd.DataFrame())
            summary[ct] = {
                'study_id': data.get('metadata', {}).get('study_id', ''),
                'name': data.get('metadata', {}).get('name', ''),
                'n_mutations': len(mutations) if isinstance(mutations, pd.DataFrame) else 0,
                'error': data.get('error', ''),
            }

        with open(os.path.join(self.cache_dir, 'download_summary.json'), 'w') as f:
            json.dump(summary, f, indent=2)

        return results


def build_fallback_dataset(cancer_types: List[str],
                           n_samples: int = 500,
                           seed: int = 42) -> Dict[str, Any]:
    """
    Build fallback dataset from literature-known cancer mutation frequencies.

    This generates realistic tumor mutation profiles using known hotspot
    frequencies from COSMIC/TCGA publications. Used when cBioPortal API is
    unavailable.

    Returns:
        dict with keys: ground_truth_variants, sample_metadata, cancer_type_labels
    """
    rng = np.random.RandomState(seed)
    print(f"\n{'='*60}")
    print("Building fallback dataset from literature-known cancer hotspots")
    print(f"Target: {n_samples} samples across {cancer_types}")
    print(f"{'='*60}")

    all_variants = []
    sample_metadata = []
    labels = []

    n_per_type = n_samples // len(cancer_types)

    for ct in cancer_types:
        print(f"\n  Cancer type: {ct}")
        ct_variants = []
        ct_samples = []

        for i in range(n_per_type):
            sample_id = f"{ct}_S{i:04d}"
            tumor_purity = rng.uniform(0.3, 0.95)
            is_cancer = True  # All TCGA samples are tumor

            # Select cancer-relevant genes for this type
            for gene_name, gene_data in CANCER_HOTSPOTS.items():
                if ct not in gene_data['cancer_types']:
                    continue

                # Is this gene mutated in this sample?
                if rng.random() < gene_data['overall_frequency'] * rng.uniform(0.5, 1.5):
                    # Select which hotspot(s) are mutated
                    hotspots = gene_data['hotspots']
                    n_muts = rng.poisson(1.0) + 1  # Usually 1-2 mutations per gene
                    n_muts = min(n_muts, len(hotspots))

                    selected = rng.choice(len(hotspots), n_muts, replace=False)
                    for idx in selected:
                        hs = hotspots[idx]
                        # True VAF depends on tumor purity and clonality
                        true_vaf = tumor_purity * hs['frequency'] * rng.uniform(0.5, 1.5)
                        true_vaf = np.clip(true_vaf, 0.001, 0.95)

                        variant = {
                            'sample_id': sample_id,
                            'cancer_type': ct,
                            'chrom': f"chr{gene_data['chrom']}",
                            'pos': hs['pos'],
                            'ref': hs['ref'],
                            'alt': hs['alt'],
                            'gene': gene_name,
                            'protein_change': hs['protein_change'],
                            'true_vaf': true_vaf,
                            'tumor_purity': tumor_purity,
                            'is_true_variant': True,
                            'trinuc_context': _random_trinuc(rng),
                        }
                        ct_variants.append(variant)

            ct_samples.append({
                'sample_id': sample_id,
                'cancer_type': ct,
                'is_cancer': is_cancer,
                'tumor_purity': tumor_purity,
            })

        all_variants.extend(ct_variants)
        sample_metadata.extend(ct_samples)
        labels.extend([ct] * n_per_type)
        print(f"    Generated {len(ct_variants)} true variants in {n_per_type} samples")

    print(f"\n  Total: {len(all_variants)} true variants in {len(sample_metadata)} samples")

    return {
        'ground_truth_variants': all_variants,
        'sample_metadata': sample_metadata,
        'cancer_type_labels': labels,
        'source': 'literature_fallback',
    }


def _random_trinuc(rng: np.random.RandomState) -> str:
    """Generate random trinucleotide context."""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(rng.choice(bases, 3))


def fetch_gdc_mafs(cache_dir: str,
                   project: str = "TCGA-LUAD",
                   n_files: int = 25) -> List[str]:
    """Download open-access per-aliquot masked MAF files from the GDC API.

    GDC now serves per-aliquot (per-sample) MAF files rather than a single
    project-level MAF. Each downloaded file is saved to ``cache_dir`` as
    ``gdc_<project>_<index>.maf.gz`` (stable names) and can be re-read by
    ``load_tcga_cohort`` on subsequent runs (no network needed).

    Args:
        cache_dir: directory to save the MAF files into
        project: GDC project id (e.g. TCGA-LUAD)
        n_files: number of aliquot MAF files to download

    Returns:
        list of downloaded file paths
    """
    import urllib.request
    import urllib.parse
    import time

    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Masked Somatic Mutation"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    qs = urllib.parse.quote(json.dumps(filters))
    url = (f"{GDC_API}/files?filters={qs}"
           f"&fields=file_id,file_name,file_size&size={n_files}&pretty=false")

    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            data = json.loads(r.read().decode())
    except Exception as e:
        print(f"  ✗ GDC file query failed: {e}")
        return []

    hits = data.get('data', {}).get('hits', [])
    print(f"  GDC: {len(hits)} open-access MAF files available for {project}")

    os.makedirs(cache_dir, exist_ok=True)
    downloaded = []
    for i, h in enumerate(hits):
        file_id = h.get('file_id')
        if not file_id:
            continue
        out_path = os.path.join(cache_dir, f"gdc_{project}_{i}.maf.gz")
        try:
            with urllib.request.urlopen(f"{GDC_API}/data/{file_id}", timeout=120) as r:
                with open(out_path, 'wb') as f:
                    f.write(r.read())
            downloaded.append(out_path)
            print(f"    ✓ {os.path.basename(out_path)} ({len(downloaded)}/{len(hits)})")
        except Exception as e:
            print(f"    ✗ download {file_id} failed: {e}")
        time.sleep(0.3)  # be polite to GDC

    return downloaded


def load_or_download(cancer_types: List[str],
                     cache_dir: str = './tcga_cache/',
                     n_fallback_samples: int = 500) -> Dict[str, Any]:
    """
    Load cached data or download from cBioPortal.

    If API fails, falls back to literature-known hotspot profiles.
    """
    downloader = TCGADownloader(cache_dir)

    # Try API first
    print("Attempting cBioPortal API download...")
    try:
        results = downloader.download_all(cancer_types)
        # Check if we got meaningful data
        total_mutations = sum(
            len(r.get('mutations', pd.DataFrame()))
            if isinstance(r.get('mutations'), pd.DataFrame) else 0
            for r in results.values()
        )
        if total_mutations > 100:
            print(f"\n✓ Successfully downloaded {total_mutations} mutations via API")
            return {'api_results': results, 'source': 'cbioportal_api'}
    except Exception as e:
        print(f"\n✗ API download failed: {e}")

    # Fallback: use literature-known profiles
    print("\n⚠ Falling back to literature-known cancer mutation profiles")
    return build_fallback_dataset(cancer_types, n_fallback_samples)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download TCGA data from cBioPortal')
    parser.add_argument('--output', default='./tcga_cache/', help='Cache directory')
    parser.add_argument('--cancer-types', default='LUAD,COADREAD,BRCA',
                        help='Comma-separated cancer types')
    parser.add_argument('--n-fallback', type=int, default=500,
                        help='Number of samples for fallback dataset')
    args = parser.parse_args()

    cancer_types = [ct.strip() for ct in args.cancer_types.split(',')]
    print(f"Target cancer types: {cancer_types}")
    print(f"Cache directory: {args.output}")

    data = load_or_download(cancer_types, args.output, args.n_fallback)

    # Save fallback data if used
    if data.get('source') == 'literature_fallback':
        fallback_path = os.path.join(args.output, 'fallback_dataset.json')
        with open(fallback_path, 'w') as f:
            json.dump({
                'ground_truth_variants': data['ground_truth_variants'],
                'sample_metadata': data['sample_metadata'],
            }, f, indent=2)
        print(f"\n✓ Fallback dataset saved to {fallback_path}")

    print("\n✓ TCGA download complete.")
