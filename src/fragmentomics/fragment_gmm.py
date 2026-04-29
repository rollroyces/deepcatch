"""
Gaussian Mixture Model for cfDNA Fragment Length Distribution (FragmentoSign)

Models cfDNA fragment size distribution as a GMM with 4 biologically meaningful
components reflecting nucleosomal protection patterns.

Components (based on Snyder et al. 2016, Cell 164:57-68):
1. Sub-nucleosomal peak: ~60-100bp  (degraded fragments, ↑ in cancer)
2. Mono-nucleosomal peak: ~160-180bp (single nucleosome protection)
3. Di-nucleosomal peak: ~320-360bp (two nucleosomes)
4. Tri-nucleosomal peak: ~480-540bp (three nucleosomes)

Cancer-specific signal: decreased mono-nucleosomal peak, increased sub-nucleosomal
and shorter fragments (<150bp).

DEPENDENCIES: numpy, scipy, sklearn
"""

import numpy as np
from typing import Dict, Tuple, List, Optional
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
import warnings


class FragmentLengthGMM:
    """
    Gaussian Mixture Model for cfDNA fragment length distributions.
    
    Models the fragment length distribution as a weighted sum of
    Gaussian components, each representing a nucleosomal protection pattern.
    
    Cancer classification based on:
    - Sub-nucleosomal fraction (higher in cancer)
    - Peak position shifts (altered nucleosome spacing)
    - Peak width changes (increased heterogeneity)
    """
    
    # Prior means for nucleosomal peaks (from Snyder et al. 2016)
    PEAK_PRIORS = {
        'sub_nucleosomal': {'mu': 80.0, 'sigma': 15.0, 'label': 'Sub-nucleosomal'},
        'mono_nucleosomal': {'mu': 167.0, 'sigma': 12.0, 'label': 'Mono-nucleosomal'},
        'di_nucleosomal': {'mu': 334.0, 'sigma': 15.0, 'label': 'Di-nucleosomal'},
        'tri_nucleosomal': {'mu': 501.0, 'sigma': 20.0, 'label': 'Tri-nucleosomal'},
    }
    
    def __init__(
        self,
        n_components: int = 4,
        use_priors: bool = True,
        random_state: int = 42
    ):
        """
        Args:
            n_components: Number of Gaussian components (default 4)
            use_priors: Initialize means from known nucleosomal peaks
            random_state: Random seed for reproducibility
        """
        self.n_components = n_components
        self.use_priors = use_priors
        self.random_state = random_state
        self.gmm = None
        self.component_labels = None
        self._fitted = False
    
    def fit(self, fragment_lengths: np.ndarray) -> 'FragmentLengthGMM':
        """
        Fit GMM to fragment length distribution.
        
        Args:
            fragment_lengths: Array of fragment lengths (bp)
        
        Returns:
            self (fitted model)
        """
        X = fragment_lengths.reshape(-1, 1)
        
        # Initialize GMM
        init_params = 'kmeans'
        means_init = None
        
        if self.use_priors and self.n_components == 4:
            # Initialize means from known nucleosomal peaks
            means_init = np.array([
                [self.PEAK_PRIORS['sub_nucleosomal']['mu']],
                [self.PEAK_PRIORS['mono_nucleosomal']['mu']],
                [self.PEAK_PRIORS['di_nucleosomal']['mu']],
                [self.PEAK_PRIORS['tri_nucleosomal']['mu']],
            ])
            init_params = 'k-means++'
        
        self.gmm = GaussianMixture(
            n_components=self.n_components,
            means_init=means_init,
            covariance_type='full',
            random_state=self.random_state,
            max_iter=200,
            n_init=5,
            init_params=init_params,
        )
        
        try:
            self.gmm.fit(X)
            self._fitted = True
            
            # Label components by mean position
            means = self.gmm.means_.flatten()
            peak_order = ['sub_nucleosomal', 'mono_nucleosomal', 
                         'di_nucleosomal', 'tri_nucleosomal']
            component_indices = np.argsort(means)
            self.component_labels = {
                peak_order[i]: component_indices[i]
                for i in range(self.n_components)
            }
            
        except Exception as e:
            warnings.warn(f"GMM fitting failed: {e}. Using prior-only estimates.")
            self._fitted = False
        
        return self
    
    def get_component_stats(self) -> Dict:
        """
        Extract per-component statistics.
        
        Returns:
            Dict with per-peak: weight, mean (bp), std (bp), and 
            derived features for cancer classification
        """
        if not self._fitted or self.gmm is None:
            return self._get_prior_stats()
        
        weights = self.gmm.weights_
        means = self.gmm.means_.flatten()
        stds = np.sqrt(self.gmm.covariances_.flatten())
        
        # Sort by mean
        order = np.argsort(means)
        weights = weights[order]
        means = means[order]
        stds = stds[order]
        
        # Derived features for cancer detection
        stats = {}
        for i, label in enumerate(['sub_nucleosomal', 'mono_nucleosomal',
                                    'di_nucleosomal', 'tri_nucleosomal']):
            if i < len(weights):
                stats[label] = {
                    'weight': float(weights[i]),
                    'mean_bp': float(means[i]),
                    'std_bp': float(stds[i]),
                }
        
        # Cancer-relevant features
        stats['sub_nucleosomal_fraction'] = float(weights[0]) if len(weights) > 0 else 0.0
        stats['mono_nucleosomal_fraction'] = float(weights[1]) if len(weights) > 1 else 0.0
        stats['nucleosome_ratio'] = (
            stats['mono_nucleosomal_fraction'] / 
            (stats['sub_nucleosomal_fraction'] + 1e-6)
        )
        stats['peak_periodicity'] = (
            float(means[2] - means[1]) if len(weights) > 2 
            else float(self.PEAK_PRIORS['di_nucleosomal']['mu'] - self.PEAK_PRIORS['mono_nucleosomal']['mu'])
        )
        
        return stats
    
    def _get_prior_stats(self) -> Dict:
        """Fallback when GMM fitting fails: use literature priors."""
        return {
            'sub_nucleosomal': {'weight': 0.15, 'mean_bp': 80.0, 'std_bp': 15.0},
            'mono_nucleosomal': {'weight': 0.60, 'mean_bp': 167.0, 'std_bp': 12.0},
            'di_nucleosomal': {'weight': 0.20, 'mean_bp': 334.0, 'std_bp': 15.0},
            'tri_nucleosomal': {'weight': 0.05, 'mean_bp': 501.0, 'std_bp': 20.0},
            'sub_nucleosomal_fraction': 0.15,
            'mono_nucleosomal_fraction': 0.60,
            'nucleosome_ratio': 4.0,
            'peak_periodicity': 167.0,
        }
    
    def compute_delfi_score(self, fragment_lengths: np.ndarray) -> float:
        """
        Compute DELFI-style fragmentomics score.
        
        Based on Cristiano et al. 2019: genome-wide fragmentation profiles
        differ between cancer and healthy cfDNA. Key features:
        - Ratio of short (<150bp) to long (>250bp) fragments
        - GMM sub-nucleosomal fraction
        
        Higher score → more cancer-like fragment profile.
        """
        n_short = np.sum(fragment_lengths < 150)
        n_long = np.sum(fragment_lengths > 250)
        n_total = len(fragment_lengths)
        
        if n_total == 0:
            return 0.0
        
        short_long_ratio = n_short / (n_long + 1)
        short_fraction = n_short / n_total
        
        # Combine with GMM if fitted
        if self._fitted:
            stats = self.get_component_stats()
            sub_nuc = stats['sub_nucleosomal_fraction']
            return 0.5 * short_fraction + 0.5 * sub_nuc
        
        return short_fraction


def compute_fragmentomics_features(
    fragment_lengths: np.ndarray,
    fit_gmm: bool = True
) -> Dict:
    """
    Complete fragmentomics feature extraction pipeline (FragmentoSign).
    
    Combines GMM decomposition + DELFI features + MDS into a single
    feature vector for downstream classification.
    
    Returns:
        dict with fragmentomics features for the DeepCatch fusion model
    """
    features = {}
    
    # Basic statistics
    features['fragment_length_mean'] = float(np.mean(fragment_lengths))
    features['fragment_length_median'] = float(np.median(fragment_lengths))
    features['fragment_length_std'] = float(np.std(fragment_lengths))
    features['fragment_length_skew'] = float(
        (np.mean(fragment_lengths) - np.median(fragment_lengths)) / 
        (np.std(fragment_lengths) + 1e-6)
    )
    
    # Short/long fragment ratios (DELFI)
    features['short_fraction_150bp'] = float(np.mean(fragment_lengths < 150))
    features['long_fraction_250bp'] = float(np.mean(fragment_lengths > 250))
    features['short_long_ratio'] = features['short_fraction_150bp'] / (features['long_fraction_250bp'] + 1e-6)
    
    # GMM features
    if fit_gmm:
        gmm = FragmentLengthGMM(n_components=4, use_priors=True, random_state=42)
        gmm.fit(fragment_lengths)
        gmm_stats = gmm.get_component_stats()
        for key in ['sub_nucleosomal_fraction', 'mono_nucleosomal_fraction',
                     'nucleosome_ratio', 'peak_periodicity']:
            features[key] = gmm_stats.get(key, 0.0)
        features['delfi_score'] = gmm.compute_delfi_score(fragment_lengths)
    
    return features
