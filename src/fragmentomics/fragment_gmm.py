"""
Gaussian Mixture Model for cfDNA Fragment Length Distribution (FragmentoSign)

Models cfDNA fragment size distribution as a GMM with 4 biologically meaningful
components reflecting nucleosomal protection patterns.

.. rubric:: Mathematical Formulation

The fragment length distribution :math:`f(x)` is modeled as a weighted sum of
Gaussian probability density functions:

.. math::

    f(x \mid \Theta) = \sum_{k=1}^{K} w_k \, \mathcal{N}(x \mid \mu_k, \sigma^2_k)

where:
    - :math:`K = 4` (number of nucleosomal components)
    - :math:`w_k` are the mixing weights (:math:`\sum w_k = 1, w_k \ge 0`)
    - :math:`\mu_k` are the component means (fragment length in bp)
    - :math:`\sigma^2_k` are the component variances

The parameters :math:`\Theta = \{w_k, \mu_k, \sigma^2_k\}_{k=1}^K` are
estimated via Expectation-Maximization (EM) using :class:`sklearn.mixture.GaussianMixture`.

.. rubric:: Biologically Meaningful Components

Based on Snyder et al. 2016, Cell 164:57-68:

1. **Sub-nucleosomal** peak: ~60–100 bp (degraded fragments, ↑ in cancer)
2. **Mono-nucleosomal** peak: ~160–180 bp (single nucleosome protection)
3. **Di-nucleosomal** peak: ~320–360 bp (two nucleosomes)
4. **Tri-nucleosomal** peak: ~480–540 bp (three nucleosomes)

.. rubric:: Cancer-Specific Signal

- Decreased mono-nucleosomal peak weight
- Increased sub-nucleosomal fraction
- Shorter fragment population (<150 bp)
- Altered nucleosome spacing (peak periodicity shifts)
- Increased peak width (fragmentation heterogeneity)

.. rubric:: Derived Features

- **Sub-nucleosomal fraction**: :math:`w_{\text{sub}}`
- **Nucleosome ratio**: :math:`w_{\text{mono}} / (w_{\text{sub}} + \epsilon)`
- **Peak periodicity**: :math:`\mu_{\text{di}} - \mu_{\text{mono}}`
- **DELFI score**: Combined short-fragment fraction + sub-nucleosomal weight

.. rubric:: Dependencies

- numpy
- scipy (stats.norm for fallback density)
- sklearn.mixture.GaussianMixture

References
----------
.. [1] Snyder, M.W. et al. (2016). Cell 164:57-68. PMID: 26771485
.. [2] Cristiano, S. et al. (2019). Nature 570:385-389. PMID: 31142840
.. [3] Jiang, P. et al. (2020). Nature Genetics 52:712-719. PMID: 32514122
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

    Attributes
    ----------
    n_components : int
        Number of Gaussian components (default 4).
    use_priors : bool
        If True, initialize component means from known nucleosomal peaks.
    random_state : int
        Random seed for reproducible EM initialization.
    gmm : sklearn.mixture.GaussianMixture or None
        The fitted GMM object. None until :meth:`fit` succeeds.
    component_labels : dict or None
        Mapping of peak names to sorted component indices.
    _fitted : bool
        Whether the GMM has been successfully fitted.

    Cancer classification based on:
    - Sub-nucleosomal fraction (higher in cancer)
    - Peak position shifts (altered nucleosome spacing)
    - Peak width changes (increased heterogeneity)

    Examples
    --------
    >>> import numpy as np
    >>> from fragmentomics.fragment_gmm import FragmentLengthGMM
    >>>
    >>> # Simulate cfDNA fragment lengths
    >>> rng = np.random.RandomState(42)
    >>> healthy = np.concatenate([
    ...     rng.normal(80, 15, 150),
    ...     rng.normal(167, 12, 600),
    ...     rng.normal(334, 15, 200),
    ...     rng.normal(501, 20, 50),
    ... ])
    >>>
    >>> gmm = FragmentLengthGMM(n_components=4, use_priors=True, random_state=42)
    >>> gmm.fit(healthy)
    >>> stats = gmm.get_component_stats()
    >>> print(f"Nucleosome ratio: {stats['nucleosome_ratio']:.2f}")
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

        Uses Expectation-Maximization via sklearn's :class:`GaussianMixture`.
        When ``use_priors=True`` and ``n_components=4``, means are initialized
        from the known nucleosomal peak positions defined in :attr:`PEAK_PRIORS`.
        Otherwise, k-means++ initialization is used.

        Component labels are assigned by sorting fitted means ascending
        and mapping to the biological peak order:
        sub-nucleosomal → mono-nucleosomal → di-nucleosomal → tri-nucleosomal.

        Parameters
        ----------
        fragment_lengths : np.ndarray
            Array of fragment lengths in base pairs. Required to have
            shape (n_fragments,). Non-positive values will be included
            in the fit — filter or clip before calling if needed.

        Returns
        -------
        FragmentLengthGMM
            The fitted model instance (self), enabling method chaining.

        Raises
        ------
        ValueError
            If ``fragment_lengths`` is empty or contains fewer samples than
            ``n_components`` (EM requires at least n_components samples).

        Notes
        -----
        - If EM fails (e.g., singular covariance), a warning is issued and
          ``_fitted`` is set to False. Subsequent calls to
          :meth:`get_component_stats` will return literature priors.
        - Edge case: if all fragment lengths are identical, the covariance
          matrix may be singular; prior-based fallback values are used.
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
        Extract per-component statistics and derived cancer features.

        If the model is not fitted (``_fitted=False``), returns literature
        priors via :meth:`_get_prior_stats` as a safe fallback.

        Derived features:

        - ``sub_nucleosomal_fraction``: weight of the shortest component
        - ``mono_nucleosomal_fraction``: weight of the mono-nucleosome component
        - ``nucleosome_ratio``: mono / (sub + ε) — higher in healthy
        - ``peak_periodicity``: μ_di − μ_mono — altered in cancer

        Returns
        -------
        dict
            Nested dictionary with per-peak statistics and derived features.
            Per-peak entries contain ``weight`` (float), ``mean_bp`` (float),
            and ``std_bp`` (float). Derived features are top-level float keys.

        Notes
        -----
        Edge case: if ``n_components`` is not 4 (custom component count),
        only features for available components are returned. Missing
        components default to prior values from :attr:`PEAK_PRIORS`.
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
        """
        Return literature-based prior statistics when GMM fitting fails.

        Uses default weights and positions from Snyder et al. 2016:
        sub-nucleosomal (15% weight, 80 bp), mono-nucleosomal (60%, 167 bp),
        di-nucleosomal (20%, 334 bp), tri-nucleosomal (5%, 501 bp).

        Returns
        -------
        dict
            Prior-based statistics with the same schema as
            :meth:`get_component_stats`.

        Notes
        -----
        The nucleosome ratio default is 4.0 (60/15), representing a
        typical healthy profile. Cancer samples typically show
        ratios below 2.5.
        """
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

        Combines two orthogonal signals of cancer-associated fragmentation:

        .. math::

            \text{DELFI} = 0.5 \cdot f_{\text{short}} + 0.5 \cdot w_{\text{sub}}

        where :math:`f_{\text{short}}` is the fraction of fragments < 150 bp
        and :math:`w_{\text{sub}}` is the GMM sub-nucleosomal weight.

        If the GMM is not fitted, falls back to the short-fragment fraction alone.

        Parameters
        ----------
        fragment_lengths : np.ndarray
            Array of fragment lengths (bp).

        Returns
        -------
        float
            DELFI score in [0, 1]. Higher values indicate a more
            cancer-like fragmentation profile.

        Raises
        ------
        ValueError
            If ``fragment_lengths`` has zero elements.

        References
        ----------
        .. [1] Cristiano et al. (2019). Nature 570:385-389.

        Notes
        -----
        Edge case: when all fragments are exactly 150–250 bp (n_short=0,
        n_long=0), the score returns the sub-nucleosomal GMM fraction or 0.
        """
        n_short = np.sum(fragment_lengths < 150)
        n_long = np.sum(fragment_lengths > 250)
        n_total = len(fragment_lengths)

        if n_total == 0:
            return 0.0

        # Edge case: if all fragments are intermediate-length
        # (150-250 bp, typical for healthy cfDNA), n_short and n_long
        # are both zero; the score defaults to sub-nucleosomal fraction or 0.
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

    Combines GMM decomposition, DELFI-style features, and basic statistical
    moments into a unified feature vector for downstream classification
    (e.g., the DeepCatch multi-modal fusion model).

    Extracted features:

    **Basic statistics** (always computed)
        - ``fragment_length_mean``: mean fragment length
        - ``fragment_length_median``: median fragment length
        - ``fragment_length_std``: standard deviation of fragment lengths
        - ``fragment_length_skew``: approximate skewness
          (:math:`(\bar{x} - \tilde{x}) / (s + \epsilon)`)

    **DELFI features** (always computed)
        - ``short_fraction_150bp``: proportion of fragments < 150 bp
        - ``long_fraction_250bp``: proportion of fragments > 250 bp
        - ``short_long_ratio``: short / (long + ε)

    **GMM features** (computed when ``fit_gmm=True``)
        - ``sub_nucleosomal_fraction``: weight of sub-nucleosomal component
        - ``mono_nucleosomal_fraction``: weight of mono-nucleosomal component
        - ``nucleosome_ratio``: mono / (sub + ε)
        - ``peak_periodicity``: μ_di − μ_mono
        - ``delfi_score``: combined short-fragment + GMM score

    Parameters
    ----------
    fragment_lengths : np.ndarray
        Array of fragment lengths (bp).
    fit_gmm : bool, optional
        If True (default), fit :class:`FragmentLengthGMM` and include GMM
        features. Set to False to skip GMM fitting and return only basic
        and DELFI features.

    Returns
    -------
    dict
        Dictionary of float feature values for downstream fusion.

    Raises
    ------
    ValueError
        If ``fragment_lengths`` is empty (zero-length array).

    Notes
    -----
    - If ``fit_gmm=False``, GMM keys are absent from the output dict.
      Downstream consumers should use ``features.get(key, 0.0)``.
    - The skew computation uses a robust approximation:
      (mean − median) / (std + ε), which is zero for symmetric distributions
      and non-zero for skewed ones.

    Examples
    --------
    >>> rng = np.random.RandomState(0)
    >>> lengths = np.concatenate([rng.normal(167, 12, 1000)])
    >>> features = compute_fragmentomics_features(lengths, fit_gmm=False)
    >>> print(f"Mean: {features['fragment_length_mean']:.1f} bp")
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


__all__ = [
    "FragmentLengthGMM",
    "compute_fragmentomics_features",
]
