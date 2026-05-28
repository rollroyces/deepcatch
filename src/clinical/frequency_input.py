#!/usr/bin/env python3
"""
Frequency Vector Input Module — P0-B1
======================================

Accepts pre-computed motif frequency vectors (e.g., Jiang 4-mer data:
samples × 256 motifs) and feeds them into the CET analysis pipeline.

Designed to bridge the gap between simulation-generated data and
real-world frequency matrices from external tools (BAM extractors,
k-mer counters, published datasets).

Classes
-------
FrequencyDataset
    Load, validate, and describe frequency data from Excel/CSV/NumPy.
PlotGenerator
    Publication-quality visualisation: volcano, heatmap, ROC, feature importance.

References
----------
Jiang et al. (2018) Nat Commun 9:XXXX.  Plasma DNA end-motif profiling
for cancer detection.  256 × 4-mer motif frequencies.
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy.stats import rankdata

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# FrequencyDataset
# ═══════════════════════════════════════════════════════════════════════════

class FrequencyDataset:
    """Load and validate pre-computed motif frequency matrices.

    Supports Excel (.xlsx/.xls), CSV, and NumPy (.npy/.npz) formats.
    Auto-detects whether the first column contains sample IDs or
    feature names, and handles label assignment from the same file
    or a separate label file.

    Parameters
    ----------
    data_path : str or Path
        Path to the frequency data file.
    label_path : str, Path, or None
        Path to a separate label file.  If None, labels are expected
        inside `data_path` (column named 'label' or last column).
    format : str
        One of 'excel', 'csv', 'npy', 'npz', or 'auto' (infer from
        extension).

    Examples
    --------
    >>> ds = FrequencyDataset('jiang_4mer.xlsx')
    >>> X, y, names = ds.load()
    >>> ds.validate()
    >>> ds.describe()
    """

    _SUPPORTED_EXTENSIONS = {
        '.xlsx': 'excel', '.xls': 'excel',
        '.csv':  'csv',   '.tsv':  'csv',
        '.npy':  'npy',   '.npz':  'npz',
    }

    def __init__(
        self,
        data_path: Union[str, Path],
        label_path: Optional[Union[str, Path]] = None,
        format: str = 'auto',
    ):
        self.data_path = Path(data_path)
        self.label_path = Path(label_path) if label_path else None
        self._fmt = format if format != 'auto' else self._infer_format()
        self._validate_path()

        # Populated by load()
        self.X: Optional[np.ndarray] = None
        self.y: Optional[np.ndarray] = None
        self.feature_names: Optional[List[str]] = None
        self.sample_ids: Optional[List[str]] = None
        self._loaded = False

    # ── public API ────────────────────────────────────────────────────────

    def load(self) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """Load data and return (X, y, feature_names).

        Returns
        -------
        X : np.ndarray  shape (n_samples, n_features)
        y : np.ndarray  shape (n_samples,) — binary 0/1 labels
        feature_names : list of str
        """
        if self._loaded:
            return self.X, self.y, self.feature_names

        raw = self._read_file()

        # ── Separate labels if embedded ────────────────────────────────
        if self.label_path is not None:
            label_raw = self._read_label_file()
            y, self._label_mapping = self._parse_labels_multiclass(label_raw)
            self._raw_labels_str = label_raw.values
            X = raw.values.astype(np.float64)
            feature_names = list(raw.columns)
            sample_ids = list(raw.index.astype(str))
        else:
            # Try to find label column in the data
            raw_df = raw.copy()  # snapshot before _split_embedded_labels transforms it
            X, y, feature_names, sample_ids = self._split_embedded_labels(raw)
            # Store original string labels (extract from the label column before it was dropped)
            if 'Group' in raw_df.columns:
                self._raw_labels_str = raw_df['Group'].astype(str).values
                _, self._label_mapping = self._parse_labels_multiclass(raw_df['Group'])
            elif 'group' in raw_df.columns:
                self._raw_labels_str = raw_df['group'].astype(str).values
                _, self._label_mapping = self._parse_labels_multiclass(raw_df['group'])
            elif 'label' in raw_df.columns:
                self._raw_labels_str = raw_df['label'].astype(str).values
                _, self._label_mapping = self._parse_labels_multiclass(raw_df['label'])

        self.X = X
        self.y = y
        self.feature_names = feature_names
        self.sample_ids = sample_ids
        self._raw_X = X.copy()  # preserve for ratio computation
        self._raw_feature_names = list(feature_names)
        self._loaded = True

        logger.info(
            "FrequencyDataset loaded: %d samples × %d features, "
            "classes: %s",
            X.shape[0], X.shape[1],
            dict(zip(*np.unique(y, return_counts=True))),
        )
        return X, y, feature_names

    def validate(self) -> Dict[str, Union[bool, float, int]]:
        """Run quality checks and return a report dict.

        Checks performed:
        - NaN/infinite values
        - Constant (zero-variance) features
        - Class balance (warns if minority < 10%)
        - Negative frequency values
        - Row sums far from 1.0 (for compositional data)
        """
        if not self._loaded:
            self.load()

        report: Dict[str, Union[bool, float, int]] = {
            'passed': True,
            'n_samples': self.X.shape[0],
            'n_features': self.X.shape[1],
        }

        # 1. NaN / infinity
        nan_count = int(np.sum(np.isnan(self.X)))
        inf_count = int(np.sum(np.isinf(self.X)))
        report['n_nan'] = nan_count
        report['n_inf'] = inf_count
        if nan_count > 0 or inf_count > 0:
            report['passed'] = False
            logger.warning("Found %d NaN and %d inf values", nan_count, inf_count)

        # 2. Constant features
        stds = np.std(self.X, axis=0)
        n_constant = int(np.sum(stds == 0))
        report['n_constant_features'] = n_constant
        if n_constant > 0:
            report['passed'] = False
            logger.warning("Found %d constant (zero-variance) features", n_constant)

        # 3. Class balance
        n_unique = len(np.unique(self.y))
        report['n_classes'] = n_unique
        n_pos = int(np.sum(self.y == 1))
        n_neg = int(np.sum(self.y == 0))
        total = len(self.y)
        min_pct = min(n_pos, n_neg) / max(1, total)
        report['n_positive'] = n_pos
        report['n_negative'] = n_neg
        report['minority_pct'] = round(min_pct * 100, 1)
        if min_pct < 0.10:
            logger.warning(
                "Severe class imbalance: minority class = %.1f%% (%d samples). "
                "Consider stratified sampling or SMOTE.",
                min_pct * 100, min(n_pos, n_neg),
            )
        elif min_pct < 0.20:
            logger.info("Moderate class imbalance: minority = %.1f%%", min_pct * 100)

        # 4. Negative values (frequency data should be ≥ 0)
        n_negative = int(np.sum(self.X < 0))
        report['n_negative_values'] = n_negative
        if n_negative > 0:
            logger.warning("Found %d negative values — are these really frequencies?", n_negative)

        # 5. Row-sum check (motif frequencies should sum to ~1)
        if self.X.shape[1] >= 10:
            row_sums = np.sum(self.X, axis=1)
            far_from_one = np.sum(np.abs(row_sums - 1.0) > 0.1)
            report['rows_far_from_unity'] = int(far_from_one)
            if far_from_one > self.X.shape[0] * 0.5:
                logger.info(
                    "Note: %d/%d rows don't sum to ~1.0 — data may not be "
                    "compositional frequencies.",
                    far_from_one, self.X.shape[0],
                )

        if report['passed']:
            logger.info("✅ FrequencyDataset validation PASSED")
        else:
            logger.warning("⚠️  FrequencyDataset validation found issues (see report)")

        return report

    def to_cet_format(self) -> Dict:
        """Convert to a dict consumable by the CET analysis pipeline.

        Returns a dict with keys matching what `cet_cross_validator`
        expects: ``X``, ``y``, ``feature_names``, ``sample_ids``.
        """
        if not self._loaded:
            self.load()

        return {
            'X': self.X,
            'y': self.y,
            'feature_names': self.feature_names,
            'sample_ids': self.sample_ids,
            'n_samples': self.X.shape[0],
            'n_features': self.X.shape[1],
            'class_distribution': {
                'positive': int(np.sum(self.y == 1)),
                'negative': int(np.sum(self.y == 0)),
            },
        }

    def describe(self) -> str:
        """Return a human-readable summary string."""
        if not self._loaded:
            self.load()

        lines = [
            "═" * 50,
            " FrequencyDataset Summary",
            "═" * 50,
            f"  Source:       {self.data_path}",
            f"  Samples:      {self.X.shape[0]}",
            f"  Features:     {self.X.shape[1]}",
            f"  Feature range: [{self.X.min():.6f}, {self.X.max():.6f}]",
            f"  Mean:         {self.X.mean():.6f}",
            f"  Std:          {self.X.std():.6f}",
        ]

        if self.y is not None:
            n_pos = int(np.sum(self.y == 1))
            n_neg = int(np.sum(self.y == 0))
            lines.extend([
                f"  Class +:      {n_pos} ({n_pos / len(self.y) * 100:.1f}%)",
                f"  Class −:      {n_neg} ({n_neg / len(self.y) * 100:.1f}%)",
            ])

        n_nan = int(np.sum(np.isnan(self.X)))
        lines.append(f"  NaN values:   {n_nan}")

        if self.feature_names:
            lines.append(f"  First 5 features: {self.feature_names[:5]}")
            lines.append(f"  Last 5 features:  {self.feature_names[-5:]}")

        lines.append("═" * 50)
        return "\n".join(lines)

    def add_rank_features(self) -> None:
        """Convert motif frequencies to ranks per sample.

        For each sample, replaces raw motif frequencies with their
        rank order (1 = highest frequency, 256 = lowest). This
        eliminates amplitude-based batch effects while preserving
        the relative motif abundance pattern — which carries the
        cancer signal.

        Rank transformation has been shown to improve cross-platform
        reproducibility in cfDNA analysis and provides robust
        batch-effect normalisation without external reference.
        """
        if not self._loaded or self.X is None:
            self.load()

        original_shape = self.X.shape
        self.X = np.array([
            rankdata(-self.X[i, :]) for i in range(self.X.shape[0])
        ])
        logger.info(
            "Applied rank transformation: %s (preserves %s shape)",
            list(self.X.shape), list(original_shape),
        )

    def add_composition_ratios(self) -> None:

        feature_names = self._raw_feature_names if hasattr(self, '_raw_feature_names') else self.feature_names
        X_raw = self._raw_X if hasattr(self, '_raw_X') else self.X

        # Identify pure CG-rich (C+G, no A,T) and pure AT-rich (A/T, no C,G) motifs
        cg_indices = [
            i for i, m in enumerate(feature_names)
            if 'C' in m and 'G' in m and 'A' not in m and 'T' not in m
        ]
        at_indices = [
            i for i, m in enumerate(feature_names)
            if ('A' in m or 'T' in m) and 'C' not in m and 'G' not in m
        ]

        if not cg_indices or not at_indices:
            logger.warning("Could not identify CG/AT motifs — skipping ratio features")
            return

        # Compute from raw X (before filtering) to preserve global statistics
        cg_total = X_raw[:, cg_indices].sum(axis=1)
        at_total = X_raw[:, at_indices].sum(axis=1)
        ratio = np.where(at_total > 1e-12, cg_total / at_total, 0.0)

        # If already filtered, subset
        if self.X.shape[0] < X_raw.shape[0]:
            # self.X already filtered — need to re-index
            if hasattr(self, '_filter_mask') and self._filter_mask is not None:
                cg_total = cg_total[self._filter_mask]
                at_total = at_total[self._filter_mask]
                ratio = ratio[self._filter_mask]

        ratio_features = np.column_stack([cg_total, at_total, ratio])
        ratio_names = ['cg_total', 'at_total', 'cg_at_ratio']

        self.X = np.column_stack([self.X, ratio_features])
        self.feature_names = list(self.feature_names) + ratio_names

        logger.info(
            "Added %d ratio features: %s (X now %s)",
            len(ratio_names), ratio_names, list(self.X.shape),
        )

    def filter_by_label(self, target_label: str,
                         control_label: Optional[str] = None) -> None:
        """Filter dataset to binary comparison (target vs control)."""
        if not self._loaded:
            self.load()

        # Use raw string labels for matching
        if hasattr(self, '_raw_labels_str') and self._raw_labels_str is not None:
            raw_labels = self._raw_labels_str
        else:
            raw_labels = np.array([str(v) for v in self.y])

        # Create boolean mask
        target_str = str(target_label).strip()
        mask_pos = np.array([str(v).strip() == target_str for v in raw_labels])

        if control_label is not None:
            control_str = str(control_label).strip()
            mask_neg = np.array([str(v).strip() == control_str for v in raw_labels])
        else:
            mask_neg = ~mask_pos

        mask = mask_pos | mask_neg
        y_new = np.where(mask_pos[mask], 1, 0)

        # Store filter mask for ratio computation
        self._filter_mask = mask

        self.X = self.X[mask]
        self.y = y_new
        if hasattr(self, 'sample_ids') and self.sample_ids:
            self.sample_ids = [s for i, s in enumerate(self.sample_ids) if mask[i]]
        n_pos = int(y_new.sum())
        n_neg = len(y_new) - n_pos
        logger.info("filter_by_label(%s): %d samples (%d pos, %d neg)",
                     target_label, len(self.y), n_pos, n_neg)

    # ── internal helpers ──────────────────────────────────────────────────

    def _infer_format(self) -> str:
        ext = self.data_path.suffix.lower()
        fmt = self._SUPPORTED_EXTENSIONS.get(ext)
        if fmt is None:
            raise ValueError(
                f"Unsupported file extension '{ext}'. "
                f"Supported: {list(self._SUPPORTED_EXTENSIONS)}"
            )
        return fmt

    def _validate_path(self) -> None:
        if not self.data_path.exists():
            raise FileNotFoundError(f"Data file not found: {self.data_path}")
        if self.label_path is not None and not self.label_path.exists():
            raise FileNotFoundError(f"Label file not found: {self.label_path}")

    def _infer_and_read_excel(self) -> pd.DataFrame:
        """Read Excel file with smart header detection.

        Handles common layouts:
        1. Standard: header=0, data starts at row 1
        2. Title row: row 0 is a title, header at row 1 (prof_jiang format)
        3. Column names in row 0, data row 1

        Returns DataFrame with proper column names and no title rows.
        """
        # Quick peek to detect multi-row header
        peek = pd.read_excel(self.data_path, nrows=3, header=None)
        first_val = peek.iloc[0, 0]

        # If first cell of first row is a long string (title), skip it
        if isinstance(first_val, str) and len(str(first_val)) > 20:
            logger.info("Detected title row — using header=1")
            return pd.read_excel(self.data_path, header=1)

        # Standard case
        return pd.read_excel(self.data_path, index_col=None)

    def _read_file(self) -> pd.DataFrame:
        if self._fmt == 'excel':
            return self._infer_and_read_excel()
        elif self._fmt == 'csv':
            return pd.read_csv(self.data_path, index_col=None)
        elif self._fmt == 'npy':
            arr = np.load(self.data_path)
            return pd.DataFrame(arr)
        elif self._fmt == 'npz':
            data = np.load(self.data_path)
            # Try common keys
            for key in ['data', 'X', 'arr_0', 'features']:
                if key in data:
                    return pd.DataFrame(data[key])
            raise KeyError(
                f"Unknown keys in .npz: {list(data.keys())}. "
                f"Expected one of: data, X, arr_0, features"
            )
        else:
            raise ValueError(f"Unknown format: {self._fmt}")

    def _read_label_file(self) -> pd.Series:
        label_ext = self.label_path.suffix.lower()
        if label_ext in ('.xlsx', '.xls'):
            return pd.read_excel(self.label_path, index_col=None).iloc[:, 0]
        elif label_ext == '.csv':
            return pd.read_csv(self.label_path, index_col=None).iloc[:, 0]
        elif label_ext in ('.npy', '.npz'):
            if label_ext == '.npy':
                arr = np.load(self.label_path)
            else:
                data = np.load(self.label_path)
                arr = data[list(data.keys())[0]]
            return pd.Series(arr.flatten())
        else:
            raise ValueError(f"Unsupported label format: {label_ext}")

    @staticmethod
    def _parse_labels(series: pd.Series) -> np.ndarray:
        """Convert a label series to binary 0/1 (kept for backward compat)."""
        vals = series.dropna()
        if vals.dtype == bool:
            return vals.astype(int).values
        if vals.dtype in (int, np.int64, np.int32):
            unique = set(vals.unique())
            if unique <= {0, 1}:
                return vals.values.astype(int)
            logger.warning("Labels have %d unique values — expected binary", len(unique))
            return vals.values.astype(int)

        unique_vals = vals.unique()
        if len(unique_vals) == 2:
            mapping = {str(v): i for i, v in enumerate(sorted(unique_vals))}
            logger.info("Label mapping: %s", mapping)
            return np.array([mapping[str(v)] for v in vals])

        logger.warning(
            "Non-binary string labels with %d categories. "
            "Mapping alphabetically.", len(unique_vals)
        )
        mapping = {str(v): i for i, v in enumerate(sorted(unique_vals))}
        return np.array([mapping[str(v)] for v in vals])

    @staticmethod
    def _parse_labels_multiclass(series: pd.Series) -> Tuple[np.ndarray, dict]:
        """Convert a label series to numeric labels, returning mapping."""
        vals = series.dropna()
        unique_vals = sorted(vals.unique(), key=str)
        mapping = {str(v): i for i, v in enumerate(unique_vals)}
        inv_mapping = {i: str(v) for v, i in mapping.items()}
        logger.info("Label mapping (%d classes): %s", len(mapping), dict(zip(mapping.values(), mapping.keys())))
        labels = np.array([mapping[str(v)] for v in vals])
        return labels, inv_mapping

    @staticmethod
    def _split_embedded_labels(
        df: pd.DataFrame,
    ) -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
        """Try to detect label column and sample IDs within the DataFrame.

        Heuristics (in order):
        1. Column named 'label', 'class', 'target', 'y' → labels
        2. Last column is binary/categorical with few unique values → labels
        3. First column is string/object → sample IDs
        4. First column is numeric but looks like IDs → sample IDs
        """
        # Detect label column
        label_col = None
        for col_name in ['label', 'class', 'target', 'y', 'group', 'Label', 'Class', 'Target', 'Group', 'diagnosis']:
            if col_name in df.columns:
                label_col = col_name
                break

        if label_col is None:
            # Guess: last column if it has ≤ 10 unique values
            last_col = df.columns[-1]
            if df[last_col].nunique() <= 10:
                label_col = last_col

        # Extract labels
        if label_col is not None:
            y_series = df[label_col]
            df_feat = df.drop(columns=[label_col])
            y, _ = FrequencyDataset._parse_labels_multiclass(y_series)
        else:
            # No labels found — assume all columns are features
            y = np.zeros(len(df), dtype=int)
            df_feat = df.copy()
            logger.warning(
                "No label column detected. All labels set to 0. "
                "Provide --labels to specify labels separately."
            )

        # Detect sample IDs (first column if string/object)
        sample_ids = None
        first_col = df_feat.columns[0]
        if pd.api.types.is_object_dtype(df_feat[first_col]) or pd.api.types.is_string_dtype(df_feat[first_col]):
            sample_ids = df_feat[first_col].astype(str).tolist()
            df_feat = df_feat.drop(columns=[first_col])
        elif pd.api.types.is_object_dtype(df_feat.index) or pd.api.types.is_string_dtype(df_feat.index):
            sample_ids = df_feat.index.astype(str).tolist()

        # Auto-drop any remaining non-numeric columns (e.g. 'group', 'sample_id')
        numeric_cols = df_feat.select_dtypes(include=['number']).columns.tolist()
        if len(numeric_cols) < len(df_feat.columns):
            dropped = [c for c in df_feat.columns if c not in numeric_cols]
            logger.info("Dropped non-numeric columns: %s", dropped)
            df_feat = df_feat[numeric_cols]

        if sample_ids is None:
            sample_ids = [f"sample_{i:04d}" for i in range(len(df_feat))]

        # Ensure numeric features
        X = df_feat.values.astype(np.float64)
        feature_names = list(df_feat.columns)

        return X, y, feature_names, sample_ids


# ═══════════════════════════════════════════════════════════════════════════
# PlotGenerator
# ═══════════════════════════════════════════════════════════════════════════

class PlotGenerator:
    """Publication-quality plots for CET / frequency-vector analysis.

    All methods return the ``matplotlib.figure.Figure`` and
    simultaneously save to disk.

    Parameters
    ----------
    style : str
        Matplotlib style name.  Default ``'seaborn-v0_8-whitegrid'``.
    dpi : int
        Output resolution.
    figsize : tuple
        Default figure size (width, height) in inches.
    """

    def __init__(
        self,
        style: str = 'seaborn-v0_8-whitegrid',
        dpi: int = 150,
        figsize: Tuple[float, float] = (10, 6),
    ):
        import matplotlib
        matplotlib.use('Agg')  # non-interactive backend
        import matplotlib.pyplot as plt

        self.plt = plt
        try:
            plt.style.use(style)
        except Exception:
            plt.style.use('default')
        self.dpi = dpi
        self.figsize = figsize

    # ── volcano plot ──────────────────────────────────────────────────────

    def volcano_plot(
        self,
        p_values: np.ndarray,
        effect_sizes: np.ndarray,
        labels: Optional[List[str]] = None,
        title: str = "Volcano Plot",
        save_path: Optional[Union[str, Path]] = None,
        alpha: float = 0.05,
        effect_threshold: float = 0.2,
        annotate_top: int = 10,
    ):
        """Volcano plot: −log₁₀(p) vs effect size.

        Parameters
        ----------
        p_values : (n_features,)
        effect_sizes : (n_features,)
        labels : list of str or None
        title : str
        save_path : str, Path, or None
        alpha : float — significance threshold (horizontal line)
        effect_threshold : float — practical significance (vertical lines)
        annotate_top : int — label top N most significant points
        """
        p = np.asarray(p_values, dtype=float)
        es = np.asarray(effect_sizes, dtype=float)
        p_safe = np.maximum(p, 1e-300)  # avoid log(0)
        neg_log_p = -np.log10(p_safe)

        sig = p < alpha
        large_effect = np.abs(es) >= effect_threshold

        # Categories for colouring
        upregulated = sig & (es > 0) & large_effect
        downregulated = sig & (es < 0) & large_effect
        significant_only = sig & ~large_effect
        nonsig = ~sig

        fig, ax = self.plt.subplots(figsize=self.figsize, dpi=self.dpi)

        ax.scatter(es[nonsig], neg_log_p[nonsig], c='grey', alpha=0.4, s=12,
                   label='NS', edgecolors='none')
        ax.scatter(es[significant_only], neg_log_p[significant_only],
                   c='orange', alpha=0.6, s=14, label=f'p < {alpha}',
                   edgecolors='none')
        ax.scatter(es[upregulated], neg_log_p[upregulated],
                   c='#d62728', alpha=0.8, s=18,
                   label=f'Up (|δ| ≥ {effect_threshold})', edgecolors='none')
        ax.scatter(es[downregulated], neg_log_p[downregulated],
                   c='#1f77b4', alpha=0.8, s=18,
                   label=f'Down (|δ| ≥ {effect_threshold})', edgecolors='none')

        # Threshold lines
        ax.axhline(-np.log10(alpha), color='red', linestyle='--', linewidth=0.8,
                   label=f'α = {alpha}')
        ax.axvline(effect_threshold, color='grey', linestyle=':', linewidth=0.6)
        ax.axvline(-effect_threshold, color='grey', linestyle=':', linewidth=0.6)

        # Annotate top hits
        if labels is not None and annotate_top > 0:
            order = np.argsort(p_safe)
            for i in order[:annotate_top]:
                ax.annotate(
                    str(labels[i]),
                    (es[i], neg_log_p[i]),
                    fontsize=7, alpha=0.9,
                    xytext=(5, 5), textcoords='offset points',
                )

        ax.set_xlabel("Effect Size (Cliff's δ)")
        ax.set_ylabel("−log₁₀(p-value)")
        ax.set_title(title)
        ax.legend(loc='best', fontsize=8, framealpha=0.7)

        fig.tight_layout()
        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            logger.info("Volcano plot saved → %s", save_path)

        return fig

    # ── heatmap ───────────────────────────────────────────────────────────

    def heatmap(
        self,
        data: np.ndarray,
        feature_names: Optional[List[str]] = None,
        sample_labels: Optional[List[str]] = None,
        n_top: int = 30,
        save_path: Optional[Union[str, Path]] = None,
        cmap: str = 'RdBu_r',
        cluster: bool = True,
    ):
        """Heatmap of the top *n_top* most variable features.

        Parameters
        ----------
        data : (n_samples, n_features)
        feature_names : list of str or None
        sample_labels : list of str or None — used for colour bar
        n_top : int
        save_path : str, Path, or None
        cmap : str — matplotlib colormap
        cluster : bool — hierarchical clustering on rows/cols
        """
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import pdist

        X = np.asarray(data, dtype=float)

        # Select top N most variable features
        var = np.var(X, axis=0)
        top_idx = np.argsort(-var)[:min(n_top, X.shape[1])]
        X_sub = X[:, top_idx]

        if feature_names is not None:
            top_names = [feature_names[i] for i in top_idx]
        else:
            top_names = [f"F{i}" for i in top_idx]

        # Optional clustering
        if cluster and X_sub.shape[0] > 2:
            row_linkage = linkage(pdist(X_sub), method='average')
            row_order = leaves_list(row_linkage)
            X_sub = X_sub[row_order, :]
        if cluster and X_sub.shape[1] > 2:
            col_linkage = linkage(pdist(X_sub.T), method='average')
            col_order = leaves_list(col_linkage)
            X_sub = X_sub[:, col_order]
            top_names = [top_names[i] for i in col_order]

        # Z-score normalise rows
        row_means = X_sub.mean(axis=1, keepdims=True)
        row_stds = X_sub.std(axis=1, keepdims=True)
        row_stds = np.where(row_stds == 0, 1, row_stds)
        X_norm = (X_sub - row_means) / row_stds

        fig_height = max(6, X_norm.shape[0] * 0.25)
        fig_width = max(10, X_norm.shape[1] * 0.2)
        fig, ax = self.plt.subplots(figsize=(fig_width, fig_height), dpi=self.dpi)

        im = ax.imshow(X_norm, aspect='auto', cmap=cmap, interpolation='nearest')

        ax.set_xticks(range(len(top_names)))
        ax.set_xticklabels(top_names, rotation=90, fontsize=6, ha='center')
        ax.set_yticks(range(X_norm.shape[0]))
        if sample_labels is not None:
            ax.set_yticklabels(
                [sample_labels[i] for i in (row_order if cluster else range(len(sample_labels)))],
                fontsize=5,
            )
        else:
            ax.set_yticklabels([])

        ax.set_title(f"Top {n_top} Most Variable Features")
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('Z-score')

        fig.tight_layout()
        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            logger.info("Heatmap saved → %s", save_path)

        return fig

    # ── ROC curve ─────────────────────────────────────────────────────────

    def roc_curve(
        self,
        y_true: np.ndarray,
        y_scores: np.ndarray,
        title: str = "ROC Curve",
        save_path: Optional[Union[str, Path]] = None,
    ):
        """ROC curve with AUC.

        Parameters
        ----------
        y_true : (n_samples,) binary
        y_scores : (n_samples,) continuous probability scores
        title : str
        save_path : str, Path, or None
        """
        from sklearn.metrics import roc_curve, auc as _auc

        yt = np.asarray(y_true, dtype=int)
        ys = np.asarray(y_scores, dtype=float)

        # Handle edge case: only one class
        unique_cls = np.unique(yt)
        if len(unique_cls) < 2:
            fig, ax = self.plt.subplots(figsize=self.figsize, dpi=self.dpi)
            ax.text(0.5, 0.5, "ROC undefined (single class)", ha='center', va='center',
                    transform=ax.transAxes, fontsize=14)
            ax.set_title(title)
            if save_path:
                Path(save_path).parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            return fig

        fpr, tpr, _ = roc_curve(yt, ys)
        auc_val = _auc(fpr, tpr)

        fig, ax = self.plt.subplots(figsize=self.figsize, dpi=self.dpi)

        ax.plot(fpr, tpr, color='#d62728', linewidth=2,
                label=f'AUC = {auc_val:.3f}')
        ax.plot([0, 1], [0, 1], color='grey', linestyle='--', linewidth=0.8,
                label='Random (0.5)')

        ax.set_xlabel('False Positive Rate (1 − Specificity)')
        ax.set_ylabel('True Positive Rate (Sensitivity)')
        ax.set_title(title)
        ax.legend(loc='lower right')
        ax.set_xlim([-0.02, 1.02])
        ax.set_ylim([-0.02, 1.02])

        fig.tight_layout()
        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            logger.info("ROC curve saved → %s", save_path)

        return fig

    # ── feature importance ────────────────────────────────────────────────

    def feature_importance(
        self,
        coefficients: np.ndarray,
        feature_names: Optional[List[str]] = None,
        n_top: int = 20,
        save_path: Optional[Union[str, Path]] = None,
        title: str = "Feature Importance",
    ):
        """Horizontal bar chart of top feature coefficients.

        Parameters
        ----------
        coefficients : (n_features,) — e.g., logistic regression coefs
        feature_names : list of str or None
        n_top : int — show top N by absolute value
        save_path : str, Path, or None
        title : str
        """
        coef = np.asarray(coefficients, dtype=float).flatten()
        n = len(coef)

        if feature_names is None:
            feature_names = [f"F{i}" for i in range(n)]

        # Top N by absolute value
        abs_coef = np.abs(coef)
        top_idx = np.argsort(-abs_coef)[:min(n_top, n)]

        top_coef = coef[top_idx]
        top_names = [feature_names[i] for i in top_idx]

        # Sort by coefficient value for display
        sort_idx = np.argsort(top_coef)
        top_coef = top_coef[sort_idx]
        top_names = [top_names[i] for i in sort_idx]

        colours = ['#1f77b4' if c < 0 else '#d62728' for c in top_coef]

        fig_height = max(5, n_top * 0.3)
        fig, ax = self.plt.subplots(figsize=(8, fig_height), dpi=self.dpi)

        y_pos = range(len(top_coef))
        ax.barh(y_pos, top_coef, color=colours, edgecolor='white', linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_names, fontsize=8)
        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_xlabel('Coefficient')
        ax.set_title(title)
        ax.invert_yaxis()  # largest at top

        fig.tight_layout()
        if save_path:
            Path(save_path).parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path, dpi=self.dpi, bbox_inches='tight')
            logger.info("Feature importance plot saved → %s", save_path)

        return fig


# ═══════════════════════════════════════════════════════════════════════════
# Demo / self-test
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    print("=" * 60)
    print("FrequencyDataset + PlotGenerator — self-test")
    print("=" * 60)

    # ── Create synthetic frequency data ───────────────────────────────
    rng = np.random.RandomState(42)
    n_samples, n_features = 100, 64
    feature_names = [f"motif_{i:03d}" for i in range(n_features)]

    # Generate realistic frequency vectors (sum ~ 1)
    X = rng.dirichlet(np.ones(n_features) * 2, size=n_samples)
    # Add signal: cancer samples have elevated first 5 motifs
    cancer_idx = rng.choice(n_samples, size=40, replace=False)
    X[cancer_idx, :5] += 0.02
    X = X / X.sum(axis=1, keepdims=True)  # re-normalise

    y = np.zeros(n_samples, dtype=int)
    y[cancer_idx] = 1

    # ── Save to temp files and load ──────────────────────────────────
    import tempfile, os

    with tempfile.TemporaryDirectory() as tmpdir:
        xlsx_path = os.path.join(tmpdir, 'test_frequencies.xlsx')
        df = pd.DataFrame(X, columns=feature_names)
        df.insert(0, 'sample_id', [f"S{i:04d}" for i in range(n_samples)])
        df['label'] = y
        df.to_excel(xlsx_path, index=False)

        # FrequencyDataset
        ds = FrequencyDataset(xlsx_path)
        X_ld, y_ld, names = ds.load()
        print(ds.describe())

        report = ds.validate()
        print(f"\nValidation: {'PASSED' if report['passed'] else 'ISSUES'}")
        for k, v in report.items():
            print(f"  {k}: {v}")

        cet_data = ds.to_cet_format()
        print(f"\nCET format keys: {list(cet_data)}")

        # ── PlotGenerator ────────────────────────────────────────────
        pg = PlotGenerator()

        # Volcano
        from scipy.stats import mannwhitneyu

        p_values = []
        effect_sizes = []
        for j in range(n_features):
            pos_vals = X[cancer_idx, j]
            neg_vals = X[~np.isin(range(n_samples), cancer_idx), j]
            stat, p = mannwhitneyu(pos_vals, neg_vals, alternative='two-sided')
            p_values.append(p)
            # Cliff's delta approximation
            diff = np.mean(pos_vals) - np.mean(neg_vals)
            pooled_std = np.sqrt((np.var(pos_vals) + np.var(neg_vals)) / 2)
            d = diff / max(pooled_std, 1e-12)
            effect_sizes.append(d)

        fig_v = pg.volcano_plot(
            np.array(p_values), np.array(effect_sizes),
            labels=feature_names, title="Motif Volcano Plot",
            save_path=os.path.join(tmpdir, 'volcano.png'),
        )
        print(f"\nVolcano: {fig_v.get_size_inches()}")

        # Heatmap
        sample_lbl = ['Cancer' if yi else 'Healthy' for yi in y]
        fig_h = pg.heatmap(
            X, feature_names=feature_names,
            sample_labels=sample_lbl, n_top=20,
            save_path=os.path.join(tmpdir, 'heatmap.png'),
        )
        print(f"Heatmap: {fig_h.get_size_inches()}")

        # ROC
        scores = X[:, 0]  # dummy scores
        fig_r = pg.roc_curve(
            y, scores, title="Motif-1 ROC",
            save_path=os.path.join(tmpdir, 'roc.png'),
        )
        print(f"ROC: {fig_r.get_size_inches()}")

        # Feature importance
        from sklearn.linear_model import LogisticRegression
        lr = LogisticRegression(max_iter=500).fit(X, y)
        fig_fi = pg.feature_importance(
            lr.coef_.flatten(), feature_names=feature_names,
            n_top=15, title="LR Coefficients",
            save_path=os.path.join(tmpdir, 'feature_importance.png'),
        )
        print(f"Feature importance: {fig_fi.get_size_inches()}")

    print("\n✅ Self-test complete — all plots generated successfully.")
    import matplotlib.pyplot as plt
    plt.close('all')
