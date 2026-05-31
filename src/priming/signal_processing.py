#!/usr/bin/env python3
"""
Post-Priming Signal Processing
================================

Signal denoising and enhancement for cfDNA data collected after
priming agent administration. Removes priming-related artifacts
while amplifying true ctDNA signal.

Classes:
- PostPrimingDenoiser: Removes priming artifacts
- SignalEnhancer: Amplifies ctDNA signal
- BaselineCorrector: Corrects baseline drift
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple

import numpy as np
from scipy import signal as scipy_signal
from scipy import stats as scipy_stats

logger = logging.getLogger(__name__)


class PostPrimingDenoiser:
    """Denoise cfDNA signal after priming agent administration.

    Removes priming-induced background noise while preserving true
    ctDNA signal using multi-stage processing:
    1. Moving average smoothing
    2. Outlier detection and removal
    3. Trend correction
    """

    def __init__(self, window: int = 3):
        """
        Parameters
        ----------
        window : int
            Moving average window size in samples.
        """
        self.window = max(3, window)  # Minimum window of 3

    def denoise(
        self,
        cfDNA_signal: np.ndarray,
        priming_timing: Optional[float] = None,
    ) -> np.ndarray:
        """Denoise a cfDNA signal array.

        Parameters
        ----------
        cfDNA_signal : np.ndarray, shape (n_samples,) or (n_samples, n_features)
            Raw cfDNA signal (1D or 2D).
        priming_timing : float, optional
            Time (in samples) when priming agent was administered.
            Samples before this are treated as pre-priming baseline.

        Returns
        -------
        np.ndarray : Denoised signal, same shape as input.
        """
        signal = np.asarray(cfDNA_signal, dtype=np.float64)
        if signal.ndim > 2:
            raise ValueError(f"Expected 1D or 2D signal, got {signal.ndim}D")

        if signal.size == 0:
            return signal

        if signal.ndim == 1:
            return self._denoise_1d(signal, priming_timing)
        else:
            result = np.zeros_like(signal)
            for i in range(signal.shape[1]):
                result[:, i] = self._denoise_1d(signal[:, i], priming_timing)
            return result

    def _denoise_1d(
        self, signal: np.ndarray, priming_timing: Optional[float] = None
    ) -> np.ndarray:
        """1D denoising pipeline."""
        n = len(signal)

        # Stage 1: Moving average smoothing
        if n >= self.window:
            kernel = np.ones(self.window) / self.window
            smoothed = np.convolve(signal, kernel, mode="same")
            # Fix boundary effects
            smoothed[: self.window // 2] = signal[: self.window // 2]
            smoothed[-self.window // 2 :] = signal[-self.window // 2 :]
        else:
            smoothed = signal.copy()

        # Stage 2: Outlier detection using modified Z-score
        if n > 5:
            median = np.median(smoothed)
            mad = np.median(np.abs(smoothed - median))
            if mad > 1e-10:
                modified_z = 0.6745 * (smoothed - median) / mad
                outliers = np.abs(modified_z) > 3.5
                if np.any(outliers):
                    # Replace outliers with local median
                    cleaned = smoothed.copy()
                    for idx in np.where(outliers)[0]:
                        start = max(0, idx - 2)
                        end = min(n, idx + 3)
                        cleaned[idx] = np.median(smoothed[start:end])
                    smoothed = cleaned

        # Stage 3: Trend correction (remove linear drift)
        if n > 3:
            x = np.arange(n)
            try:
                slope, intercept, _, _, _ = scipy_stats.linregress(x, smoothed)
                trend = slope * x + intercept
                detrended = smoothed - trend + np.mean(smoothed)
                smoothed = detrended
            except (ValueError, scipy_stats.LinAlgError):
                pass

        # Stage 4: If priming timing provided, preserve pre-priming baseline
        if priming_timing is not None and priming_timing > 0:
            pre_region = slice(0, max(1, int(priming_timing)))
            pre_mean = np.mean(signal[pre_region])
            smoothed_mean = np.mean(smoothed[pre_region])
            if abs(smoothed_mean - pre_mean) > 1e-10:
                offset = pre_mean - smoothed_mean
                smoothed = smoothed + offset

        return smoothed

    def estimate_background_noise(self, signal: np.ndarray, window: Optional[int] = None) -> float:
        """Estimate background noise level.

        Uses median absolute deviation (MAD) as robust estimator.

        Parameters
        ----------
        signal : np.ndarray
            Input signal.
        window : int, optional
            Rolling window size; if None, uses global estimate.

        Returns
        -------
        float : Noise estimate.
        """
        signal = np.asarray(signal, dtype=np.float64)
        if signal.size == 0:
            return 0.0

        if signal.ndim > 1:
            signal = signal.flatten()

        if window is not None and len(signal) >= window:
            # Rolling MAD
            noises = []
            for i in range(0, len(signal) - window + 1, window // 2):
                chunk = signal[i : i + window]
                mad = np.median(np.abs(chunk - np.median(chunk)))
                noises.append(mad)
            return float(np.median(noises)) if noises else 0.0

        mad = np.median(np.abs(signal - np.median(signal)))
        return float(mad * 1.4826)  # Scale factor for normal distribution


class SignalEnhancer:
    """Enhance ctDNA signal by leveraging priming boost information.

    Amplifies true ctDNA features while suppressing background noise
    using adaptive thresholding and signal-to-noise weighting.
    """

    def __init__(self, boost_factor: float = 10.0):
        """
        Parameters
        ----------
        boost_factor : float
            Maximum expected ctDNA boost from priming.
        """
        self.boost_factor = boost_factor

    def enhance(
        self,
        features: np.ndarray,
        signal_to_noise: Optional[float] = None,
        priming_boost: Optional[float] = None,
        threshold_percentile: float = 50.0,
    ) -> np.ndarray:
        """Enhance ctDNA signal in feature array.

        Parameters
        ----------
        features : np.ndarray, shape (n_features,) or (n_samples, n_features)
            Feature array to enhance.
        signal_to_noise : float, optional
            Estimated signal-to-noise ratio. Higher → less suppression.
        priming_boost : float, optional
            Actual priming boost factor. Defaults to self.boost_factor.
        threshold_percentile : float
            Percentile for adaptive threshold (0-100).

        Returns
        -------
        np.ndarray : Enhanced features, same shape as input.
        """
        features = np.asarray(features, dtype=np.float64)
        if features.size == 0:
            return features

        boost = priming_boost if priming_boost is not None else self.boost_factor
        snr = signal_to_noise if signal_to_noise is not None else 1.0
        snr = max(0.01, snr)

        was_1d = features.ndim == 1
        if was_1d:
            features = features.reshape(1, -1)

        result = np.zeros_like(features)
        for i in range(features.shape[0]):
            row = features[i].copy()

            # Adaptive threshold: keep features above noise floor
            threshold = np.percentile(np.abs(row), threshold_percentile)
            mask = np.abs(row) >= threshold
            background = np.abs(row) < threshold

            # Suppress background: scale down noise components
            # SNR of 1 → 0.5x suppression; SNR >> 1 → minimal suppression
            suppression_factor = 1.0 / (1.0 + snr)
            row[background] *= suppression_factor

            # Amplify signal: boost true features proportional to priming effect
            # Cap amplification at boost_factor
            amp_factor = min(boost, self.boost_factor)
            row[mask] *= amp_factor

            # Weight by signal-to-noise
            snr_weight = snr / (snr + 1.0)
            result[i] = row * snr_weight + features[i] * (1 - snr_weight)

        if was_1d:
            result = result.flatten()

        return result


class BaselineCorrector:
    """Correct for baseline drift in repeated cfDNA sampling.

    When multiple blood draws are performed (pre-priming, post-priming,
    follow-up), baseline cfDNA levels may drift. This class corrects
    for that drift using pre-priming measurements as reference.
    """

    def __init__(self, correction_method: str = "subtractive"):
        """
        Parameters
        ----------
        correction_method : str
            "subtractive" or "ratio". Subtractive removes absolute drift;
            ratio normalizes by baseline.
        """
        if correction_method not in ("subtractive", "ratio"):
            raise ValueError(f"Unknown correction method: {correction_method}")
        self.correction_method = correction_method

    def correct(
        self,
        post_priming: np.ndarray,
        pre_priming_baseline: np.ndarray,
    ) -> np.ndarray:
        """Correct post-priming signal for baseline drift.

        Parameters
        ----------
        post_priming : np.ndarray
            Post-priming measurements.
        pre_priming_baseline : np.ndarray
            Pre-priming baseline measurements. Must be same shape as
            post_priming or broadcastable.

        Returns
        -------
        np.ndarray : Baseline-corrected signal.
        """
        post = np.asarray(post_priming, dtype=np.float64)
        pre = np.asarray(pre_priming_baseline, dtype=np.float64)

        if self.correction_method == "subtractive":
            # Remove absolute baseline drift
            if pre.ndim == 0 or (pre.ndim == 1 and pre.shape[0] == 1):
                corrected = post - float(pre.flat[0])
            else:
                corrected = post - pre
        else:
            # Ratio correction: normalize by baseline
            # Avoid division by zero
            pre_safe = np.where(np.abs(pre) < 1e-10, 1.0, pre)
            if pre.ndim == 0 or (pre.ndim == 1 and pre.shape[0] == 1):
                corrected = post / float(pre_safe.flat[0])
            else:
                corrected = post / pre_safe

        return corrected

    def estimate_drift(
        self,
        measurements: np.ndarray,
        times: Optional[np.ndarray] = None,
    ) -> float:
        """Estimate baseline drift rate.

        Parameters
        ----------
        measurements : np.ndarray
            Sequential measurements.
        times : np.ndarray, optional
            Corresponding time points. If None, assumes uniform spacing.

        Returns
        -------
        float : Drift rate per time unit.
        """
        measurements = np.asarray(measurements, dtype=np.float64).flatten()
        n = len(measurements)

        if n < 2:
            return 0.0

        if times is None:
            times = np.arange(n)

        try:
            slope, _, _, _, _ = scipy_stats.linregress(times, measurements)
            return float(slope)
        except (ValueError, scipy_stats.LinAlgError):
            return 0.0
