"""
Variational Autoencoder (beta-VAE) for unsupervised cancer anomaly detection.

Ported from ``cancer-screening/src/ml_models/vae_anomaly.py`` (SHA 697712be).

Trained primarily on healthy samples, this model learns the normal cfDNA
manifold. Samples far from the healthy manifold are flagged as anomalous --
potentially cancer, even if the cancer subtype wasn't in the training set.

Complements the supervised models by catching:
  - Rare cancer subtypes not represented in training data
  - Novel mutation patterns
  - Biological variation that looks like cancer but isn't
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
import logging
from typing import Optional

logger = logging.getLogger(__name__)


class AnomalyVAE(nn.Module):
    """
    beta-VAE for unsupervised anomaly detection in multi-omics cancer data.

    Detection: Samples with high reconstruction error + high KL divergence
    from the healthy prior are flagged as anomalous.

    beta > 1 forces a more disentangled latent space, sharpening anomaly scores.
    """

    def __init__(
        self,
        input_dim: int = 736,
        latent_dim: int = 32,
        hidden_dims: Optional[list] = None,
        beta: float = 4.0,
    ):
        if hidden_dims is None:
            hidden_dims = [256, 128]
        super().__init__()
        self.latent_dim = latent_dim
        self.beta = beta

        encoder_layers = []
        prev_dim = input_dim
        for h_dim in hidden_dims:
            encoder_layers.extend([
                nn.Linear(prev_dim, h_dim),
                nn.LayerNorm(h_dim),
                nn.GELU(),
                nn.Dropout(0.1),
            ])
            prev_dim = h_dim
        encoder_layers.append(nn.Linear(prev_dim, latent_dim * 2))
        self.encoder = nn.Sequential(*encoder_layers)

        decoder_layers = []
        prev_dim = latent_dim
        for h_dim in reversed(hidden_dims):
            decoder_layers.extend([
                nn.Linear(prev_dim, h_dim),
                nn.LayerNorm(h_dim),
                nn.GELU(),
                nn.Dropout(0.1),
            ])
            prev_dim = h_dim
        decoder_layers.append(nn.Linear(prev_dim, input_dim))
        self.decoder = nn.Sequential(*decoder_layers)

    def encode(self, x: torch.Tensor) -> tuple:
        params = self.encoder(x)
        mu = params[:, :self.latent_dim]
        logvar = params[:, self.latent_dim:]
        return mu, logvar

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def forward(self, x: torch.Tensor) -> dict:
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decode(z)
        return {'reconstruction': x_recon, 'mu': mu, 'logvar': logvar, 'z': z}

    def anomaly_score(self, x: torch.Tensor) -> torch.Tensor:
        """Per-sample anomaly score = beta * KL + reconstruction error."""
        mu, logvar = self.encode(x)
        z = self.reparameterize(mu, logvar)
        x_recon = self.decode(z)

        kl_div = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=-1)
        recon_error = F.mse_loss(x_recon, x, reduction='none').sum(dim=-1)

        anomaly = self.beta * kl_div + recon_error
        return anomaly


class AnomalyDetector:
    """
    Unsupervised cancer anomaly detector using beta-VAE.

    Usage:
        detector = AnomalyDetector()
        detector.train_on_healthy(healthy_features_df)  # No cancer labels needed
        scores = detector.score(samples_df)  # Higher = more anomalous
    """

    def __init__(
        self,
        latent_dim: int = 32,
        beta: float = 4.0,
        device: str = 'cuda' if torch.cuda.is_available() else 'cpu',
    ):
        self.device = device
        self.vae = AnomalyVAE(latent_dim=latent_dim, beta=beta).to(device)
        self.optimizer = torch.optim.Adam(self.vae.parameters(), lr=1e-3)
        self.is_trained = False
        self.threshold: Optional[float] = None

    def train(
        self,
        healthy_data: np.ndarray,
        epochs: int = 100,
        batch_size: int = 64,
    ) -> dict:
        """Train the VAE on healthy (normal) samples only."""
        dataset = torch.FloatTensor(healthy_data).to(self.device)
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        losses = []
        for epoch in range(epochs):
            epoch_loss = 0.0
            for batch in loader:
                self.optimizer.zero_grad()

                out = self.vae(batch)
                mu, logvar = out['mu'], out['logvar']

                recon_loss = F.mse_loss(out['reconstruction'], batch, reduction='mean')
                kl_loss = -0.5 * torch.mean(
                    1 + logvar - mu.pow(2) - logvar.exp()
                )

                loss = recon_loss + self.vae.beta * kl_loss
                loss.backward()
                self.optimizer.step()

                epoch_loss += loss.item()

            avg_loss = epoch_loss / len(loader)
            losses.append(avg_loss)

            if (epoch + 1) % 20 == 0:
                logger.info(f"VAE Epoch {epoch+1}: loss={avg_loss:.4f}")

        self.is_trained = True
        self.healthy_data = torch.FloatTensor(healthy_data).to(self.device)

        with torch.no_grad():
            healthy_scores = self.vae.anomaly_score(self.healthy_data)
            self.threshold = float(torch.quantile(healthy_scores, 0.99))

        logger.info(
            f"Training complete. Healthy anomaly threshold (99th pct): {self.threshold:.4f}"
        )

        return {'losses': losses, 'threshold': self.threshold}

    def score(self, data: np.ndarray) -> np.ndarray:
        """Score samples for anomaly. Higher = more anomalous."""
        if not self.is_trained:
            raise ValueError("Model not trained. Call train() first.")

        self.vae.eval()
        with torch.no_grad():
            batch = torch.FloatTensor(data).to(self.device)
            scores = self.vae.anomaly_score(batch)

        return scores.cpu().numpy()

    def predict(
        self,
        data: np.ndarray,
        threshold: Optional[float] = None,
    ) -> tuple:
        """Predict anomalies (binary) and return scores."""
        scores = self.score(data)
        thresh = threshold if threshold is not None else self.threshold
        is_anomaly = (scores > thresh).astype(int)
        return is_anomaly, scores
