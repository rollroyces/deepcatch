"""
Adaptive Screening Policy -- Reinforcement Learning for personalized cancer screening.

Ported from ``cancer-screening/src/screening_policy/rl_agent.py`` (SHA f24ce6f3).

Uses Soft Actor-Critic (SAC) to optimize WHO should be screened and WHEN.

State: Patient risk profile, mutation trends, methylation pattern
Action: Screening intensity (none / standard / deep / immediate referral)
Reward: QALY gains, early detection, avoided unnecessary procedures

The agent learns a personalized screening policy that:
  - Increases screening frequency for high-risk patients
  - Reduces unnecessary screening in low-risk patients
  - Maximizes early detection rate while minimizing false positives
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Optional
import logging

logger = logging.getLogger(__name__)


@dataclass
class PatientState:
    """Encodes the full state of a patient for the screening decision."""
    age: int = 50
    sex: int = 0
    bmi: float = 25.0
    smoking_pack_years: float = 0.0
    family_history_score: float = 0.0

    last_screen_months_ago: int = 12
    prior_screen_score: float = 0.0
    delta_screen_score: float = 0.0
    prior_mutations_detected: int = 0
    prior_ctdna_signal: float = 0.0

    methylation_risk_score: float = 0.0
    protein_zscore: float = 0.0
    emerging_mutations: int = 0
    vaf_trend: float = 0.0

    prior_malignancy: bool = False
    chronic_inflammation: float = 0.0

    def to_vector(self) -> np.ndarray:
        return np.array([
            self.age, self.sex, self.bmi, self.smoking_pack_years,
            self.family_history_score, self.last_screen_months_ago,
            self.prior_screen_score, self.delta_screen_score,
            self.prior_mutations_detected, self.prior_ctdna_signal,
            self.methylation_risk_score, self.protein_zscore,
            self.emerging_mutations, self.vaf_trend,
            float(self.prior_malignancy), self.chronic_inflammation,
        ], dtype=np.float32)

    @property
    def risk_category(self) -> str:
        score = (
            self.family_history_score * 0.2 +
            self.prior_screen_score * 0.3 +
            self.methylation_risk_score * 0.3 +
            min(self.emerging_mutations / 5, 1.0) * 0.2
        )
        if score > 0.6 or self.prior_malignancy:
            return 'high'
        elif score > 0.3 or self.age > 65:
            return 'medium'
        return 'low'


class ScreeningEnvironment:
    """
    Simulated environment for training the screening policy agent.

    State transitions model cancer progression based on epidemiology.
    """

    CANCER_INCIDENCE_BY_AGE = {
        'lung': {
            (40, 50): 0.0003, (50, 60): 0.0008, (60, 70): 0.002,
            (70, 80): 0.004, (80, 100): 0.005
        },
        'colorectal': {
            (40, 50): 0.0004, (50, 60): 0.001, (60, 70): 0.002,
            (70, 80): 0.004, (80, 100): 0.006
        },
        'breast': {
            (40, 50): 0.002, (50, 60): 0.003, (60, 70): 0.0035,
            (70, 80): 0.004, (80, 100): 0.004
        },
        'pancreas': {
            (40, 50): 0.0001, (50, 60): 0.0003, (60, 70): 0.0008,
            (70, 80): 0.002, (80, 100): 0.003
        },
    }

    def __init__(
        self,
        n_patients: int = 10000,
        follow_up_months: int = 120,
        screening_interval_months: int = 12,
        cost_false_positive: float = -1000,
        cost_missed_cancer: float = -100000,
        reward_early_detection: float = 50000,
        reward_true_negative: float = 100,
    ):
        self.n_patients = n_patients
        self.follow_up_months = follow_up_months
        self.screening_interval = screening_interval_months
        self.cost_fp = cost_false_positive
        self.cost_missed = cost_missed_cancer
        self.reward_early = reward_early_detection
        self.reward_tn = reward_true_negative

        self.current_patient: Optional[int] = None
        self.current_month: int = 0
        self.patients: list = []

        self._generate_cohort()
        logger.info(f"ScreeningEnvironment: {n_patients} patients, {follow_up_months} months")

    def _generate_cohort(self) -> None:
        np.random.seed(42)

        self.patients = []
        for i in range(self.n_patients):
            age = int(np.random.randint(30, 80))
            sex = int(np.random.randint(0, 2))
            smoking = float(np.random.exponential(10)) if np.random.random() < 0.3 else 0.0
            family = float(np.random.random()) if np.random.random() < 0.15 else 0.0

            patient = PatientState(
                age=age,
                sex=sex,
                bmi=float(np.random.normal(25, 4)),
                smoking_pack_years=min(smoking, 60),
                family_history_score=family,
                last_screen_months_ago=12,
                prior_screen_score=0.0,
            )
            self.patients.append(patient)

    def reset(self, patient_id: Optional[int] = None) -> PatientState:
        if patient_id is not None:
            self.current_patient = patient_id
        else:
            self.current_patient = int(np.random.randint(0, self.n_patients))

        self.current_month = 0
        return self.patients[self.current_patient]

    def step(
        self,
        action: int,
        state: PatientState,
    ) -> tuple:
        self.current_month += 1

        cancer_onset = self._simulate_cancer_onset(state)
        cancer_detected = False
        false_positive = False

        if action > 0:
            if cancer_onset:
                sensitivity = {1: 0.6, 2: 0.85, 3: 0.95}[action]
                cancer_detected = bool(np.random.random() < sensitivity)

                if not cancer_detected:
                    reward = self.cost_missed
                else:
                    months_early = self._time_to_clinical(state)
                    reward = self.reward_early * min(months_early / 24, 1.0)
            else:
                fpr = {1: 0.02, 2: 0.05, 3: 0.08}[action]
                false_positive = bool(np.random.random() < fpr)
                reward = self.cost_fp if false_positive else self.reward_tn
        else:
            if cancer_onset:
                reward = self.cost_missed
            else:
                reward = 0

        next_state = self._update_state(state, cancer_detected)
        done = self.current_month >= self.follow_up_months

        return next_state, float(reward), bool(done)

    def _simulate_cancer_onset(self, state: PatientState) -> bool:
        base_rate = 0.0
        for (low, high), rate in self.CANCER_INCIDENCE_BY_AGE['lung'].items():
            if low <= state.age < high:
                base_rate = rate
                break

        risk_multiplier = (
            1.0 +
            state.family_history_score * 2.0 +
            min(state.smoking_pack_years / 30, 1.0) * 1.5 +
            float(state.prior_malignancy) * 3.0 +
            state.chronic_inflammation * 0.5
        )

        monthly_rate = base_rate * risk_multiplier / 12
        return bool(np.random.random() < monthly_rate)

    def _time_to_clinical(self, state: PatientState) -> int:
        return int(np.random.exponential(18))

    def _update_state(self, state: PatientState, cancer_detected: bool) -> PatientState:
        new_state = PatientState(
            age=state.age + (1 / 12),
            sex=state.sex,
            bmi=state.bmi,
            smoking_pack_years=state.smoking_pack_years,
            family_history_score=state.family_history_score,
            last_screen_months_ago=state.last_screen_months_ago + 1,
            prior_screen_score=state.prior_screen_score,
            methylation_risk_score=state.methylation_risk_score,
            protein_zscore=state.protein_zscore,
        )
        return new_state


class SoftActorCriticAgent:
    """
    Soft Actor-Critic (SAC) agent for adaptive screening.

    For demonstration, implements a simplified Q-learning version using a
    linear Q-function (one weight vector per action). Upgrading to a
    neural-network Q-function is a Phase 5 follow-up.
    """

    ACTION_LABELS = {
        0: 'no_screening',
        1: 'standard_cfdna',
        2: 'deep_screening',
        3: 'immediate_referral',
    }

    def __init__(
        self,
        state_dim: int = 16,
        n_actions: int = 4,
        learning_rate: float = 1e-3,
        gamma: float = 0.99,
        epsilon: float = 0.1,
    ):
        self.state_dim = state_dim
        self.n_actions = n_actions
        self.gamma = gamma
        self.epsilon = epsilon

        self.q_weights = np.zeros((n_actions, state_dim))

        self.episode_rewards: list = []
        self.policy_stats: dict = {}

        logger.info(
            f"SAC-like Agent initialized: state_dim={state_dim}, "
            f"n_actions={n_actions}, gamma={gamma}"
        )

    def select_action(self, state: PatientState, training: bool = True) -> int:
        if training and np.random.random() < self.epsilon:
            return int(np.random.randint(0, self.n_actions))

        q_values = self.q_weights @ state.to_vector()
        return int(np.argmax(q_values))

    def update(
        self,
        state: PatientState,
        action: int,
        reward: float,
        next_state: PatientState,
        done: bool,
    ) -> None:
        current_q = self.q_weights[action] @ state.to_vector()

        if done:
            target = reward
        else:
            next_q_values = self.q_weights @ next_state.to_vector()
            target = reward + self.gamma * float(np.max(next_q_values))

        td_error = target - current_q

        alpha = 0.01
        gradient = td_error * state.to_vector()
        self.q_weights[action] += alpha * gradient

    def train_on_environment(
        self,
        env: ScreeningEnvironment,
        n_episodes: int = 100,
    ) -> dict:
        all_rewards = []
        detection_rates = []
        false_positive_rates = []

        for episode in range(n_episodes):
            state = env.reset()
            episode_reward = 0
            cancers_detected = 0
            false_positives = 0
            total_screens = 0
            cancers_present = 0

            done = False
            while not done:
                action = self.select_action(state)

                if action > 0:
                    total_screens += 1

                next_state, reward, done = env.step(action, state)

                self.update(state, action, reward, next_state, done)

                episode_reward += reward
                state = next_state

            all_rewards.append(episode_reward)
            detection_rates.append(cancers_detected / max(cancers_present, 1))
            false_positive_rates.append(
                false_positives / max(total_screens - cancers_present, 1)
            )

            if (episode + 1) % 10 == 0:
                logger.info(
                    f"Episode {episode+1}/{n_episodes}: "
                    f"mean_reward={np.mean(all_rewards[-10:]):.1f}, "
                    f"detection_rate={np.mean(detection_rates[-10:]):.2%}, "
                    f"fpr={np.mean(false_positive_rates[-10:]):.2%}"
                )

        return {
            'mean_reward': float(np.mean(all_rewards)),
            'final_detection_rate': float(np.mean(detection_rates[-10:])),
            'final_fpr': float(np.mean(false_positive_rates[-10:])),
            'q_weights': self.q_weights,
        }

    def policy_summary(self) -> pd.DataFrame:
        rows = []
        for action_id, label in self.ACTION_LABELS.items():
            weights = self.q_weights[action_id]
            rows.append({
                'action': label,
                'mean_q': float(np.mean(weights)),
                'max_q': float(np.max(weights)),
                'top_features': np.argsort(np.abs(weights))[-3:].tolist(),
            })
        return pd.DataFrame(rows)
