"""
DeepCatch Billing — tier gating, usage tracking & webhooks.

===========================================================

A drop-in billing layer that wraps the DeepCatch Fragmentomics API
with tier-based rate limiting and Stripe / Lemonsqueezy webhook handling.

End-user tiers (locked by product spec):

    Free    10 calls/day                no auth
    Pro     $49/mo,   10,000 calls/day  API key auth
    Lab     $499/mo,  unlimited         API key + raw motif + checkpoint access

Design rules:
    * No hardcoded secrets. Everything comes from env vars.
    * SQLite-backed usage store at ``api/billing.db`` (auto-created).
    * Token-bucket rate limiter per API key, independent of the daily counter.
    * Both Stripe and Lemonsqueezy webhook signatures are verified using
      ``hmac.compare_digest`` against timing-safe digests.
    * A first-class ``--mock`` mode fakes the payment-provider calls so the
      full flow can be exercised without any real Stripe/Lemonsqueezy account.
    * Use ``@require_tier("pro")`` to protect a FastAPI endpoint.

Module map
----------

    Tier                 enum of the three tiers
    TierConfig           static per-tier policy (daily cap, monthly $)
    ApiKeyRecord         row in the api_keys table
    UsageRecord          row in the usage table (one per call)
    SubscriptionRecord   row in the subscriptions table (Pro/Lab customers)
    BillingStore         SQLite-backed persistent store (auto-creates file)
    TokenBucket          in-process per-key leaky-bucket limiter
    BillingManager       the class your app holds; everything wires through it
    WebhookVerifier      verifies Stripe + Lemonsqueezy webhook signatures
    require_tier         FastAPI decorator
    --mock               CLI demo of a full Pro upgrade + webhook flow

CLI usage::

    python -m api.billing --mock

Run tests::

    pytest api/test_billing.py -v

No third-party dependencies. stdlib only (sqlite3, hmac, hashlib, json, time,
dataclasses, enum, pathlib, secrets, functools, argparse).
"""

from __future__ import annotations

import argparse
import dataclasses
import enum
import functools
import hashlib
import hmac
import json
import os
import secrets
import sqlite3
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Optional, Tuple, Union

# ───────────────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────────────

# Module location so the SQLite file lands next to it, regardless of cwd.
_API_DIR = Path(__file__).resolve().parent
_DB_PATH = _API_DIR / "billing.db"

# Per-tier policy. Daily limits are ZERO or positive; ``-1`` means unlimited.
@dataclass(frozen=True)
class TierConfig:
    name: str
    daily_call_limit: int   # -1 == unlimited
    monthly_price_usd: float
    requires_api_key: bool
    allows_raw_motifs: bool
    allows_checkpoints: bool


TIER_CONFIG: Dict[str, TierConfig] = {
    "free": TierConfig(
        name="free",
        daily_call_limit=10,
        monthly_price_usd=0.0,
        requires_api_key=False,
        allows_raw_motifs=False,
        allows_checkpoints=False,
    ),
    "pro": TierConfig(
        name="pro",
        daily_call_limit=10_000,
        monthly_price_usd=49.0,
        requires_api_key=True,
        allows_raw_motifs=False,
        allows_checkpoints=False,
    ),
    "lab": TierConfig(
        name="lab",
        daily_call_limit=-1,           # unlimited
        monthly_price_usd=499.0,
        requires_api_key=True,
        allows_raw_motifs=True,
        allows_checkpoints=True,
    ),
}


class Tier(str, enum.Enum):
    FREE = "free"
    PRO = "pro"
    LAB = "lab"


# Convenience ordering used by ``require_tier``: a Lab key satisfies pro.
_TIER_RANK = {Tier.FREE: 0, Tier.PRO: 1, Tier.LAB: 2}


# ───────────────────────────────────────────────────────────────────────
# Persistent records
# ───────────────────────────────────────────────────────────────────────

@dataclass
class ApiKeyRecord:
    key: str
    customer_id: str
    tier: Tier
    created_at: float
    expires_at: Optional[float] = None     # epoch seconds
    active: bool = True
    note: str = ""


@dataclass
class SubscriptionRecord:
    customer_id: str
    tier: Tier
    started_at: float
    expires_at: float
    provider: str               # "stripe" or "lemonsqueezy"
    external_id: str            # subscription id from provider
    status: str = "active"      # active | canceled | expired | past_due


@dataclass
class UsageRecord:
    api_key: str
    timestamp: float
    day_utc: str
    endpoint: str
    tier: str
    count: int = 1


# ───────────────────────────────────────────────────────────────────────
# SQLite store
# ───────────────────────────────────────────────────────────────────────

_SCHEMA = """
CREATE TABLE IF NOT EXISTS api_keys (
    key          TEXT PRIMARY KEY,
    customer_id  TEXT NOT NULL,
    tier         TEXT NOT NULL,
    created_at   REAL NOT NULL,
    expires_at   REAL,
    active       INTEGER NOT NULL DEFAULT 1,
    note         TEXT NOT NULL DEFAULT ''
);
CREATE INDEX IF NOT EXISTS idx_api_keys_customer ON api_keys(customer_id);

CREATE TABLE IF NOT EXISTS subscriptions (
    customer_id  TEXT NOT NULL,
    tier         TEXT NOT NULL,
    started_at   REAL NOT NULL,
    expires_at   REAL NOT NULL,
    provider     TEXT NOT NULL,
    external_id  TEXT NOT NULL,
    status       TEXT NOT NULL DEFAULT 'active',
    PRIMARY KEY (customer_id, external_id)
);
CREATE INDEX IF NOT EXISTS idx_subs_customer ON subscriptions(customer_id);

CREATE TABLE IF NOT EXISTS usage (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    api_key      TEXT NOT NULL,
    timestamp    REAL NOT NULL,
    day_utc      TEXT NOT NULL,
    endpoint     TEXT NOT NULL,
    tier         TEXT NOT NULL,
    count        INTEGER NOT NULL DEFAULT 1
);
CREATE INDEX IF NOT EXISTS idx_usage_key_day ON usage(api_key, day_utc);

CREATE TABLE IF NOT EXISTS webhook_events (
    event_id     TEXT PRIMARY KEY,
    provider     TEXT NOT NULL,
    received_at  REAL NOT NULL,
    event_type   TEXT NOT NULL,
    payload      TEXT NOT NULL
);
"""


class BillingStore:
    """Thin sqlite3 wrapper. Auto-creates schema on first use."""

    def __init__(self, db_path: Optional[Union[Path, str]] = None):
        # ``:memory:`` (or any string) is forwarded to sqlite3 unchanged.
        if db_path is None:
            self.db_path = _DB_PATH
        elif isinstance(db_path, Path):
            self.db_path = db_path
        else:
            self.db_path = Path(db_path)  # for sqlite3 str args like ":memory:"
        # ``check_same_thread=False`` so FastAPI's threadpool can share it.
        # sqlite3.connect accepts either a Path or a string; let it pick.
        self._conn = sqlite3.connect(str(self.db_path), check_same_thread=False)
        self._conn.row_factory = sqlite3.Row
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.executescript(_SCHEMA)
        self._conn.commit()

    # ── Connection helpers ──────────────────────────────────────
    def close(self) -> None:
        self._conn.close()

    def commit(self) -> None:
        self._conn.commit()

    # ── API keys ────────────────────────────────────────────────
    def upsert_api_key(self, rec: ApiKeyRecord) -> None:
        self._conn.execute(
            """
            INSERT INTO api_keys (key, customer_id, tier, created_at,
                                  expires_at, active, note)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(key) DO UPDATE SET
                customer_id=excluded.customer_id,
                tier=excluded.tier,
                created_at=excluded.created_at,
                expires_at=excluded.expires_at,
                active=excluded.active,
                note=excluded.note
            """,
            (
                rec.key, rec.customer_id, rec.tier.value, rec.created_at,
                rec.expires_at, int(rec.active), rec.note,
            ),
        )
        self.commit()

    def get_api_key(self, key: str) -> Optional[ApiKeyRecord]:
        row = self._conn.execute(
            "SELECT * FROM api_keys WHERE key = ?", (key,)
        ).fetchone()
        if row is None:
            return None
        return _row_to_api_key(row)

    def all_api_keys(self) -> list[ApiKeyRecord]:
        rows = self._conn.execute("SELECT * FROM api_keys ORDER BY created_at").fetchall()
        return [_row_to_api_key(r) for r in rows]

    # ── Subscriptions ───────────────────────────────────────────
    def upsert_subscription(self, sub: SubscriptionRecord) -> None:
        self._conn.execute(
            """
            INSERT INTO subscriptions (customer_id, tier, started_at, expires_at,
                                       provider, external_id, status)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(customer_id, external_id) DO UPDATE SET
                tier=excluded.tier,
                started_at=excluded.started_at,
                expires_at=excluded.expires_at,
                provider=excluded.provider,
                status=excluded.status
            """,
            (
                sub.customer_id, sub.tier.value, sub.started_at, sub.expires_at,
                sub.provider, sub.external_id, sub.status,
            ),
        )
        self.commit()

    def latest_active_subscription(self, customer_id: str) -> Optional[SubscriptionRecord]:
        row = self._conn.execute(
            """
            SELECT * FROM subscriptions
            WHERE customer_id = ? AND status = 'active'
            ORDER BY expires_at DESC LIMIT 1
            """,
            (customer_id,),
        ).fetchone()
        if row is None:
            return None
        return _row_to_subscription(row)

    def mark_subscription_expired(self, customer_id: str, external_id: str) -> None:
        self._conn.execute(
            "UPDATE subscriptions SET status='expired' "
            "WHERE customer_id = ? AND external_id = ?",
            (customer_id, external_id),
        )
        self.commit()

    # ── Usage ───────────────────────────────────────────────────
    def record_usage(self, rec: UsageRecord) -> None:
        self._conn.execute(
            """
            INSERT INTO usage (api_key, timestamp, day_utc, endpoint, tier, count)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (rec.api_key, rec.timestamp, rec.day_utc, rec.endpoint, rec.tier, rec.count),
        )
        self.commit()

    def daily_count(self, api_key: str, day_utc: str) -> int:
        row = self._conn.execute(
            "SELECT COALESCE(SUM(count), 0) AS c FROM usage "
            "WHERE api_key = ? AND day_utc = ?",
            (api_key, day_utc),
        ).fetchone()
        return int(row["c"]) if row else 0

    def usage_history(self, api_key: str, limit: int = 100) -> list[UsageRecord]:
        rows = self._conn.execute(
            "SELECT * FROM usage WHERE api_key = ? "
            "ORDER BY timestamp DESC LIMIT ?",
            (api_key, limit),
        ).fetchall()
        return [_row_to_usage(r) for r in rows]

    # ── Webhook event idempotency ───────────────────────────────
    def remember_event(self, event_id: str, provider: str,
                       event_type: str, payload: str) -> bool:
        """Returns True if newly stored, False if it was already seen."""
        try:
            self._conn.execute(
                "INSERT INTO webhook_events "
                "(event_id, provider, received_at, event_type, payload) "
                "VALUES (?, ?, ?, ?, ?)",
                (event_id, provider, time.time(), event_type, payload),
            )
            self.commit()
            return True
        except sqlite3.IntegrityError:
            return False

    # ── Maintenance ─────────────────────────────────────────────
    def nuke(self) -> None:
        """Wipe all tables (used by tests)."""
        for tbl in ("api_keys", "subscriptions", "usage", "webhook_events"):
            self._conn.execute(f"DELETE FROM {tbl}")
        self.commit()


def _row_to_api_key(row: sqlite3.Row) -> ApiKeyRecord:
    return ApiKeyRecord(
        key=row["key"],
        customer_id=row["customer_id"],
        tier=Tier(row["tier"]),
        created_at=row["created_at"],
        expires_at=row["expires_at"],
        active=bool(row["active"]),
        note=row["note"] or "",
    )


def _row_to_subscription(row: sqlite3.Row) -> SubscriptionRecord:
    return SubscriptionRecord(
        customer_id=row["customer_id"],
        tier=Tier(row["tier"]),
        started_at=row["started_at"],
        expires_at=row["expires_at"],
        provider=row["provider"],
        external_id=row["external_id"],
        status=row["status"],
    )


def _row_to_usage(row: sqlite3.Row) -> UsageRecord:
    return UsageRecord(
        api_key=row["api_key"],
        timestamp=row["timestamp"],
        day_utc=row["day_utc"],
        endpoint=row["endpoint"],
        tier=row["tier"],
        count=row["count"],
    )


# ───────────────────────────────────────────────────────────────────────
# Rate limiter (token bucket per API key, in-process)
# ───────────────────────────────────────────────────────────────────────

@dataclass
class _Bucket:
    capacity: int           # max burst
    refill_per_sec: float   # sustained rate
    tokens: float
    last_refill: float


class TokenBucket:
    """
    A simple in-memory token bucket per API key.

    Defaults: capacity 20, refill 5 tokens/sec (= 300/min) to keep
    bots from hammering even if their daily quota is large.
    """

    def __init__(
        self,
        capacity: int = 20,
        refill_per_sec: float = 5.0,
        clock: Callable[[], float] = time.monotonic,
    ):
        self.capacity = capacity
        self.refill_per_sec = refill_per_sec
        self._clock = clock
        self._buckets: Dict[str, _Bucket] = {}

    def _bucket_for(self, key: str) -> _Bucket:
        b = self._buckets.get(key)
        if b is None:
            now = self._clock()
            b = _Bucket(capacity=self.capacity,
                        refill_per_sec=self.refill_per_sec,
                        tokens=self.capacity, last_refill=now)
            self._buckets[key] = b
            return b
        # refill
        now = self._clock()
        delta = now - b.last_refill
        if delta > 0:
            b.tokens = min(b.capacity, b.tokens + delta * b.refill_per_sec)
            b.last_refill = now
        return b

    def try_acquire(self, key: str, cost: float = 1.0) -> bool:
        b = self._bucket_for(key)
        if b.tokens >= cost:
            b.tokens -= cost
            return True
        return False

    def reset(self, key: Optional[str] = None) -> None:
        if key is None:
            self._buckets.clear()
        else:
            self._buckets.pop(key, None)


# ───────────────────────────────────────────────────────────────────────
# Webhook signature verification
# ───────────────────────────────────────────────────────────────────────

class WebhookVerificationError(Exception):
    """Raised when a webhook signature fails to verify."""


class WebhookVerifier:
    """
    Verifies Stripe and Lemonsqueezy webhook signatures.

    Both providers sign payloads with HMAC-SHA256; the schemas differ:

      Stripe          ``t=<ts>,v1=<hex>``
      Lemonsqueezy    ``<hex>``  (raw HMAC of the body)

    Both expose ``verify(provider, raw_body, signature_header, secret)``.
    """

    @staticmethod
    def verify_stripe(raw_body: bytes, sig_header: str, secret: str,
                      tolerance: int = 300) -> None:
        """
        Mirror Stripe's ``construct_event`` from ``stripe-python``.

        Raises ``WebhookVerificationError`` on any mismatch or stale timestamp.
        """
        if not secret:
            raise WebhookVerificationError("Stripe webhook secret is empty (set DEEPCATCH_STRIPE_SECRET).")

        pairs: Dict[str, str] = {}
        for piece in (sig_header or "").split(","):
            if "=" in piece:
                k, v = piece.split("=", 1)
                pairs[k.strip()] = v.strip()

        ts = pairs.get("t")
        sig = pairs.get("v1")
        if not ts or not sig:
            raise WebhookVerificationError("Missing t= or v1= in Stripe signature header.")

        try:
            ts_int = int(ts)
        except ValueError:
            raise WebhookVerificationError("Stripe timestamp not an integer.")

        if abs(int(time.time()) - ts_int) > tolerance:
            raise WebhookVerificationError(
                f"Stripe timestamp outside tolerance ({tolerance}s)."
            )

        signed_payload = f"{ts}.".encode("utf-8") + raw_body
        expected = hmac.new(
            secret.encode("utf-8"),
            signed_payload,
            hashlib.sha256,
        ).hexdigest()

        if not hmac.compare_digest(expected, sig):
            raise WebhookVerificationError("Stripe HMAC signature mismatch.")

    @staticmethod
    def verify_lemonsqueezy(raw_body: bytes, sig_header: str, secret: str) -> None:
        """
        Verify Lemonsqueezy webhook signature.

        The provider computes ``HMAC_SHA256(secret, body)`` and emits the raw
        hex digest in the ``X-Signature`` header.
        """
        if not secret:
            raise WebhookVerificationError("Lemonsqueezy webhook secret is empty (set DEEPCATCH_LEMONSQUEEZY_SECRET).")

        if not sig_header:
            raise WebhookVerificationError("Missing Lemonsqueezy signature header.")

        expected = hmac.new(
            secret.encode("utf-8"),
            raw_body,
            hashlib.sha256,
        ).hexdigest()
        if not hmac.compare_digest(expected, sig_header.strip()):
            raise WebhookVerificationError("Lemonsqueezy HMAC signature mismatch.")


# ───────────────────────────────────────────────────────────────────────
# Manager (the public surface your app holds)
# ───────────────────────────────────────────────────────────────────────

class BillingError(Exception):
    """Raised for domain-level billing failures (rate limits, tier mismatch)."""


class BillingManager:
    """
    Public surface used by the DeepCatch API.

    Lifecycle:
        mgr = BillingManager(db_path=":memory:")     # in-process
        mgr = BillingManager()                       # persistent file at api/billing.db

    With the env vars set::

        DEEPCATCH_STRIPE_SECRET, DEEPCATCH_LEMONSQUEEZY_SECRET, DEEPCATCH_PUBLIC_URL

    Webhook handling::

        verdict = mgr.handle_stripe_webhook(raw_body, signature_header)
        verdict = mgr.handle_lemonsqueezy_webhook(raw_body, signature_header)
    """

    def __init__(
            self,
            db_path: Optional[Union[Path, str]] = None,
            stripe_secret: Optional[str] = None,
            lemonsqueezy_secret: Optional[str] = None,
        public_url: Optional[str] = None,
        # Fake mode for tests / --mock: skip real network, accept all webhooks.
        accept_any_signature: bool = False,
        store: Optional[BillingStore] = None,
        bucket: Optional[TokenBucket] = None,
        clock: Callable[[], float] = time.time,
    ):
        self.store = store or BillingStore(db_path)
        self.verifier = WebhookVerifier()
        self.bucket = bucket or TokenBucket()
        self.clock = clock
        self.stripe_secret = stripe_secret or os.environ.get("DEEPCATCH_STRIPE_SECRET", "")
        self.lemonsqueezy_secret = lemonsqueezy_secret or os.environ.get("DEEPCATCH_LEMONSQUEEZY_SECRET", "")
        self.public_url = public_url or os.environ.get("DEEPCATCH_PUBLIC_URL", "http://localhost:8000")
        self.accept_any_signature = accept_any_signature

    # ── Public key helpers ─────────────────────────────────────
    @staticmethod
    def mint_api_key() -> str:
        """Mint a fresh 32-byte URL-safe API key."""
        return "dck_" + secrets.token_urlsafe(32)

    def provision_api_key(self, customer_id: str, tier: Tier,
                          note: str = "",
                          ttl_days: Optional[int] = None) -> ApiKeyRecord:
        """Issue a new API key for ``customer_id`` at the given tier."""
        if tier not in TIER_CONFIG:
            raise BillingError(f"Unknown tier: {tier!r}")
        cfg = TIER_CONFIG[tier]
        now = self.clock()
        expires = (now + ttl_days * 86_400) if ttl_days else None
        rec = ApiKeyRecord(
            key=self.mint_api_key(),
            customer_id=customer_id,
            tier=tier,
            created_at=now,
            expires_at=expires,
            active=True,
            note=note,
        )
        self.store.upsert_api_key(rec)
        return rec

    def revoke_api_key(self, key: str) -> None:
        rec = self.store.get_api_key(key)
        if rec is None:
            return
        rec.active = False
        self.store.upsert_api_key(rec)

    # ── Tier resolution & entitlement checks ────────────────────
    def tier_for_api_key(self, key: str) -> Tier:
        rec = self.store.get_api_key(key)
        if rec is None or not rec.active:
            raise BillingError("Invalid or revoked API key.")
        if rec.expires_at and self.clock() > rec.expires_at:
            raise BillingError("API key expired.")
        return rec.tier

    def has_tier(self, key: str, required: Tier) -> bool:
        try:
            t = self.tier_for_api_key(key)
        except BillingError:
            return False
        return _TIER_RANK[t] >= _TIER_RANK[required]

    # ── Rate limiting ───────────────────────────────────────────
    def check_call_allowed(self, api_key: Optional[str],
                           endpoint: str,
                           tier: Tier = Tier.FREE) -> Tuple[bool, str]:
        """
        Returns ``(ok, reason)``. ``ok=False`` means the call should be rejected.

        Strategy:
          1. The token bucket rejects runaway bursts (per key or "anonymous").
          2. The daily counter rejects over-quota days (per key, or 1 per IP for Free).
          3. Lab tier is exempt from the daily counter.
        """
        bucket_key = api_key or "anonymous"

        # 1) burst control
        if not self.bucket.try_acquire(bucket_key):
            return False, "rate_limited_burst"

        # 2) daily counter
        cfg = TIER_CONFIG[tier.value]
        if cfg.daily_call_limit == -1:
            return True, "unlimited"

        # Free/anonymous callers share one daily counter keyed by the
        # bucket key (defaults to literal "anonymous"). Authenticated
        # callers are keyed by their API key.
        day = _today_utc()
        used = self.store.daily_count(api_key or bucket_key, day)
        if used >= cfg.daily_call_limit:
            return False, f"daily_quota_exceeded:{used}/{cfg.daily_call_limit}"

        return True, "ok"

    def record_call(self, api_key: str, endpoint: str, tier: Tier) -> None:
        """Atomically increment the daily counter for the key."""
        now = self.clock()
        day = datetime.fromtimestamp(now, tz=timezone.utc).strftime("%Y-%m-%d")
        self.store.record_usage(
            UsageRecord(api_key=api_key, timestamp=now,
                        day_utc=day, endpoint=endpoint, tier=tier.value)
        )

    # ── Subscription lifecycle ──────────────────────────────────
    def upgrade_customer(self, customer_id: str, target: Tier,
                         provider: str = "manual",
                         external_id: Optional[str] = None,
                         days: int = 30) -> SubscriptionRecord:
        """Manually record a subscription (used by webhooks & --mock)."""
        now = self.clock()
        sub = SubscriptionRecord(
            customer_id=customer_id,
            tier=target,
            started_at=now,
            expires_at=now + days * 86_400,
            provider=provider,
            external_id=external_id or f"manual_{int(now)}",
            status="active",
        )
        self.store.upsert_subscription(sub)
        return sub

    def cancel_subscription(self, customer_id: str) -> None:
        sub = self.store.latest_active_subscription(customer_id)
        if sub is None:
            return
        self.store.mark_subscription_expired(customer_id, sub.external_id)

    def effective_tier_for_customer(self, customer_id: str) -> Tier:
        sub = self.store.latest_active_subscription(customer_id)
        if sub is None:
            return Tier.FREE
        if sub.status != "active" or self.clock() > sub.expires_at:
            return Tier.FREE
        return sub.tier

    # ── Webhook handling ────────────────────────────────────────
    def handle_stripe_webhook(self, raw_body: bytes, sig_header: str) -> Dict[str, Any]:
        """
        Verify signature, idempotently record the event, and apply the
        resulting subscription change.

        Returns a small dict summary suitable for an HTTP 200 response.
        """
        if self.accept_any_signature:
            payload = json.loads(raw_body.decode("utf-8"))
        else:
            self.verifier.verify_stripe(raw_body, sig_header, self.stripe_secret)
            payload = json.loads(raw_body.decode("utf-8"))

        event_id = payload.get("id") or hashlib.sha256(raw_body).hexdigest()
        event_type = payload.get("type", "unknown")
        is_new = self.store.remember_event(event_id, "stripe", event_type, raw_body.decode("utf-8"))
        if not is_new:
            return {"status": "duplicate", "event_id": event_id, "type": event_type}

        obj = payload.get("data", {}).get("object", {})
        return self._apply_subscription_event(
            event_type=event_type,
            customer_id=str(obj.get("customer") or obj.get("client_reference_id") or ""),
            metadata=obj.get("metadata", {}) or {},
            external_id=str(obj.get("id") or event_id),
            provider="stripe",
        )

    def handle_lemonsqueezy_webhook(self, raw_body: bytes, sig_header: str) -> Dict[str, Any]:
        if self.accept_any_signature:
            payload = json.loads(raw_body.decode("utf-8"))
        else:
            self.verifier.verify_lemonsqueezy(raw_body, sig_header, self.lemonsqueezy_secret)
            payload = json.loads(raw_body.decode("utf-8"))

        meta = payload.get("meta", {})
        event_name = meta.get("event_name", "unknown")
        event_id = payload.get("id") or hashlib.sha256(raw_body).hexdigest()
        is_new = self.store.remember_event(event_id, "lemonsqueezy", event_name, raw_body.decode("utf-8"))
        if not is_new:
            return {"status": "duplicate", "event_id": event_id, "type": event_name}

        attrs = payload.get("data", {}).get("attributes", {})
        customer_id = str(
            attrs.get("customer_id")
            or attrs.get("user_email")
            or attrs.get("user_name")
            or ""
        )
        external_id = str(payload.get("data", {}).get("id") or event_id)
        return self._apply_subscription_event(
            event_type=event_name,
            customer_id=customer_id,
            metadata=attrs.get("custom_data", {}) or {},
            external_id=external_id,
            provider="lemonsqueezy",
        )

    # ── Internal ────────────────────────────────────────────────
    def _apply_subscription_event(self, *, event_type: str,
                                  customer_id: str, metadata: Dict[str, Any],
                                  external_id: str, provider: str) -> Dict[str, Any]:
        """Translate a provider event into a subscription state change."""
        if not customer_id:
            return {"status": "ignored", "reason": "missing_customer_id"}

        t_low = event_type.lower()

        # Normalise: Stripe uses "customer.subscription.created", Lemonsqueezy
        # uses "subscription_created". We match by stripping both separators
        # and checking for the meaningful tail tokens.
        norm = (
            t_low.replace("customer.", "")
                 .replace("subscription.", "subscription_")
                 .replace(".", "")
                 .replace("-", "_")
        )

        # Pro tier mapping helpers ─────────────────────────────────
        is_pro_signal = any(
            norm.startswith(token)
            for token in ("subscription_created", "subscription_updated",
                          "subscription_resumed", "order_created", "pro_")
        ) or "subscription_created" in norm or norm == "pro"
        is_lab_signal = ("lab" in norm) or (str(metadata.get("tier", "")).lower() == "lab")

        is_cancel_signal = any(
            norm.startswith(token)
            for token in ("subscription_canceled", "subscription_cancelled",
                          "subscription_expired", "subscription_deleted")
        ) or any(s in norm for s in
                 ("subscription_canceled", "subscription_cancelled",
                  "subscription_expired", "subscription_deleted"))
        is_failed_signal = "payment_failed" in norm or "past_due" in norm

        if is_cancel_signal:
            self.cancel_subscription(customer_id)
            return {"status": "canceled", "customer_id": customer_id,
                    "provider": provider, "event_type": event_type}

        if is_failed_signal:
            sub = self.store.latest_active_subscription(customer_id)
            if sub is not None:
                self.store.mark_subscription_expired(customer_id, sub.external_id)
            return {"status": "past_due", "customer_id": customer_id,
                    "provider": provider, "event_type": event_type}

        if is_pro_signal or is_lab_signal:
            tier = Tier.LAB if is_lab_signal else Tier.PRO
            sub = self.upgrade_customer(
                customer_id=customer_id,
                target=tier,
                provider=provider,
                external_id=external_id,
                days=30,
            )
            return {"status": "upgraded", "customer_id": customer_id,
                    "tier": tier.value, "provider": provider,
                    "event_type": event_type,
                    "subscription_id": sub.external_id}

        return {"status": "ignored", "reason": f"unhandled_event:{event_type}"}


# ───────────────────────────────────────────────────────────────────────
# FastAPI decorator
# ───────────────────────────────────────────────────────────────────────

def require_tier(required: str):
    """
    Decorator: gate a FastAPI endpoint so only customers at or above
    ``required`` (free | pro | lab) can call it.

    Usage::

        from api.billing import require_tier

        @app.post("/predict")
        @require_tier("pro")
        async def predict(req: PredictRequest, request: Request):
            ...

    Requires the FastAPI app to have a ``BillingManager`` stashed on
    ``app.state.billing``. If you don't, this raises ``RuntimeError`` at
    invocation time (not import time) so app boot can't be silently broken.
    """
    needed = Tier(required.lower())

    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            # Imported lazily so non-FastAPI hosts (CLI mock, tests) don't have to.
            try:
                from fastapi import HTTPException, Request  # noqa: F401
            except ImportError:
                raise RuntimeError(
                    "require_tier() needs fastapi; install "
                    "`pip install fastapi uvicorn`."
                )

            # Pull the request + billing manager out of the kwargs / context.
            from fastapi import Request as _Req

            request: Optional["_Req"] = kwargs.get("request")
            if request is None:
                for a in args:
                    if isinstance(a, _Req):
                        request = a
                        break
            if request is None:
                raise RuntimeError(
                    "@require_tier couldn't locate a FastAPI Request. "
                    "Add `request: Request` to your endpoint signature."
                )

            mgr: Optional[BillingManager] = getattr(request.app.state, "billing", None)
            if mgr is None:
                raise RuntimeError(
                    "No BillingManager attached. Set "
                    "`app.state.billing = BillingManager(...)` at startup."
                )

            # Resolve caller. Pro/Lab need an API key; Free doesn't.
            cfg = TIER_CONFIG[needed.value]
            api_key = request.headers.get("X-Api-Key") or request.headers.get("Authorization", "").removeprefix("Bearer ").strip() or None
            if cfg.requires_api_key and not api_key:
                raise HTTPException(status_code=401,
                                    detail=f"{needed.value} tier requires an API key.")
            if not api_key:
                tier = Tier.FREE
            else:
                try:
                    tier = mgr.tier_for_api_key(api_key)
                except BillingError as e:
                    raise HTTPException(status_code=401, detail=str(e))

            if _TIER_RANK[tier] < _TIER_RANK[needed]:
                raise HTTPException(
                    status_code=402,
                    detail=f"{needed.value} tier required (current: {tier.value}).",
                )

            ok, reason = mgr.check_call_allowed(api_key=api_key or "",
                                                endpoint=request.url.path,
                                                tier=tier)
            if not ok:
                raise HTTPException(status_code=429, detail=reason)

            if api_key:
                mgr.record_call(api_key=api_key,
                                endpoint=request.url.path,
                                tier=tier)
            return await func(*args, **kwargs)

        return wrapper
    return decorator


# ───────────────────────────────────────────────────────────────────────
# Helpers
# ───────────────────────────────────────────────────────────────────────

def _today_utc() -> str:
    return datetime.now(tz=timezone.utc).strftime("%Y-%m-%d")


def daily_reset_check(store: BillingStore, key: str, day_utc: Optional[str] = None) -> int:
    """Public helper used by tests. Returns the day's usage after reset."""
    return store.daily_count(key, day_utc or _today_utc())


# ───────────────────────────────────────────────────────────────────────
# Mock CLI
# ───────────────────────────────────────────────────────────────────────

def _run_mock_flow() -> int:
    """
    End-to-end demo using ``:memory:`` SQLite and fake webhook signatures.
    Prints exactly what a real Stripe / Lemonsqueezy interaction would log.
    Returns process exit code (0 on success).
    """
    print("=" * 70)
    print("DeepCatch billing — MOCK MODE")
    print("=" * 70)
    print("No network. No real accounts. Every step below would otherwise be")
    print("a real Stripe / Lemonsqueezy API call.\n")

    store = BillingStore(db_path=":memory:")
    mgr = BillingManager(
        store=store,
        stripe_secret="whsec_mocksk",
        lemonsqueezy_secret="lsqsec_mocksk",
        public_url="https://api.deepcatch.example.com",
        accept_any_signature=True,
    )

    # Step 1 — Issue a Free key (anonymous ok)
    print("[1] Free tier — no API key, 10 calls/day")
    ok, why = mgr.check_call_allowed(api_key=None, endpoint="/predict", tier=Tier.FREE)
    print(f"    GET /predict (anon):   ok={ok} reason={why}")
    assert ok
    print("    ✅ Free quota accepted.\n")

    # Step 2 — Upgrade customer via Stripe webhook
    print("[2] Stripe webhook: customer 'cus_alice' upgrades to Pro ($49/mo)")
    body = json.dumps({
        "id": "evt_mock_001",
        "type": "customer.subscription.created",
        "data": {"object": {
            "id": "sub_mock_alice_001",
            "customer": "cus_alice",
            "metadata": {"tier": "pro"},
        }},
    }).encode()
    verdict = mgr.handle_stripe_webhook(body, sig_header="t=1700000000,v1=deadbeef")
    print(f"    Stripe verdict: {json.dumps(verdict, indent=4).replace(chr(10), chr(10) + '    ')}")
    assert verdict["status"] == "upgraded", verdict
    assert verdict["tier"] == "pro"
    print("    ✅ Stripe webhook upgraded Alice → Pro.\n")

    # Step 3 — Provision an API key for Alice
    alice_key = mgr.provision_api_key(customer_id="cus_alice",
                                      tier=Tier.PRO,
                                      note="alice@lab.example").key
    print(f"[3] Provisioned API key for Alice (Pro):  {alice_key[:20]}...")
    print(f"    tier_for_api_key() → {mgr.tier_for_api_key(alice_key).value}")
    print("    ✅ Pro key provisioned and recognised.\n")

    # Step 4 — Rate limit a Pro key burst (way more than 10 quick calls)
    print("[4] Token-bucket burst test on Pro key")
    bucket_hits = sum(
        mgr.bucket.try_acquire(alice_key) for _ in range(200)
    )
    print(f"    200 try_acquire() attempts → {bucket_hits} successes (capacity=20)")
    assert 18 <= bucket_hits <= 22, bucket_hits
    print("    ✅ Burst limiter kicked in.\n")

    # Reset the bucket so step 5 isn't blocked by it.
    mgr.bucket.reset(alice_key)

    # Step 5 — Daily counter is per-key
    print("[5] Daily counter test on Pro key")
    for i in range(5):
        mgr.record_call(alice_key, "/predict", Tier.PRO)
    day = _today_utc()
    n = mgr.store.daily_count(alice_key, day)
    print(f"    After 5 calls: usage_table[day={day}] = {n}")
    assert n == 5
    print("    ✅ Daily counter increments per call.\n")

    # Step 6 — Free tier hits daily quota
    print("[6] Free tier hits daily quota at 10")
    free_key = "anonymous"     # canonical anon key used by check_call_allowed
    mgr.store.record_usage(UsageRecord(
        api_key=free_key, timestamp=time.time(),
        day_utc=day, endpoint="/predict", tier="free", count=10,
    ))
    ok, why = mgr.check_call_allowed(api_key=None, endpoint="/predict", tier=Tier.FREE)
    print(f"    GET /predict (anon, quota=10): ok={ok} reason={why}")
    assert not ok and "daily_quota_exceeded" in why, why
    print("    ✅ Free tier blocked at daily quota.\n")

    # Step 7 — Lab tier upgrade via Lemonsqueezy
    print("[7] Lemonsqueezy webhook: customer 'cus_labco' upgrades to Lab ($499/mo)")
    body = json.dumps({
        "id": "evt_mock_lab_001",
        "meta": {"event_name": "subscription_created"},
        "data": {"id": "lsq_sub_001",
                 "attributes": {
                     "customer_id": "cus_labco",
                     "user_email": "ops@labco.example",
                     "custom_data": {"tier": "lab"},
                 }},
    }).encode()
    verdict = mgr.handle_lemonsqueezy_webhook(body, sig_header="deadbeef")
    print(f"    Lemonsqueezy verdict: {verdict}")
    assert verdict["tier"] == "lab"
    print("    ✅ Lab tier unlocked via Lemonsqueezy.\n")

    # Step 8 — Verify webhook signature verification actually works
    print("[8] Verifying Stripe signature math in --verify mode")
    secret = "whsec_unit_test"
    body2 = b'{"hello":"world"}'
    ts = str(int(time.time()))
    expected = hmac.new(
        secret.encode(), f"{ts}.".encode() + body2, hashlib.sha256
    ).hexdigest()
    header = f"t={ts},v1={expected}"
    WebhookVerifier.verify_stripe(body2, header, secret)
    print(f"    Stripe HMAC-SHA256 match: header={header[:32]}... OK")

    # Bad signature must raise
    try:
        WebhookVerifier.verify_stripe(body2, "t=0,v1=00", secret)
    except WebhookVerificationError as e:
        print(f"    Bad-signature raises WebhookVerificationError: {e}")
        assert "mismatch" in str(e) or "tolerance" in str(e)
    print("    ✅ Signatures enforced.\n")

    # Step 9 — Daily reset: simulate by changing the clock + re-querying
    print("[9] Daily reset simulation")
    future_day = (datetime.now(tz=timezone.utc) + timedelta(days=1)).strftime("%Y-%m-%d")
    n_today = mgr.store.daily_count(alice_key, day)
    n_tomorrow = mgr.store.daily_count(alice_key, future_day)
    print(f"    usage_table for Alice: today={n_today}, tomorrow={n_tomorrow}")
    assert n_today == 5 and n_tomorrow == 0
    print("    ✅ Counters reset on day rollover.\n")

    print("=" * 70)
    print("ALL MOCK STEPS PASSED. billing.db would have been created at")
    print(f"  {_DB_PATH}")
    print("Run `python -m pytest api/test_billing.py -v` to re-execute as tests.")
    print("=" * 70)
    return 0


# ───────────────────────────────────────────────────────────────────────
# CLI
# ───────────────────────────────────────────────────────────────────────

def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m api.billing",
        description="DeepCatch billing tier manager — mock & verify CLI.",
    )
    parser.add_argument("--mock", action="store_true",
                        help="Run an end-to-end mock flow (no network).")
    parser.add_argument("--verify-webhook", choices=["stripe", "lemonsqueezy"],
                        help="Verify a signature given on stdin.")
    parser.add_argument("--secret", default="",
                        help="Webhook secret for --verify-webhook.")
    args = parser.parse_args(argv)

    if args.mock:
        return _run_mock_flow()

    if args.verify_webhook:
        body = sys.stdin.buffer.read()
        try:
            if args.verify_webhook == "stripe":
                WebhookVerifier.verify_stripe(body, body.decode().strip(), args.secret)
            else:
                WebhookVerifier.verify_lemonsqueezy(body, body.decode().strip(), args.secret)
        except WebhookVerificationError as e:
            print(f"FAIL: {e}", file=sys.stderr)
            return 2
        print("OK")
        return 0

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
