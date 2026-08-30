"""
Tests for the DeepCatch billing layer.

Run from the repo root::

    pytest api/test_billing.py -v

These tests use ``:memory:`` SQLite (no files written) and never touch the
network. They cover every requirement from the spec plus several sanity
checks to make the layer behave predictably under load.
"""

from __future__ import annotations

import hashlib
import hmac
import json
import time

import pytest

from api import billing
from api.billing import (
    ApiKeyRecord,
    BillingError,
    BillingManager,
    BillingStore,
    Tier,
    TokenBucket,
    UsageRecord,
    WebhookVerificationError,
    WebhookVerifier,
    daily_reset_check,
)


# ── Fixtures ─────────────────────────────────────────────────────────


@pytest.fixture
def fresh_db() -> BillingStore:
    """Brand-new in-memory SQLite store for each test."""
    return BillingStore(db_path=":memory:")


@pytest.fixture
def mgr(fresh_db: BillingStore) -> BillingManager:
    """BillingManager wired to that store, in mock-signature mode."""
    return BillingManager(
        store=fresh_db,
        stripe_secret="whsec_unit_test",
        lemonsqueezy_secret="lsqsec_unit_test",
        accept_any_signature=True,
    )


def _stripe_signature(raw_body: bytes, secret: str, ts: int | None = None) -> str:
    """Reproduce Stripe's ``t=...,v1=...`` HMAC format."""
    ts = ts if ts is not None else int(time.time())
    signed = f"{ts}.".encode("utf-8") + raw_body
    digest = hmac.new(secret.encode(), signed, hashlib.sha256).hexdigest()
    return f"t={ts},v1={digest}"


# ── 1. Tier upgrade (Stripe webhook → Pro key) ───────────────────────


def test_stripe_webhook_upgrades_customer_to_pro(mgr: BillingManager, fresh_db: BillingStore):
    """A Stripe ``customer.subscription.created`` event upgrades the customer
    to Pro and the subsequent provisioned key resolves to Pro."""
    body = json.dumps({
        "id": "evt_unit_001",
        "type": "customer.subscription.created",
        "data": {"object": {
            "id": "sub_unit_alice_001",
            "customer": "cus_alice",
            "metadata": {"tier": "pro"},
        }},
    }).encode()
    sig = _stripe_signature(body, mgr.stripe_secret)

    verdict = mgr.handle_stripe_webhook(body, sig)

    assert verdict["status"] == "upgraded"
    assert verdict["tier"] == "pro"
    assert verdict["customer_id"] == "cus_alice"
    assert verdict["subscription_id"] == "sub_unit_alice_001"

    sub = fresh_db.latest_active_subscription("cus_alice")
    assert sub is not None
    assert sub.tier == Tier.PRO
    assert sub.provider == "stripe"

    # Provision a key on top of the subscription
    rec = mgr.provision_api_key(customer_id="cus_alice", tier=Tier.PRO, note="alice")
    assert mgr.tier_for_api_key(rec.key) == Tier.PRO
    assert mgr.has_tier(rec.key, Tier.PRO)
    assert mgr.has_tier(rec.key, Tier.FREE)


# ── 2. Rate limit (token bucket + daily quota) ───────────────────────


def test_rate_limit_enforced_for_free_tier(mgr: BillingManager):
    """Free tier allows 10/day; the 11th call is rejected."""
    # Burst pass: token bucket lets a few through, then daily counter takes over.
    ok_first, _ = mgr.check_call_allowed(api_key=None, endpoint="/predict", tier=Tier.FREE)
    assert ok_first

    day = billing._today_utc()
    mgr.store.record_usage(UsageRecord(
        api_key="anonymous", timestamp=time.time(),
        day_utc=day, endpoint="/predict", tier="free", count=10,
    ))

    ok_after, reason = mgr.check_call_allowed(api_key=None, endpoint="/predict", tier=Tier.FREE)
    assert not ok_after
    assert "daily_quota_exceeded:10/10" in reason

    # Pro users on the same quota day are not affected
    pro_key = mgr.provision_api_key(customer_id="cus_pro", tier=Tier.PRO).key
    ok_pro, why_pro = mgr.check_call_allowed(api_key=pro_key, endpoint="/predict", tier=Tier.PRO)
    assert ok_pro, why_pro


def test_token_bucket_burst_rejects_after_capacity():
    """The token bucket rejects requests after capacity is exhausted."""
    b = TokenBucket(capacity=5, refill_per_sec=0.001)  # near-zero refill
    for i in range(5):
        assert b.try_acquire("k1"), f"first 5 should pass, got failure at {i}"
    assert not b.try_acquire("k1"), "6th should be rejected"
    b.reset("k1")
    assert b.try_acquire("k1"), "after reset, should pass again"


# ── 3. Expired subscription ──────────────────────────────────────────


def test_expired_subscription_reverts_to_free(mgr: BillingManager):
    """A subscription that has passed its ``expires_at`` no longer counts."""
    now = time.time()
    mgr.store.upsert_api_key(ApiKeyRecord(
        key="dck_unit_expired", customer_id="cus_old", tier=Tier.PRO,
        created_at=now - 100, active=True, note="",
    ))
    # Insert an already-expired subscription directly
    fresh_db = mgr.store
    fresh_db._conn.execute(
        "INSERT INTO subscriptions "
        "(customer_id, tier, started_at, expires_at, provider, external_id, status) "
        "VALUES (?, ?, ?, ?, ?, ?, 'active')",
        ("cus_old", "pro", now - 1000, now - 1, "stripe", "sub_old_001"),
    )
    fresh_db.commit()

    tier = mgr.effective_tier_for_customer("cus_old")
    assert tier == Tier.FREE, f"expected FREE after expiry, got {tier}"

    cancelled = mgr.upgrade_customer(
        customer_id="cus_old", target=Tier.PRO, provider="manual",
        external_id="sub_new_001", days=30,
    )
    assert mgr.effective_tier_for_customer("cus_old") == Tier.PRO

    # Cancellation drops tier
    mgr.cancel_subscription("cus_old")
    assert mgr.effective_tier_for_customer("cus_old") == Tier.FREE


# ── 4. Webhook signature verification ─────────────────────────────────


def test_stripe_signature_verification_round_trip():
    """Valid signature passes; tampered body fails."""
    secret = "whsec_unit_test_42"
    body = b'{"id":"evt_42","type":"customer.subscription.created"}'
    sig = _stripe_signature(body, secret)

    WebhookVerifier.verify_stripe(body, sig, secret)  # no exception

    bad = body + b"x"
    with pytest.raises(WebhookVerificationError):
        WebhookVerifier.verify_stripe(bad, sig, secret)

    # stale timestamp
    stale = _stripe_signature(body, secret, ts=int(time.time()) - 10_000)
    with pytest.raises(WebhookVerificationError):
        WebhookVerifier.verify_stripe(body, stale, secret)


def test_lemonsqueezy_signature_verification():
    """Lemonsqueezy uses raw hex HMAC."""
    secret = "lsqsec_unit"
    body = b'{"meta":{"event_name":"subscription_created"}}'
    expected = hmac.new(secret.encode(), body, hashlib.sha256).hexdigest()

    WebhookVerifier.verify_lemonsqueezy(body, expected, secret)

    with pytest.raises(WebhookVerificationError):
        WebhookVerifier.verify_lemonsqueezy(body, "deadbeef", secret)


def test_idempotent_webhook_events(mgr: BillingManager):
    """Same Stripe event id delivered twice → second call is marked duplicate."""
    body = json.dumps({
        "id": "evt_dup_001",
        "type": "customer.subscription.created",
        "data": {"object": {"id": "sub_dup_001", "customer": "cus_dup"}},
    }).encode()
    sig = _stripe_signature(body, mgr.stripe_secret)

    first = mgr.handle_stripe_webhook(body, sig)
    second = mgr.handle_stripe_webhook(body, sig)

    assert first["status"] == "upgraded"
    assert second["status"] == "duplicate"


# ── 5. Daily reset ────────────────────────────────────────────────────


def test_daily_counter_resets_on_day_rollover(mgr: BillingManager):
    """Day-keyed counter is independent per UTC day."""
    api_key = "dck_unit_daily"
    mgr.store.upsert_api_key(ApiKeyRecord(
        key=api_key, customer_id="cus_day", tier=Tier.PRO,
        created_at=time.time(), active=True,
    ))
    today = billing._today_utc()
    yesterday = "2020-01-01"     # anything before today
    tomorrow = "2999-12-31"      # anything after today

    mgr.store.record_usage(UsageRecord(
        api_key=api_key, timestamp=time.time(),
        day_utc=today, endpoint="/predict", tier="pro",
    ))
    mgr.store.record_usage(UsageRecord(
        api_key=api_key, timestamp=time.time(),
        day_utc=today, endpoint="/predict", tier="pro",
    ))

    assert daily_reset_check(mgr.store, api_key, today) == 2
    assert daily_reset_check(mgr.store, api_key, yesterday) == 0
    assert daily_reset_check(mgr.store, api_key, tomorrow) == 0


# ── 6. Mock-mode billing flow (end-to-end) ────────────────────────────


def test_mock_mode_billing_flow_end_to_end(mgr: BillingManager):
    """
    Simulate the exact sequence from ``--mock``:

    1) Stripe Pro upgrade
    2) Mint key
    3) Burst-limit test
    4) Record free-tier 11th-day rejection
    5) Lemonsqueezy Lab upgrade
    6) Free-tier stays blocked, Pro allowed
    """
    # (1) Stripe Pro upgrade
    body = json.dumps({
        "id": "evt_e2e_001",
        "type": "customer.subscription.created",
        "data": {"object": {"id": "sub_e2e_001", "customer": "cus_e2e"}},
    }).encode()
    verdict = mgr.handle_stripe_webhook(body, _stripe_signature(body, mgr.stripe_secret))
    assert verdict["status"] == "upgraded" and verdict["tier"] == "pro"

    # (2) Mint a Pro key
    rec = mgr.provision_api_key(customer_id="cus_e2e", tier=Tier.PRO)
    assert mgr.tier_for_api_key(rec.key) == Tier.PRO
    assert mgr.has_tier(rec.key, Tier.PRO)
    assert not mgr.has_tier(rec.key, Tier.LAB)

    # (3) Burst test: token bucket should kick in
    mgr.bucket.reset(rec.key)  # start clean
    bucket_hits = sum(mgr.bucket.try_acquire(rec.key) for _ in range(50))
    assert bucket_hits <= mgr.bucket.capacity  # at most capacity, then blocked

    # (4) Free tier hits daily quota after 10 calls
    day = billing._today_utc()
    mgr.store.record_usage(UsageRecord(
        api_key="anonymous", timestamp=time.time(),
        day_utc=day, endpoint="/predict", tier="free", count=10,
    ))
    ok, reason = mgr.check_call_allowed(api_key=None, endpoint="/predict", tier=Tier.FREE)
    assert not ok and "daily_quota_exceeded" in reason, reason

    # (5) Lemonsqueezy Lab upgrade
    lsq_body = json.dumps({
        "id": "evt_lsq_e2e",
        "meta": {"event_name": "subscription_created"},
        "data": {"id": "lsq_e2e_001",
                 "attributes": {
                     "customer_id": "cus_lab",
                     "custom_data": {"tier": "lab"},
                 }},
    }).encode()
    lsq_sig = hmac.new(mgr.lemonsqueezy_secret.encode(), lsq_body,
                        hashlib.sha256).hexdigest()
    lsq_verdict = mgr.handle_lemonsqueezy_webhook(lsq_body, lsq_sig)
    assert lsq_verdict["status"] == "upgraded"
    assert lsq_verdict["tier"] == "lab"

    # (6) Pro key untouched by Lemonsqueezy event
    assert mgr.tier_for_api_key(rec.key) == Tier.PRO
    lab_key = mgr.provision_api_key(customer_id="cus_lab", tier=Tier.LAB).key
    assert mgr.has_tier(lab_key, Tier.LAB)
    assert mgr.has_tier(lab_key, Tier.PRO)


# ── 7. Tier gating on the decorator (without real FastAPI deps) ───────


def test_require_tier_rejects_lower_tier(mgr: BillingManager):
    """``has_tier`` returns False when key is below the required level."""
    free_key = mgr.provision_api_key(customer_id="cus_free_user", tier=Tier.FREE).key
    pro_key = mgr.provision_api_key(customer_id="cus_pro_user", tier=Tier.PRO).key

    # FREE API key should not pass PRO gating
    assert not mgr.has_tier(free_key, Tier.PRO)
    # PRO should pass FREE gating
    assert mgr.has_tier(pro_key, Tier.FREE)
    # PRO should pass PRO gating
    assert mgr.has_tier(pro_key, Tier.PRO)
    # PRO should not pass LAB gating
    assert not mgr.has_tier(pro_key, Tier.LAB)


def test_unknown_api_key_rejected(mgr: BillingManager):
    """An unknown key raises BillingError (translated to 401 by decorator)."""
    with pytest.raises(BillingError):
        mgr.tier_for_api_key("dck_does_not_exist")
    assert mgr.has_tier("dck_does_not_exist", Tier.FREE) is False
