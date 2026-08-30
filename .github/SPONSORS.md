# Sponsor DeepCatch

> **Independent, reproducible, open-source cancer detection.** Built in Hong Kong. No institution. No paywall. No lock-in.

If DeepCatch helps your research or pipeline, please consider sponsoring its continued development.

## Why sponsor

DeepCatch is a solo project. Every line of code, every validation, every published result is built on my own time without institutional support. Sponsorship funds:

- **Compute** for the 627-sample cross-study validation (real FinaleDB WGS data)
- **Data licensing** where public access isn't enough
- **Time** to respond to issues, review PRs, and maintain reproducibility

100% of sponsor funds go to keeping DeepCatch open, reproducible, and at the state of the art.

## Tiers

All tiers are recurring monthly. Cancel anytime. GitHub takes 0% for maintainers; only Stripe processing (~2.9% + 30¢) is deducted.

### 🟢 Open Backer — **$5 / month**

A recurring thank-you. Your name listed in [`SPONSORS.md`](SPONSORS.md#current-backers) and in the monthly sponsor-update email.

Best for: students, indie researchers, people who find the project useful but don't need API access.

### 🔵 Pro Sponsor — **$49 / month**

Get a **single Pro API key** with:
- **10,000 calls/day** to the public DeepCatch API (vs. 10/day free tier)
- All 7 input modalities (5-mer motifs, FSD, WPS, coverage, IG, deconv, panel)
- Integrated-gradients explanations for every prediction
- Priority issue response (≤48h)

Use case: researchers running batch screening, lab pipelines needing higher rate limits, biotech consultants prototyping cfDNA assays.

**Get the key:** Subscribe at https://github.com/sponsors/rollroyces → "Manage subscription" → copy your Pro API key.

### 🟣 Lab Sponsor — **$499 / month**

Everything in Pro, plus:
- **Unlimited calls**
- **Raw motif dumps** (per-sample 4-mer frequency tables for offline analysis)
- **Model checkpoint access** (`gnn_pretrained.pt`, `deconv.pt`) on request
- Email support, custom integration help
- Quarterly research sync call (30 min)

Best for: biotech labs, contract research orgs, ML teams building on top of DeepCatch.

**Activation:** Email `deepcatch+sponsors@rollroyces.dev` after subscribing for the S3 download link and onboarding call.

## What you get (every tier)

- **Monthly sponsor-update email** — what shipped, what's validated, what's next. Written for humans, not marketing.
- **Early access** to major releases (v3 → v4 etc.)
- **Direct line** to the maintainer for questions, bug reports, or design feedback
- Your name on this page if you opt in

## What you DON'T get

- No private features. Everything in DeepCatch is open-source.
- No proprietary data access. The datasets used are all public (TCGA GDC, FinaleDB, GEO).
- No "guaranteed outcomes" or accuracy claims. Every result in the repo is 5-seed mean ± std with bootstrap CI — sponsors get the same honest numbers everyone else does.

## Current backers

<!--SPONSORS-LIST-START-->
_No backers yet — be the first._ 🟢
<!--SPONSORS-LIST-END-->

This list is auto-generated monthly from GitHub Sponsors data. It updates on the 1st of each month.

## Refund policy

If you sponsor and decide DeepCatch isn't useful, email me within 30 days for a full refund. No questions, no friction.

## Contact

- **Issues / questions:** https://github.com/rollroyces/deepcatch/issues
- **Email:** `deepcatch+sponsors@rollroyces.dev`
- **Twitter:** _(placeholder — update when you set up xurl)_

---

## Why I built this

I'm an independent researcher in Hong Kong. I don't have an institutional affiliation, an academic email, or lab access. I built DeepCatch because the problem is important and the open data is good enough.

The tradeoff: I can publish on bioRxiv (with ORCID + "Independent Researcher" affiliation), I can build on GitHub, I can use only public datasets. I can't access paywalled supplements, run wet-lab validations, or hire postdocs.

Sponsorship closes that gap. It buys time and compute, not infrastructure I don't have.

If you're a researcher, lab, or just someone who thinks this work matters — thank you.

— Royce (rollroyces)
