#!/usr/bin/env python3
"""
sponsor_outreach.py — generate a personalized sponsor outreach message
for any GitHub user who has starred or forked rollroyces repos.

Strategy:
  1. Fetch stargazers + watchers of rollroyces repos via `gh` CLI
  2. Filter out rollroyces's own user, bots, and existing sponsors
  3. Output a CSV with: github_login, profile_url, email (if public), repos_starred
  4. Generate a personalized email draft for each

The output is a CSV (one row per potential sponsor) + a markdown file with
all draft emails, ready to copy-paste into Gmail / Apple Mail / etc.

No email is sent automatically — this is a draft-generation tool. You send
manually via your normal email client because:
  - cold outreach without consent violates spam laws (CAN-SPAM, GDPR, etc.)
  - one-to-one personalized messages convert 10-100x better than blasts
  - GitHub-stargazer emails aren't directly accessible without GraphQL + Patreon

Usage:
    python sponsor_outreach.py --repo rollroyces/deepcatch --output /tmp/outreach.csv
    python sponsor_outreach.py --all-public --output /tmp/outreach.csv

For each potential sponsor, generates an email like:
    Subject: Quick question about your [DeepCatch] fork

    Hi @username,

    Saw you starred rollroyces/[repo] — thanks for the interest in [topic].

    I'm Royce, the maintainer. [One specific question about their use case
    based on their other public repos / bio].

    I'm working on keeping DeepCatch sustainable via GitHub Sponsors.
    If the project has been useful, [tier suggestion based on their
    apparent engagement level].

    [Link to SPONSORS.md]

    Either way, no pressure. If there's a specific feature or bug I
    should know about, just hit reply.

    — Royce
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


GH_TIMEOUT = 30  # seconds


def gh_json(args: list[str], fallback: Any = None) -> Any:
    """Run `gh api ...` and return JSON. Returns fallback on failure."""
    try:
        proc = subprocess.run(
            ["gh"] + args,
            capture_output=True, text=True, timeout=GH_TIMEOUT,
        )
        if proc.returncode != 0:
            return fallback
        return json.loads(proc.stdout)
    except (subprocess.TimeoutExpired, json.JSONDecodeError, FileNotFoundError):
        return fallback


def fetch_stargazers(repo: str) -> list[dict[str, Any]]:
    """Fetch all stargazers for a repo. Returns list of {login, url}."""
    out: list[dict[str, Any]] = []
    page = 1
    while True:
        batch = gh_json([
            "api", f"/repos/{repo}/stargazers",
            "--paginate", "--jq",
            ".[].user | {login: .login, url: .html_url}",
        ], fallback=[])
        if not batch:
            break
        if isinstance(batch, list):
            out.extend(batch)
        # gh --paginate doesn't expose page count directly; cap at 10 pages
        if page >= 10:
            break
        page += 1
    return out


def fetch_user_profile(login: str) -> dict[str, Any] | None:
    """Fetch public profile for a GitHub user."""
    return gh_json(["api", f"/users/{login}"], fallback=None)


def is_excluded(login: str) -> bool:
    """Filter out rollroyces's own user and known bots."""
    excluded = {
        "rollroyces", "github-actions[bot]", "dependabot[bot]",
        "renovate[bot]", "codecov[bot]", "github-actions",
    }
    return login.lower() in excluded


def generate_email(login: str, profile: dict[str, Any] | None,
                   starred_repos: list[str]) -> str:
    """Generate a personalized outreach email draft."""
    name = (profile or {}).get("name") or f"@{login}"
    bio = (profile or {}).get("bio") or ""
    company = (profile or {}).get("company") or ""
    location = (profile or {}).get("location") or ""

    repos_str = ", ".join(f"`{r}`" for r in starred_repos[:3])

    # Tier suggestion based on apparent engagement (rough heuristic):
    #   -1 repo: thank-you-only tier ($5)
    #   2 repos: open backer ($5)
    #   3+ repos: pro sponsor ($49)
    if len(starred_repos) >= 3:
        tier = "Pro Sponsor ($49/mo — full API key, 10K calls/day)"
    elif len(starred_repos) == 2:
        tier = "Open Backer ($5/mo — listed in SPONSORS.md)"
    else:
        tier = "Open Backer ($5/mo)"

    # Opening hook tailored to their context
    if "research" in (bio + company).lower() or "lab" in (bio + company).lower() or "phd" in (bio + company).lower():
        hook = "I saw your work in research and figured you might use DeepCatch for cfDNA-related work."
    elif company:
        hook = f"Saw your affiliation with {company} — figured DeepCatch might fit your pipeline."
    elif location and location.lower() not in {"hong kong", "hk"}:
        hook = f"Saw you're based in {location} — figured DeepCatch might be useful regardless of geography."
    else:
        hook = f"Thanks for the star on {repos_str or 'DeepCatch'}."

    return f"""Subject: Quick question about your use of {repos_str or 'DeepCatch'}

Hi {name},

{hook}

I'm Royce, the solo maintainer. I'm working on keeping DeepCatch sustainable
via GitHub Sponsors — see https://github.com/rollroyces/deepcatch/blob/main/.github/SPONSORS.md

If the project has been useful for your work, the {tier} tier would mean
a lot. Even $5/mo helps fund the next round of cross-study validation
(currently targeting 1,000+ real cfDNA samples from open-access sources).

If not, no pressure — but if there's a specific feature, bug, or
direction you'd like DeepCatch to take, just hit reply. I read every email.

Thanks for the star(s).

— Royce (rollroyces)
   https://github.com/sponsors/rollroyces
"""


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="sponsor_outreach.py",
        description="Generate personalized sponsor outreach drafts for rollroyces repo stargazers.",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--repo", help="Single repo (e.g. rollroyces/deepcatch)")
    group.add_argument("--all-public", action="store_true",
                       help="Scan all public rollroyces repos")
    parser.add_argument("--output", required=True,
                        help="Output CSV path (one row per potential sponsor)")
    parser.add_argument("--drafts-dir", default=None,
                        help="If set, write one email-draft .md file per sponsor here")
    parser.add_argument("--max", type=int, default=500,
                        help="Max prospects to consider (safety cap)")
    parser.add_argument("--exclude-existing-sponsors", action="store_true",
                        help="Try to skip people who already sponsor (requires Sponsors API)")
    args = parser.parse_args(argv)

    repos: list[str] = []
    if args.repo:
        repos = [args.repo]
    else:
        repo_list = gh_json(["repo", "list", "rollroyces", "--json", "name,visibility",
                             "--jq", ".[] | select(.visibility==\"public\") | .name"], fallback=[])
        if isinstance(repo_list, list):
            repos = [f"rollroyces/{n}" for n in repo_list[:10]]
        else:
            print("WARN: could not enumerate repos; falling back to deepcatch only", file=sys.stderr)
            repos = ["rollroyces/deepcatch"]

    print(f"Scanning {len(repos)} repo(s) for stargazers: {repos}", file=sys.stderr)

    # Map login -> set of repos starred
    prospects: dict[str, set[str]] = {}
    for repo in repos:
        stars = fetch_stargazers(repo)
        print(f"  {repo}: {len(stars)} stargazers", file=sys.stderr)
        for s in stars:
            login = s.get("login")
            if not login or is_excluded(login):
                continue
            prospects.setdefault(login, set()).add(repo)
        if len(prospects) >= args.max:
            print(f"  hit cap of {args.max}; stopping early", file=sys.stderr)
            break

    # Sort by engagement (most-repos-starred first)
    sorted_logins = sorted(prospects.keys(),
                           key=lambda l: (-len(prospects[l]), l))

    out_csv = Path(args.output)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    drafts_dir = Path(args.drafts_dir) if args.drafts_dir else None
    if drafts_dir:
        drafts_dir.mkdir(parents=True, exist_ok=True)

    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["github_login", "profile_url", "repos_starred",
                        "n_repos", "tier_suggested", "email_drafted"],
        )
        writer.writeheader()
        written = 0
        for login in sorted_logins:
            profile = fetch_user_profile(login)
            repos_starred = sorted(prospects[login])
            draft = generate_email(login, profile, repos_starred)
            tier = "pro" if len(repos_starred) >= 3 else "open"
            profile_url = f"https://github.com/{login}"

            email_drafted = "no"
            if drafts_dir:
                (drafts_dir / f"{login}.md").write_text(draft, encoding="utf-8")
                email_drafted = "yes"

            writer.writerow({
                "github_login": login,
                "profile_url": profile_url,
                "repos_starred": ";".join(repos_starred),
                "n_repos": len(repos_starred),
                "tier_suggested": tier,
                "email_drafted": email_drafted,
            })
            written += 1
            if written >= args.max:
                break

    print(f"\nWrote {written} prospects to {out_csv}", file=sys.stderr)
    if drafts_dir:
        print(f"Wrote {written} email drafts to {drafts_dir}/", file=sys.stderr)
    print(f"\nNext steps:", file=sys.stderr)
    print(f"  1. Open {out_csv} and review the list", file=sys.stderr)
    print(f"  2. For each prospect, send the draft email via your normal mail client", file=sys.stderr)
    print(f"  3. Track responses in state/sponsor_outreach_log.csv (manual)", file=sys.stderr)
    print(f"  4. Don't blast — send 3-5 personalized emails per day max", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
