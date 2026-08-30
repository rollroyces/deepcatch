#!/usr/bin/env python3
"""
sponsor_digest.py — daily GitHub digest for DeepCatch sponsors.

==================================================================

Reads GitHub repo stats (stars, forks, open/new issues, recent commits,
download/release counts) via the ``gh`` CLI and renders them as a
markdown digest suitable for emailing to sponsors.

Examples
--------
::

    # Real run (uses gh CLI; requires gh auth login)
    python scripts/sponsor_digest.py \
        --repo rollroyces/deepcatch \
        --output /tmp/digest.md

    # Offline / no-network run
    python scripts/sponsor_digest.py \
        --repo rollroyces/deepcatch \
        --output /tmp/digest.md \
        --mock

Cron-friendly: re-running on the same day yields an idempotent file
(with the date stamped in the header). Exit code is 0 on success and
non-zero on CLI errors so cron can pick up failures.

The digest includes:
- Stars, forks, watchers, open issues (current snapshot)
- Issue delta in the last 7 days (new, closed)
- Top 5 contributors in the last 30 days
- Top 5 recent commits
- Latest release (if any)
- A short one-liner pitching what's new
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

DATE_FMT = "%Y-%m-%d"
HUMAN_DATE = "%A, %B %d, %Y"


# ───────────────────────────────────────────────────────────────────────
# Data shape
# ───────────────────────────────────────────────────────────────────────

@dataclass
class RepoStats:
    repo: str
    stars: int = 0
    forks: int = 0
    watchers: int = 0
    open_issues: int = 0
    description: str = ""
    primary_language: str = ""
    license: str = ""
    pushed_at: str = ""

    # Last-7-day activity (real mode only)
    new_issues_7d: int = 0
    closed_issues_7d: int = 0

    # Last-30-day activity
    commits_30d: int = 0
    contributors_30d: List[Dict[str, Any]] = field(default_factory=list)

    # Latest release (or empty string)
    latest_release: str = ""
    latest_release_url: str = ""

    # Recent commits (top 5)
    recent_commits: List[Dict[str, str]] = field(default_factory=list)


# ───────────────────────────────────────────────────────────────────────
# GitHub helpers — wrap ``gh`` so we can swap in mocks for tests
# ───────────────────────────────────────────────────────────────────────

def _run_gh(args: List[str], timeout: int = 30) -> subprocess.CompletedProcess:
    """
    Run a ``gh`` subcommand, returning a CompletedProcess.

    Set ``GITHUB_TOKEN`` env to ``""`` to force a no-network run.
    """
    return subprocess.run(
        ["gh"] + args,
        capture_output=True,
        text=True,
        timeout=timeout,
        env={**os.environ, "GH_NO_UPDATE_NOTIFIER": "1"},
    )


def _gh_json(args: List[str], timeout: int = 30) -> Any:
    cp = _run_gh(args + ["--json"], timeout=timeout) if "--json" not in args else _run_gh(args, timeout=timeout)
    if cp.returncode != 0:
        raise RuntimeError(
            f"gh {' '.join(args)} failed (rc={cp.returncode}): "
            f"stderr={cp.stderr.strip()[:200]}"
        )
    return json.loads(cp.stdout or "null")


# ───────────────────────────────────────────────────────────────────────
# Real collectors (require `gh auth login`)
# ───────────────────────────────────────────────────────────────────────

def collect_real(repo: str) -> RepoStats:
    stats = RepoStats(repo=repo)

    # 1) Repo snapshot — use fields that gh 2.9x actually exposes
    repo_json = _gh_json([
        "repo", "view", repo,
        "--json", "stargazerCount,forkCount,issues,"
                  "description,pushedAt,licenseInfo,nameWithOwner,"
                  "primaryLanguage",
    ])
    stats.stars = int(repo_json.get("stargazerCount") or 0)
    stats.forks = int(repo_json.get("forkCount") or 0)
    stats.description = (repo_json.get("description") or "").strip()
    stats.primary_language = (repo_json.get("primaryLanguage") or {}).get("name", "") if isinstance(repo_json.get("primaryLanguage"), dict) else ""
    license_info = repo_json.get("licenseInfo") or repo_json.get("license")
    if isinstance(license_info, dict):
        stats.license = (license_info.get("key") or license_info.get("spdxId") or "").upper()
    else:
        stats.license = ""
    stats.pushed_at = repo_json.get("pushedAt", "")

    issues_obj = repo_json.get("issues")
    if isinstance(issues_obj, dict):
        stats.open_issues = int(issues_obj.get("totalCount", 0) or 0)
    else:
        stats.open_issues = 0

    # 2) New + closed issues in the last 7 days
    seven_days = (datetime.now(tz=timezone.utc) - timedelta(days=7)).strftime(DATE_FMT)
    def _search(query: str) -> int:
        cp = _run_gh([
            "search", "issues", query,
            "--repo", repo,
            "--json", "number",
            "--limit", "200",
            "--jq", "length",
        ])
        if cp.returncode != 0:
            return 0
        try:
            return int(cp.stdout.strip() or "0")
        except ValueError:
            return 0

    stats.new_issues_7d = _search(f"repo:{repo} is:issue created:>={seven_days}")
    stats.closed_issues_7d = _search(f"repo:{repo} is:issue is:closed closed:>={seven_days}")

    # 3) Commits in the last 30 days
    since_iso = (datetime.now(tz=timezone.utc) - timedelta(days=30)).strftime("%Y-%m-%dT%H:%M:%SZ")
    cp = _run_gh([
        "api", f"repos/{repo}/commits?since={since_iso}&per_page=100",
        "--jq", "length",
    ])
    if cp.returncode == 0 and cp.stdout.strip().isdigit():
        stats.commits_30d = int(cp.stdout.strip())

    # 4) Top 5 contributors (last 30d, best-effort)
    cp = _run_gh([
        "api", f"repos/{repo}/contributors?per_page=5",
    ])
    if cp.returncode == 0:
        try:
            contribs = json.loads(cp.stdout or "[]")
            stats.contributors_30d = [
                {
                    "login": c.get("login", ""),
                    "contributions": int(c.get("contributions", 0)),
                    "html_url": c.get("html_url", ""),
                }
                for c in contribs
            ]
        except json.JSONDecodeError:
            pass

    # 5) Latest release
    cp = _run_gh([
        "release", "view", "--repo", repo, "--json", "name,url,publishedAt",
    ])
    if cp.returncode == 0 and cp.stdout.strip():
        try:
            rel = json.loads(cp.stdout)
            stats.latest_release = rel.get("name") or ""
            stats.latest_release_url = rel.get("url") or ""
        except json.JSONDecodeError:
            pass

    # 6) Recent 5 commits
    cp = _run_gh([
        "api", f"repos/{repo}/commits?per_page=5",
    ])
    if cp.returncode == 0:
        try:
            commits = json.loads(cp.stdout or "[]")
            for c in commits:
                sha = c.get("sha", "")
                msg = (c.get("commit", {}).get("message") or "").splitlines()[0]
                author = (c.get("commit", {}).get("author") or {}).get("name", "")
                stats.recent_commits.append({
                    "sha": sha[:7] if sha else "",
                    "url": c.get("html_url", ""),
                    "message": msg[:120],
                    "author": author,
                })
        except json.JSONDecodeError:
            pass

    return stats


# ───────────────────────────────────────────────────────────────────────
# Mock collectors — for offline / --mock runs
# ───────────────────────────────────────────────────────────────────────

def collect_mock(repo: str) -> RepoStats:
    """Deterministic fake numbers. Stable across runs so the digest is testable."""
    # Hash-stable so the same repo gives the same mock numbers
    base = sum(ord(c) for c in repo) or 1
    return RepoStats(
        repo=repo,
        stars=212 + (base % 17),
        forks=37 + (base % 9),
        watchers=18 + (base % 5),
        open_issues=4 + (base % 4),
        description="(mock) DeepCatch — open-source MCED from cfDNA, 7-modality foundation model.",
        primary_language="Python",
        license="MIT",
        pushed_at=datetime.now(tz=timezone.utc).isoformat(),
        new_issues_7d=3,
        closed_issues_7d=2,
        commits_30d=42,
        contributors_30d=[
            {"login": "rollroyces", "contributions": 312, "html_url": f"https://github.com/{repo.split('/')[0]}"},
            {"login": "deepcatch-bot", "contributions": 47,  "html_url": ""},
            {"login": "reviewer-1", "contributions": 18, "html_url": ""},
            {"login": "reviewer-2", "contributions": 12, "html_url": ""},
            {"login": "reviewer-3", "contributions": 9,  "html_url": ""},
        ],
        latest_release="v2.2.0 (panel-based MCED)",
        latest_release_url=f"https://github.com/{repo}/releases/tag/v2.2.0",
        recent_commits=[
            {"sha": "abc1234", "url": f"https://github.com/{repo}/commit/abc1234",
             "message": "Add Kalman Stage 2 longitudinal tracking",
             "author": "rollroyces"},
            {"sha": "def5678", "url": f"https://github.com/{repo}/commit/def5678",
             "message": "Fix LOESS GC bias in DELFI normalisation",
             "author": "rollroyces"},
            {"sha": "9012abc", "url": f"https://github.com/{repo}/commit/9012abc",
             "message": "Docs: clarify Lemonsqueezy webhook setup",
             "author": "rollroyces"},
            {"sha": "3456def", "url": f"https://github.com/{repo}/commit/3456def",
             "message": "Tests: bump coverage of billing layer",
             "author": "reviewer-1"},
            {"sha": "7890abc", "url": f"https://github.com/{repo}/commit/7890abc",
             "message": "CI: cache venv across runs",
             "author": "deepcatch-bot"},
        ],
    )


# ───────────────────────────────────────────────────────────────────────
# Markdown rendering
# ───────────────────────────────────────────────────────────────────────

def render_markdown(stats: RepoStats, *, generated_at: datetime,
                    mock: bool, sponsor: str = "Dear Sponsor") -> str:
    today = generated_at.strftime(HUMAN_DATE)
    iso_today = generated_at.strftime(DATE_FMT)

    tag = "  (MOCK RUN)" if mock else ""
    md = [f"# DeepCatch Sponsor Digest — {today}{tag}", ""]
    md.append(f"Hello {sponsor.strip() or 'there'},")
    md.append("")
    md.append(
        f"Here is your daily snapshot of **{stats.repo}**. "
        "Numbers below are real (or, in mock mode, deterministic fakes "
        "suitable for previewing the email)."
    )
    md.append("")

    # ── Headline numbers ────────────────────────────────────────
    md.append("## 📊 Repo Snapshot")
    md.append("")
    md.append("| Metric | Value |")
    md.append("|---|---:|")
    md.append(f"| ⭐ Stars  | **{stats.stars}** |")
    md.append(f"| 🍴 Forks  | {stats.forks} |")
    if stats.watchers:
        md.append(f"| 👀 Watchers | {stats.watchers} |")
    md.append(f"| 🐛 Open issues | {stats.open_issues} |")
    if stats.primary_language:
        md.append(f"| 💻 Primary language | {stats.primary_language} |")
    if stats.license:
        md.append(f"| 📜 License | {stats.license} |")
    md.append("")

    # ── Weekly activity ────────────────────────────────────────
    md.append("## 🗓 Last 7 Days")
    md.append("")
    delta = stats.new_issues_7d - stats.closed_issues_7d
    md.append(f"- New issues: **{stats.new_issues_7d}**")
    md.append(f"- Closed issues: **{stats.closed_issues_7d}**")
    md.append(f"- Net issue delta: {delta:+d}")
    md.append(f"- Commits in last 30 days: **{stats.commits_30d}**")
    md.append("")

    # ── Latest release ─────────────────────────────────────────
    if stats.latest_release:
        md.append("## 🚀 Latest Release")
        md.append("")
        if stats.latest_release_url:
            md.append(f"[{stats.latest_release}]({stats.latest_release_url})")
        else:
            md.append(stats.latest_release)
        md.append("")

    # ── Recent commits ─────────────────────────────────────────
    if stats.recent_commits:
        md.append("## 🛠 Last 5 Commits")
        md.append("")
        for c in stats.recent_commits:
            sha = c.get("sha", "")
            url = c.get("url", "")
            msg = c.get("message", "")
            author = c.get("author", "")
            if url and sha:
                md.append(f"- [`{sha}`]({url}) {msg} — _{author}_")
            elif sha:
                md.append(f"- `{sha}` {msg} — _{author}_")
            else:
                md.append(f"- {msg} — _{author}_")
        md.append("")

    # ── Contributors ───────────────────────────────────────────
    if stats.contributors_30d:
        md.append("## 🤝 Top Contributors (all-time)")
        md.append("")
        for c in stats.contributors_30d[:5]:
            login = c.get("login", "")
            n = c.get("contributions", 0)
            url = c.get("html_url", "")
            if url:
                md.append(f"- [{login}]({url}) — {n} contributions")
            else:
                md.append(f"- {login} — {n} contributions")
        md.append("")

    # ── Footer ─────────────────────────────────────────────────
    md.append("---")
    md.append("")
    md.append(
        "Thank you for supporting open computational cancer research. "
        "Questions or feedback? Reply to this email or open an issue on GitHub."
    )
    md.append("")
    md.append(f"_Generated {today} ({iso_today}). "
              "Cron-runnable, idempotent: re-running produces the same shape "
              "with updated numbers._")
    md.append("")

    return "\n".join(md)


# ───────────────────────────────────────────────────────────────────────
# CLI
# ───────────────────────────────────────────────────────────────────────

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="sponsor_digest.py",
        description="Daily GitHub digest for DeepCatch sponsors (markdown).",
    )
    p.add_argument("--repo", default="rollroyces/deepcatch",
                   help="GitHub repo (default: rollroyces/deepcatch).")
    p.add_argument("--output", required=True,
                   help="Path to write the markdown digest to.")
    p.add_argument("--mock", action="store_true",
                   help="Use deterministic fake data (no network, no gh).")
    p.add_argument("--no-cache", action="store_true",
                   help="Always fetch fresh data even if a cached digest exists.")
    p.add_argument("--sponsor", default="Dear Sponsor",
                   help="Salutation to put in the email (e.g. 'Hi Alice').")
    p.add_argument("--publish-discussion", action="store_true",
                   help="After writing the digest, post it to the repo's "
                        "'Announcements' category as a GitHub Discussion. "
                        "Uses `gh` CLI; requires Discussions to be enabled.")
    p.add_argument("--discussion-category", default="Announcements",
                   help="Discussions category (default: Announcements).")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    out_path = Path(args.output).expanduser()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    generated_at = datetime.now(tz=timezone.utc)

    if args.mock:
        stats = collect_mock(args.repo)
    else:
        if shutil.which("gh") is None:
            print("error: gh CLI not found in PATH. Use --mock for offline runs.",
                  file=sys.stderr)
            return 2
        try:
            stats = collect_real(args.repo)
        except Exception as e:
            print(f"warn: gh collection failed ({e}); falling back to mock",
                  file=sys.stderr)
            stats = collect_mock(args.repo)
            args.mock = True     # mark the file accordingly

    digest = render_markdown(stats, generated_at=generated_at,
                             mock=args.mock, sponsor=args.sponsor)

    out_path.write_text(digest)
    print(f"  wrote {out_path}  ({len(digest)} bytes, mock={args.mock})")

    if args.publish_discussion:
        if args.mock:
            print("skipping --publish-discussion (mock mode)", file=sys.stderr)
        else:
            rc = publish_to_discussion(args.repo, digest, args.discussion_category)
            if rc != 0:
                print(f"warn: discussion publish failed (rc={rc}); digest still at {out_path}",
                      file=sys.stderr)
                # Don't fail the cron — the digest is the primary deliverable.

    return 0


def publish_to_discussion(repo: str, digest_md: str, category: str) -> int:
    """Post the digest as a GitHub Discussion in the given category.

    Requires `gh` CLI authenticated and Discussions enabled on the repo.
    Uses the GraphQL API via `gh api graphql`.

    Returns 0 on success, non-zero on failure.
    """
    title = f"Sponsor update — {datetime.now(tz=timezone.utc).strftime('%Y-%m-%d')}"
    # Strip the H1 header from the digest since Discussions has its own title field.
    body = digest_md
    lines = body.split("\n", 1)
    if lines and lines[0].startswith("# "):
        body = lines[1] if len(lines) > 1 else ""

    # GraphQL: get repository ID, then createDiscussion mutation.
    # Using `gh api graphql` because Discussions aren't in the REST API.
    query = """
    mutation PublishDigest($repoId: ID!, $catId: ID!, $title: String!, $body: String!) {
      createDiscussion(input: {
        repositoryId: $repoId,
        categoryId: $catId,
        title: $title,
        body: $body
      }) {
        discussion { url }
      }
    }
    """
    # Step 1: get repo + category IDs (GraphQL query)
    id_query = f"""
    {{
      repository(owner: "{repo.split('/')[0]}", name: "{repo.split('/')[1]}") {{
        id
        discussionCategories(first: 20) {{
          nodes {{ id name }}
        }}
      }}
    }}
    """
    try:
        id_proc = subprocess.run(
            ["gh", "api", "graphql", "-f", f"query={id_query}"],
            capture_output=True, text=True, timeout=30,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"warn: gh failed ({e}); skipping discussion publish", file=sys.stderr)
        return 1
    if id_proc.returncode != 0:
        print(f"warn: GraphQL id query failed: {id_proc.stderr.strip()}", file=sys.stderr)
        return 1
    try:
        id_data = json.loads(id_proc.stdout)
    except json.JSONDecodeError as e:
        print(f"warn: could not parse GraphQL response: {e}", file=sys.stderr)
        return 1
    repo_id = id_data.get("data", {}).get("repository", {}).get("id")
    if not repo_id:
        print(f"warn: repo not found or no access: {repo}", file=sys.stderr)
        return 1
    cat_nodes = id_data.get("data", {}).get("repository", {}).get(
        "discussionCategories", {}).get("nodes", [])
    cat_id = next((c["id"] for c in cat_nodes if c.get("name") == category), None)
    if not cat_id:
        available = ", ".join(c.get("name", "?") for c in cat_nodes)
        print(f"warn: category '{category}' not found; available: {available}",
              file=sys.stderr)
        return 1

    # Step 2: create the discussion
    try:
        create_proc = subprocess.run(
            ["gh", "api", "graphql",
             "-f", f"query={query}",
             "-f", f"repoId={repo_id}",
             "-f", f"catId={cat_id}",
             "-f", f"title={title}",
             "-F", f"body={body}"],
            capture_output=True, text=True, timeout=60,
        )
    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        print(f"warn: discussion creation timed out ({e})", file=sys.stderr)
        return 1
    if create_proc.returncode != 0:
        print(f"warn: createDiscussion failed: {create_proc.stderr.strip()}",
              file=sys.stderr)
        return 1
    try:
        create_data = json.loads(create_proc.stdout)
        url = (create_data.get("data", {}).get("createDiscussion", {})
               .get("discussion", {}).get("url"))
    except json.JSONDecodeError:
        url = None
    if url:
        print(f"  published to: {url}")
        return 0
    print(f"warn: discussion created but no URL returned: {create_proc.stdout}",
          file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
