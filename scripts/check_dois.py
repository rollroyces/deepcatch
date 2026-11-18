#!/usr/bin/env python3
"""Check every DOI in the tracked tree against the registry that issued it.

Why this exists
---------------
A DOI that does not resolve is a dead link in a document a reviewer reads.
Three of them were found in this repository's bibliographies, and one had two
digits transposed (`...-21394-w` for `...-24994-w`) -- a typo that reads
like a DOI and resolves to nothing.

Two registries, because two kinds of DOI are in use here:

* Crossref indexes the journal literature (`api.crossref.org/works/<doi>`).
  A 404 from it means the identifier is not registered there.
* DataCite indexes Zenodo, which is where this repository's own releases are
  archived (`api.datacite.org/dois/<doi>`). A Zenodo DOI is *not* in
  Crossref, so checking it there would report a false failure.

Usage
-----
    python3 scripts/check_dois.py            # check and report
    python3 scripts/check_dois.py --quiet    # only the failures

Exit status is 1 when an identifier is not registered, and 2 when the
registry could not be reached (an offline machine is not a failure).
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent

# A DOI in running text: the trailing characters of a markdown link, a
# BibTeX field or an escaped JSON quote are not part of it. Parentheses
# are allowed because Elsevier's DOIs contain them.
_DOI = re.compile(
    r"10\.\d{4,9}/[^\s\"'<>,;\\\(\)\]]*"
    r"(?:\([^\s\"'<>,;\\\(\)\]]*\)[^\s\"'<>,;\\\(\)\]]*)*")


def _trim_unbalanced(text: str) -> str:
    """Drop trailing parentheses a sentence or a markdown link added.

    Elsevier's DOIs contain parentheses, so they cannot simply be excluded
    from the character class -- but a DOI that ends a markdown link picks up
    the link's closing parenthesis, which is not part of it.
    """
    while text.endswith(")") and text.count(")") > text.count("("):
        text = text[:-1]
    return text

# Files that are not text.
_BINARY_SUFFIXES = {
    ".pdf", ".png", ".jpg", ".webm", ".mp4", ".gif", ".pt", ".npz", ".npy",
    ".pyc", ".gz", ".zip", ".woff", ".woff2",
}

_USER_AGENT = "DeepCatch DOI check (mailto:deepcatch@example.org)"

# Identifiers the repository mentions on purpose, because it documents them
# as placeholders or as dead. Each is named where it appears, so the
# mention is the evidence; a checker that flagged them would train its
# reader to ignore it.
_INTENTIONALLY_UNRESOLVED = {
    "10.5281/zenodo.to-be-assigned": (
        "the placeholder string the citation guard uses as a fixture"),
    "10.5281/zenodo.to-be-assigned-on-first-release": (
        "the placeholder the repository removed; the CHANGELOG and "
        "ZENODO_SETUP.md explain why"),
    "10.1101/2023.10.13.23296542": (
        "the Finaletoolkit identifier, which does not resolve; the "
        "bibliography file says so where it would have been"),
}


def _tracked_text_files() -> list[Path]:
    """Every tracked file worth reading as text."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), "ls-files"],
            capture_output=True, text=True, timeout=60, check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        raise SystemExit("not a git checkout; the tracked tree cannot be listed")
    out = []
    for name in proc.stdout.splitlines():
        path = _REPO_ROOT / name
        if path.suffix.lower() in _BINARY_SUFFIXES or not path.is_file():
            continue
        out.append(path)
    return out


def _dois_in_tree() -> dict[str, set[str]]:
    """Every DOI in the tracked tree, and where it is mentioned."""
    found: dict[str, set[str]] = {}
    for path in _tracked_text_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for raw in _DOI.findall(text):
            doi = _trim_unbalanced(raw)
            doi = _trim_unbalanced(doi.rstrip("`}.,;:\\，；：。（）、"))
            if doi:
                found.setdefault(doi, set()).add(str(path.relative_to(_REPO_ROOT)))
    return found


def _fetch(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(request, timeout=40) as response:
        return json.load(response)


def _is_registered(doi: str) -> tuple[bool | None, str]:
    """(registered, note). None means the registry could not be reached."""
    zenodo = doi.startswith("10.5281/")
    if zenodo:
        url = "https://api.datacite.org/dois/" + urllib.parse.quote(doi, safe="")
        try:
            data = _fetch(url)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return False, "not in DataCite"
            return None, f"DataCite error {exc.code}"
        except Exception as exc:  # network, timeout, DNS
            return None, f"DataCite unreachable ({type(exc).__name__})"
        title = data["data"]["attributes"].get("titles", [{}])[0].get("title", "")
        return True, title

    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    try:
        data = _fetch(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return False, "not in Crossref"
        return None, f"Crossref error {exc.code}"
    except Exception as exc:
        return None, f"Crossref unreachable ({type(exc).__name__})"
    title = (data["message"].get("title") or [""])[0]
    return True, title


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true",
                        help="print only the identifiers that fail")
    parser.add_argument("--delay", type=float, default=1.0,
                        help="seconds to wait between registry queries")
    args = parser.parse_args()

    dois = _dois_in_tree()
    if not dois:
        print("no DOIs found in the tracked tree")
        return 0

    failures: list[str] = []
    unreachable = False
    for doi in sorted(dois):
        if doi in _INTENTIONALLY_UNRESOLVED:
            if not args.quiet:
                print(f"SKIP {doi}  (mentioned on purpose: "
                      f"{_INTENTIONALLY_UNRESOLVED[doi]})")
            continue
        registered, note = _is_registered(doi)
        if registered is False:
            failures.append(f"{doi}  ({note}; mentioned in "
                            f"{', '.join(sorted(dois[doi])[:3])})")
        elif registered is None:
            unreachable = True
            if not args.quiet:
                print(f"SKIP {doi}  ({note})")
        elif not args.quiet:
            print(f"OK   {doi}  {note[:70]}")
        time.sleep(args.delay)

    print()
    print(f"distinct DOIs checked: {len(dois)}")
    if failures:
        print(f"not registered: {len(failures)}")
        for line in failures:
            print("  ", line)
        return 1
    print("not registered: 0")
    if unreachable:
        print("note: at least one registry was unreachable, so some "
              "identifiers were not checked")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
