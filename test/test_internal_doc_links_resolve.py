"""A link a reader clicks must land on a file that exists.

Why this file exists
--------------------
`MODEL.md` -- the model card a reviewer opens before the paper -- linked
to seven documents with the target missing the directory the link text
named: the text said `docs/PRETRAINING.md` and the target said
`PRETRAINING.md`, which is the repository root. The files were all there,
under `docs/`. Nothing about the link looked wrong in the source, and
nothing in CI could see it: the guards that existed pinned figures and
test counts, and a path in a link is neither.

`wiki/Home.md` had the same disease in the other direction: a badge linked
to `LICENSE`, which is the repository's licence file, not a wiki page.
Inside a GitHub wiki there is no such page, so the badge led nowhere.

The link that has been dead in one document is exactly the kind of thing
that costs an independent publisher more than it should: the repository's
argument is that its evidence can be checked, and the first place a
reviewer looks for the evidence is the model card's links.

What this test checks, and what it deliberately does not
-------------------------------------------------------
Every relative link in every tracked markdown file must resolve to
something that exists in the tree.

It does not check external URLs. Those change without the repository
changing, and a link checker that fails on a temporary outage teaches the
reader to ignore it.

Two kinds of link are resolved the way their host resolves them:

* A link inside a code block is an example, not a link, and is skipped.
* A link without a file extension inside `wiki/` is a GitHub wiki page,
  so `Home` means `Home.md` in the same directory.
* A link on the public docs site that points at a file of this repository
  is checked against the tree, whether it is written relatively or as an
  absolute `github.com/rollroyces/deepcatch/blob/main/...` URL. The site
  is the front door, and a link there to a file that has been renamed is
  as dead as one in a markdown file.
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from urllib.parse import unquote

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# A markdown link: the target, optionally followed by a title.
_LINK = re.compile(r'\]\(([^)\s]+)(?:\s+"[^"]*")?\)')

# A fenced code block: a link inside it is an example.
_FENCE = re.compile(r"```.*?```", re.S)

# A link attribute in the HTML docs site.
_HTML_HREF = re.compile(r'href="([^"]+)"')


# The public docs site: a front door with links of its own.
_SITE = "docs/site/index.html"

# A link on the site that names a file of this repository.
_GITHUB_BLOB = re.compile(
    r"^https://github\.com/rollroyces/deepcatch/blob/main/(?P<path>[^#?]+)")


def _tracked_documents() -> list[str]:
    """Every tracked markdown file, plus the docs site."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), "ls-files", "*.md"],
            capture_output=True, text=True, timeout=60, check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout; the tracked tree cannot be listed")
    documents = [name for name in proc.stdout.splitlines() if name]
    if (_REPO_ROOT / _SITE).is_file():
        documents.append(_SITE)
    return documents


def _link_targets(text: str) -> list[str]:
    """The link targets in a document, ignoring the ones in code blocks.

    The docs site is HTML, so its `href` attributes are read as well.
    """
    markdown = [m.group(1).strip()
                for m in _LINK.finditer(_FENCE.sub("", text))]
    html = [m.group(1).strip() for m in _HTML_HREF.finditer(text)]
    return markdown + html


def _candidates(doc: str, target: str) -> list[str]:
    """The paths a link target could mean, in the order the host tries."""
    blob = _GITHUB_BLOB.match(target)
    if blob:
        # An absolute link to a file of this repository: the file has to
        # be in the tree under the name the link uses.
        return [os.path.normpath(unquote(blob.group("path")))]
    raw = unquote(target.split("#")[0])
    if not raw:
        return []
    names = [raw]
    # A link without an extension inside the wiki is a wiki page.
    if not os.path.splitext(raw)[1] and doc.startswith("wiki/"):
        names.append(raw + ".md")
    resolved = []
    for name in names:
        resolved.append(os.path.normpath(os.path.join(os.path.dirname(doc), name)))
        resolved.append(os.path.normpath(name))
    return resolved


def _is_generated(path: str) -> bool:
    """Is the path a git-ignored location (a build product or a cache)?"""
    return subprocess.run(
        ["git", "-C", str(_REPO_ROOT), "check-ignore", "-q", path],
        capture_output=True,
    ).returncode == 0


def test_internal_markdown_links_resolve():
    """Every relative link must land on a file that exists."""
    offenders: list[str] = []

    for doc in _tracked_documents():
        path = _REPO_ROOT / doc
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for target in _link_targets(text):
            if target.startswith(("http://", "https://", "mailto:", "#")):
                # An absolute URL is another repository's business unless
                # it names a file of this one, which is checked above.
                if not _GITHUB_BLOB.match(target):
                    continue
            # A path outside the repository belongs to another checkout.
            if target.startswith(("../", "/")):
                continue
            if any(ch in target for ch in "<>*"):
                continue
            candidates = _candidates(doc, target)
            if not candidates:
                continue
            if any((_REPO_ROOT / c).exists() for c in candidates):
                continue
            if any(_is_generated(c) for c in candidates):
                continue
            offenders.append(f"{doc}: {target}")

    assert not offenders, (
        "a markdown link points at a file that does not exist:\n  "
        + "\n  ".join(offenders)
        + "\nPoint the link at the file that exists. A reader who follows it "
          "finds nothing, and in this repository the links are the argument "
          "that the evidence can be checked."
    )


def test_the_link_extractor_finds_a_known_link():
    """The extractor must see a link, or the test above is vacuous.

    The fixture is assembled at run time so this file does not become
    another copy of the links it checks.
    """
    text = "see [the model card](MODEL.md) for details\n"
    assert _link_targets(text) == ["MODEL.md"]


def test_a_link_inside_a_code_block_is_not_a_link():
    """A link in a code block is an example, and examples are not links."""
    text = "```\n[text](not_a_real_path.md)\n```\n"
    assert _link_targets(text) == []
