"""The public tree must not carry a personal mailbox address.

Why this file exists
--------------------
The April 2026 v1.0 draft (`paper/deepcatch_final.tex`, deprecated since
and marked "do not use" in `paper/README.md`) carried the author's
personal email address in two places: the author block and the
abstract's contact line. It was committed, so it is public, and a
public git history cannot be taken back -- the address is still in this
repository's history even though the current tree no longer carries it.

Nothing would have put it back on its own. A draft is written, an
earlier version is kept "for reference", and the contact line comes
along for the ride. That is what happened once, and it will happen
again unless something checks.

What this test checks, and what it deliberately does not
-------------------------------------------------------
It does NOT contain the author's address. Writing the address into
this test file would put it back into the public tree, which is the
exact defect the test exists to prevent -- a guard that reproduces the
address in order to search for it.

So the rule is about the *kind* of address instead: an address at a
consumer mailbox provider is a personal mailbox, and a personal
mailbox does not belong in a public repository. Institutional addresses
are a different matter and are left alone: the FLARE data provider's
contact is functional -- it is how a collaborator requests the data --
and the project's own `deepcatch+sponsors@` alias is a role address.
Both are already public in the upstream records they come from.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parent.parent

# Consumer mailbox providers. An address at one of these is a personal
# mailbox. This list is a list of DOMAINS, not of addresses: no address
# appears in this file.
_PERSONAL_DOMAINS = frozenset({
    "me.com",
    "gmx.com",
    "gmail.com",
    "googlemail.com",
    "hotmail.com",
    "outlook.com",
    "live.com",
    "yahoo.com",
    "icloud.com",
    "me.net",
    "aol.com",
    "proton.me",
    "protonmail.com",
    "mail.com",
    "zoho.com",
})

# A local part, an @, and a dotted domain.
_ADDRESS = re.compile(r"[A-Za-z0-9._%+-]+@([A-Za-z0-9.-]+\.[A-Za-z]{2,})")

# Files worth reading as text. A binary is skipped rather than guessed at:
# a false positive inside a PNG would be noise, and the addresses that
# matter are in the text of manuscripts and documentation.
_TEXT_SUFFIXES = {
    ".md", ".py", ".tex", ".txt", ".yaml", ".yml", ".cff", ".toml",
    ".cfg", ".ini", ".json", ".sh", ".bib", ".html", ".csv", ".tsv",
}
_TEXT_NAMES = {"AUTHORS", "LICENSE", "CITATION.cff", "NOTICE"}


def _personal_domains(text: str) -> list[str]:
    """The personal-mailbox domains used in `text`."""
    return sorted(
        domain.lower()
        for domain in _ADDRESS.findall(text)
        if domain.lower() in _PERSONAL_DOMAINS
    )


def _tracked_text_files() -> list[Path]:
    """Every tracked file that is worth reading as text."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(_REPO_ROOT), "ls-files", "-z"],
            capture_output=True, text=True, timeout=60, check=True,
        )
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout; the tracked tree cannot be listed")
    out = []
    for name in proc.stdout.split("\0"):
        if not name:
            continue
        path = _REPO_ROOT / name
        if path.suffix in _TEXT_SUFFIXES or path.name in _TEXT_NAMES:
            out.append(path)
    return out


def test_no_personal_mailbox_address_is_in_the_tracked_tree():
    """The tracked tree must not carry an address at a consumer mailbox."""
    offenders = []
    for path in _tracked_text_files():
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        domains = _personal_domains(text)
        if domains:
            offenders.append(
                f"{path.relative_to(_REPO_ROOT)}: {', '.join(domains)}")
    assert not offenders, (
        "a personal mailbox address is committed in the public tree:\n  "
        + "\n  ".join(offenders)
        + "\nKeep the author's contact address out of the repository; it "
          "belongs in the local, gitignored submission checklist."
    )


def test_the_detector_sees_a_personal_address():
    """The detector must actually see one, or the test above is vacuous."""
    assert _personal_domains("write to someone@example.org") == []
    assert _personal_domains("write to someone@gmail.com") == ["gmail.com"]
    assert _personal_domains("a@gmail.com b@yahoo.com") == [
        "gmail.com", "yahoo.com"]


def test_institutional_and_role_addresses_are_not_flagged():
    """A data provider's contact and the project's alias are functional.

    They are how a collaborator requests data and how a sponsor reaches
    the project. Flagging them would make the guard cry wolf, and a guard
    that cries wolf gets deleted.
    """
    assert _personal_domains(
        "Loris De Cecco (loris.dececco@istitutotumori.mi.it)") == []
    assert _personal_domains(
        "deepcatch+sponsors@rollroyces.dev") == []
    assert _personal_domains("lab-pi@example.org") == []
