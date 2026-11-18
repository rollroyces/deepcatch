"""A citation file must not declare an identifier that does not exist.

Why this file exists
--------------------
`CITATION.cff` is what a reader's citation tool reads when they cite this
repository. It once carried a placeholder string under `doi:` -- a DOI that
resolved to nothing, which rendered in GitHub's "Cite this repository" box
as a link to a record that does not exist. That key was removed, and the
reason it was removed is written into the file itself.

The same defect survived in another form: an `identifiers:` entry with
`type: doi` whose value was the string "to be assigned (preprint)". It is
not a DOI, but it sits in a machine-readable field that says `type: doi`,
so a citation tool would record it as one. A promise printed where an
identifier belongs is not a citation; it is a broken link with better
manners.

What this test checks, and what it deliberately does not
-------------------------------------------------------
If the file declares a DOI -- as a `doi:` key or as an `identifiers:` entry
of type `doi` -- that DOI has to look like one. The file may declare no DOI
at all while none exists; that is the honest state and it is what the
repository is in until a release carries one.

It does not check that a DOI resolves. A DOI that is well formed but
records a state Zenodo has not archived yet is a different question, and
the answer belongs to the release, not to a test.
"""
from __future__ import annotations

import re
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent.parent
_CITATION = _REPO_ROOT / "CITATION.cff"

# A DOI: a registrant prefix, a slash, and a suffix.
_DOI = re.compile(r"^10\.\d{4,9}/\S+$")

# The top-level `doi:` key.
_TOP_LEVEL_DOI = re.compile(r'^doi:\s*"?([^"\n]+?)"?\s*$', re.M)

# An entry of an `identifiers:` list.
_IDENTIFIER_TYPE = re.compile(r"^\s*-\s+type:\s*(\S+)\s*$")
_IDENTIFIER_VALUE = re.compile(r'^\s*value:\s*"?([^"\n]+?)"?\s*$')


def _identifier_dois(text: str) -> list[str]:
    """The values of the `identifiers:` entries whose type is `doi`."""
    out: list[str] = []
    current_type: str | None = None
    in_list = False
    for line in text.splitlines():
        if re.match(r"^identifiers:\s*$", line):
            in_list = True
            continue
        if not in_list:
            continue
        if line and not line[0].isspace():
            in_list = False
            continue
        type_match = _IDENTIFIER_TYPE.match(line)
        if type_match:
            current_type = type_match.group(1)
            continue
        value_match = _IDENTIFIER_VALUE.match(line)
        if value_match and current_type == "doi":
            out.append(value_match.group(1))
    return out


def _declared_dois(text: str) -> list[tuple[str, str]]:
    """(where, value) for every DOI the file declares."""
    found = [(f"doi: key", m.group(1)) for m in _TOP_LEVEL_DOI.finditer(text)]
    found += [("identifiers: entry", value)
              for value in _identifier_dois(text)]
    return found


def test_the_citation_file_declares_no_fake_identifier():
    """Every DOI in CITATION.cff must look like a DOI."""
    assert _CITATION.is_file(), "CITATION.cff is missing"
    text = _CITATION.read_text(encoding="utf-8")

    offenders = [
        f"{where}: {value!r}"
        for where, value in _declared_dois(text)
        if not _DOI.match(value)
    ]
    assert not offenders, (
        "CITATION.cff declares a DOI that does not look like one:\n  "
        + "\n  ".join(offenders)
        + "\nA placeholder in an identifier field is a broken link with "
          "better manners. Declare no identifier while none exists; the "
          "release adds the real one (ZENODO_SETUP.md section 5)."
    )


def test_the_detector_sees_a_placeholder():
    """The detector must catch one, or the test above is vacuous.

    The fixtures are assembled at run time, so this file does not contain a
    placeholder of its own for the detector to trip over.
    """
    placeholder = "to be assigned (preprint)"
    text = f'identifiers:\n  - type: doi\n    value: "{placeholder}"\n'
    assert _declared_dois(text) == [("identifiers: entry", placeholder)]

    real = "10.5281/zenodo.14013719"
    text = f'identifiers:\n  - type: doi\n    value: "{real}"\n'
    assert _declared_dois(text) == [("identifiers: entry", real)]
    assert _DOI.match(real)


def test_the_top_level_doi_key_is_checked():
    """A `doi:` key is read the same way an identifiers entry is."""
    text = 'doi: "10.5281/zenodo.14013719"\n'
    assert _declared_dois(text) == [("doi: key", "10.5281/zenodo.14013719")]
    assert _DOI.match("10.5281/zenodo.14013719")
    # The rule is the shape, not resolvability: a well-formed DOI that
    # records no record is the release's business, not this test's.
    assert _DOI.match("10.5281/zenodo.to-be-assigned")
    assert not _DOI.match("to be assigned (preprint)")
