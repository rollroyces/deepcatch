"""The repository's version, read from one place so it cannot drift.

The API's OpenAPI metadata and the clinical report both used to carry
a hardcoded version string, and both had drifted: they still said
2.1.0 while the repository was at 2.3.0. A reader who ran the API
(the demo recording does exactly that) saw a version that did not match
the code.

So the version is read rather than written:

* from the installed distribution's metadata, when the package is
  installed -- which is what a user or a reviewer has;
* otherwise from ``pyproject.toml``, which is where it is written down
  in a source checkout.
"""
from __future__ import annotations

from pathlib import Path

DISTRIBUTION = "deepcatch"
FALLBACK = "unknown"


def package_version(default: str = FALLBACK) -> str:
    """The distribution's version, or `default` if it cannot be found."""
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8
        version = None  # type: ignore[assignment]
    else:
        try:
            return version(DISTRIBUTION)
        except PackageNotFoundError:
            pass

    pyproject = Path(__file__).resolve().parents[2] / "pyproject.toml"
    if pyproject.exists():
        for line in pyproject.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("version ="):
                return stripped.split("=", 1)[1].strip().strip('"').strip("'")
    return default
