"""The canonical label vocabulary.

A label cell says either "this sample is a cancer case" or "this sample
is a non-cancer control", and the words that say so are a closed
vocabulary. Reading a label by exact string equality --
``1 if cell == "cancer" else 0`` -- turns every other spelling
("Cancer", "tumor", "control", "CRC") into a **silent negative**: the
sample is filed as a healthy donor, and the benchmark quietly learns
the wrong thing while reporting a plausible number.

Two readers of the same labels file did exactly that, and produced
different label vectors from the same file -- one of them wrote the
headline cross-study AUC. This module is the single place the
vocabulary lives; every reader resolves through it, and an unknown
label is an error to be fixed rather than a control to be invented.
"""
from __future__ import annotations

#: Label cells that mean "cancer case" (compared case-insensitively).
POSITIVE = frozenset({"cancer", "positive", "tumor", "y", "1", "true", "case"})

#: Label cells that mean "non-cancer control".
NEGATIVE = frozenset({"healthy", "control", "normal", "n", "0", "false"})


class LabelVocabularyError(ValueError):
    """A label cell that is neither a case nor a control."""


def resolve_label(raw: str, sample: str = "") -> int:
    """Map a label cell to 1 (cancer) or 0 (non-cancer control).

    Raises :class:`LabelVocabularyError` on anything else. The sample
    id, when given, is included in the message so the offending row can
    be found without re-reading the file.
    """
    value = (raw or "").strip().lower()
    if value in POSITIVE:
        return 1
    if value in NEGATIVE:
        return 0
    where = f" for sample {sample!r}" if sample else ""
    raise LabelVocabularyError(
        f"Unrecognised label {raw!r}{where}; expected one of "
        f"{sorted(POSITIVE | NEGATIVE)}")
