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


#: Column names that identify the sample / label / study columns.
SAMPLE_NAMES = ("sample", "sample_id", "id")
LABEL_NAMES = ("label", "class", "y", "target")
STUDY_NAMES = ("study", "study_id", "publication", "cohort", "dataset")
DISEASE_NAMES = ("disease_class", "disease", "cancer_type")


def looks_like_header(parts) -> bool:
    """Is this row a header row rather than a data row?

    A header names its columns. The check is on the NAMES, not on the
    row's position: a data row whose first cell happens to be "sample"
    (a sample id) must not be mistaken for a header and dropped.
    """
    names = [p.strip().lower() for p in parts]
    if not names:
        return False
    if names[0] in SAMPLE_NAMES:
        return True
    return any(n in LABEL_NAMES or n in STUDY_NAMES or n in DISEASE_NAMES
               for n in names)


def resolve_columns(parts, headered: bool) -> tuple[int, int | None]:
    """Return (label column, study column) for one labels row.

    With a header, the columns are found by NAME. Without one, the
    layout is decided by the WIDTH, because the two documented
    headerless layouts are different:

        sample \\t label                       -> label 1, study -
        sample \\t label \\t study              -> label 1, study 2
        sample \\t disease_class \\t label \\t study
                                               -> label 2, study 3

    Reading the label positionally in the third layout picks up the
    disease class ("CRC"), and reading it in the second picks up the
    label — which is how two readers of the same file ended up with
    different label vectors.
    """
    if headered:
        label_col, study_col = 1, None
        for index, name in enumerate(parts):
            lowered = name.strip().lower()
            if lowered in LABEL_NAMES:
                label_col = index
            elif lowered in STUDY_NAMES:
                study_col = index
        return label_col, study_col
    width = len(parts)
    if width >= 4:
        return 2, 3
    if width == 3:
        return 1, 2
    return 1, None


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
        f"Unrecognized label {raw!r}{where}; expected one of "
        f"{sorted(POSITIVE | NEGATIVE)}")
