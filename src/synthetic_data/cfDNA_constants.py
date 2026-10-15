"""Stub for `cfDNA_constants` (replaces broken symlink).

Original repository layout used symlinks from `src/synthetic_data/`
into per-agent development directories under `src/agent5-synthetic-data/`. Those agent
directories were never committed to git, leaving these symlinks pointing to
non-existent paths on every clone (broken at HEAD and in working tree).

This stub preserves the module path so existing imports do not break, but raises
`NotImplementedError` if any callable is actually invoked. To restore real
functionality, port the original agent5 (synthetic data) implementation into this file
(or restore it from an external backup) and remove this stub.

Bug surfaced during cross-repo shallow review (commit history shows the
symlinks were committed in the initial import — the target directories never
existed in any branch).
"""


class CfdnaConstantsPlaceholder:
    """Placeholder class — original implementation was in `src/agent5-synthetic-data/`."""

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "cfDNA_constants is a placeholder stub. The original implementation "
            "lived in `src/agent5-synthetic-data/` and was never committed. Restore the "
            "source to enable this module."
        )


def _not_implemented(*args, **kwargs):  # pragma: no cover
    raise NotImplementedError(
        "cfDNA_constants is a placeholder stub. The original implementation "
        "lived in `src/agent5-synthetic-data/` and was never committed. Restore the "
        "source to enable this module."
    )


# Convenience aliases — replaced callers can import the original names without
# ImportError. They will fail loudly at call time.
__all__ = ["CfdnaConstantsPlaceholder"]
