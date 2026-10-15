"""Stub for `temporal_transformer_fixed` (replaces broken symlink).

Original repository layout used symlinks from `src/longitudinal/`
into per-agent development directories under `src/agent3-longitudinal/`. Those agent
directories were never committed to git, leaving these symlinks pointing to
non-existent paths on every clone (broken at HEAD and in working tree).

This stub preserves the module path so existing imports do not break, but raises
`NotImplementedError` if any callable is actually invoked. To restore real
functionality, port the original agent3 (longitudinal) implementation into this file
(or restore it from an external backup) and remove this stub.

Bug surfaced during cross-repo shallow review (commit history shows the
symlinks were committed in the initial import — the target directories never
existed in any branch).
"""


class TemporalTransformerPlaceholder:
    """Placeholder class — original implementation was in `src/agent3-longitudinal/`."""

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "temporal_transformer_fixed is a placeholder stub. The original implementation "
            "lived in `src/agent3-longitudinal/` and was never committed. Restore the "
            "source to enable this module."
        )


def _not_implemented(*args, **kwargs):  # pragma: no cover
    raise NotImplementedError(
        "temporal_transformer_fixed is a placeholder stub. The original implementation "
        "lived in `src/agent3-longitudinal/` and was never committed. Restore the "
        "source to enable this module."
    )


# Convenience aliases — replaced callers can import the original names without
# ImportError. They will fail loudly at call time.
__all__ = ["TemporalTransformerPlaceholder"]
