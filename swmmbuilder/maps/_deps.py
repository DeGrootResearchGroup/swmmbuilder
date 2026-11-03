"""Optional dependencies (the ``[maps]`` extra)."""

import importlib


def _optional(module):
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(
            f"This needs {module}, an optional dependency: pip install 'swmmbuilder[maps]'"
        ) from e
