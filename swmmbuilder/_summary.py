"""The flat dict that every ``summary()`` returns, with its former keys still readable."""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from typing import Any


class Summary(dict[str, Any]):
    """A ``summary()`` result: a plain dict of numbers under the current keys.

    A key renamed by #161 still reads (``s['old']``) with a ``DeprecationWarning`` naming the
    new key, but is not listed by ``keys()``, iteration, ``in``, or ``get``.
    """

    def __init__(self, values: Mapping[str, Any], renamed: Mapping[str, str] | None = None) -> None:
        super().__init__(values)
        self._renamed = {old: new for old, new in (renamed or {}).items() if new in self}

    def __missing__(self, key: str) -> Any:
        new = self._renamed.get(key)
        if new is None:
            raise KeyError(key)
        warnings.warn(
            f"The summary key {key!r} is now {new!r}; the old key will be removed.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self[new]
