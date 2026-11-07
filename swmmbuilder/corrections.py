"""Documented corrections to pipes, applied by :class:`~swmmbuilder.Network` on loading.

Every municipal export needs a few corrections that a model should record: pipes to drop
(loops, self-loops, force mains to ignore), statuses to override (built pipes still marked
proposed), node IDs recorded the wrong way round, and single attribute fixes. Pass them as
``Network(..., corrections=...)``, so that they are logged with their reasons and listed in
``Network.corrections_report``, instead of editing the data before loading it.
"""

from __future__ import annotations

import os
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .schema import CANONICAL_COLUMNS

if TYPE_CHECKING:
    from ._types import StrPath

# Canonical columns that ``set`` may not change: the pipe ID is the key of the correction, and
# ``abandoned`` follows ``status``.
FIXED_COLUMNS = frozenset({"pipe_id", "abandoned"})

# Columns of Network.corrections_report.
REPORT_COLUMNS = ["pipe_id", "action", "column", "old", "new", "reason"]

ACTIONS = ("drop", "set", "swap_nodes")


@dataclass(frozen=True)
class PipeCorrections:
    """Corrections to pipes, keyed by the pipe ID in the source data.

    Pipe IDs match the data's IDs exactly or, failing that, by their text, so that ``112216``
    and ``'112216'`` name the same pipe (and the keys of a TOML file, which are always text,
    match numeric IDs).

    Parameters
    ----------
    drop : dict
        {pipe ID: reason}. The pipe is removed before any other filter, and counted in
        ``Network.removed['corrections']['dropped']``.

    set : dict
        {pipe ID: {canonical column: value, ..., 'reason': str}}. Each value replaces the
        pipe's value of a canonical column (see :mod:`swmmbuilder.schema`), so it is in the
        canonical units and codes: ``{'status': 'built'}``, not the source code ``'Built'``;
        inverts and sizes in m. ``pipe_id`` and ``abandoned`` cannot be set (``abandoned``
        follows ``status``). A ``from_node`` or ``to_node`` value goes through the network's
        ``normalize_ids``. ``'reason'`` is optional.

    swap_nodes : dict
        {pipe ID: reason}. The pipe's ``from_node`` and ``to_node`` are swapped, with their
        node types (``us_node_type`` and ``ds_node_type``, which describe the node of each ID),
        and nothing else: its inverts and geometry stay as recorded. This is for a pipe whose
        IDs were recorded the wrong way round while its drawing and inverts agree with each
        other (the geometry runs from the higher invert to the lower one, and the IDs say the
        opposite). ``us_invert`` is the invert at the pipe's upstream end, which is the end the
        drawing starts from; after the swap ``from_node`` is the node at that end, so the
        upstream invert stays with the pipe's upstream end. To reverse a pipe whose inverts are
        wrong instead, ``set`` its ``us_invert`` and ``ds_invert``.

    Raises
    ------
    ValueError
        If a reason is not a string, a ``set`` entry is not a mapping or names a column that is
        not canonical (or cannot be set), or the same pipe is both dropped and set or swapped,
        or both swapped and given a new ``from_node`` or ``to_node``.
    """

    drop: Mapping[Any, str] = field(default_factory=dict)
    set: Mapping[Any, Mapping[str, Any]] = field(default_factory=dict)
    swap_nodes: Mapping[Any, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for action in ("drop", "swap_nodes"):
            for pipe, reason in getattr(self, action).items():
                if not isinstance(reason, str):
                    raise ValueError(
                        f"corrections {action}: the value for pipe {pipe} must be its reason (a "
                        f"string), not {reason!r}."
                    )
        settable = set(CANONICAL_COLUMNS) - FIXED_COLUMNS
        for pipe, values in self.set.items():
            if not isinstance(values, Mapping):
                raise ValueError(
                    f"corrections set: the value for pipe {pipe} must be a mapping of canonical "
                    f"columns to values, not {values!r}."
                )
            unknown = sorted(set(values) - settable - {"reason"})
            if unknown:
                raise ValueError(
                    f"corrections set: pipe {pipe} names columns that cannot be set: {unknown}. "
                    f"Settable columns: {sorted(settable)}."
                )
            reason = values.get("reason")
            if reason is not None and not isinstance(reason, str):
                raise ValueError(f"corrections set: the reason for pipe {pipe} must be a string.")
        conflicts = [
            *(f"{p} (drop and set)" for p in _common(self.drop, self.set)),
            *(f"{p} (drop and swap_nodes)" for p in _common(self.drop, self.swap_nodes)),
            *(
                f"{p} (swap_nodes and set {'/'.join(c for c in ('from_node', 'to_node') if c in v)})"
                for p, v in self.set.items()
                if str(p) in {str(s) for s in self.swap_nodes}
                and ({"from_node", "to_node"} & set(v))
            ),
        ]
        if conflicts:
            raise ValueError(f"corrections: conflicting entries for pipes {', '.join(conflicts)}.")

    def __len__(self) -> int:
        return len(self.drop) + len(self.set) + len(self.swap_nodes)

    @classmethod
    def from_toml(cls, path: StrPath) -> PipeCorrections:
        """Read corrections from a TOML file with the tables ``[drop]``, ``[set.<pipe ID>]``,
        and ``[swap_nodes]``, e.g.

        .. code-block:: toml

            [drop]
            112216 = "force main taken as gravity pass-through"

            [set.328065]
            status = "built"
            reason = "Proposed in GIS, built"

            [swap_nodes]
            144019 = "IDs recorded the wrong way round; drawing and inverts agree"

        Raises
        ------
        ValueError
            If the file has other top-level keys, or as :class:`PipeCorrections`.
        """
        with open(path, "rb") as file:
            values = tomllib.load(file)
        return cls.coerce(values, source=f"corrections file {path}")

    @classmethod
    def coerce(
        cls,
        corrections: PipeCorrections | Mapping[str, Mapping[Any, Any]] | StrPath,
        source: str = "corrections",
    ) -> PipeCorrections:
        """``corrections`` as a :class:`PipeCorrections`: itself, built from a mapping with the
        keys ``drop``, ``set``, and ``swap_nodes`` (each optional), or read from a TOML file
        (:meth:`from_toml`).

        Raises
        ------
        ValueError
            If a mapping has other keys, or as :class:`PipeCorrections`.
        """
        if isinstance(corrections, PipeCorrections):
            return corrections
        if isinstance(corrections, str | os.PathLike):
            return cls.from_toml(Path(corrections))
        unknown = sorted(set(corrections) - set(ACTIONS))
        if unknown:
            raise ValueError(f"Unknown keys in {source}: {unknown}. Known keys: {list(ACTIONS)}.")
        return cls(**{action: dict(corrections.get(action) or {}) for action in ACTIONS})


def _common(first, second):
    """The keys of ``first`` that name (by text) a key of ``second``."""
    texts = {str(key) for key in second}
    return [key for key in first if str(key) in texts]
