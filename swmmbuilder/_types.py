"""Type aliases of the public API (for annotations only)."""

import os
from collections.abc import Iterable, Mapping
from typing import Any, TypeAlias

import pandas as pd

# A node or pipe ID, as in the data: text, or a number for cities that number them (Kitchener).
NodeId: TypeAlias = str | int
PipeId: TypeAlias = str | int

# A file path.
StrPath: TypeAlias = str | os.PathLike[str]

# The outfalls of a model: {node: outfall type}, the outfall nodes, or one node.
Outfalls: TypeAlias = Mapping[NodeId, str] | Iterable[NodeId] | NodeId

# Values per pipe, node or group ID, for plots and exports.
Values: TypeAlias = pd.Series | pd.DataFrame | Mapping[Any, float]

# An option given as a dataclass or as a dict of its fields (e.g. a Pump or {'curve': '*'}).
Fields: TypeAlias = Mapping[str, Any]
