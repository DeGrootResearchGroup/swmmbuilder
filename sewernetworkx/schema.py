"""Mapping from a city's GIS sewer schema to the canonical columns used by SewerNetworkx.

:class:`~sewernetworkx.Network` reads a GeoJSON file and passes it through
:meth:`Schema.normalize`, which renames, recodes and converts the columns. Everything downstream
reads only the canonical columns:

=================  ===========================================================================
``pipe_id``        unique pipe ID (the SWMM conduit name and the graph edge key)
``from_node``      upstream node ID (NaN if missing)
``to_node``        downstream node ID (NaN if missing)
``length``         pipe length, in m
``us_invert``      upstream invert, in m (NaN if missing)
``ds_invert``      downstream invert, in m (NaN if missing)
``drop_invert``    drop-structure invert at the downstream node, in m (NaN if missing)
``shape``          SWMM cross-section name (NaN if the shape code is unknown or missing)
``shape_code``     the source shape code
``diameter``       in m
``height``         full height, in m
``width``          maximum width, in m
``flow_type``      e.g. ``'sanitary'``, ``'storm'``, ``'combined'``; unmapped codes pass through
``material``       passed through unchanged (the key of ``manning_n_by_material``)
``abandoned``      bool
=================  ===========================================================================

Columns that are not mapped keep their source names.
"""

import warnings
from dataclasses import dataclass, field

import geopandas as gpd
import pandas as pd

# Canonical columns, in the order they are written.
CANONICAL_COLUMNS = (
    "pipe_id",
    "from_node",
    "to_node",
    "length",
    "us_invert",
    "ds_invert",
    "drop_invert",
    "shape",
    "shape_code",
    "diameter",
    "height",
    "width",
    "flow_type",
    "material",
    "abandoned",
)

# Canonical columns a schema must map, and that must be present in the source data.
REQUIRED_COLUMNS = ("pipe_id", "from_node", "to_node", "length")

# Units of length -> units per metre.
UNITS = {"m": 1, "cm": 100, "mm": 1000, "ft": 1 / 0.3048, "in": 1 / 0.0254}

INVERT_COLUMNS = ("us_invert", "ds_invert", "drop_invert")
SIZE_COLUMNS = ("diameter", "height", "width")


@dataclass(frozen=True)
class Schema:
    """How a city's GIS sewer data maps to the canonical columns.

    Parameters
    ----------
    name : str
        Name of the schema, e.g. the city.
    columns : dict
        Canonical column -> source column. ``pipe_id``, ``from_node``, ``to_node`` and ``length``
        are required; ``shape_code`` and ``abandoned`` are derived and cannot be mapped.
        ``columns['shape']`` is the source column holding the shape code.
    shapes : dict
        Source shape code -> SWMM cross-section name (CIRCULAR uses ``diameter``, EGG ``height``,
        the others ``height`` and ``width``).
    flow_types : dict
        Source flow-type code -> canonical flow type. Codes not in the mapping pass through.
    abandoned : str or None
        Source column that flags abandoned pipes.
    abandoned_values : frozenset
        Values of ``abandoned`` that mean the pipe is abandoned.
    null_nodes : frozenset
        Node IDs that mean "no node" (set to NaN).
    nonpositive_inverts_missing : bool
        Treat inverts at or below 0 as missing (some GIS exports use 0 for none).
    length_unit, invert_unit, size_unit : str
        Units of the source lengths, inverts and sizes: one of ``'m'``, ``'cm'``, ``'mm'``,
        ``'ft'``, ``'in'``.
    crs : str or None
        If set, the geometry is reprojected to this CRS on loading.
    """

    name: str
    columns: dict
    shapes: dict = field(default_factory=dict)
    flow_types: dict = field(default_factory=dict)
    abandoned: str | None = None
    abandoned_values: frozenset = frozenset()
    null_nodes: frozenset = frozenset()
    nonpositive_inverts_missing: bool = False
    length_unit: str = "m"
    invert_unit: str = "m"
    size_unit: str = "m"
    crs: str | None = None

    def __post_init__(self):
        unknown = set(self.columns) - set(CANONICAL_COLUMNS)
        derived = set(self.columns) & {"shape_code", "abandoned"}
        if unknown or derived:
            raise ValueError(
                f"Schema {self.name!r} maps columns that cannot be mapped: "
                f"{sorted(unknown | derived)}. Canonical columns: {list(CANONICAL_COLUMNS)}."
            )
        missing = [c for c in REQUIRED_COLUMNS if c not in self.columns]
        if missing:
            raise ValueError(f"Schema {self.name!r} does not map required columns: {missing}.")
        for unit in (self.length_unit, self.invert_unit, self.size_unit):
            if unit not in UNITS:
                raise ValueError(f"Unknown unit {unit!r}; use one of {list(UNITS)}.")

    def flow_type(self, code):
        """Canonical flow type of a source code (codes not in ``flow_types`` pass through)."""
        return self.flow_types.get(code, code)

    def normalize(self, data):
        """Return a copy of ``data`` with the canonical columns.

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Pipes as read from the GIS file.

        Returns
        ----------
        geopandas.GeoDataFrame
            The canonical columns, followed by the unmapped source columns and the geometry.

        Raises
        ----------
        ValueError
            If a required source column is missing or ``pipe_id`` is not unique.

        Warns
        ----------
        UserWarning
            For shape codes that are not in ``shapes`` (those pipes are rejected when written).
        """
        sources = dict(self.columns)
        if self.abandoned is not None:
            sources["abandoned"] = self.abandoned

        problems = []
        missing = [sources[c] for c in REQUIRED_COLUMNS if sources[c] not in data.columns]
        if missing:
            problems.append(f"missing required columns {missing}")
        elif data[sources["pipe_id"]].duplicated().any():
            duplicated = data.loc[data[sources["pipe_id"]].duplicated(), sources["pipe_id"]]
            problems.append(f"duplicate {sources['pipe_id']} values {sorted(set(duplicated))}")
        if problems:
            raise ValueError(
                f"Cannot load the data with schema {self.name!r}: {'; '.join(problems)}."
            )

        if self.crs is not None:
            data = data.to_crs(self.crs)

        def source(canonical):
            column = sources.get(canonical)
            if column is None or column not in data.columns:
                return pd.Series(pd.NA, index=data.index, dtype=object)
            return data[column]

        out = {}
        out["pipe_id"] = source("pipe_id")
        for node in ("from_node", "to_node"):
            out[node] = source(node).where(~source(node).isin(self.null_nodes))
        out["length"] = source("length") / UNITS[self.length_unit]
        for column in INVERT_COLUMNS:
            inverts = pd.to_numeric(source(column), errors="coerce") / UNITS[self.invert_unit]
            if self.nonpositive_inverts_missing:
                inverts = inverts.where(inverts > 0)
            out[column] = inverts
        codes = source("shape")
        out["shape"] = codes.map(self.shapes)
        out["shape_code"] = codes
        for column in SIZE_COLUMNS:
            out[column] = pd.to_numeric(source(column), errors="coerce") / UNITS[self.size_unit]
        out["flow_type"] = source("flow_type").map(self.flow_type)
        out["material"] = source("material")
        if self.abandoned is not None and self.abandoned in data.columns:
            flag = data[self.abandoned]
            numeric = pd.to_numeric(flag, errors="coerce")
            out["abandoned"] = flag.isin(self.abandoned_values) | numeric.isin(
                self.abandoned_values
            )
        else:
            out["abandoned"] = pd.Series(False, index=data.index)

        unknown = codes[codes.notna() & (codes.astype(str).str.strip() != "")]
        unknown = unknown[~unknown.isin(self.shapes)]
        if len(unknown):
            counts = unknown.astype(str).value_counts().sort_index().to_dict()
            warnings.warn(
                f"Shape codes not in schema {self.name!r} (pipes with them cannot be written to "
                f".inp): {counts}. Known codes: {list(self.shapes)}.",
                stacklevel=3,
            )

        geometry = data.geometry.name
        mapped = set(sources.values()) | set(CANONICAL_COLUMNS) | {geometry}
        rest = [c for c in data.columns if c not in mapped]
        frame = pd.concat([pd.DataFrame(out, index=data.index), data[[*rest, geometry]]], axis=1)
        return gpd.GeoDataFrame(frame, geometry=geometry, crs=data.crs)


# London, ON. In the London GIS, PipeShape 'R' means round; rectangular pipes are 'REC'. Sizes
# are in mm, an invert of 0 means none, and a node ID of 'Abandoned' means no node.
LONDON_ON = Schema(
    name="London, ON",
    columns={
        "pipe_id": "OBJECTID",
        "from_node": "UpstreamID",
        "to_node": "DownstreamID",
        "length": "Length",
        "us_invert": "UpstreamInvert",
        "ds_invert": "DownstreamInvert",
        "drop_invert": "DropStructureInvert",
        "shape": "PipeShape",
        "diameter": "Diameter",
        "height": "Height",
        "width": "Width",
        "flow_type": "FlowType",
        "material": "Material",
    },
    shapes={
        "R": "CIRCULAR",
        "Round": "CIRCULAR",
        "REC": "RECT_CLOSED",
        "EGG": "EGG",
        "ARC": "ARCH",
        "ELH": "HORIZ_ELLIPSE",
        "ELV": "VERT_ELLIPSE",
    },
    flow_types={"SAN": "sanitary", "STM": "storm", "COM": "combined"},
    abandoned="Abandoned",
    abandoned_values=frozenset({1}),
    null_nodes=frozenset({"Abandoned"}),
    nonpositive_inverts_missing=True,
    size_unit="mm",
)
