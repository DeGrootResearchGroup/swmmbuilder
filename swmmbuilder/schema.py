"""Mapping from a city's GIS sewer schema to the canonical columns used by swmmbuilder.

:class:`~swmmbuilder.Network` reads a GeoJSON file and passes it through
:meth:`Schema.normalize`, which renames, recodes and converts the columns. Everything downstream
reads only the canonical columns:

=================  ===========================================================================
``pipe_id``        unique pipe ID (the SWMM conduit name and the graph edge key)
``from_node``      upstream node ID (NaN if missing)
``to_node``        downstream node ID (NaN if missing)
``length``         pipe length, in m
``us_invert``      upstream invert, in m (NaN if missing or implausible)
``ds_invert``      downstream invert, in m (NaN if missing or implausible)
``drop_invert``    drop-structure invert at the downstream node, in m (NaN if missing or
                   implausible: see ``invert_range``)
``shape``          SWMM cross-section name (NaN if the shape code is unknown or missing)
``shape_code``     the source shape code
``diameter``       in m
``height``         full height, in m
``width``          maximum width, in m
``flow_type``      e.g. ``'sanitary'``, ``'storm'``, ``'combined'``; unmapped codes pass through
``material``       passed through unchanged (the key of ``manning_n_by_material``)
``status``         ``'built'``, ``'abandoned'``, ``'proposed'`` or ``'removed'``; unmapped codes
                   pass through. Missing: ``'abandoned'`` if the ``abandoned`` flag is set, else
                   ``'built'``
``abandoned``      bool: ``status == 'abandoned'``
``structure_type`` e.g. ``'gravity'``, ``'overflow'``, ``'force_main'``; unmapped codes pass
                   through
``us_node_type``   type of the upstream node (see ``node_type`` below), from the pipe's
                   inventory type
``ds_node_type``   type of the downstream node
=================  ===========================================================================

Columns that are not mapped keep their source names.

Node layers (maintenance holes, pumping stations, outfalls and other nodes) are read with
:meth:`Schema.normalize_nodes` into the canonical node columns; source columns that are not mapped
are dropped:

=================  ===========================================================================
``node_id``        node ID, as in the pipes' ``from_node`` and ``to_node``
``node_type``      ``'maintenance_hole'``, ``'pump_station'``, ``'outfall'``, ``'cap'``,
                   ``'tee'``, ``'chamber'``, ...; unmapped codes pass through, missing ones are
                   ``'other'``
``name``           e.g. the pumping station's name
``invert``         in m (NaN if missing)
``rim_elevation``  lid or top elevation, in m; NaN if missing or outside ``rim_range``
``depth``          rim to invert, in m
``has_drop``       bool (NA if missing): the node has a drop structure
``overflow``       bool: the node has an overflow
``weir``           bool: the node has a weir
``pipes_in``       number of pipes entering the node, as recorded in the node layer
``pipes_out``      number of pipes leaving the node
``flow_monitor``   bool: a flow monitor is installed
``access``         passed through, e.g. London's ``Location`` (road or easement)
``verified``       bool: the record has been verified
=================  ===========================================================================

Schemas can be written in TOML (see :meth:`Schema.from_toml`); the presets (:data:`LONDON_ON`,
:data:`TORONTO_ON`, :data:`BARRIE_ON`, :data:`KITCHENER_ON`, :data:`PETERBOROUGH_ON`,
:data:`SARNIA_ON`, :data:`WATERLOO_ON` and :data:`WINDSOR_ON`) are loaded from
``swmmbuilder/schemas``.
"""

from __future__ import annotations

import tomllib
import warnings
from dataclasses import dataclass, field, fields
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import pandas as pd

from .design import check_design_rule
from .diagnostics import DataWarning

if TYPE_CHECKING:
    from ._types import StrPath

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
    "status",
    "abandoned",
    "structure_type",
    "us_node_type",
    "ds_node_type",
)

# Canonical node columns, in the order they are written.
NODE_COLUMNS = (
    "node_id",
    "node_type",
    "name",
    "invert",
    "rim_elevation",
    "depth",
    "has_drop",
    "overflow",
    "weir",
    "pipes_in",
    "pipes_out",
    "flow_monitor",
    "access",
    "verified",
)
NODE_BOOL_COLUMNS = ("has_drop", "overflow", "weir", "flow_monitor", "verified")
TRUE_STRINGS = frozenset({"1", "Y", "YES", "T", "TRUE"})
FALSE_STRINGS = frozenset({"0", "N", "NO", "F", "FALSE"})

# Canonical columns a schema must map, and that must be present in the source data.
REQUIRED_COLUMNS = ("pipe_id", "from_node", "to_node")

# Units of length -> units per metre.
UNITS = {"m": 1, "cm": 100, "mm": 1000, "ft": 1 / 0.3048, "in": 1 / 0.0254}

INVERT_COLUMNS = ("us_invert", "ds_invert", "drop_invert")

# Design rules a schema's ``design`` table may hold (see :func:`swmmbuilder.gaps.fill_gaps` and
# the ``force_main_roughness`` option of :func:`swmmbuilder.swmm.build_model`).
DESIGN_RULES = (
    "force_main_discharge",
    "force_main_roughness",
    "wet_well_defaults",
    "force_main_checks",
    "default_pump",
)
SIZE_COLUMNS = ("diameter", "height", "width")


@dataclass(frozen=True)
class Schema:
    """How a city's GIS sewer data maps to the canonical columns.

    Parameters
    ----------
    name : str
        Name of the schema, e.g. the city.
    columns : dict
        Canonical column -> source column. ``pipe_id``, ``from_node`` and ``to_node`` are
        required; ``shape_code`` and ``abandoned`` are derived and cannot be mapped.
        ``columns['shape']`` is the source column holding the shape code. Missing lengths
        (including all of them, if ``length`` is not mapped) are computed from the geometry, which
        needs a projected CRS (see ``crs``).
    shapes : dict
        Source shape code -> SWMM cross-section name (CIRCULAR uses ``diameter``, EGG ``height``,
        the others ``height`` and ``width``).
    flow_types : dict
        Source flow-type code -> canonical flow type. Codes not in the mapping pass through.
    statuses : dict
        Source status code -> canonical status (``'built'``, ``'abandoned'``, ``'proposed'``,
        ``'removed'``). Codes not in the mapping pass through.
    structure_types : dict
        Source pipe structure-type code -> canonical structure type (e.g. ``'gravity'``,
        ``'overflow'``, ``'force_main'``). Codes not in the mapping pass through.
    node_types : dict
        Source node-type code -> canonical node type (e.g. ``'maintenance_hole'``,
        ``'pump_station'``, ``'outfall'``, ``'cap'``), for the pipes' inventory types and the
        ``node_type`` column of the other-nodes layer. Codes not in the mapping pass through.
    node_columns : dict
        Canonical node column -> source column of the node layers, or a list of source columns of
        which the first one present is used (the layers have different columns). ``node_id`` is
        required to read node layers.
    rim_range : tuple of float or None
        (lowest, highest) plausible rim elevation, in m; rims outside it are treated as missing.
    invert_range : tuple of float or None
        (lowest, highest) plausible pipe and drop-structure invert, in m; inverts outside it
        (e.g. 25089.0 for 250.89, or a drop height in the invert column) are treated as missing
        and listed in :attr:`Network.invert_report`. Not applied unless given (it is independent of
        ``rim_range``), e.g. ``dataclasses.replace(LONDON_ON, invert_range=LONDON_ON.rim_range)``.
    default_shape : str or None
        SWMM cross-section of every pipe, for data with no shape column (``shape`` must then not
        be mapped); e.g. ``'CIRCULAR'`` where only a diameter is published.
    abandoned : str or None
        Source column that flags abandoned pipes (older exports; ``status`` takes precedence).
    abandoned_values : frozenset
        Values of ``abandoned`` that mean the pipe is abandoned.
    null_values : frozenset
        Values that mean "missing" in every mapped column (e.g. the string ``'None'``).
    null_nodes : frozenset
        Node IDs that mean "no node" (set to NaN).
    nonpositive_inverts_missing : bool
        Treat inverts at or below 0 as missing (some GIS exports use 0 for none).
    length_unit, invert_unit, size_unit : str
        Units of the source lengths, inverts and sizes: one of ``'m'``, ``'cm'``, ``'mm'``,
        ``'ft'``, ``'in'``.
    crs : str or None
        If set, the geometry is reprojected to this CRS on loading. Data with no CRS are taken
        to be in this CRS already.
    design : dict
        Design rule -> its parameters, from a design standard (see :data:`DESIGN_RULES`), e.g.
        ``{'force_main_discharge': {'max_height_above_receiving_invert': 0.2}}``. They are used
        only when passed to the .inp writer, e.g.
        ``to_inp(..., force_main_discharge=schema.design['force_main_discharge'])`` or
        ``to_inp(..., force_main_roughness=schema.design['force_main_roughness'])``; likewise
        ``wet_well_defaults``, ``force_main_checks`` and ``default_pump``, or several at once with
        ``to_inp(..., design='schema')`` (see :func:`swmmbuilder.swmm.build_model`). So a schema
        can map the columns without adopting the standard.
    """

    name: str
    columns: dict
    shapes: dict = field(default_factory=dict)
    flow_types: dict = field(default_factory=dict)
    statuses: dict = field(default_factory=dict)
    structure_types: dict = field(default_factory=dict)
    node_types: dict = field(default_factory=dict)
    node_columns: dict = field(default_factory=dict)
    rim_range: tuple | None = None
    invert_range: tuple | None = None
    default_shape: str | None = None
    abandoned: str | None = None
    abandoned_values: frozenset = frozenset()
    null_values: frozenset = frozenset()
    null_nodes: frozenset = frozenset()
    nonpositive_inverts_missing: bool = False
    length_unit: str = "m"
    invert_unit: str = "m"
    size_unit: str = "m"
    crs: str | None = None
    design: dict = field(default_factory=dict)

    def _reprojected(self, data):
        """``data`` in :attr:`crs`, if set (data with no CRS are taken to be in it already)."""
        if self.crs is None:
            return data
        if data.crs is None:
            return data.set_crs(self.crs)
        return data.to_crs(self.crs)

    def __post_init__(self):
        unknown = set(self.columns) - set(CANONICAL_COLUMNS)
        derived = set(self.columns) & {"shape_code", "abandoned"}
        if unknown or derived:
            raise ValueError(
                f"Schema {self.name!r} maps columns that cannot be mapped: "
                f"{sorted(unknown | derived)}. Canonical columns: {list(CANONICAL_COLUMNS)}."
            )
        if self.default_shape is not None and "shape" in self.columns:
            raise ValueError(f"Schema {self.name!r} maps shape and also sets default_shape.")
        missing = [c for c in REQUIRED_COLUMNS if c not in self.columns]
        if missing:
            raise ValueError(f"Schema {self.name!r} does not map required columns: {missing}.")
        unknown = set(self.node_columns) - set(NODE_COLUMNS)
        if unknown:
            raise ValueError(
                f"Schema {self.name!r} maps node columns that cannot be mapped: "
                f"{sorted(unknown)}. Canonical node columns: {list(NODE_COLUMNS)}."
            )
        unknown = set(self.design) - set(DESIGN_RULES)
        if unknown:
            raise ValueError(
                f"Schema {self.name!r} has unknown design rules: {sorted(unknown)}. "
                f"Known rules: {list(DESIGN_RULES)}."
            )
        for rule, params in self.design.items():
            try:
                check_design_rule(rule, params)
            except (TypeError, ValueError) as error:
                raise ValueError(f"Schema {self.name!r}: design rule {rule}: {error}") from None
        for unit in (self.length_unit, self.invert_unit, self.size_unit):
            if unit not in UNITS:
                raise ValueError(f"Unknown unit {unit!r}; use one of {list(UNITS)}.")

    def __hash__(self):
        # The generated hash would hash the dict fields and fail. Hash the scalar fields only:
        # equal schemas have equal scalars, so this stays consistent with the generated __eq__.
        return hash(
            (
                self.name,
                self.default_shape,
                self.abandoned,
                self.length_unit,
                self.invert_unit,
                self.size_unit,
                self.crs,
            )
        )

    @classmethod
    def from_toml(cls, path):
        """Read a schema from a TOML file.

        The top-level keys are the parameters of :class:`Schema`, and ``columns``, ``shapes`` and
        ``flow_types`` are tables; see ``swmmbuilder/schemas/*.toml`` for examples.

        Parameters
        ----------
        path : str, Path or importlib.resources.abc.Traversable
            The TOML file.
        """
        path = Path(path) if isinstance(path, str) else path
        with path.open("rb") as file:
            values = tomllib.load(file)
        unknown = set(values) - {f.name for f in fields(cls)}
        if unknown:
            raise ValueError(f"Unknown keys in schema file {path}: {sorted(unknown)}.")
        for key in ("abandoned_values", "null_values", "null_nodes"):
            if key in values:
                values[key] = frozenset(values[key])
        for key in ("rim_range", "invert_range"):
            if key in values:
                values[key] = tuple(values[key])
        return cls(**values)

    def _lengths(self, data, lengths):
        """Pipe lengths in m; missing ones are computed from the (projected) geometry."""
        lengths = pd.to_numeric(lengths, errors="coerce") / UNITS[self.length_unit]
        missing = lengths.isna()
        if not missing.any():
            return lengths
        if data.crs is None or not data.crs.is_projected:
            raise ValueError(
                f"{int(missing.sum())} pipes have no length, and it cannot be computed from the "
                f"geometry, whose CRS ({data.crs}) is not projected. Set the schema's crs to a "
                "projected CRS in metres (e.g. the local UTM or MTM zone)."
            )
        metres = data.crs.axis_info[0].unit_conversion_factor
        return lengths.where(~missing, data.geometry.length * metres)

    def flow_type(self, code):
        """Canonical flow type of a source code (codes not in ``flow_types`` pass through)."""
        return self.flow_types.get(code, code)

    def node_type(self, code):
        """Canonical node type of a source code (codes not in ``node_types`` pass through)."""
        return self.node_types.get(code, code)

    def normalize(self, data):
        """Return a copy of ``data`` with the canonical columns.

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Pipes as read from the GIS file.

        Returns
        -------
        geopandas.GeoDataFrame
            The canonical columns, followed by the unmapped source columns and the geometry.

        Raises
        ------
        ValueError
            If a required source column is missing, ``pipe_id`` is not unique, or lengths must be
            computed from a geometry that is not in a projected CRS.

        Warns
        -----
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
        else:
            # Missing IDs (and null values) are not duplicates: Network drops those pipes and
            # counts them under removed['invalid']['no pipe_id'].
            ids = data[sources["pipe_id"]]
            if self.null_values:
                ids = ids.where(~ids.isin(self.null_values))
            ids = ids[ids.notna()]
            duplicated = ids[ids.duplicated()]
            if len(duplicated):
                values = sorted(set(duplicated), key=str)
                problems.append(f"duplicate {sources['pipe_id']} values {values}")
        if problems:
            raise ValueError(
                f"Cannot load the data with schema {self.name!r}: {'; '.join(problems)}."
            )

        data = self._reprojected(data)

        def source(canonical):
            column = sources.get(canonical)
            if column is None or column not in data.columns:
                return pd.Series(pd.NA, index=data.index, dtype=object)
            values = data[column]
            return values.where(~values.isin(self.null_values)) if self.null_values else values

        out = {}
        out["pipe_id"] = source("pipe_id")
        for node in ("from_node", "to_node"):
            out[node] = source(node).where(~source(node).isin(self.null_nodes))
        out["length"] = self._lengths(data, source("length"))
        bounds = self.invert_range
        for column in INVERT_COLUMNS:
            inverts = self._inverts(source(column))
            if bounds is not None:
                inverts = inverts.where(inverts.isna() | inverts.between(*bounds))
            out[column] = inverts
        codes = source("shape")
        out["shape"] = codes.map(self.shapes)
        if self.default_shape is not None:
            out["shape"] = pd.Series(self.default_shape, index=data.index, dtype=object)
        out["shape_code"] = codes
        for column in SIZE_COLUMNS:
            out[column] = pd.to_numeric(source(column), errors="coerce") / UNITS[self.size_unit]
        out["flow_type"] = source("flow_type").map(self.flow_type)
        out["material"] = source("material")
        if self.abandoned is not None and self.abandoned in data.columns:
            flag = data[self.abandoned]
            numeric = pd.to_numeric(flag, errors="coerce")
            flagged = flag.isin(self.abandoned_values) | numeric.isin(self.abandoned_values)
            fallback = flagged.map({True: "abandoned", False: "built"})
        else:
            fallback = "built"
        status = source("status").map(lambda code: self.statuses.get(code, code))
        out["status"] = status.where(status.notna(), fallback)
        out["abandoned"] = out["status"] == "abandoned"
        out["structure_type"] = source("structure_type").map(
            lambda code: self.structure_types.get(code, code)
        )
        for column in ("us_node_type", "ds_node_type"):
            out[column] = source(column).map(self.node_type)

        unknown = codes[codes.notna() & (codes.astype(str).str.strip() != "")]
        unknown = unknown[~unknown.isin(self.shapes)]
        if len(unknown):
            counts = unknown.astype(str).value_counts().sort_index().to_dict()
            warnings.warn(
                f"Shape codes not in schema {self.name!r} (pipes with them cannot be written to "
                f".inp): {counts}. Known codes: {list(self.shapes)}.",
                DataWarning,
                stacklevel=3,
            )

        geometry = data.geometry.name
        mapped = set(sources.values()) | set(CANONICAL_COLUMNS) | {geometry}
        rest = [c for c in data.columns if c not in mapped]
        frame = pd.concat([pd.DataFrame(out, index=data.index), data[[*rest, geometry]]], axis=1)
        return gpd.GeoDataFrame(frame, geometry=geometry, crs=data.crs)

    def _inverts(self, values):
        """Inverts in m from source values; non-numeric (and, with
        ``nonpositive_inverts_missing``, non-positive) values are NaN."""
        inverts = pd.to_numeric(values, errors="coerce").astype(float) / UNITS[self.invert_unit]
        if self.nonpositive_inverts_missing:
            inverts = inverts.where(inverts > 0)
        return inverts

    def implausible_inverts(self, data) -> pd.DataFrame:
        """The inverts of pipe data (as read, before :meth:`normalize`) outside
        ``invert_range``, which :meth:`normalize` treats as missing: one row per value,
        with ``pipe_id``, ``column`` (canonical), ``source_column``, and ``value`` (in m)."""
        columns = ["pipe_id", "column", "source_column", "value"]
        bounds = self.invert_range
        if bounds is None:
            return pd.DataFrame(columns=columns)
        rows = []
        for column in INVERT_COLUMNS:
            source = self.columns.get(column)
            if source is None or source not in data.columns:
                continue
            values = data[source]
            if self.null_values:
                values = values.where(~values.isin(self.null_values))
            inverts = self._inverts(values)
            bad = inverts.notna() & ~inverts.between(*bounds)
            for pipe, value in zip(
                data.loc[bad, self.columns["pipe_id"]], inverts[bad], strict=True
            ):
                rows.append((pipe, column, source, float(value)))
        return pd.DataFrame(rows, columns=columns)

    def normalize_nodes(self, data, node_type=None):
        """Return the canonical node columns of a node layer.

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            Node points as read from the GIS file.
        node_type : str, optional
            Canonical type of every node of the layer (e.g. ``'maintenance_hole'``); by default
            it is read from the ``node_type`` source column through ``node_types``.

        Returns
        -------
        geopandas.GeoDataFrame
            The columns of :data:`NODE_COLUMNS` and the geometry.

        Raises
        ------
        ValueError
            If the schema maps no ``node_id`` column or it is not in the data.
        """
        candidates = {
            c: [s] if isinstance(s, str) else list(s) for c, s in self.node_columns.items()
        }
        id_column = next((s for s in candidates.get("node_id", []) if s in data.columns), None)
        if id_column is None:
            raise ValueError(
                f"Cannot load the node layer with schema {self.name!r}: its node_id column "
                f"{candidates.get('node_id')} is not in the data."
            )
        data = self._reprojected(data)

        def source(canonical):
            column = next((s for s in candidates.get(canonical, []) if s in data.columns), None)
            if column is None:
                return pd.Series(pd.NA, index=data.index, dtype=object)
            values = data[column]
            return values.where(~values.isin(self.null_values)) if self.null_values else values

        def elevation(canonical):
            values = pd.to_numeric(source(canonical), errors="coerce") / UNITS[self.invert_unit]
            return values.where(values > 0) if self.nonpositive_inverts_missing else values

        out = {}
        out["node_id"] = source("node_id").where(~source("node_id").isin(self.null_nodes))
        if node_type is None:
            out["node_type"] = source("node_type").map(self.node_type).fillna("other")
        else:
            out["node_type"] = pd.Series(node_type, index=data.index)
        out["name"] = source("name")
        out["invert"] = elevation("invert")
        rim = elevation("rim_elevation")
        if self.rim_range is not None:
            low, high = self.rim_range
            rim = rim.where(rim.between(low, high))
        out["rim_elevation"] = rim
        out["depth"] = pd.to_numeric(source("depth"), errors="coerce") / UNITS[self.invert_unit]
        for column in NODE_BOOL_COLUMNS:
            out[column] = _booleans(source(column))
        for column in ("pipes_in", "pipes_out"):
            out[column] = pd.to_numeric(source(column), errors="coerce").astype("Int64")
        out["access"] = source("access")
        geometry = data.geometry.name
        out = pd.DataFrame(out, index=data.index)[list(NODE_COLUMNS)]
        frame = pd.concat([out, data[[geometry]]], axis=1)
        return gpd.GeoDataFrame(frame, geometry=geometry, crs=data.crs)


def _booleans(values):
    """Nullable booleans from flags such as 1/0, 'Y'/'N' or True/False (NA if unrecognised)."""

    def flag(value):
        if pd.isna(value):
            return pd.NA
        if not isinstance(value, str):
            return bool(value)
        text = value.strip().upper()
        return True if text in TRUE_STRINGS else False if text in FALSE_STRINGS else pd.NA

    return values.map(flag).astype("boolean")


def _preset(name):
    return Schema.from_toml(resources.files("swmmbuilder") / "schemas" / f"{name}.toml")


def as_schema(schema: Schema | StrPath) -> Schema:
    """A :class:`Schema`, or one read from a TOML file path."""
    return schema if isinstance(schema, Schema) else Schema.from_toml(schema)


# London, ON (the default). See swmmbuilder/schemas/london_on.toml.
LONDON_ON = _preset("london_on")

# The other presets are read on first use (module __getattr__), so that importing swmmbuilder does
# not read and validate every city's schema: {name: (TOML file stem, description)}.
LAZY_PRESETS = {
    "TORONTO_ON": ("toronto_on", "Toronto, ON, Sewer Gravity Mains"),
    "BARRIE_ON": ("barrie_on", "Barrie, ON, Sanitary Pipes"),
    "KITCHENER_ON": ("kitchener_on", "Kitchener, ON, Sanitary Pipes"),
    "PETERBOROUGH_ON": ("peterborough_on", "Peterborough, ON, Sanitary and Storm Sewer Pipes"),
    "SARNIA_ON": ("sarnia_on", "Sarnia, ON, Sanitary and Storm Sewers"),
    "WATERLOO_ON": ("waterloo_on", "Waterloo, ON, Sanitary Sewer Gravity Mains"),
    "WINDSOR_ON": ("windsor_on", "Windsor, ON, Sanitary, Storm and Combined Sewers"),
}


def __getattr__(name: str) -> Schema:
    """The schema presets other than LONDON_ON, read from swmmbuilder/schemas on first use."""
    if name in LAZY_PRESETS:
        schema = globals()[name] = _preset(LAZY_PRESETS[name][0])
        return schema
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted([*globals(), *LAZY_PRESETS])
