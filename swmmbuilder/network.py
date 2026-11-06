"""A city's sewer network: its pipes and node layers, loaded through a
:class:`~swmmbuilder.Schema`, and the flow graph between their nodes."""

from __future__ import annotations

import logging
import re
import warnings
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely import Point

from .corrections import REPORT_COLUMNS, PipeCorrections
from .diagnostics import DataWarning
from .elevation import RimEstimate, dem_name, sample_dem
from .gaps import lowest_inverts
from .geometry import is_disjoint, line_coords, node_positions, reverse_line, to_geodataframe
from .schema import LONDON_ON, NODE_COLUMNS, as_schema
from .subnetwork import Subnetwork, metric_positions
from .units import _is_real

if TYPE_CHECKING:
    from collections.abc import Callable, Collection, Mapping
    from typing import Any

    from ._types import NodeId, StrPath, Values
    from .schema import Schema

logger = logging.getLogger(__name__)

# Canonical flow types that Network keeps by default: sanitary sewers only.
DEFAULT_FLOW_TYPES = frozenset({"sanitary"})

# Canonical pipe statuses that Network keeps by default: built pipes only.
DEFAULT_STATUSES = frozenset({"built"})

# Node layers, in order of precedence when a node ID is in more than one, and the node type of
# their nodes (None: read from the layer's node_type column).
NODE_LAYERS = {
    "maintenance_holes": "maintenance_hole",
    "pump_stations": "pump_station",
    "outfalls": "outfall",
    "other_nodes": None,
}

# The columns of Network.split_report.
SPLIT_REPORT_COLUMNS = (
    "node_id",
    "new_id",
    "pipes",
    "separation",
    "layer",
    "point_distance",
    "attributes",
)


def strip_leading_zeros(node_id: NodeId) -> NodeId:
    """Node ID with the leading zeros of its number removed, e.g. 'OK0329' -> 'OK329'.

    A ``normalize_ids`` hook for :class:`Network`: newer London exports drop the zeros that older
    ones (and node lists saved from them) have. IDs that are not strings pass through.
    """
    if not isinstance(node_id, str):
        return node_id
    return re.sub(r"^([A-Za-z]+)0+(?=\d)", r"\1", node_id)


def _pipe_records(data):
    """Yield (pipe_id, from_node, to_node, length, structure_type) for every pipe with two node IDs.

    Pipes with a missing node ID are skipped.
    """
    for oid, up, down, length, structure in zip(
        data["pipe_id"],
        data["from_node"],
        data["to_node"],
        data["length"],
        data["structure_type"],
        strict=True,
    ):
        if pd.isna(up) or pd.isna(down):
            continue
        yield oid, up, down, length, None if pd.isna(structure) else structure


def _ids(ids, limit=20):
    """Text listing up to ``limit`` IDs."""
    ids = [str(i) for i in ids]
    more = f" and {len(ids) - limit} more" if len(ids) > limit else ""
    return ", ".join(ids[:limit]) + more


def _nullable(dtype):
    """A dtype that holds ``dtype``'s values and missing values: NumPy ints and bools become
    pandas' nullable Int64 and boolean (a node layer's column, e.g. an all-numeric name, is
    missing for nodes that are not in the layer)."""
    if pd.api.types.is_bool_dtype(dtype):
        return "boolean"
    if pd.api.types.is_integer_dtype(dtype):
        return "Int64"
    return dtype


def _read(layer):
    return layer if isinstance(layer, gpd.GeoDataFrame) else gpd.read_file(Path(layer))


def _letters(i):
    """The suffix of the ``i``-th structure of a reused ID: 'a', 'b', ..., 'z', 'aa', 'ab', ..."""
    suffix = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        suffix = chr(ord("a") + r) + suffix
    return suffix


def _clusters(points, tolerance):
    """Single-linkage clusters of ``points`` (an (n, 2) array, in m): lists of point indices, each
    cluster's points joined by a chain of steps of at most ``tolerance``."""
    parent = list(range(len(points)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    distance = np.hypot(*(points[:, None, :] - points[None, :, :]).transpose(2, 0, 1))
    for i, j in zip(*np.nonzero(np.triu(distance <= tolerance, 1)), strict=True):
        parent[root(int(i))] = root(int(j))
    clusters: dict[int, list[int]] = {}
    for i in range(len(points)):
        clusters.setdefault(root(i), []).append(i)
    return list(clusters.values())


class Network:
    """A city's sewer network.

    Parameters
    ----------
    source : str, Path or geopandas.GeoDataFrame
        GIS file of the pipes (anything :func:`geopandas.read_file` reads), or the pipes
        themselves, in the source schema.

    flow_types : iterable of str or None
        Flow types to keep (default ``{'sanitary'}``). Source codes are translated by the schema,
        so for London ``{'SAN'}`` is the same. Every other pipe, including those with no flow type,
        is dropped and the counts are logged. ``None`` keeps every pipe.

    keep_abandoned : bool
        Also keep abandoned pipes (adds ``'abandoned'`` to ``status``). Pipes with a node ID the
        schema treats as no node (for London, 'Abandoned') are dropped from the graph either way.

    schema : Schema, str or Path
        How the GIS columns map to the canonical columns: a :class:`Schema` (default
        :data:`LONDON_ON`; see :mod:`swmmbuilder.schema` for the other presets) or the path to a
        schema TOML file.

    status : iterable of str or None
        Canonical pipe statuses to keep (default ``{'built'}``; the others are ``'abandoned'``,
        ``'proposed'`` and ``'removed'``). ``None`` keeps every pipe.

    include_overflows : bool
        Keep pipes whose structure type is ``'overflow'`` (default True). Overflow pipes carry no
        dry-weather flow and make spurious flow splits and dead ends, so pass False for a
        dry-weather network.

    maintenance_holes, pump_stations, outfalls, other_nodes : str, Path or GeoDataFrame, optional
        Node layers (points). Their nodes get the node type of the layer (``other_nodes`` reads it
        from the schema's ``node_type`` column, e.g. caps, tees and chambers); a node ID in more
        than one layer takes the first of them in this order. Nodes in no layer take the type of
        the pipes' inventory types, if the schema maps them.

    normalize_ids : callable, optional
        Applied to every node ID of the pipes and node layers (missing IDs are skipped), and to
        the node IDs passed to :meth:`catchment`, :meth:`subnetwork`,
        :meth:`Subnetwork.from_node_list` and the node-keyed options of :meth:`Subnetwork.to_inp`
        (see :meth:`Subnetwork.match_node`), e.g. :func:`strip_leading_zeros`.

    snap_tolerance : float, optional
        If set, :meth:`snap_nodes` is called with it (and ``max_invert_difference``) on loading,
        and its report is kept as ``snap_report``.

    max_invert_difference : float, optional
        See :meth:`snap_nodes`.

    split_tolerance : float, optional
        If set, node IDs used for more than one structure are split on loading, before any
        snapping: the pipe ends that carry one ID are grouped into structures, and two ends
        more than ``split_tolerance`` apart (in m) belong to different structures unless a chain
        of ends at most ``split_tolerance`` apart joins them. An ID with more than one structure
        becomes one node per structure, named by appending 'a', 'b', 'c', ... to the ID in the
        order of each structure's smallest pipe ID, so that rebuilds give the same names. A
        pipe whose geometry runs from its downstream node to its upstream node (each end is
        more than ``split_tolerance`` from the other points of its own ID, and within it of the
        other ID's) has its ends matched the other way round, so that it does not split its
        nodes, and its geometry is reversed (it is listed in ``reversed_pipes``), so that node
        positions and the .inp vertices run from its upstream node. Each node-layer point with a split ID goes to the nearest of its structures,
        which takes the point's position and attributes (node type, invert, rim); a structure
        with no point takes them from its pipes, as a node in no layer does. An ID whose pipes
        include one with no geometry is not split. Every split is logged as a warning and
        listed in ``split_report``. Must exceed ``snap_tolerance``, so that snapping cannot
        merge the structures again.

    corrections : PipeCorrections, dict, str or Path, optional
        Documented corrections to pipes: a :class:`PipeCorrections`, a mapping with its fields
        as keys (``drop``: {pipe ID: reason}; ``set``: {pipe ID: {canonical column: value, ...,
        'reason': str}}; ``swap_nodes``: {pipe ID: reason}), or the path to a TOML file of them
        (see :meth:`PipeCorrections.from_toml`). Pipe IDs are those of the source data, matched
        exactly or by text. They are applied after the schema has normalized the data and
        ``normalize_ids`` has been applied to the pipes' node IDs, and before any pipe is
        dropped by the flow-type, status, and overflow filters (or as invalid), so that e.g.
        ``set={pid: {'status': 'built'}}`` keeps a pipe the data mark as proposed. Each one is
        logged at INFO with its reason and listed in ``corrections_report``. ``swap_nodes``
        swaps only the node IDs (inverts and geometry stay as recorded; see
        :class:`PipeCorrections`). Raises ValueError naming every pipe ID that is not in the
        data.

    Attributes
    ----------
    full_data : geopandas.GeoDataFrame
        The pipes, with the canonical columns of :mod:`swmmbuilder.schema`.

    schema : Schema
        The schema the data were loaded with.

    graph : networkx.MultiDiGraph
        One edge per pipe, in the direction of flow (from_node -> to_node), keyed by pipe_id so
        that parallel pipes are kept. Each edge has the attributes ``pipe`` (pipe_id), ``length``
        and ``structure_type``; each node has the non-missing columns of ``nodes``. Nodes and
        edges are inserted in sorted order, so traversals do not depend on the order of the file.

    nodes : geopandas.GeoDataFrame
        One row per node of the graph, indexed by node ID, with the canonical node columns (see
        :mod:`swmmbuilder.schema`) and the node's point (from its node layer, else the end of a
        pipe). ``invert`` is the node layer's, else the lowest invert of the pipes connected to
        the node; ``rim_elevation`` is the node layer's, else ``invert + depth``.

    aliases : dict
        {node ID: node ID it was merged into by :meth:`snap_nodes`}. Every method that takes node
        IDs accepts either.

    removed : dict
        What was dropped on loading: ``{'flow_type': {flow_type: count}, 'status': {status:
        count}, 'structure_type': {structure_type: count}, 'invalid': {problem: count}}``. Pipes
        with no flow type are counted under ``'None'``; ``'invalid'`` counts the pipes with
        ``'no pipe_id'`` or ``'no length'`` (none in the data, and no geometry to compute it
        from), which cannot be used. Pipes with no geometry but a length are kept.
        ``'corrections'`` is ``{'dropped': count}`` for the pipes dropped by ``corrections``
        (empty if none were).

    corrections_report : pandas.DataFrame
        One row per correction applied: ``pipe_id`` (as in the data), ``action`` (``'drop'``,
        ``'set'`` or ``'swap_nodes'``), ``column`` (the column set, else missing), ``old`` and
        ``new`` (the value before and after; for ``'swap_nodes'``, ``'from_node -> to_node'``;
        missing for ``'drop'``), and ``reason``. Empty, with these columns, if there were none.

    snap_report : pandas.DataFrame or None
        The report of :meth:`snap_nodes` if ``snap_tolerance`` was given.

    splits : dict
        {node ID: [the node IDs it was split into]} by ``split_tolerance``. The original IDs are
        no longer nodes: :meth:`check_nodes` names their split IDs.

    split_report : pandas.DataFrame or None
        If ``split_tolerance`` was given, one row per node made by splitting a reused ID:
        ``node_id`` (the ID in the data), ``new_id``, ``pipes`` (the pipe IDs of the structure),
        ``separation`` (distance to the nearest end of another structure of the ID, in m),
        ``layer`` (the node layer whose point it took, or None), ``point_distance`` (from that
        point to the nearest of its pipe ends, in m) and ``attributes`` ('node layer', or
        'pipes' for a structure with no point).

    reversed_pipes : list
        IDs of the pipes drawn from their downstream node to their upstream node, whose geometry
        was reversed (found with ``split_tolerance``; empty without it).

    invert_report : pandas.DataFrame
        The pipe and drop-structure inverts of the pipes kept that are outside the schema's plausible range
        (``Schema.invert_range``), treated as missing: ``pipe_id``, ``column``,
        ``source_column``, and ``value`` (in m), as recorded.

    version : int
        Bumped each time the graph is rebuilt (by :meth:`snap_nodes`); subnetworks made from
        an older version raise when used.
    """

    def __init__(
        self,
        source: StrPath | gpd.GeoDataFrame,
        flow_types: Iterable[str] | None = DEFAULT_FLOW_TYPES,
        keep_abandoned: bool = False,
        schema: Schema | StrPath = LONDON_ON,
        status: Iterable[str] | None = DEFAULT_STATUSES,
        include_overflows: bool = True,
        maintenance_holes: StrPath | gpd.GeoDataFrame | None = None,
        pump_stations: StrPath | gpd.GeoDataFrame | None = None,
        outfalls: StrPath | gpd.GeoDataFrame | None = None,
        other_nodes: StrPath | gpd.GeoDataFrame | None = None,
        normalize_ids: Callable[[Any], Any] | None = None,
        snap_tolerance: float | None = None,
        max_invert_difference: float | None = None,
        split_tolerance: float | None = None,
        corrections: PipeCorrections | Mapping[str, Mapping[Any, Any]] | StrPath | None = None,
    ) -> None:
        self.schema = schema = as_schema(schema)
        self.normalize_ids = normalize_ids
        self.aliases: dict[NodeId, NodeId] = {}
        self.snap_report: pd.DataFrame | None = None
        self.splits: dict[NodeId, list[NodeId]] = {}
        self._split_points: dict[NodeId, Point] = {}
        self.split_report: pd.DataFrame | None = None
        self.reversed_pipes: list = []
        if (
            split_tolerance is not None
            and snap_tolerance is not None
            and not split_tolerance > snap_tolerance
        ):
            raise ValueError(
                f"split_tolerance ({split_tolerance} m) must exceed snap_tolerance "
                f"({snap_tolerance} m)."
            )
        raw = _read(source)
        self.full_data = schema.normalize(raw)
        implausible = schema.implausible_inverts(raw)
        if normalize_ids is not None:
            for column in ("from_node", "to_node"):
                self.full_data[column] = self.full_data[column].map(
                    normalize_ids, na_action="ignore"
                )
        self.removed: dict[str, dict[str, int]] = {
            "flow_type": {},
            "status": {},
            "structure_type": {},
            "invalid": {},
            "corrections": {},
        }
        self.corrections_report = self._apply_corrections(corrections)
        self._drop_invalid()

        if flow_types is not None and "flow_type" in schema.columns:
            flow_types = {schema.flow_type(f) for f in flow_types}
            keep = self.full_data["flow_type"].isin(flow_types)
            self._drop(keep, "flow_type", f"keeping {sorted(flow_types)}")

        if status is not None:
            status = set(status) | ({"abandoned"} if keep_abandoned else set())
            keep = self.full_data["status"].isin(status)
            self._drop(keep, "status", f"keeping {sorted(status)}")

        if not include_overflows:
            keep = self.full_data["structure_type"] != "overflow"
            self._drop(keep, "structure_type", "include_overflows=False")

        # Report only the pipes kept (not, say, the storm sewers filtered out above).
        kept = set(self.full_data["pipe_id"].astype(str))
        self.invert_report = implausible[implausible["pipe_id"].astype(str).isin(kept)]
        # ... and not the inverts that corrections replaced.
        corrected = self.corrections_report[self.corrections_report["action"] == "set"]
        corrected_cells = set(
            zip(corrected["pipe_id"].astype(str), corrected["column"], strict=True)
        )
        self.invert_report = self.invert_report[
            [
                (str(p), c) not in corrected_cells
                for p, c in zip(
                    self.invert_report["pipe_id"], self.invert_report["column"], strict=True
                )
            ]
        ]
        self.invert_report = self.invert_report.reset_index(drop=True)
        if len(self.invert_report) and schema.invert_range is not None:
            low, high = schema.invert_range
            for row in self.invert_report.itertuples():
                logger.info(
                    "Pipe %s: %s %s outside %g-%g m, treated as missing",
                    row.pipe_id,
                    row.source_column,
                    row.value,
                    low,
                    high,
                )
            warnings.warn(
                f"{len(self.invert_report)} pipe inverts outside {low:g}-{high:g} m treated as "
                "missing (see Network.invert_report)",
                DataWarning,
                stacklevel=2,
            )

        self._layers = self._read_layers(
            dict(
                maintenance_holes=maintenance_holes,
                pump_stations=pump_stations,
                outfalls=outfalls,
                other_nodes=other_nodes,
            )
        )
        if split_tolerance is not None:
            self.split_report = self._split_reused_ids(split_tolerance)
        self._build()
        if snap_tolerance is not None:
            self.snap_report = self.snap_nodes(snap_tolerance, max_invert_difference)

    def _drop_invalid(self):
        """Drop the pipes that cannot be used (no pipe_id, or no length from the data or the
        geometry), counting them in ``removed['invalid']``, and warn about pipes with no
        geometry or whose MultiLineString parts do not join."""
        data = self.full_data
        problems = {
            "no pipe_id": data["pipe_id"].isna(),
            "no length": data["pipe_id"].notna() & data["length"].isna(),
        }
        for problem, mask in problems.items():
            if mask.any():
                self.removed["invalid"][problem] = int(mask.sum())
                logger.warning(
                    "Removed %d pipes with %s%s",
                    int(mask.sum()),
                    problem,
                    "" if problem == "no pipe_id" else f": {_ids(data.loc[mask, 'pipe_id'])}",
                )
        self.full_data = data[~(problems["no pipe_id"] | problems["no length"])]
        geometry = self.full_data.geometry
        missing = geometry.isna() | geometry.is_empty
        if missing.any():
            logger.warning(
                "%d pipes have no geometry; they are kept (their length is from the data), but "
                "they give no node positions or vertices: %s",
                int(missing.sum()),
                _ids(self.full_data.loc[missing, "pipe_id"]),
            )
        disjoint = geometry[~missing].map(is_disjoint)
        if disjoint.any():
            logger.warning(
                "%d pipes are MultiLineStrings whose parts do not join; their parts are joined "
                "in order, with a straight segment across each gap: %s",
                int(disjoint.sum()),
                _ids(self.full_data.loc[disjoint[disjoint].index, "pipe_id"]),
            )

    def _drop(self, keep, column, reason):
        """Drop the pipes not in ``keep``, counting them by ``column`` in ``removed``."""
        dropped = self.full_data.loc[~keep, column].fillna("None").astype(str)
        self.removed[column] = dropped.value_counts().sort_index().to_dict()
        self.full_data = self.full_data[keep]
        for value, count in self.removed[column].items():
            logger.info(
                "Removed %d pipes with %s %s (%s)", count, column.replace("_", " "), value, reason
            )

    def _apply_corrections(self, corrections):
        """Apply ``corrections`` (see :class:`PipeCorrections`) to ``full_data``, logging each
        one, and return the corrections report."""
        if corrections is None:
            return pd.DataFrame(columns=REPORT_COLUMNS)
        corrections = PipeCorrections.coerce(corrections)
        data = self.full_data.copy()
        positions = {p: i for i, p in enumerate(data["pipe_id"]) if pd.notna(p)}
        by_text = {str(p): p for p in positions}

        missing = set()
        repeated = set()
        matched = {}
        for action in ("drop", "set", "swap_nodes"):
            entries = {}
            for key, value in getattr(corrections, action).items():
                pipe = key if key in positions else by_text.get(str(key))
                if pipe is None:
                    missing.add(str(key))
                elif pipe in entries:
                    repeated.add(str(pipe))
                else:
                    entries[pipe] = value
            matched[action] = entries
        if missing:
            raise ValueError(
                f"corrections: pipes not in the data: {', '.join(sorted(missing))}. Pipe IDs are "
                f"those of the source data ({self.schema.columns['pipe_id']!r})."
            )
        if repeated:
            raise ValueError(
                f"corrections: pipes given more than once: {', '.join(sorted(repeated))}."
            )

        rows = []
        for pipe, reason in matched["drop"].items():
            logger.info("Pipe %s dropped by corrections: %s", pipe, reason)
            rows.append((pipe, "drop", None, None, None, reason))
        column_of = data.columns.get_loc
        for pipe, values in matched["set"].items():
            i = positions[pipe]
            reason = values.get("reason")
            for column, value in values.items():
                if column == "reason":
                    continue
                if column in ("from_node", "to_node") and self.normalize_ids and pd.notna(value):
                    value = self.normalize_ids(value)
                old = data.iat[i, column_of(column)]
                data.iat[i, column_of(column)] = value
                if column == "status":
                    data.iat[i, column_of("abandoned")] = value == "abandoned"
                logger.info(
                    "Pipe %s: %s set from %r to %r by corrections%s",
                    pipe,
                    column,
                    old,
                    value,
                    f": {reason}" if reason else "",
                )
                rows.append((pipe, "set", column, old, value, reason))
        ends = [column_of("from_node"), column_of("to_node")]
        for pipe, reason in matched["swap_nodes"].items():
            i = positions[pipe]
            up, down = data.iat[i, ends[0]], data.iat[i, ends[1]]
            data.iat[i, ends[0]], data.iat[i, ends[1]] = down, up
            logger.info(
                "Pipe %s: node IDs swapped by corrections, %s -> %s is now %s -> %s: %s",
                pipe,
                up,
                down,
                down,
                up,
                reason,
            )
            rows.append((pipe, "swap_nodes", None, f"{up} -> {down}", f"{down} -> {up}", reason))

        if matched["drop"]:
            dropped = data["pipe_id"].isin(list(matched["drop"]))
            data = data[~dropped]
            self.removed["corrections"] = {"dropped": int(dropped.sum())}
        self.full_data = data
        return pd.DataFrame(rows, columns=REPORT_COLUMNS)

    def _read_layers(self, layers):
        """The canonical node layers, in order of precedence, each with a ``layer`` column."""
        frames = []
        for name, layer in layers.items():
            if layer is None:
                continue
            nodes = self.schema.normalize_nodes(_read(layer), NODE_LAYERS[name])
            if self.normalize_ids is not None:
                nodes["node_id"] = nodes["node_id"].map(self.normalize_ids, na_action="ignore")
            if self.full_data.crs is not None and nodes.crs != self.full_data.crs:
                nodes = nodes.to_crs(self.full_data.crs)
            frames.append(nodes.assign(layer=name))
            logger.info("Read %d nodes from the %s layer", len(nodes), name)
        if not frames:
            return None
        return pd.concat(frames, ignore_index=True)

    def _split_reused_ids(self, tolerance):
        """Split the node IDs used for more than one structure (see ``split_tolerance``),
        renaming them in ``full_data`` and the node layers, and return the report."""
        data = self.full_data
        if data.crs is None:
            raise ValueError("Cannot split reused node IDs: the data have no CRS.")
        geometry = data.geometry
        layers = self._layers
        if not data.crs.is_projected:
            utm = geometry.estimate_utm_crs()
            geometry = geometry.to_crs(utm)
            if layers is not None:
                layers = layers.to_crs(utm)
        metres = geometry.crs.axis_info[0].unit_conversion_factor

        ids = (data["from_node"].to_numpy(), data["to_node"].to_numpy())
        pipe_ids = data["pipe_id"].to_numpy()
        # ends[i] = [upstream end, downstream end] of pipe i (by position), in m, from its
        # geometry; ends_of[node] = [(i, side)], side 0 for a pipe leaving the node, 1 entering.
        ends: dict[int, list[np.ndarray]] = {}
        drawn: dict[int, list[tuple]] = {}  # the same ends in the data's coordinates
        ends_of: dict[NodeId, list[tuple[int, int]]] = {}
        no_geometry: set[NodeId] = set()
        for i, (g, g_data) in enumerate(zip(geometry, data.geometry, strict=True)):
            coords = line_coords(g)
            for side in (0, 1):
                node = ids[side][i]
                if pd.isna(node):
                    continue
                if not coords:
                    no_geometry.add(node)
                    continue
                ends_of.setdefault(node, []).append((i, side))
            if coords:
                ends[i] = [np.array(coords[0]) * metres, np.array(coords[-1]) * metres]
                coords = line_coords(g_data)
                drawn[i] = [coords[0], coords[-1]]
        points_of: dict[NodeId, list[tuple[int, np.ndarray]]] = {}
        if layers is not None:
            for j, (node, point) in enumerate(zip(layers["node_id"], layers.geometry, strict=True)):
                if pd.notna(node) and point is not None and not point.is_empty:
                    points_of.setdefault(node, []).append((j, np.array(point.coords[0]) * metres))

        def extent(node):
            xy = np.array([ends[i][side] for i, side in ends_of[node]])
            return float(np.hypot(*(xy.max(axis=0) - xy.min(axis=0))))

        candidates = {n for n in ends_of if len(ends_of[n]) > 1 and extent(n) > tolerance}

        def near(node, xy, pipe):
            """True if ``xy`` is within tolerance of a point of ``node`` other than ``pipe``'s."""
            others = [ends[i][side] for i, side in ends_of.get(node, ()) if i != pipe]
            others += [p for _, p in points_of.get(node, ())]
            return any(np.hypot(*(xy - p)) <= tolerance for p in others)

        # Pipes drawn from their downstream node to their upstream node: match their ends the
        # other way round, so that the drawing direction does not split their nodes.
        backwards = []
        for i in sorted({i for n in candidates for i, _ in ends_of[n]}):
            up, down = ids[0][i], ids[1][i]
            first, last = ends[i]
            if (
                up != down
                and not near(up, first, i)
                and not near(down, last, i)
                and near(up, last, i)
                and near(down, first, i)
            ):
                ends[i] = [last, first]
                drawn[i] = drawn[i][::-1]
                backwards.append(pipe_ids[i])
        if backwards:
            logger.info(
                "%d pipes are drawn from their downstream node to their upstream node; their ends "
                "were matched to their node IDs the other way round: %s",
                len(backwards),
                _ids(sorted(backwards)),
            )

        taken = {n for column in ids for n in column if pd.notna(n)} | set(points_of)
        renamed = [ids[0].copy(), ids[1].copy()]
        layer_ids = None if layers is None else self._layers["node_id"].to_numpy().copy()
        rows, clashes, unsplit = [], [], []
        for node in sorted(candidates, key=str):
            members = ends_of[node]
            xy = np.array([ends[i][side] for i, side in members])
            clusters = _clusters(xy, tolerance)
            if len(clusters) == 1:
                continue
            if node in no_geometry:
                unsplit.append(node)
                continue
            clusters.sort(
                key=lambda c: (
                    sorted(pipe_ids[members[k][0]] for k in c),
                    *xy[c].mean(axis=0).round(3),
                )
            )
            new_ids = [f"{node}{_letters(k)}" for k in range(len(clusters))]
            clashes += [n for n in new_ids if n in taken]
            # Each node-layer point goes to its nearest structure; the first in layer order wins.
            point: dict[int, tuple[str, float]] = {}
            for j, p in points_of.get(node, ()):
                gaps = [float(np.hypot(*(xy[c] - p).T).min()) for c in clusters]
                k = int(np.argmin(gaps))
                layer_ids[j] = new_ids[k]
                point.setdefault(k, (layers["layer"].iloc[j], gaps[k]))
            for k, cluster in enumerate(clusters):
                for m in cluster:
                    i, side = members[m]
                    renamed[side][i] = new_ids[k]
                others = [m for m in range(len(members)) if m not in cluster]
                separation = float(
                    np.hypot(*(xy[cluster][:, None, :] - xy[others][None, :, :]).T).min()
                )
                layer, distance = point.get(k, (None, np.nan))
                if layer is None:
                    # No point: the node is at its pipe end nearest the structure's centre.
                    centre = xy[cluster].mean(axis=0)
                    m = cluster[int(np.argmin(np.hypot(*(xy[cluster] - centre).T)))]
                    i, side = members[m]
                    self._split_points[new_ids[k]] = Point(drawn[i][side])
                rows.append(
                    (
                        node,
                        new_ids[k],
                        tuple(sorted({pipe_ids[members[m][0]] for m in cluster})),
                        separation,
                        layer,
                        distance,
                        "pipes" if layer is None else "node layer",
                    )
                )
            self.splits[node] = new_ids
            logger.warning(
                "Node ID %s is used for %d structures, at least %.0f m apart; split into %s",
                node,
                len(clusters),
                min(r[3] for r in rows[-len(clusters) :]),
                "; ".join(
                    f"{r[1]} (pipes {_ids(r[2])}; "
                    + (
                        f"{r[4]} point {r[5]:.1f} m away)"
                        if r[4] is not None
                        else "no node-layer point: attributes from its pipes)"
                    )
                    for r in rows[-len(clusters) :]
                ),
            )
        if clashes:
            raise ValueError(
                f"Cannot split reused node IDs: the new IDs are already in use: {_ids(clashes)}"
            )
        if unsplit:
            logger.warning(
                "%d node IDs have pipe ends more than %s m apart but a pipe with no geometry, "
                "so they are not split: %s",
                len(unsplit),
                tolerance,
                _ids(unsplit),
            )
        if self.splits:
            self.full_data = data.assign(from_node=renamed[0], to_node=renamed[1])
            if layer_ids is not None:
                self._layers = self._layers.assign(node_id=layer_ids)
        if backwards:
            # Reverse the geometry itself, so that everything reading it (node positions, the
            # .inp coordinates and vertices) starts at the upstream node (#121).
            flip = self.full_data["pipe_id"].isin(backwards).to_numpy()
            geometry = self.full_data.geometry.copy()
            geometry[flip] = [reverse_line(g) for g in geometry[flip]]
            self.full_data = self.full_data.set_geometry(geometry)
            self.reversed_pipes = sorted(backwards, key=str)
        logger.info(
            "Split %d reused node IDs into %d nodes (structures more than %s m apart)",
            len(self.splits),
            len(rows),
            tolerance,
        )
        return pd.DataFrame(rows, columns=list(SPLIT_REPORT_COLUMNS))

    def _build(self):
        """(Re)build ``graph`` and ``nodes`` from ``full_data`` and the node layers.

        Each rebuild bumps ``version``, so that subnetworks made before it can tell they are
        out of date (see :attr:`Subnetwork.network`).
        """
        self.version = getattr(self, "version", 0) + 1
        self.graph = self._build_graph(self.full_data)
        self.nodes = self._build_nodes()
        attrs = self.nodes.drop(columns="geometry").to_dict("index")
        nx.set_node_attributes(
            self.graph,
            {n: {k: v for k, v in a.items() if not pd.isna(v)} for n, a in attrs.items()},
        )

    @staticmethod
    def _build_graph(data):
        """Build the flow-direction MultiDiGraph of the pipes in ``data``."""
        records = sorted(_pipe_records(data), key=lambda r: (r[1], r[2], r[0]))
        G = nx.MultiDiGraph()
        G.add_nodes_from(sorted({n for _, up, down, _, _ in records for n in (up, down)}))
        for oid, up, down, length, structure in records:
            G.add_edge(up, down, key=oid, pipe=oid, length=length, structure_type=structure)
        return G

    def _build_nodes(self):
        """The ``nodes`` table of the graph nodes."""
        ids = list(self.graph.nodes)
        data = self.full_data
        columns = [c for c in NODE_COLUMNS if c != "node_id"]
        nodes = pd.DataFrame(index=pd.Index(ids, name="node_id"), columns=[*columns, "layer"])
        geometry = pd.Series(None, index=nodes.index, dtype=object)

        if self._layers is not None:
            layers = self._layers.assign(node_id=self._layers["node_id"].map(self.resolve))
            layers = layers[layers["node_id"].isin(set(ids))]
            duplicated = layers["node_id"].duplicated()
            if duplicated.any():
                logger.info(
                    "%d node IDs are in more than one node layer (or were merged); the first "
                    "layer's record is kept",
                    int(duplicated.sum()),
                )
            layers = layers[~duplicated].set_index("node_id")
            nodes.update(layers[[*columns, "layer"]])
            for column in columns:
                nodes[column] = nodes[column].astype(_nullable(layers[column].dtype))
            geometry.update(layers.geometry)

        inventory = pd.concat(
            [
                data[["from_node", "us_node_type"]].set_axis(["node", "type"], axis=1),
                data[["to_node", "ds_node_type"]].set_axis(["node", "type"], axis=1),
            ]
        ).dropna()
        inventory = inventory.drop_duplicates("node").set_index("node")["type"]
        nodes["node_type"] = nodes["node_type"].fillna(inventory.reindex(nodes.index))

        inverts = pd.Series(lowest_inverts(data, set(ids)), dtype=float)
        nodes["invert"] = nodes["invert"].astype(float).fillna(inverts.reindex(nodes.index))
        nodes["rim_elevation"] = (
            nodes["rim_elevation"]
            .astype(float)
            .fillna(nodes["invert"] + nodes["depth"].astype(float))
        )
        split_points = pd.Series(self._split_points, dtype=object).reindex(nodes.index)
        geometry = geometry.where(geometry.notna(), split_points)
        missing = geometry.isna()
        if missing.any():
            ends = pipe_end_points(data, set(nodes.index[missing]))
            geometry = geometry.where(~missing, pd.Series(ends).reindex(nodes.index))
        geometry = gpd.GeoSeries(geometry.where(geometry.notna(), None), crs=data.crs)
        return gpd.GeoDataFrame(nodes, geometry=geometry)

    def resolve(self, node: NodeId) -> NodeId:
        """The node ID that ``node`` refers to: normalized and with merged IDs replaced."""
        if self.normalize_ids is not None:
            node = self.normalize_ids(node)
        return self.aliases.get(node, node)

    def snap_nodes(
        self, tolerance: float = 0.5, max_invert_difference: float | None = None
    ) -> pd.DataFrame:
        """Merge node IDs whose points are within ``tolerance`` of each other.

        The same structure sometimes has two IDs at one point, which disconnects the network
        there. Each group of nodes within ``tolerance`` (in m; geographic coordinates are
        projected to the local UTM zone) is merged into one ID, preferring an ID in the
        maintenance-hole layer, then one in another node layer, then the smallest. Pairs whose
        inverts differ by more than ``max_invert_difference`` (pipes that cross rather than join),
        pairs already joined by a pipe, and pairs split from one reused ID (see
        ``split_tolerance``) are not merged. Since merging is transitive, both
        rules hold for each merged group as a whole: a pair is not merged if that would put two
        nodes joined by a pipe, or two inverts more than ``max_invert_difference`` apart, into
        one group. Merged IDs are kept in ``aliases``, so either ID can be used.

        Snapping changes the network in place: subnetworks made before a snap that merges
        nodes are out of date and raise when used (rebuild them with :meth:`subnetwork`). Prefer
        ``snap_tolerance`` when loading, which snaps before any subnetwork exists.

        Parameters
        ----------
        tolerance : float
            Largest distance between merged nodes, in m.

        max_invert_difference : float, optional
            Largest difference between the inverts of merged nodes, in m.

        Returns
        -------
        pandas.DataFrame
            One row per pair of nodes within ``tolerance``: ``node_1``, ``node_2``, ``distance``,
            ``invert_1``, ``invert_2``, ``snapped`` (bool), ``kept`` (the ID they were merged
            into) and ``reason`` (why a pair was not merged).

        Raises
        ------
        ValueError
            If the data have no CRS.
        """
        points = self.nodes.geometry
        if points.crs is None:
            raise ValueError("Cannot snap nodes: the data have no CRS.")
        if not points.crs.is_projected:
            points = points.to_crs(points.estimate_utm_crs())
        metres = points.crs.axis_info[0].unit_conversion_factor
        present = points[points.notna()]
        left, right = present.sindex.query(
            present, predicate="dwithin", distance=tolerance / metres
        )
        ids = list(present.index)
        inverts = self.nodes["invert"]
        pairs: list[tuple] = []
        parent: dict[NodeId, NodeId] = {}
        members: dict[NodeId, set[NodeId]] = {}

        def root(node):
            while parent.get(node, node) != node:
                node = parent[node]
            return node

        def group(node):
            return members.get(root(node), {node})

        split_from = {new: node for node, new_ids in self.splits.items() for new in new_ids}

        def neighbours(node):
            return set(self.graph.successors(node)) | set(self.graph.predecessors(node))

        # The checks apply to the whole merged group, not only to the pair: merging is
        # transitive, so A~B and B~C would otherwise merge A and C even if a pipe joins them (a
        # self-loop) or their inverts differ by more than max_invert_difference.
        for i, j in sorted({(i, j) for i, j in zip(left, right, strict=True) if i < j}):
            a, b = sorted((ids[i], ids[j]))
            difference = abs(inverts[a] - inverts[b])
            reason = None
            if self.graph.has_edge(a, b) or self.graph.has_edge(b, a):
                reason = "joined by a pipe"
            elif a in split_from and split_from[a] == split_from.get(b):
                reason = "split from one node ID"
            elif max_invert_difference is not None and not difference <= max_invert_difference:
                reason = f"inverts differ by more than {max_invert_difference} m"
            elif root(a) != root(b):
                union = group(a) | group(b)
                spread = inverts[sorted(union)]
                if any(neighbours(n) & union for n in group(a)):
                    reason = "would merge nodes joined by a pipe"
                elif max_invert_difference is not None and (
                    spread.isna().any() or not spread.max() - spread.min() <= max_invert_difference
                ):
                    reason = (
                        "would merge nodes whose inverts differ by more than "
                        f"{max_invert_difference} m"
                    )
                else:
                    ra, rb = root(a), root(b)
                    parent[rb] = ra
                    members[ra] = union
                    members.pop(rb, None)
            distance = present[a].distance(present[b]) * metres
            pairs.append((a, b, distance, inverts[a], inverts[b], reason is None, reason))

        groups: dict[NodeId, set[NodeId]] = {}
        for node in set(parent) | set(parent.values()):
            groups.setdefault(root(node), set()).add(node)
        rank = {"maintenance_holes": 0}
        merged: dict[NodeId, NodeId] = {}
        for snapped in groups.values():
            layers = self.nodes.loc[sorted(snapped), "layer"]
            kept = min(
                snapped, key=lambda n: (rank.get(layers[n], 1 if pd.notna(layers[n]) else 2), n)
            )
            merged.update({n: kept for n in snapped if n != kept})

        report = pd.DataFrame(
            pairs,
            columns=["node_1", "node_2", "distance", "invert_1", "invert_2", "snapped", "reason"],
        )
        report.insert(6, "kept", report["node_1"].map(lambda n: merged.get(n, n)))
        report["kept"] = report["kept"].where(report["snapped"])
        if merged:
            self.full_data = self.full_data.assign(
                from_node=self.full_data["from_node"].replace(merged),
                to_node=self.full_data["to_node"].replace(merged),
            )
            self.aliases = {k: merged.get(v, v) for k, v in self.aliases.items()} | merged
            self._build()
        logger.info(
            "Snapped %d node IDs within %s m into %d nodes", len(merged), tolerance, len(groups)
        )
        return report

    def subnetwork(self, nodes: Iterable[NodeId]) -> Subnetwork:
        """The :class:`Subnetwork` made of ``nodes`` and every pipe between two of them.

        Raises
        ------
        ValueError
            If a node is not in the network.
        """
        return Subnetwork(self, nodes)

    def catchment(
        self,
        node: NodeId,
        max_hops: int | None = None,
        max_length: float | None = None,
        outlet: NodeId | Iterable[NodeId] | None = None,
        near: float | None = 100.0,
    ) -> Subnetwork:
        """The :class:`Subnetwork` of ``node``'s catchment: every node that drains to it.

        With ``max_hops`` or ``max_length`` the catchment is truncated, and with ``outlet`` the
        path from ``node`` down to the outlet(s) is added, so that the subnetwork discharges there.

        A node with two outlets, one draining to ``node`` and one elsewhere, is in the catchment,
        but its other outlet pipe is not, so the written model sends its flow to ``node``. Pipes
        leaving the catchment from nodes with area upstream, and force mains leaving it, are
        warned about (see :meth:`Subnetwork.boundary_pipes`). So is any other network whose end
        lies within ``near`` of ``node`` without draining to it, which may be a second inlet of
        the same plant (see :meth:`ends`).

        Parameters
        ----------
        node : str
            Node whose catchment is selected. Pass the outlet itself for its whole catchment.

        max_hops : int, optional
            Keep only nodes at most this many pipes upstream of ``node``.

        max_length : float, optional
            Keep only nodes whose shortest upstream path to ``node`` is at most this long, in m.

        outlet : node ID or list of node IDs, optional
            Outlet(s) downstream of ``node``. A shortest path (by length) is added to every listed
            outlet that ``node`` drains to.

        near : float, optional
            Distance (m) within which another network's end is reported (default 100 m); None
            skips the check.

        Warns
        -----
        DataWarning
            If pipes leave the catchment from nodes with area upstream or as force mains, or a
            network ends within ``near`` of ``node`` without draining to it.

        Raises
        ------
        ValueError
            If a node is not in the network, or ``node`` drains to none of the outlets.
        """
        node = self.resolve(node)
        self.check_nodes([node])
        upstream = self.graph.reverse(copy=False)

        if max_hops is None:
            selected = set(nx.descendants(upstream, node))
        else:
            selected = set(nx.single_source_shortest_path_length(upstream, node, cutoff=max_hops))
        if max_length is not None:
            within = nx.single_source_dijkstra_path_length(
                upstream, node, cutoff=max_length, weight=self._pipe_length
            )
            selected &= set(within)
        selected.add(node)

        if outlet is not None:
            if isinstance(outlet, str) or not isinstance(outlet, Iterable):
                outlets = [self.resolve(outlet)]
            else:
                outlets = [self.resolve(o) for o in outlet]
            self.check_nodes(outlets)
            reachable = [o for o in sorted(outlets) if nx.has_path(self.graph, node, o)]
            if not reachable:
                raise ValueError(f"Node {node} does not drain to any of the outlets {outlets}.")
            for o in reachable:
                selected.update(nx.dijkstra_path(self.graph, node, o, weight=self._pipe_length))

        sub = self.subnetwork(selected)
        self._check_catchment(sub, node, near)
        return sub

    def _check_catchment(self, sub, node, near):
        """Warn about the pipes leaving ``sub`` with flow to take, and the networks ending near
        ``node`` that do not drain to it."""
        leaving = sub.boundary_pipes(ends=False)
        leaving = leaving[leaving["direction"] == "out"]
        serious = leaving[leaving["kind"].isin(("split", "force_main"))]
        summits = int((leaving["kind"] == "summit").sum())
        if summits:
            logger.info(
                "%d pipes leave the catchment of %s from summit nodes (nothing upstream)",
                summits,
                node,
            )
        if len(serious):
            listed = ", ".join(
                f"{r.pipe_id} ({r.kind}, from {r.node}, {r.upstream_nodes} nodes upstream)"
                for r in serious.head(5).itertuples()
            )
            warnings.warn(
                f"{len(serious)} pipes leave the catchment of {node} from nodes with area "
                f"upstream ({int((serious['kind'] == 'force_main').sum())} force mains), so the "
                f"written model sends their flow down the nodes' other outlets: {listed}"
                f"{', ...' if len(serious) > 5 else ''}. See Subnetwork.boundary_pipes().",
                DataWarning,
                stacklevel=3,
            )
        if near is None:
            return
        for end, upstream, distance in self._ends_near(node, near, exclude=set(sub.graph)):
            warnings.warn(
                f"{end} ends a network of {upstream} nodes {distance:.0f} m from {node} without "
                "draining to it: if both reach the same plant, it is a second inlet to join.",
                DataWarning,
                stacklevel=3,
            )

    def _ends_near(self, node, distance, exclude=()):
        """(end, nodes draining to it, distance in m) of each node with no outlet pipe, not in
        ``exclude``, within ``distance`` of ``node``."""
        ends = [n for n in self.graph if self.graph.out_degree(n) == 0 and n not in exclude]
        pos = metric_positions(self, [*ends, node])
        if node not in pos.index:
            return []
        x0, y0 = pos.loc[node, "x"], pos.loc[node, "y"]
        found = []
        for end in ends:
            if end in pos.index:
                d = float(np.hypot(pos.loc[end, "x"] - x0, pos.loc[end, "y"] - y0))
                if d <= distance:
                    found.append((end, len(nx.ancestors(self.graph, end)), d))
        return sorted(found, key=lambda f: f[2])

    def ends(self) -> pd.DataFrame:
        """Every node with no outlet pipe: where the network's flow ends.

        The treatment plants' inlets are the few ends with many nodes draining to them, so
        listing them checks that the outlets being modelled account for the city; the many
        small ends are caps, outfalls, and ends of the mapped data.

        Returns
        -------
        pandas.DataFrame
            Indexed by node, with ``node_type`` and ``upstream_nodes`` (the nodes draining to
            it, itself excluded), the largest first.
        """
        ends = [n for n in self.graph if self.graph.out_degree(n) == 0]
        frame = pd.DataFrame(
            {
                "node_type": [self.graph.nodes[n].get("node_type") for n in ends],
                "upstream_nodes": [len(nx.ancestors(self.graph, n)) for n in ends],
            },
            index=pd.Index(ends, name="node"),
        )
        return frame.sort_values("upstream_nodes", ascending=False, kind="stable")

    def to_geodataframe(
        self, kind: str = "pipes", values: Values | None = None
    ) -> gpd.GeoDataFrame:
        """The pipes or nodes as a GeoDataFrame, in the network's CRS, for plotting or export.

        Parameters
        ----------
        kind : {'pipes', 'nodes'}
            'pipes': the rows of ``full_data`` (geometry and canonical columns), indexed by
            ``pipe_id``. 'nodes': the rows of :attr:`nodes` (point and canonical node
            columns), indexed by ``node_id``.

        values : pandas.Series, pandas.DataFrame or dict, optional
            Values indexed by pipe or node ID, left-joined as extra columns (a Series by its
            name, or 'value' if it has none). Pipes or nodes with no value get NaN.
        """
        return to_geodataframe(self.full_data, self.nodes, kind, values)

    def check_nodes(self, nodes: Iterable[NodeId]) -> None:
        """Raise ValueError listing the ``nodes`` that are not in the network."""
        missing = sorted(n for n in nodes if n not in self.graph)
        if missing:
            split = [
                f"{n} (split into {', '.join(map(str, self.splits[n]))})"
                for n in missing
                if n in self.splits
            ]
            raise ValueError(
                f"Nodes not in the network: {', '.join(map(str, missing))}"
                + (f"; reused IDs that were split: {'; '.join(split)}" if split else "")
            )

    _check_nodes = check_nodes  # the former name

    def ids(self) -> set[NodeId]:
        """Every node ID and pipe ID of the network, e.g. to check that a new name is free."""
        return set(self.graph.nodes) | set(self.full_data["pipe_id"])

    def lowest_inverts(self, nodes: Collection[NodeId]) -> dict[NodeId, float]:
        """{node: lowest invert of the network's pipes connected to it} for ``nodes`` (see
        :func:`swmmbuilder.gaps.lowest_inverts`)."""
        return lowest_inverts(self.full_data, nodes)

    def estimate_rims(
        self,
        dem: StrPath | Any | Iterable[StrPath | Any],
        window: float = 0.0,
        statistic: str = "median",
        offset: float | str = 0.0,
        nodes: Iterable[NodeId] | None = None,
        flag: float | None = None,
        replace_above: float | None = None,
        replace_below: float | None = None,
    ) -> RimEstimate:
        """Compare the recorded rims with a ground-elevation model (DEM), and estimate the
        missing ones from it.

        A maintenance-hole lid in a street is generally close to flush with the pavement, so the
        ground elevation at a node estimates its rim (see :mod:`swmmbuilder.elevation` for the
        caveats). The DEM is sampled at each node's point. Nodes with a usable rim (one that is
        recorded and above the invert, or with no invert to compare) are compared with it: see
        :attr:`RimEstimate.validation` and :meth:`RimEstimate.largest_differences`. Nodes with
        no usable rim are estimated as the DEM value plus ``offset``, unless that is at or below
        their invert.

        A usable recorded rim far from the DEM is probably wrong (a typo, a rim from another
        structure) or a lid not flush with the ground (raised in a field, a ravine, or an
        easement). Its ``deviation`` is recorded rim - DEM minus the median of that difference,
        so that a datum shift does not count. With ``flag``, the rims deviating by more than
        ``flag`` metres either way are flagged as implausible; with ``replace_above`` or
        ``replace_below``, those recorded more than that far above or below the ground are
        replaced by the DEM estimate. Lids are more often raised than sunk, so the threshold
        above can be the larger. Otherwise usable recorded rims are kept.

        The network is not changed. Pass :attr:`RimEstimate.rims` to the ``rims`` option of
        :meth:`Subnetwork.to_inp` (or :meth:`Subnetwork.data_gaps`, or
        :func:`swmmbuilder.gaps.fill_gaps`) to fill those nodes: each is listed in the data-gaps
        report (gap ``no_rim``, filled from ``estimate_rims``, with the DEM as its source), and
        its maximum depth is rim - invert. The flagged rims are listed too (gap
        ``implausible_rim``), filled where replaced.

        Parameters
        ----------
        dem : path, rasterio dataset, or list of them
            A single-band elevation raster (e.g. a GeoTIFF of a lidar bare-earth model) or a list
            of tiles, in any CRS, with elevations in m in the vertical datum of the network's
            data (see :func:`swmmbuilder.elevation.sample_dem`). Needs rasterio:
            ``pip install 'swmmbuilder[elevation]'``.

        window : float
            0 (default): the cell at the node. Otherwise the ``statistic`` of the cells within
            ``window`` metres of the node, to soften steep edges such as curbs and ditches.

        statistic : {'median', 'mean'}
            How the cells of a window are combined.

        offset : float or 'median'
            Added to the DEM value to estimate a rim, in m: a known datum shift or lid height,
            or ``'median'`` for the median recorded rim - DEM over the nodes compared.

        nodes : iterable, optional
            The nodes to sample (e.g. a catchment's); by default every node. With
            ``offset='median'``, the median is over these nodes.

        flag : float, optional
            Flag the recorded rims deviating from the DEM by more than this, in m (e.g. 1.0).

        replace_above, replace_below : float, optional
            Replace by the DEM estimate the recorded rims more than this far above, or below,
            the ground (``deviation`` above ``replace_above`` or below ``-replace_below``), in
            m. A replaced rim is flagged too. A DEM estimate at or below the invert does not
            replace a rim.

        Returns
        -------
        RimEstimate
            The DEM value, the comparison, and the estimate at each node.

        Raises
        ------
        ImportError
            If rasterio is not installed.

        ValueError
            If a node is not in the network, ``offset`` or a threshold is invalid, or
            ``offset='median'`` (or a threshold) with no node to compare (see also
            :func:`swmmbuilder.elevation.sample_dem`).
        """
        if isinstance(offset, str):
            if offset != "median":
                raise ValueError(f"offset must be a number or 'median'; got {offset!r}.")
        elif not np.isfinite(offset):
            raise ValueError(f"offset must be a finite number or 'median'; got {offset}.")
        thresholds = {"flag": flag, "replace_above": replace_above, "replace_below": replace_below}
        for name, value in thresholds.items():
            if value is not None and not (_is_real(value) and value >= 0):
                raise ValueError(f"{name} must be a number of metres >= 0; got {value!r}.")
        ids = list(self.graph) if nodes is None else sorted({self.resolve(n) for n in nodes})
        self.check_nodes(ids)
        table = self.nodes.loc[ids]
        recorded = table["rim_elevation"].astype(float)
        invert = table["invert"].astype(float)
        dem_values = sample_dem(dem, table.geometry, window, statistic)

        usable = recorded.notna() & ~(recorded <= invert)
        has_dem = dem_values.notna()
        compared = usable & has_dem
        difference = (recorded - dem_values).where(compared)
        if offset == "median":
            if not compared.any():
                raise ValueError(
                    "offset='median': no node has both a usable recorded rim and a DEM value."
                )
            offset = float(difference.median())
        offset = float(offset)
        estimate = (dem_values + offset).round(4)
        below = has_dem & ~usable & (estimate <= invert)
        if any(v is not None for v in thresholds.values()) and not compared.any():
            raise ValueError(
                "flag/replace: no node has both a usable recorded rim and a DEM value."
            )
        deviation = difference - difference.median()
        never = pd.Series(False, index=table.index)
        replace = never.copy()
        if replace_above is not None:
            replace |= deviation > replace_above
        if replace_below is not None:
            replace |= deviation < -replace_below
        replace &= estimate > invert.fillna(-np.inf)
        flagged = (deviation.abs() > flag if flag is not None else never) | replace
        status = pd.Series("no_dem_value", index=table.index, dtype=object)
        status[compared] = "compared"
        status[flagged] = "implausible"
        status[replace] = "replaced"
        status[has_dem & ~usable] = "estimated"
        status[below] = "dem_not_above_invert"
        samples = pd.DataFrame(
            {
                "recorded_rim": recorded,
                "invert": invert,
                "dem": dem_values,
                "difference": difference,
                "deviation": deviation,
                "estimate": estimate,
                "status": status,
            }
        )
        source = f"DEM {dem_name(dem)}"
        if window:
            source += f", {statistic} within {window:g} m"
        if offset:
            source += f", offset {offset:+.3f} m"
        result = RimEstimate(samples, offset, source, flag, replace_above, replace_below)

        logger.info("%s", result.summary())
        if below.any():
            warnings.warn(
                f"{int(below.sum())} rims from the DEM are at or below the invert and are not "
                "used (status 'dem_not_above_invert' in RimEstimate.samples)",
                DataWarning,
                stacklevel=2,
            )
        if flagged.any():
            warnings.warn(
                f"{int(flagged.sum())} recorded rims are implausible against the DEM "
                f"({int(replace.sum())} replaced); see RimEstimate.implausible",
                DataWarning,
                stacklevel=2,
            )
        missing = int((~usable & ~has_dem).sum())
        if missing:
            logger.warning("%d nodes with no usable rim have no DEM value", missing)
        return result

    @staticmethod
    def _pipe_length(u, v, pipes):
        """Dijkstra weight: the shortest of the (possibly parallel) pipes from u to v."""
        return min(attrs["length"] for attrs in pipes.values())


def pipe_end_points(data: pd.DataFrame, nodes: Iterable[NodeId]) -> dict[NodeId, Point]:
    """{node: shapely Point} for the nodes in ``nodes``, from the pipe geometry.

    A node is placed at the upstream end of a pipe leaving it or, if no pipe leaves it, at the
    downstream end of a pipe entering it.
    """
    return {node: Point(xy) for node, xy in node_positions(data, nodes).items()}
