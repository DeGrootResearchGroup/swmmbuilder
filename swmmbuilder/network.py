"""A city's sewer network: its pipes and node layers, loaded through a
:class:`~swmmbuilder.Schema`, and the flow graph between their nodes."""

import logging
import re
from collections.abc import Iterable
from pathlib import Path

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely import Point

from .gaps import lowest_inverts
from .geometry import is_disjoint, node_positions, to_geodataframe
from .schema import LONDON_ON, NODE_COLUMNS, as_schema
from .subnetwork import Subnetwork

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


def strip_leading_zeros(node_id):
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
        Applied to every node ID of the pipes and node layers, and to the node IDs passed to
        :meth:`catchment`, :meth:`subnetwork` and :meth:`Subnetwork.from_node_list`, e.g.
        :func:`strip_leading_zeros`.

    snap_tolerance : float, optional
        If set, :meth:`snap_nodes` is called with it (and ``max_invert_difference``) on loading,
        and its report is kept as ``snap_report``.

    max_invert_difference : float, optional
        See :meth:`snap_nodes`.

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

    snap_report : pandas.DataFrame or None
        The report of :meth:`snap_nodes` if ``snap_tolerance`` was given.
    """

    def __init__(
        self,
        source,
        flow_types=DEFAULT_FLOW_TYPES,
        keep_abandoned=False,
        schema=LONDON_ON,
        status=DEFAULT_STATUSES,
        include_overflows=True,
        maintenance_holes=None,
        pump_stations=None,
        outfalls=None,
        other_nodes=None,
        normalize_ids=None,
        snap_tolerance=None,
        max_invert_difference=None,
    ):
        self.schema = schema = as_schema(schema)
        self.normalize_ids = normalize_ids
        self.aliases = {}
        self.snap_report = None
        self.full_data = schema.normalize(_read(source))
        if normalize_ids is not None:
            for column in ("from_node", "to_node"):
                self.full_data[column] = self.full_data[column].map(normalize_ids)
        self.removed = {"flow_type": {}, "status": {}, "structure_type": {}, "invalid": {}}
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

        self._layers = self._read_layers(
            dict(
                maintenance_holes=maintenance_holes,
                pump_stations=pump_stations,
                outfalls=outfalls,
                other_nodes=other_nodes,
            )
        )
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

    def _read_layers(self, layers):
        """The canonical node layers, in order of precedence, each with a ``layer`` column."""
        frames = []
        for name, layer in layers.items():
            if layer is None:
                continue
            nodes = self.schema.normalize_nodes(_read(layer), NODE_LAYERS[name])
            if self.normalize_ids is not None:
                nodes["node_id"] = nodes["node_id"].map(self.normalize_ids)
            if self.full_data.crs is not None and nodes.crs != self.full_data.crs:
                nodes = nodes.to_crs(self.full_data.crs)
            frames.append(nodes.assign(layer=name))
            logger.info("Read %d nodes from the %s layer", len(nodes), name)
        if not frames:
            return None
        return pd.concat(frames, ignore_index=True)

    def _build(self):
        """(Re)build ``graph`` and ``nodes`` from ``full_data`` and the node layers."""
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
        missing = geometry.isna()
        if missing.any():
            ends = pipe_end_points(data, set(nodes.index[missing]))
            geometry = geometry.where(~missing, pd.Series(ends).reindex(nodes.index))
        geometry = gpd.GeoSeries(geometry.where(geometry.notna(), None), crs=data.crs)
        return gpd.GeoDataFrame(nodes, geometry=geometry)

    def resolve(self, node):
        """The node ID that ``node`` refers to: normalized and with merged IDs replaced."""
        if self.normalize_ids is not None:
            node = self.normalize_ids(node)
        return self.aliases.get(node, node)

    def snap_nodes(self, tolerance=0.5, max_invert_difference=None):
        """Merge node IDs whose points are within ``tolerance`` of each other.

        The same structure sometimes has two IDs at one point, which disconnects the network
        there. Each group of nodes within ``tolerance`` (in m; geographic coordinates are
        projected to the local UTM zone) is merged into one ID, preferring an ID in the
        maintenance-hole layer, then one in another node layer, then the smallest. Pairs whose
        inverts differ by more than ``max_invert_difference`` (pipes that cross rather than join)
        and pairs already joined by a pipe are not merged. Since merging is transitive, both
        rules hold for each merged group as a whole: a pair is not merged if that would put two
        nodes joined by a pipe, or two inverts more than ``max_invert_difference`` apart, into
        one group. Merged IDs are kept in ``aliases``, so either ID can be used.

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
        pairs, parent, members = [], {}, {}

        def root(node):
            while parent.get(node, node) != node:
                node = parent[node]
            return node

        def group(node):
            return members.get(root(node), {node})

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
            elif max_invert_difference is not None and not difference <= max_invert_difference:
                reason = f"inverts differ by more than {max_invert_difference} m"
            elif root(a) != root(b):
                merged = group(a) | group(b)
                spread = inverts[sorted(merged)]
                if any(neighbours(n) & merged for n in group(a)):
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
                    members[ra] = merged
                    members.pop(rb, None)
            distance = present[a].distance(present[b]) * metres
            pairs.append((a, b, distance, inverts[a], inverts[b], reason is None, reason))

        groups = {}
        for node in set(parent) | set(parent.values()):
            groups.setdefault(root(node), set()).add(node)
        rank = {"maintenance_holes": 0}
        merged = {}
        for group in groups.values():
            layers = self.nodes.loc[sorted(group), "layer"]
            kept = min(
                group, key=lambda n: (rank.get(layers[n], 1 if pd.notna(layers[n]) else 2), n)
            )
            merged.update({n: kept for n in group if n != kept})

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

    def subnetwork(self, nodes):
        """The :class:`Subnetwork` made of ``nodes`` and every pipe between two of them.

        Raises
        ------
        ValueError
            If a node is not in the network.
        """
        return Subnetwork(self, nodes)

    def catchment(self, node, max_hops=None, max_length=None, outlet=None):
        """The :class:`Subnetwork` of ``node``'s catchment: every node that drains to it.

        With ``max_hops`` or ``max_length`` the catchment is truncated, and with ``outlet`` the
        path from ``node`` down to the outlet(s) is added, so that the subnetwork discharges there.

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

        Raises
        ------
        ValueError
            If a node is not in the network, or ``node`` drains to none of the outlets.
        """
        node = self.resolve(node)
        self._check_nodes([node])
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
            single = isinstance(outlet, str) or not isinstance(outlet, Iterable)
            outlets = [outlet] if single else list(outlet)
            outlets = [self.resolve(o) for o in outlets]
            self._check_nodes(outlets)
            reachable = [o for o in sorted(outlets) if nx.has_path(self.graph, node, o)]
            if not reachable:
                raise ValueError(f"Node {node} does not drain to any of the outlets {outlets}.")
            for o in reachable:
                selected.update(nx.dijkstra_path(self.graph, node, o, weight=self._pipe_length))

        return self.subnetwork(selected)

    def to_geodataframe(self, kind="pipes", values=None):
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

    def _check_nodes(self, nodes):
        missing = sorted(n for n in nodes if n not in self.graph)
        if missing:
            raise ValueError(f"Nodes not in the network: {', '.join(map(str, missing))}")

    @staticmethod
    def _pipe_length(u, v, pipes):
        """Dijkstra weight: the shortest of the (possibly parallel) pipes from u to v."""
        return min(attrs["length"] for attrs in pipes.values())


def pipe_end_points(data, nodes):
    """{node: shapely Point} for the nodes in ``nodes``, from the pipe geometry.

    A node is placed at the upstream end of a pipe leaving it or, if no pipe leaves it, at the
    downstream end of a pipe entering it.
    """
    return {node: Point(xy) for node, xy in node_positions(data, nodes).items()}
