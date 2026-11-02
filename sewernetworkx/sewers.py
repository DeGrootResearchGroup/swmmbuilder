"""Module for manipulating GeoJSON sewer data.

The module takes full GeoJSON data sets of a city sewer network,
generates a subnetwork given an initial starting node and number
of levels, and generates a .inp file to be used by pyswmm.
"""

import os
import re
import warnings
from dataclasses import dataclass, field, replace
from importlib import resources
from pathlib import Path

import geopandas as gpd
import gmplot
import networkx as nx
import pandas as pd

from .schema import CANONICAL_COLUMNS, LONDON_ON, as_schema

# SWMM's SI flow units. With any of these, lengths, elevations and section sizes are in metres,
# matching the GIS data; with a US unit (CFS, GPM, MGD) SWMM would read them as feet.
SI_FLOW_UNITS = ("CMS", "LPS", "MLD")

# Canonical flow types that Network keeps by default: sanitary sewers only.
DEFAULT_FLOW_TYPES = frozenset({"sanitary"})

# Placeholders that write_pyswmm fills in the template.
TEMPLATE_PLACEHOLDERS = (
    "nodes",
    "outfalls",
    "storage",
    "conduits",
    "pumps",
    "xsections",
    "curves",
    "coordinates",
    "vertices",
)

# SWMM pump curve types, and what the x values of each are (y is flow, in the FLOW_UNITS).
PUMP_CURVE_TYPES = {
    "PUMP1": "wet well volume",
    "PUMP2": "inlet node depth",
    "PUMP3": "head difference",
    "PUMP4": "inlet node depth",
    "PUMP5": "head difference",
}


def default_template_path():
    """Path-like (Traversable) to the SI template shipped with the package."""
    return resources.files("sewernetworkx") / "templates" / "template.inp"


def check_template(content):
    """Check that a .inp template is one the writer can fill.

    Parameters
    ----------
    content : str
        Text of the .inp template.

    Raises
    ----------
    ValueError
        If FLOW_UNITS is missing or not an SI unit (the GIS data are metric), or if the template
        contains a placeholder the writer does not fill (see ``TEMPLATE_PLACEHOLDERS``).
    """
    match = re.search(r"^\s*FLOW_UNITS\s+(\S+)", content, re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise ValueError(
            "The template has no FLOW_UNITS option (SWMM would default to CFS, US units). "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )
    flow_units = match.group(1).upper()
    if flow_units not in SI_FLOW_UNITS:
        raise ValueError(
            f"The template sets FLOW_UNITS {flow_units}, a US unit, so SWMM would read the metric "
            "GIS lengths, elevations and sizes as feet. "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )

    unknown = sorted(set(re.findall(r"\{(\w+)\}", content)) - set(TEMPLATE_PLACEHOLDERS))
    if unknown:
        raise ValueError(
            f"The template contains placeholders that are not filled: {unknown}. "
            f"Supported placeholders: {list(TEMPLATE_PLACEHOLDERS)}."
        )


def line_coords(geometry):
    """(x, y) vertices of a pipe geometry, in order from upstream to downstream.

    The parts of a MultiLineString are joined in order. Any z values are dropped.
    """
    parts = geometry.geoms if geometry.geom_type == "MultiLineString" else [geometry]
    return [(c[0], c[1]) for part in parts for c in part.coords]


def xsection(pipe):
    """SWMM cross-section of a pipe, from its canonical shape and sizes (in m).

    CIRCULAR uses the diameter and EGG the full height; the other shapes use the full height
    and maximum width.

    Parameters
    ----------
    pipe : mapping
        A row of the canonical data (see :mod:`sewernetworkx.schema`).

    Returns
    ----------
    tuple or None
        (shape, geom1, geom2), or None if the shape is unknown or missing, or a size it needs is
        missing or not positive.
    """
    shape = pipe.get("shape")
    if pd.isna(shape):
        return None
    if shape == "CIRCULAR":
        sizes = [pipe.get("diameter")]
    elif shape == "EGG":
        sizes = [pipe.get("height")]
    else:
        sizes = [pipe.get("height"), pipe.get("width")]
    if any(pd.isna(s) or s <= 0 for s in sizes):
        return None
    return shape, *sizes, *[0.0] * (2 - len(sizes))


@dataclass
class Pump:
    """A pump that replaces a pipe in the .inp file, as passed to :meth:`SubNetwork.write_pyswmm`.

    Attributes
    ----------
    curve : list of (float, float) or str
        Pump curve points (x, y), y being the flow in the template's FLOW_UNITS (L/s with the
        shipped template) and x as given by ``curve_type``; or ``'*'`` for an ideal pump, which
        pumps whatever flows into its inlet node.

    curve_type : str
        PUMP1 (x: wet well volume, m3; the inlet must be a wet well), PUMP2 (x: inlet depth, m,
        stepwise), PUMP3 (x: head difference, m), PUMP4 (x: inlet depth, m, continuous) or PUMP5
        (variable-speed PUMP3). Ignored for an ideal pump.

    status : str
        Initial status, 'ON' or 'OFF'.

    startup, shutoff : float
        Inlet node depths (m) at which the pump turns on and off; 0 for none.
    """

    curve: object
    curve_type: str = "PUMP3"
    status: str = "ON"
    startup: float = 0.0
    shutoff: float = 0.0

    @classmethod
    def coerce(cls, pump_id, value):
        """A Pump from a Pump or a dict of its fields, checked; raise ValueError if it is invalid."""
        pump = replace(value) if isinstance(value, cls) else cls(**value)
        problems = []
        if isinstance(pump.curve, str):
            if pump.curve != "*":
                problems.append(f"curve must be a list of (x, y) points or '*', not {pump.curve!r}")
        else:
            points = [tuple(pt) for pt in pump.curve]
            if not points or any(len(pt) != 2 for pt in points):
                problems.append("curve must be a non-empty list of (x, y) points")
            else:
                pump.curve = [(float(x), float(y)) for x, y in points]
            if pump.curve_type.upper() not in PUMP_CURVE_TYPES:
                problems.append(
                    f"curve_type {pump.curve_type!r} is not one of {list(PUMP_CURVE_TYPES)}"
                )
            pump.curve_type = pump.curve_type.upper()
        pump.status = pump.status.upper()
        if pump.status not in ("ON", "OFF"):
            problems.append(f"status must be 'ON' or 'OFF', not {pump.status!r}")
        if problems:
            raise ValueError(f"Pump {pump_id}: {'; '.join(problems)}.")
        return pump


@dataclass
class WetWell:
    """A node written as a [STORAGE] node of constant plan area, e.g. a pumping station wet well.

    Attributes
    ----------
    max_depth : float
        Depth (m) from the node invert to the top of the wet well.

    area : float
        Plan area (m2).

    init_depth : float
        Initial water depth (m).
    """

    max_depth: float
    area: float
    init_depth: float = 0.0

    @classmethod
    def coerce(cls, node, value):
        """A WetWell from a WetWell or a dict of its fields; raise ValueError if it is invalid."""
        well = replace(value) if isinstance(value, cls) else cls(**value)
        if not (well.max_depth > 0 and well.area > 0 and well.init_depth >= 0):
            raise ValueError(
                f"Wet well {node}: max_depth and area must be positive and init_depth not negative."
            )
        return well


def _invert(value):
    """An invert elevation, or None if it is missing."""
    return None if pd.isna(value) else float(value)


def _pipe_records(data):
    """Yield (pipe_id, from_node, to_node, length) for every pipe with two node IDs.

    Pipes with a missing node ID are skipped.
    """
    for oid, up, down, length in zip(
        data["pipe_id"], data["from_node"], data["to_node"], data["length"], strict=True
    ):
        if pd.isna(up) or pd.isna(down):
            continue
        yield oid, up, down, length


@dataclass
class TopologyReport:
    """Topology problems of a subnetwork, as returned by :meth:`SubNetwork.validate`.

    Attributes
    ----------
    splits : dict
        Node -> sorted pipe IDs of its outlet pipes, for every node with more than one outlet
        pipe. Parallel pipes count separately, as they do in SWMM.

    cycles : list of list
        Each directed cycle as a list of nodes in flow order, starting from its smallest node.

    dead_ends : list
        Nodes with no outlet pipe that are not outfalls.

    parallel_pipes : dict
        (upstream, downstream) -> sorted pipe IDs, for every node pair joined by more than one pipe.
    """

    splits: dict = field(default_factory=dict)
    cycles: list = field(default_factory=list)
    dead_ends: list = field(default_factory=list)
    parallel_pipes: dict = field(default_factory=dict)

    @property
    def is_dendritic(self):
        """True when the network has no splits and no cycles, as kinematic-wave routing requires."""
        return not self.splits and not self.cycles


class Network:
    """Class used to represent an entire cities sewer network

    Parameters
    ----------
    file_path : str, Path
        Name of GeoJSON file

    Attributes
    ----------
    full_data : geopandas.GeoDataFrame
        Pipe data read from the GeoJSON file, with the canonical columns of
        :mod:`sewernetworkx.schema`.

    schema : Schema
        The schema the data were loaded with.

    graph : networkx.MultiDiGraph
        One edge per pipe, in the direction of flow (from_node -> to_node), keyed by pipe_id so
        that parallel pipes are kept. Each edge has the attributes ``pipe`` (pipe_id) and
        ``length``. Nodes and edges are inserted in sorted order, so traversals do not depend
        on the order of the file.

    remove_storm : bool, optional
        Deprecated. ``False`` is the same as ``flow_types=None``; ``True`` keeps the default.

    flow_types : iterable of str or None
        Flow types to keep (default ``{'sanitary'}``). Source codes are translated by the schema,
        so for London ``{'SAN'}`` is the same. Every other pipe, including those with no flow type,
        is dropped and the counts are printed. ``None`` keeps every pipe.

    keep_abandoned : bool
        Keep pipes flagged abandoned (default False: they are dropped). Pipes with a node ID the
        schema treats as no node (for London, 'Abandoned') are dropped from the graph either way.

    schema : Schema, str or Path
        How the GIS columns map to the canonical columns: a :class:`Schema` (default
        :data:`LONDON_ON`; see also :data:`TORONTO_ON`) or the path to a schema TOML file.

    Attributes
    ----------
    removed : dict
        What was dropped on loading: ``{'flow_type': {flow_type: count}, 'abandoned': count}``.
        Pipes with no flow type are counted under ``'None'``.
    """

    def __init__(
        self,
        file_path,
        remove_storm=None,
        flow_types=DEFAULT_FLOW_TYPES,
        keep_abandoned=False,
        schema=LONDON_ON,
    ):
        """Class initialization"""

        self.file_path = file_path
        self.schema = schema = as_schema(schema)
        self.full_data = schema.normalize(gpd.read_file(self.file_path))

        if remove_storm is not None:
            warnings.warn(
                "remove_storm is deprecated; use flow_types (flow_types=None keeps every pipe).",
                DeprecationWarning,
                stacklevel=2,
            )
            if not remove_storm:
                flow_types = None

        self.removed = {"flow_type": {}, "abandoned": 0}

        if flow_types is not None and "flow_type" in schema.columns:
            flow_types = {schema.flow_type(f) for f in flow_types}
            keep = self.full_data["flow_type"].isin(flow_types)
            dropped = self.full_data.loc[~keep, "flow_type"].fillna("None").astype(str)
            self.removed["flow_type"] = dropped.value_counts().sort_index().to_dict()
            self.full_data = self.full_data[keep]
            for flow_type, count in self.removed["flow_type"].items():
                print(
                    f"Removed {count} pipes with flow type {flow_type} "
                    f"(keeping {sorted(flow_types)})"
                )

        if not keep_abandoned and schema.abandoned is not None:
            abandoned = self.full_data["abandoned"]
            self.removed["abandoned"] = int(abandoned.sum())
            self.full_data = self.full_data[~abandoned]
            print(f"Removed {self.removed['abandoned']} pipes flagged abandoned")

        self.graph = self._build_graph(self.full_data)

    @staticmethod
    def _build_graph(data):
        """Build the flow-direction MultiDiGraph of the pipes in ``data``."""
        records = sorted(_pipe_records(data), key=lambda r: (r[1], r[2], r[0]))
        G = nx.MultiDiGraph()
        G.add_nodes_from(sorted({n for _, up, down, _ in records for n in (up, down)}))
        for oid, up, down, length in records:
            G.add_edge(up, down, key=oid, pipe=oid, length=length)
        return G

    """
    # Needs to be updated
    def filter_data_gisid(self, network_data, ppty):
        Filter Network class data.

                Parameters
                ----------
                network_data : class '__main__.Network'
                    Sewer network class

                ppty : str
                    Property to filter data by

                Returns
                -------
                str
                    The GIS_ID of the filtered data

        if not isinstance(ppty, str):
            raise TypeError("Input must be a string.")

        return network_data[ppty].iloc[0]
    """

    def get_upstream_invert(self, upstream_id):
        """Find the upstream invert of a pipe leaving a node (or the downstream one of a pipe
        entering it, if none leaves it).

        Parameters
        ----------
        upstream_id : str
            The node to search for.

        Returns
        ----------
        float or None
            The corresponding invert.
        """
        if not isinstance(upstream_id, str):
            raise TypeError("Input must be a string")

        inverts = self._first_inverts()
        if upstream_id not in inverts:
            raise TypeError("Cant find either upstreaminvert of downstreaminvert")
        return inverts[upstream_id]

    def _first_inverts(self):
        """{node: us_invert of the first pipe leaving it or, if no pipe leaves it, ds_invert of
        the first pipe entering it}, over ``full_data``.

        Built on first use and rebuilt if ``full_data`` is replaced.
        """
        cached = getattr(self, "_inverts", None)
        if cached is not None and cached[0] is self.full_data:
            return cached[1]

        data = self.full_data
        upstream, downstream = {}, {}
        for up, down, up_inv, down_inv in zip(
            data["from_node"],
            data["to_node"],
            data["us_invert"],
            data["ds_invert"],
            strict=True,
        ):
            upstream.setdefault(up, up_inv)
            downstream.setdefault(down, down_inv)

        inverts = {**downstream, **upstream}
        self._inverts = (data, inverts)
        return inverts


class SubNetwork(Network):
    """Class used to represent an entire cities sewer network

    Parameters
    ----------
    subnetwork_data : dict
        dict containing all required information to generate subnetwork.
        Required entries:

        network_data : geoDataFrame
        starting_node : str
        outfall_node : list, str
        outfall_node_type : list, str
        depth : int
        scale : int

    file_path : str, Path
        Path to GeoJSON file. If given, the pipes are read from this file instead of
        ``network_data``.

    schema : Schema, str or Path
        Schema of the file at ``file_path`` (default :data:`LONDON_ON`), or the path to a schema
        TOML file.

    Notes
    -----
    Select the subnetwork with one of :meth:`catchment`, :meth:`from_nodes` or
    :meth:`from_node_list`; each sets ``subgraph`` deterministically. ``subgraph`` is a
    ``networkx.DiGraph`` whose edges are stored as (downstream, upstream), the orientation the
    .inp writer expects, with the attribute ``pipes``: the sorted pipe IDs of the (possibly
    parallel) pipes between the two nodes, by pipe_id.
    """

    def __init__(self, subnetwork_data, file_path=None, schema=LONDON_ON):
        if file_path:
            super().__init__(file_path, schema=schema)
            source = self
        else:
            source = subnetwork_data.get("network_data")

        required_keys = [
            "network_data",
            "starting_node",
            "outfall_node",
            "outfall_node_type",
            "depth",
            "scale",
        ]
        # Ensure all keys are being specified
        missing_keys = [key for key in required_keys if key not in subnetwork_data]
        if missing_keys:
            raise ValueError(f"Missing the following keys: {', '.join(missing_keys)}")

        self.network_data = {}
        self.start = subnetwork_data["starting_node"]
        self.depth = subnetwork_data["depth"]
        self.outfall_node = subnetwork_data["outfall_node"]
        self.outfall_node_type = subnetwork_data["outfall_node_type"]
        self.scale = subnetwork_data["scale"]
        self.subgraph = nx.DiGraph()
        self.node_colour = []
        self.coords = {}

        if len(self.outfall_node) != len(self.outfall_node_type):
            raise ValueError(
                f"Length of {self.outfall_node} and {self.outfall_node_type} are not the same."
            )

        self.schema = source.schema
        self.full_data = source.full_data
        self.graph = source.graph
        self.reduced_keys_data = source.full_data[[*CANONICAL_COLUMNS, "geometry"]].copy()

        print("Creating Pipes")
        for oid, up, down, _ in _pipe_records(self.reduced_keys_data):
            self.network_data[oid] = (down, up)

    def catchment(self, node, max_hops=None, max_length=None, outlet=None):
        """Select the catchment of ``node``: every node that drains to it.

        With ``max_hops`` or ``max_length`` the catchment is truncated, and the path from ``node``
        down to the outlet is added so that the subnetwork still discharges there.

        Parameters
        ----------
        node : str
            Node whose catchment is selected. Pass the outlet itself for its whole catchment.

        max_hops : int, optional
            Keep only nodes at most this many pipes upstream of ``node``.

        max_length : float, optional
            Keep only nodes whose shortest upstream path to ``node`` is at most this long, in the
            units of the ``Length`` column.

        outlet : str or list of str, optional
            Outlet(s) to connect ``node`` to. Defaults to ``outfall_node``. A shortest path (by
            length) is added to every listed outlet that ``node`` drains to.

        Returns
        ----------
        networkx.DiGraph
            The new ``subgraph``.
        """
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

        if max_hops is not None or max_length is not None:
            outlets = self.outfall_node if outlet is None else outlet
            outlets = [outlets] if isinstance(outlets, str) else list(outlets)
            self._check_nodes(outlets)
            reachable = [o for o in sorted(outlets) if nx.has_path(self.graph, node, o)]
            if not reachable:
                raise ValueError(f"Node {node} does not drain to any of the outlets {outlets}.")
            for o in reachable:
                selected.update(nx.dijkstra_path(self.graph, node, o, weight=self._pipe_length))

        return self.from_nodes(selected)

    def from_nodes(self, nodes):
        """Select the subnetwork made of ``nodes`` and every pipe between two of them.

        Parameters
        ----------
        nodes : iterable of str
            Node IDs.

        Returns
        ----------
        networkx.DiGraph
            The new ``subgraph``.
        """
        nodes = sorted(set(nodes))
        self._check_nodes(nodes)
        induced = self.graph.subgraph(nodes)

        subgraph = nx.DiGraph()
        subgraph.add_nodes_from(nodes)
        pipes = {}
        for up, down, oid in induced.edges(keys=True):
            pipes.setdefault((down, up), []).append(oid)
        for edge in sorted(pipes):
            subgraph.add_edge(*edge, pipes=sorted(pipes[edge]))

        self.subgraph = subgraph
        return subgraph

    def to_node_list(self, path):
        """Write the nodes of ``subgraph`` to a text file, one per line, in sorted order.

        :meth:`from_node_list` rebuilds the same subnetwork from this file.

        Parameters
        ----------
        path : str, Path
            File to write.
        """
        lines = ["# SewerNetworkx node list: one node ID per line"]
        lines += sorted(self.subgraph.nodes)
        Path(path).write_text("\n".join(lines) + "\n")

    def from_node_list(self, path):
        """Select the subnetwork whose nodes are listed in a file written by :meth:`to_node_list`.

        Blank lines and lines starting with '#' are ignored.

        Parameters
        ----------
        path : str, Path
            File to read.

        Returns
        ----------
        networkx.DiGraph
            The new ``subgraph``.
        """
        lines = (line.strip() for line in Path(path).read_text().splitlines())
        return self.from_nodes(line for line in lines if line and not line.startswith("#"))

    def _check_nodes(self, nodes):
        missing = sorted(n for n in nodes if n not in self.graph)
        if missing:
            raise ValueError(f"Nodes not in the network: {', '.join(map(str, missing))}")

    @staticmethod
    def _pipe_length(u, v, pipes):
        """Dijkstra weight: the shortest of the (possibly parallel) pipes from u to v."""
        return min(attrs["length"] for attrs in pipes.values())

    def _outfalls(self):
        return (
            [self.outfall_node] if isinstance(self.outfall_node, str) else list(self.outfall_node)
        )

    def _outfall_types(self):
        types = self.outfall_node_type
        return [types] if isinstance(types, str) else list(types)

    def _edge_pipes(self, up, down):
        """pipe IDs of the pipes from ``up`` to ``down`` in ``subgraph``.

        Uses the ``pipes`` edge attribute, or the full network's pipes when a subgraph was built
        without it.
        """
        pipes = self.subgraph.edges[down, up].get("pipes")
        if pipes is None:
            pipes = sorted(self.graph[up][down]) if self.graph.has_edge(up, down) else []
        return pipes

    def _outlet_pipes(self, node):
        """{downstream node: pipe IDs} of the pipes leaving ``node`` in ``subgraph``."""
        return {down: self._edge_pipes(node, down) for down in self.subgraph.predecessors(node)}

    def validate(self):
        """Check the topology of ``subgraph`` before it is written to an .inp file.

        Returns
        ----------
        TopologyReport
            Flow splits, cycles, dead ends and parallel pipes. Splits and cycles are only allowed
            with dynamic-wave routing; dead ends are never valid.
        """
        outfalls = set(self._outfalls())
        report = TopologyReport()

        for node in sorted(self.subgraph.nodes):
            outlets = self._outlet_pipes(node)
            pipes = sorted(p for down in outlets for p in outlets[down])
            # An edge with no matching pipe still counts as one outlet link.
            n_links = sum(len(p) or 1 for p in outlets.values())
            if n_links > 1:
                report.splits[node] = pipes
            if n_links == 0 and node not in outfalls:
                report.dead_ends.append(node)
            for down in sorted(outlets):
                if len(outlets[down]) > 1:
                    report.parallel_pipes[(node, down)] = outlets[down]

        # subgraph edges are (downstream, upstream), so reverse cycles into flow order.
        for cycle in nx.simple_cycles(self.subgraph):
            cycle = cycle[::-1]
            i = cycle.index(min(cycle))
            report.cycles.append(cycle[i:] + cycle[:i])
        report.cycles.sort()

        return report

    def create_sub_network(self, ignore_successors=False, ignore_predecessors=False):
        """Create subnetwork, generate node colouring, and generate node locations.

        .. deprecated::
            The selected nodes are not a well-defined function of ``depth``. Use
            :meth:`catchment`, :meth:`from_nodes` or :meth:`from_node_list` instead.

        Returns
        ----------
        ignore_successors : bool
            ignore successor branch from starting node

        ignore_predecessors : bool
            ignore predecessors from starting node
        """
        warnings.warn(
            "create_sub_network is deprecated; use catchment, from_nodes or from_node_list.",
            DeprecationWarning,
            stacklevel=2,
        )
        print("Creating subnetwork for plotting")

        G = self._initialize_digraph()

        # Generate subnetwork starting with successors of starting node
        print("Checking successor nodes")
        successors_subgraph = self._explore_nodes(G, direction="successors")

        # Generate subnetwork with predecessors of starting node
        print("Checking predecessor nodes")
        predecessors_subgraph = self._explore_nodes(G, direction="predecessors")

        if ignore_successors:
            connected_subgraph = predecessors_subgraph
        elif ignore_predecessors:
            connected_subgraph = successors_subgraph
        else:
            connected_subgraph = successors_subgraph | predecessors_subgraph

        # Create a subgraph from the connected nodes
        self.subgraph = G.subgraph(connected_subgraph)

    def _initialize_digraph(self):
        # Initialize Networkx graph
        G = nx.DiGraph()

        # Add edges from network_data
        G.add_edges_from(self.network_data.values())

        return G

    def add_colouring(self):
        """Apply colouring for subnetwork nodes.

        Note: This method is intended for internal use and should not be called directly
        by external code.

        Returns
        ----------
        None
        """

        for node in self.subgraph.nodes:
            if node == self.start or node == self.outfall_node:
                self.node_colour.append("green")
            elif len(list(self.subgraph.successors(node))) > 1:
                self.node_colour.append("red")
            elif len(list(self.subgraph.successors(node))) == 0:
                self.node_colour.append("orange")
            else:
                self.node_colour.append("blue")

    def add_coordinates(self):
        """Apply node coordinates for subnetwork nodes.

        Note: This method is intended for internal use and should not be called directly
        by external code.

        Returns
        ----------
        None
        """

        for node, (x, y) in self.node_positions(self.subgraph.nodes()).items():
            self.coords[node] = (x * self.scale, y * self.scale)

        missing_nodes = set(self.subgraph.nodes()) - set(self.coords.keys())
        if missing_nodes:
            print(f"Error: Positions for nodes {missing_nodes} could not be found.")

    def node_positions(self, nodes):
        """Unscaled (x, y) position of each node, from the GIS pipe geometry.

        A node is placed at the upstream end of a pipe leaving it or, if no pipe leaves it, at the
        downstream end of a pipe entering it.

        Parameters
        ----------
        nodes : iterable
            Nodes to place.

        Returns
        ----------
        dict
            {node: (x, y)} for the nodes that were found; nodes with no pipe are left out.
        """
        nodes = set(nodes)
        upstream_ends, downstream_ends = {}, {}
        data = self.reduced_keys_data
        for upstream_id, downstream_id, geometry in zip(
            data["from_node"], data["to_node"], data["geometry"], strict=True
        ):
            if upstream_id in nodes and upstream_id not in upstream_ends:
                upstream_ends[upstream_id] = line_coords(geometry)[0]
            if downstream_id in nodes and downstream_id not in downstream_ends:
                downstream_ends[downstream_id] = line_coords(geometry)[-1]
        return {**downstream_ends, **upstream_ends}

    def explore_nodes(self, G, direction):
        """Generate sewer network based on supplied direction

        .. deprecated::
            Used by the deprecated :meth:`create_sub_network`.
        """
        warnings.warn(
            "explore_nodes is deprecated; use catchment, from_nodes or from_node_list.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._explore_nodes(G, direction)

    def _explore_nodes(self, G, direction):
        """Generate sewer network based on supplied direction

        Parameters
        ----------
        G : networkx.digraph(),
            data from Network class

        direction : str
            specifying successor or predecessor

        Return
        ----------
        None
        """
        connected_subgraph = {self.start}
        current_depth = 0
        opposite_direction = "successors" if direction == "predecessors" else "predecessors"

        while current_depth < self.depth:
            new_nodes = set()
            for node in sorted(connected_subgraph):
                # Find the next node in the sewer network
                next_nodes = set(getattr(G, direction)(node))
                if all(
                    not any(n in node_set for n in next_nodes) for node_set in connected_subgraph
                ):
                    if len(next_nodes) > 1:
                        print(f"Branch Detected at node: {node}")
                        branch_nodes = self._explore_branch(
                            next_nodes, new_nodes, G, direction, opposite_direction
                        )
                        new_nodes.update(branch_nodes)
                        current_depth += 1
                    elif len(next_nodes) == 1:
                        new_nodes.update(next_nodes)
                        opposite_nodes = set(getattr(G, opposite_direction)(next(iter(next_nodes))))
                        print(
                            f"Checking {opposite_direction.capitalize()} of {next_nodes}: {opposite_nodes}"
                        )
                        if next(iter(opposite_nodes)) == node and len(opposite_nodes) == 1:
                            break
                        else:
                            print(f"Possible branch detected at: {next_nodes}")
                            if next_nodes in connected_subgraph:
                                break
                            else:
                                current_depth += 1
                                inner_depth = current_depth
                                while inner_depth < self.depth:
                                    branch_nodes = self._explore_branch(
                                        opposite_nodes, new_nodes, G, opposite_direction, direction
                                    )
                                    new_nodes.update(branch_nodes)
                                    inner_depth += 1

                    elif len(next_nodes) == 0:
                        print("Endpoint reached")
                        new_nodes.update(next_nodes)
                        current_depth += 1
                        break

            connected_subgraph.update(new_nodes)

        return connected_subgraph

    @staticmethod
    def explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction):
        """Exploring branches found in explore_nodes function

        .. deprecated::
            Used by the deprecated :meth:`create_sub_network`.
        """
        warnings.warn(
            "explore_branch is deprecated; use catchment, from_nodes or from_node_list.",
            DeprecationWarning,
            stacklevel=2,
        )
        return SubNetwork._explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction)

    @staticmethod
    def _explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction):
        """Exploring branches found in explore_nodes function

        Parameters
        ----------
        branch_nodes : set
            nodes where a branch is found

        new_nodes : set
            nodes that will be added to subnetwork

        G : networkx.digraph(),
            data from Network class

        direction : str
            specifying successor or predecessor

        opposite_direction : str
            value opposite to direction

        Return
        ----------
        new_nodes : set
            nodes to be added to the subnetwork from branches
        """
        for branch in sorted(branch_nodes):
            print(f"Exploring branch from node: {branch}")
            new_nodes.add(branch)
            # print(f'{direction.capitalize()} of {branch} is: {list(getattr(G, direction)(branch))}')
            branch_next = set(getattr(G, direction)(branch))

            new_nodes.update(branch_next)
            while len(branch_next) == 1:
                next_node = next(iter(branch_next))
                branch_next = set(getattr(G, direction)(next_node))
                if len(branch_next) >= 2:
                    print(f"New branch detected at {next_node}")
                    break
                elif len(branch_next) == 0:
                    print(f"Endpoint reached at {next_node}")
                    new_nodes.add(next_node)
                    new_nodes.update(branch_next)

                    # Check for successors or predecessors based on the current direction
                    opposite_nodes = set(getattr(G, opposite_direction)(next_node))
                    print(
                        f"Checking {opposite_direction.capitalize()} of {next_node}: "
                        f"{opposite_nodes}"
                    )
                    break
                else:
                    new_nodes.update(branch_next)
        return new_nodes

    def manipulate_network(self, operation, direction, start_nodes, copy=False):
        """Manipulate graph based on specified operation

        Parameters
        ----------
        operation : str
            operation to perform.
            options: trim, extend

        direction : str
            direction to trim or extend.

        start_nodes: list
            node(s) to start trimming or extending direction.

        copy : bool, optional
            manipulate existing graph or return a new one.

        Returns
        ----------
        nx.Digraph
        """

        G = self._initialize_digraph()

        if len(start_nodes) == 0:
            raise ValueError("No nodes listed to trim. At least 1 node needs to be provided.")

        self.validate_direction(direction)

        # if tuple(str(start_nodes)) not in self.subgraph.nodes():
        #    print(f"Node {start_nodes} is not in the current graph. Exiting manipulate_subnetwork.")
        #    return

        if operation == "trim":
            if copy:
                trimmed_graph = self.subgraph.copy()
                for nodes in start_nodes:
                    remove_nodes = set(getattr(trimmed_graph, direction)(nodes))
                    remove_nodes.update(start_nodes)
                    trimmed_graph.remove_nodes_from(remove_nodes)
                return trimmed_graph
            else:
                self.subgraph = nx.DiGraph(self.subgraph)
                remove_nodes = set()
                for nodes in start_nodes:
                    remove_nodes = set(getattr(self.subgraph, direction)(nodes))
                    remove_nodes.update(start_nodes)
                    self.subgraph.remove_nodes_from(remove_nodes)

        elif operation == "extend":
            if copy:
                extended_graph = self.subgraph.copy()
                for nodes in start_nodes:
                    add_nodes = set(getattr(G, direction)(nodes))
                    add_nodes.update(start_nodes)
                    extended_graph.update(add_nodes)
                return extended_graph
            else:
                self.subgraph = nx.DiGraph(self.subgraph)
                for nodes in start_nodes:
                    neighbour = getattr(G, direction)(nodes)
                    add_nodes = set(neighbour)
                    add_nodes.update(start_nodes)
                    print(add_nodes)
                    self.subgraph.add_nodes_from(add_nodes)
        else:
            raise ValueError("Invalid operation specified. Must be: extend or trim")

        # Add the ability to trim / extend branches now

    @staticmethod
    def validate_direction(direction):
        """Validate the direction parameter.

        Parameters:
        ----------
        direction : str
            The direction parameter to be validated.

        Raises
        ----------
        ValueError
            If the direction is not 'successors' or 'predecessors'.
        """

        if direction not in ["successors", "predecessors"]:
            raise ValueError("Invalid direction. Use 'successors' or 'predecessors'.")

    def subnetwork_pipes(self):
        """GIS rows of the pipes in ``subgraph`` (including parallel pipes), sorted by pipe_id."""
        pipes = {p for down, up in self.subgraph.edges for p in self._edge_pipes(up, down)}
        return self.full_data[self.full_data["pipe_id"].isin(pipes)].sort_values("pipe_id")

    def node_inverts(self, height_correction=None):
        """Invert elevation of every node in ``subgraph``.

        A node's invert is the lowest invert of the subnetwork pipes connected to it, including
        the drop-structure inverts of the pipes entering it, so that every pipe sits at or above
        its nodes (the pipes are then attached with offsets). A node with no subnetwork pipe takes
        the lowest invert of the network pipes connected to it.

        Parameters
        ----------
        height_correction : dict, optional
            {node: invert} overriding the computed inverts.

        Returns
        ----------
        dict
            {node: invert}

        Raises
        ----------
        ValueError
            If a node has no pipe with an invert.
        """
        nodes = set(self.subgraph.nodes)
        inverts = self._lowest_inverts(self.subnetwork_pipes(), nodes)
        if nodes - set(inverts):
            inverts.update(self._lowest_inverts(self.full_data, nodes - set(inverts)))
        missing = sorted(nodes - set(inverts))
        if missing:
            raise ValueError(f"No pipe invert found for nodes: {', '.join(map(str, missing))}")
        for node, invert in (height_correction or {}).items():
            if node in nodes:
                inverts[node] = float(invert)
        return inverts

    @staticmethod
    def _lowest_inverts(pipes, nodes):
        """{node: lowest invert of the pipes connected to it} for the nodes in ``nodes``."""
        found = {}
        for up, down, inv_up, inv_down, drop in zip(
            pipes["from_node"],
            pipes["to_node"],
            pipes["us_invert"],
            pipes["ds_invert"],
            pipes["drop_invert"],
            strict=True,
        ):
            for node, value in ((up, inv_up), (down, inv_down), (down, drop)):
                value = _invert(value)
                if node in nodes and value is not None:
                    found[node] = min(found.get(node, value), value)
        return found

    def pumps(self, pumps):
        """Check and resolve the ``pumps`` argument of :meth:`write_pyswmm`.

        Parameters
        ----------
        pumps : dict or None
            {pipe_id: Pump or dict of its fields}. The keys are matched to the subnetwork pipe IDs
            as strings, so ``1010`` and ``'1010'`` are the same pipe.

        Returns
        ----------
        dict
            {pipe_id: Pump}, keyed by the pipe IDs of the data.

        Raises
        ----------
        ValueError
            If a pipe is not in the subnetwork or a pump is invalid.
        """
        if not pumps:
            return {}
        ids = {str(p): p for p in self.subnetwork_pipes()["pipe_id"]}
        missing = sorted(str(k) for k in pumps if str(k) not in ids)
        if missing:
            raise ValueError(f"Pumped pipes not in the subnetwork: {', '.join(missing)}")
        return {ids[str(k)]: Pump.coerce(k, v) for k, v in pumps.items()}

    def wet_wells(self, wet_wells):
        """Check the ``wet_wells`` argument of :meth:`write_pyswmm`.

        Parameters
        ----------
        wet_wells : dict or None
            {node: WetWell or dict of its fields}.

        Returns
        ----------
        dict
            {node: WetWell}

        Raises
        ----------
        ValueError
            If a node is not in the subnetwork or is an outfall, or a wet well is invalid.
        """
        if not wet_wells:
            return {}
        missing = sorted(str(n) for n in wet_wells if n not in self.subgraph)
        if missing:
            raise ValueError(f"Wet wells not in the subnetwork: {', '.join(missing)}")
        outfalls = sorted(str(n) for n in wet_wells if n in set(self._outfalls()))
        if outfalls:
            raise ValueError(f"Wet wells cannot be outfalls: {', '.join(outfalls)}")
        return {node: WetWell.coerce(node, v) for node, v in wet_wells.items()}

    def xsections(self, pipes):
        """{pipe_id: (shape, geom1, geom2)} for the rows in ``pipes`` (see :func:`xsection`).

        Raises
        ----------
        ValueError
            If any pipe has an unknown or missing shape code, or is missing a size; the message
            lists every such pipe.
        """
        sections, bad = {}, []
        for _, pipe in pipes.iterrows():
            section = xsection(pipe)
            if section is None:
                bad.append(f"{pipe['pipe_id']} (shape code {pipe.get('shape_code')!r})")
            sections[pipe["pipe_id"]] = section
        if bad:
            raise ValueError(
                "Cannot build a cross-section for pipes with an unknown or missing shape code or "
                f"size: {', '.join(bad)}. Known shape codes in schema {self.schema.name!r}: "
                f"{list(self.schema.shapes)}."
            )
        return sections

    def write_nodes_pyswmm(
        self,
        pyswmm_file,
        max_depth,
        init_depth,
        surcharge_depth,
        ponding_area,
        height_correction,
        extra_outfalls=(),
        wet_wells=None,
    ):
        """Add subnetwork data into the Nodes section of a .inp file.

        Each outfall is written once; wet wells are written to [STORAGE] and every other node is a
        junction. Node inverts come from :meth:`node_inverts`.

        Parameters
        ----------
        pyswmm_file : str
            Path to the write out PySWMM input file (.inp).

        max_depth : float)
            Maximum depth for node (default is 0.0).

        init_depth : float
            Initial depth for node (default is 0.0).

        surcharge_depth : float
            Surcharge depth for node (default is 0.0).

        ponding_area : float
            Ponding area for node (default is 0.0).

        height_correction : dict
            User defined node inverts, overriding those computed by :meth:`node_inverts`.

        extra_outfalls : iterable of str
            Further nodes to write as FREE outfalls, e.g. dead ends.

        wet_wells : dict, optional
            {node: WetWell or dict}, written to the {storage} placeholder (see :meth:`wet_wells`).

        Notes
        -----
        Also fills the {coordinates} placeholder with each node's position in the units of the
        GIS geometry (unscaled; see node_positions).
        """

        with open(pyswmm_file) as file:
            content = file.read()

        outfall_types = dict(zip(self._outfalls(), self._outfall_types(), strict=True))
        outfall_types.update({node: "FREE" for node in extra_outfalls})
        inverts = self.node_inverts(height_correction)
        wet_wells = self.wet_wells(wet_wells)

        nodes_content = ""
        outfall_content = ""
        storage_content = ""
        # Iterate through the nodes in the graph and build the content
        for node in self.subgraph.nodes:
            invert = inverts[node]
            if height_correction is not None and node in height_correction:
                print(f"Node {node}: Updating invert value to: {invert}")

            if node in outfall_types:
                outfall_content += (
                    f"{node:<16} {invert:<11}  {outfall_types[node]:<28} {'NO':<10}\n"
                )
            elif node in wet_wells:
                well = wet_wells[node]
                # FUNCTIONAL area = 0 * depth^0 + area: a constant plan area.
                storage_content += (
                    f"{node:<16} {invert:<11} {well.max_depth:<11} {well.init_depth:<11} "
                    f"{'FUNCTIONAL':<12} {0:<11} {0:<11} {well.area:<11} {0:<11} {0:<11}\n"
                )
            else:
                nodes_content += (
                    f"{node:<16} {invert:<11} {max_depth:<11} {init_depth:<11}"
                    f"{surcharge_depth:<11} {ponding_area:<11}\n"
                )

        coordinates_content = ""
        positions = self.node_positions(self.subgraph.nodes)
        for node in self.subgraph.nodes:
            if node in positions:
                x, y = positions[node]
                coordinates_content += f"{node:<16} {x:<18} {y:<18}\n"
            else:
                print(f"Warning: Node {node} has no position; it is left out of [COORDINATES]")

        # Replace {nodes} placeholder with the built nodes content
        modified_content = (
            content.replace("{nodes}", nodes_content)
            .replace("{outfalls}", outfall_content)
            .replace("{storage}", storage_content)
            .replace("{coordinates}", coordinates_content)
        )

        # Write the modified content back to the file
        with open(pyswmm_file, "w") as file:
            file.write(modified_content)

    def write_conduits_pyswmm(
        self,
        pyswmm_file,
        manning_n,
        init_flow,
        max_flow,
        barrels,
        manning_n_by_material=None,
        height_correction=None,
        pumps=None,
    ):
        """Add subnetwork data into the Conduits and XSection section of a .inp file.

        Parameters
        ----------

        pyswmm_file : str
            Path to the write out PySWMM input file (.inp).

        manning_n : float
            Manning N coefficient for conduit (default is 0.01).

        init_flow : float
            Initial flow value for conduit (default is 0.0).

        max_flow : float
            Maximum flow value for conduit (default is 0.0).

        barrels : float
            Barrel value for conduit (default is 0.0).

        manning_n_by_material : dict, optional
            material -> Manning N. Pipes whose material is not in the mapping use ``manning_n``.

        height_correction : dict, optional
            User defined node inverts, overriding those computed by :meth:`node_inverts`.

        pumps : dict, optional
            {pipe_id: Pump or dict} (see :meth:`pumps`). These pipes are written to the {pumps}
            placeholder, and their curves to {curves}, instead of as conduits.

        Notes
        -----
        The inlet and outlet offsets are the pipe's GIS inverts minus its nodes' inverts
        (LINK_OFFSETS DEPTH), so each conduit has the GIS pipe inverts. Conduits that slope uphill
        are reported with a warning. Also fills the {vertices} placeholder with the interior
        vertices of multi-vertex pipes, in the units of the GIS geometry (unscaled).
        """

        with open(pyswmm_file) as file:
            content = file.read()

        pipes = self.subnetwork_pipes()
        pumps = self.pumps(pumps)
        sections = self.xsections(pipes[~pipes["pipe_id"].isin(list(pumps))])
        inverts = self.node_inverts(height_correction)

        conduits_content = ""
        xsection_content = ""
        pumps_content = ""
        curves_content = ""
        vertices_content = ""
        adverse = []

        for _, row in pipes.iterrows():
            name = row["pipe_id"]
            inlet = row["from_node"]
            outlet = row["to_node"]
            length = row["length"]

            for x, y in line_coords(row["geometry"])[1:-1]:
                vertices_content += f"{name:<16} {x:<18} {y:<18}\n"

            if name in pumps:
                pump = pumps[name]
                curve = "*" if pump.curve == "*" else f"{name}_curve"
                pumps_content += (
                    f"{name:<16} {inlet:<16} {outlet:<16} {curve:<16} {pump.status:<6} "
                    f"{pump.startup:<11} {pump.shutoff:<11}\n"
                )
                if curve != "*":
                    for i, (x, y) in enumerate(pump.curve):
                        kind = pump.curve_type if i == 0 else ""
                        curves_content += f"{curve:<16} {kind:<9} {x:<11} {y:<11}\n"
                    curves_content += "\n"
                continue

            roughness = (manning_n_by_material or {}).get(row["material"], manning_n)
            shape, geom1, geom2 = sections[name]

            ends = []
            for node, pipe_invert in (
                (inlet, row["us_invert"]),
                (outlet, row["ds_invert"]),
            ):
                pipe_invert = _invert(pipe_invert)
                offset = 0.0 if pipe_invert is None else round(pipe_invert - inverts[node], 4)
                if offset < 0:
                    warnings.warn(
                        f"Pipe {name}: invert {pipe_invert} is below the invert of node {node} "
                        f"({inverts[node]}); its offset there is set to 0.",
                        stacklevel=2,
                    )
                    offset = 0.0
                ends.append((offset, inverts[node] + offset))
            (inlet_offset, inlet_elev), (outlet_offset, outlet_elev) = ends
            if outlet_elev > inlet_elev:
                adverse.append(f"{name} ({inlet} {inlet_elev:.3f} -> {outlet} {outlet_elev:.3f})")

            conduits_content += (
                f"{name:<16} {inlet:<16} {outlet:<16}"
                f"{length:<11} {roughness:<11} "
                f"{inlet_offset:<11} {outlet_offset:<11}"
                f"{init_flow:<11} {max_flow:<11}\n"
            )

            xsection_content += (
                f"{name:<16} {shape:<16} {geom1:<16}{geom2:<11} {0.0:<11} {0.0:<11} {barrels:<11}\n"
            )

        if adverse:
            warnings.warn(
                f"Conduits with an adverse slope (outlet invert above inlet invert): "
                f"{', '.join(adverse)}.",
                stacklevel=2,
            )

        # Replace variable in template file with conduit information
        modified_content = (
            content.replace("{conduits}", conduits_content)
            .replace("{xsections}", xsection_content)
            .replace("{pumps}", pumps_content)
            .replace("{curves}", curves_content)
            .replace("{vertices}", vertices_content)
        )

        with open(pyswmm_file, "w") as file:
            file.write(modified_content)

    def write_pyswmm(
        self,
        output_file=None,
        max_depth=0.0,
        init_depth=0.0,
        surcharge_depth=0.0,
        ponding_area=0.0,
        manning_n=0.01,
        manning_n_by_material=None,
        inlet_offset=None,
        outlet_offset=None,
        init_flow=0.0,
        max_flow=0.0,
        barrels=1.0,
        template_path=None,
        output_folder=None,
        height_correction=None,
        dead_ends="raise",
        pumps=None,
        wet_wells=None,
    ):
        """Generate a PySWMM input file based on a networkx digraph and GeoPandas DataFrame.

        The subnetwork is checked with :meth:`validate` first. Flow splits (including parallel
        pipes) and cycles are refused unless the template uses ``FLOW_ROUTING DYNWAVE``.

        Node inverts are the lowest invert of the pipes connected to each node (see
        :meth:`node_inverts`), and conduits are attached with offsets so that they keep their GIS
        inverts. Pipe cross-sections come from the shape codes, mapped by the schema's ``shapes``.

        Parameters
        ----------
        output_file : str
            Name of new pyswmm input file (.inp).

        max_depth : float
            Maximum depth for node (default is 0.0).

        init_depth : float
            Initial depth for node (default is 0.0).

        surcharge_depth : float
            Surcharge depth for node (default is 0.0).

        ponding_area : float
            Ponding area for node (default is 0.0).

        manning_n : float
            Manning N coefficient for conduit (default is 0.01). With ``manning_n_by_material``,
            the fallback for pipes whose material is not in the mapping.

        manning_n_by_material : dict, optional
            material -> Manning N, e.g. ``{'PVC': 0.011, 'CONC': 0.013, 'VIT': 0.013}``.

        inlet_offset, outlet_offset : float
            Deprecated and ignored: the offsets are computed from the GIS inverts.

        init_flow : float
            Initial flow value for conduit (default is 0.0).

        max_flow : float
            Maximum flow value for conduit (default is 0.0).

        barrels : float
            Barrel value for conduit (default is 0.0).

        template_path : str
            Path to the pyswmm model template file. The default is the template shipped with the
            package (``sewernetworkx/templates/template.inp``), which uses SI units (FLOW_UNITS
            LPS). A custom template must use an SI FLOW_UNITS (CMS, LPS or MLD), since the GIS
            lengths, inverts and sizes are in metres, and may only contain the placeholders
            {nodes}, {outfalls}, {storage}, {conduits}, {pumps}, {xsections}, {curves},
            {coordinates} and {vertices}; {storage}, {pumps} and {curves} are needed only with
            ``wet_wells`` or ``pumps``. The offsets are depths, so the template should set LINK_OFFSETS DEPTH (the default).

        output_folder : str
            Path to were the new .inp file will be written.

        height_correction : dict
            {node: invert} overriding the computed node inverts, e.g. where the GIS inverts are
            wrong. Conduit offsets are measured from the corrected inverts.

        dead_ends : {'raise', 'outfall'}
            What to do with nodes that have no outlet pipe and are not outfalls: raise a
            ``ValueError`` (default) or write them as FREE outfalls.

        pumps : dict, optional
            {pipe_id: :class:`Pump` or dict of its fields}, e.g.
            ``{1010: dict(curve=[(0, 5.0), (5, 3.0), (10, 0.0)], curve_type='PUMP3')}``. Each pipe
            is written as a pump from its upstream to its downstream node, replacing its conduit,
            and its curve to [CURVES]. Curve flows are in the template's FLOW_UNITS (L/s with the
            shipped template). The GIS data do not identify pumps, so they must be given here.

        wet_wells : dict, optional
            {node: :class:`WetWell` or dict of its fields}, e.g. ``{'D1': dict(max_depth=3.0,
            area=4.0)}``. These nodes are written as [STORAGE] nodes of constant plan area
            instead of junctions. The inlet of a PUMP1 pump must be a wet well.

        Raises
        ----------
        ValueError
            If the topology is invalid for the template's flow routing, there are dead ends and
            ``dead_ends`` is 'raise', a pipe has an unknown or missing shape code or size, a
            node has no invert, a pump or wet well is invalid, or the template lacks a placeholder
            that ``pumps`` or ``wet_wells`` needs.
        """

        if output_file is None:
            raise ValueError("The 'output_path' must be specified.")

        if inlet_offset is not None or outlet_offset is not None:
            warnings.warn(
                "inlet_offset and outlet_offset are ignored: conduit offsets are computed from the "
                "GIS pipe inverts.",
                DeprecationWarning,
                stacklevel=2,
            )

        if template_path is None:
            template = default_template_path().read_text()
        else:
            if not os.path.exists(template_path):
                raise FileNotFoundError(
                    f"The template file does not exist at the specified path: {template_path}"
                )
            with open(template_path) as file:
                template = file.read()

        check_template(template)

        if dead_ends not in ("raise", "outfall"):
            raise ValueError("dead_ends must be 'raise' or 'outfall'.")
        report = self.validate()
        self._check_topology(report, self._template_routing(template), dead_ends)

        # Fail before writing anything if a pipe or node cannot be represented
        resolved_pumps = self.pumps(pumps)
        resolved_wells = self.wet_wells(wet_wells)
        self._check_pumps(template, resolved_pumps, resolved_wells)
        pipes = self.subnetwork_pipes()
        self.xsections(pipes[~pipes["pipe_id"].isin(list(resolved_pumps))])
        self.node_inverts(height_correction)

        # Always create a copy of the template file
        if output_folder is not None:
            output_file_path = os.path.join(output_folder, output_file)
        else:
            # If output_folder is None, create a copy in the current working directory
            output_file_path = os.path.join(os.getcwd(), output_file)

        with open(output_file_path, "w") as file:
            file.write(template)
        print(f"Template file copied to: {output_file_path}")

        # Write nodes to pyswmm .inp file
        self.write_nodes_pyswmm(
            output_file_path,
            max_depth,
            init_depth,
            surcharge_depth,
            ponding_area,
            height_correction,
            extra_outfalls=report.dead_ends if dead_ends == "outfall" else (),
            wet_wells=resolved_wells,
        )

        # Write conduits to pyswmm inp file
        self.write_conduits_pyswmm(
            output_file_path,
            manning_n,
            init_flow,
            max_flow,
            barrels,
            manning_n_by_material,
            height_correction,
            pumps=resolved_pumps,
        )

    def _check_pumps(self, template, pumps, wet_wells):
        """Raise ValueError if the template cannot hold the pumps and wet wells, or a PUMP1 pump
        has no wet well at its inlet."""
        needed = []
        if pumps:
            needed.append("pumps")
            if any(p.curve != "*" for p in pumps.values()):
                needed.append("curves")
        if wet_wells:
            needed.append("storage")
        missing = [f"{{{name}}}" for name in needed if f"{{{name}}}" not in template]
        if missing:
            raise ValueError(
                f"The template has no {', '.join(missing)} placeholder, so the pumps or wet wells "
                "would not be written."
            )

        inlets = self.subnetwork_pipes().set_index("pipe_id")["from_node"]
        dry = sorted(
            f"{name} (inlet {inlets[name]})"
            for name, pump in pumps.items()
            if pump.curve != "*" and pump.curve_type == "PUMP1" and inlets[name] not in wet_wells
        )
        if dry:
            raise ValueError(
                f"PUMP1 pumps need a wet well at their inlet node (pass wet_wells): {', '.join(dry)}"
            )

    @staticmethod
    def _template_routing(template):
        """FLOW_ROUTING of the text of an .inp template; SWMM's default, KINWAVE, if it is not set."""
        match = re.search(r"^\s*FLOW_ROUTING\s+(\w+)", template, flags=re.MULTILINE | re.IGNORECASE)
        return match.group(1).upper() if match else "KINWAVE"

    @staticmethod
    def _check_topology(report, routing, dead_ends):
        problems = []
        if routing != "DYNWAVE":
            for node, pipes in report.splits.items():
                problems.append(
                    f"node {node} has more than one outlet pipe ({', '.join(map(str, pipes))})"
                )
            for cycle in report.cycles:
                problems.append(f"cycle {' -> '.join(cycle + cycle[:1])}")
        if dead_ends == "raise" and report.dead_ends:
            problems.append(
                f"dead ends with no outlet that are not outfalls: {', '.join(report.dead_ends)} "
                "(pass dead_ends='outfall' to write them as outfalls)"
            )
        if not problems:
            return
        message = "Invalid subnetwork topology:\n  " + "\n  ".join(problems)
        if routing != "DYNWAVE" and not report.is_dendritic:
            message += (
                f"\n{routing} routing requires a dendritic network; use a template with "
                "FLOW_ROUTING DYNWAVE for flow splits, parallel pipes and loops."
            )
        raise ValueError(message)

    def generate_gmap(self, file_name, scale=1, origin=None):
        """Generate an hmtl file to show the sewer network overlaid on Google Maps.

        Parameters
        ----------
        file_name : str
            Name of the html file

        scale : int
            Value to scale coordinate data by. Default value is 1.

        origin : list
            Location of centre for the overlaid map.
        """

        self.validate_file_extension(file_name, ".html")

        if origin is not None:
            gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], zoom=18)
        else:
            if not self.coords:
                raise ValueError(
                    "Coordinates have not been generated for this graph. "
                    "You must run the function '.add_coordinates' prior to "
                    "running 'generate.gmap'."
                )

            coords = next(iter(self.coords.values()))
            gmap = gmplot.GoogleMapPlotter(coords[1] / scale, coords[0] / scale, zoom=18)

        for node, coords in self.coords.items():
            gmap.marker(coords[1] / scale, coords[0] / scale, "cornflowerblue", label=str(node))

        for edge in self.subgraph.edges():
            start_node_pos = self.coords[edge[0]]
            end_node_pos = self.coords[edge[1]]
            gmap.plot(
                [start_node_pos[1] / scale, end_node_pos[1] / scale],
                [start_node_pos[0] / scale, end_node_pos[0] / scale],
                "red",
                edge_width=2,
            )

        gmap.draw(file_name)

    @staticmethod
    def validate_file_extension(file_name, expected_extension):
        _, file_extension = os.path.splitext(file_name)

        if file_extension != expected_extension:
            raise ValueError(
                f"Wrong file extension used. The proper extension is: {expected_extension}."
            )
