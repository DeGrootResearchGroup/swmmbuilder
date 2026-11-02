"""Module for manipulating GeoJSON sewer data.

The module takes full GeoJSON data sets of a city sewer network,
generates a subnetwork given an initial starting node and number
of levels, and generates a .inp file to be used by pyswmm.
"""

import logging
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import geopandas as gpd
import gmplot
import networkx as nx
import pandas as pd

from . import swmm
from .geometry import line_coords
from .schema import CANONICAL_COLUMNS, LONDON_ON, as_schema

logger = logging.getLogger(__name__)

# Canonical flow types that Network keeps by default: sanitary sewers only.
DEFAULT_FLOW_TYPES = frozenset({"sanitary"})


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
                logger.info(
                    "Removed %d pipes with flow type %s (keeping %s)",
                    count,
                    flow_type,
                    sorted(flow_types),
                )

        if not keep_abandoned and schema.abandoned is not None:
            abandoned = self.full_data["abandoned"]
            self.removed["abandoned"] = int(abandoned.sum())
            self.full_data = self.full_data[~abandoned]
            logger.info("Removed %d pipes flagged abandoned", self.removed["abandoned"])

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
        self.coords = {}

        if len(self.outfall_node) != len(self.outfall_node_type):
            raise ValueError(
                f"Length of {self.outfall_node} and {self.outfall_node_type} are not the same."
            )

        self.schema = source.schema
        self.full_data = source.full_data
        self.graph = source.graph
        self.reduced_keys_data = source.full_data[[*CANONICAL_COLUMNS, "geometry"]].copy()

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
        G = self._initialize_digraph()

        # Generate subnetwork starting with successors of starting node
        successors_subgraph = self._explore_nodes(G, direction="successors")

        # Generate subnetwork with predecessors of starting node
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
            logger.warning("Positions for nodes %s could not be found.", sorted(missing_nodes))

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
                        logger.debug("Branch detected at node %s", node)
                        branch_nodes = self._explore_branch(
                            next_nodes, new_nodes, G, direction, opposite_direction
                        )
                        new_nodes.update(branch_nodes)
                        current_depth += 1
                    elif len(next_nodes) == 1:
                        new_nodes.update(next_nodes)
                        opposite_nodes = set(getattr(G, opposite_direction)(next(iter(next_nodes))))
                        logger.debug("%s of %s: %s", opposite_direction, next_nodes, opposite_nodes)
                        if next(iter(opposite_nodes)) == node and len(opposite_nodes) == 1:
                            break
                        else:
                            logger.debug("Possible branch detected at %s", next_nodes)
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
                        logger.debug("Endpoint reached at %s", node)
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
            logger.debug("Exploring branch from node %s", branch)
            new_nodes.add(branch)
            branch_next = set(getattr(G, direction)(branch))

            new_nodes.update(branch_next)
            while len(branch_next) == 1:
                next_node = next(iter(branch_next))
                branch_next = set(getattr(G, direction)(next_node))
                if len(branch_next) >= 2:
                    logger.debug("New branch detected at %s", next_node)
                    break
                elif len(branch_next) == 0:
                    logger.debug("Endpoint reached at %s", next_node)
                    new_nodes.add(next_node)
                    new_nodes.update(branch_next)

                    # Check for successors or predecessors based on the current direction
                    opposite_nodes = set(getattr(G, opposite_direction)(next_node))
                    logger.debug("%s of %s: %s", opposite_direction, next_node, opposite_nodes)
                    break
                else:
                    new_nodes.update(branch_next)
        return new_nodes

    def subnetwork_pipes(self):
        """GIS rows of the pipes in ``subgraph`` (including parallel pipes), sorted by pipe_id."""
        pipes = {p for down, up in self.subgraph.edges for p in self._edge_pipes(up, down)}
        return self.full_data[self.full_data["pipe_id"].isin(pipes)].sort_values("pipe_id")

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
    ):
        """Write the subnetwork to ``output_folder/output_file`` with :func:`swmm.write_inp`.

        The subnetwork is checked with :meth:`validate` first. Flow splits (including parallel
        pipes) and cycles are refused unless the template uses ``FLOW_ROUTING DYNWAVE``.

        Node inverts are the lowest invert of the pipes connected to each node (see
        :func:`swmm.node_inverts`), and conduits are attached with offsets so that they keep their
        GIS inverts. Pipe cross-sections come from the shape codes, mapped by the schema's
        ``shapes``. The outfalls are ``outfall_node``, with the types ``outfall_node_type``.

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
            {nodes}, {outfalls}, {conduits}, {xsections}, {coordinates} and {vertices}.
            Pumps are not supported yet, so a template with a {pumps} placeholder is rejected.
            The offsets are depths, so the template should set LINK_OFFSETS DEPTH (the default).

        output_folder : str
            Path to were the new .inp file will be written.

        height_correction : dict
            {node: invert} overriding the computed node inverts, e.g. where the GIS inverts are
            wrong. Conduit offsets are measured from the corrected inverts.

        dead_ends : {'raise', 'outfall'}
            What to do with nodes that have no outlet pipe and are not outfalls: raise a
            ``ValueError`` (default) or write them as FREE outfalls.

        Raises
        ----------
        ValueError
            If the topology is invalid for the template's flow routing, there are dead ends and
            ``dead_ends`` is 'raise', a pipe has an unknown or missing shape code or size, or a
            node has no invert.
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

        folder = os.getcwd() if output_folder is None else output_folder
        swmm.write_inp(
            self,
            os.path.join(folder, output_file),
            template_path=template_path,
            dead_ends=dead_ends,
            max_depth=max_depth,
            init_depth=init_depth,
            surcharge_depth=surcharge_depth,
            ponding_area=ponding_area,
            manning_n=manning_n,
            manning_n_by_material=manning_n_by_material,
            init_flow=init_flow,
            max_flow=max_flow,
            barrels=barrels,
            height_correction=height_correction,
        )

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
