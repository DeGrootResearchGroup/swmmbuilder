"""A subnetwork: a selection of nodes of a :class:`~swmmbuilder.Network` and the pipes between
them."""

from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx

from . import gaps, swmm
from .geometry import node_positions, to_geodataframe


@dataclass
class TopologyReport:
    """Topology problems of a subnetwork, as returned by :meth:`Subnetwork.validate`.

    Attributes
    ----------
    splits : dict
        Node -> sorted pipe IDs of its outlet pipes, for every node with more than one outlet
        pipe. Parallel pipes count separately, as they do in SWMM.

    cycles : list of list
        Each directed cycle as a list of nodes in flow order, starting from its smallest node.

    dead_ends : list
        Nodes with no outlet pipe that are not outfalls.

    dead_end_types : dict
        Dead end -> ``'capped'``, ``'outfall'`` (an outfall node that was not passed as an
        outfall), ``'pump_station'`` or ``'unexplained'`` (a likely data gap), from the node type.

    parallel_pipes : dict
        (upstream, downstream) -> sorted pipe IDs, for every node pair joined by more than one pipe.

    force_mains : list
        Sorted pipe IDs of the force mains (pressurised pipes downstream of a pump).

    stations_without_force_main : list
        Pump-station nodes with no force-main outlet pipe.

    orphan_force_mains : list
        Sorted pipe IDs of the force mains that start a chain (see
        :meth:`Subnetwork.force_main_chain`) at a node that is not a pump station.

    force_main_ends : dict
        First pipe ID of a force-main chain -> (node, ``'force_main'`` or ``'dead_end'``), for the
        chains that join another force main or end at a dead end.

    pipe_count_mismatches : dict
        Node -> {'pipes_in': (node layer, graph), 'pipes_out': (node layer, graph)}, for every node
        whose pipe counts in its node layer differ from those of the network's graph (only the
        counts that differ are listed).

    data_gaps : dict
        Gap kind -> number of gaps of that kind in the data, with no fills (see
        :meth:`Subnetwork.data_gaps`), for the kinds that occur.
    """

    splits: dict = field(default_factory=dict)
    cycles: list = field(default_factory=list)
    dead_ends: list = field(default_factory=list)
    dead_end_types: dict = field(default_factory=dict)
    parallel_pipes: dict = field(default_factory=dict)
    force_mains: list = field(default_factory=list)
    stations_without_force_main: list = field(default_factory=list)
    orphan_force_mains: list = field(default_factory=list)
    force_main_ends: dict = field(default_factory=dict)
    pipe_count_mismatches: dict = field(default_factory=dict)
    data_gaps: dict = field(default_factory=dict)

    @property
    def is_dendritic(self):
        """True when the network has no splits and no cycles, as kinematic-wave routing requires."""
        return not self.splits and not self.cycles


class Subnetwork:
    """The nodes ``nodes`` of ``network`` and every pipe between two of them.

    Usually made with :meth:`Network.catchment` or :meth:`Network.subnetwork`.

    Parameters
    ----------
    network : Network
        The network.

    nodes : iterable of str
        Node IDs (see :meth:`Network.resolve`).

    Attributes
    ----------
    network : Network
        The network the subnetwork belongs to.

    graph : networkx.MultiDiGraph
        The induced subgraph of ``network.graph``: edges in the direction of flow, keyed by
        pipe_id (parallel pipes are separate edges), nodes and edges in sorted order, with the
        attributes of ``network.graph``.

    Raises
    ----------
    ValueError
        If a node is not in the network.
    """

    def __init__(self, network, nodes):
        nodes = sorted({network.resolve(n) for n in nodes})
        network._check_nodes(nodes)
        self.network = network
        # Built explicitly: networkx's subgraph views iterate their nodes in set order.
        selected = set(nodes)
        self.graph = nx.MultiDiGraph()
        self.graph.add_nodes_from((n, network.graph.nodes[n]) for n in nodes)
        self.graph.add_edges_from(
            (up, down, key, attrs)
            for up, down, key, attrs in network.graph.edges(nodes, keys=True, data=True)
            if down in selected
        )
        self._ends = {pipe: (up, down) for up, down, pipe in self.graph.edges(keys=True)}

    @classmethod
    def from_node_list(cls, network, path):
        """The subnetwork of ``network`` whose nodes are listed in a file written by
        :meth:`to_node_list`. Blank lines and lines starting with '#' are ignored."""
        lines = (line.strip() for line in Path(path).read_text().splitlines())
        return cls(network, (line for line in lines if line and not line.startswith("#")))

    def to_node_list(self, path):
        """Write the nodes to a text file, one per line, in sorted order.

        :meth:`from_node_list` rebuilds the same subnetwork from this file.
        """
        lines = ["# swmmbuilder node list: one node ID per line", *self.nodes]
        Path(path).write_text("\n".join(lines) + "\n")

    @property
    def nodes(self):
        """The node IDs, sorted."""
        return list(self.graph.nodes)

    @property
    def schema(self):
        return self.network.schema

    @property
    def outfalls(self):
        """The nodes whose node type is ``'outfall'``, sorted: the default outfalls of
        :meth:`validate` and :meth:`to_inp`."""
        return [n for n, t in self.graph.nodes(data="node_type") if t == "outfall"]

    def pipes(self):
        """The rows of ``network.full_data`` of the pipes in the subnetwork, sorted by pipe_id."""
        data = self.network.full_data
        ids = {key for _, _, key in self.graph.edges(keys=True)}
        return data[data["pipe_id"].isin(ids)].sort_values("pipe_id")

    def node_positions(self, nodes=None):
        """Unscaled (x, y) position of each node (default: every node), from the pipe geometry.

        A node is placed at the upstream end of a pipe of the network leaving it or, if no pipe
        leaves it, at the downstream end of a pipe entering it. Nodes with no pipe are left out.
        """
        return node_positions(self.network.full_data, self.nodes if nodes is None else nodes)

    def to_geodataframe(self, kind="pipes", values=None):
        """The subnetwork's pipes or nodes as a GeoDataFrame; see
        :meth:`Network.to_geodataframe`."""
        return to_geodataframe(self.pipes(), self.network.nodes.loc[self.nodes], kind, values)

    def _is_force_main(self, pipe):
        up, down = self._ends[pipe]
        return self.graph.edges[up, down, pipe].get("structure_type") == "force_main"

    def force_main_outlets(self, node):
        """Sorted pipe IDs of the force mains leaving ``node``."""
        return sorted(
            pipe
            for _, _, pipe, structure in self.graph.out_edges(
                node, keys=True, data="structure_type"
            )
            if structure == "force_main"
        )

    def force_main_chain(self, pipe, outfalls=None):
        """The force main that starts with ``pipe``, followed downstream.

        A force main is often a chain of several GIS pipes. The chain continues through every node
        whose only outlet pipe is a force main, and stops where it joins another force main.

        Parameters
        ----------
        pipe : pipe ID
            The first pipe of the chain.

        outfalls : iterable of str, optional
            The outfall nodes (default: :attr:`outfalls`).

        Returns
        ----------
        pipes : list
            The pipe IDs of the chain, in flow order.

        end : str
            The node where the chain ends.

        problem : str or None
            ``'force_main'`` if the chain joins another force main (a node with more than one
            force main entering it), ``'dead_end'`` if it ends at a node with no outlet that is not
            an outfall, else None.
        """
        outfalls = set(self.outfalls if outfalls is None else outfalls)
        pipes = [pipe]
        node = self._ends[pipe][1]
        while True:
            entering = [p for _, _, p in self.graph.in_edges(node, keys=True)]
            if sum(self._is_force_main(p) for p in entering) > 1:
                return pipes, node, "force_main"
            outlets = [p for _, _, p in self.graph.out_edges(node, keys=True)]
            if not outlets and node not in outfalls:
                return pipes, node, "dead_end"
            if len(outlets) != 1 or not self._is_force_main(outlets[0]) or outlets[0] in pipes:
                return pipes, node, None
            pipes.append(outlets[0])
            node = self._ends[outlets[0]][1]

    def validate(self, outfalls=None):
        """Check the topology before the subnetwork is written to an .inp file.

        Parameters
        ----------
        outfalls : iterable of str, optional
            The outfall nodes, which need no outlet pipe (default: :attr:`outfalls`).

        Returns
        ----------
        TopologyReport
            Flow splits, cycles, dead ends (classified by node type) and parallel pipes, plus the
            force mains (with the pump stations that have none, the force mains that start at no
            pump station and the chains that end badly), the nodes whose recorded pipe counts
            differ from the graph and the counts of gaps in the data. Splits and cycles are only
            allowed with dynamic-wave routing; dead ends are never valid.
        """
        outfalls = set(self.outfalls if outfalls is None else outfalls)
        report = TopologyReport()

        for node in self.nodes:
            outlets = {}
            for _, down, pipe in self.graph.out_edges(node, keys=True):
                outlets.setdefault(down, []).append(pipe)
            pipes = sorted(p for down in outlets for p in outlets[down])
            if len(pipes) > 1:
                report.splits[node] = pipes
            if not pipes and node not in outfalls:
                report.dead_ends.append(node)
            for down in sorted(outlets):
                if len(outlets[down]) > 1:
                    report.parallel_pipes[(node, down)] = sorted(outlets[down])

        for cycle in nx.simple_cycles(nx.DiGraph(self.graph)):
            i = cycle.index(min(cycle))
            report.cycles.append(cycle[i:] + cycle[:i])
        report.cycles.sort()

        dead_end_types = {"cap": "capped", "outfall": "outfall", "pump_station": "pump_station"}
        for node in report.dead_ends:
            node_type = self.graph.nodes[node].get("node_type")
            report.dead_end_types[node] = dead_end_types.get(node_type, "unexplained")

        report.force_mains = sorted(
            pipe
            for _, _, pipe, structure in self.graph.edges(keys=True, data="structure_type")
            if structure == "force_main"
        )
        stations = {n for n, t in self.graph.nodes(data="node_type") if t == "pump_station"}
        report.stations_without_force_main = sorted(
            n for n in stations if not self.force_main_outlets(n)
        )
        for pipe in report.force_mains:
            up = self._ends[pipe][0]
            entering = [p for _, _, p in self.graph.in_edges(up, keys=True)]
            if up not in stations and any(self._is_force_main(p) for p in entering):
                continue  # inside a chain, not its start
            if up not in stations:
                report.orphan_force_mains.append(pipe)
            _, end, problem = self.force_main_chain(pipe, outfalls)
            if problem is not None:
                report.force_main_ends[pipe] = (end, problem)

        full = self.network.graph
        for node, attrs in self.graph.nodes(data=True):
            counts = {
                column: (int(attrs[column]), degree)
                for column, degree in (
                    ("pipes_in", full.in_degree(node)),
                    ("pipes_out", full.out_degree(node)),
                )
                if column in attrs and attrs[column] != degree
            }
            if counts:
                report.pipe_count_mismatches[node] = counts

        counts = self.data_gaps()["gap"].value_counts()
        report.data_gaps = {gap: int(counts[gap]) for gap in gaps.GAP_KINDS if gap in counts}
        return report

    def data_gaps(self, **fills):
        """Every gap in the data of the subnetwork, and how the fill rules would fill it.

        A gap is a pipe with an unknown or missing shape code or a missing size, a pipe end with
        no invert, or a node with no invert from any pipe (see :mod:`swmmbuilder.gaps`). With no
        options nothing is filled, except that a round pipe with no diameter whose height equals
        its width takes that value as its diameter.

        Parameters
        ----------
        **fills
            The fill options of :meth:`to_inp` (``height_correction``, ``size_overrides``,
            ``shape_codes``, ``size_fill``, ``invert_fill``, ``headwater_slope``,
            ``force_main_discharge``, and ``pumps``, whose pipes need no cross-section); see
            :func:`swmmbuilder.gaps.fill_gaps`.

        Returns
        ----------
        pandas.DataFrame
            One row per gap: ``element`` ('pipe' or 'node'), ``id``, ``gap``, ``detail``,
            ``filled_from`` (the rule that filled it; NaN if it is not filled), ``value``,
            ``source`` and ``order`` (the order of the fills).
        """
        if "pumps" in fills:
            fills["pumps"] = swmm.resolve_pumps(self, fills["pumps"])
        return gaps.fill_gaps(self, **fills).report

    def to_inp(self, path, outfalls=None, **options):
        """Write the subnetwork to a SWMM .inp file; see :func:`swmmbuilder.swmm.write_inp`.

        Returns
        ----------
        SwmmModel
            The model that was written.
        """
        return swmm.write_inp(self, path, outfalls, **options)
