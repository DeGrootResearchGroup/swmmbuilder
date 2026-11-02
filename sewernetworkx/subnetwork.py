"""A subnetwork: a selection of nodes of a :class:`~sewernetworkx.Network` and the pipes between
them."""

from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx

from . import swmm
from .geometry import line_coords


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


class Subnetwork:
    """The nodes ``nodes`` of ``network`` and every pipe between two of them.

    Usually made with :meth:`Network.catchment` or :meth:`Network.subnetwork`.

    Parameters
    ----------
    network : Network
        The network.

    nodes : iterable of str
        Node IDs.

    Attributes
    ----------
    network : Network
        The network the subnetwork belongs to.

    graph : networkx.MultiDiGraph
        The induced subgraph of ``network.graph``: edges in the direction of flow, keyed by
        pipe_id (parallel pipes are separate edges), nodes and edges in sorted order.

    Raises
    ----------
    ValueError
        If a node is not in the network.
    """

    def __init__(self, network, nodes):
        nodes = sorted(set(nodes))
        network._check_nodes(nodes)
        self.network = network
        # Built explicitly: networkx's subgraph views iterate their nodes in set order.
        selected = set(nodes)
        self.graph = nx.MultiDiGraph()
        self.graph.add_nodes_from(nodes)
        self.graph.add_edges_from(
            (up, down, key, attrs)
            for up, down, key, attrs in network.graph.edges(nodes, keys=True, data=True)
            if down in selected
        )

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
        lines = ["# SewerNetworkx node list: one node ID per line", *self.nodes]
        Path(path).write_text("\n".join(lines) + "\n")

    @property
    def nodes(self):
        """The node IDs, sorted."""
        return list(self.graph.nodes)

    @property
    def schema(self):
        return self.network.schema

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
        nodes = set(self.nodes if nodes is None else nodes)
        upstream_ends, downstream_ends = {}, {}
        data = self.network.full_data
        for upstream_id, downstream_id, geometry in zip(
            data["from_node"], data["to_node"], data["geometry"], strict=True
        ):
            if upstream_id in nodes and upstream_id not in upstream_ends:
                upstream_ends[upstream_id] = line_coords(geometry)[0]
            if downstream_id in nodes and downstream_id not in downstream_ends:
                downstream_ends[downstream_id] = line_coords(geometry)[-1]
        return {**downstream_ends, **upstream_ends}

    def validate(self, outfalls=()):
        """Check the topology before the subnetwork is written to an .inp file.

        Parameters
        ----------
        outfalls : iterable of str
            The outfall nodes, which need no outlet pipe.

        Returns
        ----------
        TopologyReport
            Flow splits, cycles, dead ends and parallel pipes. Splits and cycles are only allowed
            with dynamic-wave routing; dead ends are never valid.
        """
        outfalls = set(outfalls)
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

        return report

    def to_inp(self, path, outfalls, **options):
        """Write the subnetwork to a SWMM .inp file; see :func:`sewernetworkx.swmm.write_inp`.

        Returns
        ----------
        SwmmModel
            The model that was written.
        """
        return swmm.write_inp(self, path, outfalls, **options)
