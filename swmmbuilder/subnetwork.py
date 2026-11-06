"""A subnetwork: a selection of nodes of a :class:`~swmmbuilder.Network` and the pipes between
them."""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np
import pandas as pd

from . import gaps, swmm
from .diagnostics import DataWarning, check_options
from .geometry import node_positions, to_geodataframe

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from typing import Any, TypeVar

    import geopandas as gpd

    from ._types import NodeId, Outfalls, PipeId, StrPath, Values
    from .network import Network
    from .schema import Schema
    from .swmm.records import SwmmModel

    V = TypeVar("V")


logger = logging.getLogger(__name__)


SPLITS = ("equal", "full")

# The share of a node's total below which it counts as not leaving a loop (see _trapped).
EXIT_TOLERANCE = 1e-12


def _trapped(order: list[Any], a: np.ndarray) -> list[Any]:
    """The nodes of a loop (``order``, with ``a[v, u]`` the share of ``u``'s total that goes to
    ``v`` within it) from which nothing leaves the loop: those that cannot reach, along pipes
    carrying a share, a node that sends part of its total out of the loop (or nowhere).

    The totals within the loop solve ``(I - a) x = v``, which has one solution exactly when every
    node can reach such an exit (the shares are then substochastic with spectral radius below 1).
    Unlike a test of the determinant, which shrinks geometrically with the size of the loop, this
    does not depend on the loop's size.
    """
    exits = [i for i in range(len(order)) if a[:, i].sum() < 1.0 - EXIT_TOLERANCE]
    reached = set(exits)
    stack = list(exits)
    while stack:
        v = stack.pop()
        for u in np.flatnonzero(a[v] > 0):
            if int(u) not in reached:
                reached.add(int(u))
                stack.append(int(u))
    return [n for i, n in enumerate(order) if i not in reached]


@dataclass
class UpstreamTotals:
    """The result of :meth:`Subnetwork.accumulate`.

    Attributes
    ----------
    nodes : pandas.Series or pandas.DataFrame
        The total at each node: its own value plus what reaches it from upstream, indexed by
        node (a DataFrame, one column per column of the values, if the values were one).

    pipes : pandas.Series or pandas.DataFrame
        The total carried by each pipe, indexed by pipe ID: its share of the total at its
        upstream node (the whole total with ``split='full'``).

    fractions : pandas.Series
        The share of its upstream node's total that each pipe carries (1 for every pipe with
        ``split='full'``).

    split : str
        ``'equal'``, ``'full'``, or ``'given'``.
    """

    nodes: pd.Series | pd.DataFrame
    pipes: pd.Series | pd.DataFrame
    fractions: pd.Series
    split: str


@dataclass
class PumpStationReport:
    """The result of :meth:`Subnetwork.find_pump_stations`.

    Attributes
    ----------
    stations : list
        Nodes of type 'pump_station' with a force main leaving them: the stations to pump, e.g.
        ``{n: PumpStation(pump=dict(curve='*')) for n in report.stations}``.

    no_force_main : list
        Nodes of type 'pump_station' with no force main leaving them, which flow on by gravity
        (or whose force main is missing from the data); they are written as junctions unless
        given in ``pump_stations``.

    unlisted : pandas.DataFrame
        Likely lift stations missing from the station layer: one row per force main that starts
        at a node not of type 'pump_station', with no force main entering it, and climbs by more
        than ``min_rise``. Columns ``node``, ``pipes`` (the chain, in flow order), ``end``,
        ``us_invert``, ``ds_invert``, and ``rise`` (m).
    """

    stations: list = field(default_factory=list)
    no_force_main: list = field(default_factory=list)
    unlisted: pd.DataFrame = field(default_factory=pd.DataFrame)


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
        chains that join another force main or end at a dead end. The .inp writer keeps the shared
        main downstream of a join as a force main of both stations; a dead end is refused like any
        other.

    pipe_count_mismatches : dict
        Node -> {'pipes_in': (node layer, graph), 'pipes_out': (node layer, graph)}, for every node
        whose pipe counts in its node layer differ from those of the network's graph (only the
        counts that differ are listed).

    data_gaps : dict
        Gap kind -> number of gaps of that kind in the data, with no fills (see
        :meth:`Subnetwork.data_gaps`), for the kinds that occur.

    nodes, pipes : int
        Number of nodes and pipes of the subnetwork.
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
    nodes: int = 0
    pipes: int = 0

    @property
    def is_dendritic(self) -> bool:
        """True when the network has no splits and no cycles, as kinematic-wave routing requires."""
        return not self.splits and not self.cycles

    def summary(self) -> dict[str, int]:
        """Counts of the report, to cite in a report.

        Returns
        -------
        dict
            Plain Python ints, under these keys:

            - ``nodes``, ``pipes``: nodes and pipes of the subnetwork.
            - ``splits``: nodes with more than one outlet pipe (:attr:`splits`).
            - ``cycles``: directed cycles (:attr:`cycles`).
            - ``dead_ends``: nodes with no outlet pipe that are not outfalls.
            - ``parallel_pipes``: node pairs joined by more than one pipe.
            - ``force_mains``: force-main pipes (:attr:`force_mains`).
            - ``stations_without_force_main``: pump stations with no force-main outlet.
            - ``orphan_force_mains``: force mains that start at no pump station.
            - ``force_main_ends``: force-main chains that join another or end at a dead end.
            - ``pipe_count_mismatches``: nodes whose recorded pipe counts differ from the graph.
            - ``data_gaps``: gaps in the data, all kinds (:attr:`data_gaps`).

        Examples
        --------
        ::

            s = sub.validate().summary()
            print(f"{s['nodes']} nodes, {s['pipes']} pipes, {s['splits']} splits")
        """
        return {
            "nodes": int(self.nodes),
            "pipes": int(self.pipes),
            "splits": len(self.splits),
            "cycles": len(self.cycles),
            "dead_ends": len(self.dead_ends),
            "parallel_pipes": len(self.parallel_pipes),
            "force_mains": len(self.force_mains),
            "stations_without_force_main": len(self.stations_without_force_main),
            "orphan_force_mains": len(self.orphan_force_mains),
            "force_main_ends": len(self.force_main_ends),
            "pipe_count_mismatches": len(self.pipe_count_mismatches),
            "data_gaps": sum(self.data_gaps.values()),
        }


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
        The network the subnetwork belongs to. Reading it raises ``RuntimeError`` once the
        network has changed since the subnetwork was made (e.g. by :meth:`Network.snap_nodes`):
        the subnetwork's nodes and pipes would no longer match the network's.

    graph : networkx.MultiDiGraph
        The induced subgraph of ``network.graph``: edges in the direction of flow, keyed by
        pipe_id (parallel pipes are separate edges), nodes and edges in sorted order, with the
        attributes of ``network.graph``.

    Raises
    ------
    ValueError
        If a node is not in the network.
    """

    def __init__(self, network: Network, nodes: Iterable[NodeId]) -> None:
        nodes = sorted({network.resolve(n) for n in nodes})
        network.check_nodes(nodes)
        self._network = network
        self._version = network.version
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
        self._by_text: dict[str, NodeId] | None = None
        self._pipes_by_text: dict[str, PipeId] | None = None
        # The text of every node of the whole network, for in_network.
        self._network_text: set[str] | None = None

    @classmethod
    def from_node_list(cls, network: Network, path: StrPath) -> Subnetwork:
        """The subnetwork of ``network`` whose nodes are listed in a file written by
        :meth:`to_node_list`. Blank lines and lines starting with '#' are ignored."""
        lines = (line.strip() for line in Path(path).read_text().splitlines())
        # The file holds text; map each line back to the node ID it was written from, so that
        # numeric IDs (e.g. Kitchener's) round-trip.
        known = {str(n): n for n in [*network.aliases, *network.graph]}
        return cls(
            network,
            (known.get(line, line) for line in lines if line and not line.startswith("#")),
        )

    def to_node_list(self, path: StrPath) -> None:
        """Write the nodes to a text file, one per line, in sorted order.

        :meth:`from_node_list` rebuilds the same subnetwork from this file.
        """
        lines = ["# swmmbuilder node list: one node ID per line", *map(str, self.nodes)]
        Path(path).write_text("\n".join(lines) + "\n")

    @property
    def network(self) -> Network:
        if self._network.version != self._version:
            raise RuntimeError(
                "This subnetwork is out of date: its network changed since it was made (e.g. "
                "nodes were merged by snap_nodes). Rebuild it with "
                "network.subnetwork(subnetwork.nodes); merged IDs resolve to the IDs kept."
            )
        return self._network

    @property
    def nodes(self) -> list[NodeId]:
        """The node IDs, sorted."""
        return list(self.graph.nodes)

    def match_node(self, key: Any) -> Any:
        """The node ID that an option key names, tried in this order: ``key`` itself if it is a
        node; the node it resolves to through the network's ``normalize_ids`` and snap aliases
        (see :meth:`Network.resolve`); the node whose ID has the same text as either (so
        ``2151735`` and ``'2151735'`` name the same node); else ``key`` unchanged, for the
        caller to report as not in the subnetwork."""
        if key in self.graph:
            return key
        resolved = self.network.resolve(key)
        if resolved in self.graph:
            return resolved
        if self._by_text is None:
            self._by_text = {str(n): n for n in self.graph}
        by_text = self._by_text
        for candidate in (key, resolved):
            if str(candidate) in by_text:
                return by_text[str(candidate)]
        return key

    def in_network(self, key: Any) -> bool:
        """True if ``key`` names a node of the whole network (see :meth:`match_node`)."""
        network = self.network
        if key in network.graph or network.resolve(key) in network.graph:
            return True
        if self._network_text is None:
            self._network_text = {str(n) for n in network.graph}
        return str(key) in self._network_text

    def match_nodes(self, mapping: Mapping[Any, V], option: str) -> dict[Any, V]:
        """``mapping`` with its keys replaced by :meth:`match_node`.

        Raises
        ------
        ValueError
            If two keys name the same node with different values.
        """
        matched: dict[Any, V] = {}
        clashes: list[str] = []
        for key, value in mapping.items():
            node = self.match_node(key)
            if node in matched and matched[node] != value:
                clashes.append(str(node))
            matched[node] = value
        if clashes:
            raise ValueError(f"{option}: nodes given more than once: {', '.join(clashes)}")
        return matched

    def match_pipe(self, key: Any) -> Any:
        """The pipe ID that an option key names: ``key`` itself if it is a pipe of the
        subnetwork, else the pipe whose ID has the same text (so ``1010`` and ``'1010'`` name the
        same pipe), else ``key`` unchanged."""
        if key in self._ends:
            return key
        if self._pipes_by_text is None:
            self._pipes_by_text = {str(p): p for p in self._ends}
        by_text = self._pipes_by_text
        return by_text.get(str(key), key)

    def match_pipes(
        self, mapping: Mapping[Any, V], option: str, unknown: str = "raise"
    ) -> dict[Any, V]:
        """``mapping`` with its keys replaced by :meth:`match_pipe`.

        Parameters
        ----------
        mapping : dict
            {pipe ID: value}, e.g. an option of :func:`~swmmbuilder.swmm.build_model`.
        option : str
            The option's name, for the error messages.
        unknown : {'raise', 'ignore'}
            Keys that name no pipe of the subnetwork: raise ValueError, or leave them out (for
            options that may list pipes of other subnetworks too).

        Raises
        ------
        ValueError
            If a key names no pipe of the subnetwork (with ``unknown='raise'``), or two keys
            name the same pipe with different values.
        """
        matched: dict[Any, V] = {}
        clashes: list[str] = []
        missing: list[str] = []
        for key, value in mapping.items():
            pipe = self.match_pipe(key)
            if pipe not in self._ends:
                missing.append(str(key))
                continue
            if pipe in matched and matched[pipe] != value:
                clashes.append(str(pipe))
            matched[pipe] = value
        if missing and unknown == "raise":
            raise ValueError(f"{option}: pipes not in the subnetwork: {', '.join(sorted(missing))}")
        if clashes:
            raise ValueError(f"{option}: pipes given more than once: {', '.join(clashes)}")
        return matched

    @property
    def schema(self) -> Schema:
        return self.network.schema

    @property
    def outfalls(self) -> list[NodeId]:
        """The nodes whose node type is ``'outfall'``, sorted: the default outfalls of
        :meth:`validate` and :meth:`to_inp`."""
        return [n for n, t in self.graph.nodes(data="node_type") if t == "outfall"]

    def pipes(self) -> gpd.GeoDataFrame:
        """The rows of ``network.full_data`` of the pipes in the subnetwork, sorted by pipe_id."""
        data = self.network.full_data
        ids = {key for _, _, key in self.graph.edges(keys=True)}
        return data[data["pipe_id"].isin(ids)].sort_values("pipe_id")

    def node_positions(
        self, nodes: Iterable[NodeId] | None = None
    ) -> dict[NodeId, tuple[float, float]]:
        """Unscaled (x, y) position of each node (default: every node), from the pipe geometry.

        A node is placed at the upstream end of a pipe of the network leaving it or, if no pipe
        leaves it, at the downstream end of a pipe entering it. Nodes with no pipe are left out.
        """
        return node_positions(self.network.full_data, self.nodes if nodes is None else nodes)

    def to_geodataframe(
        self, kind: str = "pipes", values: Values | None = None
    ) -> gpd.GeoDataFrame:
        """The subnetwork's pipes or nodes as a GeoDataFrame; see
        :meth:`Network.to_geodataframe`."""
        return to_geodataframe(self.pipes(), self.network.nodes.loc[self.nodes], kind, values)

    def _is_force_main(self, pipe):
        up, down = self._ends[pipe]
        return self.graph.edges[up, down, pipe].get("structure_type") == "force_main"

    def force_main_outlets(self, node: NodeId) -> list[PipeId]:
        """Sorted pipe IDs of the force mains leaving ``node``."""
        return sorted(
            pipe
            for _, _, pipe, structure in self.graph.out_edges(
                node, keys=True, data="structure_type"
            )
            if structure == "force_main"
        )

    def force_main_chain(
        self,
        pipe: PipeId,
        outfalls: Iterable[NodeId] | None = None,
        through_joins: bool = False,
        stations: Iterable[NodeId] | None = None,
    ) -> tuple[list[PipeId], NodeId, str | None]:
        """The force main that starts with ``pipe``, followed downstream.

        A force main is often a chain of several GIS pipes. The chain continues through every node
        whose only outlet pipe is a force main, and stops at a pump station (whose own pump starts
        the next force main) and where it joins another force main.

        Parameters
        ----------
        pipe : pipe ID
            The first pipe of the chain.

        outfalls : iterable of str, optional
            The outfall nodes (default: :attr:`outfalls`).

        through_joins : bool
            Continue through a join with another force main into the shared main downstream
            (default False). The join is still reported as the ``problem``.

        stations : iterable of str, optional
            The pump-station nodes where the chain stops (default: the nodes of type
            'pump_station').

        Returns
        -------
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
        if stations is None:
            stations = (n for n, t in self.graph.nodes(data="node_type") if t == "pump_station")
        stations = set(stations)
        pipes, problem = [pipe], None
        node = self._ends[pipe][1]
        while True:
            if node in stations:
                return pipes, node, problem
            entering = [p for _, _, p in self.graph.in_edges(node, keys=True)]
            if sum(self._is_force_main(p) for p in entering) > 1:
                if not through_joins:
                    return pipes, node, "force_main"
                problem = "force_main"
            outlets = [p for _, _, p in self.graph.out_edges(node, keys=True)]
            if not outlets and node not in outfalls:
                return pipes, node, "dead_end"
            if len(outlets) != 1 or not self._is_force_main(outlets[0]) or outlets[0] in pipes:
                return pipes, node, problem
            pipes.append(outlets[0])
            node = self._ends[outlets[0]][1]

    def find_pump_stations(self, min_rise: float = 0.0) -> PumpStationReport:
        """The pump stations of the subnetwork, and force mains that climb from no station.

        A node of type 'pump_station' (from the node layers) is a station to pump if a force
        main leaves it; one with no force main is reported apart. A force main that starts at
        any other node with no force main entering it (so not inside another station's force
        main) and whose chain (:meth:`force_main_chain`) climbs by more than ``min_rise`` metres
        from its first upstream invert to its last downstream invert is reported in
        ``unlisted``: water does not climb a pressurised main without a pump, so the node is
        likely a lift station missing from the station layer.

        Parameters
        ----------
        min_rise : float
            The rise, in m, above which a force main from a node that is not a station is
            reported (default 0).

        Returns
        -------
        PumpStationReport
        """
        types = dict(self.graph.nodes(data="node_type"))
        stations: list[NodeId] = []
        no_force_main: list[NodeId] = []
        rows: list[dict[str, Any]] = []
        pipes = self.pipes().set_index("pipe_id")
        for node in self.nodes:
            outlets = self.force_main_outlets(node)
            if types.get(node) == "pump_station":
                (stations if outlets else no_force_main).append(node)
                continue
            entering = [p for _, _, p in self.graph.in_edges(node, keys=True)]
            if not outlets or any(self._is_force_main(p) for p in entering):
                continue
            for first in outlets:
                chain, end, _ = self.force_main_chain(first)
                us = gaps._invert(pipes.loc[first, "us_invert"])
                ds = gaps._invert(pipes.loc[chain[-1], "ds_invert"])
                if us is None or ds is None or ds - us <= min_rise:
                    continue
                rows.append(
                    {
                        "node": node,
                        "pipes": chain,
                        "end": end,
                        "us_invert": us,
                        "ds_invert": ds,
                        "rise": round(ds - us, 4),
                    }
                )
        columns = ["node", "pipes", "end", "us_invert", "ds_invert", "rise"]
        unlisted = pd.DataFrame(rows, columns=columns)
        if len(unlisted):
            logger.warning(
                "%d force mains climb from nodes that are not pump stations (likely lift stations "
                "missing from the station layer): %s",
                len(unlisted),
                ", ".join(f"{r.node} ({r.rise:.1f} m)" for r in unlisted.itertuples()),
            )
        return PumpStationReport(stations, no_force_main, unlisted)

    def accumulate(
        self,
        values: Mapping[NodeId, float] | pd.Series | pd.DataFrame,
        split: str | Mapping[PipeId, float] = "equal",
    ) -> UpstreamTotals:
        """Totals of a node value (population, flow, load, ...) upstream of every node and pipe.

        A node's total is its own value plus what its inlet pipes carry; a pipe carries a share of
        the total at its upstream node.

        Parameters
        ----------
        values : dict, pandas.Series, or pandas.DataFrame
            Value of each node, keyed by node ID (missing nodes count 0). A DataFrame is
            accumulated column by column (e.g. population by zone).

        split : 'equal', 'full', or dict
            How a node's total is shared among its outlet pipes:

            - ``'equal'`` (default): equally, so totals are conserved and the outlet's total is
              the sum of all values;
            - ``'full'``: every outlet pipe carries the whole total, i.e. each node counts every
              node upstream of it once, whatever splits lie between (totals are not conserved
              at a split);
            - {pipe ID: fraction}: given fractions, e.g. from simulated flows. Outlet pipes not
              given share what is left at their node equally. The fractions at a node must not
              be negative or add up to more than 1. Where every outlet of a node is given and
              they add up to less than 1, the rest leaves the totals (with a warning).

            Loops are handled exactly: within each strongly connected part of the network the
            totals are solved as a linear system, which needs every loop to have an exit (it
            raises otherwise).

        Returns
        -------
        UpstreamTotals
            ``nodes`` and ``pipes`` totals (e.g. for ``maps.plot_pipes(width=totals.pipes)``),
            and the ``fractions`` used.
        """
        graph = self.graph
        nodes = list(graph.nodes)
        frame = values if isinstance(values, pd.DataFrame) else pd.Series(values, dtype=float)
        frame = frame.rename(index=self.match_node) if len(frame) else frame
        frame = frame.groupby(level=0).sum().reindex(nodes).fillna(0.0).astype(float)
        matrix = frame.to_numpy(dtype=float).reshape(len(nodes), -1)
        edges = list(graph.edges(keys=True))
        if isinstance(split, str):
            if split not in SPLITS:
                raise ValueError(
                    f"split must be one of {', '.join(SPLITS)}, or a dict, not {split!r}."
                )
            if split == "full":
                fractions = {key: 1.0 for _, _, key in edges}
            else:
                fractions = {key: 1.0 / graph.out_degree(u) for u, _, key in edges}
            name = split
        else:
            fractions = self._given_fractions(split, edges)
            name = "given"
        index = {n: i for i, n in enumerate(nodes)}
        if name == "full":
            totals = self._full_totals(nodes, index, matrix)
        else:
            totals = self._shared_totals(nodes, index, matrix, edges, fractions)
        pipe_rows = [totals[index[u]] * fractions[key] for u, _, key in edges]
        pipe_index = pd.Index([key for _, _, key in edges], name="pipe_id")
        node_index = pd.Index(nodes, name="node")
        if isinstance(frame, pd.DataFrame):
            node_totals = pd.DataFrame(totals, index=node_index, columns=frame.columns)
            pipe_totals = pd.DataFrame(
                np.reshape(pipe_rows, (len(edges), matrix.shape[1])),
                index=pipe_index,
                columns=frame.columns,
            )
        else:
            node_totals = pd.Series(totals[:, 0], index=node_index, name=frame.name)
            pipe_totals = pd.Series(
                [row[0] for row in pipe_rows], index=pipe_index, name=frame.name, dtype=float
            )
        return UpstreamTotals(
            node_totals,
            pipe_totals,
            pd.Series([fractions[key] for _, _, key in edges], index=pipe_index, dtype=float),
            name,
        )

    def _given_fractions(self, split, edges):
        """{pipe: fraction} from given fractions, the rest of each node shared equally."""
        given = {p: float(f) for p, f in self.match_pipes(split, "split").items()}
        problems = [f"{p} ({f})" for p, f in given.items() if not (np.isfinite(f) and f >= 0)]
        if problems:
            raise ValueError(f"split: fractions must be finite and >= 0: {', '.join(problems)}")
        fractions, over, under = {}, [], []
        for node in self.graph.nodes:
            outlets = [key for _, _, key in self.graph.out_edges(node, keys=True)]
            known = sum(given[p] for p in outlets if p in given)
            if known > 1 + 1e-9:
                over.append(f"{node} ({known:.6g})")
            rest = [p for p in outlets if p not in given]
            if outlets and not rest and known < 1 - 1e-9:
                under.append(f"{node} ({known:.6g})")
            for p in outlets:
                fractions[p] = given[p] if p in given else max(1.0 - known, 0.0) / len(rest)
        if over:
            raise ValueError(f"split: fractions add up to more than 1 at {', '.join(over)}")
        if under:
            warnings.warn(
                f"split: the fractions of every outlet add up to less than 1 at {len(under)} "
                f"nodes, so the rest leaves the totals: {', '.join(under[:20])}",
                DataWarning,
                stacklevel=3,
            )
        return fractions

    def _shared_totals(self, nodes, index, matrix, edges, fractions):
        """Totals with shared outflow: x = v + A x, solved component by component in flow
        order (a dense solve within each loop)."""
        graph = self.graph
        totals = matrix.copy()
        condensed = nx.condensation(graph)
        for component in nx.topological_sort(condensed):
            members = condensed.nodes[component]["members"]
            if len(members) == 1:
                (node,) = members
                loops = sum(
                    fractions[key] for _, v, key in graph.out_edges(node, keys=True) if v == node
                )
                if loops:
                    if loops >= 1.0 - EXIT_TOLERANCE:
                        raise ValueError(f"accumulate: the loop at {node} has no exit.")
                    totals[index[node]] /= 1.0 - loops
            else:
                order = sorted(members, key=index.get)
                local = {n: i for i, n in enumerate(order)}
                a = np.zeros((len(order), len(order)))
                for u in order:
                    for _, v, key in graph.out_edges(u, keys=True):
                        if v in local:
                            a[local[v], local[u]] += fractions[key]
                trapped = _trapped(order, a)
                if trapped:
                    raise ValueError(
                        f"accumulate: a loop has no exit: {', '.join(map(str, trapped))}"
                    )
                system = np.eye(len(order)) - a
                rows = [index[n] for n in order]
                totals[rows] = np.linalg.solve(system, totals[rows])
            for u in members:
                for _, v, key in graph.out_edges(u, keys=True):
                    if v not in members:
                        totals[index[v]] += fractions[key] * totals[index[u]]
        return totals

    def _full_totals(self, nodes, index, matrix):
        """Totals with every upstream node counted once: the sum over each node's ancestors and
        itself, from ancestor sets (bit masks) built over the loops' condensation."""
        condensed = nx.condensation(self.graph)
        masks: dict[int, int] = {}
        totals = np.zeros_like(matrix)
        width = (len(nodes) + 7) // 8
        for component in nx.topological_sort(condensed):
            mask = 0
            for n in condensed.nodes[component]["members"]:
                mask |= 1 << index[n]
            for parent in condensed.predecessors(component):
                mask |= masks[parent]
            masks[component] = mask
            bits = np.unpackbits(
                np.frombuffer(mask.to_bytes(width, "little"), dtype=np.uint8), bitorder="little"
            )[: len(nodes)].astype(bool)
            total = matrix[bits].sum(axis=0)
            for n in condensed.nodes[component]["members"]:
                totals[index[n]] = total
        return totals

    def validate(self, outfalls: Iterable[NodeId] | None = None) -> TopologyReport:
        """Check the topology before the subnetwork is written to an .inp file.

        Parameters
        ----------
        outfalls : iterable of str, optional
            The outfall nodes, which need no outlet pipe (default: :attr:`outfalls`).

        Returns
        -------
        TopologyReport
            Flow splits, cycles, dead ends (classified by node type) and parallel pipes, plus the
            force mains (with the pump stations that have none, the force mains that start at no
            pump station and the chains that end badly), the nodes whose recorded pipe counts
            differ from the graph and the counts of gaps in the data. Splits and cycles are only
            allowed with dynamic-wave routing; dead ends are never valid.
        """
        outfalls = set(self.outfalls if outfalls is None else outfalls)
        report = TopologyReport(
            nodes=self.graph.number_of_nodes(), pipes=self.graph.number_of_edges()
        )

        for node in self.nodes:
            outlets: dict[NodeId, list[PipeId]] = {}
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

    def data_gaps(self, **fills: Any) -> pd.DataFrame:
        """Every gap in the data of the subnetwork, and how the fill rules would fill it.

        A gap is a pipe with an unknown or missing shape code or a missing size, a pipe end with
        no invert, a node with no invert from any pipe, or a node with no usable rim (see
        :mod:`swmmbuilder.gaps`). With no options nothing is filled, except that a round pipe
        with no diameter whose height equals its width takes that value as its diameter.

        Parameters
        ----------
        **fills
            The fill options of :meth:`to_inp` (``height_correction``, ``size_overrides``,
            ``shape_codes``, ``size_fill``, ``invert_fill``, ``headwater_slope``,
            ``force_main_discharge``, ``rims``, and ``pumps``, whose pipes need no
            cross-section); see
            :func:`swmmbuilder.gaps.fill_gaps`.

        Returns
        -------
        pandas.DataFrame
            One row per gap: ``element`` ('pipe' or 'node'), ``id``, ``gap``, ``detail``,
            ``filled_from`` (the rule that filled it; NaN if it is not filled), ``value``,
            ``source`` and ``order`` (the order of the fills).
        """
        check_options(
            fills, [*swmm.FILL_OPTIONS, "pumps", "height_correction"], "Subnetwork.data_gaps"
        )
        if "pumps" in fills:
            fills["pumps"] = swmm.resolve_pumps(self, fills["pumps"])
        return gaps.fill_gaps(self, **fills).report

    def to_inp(self, path: StrPath, outfalls: Outfalls | None = None, **options: Any) -> SwmmModel:
        """Write the subnetwork to a SWMM .inp file; see :func:`swmmbuilder.swmm.write_inp`.

        Returns
        -------
        SwmmModel
            The model that was written.
        """
        return swmm.write_inp(self, path, outfalls, **options)
