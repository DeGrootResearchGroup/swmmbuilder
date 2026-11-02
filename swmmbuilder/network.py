"""A city's sewer network: its pipes, loaded through a :class:`~swmmbuilder.Schema`, and the
flow graph between their nodes."""

import logging
from pathlib import Path

import geopandas as gpd
import networkx as nx
import pandas as pd

from .schema import LONDON_ON, as_schema
from .subnetwork import Subnetwork

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
        Keep pipes flagged abandoned (default False: they are dropped). Pipes with a node ID the
        schema treats as no node (for London, 'Abandoned') are dropped from the graph either way.

    schema : Schema, str or Path
        How the GIS columns map to the canonical columns: a :class:`Schema` (default
        :data:`LONDON_ON`; see :mod:`swmmbuilder.schema` for the other presets) or the path to a
        schema TOML file.

    Attributes
    ----------
    full_data : geopandas.GeoDataFrame
        The pipes, with the canonical columns of :mod:`swmmbuilder.schema`.

    schema : Schema
        The schema the data were loaded with.

    graph : networkx.MultiDiGraph
        One edge per pipe, in the direction of flow (from_node -> to_node), keyed by pipe_id so
        that parallel pipes are kept. Each edge has the attributes ``pipe`` (pipe_id) and
        ``length``. Nodes and edges are inserted in sorted order, so traversals do not depend
        on the order of the file.

    removed : dict
        What was dropped on loading: ``{'flow_type': {flow_type: count}, 'abandoned': count}``.
        Pipes with no flow type are counted under ``'None'``.
    """

    def __init__(
        self, source, flow_types=DEFAULT_FLOW_TYPES, keep_abandoned=False, schema=LONDON_ON
    ):
        self.schema = schema = as_schema(schema)
        data = source if isinstance(source, gpd.GeoDataFrame) else gpd.read_file(Path(source))
        self.full_data = schema.normalize(data)
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

    def subnetwork(self, nodes):
        """The :class:`Subnetwork` made of ``nodes`` and every pipe between two of them.

        Raises
        ----------
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

        outlet : str or list of str, optional
            Outlet(s) downstream of ``node``. A shortest path (by length) is added to every listed
            outlet that ``node`` drains to.

        Raises
        ----------
        ValueError
            If a node is not in the network, or ``node`` drains to none of the outlets.
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

        if outlet is not None:
            outlets = [outlet] if isinstance(outlet, str) else list(outlet)
            self._check_nodes(outlets)
            reachable = [o for o in sorted(outlets) if nx.has_path(self.graph, node, o)]
            if not reachable:
                raise ValueError(f"Node {node} does not drain to any of the outlets {outlets}.")
            for o in reachable:
                selected.update(nx.dijkstra_path(self.graph, node, o, weight=self._pipe_length))

        return self.subnetwork(selected)

    def _check_nodes(self, nodes):
        missing = sorted(n for n in nodes if n not in self.graph)
        if missing:
            raise ValueError(f"Nodes not in the network: {', '.join(map(str, missing))}")

    @staticmethod
    def _pipe_length(u, v, pipes):
        """Dijkstra weight: the shortest of the (possibly parallel) pipes from u to v."""
        return min(attrs["length"] for attrs in pipes.values())
