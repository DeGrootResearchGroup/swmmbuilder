"""Helpers for pipe geometries."""

import geopandas as gpd
import pandas as pd
from shapely import Point


def line_coords(geometry):
    """(x, y) vertices of a pipe geometry, in order from upstream to downstream.

    The parts of a MultiLineString are joined in order. Any z values are dropped.
    """
    parts = geometry.geoms if geometry.geom_type == "MultiLineString" else [geometry]
    return [(c[0], c[1]) for part in parts for c in part.coords]


def node_positions(data, nodes):
    """(x, y) position of each of ``nodes`` from the pipes in ``data``.

    A node is placed at the upstream end of a pipe leaving it or, if no pipe leaves it, at the
    downstream end of a pipe entering it. Nodes with no pipe are left out.
    """
    nodes = set(nodes)
    upstream_ends, downstream_ends = {}, {}
    for upstream_id, downstream_id, geometry in zip(
        data["from_node"], data["to_node"], data["geometry"], strict=True
    ):
        if upstream_id in nodes and upstream_id not in upstream_ends:
            upstream_ends[upstream_id] = line_coords(geometry)[0]
        if downstream_id in nodes and downstream_id not in downstream_ends:
            downstream_ends[downstream_id] = line_coords(geometry)[-1]
    return {**downstream_ends, **upstream_ends}


def to_geodataframe(pipes, positions, kind="pipes", values=None):
    """GeoDataFrame of ``pipes`` (kind 'pipes') or of the nodes at ``positions`` (kind 'nodes'),
    optionally joined with ``values``; see :meth:`swmmbuilder.Network.to_geodataframe`."""
    if kind == "pipes":
        data = pipes.set_index("pipe_id")
    elif kind == "nodes":
        ids = sorted(positions)
        data = gpd.GeoDataFrame(
            {"geometry": [Point(positions[n]) for n in ids]},
            index=pd.Index(ids, name="node_id"),
            crs=pipes.crs,
        )
    else:
        raise ValueError(f"kind must be 'pipes' or 'nodes', not {kind!r}.")
    if values is None:
        return data
    if isinstance(values, dict):
        values = pd.Series(values)
    if isinstance(values, pd.Series):
        values = values.to_frame("value" if values.name is None else values.name)
    return data.join(values)
