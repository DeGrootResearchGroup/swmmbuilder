"""Helpers for pipe geometries."""

import pandas as pd
import shapely


def _parts(geometry):
    """The lines of a pipe geometry: a MultiLineString's parts are first merged where they
    touch end to start, so a pipe split into connected pieces is one line."""
    if geometry is None or geometry.is_empty:
        return []
    if geometry.geom_type != "MultiLineString":
        return [geometry]
    merged = shapely.line_merge(geometry, directed=True)
    if merged.is_empty:  # e.g. degenerate zero-length parts: keep the parts as they are
        return list(geometry.geoms)
    return list(merged.geoms) if hasattr(merged, "geoms") else [merged]


def is_disjoint(geometry):
    """True if a pipe geometry is a MultiLineString whose parts do not join into one line."""
    return len(_parts(geometry)) > 1


def line_coords(geometry):
    """(x, y) vertices of a pipe geometry, in order from upstream to downstream; [] if the
    geometry is missing or empty.

    The parts of a MultiLineString are merged where they join; parts that do not join (see
    :func:`is_disjoint`) are concatenated in order. Any z values are dropped.
    """
    return [(c[0], c[1]) for part in _parts(geometry) for c in part.coords]


def node_positions(data, nodes):
    """(x, y) position of each of ``nodes`` from the pipes in ``data``.

    A node is placed at the upstream end of a pipe leaving it or, if no pipe leaves it, at the
    downstream end of a pipe entering it. Pipes with no geometry are skipped, and nodes with no
    pipe geometry are left out.
    """
    nodes = set(nodes)
    upstream_ends, downstream_ends = {}, {}
    for upstream_id, downstream_id, geometry in zip(
        data["from_node"], data["to_node"], data["geometry"], strict=True
    ):
        coords = line_coords(geometry)
        if not coords:
            continue
        if upstream_id in nodes and upstream_id not in upstream_ends:
            upstream_ends[upstream_id] = coords[0]
        if downstream_id in nodes and downstream_id not in downstream_ends:
            downstream_ends[downstream_id] = coords[-1]
    return {**downstream_ends, **upstream_ends}


def to_geodataframe(pipes, nodes, kind="pipes", values=None):
    """``pipes`` indexed by pipe_id (kind 'pipes') or the ``nodes`` table (kind 'nodes'),
    optionally joined with ``values``; see :meth:`swmmbuilder.Network.to_geodataframe`."""
    if kind == "pipes":
        data = pipes.set_index("pipe_id")
    elif kind == "nodes":
        data = nodes
    else:
        raise ValueError(f"kind must be 'pipes' or 'nodes', not {kind!r}.")
    if values is None:
        return data
    if isinstance(values, dict):
        values = pd.Series(values)
    if isinstance(values, pd.Series):
        values = values.to_frame("value" if values.name is None else values.name)
    return data.join(values)
