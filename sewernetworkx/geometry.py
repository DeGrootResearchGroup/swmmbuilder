"""Helpers for pipe geometries."""


def line_coords(geometry):
    """(x, y) vertices of a pipe geometry, in order from upstream to downstream.

    The parts of a MultiLineString are joined in order. Any z values are dropped.
    """
    parts = geometry.geoms if geometry.geom_type == "MultiLineString" else [geometry]
    return [(c[0], c[1]) for part in parts for c in part.coords]
