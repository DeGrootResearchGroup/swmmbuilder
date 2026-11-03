"""Coordinates: the plot CRS (always projected in metres), WGS 84 positions, and the
extent and footprint of a network."""

import math

import numpy as np
import shapely
from pyproj import CRS, Transformer

WGS84 = CRS.from_epsg(4326)


def lonlat_positions(subnetwork):
    """{node: (longitude, latitude)} of the subnetwork's nodes, in WGS 84.

    Positions in another CRS (e.g. a projected one such as :data:`TORONTO_ON`'s EPSG:2952) are
    reprojected. A network with no CRS is taken to be in longitude/latitude already.
    """
    positions = subnetwork.node_positions()
    crs = subnetwork.network.full_data.crs
    if crs is None or CRS.from_user_input(crs) == WGS84:
        return positions
    transformer = Transformer.from_crs(crs, WGS84, always_xy=True)
    return {node: transformer.transform(x, y) for node, (x, y) in positions.items()}


def plot_crs(net, crs=None):
    """The CRS the static plots of ``net`` are drawn in, always projected in metres, so that the
    scale bar, buffers, footprints and ribbon widths are in metres.

    ``crs`` if given; otherwise the network's CRS if it is projected in metres, or else (a
    geographic CRS, or one in feet) the UTM zone of the whole network. A network with no CRS is
    taken to be in longitude/latitude.

    Raises
    ------
    ValueError
        If ``crs`` is not projected in metres, or the network has no CRS and its coordinates
        cannot be longitude/latitude.
    """
    if crs is not None:
        crs = CRS.from_user_input(crs)
        if not _in_metres(crs):
            raise ValueError(
                f"The plot CRS must be projected in metres, not {crs.name} "
                f"({'geographic' if crs.is_geographic else crs.axis_info[0].unit_name})."
            )
        return crs
    data = _network(net).full_data
    if data.crs is not None and _in_metres(CRS.from_user_input(data.crs)):
        return CRS.from_user_input(data.crs)
    if data.crs is None:
        xmin, ymin, xmax, ymax = data.total_bounds
        if not (-180 <= xmin <= xmax <= 180 and -90 <= ymin <= ymax <= 90):
            raise ValueError(
                "The network has no CRS, and its coordinates are not longitude/latitude; set "
                "the schema's crs (or the GIS file's) to the CRS of the data."
            )
    return _with_crs(data).estimate_utm_crs()


def _in_metres(crs):
    """True if ``crs`` is projected with axes in metres."""
    return crs.is_projected and math.isclose(crs.axis_info[0].unit_conversion_factor, 1.0)


def _network(net):
    """The Network of a Network or Subnetwork."""
    return getattr(net, "network", net)


def _with_crs(data):
    return data.set_crs(WGS84) if data.crs is None else data


def _geodata(net, kind, crs=None, values=None):
    """``net.to_geodataframe(kind, values)`` in ``plot_crs(net, crs)``, without empty geometries."""
    data = _with_crs(net.to_geodataframe(kind, values)).to_crs(plot_crs(net, crs))
    return data[data.geometry.notna() & ~data.geometry.is_empty]


def _check_nodes(net, nodes):
    missing = sorted({str(n) for n in nodes if n not in net.graph})
    if missing:
        raise ValueError(f"Nodes not in the network: {', '.join(missing)}")


FOOTPRINT_PIPES = 2.0


FOOTPRINT_EXTENT = 0.05


def network_extent(net, crs=None):
    """The larger side of the bounding box of the pipes of ``net``, in metres of the plot CRS."""
    xmin, ymin, xmax, ymax = _geodata(net, "pipes", crs).total_bounds
    return float(max(xmax - xmin, ymax - ymin))


def footprint_distance(net, crs=None):
    """The default :func:`network_footprint` distance of ``net``: twice the median pipe length,
    at most 5 % of :func:`network_extent`, in metres."""
    pipes = _geodata(net, "pipes", crs)
    median = float(np.median(pipes.geometry.length)) if len(pipes) else 0.0
    return min(FOOTPRINT_PIPES * median, FOOTPRINT_EXTENT * network_extent(net, crs))


def network_footprint(net, distance=None, crs=None):
    """The area served by the pipes of ``net``: within about ``distance`` metres of a pipe
    (default :func:`footprint_distance`), with narrow gaps between pipes filled in, as a shapely
    geometry in the plot CRS."""
    if distance is None:
        distance = footprint_distance(net, crs)
    pipes = _geodata(net, "pipes", crs).geometry.to_numpy()
    smooth = distance / 3
    return shapely.union_all(shapely.buffer(pipes, distance + smooth)).buffer(-smooth)
