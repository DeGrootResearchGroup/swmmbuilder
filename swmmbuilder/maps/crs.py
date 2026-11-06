"""Coordinates: the plot CRS (always projected in metres), WGS 84 positions, and the
extent and footprint of a network."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import shapely
from pyproj import CRS, Transformer

if TYPE_CHECKING:
    from typing import Any

    from .._types import NodeId
    from ..network import Network
    from ..subnetwork import Subnetwork

WGS84 = CRS.from_epsg(4326)


def lonlat_positions(subnetwork: Subnetwork) -> dict[NodeId, tuple[float, float]]:
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


def plot_crs(net: Network | Subnetwork, crs: Any = None) -> CRS:
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


def _pipes(net, crs=None, pipes=None):
    """The pipes of ``net`` in the plot CRS, keeping only those of ``pipes`` (see
    :func:`network_footprint`)."""
    data = _geodata(net, "pipes", crs)
    if pipes is None:
        return data
    if isinstance(pipes, str):
        if pipes != "gravity":
            raise ValueError(
                f"pipes must be 'gravity', a list of pipe IDs, or a function, not {pipes!r}."
            )
        keep = data["structure_type"] != "force_main"
    elif callable(pipes):
        keep = np.asarray(pipes(data), dtype=bool)
    else:
        wanted = {str(p) for p in pipes}
        keep = data.index.map(lambda p: str(p) in wanted).to_numpy(dtype=bool)
    return data[keep]


def network_extent(net: Network | Subnetwork, crs: Any = None, pipes: Any = None) -> float:
    """The larger side of the bounding box of the pipes of ``net`` (or of ``pipes``; see
    :func:`network_footprint`), in metres of the plot CRS."""
    xmin, ymin, xmax, ymax = _pipes(net, crs, pipes).total_bounds
    return float(max(xmax - xmin, ymax - ymin))


def footprint_distance(net: Network | Subnetwork, crs: Any = None, pipes: Any = None) -> float:
    """The default :func:`network_footprint` distance of ``net``: twice the median pipe length,
    at most 5 % of :func:`network_extent`, in metres, over the pipes of ``pipes`` (see
    :func:`network_footprint`)."""
    data = _pipes(net, crs, pipes)
    median = float(np.median(data.geometry.length)) if len(data) else 0.0
    return min(FOOTPRINT_PIPES * median, FOOTPRINT_EXTENT * network_extent(net, crs, pipes))


def network_footprint(
    net: Network | Subnetwork,
    distance: float | None = None,
    crs: Any = None,
    pipes: Any = None,
) -> shapely.Geometry:
    """The area served by the pipes of ``net``: within about ``distance`` metres of a pipe
    (default :func:`footprint_distance` of the same pipes), with narrow gaps between pipes
    filled in, as a shapely geometry in the plot CRS.

    Parameters
    ----------
    net : Network or Subnetwork
        Whose pipes are used.

    distance : float, optional
        Distance in metres (default :func:`footprint_distance`).

    crs : optional
        The plot CRS (see :func:`plot_crs`).

    pipes : 'gravity', list of pipe IDs, or callable, optional
        The pipes to use: every pipe of ``net`` (default None), its gravity pipes
        (``'gravity'``: every pipe but the force mains, e.g. for a service area, since a force
        main crossing open land serves no one), the pipes with these IDs, or the pipes for
        which a function of the pipe rows (a GeoDataFrame indexed by pipe ID) returns True.

    Returns
    -------
    shapely.Geometry
    """
    if distance is None:
        distance = footprint_distance(net, crs, pipes)
    lines = _pipes(net, crs, pipes).geometry.to_numpy()
    smooth = distance / 3
    return shapely.union_all(shapely.buffer(lines, distance + smooth)).buffer(-smooth)
