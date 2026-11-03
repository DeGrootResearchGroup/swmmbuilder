"""Polygons covering groups of nodes (e.g. drainage zones)."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from matplotlib.path import Path
from shapely import MultiPoint

from .crs import _check_nodes, _geodata, _with_crs, network_footprint, plot_crs

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping
    from typing import Any

    from .._types import NodeId
    from ..network import Network
    from ..subnetwork import Subnetwork

logger = logging.getLogger(__name__)


def _polygons(geometry):
    """The Polygons of a (Multi)Polygon or GeometryCollection."""
    if geometry.geom_type == "Polygon":
        return [] if geometry.is_empty else [geometry]
    return [p for part in getattr(geometry, "geoms", []) for p in _polygons(part)]


def _path(geometry):
    """A matplotlib Path of the polygons of ``geometry`` (holes included), or None if it has
    none."""
    rings = [r for p in _polygons(geometry) for r in [p.exterior, *p.interiors]]
    if not rings:
        return None
    return Path.make_compound_path(*[Path(np.asarray(r.coords)[:, :2], closed=True) for r in rings])


def node_group_polygons(
    net: Network | Subnetwork,
    groups: Mapping[Any, Iterable[NodeId]],
    method: str = "voronoi",
    *,
    clip_to: Any = None,
    footprint: float | bool | None = None,
    buffer: float = 10.0,
    ratio: float = 0.3,
    crs: Any = None,
) -> gpd.GeoDataFrame:
    """Polygons covering groups of nodes (e.g. drainage zones).

    Parameters
    ----------
    net : Network or Subnetwork
        The nodes.

    groups : dict
        Group name -> node IDs.

    method : {'voronoi', 'concave_hull'}
        'voronoi': the union of the Voronoi tiles of each group's nodes (computed over the
        grouped nodes only); a node may then be in only one group. 'concave_hull': the concave
        hull of each group's nodes (see ``ratio``) grown by ``buffer``; hulls may overlap.

    clip_to : shapely geometry, GeoSeries or GeoDataFrame, optional
        Area to clip to. A shapely geometry must be in the plot CRS; a GeoSeries or GeoDataFrame
        is reprojected and its geometries merged. Default: the network footprint (see
        ``footprint``), so that the polygons cover the serviced area.

    footprint : float, None or False
        Default clip area: the :func:`network_footprint` of ``net`` at this distance in metres,
        or by default at :func:`footprint_distance` (twice the median pipe length, at most 5 %
        of the extent). False: the convex hull of all the nodes of ``net`` grown by ``buffer``
        for 'voronoi', and no clipping for 'concave_hull'.

    buffer : float
        Distance in metres by which hulls are grown, so that groups of one or two nodes, and
        nodes on the edge, still get an area.

    ratio : float
        Concaveness of 'concave_hull' (0: most concave, 1: convex hull); see
        :func:`shapely.concave_hull`.

    crs : optional
        CRS of the polygons; default :func:`plot_crs`.

    Returns
    -------
    geopandas.GeoDataFrame
        A (Multi)Polygon per group, in the order of ``groups``, indexed by group name.

    Raises
    ------
    ValueError
        If a node is not in ``net``, a group has no node with a position, or, for 'voronoi', a
        node position is in two groups.
    """
    crs = plot_crs(net, crs)
    positions = _geodata(net, "nodes", crs).geometry
    _check_nodes(net, [n for nodes in groups.values() for n in nodes])
    points = {}
    for group, nodes in groups.items():
        points[group] = [positions[n] for n in nodes if n in positions.index]
        if not points[group]:
            raise ValueError(f"Group {group!r} has no node with a position.")

    if clip_to is not None and not isinstance(clip_to, shapely.Geometry):
        clip_to = _with_crs(clip_to).to_crs(crs).geometry.union_all()
    if clip_to is None and footprint is not False:
        clip_to = network_footprint(net, footprint, crs)

    if method == "concave_hull":
        polygons = {
            g: shapely.concave_hull(MultiPoint(p), ratio=ratio).buffer(buffer)
            for g, p in points.items()
        }
        if clip_to is not None:
            polygons = {g: p.intersection(clip_to) for g, p in polygons.items()}
    elif method == "voronoi":
        if clip_to is None:
            clip_to = MultiPoint(list(positions)).convex_hull.buffer(buffer)
        owner: dict[tuple[float, float], Any] = {}
        for group, ps in points.items():
            for p in ps:
                if owner.setdefault((p.x, p.y), group) != group:
                    raise ValueError(
                        f"Node at {(p.x, p.y)} is in groups {owner[(p.x, p.y)]!r} and {group!r}."
                    )
        seeds = [shapely.Point(xy) for xy in owner]
        groups_of_seeds = list(owner.values())
        tiles: dict[Any, list] = {g: [] for g in points}
        if len(seeds) == 1:
            # The whole area, even if the one node is outside it.
            tiles[groups_of_seeds[0]].append(clip_to)
        else:
            cells = list(shapely.voronoi_polygons(MultiPoint(seeds), extend_to=clip_to).geoms)
            tree = shapely.STRtree(seeds)
            cell_index, seed_index = tree.query(cells, predicate="contains")
            for c, i in zip(cell_index, seed_index, strict=True):
                tiles[groups_of_seeds[i]].append(cells[c])
        polygons = {g: shapely.union_all(t).intersection(clip_to) for g, t in tiles.items()}
    else:
        raise ValueError(f"method must be 'voronoi' or 'concave_hull', not {method!r}.")

    empty = [str(g) for g, polygon in polygons.items() if polygon.is_empty]
    if empty:
        logger.warning(
            "Groups with an empty polygon (their nodes are outside clip_to or the footprint): %s",
            ", ".join(empty),
        )
    return gpd.GeoDataFrame(
        {"geometry": list(polygons.values())},
        index=pd.Index(list(polygons), name="group"),
        crs=crs,
    )
