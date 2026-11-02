"""Maps of a network or subnetwork: static plots with matplotlib and interactive Leaflet maps.

The functions here only draw: they take per-pipe, per-node or per-group values computed elsewhere
(concentrations, flows, the share of each upstream area in a pipe's load, ...), as pandas Series
or DataFrames indexed by pipe ID, node ID or group name.

Static plots take a :class:`~swmmbuilder.Network` or :class:`~swmmbuilder.Subnetwork`
(``net``) and draw in a metric CRS (see :func:`plot_crs`), so ribbon widths and polygon buffers
are in metres. Functions drawing on the same axes must use the same ``crs``; the default is the
same for a network and all its subnetworks.

``contextily`` (web-tile basemaps) and ``folium`` (:func:`to_leaflet`) are optional: install them
with ``pip install 'swmmbuilder[maps]'``.
"""

import importlib
import logging
import os
import warnings
from dataclasses import dataclass, field

import geopandas as gpd
import gmplot
import matplotlib as mpl
import numpy as np
import pandas as pd
import shapely
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import LogNorm, Normalize, to_hex
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
from pyproj import CRS, Transformer
from shapely import LineString, MultiPoint

from .geometry import line_coords

logger = logging.getLogger(__name__)

WGS84 = CRS.from_epsg(4326)

# Name of the strand that lumps the categories dropped by plot_ribbons' top_k and min_share.
OTHER = "other"

DEFAULT_HATCHES = ("//", "\\\\", "xx", "..", "++", "oo")


# ---- Coordinates ----------------------------------------------------------------------------


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
    """The CRS the static plots of ``net`` are drawn in.

    ``crs`` if given; otherwise the network's CRS if it is projected, or else the UTM zone of the
    whole network (a network with no CRS is taken to be in longitude/latitude).
    """
    if crs is not None:
        return CRS.from_user_input(crs)
    data = _network(net).full_data
    if data.crs is not None and CRS.from_user_input(data.crs).is_projected:
        return CRS.from_user_input(data.crs)
    return _with_crs(data).estimate_utm_crs()


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


# ---- Values and colours ---------------------------------------------------------------------


def _aligned(values, index, what="values"):
    """``values`` (Series or dict) as floats in the order of ``index``, NaN where missing."""
    series = pd.Series(values, dtype=float) if isinstance(values, dict) else values.astype(float)
    if len(series) and not series.index.isin(index).any():
        logger.warning("None of the IDs of the %s are in the network", what)
    return series.reindex(index).to_numpy()


def _scaled(values, out_range):
    """``values`` mapped linearly from their range onto ``out_range``; NaN maps to the low end
    and constant values to the middle."""
    lo, hi = out_range
    finite = values[np.isfinite(values)]
    if not len(finite):
        return np.full(len(values), float(lo))
    vmin, vmax = finite.min(), finite.max()
    if vmax == vmin:
        out = np.full(len(values), (lo + hi) / 2)
    else:
        out = lo + (values - vmin) / (vmax - vmin) * (hi - lo)
    return np.where(np.isfinite(values), out, lo)


def _norm(norm, vmin, vmax):
    if isinstance(norm, Normalize):
        return norm
    if norm == "linear":
        return Normalize(vmin, vmax)
    if norm == "log":
        return LogNorm(vmin, vmax)
    raise ValueError(f"norm must be 'linear', 'log' or a matplotlib Normalize, not {norm!r}.")


def _colour_by(collection, values, cmap, norm, vmin, vmax, missing_color):
    """Colour ``collection`` by ``values`` through ``cmap``, with NaN (and, for a log norm,
    non-positive values) in ``missing_color``."""
    collection.set_array(np.ma.masked_invalid(values))
    collection.set_cmap(mpl.colormaps.get_cmap(cmap).with_extremes(bad=missing_color))
    collection.set_norm(_norm(norm, vmin, vmax))
    collection.autoscale_None()


def _hex_colours(values, cmap, norm, missing_color):
    """Hex colour of each value through ``cmap`` and ``norm`` (autoscaled to the values), with
    NaN (and, for a log norm, non-positive values) in ``missing_color``."""
    masked = np.ma.masked_invalid(values)
    norm.autoscale_None(masked)
    cmap = mpl.colormaps.get_cmap(cmap).with_extremes(bad=missing_color)
    return [to_hex(c) for c in cmap(np.ma.masked_invalid(norm(masked)))]


def _finish(ax):
    ax.set_aspect("equal")
    ax.autoscale_view()


# ---- Shapes ---------------------------------------------------------------------------------


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


def _line(geometry):
    """A pipe geometry as one LineString (the parts of a MultiLineString joined in order)."""
    return LineString(line_coords(geometry))


def ribbon_strands(line, shares, width):
    """Polygons of the strands of a ribbon along ``line``.

    Strand ``i`` is ``shares[i] / sum(shares) * width`` wide. The strands lie side by side, from
    the right of the line to its left looking along it (downstream for a pipe), and the ribbon
    is centred on the line. Units are those of the line's coordinates.

    Returns
    ----------
    list of shapely.Polygon
        One per share; empty for a zero share.
    """
    shares = np.asarray(shares, dtype=float)
    widths = width * shares / shares.sum()
    edges = np.concatenate([[0.0], np.cumsum(widths)]) - width / 2
    strands = []
    for right, w in zip(edges[:-1], widths, strict=True):
        if w <= 0:
            strands.append(shapely.Polygon())
            continue
        base = line.offset_curve(right, join_style="mitre")
        strands.append(base.buffer(w, cap_style="flat", join_style="mitre", single_sided=True))
    return strands


def _shares_table(shares, top_k, min_share):
    """``shares`` with rows normalised to sum to 1 and the categories below ``min_share`` or
    outside each row's ``top_k`` lumped into an 'other' column."""
    table = shares.astype(float).fillna(0.0)
    if (table.to_numpy() < 0).any():
        raise ValueError("Shares must not be negative.")
    if OTHER in table.columns:
        raise ValueError(f"{OTHER!r} is reserved for the lumped categories.")
    totals = table.sum(axis=1)
    table = table[totals > 0].div(totals[totals > 0], axis=0)
    keep = (table > 0) & (table >= min_share)
    if top_k is not None:
        keep &= table.rank(axis=1, ascending=False, method="first") <= top_k
    other = table.where(~keep, 0.0).sum(axis=1)
    table = table.where(keep, 0.0)
    if (other > 0).any():
        table[OTHER] = other
    return table.loc[:, (table > 0).any()]


def _category_colours(categories, colors, other_color):
    named = [c for c in categories if c != OTHER]
    if colors is None:
        cmap = mpl.colormaps["tab10" if len(named) <= 10 else "tab20"]
        colors = {c: cmap(i % cmap.N) for i, c in enumerate(named)}
    elif not isinstance(colors, dict):
        colors = dict(zip(named, colors, strict=False))
    missing = [c for c in named if c not in colors]
    if missing:
        raise ValueError(f"No colour for categories {missing}.")
    return {**{c: colors[c] for c in named}, OTHER: other_color}


# ---- Static plots ---------------------------------------------------------------------------


def plot_pipes(
    ax,
    net,
    values=None,
    *,
    cmap="viridis",
    norm="linear",
    vmin=None,
    vmax=None,
    color="0.3",
    missing_color="0.85",
    width=None,
    width_range=(0.5, 6.0),
    colorbar=True,
    label=None,
    crs=None,
    zorder=2,
):
    """Draw the pipes of ``net``, coloured by a per-pipe value.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    net : Network or Subnetwork
        Whose pipes are drawn.

    values : pandas.Series or dict, optional
        Value per pipe ID. Pipes with no value (or NaN, or a non-positive value with a log norm)
        are drawn in ``missing_color``. Without values every pipe is drawn in ``color``.

    cmap, norm, vmin, vmax
        Colour map, ``'linear'``, ``'log'`` or a matplotlib Normalize, and its limits (default:
        the range of the values).

    width : float, pandas.Series or dict, optional
        Line width in points: one for every pipe, or a value per pipe (e.g. flow) mapped
        linearly onto ``width_range``. Default 1.5.

    colorbar : bool
        Add a colour bar (labelled ``label``) when ``values`` are given.

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    ----------
    matplotlib.collections.LineCollection
        The pipes, one line each, in the order of ``net.to_geodataframe('pipes')``. Its colour
        bar, if any, is its ``colorbar`` attribute.
    """
    pipes = _geodata(net, "pipes", crs)
    if width is None:
        linewidths = 1.5
    elif np.isscalar(width):
        linewidths = float(width)
    else:
        linewidths = _scaled(_aligned(width, pipes.index, "widths"), width_range)
    segments = [np.asarray(line_coords(g)) for g in pipes.geometry]
    lines = LineCollection(segments, linewidths=linewidths, capstyle="round", zorder=zorder)
    if values is None:
        lines.set_color(color)
    else:
        values = _aligned(values, pipes.index)
        _colour_by(lines, values, cmap, norm, vmin, vmax, missing_color)
    ax.add_collection(lines)
    if values is not None and colorbar:
        ax.figure.colorbar(lines, ax=ax, label=label)
    _finish(ax)
    return lines


def plot_ribbons(
    ax,
    net,
    shares,
    colors=None,
    *,
    width=10.0,
    top_k=None,
    other_color="0.7",
    min_share=0.0,
    order="category",
    legend=True,
    crs=None,
    zorder=2,
):
    """Draw each pipe as a ribbon of parallel strands, one per category, each as wide as its
    share of the pipe (a mass-flow map).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    net : Network or Subnetwork
        The pipes. Pipes with no shares are not drawn.

    shares : pandas.DataFrame
        Pipes (index: pipe ID) x categories (e.g. upstream areas). Each row is normalised to sum
        to 1, so loads can be given directly; NaN counts as 0 and rows summing to 0 are skipped.

    colors : dict or list, optional
        Colour per category (a list in the order of the columns). Default: tab10, or tab20 for
        more than 10 categories.

    width : float, pandas.Series or dict
        Total ribbon width in metres (of the plot CRS): one for every pipe, or per pipe ID.

    top_k : int, optional
        Keep only each pipe's ``top_k`` largest categories and lump the rest into 'other'.

    other_color
        Colour of the 'other' strand.

    min_share : float
        Lump categories with a smaller share of the pipe into 'other'.

    order : {'category', 'share'}
        Order of the strands from the right of the pipe to its left (looking downstream): the
        order of the columns (with 'other' last), or by decreasing share.

    legend : bool
        Add a legend of the categories drawn.

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    ----------
    matplotlib.collections.PatchCollection
        One patch per strand drawn.
    """
    if order not in ("category", "share"):
        raise ValueError(f"order must be 'category' or 'share', not {order!r}.")
    pipes = _geodata(net, "pipes", crs)
    table = _shares_table(shares, top_k, min_share)
    table = table.loc[table.index.isin(pipes.index)]
    widths = (
        pd.Series(float(width), index=table.index)
        if np.isscalar(width)
        else pd.Series(_aligned(width, table.index, "widths"), index=table.index)
    )
    palette = _category_colours(table.columns, colors, other_color)

    patches, facecolors = [], []
    for pipe_id, row in table.iterrows():
        row = row[row > 0]
        if order == "share":
            row = row.sort_values(ascending=False, kind="stable")
        if not np.isfinite(widths[pipe_id]) or widths[pipe_id] <= 0:
            continue
        strands = ribbon_strands(_line(pipes.geometry[pipe_id]), row.to_numpy(), widths[pipe_id])
        for category, strand in zip(row.index, strands, strict=True):
            path = _path(strand)
            if path is not None:
                patches.append(PathPatch(path))
                facecolors.append(palette[category])

    ribbons = PatchCollection(
        patches, facecolors=facecolors, edgecolors="none", linewidths=0, zorder=zorder
    )
    ax.add_collection(ribbons)
    if legend:
        handles = [Patch(color=palette[c], label=str(c)) for c in table.columns]
        ax.legend(handles=handles)
    _finish(ax)
    return ribbons


def node_group_polygons(
    net, groups, method="voronoi", *, clip_to=None, buffer=10.0, ratio=0.3, crs=None
):
    """Polygons covering groups of nodes (e.g. drainage zones).

    Parameters
    ----------
    net : Network or Subnetwork
        The nodes.

    groups : dict
        Group name -> node IDs.

    method : {'voronoi', 'concave_hull'}
        'voronoi': the union of the Voronoi tiles of each group's nodes (computed over the
        grouped nodes only), clipped to ``clip_to``; a node may then be in only one group.
        'concave_hull': the concave hull of each group's nodes (see ``ratio``) grown by
        ``buffer``, clipped to ``clip_to`` if given; hulls may overlap.

    clip_to : shapely geometry, GeoSeries or GeoDataFrame, optional
        Area to clip to. A shapely geometry must be in the plot CRS; a GeoSeries or GeoDataFrame
        is reprojected and its geometries merged. Default for 'voronoi': the convex hull of all
        the nodes of ``net``, grown by ``buffer``.

    buffer : float
        Distance in metres by which hulls are grown, so that groups of one or two nodes, and
        nodes on the edge, still get an area.

    ratio : float
        Concaveness of 'concave_hull' (0: most concave, 1: convex hull); see
        :func:`shapely.concave_hull`.

    crs : optional
        CRS of the polygons; default :func:`plot_crs`.

    Returns
    ----------
    geopandas.GeoDataFrame
        A (Multi)Polygon per group, in the order of ``groups``, indexed by group name.

    Raises
    ----------
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
        owner = {}
        for group, ps in points.items():
            for p in ps:
                if owner.setdefault((p.x, p.y), group) != group:
                    raise ValueError(
                        f"Node at {(p.x, p.y)} is in groups {owner[(p.x, p.y)]!r} and {group!r}."
                    )
        seeds = [shapely.Point(xy) for xy in owner]
        groups_of_seeds = list(owner.values())
        if len(seeds) == 1:
            cells = [clip_to]
        else:
            cells = list(shapely.voronoi_polygons(MultiPoint(seeds), extend_to=clip_to).geoms)
        tree = shapely.STRtree(seeds)
        cell_index, seed_index = tree.query(cells, predicate="contains")
        tiles = {g: [] for g in points}
        for c, s in zip(cell_index, seed_index, strict=True):
            tiles[groups_of_seeds[s]].append(cells[c])
        polygons = {g: shapely.union_all(t).intersection(clip_to) for g, t in tiles.items()}
    else:
        raise ValueError(f"method must be 'voronoi' or 'concave_hull', not {method!r}.")

    return gpd.GeoDataFrame(
        {"geometry": list(polygons.values())},
        index=pd.Index(list(polygons), name="group"),
        crs=crs,
    )


def plot_groups(
    ax,
    polygons,
    values=None,
    *,
    cmap="viridis",
    norm="linear",
    vmin=None,
    vmax=None,
    color="0.85",
    missing_color="0.85",
    edgecolor="0.3",
    linewidth=0.8,
    alpha=0.6,
    hatch_groups=None,
    hatches=DEFAULT_HATCHES,
    colorbar=True,
    label=None,
    zorder=1,
):
    """Draw group polygons, shaded by a per-group value.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    polygons : geopandas.GeoDataFrame
        Polygons indexed by group name, in the CRS of the axes, e.g. from
        :func:`node_group_polygons`.

    values : pandas.Series or dict, optional
        Value per group, coloured through ``cmap``/``norm`` (groups with no value in
        ``missing_color``). Without values every polygon is filled with ``color``.

    hatch_groups : list of list, optional
        Sets of groups that belong together: the outline of each set is drawn hatched, with
        the hatches of ``hatches`` in turn.

    colorbar : bool
        Add a colour bar (labelled ``label``) when ``values`` are given.

    Returns
    ----------
    matplotlib.collections.PatchCollection
        The polygons, one patch per group that has an area. The hatched sets are separate
        collections, labelled 'set 1', 'set 2', ... for a legend.
    """
    paths = [(g, _path(geometry)) for g, geometry in polygons.geometry.items()]
    paths = [(g, p) for g, p in paths if p is not None]
    shaded = PatchCollection(
        [PathPatch(p) for _, p in paths],
        edgecolors=edgecolor,
        linewidths=linewidth,
        alpha=alpha,
        zorder=zorder,
    )
    if values is None:
        shaded.set_facecolor(color)
    else:
        values = _aligned(values, pd.Index([g for g, _ in paths]), "group values")
        _colour_by(shaded, values, cmap, norm, vmin, vmax, missing_color)
    ax.add_collection(shaded)

    for i, members in enumerate(hatch_groups or []):
        path = _path(shapely.union_all(list(polygons.geometry.loc[list(members)])))
        if path is None:
            continue
        ax.add_collection(
            PatchCollection(
                [PathPatch(path)],
                facecolors="none",
                edgecolors=edgecolor,
                linewidths=1.5 * linewidth,
                hatch=hatches[i % len(hatches)],
                zorder=zorder + 0.1,
                label=f"set {i + 1}",
            )
        )

    if values is not None and colorbar:
        ax.figure.colorbar(shaded, ax=ax, label=label)
    _finish(ax)
    return shaded


def plot_markers(
    ax,
    net,
    node_ids,
    sizes=None,
    labels=True,
    *,
    color="crimson",
    marker="o",
    size=60.0,
    size_range=(20.0, 200.0),
    edgecolor="white",
    fontsize=8,
    crs=None,
    zorder=4,
):
    """Mark nodes, e.g. sampler or site locations.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    net : Network or Subnetwork
        The nodes.

    node_ids : iterable
        Nodes to mark.

    sizes : pandas.Series or dict, optional
        Value per node, mapped linearly onto marker areas ``size_range`` (points²). Default:
        every marker of area ``size``.

    labels : bool, dict or pandas.Series
        Label each marker with its node ID (True), with the given text per node, or not at all.

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    ----------
    matplotlib.collections.PathCollection
        The markers, in the order of ``node_ids``.

    Raises
    ----------
    ValueError
        If a node is not in ``net`` or has no position.
    """
    node_ids = list(node_ids)
    _check_nodes(net, node_ids)
    positions = _geodata(net, "nodes", crs).geometry
    unplaced = [str(n) for n in node_ids if n not in positions.index]
    if unplaced:
        raise ValueError(f"Nodes with no position: {', '.join(unplaced)}")
    points = positions.loc[node_ids]
    x, y = points.x.to_numpy(), points.y.to_numpy()
    s = size if sizes is None else _scaled(_aligned(sizes, points.index, "sizes"), size_range)
    markers = ax.scatter(x, y, s=s, c=color, marker=marker, edgecolors=edgecolor, zorder=zorder)
    if labels is not False and labels is not None:
        texts = {n: n for n in node_ids} if labels is True else dict(labels)
        for node, xi, yi in zip(node_ids, x, y, strict=True):
            if node in texts:
                ax.annotate(
                    str(texts[node]),
                    (xi, yi),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=fontsize,
                    zorder=zorder,
                )
    _finish(ax)
    return markers


def _optional(module):
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(
            f"This needs {module}, an optional dependency: pip install 'swmmbuilder[maps]'"
        ) from e


def add_basemap(
    ax,
    net,
    source=None,
    *,
    layer=None,
    attribution=None,
    color="0.75",
    linewidth=0.6,
    crs=None,
    zorder=0,
    **kwargs,
):
    """Draw a basemap beneath what is on ``ax``: web tiles, or a vector layer for offline and
    publication use.

    Call it after plotting the network, so that the tiles cover the axes' extent.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes drawn by the plots of ``net``.

    net : Network or Subnetwork
        Gives the CRS of the axes (with ``crs``; see :func:`plot_crs`).

    source : optional
        Tile provider for ``contextily`` (an ``xyzservices`` provider or a URL); default
        CartoDB Positron. The provider's attribution is written on the plot unless
        ``attribution`` is given (False for none). Needs the ``[maps]`` extra.

    layer : geopandas.GeoDataFrame, GeoSeries or path, optional
        A vector layer (roads, rivers, ...) drawn instead of tiles, in ``color`` and
        ``linewidth``. It must have a CRS.

    **kwargs
        Passed to :func:`contextily.add_basemap` or to the layer's ``plot``.
    """
    crs = plot_crs(net, crs)
    if layer is not None:
        if not isinstance(layer, gpd.GeoDataFrame | gpd.GeoSeries):
            layer = gpd.read_file(layer)
        if layer.crs is None:
            raise ValueError("The basemap layer has no CRS.")
        layer.to_crs(crs).plot(ax=ax, color=color, linewidth=linewidth, zorder=zorder, **kwargs)
        _finish(ax)
        return
    cx = _optional("contextily")
    if source is None:
        source = cx.providers.CartoDB.Positron
    cx.add_basemap(
        ax, crs=crs.to_string(), source=source, attribution=attribution, zorder=zorder, **kwargs
    )


# ---- Interactive maps and export ------------------------------------------------------------


@dataclass
class PipeLayer:
    """Pipes on a :func:`to_leaflet` map, coloured by ``values`` (per pipe ID) like
    :func:`plot_pipes`, or all in ``color``. Tooltips show the pipe ID and value."""

    values: object = None
    name: str = "Pipes"
    cmap: object = "viridis"
    norm: object = "linear"
    vmin: float | None = None
    vmax: float | None = None
    color: str = "#3366cc"
    missing_color: str = "#bbbbbb"
    weight: float = 3.0
    show: bool = True


@dataclass
class GroupLayer:
    """Group polygons on a :func:`to_leaflet` map (e.g. from :func:`node_group_polygons`),
    shaded by ``values`` (per group) or all in ``color``. Tooltips show the group and value."""

    polygons: gpd.GeoDataFrame
    values: object = None
    name: str = "Groups"
    cmap: object = "viridis"
    norm: object = "linear"
    vmin: float | None = None
    vmax: float | None = None
    color: str = "#999999"
    missing_color: str = "#dddddd"
    opacity: float = 0.5
    show: bool = True


@dataclass
class MarkerLayer:
    """Node markers on a :func:`to_leaflet` map. Tooltips show the node ID, or ``labels[node]``."""

    node_ids: list
    name: str = "Markers"
    labels: dict = field(default_factory=dict)
    color: str = "crimson"
    radius: float = 6.0
    show: bool = True


def _colour_legend(m, cmap, norm, caption):
    branca = importlib.import_module("branca.colormap")
    stops = np.linspace(0, 1, 11)
    index = [float(norm.inverse(s)) for s in stops]
    colours = [to_hex(c) for c in mpl.colormaps.get_cmap(cmap)(stops)]
    branca.LinearColormap(
        colours, index=index, vmin=index[0], vmax=index[-1], caption=caption
    ).add_to(m)


def _value_colours(values, layer, default_colour):
    """(hex colours, autoscaled norm) of ``values``, or (``default_colour`` each, None)."""
    if layer.values is None:
        return [default_colour] * len(values), None
    norm = _norm(layer.norm, layer.vmin, layer.vmax)
    return _hex_colours(values, layer.cmap, norm, layer.missing_color), norm


def _add_pipes(m, folium, net, layer):
    pipes = _with_crs(net.to_geodataframe("pipes")).to_crs(WGS84)
    data = gpd.GeoDataFrame({"geometry": pipes.geometry.to_numpy()}, crs=WGS84)
    data.insert(0, "pipe_id", pipes.index.astype(str))
    fields, aliases = ["pipe_id"], ["Pipe"]
    values = np.full(len(pipes), np.nan)
    if layer.values is not None:
        values = _aligned(layer.values, pipes.index)
        data.insert(1, "value", values)
        fields, aliases = [*fields, "value"], [*aliases, layer.name]
    colours, norm = _value_colours(values, layer, layer.color)
    data["_colour"] = colours
    group = folium.FeatureGroup(name=layer.name, show=layer.show)
    folium.GeoJson(
        data.to_json(na="null"),
        style_function=lambda f: {"color": f["properties"]["_colour"], "weight": layer.weight},
        tooltip=folium.GeoJsonTooltip(fields=fields, aliases=aliases),
    ).add_to(group)
    group.add_to(m)
    if norm is not None:
        _colour_legend(m, layer.cmap, norm, layer.name)


def _add_groups(m, folium, layer):
    polygons = _with_crs(layer.polygons).to_crs(WGS84)
    data = gpd.GeoDataFrame({"geometry": polygons.geometry.to_numpy()}, crs=WGS84)
    data.insert(0, "group", polygons.index.astype(str))
    fields, aliases = ["group"], ["Group"]
    values = np.full(len(polygons), np.nan)
    if layer.values is not None:
        values = _aligned(layer.values, polygons.index, "group values")
        data.insert(1, "value", values)
        fields, aliases = [*fields, "value"], [*aliases, layer.name]
    colours, norm = _value_colours(values, layer, layer.color)
    data["_colour"] = colours
    group = folium.FeatureGroup(name=layer.name, show=layer.show)
    folium.GeoJson(
        data.to_json(na="null"),
        style_function=lambda f: {
            "fillColor": f["properties"]["_colour"],
            "color": "#444444",
            "weight": 1,
            "fillOpacity": layer.opacity,
        },
        tooltip=folium.GeoJsonTooltip(fields=fields, aliases=aliases),
    ).add_to(group)
    group.add_to(m)
    if norm is not None:
        _colour_legend(m, layer.cmap, norm, layer.name)


def _add_markers(m, folium, net, layer):
    _check_nodes(net, layer.node_ids)
    nodes = _with_crs(net.to_geodataframe("nodes")).to_crs(WGS84).geometry
    group = folium.FeatureGroup(name=layer.name, show=layer.show)
    for node in layer.node_ids:
        if node not in nodes.index:
            continue
        folium.CircleMarker(
            location=[nodes[node].y, nodes[node].x],
            radius=layer.radius,
            color=layer.color,
            fill=True,
            fill_opacity=0.9,
            tooltip=str(layer.labels.get(node, node)),
        ).add_to(group)
    group.add_to(m)


def to_leaflet(net, path, layers=None, *, tiles="OpenStreetMap"):
    """Write an interactive Leaflet map of ``net`` to an HTML file, with ``folium``.

    The data are embedded in the file; Leaflet and the basemap tiles are loaded from the web
    when it is opened. Needs the ``[maps]`` extra.

    Parameters
    ----------
    net : Network or Subnetwork
        The network.

    path : str or Path
        HTML file to write.

    layers : list of PipeLayer, GroupLayer and MarkerLayer, optional
        Layers, drawn in order and listed in a layer switcher. Default: the pipes.

    tiles : str
        Basemap tiles, as for :class:`folium.Map` (None for none).

    Returns
    ----------
    folium.Map
        The map that was written.
    """
    folium = _optional("folium")
    layers = [PipeLayer()] if layers is None else layers
    m = folium.Map(tiles=tiles, control_scale=True)
    xmin, ymin, xmax, ymax = _with_crs(net.to_geodataframe("pipes")).to_crs(WGS84).total_bounds
    m.fit_bounds([[ymin, xmin], [ymax, xmax]])
    for layer in layers:
        if isinstance(layer, PipeLayer):
            _add_pipes(m, folium, net, layer)
        elif isinstance(layer, GroupLayer):
            _add_groups(m, folium, layer)
        elif isinstance(layer, MarkerLayer):
            _add_markers(m, folium, net, layer)
        else:
            raise TypeError(f"Not a map layer: {layer!r}")
    folium.LayerControl(collapsed=False).add_to(m)
    m.save(str(path))
    return m


def to_geojson(net, path, kind="pipes", values=None):
    """Write the pipes or nodes of ``net`` (with ``values``; see
    :meth:`~swmmbuilder.Network.to_geodataframe`) to a GeoJSON file in WGS 84."""
    data = _with_crs(net.to_geodataframe(kind, values)).to_crs(WGS84)
    data.to_file(path, driver="GeoJSON")


def generate_gmap(subnetwork, file_name, origin=None):
    """Write an html file showing the subnetwork overlaid on Google Maps.

    .. deprecated::
        Use :func:`to_leaflet`, which needs no API key and colours pipes by value.

    The node positions are reprojected to WGS 84 (see :func:`lonlat_positions`), so networks in
    a projected CRS are drawn in the right place.

    Parameters
    ----------
    subnetwork : Subnetwork
        The subnetwork.

    file_name : str
        Name of the html file.

    origin : (float, float), optional
        (latitude, longitude) of the centre of the map. Defaults to the first node.
    """
    warnings.warn(
        "generate_gmap is deprecated; use swmmbuilder.maps.to_leaflet instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    if os.path.splitext(file_name)[1] != ".html":
        raise ValueError("Wrong file extension used. The proper extension is: .html.")

    positions = lonlat_positions(subnetwork)
    if origin is None:
        if not positions:
            raise ValueError("The subnetwork has no node positions to centre the map on.")
        lon, lat = next(iter(positions.values()))
        origin = (lat, lon)
    gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], zoom=18)

    for node, (lon, lat) in positions.items():
        gmap.marker(lat, lon, "cornflowerblue", label=str(node))

    for up, down in sorted(set(subnetwork.graph.edges())):
        (lon0, lat0), (lon1, lat1) = positions[up], positions[down]
        gmap.plot([lat0, lat1], [lon0, lon1], "red", edge_width=2)

    gmap.draw(file_name)
