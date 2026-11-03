"""Static plots of pipes, group polygons and node markers."""

import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import shapely
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch

from ..geometry import line_coords
from .colours import _aligned, _colour_by, _scaled, category_colors
from .crs import _check_nodes, _geodata
from .polygons import _path
from .style import DrawnLayer, _finish, add_legends

# Pipes drawn as context under other layers (plot_pipes(context=True)): colour and width factor.
CONTEXT_COLOR = "0.6"


CONTEXT_WIDTH = 0.6


DEFAULT_HATCHES = ("//", "\\\\", "xx", "..", "++", "oo")


def plot_pipes(
    ax,
    net,
    values=None,
    *,
    cmap="viridis",
    norm="linear",
    vmin=None,
    vmax=None,
    color=None,
    missing_color="0.75",
    width=None,
    width_range=(0.4, 4.0),
    width_scale="sqrt",
    context=False,
    legend=True,
    label=None,
    style=True,
    crs=None,
    zorder=None,
):
    """Draw the pipes of ``net``, coloured by a per-pipe value, or as light context under other
    layers.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    net : Network or Subnetwork
        Whose pipes are drawn.

    values : pandas.Series or dict, optional
        Value per pipe ID. Pipes with no value (or NaN, or a non-positive value with a log norm)
        are drawn in ``missing_color``. Without values every pipe is drawn in ``color``
        (default ``0.35`` grey).

    cmap, norm, vmin, vmax
        Colour map, ``'linear'``, ``'log'`` or a matplotlib Normalize, and its limits (default:
        the range of the values).

    width : float, pandas.Series or dict, optional
        Line width in points (constant on screen when zooming): one for every pipe, or a value
        per pipe (e.g. flow) mapped onto ``width_range`` by ``width_scale``. Default 1.

    width_scale : {'sqrt', 'linear', 'log'}
        How per-pipe widths are scaled: proportional to the square root of the value (the
        default, readable over several orders of magnitude) or to the value, from 0 at the low
        end of ``width_range``; or the range of the logarithms onto ``width_range``.

    context : bool
        Draw the pipes as context beneath other layers (e.g. zones or ribbons): lighter (default
        ``color`` ``0.6`` grey), 0.6 times as wide, and between group fills (zorder 1) and
        ribbons (zorder 2). ``values`` cannot be given.

    legend : bool
        Add a colour bar (labelled ``label``) next to the axes when ``values`` are given.

    style : bool
        Style the axes as a map (:func:`map_axes`).

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    -------
    DrawnLayer
        ``artist``: a LineCollection with one line per pipe, thinnest first so that thick pipes
        are drawn on top; ``ids``: the pipe IDs in that order.
    """
    if context and values is not None:
        raise ValueError("Context pipes are drawn in one colour; values cannot be given.")
    if color is None:
        color = CONTEXT_COLOR if context else "0.35"
    if zorder is None:
        zorder = 1.5 if context else 2
    pipes = _geodata(net, "pipes", crs)
    if width is None:
        linewidths = np.full(len(pipes), 1.0)
    elif np.isscalar(width):
        linewidths = np.full(len(pipes), float(width))
    else:
        linewidths = _scaled(_aligned(width, pipes.index, "widths"), width_range, width_scale)
    if context:
        linewidths = CONTEXT_WIDTH * linewidths
    order = np.argsort(linewidths, kind="stable")
    pipes, linewidths = pipes.iloc[order], linewidths[order]

    segments = [np.asarray(line_coords(g)) for g in pipes.geometry]
    lines = LineCollection(segments, linewidths=linewidths, capstyle="round", zorder=zorder)
    coloured = False
    if values is None:
        lines.set_color(color)
    else:
        aligned = _aligned(values, pipes.index)
        coloured = _colour_by(lines, aligned, cmap, norm, vmin, vmax, missing_color)
    ax.add_collection(lines)

    layer = DrawnLayer(lines, list(pipes.index), mappable=lines if coloured else None, title=label)
    if legend:
        add_legends(ax, layer)
    _finish(ax, style)
    return layer


def plot_groups(
    ax,
    polygons,
    values=None,
    colors=None,
    *,
    cmap="viridis",
    norm="linear",
    vmin=None,
    vmax=None,
    color=None,
    missing_color="0.85",
    edgecolor="white",
    linewidth=1.2,
    alpha=0.3,
    hatch_groups=None,
    hatches=DEFAULT_HATCHES,
    hatch_color="0.55",
    legend=True,
    label=None,
    style=True,
    zorder=1,
):
    """Draw group polygons as light fills: by category (the default), by a per-group value, or
    in one colour.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    polygons : geopandas.GeoDataFrame
        Polygons indexed by group name, in the CRS of the axes, e.g. from
        :func:`node_group_polygons`.

    values : pandas.Series or dict, optional
        Value per group, coloured through ``cmap``/``norm`` with a colour bar (groups with no
        value in ``missing_color``).

    colors : dict or list, optional
        Without ``values``: colour per group (default :func:`category_colors` of the groups),
        with a legend of the groups.

    color : optional
        Without ``values``: one colour for every group, and no legend.

    edgecolor, linewidth, alpha
        Outlines and opacity of the fills; light by default, so that pipes drawn on top stand
        out.

    hatch_groups : list of list, or dict, optional
        Sets of groups that belong together: the outline of each set is hatched in
        ``hatch_color``, with the hatches of ``hatches`` in turn. A dict names the sets for the
        legend (default 'set 1', 'set 2', ...).

    legend : bool
        Add the colour bar (labelled ``label``) or the legend (titled ``label``) next to the
        axes, with the hatched sets.

    style : bool
        Style the axes as a map (:func:`map_axes`).

    Returns
    -------
    DrawnLayer
        ``artist``: a PatchCollection with one patch per group that has an area; ``ids``: those
        groups.
    """
    paths = [(g, _path(geometry)) for g, geometry in polygons.geometry.items()]
    paths = [(g, p) for g, p in paths if p is not None]
    ids = [g for g, _ in paths]
    shaded = PatchCollection(
        [PathPatch(p) for _, p in paths],
        edgecolors=edgecolor,
        linewidths=linewidth,
        alpha=alpha,
        zorder=zorder,
    )
    handles, mappable = [], None
    if values is not None:
        values = _aligned(values, pd.Index(ids), "group values")
        if _colour_by(shaded, values, cmap, norm, vmin, vmax, missing_color):
            mappable = shaded
    elif color is not None:
        shaded.set_facecolor(color)
    else:
        # Over every group, drawn or not, so that a list of colours keeps its order.
        palette = category_colors(list(polygons.index), colors)
        shaded.set_facecolor([palette[g] for g in ids])
        handles = [
            Patch(facecolor=palette[g], alpha=min(1.0, 2 * alpha), label=str(g)) for g in ids
        ]
    ax.add_collection(shaded)

    sets = hatch_groups or []
    if not isinstance(sets, dict):
        sets = {f"set {i + 1}": members for i, members in enumerate(sets)}
    for i, (name, members) in enumerate(sets.items()):
        path = _path(shapely.union_all(list(polygons.geometry.loc[list(members)])))
        if path is None:
            continue
        hatch = hatches[i % len(hatches)]
        ax.add_collection(
            PatchCollection(
                [PathPatch(path)],
                facecolors="none",
                edgecolors=hatch_color,
                linewidths=linewidth,
                hatch=hatch,
                zorder=zorder + 0.1,
                label=name,
            )
        )
        handles.append(Patch(facecolor="none", edgecolor=hatch_color, hatch=hatch, label=name))

    layer = DrawnLayer(shaded, ids, mappable=mappable, handles=handles, title=label)
    if legend:
        add_legends(ax, layer)
    _finish(ax, style)
    return layer


def plot_markers(
    ax,
    net,
    node_ids,
    sizes=None,
    labels=False,
    *,
    name=None,
    color="#c0392b",
    marker="o",
    size=60.0,
    size_range=(20.0, 200.0),
    edgecolor="white",
    fontsize=8,
    legend=True,
    style=True,
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
        Value per node; marker areas (points²) are proportional to it, up to ``size_range[1]``
        and from ``size_range[0]``. Default: every marker of area ``size``.

    labels : bool, dict or pandas.Series
        No labels (the default), the node IDs (True), or the given text per node, e.g.
        ``{"WM6002": "Medway P.S."}``. Labels have a white halo.

    name : str, optional
        Legend entry for the markers.

    legend : bool
        Add the legend entry ``name`` next to the axes, if given.

    style : bool
        Style the axes as a map (:func:`map_axes`).

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    -------
    DrawnLayer
        ``artist``: the markers (a PathCollection), in the order of ``node_ids``.

    Raises
    ------
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
    markers = ax.scatter(
        x, y, s=s, c=color, marker=marker, edgecolors=edgecolor, linewidths=1.2, zorder=zorder
    )
    if labels is not False and labels is not None:
        texts = {n: n for n in node_ids} if labels is True else dict(labels)
        halo = [pe.withStroke(linewidth=3, foreground="white")]
        for node, xi, yi in zip(node_ids, x, y, strict=True):
            if node in texts:
                ax.annotate(
                    str(texts[node]),
                    (xi, yi),
                    xytext=(6, 4),
                    textcoords="offset points",
                    fontsize=fontsize,
                    color="0.1",
                    path_effects=halo,
                    zorder=zorder + 0.1,
                )
    handles = []
    if name is not None:
        handles = [
            Line2D(
                [],
                [],
                linestyle="none",
                marker=marker,
                markerfacecolor=color,
                markeredgecolor=edgecolor,
                markersize=np.sqrt(size),
                label=name,
            )
        ]
    layer = DrawnLayer(markers, node_ids, handles=handles)
    if legend:
        add_legends(ax, layer)
    _finish(ax, style)
    return layer
