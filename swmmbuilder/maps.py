"""Maps of a network or subnetwork: static plots with matplotlib and interactive Leaflet maps.

The functions here only draw: they take per-pipe, per-node or per-group values computed elsewhere
(concentrations, flows, the share of each upstream area in a pipe's load, ...), as pandas Series
or DataFrames indexed by pipe ID, node ID or group name.

Static plots take a :class:`~swmmbuilder.Network` or :class:`~swmmbuilder.Subnetwork`
(``net``) and draw in a metric CRS (see :func:`plot_crs`). Functions drawing on the same axes must
use the same ``crs``; the default is the same for a network and all its subnetworks.

Units: pipe line widths and marker sizes are in points, so they stay the same on screen when you
zoom; ribbon widths, footprints and buffers are in metres of the plot CRS, so they scale with the
map.

Each static plot styles the axes as a map (:func:`map_axes`) and returns a :class:`DrawnLayer`.
With ``legend=True`` (the default) its colour bar or legend is stacked to the right of the axes
(see :func:`add_legends`), below those of earlier plots on the same axes.

``contextily`` (web-tile basemaps) and ``folium`` (:func:`to_leaflet`) are optional: install them
with ``pip install 'swmmbuilder[maps]'``.
"""

import importlib
import logging
import os
import warnings
import weakref
from dataclasses import dataclass, field

import geopandas as gpd
import gmplot
import matplotlib as mpl
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import shapely
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.colors import LogNorm, Normalize, to_hex
from matplotlib.legend import Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, NullFormatter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from pyproj import CRS, Transformer
from shapely import LineString, MultiPoint

from .geometry import line_coords

logger = logging.getLogger(__name__)

WGS84 = CRS.from_epsg(4326)

# Name of the strand that lumps the categories dropped by plot_ribbons' top_k and min_share.
OTHER = "other"

# Categories drawn in NEUTRAL_COLOR by category_colors.
NEUTRAL_CATEGORIES = frozenset({"other", "rest"})

# Tableau 10 (2016), designed to stay distinguishable for colour-blind readers, without its grey,
# which is kept for the neutral categories.
TABLEAU10 = (
    "#4e79a7",
    "#f28e2b",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
)
NEUTRAL_COLOR = "#bab0ac"

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


def network_footprint(net, distance=100.0, crs=None):
    """The area served by the pipes of ``net``: within about ``distance`` metres of a pipe, with
    narrow gaps between pipes filled in, as a shapely geometry in the plot CRS."""
    pipes = _geodata(net, "pipes", crs).geometry.to_numpy()
    smooth = distance / 3
    return shapely.union_all(shapely.buffer(pipes, distance + smooth)).buffer(-smooth)


# ---- Values and colours ---------------------------------------------------------------------


def _aligned(values, index, what="values"):
    """``values`` (Series or dict) as floats in the order of ``index``, NaN where missing."""
    series = pd.Series(values, dtype=float) if isinstance(values, dict) else values.astype(float)
    if len(series) and not series.index.isin(index).any():
        logger.warning("None of the IDs of the %s are in the network", what)
    return series.reindex(index).to_numpy()


SCALES = ("sqrt", "linear", "log")


def _scaled(values, out_range, scale="linear"):
    """``values`` mapped onto ``out_range``.

    'linear' and 'sqrt' are proportional to the value (or its square root) from 0 at the low end
    to the largest value at the high end; 'log' maps the range of the logarithms. NaN, and values
    that the scale cannot take (negative, or non-positive for 'log'), map to the low end.
    """
    if scale not in SCALES:
        raise ValueError(f"scale must be one of {SCALES}, not {scale!r}.")
    lo, hi = out_range
    values = np.asarray(values, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        if scale == "log":
            t = np.where(values > 0, np.log10(values), np.nan)
        else:
            t = np.where(values >= 0, values, np.nan)
            t = np.sqrt(t) if scale == "sqrt" else t
    finite = t[np.isfinite(t)]
    if not len(finite):
        return np.full(len(t), float(lo))
    if scale == "log":
        tmin, tmax = finite.min(), finite.max()
    else:
        tmin, tmax = 0.0, finite.max()
    if tmax == tmin:
        out = np.full(len(t), float(hi))
    else:
        out = lo + (t - tmin) / (tmax - tmin) * (hi - lo)
    return np.where(np.isfinite(t), out, lo)


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


def category_colors(categories, colors=None, other=NEUTRAL_COLOR):
    """{category: colour} for categorical maps, e.g. zones and the ribbons of their loads.

    ``colors`` (a dict, or a list in the order of the categories) where given; otherwise Tableau 10
    in order, repeating after 9 categories. The neutral categories 'other' and 'rest' get
    ``other`` unless ``colors`` names them. Pass the result as ``colors`` to both
    :func:`plot_groups` and :func:`plot_ribbons` so that a group and its strands match.
    """
    categories = list(dict.fromkeys(categories))
    named = [c for c in categories if c not in NEUTRAL_CATEGORIES]
    if colors is None:
        if len(named) > len(TABLEAU10):
            logger.warning(
                "%d categories for %d colours: colours repeat", len(named), len(TABLEAU10)
            )
        colors = {c: TABLEAU10[i % len(TABLEAU10)] for i, c in enumerate(named)}
    elif not isinstance(colors, dict):
        colors = dict(zip(named, colors, strict=False))
    missing = [c for c in named if c not in colors]
    if missing:
        raise ValueError(f"No colour for categories {missing}.")
    return {c: colors.get(c, other) if c in NEUTRAL_CATEGORIES else colors[c] for c in categories}


# ---- Map axes and keys ----------------------------------------------------------------------


@dataclass
class DrawnLayer:
    """What a static plot drew, and its key for :func:`add_legends`.

    Attributes
    ----------
    artist : matplotlib.artist.Artist
        The main collection drawn.

    ids : list
        Pipe IDs, node IDs or group names of the elements of ``artist``, in its order.

    mappable : matplotlib.cm.ScalarMappable or None
        Source of a colour bar, for plots coloured by value.

    handles : list
        Legend handles, for categories, hatched sets or named markers.

    title : str or None
        Colour bar label or legend title.

    colorbar : matplotlib.colorbar.Colorbar or None
        The colour bar, once :func:`add_legends` has drawn it.

    legend : matplotlib.legend.Legend or None
        The legend, once :func:`add_legends` has drawn it.
    """

    artist: object
    ids: list = field(default_factory=list)
    mappable: object = None
    handles: list = field(default_factory=list)
    title: str | None = None
    colorbar: object = None
    legend: object = None


# Per-axes state: the scale bar and north arrow (replaced on restyling) and the top of the free
# space in the column of keys to the right of the axes, in axes fraction.
_AXES_STATE = weakref.WeakKeyDictionary()


def _state(ax):
    return _AXES_STATE.setdefault(ax, {"key_top": 1.0, "decorations": []})


def _nice_length(length):
    """The largest 1, 2 or 5 x 10^n not above ``length``."""
    power = 10 ** np.floor(np.log10(length))
    return max(m * power for m in (1, 2, 5) if m * power <= length)


def map_axes(ax, scalebar=True, north_arrow=False, units="m"):
    """Style ``ax`` as a map: equal aspect, no ticks or frame, a scale bar about a fifth of the
    width (rounded to 1, 2 or 5 x 10^n, in km from 1000 m) and optionally a north arrow.

    Every static plot calls it unless ``style=False``; calling it again (e.g. after changing the
    limits) replaces the scale bar and arrow. ``units`` are those of the plot CRS.
    """
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set_axis_off()
    state = _state(ax)
    for artist in state["decorations"]:
        artist.remove()
    state["decorations"] = []

    if scalebar:
        (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
        if x1 > x0:
            length = _nice_length((x1 - x0) / 5)
            text = (
                f"{length / 1000:g} km"
                if units == "m" and length >= 1000
                else f"{length:g} {units}"
            )
            bar = AnchoredSizeBar(
                ax.transData,
                length,
                text,
                "lower left",
                pad=0.4,
                sep=3,
                frameon=False,
                size_vertical=abs(y1 - y0) / 200,
                fontproperties={"size": 8},
            )
            ax.add_artist(bar)
            state["decorations"].append(bar)

    if north_arrow:
        arrow = ax.annotate(
            "",
            xy=(0.96, 0.93),
            xytext=(0.96, 0.83),
            xycoords="axes fraction",
            arrowprops={"arrowstyle": "-|>", "color": "0.2", "linewidth": 1.2},
        )
        letter = ax.text(
            0.96, 0.94, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=9
        )
        state["decorations"] += [arrow, letter]


def _finish(ax, style):
    if style:
        map_axes(ax)
    else:
        ax.set_aspect("equal")
        ax.autoscale_view()


def _renderer(fig):
    try:
        return fig.canvas.get_renderer()
    except AttributeError:
        return None


def add_legends(ax, *layers, x=1.02, colorbar_height=0.3):
    """Draw the colour bars and legends of ``layers`` in a column to the right of ``ax``, from the
    top, below the keys already added to ``ax``.

    Colour bars are compact (``colorbar_height`` of the axes) with small type, and only major
    ticks on a log norm; legends are frameless, titled and left-aligned. The colour bar and
    legend are stored on each layer.
    """
    state = _state(ax)
    fig = ax.figure
    for layer in layers:
        if layer.mappable is not None:
            top = state["key_top"]
            cax = ax.inset_axes([x, top - colorbar_height, 0.03, colorbar_height])
            layer.colorbar = fig.colorbar(layer.mappable, cax=cax)
            layer.colorbar.set_label(layer.title or "", size=8)
            layer.colorbar.ax.tick_params(labelsize=7)
            norm = layer.mappable.norm
            if isinstance(norm, LogNorm) and norm.vmin and norm.vmax:
                # Label the decades only, or 1, 2 and 5 x 10^n over less than two decades.
                decades = np.log10(norm.vmax / norm.vmin)
                subs = (1.0,) if decades >= 2 else (1.0, 2.0, 5.0)
                layer.colorbar.ax.yaxis.set_major_locator(LogLocator(subs=subs))
                layer.colorbar.ax.yaxis.set_major_formatter(
                    LogFormatterSciNotation(labelOnlyBase=False, minor_thresholds=(np.inf, np.inf))
                )
                layer.colorbar.ax.yaxis.set_minor_formatter(NullFormatter())
            state["key_top"] = top - colorbar_height - 0.08
        if layer.handles:
            top = state["key_top"]
            legend = Legend(
                ax,
                layer.handles,
                [h.get_label() for h in layer.handles],
                loc="upper left",
                bbox_to_anchor=(x - 0.01, top),
                bbox_transform=ax.transAxes,
                frameon=False,
                fontsize=8,
                title=None if layer.mappable is not None else layer.title,
                title_fontsize=8,
                alignment="left",
            )
            ax.add_artist(legend)
            legend.set_clip_on(False)  # outside the axes, and in a tight bounding box
            layer.legend = legend
            renderer = _renderer(fig)
            if renderer is not None:
                height = (
                    legend.get_window_extent(renderer).height
                    / ax.get_window_extent(renderer).height
                )
            else:
                height = 0.06 * (len(layer.handles) + 1)
            state["key_top"] = top - height - 0.04


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


def _joined(geometries):
    """One LineString through pipe geometries in flow order, without repeating shared ends."""
    coords = []
    for geometry in geometries:
        part = line_coords(geometry)
        if coords and np.allclose(coords[-1], part[0]):
            part = part[1:]
        coords += part
    return LineString(coords)


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


def _reaches(graph, pipe_ids):
    """``pipe_ids`` as runs in flow order between junctions: a run continues through a node with
    exactly one pipe in and one out (in ``graph``) when both pipes are in ``pipe_ids``. Pipes not
    in ``graph`` are runs of their own."""
    wanted = set(pipe_ids)
    ends = {key: (up, down) for up, down, key in graph.edges(keys=True) if key in wanted}

    def through(node):
        return graph.in_degree(node) == 1 and graph.out_degree(node) == 1

    def pipe_in(node):
        return next(iter(graph.in_edges(node, keys=True)))[2]

    def pipe_out(node):
        return next(iter(graph.out_edges(node, keys=True)))[2]

    runs, seen = [], set()
    for key in pipe_ids:
        if key not in ends or key in seen:
            continue
        up, _ = ends[key]
        if through(up) and pipe_in(up) in ends:
            continue  # not the start of a run
        run = [key]
        while through(ends[run[-1]][1]) and pipe_out(ends[run[-1]][1]) in ends:
            following = pipe_out(ends[run[-1]][1])
            if following in seen or following in run:
                break
            run.append(following)
        seen.update(run)
        runs.append(run)
    # Pipes on cycles of through nodes, and pipes not in the graph.
    runs += [[key] for key in pipe_ids if key not in seen]
    return runs


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
    color="0.35",
    missing_color="0.75",
    width=None,
    width_range=(0.4, 4.0),
    width_scale="sqrt",
    legend=True,
    label=None,
    style=True,
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
        Line width in points (constant on screen when zooming): one for every pipe, or a value
        per pipe (e.g. flow) mapped onto ``width_range`` by ``width_scale``. Default 1.

    width_scale : {'sqrt', 'linear', 'log'}
        How per-pipe widths are scaled: proportional to the square root of the value (the
        default, readable over several orders of magnitude) or to the value, from 0 at the low
        end of ``width_range``; or the range of the logarithms onto ``width_range``.

    legend : bool
        Add a colour bar (labelled ``label``) next to the axes when ``values`` are given.

    style : bool
        Style the axes as a map (:func:`map_axes`).

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    ----------
    DrawnLayer
        ``artist``: a LineCollection with one line per pipe, thinnest first so that thick pipes
        are drawn on top; ``ids``: the pipe IDs in that order.
    """
    pipes = _geodata(net, "pipes", crs)
    if width is None:
        linewidths = np.full(len(pipes), 1.0)
    elif np.isscalar(width):
        linewidths = np.full(len(pipes), float(width))
    else:
        linewidths = _scaled(_aligned(width, pipes.index, "widths"), width_range, width_scale)
    order = np.argsort(linewidths, kind="stable")
    pipes, linewidths = pipes.iloc[order], linewidths[order]

    segments = [np.asarray(line_coords(g)) for g in pipes.geometry]
    lines = LineCollection(segments, linewidths=linewidths, capstyle="round", zorder=zorder)
    if values is None:
        lines.set_color(color)
    else:
        _colour_by(lines, _aligned(values, pipes.index), cmap, norm, vmin, vmax, missing_color)
    ax.add_collection(lines)

    layer = DrawnLayer(
        lines, list(pipes.index), mappable=None if values is None else lines, title=label
    )
    if legend:
        add_legends(ax, layer)
    _finish(ax, style)
    return layer


def plot_ribbons(
    ax,
    net,
    shares,
    colors=None,
    *,
    width="total",
    max_width=30.0,
    min_width=0.0,
    top_k=None,
    min_share=0.0,
    other_color=NEUTRAL_COLOR,
    order="category",
    context=True,
    legend=True,
    title=None,
    style=True,
    crs=None,
    zorder=2,
):
    """Draw pipes as ribbons of parallel strands, one per category, each as wide as its share of
    the pipe's load (a mass-flow map).

    Consecutive pipes between junctions with the same shares and width are drawn as one ribbon,
    so there are no notches where they meet.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to draw on.

    net : Network or Subnetwork
        The pipes. Pipes with no shares are not drawn as ribbons.

    shares : pandas.DataFrame
        Pipes (index: pipe ID) x categories (e.g. upstream areas): the load of each category in
        each pipe, or its share. NaN counts as 0 and rows summing to 0 are skipped.

    colors : dict or list, optional
        Colour per category; default :func:`category_colors` of the columns.

    width : 'total', float, pandas.Series or dict
        Ribbon width in metres (it scales with the map). 'total' (the default): proportional to
        the row total, ``max_width`` for the largest, so that the ribbons show the load. A
        number: the same for every pipe; a Series or dict: per pipe ID.

    max_width : float
        Width in metres of the pipe with the largest total, for ``width='total'``.

    min_width : float
        Pipes narrower than this (in metres) are not drawn as ribbons; with ``context`` they
        are still drawn as thin grey lines.

    top_k : int, optional
        Keep only each pipe's ``top_k`` largest categories and lump the rest into 'other'.

    min_share : float
        Lump categories with a smaller share of the pipe into 'other'.

    other_color
        Colour of the 'other' strand.

    order : {'category', 'share'}
        Order of the strands from the right of the pipe to its left (looking downstream): the
        order of the columns (with 'other' last), or by decreasing share.

    context : bool
        Draw every pipe as a thin grey line beneath the ribbons.

    legend : bool
        Add a legend (titled ``title``) of the categories drawn.

    style : bool
        Style the axes as a map (:func:`map_axes`).

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    ----------
    DrawnLayer
        ``artist``: a PatchCollection with one patch per strand; ``ids``: the pipe IDs drawn as
        ribbons; ``handles``: one per category drawn.
    """
    if order not in ("category", "share"):
        raise ValueError(f"order must be 'category' or 'share', not {order!r}.")
    pipes = _geodata(net, "pipes", crs)
    raw = shares.astype(float).fillna(0.0)
    raw = raw.loc[raw.index.isin(pipes.index)]
    if isinstance(width, str):
        if width != "total":
            raise ValueError(f"width must be 'total', a number or per-pipe values, not {width!r}.")
        totals = raw.sum(axis=1)
        widths = totals / totals.max() * max_width if len(totals) else totals
    elif np.isscalar(width):
        widths = pd.Series(float(width), index=raw.index)
    else:
        widths = pd.Series(_aligned(width, raw.index, "widths"), index=raw.index)
    widths = widths[np.isfinite(widths) & (widths > 0) & (widths >= min_width)]

    table = _shares_table(raw.loc[widths.index], top_k, min_share)
    widths = widths.loc[table.index]
    palette = category_colors(table.columns, colors, other_color)

    if context:
        ax.add_collection(
            LineCollection(
                [np.asarray(line_coords(g)) for g in pipes.geometry],
                linewidths=0.3,
                color="0.8",
                zorder=zorder - 0.5,
            )
        )

    patches, facecolors, drawn = [], [], set()
    for reach in _reaches(net.graph, list(table.index)):
        # Split the reach where the shares or the width change.
        runs = [[reach[0]]]
        for pipe in reach[1:]:
            previous = runs[-1][-1]
            same = np.allclose(table.loc[pipe], table.loc[previous], rtol=1e-6) and np.isclose(
                widths[pipe], widths[previous], rtol=1e-6
            )
            if same:
                runs[-1].append(pipe)
            else:
                runs.append([pipe])
        for run in runs:
            row = table.loc[run[0]]
            row = row[row > 0]
            if order == "share":
                row = row.sort_values(ascending=False, kind="stable")
            line = _joined(pipes.geometry[p] for p in run)
            for category, strand in zip(
                row.index, ribbon_strands(line, row.to_numpy(), widths[run[0]]), strict=True
            ):
                path = _path(strand)
                if path is not None:
                    patches.append(PathPatch(path))
                    facecolors.append(palette[category])
                    drawn.add(category)

    ribbons = PatchCollection(
        patches, facecolors=facecolors, edgecolors="none", linewidths=0, zorder=zorder
    )
    ax.add_collection(ribbons)
    handles = [Patch(color=palette[c], label=str(c)) for c in table.columns if c in drawn]
    layer = DrawnLayer(ribbons, list(table.index), handles=handles, title=title)
    if legend:
        add_legends(ax, layer)
    _finish(ax, style)
    return layer


def node_group_polygons(
    net,
    groups,
    method="voronoi",
    *,
    clip_to=None,
    footprint=100.0,
    buffer=10.0,
    ratio=0.3,
    crs=None,
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
        grouped nodes only); a node may then be in only one group. 'concave_hull': the concave
        hull of each group's nodes (see ``ratio``) grown by ``buffer``; hulls may overlap.

    clip_to : shapely geometry, GeoSeries or GeoDataFrame, optional
        Area to clip to. A shapely geometry must be in the plot CRS; a GeoSeries or GeoDataFrame
        is reprojected and its geometries merged. Default: the network footprint (see
        ``footprint``), so that the polygons cover the serviced area.

    footprint : float or None
        Default clip area: the :func:`network_footprint` of ``net`` at this distance in metres.
        None: the convex hull of all the nodes of ``net`` grown by ``buffer`` for 'voronoi', and
        no clipping for 'concave_hull'.

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
    if clip_to is None and footprint is not None:
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
    ----------
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
        _colour_by(shaded, values, cmap, norm, vmin, vmax, missing_color)
        mappable = shaded
    elif color is not None:
        shaded.set_facecolor(color)
    else:
        palette = category_colors(ids, colors)
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
    ----------
    DrawnLayer
        ``artist``: the markers (a PathCollection), in the order of ``node_ids``.

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
    publication use. The extent of the axes is kept.

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
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        layer.to_crs(crs).plot(ax=ax, color=color, linewidth=linewidth, zorder=zorder, **kwargs)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
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
