"""Ribbons: pipes drawn as bands of strands, one per share of the flow, along main
stems with mitred joins."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import shapely
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Patch, PathPatch
from shapely import LineString

from ..geometry import line_coords
from .colours import NEUTRAL_COLOR, _aligned, category_colors
from .crs import _geodata, network_extent
from .polygons import _path
from .style import DrawnLayer, _finish, add_legends

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from ..network import Network
    from ..subnetwork import Subnetwork

# Name of the strand that lumps the categories dropped by plot_ribbons' top_k and min_share.
OTHER = "other"


# Defaults relative to the size of the network drawn, so that they suit a 30-pipe catchment and a
# city alike: ribbons at most RIBBON_WIDTH x the extent wide, ribbons narrower than
# MIN_RIBBON_WIDTH x the widest left as context, and footprints of FOOTPRINT_PIPES x the median
# pipe length, at most FOOTPRINT_EXTENT x the extent.
RIBBON_WIDTH = 0.025


MIN_RIBBON_WIDTH = 0.05


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
    -------
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


def _stems(graph, pipe_ids, widths=None):
    """``pipe_ids`` as chains in flow order along main stems. At each node, the widest incoming
    pipe in ``pipe_ids`` (by ``widths``, a Series by pipe ID; else the first) continues into the
    widest outgoing one; other pipes start or end a chain there. Pipes not in ``graph`` are
    chains of their own."""
    wanted = set(pipe_ids)
    ends = {key: (up, down) for up, down, key in graph.edges(keys=True) if key in wanted}
    rank = {key: i for i, key in enumerate(pipe_ids)}

    def widest(keys):
        keys = [k for k in keys if k in ends]
        if not keys:
            return None
        if widths is None:
            return min(keys, key=rank.get)
        return max(keys, key=lambda k: (widths[k], -rank[k]))

    following = {}
    for node in {n for ends_ in ends.values() for n in ends_}:
        main_in = widest(k for _, _, k in graph.in_edges(node, keys=True))
        main_out = widest(k for _, _, k in graph.out_edges(node, keys=True))
        if main_in is not None and main_out is not None and main_in != main_out:
            following[main_in] = main_out
    preceding = {b: a for a, b in following.items()}

    chains, seen = [], set()
    # Chain heads first, then what is left: pipes on cycles and pipes not in the graph.
    starts = [k for k in pipe_ids if k in ends and k not in preceding] + list(pipe_ids)
    for key in starts:
        if key in seen:
            continue
        chain, members = [key], {key}
        while (after := following.get(chain[-1])) is not None and not (
            after in seen or after in members
        ):
            chain.append(after)
            members.add(after)
        seen.update(chain)
        chains.append(chain)
    return chains


def _merged(chain, table, widths, width_tol, share_tol):
    """``chain`` split into runs, each drawn with the shares and width of its downstream pipe:
    walking upstream, a run continues while a pipe's width is within ``width_tol`` (metres) and
    each of its shares within ``share_tol`` of the run's downstream pipe."""
    runs = [[chain[-1]]]
    for pipe in reversed(chain[:-1]):
        last = runs[-1][-1]
        if abs(widths[pipe] - widths[last]) <= width_tol + 1e-9 and (
            np.abs(table.loc[pipe].to_numpy() - table.loc[last].to_numpy()).max()
            <= share_tol + 1e-9
        ):
            runs[-1].insert(0, pipe)
        else:
            runs.append([pipe])
    return runs[::-1]


def _direction(a, b):
    d = np.asarray(b, dtype=float)[:2] - np.asarray(a, dtype=float)[:2]
    n = np.hypot(*d)
    return d / n if n > 0 else None


def _half_plane(point, normal, size):
    """The square of side ``2 * size`` on the side of the line through ``point`` that ``normal``
    points to."""
    p = np.asarray(point, dtype=float)[:2]
    t = np.array([-normal[1], normal[0]])
    return shapely.Polygon(
        [p - t * size, p + t * size, p + (t + normal) * size, p + (normal - t) * size]
    )


def mitred_strands(line, shares, width, before=None, after=None):
    """:func:`ribbon_strands` along ``line``, with the ends cut on the bisector of the angle
    with the neighbouring pieces, so that consecutive pieces meet without a gap or an overlap.

    ``before`` is the direction (a 2-vector) in which the line arrives at its start from the
    piece upstream, ``after`` the direction in which it leaves its end into the piece
    downstream; None leaves that end cut square, as in :func:`ribbon_strands`.
    """
    coords = [tuple(c[:2]) for c in line.coords]
    first = _direction(coords[0], coords[1])
    last = _direction(coords[-2], coords[-1])
    reach = 2 * width + 1.0
    cuts = []
    if before is not None and first is not None:
        normal = before + first
        if np.hypot(*normal) > 1e-6:
            coords.insert(0, tuple(np.asarray(coords[0]) - before * reach))
            cuts.append((line.coords[0], normal / np.hypot(*normal)))
    if after is not None and last is not None:
        normal = last + after
        if np.hypot(*normal) > 1e-6:
            coords.append(tuple(np.asarray(coords[-1]) + after * reach))
            cuts.append((line.coords[-1], -normal / np.hypot(*normal)))
    strands = ribbon_strands(LineString(coords), shares, width)
    size = line.length + 4 * reach
    for point, normal in cuts:
        keep = _half_plane(point, normal, size)
        strands = [s.intersection(keep) if not s.is_empty else s for s in strands]
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


def plot_ribbons(
    ax: Axes,
    net: Network | Subnetwork,
    shares: pd.DataFrame,
    colors=None,
    *,
    width="total",
    max_width=None,
    min_width=None,
    top_k=None,
    min_share=0.0,
    other_color=NEUTRAL_COLOR,
    order="category",
    merge_tol=(0.05, 0.02),
    context=True,
    legend=True,
    title=None,
    style=True,
    crs=None,
    zorder=2,
) -> DrawnLayer:
    """Draw pipes as ribbons of parallel strands, one per category, each as wide as its share of
    the pipe's load (a mass-flow map).

    Pipes are joined into runs along main stems (at each node, the widest incoming ribbon
    continues into the widest outgoing one), and a run is split only where its ribbon visibly
    changes (see ``merge_tol``). Each run is drawn as one ribbon with the shares and width of its
    downstream pipe, and consecutive runs are cut on the bisector of the angle between them, so
    a trunk reads as one band with no gaps or overlaps at its bends.

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

    max_width : float, optional
        Width in metres of the pipe with the largest total, for ``width='total'``. Default: 2.5 %
        of :func:`network_extent`, so that the ribbons are as visible on a city as on a street.

    min_width : float, optional
        Pipes narrower than this (in metres) are not drawn as ribbons; with ``context`` they
        are still drawn as thin grey lines. Default: 5 % of the widest ribbon.

    top_k : int, optional
        Keep only each pipe's ``top_k`` largest categories and lump the rest into 'other'.

    min_share : float
        Lump categories with a smaller share of the pipe into 'other'.

    other_color
        Colour of the 'other' strand.

    order : {'category', 'share'}
        Order of the strands from the right of the pipe to its left (looking downstream): the
        order of the columns (with 'other' last), or by decreasing share.

    merge_tol : float or (float, float)
        ``(width, share)``: a run continues while the width differs from that of its downstream
        pipe by at most ``width`` times the widest ribbon, and every share by at most ``share``.
        A number sets both. 0 joins only pipes with the same width and shares.

    context : bool
        Draw every pipe as a thin grey line beneath the ribbons.

    legend : bool
        Add a legend (titled ``title``) of the categories drawn.

    style : bool
        Style the axes as a map (:func:`map_axes`).

    crs : optional
        CRS to draw in; default :func:`plot_crs`.

    Returns
    -------
    DrawnLayer
        ``artist``: a PatchCollection with one patch per strand; ``ids``: the pipe IDs drawn as
        ribbons; ``handles``: one per category drawn.
    """
    if order not in ("category", "share"):
        raise ValueError(f"order must be 'category' or 'share', not {order!r}.")
    width_tol, share_tol = (merge_tol, merge_tol) if np.isscalar(merge_tol) else merge_tol
    pipes = _geodata(net, "pipes", crs)
    raw = shares.astype(float).fillna(0.0)
    raw = raw.loc[raw.index.isin(pipes.index)]
    if isinstance(width, str):
        if width != "total":
            raise ValueError(f"width must be 'total', a number or per-pipe values, not {width!r}.")
        if max_width is None:
            max_width = RIBBON_WIDTH * network_extent(net, crs)
        totals = raw.sum(axis=1)
        widths = totals / totals.max() * max_width if len(totals) else totals
    elif np.isscalar(width):
        widths = pd.Series(float(width), index=raw.index)
    else:
        widths = pd.Series(_aligned(width, raw.index, "widths"), index=raw.index)
    widths = widths[np.isfinite(widths) & (widths > 0)]
    if min_width is None:
        min_width = MIN_RIBBON_WIDTH * widths.max() if len(widths) else 0.0
    widths = widths[widths >= min_width]

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

    width_tol *= widths.max() if len(widths) else 0.0
    runs = [
        _merged(chain, table, widths, width_tol, share_tol)
        for chain in _stems(net.graph, list(table.index), widths)
    ]
    lines = {run[-1]: _joined(pipes.geometry[p] for p in run) for chain in runs for run in chain}

    def heading(line, end):
        coords = list(line.coords)
        return _direction(*coords[-2:]) if end else _direction(*coords[:2])

    patches, facecolors, drawn = [], [], set()
    for chain in runs:
        for i, run in enumerate(chain):
            row = table.loc[run[-1]]
            row = row[row > 0]
            if order == "share":
                row = row.sort_values(ascending=False, kind="stable")
            before = heading(lines[chain[i - 1][-1]], end=True) if i > 0 else None
            after = heading(lines[chain[i + 1][-1]], end=False) if i + 1 < len(chain) else None
            strands = mitred_strands(lines[run[-1]], row.to_numpy(), widths[run[-1]], before, after)
            for category, strand in zip(row.index, strands, strict=True):
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
