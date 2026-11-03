"""Map axes (scale bar, north arrow) and the column of keys (colour bars and legends)
to the right of them."""

import weakref
from dataclasses import dataclass, field

import matplotlib as mpl
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.legend import Legend
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, NullFormatter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


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
    """The state of ``ax``, started afresh if ``ax`` was cleared since (its decorations and
    keys are then gone, and the column of keys is empty)."""
    state = _AXES_STATE.get(ax)
    if state is None or _cleared(ax, state):
        state = _AXES_STATE[ax] = {"key_top": 1.0, "decorations": [], "keys": []}
    return state


def _cleared(ax, state):
    """True if an artist or key axes recorded in ``state`` is no longer on ``ax`` (e.g. after
    ``ax.clear()``, which detaches artists and drops inset axes)."""
    for artist in [*state["decorations"], *state["keys"]]:
        if artist in ax.figure.axes:  # a colour bar's inset axes
            continue
        if isinstance(artist, mpl.axes.Axes) or artist.axes is not ax:
            return True
    return False


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
            state["keys"].append(cax)
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
            state["keys"].append(legend)
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
