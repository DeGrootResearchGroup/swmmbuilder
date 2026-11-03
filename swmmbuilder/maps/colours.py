"""Values to colours and sizes: scaling, norms, colour maps and category palettes."""

import copy
import logging

import matplotlib as mpl
import matplotlib.cbook as cbook
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, Normalize, to_hex

logger = logging.getLogger(__name__)


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
    to the largest value at the high end (so all zeros map to the low end); 'log' maps the range
    of the logarithms (equal values map to the high end). NaN, and values that the scale cannot
    take (negative, or non-positive for 'log'), map to the low end.
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
        # Equal values: all at the largest value (the high end), except that zeros on a
        # proportional scale are 0 (the low end).
        out = np.full(len(t), float(lo if scale != "log" and tmax == 0 else hi))
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


def _fitted_norm(norm, vmin, vmax, values):
    """(norm, has_data): ``norm`` (a name, or a Normalize, which is copied, never changed)
    with its unset limits scaled to the values it can show: the finite ones, and only the
    positive ones for a log norm.

    ``has_data`` is False when there are none (every value is NaN, or non-positive under a log
    norm); the limits are then a placeholder, every value is drawn in the missing colour, and
    callers draw no colour bar or legend. A zero-width range (all values equal) is widened
    around the value, so a colour bar or legend can still be drawn.
    """
    if isinstance(norm, Normalize):
        # A shallow copy with its own callbacks: a deep copy would also copy every artist the
        # caller's norm is already attached to.
        norm = copy.copy(norm)
        norm.callbacks = cbook.CallbackRegistry(signals=["changed"])
    else:
        norm = _norm(norm, vmin, vmax)
    log = isinstance(norm, LogNorm)
    values = np.asarray(values, dtype=float)
    valid = values[np.isfinite(values) & ((values > 0) if log else True)]
    has_data = len(valid) > 0
    if norm.vmin is None:
        norm.vmin = float(valid.min()) if has_data else (1.0 if log else 0.0)
    if norm.vmax is None:
        norm.vmax = float(valid.max()) if has_data else (10.0 if log else 1.0)
    if norm.vmin == norm.vmax:
        if log:
            norm.vmin, norm.vmax = norm.vmin / 10**0.5, norm.vmax * 10**0.5
        else:
            pad = 0.05 * abs(norm.vmin) or 0.05
            norm.vmin, norm.vmax = norm.vmin - pad, norm.vmax + pad
    return norm, has_data


def _colour_by(collection, values, cmap, norm, vmin, vmax, missing_color):
    """Colour ``collection`` by ``values`` through ``cmap``, with NaN (and, for a log norm,
    non-positive values) in ``missing_color``. Return False if no value could be coloured (see
    :func:`_fitted_norm`), in which case there is nothing for a colour bar to show."""
    norm, has_data = _fitted_norm(norm, vmin, vmax, values)
    collection.set_array(np.ma.masked_invalid(values))
    collection.set_cmap(mpl.colormaps.get_cmap(cmap).with_extremes(bad=missing_color))
    collection.set_norm(norm)
    if not has_data:
        logger.info("No value can be coloured (all missing, or non-positive on a log norm)")
    return has_data


def _hex_colours(values, cmap, norm, missing_color):
    """Hex colour of each value through ``cmap`` and a fitted ``norm`` (see
    :func:`_fitted_norm`), with NaN (and, for a log norm, non-positive values) in
    ``missing_color``."""
    masked = np.ma.masked_invalid(values)
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
