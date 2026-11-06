"""Ground elevations from a digital elevation model (DEM), and rims estimated from them.

A maintenance-hole lid in a street is generally close to flush with the pavement, so the ground
elevation at a node estimates its rim. :meth:`Network.estimate_rims` samples a bare-earth DEM
(e.g. one derived from airborne lidar) at the nodes, compares it with the recorded rims, and
estimates the rims that are missing or not above the invert; passing its :attr:`RimEstimate.rims`
to the ``rims`` option of :meth:`Subnetwork.to_inp` uses them. Usable recorded rims are never
replaced. This module does the sampling (:func:`sample_dem`) and holds the result
(:class:`RimEstimate`).

Reading a raster needs rasterio, an optional dependency: ``pip install 'swmmbuilder[elevation]'``.

Caveats:

- Lids in fields, ravines, or easements are often raised above grade, and structures under
  bridges, at creek crossings, or inside buildings cannot be read from the ground surface.
- A bare-earth model removes buildings and vegetation but may smooth sharp features such as curbs
  and ditches. A window statistic (``window``) softens steep edges, and the comparison with the
  recorded rims (:attr:`RimEstimate.validation`) tells whether the DEM is good enough.
- The DEM and the recorded elevations must share a vertical datum (in Ontario, CGVD28 and
  CGVD2013 differ by decimetres). The CRS of the DEM gives its horizontal coordinates only:
  check the datum, or convert, before comparing. A constant ``offset`` can apply a known shift.
"""

from __future__ import annotations

import importlib
import logging
import math
import os
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from typing import Any

    import geopandas as gpd

    from ._types import StrPath

logger = logging.getLogger(__name__)

STATISTICS = ("median", "mean")

# The status of each node in RimEstimate.samples.
STATUSES = {
    "compared": "usable recorded rim, compared with the DEM",
    "no_dem_value": "no DEM value at the node (outside the DEM, no data, or no point)",
    "estimated": "no usable rim; estimated from the DEM (in RimEstimate.rims)",
    "dem_not_above_invert": "no usable rim; the DEM estimate is not above the invert, not used",
    "implausible": "usable recorded rim, but farther from the DEM than ``flag``; kept",
    "replaced": "usable recorded rim beyond ``replace_above`` or ``replace_below``; replaced by "
    "the DEM estimate (in RimEstimate.rims)",
}

# Metres per degree of latitude (and of longitude at the equator), for a window on a DEM in
# geographic coordinates.
_METRES_PER_DEGREE = 111_320.0


def _rasterio():
    try:
        return importlib.import_module("rasterio")
    except ImportError as e:
        raise ImportError(
            "Reading a DEM needs rasterio, an optional dependency: "
            "pip install 'swmmbuilder[elevation]'"
        ) from e


def _sources(dem):
    """The DEM as a list of rasters: paths or open rasterio datasets."""
    if isinstance(dem, (str, os.PathLike)) or hasattr(dem, "read"):
        return [dem]
    if isinstance(dem, Iterable):
        sources = list(dem)
        if not sources:
            raise ValueError("dem: no rasters given.")
        return sources
    raise TypeError(
        f"dem must be a raster path, an open rasterio dataset, or a list of them; got "
        f"{type(dem).__name__}."
    )


def dem_name(dem: Any) -> str:
    """A short name of a DEM for reports: its file name, or the number of tiles."""
    sources = _sources(dem)
    if len(sources) > 1:
        return f"{len(sources)} tiles"
    source = sources[0]
    name = getattr(source, "name", source)
    return os.path.basename(str(name))


def sample_dem(
    dem: StrPath | Any | Iterable[StrPath | Any],
    points: gpd.GeoSeries,
    window: float = 0.0,
    statistic: str = "median",
) -> pd.Series:
    """The ground elevation of a DEM at each point.

    Parameters
    ----------
    dem : path, rasterio dataset, or list of them
        A single-band elevation raster (e.g. a GeoTIFF), or a list of tiles. Each raster may be
        in any CRS: the points are transformed to it, so the raster is not resampled. Cells equal
        to the raster's nodata value, and NaN cells, are left out.

    points : geopandas.GeoSeries
        The points. Missing points get NaN.

    window : float
        0 (default): the value of the cell containing the point. Otherwise the ``statistic`` of
        the cells whose centres are within ``window`` metres of the point (and of the cell
        containing it), over every tile the window reaches.

    statistic : {'median', 'mean'}
        How the cells of a window are combined.

    Returns
    -------
    pandas.Series
        The elevations, indexed like ``points``; NaN where no tile has a value at the point.

    Raises
    ------
    ImportError
        If rasterio is not installed.

    ValueError
        If ``window`` is negative, ``statistic`` is unknown, a raster has more than one band, or
        the points have a CRS and a raster has none (or the other way round).
    """
    if not window >= 0:
        raise ValueError(f"window must be 0 or more metres; got {window}.")
    if statistic not in STATISTICS:
        raise ValueError(f"statistic must be one of {list(STATISTICS)}; got {statistic!r}.")
    combine = np.median if statistic == "median" else np.mean
    rasterio = _rasterio()
    present = points[points.notna() & ~points.is_empty]
    # The points in each raster CRS, transformed once per CRS (tiles usually share one).
    projected: dict[Any, tuple[np.ndarray, np.ndarray]] = {}
    found: dict[Any, list[np.ndarray]] = {}
    for source in _sources(dem):
        opened = not hasattr(source, "read")
        dataset = rasterio.open(source) if opened else source
        try:
            if (points.crs is None) != (dataset.crs is None):
                raise ValueError(
                    f"DEM {dem_name(dataset)}: the CRS of the raster ({dataset.crs}) and of the "
                    f"network ({points.crs}) must both be known to sample it."
                )
            key = None if dataset.crs is None else dataset.crs.to_wkt()
            if key not in projected:
                moved = present if key is None else present.to_crs(key)
                projected[key] = (moved.x.to_numpy(), moved.y.to_numpy())
            xs, ys = projected[key]
            for i, values in _sample_raster(dataset, xs, ys, window):
                found.setdefault(present.index[i], []).append(values)
        finally:
            if opened:
                dataset.close()
    result = pd.Series(np.nan, index=points.index, dtype=float)
    for index, parts in found.items():
        values = np.concatenate(parts)
        if values.size:
            result[index] = float(combine(values))
    return result


def _sample_raster(dataset, xs, ys, window):
    """Yield (position in ``xs``, values of the cells sampled) for the points (in the raster's
    CRS) that the raster covers, within the window's reach."""
    if dataset.count != 1:
        raise ValueError(
            f"DEM {dem_name(dataset)} has {dataset.count} bands; give a single-band raster."
        )
    if not len(xs):
        return
    transform = dataset.transform
    cols, rows = _apply(~transform, xs, ys)
    cols, rows = np.floor(cols).astype(int), np.floor(rows).astype(int)

    # The window radius in cells, and the metres per raster unit in x and y at each point.
    geographic = dataset.crs is not None and dataset.crs.is_geographic
    if geographic:
        metres_x = _METRES_PER_DEGREE * np.cos(np.radians(ys))
        metres_y = np.full(len(ys), _METRES_PER_DEGREE)
    else:
        factor = dataset.crs.linear_units_factor[1] if dataset.crs is not None else 1.0
        metres_x = metres_y = np.full(len(ys), factor)
    if window > 0:
        cell_x = math.hypot(transform.a, transform.d)
        cell_y = math.hypot(transform.b, transform.e)
        reach = np.ceil(np.maximum(window / metres_x / cell_x, window / metres_y / cell_y)).astype(
            int
        )
    else:
        reach = np.zeros(len(ys), dtype=int)

    # Only the points whose window meets the raster (vectorized: a city's nodes against one of
    # many tiles).
    inside = np.flatnonzero(
        (rows + reach >= 0)
        & (rows - reach < dataset.height)
        & (cols + reach >= 0)
        & (cols - reach < dataset.width)
    )
    nodata = dataset.nodata
    for i in inside:
        r0, r1 = rows[i] - reach[i], rows[i] + reach[i] + 1
        c0, c1 = cols[i] - reach[i], cols[i] + reach[i] + 1
        r0, c0 = max(r0, 0), max(c0, 0)
        r1, c1 = min(r1, dataset.height), min(c1, dataset.width)
        if r0 >= r1 or c0 >= c1:
            continue
        block = dataset.read(1, window=((r0, r1), (c0, c1)), out_dtype="float64")
        valid = ~np.isnan(block)
        if nodata is not None and not math.isnan(nodata):
            valid &= block != nodata
        if reach[i] > 0:
            rr, cc = np.mgrid[r0:r1, c0:c1]
            cx, cy = _apply(transform, cc + 0.5, rr + 0.5)
            dx, dy = (cx - xs[i]) * metres_x[i], (cy - ys[i]) * metres_y[i]
            near = np.hypot(dx, dy) <= window
            near |= (rr == rows[i]) & (cc == cols[i])
            valid &= near
        yield i, block[valid]


def _apply(transform, x, y):
    """An affine transform of arrays of coordinates."""
    return (
        transform.a * x + transform.b * y + transform.c,
        transform.d * x + transform.e * y + transform.f,
    )


class Rims(dict):
    """{node: rim (m)} with the ``source`` the rims came from, written in the ``source`` column
    of the data-gaps report when they are passed as the ``rims`` option (see
    :attr:`RimEstimate.rims`), and ``implausible``: {node: what is wrong} of the recorded rims
    flagged as implausible, listed in the report as gap ``implausible_rim`` (filled where the
    node is also in the mapping, i.e. replaced)."""

    def __init__(self, rims: Any = (), source: str = "", implausible: Any = ()) -> None:
        super().__init__(rims)
        self.source = source
        self.implausible = dict(implausible)


@dataclass
class RimEstimate:
    """Rims estimated from a DEM by :meth:`Network.estimate_rims`. The network is not changed:
    pass :attr:`rims` to the ``rims`` option of :meth:`Subnetwork.to_inp` (or
    :meth:`Subnetwork.data_gaps`) to use them.

    Attributes
    ----------
    samples : pandas.DataFrame
        One row per node sampled, indexed by node ID: ``recorded_rim`` (NaN if none),
        ``invert``, ``dem`` (the DEM value, NaN if none), ``difference`` (recorded rim - DEM, for
        the nodes compared), ``deviation`` (``difference`` minus its median: how far the
        recorded rim is from the DEM once a datum shift is taken out; positive when recorded
        above the ground), ``estimate`` (DEM + ``offset``), and ``status`` (see
        :data:`STATUSES`): 'compared', 'implausible', 'replaced', 'no_dem_value', 'estimated',
        or 'dem_not_above_invert'.

    offset : float
        Added to the DEM values to estimate the rims, in m (the median difference if
        ``offset='median'`` was given).

    source : str
        What the estimates came from, as written in the ``source`` column of the data-gaps
        report.

    flag, replace_above, replace_below : float or None
        The thresholds the estimate was made with (see :meth:`Network.estimate_rims`).
    """

    samples: pd.DataFrame
    offset: float
    source: str
    flag: float | None = None
    replace_above: float | None = None
    replace_below: float | None = None

    @property
    def compared(self) -> pd.DataFrame:
        """The rows of ``samples`` with a usable recorded rim and a DEM value (including those
        flagged implausible or replaced)."""
        return self.samples[self.samples["difference"].notna()]

    @property
    def implausible(self) -> pd.DataFrame:
        """The rows of ``samples`` whose recorded rim was flagged as implausible, replaced or
        not."""
        return self.samples[self.samples["status"].isin(("implausible", "replaced"))]

    @property
    def replaced(self) -> pd.DataFrame:
        """The rows of ``samples`` whose recorded rim is replaced by the DEM estimate."""
        return self.samples[self.samples["status"] == "replaced"]

    @property
    def estimated(self) -> pd.DataFrame:
        """The rows of ``samples`` with no usable rim and an estimate above the invert."""
        return self.samples[self.samples["status"] == "estimated"]

    @property
    def rims(self) -> Rims:
        """{node: estimated rim} of the nodes with no usable rim whose estimate is above the
        invert, and of the recorded rims replaced, for the ``rims`` option of
        :meth:`Subnetwork.to_inp`. It fills the first (gap ``no_rim``) and replaces the second
        (gap ``implausible_rim``), both filled from ``estimate_rims`` with :attr:`source` as
        their source; the rims flagged but kept are listed as ``implausible_rim`` gaps not
        filled."""
        rows = self.samples[self.samples["status"].isin(("estimated", "replaced"))]
        flagged = {
            node: f"recorded rim {row.recorded_rim:g}, {row.deviation:+.2f} m from the DEM"
            for node, row in self.implausible.iterrows()
        }
        return Rims({n: float(v) for n, v in rows["estimate"].items()}, self.source, flagged)

    @property
    def validation(self) -> pd.Series:
        """Recorded rim - DEM over the nodes compared: ``count``, ``median``, ``mean``, ``std``,
        ``iqr`` (interquartile range), ``p05`` and ``p95`` (5th and 95th percentiles), and
        ``max_abs`` (the largest disagreement), in m. A median far from 0 suggests a vertical
        datum or a systematic difference; a wide spread, a DEM (or recorded rims) not good
        enough to fill from."""
        d = self.compared["difference"]
        if d.empty:
            values = {"count": 0} | dict.fromkeys(
                ("median", "mean", "std", "iqr", "p05", "p95", "max_abs"), math.nan
            )
        else:
            values = {
                "count": len(d),
                "median": d.median(),
                "mean": d.mean(),
                "std": d.std(),
                "iqr": d.quantile(0.75) - d.quantile(0.25),
                "p05": d.quantile(0.05),
                "p95": d.quantile(0.95),
                "max_abs": d.abs().max(),
            }
        return pd.Series(values, dtype=float)

    def largest_differences(self, n: int = 10) -> pd.DataFrame:
        """The ``n`` nodes compared whose recorded rim disagrees most with the DEM, measured from
        the median difference (so that a datum shift does not hide them): recorded rims that
        are probably wrong, or lids that are not flush with the ground."""
        compared = self.compared
        if compared.empty:
            return compared
        spread = (compared["difference"] - compared["difference"].median()).abs()
        return compared.loc[spread.sort_values(ascending=False, kind="stable").index[:n]]

    def summary(self) -> str:
        """A few lines describing the comparison and the fill."""
        v = self.validation
        counts = self.samples["status"].value_counts()
        lines = [f"Rims from {self.source}: {len(self.samples)} nodes sampled"]
        if v["count"]:
            lines.append(
                f"  recorded - DEM over {int(v['count'])} recorded rims: median "
                f"{v['median']:+.3f} m, IQR {v['iqr']:.3f} m, 5-95% {v['p05']:+.3f} to "
                f"{v['p95']:+.3f} m, largest {v['max_abs']:.3f} m"
            )
        else:
            lines.append("  no recorded rim to compare with")
        if self.flag is not None:
            lines.append(
                f"  implausible (more than {self.flag:g} m from the median difference): "
                f"{len(self.implausible)}, of which {len(self.replaced)} replaced"
            )
        for status in STATUSES:
            if status not in ("compared", "implausible", "replaced") and counts.get(status, 0):
                lines.append(f"  {status}: {counts[status]}")
        return "\n".join(lines)
