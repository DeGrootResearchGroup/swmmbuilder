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
- Elevations are taken to be in metres unless the raster says otherwise (a vertical CRS, or band
  units, in feet), or ``z_factor`` is given.
"""

from __future__ import annotations

import importlib
import logging
import math
import os
import warnings
from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ._summary import Summary
from .diagnostics import DataWarning

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

# Metres per unit of the elevation units a raster may declare (band units, lower case).
_Z_UNITS = {
    "m": 1.0,
    "metre": 1.0,
    "metres": 1.0,
    "meter": 1.0,
    "meters": 1.0,
    "ft": 0.3048,
    "foot": 0.3048,
    "feet": 0.3048,
    "international foot": 0.3048,
    "us survey foot": 1200 / 3937,
    "us-ft": 1200 / 3937,
    "ftus": 1200 / 3937,
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


def _vertical_factor(dataset: Any) -> tuple[float, str] | None:
    """(metres per elevation unit, where it was read) from a raster's vertical CRS or band
    units, or None if it declares neither."""
    crs = getattr(dataset, "crs", None)
    if crs is not None:
        from pyproj import CRS

        try:
            parts = CRS.from_wkt(crs.to_wkt()).sub_crs_list
        except Exception:  # an unusual CRS: fall back to the band units
            parts = []
        for part in parts:
            if part.is_vertical and part.axis_info:
                axis = part.axis_info[0]
                return float(axis.unit_conversion_factor), f"vertical CRS in {axis.unit_name}"
    units = (getattr(dataset, "units", None) or ("",))[0] or ""
    if units.strip():
        factor = _Z_UNITS.get(units.strip().lower())
        if factor is None:
            warnings.warn(
                f"DEM {dem_name(dataset)}: its elevations are in {units!r}, a unit swmmbuilder "
                "does not know, so they are taken as metres; pass z_factor (metres per unit) "
                "to convert them",
                DataWarning,
                stacklevel=2,
            )
            return None
        return factor, f"band units {units!r}"
    return None


def _z_factor(dataset: Any, z_factor: float | None) -> float:
    """Metres per elevation unit of a raster: ``z_factor`` if given, else what the raster
    declares, else 1 (metres)."""
    if z_factor is not None:
        return float(z_factor)
    declared = _vertical_factor(dataset)
    if declared is None:
        return 1.0
    factor, where = declared
    if factor != 1.0:
        logger.info(
            "DEM %s: elevations converted to metres (x %g, from its %s)",
            dem_name(dataset),
            factor,
            where,
        )
    return factor


def _bands(dataset: Any) -> int:
    """The number of bands of a raster that are not alpha (transparency) bands."""
    rasterio = _rasterio()
    alpha = rasterio.enums.ColorInterp.alpha
    interp = getattr(dataset, "colorinterp", None) or ()
    return dataset.count - sum(1 for c in interp if c == alpha)


def sample_dem(
    dem: StrPath | Any | Iterable[StrPath | Any],
    points: gpd.GeoSeries,
    window: float = 0.0,
    statistic: str = "median",
    z_factor: float | None = None,
) -> pd.Series:
    """The ground elevation of a DEM at each point.

    Parameters
    ----------
    dem : path, rasterio dataset, or list of them
        A single-band elevation raster (e.g. a GeoTIFF), or a list of tiles. Each raster may be
        in any CRS: the points are transformed to it, so the raster is not resampled. Cells the
        raster masks (its nodata value, a mask band or file, or an alpha band), and NaN cells,
        are left out.

    points : geopandas.GeoSeries
        The points. Missing points get NaN.

    window : float
        0 (default): the value of the cell containing the point. Otherwise the ``statistic`` of
        the cells whose centres are within ``window`` metres of the point (and of the cell
        containing it), over every tile the window reaches.

    statistic : {'median', 'mean'}
        How the cells of a window are combined.

    z_factor : float, optional
        Metres per elevation unit of the rasters (e.g. 0.3048 for a DEM in feet). Default: from
        each raster's vertical CRS, or else its band units, or else 1 (elevations in metres).

    Returns
    -------
    pandas.Series
        The elevations in metres, indexed like ``points``; NaN where no tile has a value at the
        point.

    Raises
    ------
    ImportError
        If rasterio is not installed.

    ValueError
        If ``window`` is negative, ``statistic`` is unknown, ``z_factor`` is not a positive
        number, a raster has more than one band (other than an alpha band), or the points
        have a CRS and a raster has none (or the other way round).

    Warns
    -----
    DataWarning
        If a raster declares elevation units swmmbuilder does not know (taken as metres).
    """
    if not window >= 0:
        raise ValueError(f"window must be 0 or more metres; got {window}.")
    if statistic not in STATISTICS:
        raise ValueError(f"statistic must be one of {list(STATISTICS)}; got {statistic!r}.")
    if z_factor is not None and not (
        isinstance(z_factor, (int, float)) and math.isfinite(z_factor) and z_factor > 0
    ):
        raise ValueError(f"z_factor must be a positive number of metres per unit; got {z_factor}.")
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
            factor = _z_factor(dataset, z_factor)
            for i, values in _sample_raster(dataset, xs, ys, window):
                found.setdefault(present.index[i], []).append(values * factor)
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
    if _bands(dataset) != 1:
        raise ValueError(
            f"DEM {dem_name(dataset)} has {dataset.count} bands; give a single-band raster "
            "(an alpha band besides it is allowed)."
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
    for i in inside:
        r0, r1 = rows[i] - reach[i], rows[i] + reach[i] + 1
        c0, c1 = cols[i] - reach[i], cols[i] + reach[i] + 1
        r0, c0 = max(r0, 0), max(c0, 0)
        r1, c1 = min(r1, dataset.height), min(c1, dataset.width)
        if r0 >= r1 or c0 >= c1:
            continue
        # masked=True leaves out the nodata value, a mask band or file, and an alpha band.
        read = dataset.read(1, window=((r0, r1), (c0, c1)), out_dtype="float64", masked=True)
        block = np.ma.getdata(read)
        valid = ~np.ma.getmaskarray(read) & ~np.isnan(block)
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

    def summary(self) -> dict[str, float | int]:
        """Counts and statistics of the comparison and the fill, to cite in a report.

        The keys follow the rules of every ``summary()`` (see the Diagnostics page of the
        documentation). :meth:`report` gives the same as text.

        Returns
        -------
        dict
            Plain Python numbers (``int`` for counts, ``float`` for values in m, NaN where there
            is nothing to measure), under these keys:

            - ``nodes``: nodes sampled (rows of :attr:`samples`).
            - ``compared``: those with a usable recorded rim and a DEM value (:attr:`compared`,
              including those flagged or replaced).
            - ``implausible``: recorded rims flagged as implausible, replaced or not
              (:attr:`implausible`).
            - ``replaced``: those replaced by the DEM estimate (:attr:`replaced`).
            - ``estimated``: nodes with no usable rim given an estimate (:attr:`estimated`).
            - ``no_dem_value``, ``dem_not_above_invert``: nodes of those statuses (see
              :data:`STATUSES`).
            - ``offset``: added to the DEM values (:attr:`offset`).
            - ``difference_<statistic>``: recorded rim - DEM over the nodes compared, for each
              statistic of :attr:`validation` but the count (``median``, ``mean``, ``std``,
              ``iqr``, ``p05``, ``p95``, and ``max_abs``).

        Examples
        --------
        ::

            s = est.summary()
            print(f"{s['estimated']} rims estimated; median difference {s['difference_median']:+.2f} m")
        """
        counts = self.samples["status"].value_counts()
        values: dict[str, float | int] = {
            "nodes": len(self.samples),
            "compared": len(self.compared),
            "implausible": len(self.implausible),
            "replaced": len(self.replaced),
            "estimated": len(self.estimated),
            "no_dem_value": int(counts.get("no_dem_value", 0)),
            "dem_not_above_invert": int(counts.get("dem_not_above_invert", 0)),
            "offset": float(self.offset),
        }
        validation = self.validation.drop("count")
        values.update({f"difference_{k}": float(v) for k, v in validation.items()})
        return Summary(values)

    def report(self) -> str:
        """A few lines describing the comparison and the fill, to print (also ``str()``)."""
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

    def __str__(self) -> str:
        return self.report()
