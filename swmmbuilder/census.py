"""Census polygons with their population, to pass to :func:`~swmmbuilder.allocate_population`.

:func:`census_polygons` gets the polygons covering a network (or any area) from a *source*, keeps
those near its pipes, and returns them with standard columns. A cache file makes a model rebuild
offline. The source is pluggable: anything with a ``name`` and a ``fetch`` method (see
:class:`PolygonSource`) will do, and a source with a ``provenance`` method has its cache checked
against its parameters. Two are provided:

- :class:`CensusMapper`, Canada's 2021 Census (or another CensusMapper dataset) through the
  CensusMapper API, with pycancensus (``pip install 'swmmbuilder[census]'``);
- :class:`PolygonFile`, any local file of polygons with a population column, e.g. another
  country's census.

Canadian census terms: a *dissemination block* (``"DB"``) is the smallest census area, usually
one city block bounded by streets; only its population, dwellings, and households are published.
A *dissemination area* (``"DA"``) is a group of neighbouring blocks with about 400 to 700 people.
Blocks suit :func:`~swmmbuilder.allocate_population`'s ``max_distance`` rule, since their edges
follow the streets the sewers run under; dissemination areas suit the midpoint rule.

**The API key and the citation are the user's.** swmmbuilder provides no API key and ships no census data. A CensusMapper key is yours to get
(free, with a CensusMapper account) and to provide, and its quota is yours. If you use the data,
citing it is your responsibility: Statistics Canada's data are under the Statistics Canada Open
Licence, which asks for a source notice, and CensusMapper asks that cancensus be cited. The
cache's provenance (source, dataset, regions, and fetch date) is recorded to help.
"""

from __future__ import annotations

import importlib
import json
import os
import warnings
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib import metadata
from pathlib import Path
from typing import Any, Literal, Protocol, runtime_checkable

import geopandas as gpd
import pandas as pd
import shapely
from pyproj import CRS

from .diagnostics import DataWarning, Problems
from .maps.crs import _in_metres, network_footprint, plot_crs
from .network import Network
from .subnetwork import Subnetwork

KEY_VARIABLES = ("CENSUSMAPPER_API_KEY", "CANCENSUS_API_KEY")
SIGN_UP = "https://censusmapper.ca/users/sign_up"
POLYGON_LAYER = "polygons"
AREA_LAYER = "area"
# A cached area must cover the area asked for to within this distance (m), for reprojection.
COVER_TOLERANCE = 1.0


@runtime_checkable
class PolygonSource(Protocol):
    """A source of polygons with a population: what :func:`census_polygons` fetches from.

    Attributes
    ----------
    name : str
        A short name, recorded in the cache, e.g. ``"CensusMapper"``.

    Notes
    -----
    A source may also have a ``provenance()`` method, returning a dict of str of the parameters
    that decide what it fetches (e.g. the dataset, or a file and its columns). When a cache is
    read with such a source, :func:`census_polygons` compares them with those recorded in the
    cache and refuses a cache fetched with other values. :class:`CensusMapper` and
    :class:`PolygonFile` have one.
    """

    name: str

    def fetch(self, area: shapely.Geometry, crs: CRS, level: str) -> gpd.GeoDataFrame:
        """The polygons of ``level`` that meet ``area`` (in ``crs``), and possibly more.

        Returns a GeoDataFrame with a CRS, indexed by polygon ID, with ``population`` and,
        optionally, ``parent`` (the ID of the enclosing area). Anything in
        ``.attrs["provenance"]`` (a dict of str) is recorded in the cache.
        """
        ...


def _optional(module):
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(
            f"This needs {module}, an optional dependency: pip install 'swmmbuilder[census]'"
        ) from e


def _version(package):
    try:
        return metadata.version(package)
    except metadata.PackageNotFoundError:
        return "unknown"


def _id(value):
    if pd.isna(value):
        return None
    if isinstance(value, float) and value.is_integer():
        value = int(value)  # a column of integer IDs with gaps is read as float
    return str(value)


def _ids(values):
    """IDs as str, keeping missing values missing."""
    return pd.Series(values, dtype=object).map(_id).to_numpy(dtype=object)


def _frame(ids, population, parent, geometry, crs):
    """The standard GeoDataFrame of a source: indexed by ``id``, with population and parent."""
    return gpd.GeoDataFrame(
        {"population": population, "parent": _ids(parent) if parent is not None else None},
        geometry=geometry,
        index=pd.Index(_ids(ids), name="id"),
        crs=crs,
    )


def _regions_text(regions: Any) -> str:
    """Regions ({level: [IDs]}) as JSON text with the levels and IDs sorted, so that the same
    regions are always recorded the same way."""
    if isinstance(regions, str):
        regions = json.loads(regions)
    return json.dumps(
        {str(level): sorted(str(i) for i in ids) for level, ids in regions.items()},
        sort_keys=True,
    )


def _missing_columns(data, columns, what):
    missing = [c for c in columns if c is not None and c not in data.columns]
    if missing:
        raise ValueError(f"{what} has no column {', '.join(map(repr, missing))}.")


@dataclass
class CensusMapper:
    """Census polygons from the CensusMapper API (censusmapper.ca), through pycancensus.

    Needs ``pip install 'swmmbuilder[census]'`` and a CensusMapper API key, which you provide:
    swmmbuilder has none. Get a free key with a CensusMapper account; its quota is yours.

    Parameters
    ----------
    api_key : str, optional
        Your CensusMapper API key. Default: the ``CENSUSMAPPER_API_KEY`` environment variable, or
        else ``CANCENSUS_API_KEY``. It is passed to each pycancensus call and never stored,
        logged, or shown (it is left out of ``repr``). It is only needed when there is no cache.

    dataset : str
        The CensusMapper dataset (default ``"CA21"``, the 2021 Census).

    regions : dict, optional
        The census regions to fetch, as pycancensus takes them, e.g. ``{"CSD": ["3539036"]}``
        (the City of London). Default: the regions of ``region_level`` meeting the area, found
        with CensusMapper's intersecting-geometries query (which has its own quota).

    region_level : str
        The level of the regions found (default ``"CSD"``, census subdivisions: municipalities).

    resolution : {'high', 'simplified'}
        The boundaries' resolution (default ``"high"``).

    Notes
    -----
    pycancensus keeps its own cache of what it downloads, in ``CANCENSUS_CACHE_PATH`` or
    ``~/.cancensus_cache``; the key is not in it. swmmbuilder never reads that cache (each fetch
    is fresh), nor pycancensus's saved settings.
    """

    api_key: str | None = field(default=None, repr=False)
    dataset: str = "CA21"
    regions: dict[str, list[str]] | None = None
    region_level: str = "CSD"
    resolution: Literal["high", "simplified"] = "high"
    name: str = field(default="CensusMapper", init=False)

    def provenance(self) -> dict[str, str]:
        """The parameters that decide what is fetched, as recorded in a cache: ``dataset``,
        ``resolution``, and ``regions`` (as JSON text) if given, or else ``region_level``.

        Returns
        -------
        dict of str
            Compared with a cache's provenance when the cache is read with this source.
        """
        found = {"dataset": self.dataset, "resolution": self.resolution}
        if self.regions:
            found["regions"] = _regions_text(self.regions)
        else:
            found["region_level"] = self.region_level
        return found

    def _key(self) -> str:
        if self.api_key:
            return self.api_key
        for variable in KEY_VARIABLES:
            if os.environ.get(variable):
                return os.environ[variable]
        raise ValueError(
            "CensusMapper needs an API key, and it is yours to provide: swmmbuilder has none. "
            f"A free key comes with a CensusMapper account ({SIGN_UP}). Pass it as "
            f"CensusMapper(api_key=...) or set {' or '.join(KEY_VARIABLES)}."
        )

    def _intersecting(self, pc, area, crs, key):
        """{level: [region IDs]} of the regions meeting ``area``, from a simplified outline of
        it (enlarged by the simplification tolerance, so it still covers the area)."""
        xmin, ymin, xmax, ymax = area.bounds
        tolerance = 0.01 * max(xmax - xmin, ymax - ymin, 1.0)
        outline = shapely.simplify(area, tolerance).buffer(tolerance)
        found = pc.get_intersecting_geometries(
            dataset=self.dataset,
            level=self.region_level,
            geometry=gpd.GeoSeries([outline], crs=crs),
            simplified=False,
            use_cache=False,
            quiet=True,
            api_key=key,
        )
        if isinstance(found, list):
            found = {self.region_level: found}
        regions = {level: [str(r) for r in ids] for level, ids in (found or {}).items() if ids}
        if not regions:
            raise ValueError(
                f"CensusMapper found no {self.region_level} regions of {self.dataset} meeting "
                "the area. Is it in Canada? Pass regions= to name them."
            )
        return regions

    def fetch(self, area: shapely.Geometry, crs: CRS, level: str) -> gpd.GeoDataFrame:
        """The polygons of ``level`` (e.g. ``"DB"`` or ``"DA"``) in the regions meeting
        ``area``, with ``population`` (CensusMapper's ``pop``) and ``parent`` (``rpid``)."""
        pc = _optional("pycancensus")
        key = self._key()
        regions = self.regions or self._intersecting(pc, area, crs, key)
        data = pc.get_census(
            dataset=self.dataset,
            regions=regions,
            level=level,
            geo_format="geopandas",
            resolution=self.resolution,
            use_cache=False,
            quiet=True,
            api_key=key,
        )
        _missing_columns(data, ["id", "pop"], "CensusMapper's response")
        out = _frame(
            data["id"],
            data["pop"].to_numpy(),
            data["rpid"] if "rpid" in data.columns else None,
            data.geometry.to_numpy(),
            data.crs,
        )
        out.attrs["provenance"] = {
            "dataset": self.dataset,
            "regions": _regions_text(regions),
            "region_level": self.region_level,
            "resolution": self.resolution,
            "pycancensus": _version("pycancensus"),
        }
        return out


@dataclass
class PolygonFile:
    """Polygons with a population from a local file (any format geopandas reads), e.g. another
    country's census or a city's own data.

    Parameters
    ----------
    path : str or path
        The file.
    id : str
        The column of unique polygon IDs.
    population : str
        The population column.
    parent : str, optional
        A column of the enclosing area's ID.
    layer : str, optional
        The layer, for a file with several.

    Notes
    -----
    The whole file is read; :func:`census_polygons` keeps the polygons meeting the area. The
    ``level`` asked for is recorded but not used. A cache records the file's name (not its
    folder, so a project can move) and the columns, and is refused for another file or
    other columns.
    """

    path: str | os.PathLike
    id: str
    population: str
    parent: str | None = None
    layer: str | None = None
    name: str = field(default="file", init=False)

    def provenance(self) -> dict[str, str]:
        """The file's name and the columns read, as recorded in a cache: ``file``, ``id``,
        ``population``, ``parent``, and ``layer`` (an empty string where not given).

        Returns
        -------
        dict of str
            Compared with a cache's provenance when the cache is read with this source.
        """
        return {
            "file": Path(self.path).name,
            "id": self.id,
            "population": self.population,
            "parent": self.parent or "",
            "layer": self.layer or "",
        }

    def fetch(self, area: shapely.Geometry, crs: CRS, level: str) -> gpd.GeoDataFrame:
        """Every polygon of the file, with ``population`` and ``parent``."""
        data = gpd.read_file(self.path, layer=self.layer)
        _missing_columns(data, [self.id, self.population, self.parent], str(self.path))
        out = _frame(
            data[self.id],
            data[self.population].to_numpy(),
            data[self.parent] if self.parent else None,
            data.geometry.to_numpy(),
            data.crs,
        )
        out.attrs["provenance"] = self.provenance()
        return out


def _region(area, crs):
    """(geometry, CRS in metres) of ``area``."""
    if isinstance(area, (Network, Subnetwork)):
        work = plot_crs(area)
        return network_footprint(area, crs=work, pipes="gravity"), work
    if isinstance(area, shapely.Geometry):
        if crs is None:
            raise ValueError("A shapely area needs crs=, the CRS of its coordinates.")
        area = gpd.GeoSeries([area], crs=crs)
    if isinstance(area, (gpd.GeoDataFrame, gpd.GeoSeries)):
        if area.crs is None:
            raise ValueError("The area has no CRS.")
        geometry = area.geometry
        work = CRS.from_user_input(area.crs)
        if not _in_metres(work):
            work = geometry.estimate_utm_crs()
        return shapely.union_all(geometry.to_crs(work).to_numpy()), work
    raise TypeError(
        "area must be a Network, a Subnetwork, a GeoDataFrame, a GeoSeries, or a shapely "
        f"geometry, not {type(area).__name__}."
    )


def _standard(data, source):
    """The polygons of a source, checked: indexed by str ID, with numeric population."""
    what = f"The polygons of {source.name}"
    if not isinstance(data, gpd.GeoDataFrame):
        raise TypeError(f"{what} must be a GeoDataFrame, not {type(data).__name__}.")
    problems = Problems(f"Invalid polygons from {source.name}")
    if data.crs is None:
        problems.add("they have no CRS")
    if "population" not in data.columns:
        problems.add("there is no population column")
    repeated = data.index[data.index.duplicated()].unique()
    if len(repeated):
        problems.add(f"repeated IDs: {', '.join(map(str, repeated[:10]))}")
    problems.raise_if_any()
    parent = data["parent"] if "parent" in data.columns else None
    population = pd.to_numeric(pd.Series(data["population"], dtype=object), errors="coerce")
    return _frame(
        data.index, population.to_numpy(dtype=float), parent, data.geometry.to_numpy(), data.crs
    )


def _write_cache(path, polygons, region, crs, provenance):
    """Write the polygons and the area they were fetched for, with the provenance, to a
    GeoPackage, replacing ``path`` only once both layers are written."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.stem}.partial.gpkg")
    partial.unlink(missing_ok=True)
    try:
        polygons.reset_index().to_file(partial, layer=POLYGON_LAYER, driver="GPKG")
        area = gpd.GeoDataFrame({k: [v] for k, v in provenance.items()}, geometry=[region], crs=crs)
        area.to_file(partial, layer=AREA_LAYER, driver="GPKG")
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


def _read_cache(path, level, source, region, crs):
    """The polygons and provenance of the cache at ``path``, checked against the request."""
    try:
        polygons = gpd.read_file(path, layer=POLYGON_LAYER)
        area = gpd.read_file(path, layer=AREA_LAYER)
    except Exception as e:
        raise ValueError(f"{path} is not a census_polygons cache: {e}") from e
    provenance = {
        k: str(v) for k, v in area.drop(columns="geometry").iloc[0].items() if pd.notna(v)
    }
    if provenance.get("level") != level:
        raise ValueError(
            f"{path} holds level {provenance.get('level')!r}, not {level!r}; use another cache "
            "file, or refresh=True to replace it."
        )
    if source is not None and provenance.get("source") != source.name:
        raise ValueError(
            f"{path} was fetched from {provenance.get('source')!r}, not {source.name!r}; use "
            "another cache file, or refresh=True to replace it."
        )
    if source is not None:
        _check_provenance(path, provenance, source)
    asked = gpd.GeoSeries([region], crs=crs).to_crs(area.crs).iloc[0]
    if not area.geometry.iloc[0].buffer(COVER_TOLERANCE).covers(asked):
        raise ValueError(
            f"{path} covers a smaller area than the one asked for (a larger network or margin?); "
            "refresh=True fetches it again."
        )
    polygons = polygons.set_index("id")
    return _standard(polygons, _Named(provenance.get("source", "the cache"))), provenance


def _check_provenance(path, cached, source):
    """Refuse a cache whose recorded parameters differ from those of ``source`` (its
    ``provenance()``, if it has one); warn about those the cache does not record."""
    method = getattr(source, "provenance", None)
    if not callable(method):
        return
    asked = {k: str(v) for k, v in method().items()}
    if "regions" in cached and "regions" in asked:
        try:
            cached = {**cached, "regions": _regions_text(cached["regions"])}
        except (ValueError, AttributeError, TypeError):
            pass  # not JSON of {level: [IDs]}: compared as text
    differ = [
        f"{k}={cached[k]!r}, not {v!r}" for k, v in asked.items() if k in cached and cached[k] != v
    ]
    if differ:
        raise ValueError(
            f"{path} was fetched with {'; '.join(differ)}; use another cache file, or "
            "refresh=True to replace it."
        )
    unrecorded = [k for k in asked if k not in cached]
    if unrecorded:
        warnings.warn(
            f"{path} does not record {', '.join(unrecorded)} (written by an older swmmbuilder?), "
            "so the cache could not be checked against them; refresh=True records them",
            DataWarning,
            stacklevel=4,
        )


@dataclass
class _Named:
    name: str


def census_polygons(
    area: Network | Subnetwork | gpd.GeoDataFrame | gpd.GeoSeries | shapely.Geometry,
    level: str = "DB",
    source: PolygonSource | None = None,
    cache: str | os.PathLike | None = None,
    refresh: bool = False,
    margin: float = 0.0,
    crs: Any = None,
) -> gpd.GeoDataFrame:
    """Census polygons meeting ``area``, with their population, for
    :func:`~swmmbuilder.allocate_population`.

    Parameters
    ----------
    area : Network, Subnetwork, GeoDataFrame, GeoSeries, or shapely geometry
        Where polygons are wanted. For a network, the footprint of its gravity pipes
        (:func:`~swmmbuilder.maps.network_footprint` with ``pipes='gravity'``, so a force main
        across open land brings in no rural polygons). Pass the whole Network, not a catchment,
        as :func:`~swmmbuilder.allocate_population` wants. A shapely geometry needs ``crs``.

    level : str
        The census level: ``"DB"`` (dissemination blocks, default) or ``"DA"`` (dissemination
        areas) for Canada, or whatever the source takes.

    source : PolygonSource, optional
        Where the polygons come from. Default: :class:`CensusMapper` (which needs your API key).
        With a cache, it is only used to fetch; when reading a cache, a source given must be the
        one the cache was fetched from.

    cache : str or path, optional
        A GeoPackage. If it exists (and ``refresh`` is False), the polygons are read from it and
        nothing is fetched, so no network access or key is needed; it must be for the same level
        and cover the area, and, if ``source`` is given, have been fetched from it with the
        same parameters (for :class:`CensusMapper`, the dataset, regions, and resolution; for
        :class:`PolygonFile`, the file name and columns). Otherwise they are fetched and written to it, with the area and the
        provenance (source, level, dataset, regions, fetch date, versions; never the key).

    refresh : bool
        Fetch even if the cache exists, and replace it.

    margin : float
        Metres around the area (default 0). Polygons meeting the area enlarged by ``margin`` are
        kept. With ``allocate_population(sewered=...)``, use at least its
        ``max_nearest_distance``, so that flagged polygons a little away from the pipes are kept.

    crs : optional
        The CRS of a shapely ``area``, and of the result. Default for the result: the CRS in
        metres the area was measured in (the network's own CRS if it is projected).

    Returns
    -------
    geopandas.GeoDataFrame
        The polygons meeting the area, whole (not clipped), indexed by ID (``id``, str), with
        ``population`` (float; missing where the source has none) and ``parent`` (the enclosing
        area's ID, e.g. a block's dissemination area). ``.attrs["provenance"]`` holds the
        provenance. Pass it as ``allocate_population(net, polygons)``.

    Raises
    ------
    ValueError
        If the area is empty or has no CRS, the polygons are invalid (no CRS, no population,
        repeated IDs), a CensusMapper key is needed and none is given, or the cache does not
        match the request (naming the parameter that differs).

    Warns
    -----
    DataWarning
        If some polygons have no population (e.g. suppressed by the census); they count as 0.
        If a cache does not record a parameter of the source, so it cannot be checked.

    Examples
    --------
    The 2021 dissemination blocks near a network's pipes, cached so later builds are offline::

        from swmmbuilder.census import CensusMapper, census_polygons

        blocks = census_polygons(
            net, "DB", CensusMapper(api_key=my_key), cache="data/blocks_2021.gpkg", margin=200
        )
        alloc = allocate_population(net, blocks, max_distance=15)

    Polygons of any other source, from a file::

        zones = census_polygons(net, source=PolygonFile("zones.gpkg", id="ZONE", population="POP"))
    """
    if margin < 0:
        raise ValueError(f"margin must be at least 0, not {margin}.")
    region, work = _region(area, crs)
    if region.is_empty:
        raise ValueError("The area is empty (a network with no gravity pipes?).")
    if margin:
        region = region.buffer(margin)
    path = Path(cache) if cache is not None else None
    if path is not None and path.exists() and not refresh:
        polygons, provenance = _read_cache(path, level, source, region, work)
    else:
        source = CensusMapper() if source is None else source
        fetched = source.fetch(region, work, level)
        polygons = _standard(fetched, source)
        provenance = {
            "source": source.name,
            "level": level,
            **{k: str(v) for k, v in fetched.attrs.get("provenance", {}).items()},
            "fetched": datetime.now(UTC).isoformat(timespec="seconds"),
            "swmmbuilder": _version("swmmbuilder"),
        }
        if path is not None:
            _write_cache(path, polygons, region, work, provenance)
    polygons = polygons.to_crs(work)
    kept = polygons.iloc[sorted(polygons.sindex.query(region, predicate="intersects"))]
    kept = kept.to_crs(crs if crs is not None else work)
    missing = int(kept["population"].isna().sum())
    if missing:
        warnings.warn(
            f"{missing} of the {len(kept)} polygons have no population (suppressed or missing in "
            "the source); allocate_population counts them as 0",
            DataWarning,
            stacklevel=2,
        )
    kept.attrs["provenance"] = provenance
    return kept
