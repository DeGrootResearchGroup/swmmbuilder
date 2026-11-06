"""Census polygons from pluggable sources (swmmbuilder.census). No test uses the network: the
sources are fakes, and pycancensus is replaced by a fake module."""

import sys
import types
import warnings

import geopandas as gpd
import pandas as pd
import pytest
from pyproj import CRS as ProjCRS
from shapely.geometry import box
from test_population import CRS, X0, Y0, make_network

from swmmbuilder import DataWarning, allocate_population
from swmmbuilder.census import CensusMapper, PolygonFile, PolygonSource, census_polygons

# The test network's pipes run along y = 50 m, x = 10 to 250 m (from X0, Y0); its footprint
# reaches about 12 m from them. "near" and "edge" meet it; "far" is 150 m away, "away" 1 km.
BOXES = {
    "near": (0, 0, 100, 100, 80.0),
    "edge": (100, 55, 200, 100, 20.0),
    "far": (300, 200, 400, 300, 30.0),
    "away": (0, 1100, 100, 1200, 5.0),
}


def boxes(spec=BOXES, crs=CRS, **columns):
    return gpd.GeoDataFrame(
        {"population": [s[4] for s in spec.values()], **columns},
        geometry=[box(X0 + a, Y0 + b, X0 + c, Y0 + d) for a, b, c, d, _ in spec.values()],
        index=pd.Index(list(spec), name="id"),
        crs=crs,
    )


class FakeSource:
    """Returns ``data`` and records each fetch."""

    def __init__(self, data=None, name="fake", provenance=None):
        self.data = boxes() if data is None else data
        self.name = name
        self.calls = []
        if provenance:
            self.data.attrs["provenance"] = provenance

    def fetch(self, area, crs, level):
        self.calls.append((area, crs, level))
        return self.data.copy()


class NoFetch:
    name = "fake"

    def fetch(self, area, crs, level):
        raise AssertionError("fetched instead of reading the cache")


@pytest.fixture
def net():
    return make_network()


def test_keeps_the_polygons_meeting_the_footprint(net):
    source = FakeSource()
    polygons = census_polygons(net, source=source)
    assert list(polygons.index) == ["near", "edge"]
    assert polygons.index.name == "id"
    assert list(polygons.columns) == ["population", "parent", "geometry"]
    assert polygons.crs == ProjCRS(CRS)
    assert polygons.loc["near", "geometry"].area == pytest.approx(100 * 100)  # whole, not clipped
    area, crs, level = source.calls[0]
    assert level == "DB" and crs == ProjCRS(CRS)
    assert area.distance(box(X0 + 10, Y0 + 50, X0 + 250, Y0 + 50)) == 0
    assert isinstance(source, PolygonSource)


def test_margin_keeps_polygons_farther_away(net):
    polygons = census_polygons(net, source=FakeSource(), margin=200)
    assert list(polygons.index) == ["near", "edge", "far"]


def test_result_crs_and_source_crs(net):
    source = FakeSource(boxes().to_crs(4326))
    polygons = census_polygons(net, source=source, crs=4326)
    assert polygons.crs == ProjCRS(4326)
    assert list(polygons.index) == ["near", "edge"]


def test_areas_other_than_a_network():
    region = box(X0, Y0, X0 + 50, Y0 + 50)
    assert list(census_polygons(region, source=FakeSource(), crs=CRS).index) == ["near"]
    frame = gpd.GeoDataFrame(geometry=[region], crs=CRS).to_crs(4326)
    polygons = census_polygons(frame, source=FakeSource())
    assert list(polygons.index) == ["near"]
    assert polygons.crs.is_projected  # measured in UTM, not degrees
    with pytest.raises(ValueError, match="needs crs="):
        census_polygons(region, source=FakeSource())
    with pytest.raises(TypeError, match="area must be"):
        census_polygons("London", source=FakeSource())
    with pytest.raises(ValueError, match="margin"):
        census_polygons(region, source=FakeSource(), crs=CRS, margin=-1)


def test_population_is_numeric_and_missing_values_warn(net):
    data = boxes()
    data["population"] = data["population"].astype(object)
    data.loc["edge", "population"] = "x"  # suppressed
    with pytest.warns(DataWarning, match="1 of the 2 polygons have no population"):
        polygons = census_polygons(net, source=FakeSource(data))
    assert polygons["population"].dtype == float
    assert pd.isna(polygons.loc["edge", "population"])


def test_ids_and_parents_are_str(net):
    data = boxes(parent=[35390001, None, 3, 4])
    data.index = pd.Index([101, 102, 103, 104], name="id")
    polygons = census_polygons(net, source=FakeSource(data))
    assert list(polygons.index) == ["101", "102"]
    assert polygons.loc["101", "parent"] == "35390001"
    assert pd.isna(polygons.loc["102", "parent"])


@pytest.mark.parametrize(
    "change, message",
    [
        (
            lambda d: gpd.GeoDataFrame(d.drop(columns="geometry"), geometry=list(d.geometry)),
            "no CRS",
        ),
        (lambda d: d.rename(columns={"population": "pop"}), "no population column"),
        (lambda d: d.set_axis(["a", "a", "b", "c"]), "repeated IDs: a"),
    ],
)
def test_invalid_polygons_are_refused(net, change, message):
    with pytest.raises(ValueError, match=message):
        census_polygons(net, source=FakeSource(change(boxes())))


def test_cache_is_written_then_read_without_fetching(net, tmp_path):
    cache = tmp_path / "census" / "blocks.gpkg"
    first = census_polygons(
        net, source=FakeSource(provenance={"dataset": "CA21"}), cache=cache, margin=50
    )
    assert cache.exists()
    assert not list(tmp_path.rglob("*.partial*"))
    provenance = first.attrs["provenance"]
    assert provenance["source"] == "fake" and provenance["level"] == "DB"
    assert provenance["dataset"] == "CA21" and "fetched" in provenance
    # Read back: same polygons, no fetch, with or without the source named.
    for source in (NoFetch(), None):
        again = census_polygons(net, source=source, cache=cache, margin=50)
        pd.testing.assert_frame_equal(
            pd.DataFrame(again.drop(columns="geometry")),
            pd.DataFrame(first.drop(columns="geometry")),
        )
        assert again.geom_equals_exact(first.geometry, tolerance=1e-6).all()
        assert again.attrs["provenance"] == provenance
    # A smaller margin is covered by the cache too.
    assert len(census_polygons(net, source=NoFetch(), cache=cache)) == 2


def test_cache_must_match_the_request(net, tmp_path):
    cache = tmp_path / "blocks.gpkg"
    census_polygons(net, source=FakeSource(), cache=cache)
    with pytest.raises(ValueError, match="holds level 'DB', not 'DA'"):
        census_polygons(net, "DA", cache=cache)
    with pytest.raises(ValueError, match="fetched from 'fake', not 'other'"):
        census_polygons(net, source=FakeSource(name="other"), cache=cache)
    with pytest.raises(ValueError, match="covers a smaller area"):
        census_polygons(net, source=NoFetch(), cache=cache, margin=200)
    refreshed = census_polygons(net, source=FakeSource(), cache=cache, margin=200, refresh=True)
    assert "far" in refreshed.index
    (tmp_path / "junk.gpkg").write_text("not a geopackage")
    with pytest.raises(ValueError, match="not a census_polygons cache"):
        census_polygons(net, cache=tmp_path / "junk.gpkg")


def test_polygons_feed_allocate_population(net):
    polygons = census_polygons(net, source=FakeSource())
    alloc = allocate_population(net, polygons, max_distance=15)
    assert alloc.total == pytest.approx(100.0)


def test_polygon_file(net, tmp_path):
    path = tmp_path / "zones.gpkg"
    data = boxes().rename(columns={"population": "POP"}).assign(DA=["d1", "d1", "d2", "d3"])
    data.reset_index().rename(columns={"id": "ZONE"}).to_file(path, layer="zones")
    source = PolygonFile(path, id="ZONE", population="POP", parent="DA", layer="zones")
    polygons = census_polygons(net, source=source)
    assert list(polygons.index) == ["near", "edge"]
    assert list(polygons["parent"]) == ["d1", "d1"]
    assert polygons.attrs["provenance"]["file"] == "zones.gpkg"
    with pytest.raises(ValueError, match="no column 'PEOPLE'"):
        census_polygons(net, source=PolygonFile(path, id="ZONE", population="PEOPLE"))


# --- CensusMapper, with a fake pycancensus --------------------------------------------------

KEY = "CensusMapper_test0123456789abcdef"  # not a real key


def censusmapper_response():
    """What pycancensus.get_census returns for blocks with geo_format='geopandas': CensusMapper's
    short column names, in WGS 84. Made up, not census data."""
    data = boxes().to_crs(4326).reset_index()
    return gpd.GeoDataFrame(
        {
            "id": ["35390001001", "35390001002", "35390002001", "35390099001"],
            "pop": pd.array([80, 20, 30, None], dtype="Int64"),
            "dw": [30, 8, 12, 2],
            "hh": [29, 8, 11, 2],
            "rpid": ["35390001", "35390001", "35390002", "35390099"],
            "name": ["", "", "", ""],
        },
        geometry=data.geometry,
        crs=4326,
    )


@pytest.fixture
def pycancensus(monkeypatch):
    """A fake pycancensus module recording its calls."""
    fake = types.ModuleType("pycancensus")
    fake.calls = []

    def get_intersecting_geometries(**kwargs):
        fake.calls.append(("intersecting", kwargs))
        return {"CSD": [3539036]}

    def get_census(**kwargs):
        fake.calls.append(("census", kwargs))
        return censusmapper_response()

    fake.get_intersecting_geometries = get_intersecting_geometries
    fake.get_census = get_census
    monkeypatch.setitem(sys.modules, "pycancensus", fake)
    for variable in ("CENSUSMAPPER_API_KEY", "CANCENSUS_API_KEY"):
        monkeypatch.delenv(variable, raising=False)
    return fake


def test_censusmapper_finds_regions_and_maps_columns(net, pycancensus):
    polygons = census_polygons(net, source=CensusMapper(api_key=KEY))
    (kind, found), (_, census) = pycancensus.calls
    assert kind == "intersecting"
    assert found["level"] == "CSD" and found["dataset"] == "CA21" and found["api_key"] == KEY
    assert found["use_cache"] is False
    assert census["regions"] == {"CSD": ["3539036"]}
    assert census["level"] == "DB" and census["geo_format"] == "geopandas"
    assert census["api_key"] == KEY and census["use_cache"] is False
    assert list(polygons.index) == ["35390001001", "35390001002"]
    assert list(polygons["population"]) == [80.0, 20.0]
    assert list(polygons["parent"]) == ["35390001", "35390001"]
    provenance = polygons.attrs["provenance"]
    assert provenance["source"] == "CensusMapper" and provenance["dataset"] == "CA21"
    assert provenance["regions"] == '{"CSD": ["3539036"]}'


def test_censusmapper_given_regions_are_not_looked_up(net, pycancensus):
    source = CensusMapper(api_key=KEY, regions={"CSD": ["3539036"]}, dataset="CA16")
    census_polygons(net, "DA", source=source)
    [(kind, census)] = pycancensus.calls
    assert kind == "census" and census["level"] == "DA" and census["dataset"] == "CA16"


def test_censusmapper_no_regions_found(net, pycancensus):
    pycancensus.get_intersecting_geometries = lambda **kwargs: {"CSD": []}
    with pytest.raises(ValueError, match="found no CSD regions"):
        census_polygons(net, source=CensusMapper(api_key=KEY))


@pytest.mark.parametrize("variable", ["CENSUSMAPPER_API_KEY", "CANCENSUS_API_KEY"])
def test_censusmapper_key_from_the_environment(net, pycancensus, monkeypatch, variable):
    monkeypatch.setenv(variable, KEY)
    census_polygons(net, source=CensusMapper())
    assert all(call["api_key"] == KEY for _, call in pycancensus.calls)


def test_censusmapper_key_is_the_users_to_provide(net, pycancensus):
    with pytest.raises(ValueError, match="yours to provide: swmmbuilder has none"):
        census_polygons(net)  # the default source is CensusMapper
    assert pycancensus.calls == []


def test_censusmapper_cache_needs_no_key(net, pycancensus, tmp_path):
    cache = tmp_path / "blocks.gpkg"
    census_polygons(net, source=CensusMapper(api_key=KEY), cache=cache)
    pycancensus.calls.clear()
    assert len(census_polygons(net, cache=cache)) == 2  # no key anywhere, and no call
    assert pycancensus.calls == []


def test_censusmapper_key_is_never_kept(net, pycancensus, tmp_path, caplog):
    cache = tmp_path / "blocks.gpkg"
    source = CensusMapper(api_key=KEY)
    caplog.set_level("DEBUG")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        polygons = census_polygons(net, source=source, cache=cache, margin=2000)
    assert KEY not in repr(source)
    assert KEY not in str(polygons.attrs)
    assert KEY not in caplog.text
    assert not any(KEY in str(w.message) for w in caught)
    assert KEY.encode() not in cache.read_bytes()


def test_censusmapper_without_pycancensus(net, monkeypatch):
    monkeypatch.setitem(sys.modules, "pycancensus", None)
    with pytest.raises(ImportError, match=r"swmmbuilder\[census\]"):
        census_polygons(net, source=CensusMapper(api_key=KEY))
