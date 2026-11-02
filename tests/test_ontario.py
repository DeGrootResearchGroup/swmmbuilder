"""The BARRIE_ON, KITCHENER_ON, PETERBOROUGH_ON, SARNIA_ON, WATERLOO_ON and WINDSOR_ON schemas, on
a sanitary catchment from each city's open data (tests/data/<city>_catchment.geojson; see
tests/data/make_ontario_fixtures.py)."""

import json
import warnings
from dataclasses import dataclass

import geopandas as gpd
import pytest
from conftest import TESTS, read_sections, run_swmm, write_inp

import swmmbuilder
from swmmbuilder import Network, maps


@dataclass
class City:
    name: str
    pipes: int
    crs: str
    upstream: str  # source column of the upstream node
    outfall: object  # downstream node of the pipe leaving the catchment's outlet
    lonlat: tuple  # (min lon, max lon, min lat, max lat) around the catchment
    routing: str | None = None  # None: the template's KINWAVE; DYNWAVE where flow splits

    @property
    def schema(self):
        return getattr(swmmbuilder, f"{self.name.upper()}_ON")

    @property
    def path(self):
        return TESTS / "data" / f"{self.name}_catchment.geojson"


CITIES = [
    City("barrie", 14, "EPSG:26917", "FROM_ID", "120264435", (-79.8, -79.6, 44.3, 44.4)),
    City(
        "kitchener",
        19,
        "EPSG:26917",
        "UP_SANMANHOLEID",
        2151735,
        (-80.6, -80.4, 43.4, 43.5),
        "DYNWAVE",
    ),
    City("peterborough", 36, "EPSG:2958", "UP_MH_ID", "83242", (-78.4, -78.3, 44.2, 44.4)),
    City("sarnia", 18, "EPSG:26917", "MH_Upstream", "MH7922", (-82.4, -82.3, 42.9, 43.0)),
    City("waterloo", 46, "EPSG:26917", "FROMMH", "W2-11C", (-80.6, -80.5, 43.4, 43.5), "DYNWAVE"),
    City("windsor", 14, "EPSG:26917", "UPSTREAM_MH", "1S4101", (-83.0, -82.9, 42.3, 42.4)),
]


@pytest.fixture(params=CITIES, ids=lambda city: city.name)
def city(request):
    return request.param


@pytest.fixture
def network(city):
    return Network(str(city.path), schema=city.schema)


def test_load(city, network):
    data = network.full_data
    assert len(data) == city.pipes and network.graph.number_of_edges() == city.pipes
    assert data.crs == city.crs
    assert set(data["shape"]) == {"CIRCULAR"}
    assert data["diameter"].between(0.15, 0.6).all()
    assert data[["from_node", "to_node", "us_invert", "ds_invert"]].notna().all().all()
    assert data["flow_type"].isin({"sanitary"}).all() or data["flow_type"].isna().all()


def test_lengths_match_geometry(network):
    """Mapped lengths are as measured, so agree with the reprojected geometry only roughly."""
    data = network.full_data
    assert ((data.geometry.length - data["length"]).abs() < 0.1 * data["length"]).all()


def test_catchment_is_whole_fixture(city, network):
    assert network.catchment(city.outfall).graph.number_of_edges() == city.pipes


def test_inp_runs_in_swmm(city, network, tmp_path):
    sub = network.catchment(city.outfall)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        path = write_inp(sub, tmp_path, f"{city.name}.inp", city.routing, outfalls=(city.outfall,))
    sections = read_sections(path)
    assert len(sections["CONDUITS"]) == city.pipes
    assert [r[0] for r in sections["OUTFALLS"]] == [str(city.outfall)]
    assert run_swmm(path, tmp_path) == []


@pytest.mark.parametrize("city", [c for c in CITIES if c.routing], ids=lambda c: c.name)
def test_flow_split_needs_dynwave(city, network, tmp_path):
    with pytest.raises(ValueError, match="more than one outlet pipe"):
        write_inp(network.catchment(city.outfall), tmp_path, outfalls=(city.outfall,))


def test_map_positions_are_reprojected_to_wgs84(city, network):
    """The network is in a projected CRS; the map gets back the WGS 84 of the GeoJSON."""
    ends = {}
    for feature in json.loads(city.path.read_text())["features"]:
        node = feature["properties"][city.upstream]
        ends.setdefault(node, tuple(feature["geometry"]["coordinates"][0]))
    positions = maps.lonlat_positions(network.catchment(city.outfall))
    assert set(ends) < set(positions)
    for node, lonlat in ends.items():
        assert positions[node] == pytest.approx(lonlat, abs=1e-7), node
    west, east, south, north = city.lonlat
    for lon, lat in positions.values():
        assert west < lon < east and south < lat < north


def edited(city, tmp_path, rows, **values):
    """The fixture with ``values`` ({column: value}) set on ``rows``, written to tmp_path."""
    data = gpd.read_file(city.path)
    for column, value in values.items():
        data.loc[rows, column] = value
    path = tmp_path / city.path.name
    data.to_file(path, driver="GeoJSON")
    return str(path)


def by_name(name):
    return next(c for c in CITIES if c.name == name)


def test_waterloo_placeholder_invert(tmp_path):
    waterloo = by_name("waterloo")
    path = edited(waterloo, tmp_path, [0], UPSTREAMINVERT=111.111)
    data = Network(path, schema=waterloo.schema).full_data
    assert data["us_invert"].isna().tolist() == [True] + [False] * (waterloo.pipes - 1)


def test_kitchener_missing_node_and_force_main(tmp_path):
    kitchener = by_name("kitchener")
    path = edited(kitchener, tmp_path, [0], UP_SANMANHOLEID=-1, CATEGORY="FORCEMAIN")
    network = Network(path, schema=kitchener.schema, flow_types=None)
    assert network.full_data["from_node"].isna().sum() == 1
    network = Network(path, schema=kitchener.schema)
    assert network.removed["flow_type"] == {"FORCEMAIN": 1}


def test_barrie_force_main_dropped(tmp_path):
    barrie = by_name("barrie")
    path = edited(barrie, tmp_path, [0, 1], TYPE="FORCE")
    assert Network(path, schema=barrie.schema).removed["flow_type"] == {"FORCE": 2}


def test_sarnia_non_circular_size(tmp_path):
    """Sarnia gives non-circular sizes as "W x H": those pipes get no diameter."""
    sarnia = by_name("sarnia")
    path = edited(sarnia, tmp_path, [0], Diam_m="600 x 750")
    data = Network(path, schema=sarnia.schema).full_data
    assert data["diameter"].isna().sum() == 1
    assert set(data["shape"]) == {"CIRCULAR"} and data["shape_code"].isna().all()


def test_windsor_ellipse(tmp_path):
    windsor = by_name("windsor")
    path = edited(windsor, tmp_path, [0], PIPE_SHAPE="EL", PIPE_SIZE=600, PIPE_HGT=900)
    pipe = Network(path, schema=windsor.schema).full_data.iloc[0]
    assert pipe["shape"] == "VERT_ELLIPSE"
    assert (pipe["height"], pipe["width"]) == pytest.approx((0.9, 0.6))
