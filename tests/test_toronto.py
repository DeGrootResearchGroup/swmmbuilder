"""The TORONTO_ON schema, on a sanitary catchment from the City of Toronto Open Data Portal
(tests/data/toronto_catchment.geojson; see tests/data/make_toronto_fixture.py)."""

import gzip
import json

import geopandas as gpd
import pytest
from conftest import TESTS, read_sections, run_swmm, write_inp

from swmmbuilder import TORONTO_ON, Network, maps

TORONTO = TESTS / "data" / "toronto_catchment.geojson"
TORONTO_LARGE = TESTS / "data" / "toronto_large_catchment.geojson.gz"
TORONTO_LARGE_OUTLET = "MH3948903491"
# The catchment drains to MH5130118747, and pipe SL963 takes it on to OUTFALL.
OUTFALL = "MH5128918783"


@pytest.fixture
def toronto():
    return Network(str(TORONTO), schema=TORONTO_ON)


@pytest.fixture
def catchment(toronto):
    return toronto.catchment(OUTFALL)


def test_load(toronto):
    data = toronto.full_data
    assert len(data) == 31 and toronto.graph.number_of_edges() == 31
    assert data.crs == "EPSG:2952"
    assert set(data["flow_type"]) == {"sanitary"} and set(data["shape"]) == {"CIRCULAR"}
    assert data["diameter"].between(0.2, 0.4).all()
    assert data[["us_invert", "ds_invert"]].notna().all().all()
    assert "Sewer Gravity Location Description" in data.columns


def test_structure_types(tmp_path):
    """Structure types use London's codes: SL is gravity, OVF an overflow pipe (#154); other codes
    pass through."""
    assert set(Network(str(TORONTO), schema=TORONTO_ON).full_data["structure_type"]) == {"gravity"}
    data = gpd.read_file(TORONTO)
    data.loc[0, "Sewer Gravity Structure Type"] = "OVF"
    data.loc[1, "Sewer Gravity Structure Type"] = "SP"
    path = tmp_path / "toronto.geojson"
    data.to_file(path, driver="GeoJSON")
    types = Network(str(path), schema=TORONTO_ON).full_data["structure_type"]
    assert types.value_counts().to_dict() == {"gravity": 29, "overflow": 1, "SP": 1}
    network = Network(str(path), schema=TORONTO_ON, include_overflows=False)
    assert network.removed["structure_type"] == {"overflow": 1}
    assert len(network.full_data) == 30


def test_lengths_match_reprojected_geometry(toronto):
    """The measured lengths agree with the lengths of the reprojected geometry (to 10 cm)."""
    data = toronto.full_data
    assert (data.geometry.length - data["length"]).abs().max() < 0.1


def test_missing_length_from_geometry(tmp_path):
    data = gpd.read_file(TORONTO)
    data.loc[0, "Sewer Gravity Measured Length"] = None
    path = tmp_path / "toronto.geojson"
    data.to_file(path, driver="GeoJSON")
    network = Network(str(path), schema=TORONTO_ON)
    pipe = network.full_data.iloc[0]
    assert pipe["length"] == pytest.approx(pipe.geometry.length)


def test_catchment_is_whole_fixture(catchment):
    assert catchment.graph.number_of_edges() == 31


def test_large_catchment():
    """The 2,692-pipe catchment (with trunk sewers) used by the map gallery at full scale; see
    tests/data/make_toronto_large_fixture.py."""
    data = gpd.read_file(gzip.decompress(TORONTO_LARGE.read_bytes()))
    with pytest.warns(UserWarning, match="Shape codes not in schema"):
        network = Network(data, schema=TORONTO_ON)
    sub = network.catchment(TORONTO_LARGE_OUTLET)
    assert sub.graph.number_of_edges() == len(data) == 2692
    assert (data["Sewer Gravity Trunk Sewer"] == "Yes").sum() == 159
    xmin, ymin, xmax, ymax = sub.to_geodataframe("pipes").total_bounds
    assert 8000 < max(xmax - xmin, ymax - ymin) < 10000  # metres, in EPSG:2952


def test_toronto_inp_runs_in_swmm(catchment, tmp_path):
    # SL939 slopes uphill in the GIS inverts; KINWAVE rejects it, so run DYNWAVE.
    with pytest.warns(UserWarning, match=r"adverse slope.*SL939"):
        path = write_inp(catchment, tmp_path, "toronto.inp", "DYNWAVE", outfalls=(OUTFALL,))
    sections = read_sections(path)
    assert len(sections["CONDUITS"]) == 31
    assert [r[0] for r in sections["OUTFALLS"]] == [OUTFALL]
    x, y = map(float, sections["COORDINATES"][0][1:])
    assert 300_000 < x < 340_000 and 4_830_000 < y < 4_870_000  # MTM zone 10, metres
    assert run_swmm(path, tmp_path) == []


def raw_upstream_ends():
    """{node: (lon, lat)} of the upstream end of each pipe leaving it, from the WGS 84 fixture."""
    ends = {}
    for feature in json.loads(TORONTO.read_text())["features"]:
        node = feature["properties"]["Sewer Gravity Upstream Maintenance Hole"]
        ends.setdefault(node, tuple(feature["geometry"]["coordinates"][0][0]))
    return ends


def test_map_positions_are_reprojected_to_wgs84(catchment):
    """The network is in EPSG:2952 (metres); the map gets back the WGS 84 of the GeoJSON."""
    positions = maps.lonlat_positions(catchment)
    ends = raw_upstream_ends()
    assert set(ends) < set(positions)
    for node, lonlat in ends.items():
        assert positions[node] == pytest.approx(lonlat, abs=1e-7), node


def test_toronto_gmap(catchment, tmp_path, monkeypatch):
    markers = {}

    class FakePlotter:
        def __init__(self, lat, lng, zoom):
            self.centre = (lat, lng)

        def marker(self, lat, lng, color, label):
            markers[label] = (lat, lng)

        def plot(self, lats, lngs, color, edge_width):
            pass

        def draw(self, file_name):
            pass

    monkeypatch.setattr("gmplot.GoogleMapPlotter", FakePlotter)
    maps.generate_gmap(catchment, str(tmp_path / "map.html"))
    assert len(markers) == len(catchment.nodes)
    for lat, lng in markers.values():
        assert 43.5 < lat < 43.9 and -79.7 < lng < -79.1  # Toronto
