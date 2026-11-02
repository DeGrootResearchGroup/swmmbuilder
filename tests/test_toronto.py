"""The TORONTO_ON schema, on a sanitary catchment from the City of Toronto Open Data Portal
(tests/data/toronto_catchment.geojson; see tests/data/make_toronto_fixture.py)."""

import json

import pytest
from conftest import TESTS, read_sections, run_swmm, write_inp

from sewernetworkx import TORONTO_ON, Network, maps

TORONTO = TESTS / "data" / "toronto_catchment.geojson"
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


def test_lengths_match_reprojected_geometry(toronto):
    """The measured lengths agree with the lengths of the reprojected geometry (to 10 cm)."""
    data = toronto.full_data
    assert (data.geometry.length - data["length"]).abs().max() < 0.1


def test_missing_length_from_geometry(tmp_path):
    import geopandas as gpd

    data = gpd.read_file(TORONTO)
    data.loc[0, "Sewer Gravity Measured Length"] = None
    path = tmp_path / "toronto.geojson"
    data.to_file(path, driver="GeoJSON")
    network = Network(str(path), schema=TORONTO_ON)
    pipe = network.full_data.iloc[0]
    assert pipe["length"] == pytest.approx(pipe.geometry.length)


def test_catchment_is_whole_fixture(catchment):
    assert catchment.graph.number_of_edges() == 31


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

    monkeypatch.setattr(maps.gmplot, "GoogleMapPlotter", FakePlotter)
    maps.generate_gmap(catchment, str(tmp_path / "map.html"))
    assert len(markers) == len(catchment.nodes)
    for lat, lng in markers.values():
        assert 43.5 < lat < 43.9 and -79.7 < lng < -79.1  # Toronto
