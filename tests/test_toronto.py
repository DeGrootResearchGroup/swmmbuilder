"""The TORONTO_ON schema, on a sanitary catchment from the City of Toronto Open Data Portal
(tests/data/toronto_catchment.geojson; see tests/data/make_toronto_fixture.py)."""

import pytest
from conftest import TESTS, quiet, read_sections, run_swmm, write_inp

from swmmbuilder import TORONTO_ON, Network, SubNetwork

TORONTO = TESTS / "data" / "toronto_catchment.geojson"
# The catchment drains to MH5130118747, and pipe SL963 takes it on to OUTFALL.
OUTFALL = "MH5128918783"


@pytest.fixture
def toronto():
    return quiet(Network, str(TORONTO), schema=TORONTO_ON)


@pytest.fixture
def catchment(toronto):
    sub = quiet(
        SubNetwork,
        {
            "network_data": toronto,
            "starting_node": OUTFALL,
            "outfall_node": [OUTFALL],
            "depth": 1,
            "outfall_node_type": ["FREE"],
            "scale": 1.0,
        },
    )
    quiet(sub.catchment, OUTFALL)
    return sub


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
    network = quiet(Network, str(path), schema=TORONTO_ON)
    pipe = network.full_data.iloc[0]
    assert pipe["length"] == pytest.approx(pipe.geometry.length)


def test_catchment_is_whole_fixture(catchment):
    assert catchment.subgraph.number_of_edges() == 31


def test_toronto_inp_runs_in_swmm(catchment, tmp_path):
    # SL939 slopes uphill in the GIS inverts; KINWAVE rejects it, so run DYNWAVE.
    with pytest.warns(UserWarning, match=r"adverse slope.*SL939"):
        path = write_inp(catchment, tmp_path, "toronto.inp", routing="DYNWAVE")
    sections = read_sections(path)
    assert len(sections["CONDUITS"]) == 31
    assert [r[0] for r in sections["OUTFALLS"]] == [OUTFALL]
    x, y = map(float, sections["COORDINATES"][0][1:])
    assert 300_000 < x < 340_000 and 4_830_000 < y < 4_870_000  # MTM zone 10, metres
    assert run_swmm(path, tmp_path) == []
