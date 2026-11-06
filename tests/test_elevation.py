"""Rims estimated from a ground-elevation model (DEM) (#112), on small synthetic rasters."""

import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from conftest import DENDRITIC, nodes_of
from pyproj import Transformer
from shapely import Point

from swmmbuilder import DataWarning, OptionWarning, swmm
from swmmbuilder.elevation import RimEstimate, sample_dem

CRS = "EPSG:2958"  # the fixture's: NAD83(CSRS) / UTM zone 17N
X0, Y0 = 475_300.0, 4_758_600.0  # lower left corner of the rasters, around the fixture's nodes
WIDTH, HEIGHT = 700, 750  # 1 m cells
# Recorded rims in the tests are this far above the DEM.
LID = 0.1


def ground(x, y):
    """The synthetic ground: a plane rising 1 m per km east and 2 m per km north from 110 m."""
    return 110.0 + 0.001 * (np.asarray(x) - X0) + 0.002 * (np.asarray(y) - Y0)


def grid(x0=X0, y0=Y0, width=WIDTH, height=HEIGHT, cell=1.0):
    """(x, y) of the cell centres of a north-up raster whose lower left corner is (x0, y0)."""
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    return x0 + (cols + 0.5) * cell, y0 + height * cell - (rows + 0.5) * cell


def write_dem(path, values, x0=X0, y0=Y0, cell=1.0, crs=CRS, nodata=None):
    """Write ``values`` (rows from north to south) as a single-band GeoTIFF."""
    rasterio = pytest.importorskip("rasterio")
    from rasterio.transform import from_origin

    height, width = values.shape
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        width=width,
        height=height,
        count=1,
        dtype="float32",
        crs=crs,
        transform=from_origin(x0, y0 + height * cell, cell, cell),
        nodata=nodata,
    ) as dst:
        dst.write(values.astype("float32"), 1)
    return path


@pytest.fixture
def dem(tmp_path):
    pytest.importorskip("rasterio")
    return write_dem(tmp_path / "dem.tif", ground(*grid()))


def record_rims(network, rims):
    """Give nodes recorded rims, as if read from a node layer."""
    for node, rim in rims.items():
        network.nodes.loc[node, "rim_elevation"] = rim
        network.graph.nodes[node]["rim_elevation"] = rim


def ground_at(network, node):
    point = network.nodes.geometry[node]
    return float(ground(point.x, point.y))


@pytest.fixture
def recorded(network):
    """A1, A2, B1, and J1 have rims LID above the ground; the other nodes have none."""
    rims = {n: round(ground_at(network, n) + LID, 3) for n in ("A1", "A2", "B1", "J1")}
    record_rims(network, rims)
    return rims


# --- sampling ---------------------------------------------------------------------------------


def test_sample_the_cell_at_each_point(network, dem):
    points = network.nodes.geometry
    sampled = sample_dem(dem, points)
    expected = ground(points.x, points.y)
    # The cell centre is within 0.5 m of the point in each direction.
    assert np.allclose(sampled, expected, atol=0.0016)


def test_a_window_softens_a_spike(network, tmp_path):
    point = network.nodes.geometry["A1"]
    xs, ys = grid()
    values = ground(xs, ys)
    curb = (np.abs(xs - point.x) <= 0.6) & (np.abs(ys - point.y) <= 0.6)
    values[curb] += 0.5  # a curb or wall at the node
    path = write_dem(tmp_path / "curb.tif", values)
    points = network.nodes.geometry[["A1"]]
    assert sample_dem(path, points)["A1"] == pytest.approx(ground_at(network, "A1") + 0.5, abs=2e-3)
    softened = sample_dem(path, points, window=3.0)["A1"]
    assert softened == pytest.approx(ground_at(network, "A1"), abs=5e-3)
    mean = sample_dem(path, points, window=3.0, statistic="mean")["A1"]
    assert softened < mean < ground_at(network, "A1") + 0.5


def test_a_dem_in_another_crs_is_sampled_at_the_transformed_points(network, tmp_path):
    # A DEM in geographic coordinates (EPSG:4326) whose cells hold the ground at their centres.
    to_wgs84 = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
    lon0, lat0 = to_wgs84.transform(X0 - 50, Y0 - 50)
    lon1, lat1 = to_wgs84.transform(X0 + WIDTH + 50, Y0 + HEIGHT + 50)
    cell = 1e-5  # about 1 m
    width, height = int((lon1 - lon0) / cell) + 1, int((lat1 - lat0) / cell) + 1
    lons, lats = grid(lon0, lat0, width, height, cell)
    xs, ys = Transformer.from_crs("EPSG:4326", CRS, always_xy=True).transform(lons, lats)
    path = write_dem(tmp_path / "wgs84.tif", ground(xs, ys), lon0, lat0, cell, "EPSG:4326")

    points = network.nodes.geometry
    assert np.allclose(sample_dem(path, points), ground(points.x, points.y), atol=0.002)
    # A window in metres on a geographic raster.
    windowed = sample_dem(path, points, window=2.0)
    assert np.allclose(windowed, ground(points.x, points.y), atol=0.002)


def test_tiles_give_the_same_values_as_one_raster(network, dem, tmp_path):
    xs, ys = grid()
    values = ground(xs, ys)
    half = WIDTH // 2
    west = write_dem(tmp_path / "west.tif", values[:, :half])
    east = write_dem(tmp_path / "east.tif", values[:, half:], x0=X0 + half)
    points = network.nodes.geometry
    whole = sample_dem(dem, points, window=2.0)
    assert np.allclose(sample_dem([west, east], points, window=2.0), whole, atol=1e-6)


def test_points_outside_the_dem_or_on_nodata_get_nan(network, tmp_path):
    xs, ys = grid()
    values = ground(xs, ys)
    point = network.nodes.geometry["A1"]
    values[(np.abs(xs - point.x) <= 3) & (np.abs(ys - point.y) <= 3)] = -9999.0
    path = write_dem(tmp_path / "holes.tif", values, nodata=-9999.0)
    points = gpd.GeoSeries(
        [network.nodes.geometry["A1"], Point(0, 0), network.nodes.geometry["A2"], None],
        index=["A1", "far", "A2", "none"],
        crs=CRS,
    )
    sampled = sample_dem(path, points)
    assert sampled[["A1", "far", "none"]].isna().all()
    assert sampled["A2"] == pytest.approx(ground_at(network, "A2"), abs=2e-3)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"window": -1.0}, "window must be 0 or more"),
        ({"statistic": "max"}, "statistic must be one of"),
    ],
)
def test_invalid_sampling_options(network, dem, kwargs, message):
    with pytest.raises(ValueError, match=message):
        sample_dem(dem, network.nodes.geometry, **kwargs)


def test_a_dem_with_no_crs_is_refused(network, tmp_path):
    path = write_dem(tmp_path / "nocrs.tif", ground(*grid()), crs=None)
    with pytest.raises(ValueError, match="must both be known"):
        sample_dem(path, network.nodes.geometry)


def test_rasterio_is_needed(network, monkeypatch):
    monkeypatch.setitem(sys.modules, "rasterio", None)
    with pytest.raises(ImportError, match=r"swmmbuilder\[elevation\]"):
        sample_dem("dem.tif", network.nodes.geometry)


# --- validation -------------------------------------------------------------------------------


def test_the_validation_report(network, recorded, dem):
    record_rims(network, {"J1": recorded["J1"] + 1.5})  # a probably wrong recorded rim
    estimate = network.estimate_rims(dem)
    assert isinstance(estimate, RimEstimate)
    v = estimate.validation
    assert v["count"] == 4
    assert v["median"] == pytest.approx(LID, abs=2e-3)
    assert v["max_abs"] == pytest.approx(LID + 1.5, abs=2e-3)
    assert list(estimate.largest_differences(1).index) == ["J1"]
    assert set(estimate.compared.index) == set(recorded)
    assert "median +0.10" in estimate.summary()


def test_estimating_changes_nothing(network, recorded, dem):
    version, nodes = network.version, network.nodes.copy()
    network.estimate_rims(dem)
    assert network.version == version
    pd.testing.assert_frame_equal(network.nodes, nodes)
    assert "rim_elevation" not in network.graph.nodes["C1"]


# --- estimates --------------------------------------------------------------------------------


def test_only_nodes_with_no_usable_rim_are_estimated(network, recorded, dem):
    estimate = network.estimate_rims(dem)
    missing = set(network.nodes.index) - set(recorded)
    assert set(estimate.estimated.index) == missing
    rims = estimate.rims
    assert set(rims) == missing
    for node in missing:
        assert rims[node] == pytest.approx(ground_at(network, node), abs=2e-3)
    assert rims.source == "DEM dem.tif"


def test_median_offset(network, recorded, dem):
    estimate = network.estimate_rims(dem, offset="median")
    assert estimate.offset == pytest.approx(LID, abs=2e-3)
    assert estimate.rims["C1"] == pytest.approx(ground_at(network, "C1") + LID, abs=3e-3)


def test_a_rim_below_the_invert_is_not_usable(network, recorded, dem):
    invert = network.nodes.loc["C1", "invert"]
    record_rims(network, {"C1": invert - 1.0})
    estimate = network.estimate_rims(dem)
    assert estimate.samples.loc["C1", "status"] == "estimated"
    assert estimate.rims["C1"] == pytest.approx(ground_at(network, "C1"), abs=2e-3)


def test_a_dem_rim_at_or_below_the_invert_is_refused(network, recorded, dem):
    network.nodes.loc["C1", "invert"] = ground_at(network, "C1") + 1.0  # a structure under a bridge
    with pytest.warns(DataWarning, match="1 rims from the DEM are at or below the invert"):
        estimate = network.estimate_rims(dem)
    assert estimate.samples.loc["C1", "status"] == "dem_not_above_invert"
    assert "C1" not in estimate.rims


def test_a_node_with_no_dem_value_is_left(network, recorded, tmp_path):
    # A DEM covering only the western half: C1 (east) has no value.
    xs, ys = grid(width=300)
    path = write_dem(tmp_path / "west.tif", ground(xs, ys))
    estimate = network.estimate_rims(path)
    assert estimate.samples.loc["C1", "status"] == "no_dem_value"
    assert "C1" not in estimate.rims


def test_some_nodes_only(network, recorded, dem):
    estimate = network.estimate_rims(dem, nodes=["C1", "A1"])
    assert sorted(estimate.samples.index) == ["A1", "C1"]
    assert list(estimate.rims) == ["C1"]
    with pytest.raises(ValueError, match="not in the network"):
        network.estimate_rims(dem, nodes=["nowhere"])


# --- the rims option --------------------------------------------------------------------------


def no_rim_rows(report):
    return report[report["gap"] == "no_rim"].set_index("id")


def test_every_node_with_no_usable_rim_is_a_gap(network, recorded):
    invert = network.nodes.loc["C1", "invert"]
    record_rims(network, {"C1": invert - 1.0})
    sub = network.subnetwork(nodes_of(DENDRITIC))
    rows = no_rim_rows(sub.data_gaps())
    assert set(rows.index) == {"B2", "C1", "J2", "OUT"}
    assert rows["filled_from"].isna().all()
    assert rows.loc["B2", "detail"] == "no rim recorded"
    assert (
        rows.loc["C1", "detail"] == f"recorded rim {invert - 1.0:g} not above the invert {invert:g}"
    )
    assert sub.validate().data_gaps["no_rim"] == 4


def test_estimated_rims_are_flagged_in_the_data_gaps(network, recorded, dem):
    estimate = network.estimate_rims(dem, window=1.0, offset=0.05)
    sub = network.subnetwork(nodes_of(DENDRITIC))
    rows = no_rim_rows(sub.data_gaps(rims=estimate.rims))
    assert set(rows.index) == {"B2", "C1", "J2", "OUT"}
    assert set(rows["filled_from"]) == {"estimate_rims"}
    assert set(rows["source"]) == {"DEM dem.tif, median within 1 m, offset +0.050 m"}
    assert rows.loc["C1", "value"] == estimate.rims["C1"]
    assert rows.loc["C1", "detail"] == "no rim recorded"


def test_a_plain_mapping_is_filled_from_rims(network, recorded):
    sub = network.subnetwork(nodes_of(DENDRITIC))
    rows = no_rim_rows(sub.data_gaps(rims={"C1": 112.0}))
    assert rows.loc["C1", "filled_from"] == "rims" and rows.loc["C1", "value"] == 112.0
    assert rows["filled_from"].notna().sum() == 1


def test_the_model_takes_rim_minus_invert(network, recorded, dem):
    estimate = network.estimate_rims(dem)
    sub = network.subnetwork(nodes_of(DENDRITIC))
    model = swmm.build_model(sub, ["OUT"], rims=estimate.rims)
    depths = {j.name: j.max_depth for j in model.junctions}
    invert = network.nodes.loc["C1", "invert"]
    assert depths["C1"] == pytest.approx(round(estimate.rims["C1"] - invert, 4))
    assert depths["A1"] == pytest.approx(
        round(recorded["A1"] - network.nodes.loc["A1", "invert"], 4)
    )
    assert model.data_gaps["filled_from"].eq("estimate_rims").sum() == 4  # B2, C1, J2, OUT


def test_estimated_rims_count_for_the_median_depth(network, recorded, dem):
    sub = network.subnetwork(nodes_of(DENDRITIC))
    before = swmm.build_model(sub, ["OUT"], max_depth="median")
    assert {d.id for d in before.derived} == {"B2", "C1", "J2"}
    estimate = network.estimate_rims(dem)
    after = swmm.build_model(sub, ["OUT"], max_depth="median", rims=estimate.rims)
    assert after.derived == []
    inverts = swmm.node_inverts(sub)
    junctions = [n for n in sub.nodes if n != "OUT"]
    rims = {**recorded, **estimate.rims}
    expected = np.median([round(rims[n] - inverts[n], 4) for n in junctions])
    written = {j.name: j.max_depth for j in after.junctions}
    assert np.median([written[n] for n in junctions]) == pytest.approx(expected)


def test_a_given_rim_at_or_below_the_invert_is_refused_at_build_time(network, recorded):
    sub = network.subnetwork(nodes_of(DENDRITIC))
    invert = network.nodes.loc["C1", "invert"]
    with pytest.warns(DataWarning, match="rims at or below the node invert, not used: C1"):
        model = swmm.build_model(sub, ["OUT"], max_depth=1.5, rims={"C1": invert})
    assert {j.name: j.max_depth for j in model.junctions}["C1"] == 1.5
    rows = no_rim_rows(model.data_gaps)
    assert pd.isna(rows.loc["C1", "filled_from"])


def test_rims_for_unknown_nodes_are_named(network, recorded):
    sub = network.subnetwork(nodes_of(DENDRITIC))
    with pytest.warns(OptionWarning, match="rims: nodes not in the network, ignored: nowhere"):
        sub.data_gaps(rims={"nowhere": 100.0})


@pytest.mark.parametrize(
    "offset, message",
    [("mean", "offset must be a number or 'median'"), (float("nan"), "finite number")],
)
def test_invalid_offset(network, dem, offset, message):
    with pytest.raises(ValueError, match=message):
        network.estimate_rims(dem, offset=offset)


def test_median_offset_needs_recorded_rims(network, dem):
    with pytest.raises(ValueError, match="no node has both a usable recorded rim"):
        network.estimate_rims(dem, offset="median")
