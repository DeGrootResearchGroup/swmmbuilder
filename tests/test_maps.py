"""Map plotting (#26): to_geodataframe, static plots, group polygons, basemaps and export."""

import json
import sys
import types

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from conftest import DENDRITIC, FIXTURE, nodes_of, subnetwork
from matplotlib.collections import LineCollection, PatchCollection
from pyproj import CRS
from shapely import LineString, Point

from swmmbuilder import Network, maps

matplotlib.use("Agg")

GROUPS = {"A": ["A1", "A2"], "B": ["B1", "B2"], "C": ["C1"], "trunk": ["J1", "J2", "OUT"]}


@pytest.fixture
def ax():
    fig, ax = plt.subplots()
    yield ax
    plt.close(fig)


@pytest.fixture
def sub(network):
    return subnetwork(network, DENDRITIC)


# ---- to_geodataframe ------------------------------------------------------------------------


def test_pipes_geodataframe_joins_values(network):
    values = pd.Series({1000: 2.5, 1001: 4.0, 99999: 1.0}, name="conc")
    pipes = network.to_geodataframe("pipes", values)
    assert pipes.index.name == "pipe_id"
    assert len(pipes) == len(network.full_data)
    assert pipes.loc[1000, "conc"] == 2.5
    assert pipes["conc"].notna().sum() == 2  # 99999 is not a pipe
    assert {"from_node", "to_node", "length", "geometry"} <= set(pipes.columns)
    assert pipes.crs == network.full_data.crs


def test_unnamed_values_are_called_value(network):
    assert "value" in network.to_geodataframe("pipes", pd.Series({1000: 1.0})).columns
    assert "value" in network.to_geodataframe("nodes", {"A1": 1.0}).columns


def test_nodes_geodataframe(network, sub):
    nodes = sub.to_geodataframe("nodes")
    assert nodes.index.name == "node_id"
    assert list(nodes.index) == sorted(nodes_of(DENDRITIC))
    positions = sub.node_positions()
    for node, point in nodes.geometry.items():
        assert (point.x, point.y) == pytest.approx(positions[node])
    assert set(network.to_geodataframe("nodes").index) == set(network.graph.nodes)


def test_subnetwork_pipes_geodataframe(sub):
    pipes = sub.to_geodataframe("pipes")
    assert set(pipes.index) == {key for _, _, key in sub.graph.edges(keys=True)}


def test_bad_kind(network):
    with pytest.raises(ValueError, match="kind"):
        network.to_geodataframe("manholes")


# ---- plot_crs -------------------------------------------------------------------------------


def test_plot_crs_is_utm_for_lonlat(network, sub):
    assert maps.plot_crs(network) == CRS.from_epsg(32617)  # London, ON: UTM 17N
    assert maps.plot_crs(sub) == maps.plot_crs(network)
    assert maps.plot_crs(network, "EPSG:3857") == CRS.from_epsg(3857)


def test_plot_crs_keeps_projected_crs():
    data = gpd.read_file(FIXTURE).to_crs("EPSG:2952")
    assert maps.plot_crs(Network(data)) == CRS.from_epsg(2952)


def test_plot_crs_without_crs_is_utm():
    data = gpd.read_file(FIXTURE).set_crs(None, allow_override=True)
    assert maps.plot_crs(Network(data)) == CRS.from_epsg(32617)


# ---- plot_pipes -----------------------------------------------------------------------------


def test_plot_pipes_colours_by_value(ax, network):
    pipes = network.to_geodataframe("pipes")
    values = pd.Series(np.arange(len(pipes), dtype=float), index=pipes.index)
    values.iloc[0] = np.nan
    lines = maps.plot_pipes(ax, network, values, label="conc")
    assert isinstance(lines, LineCollection)
    assert len(lines.get_paths()) == len(pipes)
    array = lines.get_array()
    assert array.mask[0] and not array.mask[1:].any()  # NaN drawn as missing
    assert lines.colorbar is not None
    assert lines.get_cmap().get_bad() == pytest.approx(matplotlib.colors.to_rgba("0.85"))


def test_plot_pipes_missing_pipes_are_masked(ax, sub):
    first = sub.to_geodataframe("pipes").index[0]
    lines = maps.plot_pipes(ax, sub, {first: 3.0}, colorbar=False)
    assert lines.get_array().count() == 1
    assert lines.colorbar is None


def test_plot_pipes_log_norm_masks_nonpositive(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    values = pd.Series(np.linspace(0, 10, len(ids)), index=ids)
    lines = maps.plot_pipes(ax, sub, values, norm="log")
    ax.figure.canvas.draw()
    assert isinstance(lines.norm, matplotlib.colors.LogNorm)
    assert lines.norm.vmin == pytest.approx(values[values > 0].min())


def test_plot_pipes_widths(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    flow = pd.Series(np.linspace(1, 3, len(ids)), index=ids)
    lines = maps.plot_pipes(ax, sub, width=flow, width_range=(1, 5))
    widths = np.asarray(lines.get_linewidths())
    assert widths.min() == pytest.approx(1) and widths.max() == pytest.approx(5)
    assert maps.plot_pipes(ax, sub, width=2).get_linewidths() == pytest.approx([2])


def test_plot_pipes_in_metric_crs(ax, network):
    maps.plot_pipes(ax, network)
    xmin, xmax = ax.get_xlim()
    assert 1e5 < xmin < xmax < 1e6  # UTM eastings, not longitudes


# ---- Ribbons --------------------------------------------------------------------------------


def test_ribbon_strand_widths_and_no_overlap():
    line = LineString([(0, 0), (100, 0), (150, 50)])
    shares = (0.5, 0.3, 0.2)
    strands = maps.ribbon_strands(line, shares, width=10)
    assert all(s.is_valid and not s.is_empty for s in strands)
    widths = [s.area / line.length for s in strands]
    assert np.array(widths) / sum(widths) == pytest.approx(shares, abs=0.01)
    assert sum(widths) == pytest.approx(10, rel=0.02)
    for i in range(3):
        for j in range(i + 1, 3):
            assert strands[i].intersection(strands[j]).area == pytest.approx(0, abs=1e-6)
    # Right to left looking along the line: on its first (eastward) leg, increasing y.
    first_leg = LineString([(0, -20), (100, -20), (100, 20), (0, 20)]).envelope
    ys = [s.intersection(first_leg).centroid.y for s in strands]
    assert ys == sorted(ys)
    assert ys[0] < 0 < ys[-1]


def test_ribbon_is_centred():
    line = LineString([(0, 0), (100, 0)])
    union = maps.ribbon_strands(line, (1, 1), width=8)
    ymin = min(s.bounds[1] for s in union)
    ymax = max(s.bounds[3] for s in union)
    assert (ymin, ymax) == pytest.approx((-4, 4))


def test_shares_table_top_k_lumps_into_other():
    shares = pd.DataFrame({"a": [0.5, 0.1], "b": [0.3, 0.6], "c": [0.2, 0.3]}, index=[1000, 1001])
    table = maps._shares_table(shares, top_k=2, min_share=0.0)
    assert list(table.columns) == ["a", "b", "c", "other"]
    assert table.loc[1000].to_dict() == pytest.approx({"a": 0.5, "b": 0.3, "c": 0, "other": 0.2})
    assert table.loc[1001].to_dict() == pytest.approx({"a": 0, "b": 0.6, "c": 0.3, "other": 0.1})


def test_shares_table_normalises_and_applies_min_share():
    shares = pd.DataFrame({"a": [6.0], "b": [3.0], "c": [1.0]}, index=[1000])
    table = maps._shares_table(shares, top_k=None, min_share=0.15)
    assert table.loc[1000].to_dict() == pytest.approx({"a": 0.6, "b": 0.3, "other": 0.1})


def test_shares_must_not_be_negative():
    with pytest.raises(ValueError, match="negative"):
        maps._shares_table(pd.DataFrame({"a": [-1.0]}), None, 0.0)


def test_plot_ribbons(ax, sub):
    ids = list(sub.to_geodataframe("pipes").index)
    rng = np.random.default_rng(0)
    shares = pd.DataFrame(rng.random((len(ids), 4)), index=ids, columns=list("abcd"))
    ribbons = maps.plot_ribbons(ax, sub, shares, width=5, top_k=2)
    assert isinstance(ribbons, PatchCollection)
    assert len(ribbons.get_paths()) == 3 * len(ids)  # top 2 + other on every pipe
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels[-1] == "other"


def test_plot_ribbons_colours(ax, sub):
    pipe = sub.to_geodataframe("pipes").index[0]
    shares = pd.DataFrame({"x": [0.25], "y": [0.75]}, index=[pipe])
    ribbons = maps.plot_ribbons(ax, sub, shares, {"x": "red", "y": "blue"}, order="share")
    faces = [matplotlib.colors.to_hex(c) for c in ribbons.get_facecolors()]
    assert faces == ["#0000ff", "#ff0000"]  # largest share first
    with pytest.raises(ValueError, match="colour"):
        maps.plot_ribbons(ax, sub, shares, {"x": "red"})


# ---- Group polygons -------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["voronoi", "concave_hull"])
def test_node_group_polygons(sub, method):
    polygons = maps.node_group_polygons(sub, GROUPS, method)
    assert list(polygons.index) == list(GROUPS)
    assert polygons.crs == maps.plot_crs(sub)
    nodes = maps._geodata(sub, "nodes").geometry
    for group, geometry in polygons.geometry.items():
        assert geometry.is_valid and geometry.area > 0, group
        for node in GROUPS[group]:
            assert geometry.covers(nodes[node]), (group, node)


def test_voronoi_tiles_partition_the_hull(sub):
    polygons = maps.node_group_polygons(sub, GROUPS)
    nodes = maps._geodata(sub, "nodes").geometry
    hull = nodes.union_all().convex_hull.buffer(10.0)
    union = polygons.geometry.union_all()
    assert union.area == pytest.approx(hull.area, rel=1e-6)
    assert polygons.geometry.area.sum() == pytest.approx(hull.area, rel=1e-6)  # no overlap


def test_voronoi_clip_to(sub):
    nodes = maps._geodata(sub, "nodes").geometry
    clip = nodes.union_all().convex_hull.buffer(-1.0)
    polygons = maps.node_group_polygons(sub, GROUPS, clip_to=clip)
    assert polygons.geometry.union_all().difference(clip).area == pytest.approx(0, abs=1e-6)
    clip_gdf = gpd.GeoDataFrame(geometry=[clip], crs=polygons.crs).to_crs(4326)
    again = maps.node_group_polygons(sub, GROUPS, clip_to=clip_gdf)
    assert again.geometry.area.sum() == pytest.approx(polygons.geometry.area.sum(), rel=1e-4)


def test_node_group_polygons_errors(sub):
    with pytest.raises(ValueError, match="not in the network"):
        maps.node_group_polygons(sub, {"x": ["NOPE"]})
    with pytest.raises(ValueError, match="in groups"):
        maps.node_group_polygons(sub, {"x": ["A1"], "y": ["A1"]})
    with pytest.raises(ValueError, match="method"):
        maps.node_group_polygons(sub, GROUPS, "kriging")


def test_plot_groups(ax, sub):
    polygons = maps.node_group_polygons(sub, GROUPS)
    shaded = maps.plot_groups(
        ax, polygons, {"A": 1.0, "B": 2.0, "trunk": 3.0}, hatch_groups=[["A", "B"], ["C"]]
    )
    assert isinstance(shaded, PatchCollection)
    assert len(shaded.get_paths()) == len(GROUPS)
    assert shaded.get_array().mask.tolist() == [False, False, True, False]  # C has no value
    hatched = [c for c in ax.collections if c.get_hatch()]
    assert [c.get_label() for c in hatched] == ["set 1", "set 2"]
    assert maps.plot_groups(ax, polygons).get_array() is None


# ---- Markers and basemaps -------------------------------------------------------------------


def test_plot_markers(ax, sub):
    markers = maps.plot_markers(ax, sub, ["A1", "OUT"], sizes={"A1": 1.0, "OUT": 3.0})
    nodes = maps._geodata(sub, "nodes").geometry
    expected = [[nodes[n].x, nodes[n].y] for n in ["A1", "OUT"]]
    np.testing.assert_allclose(markers.get_offsets(), expected)
    assert markers.get_sizes().tolist() == pytest.approx([20, 200])
    assert sorted(t.get_text() for t in ax.texts) == ["A1", "OUT"]


def test_plot_markers_custom_labels_and_errors(ax, sub):
    maps.plot_markers(ax, sub, ["A1", "B1"], labels={"A1": "Sampler 1"})
    assert [t.get_text() for t in ax.texts] == ["Sampler 1"]
    with pytest.raises(ValueError, match="not in the network"):
        maps.plot_markers(ax, sub, ["D1"])


def test_add_basemap_vector_layer(ax, sub, tmp_path):
    maps.plot_pipes(ax, sub)
    road = gpd.GeoDataFrame(geometry=[LineString([(-81.31, 42.98), (-81.28, 42.99)])], crs=4326)
    road.to_file(tmp_path / "roads.geojson")
    before = len(ax.collections)
    maps.add_basemap(ax, sub, layer=tmp_path / "roads.geojson")
    assert len(ax.collections) == before + 1
    assert ax.collections[-1].get_zorder() == 0
    with pytest.raises(ValueError, match="CRS"):
        maps.add_basemap(ax, sub, layer=road.set_crs(None, allow_override=True))


def test_add_basemap_tiles(ax, sub, monkeypatch):
    calls = {}

    def add_basemap(ax, **kwargs):
        calls.update(kwargs)

    fake = types.SimpleNamespace(
        add_basemap=add_basemap,
        providers=types.SimpleNamespace(CartoDB=types.SimpleNamespace(Positron="positron")),
    )
    monkeypatch.setitem(sys.modules, "contextily", fake)
    maps.add_basemap(ax, sub)
    assert calls["source"] == "positron"
    assert CRS.from_user_input(calls["crs"]) == maps.plot_crs(sub)


def test_add_basemap_without_contextily(ax, sub, monkeypatch):
    monkeypatch.setitem(sys.modules, "contextily", None)
    with pytest.raises(ImportError, match=r"\[maps\]"):
        maps.add_basemap(ax, sub)


# ---- Interactive export ---------------------------------------------------------------------


def test_to_leaflet(sub, tmp_path):
    pytest.importorskip("folium")
    ids = sub.to_geodataframe("pipes").index
    path = tmp_path / "map.html"
    maps.to_leaflet(
        sub,
        path,
        [
            maps.PipeLayer(pd.Series(np.arange(len(ids)), index=ids), name="Concentration"),
            maps.GroupLayer(maps.node_group_polygons(sub, GROUPS), {"A": 1.0}, name="Zones"),
            maps.MarkerLayer(["A1", "OUT"], name="Samplers", labels={"OUT": "Plant"}),
        ],
        tiles=None,
    )
    html = path.read_text()
    for name in ["Concentration", "Zones", "Samplers", "Plant", str(ids[0])]:
        assert name in html, name


def test_to_leaflet_default_layer(network, tmp_path):
    pytest.importorskip("folium")
    maps.to_leaflet(network, tmp_path / "map.html", tiles=None)
    assert "Pipes" in (tmp_path / "map.html").read_text()


def test_to_geojson(sub, tmp_path):
    values = pd.Series({"A1": 7.0}, name="samples")
    maps.to_geojson(sub, tmp_path / "nodes.geojson", kind="nodes", values=values)
    features = json.loads((tmp_path / "nodes.geojson").read_text())["features"]
    by_id = {f["properties"]["node_id"]: f for f in features}
    assert set(by_id) == nodes_of(DENDRITIC)
    assert by_id["A1"]["properties"]["samples"] == 7.0
    lon, lat = by_id["A1"]["geometry"]["coordinates"]
    assert Point(lon, lat).distance(Point(sub.node_positions()["A1"])) < 1e-9


def test_to_geojson_reprojects_to_wgs84(tmp_path):
    projected = Network(gpd.read_file(FIXTURE).to_crs("EPSG:2952"))
    maps.to_geojson(projected, tmp_path / "pipes.geojson")
    written = gpd.read_file(tmp_path / "pipes.geojson")
    assert written.crs == CRS.from_epsg(4326)
    assert written.total_bounds[0] == pytest.approx(-81.3, abs=0.01)
