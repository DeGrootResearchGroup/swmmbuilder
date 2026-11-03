"""Map plotting (#26): to_geodataframe, map styling and keys, static plots, group polygons,
basemaps and export."""

import dataclasses
import json
import sys
import types

import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from conftest import DENDRITIC, FIXTURE, nodes_of, subnetwork, with_pipes
from matplotlib.collections import LineCollection, PatchCollection
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from pyproj import CRS
from shapely import LineString, Point, Polygon

from swmmbuilder import LONDON_ON, Network, maps

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


LONLAT_SCHEMA = dataclasses.replace(LONDON_ON, crs=None)  # keeps the fixture in lon/lat


def test_plot_crs_is_utm_for_lonlat():
    network = Network(FIXTURE, schema=LONLAT_SCHEMA)
    sub = subnetwork(network, DENDRITIC)
    assert maps.plot_crs(network) == CRS.from_epsg(32617)  # London, ON: UTM 17N
    assert maps.plot_crs(sub) == maps.plot_crs(network)
    assert maps.plot_crs(network, "EPSG:3857") == CRS.from_epsg(3857)


def test_plot_crs_keeps_projected_crs(network):
    assert maps.plot_crs(network) == CRS.from_epsg(2958)  # the London schema's CRS
    data = gpd.read_file(FIXTURE).to_crs("EPSG:2952")
    assert maps.plot_crs(Network(data, schema=LONLAT_SCHEMA)) == CRS.from_epsg(2952)


def test_plot_crs_without_crs_is_utm():
    data = gpd.read_file(FIXTURE).set_crs(None, allow_override=True)
    assert maps.plot_crs(Network(data, schema=LONLAT_SCHEMA)) == CRS.from_epsg(32617)


FEET = "EPSG:2263"  # NAD83 / New York Long Island, US survey feet


def _scale_label(ax):
    (bar,) = [a for a in ax.artists if isinstance(a, AnchoredSizeBar)]
    return bar.txt_label.get_text()


def test_plot_crs_in_feet_is_reprojected_to_metres(network):
    feet = Network(gpd.read_file(FIXTURE).to_crs(FEET), schema=LONLAT_SCHEMA)
    assert feet.full_data.crs == CRS.from_user_input(FEET)
    assert maps.plot_crs(feet) == CRS.from_epsg(32617)
    metric, imperial = subnetwork(network, DENDRITIC), subnetwork(feet, DENDRITIC)
    assert maps.network_extent(imperial) == pytest.approx(maps.network_extent(metric), rel=0.01)
    labels = []
    for net in (metric, imperial):
        fig, ax = plt.subplots()
        maps.plot_pipes(ax, net)
        labels.append(_scale_label(ax))
        plt.close(fig)
    assert labels[0] == labels[1] and labels[0].endswith(" m")


@pytest.mark.parametrize("crs", ["EPSG:4326", FEET])
def test_plot_crs_must_be_metres(network, crs):
    with pytest.raises(ValueError, match="projected in metres"):
        maps.plot_crs(network, crs)


def test_plot_crs_without_crs_needs_lonlat():
    data = gpd.read_file(FIXTURE).to_crs("EPSG:2958").set_crs(None, allow_override=True)
    with pytest.raises(ValueError, match="no CRS, and its coordinates are not longitude"):
        maps.plot_crs(Network(data, schema=LONLAT_SCHEMA))


# ---- Scaling and colours --------------------------------------------------------------------


def test_scaled_sqrt_is_proportional_from_zero():
    out = maps._scaled(np.array([0.0, 1.0, 4.0, np.nan, -1.0]), (1, 5), "sqrt")
    assert out == pytest.approx([1, 3, 5, 1, 1])


def test_scaled_linear_and_log():
    assert maps._scaled(np.array([1.0, 2.0, 4.0]), (0, 4), "linear") == pytest.approx([1, 2, 4])
    assert maps._scaled(np.array([1.0, 10.0, 100.0, 0.0]), (1, 3), "log") == pytest.approx(
        [1, 2, 3, 1]
    )
    with pytest.raises(ValueError, match="scale"):
        maps._scaled(np.array([1.0]), (0, 1), "cube")


def test_category_colors():
    colours = maps.category_colors(["Z1", "Z2", "rest", "other"])
    assert colours == {
        "Z1": maps.TABLEAU10[0],
        "Z2": maps.TABLEAU10[1],
        "rest": maps.NEUTRAL_COLOR,
        "other": maps.NEUTRAL_COLOR,
    }
    assert maps.category_colors(["a", "b"], ["red", "blue"]) == {"a": "red", "b": "blue"}
    assert maps.category_colors(["a", "rest"], {"a": "red", "rest": "k"})["rest"] == "k"
    with pytest.raises(ValueError, match="colour"):
        maps.category_colors(["a", "b"], {"a": "red"})


# ---- Map axes and keys ----------------------------------------------------------------------


def test_map_axes_style(ax, sub):
    maps.plot_pipes(ax, sub)
    maps.plot_markers(ax, sub, ["OUT"])  # restyling replaces the scale bar
    assert not ax.axison
    bars = [a for a in ax.artists if isinstance(a, AnchoredSizeBar)]
    assert len(bars) == 1
    length = bars[0].size_bar.get_children()[0].get_width()
    x0, x1 = ax.get_xlim()
    assert (x1 - x0) / 10 < length <= (x1 - x0) / 5
    assert length / 10 ** np.floor(np.log10(length)) in (1, 2, 5)


def test_map_axes_north_arrow_and_no_style(ax, sub):
    maps.map_axes(ax, north_arrow=True)
    assert "N" in [t.get_text() for t in ax.texts]
    fig, plain = plt.subplots()
    maps.plot_pipes(plain, sub, style=False)
    assert plain.axison and not plain.artists
    plt.close(fig)


def test_nice_length():
    assert [maps._nice_length(x) for x in (0.9, 1.7, 3.0, 740.0, 1999.0)] == [0.5, 1, 2, 500, 1000]


def test_add_legends_stacks_keys_outside(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    values = pd.Series(np.arange(1, len(ids) + 1, dtype=float), index=ids)
    first = maps.plot_pipes(ax, sub, values, label="conc", legend=False)
    second = maps.plot_pipes(ax, sub, values * 2, label="flow", legend=False)
    groups = maps.plot_groups(ax, maps.node_group_polygons(sub, GROUPS), legend=False)
    assert first.colorbar is None and groups.legend is None
    maps.add_legends(ax, first, second, groups)
    tops = [layer.colorbar.ax.get_position().y1 for layer in (first, second)]
    assert tops[0] > tops[1]
    for layer in (first, second):
        assert layer.colorbar.ax.get_position().x0 > ax.get_position().x1  # outside the map
    assert [t.get_text() for t in groups.legend.get_texts()] == list(GROUPS)


def test_log_colorbar_has_no_minor_labels(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    values = pd.Series(np.logspace(0, 6, len(ids)), index=ids)
    layer = maps.plot_pipes(ax, sub, values, norm="log", label="conc")
    assert isinstance(
        layer.colorbar.ax.yaxis.get_minor_formatter(), matplotlib.ticker.NullFormatter
    )


# ---- plot_pipes -----------------------------------------------------------------------------


def test_plot_pipes_colours_by_value(ax, network):
    pipes = network.to_geodataframe("pipes")
    values = pd.Series(np.arange(len(pipes), dtype=float), index=pipes.index)
    values.iloc[0] = np.nan
    layer = maps.plot_pipes(ax, network, values, label="conc")
    lines = layer.artist
    assert isinstance(lines, LineCollection)
    assert len(lines.get_paths()) == len(pipes) == len(layer.ids)
    array = lines.get_array()
    assert array.mask[0] and not array.mask[1:].any()  # NaN drawn as missing
    assert layer.mappable is lines and layer.colorbar is not None
    assert layer.colorbar.ax.yaxis.label.get_text() == "conc"
    assert lines.get_cmap().get_bad() == pytest.approx(matplotlib.colors.to_rgba("0.75"))


def test_plot_pipes_missing_pipes_are_masked(ax, sub):
    first = sub.to_geodataframe("pipes").index[0]
    layer = maps.plot_pipes(ax, sub, {first: 3.0}, legend=False)
    assert layer.artist.get_array().count() == 1
    assert layer.colorbar is None


def test_plot_pipes_log_norm_masks_nonpositive(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    values = pd.Series(np.linspace(0, 10, len(ids)), index=ids)
    lines = maps.plot_pipes(ax, sub, values, norm="log").artist
    ax.figure.canvas.draw()
    assert isinstance(lines.norm, matplotlib.colors.LogNorm)
    assert lines.norm.vmin == pytest.approx(values[values > 0].min())


def test_plot_pipes_sqrt_widths_thin_first(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    flow = pd.Series(np.linspace(1, 1000, len(ids))[::-1], index=ids)
    layer = maps.plot_pipes(ax, sub, width=flow, width_range=(0, 4))
    widths = np.asarray(layer.artist.get_linewidths())
    assert list(widths) == sorted(widths)  # thick pipes drawn last
    expected = 4 * np.sqrt(flow[layer.ids] / flow.max())
    assert widths == pytest.approx(expected.to_numpy())
    lines = maps.plot_pipes(ax, sub, width=flow, width_scale="linear", width_range=(0, 4)).artist
    assert max(lines.get_linewidths()) == pytest.approx(4)
    assert maps.plot_pipes(ax, sub, width=2).artist.get_linewidths() == pytest.approx(
        [2] * len(ids)
    )


def test_plot_pipes_as_context(ax, sub):
    ids = sub.to_geodataframe("pipes").index
    flow = pd.Series(np.linspace(1, 100, len(ids)), index=ids)
    full = maps.plot_pipes(ax, sub, width=flow, legend=False).artist
    layer = maps.plot_pipes(ax, sub, width=flow, context=True)
    assert np.asarray(layer.artist.get_linewidths()) == pytest.approx(
        0.6 * np.asarray(full.get_linewidths())
    )
    assert matplotlib.colors.to_hex(layer.artist.get_color()[0]) == matplotlib.colors.to_hex("0.6")
    assert layer.artist.get_zorder() == 1.5  # above group fills, below ribbons
    assert layer.mappable is None
    with pytest.raises(ValueError, match="Context"):
        maps.plot_pipes(ax, sub, flow, context=True)


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


def test_stems_follow_the_widest_pipe(network):
    sub = network.subnetwork(nodes_of(DENDRITIC))
    pipe = {(u, v): k for u, v, k in sub.graph.edges(keys=True)}
    ids = sorted(pipe.values())
    widths = pd.Series(1.0, index=ids)
    widths[pipe["B2", "J1"]] = 2.0  # B is the main stem at J1
    chains = maps._stems(sub.graph, ids, widths)
    assert sorted(p for chain in chains for p in chain) == ids
    trunk = [pipe["J1", "J2"], pipe["J2", "OUT"]]
    assert [pipe["B1", "B2"], pipe["B2", "J1"], *trunk] in chains
    assert [pipe["A1", "A2"], pipe["A2", "J1"]] in chains  # ends at J1
    assert [pipe["C1", "J1"]] in chains


def test_stems_survive_a_cycle():
    net = with_pipes(("A2", "A1"))  # A1 <-> A2
    pipe = {(u, v): k for u, v, k in net.graph.edges(keys=True)}
    ids = [pipe["A1", "A2"], 2000]
    chains = maps._stems(net.graph, ids)
    assert sorted(p for chain in chains for p in chain) == sorted(ids)


def test_merged_splits_only_visible_changes():
    chain = [1, 2, 3, 4]
    table = pd.DataFrame({"x": [0.50, 0.59, 0.60, 0.61]}, chain).assign(y=lambda t: 1 - t.x)
    widths = pd.Series([10.0, 10.2, 10.4, 15.0], index=chain)
    # 4 is much wider than 3; 2 is within 0.5 m and 0.02 of 3; 1 differs from 3 by 0.1.
    assert maps._merged(chain, table, widths, 0.5, 0.02) == [[1], [2, 3], [4]]
    assert maps._merged(chain, table, widths, 0.0, 0.0) == [[1], [2], [3], [4]]
    assert maps._merged(chain, table, widths, 10.0, 1.0) == [chain]


@pytest.mark.parametrize("turn", [90, 30, -60])
def test_mitred_bend_has_no_gap(turn):
    angle = np.radians(turn)
    joint = np.array([100.0, 0.0])
    end = joint + 80 * np.array([np.cos(angle), np.sin(angle)])
    upstream, downstream = LineString([(0, 0), joint]), LineString([joint, end])
    before, after = np.array([1.0, 0.0]), (end - joint) / 80
    shares = (0.6, 0.4)
    first = maps.mitred_strands(upstream, shares, 12, after=after)
    second = maps.mitred_strands(downstream, shares, 8, before=before)
    # Strand offsets: (-6, 1.2) and (1.2, 6) upstream, (-4, 0.8) and (0.8, 4) downstream.
    for a, b, overlap in zip(first, second, (4.8, 2.8), strict=True):
        assert_meet_along_cut(a, b, overlap)
    # Square-cut ribbons overlap on the inside of the bend and leave a wedge on the outside.
    square = maps.ribbon_strands(upstream, shares, 12), maps.ribbon_strands(downstream, shares, 8)
    inner = 1 if turn > 0 else 0  # the left of the line is the inside of a left turn
    assert square[0][inner].intersection(square[1][inner]).area > 1


def assert_meet_along_cut(a, b, width):
    """Strands ``a`` and ``b`` don't overlap and share an edge at least ``width`` long."""
    assert a.intersection(b).area == pytest.approx(0, abs=1e-6)
    eps = 1e-6
    shared = a.buffer(eps).intersection(b.buffer(eps)).area / (2 * eps)
    assert shared >= 0.99 * width


def test_plot_ribbons_trunk_is_one_band(ax, sub):
    pipe = {(u, v): k for u, v, k in sub.graph.edges(keys=True)}
    stem = [pipe["A1", "A2"], pipe["A2", "J1"], pipe["J1", "J2"], pipe["J2", "OUT"]]
    shares = pd.DataFrame({"x": [1.0, 1.01, 2.0, 2.02], "y": [1.0, 1.0, 2.0, 2.0]}, index=stem)
    layer = maps.plot_ribbons(ax, sub, shares, context=False, max_width=10)
    assert len(layer.artist.get_paths()) == 4  # A1 -> J1 and J1 -> OUT, two strands each
    pieces = [Polygon(p.to_polygons()[0]) for p in layer.artist.get_paths()]
    assert_meet_along_cut(pieces[0], pieces[2], 0.99 * 5 / 2)  # strand x: offsets (-2.5, 0)
    layer = maps.plot_ribbons(ax, sub, shares, context=False, max_width=10, merge_tol=0)
    assert len(layer.artist.get_paths()) == 8


def test_plot_ribbons(ax, sub):
    ids = list(sub.to_geodataframe("pipes").index)
    rng = np.random.default_rng(0)
    shares = pd.DataFrame(rng.random((len(ids), 4)), index=ids, columns=list("abcd"))
    layer = maps.plot_ribbons(ax, sub, shares, width=5, top_k=2)
    assert isinstance(layer.artist, PatchCollection)
    assert len(layer.artist.get_paths()) == 3 * len(ids)  # top 2 + other; every pipe differs
    assert [t.get_text() for t in layer.legend.get_texts()][-1] == "other"
    assert len(ax.collections) == 2  # grey context beneath the ribbons


def test_plot_ribbons_merges_reaches(ax, sub):
    pipe = {(u, v): k for u, v, k in sub.graph.edges(keys=True)}
    chain = [pipe["A1", "A2"], pipe["A2", "J1"]]
    shares = pd.DataFrame({"x": [1.0, 1.0], "y": [3.0, 3.0]}, index=chain)
    layer = maps.plot_ribbons(ax, sub, shares, width=6, context=False)
    assert len(layer.artist.get_paths()) == 2  # one ribbon of two strands for the reach
    assert [t.get_text() for t in layer.legend.get_texts()] == ["x", "y"]  # no 'other'


def test_plot_ribbons_width_proportional_to_load(ax, sub):
    pipe = {(u, v): k for u, v, k in sub.graph.edges(keys=True)}
    big, small, tiny = pipe["J2", "OUT"], pipe["B1", "B2"], pipe["C1", "J1"]
    shares = pd.DataFrame({"x": [6.0, 3.0, 0.3], "y": [4.0, 0.0, 0.0]}, index=[big, small, tiny])
    layer = maps.plot_ribbons(ax, sub, shares, max_width=20, min_width=3, context=False)
    assert layer.ids == [big, small]  # tiny is 0.6 m wide, below min_width
    lengths = maps._geodata(sub, "pipes").geometry.length
    areas = [Polygon(p.to_polygons()[0]).area for p in layer.artist.get_paths()]
    assert len(areas) == 3  # x and y on the big pipe, x on the small one
    assert sum(areas[:2]) / lengths[big] == pytest.approx(20, rel=0.05)
    assert areas[2] / lengths[small] == pytest.approx(6, rel=0.05)  # 3 / 10 of the largest load


def test_ribbon_widths_default_to_the_extent(ax, sub):
    pipe = {(u, v): k for u, v, k in sub.graph.edges(keys=True)}
    big, small, tiny = pipe["J2", "OUT"], pipe["B1", "B2"], pipe["C1", "J1"]
    shares = pd.DataFrame({"x": [100.0, 10.0, 4.0]}, index=[big, small, tiny])
    layer = maps.plot_ribbons(ax, sub, shares, context=False)
    max_width = 0.025 * maps.network_extent(sub)
    assert layer.ids == [big, small]  # tiny is 4 % of the widest, below the 5 % default
    lengths = maps._geodata(sub, "pipes").geometry.length
    areas = [Polygon(p.to_polygons()[0]).area for p in layer.artist.get_paths()]
    assert areas[0] / lengths[big] == pytest.approx(max_width, rel=0.05)


def test_network_extent(sub):
    xmin, ymin, xmax, ymax = maps._geodata(sub, "pipes").total_bounds
    assert maps.network_extent(sub) == pytest.approx(max(xmax - xmin, ymax - ymin))


def test_plot_ribbons_colours(ax, sub):
    pipe = sub.to_geodataframe("pipes").index[0]
    shares = pd.DataFrame({"x": [0.25], "y": [0.75]}, index=[pipe])
    layer = maps.plot_ribbons(ax, sub, shares, {"x": "red", "y": "blue"}, order="share")
    faces = [matplotlib.colors.to_hex(c) for c in layer.artist.get_facecolors()]
    assert faces == ["#0000ff", "#ff0000"]  # largest share first
    with pytest.raises(ValueError, match="colour"):
        maps.plot_ribbons(ax, sub, shares, {"x": "red"})


def test_ribbons_and_groups_share_colours(ax, sub):
    polygons = maps.node_group_polygons(sub, GROUPS)
    groups = maps.plot_groups(ax, polygons)
    ids = list(sub.to_geodataframe("pipes").index)
    shares = pd.DataFrame(1.0, index=ids, columns=list(GROUPS))
    ribbons = maps.plot_ribbons(ax, sub, shares)
    fills = {h.get_label(): matplotlib.colors.to_hex(h.get_facecolor()) for h in groups.handles}
    strands = {h.get_label(): matplotlib.colors.to_hex(h.get_facecolor()) for h in ribbons.handles}
    assert fills == strands


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


def test_voronoi_tiles_partition_the_footprint(sub):
    polygons = maps.node_group_polygons(sub, GROUPS, footprint=50)
    footprint = maps.network_footprint(sub, 50)
    assert polygons.geometry.union_all().area == pytest.approx(footprint.area, rel=1e-6)
    assert polygons.geometry.area.sum() == pytest.approx(footprint.area, rel=1e-6)  # no overlap
    pipes = maps._geodata(sub, "pipes").geometry
    assert footprint.contains(pipes.union_all())
    assert footprint.area < pipes.union_all().convex_hull.buffer(50).area


def test_footprint_distance_scales_with_the_network(sub):
    lengths = maps._geodata(sub, "pipes").geometry.length
    expected = min(2 * lengths.median(), 0.05 * maps.network_extent(sub))
    assert maps.footprint_distance(sub) == pytest.approx(expected)
    default = maps.network_footprint(sub)
    assert default.area == pytest.approx(maps.network_footprint(sub, expected).area)
    polygons = maps.node_group_polygons(sub, GROUPS)
    assert polygons.geometry.union_all().area == pytest.approx(default.area, rel=1e-6)


def test_voronoi_without_footprint_fills_the_hull(sub):
    polygons = maps.node_group_polygons(sub, GROUPS, footprint=False)
    nodes = maps._geodata(sub, "nodes").geometry
    hull = nodes.union_all().convex_hull.buffer(10.0)
    assert polygons.geometry.area.sum() == pytest.approx(hull.area, rel=1e-6)


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


def test_plot_groups_by_category(ax, sub):
    polygons = maps.node_group_polygons(sub, GROUPS)
    layer = maps.plot_groups(ax, polygons, hatch_groups={"North": ["A", "B"]})
    assert layer.mappable is None and layer.colorbar is None
    faces = [matplotlib.colors.to_hex(c) for c in layer.artist.get_facecolors()]
    assert faces == list(maps.TABLEAU10[:4])
    labels = [t.get_text() for t in layer.legend.get_texts()]
    assert labels == [*GROUPS, "North"]
    hatched = [c for c in ax.collections if c.get_hatch()]
    assert matplotlib.colors.to_hex(hatched[0].get_edgecolor()[0]) == matplotlib.colors.to_hex(
        "0.55"
    )


def test_plot_groups_by_value(ax, sub):
    polygons = maps.node_group_polygons(sub, GROUPS)
    layer = maps.plot_groups(
        ax, polygons, {"A": 1.0, "B": 2.0, "trunk": 3.0}, hatch_groups=[["A", "B"], ["C"]]
    )
    assert isinstance(layer.artist, PatchCollection)
    assert len(layer.artist.get_paths()) == len(GROUPS)
    assert layer.artist.get_array().mask.tolist() == [False, False, True, False]  # C: no value
    assert layer.colorbar is not None
    assert [h.get_label() for h in layer.handles] == ["set 1", "set 2"]
    single = maps.plot_groups(ax, polygons, color="tan")
    assert single.artist.get_array() is None and not single.handles


# ---- Markers and basemaps -------------------------------------------------------------------


def test_plot_markers(ax, sub):
    layer = maps.plot_markers(ax, sub, ["A1", "OUT"], sizes={"A1": 1.0, "OUT": 4.0}, labels=True)
    nodes = maps._geodata(sub, "nodes").geometry
    expected = [[nodes[n].x, nodes[n].y] for n in ["A1", "OUT"]]
    np.testing.assert_allclose(layer.artist.get_offsets(), expected)
    assert layer.artist.get_sizes().tolist() == pytest.approx([20 + 180 / 4, 200])
    assert sorted(t.get_text() for t in ax.texts) == ["A1", "OUT"]
    assert all(t.get_path_effects() for t in ax.texts)  # white halo


def test_plot_markers_labels_and_errors(ax, sub):
    maps.plot_markers(ax, sub, ["A1", "B1"])
    assert not ax.texts  # no labels by default
    layer = maps.plot_markers(ax, sub, ["A1", "B1"], labels={"A1": "Sampler 1"}, name="Samplers")
    assert [t.get_text() for t in ax.texts] == ["Sampler 1"]
    assert [t.get_text() for t in layer.legend.get_texts()] == ["Samplers"]
    with pytest.raises(ValueError, match="not in the network"):
        maps.plot_markers(ax, sub, ["D1"])


def test_add_basemap_vector_layer(ax, sub, tmp_path):
    maps.plot_pipes(ax, sub)
    limits = ax.get_xlim(), ax.get_ylim()
    road = gpd.GeoDataFrame(geometry=[LineString([(-81.4, 42.9), (-81.2, 43.1)])], crs=4326)
    road.to_file(tmp_path / "roads.geojson")
    before = len(ax.collections)
    maps.add_basemap(ax, sub, layer=tmp_path / "roads.geojson")
    assert len(ax.collections) == before + 1
    assert ax.collections[-1].get_zorder() == 0
    assert (ax.get_xlim(), ax.get_ylim()) == limits  # the long road does not zoom out
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
    assert Point(lon, lat).distance(Point(maps.lonlat_positions(sub)["A1"])) < 1e-9


def test_to_geojson_reprojects_to_wgs84(tmp_path):
    projected = Network(gpd.read_file(FIXTURE).to_crs("EPSG:2952"), schema=LONLAT_SCHEMA)
    maps.to_geojson(projected, tmp_path / "pipes.geojson")
    written = gpd.read_file(tmp_path / "pipes.geojson")
    assert written.crs == CRS.from_epsg(4326)
    assert written.total_bounds[0] == pytest.approx(-81.3, abs=0.01)
