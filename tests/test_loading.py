"""Messy GIS input on loading (#61): pipes with no geometry or no ID, numeric node-layer
columns, multi-part pipe geometry, and zero-length invert interpolation."""

import dataclasses
import logging

import geopandas as gpd
import networkx as nx
import pandas as pd
import pytest
from conftest import FIXTURE, read_sections, write_inp
from shapely.geometry import LineString, MultiLineString, Point

from swmmbuilder import LONDON_ON, Network
from swmmbuilder.gaps import _fill_node_inverts
from swmmbuilder.geometry import is_disjoint, line_coords, reverse_line


@pytest.fixture
def data():
    return gpd.read_file(FIXTURE)


def row(data, pipe_id):
    return data.index[data["OBJECTID"] == pipe_id][0]


def test_pipe_with_no_geometry_is_kept_by_its_length(data, tmp_path, caplog):
    data.loc[row(data, 1004), "geometry"] = None  # C1 -> J1, Length known; C1 has no other pipe
    with caplog.at_level(logging.WARNING, logger="swmmbuilder"):
        net = Network(data)
    assert net.graph.has_edge("C1", "J1", 1004)
    assert "1 pipes have no geometry" in caplog.text
    assert net.nodes.geometry["C1"] is None and net.nodes.geometry["J1"] is not None
    sub = net.subnetwork(["C1", "J1", "J2", "OUT"])
    sections = read_sections(
        write_inp(sub, tmp_path, routing="DYNWAVE", outfalls=("OUT",), dead_ends="outfall")
    )
    assert "1004" in {r[0] for r in sections["CONDUITS"]}
    assert "C1" not in {r[0] for r in sections["COORDINATES"]}


def test_pipe_with_no_geometry_and_no_length_is_removed(data):
    i = row(data, 1000)
    data.loc[i, "geometry"] = None
    data.loc[i, "Length"] = None
    net = Network(data)
    assert not net.graph.has_edge("A1", "A2")
    assert net.removed["invalid"] == {"no length": 1}


def test_pipe_with_no_id_is_removed(data):
    data["OBJECTID"] = data["OBJECTID"].astype(float)
    data.loc[row(data, 1000), "OBJECTID"] = float("nan")
    net = Network(data)
    assert not net.graph.has_edge("A1", "A2")
    assert all(key == key for _, _, key in net.graph.edges(keys=True))  # no NaN key
    assert net.removed["invalid"] == {"no pipe_id": 1}


def test_several_pipes_with_no_id_are_removed_not_duplicates(data):
    # Missing IDs are not duplicates of each other (#120).
    data["OBJECTID"] = data["OBJECTID"].astype(float)
    for pipe in (1000, 1004):
        data.loc[row(data, pipe), "OBJECTID"] = float("nan")
    net = Network(data)
    assert not net.graph.has_edge("A1", "A2")
    assert net.removed["invalid"] == {"no pipe_id": 2}


def test_pipe_ids_that_are_null_values_are_removed_not_duplicates(data):
    schema = dataclasses.replace(LONDON_ON, null_values=frozenset({"None"}))
    data["OBJECTID"] = data["OBJECTID"].astype(object)
    for pipe in (1000, 1004):
        data.loc[row(data, pipe), "OBJECTID"] = "None"
    net = Network(data, schema=schema)
    assert net.removed["invalid"] == {"no pipe_id": 2}


def test_real_duplicate_ids_are_still_refused_beside_missing_ones(data):
    data["OBJECTID"] = data["OBJECTID"].astype(float)
    data.loc[row(data, 1000), "OBJECTID"] = float("nan")
    data.loc[row(data, 1004), "OBJECTID"] = float("nan")
    data.loc[row(data, 1001), "OBJECTID"] = 1002.0
    with pytest.raises(ValueError, match=r"duplicate OBJECTID values \[1002\.0\]"):
        Network(data)


def test_numeric_node_layer_column_covering_some_nodes(data):
    holes = gpd.GeoDataFrame(
        {"GIS_FeatureKey": ["A1", "B1"], "Fullname": [5, 6], "PipesIn": [0, 0]},
        geometry=[Point(0, 0), Point(1, 1)],
        crs=data.crs,
    )
    net = Network(FIXTURE, maintenance_holes=holes)
    assert net.nodes.loc["A1", "name"] == 5
    assert net.nodes.loc["J1", "name"] is pd.NA
    assert str(net.nodes["name"].dtype) == "Int64"


def test_connected_multilinestring_is_one_line():
    line = MultiLineString([[(0, 0), (1, 0)], [(1, 0), (2, 1)]])
    assert line_coords(line) == [(0, 0), (1, 0), (2, 1)]  # no repeated joint vertex
    assert not is_disjoint(line)


def test_disjoint_multilinestring_is_warned(data, caplog):
    line = MultiLineString([[(0, 0), (1, 0)], [(5, 5), (6, 5)]])
    assert is_disjoint(line)
    assert line_coords(line) == [(0, 0), (1, 0), (5, 5), (6, 5)]
    i = row(data, 1000)
    start, end = data.geometry[i].coords[0], data.geometry[i].coords[-1]
    data.loc[i, "geometry"] = MultiLineString([[start, start], [end, end]])
    with caplog.at_level(logging.WARNING, logger="swmmbuilder"):
        Network(data)
    assert "parts do not join" in caplog.text and "1000" in caplog.text


def test_missing_geometry_gives_no_coords():
    assert line_coords(None) == []
    assert line_coords(LineString()) == []


def test_interpolation_between_coincident_known_inverts():
    graph = nx.MultiDiGraph()
    graph.add_edge("U", "X", key=1, length=0.0)
    graph.add_edge("X", "D", key=2, length=0.0)
    inverts, records = {"U": 10.0, "D": 9.0}, []
    _fill_node_inverts(graph, ["X"], inverts, "interpolate", None, lambda *a: records.append(a))
    assert inverts["X"] == 9.0  # no ZeroDivisionError; takes the downstream invert


def test_reverse_line_reverses_parts_and_their_order():
    line = LineString([(0, 0), (1, 0), (2, 0)])
    assert line_coords(reverse_line(line)) == [(2, 0), (1, 0), (0, 0)]
    multi = MultiLineString([[(0, 0), (1, 0)], [(5, 0), (6, 0)]])  # disjoint parts
    assert line_coords(reverse_line(multi)) == [(6, 0), (5, 0), (1, 0), (0, 0)]
    assert reverse_line(None) is None
