"""Node IDs reused for structures far apart, split by ``Network(split_tolerance=...)``."""

import logging

import geopandas as gpd
import networkx as nx
import pytest
from conftest import FIXTURE
from shapely.geometry import LineString, Point

from swmmbuilder import Network

CRS = "EPSG:2958"
X0, Y0 = 480_000.0, 4_760_000.0

# Two subdivisions 900 m apart that reuse the IDs X and Y. In the west one X drains to Y, in the
# east one Y drains to X, so read by ID alone the two make a cycle X -> Y -> X. Pipes are
# (pipe, up, down, (x, y) of each vertex in m from (X0, Y0)).
WEST = [
    (1, "X", "Y", [(0, 0), (100, 0)]),
    (2, "Y", "OUT", [(100, 0), (200, 0)]),
]
EAST = [
    (11, "Y", "X", [(1000, 0), (1100, 0)]),
    (12, "X", "E2", [(1100, 0), (1200, 0)]),
]
# Maintenance-hole points of the west subdivision only.
HOLES = {"X": (0, 0), "Y": (100, 0), "OUT": (200, 0)}


def pipes(spec):
    """The pipes of ``spec`` in the London schema."""
    template = gpd.read_file(FIXTURE).iloc[0].drop("geometry").to_dict()
    rows = []
    for pipe, up, down, xy in spec:
        line = None if xy is None else LineString([(X0 + x, Y0 + y) for x, y in xy])
        rows.append(
            {
                **template,
                "OBJECTID": pipe,
                "UpstreamID": up,
                "DownstreamID": down,
                "Length": 100.0 if line is None else line.length,
                "geometry": line,
            }
        )
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=CRS)


def holes(points):
    """Maintenance-hole points {ID: (x, y)}, with a rim elevation."""
    return gpd.GeoDataFrame(
        {"GIS_FeatureKey": list(points), "LidElevation": 260.0, "HeadwallElevation": None},
        geometry=[Point(X0 + x, Y0 + y) for x, y in points.values()],
        crs=CRS,
    )


def load(spec=WEST + EAST, points=HOLES, **options):
    options = {"split_tolerance": 50.0, **options}
    return Network(pipes(spec), maintenance_holes=holes(points), **options)


def ends(net, pipe):
    """(from_node, to_node) of ``pipe`` in the network."""
    row = net.full_data.set_index("pipe_id").loc[pipe]
    return row["from_node"], row["to_node"]


def test_reused_ids_are_split_into_one_node_per_structure():
    net = load()
    assert net.splits == {"X": ["Xa", "Xb"], "Y": ["Ya", "Yb"]}
    assert ends(net, 1) == ("Xa", "Ya")
    assert ends(net, 2) == ("Ya", "OUT")
    assert ends(net, 11) == ("Yb", "Xb")
    assert ends(net, 12) == ("Xb", "E2")
    assert "X" not in net.graph and "Y" not in net.graph


def test_splitting_removes_the_false_cycle():
    assert not nx.is_directed_acyclic_graph(load(split_tolerance=None).graph)
    assert nx.is_directed_acyclic_graph(load().graph)


def test_off_by_default():
    net = load(split_tolerance=None)
    assert net.splits == {} and net.split_report is None
    assert ends(net, 11) == ("Y", "X")


def test_report_flags_every_split():
    report = load().split_report.set_index("new_id")
    assert list(report.index) == ["Xa", "Xb", "Ya", "Yb"]
    assert report.loc["Xb", "node_id"] == "X"
    assert report.loc["Xa", "pipes"] == (1,)
    assert report.loc["Xb", "pipes"] == (11, 12)
    assert report["separation"].to_dict() == pytest.approx(
        {"Xa": 1100.0, "Xb": 1100.0, "Ya": 900.0, "Yb": 900.0}
    )


def test_split_is_logged_as_a_warning(caplog):
    with caplog.at_level(logging.WARNING, logger="swmmbuilder.network"):
        load()
    text = caplog.text
    assert "Node ID X is used for 2 structures, at least 1100 m apart" in text
    assert "Xb (pipes 11, 12; no node-layer point: attributes from its pipes)" in text


def test_the_point_goes_to_the_nearest_structure():
    net = load()
    report = net.split_report.set_index("new_id")
    # The west structures have the maintenance-hole points; the east ones take their pipes'.
    assert report.loc["Xa", "attributes"] == "node layer"
    assert report.loc["Xa", "layer"] == "maintenance_holes"
    assert report.loc["Xa", "point_distance"] == pytest.approx(0.0)
    assert report.loc["Xb", "attributes"] == "pipes"
    assert net.nodes.loc["Xa", "rim_elevation"] == 260.0
    assert net.nodes.loc["Xa", "layer"] == "maintenance_holes"
    assert net.nodes["layer"].isna()["Xb"]
    # Attributes from the pipes: the lowest invert of the pipes connected to it.
    assert net.nodes.loc["Xb", "invert"] == pytest.approx(104.5)
    # Each node is at its own structure.
    assert net.nodes.geometry["Xa"].distance(Point(X0, Y0)) == pytest.approx(0.0)
    assert net.nodes.geometry["Xb"].distance(Point(X0 + 1100, Y0)) == pytest.approx(0.0)


def test_a_point_at_the_other_structure_goes_there():
    net = load(points={"X": (1100, 0)})
    report = net.split_report.set_index("new_id")
    assert report.loc["Xa", "attributes"] == "pipes"
    assert report.loc["Xb", "attributes"] == "node layer"
    assert net.nodes.loc["Xb", "rim_elevation"] == 260.0


def test_letters_follow_the_smallest_pipe_id_whatever_the_row_order():
    # Three structures of Z; the one at x = 2000 has the smallest pipe ID.
    spec = [
        (32, "Z", "A", [(0, 300), (100, 300)]),
        (33, "Z", "B", [(1000, 300), (1100, 300)]),
        (31, "Z", "C", [(2000, 300), (2100, 300)]),
    ]
    for order in (spec, spec[::-1], [spec[1], spec[2], spec[0]]):
        net = load(order, points={})
        assert net.splits == {"Z": ["Za", "Zb", "Zc"]}
        assert ends(net, 31)[0] == "Za"
        assert ends(net, 32)[0] == "Zb"
        assert ends(net, 33)[0] == "Zc"


def test_ends_within_the_tolerance_are_one_structure():
    # Pipe ends of one maintenance hole drawn up to 30 m apart, and chained: 0, 30 and 60 m.
    spec = [
        (1, "A", "M", [(-100, 0), (0, 0)]),
        (2, "B", "M", [(30, -100), (30, 0)]),
        (3, "M", "C", [(60, 0), (160, 0)]),
    ]
    net = load(spec, points={})
    assert net.splits == {}
    assert ends(net, 3) == ("M", "C")


def test_a_pipe_drawn_backwards_does_not_split_its_nodes(caplog):
    # Pipe 21 runs P -> Q but is drawn from Q (x = 560) to P (x = 500).
    spec = [
        (20, "O", "P", [(400, 0), (500, 0)]),
        (21, "P", "Q", [(560, 0), (500, 0)]),
        (22, "Q", "R", [(560, 0), (660, 0)]),
    ]
    with caplog.at_level(logging.INFO, logger="swmmbuilder.network"):
        net = load(spec, points={}, split_tolerance=30.0)
    assert net.splits == {}
    assert ends(net, 21) == ("P", "Q")
    assert "1 pipes are drawn from their downstream node to their upstream node" in caplog.text
    # Its geometry is reversed, so P sits at its own end, not at Q (#121).
    assert net.reversed_pipes == [21]
    line = net.full_data.set_index("pipe_id").geometry[21]
    assert line.coords[0] == (X0 + 500, Y0) and line.coords[-1] == (X0 + 560, Y0)
    positions = net.subnetwork(net.graph.nodes).node_positions()
    assert positions["P"] == (X0 + 500, Y0)
    assert positions["Q"] == (X0 + 560, Y0)
    assert net.nodes.geometry["P"] == Point(X0 + 500, Y0)


def test_no_pipe_is_reversed_without_split_tolerance():
    spec = [(20, "O", "P", [(400, 0), (500, 0)]), (21, "P", "Q", [(560, 0), (500, 0)])]
    net = Network(pipes(spec))
    assert net.reversed_pipes == []


def test_a_self_loop_across_two_structures_is_split():
    # Pipe 5 has the same ID at both ends, 300 m apart.
    spec = [(5, "S", "S", [(0, 500), (300, 500)]), (6, "S", "T", [(300, 500), (400, 500)])]
    net = load(spec, points={})
    assert ends(net, 5) == ("Sa", "Sb")
    assert ends(net, 6) == ("Sb", "T")


def test_a_mistyped_end_snaps_to_the_structure_it_is_at():
    # Pipe 43's downstream end is at V but carries W, the ID of a structure 500 m away.
    spec = [
        (40, "W1", "W", [(-100, 200), (0, 200)]),
        (41, "W", "W2", [(0, 200), (100, 200)]),
        (42, "V", "V2", [(500, 200), (600, 200)]),
        (43, "V1", "W", [(400, 200), (500, 200)]),
    ]
    net = load(spec, points={"W": (0, 200), "V": (500, 200)}, snap_tolerance=0.5)
    assert net.splits == {"W": ["Wa", "Wb"]}
    assert net.resolve("Wb") == "V"
    assert ends(net, 43) == ("V1", "V")
    assert ends(net, 41) == ("Wa", "W2")


def test_split_tolerance_must_exceed_snap_tolerance():
    with pytest.raises(ValueError, match="must exceed snap_tolerance"):
        load(split_tolerance=0.5, snap_tolerance=0.5)


def test_snapping_never_merges_the_structures_of_one_id():
    net = load(snap_tolerance=0.5)
    # Move Xb onto Xa and snap again: they stay apart.
    net.nodes.loc["Xb", "geometry"] = net.nodes.geometry["Xa"]
    report = net.snap_nodes(0.5)
    row = report[(report["node_1"] == "Xa") & (report["node_2"] == "Xb")].iloc[0]
    assert not row["snapped"] and row["reason"] == "split from one node ID"


def test_the_original_id_names_its_split_ids():
    net = load()
    with pytest.raises(ValueError, match=r"X \(split into Xa, Xb\)"):
        net.catchment("X")
    assert net.catchment("Ya").nodes == ["Xa", "Ya"]


def test_an_id_with_a_pipe_without_geometry_is_not_split(caplog):
    spec = [*WEST, *EAST, (13, "E3", "X", None)]
    with caplog.at_level(logging.WARNING, logger="swmmbuilder.network"):
        net = load(spec)
    assert net.splits == {"Y": ["Ya", "Yb"]}
    assert ends(net, 12) == ("X", "E2")
    assert "a pipe with no geometry, so they are not split: X" in caplog.text


def test_a_new_id_already_in_use_raises():
    with pytest.raises(ValueError, match="already in use: Xa"):
        load([*WEST, *EAST, (99, "Xa", "OUT", [(5000, 0), (5100, 0)])])


def test_geographic_coordinates_are_measured_in_metres():
    data = pipes(WEST + EAST).to_crs("EPSG:4326")
    net = Network(data, maintenance_holes=holes(HOLES).to_crs("EPSG:4326"), split_tolerance=50.0)
    assert net.splits == {"X": ["Xa", "Xb"], "Y": ["Ya", "Yb"]}
    assert net.split_report["separation"].min() == pytest.approx(900.0, rel=1e-3)
