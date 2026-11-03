"""Node layers, pipe status and structure type, coincident-node snapping, ID normalisation and
maintenance-hole attributes (#19), on the node-layer fixture (see tests/data/make_fixture.py)."""

import math

import geopandas as gpd
import pandas as pd
import pytest
from conftest import TESTS, read_sections, run_swmm, write_inp
from pyproj import Transformer

from swmmbuilder import Network, Schema, Subnetwork, strip_leading_zeros
from swmmbuilder.swmm import build_model

DATA = TESTS / "data"
PIPES = DATA / "fixture_nodes_pipes.geojson"
LAYERS = {
    name: DATA / f"fixture_nodes_{name}.geojson"
    for name in ("maintenance_holes", "pump_stations", "outfalls", "other_nodes")
}
FORCE_MAIN = 5002  # PS1 -> M3
OVERFLOW = 5005  # M3 -> OF2


def load(**kwargs):
    return Network(PIPES, **{**LAYERS, **kwargs})


@pytest.fixture(scope="module")
def net():
    return load()


@pytest.fixture(scope="module")
def snapped():
    return load(include_overflows=False, snap_tolerance=0.5, max_invert_difference=1.0)


# --- Node types ------------------------------------------------------------------------------


def test_node_types_from_layers(net):
    types = net.nodes["node_type"]
    assert types["M1"] == "maintenance_hole"
    assert types["PS1"] == "pump_station"
    assert types["OF1"] == types["OF2"] == "outfall"
    assert types["CP1"] == "cap"
    assert types["DUP2"] == "DATE"  # unmapped codes pass through
    assert net.nodes.loc["PS1", "name"] == "Test Low Level P.S."


def test_node_types_fall_back_to_inventory_types():
    net = Network(PIPES)
    assert net.nodes["node_type"].to_dict() == {
        **{n: "maintenance_hole" for n in net.graph},
        "PS1": "pump_station",
        "OF1": "outfall",
        "OF2": "outfall",
        "CP1": "cap",
        "DUP2": "DATE",
    }
    assert net.nodes.loc["X1", "layer"] != net.nodes.loc["X1", "layer"]  # NaN: in no layer


def test_nodes_table_and_graph_attributes(net):
    assert set(net.nodes.index) == set(net.graph)
    assert net.nodes.crs == "EPSG:2958"  # the London schema reprojects the UTM 17N layers
    expected = Transformer.from_crs(32617, 2958, always_xy=True).transform(480100.0, 4760100.0)
    assert net.nodes.loc["PS1", "geometry"].coords[0] == pytest.approx(expected)
    assert net.graph.nodes["PS1"]["node_type"] == "pump_station"
    assert net.graph.nodes["M4"]["flow_monitor"]
    assert "rim_elevation" not in net.graph.nodes["M3"]  # missing values are left out
    assert net.graph.edges["PS1", "M3", FORCE_MAIN]["structure_type"] == "force_main"


def test_layer_nodes_without_pipes_are_left_out(net):
    assert "M9" not in net.nodes.index  # only an abandoned pipe


# --- Pipe status and structure type ----------------------------------------------------------


def test_status_filter(net):
    assert not {"M9", "M10", "M11"} & set(net.graph)
    assert net.removed["status"] == {"abandoned": 1, "proposed": 1, "removed": 1}
    assert set(net.full_data["status"]) == {"built"}


@pytest.mark.parametrize(
    "kwargs, kept",
    [
        (dict(status=None), {"M9", "M10", "M11"}),
        (dict(status={"built", "proposed"}), {"M10"}),
        (dict(keep_abandoned=True), {"M9"}),
    ],
)
def test_status_options(kwargs, kept):
    assert {"M9", "M10", "M11"} & set(Network(PIPES, **kwargs).graph) == kept


def test_overflows_included_by_default(net):
    assert net.graph.has_edge("M3", "OF2")
    assert net.subnetwork(net.graph).validate().splits == {"M3": [5003, OVERFLOW]}


def test_exclude_overflows(snapped):
    assert "OF2" not in snapped.graph
    assert snapped.removed["structure_type"] == {"overflow": 1}


# --- Outfalls and validation -----------------------------------------------------------------


def test_outfalls_from_data(net, snapped):
    assert net.subnetwork(net.graph).outfalls == ["OF1", "OF2"]
    assert snapped.catchment("OF1").outfalls == ["OF1"]


def test_dead_ends_classified(net):
    report = net.subnetwork(net.graph).validate(["OF1"])
    assert report.dead_end_types == {
        "CP1": "capped",
        "DUP1": "unexplained",  # until it is snapped to DUP2
        "M7": "unexplained",
        "OF2": "outfall",
    }
    assert report.dead_ends == ["CP1", "DUP1", "M7", "OF2"]


def test_dead_end_at_pump_station(net):
    report = net.subnetwork(["M1", "M2", "PS1"]).validate()
    assert report.dead_end_types == {"PS1": "pump_station"}


def test_force_mains_reported(net):
    assert net.catchment("M4").validate().force_mains == [FORCE_MAIN]


def test_pipe_count_mismatch(net):
    report = net.subnetwork(net.graph).validate()
    assert report.pipe_count_mismatches == {"M4": {"pipes_out": (2, 1)}}


# --- Coincident-node snapping ----------------------------------------------------------------


def test_duplicate_disconnects_catchment(net):
    assert "M8" not in net.catchment("M3").nodes


def test_snap_report(snapped):
    report = snapped.snap_report.set_index(["node_1", "node_2"])
    assert list(report.index) == [("DUP1", "DUP2"), ("X2", "Y1")]
    dup = report.loc[("DUP1", "DUP2")]
    assert dup["distance"] == pytest.approx(0.2, abs=1e-6)
    assert (dup["invert_1"], dup["invert_2"]) == (246.0, 245.95)
    assert dup["snapped"] and dup["kept"] == "DUP1"
    crossing = report.loc[("X2", "Y1")]
    assert crossing["distance"] == pytest.approx(0.3, abs=1e-6)
    assert not crossing["snapped"] and "inverts differ" in crossing["reason"]


def test_snapping_prefers_the_maintenance_hole_id(snapped):
    assert "DUP2" not in snapped.graph
    assert snapped.graph.has_edge("DUP1", "M3")
    assert snapped.aliases == {"DUP2": "DUP1"}
    assert snapped.nodes.loc["DUP1", "node_type"] == "maintenance_hole"


def test_snapping_joins_the_catchment(snapped):
    assert "M8" in snapped.catchment("M3").nodes


def test_merged_id_resolves(snapped):
    assert snapped.subnetwork(["DUP2", "M3"]).nodes == ["DUP1", "M3"]
    assert "M8" in snapped.catchment("DUP2").nodes


def test_crossing_pipes_merge_without_invert_limit():
    net = load(snap_tolerance=0.5)
    assert net.snap_report["snapped"].all()
    assert net.aliases == {"DUP2": "DUP1", "Y1": "X2"}
    assert net.graph.has_edge("X2", "M4", 5010)  # Y1 -> M4 now leaves X2


def test_pairs_joined_by_a_pipe_are_not_merged():
    data = gpd.read_file(PIPES)
    short = data[data["UpstreamID"] == "M8"].assign(
        OBJECTID=9000, UpstreamID="DUP1", DownstreamID="DUP2", Length=0.2
    )
    joined = Network(pd.concat([data, short]), snap_tolerance=0.5, max_invert_difference=1.0)
    assert joined.snap_report.loc[0, "reason"] == "joined by a pipe"
    assert not joined.aliases


def _chain(inverts, pipe_a_to_c):
    """Three nodes A, B, C 0.4 m apart in a row (A~B and B~C within 0.5 m, A~C 0.8 m), each the
    end of a pipe from far away, plus a pipe A -> C if ``pipe_a_to_c`` (#60)."""
    from shapely.geometry import LineString

    x0, y0 = 480000.0, 4760000.0
    xs = {"A": x0, "B": x0 + 0.4, "C": x0 + 0.8}
    rows = [
        (
            i,
            f"U{n}",
            n,
            LineString([(x0 + 50 * i, y0 + 50), (xs[n], y0)]),
            inverts[n] + 0.5,
            inverts[n],
        )
        for i, n in enumerate("ABC", 1)
    ]
    if pipe_a_to_c:
        rows.append((9, "A", "C", LineString([(xs["A"], y0), (xs["C"], y0)]), 10.0, 9.9))
    ids, up, down, geometry, us, ds = zip(*rows, strict=True)
    pipes = gpd.GeoDataFrame(
        {
            "OBJECTID": ids,
            "UpstreamID": up,
            "DownstreamID": down,
            "UpstreamInvert": us,
            "DownstreamInvert": ds,
            "PipeShape": "R",
            "Diameter": 300,
            "FlowType": "SAN",
        },
        geometry=list(geometry),
        crs="EPSG:2958",
    )
    return Network(pipes)


def test_snapping_does_not_merge_nodes_joined_through_a_chain():
    # A~B and B~C are each within tolerance and not joined by a pipe, but A -> C is: merging
    # both pairs would turn pipe 9 into a self-loop A -> A.
    net = _chain({"A": 10.0, "B": 10.0, "C": 9.9}, pipe_a_to_c=True)
    report = net.snap_nodes(0.5)
    assert net.graph.has_edge("A", "C", 9)
    assert not any(u == v for u, v in net.graph.edges())
    assert report["snapped"].tolist() == [True, False]
    assert report.loc[1, "reason"] == "would merge nodes joined by a pipe"


def test_snapping_limits_the_invert_spread_of_a_group():
    # A~B and B~C differ by 0.3 m each, within the 0.5 m limit, but A and C by 0.6 m.
    net = _chain({"A": 10.0, "B": 10.3, "C": 10.6}, pipe_a_to_c=False)
    report = net.snap_nodes(0.5, max_invert_difference=0.5)
    assert report["snapped"].tolist() == [True, False]
    assert "would merge nodes whose inverts differ" in report.loc[1, "reason"]
    assert net.aliases == {"B": "A"}
    # Without the limit, all three merge.
    net = _chain({"A": 10.0, "B": 10.3, "C": 10.6}, pipe_a_to_c=False)
    assert net.snap_nodes(0.5)["snapped"].all()
    assert net.aliases == {"B": "A", "C": "A"}


def test_snapping_in_geographic_crs():
    pipes = gpd.read_file(PIPES).to_crs("EPSG:4326")
    net = Network(pipes, snap_tolerance=0.5, max_invert_difference=1.0)
    assert net.aliases == {"DUP2": "DUP1"}  # no node layers: the smaller ID is kept


# --- ID normalisation ------------------------------------------------------------------------


@pytest.mark.parametrize(
    "node, expected",
    [("OK0329", "OK329"), ("OK329", "OK329"), ("OK0", "OK0"), ("WT6011", "WT6011"), (7, 7)],
)
def test_strip_leading_zeros(node, expected):
    assert strip_leading_zeros(node) == expected


def test_normalize_ids(tmp_path):
    holes = gpd.read_file(LAYERS["maintenance_holes"])
    holes["GIS_FeatureKey"] = holes["GIS_FeatureKey"].str.replace("M", "M0")
    net = Network(PIPES, maintenance_holes=holes, normalize_ids=strip_leading_zeros)
    assert net.nodes.loc["M1", "rim_elevation"] == 252.0
    path = tmp_path / "nodes.txt"
    path.write_text("M001\nM02\n")
    assert Subnetwork.from_node_list(net, path).nodes == ["M1", "M2"]
    assert net.catchment("M02", outlet="M004").nodes[-1] == "PS1"


def test_unnormalized_layer_ids_do_not_match():
    holes = gpd.read_file(LAYERS["maintenance_holes"])
    holes["GIS_FeatureKey"] = holes["GIS_FeatureKey"].str.replace("M", "M0")
    net = Network(PIPES, maintenance_holes=holes)
    assert math.isnan(net.nodes.loc["M1", "rim_elevation"])


# --- Maintenance-hole attributes -------------------------------------------------------------


def test_rim_elevation(net):
    rims = net.nodes["rim_elevation"]
    assert rims["M1"] == 252.0  # LidElevation
    assert rims["M2"] == pytest.approx(248.7 + 3.0)  # invert + Depth
    assert math.isnan(rims["M3"])  # 12 m is outside London's 150-350 m
    assert rims["OF1"] == 243.0  # HeadwallElevation
    assert net.nodes.loc["OF1", "invert"] == 240.0  # PipeInvert, not the pipe's 240.5


def test_maintenance_hole_fields(net):
    m4 = net.nodes.loc["M4"]
    assert (m4["pipes_in"], m4["pipes_out"]) == (3, 2)
    assert m4["flow_monitor"] and m4["verified"] and not m4["has_drop"] and not m4["weir"]
    assert m4["access"] == "Road"


def test_junction_max_depth_from_rim(snapped):
    model = build_model(snapped.catchment("OF1"), max_depth=1.5)
    depths = {j.name: j.max_depth for j in model.junctions}
    assert depths["M1"] == pytest.approx(3.0)
    assert depths["M2"] == pytest.approx(3.0)
    assert depths["M4"] == pytest.approx(251.0 - 244.9)
    assert depths["M3"] == 1.5  # no plausible rim: the default
    assert [o.name for o in model.outfalls] == ["OF1"]


def test_node_layer_needs_node_id_column():
    schema = Schema(
        name="no nodes",
        columns={
            "pipe_id": "OBJECTID",
            "from_node": "UpstreamID",
            "to_node": "DownstreamID",
            "length": "Length",
        },
    )
    with pytest.raises(ValueError, match="node_id column"):
        Network(PIPES, schema=schema, maintenance_holes=LAYERS["maintenance_holes"])


# --- SWMM ------------------------------------------------------------------------------------


def test_written_inp_runs_in_swmm(snapped, tmp_path):
    sub = snapped.catchment("OF1")
    assert sub.validate().is_dendritic and not sub.validate().dead_ends
    path = write_inp(
        sub, tmp_path, routing="DYNWAVE", outfalls=None, pumps={FORCE_MAIN: {"curve": "*"}}
    )
    sections = read_sections(path)
    assert [row[0] for row in sections["OUTFALLS"]] == ["OF1"]
    assert [row[0] for row in sections["PUMPS"]] == [str(FORCE_MAIN)]
    assert run_swmm(path, tmp_path) == []


# --- Node IDs of the .inp options (#59) -----------------------------------------------------


def test_normalize_ids_skips_missing_ids():
    # A plain str method as the hook: the pipes with no node ID ('Abandoned' -> NaN) used to
    # reach it and raise TypeError.
    net = Network(TESTS / "data" / "fixture_network.geojson", normalize_ids=str.upper)
    assert "A1" in net.graph


def test_options_resolve_normalized_ids(tmp_path):
    def upper(node):
        return node.upper() if isinstance(node, str) else node

    net = Network(TESTS / "data" / "fixture_network.geojson", normalize_ids=upper)
    sub = net.subnetwork(["a1", "a2", "j1"])

    def written(case):
        ids = {n: case(n) for n in ("A1", "A2", "J1")}
        path = write_inp(
            sub,
            tmp_path,
            name=f"{case.__name__}.inp",
            outfalls={ids["J1"]: "FREE"},
            dwf={ids["A1"]: 1.0},
            wet_wells={ids["A2"]: dict(max_depth=2.0, area=3.0)},
            height_correction={ids["A1"]: 104.9},
        )
        return read_sections(path)

    assert written(str.lower) == written(str.upper)
    sections = written(str.lower)
    assert sections["OUTFALLS"][0][0] == "J1"
    assert sections["DWF"][0][0] == "A1"
    assert sections["STORAGE"][0][0] == "A2"
    assert {r[0]: r[1] for r in sections["JUNCTIONS"]}["A1"] == "104.9"


def test_options_resolve_snapped_aliases(snapped):
    sub = snapped.subnetwork(["M8", "DUP1", "M3"])
    model = build_model(sub, outfalls=["M3"], dwf={"DUP2": 1.0}, height_correction={"DUP2": 1.0})
    assert [d.node for d in model.dwf] == ["DUP1"]
    assert {j.name: j.invert for j in model.junctions}["DUP1"] == 1.0


def test_height_correction_outside_the_subnetwork(snapped, caplog):
    sub = snapped.subnetwork(["M8", "DUP1", "M3"])
    with caplog.at_level("WARNING", logger="swmmbuilder"):
        model = build_model(sub, outfalls=["M3"], height_correction={"M1": 1.0, "NOPE": 1.0})
    # M1 is in the network but not in the subnetwork: skipped quietly; NOPE is in neither.
    assert "nodes not in the network, ignored: NOPE" in caplog.text
    assert "M1" not in caplog.text
    assert "M1" not in {j.name for j in model.junctions}


# --- Subnetworks made before a snap (#75) ---------------------------------------------------


def test_subnetwork_made_before_a_snap_is_out_of_date(tmp_path):
    net = load()
    sub = net.subnetwork(["M8", "DUP1", "DUP2", "M3"])
    sub.validate()  # fine before the snap
    net.snap_nodes(0.5, max_invert_difference=1.0)  # DUP2 -> DUP1
    for use in (sub.pipes, sub.validate, sub.data_gaps, lambda: sub.to_inp(tmp_path / "x.inp")):
        with pytest.raises(RuntimeError, match=r"out of date.*network\.subnetwork"):
            use()
    rebuilt = net.subnetwork(sub.nodes)  # DUP2 resolves to DUP1
    assert rebuilt.nodes == ["DUP1", "M3", "M8"]
    rebuilt.validate()


def test_snap_that_merges_nothing_keeps_subnetworks():
    net = load()
    sub = net.subnetwork(["M8", "DUP1"])
    net.snap_nodes(0.001)  # nothing that close
    assert not net.aliases
    sub.validate()


def test_subnetwork_made_after_a_snap_is_current():
    net = load(snap_tolerance=0.5, max_invert_difference=1.0)
    net.subnetwork(["M8", "DUP1", "M3"]).validate()
