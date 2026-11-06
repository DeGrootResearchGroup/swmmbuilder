"""Allocation of polygon populations to network nodes (swmmbuilder.population)."""

import logging
import math

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from conftest import FIXTURE
from shapely.geometry import LineString, box

from swmmbuilder import DataWarning, Network, allocate_population, dwf_from_population

CRS = "EPSG:2958"
X0, Y0 = 480_000.0, 4_760_000.0

# A straight line of pipes along y = 50 (x in m from X0), as (pipe, up, down, x_up, x_down,
# PDC_Count). Their lengths are 50, 30, 60, 40 and 60 m; midpoints at x = 35, 75, 120, 170
# and 220.
PIPES = [
    (1, "N1", "N2", 10, 60, 2),
    (2, "N2", "N3", 60, 90, 6),
    (3, "N3", "N4", 90, 150, None),  # missing: imputed
    (4, "N4", "N5", 150, 190, 2),
    (5, "N5", "N6", 190, 250, 999),  # implausible: imputed; midpoint in no polygon
]
# P holds pipes 1 and 2, Q pipes 3 and 4, R none.
POLYGONS = {"P": (0, 100, 80.0), "Q": (100, 200, 100.0), "R": (300, 400, 30.0)}


def make_network(pipes=PIPES):
    """The network of ``pipes``: (pipe, up, down, x_up, x_down, PDC_Count[, y]), along y = 50 m
    unless y is given."""
    template = gpd.read_file(FIXTURE).iloc[[0]].drop(columns="geometry")
    rows = []
    for pipe, up, down, x_up, x_down, pdc, *rest in pipes:
        y = rest[0] if rest else 50
        structure = rest[1] if len(rest) > 1 else "SL"
        row = template.iloc[0].to_dict()
        row.update(
            OBJECTID=pipe,
            UpstreamID=up,
            DownstreamID=down,
            Length=float(x_down - x_up),
            PDC_Count=pdc,
            StructureType=structure,
            geometry=LineString([(X0 + x_up, Y0 + y), (X0 + x_down, Y0 + y)]),
        )
        rows.append(row)
    return Network(gpd.GeoDataFrame(rows, crs=CRS))


def squares(spec=POLYGONS):
    """Polygons 0-100 m in y, x from x0 to x1, indexed by name, with a population column."""
    return gpd.GeoDataFrame(
        {"population": [pop for _, _, pop in spec.values()]},
        geometry=[box(X0 + x0, Y0, X0 + x1, Y0 + 100) for x0, x1, _ in spec.values()],
        index=pd.Index(list(spec), name="id"),
        crs=CRS,
    )


@pytest.fixture
def net():
    return make_network()


def allocate(net, polygons=None, **options):
    options = {"weight": "PDC_Count", "max_weight": 100, **options}
    return allocate_population(net, squares() if polygons is None else polygons, **options)


def approx(d):
    return pytest.approx(d, abs=1e-9)


def test_shared_by_weight_at_upstream_nodes(net):
    # Q's only valid weight per metre is pipe 4's, 2/40, so pipe 3 weighs 0.05 x 60 = 3. Pipe 5
    # is in no polygon, so has no rate: weight 0. P (weights 2 and 6) and Q (3 and 2).
    alloc = allocate(net)
    assert alloc.nodes == approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})
    assert alloc.pipes["polygon"].tolist()[:4] == ["P", "P", "Q", "Q"]
    assert pd.isna(alloc.pipes.loc[5, "polygon"])
    assert alloc.pipes["weight"].tolist() == approx([2, 6, 3, 2, 0])
    assert alloc.pipes["weight_imputed"].tolist() == [False, False, True, False, True]
    assert alloc.total == pytest.approx(180)
    assert alloc.unallocated == pytest.approx(30)
    assert alloc.no_pipes == ["R"]
    assert alloc.zero_weight == []
    assert alloc.polygons.loc["Q"].to_dict() == approx(
        {"population": 100, "pipes": 2, "weight": 5, "allocated": 100}
    )
    assert alloc.polygons.loc["R", "allocated"] == 0


def test_feeds_dwf_from_population(net):
    dwf = dwf_from_population(allocate(net).nodes, per_capita=86_400, unit="L/d")
    assert dwf == approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})  # L/s


@pytest.mark.parametrize(
    "placement, expected",
    [
        ("downstream", {"N2": 20, "N3": 60, "N4": 60, "N5": 40, "N6": 0}),
        ("split", {"N1": 10, "N2": 40, "N3": 60, "N4": 50, "N5": 20, "N6": 0}),
    ],
)
def test_placement(net, placement, expected):
    assert allocate(net, placement=placement).nodes == approx(expected)


def test_length_weights_by_default(net):
    alloc = allocate_population(net, squares())
    assert alloc.nodes == approx({"N1": 50, "N2": 30, "N3": 60, "N4": 40, "N5": 0})
    assert not alloc.pipes["weight_imputed"].any()


def test_callable_weight(net):
    alloc = allocate(net, weight=lambda pipes: 10 * pipes["PDC_Count"], max_weight=1000)
    assert alloc.nodes == approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})


def test_imputed_from_the_pipes_own_polygon():
    # P: pipe 1 has 3 / 30 per metre, so pipe 2 (20 m, missing) weighs 2; the network median
    # (of 0.1 and Q's 0.05) would give it 1.5.
    pipes = [
        (1, "N1", "N2", 10, 40, 3),
        (2, "N2", "N3", 40, 60, None),
        (3, "N4", "N5", 110, 150, 2),
    ]
    alloc = allocate(make_network(pipes))
    assert alloc.pipes["weight"].tolist() == approx([3, 2, 2])
    assert alloc.nodes == approx({"N1": 48, "N2": 32, "N4": 100})


def test_polygon_with_no_valid_weight_is_reported():
    # Q's only pipe has no count: with no fallback, nothing to impute from, so Q's people are not
    # allocated.
    pipes = [(1, "N1", "N2", 10, 60, 2), (3, "N3", "N4", 110, 150, None)]
    alloc = allocate(make_network(pipes), fallback_weight=None)
    assert alloc.no_valid_weight == ["Q"]
    assert alloc.zero_weight == []
    assert alloc.pipes.loc[3, "weight"] == 0
    assert alloc.nodes == approx({"N1": 80, "N3": 0})
    assert alloc.unallocated == pytest.approx(130)


def test_impute_zero(net):
    alloc = allocate(net, impute="zero")
    assert alloc.pipes["weight"].tolist() == approx([2, 6, 0, 2, 0])
    assert alloc.nodes == approx({"N1": 20, "N2": 60, "N3": 0, "N4": 100, "N5": 0})


def test_partly_sewered_and_zero_weight_polygons():
    # S holds only pipe 4, whose weight is 0; Q is only partly sewered (its east half has no pipe).
    pipes = [(1, "N1", "N2", 110, 140, 3), (4, "N4", "N5", 210, 240, 0)]
    alloc = allocate(
        make_network(pipes),
        squares({"Q": (100, 200, 50.0), "S": (200, 300, 7.0)}),
        fallback_weight=None,
    )
    assert alloc.nodes == approx({"N1": 50, "N4": 0})
    assert alloc.zero_weight == ["S"]
    assert alloc.unallocated == pytest.approx(7)


def test_fallback_where_the_median_is_zero():
    # P: counts 0, 0 and 2 (median per metre 0), and pipe 4 missing: it gets one connection.
    pipes = [
        (1, "N1", "N2", 0, 20, 0),
        (2, "N2", "N3", 20, 40, 0),
        (3, "N3", "N4", 40, 60, 2),
        (4, "N4", "N5", 60, 90, None),
    ]
    alloc = allocate(make_network(pipes), squares({"P": (0, 100, 30.0)}))
    assert alloc.pipes["weight"].tolist() == approx([0, 0, 2, 1])
    assert alloc.pipes["weight_fallback"].tolist() == [False, False, False, True]
    assert alloc.nodes == approx({"N1": 0, "N2": 0, "N3": 20, "N4": 10})
    assert alloc.fallback == ["P"]


def test_fallback_in_a_polygon_with_no_connections():
    # Every count in P is 0: each pipe gets one connection, so P is shared equally.
    pipes = [(1, "N1", "N2", 0, 20, 0), (2, "N2", "N3", 20, 80, 0), (3, "N3", "N4", 80, 90, 0)]
    alloc = allocate(make_network(pipes), squares({"P": (0, 100, 30.0)}))
    assert alloc.pipes["weight"].tolist() == approx([1, 1, 1])
    assert alloc.nodes == approx({"N1": 10, "N2": 10, "N3": 10})
    assert alloc.zero_weight == [] and alloc.unallocated == pytest.approx(0)


def test_fallback_with_no_valid_weight():
    # Zero valid counts and two missing ones: the missing get one connection each, the zeros
    # keep 0.
    pipes = [
        (1, "N1", "N2", 0, 20, 0),
        (2, "N2", "N3", 20, 80, None),
        (3, "N3", "N4", 80, 90, None),
    ]
    alloc = allocate(make_network(pipes), squares({"P": (0, 100, 30.0)}))
    assert alloc.pipes["weight"].tolist() == approx([0, 1, 1])
    assert alloc.nodes == approx({"N1": 0, "N2": 15, "N3": 15})
    # Only missing counts: every pipe gets one.
    pipes = [(2, "N2", "N3", 20, 80, None), (3, "N3", "N4", 80, 90, None)]
    alloc = allocate(make_network(pipes), squares({"P": (0, 100, 30.0)}), fallback_weight=2)
    assert alloc.pipes["weight"].tolist() == approx([2, 2])
    assert alloc.no_valid_weight == []


def test_missing_population_counts_as_zero(net):
    polygons = squares()
    polygons.loc["P", "population"] = None
    alloc = allocate(net, polygons, population=polygons["population"].astype(object))
    assert alloc.nodes["N1"] == 0 and alloc.polygons.loc["P", "population"] == 0
    assert alloc.no_pipes == ["R"]


def test_within_counts_only_the_subnetworks_pipes(net):
    alloc = allocate(net)
    sub = net.catchment("N4")  # pipes 1-3
    assert alloc.within(sub) == approx({"N1": 20, "N2": 60, "N3": 60})
    # Allocating over the subnetwork alone (here by length) puts all of Q's people on pipe 3.
    alone = allocate_population(net.catchment("N4"), squares())
    assert alone.nodes == approx({"N1": 50, "N2": 30, "N3": 100})


def test_overlapping_polygons_take_the_first(net):
    polygons = squares({"P": (0, 100, 80.0), "Big": (0, 200, 10.0)})
    with pytest.warns(DataWarning, match="2 pipe midpoints are in more than one polygon"):
        alloc = allocate(net, polygons)
    assert alloc.pipes["polygon"].tolist()[:4] == ["P", "P", "Big", "Big"]


def test_reprojects_polygons(net):
    alloc = allocate(net, squares().to_crs("EPSG:4326"))
    assert alloc.nodes == pytest.approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})


def test_invalid_options_listed_together(net):
    with pytest.raises(ValueError) as error:
        allocate(net, weight="Connections", placement="middle", impute="mean")
    message = str(error.value)
    for name in ("weight", "placement", "impute"):
        assert name in message


def test_repeated_ids_and_missing_crs(net):
    polygons = squares()
    polygons.index = ["P", "P", "R"]
    polygons = polygons.set_crs(None, allow_override=True)
    with pytest.raises(ValueError, match="repeated") as error:
        allocate(net, polygons)
    assert "polygons has no CRS" in str(error.value)


def test_negative_population(net):
    polygons = squares()
    polygons.loc["Q", "population"] = -1
    with pytest.raises(ValueError, match=r"Negative population in polygons \['Q'\]"):
        allocate(net, polygons)


def test_shares_sum_to_each_polygon(net):
    alloc = allocate(net, placement="split")
    sewered = alloc.polygons[alloc.polygons["pipes"] > 0]
    assert sewered["allocated"].tolist() == approx(sewered["population"].tolist())
    assert math.isclose(alloc.total + alloc.unallocated, squares()["population"].sum())


def blocks(spec):
    """Rectangles {name: (x0, x1, y0, y1, population)}, in m from (X0, Y0)."""
    return gpd.GeoDataFrame(
        {"population": [b[4] for b in spec.values()]},
        geometry=[box(X0 + x0, Y0 + y0, X0 + x1, Y0 + y1) for x0, x1, y0, y1, _ in spec.values()],
        index=list(spec),
        crs=CRS,
    )


# Two blocks either side of the street the pipes run under (y = 50), and one far away.
STREET = {"N": (0, 100, 55, 100, 60.0), "S": (0, 100, 0, 45, 40.0), "E": (300, 400, 0, 100, 9.0)}
STREET_PIPES = [(1, "N1", "N2", 10, 60, 2), (2, "N2", "N3", 60, 90, 6)]


def test_distance_rule_serves_the_blocks_on_both_sides():
    net = make_network(STREET_PIPES)
    alloc = allocate(net, blocks(STREET), max_distance=10)
    # Each block shares by 2 : 6 over both pipes: N 15 and 45, S 10 and 30.
    assert alloc.nodes == approx({"N1": 25, "N2": 75})
    assert alloc.pipes.loc[1, "polygon"].tolist() == ["N", "S"]
    assert alloc.by_pipe().to_dict() == approx({1: 25, 2: 75})
    assert alloc.polygons.loc[["N", "S"], "pipes"].tolist() == [2, 2]
    assert alloc.no_pipes == ["E"]
    assert alloc.unallocated == pytest.approx(9)


def test_distance_rule_leaves_blocks_with_no_pipe_that_close():
    net = make_network(STREET_PIPES)
    alloc = allocate(net, blocks(STREET), max_distance=3)
    assert alloc.no_pipes == ["N", "S", "E"]
    assert alloc.total == 0
    # The midpoint rule finds no pipe in the blocks either.
    assert allocate(net, blocks(STREET)).no_pipes == ["N", "S", "E"]


def test_distance_rule_imputes_from_each_blocks_own_pipes():
    # Pipe 2 has no count and serves N (with pipe 1: 2 / 50 per metre) and S (with pipe 1 and pipe
    # 3, 9 / 30 per metre, at y = 20: median 0.17), so it weighs 0.04 x 30 in N and 0.17 x 30 in S.
    pipes = [
        (1, "N1", "N2", 10, 60, 2),
        (2, "N2", "N3", 60, 90, None),
        (3, "M1", "M2", 60, 90, 9, 20),
    ]
    net = make_network(pipes)
    alloc = allocate(net, blocks(STREET), max_distance=10)
    weights = alloc.pipes.reset_index().set_index(["pipe_id", "polygon"])["weight"]
    assert weights[(2, "N")] == pytest.approx(1.2)
    assert weights[(2, "S")] == pytest.approx(5.1)
    assert weights[(3, "S")] == 9 and (3, "N") not in weights.index


def test_negative_distance(net):
    with pytest.raises(ValueError, match="max_distance must be a distance"):
        allocate(net, max_distance=-1)


def test_force_mains_take_no_connections():
    # Pipe 2 is a force main with a count of 6, in the street between N and S: it gets weight 0
    # and no one, so both blocks go to pipe 1.
    pipes = [(1, "N1", "N2", 10, 60, 2), (2, "N2", "N3", 60, 90, 6, 50, "FM")]
    net = make_network(pipes)
    assert net.full_data.set_index("pipe_id").loc[2, "structure_type"] == "force_main"
    alloc = allocate(net, blocks(STREET), max_distance=10)
    assert alloc.by_pipe().to_dict() == approx({1: 100, 2: 0})
    assert alloc.pipes.loc[2, "weight"] == 0 and pd.isna(alloc.pipes.loc[2, "polygon"])
    # Nor from the fallback: a block next to a force main alone gets no one.
    alone = allocate(
        make_network([(2, "N2", "N3", 60, 90, None, 50, "FM")]), blocks(STREET), max_distance=10
    )
    assert alone.pipes["weight"].tolist() == [0] and alone.total == 0
    assert alone.no_pipes == ["N", "S", "E"]
    # Unless asked.
    both = allocate(net, blocks(STREET), max_distance=10, serve_force_mains=True)
    assert both.nodes == approx({"N1": 25, "N2": 75})


def test_sewered_polygons_go_to_the_nearest_pipe():
    # N and S have no pipe within 3 m. S is flagged sewered and goes to its nearest pipe: the
    # street pipes are 5 m away (pipe 1 on the tie), pipe 3 is 15 m below it, and the force main
    # inside it does not count. N and E are not flagged, and contribute no one.
    pipes = [
        (1, "N1", "N2", 10, 60, 2),
        (2, "N2", "N3", 60, 90, 6),
        (3, "M1", "M2", 20, 40, 1, -15),
        (4, "F1", "F2", 0, 100, 0, 44, "FM"),
    ]
    polygons = blocks(STREET).assign(residential=[False, True, False])
    alloc = allocate(make_network(pipes), polygons, max_distance=3, sewered="residential")
    assert alloc.nodes == approx({"N1": 40, "N2": 0, "M1": 0, "F1": 0})
    assert alloc.nearest.to_dict() == approx({"S": 5})
    assert alloc.no_pipes == ["N", "E"]
    assert alloc.unallocated == pytest.approx(69)


def test_nearest_pipe_within_a_cap():
    # S's nearest pipe is 5 m away, E's 210 m: with a 100 m cap E still contributes no one.
    polygons = blocks(STREET).assign(residential=[False, True, True])
    net = make_network(STREET_PIPES)
    capped = allocate(
        net, polygons, max_distance=3, sewered="residential", max_nearest_distance=100
    )
    assert capped.nearest.to_dict() == approx({"S": 5})
    assert capped.no_pipes == ["N", "E"]
    free = allocate(net, polygons, max_distance=3, sewered="residential")
    assert free.nearest.to_dict() == approx({"S": 5, "E": 210})
    assert free.no_pipes == ["N"]
    with pytest.raises(ValueError, match="max_nearest_distance must be"):
        allocate(net, polygons, max_nearest_distance=-5)


@pytest.mark.parametrize("flags", [[False, False, True], np.array([False, False, True])])
def test_sewered_rejects_unindexed_flags(net, flags):
    with pytest.raises(TypeError, match=r"sewered.*column name or pandas.Series"):
        allocate(net, sewered=flags)


def test_sewered_series_aligns_by_polygon_id(net):
    flags = pd.Series([True, False, False], index=["R", "Q", "P"])
    result = allocate(net, sewered=flags)
    assert result.nearest.index.tolist() == ["R"]
    assert result.unallocated == 0


@pytest.mark.parametrize("ids", [["P", "Q"], [0, 1, 2]])
def test_sewered_rejects_missing_polygon_ids(net, ids):
    with pytest.raises(ValueError, match=r"sewered.*missing polygon IDs.*R"):
        allocate(net, sewered=pd.Series(True, index=ids))


# Pipe 2 (N2 to nowhere, 30 m, 6 connections) is in full_data but not in the graph.
NO_NODE = [PIPES[0], (2, "N2", None, 60, 90, 6), *PIPES[2:]]


def test_a_pipe_with_a_missing_node_is_left_out(caplog):
    """Its polygon's people go to the polygon's other pipes, and the people it would have taken
    are logged (#122)."""
    with caplog.at_level(logging.WARNING, logger="swmmbuilder.population"):
        alloc = allocate(make_network(NO_NODE))
    assert 2 not in alloc.pipes.index
    assert alloc.nodes["N1"] == pytest.approx(80.0)  # all of P, on pipe 1
    assert alloc.unallocated == pytest.approx(30.0)  # R only
    assert alloc.total == pytest.approx(180.0)
    # With pipe 2, P's 80 people go 2:6 to pipes 1 and 2: 60 people to a pipe in no catchment
    # (placed at N2, its upstream node, so not unallocated).
    assert "1 pipes (30 m) have a missing from_node or to_node" in caplog.text
    assert "they would take 60 people, in no catchment" in caplog.text
    assert "(not allocated: 30 with them, 30 without): [2]" in caplog.text
    # Placed downstream, pipe 2's share has no node, so it was not allocated.
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="swmmbuilder.population"):
        down = allocate(make_network(NO_NODE), placement="downstream")
    assert down.nodes["N2"] == pytest.approx(80.0)
    assert "(not allocated: 90 with them, 30 without)" in caplog.text


def test_no_pipe_is_left_out_of_a_complete_network(net, caplog):
    with caplog.at_level(logging.WARNING, logger="swmmbuilder.population"):
        allocate(net)
    assert "left out" not in caplog.text
