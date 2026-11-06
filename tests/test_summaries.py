"""summary() of the allocation, land-use, topology, and model results (#137)."""

import math

import pandas as pd
import pytest
import test_land_use
import test_population
from conftest import DENDRITIC, subnetwork, with_pipes

from swmmbuilder import allocate_land_use, swmm
from swmmbuilder.design import Derived
from swmmbuilder.swmm.records import Conduit, Junction, Outfall, PumpLink, SwmmModel

POPULATION_KEYS = [
    "total",
    "nodes",
    "weight",
    "allocation_total",
    "unallocated",
    "polygons",
    "polygons_served",
    "nearest",
    "nearest_population",
    "nearest_max_distance",
    "no_pipes",
    "no_pipes_population",
    "no_valid_weight",
    "zero_weight",
    "fallback",
]


@pytest.fixture
def net():
    return test_population.make_network()


def test_population_summary(net):
    # Weights 2, 6, 3 (imputed), 2, and 0 (pipe 5, in no polygon); R (30 people) has no pipe.
    s = test_population.allocate(net).summary()
    assert list(s) == POPULATION_KEYS
    assert s == {
        "total": 180.0,
        "nodes": 4,  # N5 gets 0
        "weight": 13.0,
        "allocation_total": 180.0,
        "unallocated": 30.0,
        "polygons": 3,
        "polygons_served": 2,
        "nearest": 0,
        "nearest_population": 0.0,
        "nearest_max_distance": s["nearest_max_distance"],
        "no_pipes": 1,
        "no_pipes_population": 30.0,
        "no_valid_weight": 0,
        "zero_weight": 0,
        "fallback": 0,
    }
    assert math.isnan(s["nearest_max_distance"])
    assert type(s["nodes"]) is int and type(s["total"]) is float
    assert f"{s['total']:,.0f} at {s['nodes']} nodes" == "180 at 4 nodes"


def test_population_summary_within_a_subnetwork(net):
    s = test_population.allocate(net).summary(net.catchment("N4"))  # pipes 1-3
    assert (s["total"], s["nodes"], s["weight"]) == (140.0, 3, 11.0)
    # The rest is of the whole allocation.
    assert (s["allocation_total"], s["no_pipes"]) == (180.0, 1)


def test_population_summary_of_nearest_polygons(net):
    # R, flagged sewered, takes pipe 5, 50 m away, which gets the fallback weight 1.
    polygons = test_population.squares()
    s = test_population.allocate(
        net, polygons, sewered=polygons.index.to_series().eq("R")
    ).summary()
    assert (s["nearest"], s["nearest_population"], s["nearest_max_distance"]) == (1, 30.0, 50.0)
    assert (s["no_pipes"], s["no_pipes_population"], s["unallocated"]) == (0, 0.0, 0.0)
    assert (s["fallback"], s["weight"], s["total"]) == (1, 14.0, 210.0)


def test_population_summary_takes_the_mean_weight_of_a_shared_pipe(net):
    # Within 5 m, pipe 3 serves P and Q. Its missing count is imputed in each: 60 m at P's median
    # 0.12 per metre (of 2/50 and 6/30) is 7.2; at Q's 0.05 (2/40) it is 3. Pipe 5 serves Q
    # (999 is invalid: 60 x 0.05 = 3). The mean counts pipe 3 once, as 5.1.
    alloc = test_population.allocate(net, max_distance=5)
    assert alloc.pipes.loc[3, "weight"].tolist() == pytest.approx([7.2, 3.0])
    s = alloc.summary()
    assert s["weight"] == pytest.approx(2 + 6 + 5.1 + 2 + 3)
    assert s["weight"] != pytest.approx(alloc.pipes["weight"].sum())


def test_land_use_summary(net):
    land = allocate_land_use(net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES)
    s = land.summary()
    assert list(s) == [
        "COMMERCIAL_polygons",
        "COMMERCIAL_area_ha",
        "COMMERCIAL_flow",
        "COMMERCIAL_nodes",
        "INDUSTRIAL_polygons",
        "INDUSTRIAL_area_ha",
        "INDUSTRIAL_flow",
        "INDUSTRIAL_nodes",
        "total",
        "nodes",
        "weight",
        "allocation_total",
        "polygons",
        "polygons_served",
        "area_ha",
        "nearest",
        "nearest_median_distance",
        "nearest_max_distance",
        "no_pipes",
        "zero_weight",
        "unallocated_design_flow",
        "scale",
    ]
    # P (commercial, 0.1 L/s at N1 and N2) and Q (industrial, 0.2 L/s at N3 and N4) are served;
    # R (commercial) has no pipe.
    expected = {
        "COMMERCIAL_polygons": 1,
        "COMMERCIAL_area_ha": 1.0,
        "COMMERCIAL_flow": 0.1,
        "COMMERCIAL_nodes": 2,
        "INDUSTRIAL_polygons": 1,
        "INDUSTRIAL_area_ha": 1.0,
        "INDUSTRIAL_flow": 0.2,
        "INDUSTRIAL_nodes": 2,
        "total": 0.3,
        "nodes": 4,
        "allocation_total": 0.3,
        "polygons": 3,
        "polygons_served": 2,
        "area_ha": 2.0,
        "nearest": 0,
        "no_pipes": 1,
        "zero_weight": 0,
        "unallocated_design_flow": 0.1,
        "scale": 1.0,
    }
    assert {k: s[k] for k in expected} == pytest.approx(expected)
    assert math.isnan(s["nearest_median_distance"]) and math.isnan(s["nearest_max_distance"])


def test_land_use_summary_with_nearest_and_total(net):
    # R takes pipe 5, 50 m away; the 0.4 L/s placed is scaled to 0.8.
    land = allocate_land_use(
        net,
        test_land_use.zoning(),
        land_use="use",
        rates=test_land_use.RATES,
        to_nearest=True,
        total=0.8,
    )
    s = land.summary()
    assert (s["nearest"], s["nearest_median_distance"], s["nearest_max_distance"]) == (
        1,
        50.0,
        50.0,
    )
    assert (s["COMMERCIAL_polygons"], s["COMMERCIAL_nodes"]) == (2, 3)
    assert [s["COMMERCIAL_flow"], s["scale"], s["total"]] == pytest.approx([0.4, 2.0, 0.8])


def test_topology_summary(network):
    sub = network.subnetwork(network.nodes.index)
    report = sub.validate()
    s = report.summary()
    assert (s["nodes"], s["pipes"]) == (sub.graph.number_of_nodes(), sub.graph.number_of_edges())
    for key in (
        "splits",
        "cycles",
        "dead_ends",
        "parallel_pipes",
        "force_mains",
        "stations_without_force_main",
        "orphan_force_mains",
        "force_main_ends",
        "pipe_count_mismatches",
    ):
        assert s[key] == len(getattr(report, key))
    assert s["data_gaps"] == sum(report.data_gaps.values())
    assert all(type(v) is int for v in s.values())


def test_topology_summary_counts_a_cycle():
    net = with_pipes(("A2", "A1"))
    sub = net.subnetwork(["A1", "A2", "J1"])
    s = sub.validate(outfalls=["J1"]).summary()
    assert (s["nodes"], s["pipes"], s["cycles"], s["splits"]) == (3, 3, 1, 1)


def test_model_summary_of_records():
    gaps = pd.DataFrame({"gap": ["a", "b", "c"], "filled_from": ["size_fill", None, "x"]})
    model = SwmmModel(
        junctions=[Junction(n, 0.0, 2.0, 0.0, 0.0, 0.0) for n in ("A", "S", "J")],
        outfalls=[Outfall("J_OUT", -0.3, "FREE"), Outfall("O2", 0.0, "FREE")],
        conduits=[
            Conduit("C1", "A", "S", 10, 0.013, 0, 0, 0, 0),
            Conduit("J_DUMMY", "J", "J_OUT", 30, 0.01, 0, 0, 0, 0),
        ],
        pumps=[
            PumpLink("S_P1", "S", "J", "*", "ON", 0, 0),
            PumpLink("S_P2", "S", "J", "*", "ON", 0, 0),
        ],
        data_gaps=gaps,
        derived=[
            Derived("junction", "A", "max_depth", 2.0, "max_depth='median'", ""),
            Derived("junction", "J", "max_depth", 2.0, "max_depth='median'", ""),
            Derived("junction", "K", "max_depth", 1.7, "rim estimate", ""),
            Derived("station", "S", "invert", -1.0, "station_invert='lowest_inlet'", ""),
        ],
    )
    assert model.summary() == {
        "junctions": 3,
        "outfalls": 2,
        "storage": 0,
        "conduits": 2,
        "pumps": 2,
        "stations": 1,
        "stations_extra": 0,
        "stations_excluded": 0,
        "outfall_links": 1,
        "data_gaps": 3,
        "data_gaps_filled": 2,
        "derived": 4,
        "derived_junction_max_depth": 3,
        "derived_junction_max_depth[max_depth='median']": 2,
        "derived_junction_max_depth[rim estimate]": 1,
        "derived_station_invert": 1,
        "derived_station_invert[station_invert='lowest_inlet']": 1,
        "default_depth_junctions": 2,
        "default_depth": 2.0,
    }


def test_model_summary_of_an_empty_model():
    s = SwmmModel().summary()
    assert (s["data_gaps"], s["derived"], s["default_depth_junctions"]) == (0, 0, 0)
    assert math.isnan(s["default_depth"])
    assert not any(k.startswith("derived_") for k in s)


def test_model_summary_of_a_built_model(network):
    # Rims 2, 3, and 4 m above the inverts of A1, A2, and B1 (as in test_median_depth).
    sub = subnetwork(network, DENDRITIC)
    inverts = swmm.node_inverts(sub)
    for node, height in (("A1", 2.0), ("A2", 3.0), ("B1", 4.0)):
        sub.graph.nodes[node]["rim_elevation"] = inverts[node] + height
    model = swmm.build_model(sub, ["OUT"], max_depth="median")
    s = model.summary()
    # B2, C1, J1, and J2 take the median depth, 3 m, of A1, A2, and B1.
    assert (s["default_depth_junctions"], s["default_depth"]) == (4, pytest.approx(3.0))
    assert s["derived_junction_max_depth[max_depth='median']"] == 4
    assert s["data_gaps"] == len(model.data_gaps)
    assert (s["outfall_links"], s["outfalls"]) == (0, 1)


def test_model_summary_counts_outfall_links(network):
    sub = subnetwork(network, [e for e in DENDRITIC if "OUT" not in e and "J2" not in e])
    s = swmm.build_model(sub, ["J1"], outfall_link=True, max_depth=2.5).summary()
    assert (s["outfall_links"], s["outfalls"], s["junctions"]) == (1, 1, 6)


# ---- one shape and one key style for every result (#161) -----------------------------------------


def test_land_use_summary_within_a_subnetwork(net):
    land = allocate_land_use(net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES)
    sub = net.catchment("N4")  # pipes 1-3: P's, and Q's pipe 3 (0.12 L/s at N3)
    s = land.summary(sub)
    assert [s["total"], s["COMMERCIAL_flow"], s["INDUSTRIAL_flow"]] == pytest.approx(
        [0.22, 0.1, 0.12]
    )
    assert (s["nodes"], s["COMMERCIAL_nodes"], s["INDUSTRIAL_nodes"]) == (3, 2, 1)
    assert s["weight"] == pytest.approx(50 + 30 + 60)  # lengths
    assert s["allocation_total"] == pytest.approx(0.3)  # the rest is of the whole allocation
    assert (s["polygons_served"], s["no_pipes"]) == (2, 1)
    assert s["total"] == pytest.approx(land.within(sub).to_numpy().sum())


def test_land_use_by_pipe(net):
    land = allocate_land_use(net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES)
    flows = land.by_pipe()
    assert flows.index.name == "pipe_id"
    assert flows.to_dict() == pytest.approx({1: 0.0625, 2: 0.0375, 3: 0.12, 4: 0.08, 5: 0.0})
    population = test_population.allocate(net).by_pipe()
    assert list(flows.index) == list(population.index)


def test_within_returns_the_type_of_nodes(net):
    sub = net.catchment("N4")
    people = test_population.allocate(net)
    land = allocate_land_use(net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES)
    assert type(people.within(sub)) is type(people.nodes)
    assert type(land.within(sub)) is type(land.nodes)
    assert list(land.within(sub).columns) == list(land.nodes.columns)


def test_land_use_summary_counts_a_polygon_with_zero_weight_once():
    # Q's pipes (3 and 4) weigh 0, so Q is served but unallocated, and R has no pipe. Q used to
    # be counted both in polygons_served and in no_pipes.
    pipes = [
        (p, u, d, x0, x1, 0 if p in (3, 4) else 1) for p, u, d, x0, x1, _ in test_population.PIPES
    ]
    net = test_population.make_network(pipes)
    land = allocate_land_use(
        net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES, weight="PDC_Count"
    )
    assert sorted(land.no_pipes) == ["Q", "R"]
    s = land.summary()
    assert (s["polygons"], s["polygons_served"], s["no_pipes"], s["zero_weight"]) == (3, 2, 1, 1)
    assert s["unallocated_design_flow"] == pytest.approx(0.1 + 0.2)


def test_former_keys_still_read_with_a_warning(net):
    s = test_population.allocate(net).summary()
    with pytest.warns(DeprecationWarning, match="'population' is now 'total'"):
        assert s["population"] == s["total"]
    with pytest.warns(DeprecationWarning, match="'population_total' is now 'allocation_total'"):
        assert s["population_total"] == s["allocation_total"]
    assert "population" not in s and "population" not in list(s)
    land = allocate_land_use(net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES)
    with pytest.warns(DeprecationWarning, match="'unallocated' is now 'unallocated_design_flow'"):
        assert land.summary()["unallocated"] == pytest.approx(0.1)
    with pytest.raises(KeyError):
        s["no_such_key"]


def test_summaries_are_flat_dicts_of_numbers(net, network):
    sub = net.catchment("N4")
    people = test_population.allocate(net)
    results = [
        people.summary(),
        people.summary(sub),
        allocate_land_use(
            net, test_land_use.zoning(), land_use="use", rates=test_land_use.RATES
        ).summary(sub),
        sub.accumulate(people.within(sub)).summary(),
        sub.validate().summary(),
        network.subnetwork(network.nodes.index).find_pump_stations().summary(),
    ]
    for s in results:
        assert isinstance(s, dict)
        assert all(isinstance(k, str) for k in s)
        assert all(type(v) in (int, float) for v in s.values())
        assert dict(s) == s


def test_upstream_totals_summary(net):
    sub = net.catchment("N4")
    totals = sub.accumulate(test_population.allocate(net).within(sub))
    assert totals.summary() == {
        "nodes": 4,
        "pipes": 3,
        "shared_pipes": 0,
        "max": pytest.approx(140.0),
    }
    zones = pd.DataFrame({"a": {"N1": 1.0}, "b": {"N2": 2.0}})
    s = sub.accumulate(zones).summary()
    assert (s["a_max"], s["b_max"]) == (1.0, 2.0) and "max" not in s


def test_pump_station_report_summary(network):
    report = network.subnetwork(network.nodes.index).find_pump_stations()
    assert report.summary() == {
        "stations": len(report.stations),
        "no_force_main": len(report.no_force_main),
        "unlisted": len(report.unlisted),
    }
