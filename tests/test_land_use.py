"""Land-use flows (swmmbuilder.land_use) and pattern inflows (#99)."""

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box
from test_population import CRS, X0, Y0, make_network

from swmmbuilder import (
    LandUseAllocation,
    allocate_land_use,
    hourly_series,
    pattern_inflows,
)

# Pipes along y = 50 (test_population.PIPES): 1 (N1-N2, 50 m, midpoint x = 35), 2 (N2-N3, 30 m,
# x = 75), 3 (N3-N4, 60 m, x = 120), 4 (N4-N5, 40 m, x = 170), 5 (N5-N6, 60 m, x = 220).
# Zones of 1 ha (100 x 100 m): P holds pipes 1 and 2, Q pipes 3 and 4, R none.
ZONES = {"P": (0, 100, "COMMERCIAL"), "Q": (100, 200, "INDUSTRIAL"), "R": (300, 400, "COMMERCIAL")}
# 8640 L/ha/d = 0.1 L/s per ha.
RATES = {"COMMERCIAL": 8640.0, "INDUSTRIAL": 17280.0}


def zoning(spec=ZONES):
    return gpd.GeoDataFrame(
        {"use": [u for *_, u in spec.values()]},
        geometry=[box(X0 + x0, Y0, X0 + x1, Y0 + 100) for x0, x1, _ in spec.values()],
        index=pd.Index(list(spec), name="id"),
        crs=CRS,
    )


@pytest.fixture
def net():
    return make_network()


def approx(x):
    return pytest.approx(x, abs=1e-12)


def test_area_times_rate_shared_by_length(net):
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES)
    assert isinstance(flows, LandUseAllocation)
    # P: 0.1 L/s over 50 and 30 m; Q: 0.2 L/s over 60 and 40 m; at the upstream nodes.
    assert flows.nodes.to_dict("index") == {
        "N1": approx({"COMMERCIAL": 0.0625, "INDUSTRIAL": 0.0}),
        "N2": approx({"COMMERCIAL": 0.0375, "INDUSTRIAL": 0.0}),
        "N3": approx({"COMMERCIAL": 0.0, "INDUSTRIAL": 0.12}),
        "N4": approx({"COMMERCIAL": 0.0, "INDUSTRIAL": 0.08}),
    }
    assert list(flows.nodes.columns) == ["COMMERCIAL", "INDUSTRIAL"]
    assert flows.total == pytest.approx(0.3)
    assert flows.no_pipes == ["R"]
    assert flows.unallocated == pytest.approx(0.1)
    assert flows.scale == 1.0
    p = flows.polygons.loc["P"]
    assert (p["land_use"], p["pipes"]) == ("COMMERCIAL", 2)
    assert [p["area_ha"], p["rate"], p["design_flow"], p["allocated"]] == approx(
        [1.0, 8640.0, 0.1, 0.1]
    )
    assert flows.pipes.loc[5, "flow"] == 0.0 and pd.isna(flows.pipes.loc[5, "polygon"])


def test_scaled_to_total(net):
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES, total=0.6)
    assert flows.scale == pytest.approx(2.0)
    assert flows.total == pytest.approx(0.6)
    assert flows.nodes.loc["N3", "INDUSTRIAL"] == pytest.approx(0.24)
    assert flows.polygons["allocated"].to_dict() == approx({"P": 0.2, "Q": 0.4, "R": 0.0})


def test_area_within_clips_the_polygons(net):
    half_of_p = box(X0, Y0, X0 + 50, Y0 + 100).union(box(X0 + 100, Y0, X0 + 200, Y0 + 100))
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES, area_within=half_of_p)
    assert flows.polygons["area_ha"].to_dict() == approx({"P": 0.5, "Q": 1.0, "R": 0.0})
    assert flows.nodes["COMMERCIAL"].sum() == pytest.approx(0.05)
    assert flows.no_pipes == []  # R has no area left, so no flow to lose


def test_unrated_land_use_takes_nothing(net):
    spec = {**ZONES, "Q": (100, 200, "RESIDENTIAL")}
    flows = allocate_land_use(net, zoning(spec), land_use="use", rates=RATES)
    assert flows.pipes.loc[3, "polygon"] == "Q"
    assert flows.pipes.loc[3, "flow"] == 0.0
    assert "Q" not in flows.polygons.index
    assert set(flows.nodes.index) == {"N1", "N2"}


def test_units_and_placement(net):
    flows = allocate_land_use(
        net,
        zoning(),
        land_use="use",
        rates={"COMMERCIAL": 8.64, "INDUSTRIAL": 17.28},
        unit="m3/ha/d",
        flow_units="CMS",
        placement="downstream",
    )
    assert flows.nodes.loc["N2", "COMMERCIAL"] == pytest.approx(0.0625e-3)
    assert flows.nodes.loc["N5", "INDUSTRIAL"] == pytest.approx(0.08e-3)
    split = allocate_land_use(net, zoning(), land_use="use", rates=RATES, placement="split")
    assert split.nodes.loc["N2", "COMMERCIAL"] == pytest.approx(0.0625 / 2 + 0.0375 / 2)


def test_weight_column(net):
    # PDC_Count: pipe 1 has 2 and pipe 2 has 6; Q's pipe 3 is missing (0) and pipe 4 has 2.
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES, weight="PDC_Count")
    assert flows.nodes.loc["N1", "COMMERCIAL"] == pytest.approx(0.025)
    assert flows.nodes.loc["N4", "INDUSTRIAL"] == pytest.approx(0.2)
    assert "N3" not in flows.nodes.index


def test_force_mains_take_nothing():
    net = make_network([(1, "N1", "N2", 10, 60, 2), (2, "N2", "N3", 60, 90, 6, 50, "FM")])
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES)
    assert flows.nodes.loc["N1", "COMMERCIAL"] == pytest.approx(0.1)
    assert "N2" not in flows.nodes.index
    both = allocate_land_use(net, zoning(), land_use="use", rates=RATES, serve_force_mains=True)
    assert both.nodes.loc["N2", "COMMERCIAL"] == pytest.approx(0.0375)


def test_within_keeps_the_subnetwork_pipes(net):
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES)
    sub = net.catchment("N3")  # pipes 1 and 2
    assert flows.within(sub).to_dict("index") == {
        "N1": approx({"COMMERCIAL": 0.0625, "INDUSTRIAL": 0.0}),
        "N2": approx({"COMMERCIAL": 0.0375, "INDUSTRIAL": 0.0}),
    }


@pytest.mark.parametrize(
    "options, message",
    [
        (dict(land_use="nope"), "no column 'nope'"),
        (dict(rates={}), "non-empty mapping"),
        (dict(rates={"COMMERCIAL": -1}), "non-negative numbers"),
        (dict(total=-1.0), "total must be"),
        (dict(placement="middle"), "placement must be"),
        (dict(weight="nope"), "the pipes have no column 'nope'"),
        (dict(area_within="everywhere"), "area_within must be"),
    ],
)
def test_invalid_options(net, options, message):
    args = {"land_use": "use", "rates": RATES, **options}
    with pytest.raises(ValueError, match=message):
        allocate_land_use(net, zoning(), **args)


def test_unit_and_flow_units_refused(net):
    with pytest.raises(ValueError, match="unit must be"):
        allocate_land_use(net, zoning(), land_use="use", rates=RATES, unit="L/d")
    with pytest.raises(ValueError, match="flow_units"):
        allocate_land_use(net, zoning(), land_use="use", rates=RATES, flow_units="CFS")


def test_total_with_nowhere_to_go(net):
    spec = {"R": ZONES["R"]}
    with pytest.raises(ValueError, match="no polygon with a rate has a pipe"):
        allocate_land_use(net, zoning(spec), land_use="use", rates=RATES, total=1.0)


def test_geographic_crs_refused():
    net = make_network()
    net_geo = gpd.GeoDataFrame(net.full_data).to_crs("EPSG:4326")
    from swmmbuilder.land_use import _check

    with pytest.raises(ValueError, match="projected in metres"):
        _check(net_geo, zoning(), "use", RATES, None, None, None, "upstream", None)


# ---- pattern_inflows ---------------------------------------------------------------------------

FLAT = [1.0] * 24
DAY = [0.5] * 12 + [1.5] * 12


def test_pattern_inflows_one_series_per_kind_and_mixed_nodes():
    flows = pd.DataFrame(
        {"COM": {"A": 2.0, "B": 0.0, "C": 1.0}, "IND": {"A": 0.0, "B": 3.0, "C": 3.0}}
    )
    inflows, series = pattern_inflows(
        flows, {"COM": FLAT, "IND": DAY}, days=2, baseline={"A": 0.5, "D": 0.25}, prefix="NR_"
    )
    assert list(inflows) == ["A", "B", "C", "D"]
    assert inflows["A"].timeseries == "NR_COM" and inflows["A"].scale == 2.0
    assert inflows["A"].baseline == 0.5
    assert inflows["B"].timeseries == "NR_IND" and inflows["B"].baseline == 0.0
    assert inflows["C"].timeseries == "NR_C" and inflows["C"].scale == 4.0
    assert inflows["D"].timeseries is None and inflows["D"].baseline == 0.25
    assert list(series) == ["NR_COM", "NR_IND", "NR_C"]
    assert series["NR_IND"] == hourly_series(DAY, 2)
    mixed = (1.0 * np.array(FLAT) + 3.0 * np.array(DAY)) / 4.0
    assert np.array(series["NR_C"]) == pytest.approx(np.array(hourly_series(list(mixed), 2)))


def test_pattern_inflows_from_a_dict_and_a_land_use_allocation(net):
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES)
    inflows, _ = pattern_inflows(flows.nodes, {"COMMERCIAL": FLAT, "INDUSTRIAL": DAY}, days=1)
    assert {n: i.scale for n, i in inflows.items()} == approx(
        {"N1": 0.0625, "N2": 0.0375, "N3": 0.12, "N4": 0.08}
    )
    inflows, _ = pattern_inflows({"COM": {"A": 1.0}}, {"COM": FLAT}, days=1)
    assert inflows["A"].timeseries == "COM"


@pytest.mark.parametrize(
    "args, message",
    [
        (({"COM": {"A": 1.0}}, {}), "no pattern for COM"),
        (({"COM": {"A": -1.0}}, {"COM": FLAT}), "non-negative and finite"),
        (({"COM": {"A": 1.0}}, {"COM": FLAT}, {"A": float("nan")}), "baselines"),
        (({"COM": {"A": 1.0}}, {"COM": FLAT[:23]}), "24 finite"),
    ],
)
def test_pattern_inflows_refuses(args, message):
    flows, patterns, *baseline = args
    with pytest.raises(ValueError, match=message):
        pattern_inflows(flows, patterns, days=1, baseline=baseline[0] if baseline else None)


def test_pattern_inflows_name_clash():
    flows = pd.DataFrame({"COM": {"IND": 1.0}, "IND": {"IND": 1.0}})
    with pytest.raises(ValueError, match="both a pattern and a node"):
        pattern_inflows(flows, {"COM": FLAT, "IND": FLAT}, days=1)


def test_pattern_inflows_written_by_to_inp(net, tmp_path):
    """The two options go straight into the writer."""
    from conftest import TEMPLATE

    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES)
    inflows, series = pattern_inflows(flows.nodes, {"COMMERCIAL": FLAT, "INDUSTRIAL": DAY}, 3)
    sub = net.catchment("N6")
    model = sub.to_inp(
        tmp_path / "m.inp", ["N6"], template_path=TEMPLATE, inflows=inflows, timeseries=series
    )
    assert {i.node for i in model.inflows} == {"N1", "N2", "N3", "N4"}
    assert {s.name for s in model.timeseries} == {"COMMERCIAL", "INDUSTRIAL"}


def test_to_nearest_gives_a_pipeless_polygon_its_nearest_pipe(net):
    # R (x 300-400) has no pipe; the nearest is pipe 5 (x 190-250, 50 m away), upstream node N5.
    flows = allocate_land_use(net, zoning(), land_use="use", rates=RATES, to_nearest=True)
    assert flows.no_pipes == []
    assert flows.nearest.to_dict() == approx({"R": 50.0})
    assert flows.nodes.loc["N5", "COMMERCIAL"] == pytest.approx(0.1)
    assert flows.total == pytest.approx(0.4)
    assert flows.polygons.loc["R", "pipes"] == 1
    scaled = allocate_land_use(
        net, zoning(), land_use="use", rates=RATES, to_nearest=True, total=0.8
    )
    assert scaled.scale == pytest.approx(2.0)


def test_to_nearest_within_a_distance(net):
    flows = allocate_land_use(
        net, zoning(), land_use="use", rates=RATES, to_nearest=True, max_nearest_distance=40
    )
    assert flows.no_pipes == ["R"] and flows.nearest.empty
    with pytest.raises(ValueError, match="max_nearest_distance"):
        allocate_land_use(net, zoning(), land_use="use", rates=RATES, max_nearest_distance=-1)


def test_to_nearest_skips_polygons_with_no_flow(net):
    clip = box(X0, Y0, X0 + 200, Y0 + 100)  # R has no area left
    flows = allocate_land_use(
        net, zoning(), land_use="use", rates=RATES, to_nearest=True, area_within=clip
    )
    assert flows.nearest.empty and "N5" not in flows.nodes.index
