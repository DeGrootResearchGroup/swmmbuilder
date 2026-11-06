"""Spreading a flow over the pipes of a network (swmmbuilder.distribute, #100)."""

import pandas as pd
import pytest
from test_population import make_network

from swmmbuilder import DataWarning, distribute_flow
from swmmbuilder.distribute import pipe_weights

# Pipes along y = 50 (test_population.PIPES): 1 N1-N2 50 m, 2 N2-N3 30 m, 3 N3-N4 60 m,
# 4 N4-N5 40 m, 5 N5-N6 60 m; 240 m in all.


@pytest.fixture
def net():
    return make_network()


def approx(x):
    return pytest.approx(x, abs=1e-12)


def test_total_shared_by_length_at_upstream_nodes(net):
    flows = distribute_flow(net, total=24.0)
    assert flows == approx({"N1": 5.0, "N2": 3.0, "N3": 6.0, "N4": 4.0, "N5": 6.0})
    assert list(flows) == ["N1", "N2", "N3", "N4", "N5"]


def test_rate_per_metre(net):
    assert distribute_flow(net, rate=0.1)["N3"] == pytest.approx(6.0)


def test_placements(net):
    down = distribute_flow(net, total=24.0, placement="downstream")
    assert down == approx({"N2": 5.0, "N3": 3.0, "N4": 6.0, "N5": 4.0, "N6": 6.0})
    split = distribute_flow(net, total=24.0, placement="split")
    assert split["N1"] == pytest.approx(2.5) and split["N2"] == pytest.approx(4.0)
    assert sum(split.values()) == pytest.approx(24.0)


def test_subnetwork_only(net):
    sub = net.catchment("N3")  # pipes 1 and 2, 80 m
    assert distribute_flow(sub, total=8.0) == approx({"N1": 5.0, "N2": 3.0})


def test_force_mains_take_nothing():
    net = make_network([(1, "N1", "N2", 10, 60, 2), (2, "N2", "N3", 60, 90, 6, 50, "FM")])
    assert distribute_flow(net, total=1.0) == approx({"N1": 1.0})
    both = distribute_flow(net, total=8.0, serve_force_mains=True)
    assert both == approx({"N1": 5.0, "N2": 3.0})


def test_weights(net):
    pipes = net.full_data.set_index("pipe_id")
    pipes["diameter"] = [0.2, 0.2, 0.4, 0.4, None]
    w = pipe_weights(pipes, "diameter_length")
    assert w.iloc[:4].tolist() == approx([10.0, 6.0, 24.0, 16.0]) and pd.isna(w.iloc[4])
    with pytest.warns(DataWarning):  # 2, 6, missing, 2, 999
        by_column = distribute_flow(net, total=10.0, weight="PDC_Count")
    assert by_column["N2"] == pytest.approx(10.0 * 6 / 1009)
    by_callable = distribute_flow(
        net, total=2.0, weight=lambda p: (p["pipe_id"] <= 2).astype(float)
    )
    assert by_callable == approx({"N1": 1.0, "N2": 1.0})


def test_missing_weights_warn(net):
    with pytest.warns(DataWarning, match="1 pipes have no valid weight"):
        flows = distribute_flow(net, total=10.0, weight="PDC_Count")
    assert "N3" not in flows


@pytest.mark.parametrize(
    "options, message",
    [
        (dict(), "give one of total and rate"),
        (dict(total=1.0, rate=1.0), "give one of total and rate"),
        (dict(total=-1.0), "total must be"),
        (dict(rate=float("nan")), "rate must be"),
        (dict(total=1.0, placement="middle"), "placement must be"),
        (dict(total=1.0, weight="nope"), "no column 'nope'"),
        (dict(total=1.0, weight=lambda p: p["length"] * 0), "weigh nothing"),
    ],
)
def test_invalid(net, options, message):
    with pytest.raises(ValueError, match=message):
        distribute_flow(net, **options)


def test_feeds_pattern_inflows(net):
    from swmmbuilder import pattern_inflows

    infiltration = distribute_flow(net, total=24.0)
    inflows, series = pattern_inflows(pd.DataFrame(), {}, days=1, baseline=infiltration)
    assert {n: i.baseline for n, i in inflows.items()} == approx(infiltration)
    assert series == {}
