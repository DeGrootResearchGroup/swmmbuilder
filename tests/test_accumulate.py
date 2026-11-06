"""Subnetwork.accumulate: totals of a node value upstream of each node and pipe (#107)."""

import numpy as np
import pandas as pd
import pytest
from conftest import DENDRITIC, subnetwork
from test_population import make_network

from swmmbuilder import UpstreamTotals
from swmmbuilder.subnetwork import _trapped

ONES = dict.fromkeys(["A1", "A2", "B1", "B2", "C1", "J1", "J2", "OUT"], 1.0)


@pytest.fixture
def sub(network):
    """Branches A1-A2, B1-B2, and C1 meet at J1; two parallel pipes J1 -> J2; then OUT."""
    return subnetwork(network, DENDRITIC)


def parallel(sub):
    return [k for u, v, k in sub.graph.edges(keys=True) if (u, v) == ("J1", "J2")]


def test_equal_split_conserves(sub):
    totals = sub.accumulate(ONES)
    assert isinstance(totals, UpstreamTotals) and totals.split == "equal"
    assert totals.nodes["J1"] == 6 and totals.nodes["J2"] == 7 and totals.nodes["OUT"] == 8
    a, b = parallel(sub)
    assert totals.pipes[a] == totals.pipes[b] == 3
    assert totals.fractions[a] == 0.5


def test_full_counts_each_upstream_node_once(sub):
    totals = sub.accumulate(ONES, split="full")
    assert totals.nodes["J2"] == 7 and totals.nodes["OUT"] == 8
    a, b = parallel(sub)
    assert totals.pipes[a] == totals.pipes[b] == 6


def test_given_fractions_and_the_rest_shared(sub):
    a, b = parallel(sub)
    totals = sub.accumulate(ONES, split={a: 0.25})
    assert totals.split == "given"
    assert totals.pipes[a] == pytest.approx(1.5) and totals.pipes[b] == pytest.approx(4.5)
    assert totals.nodes["J2"] == pytest.approx(7)
    with pytest.raises(ValueError, match="more than 1 at J1"):
        sub.accumulate(ONES, split={a: 0.7, b: 0.5})
    with pytest.raises(ValueError, match="finite and >= 0"):
        sub.accumulate(ONES, split={a: -0.1})


def test_missing_nodes_count_zero_and_a_dataframe(sub):
    assert sub.accumulate({"A1": 2.0}).nodes["OUT"] == 2.0
    frame = pd.DataFrame({"x": {"A1": 1.0, "C1": 2.0}, "y": {"B1": 5.0}})
    totals = sub.accumulate(frame)
    assert totals.nodes.loc["J1"].to_dict() == {"x": 3.0, "y": 5.0}
    assert totals.pipes.loc[parallel(sub)[0]].to_dict() == {"x": 1.5, "y": 2.5}


def test_unknown_split(sub):
    with pytest.raises(ValueError, match="split must be one of"):
        sub.accumulate(ONES, split="flow")


# A loop with an exit: N1 -> N2 -> N3, N3 -> N2 back, and N3 -> N4.
LOOP = [
    (1, "N1", "N2", 0, 50, 1),
    (2, "N2", "N3", 50, 100, 1),
    (3, "N3", "N2", 100, 50, 1, 60),
    (4, "N3", "N4", 100, 150, 1),
]


def test_a_loop_with_an_exit():
    sub = make_network(LOOP).catchment("N4")
    values = dict.fromkeys(["N1", "N2", "N3", "N4"], 1.0)
    # N2 = 1 + N1 + N3 / 2 and N3 = 1 + N2, so N2 = 5, N3 = 6, and N4 = 1 + N3 / 2 = 4.
    equal = sub.accumulate(values)
    assert equal.nodes.to_dict() == pytest.approx({"N1": 1, "N2": 5, "N3": 6, "N4": 4})
    full = sub.accumulate(values, split="full")
    assert full.nodes.to_dict() == {"N1": 1, "N2": 3, "N3": 3, "N4": 4}


def test_a_loop_with_no_exit():
    pipes = [
        (1, "N1", "N2", 0, 50, 1),
        (2, "N2", "N3", 50, 100, 1),
        (3, "N3", "N2", 100, 50, 1, 60),
    ]
    sub = make_network(pipes).subnetwork(["N1", "N2", "N3"])
    with pytest.raises(ValueError, match="has no exit"):
        sub.accumulate(dict.fromkeys(["N1", "N2", "N3"], 1.0))


def two_way_chain(n):
    """N0 ... N<n-1> joined by pipes both ways, with exits N0 -> A and N<n-1> -> B."""
    pipes = []
    for i in range(n - 1):
        pipes.append((i + 1, f"N{i}", f"N{i + 1}", 10 * i, 10 * i + 10, 1))
        pipes.append((1000 + i, f"N{i + 1}", f"N{i}", 10 * i + 10, 10 * i, 1, 60))
    pipes.append((5000, "N0", "A", -10, 0, 1))
    pipes.append((5001, f"N{n - 1}", "B", 10 * n - 10, 10 * n, 1))
    return make_network(pipes).subnetwork([*(f"N{i}" for i in range(n)), "A", "B"])


def test_a_long_loop_with_exits_is_solved():
    # A determinant test refused this: det(I - A) = 61 / 2**60, about 5e-17 (#119).
    n = 60
    sub = two_way_chain(n)
    totals = sub.accumulate({f"N{i}": 1.0 for i in range(n)})
    assert totals.nodes["A"] == pytest.approx(n / 2)
    assert totals.nodes["B"] == pytest.approx(n / 2)
    assert totals.nodes["A"] + totals.nodes["B"] == pytest.approx(n)


def test_only_the_nodes_with_no_exit_are_named():
    # Shares within a loop: 0 -> 1 -> 2 -> 1 (nodes 1 and 2 send everything to each other), and
    # node 0 sends half its total out.
    a = np.zeros((3, 3))
    a[1, 0] = 0.5
    a[2, 1] = 1.0
    a[1, 2] = 1.0
    assert _trapped(["x", "y", "z"], a) == ["y", "z"]
    a[1, 2] = 0.9  # z now sends a tenth out: every node reaches an exit
    assert _trapped(["x", "y", "z"], a) == []
