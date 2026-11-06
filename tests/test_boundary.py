"""Pipes crossing a subnetwork's boundary, and networks ending near a catchment's outlet (#152)."""

import pytest
from test_population import make_network

from swmmbuilder import DataWarning

# OUT is the outlet. J has a second outlet, pipe 3, leaving the catchment (a split) for D2, far away. S is a summit:
# nothing drains to it, and its two runs fall to A1 (in the catchment) and to D4 (outside). P sends
# a force main, pipe 7, out of the catchment. Q1 -> E is another network that ends at E, about
# 41 m from OUT, without draining to it.
PIPES = [
    (1, "A1", "J", 0, 50, 0),
    (2, "J", "OUT", 50, 100, 0),
    (3, "J", "D2", 50, 250, 0, 200),
    (4, "S", "A1", -30, 0, 0),
    (5, "S", "D4", -30, -10, 0, 70),
    (6, "P", "J", 20, 50, 0, 30),
    (7, "P", "X", 20, 60, 0, 90, "FM"),
    (8, "Q1", "E", 120, 140, 0, 60),
]


@pytest.fixture
def net():
    return make_network(PIPES)


def catchment(net):
    with pytest.warns(DataWarning):
        return net.catchment("OUT")


def test_leaving_pipes_are_classified(net):
    sub = catchment(net)
    report = sub.boundary_pipes()
    rows = report.set_index("pipe_id")
    assert list(report["pipe_id"]) == [7, 3, 5]  # force main, split, summit
    assert rows.loc[3, ["direction", "node", "other", "kind"]].tolist() == [
        "out",
        "J",
        "D2",
        "split",
    ]
    assert rows.loc[3, "upstream_nodes"] == 3  # A1, S, and P
    assert rows.loc[5, "kind"] == "summit" and rows.loc[5, "upstream_nodes"] == 0
    assert rows.loc[7, "kind"] == "force_main" and rows.loc[7, "structure_type"] == "force_main"


def test_where_the_leaving_flow_ends(net):
    rows = catchment(net).boundary_pipes().set_index("pipe_id")
    assert rows.loc[3, "ends"] == "D2"
    assert rows.loc[3, "beyond_nodes"] == 1
    assert rows.loc[3, "end_distance"] == pytest.approx(((250 - 100) ** 2 + (200 - 50) ** 2) ** 0.5)


def test_upstream_values(net):
    sub = catchment(net)
    rows = sub.boundary_pipes(values={"A1": 10, "S": 5, "J": 1, "OUT": 100}, ends=False)
    rows = rows.set_index("pipe_id")
    assert rows.loc[3, "upstream_value"] == 16  # J, A1, and S
    assert rows.loc[5, "upstream_value"] == 5


def test_catchment_warns_about_splits_force_mains_and_nearby_ends(net):
    with pytest.warns(DataWarning) as caught:
        net.catchment("OUT")
    messages = [str(w.message) for w in caught]
    leaving = next(m for m in messages if "leave the catchment" in m)
    assert leaving.startswith("2 pipes leave the catchment of OUT from nodes with area upstream")
    assert "(1 force mains)" in leaving and "7 (force_main" in leaving and "3 (split" in leaving
    assert "5 (" not in leaving  # the summit is not warned about
    near = next(m for m in messages if "ends a network" in m)
    assert near.startswith("E ends a network of 1 nodes 41 m from OUT without draining to it")


def test_no_nearby_check_with_near_none(net):
    with pytest.warns(DataWarning) as caught:
        net.catchment("OUT", near=None)
    assert not any("ends a network" in str(w.message) for w in caught)


def test_validate_reports_the_boundary(net):
    report = catchment(net).validate()
    assert report.summary()["boundary_outflows"] == 3
    assert report.summary()["boundary_splits"] == 2
    assert report.summary()["boundary_inflows"] == 0
    assert list(report.boundary["pipe_id"]) == [7, 3, 5]


def test_inflows_of_a_subnetwork(net):
    sub = net.subnetwork(["J", "OUT"])
    rows = sub.boundary_pipes().set_index("pipe_id")
    assert set(rows.index[rows["direction"] == "in"]) == {1, 6}
    assert set(rows.loc[rows["direction"] == "in", "kind"]) == {"inflow"}


def test_network_ends(net):
    ends = net.ends()
    assert set(ends.index) == {"OUT", "D2", "D4", "X", "E"}
    assert ends["upstream_nodes"].to_dict() == {"OUT": 4, "D2": 4, "D4": 1, "X": 1, "E": 1}
    assert list(ends["upstream_nodes"]) == sorted(ends["upstream_nodes"], reverse=True)


def test_a_complete_catchment_says_nothing(recwarn):
    net = make_network([(1, "A1", "J", 0, 50, 0), (2, "J", "OUT", 50, 100, 0)])
    sub = net.catchment("OUT")
    assert sub.boundary_pipes().empty
    assert not [w for w in recwarn if issubclass(w.category, DataWarning)]
