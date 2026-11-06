"""max_depth='median': junctions with no usable rim take the median recorded depth (#102)."""

import pytest
from conftest import DENDRITIC, subnetwork

from swmmbuilder import DataWarning, swmm


@pytest.fixture
def sub(network):
    """The dendritic fixture with rims 2, 3, and 4 m above the inverts of A1, A2, and B1, a rim
    below C1's invert, and none elsewhere."""
    sub = subnetwork(network, DENDRITIC)
    inverts = swmm.node_inverts(sub)
    for node, height in (("A1", 2.0), ("A2", 3.0), ("B1", 4.0), ("C1", -0.5)):
        sub.graph.nodes[node]["rim_elevation"] = inverts[node] + height
    return sub


def depths(model):
    return {j.name: j.max_depth for j in model.junctions}


def test_median_of_the_usable_depths(sub):
    with pytest.warns(DataWarning, match="C1: rim elevation"):
        model = swmm.build_model(sub, ["OUT"], max_depth="median")
    d = depths(model)
    assert (d["A1"], d["A2"], d["B1"]) == pytest.approx((2.0, 3.0, 4.0))
    assert {n: d[n] for n in ("B2", "C1", "J1", "J2")} == dict.fromkeys(
        ("B2", "C1", "J1", "J2"), 3.0
    )
    report = model.derived_report()
    rows = report[report["element"] == "junction"]
    assert sorted(rows["id"]) == ["B2", "C1", "J1", "J2"]
    assert set(rows["value"]) == {3.0}
    assert set(rows["rule"]) == {"max_depth='median'"}
    assert rows["basis"].iloc[0] == "median rim - invert of 3 junctions"


def test_a_function_of_the_depths(sub):
    with pytest.warns(DataWarning):
        model = swmm.build_model(sub, ["OUT"], max_depth=lambda d: d.min())
    assert depths(model)["J1"] == pytest.approx(2.0)
    assert set(model.derived_report()["rule"]) == {"max_depth function"}


def test_a_number_is_unchanged(sub):
    with pytest.warns(DataWarning):
        model = swmm.build_model(sub, ["OUT"], max_depth=1.5)
    assert depths(model)["J1"] == 1.5 and depths(model)["A1"] == pytest.approx(2.0)
    assert model.derived == []


def test_outfalls_are_not_counted(sub):
    inverts = swmm.node_inverts(sub)
    sub.graph.nodes["OUT"]["rim_elevation"] = inverts["OUT"] + 100.0
    with pytest.warns(DataWarning):
        model = swmm.build_model(sub, ["OUT"], max_depth="median")
    assert depths(model)["J1"] == pytest.approx(3.0)


def test_no_rims_to_take_a_median_of(network):
    sub = subnetwork(network, DENDRITIC)
    with pytest.raises(ValueError, match="no junction has a rim"):
        swmm.build_model(sub, ["OUT"], max_depth="median")


@pytest.mark.parametrize(
    "rule, message",
    [
        ("mean", "max_depth must be a number, 'median', or a function"),
        (lambda d: -1.0, "gave -1.0"),
    ],
)
def test_invalid(sub, rule, message):
    with pytest.raises(ValueError, match=message):  # before any junction is written
        swmm.build_model(sub, ["OUT"], max_depth=rule)
