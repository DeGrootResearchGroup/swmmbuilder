"""Outfalls with more than one link under DYNWAVE: refused, or written through a dummy conduit
(the outfall_link option, #101)."""

import pytest
from conftest import read_sections, run_swmm, subnetwork, write_inp

from swmmbuilder import swmm

# Three branches meeting at J1, which is the outfall: three inlet links.
TO_J1 = [("A1", "A2"), ("A2", "J1"), ("B1", "B2"), ("B2", "J1"), ("C1", "J1")]
DWF = {"A1": 1.0, "B1": 0.5, "C1": 0.75}


@pytest.fixture
def sub(network):
    return subnetwork(network, TO_J1)


def test_dynwave_refuses_a_crowded_outfall(sub, tmp_path):
    with pytest.raises(ValueError, match=r"ERROR 141.*outfall_link=True.*J1 \("):
        write_inp(sub, tmp_path, routing="DYNWAVE", outfalls=["J1"], dwf=DWF)
    assert not (tmp_path / "model.inp").exists()


def test_kinwave_needs_no_link(sub, tmp_path):
    path = write_inp(sub, tmp_path, outfalls=["J1"], dwf=DWF)  # the shipped template: KINWAVE
    assert [r[0] for r in read_sections(path)["OUTFALLS"]] == ["J1"]


def test_outfall_link_records(sub):
    model = swmm.build_model(sub, ["J1"], outfall_link=True, max_depth=2.5)
    assert swmm.outfall_inlets(model) == {}
    assert [o.name for o in model.outfalls] == ["J1_OUT"]
    j1 = next(j for j in model.junctions if j.name == "J1")
    out = model.outfalls[0]
    assert out.type == "FREE"
    assert out.invert == pytest.approx(j1.invert - 0.3)
    dummy = next(c for c in model.conduits if c.name == "J1_DUMMY")
    assert (dummy.inlet, dummy.outlet, dummy.length, dummy.roughness) == (
        "J1",
        "J1_OUT",
        30.0,
        0.01,
    )
    xs = next(x for x in model.xsections if x.name == "J1_DUMMY")
    assert (xs.shape, xs.geom1) == ("CIRCULAR", 3.0)
    x, y = model.coordinates["J1"]
    assert model.coordinates["J1_OUT"] == (x + 1.0, y)


def test_outfall_link_options(sub):
    model = swmm.build_model(sub, {"J1": "FIXED 9.0"}, outfall_link={"length": 10, "slope": 0.02})
    out = model.outfalls[0]
    assert out.type == "FIXED 9.0"
    j1 = next(j for j in model.junctions if j.name == "J1")
    assert out.invert == pytest.approx(j1.invert - 0.2)
    with pytest.raises(TypeError, match="lenght"):
        swmm.build_model(sub, ["J1"], outfall_link={"lenght": 10})
    with pytest.raises(ValueError, match="positive numbers"):
        swmm.build_model(sub, ["J1"], outfall_link={"diameter": 0})


def test_single_inlet_outfalls_are_left_alone(network):
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1")])
    model = swmm.build_model(sub, ["J1"], outfall_link=True)
    assert [o.name for o in model.outfalls] == ["J1"]
    assert not any(str(c.name).endswith("_DUMMY") for c in model.conduits)


def test_swmm_runs_with_the_link(sub, tmp_path):
    pytest.importorskip("pyswmm")
    path = write_inp(sub, tmp_path, routing="DYNWAVE", outfalls=["J1"], dwf=DWF, outfall_link=True)
    s = read_sections(path)
    assert [r[0] for r in s["OUTFALLS"]] == ["J1_OUT"]
    assert "J1" in {r[0] for r in s["JUNCTIONS"]}
    assert run_swmm(path, tmp_path) == []


def test_a_linked_outfall_given_the_rule_depth_is_derived(sub):
    """An outfall written as a junction with no usable rim takes the max_depth rule's depth and
    is listed in derived, like the other junctions (#123)."""
    model = swmm.build_model(sub, ["J1"], outfall_link=True, max_depth=lambda d: 2.0)
    assert next(j for j in model.junctions if j.name == "J1").max_depth == 2.0
    assert "J1" in set(model.derived_report()["id"])
