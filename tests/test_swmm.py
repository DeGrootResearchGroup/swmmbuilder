"""The SWMM model of a subnetwork (swmm.build_model), without writing a file."""

import pytest
from conftest import DENDRITIC, TEMPLATE, subnetwork_with_edges

from sewernetworkx import swmm


@pytest.fixture
def dendritic(network):
    return subnetwork_with_edges(network, DENDRITIC)


@pytest.fixture
def model(dendritic):
    return swmm.build_model(dendritic, manning_n=0.02, manning_n_by_material={"PVC": 0.011})


def test_junctions_and_outfalls(model):
    assert [o.name for o in model.outfalls] == ["OUT"]
    assert model.outfalls[0] == swmm.Outfall("OUT", 103.30, "FREE")
    assert {j.name for j in model.junctions} == {n for e in DENDRITIC for n in e} - {"OUT"}


def test_node_inverts_are_lowest_connected_invert(dendritic):
    inverts = swmm.node_inverts(dendritic)
    # J1: pipes in at 104.00, 103.95 and drop 103.98; the parallel J1 -> J2 pipes leave at 103.90.
    assert inverts["J1"] == pytest.approx(103.90)
    assert inverts["A1"] == pytest.approx(105.00)
    assert swmm.node_inverts(dendritic, {"J1": 103.0})["J1"] == 103.0


def test_conduit_offsets_reproduce_gis_inverts(model, dendritic):
    inverts = swmm.node_inverts(dendritic)
    pipes = dendritic.subnetwork_pipes().set_index("pipe_id")
    assert len(model.conduits) == len(pipes) == 8
    for c in model.conduits:
        assert inverts[c.inlet] + c.inlet_offset == pytest.approx(pipes.loc[c.name, "us_invert"])
        assert inverts[c.outlet] + c.outlet_offset == pytest.approx(pipes.loc[c.name, "ds_invert"])


def test_roughness_by_material(model, dendritic):
    material = dendritic.subnetwork_pipes().set_index("pipe_id")["material"]
    for c in model.conduits:
        assert c.roughness == (0.011 if material[c.name] == "PVC" else 0.02)


def test_xsections(model):
    sections = {x.name: (x.shape, x.geom1, x.geom2, x.barrels) for x in model.xsections}
    assert sections[1003] == ("RECT_CLOSED", pytest.approx(0.3), pytest.approx(0.45), 1.0)
    assert sections[1000] == ("CIRCULAR", pytest.approx(0.2), 0.0, 1.0)


def test_outfalls_override(dendritic):
    model = swmm.build_model(dendritic, outfalls={"OUT": "FIXED"})
    assert [(o.name, o.type) for o in model.outfalls] == [("OUT", "FIXED")]


def test_render_fills_every_placeholder(model):
    template = TEMPLATE.read_text()
    text = model.render(template)
    assert set(model.sections()) == set(swmm.TEMPLATE_PLACEHOLDERS)
    for name in swmm.TEMPLATE_PLACEHOLDERS:
        assert f"{{{name}}}" not in text
    assert text.startswith(template.split("{")[0])


def test_write_inp_writes_nothing_when_invalid(network, tmp_path):
    sub = subnetwork_with_edges(network, DENDRITIC + [("E1", "J2")])  # E1-J2 has no shape
    path = tmp_path / "bad.inp"
    with pytest.raises(ValueError, match="cross-section"):
        swmm.write_inp(sub, path, template_path=_dynwave(tmp_path))
    assert not path.exists()


def test_write_inp_returns_the_model(dendritic, tmp_path):
    path = tmp_path / "model.inp"
    model = swmm.write_inp(dendritic, path, template_path=_dynwave(tmp_path))
    assert path.read_text() == model.render(_dynwave(tmp_path).read_text())


def _dynwave(tmp_path):
    """The shipped template with DYNWAVE routing (the fixture's J1 -> J2 pipes are parallel)."""
    path = tmp_path / "dynwave.inp"
    path.write_text(
        TEMPLATE.read_text().replace("FLOW_ROUTING         KINWAVE", "FLOW_ROUTING DYNWAVE")
    )
    return path
