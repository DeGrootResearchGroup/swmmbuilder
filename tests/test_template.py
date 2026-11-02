"""The packaged SI template (#5): loading, sync with the example copy, and rejection of templates
the writer cannot fill."""

import pytest
from conftest import TEMPLATE, read_sections, subnetwork

from sewernetworkx.swmm import SI_FLOW_UNITS, check_template, default_template_path


@pytest.fixture
def sub(network):
    # Without the parallel J1 -> J2 pipes, which the template's KINWAVE routing rejects (#3).
    edges = [("A1", "A2"), ("A2", "J1"), ("B1", "B2"), ("B2", "J1"), ("C1", "J1")]
    return subnetwork(network, edges)


def test_packaged_template_is_si_and_fillable():
    content = default_template_path().read_text()
    check_template(content)
    opts = {row[0]: row[1] for row in read_sections(default_template_path())["OPTIONS"]}
    assert opts["FLOW_UNITS"] == "LPS"


def test_example_template_matches_packaged():
    assert TEMPLATE.read_text() == default_template_path().read_text()


def test_default_template_used_without_template_path(sub, tmp_path):
    sub.to_inp(tmp_path / "model.inp", ["J1"])
    text = (tmp_path / "model.inp").read_text()
    assert "{" not in text
    sections = read_sections(tmp_path / "model.inp")
    opts = {row[0]: row[1] for row in sections["OPTIONS"]}
    assert opts["FLOW_UNITS"] in SI_FLOW_UNITS
    assert {"JUNCTIONS", "OUTFALLS", "CONDUITS", "XSECTIONS", "COORDINATES"} <= set(sections)


def custom_template(tmp_path, old, new):
    path = tmp_path / "custom.inp"
    path.write_text(default_template_path().read_text().replace(old, new))
    return path


@pytest.mark.parametrize("units", ["MGD", "CFS", "GPM"])
def test_us_units_template_rejected(sub, tmp_path, units):
    path = custom_template(tmp_path, "FLOW_UNITS           LPS", f"FLOW_UNITS           {units}")
    with pytest.raises(ValueError, match=units):
        sub.to_inp(tmp_path / "m.inp", ["J1"], template_path=path)


def test_missing_flow_units_rejected(sub, tmp_path):
    path = custom_template(tmp_path, "FLOW_UNITS           LPS", "")
    with pytest.raises(ValueError, match="FLOW_UNITS"):
        sub.to_inp(tmp_path / "m.inp", ["J1"], template_path=path)


def test_pumps_placeholder_rejected(sub, tmp_path):
    path = custom_template(tmp_path, "[LOSSES]", "[PUMPS]\n{pumps}\n\n[LOSSES]")
    with pytest.raises(ValueError, match="pumps"):
        sub.to_inp(tmp_path / "m.inp", ["J1"], template_path=path)
    assert not (tmp_path / "m.inp").exists()


@pytest.mark.parametrize("units", SI_FLOW_UNITS)
def test_si_units_template_accepted(sub, tmp_path, units):
    path = custom_template(tmp_path, "FLOW_UNITS           LPS", f"FLOW_UNITS           {units}")
    sub.to_inp(tmp_path / "m.inp", ["J1"], template_path=path)
    assert (tmp_path / "m.inp").exists()
