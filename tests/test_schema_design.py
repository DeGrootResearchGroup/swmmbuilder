"""design='schema': the schema's design rules as the defaults of build_model (#136)."""

from dataclasses import replace

import pytest
from test_pump_stations import make_network, template

from swmmbuilder import LONDON_ON, DataWarning, PumpStation
from swmmbuilder.swmm import build_model

STATION = ["U1", "PS", "F1", "MH", "OUT"]


@pytest.fixture(scope="module")
def sub():
    return make_network().subnetwork(STATION)


def test_schema_design_matches_explicit_rules(sub, tmp_path):
    """The build of the paper: London's rules passed one by one, and an ideal pump per station."""
    tpl = template(tmp_path)
    explicit = tmp_path / "explicit.inp"
    sub.to_inp(
        explicit,
        template_path=tpl,
        pump_stations={"PS": PumpStation(pump=dict(curve="*"))},
        force_main_discharge=LONDON_ON.design["force_main_discharge"],
        force_main_roughness=LONDON_ON.design["force_main_roughness"],
    )
    by_schema = tmp_path / "schema.inp"
    sub.to_inp(by_schema, template_path=tpl, pump_stations={"PS": PumpStation()}, design="schema")
    assert by_schema.read_text() == explicit.read_text()


def test_schema_rules_applied(sub, caplog):
    caplog.set_level("INFO", logger="swmmbuilder.swmm")
    model = build_model(sub, pump_stations={"PS": {}}, design="schema")
    assert model.pumps[0].curve == "*"
    # London's Hazen-Williams C for a 150 mm force main.
    assert {x.geom2 for x in model.xsections if x.shape == "FORCE_MAIN"} == {100.0}
    assert (
        "Design rules from schema London, ON: force_main_discharge, force_main_roughness, "
        "force_main_checks, default_pump" in caplog.text
    )


def test_explicit_option_wins(sub, caplog):
    caplog.set_level("INFO", logger="swmmbuilder.swmm")
    model = build_model(sub, pump_stations={"PS": {}}, force_main_roughness=130.0, design="schema")
    assert {x.geom2 for x in model.xsections if x.shape == "FORCE_MAIN"} == {130.0}
    assert "(given explicitly: force_main_roughness)" in caplog.text


def test_list_applies_only_those(sub):
    with pytest.raises(ValueError, match="PS has no pump"):
        build_model(sub, pump_stations={"PS": {}}, design=["force_main_roughness"])
    model = build_model(
        sub, pump_stations={"PS": {}}, design=["force_main_roughness", "default_pump"]
    )
    assert model.pumps[0].curve == "*"


def test_default_applies_nothing(sub):
    with pytest.raises(ValueError, match="PS has no pump"):
        build_model(sub, pump_stations={"PS": {}}, force_main_roughness=120.0)


def test_wet_well_defaults_only_by_name(sub):
    """'schema' leaves out wet_well_defaults, which needs a design_flow; naming it applies it."""
    build_model(sub, pump_stations={"PS": {}}, design="schema")
    with pytest.raises(ValueError, match="wet_well_defaults cannot be applied"):
        build_model(
            sub,
            pump_stations={"PS": {}},
            design=["default_pump", "wet_well_defaults"],
            force_main_roughness=120.0,
        )


@pytest.mark.parametrize(
    "design, match",
    [
        ("london", "design must be 'schema'"),
        (["force_main_rougness"], "unknown design rules force_main_rougness"),
    ],
)
def test_invalid_design_rejected(sub, design, match):
    with pytest.raises(ValueError, match=match):
        build_model(sub, design=design)


def test_rule_missing_from_schema_rejected():
    net = make_network()
    net.schema = replace(LONDON_ON, name="bare", design={})
    sub = net.subnetwork(STATION)
    with pytest.raises(ValueError, match="schema 'bare' has no default_pump rule"):
        build_model(sub, design=["default_pump"])
    # 'schema' with no design table applies nothing (the force main is then a conduit).
    with pytest.warns(DataWarning, match="adverse slope"):
        build_model(sub, design="schema")
