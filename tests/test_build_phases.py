"""The phases of build_model (#77), tested directly: resolving the options, the nodes, and the
links."""

import pytest
from conftest import DENDRITIC, subnetwork
from test_pump_stations import MAIN, make_network

from swmmbuilder import PumpStation, WetWell
from swmmbuilder.gaps import fill_gaps, xsection
from swmmbuilder.swmm import SwmmModel
from swmmbuilder.swmm.build import (
    ModelInputs,
    _add_links,
    _add_nodes,
    _apply_design_rules,
    _Defaults,
    _resolve_inputs,
)

DEFAULTS = _Defaults(
    max_depth=0.0,
    init_depth=0.0,
    surcharge_depth=0.0,
    ponding_area=0.0,
    force_main_surcharge_depth=100.0,
    manning_n=0.013,
    manning_n_by_material={},
    init_flow=0.0,
    max_flow=0.0,
    barrels=1.0,
)
STATION = PumpStation(
    pump=dict(curve=[(0.0, 10.0)], curve_type="PUMP2", startup=1.0, shutoff=0.2),
    wet_well=WetWell(max_depth=4.0, area=2.0),
)


def resolve(sub, **options):
    args = dict(
        outfalls=None,
        height_correction=None,
        pumps=None,
        wet_wells=None,
        pump_stations=None,
        default_pump=None,
        dwf=None,
        patterns=None,
        dwf_patterns=None,
        timeseries=None,
        inflows=None,
        flow_units="LPS",
    )
    return _resolve_inputs(sub, **{**args, **options})


@pytest.fixture(scope="module")
def station_sub():
    return make_network().subnetwork(MAIN)


def test_resolve_inputs(station_sub):
    inputs = resolve(station_sub, pump_stations={"PS": STATION}, dwf={"U1": 2.0}, flow_units="lps")
    assert isinstance(inputs, ModelInputs)
    assert inputs.flow_units == "LPS"
    assert inputs.wet_wells == {"PS": STATION.wet_well}  # the station's well joins wet_wells
    assert inputs.chains == {2: [2, 3]} and inputs.force_mains == {2, 3}
    assert [link.name for link in inputs.station_pumps] == ["PS_P"]
    assert inputs.dwf == {"U1": 2.0}
    assert set(inputs.outfalls) == {"OUT", "OVO"}  # the outfall nodes, by default


@pytest.mark.parametrize(
    "options, message",
    [
        (dict(outfalls=["NOPE"]), "Outfalls not in the subnetwork: NOPE"),
        (dict(flow_units="CFS"), "flow_units must be one of"),
        (
            dict(pump_stations={"PS": STATION}, wet_wells={"PS": dict(max_depth=1, area=1)}),
            "Wet wells given both in wet_wells and pump_stations: PS",
        ),
        (
            dict(pump_stations={"PS": STATION}, pumps={2: dict(curve="*")}),
            "Pipes both in pumps and in a pump station's force main: 2",
        ),
    ],
)
def test_resolve_inputs_refuses_contradictions(station_sub, options, message):
    with pytest.raises(ValueError, match=message):
        resolve(station_sub, **options)


def _built(sub, **options):
    inputs = resolve(sub, **options)
    filled = fill_gaps(sub, pumps=inputs.pumps, height_correction=inputs.height_correction)
    pipes, inverts = filled.pipes, filled.inverts
    _apply_design_rules(sub, inputs, pipes, inverts, None)
    return inputs, pipes, inverts


def test_add_nodes_adds_station_discharge_nodes(station_sub):
    inputs, pipes, inverts = _built(station_sub, pump_stations={"PS": STATION})
    model = SwmmModel()
    discharge_of = _add_nodes(model, station_sub, inputs, pipes, inverts, DEFAULTS)
    assert discharge_of == {2: "PS_D"}
    assert inverts["PS_D"] == pytest.approx(249.8)  # added for the links
    junctions = {j.name: j for j in model.junctions}
    assert junctions["PS_D"].surcharge_depth == 100.0
    assert junctions["F1"].surcharge_depth == 100.0  # inside the force main
    assert [s.name for s in model.storage] == ["PS"]
    assert {o.name for o in model.outfalls} == {"OUT", "OVO"}
    assert model.coordinates["PS_D"] == model.coordinates["PS"]


def test_add_links(network):
    sub = subnetwork(network, DENDRITIC)
    inputs, pipes, inverts = _built(sub, outfalls=["OUT"])
    sections = {
        name: xsection(row)
        for name, row in zip(pipes["pipe_id"], pipes.to_dict("records"), strict=True)
    }
    model = SwmmModel()
    adverse = _add_links(model, inputs, pipes, inverts, {}, sections, {}, DEFAULTS)
    assert adverse == []
    assert sorted(c.name for c in model.conduits) == sorted(pipes["pipe_id"])
    assert {x.name for x in model.xsections} == set(pipes["pipe_id"])
    assert all(c.roughness == 0.013 for c in model.conduits)
    assert not model.pumps and not model.curves
