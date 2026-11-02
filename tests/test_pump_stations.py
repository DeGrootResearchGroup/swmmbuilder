"""Pumping stations with force mains (#34): the pump at the station node, the force main kept as
full-flowing conduits."""

import math
import re

import geopandas as gpd
import pytest
from conftest import TEMPLATE, read_sections, run_swmm
from shapely.geometry import LineString

from swmmbuilder import LONDON_ON, Network, PumpStation, WetWell
from swmmbuilder.swmm import force_main_roughness_rule

# Layout (flow direction; => is a force main, ~> an overflow), in UTM 17N metres:
#
#     U1 -> PS => F1 => MH -> OUT     station PS, two-pipe force main rising 249.8 -> 252.0 m
#           PS ~> OVO                 overflow to an outfall
#     PS2 -> M5 -> OUT5               station with no force main
#     Q1 => Q2                        force main that starts at no station, to a dead end
#     PT => T1 -> TO1, PT => T2 -> TO2    twin force mains
#     PR => RJ, PR2 => RJ, RJ => RZ -> RO   two stations' force mains joining
POS = {
    "U1": (0, 0),
    "PS": (50, 0),
    "F1": (350, 0),
    "MH": (550, 0),
    "OUT": (600, 0),
    "OVO": (50, 50),
    "PS2": (0, 200),
    "M5": (50, 200),
    "OUT5": (100, 200),
    "Q1": (0, 300),
    "Q2": (50, 300),
    "PT": (0, 400),
    "T1": (100, 400),
    "T2": (100, 450),
    "TO1": (150, 400),
    "TO2": (150, 450),
    "PR": (0, 500),
    "PR2": (0, 550),
    "RJ": (100, 500),
    "RZ": (200, 500),
    "RO": (250, 500),
}
TYPES = {"PS": "PS", "PS2": "PS", "PT": "PS", "PR": "PS", "PR2": "PS"}
TYPES.update(dict.fromkeys(["OUT", "OVO", "OUT5", "TO1", "TO2", "RO"], "OF"))

# (id, up, down, inv_up, inv_down, length, diameter mm, structure type)
PIPES = [
    (1, "U1", "PS", 250.5, 250.0, 50.0, 300, "SL"),
    (2, "PS", "F1", 249.8, 251.0, 300.0, 150, "FM"),
    (3, "F1", "MH", 251.0, 252.0, 200.0, 150, "FM"),
    (4, "MH", "OUT", 251.9, 251.0, 50.0, 300, "SL"),
    (5, "PS", "OVO", 250.8, 250.0, 50.0, 300, "OVF"),
    (6, "PS2", "M5", 250.0, 249.5, 50.0, 300, "SL"),
    (7, "M5", "OUT5", 249.4, 249.0, 50.0, 300, "SL"),
    (8, "Q1", "Q2", 250.0, 251.0, 50.0, 150, "FM"),
    (9, "PT", "T1", 249.0, 250.0, 100.0, 200, "FM"),
    (10, "PT", "T2", 249.0, 250.0, 100.0, 650, "FM"),
    (11, "T1", "TO1", 249.9, 249.5, 50.0, 300, "SL"),
    (12, "T2", "TO2", 249.9, 249.5, 50.0, 300, "SL"),
    (13, "PR", "RJ", 249.0, 250.0, 100.0, 200, "FM"),
    (14, "PR2", "RJ", 249.0, 250.0, 100.0, 200, "FM"),
    (15, "RJ", "RZ", 250.0, 251.0, 100.0, 300, "FM"),
    (16, "RZ", "RO", 250.9, 250.0, 50.0, 300, "SL"),
]
FORCE_MAIN_VOLUME = math.pi / 4 * 0.15**2 * (300.0 + 200.0)  # m3
MAIN = ["U1", "PS", "F1", "MH", "OUT", "OVO"]
IDEAL = dict(pump=dict(curve="*"))
C = 120.0


def make_network(pipes=PIPES, types=TYPES):
    rows = [
        {
            "OBJECTID": pid,
            "FlowType": "SAN",
            "UpstreamID": up,
            "DownstreamID": down,
            "UpstreamInventoryType": types.get(up, "MH"),
            "DownstreamInventoryType": types.get(down, "MH"),
            "Length": length,
            "PipeShape": "R",
            "Diameter": diameter,
            "Height": diameter,
            "Width": diameter,
            "Material": "PVC" if diameter < 300 else "CONC",
            "UpstreamInvert": inv_up,
            "DownstreamInvert": inv_down,
            "DropStructureInvert": 0.0,
            "ConstructedStatus": "Built",
            "StructureType": structure,
        }
        for pid, up, down, inv_up, inv_down, length, diameter, structure in pipes
    ]
    lines = [
        LineString([(480000 + x, 4760000 + y) for x, y in (POS[up], POS[down])])
        for _, up, down, *_ in pipes
    ]
    return Network(gpd.GeoDataFrame(rows, geometry=lines, crs="EPSG:32617"))


@pytest.fixture(scope="module")
def net():
    return make_network()


@pytest.fixture(scope="module")
def sub(net):
    return net.subnetwork(MAIN)


def template(tmp_path, inflow=None, hours=None):
    """A DYNWAVE copy of the test template, optionally with a constant inflow (L/s) at PS and a
    shorter run."""
    text = TEMPLATE.read_text().replace("KINWAVE", "DYNWAVE")
    if inflow is not None:
        line = f'PS  FLOW  ""  FLOW  1.0  1.0  {inflow}\n'
        text = text.replace("[INFLOWS]\n", "[INFLOWS]\n" + line)
    if hours is not None:
        text = re.sub(r"(?m)^END_DATE\s+\S+", "END_DATE             10/01/2016", text)
        text = re.sub(r"(?m)^END_TIME\s+\S+", f"END_TIME             {hours:02d}:00:00", text)
        text = re.sub(r"(?m)^ROUTING_STEP\s+\S+", "ROUTING_STEP         5", text)
    path = tmp_path / "template.inp"
    path.write_text(text)
    return path


def write(sub, tmp_path, outfalls=None, tpl=None, **kwargs):
    kwargs.setdefault("force_main_roughness", C)
    path = tmp_path / "model.inp"
    sub.to_inp(path, outfalls, template_path=tpl or template(tmp_path), **kwargs)
    return path


# --- What is written ----------------------------------------------------------------------------


def test_station_with_wet_well_and_overflow(sub, tmp_path):
    station = PumpStation(
        pump=dict(curve=[(0.0, 10.0)], curve_type="PUMP2", startup=1.0, shutoff=0.2),
        wet_well=WetWell(max_depth=4.0, area=2.0),
    )
    s = read_sections(write(sub, tmp_path, pump_stations={"PS": station}))

    assert s["PUMPS"] == [["PS_P", "PS", "PS_D", "PS_curve", "ON", "1.0", "0.2"]]
    assert s["CURVES"] == [["PS_curve", "PUMP2", "0.0", "10.0"]]
    (storage,) = s["STORAGE"]
    assert storage[0] == "PS" and float(storage[2]) == 4.0
    junctions = {r[0]: r for r in s["JUNCTIONS"]}
    assert float(junctions["PS_D"][1]) == pytest.approx(249.8)
    assert float(junctions["PS_D"][4]) == 100.0  # surcharge depth
    assert float(junctions["F1"][4]) == 100.0  # interior node of the force main
    assert float(junctions["MH"][4]) == 0.0  # where it ends, at a gravity maintenance hole

    conduits = {r[0]: r for r in s["CONDUITS"]}
    assert conduits["2"][1:4] == ["PS_D", "F1", "300.0"]
    assert conduits["3"][1:4] == ["F1", "MH", "200.0"]
    assert conduits["5"][1:3] == ["PS", "OVO"]  # the overflow stays a conduit
    xsections = {r[0]: r for r in s["XSECTIONS"]}
    for pipe in ("2", "3"):
        assert xsections[pipe][1:4] == ["FORCE_MAIN", "0.15", "120.0"]
    assert xsections["5"][1] == "CIRCULAR"
    coordinates = {r[0]: r[1:] for r in s["COORDINATES"]}
    assert coordinates["PS_D"] == coordinates["PS"]


def test_force_main_is_not_reported_adverse(sub, tmp_path, recwarn):
    write(
        sub, tmp_path, pump_stations={"PS": dict(pump=dict(curve=[(0, 5.0)], curve_type="PUMP2"))}
    )
    assert not [w for w in recwarn if "adverse" in str(w.message)]


def test_station_without_wet_well_is_a_junction(net, tmp_path):
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    s = read_sections(write(sub, tmp_path, pump_stations={"PS": IDEAL}))
    assert "PS" in {r[0] for r in s["JUNCTIONS"]}
    assert s["PUMPS"] == [["PS_P", "PS", "PS_D", "*", "ON", "0.0", "0.0"]]
    assert s.get("STORAGE", []) == []


def test_twin_force_mains(net, tmp_path):
    sub = net.subnetwork(["PT", "T1", "T2", "TO1", "TO2"])
    shared = read_sections(
        write(sub, tmp_path, pump_stations={"PT": dict(IDEAL, discharge="shared")})
    )
    curve = dict(curve=[(0, 5.0)], curve_type="PUMP2")
    assert [r[:3] for r in shared["PUMPS"]] == [["PT_P", "PT", "PT_D"]]
    conduits = {r[0]: r for r in shared["CONDUITS"]}
    assert conduits["9"][1] == conduits["10"][1] == "PT_D"

    separate = read_sections(
        write(sub, tmp_path, pump_stations={"PT": dict(pump=curve, discharge="separate")})
    )
    assert [r[:3] for r in separate["PUMPS"]] == [
        ["PT_P1", "PT", "PT_D1"],
        ["PT_P2", "PT", "PT_D2"],
    ]
    conduits = {r[0]: r for r in separate["CONDUITS"]}
    assert (conduits["9"][1], conduits["10"][1]) == ("PT_D1", "PT_D2")


def test_choosing_one_of_twin_force_mains(net, tmp_path):
    sub = net.subnetwork(["PT", "T1", "T2", "TO1", "TO2"])
    station = dict(pump=dict(curve=[(0, 5.0)], curve_type="PUMP2"), force_mains=9)
    with pytest.warns(UserWarning, match="force mains 10 are not pumped"):
        s = read_sections(write(sub, tmp_path, pump_stations={"PT": station}))
    conduits = {r[0]: r for r in s["CONDUITS"]}
    assert conduits["9"][1] == "PT_D" and conduits["10"][1] == "PT"
    assert {r[0]: r[1] for r in s["XSECTIONS"]}["10"] == "CIRCULAR"


# --- Roughness --------------------------------------------------------------------------------


def test_london_roughness_by_diameter(net, tmp_path):
    sub = net.subnetwork(["PT", "T1", "T2", "TO1", "TO2"])
    rule = LONDON_ON.design["force_main_roughness"]
    s = read_sections(
        write(
            sub,
            tmp_path,
            pump_stations={"PT": dict(IDEAL, discharge="shared")},
            force_main_roughness=rule,
        )
    )
    xsections = {r[0]: r for r in s["XSECTIONS"]}
    assert float(xsections["9"][3]) == 110.0  # 200 mm
    assert float(xsections["10"][3]) == 130.0  # 650 mm


@pytest.mark.parametrize(
    "diameter, c", [(0.1, 100), (0.15, 100), (0.2, 110), (0.25, 110), (0.6, 120), (0.9, 130)]
)
def test_london_table(diameter, c):
    rule = force_main_roughness_rule(LONDON_ON.design["force_main_roughness"])
    assert rule({"diameter": diameter}) == c


def test_roughness_by_material_callable_and_override(sub, tmp_path):
    stations = {"PS": dict(pump=dict(curve=[(0, 5.0)], curve_type="PUMP2"))}
    by_material = dict(by="material", table={"PVC": 140})
    s = read_sections(
        write(sub, tmp_path, pump_stations=stations, force_main_roughness=by_material)
    )
    assert {r[0]: r[3] for r in s["XSECTIONS"]}["2"] == "140.0"

    s = read_sections(
        write(
            sub,
            tmp_path,
            pump_stations=stations,
            force_main_roughness=lambda pipe: 1000 * pipe["diameter"],
            force_main_roughness_by_pipe={"3": 90},
        )
    )
    xsections = {r[0]: r[3] for r in s["XSECTIONS"]}
    assert (xsections["2"], xsections["3"]) == ("150.0", "90.0")


def test_missing_roughness_rejected(sub, tmp_path):
    with pytest.raises(ValueError, match=r"no roughness: 2 .*3 "):
        write(sub, tmp_path, pump_stations={"PS": IDEAL}, force_main_roughness=None)
    with pytest.raises(ValueError, match=r"no roughness: 2 \(diameter 0.15 m, material 'PVC'\)"):
        write(
            sub,
            tmp_path,
            pump_stations={"PS": IDEAL},
            force_main_roughness=dict(by="material", table={"CONC": 120}),
        )


def test_roughness_equation_must_match_template(sub, tmp_path):
    tpl = template(tmp_path)
    tpl.write_text(tpl.read_text().replace("FORCE_MAIN_EQUATION  H-W", "FORCE_MAIN_EQUATION  D-W"))
    with pytest.raises(ValueError, match=r"Hazen-Williams C .*FORCE_MAIN_EQUATION D-W"):
        write(
            sub,
            tmp_path,
            tpl=tpl,
            pump_stations={"PS": IDEAL},
            force_main_roughness=LONDON_ON.design["force_main_roughness"],
        )
    assert not (tmp_path / "model.inp").exists()


def test_invalid_roughness_rule_rejected():
    with pytest.raises(ValueError, match="'by'"):
        force_main_roughness_rule(dict(by="age", table={}))


# --- Validation ---------------------------------------------------------------------------------


def test_report(net):
    report = net.subnetwork(list(POS)).validate()
    assert report.stations_without_force_main == ["PS2"]
    assert report.orphan_force_mains == [8]
    assert report.force_main_ends == {
        8: ("Q2", "dead_end"),
        13: ("RJ", "force_main"),
        14: ("RJ", "force_main"),
    }


def test_force_main_chain(sub):
    assert sub.force_main_chain(2) == ([2, 3], "MH", None)


@pytest.mark.parametrize(
    "nodes, stations, match",
    [
        (["PS2", "M5", "OUT5"], {"PS2": IDEAL}, "PS2 has no force-main outlet pipe"),
        (MAIN, {"XX": IDEAL}, "XX is not in the subnetwork"),
        (MAIN, {"OUT": IDEAL}, "OUT is an outfall"),
        (["PT", "T1", "T2", "TO1", "TO2"], {"PT": IDEAL}, "PT has 2 force mains .*discharge="),
        (MAIN, {"PS": dict(IDEAL, force_mains=[5])}, "force mains 5 do not leave it"),
        (MAIN, {"PS": dict(IDEAL, discharge="both")}, "discharge must be"),
        (MAIN, {"PS": dict(pump=dict(curve=[(0, 1)], curve_type="PUMP1"))}, "PUMP1 .*no wet well"),
        (MAIN, {"PS": dict(pump=dict(curve="x"))}, r"Pump PS: curve"),
    ],
)
def test_invalid_station_rejected(net, tmp_path, nodes, stations, match):
    with pytest.raises(ValueError, match=match):
        write(net.subnetwork(nodes), tmp_path, pump_stations=stations)
    assert not (tmp_path / "model.inp").exists()


def test_ideal_pump_must_be_only_outlet(sub, tmp_path):
    """SWMM refuses an ideal pump whose inlet node has another outlet link (ERROR 134), here the
    overflow."""
    with pytest.raises(ValueError, match=r"only outlet link .*PS_P \(inlet PS, 2 outlet links\)"):
        write(sub, tmp_path, pump_stations={"PS": IDEAL})


def test_pumps_and_station_on_same_pipe_rejected(sub, tmp_path):
    with pytest.raises(ValueError, match="both in pumps and in a pump station's force main: 3"):
        write(sub, tmp_path, pump_stations={"PS": IDEAL}, pumps={3: dict(curve="*")})


def test_wet_well_given_twice_rejected(sub, tmp_path):
    well = dict(max_depth=3.0, area=2.0)
    with pytest.raises(ValueError, match="both in wet_wells and pump_stations: PS"):
        write(
            sub,
            tmp_path,
            pump_stations={"PS": dict(IDEAL, wet_well=well)},
            wet_wells={"PS": well},
        )


def test_station_names_must_be_free(tmp_path):
    """MH -> OUT renamed MH -> PS_D: the name of PS's discharge node is taken."""
    pipes = [p if p[0] != 4 else (4, "MH", "PS_D", *p[3:]) for p in PIPES]
    POS["PS_D"] = POS["OUT"]
    try:
        clash = make_network(pipes, {**TYPES, "PS_D": "OF"})
    finally:
        del POS["PS_D"]
    with pytest.raises(ValueError, match="PS: PS_D already exist"):
        write(
            clash.subnetwork(["U1", "PS", "F1", "MH", "PS_D"]),
            tmp_path,
            pump_stations={"PS": IDEAL},
        )


# --- In SWMM --------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "station",
    [
        dict(pump=dict(curve="*")),
        dict(
            pump=dict(curve=[(0.0, 10.0)], curve_type="PUMP2", startup=1.0, shutoff=0.2),
            wet_well=dict(max_depth=4.0, area=2.0),
        ),
    ],
    ids=["ideal", "pump2"],
)
def test_force_main_delays_arrival(net, tmp_path, station):
    """A constant 5 L/s at PS first reaches MH only once the pump has displaced most of the force
    main's volume (it starts empty and rises all the way).

    SWMM takes a conduit's volume from its mid-length depth, so water spills into MH when the
    upper pipe (3, 200 m) holds about half its volume: the lag is between the time to pump the
    volume of pipe 2 plus half of pipe 3 and the time to pump the whole force main."""
    from pyswmm import Links, Simulation

    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    path = write(
        sub, tmp_path, tpl=template(tmp_path, inflow=5.0, hours=3), pump_stations={"PS": station}
    )
    assert run_swmm(path, tmp_path) == []

    half = FORCE_MAIN_VOLUME * (300.0 + 200.0 / 2) / 500.0
    pumped = 0.0  # m3
    start = half_filled = filled = arrival = None
    with Simulation(str(path)) as sim:
        pump, outlet = Links(sim)["PS_P"], Links(sim)["4"]
        last = sim.start_time
        for _ in sim:
            now = sim.current_time
            seconds = (now - sim.start_time).total_seconds()
            pumped += pump.flow / 1000 * (now - last).total_seconds()
            last = now
            if start is None and pump.flow > 0:
                start = seconds
            if half_filled is None and pumped >= half:
                half_filled = seconds
            if filled is None and pumped >= FORCE_MAIN_VOLUME:
                filled = seconds
            if arrival is None and outlet.flow > 0.5:
                arrival = seconds
        routing_error = sim.flow_routing_error

    assert None not in (start, half_filled, filled, arrival)
    assert 0.95 * (half_filled - start) <= arrival - start <= 1.05 * (filled - start)
    assert abs(routing_error) < 1.0  # continuity error, %
