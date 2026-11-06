"""Pumping stations with force mains (#34): the pump at the station node, the force main kept as
full-flowing conduits."""

import math
import re
import warnings

import geopandas as gpd
import pytest
from conftest import TEMPLATE, read_sections, run_swmm
from shapely.geometry import LineString

from swmmbuilder import LONDON_ON, Network, PumpStation, Schema, WetWell
from swmmbuilder.swmm import build_model, force_main_roughness_rule

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


# Stations in series: SA => SB => SC -> SO, with a gravity sewer SG -> SB into station SB.
SERIES_POS = {"SA": (0, 0), "SB": (100, 0), "SC": (200, 0), "SO": (250, 0), "SG": (100, 50)}
SERIES_TYPES = {"SA": "PS", "SB": "PS", "SO": "OF"}
SERIES = [
    (30, "SA", "SB", 249.0, 251.0, 100.0, 200, "FM"),
    (31, "SB", "SC", 248.0, 251.0, 100.0, 200, "FM"),
    (32, "SC", "SO", 250.9, 250.0, 50.0, 300, "SL"),
    (33, "SG", "SB", 250.0, 249.5, 50.0, 300, "SL"),
]


def make_network(pipes=PIPES, types=TYPES, pos=POS, shapes=None):
    rows = [
        {
            "OBJECTID": pid,
            "FlowType": "SAN",
            "UpstreamID": up,
            "DownstreamID": down,
            "UpstreamInventoryType": types.get(up, "MH"),
            "DownstreamInventoryType": types.get(down, "MH"),
            "Length": length,
            "PipeShape": (shapes or {}).get(pid, "R"),
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
        LineString([(480000 + x, 4760000 + y) for x, y in (pos[up], pos[down])])
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


def test_diameter_rule_without_diameter():
    rule = force_main_roughness_rule(LONDON_ON.design["force_main_roughness"])
    assert rule({"diameter": float("nan")}) is None
    assert rule({"diameter": 0.2}) == 110


def test_non_circular_force_main_rejected(tmp_path):
    sub = make_network(shapes={3: "REC"}).subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    with pytest.raises(ValueError, match=r"not circular: 3 \(RECT_CLOSED\)"):
        build(sub, pump_stations={"PS": IDEAL})


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


# --- Joins and stations in series ---------------------------------------------------------------

JOIN = ["PR", "PR2", "RJ", "RZ", "RO"]


def test_force_main_chain_through_join(net):
    sub = net.subnetwork(JOIN)
    assert sub.force_main_chain(13) == ([13], "RJ", "force_main")
    assert sub.force_main_chain(13, through_joins=True) == ([13, 15], "RZ", "force_main")


def test_joined_force_mains_share_the_main_downstream(net, tmp_path):
    """#43: the main RJ => RZ after the join is a force main of both stations, RJ pressurised."""
    sub = net.subnetwork(JOIN)
    stations = {"PR": IDEAL, "PR2": IDEAL}
    s = read_sections(write(sub, tmp_path, pump_stations=stations))
    xsections = {r[0]: r for r in s["XSECTIONS"]}
    for pipe in ("13", "14", "15"):
        assert xsections[pipe][1:3] == ["FORCE_MAIN", "0.2" if pipe != "15" else "0.3"]
    junctions = {r[0]: r for r in s["JUNCTIONS"]}
    assert float(junctions["RJ"][4]) == 100.0
    assert float(junctions["RZ"][4]) == 0.0


def test_joined_force_mains_run_in_swmm(net, tmp_path):
    """20 L/s into each station of the join all reach the outfall: no flooding at RJ."""
    from pyswmm import Simulation, SystemStats

    sub = net.subnetwork(JOIN)
    text = template(tmp_path, hours=3).read_text()
    inflows = "".join(f'{n}  FLOW  ""  FLOW  1.0  1.0  20.0\n' for n in ("PR", "PR2"))
    tpl = tmp_path / "join.inp"
    tpl.write_text(text.replace("[INFLOWS]\n", "[INFLOWS]\n" + inflows))
    path = write(sub, tmp_path, tpl=tpl, pump_stations={"PR": IDEAL, "PR2": IDEAL})
    assert run_swmm(path, tmp_path) == []
    with Simulation(str(path)) as sim:
        for _ in sim:
            pass
        routing = SystemStats(sim).routing_stats
    inflow = 2 * 20.0 * 3 * 3600 / 1000  # m3
    assert routing["flooding"] == 0.0
    assert routing["outflow"] == pytest.approx(inflow, rel=0.05)  # m3, less the mains' fill
    assert abs(routing["routing_error"]) < 1.0  # %


def test_shared_main_after_join_checked_with_both_flows(net):
    """The 300 mm main after the join carries both stations' design flows: 2 x 10 L/s is
    0.28 m/s there, 10 L/s 0.32 m/s in each 200 mm main."""
    sub = net.subnetwork(JOIN)
    stations = {"PR": dict(IDEAL, design_flow=10.0), "PR2": dict(IDEAL, design_flow=10.0)}
    with pytest.warns(
        UserWarning,
        match=r"below 0.6 m/s: 13 \(0.32 m/s\), 14 \(0.32 m/s\), "
        r"15 \(0.28 m/s\)\.$",
    ):
        build(sub, pump_stations=stations, force_main_checks=LONDON_CHECKS)


@pytest.fixture(scope="module")
def series():
    net = make_network(SERIES, SERIES_TYPES, SERIES_POS)
    return net.subnetwork(list(SERIES_POS))


def test_force_main_chain_stops_at_station(series):
    """#44: SA's force main ends at station SB, whose own pump starts the next one."""
    assert series.force_main_chain(30) == ([30], "SB", None)


def test_stations_in_series(series, tmp_path):
    curve = dict(curve=[(0.0, 5.0)], curve_type="PUMP2", startup=1.0, shutoff=0.2)
    stations = {
        "SA": dict(pump=dict(curve=[(0.0, 30.0)], curve_type="PUMP2")),
        "SB": dict(pump=curve, wet_well=dict(max_depth=4.0, area=2.0)),
    }
    s = read_sections(write(series, tmp_path, pump_stations=stations))
    conduits = {r[0]: r for r in s["CONDUITS"]}
    assert conduits["30"][1:3] == ["SA_D", "SB"]
    assert conduits["31"][1:3] == ["SB_D", "SC"]
    assert "SB" in {r[0] for r in s["STORAGE"]}
    assert [r[:3] for r in s["PUMPS"]] == [["SA_P", "SA", "SA_D"], ["SB_P", "SB", "SB_D"]]


def test_station_in_series_checked_with_its_own_pump(series):
    """SB's force main 31 is checked against SB's 5 L/s (0.16 m/s), not SA's 30 L/s (0.95 m/s)."""
    stations = {
        "SA": dict(pump=dict(curve=[(0.0, 30.0)], curve_type="PUMP2")),
        "SB": dict(
            pump=dict(curve=[(0.0, 5.0)], curve_type="PUMP2"),
            wet_well=dict(max_depth=4.0, area=2.0),
        ),
    }
    with pytest.warns(UserWarning, match=r"velocity above 0.5 m/s: 30 \(0.95 m/s\)\.$"):
        build(series, pump_stations=stations, force_main_checks=dict(velocity=(0.1, 0.5)))


def test_undeclared_station_warned(series):
    with pytest.warns(UserWarning, match="not in pump_stations: SB; their force mains"):
        build(series, pump_stations={"SA": IDEAL})


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


# --- Wet-well defaults --------------------------------------------------------------------------

LONDON_WELL = LONDON_ON.design["wet_well_defaults"]
LONDON_CHECKS = LONDON_ON.design["force_main_checks"]
# PS's floor lowered to 247.0 m: its inlet (pipe 1, 300 mm) enters at 250.0 m, so London's rule
# gives a low water level 0.3 m and a high water level 250.0 - 0.3 - 247.0 = 2.7 m above the floor.
DEEP = {"PS": 247.0}
CURVE = dict(curve=[(0.0, 20.0)], curve_type="PUMP2")


def build(sub, **kwargs):
    kwargs.setdefault("force_main_roughness", C)
    return build_model(sub, **kwargs)


def test_london_presets():
    assert LONDON_WELL == dict(
        fill_time=30,
        high_level_below_inlet=0.3,
        low_level_above_pump=0.3,
        low_level_suction_diameters=2,
    )
    assert LONDON_CHECKS == dict(min_diameter=100, velocity=[0.6, 3.0])
    assert LONDON_ON.design["default_pump"] == dict(curve="*")


@pytest.mark.parametrize(
    "design, match",
    [
        (dict(force_main_checks=dict(velocity=3.0)), "force_main_checks: velocity must be"),
        (dict(wet_well_defaults=dict(fill_time=30)), "wet_well_defaults takes the keys"),
        (dict(default_pump=dict(curve=[(0, 1)])), "a curve needs a curve_type"),
        (dict(force_main_roughness=dict(by="age")), "design rule force_main_roughness"),
    ],
)
def test_invalid_preset_rejected_when_schema_loads(design, match):
    """#48: a typo in a preset is caught when the schema is made, not when it is used."""
    with pytest.raises(ValueError, match=match):
        Schema(name="x", columns=LONDON_ON.columns, design=design)


# --- Default pump and curve types (#54) ----------------------------------------------------------


def test_default_pump(net):
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    model = build(
        sub,
        pump_stations={"PS": PumpStation(wet_well=WetWell(max_depth=4.0, area=2.0))},
        default_pump=LONDON_ON.design["default_pump"],
    )
    assert model.pumps[0].curve == "*"


def test_station_pump_wins_over_default(sub):
    model = build(sub, pump_stations={"PS": dict(pump=CURVE)}, default_pump=IDEAL["pump"])
    assert model.pumps[0].curve == "PS_curve"


def test_station_without_pump_rejected(sub):
    with pytest.raises(ValueError, match="PS has no pump: pass its pump, or a default_pump"):
        build(sub, pump_stations={"PS": {}})


def test_station_curve_needs_curve_type(sub):
    with pytest.raises(ValueError, match="Pump PS: a curve needs a curve_type"):
        build(sub, pump_stations={"PS": dict(pump=dict(curve=[(0.0, 5.0)]))})


def test_low_level_from_suction_diameter(sub):
    """Two 200 mm suction diameters, 0.4 m, are above London's 0.3 m."""
    station = dict(pump=CURVE, design_flow=10.0, suction_diameter=200)
    model = build(
        sub, pump_stations={"PS": station}, wet_well_defaults=LONDON_WELL, height_correction=DEEP
    )
    assert (model.pumps[0].startup, model.pumps[0].shutoff) == (2.7, 0.4)
    (shutoff,) = [d for d in model.derived if d.attribute == "shutoff"]
    assert shutoff.basis.startswith("2 suction diameters (200 mm) above the floor")
    # 18 m3 between depths 0.4 and 2.7 m.
    assert model.storage[0].area == pytest.approx(round(18 / 2.3, 4))


def test_wet_well_from_rule(sub, tmp_path):
    """10 L/s for 30 min = 18 m3 between depths 0.3 and 2.7 m: 7.5 m2. The well rises to the crown
    of the inlet sewer, 250.3 m."""
    station = dict(pump=CURVE, design_flow=10.0)
    path = write(
        sub,
        tmp_path,
        pump_stations={"PS": station},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    s = read_sections(path)
    (storage,) = s["STORAGE"]
    assert storage[0] == "PS"
    assert float(storage[1]) == 247.0
    assert float(storage[2]) == pytest.approx(3.3)  # max depth
    assert float(storage[7]) == pytest.approx(7.5)  # area
    assert s["PUMPS"] == [["PS_P", "PS", "PS_D", "PS_curve", "ON", "2.7", "0.3"]]

    report = build(
        sub,
        pump_stations={"PS": station},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    ).derived_report()
    assert list(report.columns) == ["element", "id", "attribute", "value", "rule", "basis"]
    assert set(report["rule"]) == {"wet_well_defaults"}
    assert dict(zip(report["attribute"], report["value"], strict=True)) == pytest.approx(
        dict(area=7.5, max_depth=3.3, startup=2.7, shutoff=0.3)
    )
    assert "inlet invert 250.0 (pipe 1)" in report["basis"][0]


def test_wet_well_rule_is_logged(sub, caplog):
    caplog.set_level("INFO", logger="swmmbuilder.swmm")
    build(
        sub,
        pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    assert "Storage PS: area = 7.5 from rule wet_well_defaults, not data" in caplog.text
    assert "Pump PS_P: startup = 2.7 from rule wet_well_defaults" in caplog.text


def test_wet_well_rule_in_flow_units(sub):
    """The design flow is in the FLOW_UNITS: 0.864 MLD = 10 L/s."""
    model = build(
        sub,
        pump_stations={"PS": dict(pump=CURVE, design_flow=0.864)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
        flow_units="MLD",
    )
    assert model.storage[0].area == pytest.approx(7.5)


def test_explicit_values_win(sub):
    station = dict(
        pump=dict(CURVE, startup=2.0),
        wet_well=WetWell(max_depth=4.0, area=2.0),
        design_flow=10.0,
    )
    model = build(
        sub, pump_stations={"PS": station}, wet_well_defaults=LONDON_WELL, height_correction=DEEP
    )
    assert (model.storage[0].max_depth, model.storage[0].area) == (4.0, 2.0)
    assert (model.pumps[0].startup, model.pumps[0].shutoff) == (2.0, 0.3)
    assert [(d.attribute, d.value) for d in model.derived] == [("shutoff", 0.3)]


def test_wet_well_from_wet_wells_option_wins(sub):
    model = build(
        sub,
        pump_stations={"PS": dict(pump=CURVE)},
        wet_wells={"PS": dict(max_depth=4.0, area=2.0)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    assert model.storage[0].area == 2.0
    assert {d.attribute for d in model.derived} == {"startup", "shutoff"}


def test_nothing_to_derive(sub):
    """A station with a wet well and both levels is left alone, even with no gravity inlet."""
    station = dict(
        pump=dict(CURVE, startup=1.0, shutoff=0.3), wet_well=WetWell(max_depth=4.0, area=2.0)
    )
    model = build(
        sub.network.subnetwork(["PS", "F1", "MH", "OUT"]),
        pump_stations={"PS": station},
        wet_well_defaults=LONDON_WELL,
    )
    assert model.derived == []


def test_no_gravity_inlet_rejected(net):
    sub = net.subnetwork(["PS", "F1", "MH", "OUT"])
    with pytest.raises(ValueError, match="PS has no gravity inlet pipe with an invert"):
        build(
            sub,
            pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)},
            wet_well_defaults=LONDON_WELL,
        )


def test_wet_well_rises_to_the_rim():
    """With a rim above the inlet, 253.0 m, the well rises to it: 253.0 - 247.0 = 6.0 m."""
    sub = make_network().subnetwork(MAIN)  # its own network: the rim is set on the graph
    sub.graph.nodes["PS"]["rim_elevation"] = 253.0
    model = build(
        sub,
        pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    assert model.storage[0].max_depth == 6.0
    assert [d.basis.split(";")[0] for d in model.derived if d.attribute == "max_depth"] == [
        "rim - floor"
    ]


def test_derived_levels_name_each_separate_pump(net):
    """#53: with discharge='separate' each pump's derived levels name that pump."""
    pipes = PIPES + [(41, "TU", "PT", 250.5, 250.0, 50.0, 300, "SL")]
    sub = make_network(pipes, pos={**POS, "TU": (-50, 400)}).subnetwork(
        ["TU", "PT", "T1", "T2", "TO1", "TO2"]
    )
    station = dict(pump=CURVE, discharge="separate", wet_well=WetWell(max_depth=4.0, area=2.0))
    model = build(
        sub,
        pump_stations={"PT": station},
        wet_well_defaults=LONDON_WELL,
        height_correction={"PT": 247.0},
    )
    assert sorted((d.id, d.attribute) for d in model.derived) == [
        ("PT_P1", "shutoff"),
        ("PT_P1", "startup"),
        ("PT_P2", "shutoff"),
        ("PT_P2", "startup"),
    ]


def test_ideal_pump_gets_no_on_off_depths(net):
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    model = build(
        sub,
        pump_stations={"PS": dict(IDEAL, design_flow=10.0)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    assert model.storage[0].area == pytest.approx(7.5)
    assert (model.pumps[0].startup, model.pumps[0].shutoff) == (0.0, 0.0)
    assert {d.attribute for d in model.derived} == {"area", "max_depth"}


def test_no_rule_no_defaults(sub):
    model = build(sub, pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)})
    assert model.storage == [] and model.derived == []
    assert (model.pumps[0].startup, model.pumps[0].shutoff) == (0.0, 0.0)


@pytest.mark.parametrize(
    "station, correction, match",
    [
        (dict(pump=CURVE), DEEP, "PS has no wet well and no design_flow"),
        # The floor is 249.8 m: the high water level, 249.7 m, is below the low one, 250.1 m.
        (dict(pump=CURVE, design_flow=10.0), None, "PS: the high water level"),
    ],
)
def test_wet_well_rule_rejected(sub, station, correction, match):
    with pytest.raises(ValueError, match=match):
        build(
            sub,
            pump_stations={"PS": station},
            wet_well_defaults=LONDON_WELL,
            height_correction=correction,
        )


@pytest.mark.parametrize(
    "pump, well, match",
    [
        # #45: the derived startup, 2.7 m, is above the given well's max depth.
        ({}, WetWell(max_depth=2.0, area=2.0), "startup depth 2.7 m is not below the wet well's"),
        # #46: a given startup below the derived shutoff, 0.3 m, and the mirror case.
        (dict(startup=0.2), None, "startup depth 0.2 m is not above the shutoff depth 0.3 m"),
        (dict(shutoff=3.0), None, "startup depth 2.7 m is not above the shutoff depth 3.0 m"),
    ],
)
def test_derived_levels_checked_against_given_values(sub, pump, well, match):
    station = dict(pump=dict(CURVE, **pump), wet_well=well, design_flow=10.0)
    with pytest.raises(ValueError, match=match):
        build(
            sub,
            pump_stations={"PS": station},
            wet_well_defaults=LONDON_WELL,
            height_correction=DEEP,
        )


@pytest.mark.parametrize(
    "pump, match",
    [
        (dict(startup=0.5, shutoff=0.5), "startup depth 0.5 m is not above the shutoff"),
        (dict(startup=4.0, shutoff=0.5), "startup depth 4.0 m is not below the wet well's max"),
    ],
)
def test_given_levels_checked(sub, pump, match):
    station = dict(pump=dict(CURVE, **pump), wet_well=WetWell(max_depth=4.0, area=2.0))
    with pytest.raises(ValueError, match=match):
        build(sub, pump_stations={"PS": station})


def test_pipe_pump_levels_checked(sub):
    pump = dict(CURVE, startup=0.2, shutoff=0.4)
    with pytest.raises(ValueError, match=r"4: startup depth 0\.2 m is not above"):
        build(sub, pumps={4: pump})


PUMP1 = dict(curve=[(0.0, 0.0), (10.0, 10.0)], curve_type="PUMP1")


@pytest.mark.parametrize(
    "options",
    [
        dict(wet_wells={"PS": WetWell(max_depth=4.0, area=2.0)}),
        dict(wet_well_defaults=LONDON_WELL, height_correction=DEEP),
    ],
    ids=["wet_wells", "wet_well_defaults"],
)
def test_pump1_station_wet_well_from_options(sub, options):
    """#47: a PUMP1 station's wet well may come from wet_wells or the rule."""
    model = build(sub, pump_stations={"PS": dict(pump=PUMP1, design_flow=10.0)}, **options)
    assert [s.name for s in model.storage] == ["PS"]


def test_inlet_invert_ignores_force_mains_and_overflows(tmp_path):
    """#50: an undeclared force main into PS at 248.5 m does not set the inlet invert; the
    gravity sewer (pipe 1, 250.0 m) does."""
    pipes = PIPES + [(40, "X", "PS", 248.0, 248.5, 50.0, 150, "FM")]
    sub = make_network(pipes, pos={**POS, "X": (50, -50)}).subnetwork(MAIN + ["X"])
    model = build(
        sub,
        pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    assert (model.pumps[0].startup, model.pumps[0].shutoff) == (2.7, 0.3)
    assert "inlet invert 250.0 (pipe 1)" in model.derived[0].basis


@pytest.mark.parametrize(
    "rule, match",
    [
        (dict(fill_time=30), "missing keys"),
        (dict(LONDON_WELL, volume=1), "unknown keys"),
        (dict(LONDON_WELL, fill_time=0), "fill_time must be positive"),
        (30, "must be a dict"),
    ],
)
def test_invalid_wet_well_rule_rejected(sub, rule, match):
    with pytest.raises(ValueError, match=match):
        build(sub, pump_stations={"PS": dict(pump=CURVE)}, wet_well_defaults=rule)


def test_invalid_design_flow_rejected(sub):
    with pytest.raises(ValueError, match="design_flow must be positive"):
        build(sub, pump_stations={"PS": dict(pump=CURVE, design_flow=0)})


def test_flow_units_case_insensitive(sub, tmp_path):
    """#52"""
    write(sub, tmp_path, pump_stations={"PS": dict(pump=CURVE)}, flow_units="lps")
    model = build(sub, pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)}, flow_units="mld")
    assert model.pumps[0].name == "PS_P"


def test_station_and_pipe_curve_names_clash(tmp_path):
    """#51: pipe 7's pump curve and station 7's would both be 7_curve."""
    pipes = [
        (7, "U", "7", 250.5, 250.0, 50.0, 300, "SL"),
        (8, "7", "F", 249.8, 251.0, 100.0, 150, "FM"),
        (9, "F", "O", 250.9, 250.0, 50.0, 300, "SL"),
    ]
    pos = {"U": (0, 0), "7": (50, 0), "F": (150, 0), "O": (200, 0)}
    sub = make_network(pipes, {"7": "PS", "O": "OF"}, pos).subnetwork(list(pos))
    with pytest.raises(ValueError, match=r"Pump curve names .*: 7_curve"):
        build(sub, pumps={7: CURVE}, pump_stations={"7": dict(pump=CURVE)})


def test_flow_units_must_match_template(sub, tmp_path):
    with pytest.raises(ValueError, match="not the template's FLOW_UNITS LPS"):
        write(sub, tmp_path, pump_stations={"PS": dict(pump=CURVE)}, flow_units="CMS")


def test_derived_wet_well_runs_in_swmm(sub, tmp_path):
    """10 L/s into the derived 7.5 m2 well: the 20 L/s pump cycles between the derived depths,
    0.3 and 2.7 m (18 m3 fills in 30 min, empties in 30 min), with no flooding."""
    from pyswmm import Links, Nodes, Simulation, SystemStats

    path = write(
        sub,
        tmp_path,
        tpl=template(tmp_path, inflow=10.0, hours=3),
        pump_stations={"PS": dict(pump=CURVE, design_flow=10.0)},
        wet_well_defaults=LONDON_WELL,
        height_correction=DEEP,
    )
    assert run_swmm(path, tmp_path) == []
    depths, starts, running = [], 0, False
    with Simulation(str(path)) as sim:
        well, pump = Nodes(sim)["PS"], Links(sim)["PS_P"]
        for _ in sim:
            depths.append(well.depth)
            if pump.flow > 0 and not running:
                starts += 1
            running = pump.flow > 0
        routing = SystemStats(sim).routing_stats
    assert starts >= 2  # it cycles, about once an hour
    assert max(depths) == pytest.approx(2.7, abs=0.05)
    assert min(depths[len(depths) // 2 :]) == pytest.approx(0.3, abs=0.05)
    assert routing["flooding"] == 0.0
    # SWMM's continuity error with an on/off pump into an empty force main is about 2%.
    assert abs(routing["routing_error"]) < 2.5  # %


# --- Force-main checks ---------------------------------------------------------------------------

# The force main PS => F1 => MH is 150 mm: 1 L/s is 0.0566 m/s.


@pytest.mark.parametrize(
    "station, match",
    [
        (dict(pump=dict(curve=[(0.0, 5.0), (1.0, 10.0)], curve_type="PUMP2")), "below 0.6 m/s"),
        (dict(pump=dict(curve=[(0.0, 60.0)], curve_type="PUMP2")), "above 3 m/s"),
        (dict(IDEAL, design_flow=10.0), "below 0.6 m/s: 2 \\(0.57 m/s\\), 3 \\(0.57 m/s\\)"),
    ],
)
def test_force_main_velocity_warned(net, station, match):
    """The curve's largest flow, or an ideal pump's design flow, over the full area."""
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    with pytest.warns(UserWarning, match=match):
        build(sub, pump_stations={"PS": station}, force_main_checks=LONDON_CHECKS)


@pytest.mark.parametrize(
    "station",
    [dict(pump=dict(curve=[(0.0, 20.0)], curve_type="PUMP2")), IDEAL],
    ids=["1.13 m/s", "ideal, no design flow"],
)
def test_force_main_velocity_not_warned(net, station):
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build(sub, pump_stations={"PS": station}, force_main_checks=LONDON_CHECKS)


def test_force_main_diameter_warned(sub):
    with pytest.warns(UserWarning, match="narrower than 200 mm: 2 \\(150 mm\\), 3 \\(150 mm\\)"):
        build(sub, pump_stations={"PS": dict(pump=CURVE)}, force_main_checks=dict(min_diameter=200))


def test_shared_force_mains_velocity_not_checked(net, caplog):
    """PT's shared pump feeds 200 and 650 mm mains: the split is unknown."""
    caplog.set_level("INFO", logger="swmmbuilder.swmm")
    sub = net.subnetwork(["PT", "T1", "T2", "TO1", "TO2"])
    station = dict(IDEAL, discharge="shared", design_flow=1.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build(sub, pump_stations={"PT": station}, force_main_checks=LONDON_CHECKS)
    assert "Force-main velocity not checked" in caplog.text and ": 9, 10" in caplog.text


def test_separate_force_mains_checked_each(net):
    sub = net.subnetwork(["PT", "T1", "T2", "TO1", "TO2"])
    station = dict(pump=dict(curve=[(0.0, 50.0)], curve_type="PUMP2"), discharge="separate")
    # 50 L/s: 1.59 m/s in the 200 mm main, 0.15 m/s in the 650 mm one.
    with pytest.warns(UserWarning, match="below 0.6 m/s: 10 \\(0.15 m/s\\)\\.$"):
        build(sub, pump_stations={"PT": station}, force_main_checks=LONDON_CHECKS)


def test_no_checks_no_warnings(net):
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build(sub, pump_stations={"PS": dict(pump=dict(curve=[(0.0, 60.0)], curve_type="PUMP2"))})


def test_head_flow_curve_uses_design_flow(net, caplog):
    """#49: a PUMP3 curve's largest flow (60 L/s, at zero head) is not its flow against the
    force main: the check uses the design flow, or is skipped without one."""
    caplog.set_level("INFO", logger="swmmbuilder.swmm")
    sub = net.subnetwork(["U1", "PS", "F1", "MH", "OUT"])
    pump = dict(curve=[(0.0, 60.0), (10.0, 0.0)], curve_type="PUMP3")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        build(sub, pump_stations={"PS": dict(pump=pump)}, force_main_checks=LONDON_CHECKS)
    assert "Force-main velocity not checked" in caplog.text
    with pytest.warns(UserWarning, match=r"below 0\.6 m/s: 2 \(0\.57 m/s\)"):
        build(
            sub,
            pump_stations={"PS": dict(pump=pump, design_flow=10.0)},
            force_main_checks=LONDON_CHECKS,
        )


@pytest.mark.parametrize(
    "rule, match",
    [
        (dict(speed=1), "unknown keys"),
        (dict(velocity=3.0), "velocity must be \\(min, max\\)"),
        (dict(velocity=[0.6]), "velocity must be"),
        (dict(velocity=["a", "b"]), "velocity must be"),
        (dict(min_diameter="100"), "min_diameter must be a positive number"),
        ([("velocity", 1)], "must be a dict"),
    ],
)
def test_invalid_checks_rejected(sub, rule, match):
    """#48"""
    with pytest.raises(ValueError, match=match):
        build(sub, pump_stations={"PS": dict(pump=CURVE)}, force_main_checks=rule)


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
    if station["pump"]["curve"] == "*":
        # The ideal pump passes the 5 L/s on at once: the same bounds from V / Q.
        assert start == pytest.approx(0.0, abs=60)
        assert 0.95 * half / 0.005 <= arrival <= 1.05 * FORCE_MAIN_VOLUME / 0.005
    assert abs(routing_error) < 1.0  # continuity error, %


def test_discharge_node_given_the_rule_depth_is_derived(sub):
    """A station's discharge node takes the max_depth rule's depth and is listed in derived
    (#123)."""
    model = build_model(
        sub,
        pump_stations={"PS": dict(pump=dict(curve=[(0, 5.0)], curve_type="PUMP2"))},
        force_main_roughness=C,
        max_depth=lambda d: 2.5,
    )
    (node,) = [j for j in model.junctions if j.name == "PS_D"]
    assert node.max_depth == 2.5
    rows = model.derived_report()
    assert set(rows.loc[rows["id"] == "PS_D", "rule"]) == {"max_depth function"}
