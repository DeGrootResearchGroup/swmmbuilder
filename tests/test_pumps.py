"""Pumps and wet wells in the .inp writer (#13)."""

import pytest
from conftest import DENDRITIC, TEMPLATE, read_sections, run_swmm, subnetwork, write_inp

from sewernetworkx import Pump, WetWell

# D1 -> D2 is adverse in the GIS inverts (103.80 -> 104.10): a pumped rising main.
PUMPED = DENDRITIC + [("D1", "D2"), ("D2", "J2")]
CURVE = [(0.0, 6.0), (2.0, 4.0), (4.0, 0.0)]  # head (m) -> flow (L/s)


@pytest.fixture
def sub(network):
    return subnetwork(network, PUMPED)


def write(sub, tmp_path, **kwargs):
    return write_inp(sub, tmp_path, routing="DYNWAVE", **kwargs)


def test_pump_replaces_conduit(sub, tmp_path):
    path = write(sub, tmp_path, pumps={1010: dict(curve=CURVE, startup=0.5, shutoff=0.1)})
    s = read_sections(path)
    assert "1010" not in {r[0] for r in s["CONDUITS"]}
    assert "1010" not in {r[0] for r in s["XSECTIONS"]}
    assert s["PUMPS"] == [["1010", "D1", "D2", "1010_curve", "ON", "0.5", "0.1"]]
    rows = s["CURVES"]
    assert rows[0] == ["1010_curve", "PUMP3", "0.0", "6.0"]
    assert [(float(r[-2]), float(r[-1])) for r in rows] == CURVE
    assert all(r[0] == "1010_curve" for r in rows)
    assert "{" not in path.read_text()


def test_pumped_pipe_is_not_reported_adverse(sub, tmp_path, recwarn):
    write(sub, tmp_path, pumps={1010: dict(curve=CURVE)})
    assert not [w for w in recwarn if "adverse" in str(w.message)]


def test_pump_keys_match_as_strings(sub, tmp_path):
    s = read_sections(write(sub, tmp_path, pumps={"1010": Pump(curve="*")}))
    assert s["PUMPS"] == [["1010", "D1", "D2", "*", "ON", "0.0", "0.0"]]
    assert s.get("CURVES", []) == []


def test_pump_without_shape_is_accepted(network, tmp_path):
    """E1 -> J2 has no PipeShape: as a conduit it is rejected, as a pump it needs none."""
    sub = subnetwork(network, DENDRITIC + [("E1", "J2")])
    s = read_sections(write(sub, tmp_path, pumps={1012: dict(curve="*")}))
    assert [r[0] for r in s["PUMPS"]] == ["1012"]


def test_wet_well_written_as_storage(sub, tmp_path):
    path = write(
        sub,
        tmp_path,
        pumps={1010: dict(curve=[(0.0, 0.0), (1.0, 5.0)], curve_type="pump1")},
        wet_wells={"D1": WetWell(max_depth=3.0, area=4.0)},
    )
    s = read_sections(path)
    assert "D1" not in {r[0] for r in s["JUNCTIONS"]}
    ((row),) = s["STORAGE"]
    assert row[0] == "D1" and float(row[1]) == pytest.approx(103.80)
    assert [float(v) for v in row[2:4]] == [3.0, 0.0]
    assert row[4] == "FUNCTIONAL" and float(row[7]) == 4.0
    assert s["CURVES"][0][1] == "PUMP1"


def test_pump1_needs_wet_well(sub, tmp_path):
    with pytest.raises(ValueError, match=r"PUMP1 .*1010 \(inlet D1\)"):
        write(sub, tmp_path, pumps={1010: dict(curve=CURVE, curve_type="PUMP1")})
    assert not (tmp_path / "model.inp").exists()


@pytest.mark.parametrize(
    "pump, match",
    [
        (dict(curve=[]), "non-empty"),
        (dict(curve="pumpcurve"), "'\\*'"),
        (dict(curve=CURVE, curve_type="PUMP9"), "curve_type"),
        (dict(curve=CURVE, status="MAYBE"), "status"),
    ],
)
def test_invalid_pump_rejected(sub, tmp_path, pump, match):
    with pytest.raises(ValueError, match=f"Pump 1010: .*{match}"):
        write(sub, tmp_path, pumps={1010: pump})


def test_pump_not_in_subnetwork_rejected(sub, tmp_path):
    with pytest.raises(ValueError, match="not in the subnetwork: 1012"):
        write(sub, tmp_path, pumps={1012: dict(curve="*")})


@pytest.mark.parametrize(
    "wells, match",
    [
        ({"E1": dict(max_depth=3.0, area=4.0)}, "not in the subnetwork: E1"),
        ({"OUT": dict(max_depth=3.0, area=4.0)}, "outfalls: OUT"),
        ({"D1": dict(max_depth=0.0, area=4.0)}, "positive"),
    ],
)
def test_invalid_wet_well_rejected(sub, tmp_path, wells, match):
    with pytest.raises(ValueError, match=match):
        write(sub, tmp_path, wet_wells=wells)


def test_template_without_pumps_placeholder_rejected(sub, tmp_path):
    template = tmp_path / "nopumps.inp"
    template.write_text(TEMPLATE.read_text().replace("{pumps}", "").replace("KINWAVE", "DYNWAVE"))
    with pytest.raises(ValueError, match=r"no \{pumps\} placeholder"):
        sub.to_inp(
            tmp_path / "m.inp", ["OUT"], template_path=template, pumps={1010: dict(curve="*")}
        )
    assert not (tmp_path / "m.inp").exists()


def test_pumped_network_runs_in_swmm(sub, tmp_path):
    """A 2 L/s baseline inflow at the wet well D1 is pumped up to D2 and reaches the outfall."""
    from pyswmm import Links, Simulation

    template = tmp_path / "inflow.inp"
    text = TEMPLATE.read_text().replace("KINWAVE", "DYNWAVE")
    inflow = 'D1  FLOW  ""  FLOW  1.0  1.0  2.0\n'
    template.write_text(text.replace("[INFLOWS]\n", "[INFLOWS]\n" + inflow))
    path = tmp_path / "model.inp"
    sub.to_inp(
        path,
        ["OUT"],
        template_path=template,
        pumps={1010: dict(curve=CURVE, startup=0.3, shutoff=0.05)},
        wet_wells={"D1": dict(max_depth=3.0, area=2.0)},
    )
    assert run_swmm(path, tmp_path) == []

    pumped = 0.0
    with Simulation(str(path)) as sim:
        pump = Links(sim)["1010"]
        for _ in sim:
            pumped = max(pumped, pump.flow)
    assert pumped > 2.0  # the pump cycles on and off, above the inflow while on
