"""Dry-weather flow loading in the .inp writer: [DWF] and [PATTERNS] (#33)."""

import logging
import math

import pytest
from conftest import DENDRITIC, TEMPLATE, read_sections, run_swmm, subnetwork, write_inp

from swmmbuilder import dwf_from_population, swmm

DWF = {"A1": 1.0, "B1": 0.5, "C1": 0.75, "J1": 0.25}
PEAK_HOUR = 8
# A smooth diurnal pattern with mean 1, peaking at PEAK_HOUR and lowest 12 h later.
DIURNAL = [round(1 + 0.5 * math.cos(2 * math.pi * (h - PEAK_HOUR) / 24), 6) for h in range(24)]


@pytest.fixture
def sub(network):
    return subnetwork(network, DENDRITIC)


def write(sub, tmp_path, **kwargs):
    return write_inp(sub, tmp_path, routing="DYNWAVE", **kwargs)


def test_dwf_and_patterns_written(sub, tmp_path):
    weekly = [1.0, 1.0, 1.0, 1.0, 1.0, 1.1, 1.1]
    path = write(
        sub,
        tmp_path,
        dwf=DWF,
        patterns={"diurnal": DIURNAL, "weekly": weekly},
        dwf_patterns={"A1": ["diurnal", "weekly"], "B1": "diurnal"},
    )
    s = read_sections(path)
    rows = {r[0]: r for r in s["DWF"]}
    assert set(rows) == set(DWF)
    assert rows["A1"] == ["A1", "FLOW", "1.0", '"diurnal"', '"weekly"']
    assert rows["B1"] == ["B1", "FLOW", "0.5", '"diurnal"']
    assert rows["C1"] == ["C1", "FLOW", "0.75"]

    patterns = {}
    for name, *values in s["PATTERNS"]:
        if values[0] in swmm.PATTERN_LENGTHS:
            kind, *values = values
            patterns[name] = (kind, [])
        patterns[name][1].extend(float(v) for v in values)
    assert patterns == {"diurnal": ("HOURLY", DIURNAL), "weekly": ("DAILY", weekly)}
    assert "{" not in path.read_text()


def test_pattern_rows_name_each_line(sub, tmp_path):
    path = write(sub, tmp_path, dwf=DWF, patterns={"diurnal": DIURNAL}, dwf_patterns="diurnal")
    rows = read_sections(path)["PATTERNS"]
    assert len(rows) == 4 and all(r[0] == "diurnal" for r in rows)
    assert rows[0][1] == "HOURLY" and all(len(r) == 7 for r in rows[1:])
    assert {r[0] for r in read_sections(path)["DWF"]} == set(DWF)


def test_one_pattern_name_applies_to_every_node(sub):
    model = swmm.build_model(
        sub, ["OUT"], dwf=DWF, patterns={"diurnal": DIURNAL}, dwf_patterns="diurnal"
    )
    assert {d.node: d.patterns for d in model.dwf} == dict.fromkeys(DWF, ("diurnal",))
    assert model.patterns == [swmm.Pattern("diurnal", "HOURLY", tuple(DIURNAL))]


def test_explicit_pattern_type(sub):
    model = swmm.build_model(sub, ["OUT"], patterns={"wknd": ("weekend", DIURNAL)})
    assert model.patterns[0].type == "WEEKEND"


def test_no_dwf_needs_no_placeholders(network, tmp_path):
    template = tmp_path / "nodwf.inp"
    template.write_text(TEMPLATE.read_text().replace("{dwf}", "").replace("{patterns}", ""))
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1")])
    sub.to_inp(tmp_path / "m.inp", ["J1"], template_path=template)
    assert (tmp_path / "m.inp").exists()


@pytest.mark.parametrize("placeholder", ["dwf", "patterns"])
def test_template_without_dwf_placeholders_rejected(sub, tmp_path, placeholder):
    template = tmp_path / "nodwf.inp"
    text = TEMPLATE.read_text().replace(f"{{{placeholder}}}", "").replace("KINWAVE", "DYNWAVE")
    template.write_text(text)
    with pytest.raises(ValueError, match=rf"no \{{{placeholder}\}} placeholder"):
        sub.to_inp(
            tmp_path / "m.inp",
            ["OUT"],
            template_path=template,
            dwf=DWF,
            patterns={"diurnal": DIURNAL},
            dwf_patterns="diurnal",
        )
    assert not (tmp_path / "m.inp").exists()


def test_total_base_flow_logged(sub, tmp_path, caplog):
    with caplog.at_level(logging.INFO, logger="swmmbuilder.swmm"):
        write(sub, tmp_path, dwf=DWF)
    assert "total base flow 2.5 LPS at 4 nodes" in caplog.text


@pytest.mark.parametrize(
    "dwf, match",
    [
        ({"E1": 1.0}, "not in the subnetwork: E1"),
        ({"OUT": 1.0}, "outfalls: OUT"),
        ({"A1": -1.0}, "non-negative and finite: A1"),
        ({"A1": float("nan")}, "non-negative and finite: A1"),
        ({"A1": float("inf")}, "non-negative and finite: A1"),
        ({"A1": "lots"}, "non-negative and finite: A1"),
    ],
)
def test_invalid_dwf_rejected(sub, tmp_path, dwf, match):
    with pytest.raises(ValueError, match=match):
        write(sub, tmp_path, dwf=dwf)
    assert not (tmp_path / "model.inp").exists()


@pytest.mark.parametrize(
    "patterns, match",
    [
        ({"p": [1.0] * 5}, "p: 5 multipliers"),
        ({"p": ("HOURLY", [1.0] * 7)}, "p: a HOURLY pattern needs 24"),
        ({"p": ("YEARLY", [1.0] * 24)}, "p: type 'YEARLY'"),
        ({"p": [0.0] + [1.0] * 6}, "p: multipliers must be positive"),
        ({"p": [-1.0] + [1.0] * 6}, "p: multipliers must be positive"),
        ({"p": [float("nan")] + [1.0] * 6}, "p: multipliers must be positive"),
    ],
)
def test_invalid_pattern_rejected(sub, tmp_path, patterns, match):
    with pytest.raises(ValueError, match=match):
        write(sub, tmp_path, dwf=DWF, patterns=patterns)


@pytest.mark.parametrize(
    "dwf_patterns, match",
    [
        ("missing", "A1: patterns not defined in patterns: missing"),
        ({"A2": "diurnal"}, "no dry-weather flow: A2"),
        ({"A1": ["diurnal", "diurnal2"]}, "A1: at most one pattern of each type"),
    ],
)
def test_invalid_dwf_patterns_rejected(sub, tmp_path, dwf_patterns, match):
    with pytest.raises(ValueError, match=match):
        write(
            sub,
            tmp_path,
            dwf=DWF,
            patterns={"diurnal": DIURNAL, "diurnal2": DIURNAL},
            dwf_patterns=dwf_patterns,
        )


def test_hourly_pattern_mean_not_one_warns(sub):
    with pytest.warns(UserWarning, match=r"Pattern p: the mean .* is 2\.000"):
        swmm.build_model(sub, ["OUT"], patterns={"p": [2.0] * 24})


def test_daily_pattern_mean_not_checked(sub, recwarn):
    swmm.build_model(sub, ["OUT"], patterns={"p": [2.0] * 7})
    assert not [w for w in recwarn if "Pattern" in str(w.message)]


@pytest.mark.parametrize(
    "per_capita, unit, flow_units, expected",
    [
        (86400.0, "L/d", "LPS", 1.0),
        (86.4, "m3/d", "LPS", 1.0),
        (1.0, "L/s", "CMS", 0.001),
        (1.0, "m3/s", "LPS", 1000.0),
        (1.0, "L/s", "MLD", 0.0864),
    ],
)
def test_dwf_from_population(per_capita, unit, flow_units, expected):
    flows = dwf_from_population({"A1": 1, "B1": 3}, per_capita, unit, flow_units)
    assert flows == {"A1": pytest.approx(expected), "B1": pytest.approx(3 * expected)}


def test_dwf_from_population_rejects_unknown_units():
    with pytest.raises(ValueError, match="unit"):
        dwf_from_population({"A1": 1}, 250, "gal/d")
    with pytest.raises(ValueError, match="flow_units"):
        dwf_from_population({"A1": 1}, 250, "L/d", "CFS")


def test_dwf_runs_in_swmm(sub, tmp_path):
    """Over the second day, the outfall carries the summed base flow, following the pattern."""
    from pyswmm import Nodes, Simulation

    path = write(sub, tmp_path, dwf=DWF, patterns={"diurnal": DIURNAL}, dwf_patterns="diurnal")
    assert run_swmm(path, tmp_path) == []

    step = 300  # s
    flows = []  # outfall inflow at the end of each step of the second day
    with Simulation(str(path)) as sim:
        sim.step_advance(step)
        outfall = Nodes(sim)["OUT"]
        start = sim.start_time
        for _ in sim:
            elapsed = (sim.current_time - start).total_seconds()
            if 86400 < elapsed <= 2 * 86400:
                flows.append(outfall.total_inflow)
            elif elapsed > 2 * 86400:
                break

    assert len(flows) == 86400 // step
    total = sum(DWF.values())
    assert sum(flows) / len(flows) == pytest.approx(total, rel=0.01)

    per_hour = 3600 // step
    hourly = [sum(flows[h * per_hour : (h + 1) * per_hour]) / per_hour for h in range(24)]
    lag = (hourly.index(max(hourly)) - PEAK_HOUR) % 24
    assert 0 <= lag <= 2
    assert (hourly.index(min(hourly)) - (PEAK_HOUR + 12)) % 24 <= 2
    # Each hour's mean flow is close to the pattern's flow, at most an hour later.
    assert all(
        min(abs(hourly[h] - total * DIURNAL[(h - lag_) % 24]) for lag_ in (0, 1)) < 0.15 * total
        for h in range(24)
    )
