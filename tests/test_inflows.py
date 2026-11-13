"""External inflows and time series in the .inp writer: [INFLOWS] and [TIMESERIES] (#98)."""

import math
from itertools import pairwise

import pytest
from conftest import DENDRITIC, TEMPLATE, read_sections, run_swmm, subnetwork, write_inp

from swmmbuilder import Inflow, hourly_series, swmm

PEAK_HOUR = 10
COMMERCIAL = [round(1 + 0.6 * math.cos(2 * math.pi * (h - PEAK_HOUR) / 24), 6) for h in range(24)]


@pytest.fixture
def sub(network):
    return subnetwork(network, DENDRITIC)


def write(sub, tmp_path, **kwargs):
    return write_inp(sub, tmp_path, routing="DYNWAVE", **kwargs)


def test_inflows_and_timeseries_written(sub, tmp_path):
    path = write(
        sub,
        tmp_path,
        patterns={"flat": [1.0] * 24},
        timeseries={"commercial": hourly_series(COMMERCIAL, 3), "step": [(0, 1.0), (2.5, 3.0)]},
        inflows={
            "A1": Inflow("commercial", scale=2.5, baseline=0.4),
            "B1": dict(baseline=0.25),
            "C1": dict(timeseries="step", baseline=0.1, baseline_pattern="flat"),
        },
    )
    s = read_sections(path)
    rows = {r[0]: r for r in s["INFLOWS"]}
    assert rows["A1"] == ["A1", "FLOW", "commercial", "FLOW", "1.0", "2.5", "0.4"]
    assert rows["B1"] == ["B1", "FLOW", '""', "FLOW", "1.0", "1.0", "0.25"]
    assert rows["C1"] == ["C1", "FLOW", "step", "FLOW", "1.0", "1.0", "0.1", "flat"]

    series = {}
    for name, t, v in s["TIMESERIES"]:
        series.setdefault(name, []).append((float(t), float(v)))
    assert series["step"] == [(0.0, 1.0), (2.5, 3.0)]
    assert len(series["commercial"]) == 2 * 24 * 3 + 1
    assert "{" not in path.read_text()


def test_hourly_series_holds_each_hour():
    points = hourly_series(COMMERCIAL, 2)
    assert points[0] == (0.0, COMMERCIAL[0])
    assert points[1] == pytest.approx((59 / 60, COMMERCIAL[0]))
    assert points[2] == (1.0, COMMERCIAL[1])
    assert points[48] == (24.0, COMMERCIAL[0])  # day 2 starts again
    assert points[-1] == (48.0, COMMERCIAL[0])  # closing point at the end
    times = [t for t, _ in points]
    assert all(b > a for a, b in pairwise(times))
    assert hourly_series(COMMERCIAL, 1, ramp_minutes=6)[1] == pytest.approx((0.9, COMMERCIAL[0]))


@pytest.mark.parametrize(
    "args, message",
    [
        (([1.0] * 23, 1), "24 finite"),
        (([1.0] * 24, 0), "positive whole"),
        (([1.0] * 24, 1.5), "positive whole"),
        (([1.0] * 24, 1, 0), "ramp_minutes"),
    ],
)
def test_hourly_series_refuses(args, message):
    with pytest.raises(ValueError, match=message):
        hourly_series(*args)


def test_model_records(sub):
    model = swmm.build_model(
        sub,
        ["OUT"],
        timeseries={"ts": [(0, 1.0), (1, 2.0)]},
        inflows={"A1": dict(timeseries="ts", scale=3, baseline=0.5)},
    )
    assert model.inflows == [swmm.ExternalInflow("A1", "ts", 3.0, 0.5, None)]
    assert model.timeseries == [swmm.TimeSeries("ts", ((0.0, 1.0), (1.0, 2.0)))]


def test_small_scale_and_baseline_are_written_in_full(sub):
    """The scale of a patterned flow (its mean, from pattern_inflows) is a flow, so a small one
    in CMS keeps its significant figures, as the baseline does."""
    model = swmm.build_model(
        sub,
        ["OUT"],
        timeseries={"ts": [(0, 1.0), (1, 2.0)]},
        inflows={"A1": dict(timeseries="ts", scale=3.4e-6, baseline=4e-7)},
    )
    assert model.sections()["inflows"].split()[5:7] == ["3.4e-06", "4e-07"]


@pytest.mark.parametrize(
    "options, message",
    [
        (dict(inflows={"NOPE": dict(baseline=1.0)}), "NOPE: not in the subnetwork"),
        (dict(inflows={"OUT": dict(baseline=1.0)}), "OUT: an outfall"),
        (dict(inflows={"A1": dict(timeseries="missing")}), "'missing' not in timeseries"),
        (dict(inflows={"A1": dict(baseline=1, baseline_pattern="p")}), "'p' not in patterns"),
        (dict(inflows={"A1": dict(baseline=-1.0)}), "baseline must be non-negative"),
        (dict(inflows={"A1": dict(scale=float("nan"))}), "scale must be non-negative"),
        (dict(inflows={"A1": dict(bseline=1.0)}), "A1: .*bseline"),
        (dict(timeseries={"ts": []}), "ts: no points"),
        (dict(timeseries={"ts": [(1, 1.0), (1, 2.0)]}), "increasing"),
        (dict(timeseries={"ts": [(-1, 1.0)]}), "at least 0 h"),
        (dict(timeseries={"ts": [(0, "x")]}), "finite numbers"),
        (dict(timeseries={"ts": [1, 2]}), "pairs"),
    ],
)
def test_invalid_options(sub, options, message):
    with pytest.raises(ValueError, match=message):
        swmm.build_model(sub, ["OUT"], **options)


def test_inflows_must_be_a_mapping(sub):
    with pytest.raises(TypeError, match="mapping"):
        swmm.build_model(sub, ["OUT"], inflows=["A1"])


def test_template_without_placeholders(network, tmp_path):
    template = tmp_path / "noinflows.inp"
    template.write_text(TEMPLATE.read_text().replace("{inflows}", "").replace("{timeseries}", ""))
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1")])
    sub.to_inp(tmp_path / "plain.inp", ["J1"], template_path=template)  # not needed: fine
    with pytest.raises(ValueError, match=r"\{inflows\}, \{timeseries\} placeholder"):
        sub.to_inp(
            tmp_path / "m.inp",
            ["J1"],
            template_path=template,
            timeseries={"ts": [(0, 1.0)]},
            inflows={"A1": dict(timeseries="ts")},
        )
    assert not (tmp_path / "m.inp").exists()


def test_swmm_runs_and_delivers_the_inflow(sub, tmp_path):
    """SWMM reads the sections, and the outflow is the steady baseline plus the series' mean
    times its scale (the run is three days; day 3 is compared)."""
    pytest.importorskip("pyswmm")
    from pyswmm import Nodes, Simulation

    scale, baseline = 4.0, 1.5
    path = write(
        sub,
        tmp_path,
        timeseries={"commercial": hourly_series(COMMERCIAL, 3)},
        inflows={"A1": dict(timeseries="commercial", scale=scale), "B1": dict(baseline=baseline)},
    )
    assert run_swmm(path, tmp_path) == []
    flows = []
    with Simulation(str(tmp_path / path.name)) as sim:
        sim.step_advance(300)
        outfall = Nodes(sim)["OUT"]
        for _ in sim:
            if sim.current_time >= sim.start_time.replace(day=sim.start_time.day + 2):
                flows.append(outfall.total_inflow)
    mean = sum(flows) / len(flows)
    assert mean == pytest.approx(baseline + scale * sum(COMMERCIAL) / 24, rel=0.02)
