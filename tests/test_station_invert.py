"""station_invert='lowest_inlet': a station's wet well is kept below the pipes arriving at it
(#103)."""

import pytest
from test_pump_stations import IDEAL, PIPES, make_network

from swmmbuilder.swmm import build_model

# Station PS: its force main leaves at 249.8 m. Pipe 1 (U1 -> PS, 50 m) arrives with its upstream
# invert at 249.6 m and no downstream invert, so at a slope of 0.005 it would arrive at 249.35 m.
# The station's overflow (pipe 5) is left out: an ideal pump must be its station's only outlet.
REST = [p for p in PIPES[1:] if p[0] != 5]
NODES = ["U1", "PS", "F1", "MH", "OUT"]
ARRIVING_NO_END = [(1, "U1", "PS", 249.6, None, 50.0, 300, "SL"), *REST]


def build(pipes=ARRIVING_NO_END, **options):
    sub = make_network(pipes).subnetwork(NODES)
    return build_model(
        sub,
        ["OUT"],
        pump_stations={"PS": IDEAL},
        force_main_roughness=120.0,
        **options,
    )


def invert(model, node="PS"):
    return next(j.invert for j in model.junctions if j.name == node)


def test_lowered_to_the_extrapolated_inlet():
    model = build(station_invert="lowest_inlet", headwater_slope=0.005)
    assert invert(model) == pytest.approx(249.35)
    rows = model.derived_report().query("element == 'station'")
    assert rows["id"].tolist() == ["PS"]
    assert rows["value"].iloc[0] == pytest.approx(249.35)
    assert "pipe 1 (upstream invert - 0.005 x length)" in rows["basis"].iloc[0]
    assert "249.800" in rows["basis"].iloc[0]


def test_not_applied_by_default_or_without_a_slope():
    assert invert(build(headwater_slope=0.005)) == pytest.approx(249.8)
    assert invert(build(station_invert="lowest_inlet")) == pytest.approx(249.8)


def test_a_recorded_arriving_invert_already_counts():
    # With the downstream invert recorded, the node invert already takes it: nothing to change.
    pipes = [(1, "U1", "PS", 249.6, 249.4, 50.0, 300, "SL"), *REST]
    model = build(pipes, station_invert="lowest_inlet", headwater_slope=0.005)
    assert invert(model) == pytest.approx(249.4)
    assert model.derived_report().query("element == 'station'").empty


def test_a_station_already_below_is_unchanged():
    pipes = [(1, "U1", "PS", 251.0, None, 50.0, 300, "SL"), *REST]  # arrives at 250.75
    model = build(pipes, station_invert="lowest_inlet", headwater_slope=0.005)
    assert invert(model) == pytest.approx(249.8)


def test_height_correction_wins():
    model = build(
        station_invert="lowest_inlet", headwater_slope=0.005, height_correction={"PS": 249.0}
    )
    assert invert(model) == pytest.approx(249.0)
    assert model.derived_report().query("element == 'station'").empty


def test_unknown_rule():
    with pytest.raises(ValueError, match="station_invert must be 'lowest_inlet'"):
        build(station_invert="lowest")
