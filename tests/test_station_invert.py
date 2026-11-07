"""station_invert='lowest_inlet': a station's wet well is kept below the pipes arriving at it
(#103)."""

import pandas as pd
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


# ---- With wet_well_defaults, in the gap report, and in data_gaps() (#160) ---------------------

WET_WELL_RULE = dict(fill_time=10.0, high_level_below_inlet=0.2, low_level_above_pump=0.3)


def test_wet_well_defaults_take_the_estimated_inlet():
    # Pipe 1 arrives at 250.75 m by the slope, above the station's 249.8 m: the station stays,
    # and its wet well's inlet is the same estimated end, not "no gravity inlet".
    pipes = [(1, "U1", "PS", 251.0, None, 50.0, 300, "SL"), *REST]
    sub = make_network(pipes).subnetwork(NODES)
    options = dict(
        pump_stations={"PS": dict(pump=dict(curve="*"), design_flow=10.0)},
        force_main_roughness=120.0,
        wet_well_defaults=WET_WELL_RULE,
        headwater_slope=0.005,
    )
    with pytest.raises(ValueError, match="PS has no gravity inlet pipe with an invert"):
        build_model(sub, ["OUT"], **options)
    model = build_model(sub, ["OUT"], station_invert="lowest_inlet", **options)
    assert [(s.name, s.invert) for s in model.storage] == [("PS", pytest.approx(249.8))]
    rows = model.derived_report().query("rule == 'wet_well_defaults'")
    assert not rows.empty
    assert all(
        "inlet invert 250.75 (pipe 1, upstream invert - 0.005 x length)" in basis
        for basis in rows["basis"]
    )


# Station PS with no invert from any pipe: pipe 1 has no downstream invert, and the force main
# no upstream one. The rule gives it 249.35 m; the force main's last pipe has no downstream
# invert, which London's force_main_discharge rule fills.
NO_INVERT = [
    (1, "U1", "PS", 249.6, None, 50.0, 300, "SL"),
    (2, "PS", "F1", None, 251.0, 300.0, 150, "FM"),
    (3, "F1", "MH", 251.0, None, 200.0, 150, "FM"),
    *REST[2:],
]
NO_INVERT_OPTIONS = dict(
    pump_stations={"PS": IDEAL},
    station_invert="lowest_inlet",
    headwater_slope=0.005,
    design="schema",
)


def test_a_station_with_no_invert_is_reported_as_filled_by_the_rule(caplog):
    caplog.set_level("INFO", logger="swmmbuilder")
    sub = make_network(NO_INVERT).subnetwork(NODES)
    model = build_model(sub, ["OUT"], **NO_INVERT_OPTIONS)
    assert invert(model) == pytest.approx(249.35)
    gap = model.data_gaps.set_index(["id", "gap"]).loc[("PS", "no_node_invert")]
    assert gap["filled_from"] == "station_invert='lowest_inlet'"
    assert "pipe 1 (upstream invert - 0.005 x length)" in gap["source"]
    assert "no_node_invert filled by station_invert='lowest_inlet'" in caplog.text
    assert "filled by height_correction" not in caplog.text


def test_data_gaps_previews_design_and_station_invert():
    sub = make_network(NO_INVERT).subnetwork(NODES)
    model = build_model(sub, ["OUT"], **NO_INVERT_OPTIONS)
    preview = sub.data_gaps(**NO_INVERT_OPTIONS)
    pd.testing.assert_frame_equal(preview, model.data_gaps)
    filled = preview.set_index(["id", "gap"])["filled_from"]
    assert filled[(3, "no_ds_invert")] == "force_main_discharge"
    assert filled[("PS", "no_node_invert")] == "station_invert='lowest_inlet'"
    # Without the options, nothing is filled.
    assert sub.data_gaps()["filled_from"].isna().all()
