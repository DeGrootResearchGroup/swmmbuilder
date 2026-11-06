"""pump_stations='found': every station with a force main, with exclusions and additions (#135)."""

import warnings

import pytest
from test_pump_stations import make_network, template

from swmmbuilder import DataWarning, OptionWarning, PumpStation
from swmmbuilder.swmm import build_model

# Station PS with a force main to OUT; station PS2 with no force main (a pass-through); Q1 => Q2,
# a force main climbing 1 m from Q1, which is not in the station layer.
NODES = ["U1", "PS", "F1", "MH", "OUT", "PS2", "M5", "OUT5", "Q1", "Q2"]
OUTFALLS = ["OUT", "OUT5", "Q2"]
NOT_PUMPED = {"PS2": "gravity pass-through"}
EXTRA_STATIONS = {"Q1": "force main climbs 1 m; no station in the layer"}
IDEAL = dict(curve="*")
IDEAL_STATION = dict(pump=IDEAL)


@pytest.fixture(scope="module")
def sub():
    return make_network().subnetwork(NODES)


def build(sub, **kwargs):
    kwargs.setdefault("force_main_roughness", 120.0)
    return build_model(sub, OUTFALLS, **kwargs)


def test_found_matches_the_explicit_loop(sub, tmp_path):
    """The loop of a build script: found stations minus exclusions, plus additions, all ideal."""
    tpl = template(tmp_path)
    found = sub.find_pump_stations()
    stations = {n: PumpStation(pump=IDEAL) for n in found.stations if n not in NOT_PUMPED}
    for n in EXTRA_STATIONS:
        stations[n] = PumpStation(pump=IDEAL)
    explicit = tmp_path / "explicit.inp"
    common = dict(template_path=tpl, force_main_roughness=120.0)
    sub.to_inp(explicit, OUTFALLS, pump_stations=stations, **common)
    by_option = tmp_path / "found.inp"
    with warnings.catch_warnings():
        warnings.filterwarnings("error", "Force mains climb", DataWarning)
        model = sub.to_inp(
            by_option,
            OUTFALLS,
            pump_stations="found",
            pump_stations_exclude=NOT_PUMPED,
            pump_stations_extra=EXTRA_STATIONS,
            default_pump=IDEAL,
            **common,
        )
    assert by_option.read_text() == explicit.read_text()
    assert model.station_choices.to_dict("records") == [
        {"node": "PS", "source": "found", "reason": ""},
        {"node": "Q1", "source": "extra", "reason": EXTRA_STATIONS["Q1"]},
        {"node": "PS2", "source": "excluded", "reason": "gravity pass-through"},
    ]
    s = model.summary()
    assert (s["stations"], s["stations_extra"], s["stations_excluded"]) == (2, 1, 1)


def test_choices_are_logged(sub, caplog):
    caplog.set_level("INFO", logger="swmmbuilder.swmm")
    build(sub, pump_stations="found", pump_stations_extra=EXTRA_STATIONS, default_pump=IDEAL)
    assert "Pump station PS: found" in caplog.text
    assert f"Pump station Q1: extra ({EXTRA_STATIONS['Q1']})" in caplog.text
    assert "Pump station PS2 has no force main: written as a junction" in caplog.text


def test_unlisted_force_main_warned(sub):
    with pytest.warns(DataWarning, match=r"pumped by gravity.*Q1 \(force main 8 climbs 1 m\)"):
        model = build(sub, pump_stations="found", default_pump=IDEAL)
    assert {p.inlet for p in model.pumps} == {"PS"}


def test_unlisted_force_main_raise_or_ignore(sub):
    with pytest.raises(ValueError, match="Q1 \\(force main 8 climbs 1 m\\)"):
        build(sub, pump_stations="found", default_pump=IDEAL, pump_stations_unlisted="raise")
    with warnings.catch_warnings():
        warnings.filterwarnings("error", "Force mains climb", DataWarning)
        build(sub, pump_stations="found", default_pump=IDEAL, pump_stations_unlisted="ignore")
        # An unlisted node excluded with a reason is acknowledged, so not warned of.
        model = build(
            sub,
            pump_stations="found",
            default_pump=IDEAL,
            pump_stations_exclude={"Q1": "checked: the main is drawn the wrong way"},
        )
    assert list(model.station_choices["source"]) == ["found", "excluded"]


def test_excluded_station_not_warned_as_undeclared(sub):
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptionWarning)
        model = build(
            sub,
            pump_stations="found",
            default_pump=IDEAL,
            pump_stations_exclude={"PS": "outlet station", **NOT_PUMPED},
            pump_stations_unlisted="ignore",
        )
    assert model.pumps == []
    assert model.summary()["stations_excluded"] == 2


def test_extra_pump_station_used_as_given(sub):
    curve = dict(curve=[(0.0, 20.0)], curve_type="PUMP2")
    wet_well = dict(max_depth=3.0, area=4.0)
    model = build(
        sub,
        pump_stations="found",
        default_pump=IDEAL,
        pump_stations_extra={"PS": PumpStation(pump=curve, wet_well=wet_well), "Q1": IDEAL_STATION},
    )
    assert {p.inlet: p.curve for p in model.pumps} == {"PS": "PS_curve", "Q1": "*"}
    assert list(model.station_choices["source"]) == ["extra", "extra"]


def test_report_accepted(sub):
    report = sub.find_pump_stations(min_rise=2.0)  # Q1's 1 m climb is not reported
    with warnings.catch_warnings():
        warnings.filterwarnings("error", "Force mains climb", DataWarning)
        model = build(sub, pump_stations=report, default_pump=IDEAL)
    assert {p.inlet for p in model.pumps} == {"PS"}


def test_found_needs_a_default_pump(sub):
    with pytest.raises(ValueError, match="stations given without a pump need a default_pump: PS"):
        build(sub, pump_stations="found", pump_stations_unlisted="ignore")


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (dict(pump_stations="all"), "pump_stations must be 'found'"),
        (dict(pump_stations={"PS": {}}, pump_stations_exclude=NOT_PUMPED), "need pump_stations="),
        (dict(pump_stations_exclude={"M5": "no"}), "M5 is excluded but is not a station"),
        (
            dict(pump_stations_exclude={"Q1": "x"}, pump_stations_extra={"Q1": "y"}),
            "Q1 is both excluded and added",
        ),
        (dict(pump_stations_unlisted="quiet"), "pump_stations_unlisted must be"),
    ],
)
def test_invalid_choices_rejected(sub, kwargs, match):
    kwargs = {"pump_stations": "found", "default_pump": IDEAL, **kwargs}
    with pytest.raises(ValueError, match=match):
        build(sub, **kwargs)


def test_dict_unchanged(sub):
    model = build(sub, pump_stations={"PS": dict(pump=IDEAL)}, pump_stations_unlisted="ignore")
    assert model.station_choices is None
    assert {p.inlet for p in model.pumps} == {"PS"}
