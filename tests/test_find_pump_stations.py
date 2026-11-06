"""Subnetwork.find_pump_stations: stations to pump, stations without a force main, and force
mains that climb from a node that is not a station (#104)."""

import pytest
from test_pump_stations import PIPES, make_network

from swmmbuilder import PumpStationReport


@pytest.fixture(scope="module")
def net():
    return make_network()


def everything(net):
    return net.subnetwork(list(net.graph.nodes))


def test_stations_and_stations_without_a_force_main(net):
    report = everything(net).find_pump_stations()
    assert isinstance(report, PumpStationReport)
    assert set(report.stations) == {"PS", "PT", "PR", "PR2"}
    assert report.no_force_main == ["PS2"]


def test_a_climbing_force_main_from_no_station(net):
    # Q1 => Q2 climbs 250.0 -> 251.0; F1 and RJ are inside force mains, so not reported.
    unlisted = everything(net).find_pump_stations().unlisted
    assert unlisted["node"].tolist() == ["Q1"]
    row = unlisted.iloc[0]
    assert row["pipes"] == [8] and row["end"] == "Q2"
    assert (row["us_invert"], row["ds_invert"], row["rise"]) == pytest.approx((250.0, 251.0, 1.0))


def test_min_rise(net):
    assert everything(net).find_pump_stations(min_rise=1.0).unlisted.empty


def test_a_falling_force_main_is_not_reported():
    pipes = [p if p[0] != 8 else (8, "Q1", "Q2", 251.0, 250.0, 50.0, 150, "FM") for p in PIPES]
    sub = make_network(pipes).subnetwork(["Q1", "Q2"])
    report = sub.find_pump_stations()
    assert report.unlisted.empty and report.stations == []


def test_only_the_subnetwork(net):
    report = net.subnetwork(["U1", "PS", "F1", "MH", "OUT", "OVO"]).find_pump_stations()
    assert report.stations == ["PS"] and report.no_force_main == []
    assert report.unlisted.empty
    assert list(report.unlisted.columns) == [
        "node",
        "pipes",
        "end",
        "us_invert",
        "ds_invert",
        "rise",
    ]
