"""Pipe inverts outside the schema's invert_range are treated as missing (#105)."""

import dataclasses
from pathlib import Path

import geopandas as gpd
import pytest
from conftest import FIXTURE

import swmmbuilder
from swmmbuilder import LONDON_ON, DataWarning, Network, Schema

RANGE = (100.0, 120.0)


@pytest.fixture
def pipes():
    """The fixture pipes with an upstream invert keyed in without its decimal point (10500.0 for
    105.0) and a drop height (2.5) in the drop-structure invert column."""
    data = gpd.read_file(FIXTURE)
    first, second = data.index[:2]
    data.loc[first, "UpstreamInvert"] = data.loc[first, "UpstreamInvert"] * 100
    data.loc[second, "DropStructureInvert"] = 2.5
    return data, data.loc[first, "OBJECTID"], data.loc[second, "OBJECTID"]


def test_implausible_inverts_are_missing_and_reported(pipes):
    data, first, second = pipes
    schema = dataclasses.replace(LONDON_ON, invert_range=RANGE)
    with pytest.warns(DataWarning, match="2 pipe inverts outside 100-120 m treated as missing"):
        net = Network(data, schema=schema)
    full = net.full_data.set_index("pipe_id")
    assert full.isna().loc[first, "us_invert"]
    assert full.isna().loc[second, "drop_invert"]
    report = net.invert_report
    assert list(report.columns) == ["pipe_id", "column", "source_column", "value"]
    rows = {(r.pipe_id, r.column): (r.source_column, r.value) for r in report.itertuples()}
    assert rows[(first, "us_invert")][0] == "UpstreamInvert"
    assert rows[(first, "us_invert")][1] > 10000
    assert rows[(second, "drop_invert")] == ("DropStructureInvert", 2.5)


def test_plausible_and_null_inverts_are_kept(pipes):
    data, *_ = pipes
    schema = dataclasses.replace(LONDON_ON, invert_range=RANGE)
    with pytest.warns(DataWarning):
        net = Network(data, schema=schema)
    # Every other recorded invert is kept; the null drop-structure inverts (0) are not reported.
    assert len(net.invert_report) == 2
    plain = Network(data)
    assert (
        net.full_data["us_invert"].notna().sum() == plain.full_data["us_invert"].notna().sum() - 1
    )
    assert net.full_data["ds_invert"].equals(plain.full_data["ds_invert"])


def test_not_applied_without_a_range(pipes):
    data, first, _ = pipes
    net = Network(data)  # the London preset sets rim_range only
    assert LONDON_ON.invert_range is None
    assert net.invert_report.empty
    assert net.full_data.set_index("pipe_id").loc[first, "us_invert"] > 10000


def test_read_from_a_schema_file(tmp_path):
    text = "invert_range = [150.0, 350.0]\n" + LONDON_ON_TOML.read_text()  # before any table
    path = tmp_path / "s.toml"
    path.write_text(text)
    assert Schema.from_toml(path).invert_range == (150.0, 350.0)


LONDON_ON_TOML = Path(swmmbuilder.__file__).parent / "schemas" / "london_on.toml"


def test_only_the_pipes_kept_are_reported(pipes):
    """A pipe filtered out (here by flow type) is not in the report or its count (#123)."""
    data, first, second = pipes
    data.loc[data["OBJECTID"] == first, "FlowType"] = "STM"
    schema = dataclasses.replace(LONDON_ON, invert_range=RANGE)
    with pytest.warns(DataWarning, match="1 pipe inverts outside"):
        net = Network(data, schema=schema, flow_types=["SAN"])
    assert list(net.invert_report["pipe_id"]) == [second]
