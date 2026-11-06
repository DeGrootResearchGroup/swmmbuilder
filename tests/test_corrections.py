import logging

import geopandas as gpd
import pytest
from conftest import FIXTURE

from swmmbuilder import Network, PipeCorrections


def proposed(*pipes):
    """The fixture, with ``pipes`` marked Proposed in ConstructedStatus."""
    data = gpd.read_file(FIXTURE)
    data["ConstructedStatus"] = "Built"
    data.loc[data["OBJECTID"].isin(pipes), "ConstructedStatus"] = "Proposed"
    return data


def test_no_corrections():
    network = Network(FIXTURE)
    assert network.corrections_report.empty
    assert list(network.corrections_report.columns) == [
        "pipe_id",
        "action",
        "column",
        "old",
        "new",
        "reason",
    ]
    assert network.removed["corrections"] == {}


def test_drop():
    network = Network(FIXTURE, corrections={"drop": {1009: "overflow tail"}})
    assert 1009 not in set(network.full_data["pipe_id"])
    assert not network.graph.has_edge("OV1", "OV2")
    assert network.removed["corrections"] == {"dropped": 1}
    row = network.corrections_report.iloc[0]
    assert (row["pipe_id"], row["action"], row["reason"]) == (1009, "drop", "overflow tail")


def test_status_override_passes_the_status_filter():
    assert not Network(proposed(1000)).graph.has_edge("A1", "A2")
    network = Network(
        proposed(1000),
        corrections={"set": {1000: {"status": "built", "reason": "built since"}}},
    )
    assert network.graph.has_edge("A1", "A2")
    assert network.removed["status"] == {}
    row = network.corrections_report.iloc[0]
    assert row[["action", "column", "old", "new", "reason"]].tolist() == [
        "set",
        "status",
        "proposed",
        "built",
        "built since",
    ]


def test_setting_abandoned_status_drops_the_pipe():
    network = Network(FIXTURE, corrections={"set": {1000: {"status": "abandoned"}}})
    assert not network.graph.has_edge("A1", "A2")
    assert network.removed["status"] == {"abandoned": 2}


def test_same_as_editing_the_source():
    data = proposed(1000)
    data.loc[data["OBJECTID"] == 1000, "ConstructedStatus"] = "Built"
    edited = Network(data)
    corrected = Network(proposed(1000), corrections={"set": {1000: {"status": "built"}}})
    assert set(edited.graph.edges(keys=True)) == set(corrected.graph.edges(keys=True))
    columns = [c for c in edited.full_data.columns if c != "ConstructedStatus"]
    assert edited.full_data[columns].equals(corrected.full_data[columns])


def test_set_invert_and_node():
    network = Network(
        FIXTURE,
        normalize_ids=str.upper,
        corrections={"set": {1012: {"us_invert": 104.5, "to_node": "j1"}}},
    )
    pipe = network.full_data.set_index("pipe_id").loc[1012]
    assert pipe["us_invert"] == 104.5
    assert pipe["to_node"] == "J1"
    assert network.graph.has_edge("E1", "J1")
    assert network.corrections_report["column"].tolist() == ["us_invert", "to_node"]


def test_swap_nodes_keeps_inverts_and_geometry():
    before = Network(FIXTURE).full_data.set_index("pipe_id").loc[1010]
    network = Network(FIXTURE, corrections={"swap_nodes": {1010: "recorded backwards"}})
    after = network.full_data.set_index("pipe_id").loc[1010]
    assert (after["from_node"], after["to_node"]) == ("D2", "D1")
    assert (after["us_invert"], after["ds_invert"]) == (before["us_invert"], before["ds_invert"])
    assert after.geometry.equals(before.geometry)
    assert network.graph.has_edge("D2", "D1")
    assert not network.graph.has_edge("D1", "D2")
    row = network.corrections_report.iloc[0]
    assert row[["action", "old", "new"]].tolist() == ["swap_nodes", "D1 -> D2", "D2 -> D1"]


def test_pipe_ids_match_by_text():
    network = Network(FIXTURE, corrections={"drop": {"1009": "as text"}})
    assert network.corrections_report["pipe_id"].tolist() == [1009]


def test_unknown_pipes_named_together():
    with pytest.raises(ValueError, match="pipes not in the data: 1, 2"):
        Network(FIXTURE, corrections={"drop": {1: "x"}, "swap_nodes": {2: "y"}})


def test_unknown_column():
    with pytest.raises(ValueError, match="cannot be set: \\['colour'\\]"):
        Network(FIXTURE, corrections={"set": {1000: {"colour": "red"}}})
    with pytest.raises(ValueError, match="cannot be set"):
        PipeCorrections(set={1000: {"abandoned": False}})


def test_unknown_key():
    with pytest.raises(ValueError, match="Unknown keys"):
        Network(FIXTURE, corrections={"remove": {1000: "x"}})


@pytest.mark.parametrize(
    "corrections",
    [
        {"drop": {1000: "x"}, "set": {1000: {"status": "built"}}},
        {"drop": {1000: "x"}, "swap_nodes": {"1000": "y"}},
        {"swap_nodes": {1000: "x"}, "set": {1000: {"from_node": "A3"}}},
    ],
)
def test_conflicting_entries(corrections):
    with pytest.raises(ValueError, match="conflicting entries for pipes 1000"):
        Network(FIXTURE, corrections=corrections)


def test_same_pipe_twice():
    with pytest.raises(ValueError, match="more than once: 1000"):
        Network(FIXTURE, corrections={"drop": {1000: "x", "1000": "y"}})


def test_reason_must_be_text():
    with pytest.raises(ValueError, match="must be its reason"):
        PipeCorrections(drop={1000: True})  # type: ignore[dict-item]


def test_logged_with_reasons(caplog):
    with caplog.at_level(logging.INFO, logger="swmmbuilder"):
        Network(
            FIXTURE,
            corrections=PipeCorrections(
                drop={1009: "overflow tail"},
                set={1000: {"material": "PVC-2", "reason": "relined"}},
                swap_nodes={1010: "recorded backwards"},
            ),
        )
    assert "Pipe 1009 dropped by corrections: overflow tail" in caplog.text
    assert "Pipe 1000: material set from 'PVC' to 'PVC-2' by corrections: relined" in caplog.text
    assert "Pipe 1010: node IDs swapped by corrections" in caplog.text
    assert "recorded backwards" in caplog.text


def test_toml_file(tmp_path):
    path = tmp_path / "corrections.toml"
    path.write_text(
        '[drop]\n1009 = "overflow tail"\n\n'
        '[set.1000]\nstatus = "built"\nreason = "built since"\n\n'
        '[swap_nodes]\n1010 = "recorded backwards"\n'
    )
    network = Network(proposed(1000), corrections=path)
    assert network.corrections_report["action"].tolist() == ["drop", "set", "swap_nodes"]
    assert network.corrections_report["pipe_id"].tolist() == [1009, 1000, 1010]
    assert network.graph.has_edge("A1", "A2")
