import logging

import geopandas as gpd
from conftest import FIXTURE, REPO

from swmmbuilder import Network


def test_storm_pipes_removed(network):
    assert "storm" not in set(network.full_data["flow_type"])


def test_pipes_abandoned_by_id_removed(network):
    assert "Abandoned" not in network.graph
    assert network.graph.has_node("A1")


def test_pipes_flagged_abandoned_removed(network):
    assert not network.graph.has_edge("X1", "A2")


def test_networks_are_independent():
    first = Network(str(FIXTURE))
    n_first, nodes_first = len(first.full_data), set(first.graph)
    Network(str(REPO / "examples" / "geojson_test_data" / "test_data_small.geojson"))
    assert len(first.full_data) == n_first
    assert set(first.graph) == nodes_first


def test_flagged_abandoned_pipes_not_in_graph(network):
    assert "X1" not in network.graph
    assert network.removed == {
        "flow_type": {"storm": 1},
        "status": {"abandoned": 1},
        "structure_type": {},
        "invalid": {},
        "corrections": {},
    }


def test_keep_abandoned():
    network = Network(str(FIXTURE), keep_abandoned=True)
    assert network.graph.has_edge("X1", "A2")
    assert network.removed["status"] == {}


def test_flow_types_none_keeps_every_pipe():
    network = Network(str(FIXTURE), flow_types=None)
    assert network.graph.has_edge("S1", "S2")
    assert network.removed["flow_type"] == {}


def test_missing_flow_type_dropped_and_logged(tmp_path, caplog):
    data = gpd.read_file(FIXTURE)
    data.loc[data["UpstreamID"] == "E1", "FlowType"] = None
    path = tmp_path / "no_flow_type.geojson"
    data.to_file(path, driver="GeoJSON")
    with caplog.at_level(logging.INFO, logger="swmmbuilder"):
        network = Network(str(path))
    assert "E1" not in network.graph
    assert network.removed["flow_type"] == {"None": 1, "storm": 1}
    assert "Removed 1 pipes with flow type None" in caplog.text


def test_network_from_geodataframe():
    network = Network(gpd.read_file(FIXTURE))
    assert set(network.graph.edges(keys=True)) == set(Network(FIXTURE).graph.edges(keys=True))


def test_network_pipes_match_a_whole_subnetwork(network):
    # Network.pipes() is Subnetwork.pipes() for the whole network (#161).
    pipes = network.pipes()
    whole = network.subnetwork(network.nodes.index).pipes()
    assert list(pipes["pipe_id"]) == list(whole["pipe_id"])
    assert len(pipes) == network.graph.number_of_edges()
    assert list(pipes.columns) == list(network.full_data.columns)
