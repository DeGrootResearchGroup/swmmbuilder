import contextlib
import io

import geopandas as gpd
import pytest
from conftest import FIXTURE, REPO, quiet, subnetwork_with_edges

from sewernetworkx import Network


def pipe_pairs(sub):
    return {(up, down) for down, up in sub.network_data.values()}


def test_storm_pipes_removed(network):
    assert "STM" not in set(network.full_data["FlowType"])


def test_pipes_abandoned_by_id_removed(network):
    sub = subnetwork_with_edges(network, [])
    assert ("Abandoned", "A1") not in pipe_pairs(sub)


def test_pipes_flagged_abandoned_removed(network):
    sub = subnetwork_with_edges(network, [])
    assert ("X1", "A2") not in pipe_pairs(sub)


def test_networks_are_independent():
    first = quiet(Network, str(FIXTURE))
    n_first, nodes_first = len(first.full_data), set(first.graph)
    quiet(Network, str(REPO / "examples" / "geojson_test_data" / "test_data_small.geojson"))
    assert len(first.full_data) == n_first
    assert set(first.graph) == nodes_first


def test_flagged_abandoned_pipes_not_in_graph(network):
    assert "X1" not in network.graph
    assert network.removed == {"flow_type": {"STM": 1}, "abandoned": 1}


def test_keep_abandoned():
    network = quiet(Network, str(FIXTURE), keep_abandoned=True)
    assert network.graph.has_edge("X1", "A2")
    assert network.removed["abandoned"] == 0


def test_flow_types_none_keeps_every_pipe():
    network = quiet(Network, str(FIXTURE), flow_types=None)
    assert network.graph.has_edge("S1", "S2")
    assert network.removed["flow_type"] == {}


def test_missing_flow_type_dropped_and_logged(tmp_path):
    data = gpd.read_file(FIXTURE)
    data.loc[data["UpstreamID"] == "E1", "FlowType"] = None
    path = tmp_path / "no_flow_type.geojson"
    data.to_file(path, driver="GeoJSON")
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        network = Network(str(path))
    assert "E1" not in network.graph
    assert network.removed["flow_type"] == {"None": 1, "STM": 1}
    assert "Removed 1 pipes with FlowType None" in out.getvalue()


def test_remove_storm_deprecated():
    with pytest.warns(DeprecationWarning, match="flow_types"):
        network = quiet(Network, str(FIXTURE), remove_storm=False)
    assert network.graph.has_edge("S1", "S2")
