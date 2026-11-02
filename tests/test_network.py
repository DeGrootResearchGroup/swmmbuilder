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


@pytest.mark.xfail(strict=True, reason="#6: the Abandoned flag column is ignored")
def test_pipes_flagged_abandoned_removed(network):
    sub = subnetwork_with_edges(network, [])
    assert ("X1", "A2") not in pipe_pairs(sub)


def test_networks_are_independent():
    first = quiet(Network, str(FIXTURE))
    n_first, nodes_first = len(first.full_data), set(first.graph)
    quiet(Network, str(REPO / "examples" / "geojson_test_data" / "test_data_small.geojson"))
    assert len(first.full_data) == n_first
    assert set(first.graph) == nodes_first
