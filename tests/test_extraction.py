import os
import subprocess
import sys

import networkx as nx
import pytest
from conftest import FIXTURE, REPO, read_sections, write_inp

from sewernetworkx import Subnetwork

# Every node that drains to OUT (flow direction in tests/data/make_fixture.py). Not OV1/OV2 (the
# overflow leaves J2), S1/S2 (storm), 'Abandoned' (abandoned by ID) or X1 (Abandoned flag).
OUT_CATCHMENT = {
    "A1", "A2", "B1", "B2", "C1", "D1", "D2", "E1", "F1", "F2", "F3", "F4",
    "G1", "G2", "G3", "G4", "J1", "J2", "OUT",
}  # fmt: skip

PROBE = """
from sewernetworkx import Network
n = Network({fixture!r})
s = {call}
print(",".join(s.nodes) + ";" + ",".join(map(str, s.graph.edges(keys=True))))
"""


def extract_with_seed(seed, call):
    env = {**os.environ, "PYTHONHASHSEED": str(seed), "PYTHONPATH": str(REPO)}
    out = subprocess.run(
        [sys.executable, "-W", "ignore", "-c", PROBE.format(fixture=str(FIXTURE), call=call)],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return out.stdout.strip().splitlines()[-1]


@pytest.mark.parametrize(
    "call",
    [
        'n.catchment("OUT")',
        'n.catchment("J1", max_hops=2, outlet="OUT")',
        'n.catchment("J1", max_length=100, outlet="OUT")',
        'n.subnetwork(["OUT", "J2", "J1", "A2"])',
    ],
)
def test_extraction_is_deterministic_across_hash_seeds(call):
    results = {extract_with_seed(seed, call) for seed in range(8)}
    assert len(results) == 1, f"{len(results)} different results: {sorted(results)}"


def flow_edges(sub):
    return set(sub.graph.edges())


def test_graph_keeps_parallel_pipes(network):
    G = network.graph
    assert isinstance(G, nx.MultiDiGraph)
    assert sorted(G["J1"]["J2"]) == [1006, 1007]
    assert {G["J1"]["J2"][k]["pipe"] for k in (1006, 1007)} == {1006, 1007}
    assert "Abandoned" not in G and "S1" not in G


def test_catchment_of_outlet(network):
    sub = network.catchment("OUT")
    assert set(sub.nodes) == OUT_CATCHMENT
    assert flow_edges(sub) == {
        (up, down) for up, down in network.graph.edges() if down in OUT_CATCHMENT
    }


def test_catchment_keeps_parallel_pipes(network):
    sub = network.catchment("OUT")
    assert sorted(sub.graph["J1"]["J2"]) == [1006, 1007]
    assert sorted(sub.graph["A2"]["J1"]) == [1001]


def test_catchment_max_hops(network):
    one = network.catchment("J1", max_hops=1)
    assert set(one.nodes) == {"A2", "B2", "C1", "J1"}
    one = network.catchment("J1", max_hops=1, outlet="OUT")
    assert set(one.nodes) == {"A2", "B2", "C1", "J1", "J2", "OUT"}
    assert flow_edges(one) == {
        ("A2", "J1"),
        ("B2", "J1"),
        ("C1", "J1"),
        ("J1", "J2"),
        ("J2", "OUT"),
    }
    two = network.catchment("J1", max_hops=2, outlet="OUT")
    assert set(two.nodes) == set(one.nodes) | {"A1", "B1"}


def test_catchment_max_length(network):
    # Upstream path lengths to J1: A2 60, B2 50, C1 40, A1 140, B1 120.
    assert set(network.catchment("J1", max_length=100, outlet="OUT").nodes) == {
        "A2", "B2", "C1", "J1", "J2", "OUT",
    }  # fmt: skip
    both = network.catchment("J1", max_hops=1, max_length=55, outlet="OUT")
    assert set(both.nodes) == {"B2", "C1", "J1", "J2", "OUT"}


def test_catchment_path_to_every_listed_outlet(network):
    sub = network.catchment("J1", max_hops=0, outlet=["OUT", "OV2"])
    assert set(sub.nodes) == {"J1", "J2", "OUT", "OV1", "OV2"}


def test_catchment_rejects_unreachable_outlet_and_unknown_node(network):
    with pytest.raises(ValueError, match="OV1"):
        network.catchment("OV1", max_hops=1, outlet="OUT")
    with pytest.raises(ValueError, match="NOPE"):
        network.catchment("NOPE")


def test_subnetwork_is_induced(network):
    sub = network.subnetwork(["OUT", "J2", "J1", "OV1"])
    assert sub.nodes == ["J1", "J2", "OUT", "OV1"]
    assert flow_edges(sub) == {("J1", "J2"), ("J2", "OUT"), ("J2", "OV1")}
    assert sorted(sub.graph["J1"]["J2"]) == [1006, 1007]
    with pytest.raises(ValueError, match="Abandoned"):
        network.subnetwork(["J1", "Abandoned"])


def test_node_list_round_trip(network, tmp_path):
    original = network.catchment("J1", max_hops=2, outlet="OUT")
    path = tmp_path / "nodes.txt"
    original.to_node_list(path)
    nodes = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    assert nodes == original.nodes

    rebuilt = Subnetwork.from_node_list(network, path)
    assert rebuilt.nodes == original.nodes
    assert list(rebuilt.graph.edges(keys=True, data=True)) == list(
        original.graph.edges(keys=True, data=True)
    )


def test_writer_consumes_catchment(network, tmp_path):
    sub = network.catchment("J1", max_hops=1, outlet="OUT")
    sections = read_sections(write_inp(sub, tmp_path, routing="DYNWAVE"))
    nodes = {row[0] for row in sections["JUNCTIONS"] + sections["OUTFALLS"]}
    assert nodes == {"A2", "B2", "C1", "J1", "J2", "OUT"}
    conduits = {row[0]: (row[1], row[2]) for row in sections["CONDUITS"]}
    assert conduits["1006"] == conduits["1007"] == ("J1", "J2")  # parallel pipes both written
    assert set(conduits.values()) == flow_edges(sub)
