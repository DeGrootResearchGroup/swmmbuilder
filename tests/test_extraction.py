import os
import subprocess
import sys

import networkx as nx
import pytest
from conftest import FIXTURE, REPO, quiet, read_sections, subnetwork_with_edges, write_inp

# Every node that drains to OUT (flow direction in tests/data/make_fixture.py). Not OV1/OV2 (the
# overflow leaves J2), S1/S2 (storm) or 'Abandoned' (abandoned by ID).
OUT_CATCHMENT = {
    "A1", "A2", "B1", "B2", "C1", "D1", "D2", "E1", "F1", "F2", "F3", "F4",
    "G1", "G2", "G3", "G4", "J1", "J2", "OUT", "X1",
}  # fmt: skip

PROBE = """
import contextlib, io, warnings
from sewernetworkx import Network, SubNetwork
warnings.simplefilter("ignore", DeprecationWarning)
with contextlib.redirect_stdout(io.StringIO()):
    n = Network({fixture!r})
    s = SubNetwork({{"network_data": n, "starting_node": "J1", "outfall_node": ["OUT"], "depth": 2,
                    "outfall_node_type": ["FREE"], "scale": 1.0}})
    {call}
print(",".join(sorted(s.subgraph.nodes)) + ";" + ",".join(map(str, s.subgraph.edges)))
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
        's.catchment("OUT")',
        's.catchment("J1", max_hops=2)',
        's.catchment("J1", max_length=100)',
        's.from_nodes(["OUT", "J2", "J1", "A2"])',
        "s.create_sub_network()",  # deprecated, but made deterministic
    ],
)
def test_extraction_is_deterministic_across_hash_seeds(call):
    results = {extract_with_seed(seed, call) for seed in range(8)}
    assert len(results) == 1, f"{len(results)} different results: {sorted(results)}"


@pytest.fixture
def sub(network):
    return subnetwork_with_edges(network, [])


def flow_edges(subgraph):
    """Subgraph edges as (upstream, downstream); the subgraph stores (downstream, upstream)."""
    return {(up, down) for down, up in subgraph.edges}


def test_graph_keeps_parallel_pipes(network):
    G = network.graph
    assert isinstance(G, nx.MultiDiGraph)
    assert sorted(G["J1"]["J2"]) == [1006, 1007]
    assert {G["J1"]["J2"][k]["pipe"] for k in (1006, 1007)} == {1006, 1007}
    assert "Abandoned" not in G and "S1" not in G


def test_catchment_of_outlet(sub):
    subgraph = sub.catchment("OUT")
    assert set(subgraph.nodes) == OUT_CATCHMENT
    assert flow_edges(subgraph) == {
        (up, down) for up, down in sub.graph.edges() if down in OUT_CATCHMENT
    }
    assert subgraph is sub.subgraph


def test_catchment_keeps_parallel_pipes_on_edge(sub):
    subgraph = sub.catchment("OUT")
    assert subgraph.edges["J2", "J1"]["pipes"] == [1006, 1007]
    assert subgraph.edges["J1", "A2"]["pipes"] == [1001]


def test_catchment_max_hops(sub):
    one = sub.catchment("J1", max_hops=1)
    assert set(one.nodes) == {"A2", "B2", "C1", "J1", "J2", "OUT"}
    assert flow_edges(one) == {
        ("A2", "J1"),
        ("B2", "J1"),
        ("C1", "J1"),
        ("J1", "J2"),
        ("J2", "OUT"),
    }
    two = sub.catchment("J1", max_hops=2)
    assert set(two.nodes) == set(one.nodes) | {"A1", "X1", "B1"}


def test_catchment_max_length(sub):
    # Upstream path lengths to J1: A2 60, B2 50, C1 40, X1 90, A1 140, B1 120.
    assert set(sub.catchment("J1", max_length=100).nodes) == {
        "A2", "B2", "C1", "X1", "J1", "J2", "OUT",
    }  # fmt: skip
    both = sub.catchment("J1", max_hops=1, max_length=55)
    assert set(both.nodes) == {"B2", "C1", "J1", "J2", "OUT"}


def test_catchment_path_to_every_listed_outlet(sub):
    subgraph = sub.catchment("J1", max_hops=0, outlet=["OUT", "OV2"])
    assert set(subgraph.nodes) == {"J1", "J2", "OUT", "OV1", "OV2"}


def test_catchment_rejects_unreachable_outlet_and_unknown_node(sub):
    with pytest.raises(ValueError, match="OV1"):
        sub.catchment("OV1", max_hops=1, outlet="OUT")
    with pytest.raises(ValueError, match="NOPE"):
        sub.catchment("NOPE")


def test_from_nodes_is_induced(sub):
    subgraph = sub.from_nodes(["OUT", "J2", "J1", "OV1"])
    assert set(subgraph.nodes) == {"J1", "J2", "OUT", "OV1"}
    assert flow_edges(subgraph) == {("J1", "J2"), ("J2", "OUT"), ("J2", "OV1")}
    assert subgraph.edges["J2", "J1"]["pipes"] == [1006, 1007]
    with pytest.raises(ValueError, match="Abandoned"):
        sub.from_nodes(["J1", "Abandoned"])


def test_node_list_round_trip(network, sub, tmp_path):
    original = sub.catchment("J1", max_hops=2)
    path = tmp_path / "nodes.txt"
    sub.to_node_list(path)
    nodes = [line for line in path.read_text().splitlines() if not line.startswith("#")]
    assert nodes == sorted(original.nodes)

    rebuilt = subnetwork_with_edges(network, []).from_node_list(path)
    assert list(rebuilt.nodes) == list(original.nodes)
    assert list(rebuilt.edges(data=True)) == list(original.edges(data=True))


def test_writer_consumes_catchment(sub, tmp_path):
    sub.catchment("J1", max_hops=1)
    sections = read_sections(write_inp(sub, tmp_path))
    nodes = {row[0] for row in sections["JUNCTIONS"] + sections["OUTFALLS"]}
    assert nodes == {"A2", "B2", "C1", "J1", "J2", "OUT"}
    conduits = {row[0]: (row[1], row[2]) for row in sections["CONDUITS"]}
    assert conduits["1006"] == conduits["1007"] == ("J1", "J2")  # parallel pipes both written
    assert set(conduits.values()) == flow_edges(sub.subgraph)


def test_deprecated_methods_warn(sub):
    with pytest.warns(DeprecationWarning, match="create_sub_network"):
        quiet(sub.create_sub_network)
    with pytest.warns(DeprecationWarning, match="explore_nodes"):
        quiet(sub.explore_nodes, sub._initialize_digraph(), "successors")
    with pytest.warns(DeprecationWarning, match="explore_branch"):
        quiet(
            sub.explore_branch,
            {"A2"},
            set(),
            sub._initialize_digraph(),
            "predecessors",
            "successors",
        )
