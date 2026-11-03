"""Numeric node IDs (#58): Kitchener's node and pipe IDs load as int, and every API that takes a
node ID must accept them, and their text, the same way."""

import pytest
from conftest import TESTS, read_sections, write_inp

from swmmbuilder import KITCHENER_ON, Network, Subnetwork

OUTLET = 2151735  # downstream node of the pipe leaving the Kitchener catchment


@pytest.fixture(scope="module")
def network():
    return Network(TESTS / "data" / "kitchener_catchment.geojson", schema=KITCHENER_ON)


@pytest.fixture(scope="module")
def catchment(network):
    return network.catchment(OUTLET)


def upstream_of(network, node):
    return sorted(network.graph.predecessors(node))[0]


def test_ids_are_int(network):
    assert all(isinstance(n, int) for n in network.graph)


@pytest.mark.parametrize("outlet", [OUTLET, [OUTLET]])
def test_catchment_with_int_outlet(network, outlet):
    sub = network.catchment(upstream_of(network, OUTLET), outlet=outlet)
    assert OUTLET in sub.nodes


def test_node_list_round_trip(catchment, tmp_path):
    path = tmp_path / "nodes.txt"
    catchment.to_node_list(path)
    again = Subnetwork.from_node_list(catchment.network, path)
    assert again.nodes == catchment.nodes
    assert all(isinstance(n, int) for n in again.nodes)


def test_dead_end_error_names_int_nodes(network, tmp_path):
    # The pipe just above the outlet, without the outlet: its downstream node is a dead end.
    up = upstream_of(network, OUTLET)
    sub = network.subnetwork([up, upstream_of(network, up)])
    with pytest.raises(ValueError, match=rf"dead ends .*{up}"):
        write_inp(sub, tmp_path, outfalls=())


def test_cycle_error_names_int_nodes(network, tmp_path):
    from swmmbuilder.subnetwork import TopologyReport
    from swmmbuilder.swmm import check_topology

    with pytest.raises(ValueError, match="cycle 1 -> 2 -> 1"):
        check_topology(TopologyReport(cycles=[[1, 2]]), "KINWAVE", "raise")


def test_options_match_int_ids_and_their_text(catchment, tmp_path):
    up = upstream_of(catchment.network, OUTLET)

    def written(key):
        convert = str if key == "str" else (lambda n: n)
        path = write_inp(
            catchment,
            tmp_path,
            name=f"{key}.inp",
            routing="DYNWAVE",
            outfalls={convert(OUTLET): "FREE"},
            dead_ends="outfall",
            dwf={convert(up): 1.5},
            wet_wells={convert(up): dict(max_depth=3.0, area=4.0)},
            height_correction={convert(up): 300.0},
        )
        return read_sections(path)

    by_int, by_str = written("int"), written("str")
    assert by_int == by_str
    assert [row[:3] for row in by_int["DWF"]] == [[str(up), "FLOW", "1.5"]]
    assert [row[:2] for row in by_int["STORAGE"]] == [[str(up), "300.0"]]
    assert str(OUTLET) in {row[0] for row in by_int["OUTFALLS"]}


def test_same_node_twice_with_different_values_is_refused(catchment):
    up = upstream_of(catchment.network, OUTLET)
    with pytest.raises(ValueError, match="dwf: nodes given more than once"):
        catchment.match_nodes({up: 1.0, str(up): 2.0}, "dwf")
    assert catchment.match_nodes({up: 1.0, str(up): 1.0}, "dwf") == {up: 1.0}
