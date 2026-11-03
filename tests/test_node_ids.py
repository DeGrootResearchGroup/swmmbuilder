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


# --- Pipe IDs (#78): one matching rule for every pipe-keyed option --------------------------


# A station with a pump curve and a wet well (PS also has an overflow, so no ideal pump).
STATION = dict(
    pump=dict(curve=[(0.0, 10.0)], curve_type="PUMP2", startup=1.0, shutoff=0.2),
    wet_well=dict(max_depth=4.0, area=2.0),
)


@pytest.fixture(scope="module")
def stations():
    from test_pump_stations import MAIN, make_network

    return make_network().subnetwork(MAIN)  # int pipe IDs 1-5


def test_match_pipe(stations):
    assert stations.match_pipe(2) == 2
    assert stations.match_pipe("2") == 2
    assert stations.match_pipe("99") == "99"  # unknown: unchanged, for the caller to report
    assert stations.match_pipes({"1": "a", 2: "b"}, "opt") == {1: "a", 2: "b"}
    with pytest.raises(ValueError, match="opt: pipes not in the subnetwork: 98, 99"):
        stations.match_pipes({99: 1, "98": 1, 1: 1}, "opt")
    assert stations.match_pipes({99: 1, 1: 1}, "opt", unknown="ignore") == {1: 1}
    with pytest.raises(ValueError, match="opt: pipes given more than once: 1"):
        stations.match_pipes({1: 1.0, "1": 2.0}, "opt")


def _xsections_and_pumps(sub, convert, tmp_path):
    from test_pump_stations import C

    path = write_inp(
        sub,
        tmp_path,
        name=f"{convert.__name__}.inp",
        routing="DYNWAVE",
        outfalls=("OUT", "OVO"),
        pump_stations={"PS": dict(STATION, force_mains=[convert(2)])},
        force_main_roughness=C,
        force_main_roughness_by_pipe={convert(3): 140.0},
        size_overrides={convert(4): 0.45},
    )
    sections = read_sections(path)
    return sections["XSECTIONS"], sections["PUMPS"]


def test_pipe_keyed_options_match_int_ids_and_their_text(stations, tmp_path):
    def as_text(pipe):
        return str(pipe)

    def as_int(pipe):
        return pipe

    by_int = _xsections_and_pumps(stations, as_int, tmp_path)
    assert by_int == _xsections_and_pumps(stations, as_text, tmp_path)
    xsections = {r[0]: r for r in by_int[0]}
    assert xsections["3"][1:4] == ["FORCE_MAIN", "0.15", "140.0"]  # by_pipe override
    assert xsections["4"][2] == "0.45"  # size override


@pytest.mark.parametrize("option", ["pumps", "size_overrides"])
def test_unknown_pipe_is_reported_alike(stations, option):
    from swmmbuilder.swmm import build_model

    value = dict(curve="*") if option == "pumps" else 0.3
    with pytest.raises(ValueError, match=f"{option}: pipes not in the subnetwork: 99"):
        build_model(stations, **{option: {99: value}})


def test_roughness_by_pipe_ignores_other_pipes(stations, tmp_path):
    # One mapping may cover the force mains of several subnetworks.
    from test_pump_stations import C

    path = write_inp(
        stations,
        tmp_path,
        name="other.inp",
        routing="DYNWAVE",
        outfalls=("OUT", "OVO"),
        pump_stations={"PS": STATION},
        force_main_roughness=C,
        force_main_roughness_by_pipe={3: 140.0, 999: 1.0},
    )
    assert {r[0]: r for r in read_sections(path)["XSECTIONS"]}["3"][3] == "140.0"


def test_network_accessors(stations):
    network = stations.network
    assert {"U1", "PS", 1, 16} <= network.ids()
    # PS: pipes 1 in at 250.0, 2 out at 249.8, 5 out at 250.8.
    assert network.lowest_inverts({"PS"}) == {"PS": pytest.approx(249.8)}
    with pytest.raises(ValueError, match="Nodes not in the network: NOPE"):
        network.check_nodes(["PS", "NOPE"])
    assert network._check_nodes == network.check_nodes  # the former name still works
