import SewerNetworkx
from pyswmm import Simulation
import pytest


@pytest.fixture
def geojson_file_path():
    return f'./Testing_file.geojson'


@pytest.fixture
def swmm_simulation():
    # Initialize the PySWMM model
    simulation = Simulation('London_branch1.inp')
    return simulation


@pytest.fixture
def level_2_model(geojson_file_path):
    levels = 2
    starting_node = 'WT1420'
    network = Network(geojson_file_path)
    subnetwork = SubNetwork(network, starting_node, levels)
    return subnetwork


# Test if a Network class will be initialized
def test_init_network_with_geojson(geojson_file_path):
    network_test = Network(geojson_file_path)
    assert True


# Test if a SubNetwork class will be initialized
def test_init_subnetwork_with_geojson(level_2_model):
    assert True


def test_num_sewers(level_2_model):
    subnet_graph, node_colours = level_2_model.create_sub_network()
    assert (subnet_graph.number_of_nodes() == 53)


def test_end_nodes(level_2_model):
    subnet_graph, node_colours = level_2_model.create_sub_network()
    num = 0
    for nodes in subnet_graph.nodes:
        if len(list(subnet_graph.successors(nodes))) == 0:
            num += 1
    assert (num == 6)


def test_pyswmm_model(swmm_simulation):
    try:
        swmm_simulation.execute()
    except Exception as e:
        pytest.fail(f"pyswmm model execution failed with error: {e}")
