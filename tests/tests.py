from sewernetworkx import Network, SubNetwork
from pyswmm import Simulation
import pytest


@pytest.fixture
def geojson_file_path():
    return f'../examples/geojson_test_data/test_data_small.geojson'

@pytest.fixture
def inp_file_path():
    return f'./examples/inp_test_file/test.inp'


@pytest.fixture
def swmm_simulation(inp_file_path):
    # Initialize the PySWMM model
    simulation = Simulation(inp_file_path)
    return simulation


@pytest.fixture
def level_2_model(geojson_file_path):
    levels = 2
    starting_node = 'WT1420'
    network = Network(geojson_file_path)
    subnetwork = SubNetwork(network, starting_node, levels)
    return subnetwork

@pytest.fixture
def create_pyswmm_file(level_2_model):
    level_2_model.write_pyswmm('pytest.inp')
    file_path = f'./pytest.inp'
    assert os.path.exists(file_path)
    yield file_path
    os.remove(file_path)

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

def test_pyswmm_write(create_pyswmm_file):
    assert os.path.exists(create_pyswmm_file)
