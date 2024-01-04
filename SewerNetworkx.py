import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import shutil
from pyswmm import Simulation, Nodes, Subcatchments

class Network:
    full_data = None

    def __init__(self, file_path, remove_storm=True):
        self.file_path = file_path
        data = gpd.read_file(self.file_path)
        Network.full_data = gpd.GeoDataFrame(data, geometry='geometry')

        if remove_storm:
            print('Removing storm sewers')
            key = "FlowType"
            condition = Network.full_data[key] != "STM"
            Network.full_data = Network.full_data[condition]

    # Needs to be updated
    def filter_data_gisid(self, network_data, ppty):
        # Check is ppty is a string
        if not isinstance(ppty, str):
            raise TypeError("Input must be a string.")

        return network_data[ppty].iloc[0]

    # Needs to be updated
    def extract_data(self, criteria):
        # Check if criteria is a string
        if not isinstance(criteria, int):
            raise TypeError("Input must be an integer.")

        # Filter the GeoDataFrame based on the given criteria
        filtered_geo_series = self.full_data[self.full_data["GIS_ID"] == criteria]
        return filtered_geo_series

class SubNetwork(Network):

    def __init__(self, node, depth=None, file_path=None):
        if file_path:
            super().__init__(file_path)
        self.network_data = {}
        self.start = node
        self.depth = depth
        self.subgraph = nx.DiGraph()

        print("Removing Unnecessary Entries")
        keys_to_remove = ["GIS_ID", "GIS_FeatureKey", "FlowType", "StreetName", "Location",
                          "Material", "MaterialClass", "Slope", "AverageDepth", "MaximumDepth",
                          "PDC_Count", "ConstructionDate", "Ownership", "Assumed", "Combined",
                          "Relief", "LastEditDate", "Abandoned", "Grouted", "CCTVInspectionDate",
                          "ShapeSTLength"]

        self.reduced_keys_data = Network.full_data.copy()

        for keys in keys_to_remove:

            if keys in self.reduced_keys_data.columns:
                self.reduced_keys_data.drop(keys, axis=1, inplace=True)

        print("Creating Pipes")

        for index, pipe in self.reduced_keys_data.iterrows():
            if pipe['UpstreamID'] is not None and pipe['DownstreamID'] is not None:
                self.network_data[pipe['OBJECTID']] = (pipe['DownstreamID'], pipe['UpstreamID'])

    def create_sub_network(self):
        print("Creating subnetwork for plotting")

        # Initialize Networkx graph
        G = nx.DiGraph()

        # Add edges from network_data
        G.add_edges_from(self.network_data.values())

        # Extract the connected subgraph starting from the specified node with limited depth
        if self.depth is not None:
            connected_subgraph = {self.start}
            current_depth = 0
            while current_depth < self.depth:
                # Create set to ensure no duplicates
                new_nodes = set()
                for nodes in list(connected_subgraph):

                    if len(connected_subgraph) == 1:
                        connected_subgraph.update(G.successors(nodes))
                    else:
                        next_node = set(G.successors(nodes))
                        # Check if next_node is already in connect_graph list
                        if all(not any(node in node_set for node in next_node) for node_set in connected_subgraph):
                            if len(next_node) > 1:
                                print(f'Branch Detected at node: {nodes}')
                                for branch in list(next_node):
                                    print(f'Exploring branch from node: {branch}')
                                    new_nodes.add(branch)
                                    # print(f'Successor of {branch} is: {list(G.successors(branch))}')
                                    branch_next = set(G.successors(branch))
                                    new_nodes.update(branch_next)
                                    while len(branch_next) == 1:
                                        next_node = next(iter(branch_next))
                                        branch_next = set(G.successors(next_node))
                                        # print(branch_next)
                                        if len(branch_next) >= 2:
                                            # print(f'New branch detected at {next_node}')
                                            break
                                        elif len(branch_next) == 0:
                                            # print(f'End point reached at {next_node}')
                                            new_nodes.update(branch_next)
                                            break
                                        new_nodes.update(branch_next)

                                current_depth += 1
                            elif len(next_node) == 1:
                                new_nodes.update(next_node)

                connected_subgraph.update(new_nodes)

        else:
            connected_subgraph = nx.descendants(G, self.start).union([self.start])

        # Create a subgraph from the connected nodes
        self.subgraph = G.subgraph(connected_subgraph)

        # Add colouring
        node_colour = []

        for node in self.subgraph.nodes:
            if node == self.start:
                node_colour.append('green')
            elif len(list(self.subgraph.successors(node))) > 1:
                node_colour.append('red')
            elif len(list(self.subgraph.successors(node))) == 0:
                node_colour.append('orange')
            else:
                node_colour.append('blue')

        return self.subgraph, node_colour

    def node_coordinates(self, network_nodes, scale, plant_node=False):
        pos = {}
        for index, row in self.reduced_keys_data.iterrows():
            upstream_id = row['UpstreamID']
            downstream_id = row['DownstreamID']
            if plant_node and upstream_id == self.start:
                # print('Plant node')
                coordinates = row['geometry'].coords[0]
                pos[downstream_id] = (coordinates[0] * scale, coordinates[1] * scale)
            elif upstream_id in network_nodes:
                # print('Found upstream node')
                coordinates = row['geometry'].coords[0]
                pos[upstream_id] = (coordinates[0] * scale, coordinates[1] * scale)
        return pos

    # Probably need to come up with better names for this... as they are used above
    def trim_subnetwork(self, network_nodes, trim_node):
        trimmed_graph = network_nodes.copy()
        for trimming in set(trim_node):
            successor_nodes_to_remove = list(trimmed_graph.successors(trimming))
            trimmed_graph.remove_nodes_from(successor_nodes_to_remove)

        return trimmed_graph

    def write_pyswmm(self, file_path, file_name):
        template_file = "pyswmm_input_file.inp"
        if not isinstance(file_path, str):
            raise TypeError("File path for .inp file must be a string.")
        output_file_path = file_path
        output_file_name = file_name
        shutil.copy(template_file, output_file_name)

        # Start a pyswmm simulation to input Node data
        with Simulation(output_file_name) as sim:





