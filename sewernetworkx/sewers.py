"""Module for manipulating GeoJSON sewer data.

The module takes full GeoJSON data sets of a city sewer network,
generates a subnetwork given an initial starting node and number
of levels, and generates a .inp file to be used by pyswmm.
"""

import geopandas as gpd
import networkx as nx
import matplotlib.pyplot as plt
import shutil
from pyswmm import Simulation, Nodes, Subcatchments
import pkg_resources
import os


class Network:
    """Class used to represent an entire cities sewer network

        Parameters
        ----------
        file_path : str, Path
            Name of GeoJSON file

        Notes
        -----
        remove_store is used as a flag to automatically remove any sewer pipe
        flagged as a storm drain.
        """
    full_data = None

    def __init__(self, file_path, remove_storm=True):
        """Class initialization"""

        self.file_path = file_path
        data = gpd.read_file(self.file_path)
        Network.full_data = gpd.GeoDataFrame(data, geometry='geometry')

        if remove_storm:
            print('Removing storm sewers')
            key = "FlowType"
            condition = Network.full_data[key] != "STM"
            Network.full_data = Network.full_data[condition]

    """
    # Needs to be updated
    def filter_data_gisid(self, network_data, ppty):
        Filter Network class data.

                Parameters
                ----------
                network_data : class '__main__.Network'
                    Sewer network class

                ppty : str
                    Property to filter data by

                Returns
                -------
                str
                    The GIS_ID of the filtered data

        if not isinstance(ppty, str):
            raise TypeError("Input must be a string.")

        return network_data[ppty].iloc[0]
    """

    def get_upstream_invert(self, upstream_id):
        """Find the UpstreamInvert value based on the UpstreamID from a GeoPandas DataFrame.

            Parameters
            ----------
            upstream_id : str
                The UpstreamID to search for.

            Returns
            ----------
            float or None
                The corresponding UpstreamInvert value if found, else None.
            """
        if not isinstance(upstream_id, str):
            raise TypeError("Input must be a string")

        row = self.full_data[self.full_data['UpstreamID'] == upstream_id]
        if not row.empty:
            return row.iloc[0]['UpstreamInvert']
        else:
            return None  # Handle the case where UpstreamID is not found


class SubNetwork(Network):
    """Class used to represent an entire cities sewer network

            Parameters
            ----------
            node : str
                Starting node of the subnetwork.

            depth : int
                Number of branches to be contained within subnetwork.

            file_path : str, Path
                Path to GeoJSON file.
            """

    def __init__(self, subnetwork_data, file_path=None):
        if file_path:
            super().__init__(file_path)

        required_keys = ['network_data', 'starting_node', 'outfall_node', 'outfall_node_type', 'depth', 'scale']
        # Ensure all keys are being specified
        missing_keys = [key for key in required_keys if key not in subnetwork_data]
        if missing_keys:
            raise ValueError(f"Missing the following keys: {', '.join(missing_keys)}")

        self.network_data = {}
        self.start = subnetwork_data['starting_node']
        self.depth = subnetwork_data['depth']
        self.outfall_node = subnetwork_data['outfall_node']
        self.outfall_node_type = subnetwork_data['outfall_node_type']
        self.scale = subnetwork_data['scale']
        self.subgraph = nx.DiGraph()
        self.node_colour = []
        self.coords = {}

        if len(self.outfall_node) != len(self.outfall_node_type):
            raise ValueError(f'Length of {self.outfall_node} and {self.outfall_node_type} are not the same.')

        print("Removing Unnecessary Entries")
        keys_to_remove = ["GIS_ID", "GIS_FeatureKey", "FlowType", "StreetName", "Location",
                          "Material", "MaterialClass", "Slope", "AverageDepth", "MaximumDepth",
                          "PDC_Count", "ConstructionDate", "Ownership", "Assumed", "Combined",
                          "Relief", "LastEditDate", "Abandoned", "Grouted", "CCTVInspectionDate",
                          "ShapeSTLength"]

        self.reduced_keys_data = full_network.full_data.copy()

        for keys in keys_to_remove:
            if keys in self.reduced_keys_data.columns:
                self.reduced_keys_data.drop(keys, axis=1, inplace=True)

        print("Creating Pipes")
        for index, pipe in self.reduced_keys_data.iterrows():
            if pipe['UpstreamID'] is not None and pipe['DownstreamID'] is not None:
                if pipe['UpstreamID'] != 'Abandoned' and pipe['DownstreamID'] != 'Abandoned':
                    self.network_data[pipe['OBJECTID']] = (pipe['DownstreamID'], pipe['UpstreamID'])

    def create_sub_network(self):
        """Create subnetwork based on starting node and depth.

        Returns
        ----------
        subgraph : networkx.DiGraph
            Graph of the subnetwork

        node_colour : str
            Colour of each node if subgraph
        """

        print("Creating subnetwork for plotting")

        # Initialize Networkx graph
        G = nx.DiGraph()

        # Add edges from network_data
        G.add_edges_from(self.network_data.values())

        # Explore nodes for graph
        print('Checking successor nodes')
        successors_subgraph = self.explore_nodes(G, direction='successors')

        print('Checking predecessor nodes')
        predecessors_subgraph = self.explore_nodes(G, direction='predecessors')

        if ignore_successors:
            connected_subgraph = predecessors_subgraph
        elif ignore_predecessors:
            connected_subgraph = successors_subgraph
        else:
            connected_subgraph = successors_subgraph | predecessors_subgraph

        # Create a subgraph from the connected nodes
        self.subgraph = G.subgraph(connected_subgraph)
        self._colouring()
        self._coordinates()

    def _colouring(self):

        for node in self.subgraph.nodes:
            if node == self.start or node == self.outfall_node:
                self.node_colour.append('green')
            elif len(list(self.subgraph.successors(node))) > 1:
                self.node_colour.append('red')
            elif len(list(self.subgraph.successors(node))) == 0:
                self.node_colour.append('orange')
            else:
                self.node_colour.append('blue')

    def _coordinates(self):
        for index, row in self.reduced_keys_data.iterrows():
            upstream_id = row['UpstreamID']
            downstream_id = row['DownstreamID']
            if upstream_id == self.start:
                # print('Plant node')
                coordinates = row['geometry'].coords[0]
                self.coords[downstream_id] = (coordinates[0] * self.scale, coordinates[1] * self.scale)
            elif upstream_id in self.subgraph.nodes():
                # print('Found upstream node')
                coordinates = row['geometry'].coords[0]
                self.coords[upstream_id] = (coordinates[0] * self.scale, coordinates[1] * self.scale)

        # Check if each node in network_nodes has a position in pos
        missing_nodes = set(self.subgraph.nodes()) - set(self.coords.keys())
        if missing_nodes:
            print(
                f"Warning: Nodes {missing_nodes} do not have positions. Attempting to find positions for missing nodes.")

            for index, row in self.reduced_keys_data.iterrows():
                downstream_id = row['DownstreamID']

                if downstream_id in missing_nodes:
                    coordinates = row['geometry'].coords[1]
                    self.coords[downstream_id] = (coordinates[0] * self.scale, coordinates[1] * self.scale)
                    missing_nodes.remove(downstream_id)

        if missing_nodes:
            print(f"Error: Positions for nodes {missing_nodes} could not be found.") 

    def explore_nodes(self, G, direction):
        connected_subgraph = {self.start}
        current_depth = 0
        opposite_direction = 'successors' if direction == 'predecessors' else 'predecessors'

        while current_depth < self.depth:
            new_nodes = set()
            for node in list(connected_subgraph):
                # Find the next node in the sewer network
                next_nodes = set(getattr(G, direction)(node))
                if all(not
                       any(n in node_set for n in next_nodes)
                       for node_set in connected_subgraph
                       ):
                    if len(next_nodes) > 1:
                        print(f'Branch Detected at node: {node}')
                        branch_nodes = self.explore_branch(next_nodes, new_nodes, G, direction,
                                                           opposite_direction)
                        new_nodes.update(branch_nodes)
                        current_depth += 1
                    elif len(next_nodes) == 1:
                        new_nodes.update(next_nodes)
                        opposite_nodes = set(getattr(G, opposite_direction)(next(iter(next_nodes))))
                        print(f'Checking {opposite_direction.capitalize()} of {next_nodes}: {opposite_nodes}')
                        if next(iter(opposite_nodes)) == node and len(opposite_nodes) == 1:
                            break
                        else:

                            print(f'Possible branch detected at: {next_nodes}')
                            if next_nodes in connected_subgraph:
                                break
                            else:
                                current_depth += 1
                                inner_depth = current_depth
                                while inner_depth < self.depth:
                                    branch_nodes = self.explore_branch(opposite_nodes, new_nodes, G, opposite_direction,
                                                                       direction)
                                    new_nodes.update(branch_nodes)
                                    inner_depth += 1

                    elif len(next_nodes) == 0:
                        print('Endpoint reached')
                        new_nodes.update(next_nodes)
                        current_depth += 1
                        break

            connected_subgraph.update(new_nodes)

        return connected_subgraph

    @staticmethod    
    def explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction):
        for branch in list(branch_nodes):
            print(f'Exploring branch from node: {branch}')
            new_nodes.add(branch)
            # print(f'{direction.capitalize()} of {branch} is: {list(getattr(G, direction)(branch))}')
            branch_next = set(getattr(G, direction)(branch))

            new_nodes.update(branch_next)
            while len(branch_next) == 1:
                next_node = next(iter(branch_next))
                branch_next = set(getattr(G, direction)(next_node))
                if len(branch_next) >= 2:
                    print(f'New branch detected at {next_node}')
                    break
                elif len(branch_next) == 0:
                    print(f'Endpoint reached at {next_node}')
                    new_nodes.add(next_node)
                    new_nodes.update(branch_next)

                    # Check for successors or predecessors based on the current direction
                    opposite_nodes = set(getattr(G, opposite_direction)(next_node))
                    print(f'Checking {opposite_direction.capitalize()} of {next_node}: '
                          f'{opposite_nodes}')
                    break
                elif branch_next == 'Abandoned':
                    print('Reached abandoned node.')
                    break
                else:
                    new_nodes.update(branch_next)
        return new_nodes

    # Probably need to come up with better names for this... as they are used above
    def trim_subnetwork(self, network_nodes, trim_node, direction, copy = True):
        if direction not in ['successors', 'predecessors']:
            raise ValueError("Invalid direction. Use 'successors' or 'predecessors'.")

        if len(trim_nodes) == 0:
            raise ValueError("No nodes listed to trim. At least 1 node needs to be provided.")

        if copy:
            trimmed_graph = self.subgraph.copy()
            for nodes in trim_nodes:
                remove_nodes = set(getattr(trimmed_graph, direction)(nodes))
                remove_nodes.update(trim_nodes)
                trimmed_graph.remove_nodes_from(remove_nodes)
            return trimmed_graph
        else:
            for nodes in trim_nodes:
                remove_nodes = set(getattr(self.subgraph, direction)(nodes))
                remove_nodes.update(trim_nodes)
                self.subgraph.remove_nodes_from(remove_nodes)

    def extend_network(self):
        return self.subgraph

    def write_nodes_pyswmm(self, pyswmm_file, max_depth, init_depth, surcharge_depth, ponding_area):
        """ Add subnetwork data into the Nodes section of a .inp file.

        Parameters
        ----------

        output_file_path : str
            Path to the write out PySWMM input file (.inp).

        max_depth : float)
            Maximum depth for node (default is 0.0).

        init_depth : float
            Initial depth for node (default is 0.0).

        surcharge_depth : float
            Surcharge depth for node (default is 0.0).

        ponding_area : float
            Ponding area for node (default is 0.0).

        """

        with open(pyswmm_file, 'r') as file:
            content = file.read()

        nodes_content = ""
        outfall_content = ""
        # Iterate through the nodes in the graph and build the content
        for node in self.subgraph.nodes:
            if node in self.outfall_node:
                upstream_invert = self.get_upstream_invert(node)
                outfall_content += f"{node:<16} {upstream_invert:<11}  {self.outfall_node_type:<28} {'NO':<10}\n"
            else:
                upstream_invert = self.get_upstream_invert(node)
                nodes_content += (f"{node:<16} {upstream_invert:<11} {max_depth:<11} {init_depth:<11}"
                                  f"{surcharge_depth:<11} {ponding_area:<11}\n")

        # Replace {nodes} placeholder with the built nodes content
        modified_content = content.replace('{nodes}', nodes_content).replace('{outfalls}', outfall_content)

        # Write the modified content back to the file
        with open(pyswmm_file, 'w') as file:
            file.write(modified_content)

    def write_conduits_pyswmm(self, pyswmm_file, manning_n, inlet_offset, outlet_offset, init_flow, max_flow, barrels):
        """Add subnetwork data into the Conduits and XSection section of a .inp file.

        Parameters
        ----------

        output_file_path : str
            Path to the write out PySWMM input file (.inp).
        
        manning_n : float
            Manning N coefficient for conduit (default is 0.01).

        inlet_offset : float
            Inlet offset height for conduit (default is 0.0).

        outlet_offset : float
            Outlet offset height for conduit (default is 0.0).

        init_flow : float
            Initial flow value for conduit (default is 0.0).

        max_flow : float
            Maximum flow value for conduit (default is 0.0).

        barrels : float
            Barrel value for conduit (default is 0.0).

        """

        with open(pyswmm_file, 'r') as file:
            content = file.read()

        conduits_content = ""
        xsection_content = ""

        # Get information from nx.digraph to fill in data
        for edges in self.subgraph.edges:
            start_node, end_node = edges

            # Find the corresponding conduit in the GeoPandas DataFrame
            # note that given the order of the nx.digraph, it can be either upstream or downstream node
            geo_pipe_up = self.full_data[self.full_data['UpstreamID'] == start_node]
            geo_pipe_down = self.full_data[self.full_data['DownstreamID'] == start_node]

            if geo_pipe_up.iloc[0]["DownstreamID"] == end_node:

                name = geo_pipe_up.iloc[0]["OBJECTID"]
                outlet = geo_pipe_up.iloc[0]["DownstreamID"]
                length = geo_pipe_up.iloc[0]["Length"]

                conduits_content += (f"{name:<16} {start_node:<16} {outlet:<16}"
                                     f"{length:<11} {manning_n:<11} "
                                     f"{inlet_offset:<11} {outlet_offset:<11} "
                                     f"{init_flow:<11} {max_flow:<11}\n")

            elif geo_pipe_down.iloc[0]["UpstreamID"] == end_node:

                name = geo_pipe_down.iloc[0]["OBJECTID"]
                inlet = geo_pipe_down.iloc[0]["UpstreamID"]
                length = geo_pipe_down.iloc[0]["Length"]

                if geo_pipe_down.iloc[0]["PipeShape"] == 'R':
                    shape = 'RECT_CLOSED'
                    height = geo_pipe_down.iloc[0]["Height"] / 1000
                    width = geo_pipe_down.iloc[0]["Width"] / 1000
                elif geo_pipe_down.iloc[0]["PipeShape"] == 'Round':
                    shape = 'CIRCULAR'
                    height = geo_pipe_down.iloc[0]["Height"] / 1000
                    width = '0.0'

                conduits_content += (f"{name:<16} {inlet:<16} {start_node:<16}"
                                     f"{length:<11} {manning_n:<11} "
                                     f"{inlet_offset:<11} {outlet_offset:<11}"
                                     f"{init_flow:<11} {max_flow:<11}\n")

                xsection_content += (f"{name:<16} {shape:<16} {height:<16}"
                                     f"{width:<11} {0.0:<11} {0.0:<11} {barrels:<11}\n")

        # Replace variable in template file with conduit information
        modified_content = content.replace('{conduits}', conduits_content).replace('{xsections}', xsection_content)

        with open(pyswmm_file, 'w') as file:
            file.write(modified_content)

    def write_pyswmm(self, output_file=None,
                     max_depth=0.0, init_depth=0.0, surcharge_depth=0.0, ponding_area=0.0,
                     manning_n=0.01, inlet_offset=0.0, outlet_offset=0.0,
                     init_flow=0.0, max_flow=0.0, barrels=1.0, template_path=None, output_folder=None):
        """Generate a PySWMM input file based on a networkx digraph and GeoPandas DataFrame.

            Parameters
            ---------_

            output_file : str
                Name of new pyswmm input file (.inp).

            max_depth : float
                Maximum depth for node (default is 0.0).

            init_depth : float
                Initial depth for node (default is 0.0).

            surcharge_depth : float
                Surcharge depth for node (default is 0.0).

            ponding_area : float
                Ponding area for node (default is 0.0).

            manning_n : float
                Manning N coefficient for conduit (default is 0.01).

            inlet_offset : float
                Inlet offset height for conduit (default is 0.0).

            outlet_offset : float
                Outlet offset height for conduit (default is 0.0).

            init_flow : float
                Initial flow value for conduit (default is 0.0).

            max_flow : float
                Maximum flow value for conduit (default is 0.0).

            barrels : float
                Barrel value for conduit (default is 0.0).

            template_path : str
                Path to the pyswmm model template file (default is "template.inp").

            output_folder : str
                Path to were the new .inp file will be written.
            """

        if output_file is None:
            raise ValueError("The 'output_path' must be specified.")

        if template_path is None:
            # Get the package directory and construct the path to the template file
            package_dir = pkg_resources.resource_filename(__name__, '')
            template_path = os.path.join(package_dir, "examples", "pyswmm_template", "template.inp")

        if not os.path.exists(template_path):
            raise FileNotFoundError(f"The template file does not exist at the specified path: {template_path}")

        # Always create a copy of the template file
        if output_folder is not None:
            output_file_path = os.path.join(output_folder, output_file)
        else:
            # If output_folder is None, create a copy in the current working directory
            output_file_path = os.path.join(os.getcwd(), output_file)

        shutil.copy(template_path, output_file_path)
        print(f"Template file copied to: {output_file_path}")

        # Write nodes to pyswmm .inp file
        self.write_nodes_pyswmm(output_file_path, max_depth, init_depth,
                                surcharge_depth, ponding_area)

        # Write conduits to pyswmm inp file
        self.write_conduits_pyswmm(output_file_path, manning_n, inlet_offset,
                                   outlet_offset, init_flow, max_flow, barrels)
