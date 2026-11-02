"""Module for manipulating GeoJSON sewer data.

The module takes full GeoJSON data sets of a city sewer network,
generates a subnetwork given an initial starting node and number
of levels, and generates a .inp file to be used by pyswmm.
"""

import os
import re
from importlib import resources

import geopandas as gpd
import gmplot
import networkx as nx

# SWMM's SI flow units. With any of these, lengths, elevations and section sizes are in metres,
# matching the GIS data; with a US unit (CFS, GPM, MGD) SWMM would read them as feet.
SI_FLOW_UNITS = ("CMS", "LPS", "MLD")

# Placeholders that write_pyswmm fills in the template.
TEMPLATE_PLACEHOLDERS = ("nodes", "outfalls", "conduits", "xsections", "coordinates", "vertices")


def default_template_path():
    """Path-like (Traversable) to the SI template shipped with the package."""
    return resources.files("sewernetworkx") / "templates" / "template.inp"


def check_template(content):
    """Check that a .inp template is one the writer can fill.

    Parameters
    ----------
    content : str
        Text of the .inp template.

    Raises
    ----------
    ValueError
        If FLOW_UNITS is missing or not an SI unit (the GIS data are metric), or if the template
        contains a placeholder the writer does not fill (e.g. ``{pumps}``: pumps are not supported).
    """
    match = re.search(r"^\s*FLOW_UNITS\s+(\S+)", content, re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise ValueError(
            "The template has no FLOW_UNITS option (SWMM would default to CFS, US units). "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )
    flow_units = match.group(1).upper()
    if flow_units not in SI_FLOW_UNITS:
        raise ValueError(
            f"The template sets FLOW_UNITS {flow_units}, a US unit, so SWMM would read the metric "
            "GIS lengths, elevations and sizes as feet. "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )

    unknown = sorted(set(re.findall(r"\{(\w+)\}", content)) - set(TEMPLATE_PLACEHOLDERS))
    if unknown:
        raise ValueError(
            f"The template contains placeholders that are not filled: {unknown}. "
            f"Supported placeholders: {list(TEMPLATE_PLACEHOLDERS)}."
        )


def line_coords(geometry):
    """(x, y) vertices of a pipe geometry, in order from upstream to downstream.

    The parts of a MultiLineString are joined in order. Any z values are dropped.
    """
    parts = geometry.geoms if geometry.geom_type == "MultiLineString" else [geometry]
    return [(c[0], c[1]) for part in parts for c in part.coords]


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
        Network.full_data = gpd.GeoDataFrame(data, geometry="geometry")

        if remove_storm:
            print("Removing storm sewers")
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

        # row = self.full_data[self.full_data['UpstreamID'] == upstream_id]
        row = self.full_data[self.full_data["UpstreamID"] == upstream_id]
        if not row.empty:
            return row.iloc[0]["UpstreamInvert"]
        else:
            row = self.full_data[self.full_data["DownstreamID"] == upstream_id]
            if not row.empty:
                return row.iloc[0]["DownstreamInvert"]
            else:
                raise TypeError("Cant find either upstreaminvert of downstreaminvert")


class SubNetwork(Network):
    """Class used to represent an entire cities sewer network

    Parameters
    ----------
    subnetwork_data : dict
        dict containing all required information to generate subnetwork.
        Required entries:

        network_data : geoDataFrame
        starting_node : str
        outfall_node : list, str
        outfall_node_type : list, str
        depth : int
        scale : int

    file_path : str, Path
        Path to GeoJSON file.
    """

    def __init__(self, subnetwork_data, file_path=None):
        if file_path:
            super().__init__(file_path)

        required_keys = [
            "network_data",
            "starting_node",
            "outfall_node",
            "outfall_node_type",
            "depth",
            "scale",
        ]
        # Ensure all keys are being specified
        missing_keys = [key for key in required_keys if key not in subnetwork_data]
        if missing_keys:
            raise ValueError(f"Missing the following keys: {', '.join(missing_keys)}")

        self.network_data = {}
        self.start = subnetwork_data["starting_node"]
        self.depth = subnetwork_data["depth"]
        self.outfall_node = subnetwork_data["outfall_node"]
        self.outfall_node_type = subnetwork_data["outfall_node_type"]
        self.scale = subnetwork_data["scale"]
        self.subgraph = nx.DiGraph()
        self.node_colour = []
        self.coords = {}

        if len(self.outfall_node) != len(self.outfall_node_type):
            raise ValueError(
                f"Length of {self.outfall_node} and {self.outfall_node_type} are not the same."
            )

        print("Removing Unnecessary Entries")
        keys_to_remove = [
            "GIS_ID",
            "GIS_FeatureKey",
            "FlowType",
            "StreetName",
            "Location",
            "Material",
            "MaterialClass",
            "Slope",
            "AverageDepth",
            "MaximumDepth",
            "PDC_Count",
            "ConstructionDate",
            "Ownership",
            "Assumed",
            "Combined",
            "Relief",
            "LastEditDate",
            "Abandoned",
            "Grouted",
            "CCTVInspectionDate",
            "ShapeSTLength",
        ]

        self.reduced_keys_data = subnetwork_data["network_data"].full_data.copy()

        for keys in keys_to_remove:
            if keys in self.reduced_keys_data.columns:
                self.reduced_keys_data.drop(keys, axis=1, inplace=True)

        print("Creating Pipes")
        for index, pipe in self.reduced_keys_data.iterrows():
            if pipe["UpstreamID"] is not None and pipe["DownstreamID"] is not None:
                if pipe["UpstreamID"] != "Abandoned" and pipe["DownstreamID"] != "Abandoned":
                    self.network_data[pipe["OBJECTID"]] = (pipe["DownstreamID"], pipe["UpstreamID"])

    def create_sub_network(self, ignore_successors=False, ignore_predecessors=False):
        """Create subnetwork, generate node colouring, and generate node locations.

        Returns
        ----------
        ignore_successors : bool
            ignore successor branch from starting node

        ignore_predecessors : bool
            ignore predecessors from starting node
        """
        print("Creating subnetwork for plotting")

        G = self._initialize_digraph()

        # Generate subnetwork starting with successors of starting node
        print("Checking successor nodes")
        successors_subgraph = self.explore_nodes(G, direction="successors")

        # Generate subnetwork with predecessors of starting node
        print("Checking predecessor nodes")
        predecessors_subgraph = self.explore_nodes(G, direction="predecessors")

        if ignore_successors:
            connected_subgraph = predecessors_subgraph
        elif ignore_predecessors:
            connected_subgraph = successors_subgraph
        else:
            connected_subgraph = successors_subgraph | predecessors_subgraph

        # Create a subgraph from the connected nodes
        self.subgraph = G.subgraph(connected_subgraph)

    def _initialize_digraph(self):
        # Initialize Networkx graph
        G = nx.DiGraph()

        # Add edges from network_data
        G.add_edges_from(self.network_data.values())

        return G

    def add_colouring(self):
        """Apply colouring for subnetwork nodes.

        Note: This method is intended for internal use and should not be called directly
        by external code.

        Returns
        ----------
        None
        """

        for node in self.subgraph.nodes:
            if node == self.start or node == self.outfall_node:
                self.node_colour.append("green")
            elif len(list(self.subgraph.successors(node))) > 1:
                self.node_colour.append("red")
            elif len(list(self.subgraph.successors(node))) == 0:
                self.node_colour.append("orange")
            else:
                self.node_colour.append("blue")

    def add_coordinates(self):
        """Apply node coordinates for subnetwork nodes.

        Note: This method is intended for internal use and should not be called directly
        by external code.

        Returns
        ----------
        None
        """

        for node, (x, y) in self.node_positions(self.subgraph.nodes()).items():
            self.coords[node] = (x * self.scale, y * self.scale)

        missing_nodes = set(self.subgraph.nodes()) - set(self.coords.keys())
        if missing_nodes:
            print(f"Error: Positions for nodes {missing_nodes} could not be found.")

    def node_positions(self, nodes):
        """Unscaled (x, y) position of each node, from the GIS pipe geometry.

        A node is placed at the upstream end of a pipe leaving it or, if no pipe leaves it, at the
        downstream end of a pipe entering it.

        Parameters
        ----------
        nodes : iterable
            Nodes to place.

        Returns
        ----------
        dict
            {node: (x, y)} for the nodes that were found; nodes with no pipe are left out.
        """
        nodes = set(nodes)
        upstream_ends, downstream_ends = {}, {}
        for _, pipe in self.reduced_keys_data.iterrows():
            upstream_id, downstream_id = pipe["UpstreamID"], pipe["DownstreamID"]
            if upstream_id in nodes and upstream_id not in upstream_ends:
                upstream_ends[upstream_id] = line_coords(pipe["geometry"])[0]
            if downstream_id in nodes and downstream_id not in downstream_ends:
                downstream_ends[downstream_id] = line_coords(pipe["geometry"])[-1]
        return {**downstream_ends, **upstream_ends}

    def explore_nodes(self, G, direction):
        """Generate sewer network based on supplied direction

        Parameters
        ----------
        G : networkx.digraph(),
            data from Network class

        direction : str
            specifying successor or predecessor

        Return
        ----------
        None
        """
        connected_subgraph = {self.start}
        current_depth = 0
        opposite_direction = "successors" if direction == "predecessors" else "predecessors"

        while current_depth < self.depth:
            new_nodes = set()
            for node in list(connected_subgraph):
                # Find the next node in the sewer network
                next_nodes = set(getattr(G, direction)(node))
                if all(
                    not any(n in node_set for n in next_nodes) for node_set in connected_subgraph
                ):
                    if len(next_nodes) > 1:
                        print(f"Branch Detected at node: {node}")
                        branch_nodes = self.explore_branch(
                            next_nodes, new_nodes, G, direction, opposite_direction
                        )
                        new_nodes.update(branch_nodes)
                        current_depth += 1
                    elif len(next_nodes) == 1:
                        new_nodes.update(next_nodes)
                        opposite_nodes = set(getattr(G, opposite_direction)(next(iter(next_nodes))))
                        print(
                            f"Checking {opposite_direction.capitalize()} of {next_nodes}: {opposite_nodes}"
                        )
                        if next(iter(opposite_nodes)) == node and len(opposite_nodes) == 1:
                            break
                        else:
                            print(f"Possible branch detected at: {next_nodes}")
                            if next_nodes in connected_subgraph:
                                break
                            else:
                                current_depth += 1
                                inner_depth = current_depth
                                while inner_depth < self.depth:
                                    branch_nodes = self.explore_branch(
                                        opposite_nodes, new_nodes, G, opposite_direction, direction
                                    )
                                    new_nodes.update(branch_nodes)
                                    inner_depth += 1

                    elif len(next_nodes) == 0:
                        print("Endpoint reached")
                        new_nodes.update(next_nodes)
                        current_depth += 1
                        break

            connected_subgraph.update(new_nodes)

        return connected_subgraph

    @staticmethod
    def explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction):
        """Exploring branches found in explore_nodes function

        Parameters
        ----------
        branch_nodes : set
            nodes where a branch is found

        new_nodes : set
            nodes that will be added to subnetwork

        G : networkx.digraph(),
            data from Network class

        direction : str
            specifying successor or predecessor

        opposite_direction : str
            value opposite to direction

        Return
        ----------
        new_nodes : set
            nodes to be added to the subnetwork from branches
        """
        for branch in list(branch_nodes):
            print(f"Exploring branch from node: {branch}")
            new_nodes.add(branch)
            # print(f'{direction.capitalize()} of {branch} is: {list(getattr(G, direction)(branch))}')
            branch_next = set(getattr(G, direction)(branch))

            new_nodes.update(branch_next)
            while len(branch_next) == 1:
                next_node = next(iter(branch_next))
                branch_next = set(getattr(G, direction)(next_node))
                if len(branch_next) >= 2:
                    print(f"New branch detected at {next_node}")
                    break
                elif len(branch_next) == 0:
                    print(f"Endpoint reached at {next_node}")
                    new_nodes.add(next_node)
                    new_nodes.update(branch_next)

                    # Check for successors or predecessors based on the current direction
                    opposite_nodes = set(getattr(G, opposite_direction)(next_node))
                    print(
                        f"Checking {opposite_direction.capitalize()} of {next_node}: "
                        f"{opposite_nodes}"
                    )
                    break
                elif branch_next == "Abandoned":
                    print("Reached abandoned node.")
                    break
                else:
                    new_nodes.update(branch_next)
        return new_nodes

    def manipulate_network(self, operation, direction, start_nodes, copy=False):
        """Manipulate graph based on specified operation

        Parameters
        ----------
        operation : str
            operation to perform.
            options: trim, extend

        direction : str
            direction to trim or extend.

        start_nodes: list
            node(s) to start trimming or extending direction.

        copy : bool, optional
            manipulate existing graph or return a new one.

        Returns
        ----------
        nx.Digraph
        """

        G = self._initialize_digraph()

        if len(start_nodes) == 0:
            raise ValueError("No nodes listed to trim. At least 1 node needs to be provided.")

        self.validate_direction(direction)

        # if tuple(str(start_nodes)) not in self.subgraph.nodes():
        #    print(f"Node {start_nodes} is not in the current graph. Exiting manipulate_subnetwork.")
        #    return

        if operation == "trim":
            if copy:
                trimmed_graph = self.subgraph.copy()
                for nodes in start_nodes:
                    remove_nodes = set(getattr(trimmed_graph, direction)(nodes))
                    remove_nodes.update(start_nodes)
                    trimmed_graph.remove_nodes_from(remove_nodes)
                return trimmed_graph
            else:
                self.subgraph = nx.DiGraph(self.subgraph)
                remove_nodes = set()
                for nodes in start_nodes:
                    remove_nodes = set(getattr(self.subgraph, direction)(nodes))
                    remove_nodes.update(start_nodes)
                    self.subgraph.remove_nodes_from(remove_nodes)

        elif operation == "extend":
            if copy:
                extended_graph = self.subgraph.copy()
                for nodes in start_nodes:
                    add_nodes = set(getattr(G, direction)(nodes))
                    add_nodes.update(start_nodes)
                    extended_graph.update(add_nodes)
                return extended_graph
            else:
                self.subgraph = nx.DiGraph(self.subgraph)
                for nodes in start_nodes:
                    neighbour = getattr(G, direction)(nodes)
                    add_nodes = set(neighbour)
                    add_nodes.update(start_nodes)
                    print(add_nodes)
                    self.subgraph.add_nodes_from(add_nodes)
        else:
            raise ValueError("Invalid operation specified. Must be: extend or trim")

        # Add the ability to trim / extend branches now

    @staticmethod
    def validate_direction(direction):
        """Validate the direction parameter.

        Parameters:
        ----------
        direction : str
            The direction parameter to be validated.

        Raises
        ----------
        ValueError
            If the direction is not 'successors' or 'predecessors'.
        """

        if direction not in ["successors", "predecessors"]:
            raise ValueError("Invalid direction. Use 'successors' or 'predecessors'.")

    def write_nodes_pyswmm(
        self, pyswmm_file, max_depth, init_depth, surcharge_depth, ponding_area, height_correction
    ):
        """Add subnetwork data into the Nodes section of a .inp file.

        Parameters
        ----------
        pyswmm_file : str
            Path to the write out PySWMM input file (.inp).

        max_depth : float)
            Maximum depth for node (default is 0.0).

        init_depth : float
            Initial depth for node (default is 0.0).

        surcharge_depth : float
            Surcharge depth for node (default is 0.0).

        ponding_area : float
            Ponding area for node (default is 0.0).

        height_correction : dict
            User defined height correction in the case that the inv. elev. is not correct in geojson file.

        Notes
        -----
        Also fills the {coordinates} placeholder with each node's position in the units of the
        GIS geometry (unscaled; see node_positions).
        """

        with open(pyswmm_file) as file:
            content = file.read()

        nodes_content = ""
        outfall_content = ""
        # Iterate through the nodes in the graph and build the content
        for node in self.subgraph.nodes:
            if height_correction is not None and node in height_correction:
                upstream_invert = height_correction[node]
                nodes_content += (
                    f"{node:<16} {upstream_invert:<11} {max_depth:<11} {init_depth:<11}"
                    f"{surcharge_depth:<11} {ponding_area:<11}\n"
                )
                print(f"Node {node}: Updating upstream_invert value to: {upstream_invert}")
            else:
                if node in self.outfall_node:
                    for idx, outfall_node in enumerate(self.outfall_node):
                        upstream_invert = self.get_upstream_invert(outfall_node)
                        outfall_node_type = self.outfall_node_type[idx]
                        outfall_content += f"{node:<16} {upstream_invert:<11}  {outfall_node_type:<28} {'NO':<10}\n"
                else:
                    upstream_invert = self.get_upstream_invert(node)
                    nodes_content += (
                        f"{node:<16} {upstream_invert:<11} {max_depth:<11} {init_depth:<11}"
                        f"{surcharge_depth:<11} {ponding_area:<11}\n"
                    )

        coordinates_content = ""
        positions = self.node_positions(self.subgraph.nodes)
        for node in self.subgraph.nodes:
            if node in positions:
                x, y = positions[node]
                coordinates_content += f"{node:<16} {x:<18} {y:<18}\n"
            else:
                print(f"Warning: Node {node} has no position; it is left out of [COORDINATES]")

        # Replace {nodes} placeholder with the built nodes content
        modified_content = (
            content.replace("{nodes}", nodes_content)
            .replace("{outfalls}", outfall_content)
            .replace("{coordinates}", coordinates_content)
        )

        # Write the modified content back to the file
        with open(pyswmm_file, "w") as file:
            file.write(modified_content)

    def write_conduits_pyswmm(
        self, pyswmm_file, manning_n, inlet_offset, outlet_offset, init_flow, max_flow, barrels
    ):
        """Add subnetwork data into the Conduits and XSection section of a .inp file.

        Parameters
        ----------

        pyswmm_file : str
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

        Notes
        -----
        Also fills the {vertices} placeholder with the interior vertices of multi-vertex pipes, in
        the units of the GIS geometry (unscaled).
        """

        with open(pyswmm_file) as file:
            content = file.read()

        conduits_content = ""
        xsection_content = ""
        vertices_content = ""

        # Get information from nx.digraph to fill in data
        for edges in self.subgraph.edges:
            start_node, end_node = edges

            for index, row in self.full_data.iterrows():
                if row["UpstreamID"] == end_node and row["DownstreamID"] == start_node:
                    print("Found pair")

                    name = row["OBJECTID"]
                    inlet = row["UpstreamID"]
                    outlet = row["DownstreamID"]
                    length = row["Length"]

                    if row["PipeShape"] == "R":
                        shape = "RECT_CLOSED"
                        height = row["Height"] / 1000
                        width = row["Width"] / 1000
                    elif row["PipeShape"] == "Round":
                        shape = "CIRCULAR"
                        height = row["Height"] / 1000
                        width = "0.0"

                    conduits_content += (
                        f"{name:<16} {inlet:<16} {outlet:<16}"
                        f"{length:<11} {manning_n:<11} "
                        f"{inlet_offset:<11} {outlet_offset:<11}"
                        f"{init_flow:<11} {max_flow:<11}\n"
                    )

                    xsection_content += (
                        f"{name:<16} {shape:<16} {height:<16}"
                        f"{width:<11} {0.0:<11} {0.0:<11} {barrels:<11}\n"
                    )

                    for x, y in line_coords(row["geometry"])[1:-1]:
                        vertices_content += f"{name:<16} {x:<18} {y:<18}\n"

        # Replace variable in template file with conduit information
        modified_content = (
            content.replace("{conduits}", conduits_content)
            .replace("{xsections}", xsection_content)
            .replace("{vertices}", vertices_content)
        )

        with open(pyswmm_file, "w") as file:
            file.write(modified_content)

    def write_pyswmm(
        self,
        output_file=None,
        max_depth=0.0,
        init_depth=0.0,
        surcharge_depth=0.0,
        ponding_area=0.0,
        manning_n=0.01,
        inlet_offset=0.0,
        outlet_offset=0.0,
        init_flow=0.0,
        max_flow=0.0,
        barrels=1.0,
        template_path=None,
        output_folder=None,
        height_correction=None,
    ):
        """Generate a PySWMM input file based on a networkx digraph and GeoPandas DataFrame.

        Parameters
        ----------
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
            Path to the pyswmm model template file. The default is the template shipped with the
            package (``sewernetworkx/templates/template.inp``), which uses SI units (FLOW_UNITS
            LPS). A custom template must use an SI FLOW_UNITS (CMS, LPS or MLD), since the GIS
            lengths, inverts and sizes are in metres, and may only contain the placeholders
            {nodes}, {outfalls}, {conduits}, {xsections}, {coordinates} and {vertices}.
            Pumps are not supported yet, so a template with a {pumps} placeholder is rejected.

        output_folder : str
            Path to were the new .inp file will be written.

        height_correction : dict
            If heights are not correct in geojson file, this overwrites it with a user defined value.
        """

        if output_file is None:
            raise ValueError("The 'output_path' must be specified.")

        if template_path is None:
            template = default_template_path().read_text()
        else:
            if not os.path.exists(template_path):
                raise FileNotFoundError(
                    f"The template file does not exist at the specified path: {template_path}"
                )
            with open(template_path) as file:
                template = file.read()

        check_template(template)

        # Always create a copy of the template file
        if output_folder is not None:
            output_file_path = os.path.join(output_folder, output_file)
        else:
            # If output_folder is None, create a copy in the current working directory
            output_file_path = os.path.join(os.getcwd(), output_file)

        with open(output_file_path, "w") as file:
            file.write(template)
        print(f"Template file copied to: {output_file_path}")

        # Write nodes to pyswmm .inp file
        self.write_nodes_pyswmm(
            output_file_path,
            max_depth,
            init_depth,
            surcharge_depth,
            ponding_area,
            height_correction,
        )

        # Write conduits to pyswmm inp file
        self.write_conduits_pyswmm(
            output_file_path, manning_n, inlet_offset, outlet_offset, init_flow, max_flow, barrels
        )

    def generate_gmap(self, file_name, scale=1, origin=None):
        """Generate an hmtl file to show the sewer network overlaid on Google Maps.

        Parameters
        ----------
        file_name : str
            Name of the html file

        scale : int
            Value to scale coordinate data by. Default value is 1.

        origin : list
            Location of centre for the overlaid map.
        """

        self.validate_file_extension(file_name, ".html")

        if origin is not None:
            gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], zoom=18)
        else:
            if not self.coords:
                raise ValueError(
                    "Coordinates have not been generated for this graph. "
                    "You must run the function '.add_coordinates' prior to "
                    "running 'generate.gmap'."
                )

            coords = next(iter(self.coords.values()))
            gmap = gmplot.GoogleMapPlotter(coords[1] / scale, coords[0] / scale, zoom=18)

        for node, coords in self.coords.items():
            gmap.marker(coords[1] / scale, coords[0] / scale, "cornflowerblue", label=str(node))

        for edge in self.subgraph.edges():
            start_node_pos = self.coords[edge[0]]
            end_node_pos = self.coords[edge[1]]
            gmap.plot(
                [start_node_pos[1] / scale, end_node_pos[1] / scale],
                [start_node_pos[0] / scale, end_node_pos[0] / scale],
                "red",
                edge_width=2,
            )

        gmap.draw(file_name)

    @staticmethod
    def validate_file_extension(file_name, expected_extension):
        _, file_extension = os.path.splitext(file_name)

        if file_extension != expected_extension:
            raise ValueError(
                f"Wrong file extension used. The proper extension is: {expected_extension}."
            )
