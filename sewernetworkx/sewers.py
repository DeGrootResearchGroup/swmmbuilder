"""Module for manipulating GeoJSON sewer data.

The module takes full GeoJSON data sets of a city sewer network,
generates a subnetwork given an initial starting node and number
of levels, and generates a .inp file to be used by pyswmm.
"""

import os
import shutil
import warnings
from pathlib import Path

import geopandas as gpd
import gmplot
import networkx as nx
import pandas as pd
import pkg_resources


def _pipe_records(data):
    """Yield (OBJECTID, UpstreamID, DownstreamID, Length) for every pipe with two node IDs.

    Pipes with a missing node ID, or a node ID of 'Abandoned', are skipped.
    """
    for oid, up, down, length in zip(
        data["OBJECTID"], data["UpstreamID"], data["DownstreamID"], data["Length"], strict=True
    ):
        if pd.isna(up) or pd.isna(down) or "Abandoned" in (up, down):
            continue
        yield oid, up, down, length


class Network:
    """Class used to represent an entire cities sewer network

    Parameters
    ----------
    file_path : str, Path
        Name of GeoJSON file

    Attributes
    ----------
    full_data : geopandas.GeoDataFrame
        Pipe data read from the GeoJSON file.

    graph : networkx.MultiDiGraph
        One edge per pipe, in the direction of flow (UpstreamID -> DownstreamID), keyed by
        OBJECTID so that parallel pipes are kept. Each edge has the attributes ``pipe`` (OBJECTID)
        and ``length``. Nodes and edges are inserted in sorted order, so traversals do not depend
        on the order of the file.

    Notes
    -----
    remove_store is used as a flag to automatically remove any sewer pipe
    flagged as a storm drain.
    """

    def __init__(self, file_path, remove_storm=True):
        """Class initialization"""

        self.file_path = file_path
        data = gpd.read_file(self.file_path)
        self.full_data = gpd.GeoDataFrame(data, geometry="geometry")

        if remove_storm:
            print("Removing storm sewers")
            key = "FlowType"
            condition = self.full_data[key] != "STM"
            self.full_data = self.full_data[condition]

        self.graph = self._build_graph(self.full_data)

    @staticmethod
    def _build_graph(data):
        """Build the flow-direction MultiDiGraph of the pipes in ``data``."""
        records = sorted(_pipe_records(data), key=lambda r: (r[1], r[2], r[0]))
        G = nx.MultiDiGraph()
        G.add_nodes_from(sorted({n for _, up, down, _ in records for n in (up, down)}))
        for oid, up, down, length in records:
            G.add_edge(up, down, key=oid, pipe=oid, length=length)
        return G

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
        Path to GeoJSON file. If given, the pipes are read from this file instead of
        ``network_data``.

    Notes
    -----
    Select the subnetwork with one of :meth:`catchment`, :meth:`from_nodes` or
    :meth:`from_node_list`; each sets ``subgraph`` deterministically. ``subgraph`` is a
    ``networkx.DiGraph`` whose edges are stored as (downstream, upstream), the orientation the
    .inp writer expects, with the attribute ``pipes``: the sorted OBJECTIDs of the (possibly
    parallel) pipes between the two nodes.
    """

    def __init__(self, subnetwork_data, file_path=None):
        if file_path:
            super().__init__(file_path)
            source = self
        else:
            source = subnetwork_data.get("network_data")

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

        self.full_data = source.full_data
        self.graph = source.graph
        self.reduced_keys_data = source.full_data.copy()

        for keys in keys_to_remove:
            if keys in self.reduced_keys_data.columns:
                self.reduced_keys_data.drop(keys, axis=1, inplace=True)

        print("Creating Pipes")
        for oid, up, down, _ in _pipe_records(self.reduced_keys_data):
            self.network_data[oid] = (down, up)

    def catchment(self, node, max_hops=None, max_length=None, outlet=None):
        """Select the catchment of ``node``: every node that drains to it.

        With ``max_hops`` or ``max_length`` the catchment is truncated, and the path from ``node``
        down to the outlet is added so that the subnetwork still discharges there.

        Parameters
        ----------
        node : str
            Node whose catchment is selected. Pass the outlet itself for its whole catchment.

        max_hops : int, optional
            Keep only nodes at most this many pipes upstream of ``node``.

        max_length : float, optional
            Keep only nodes whose shortest upstream path to ``node`` is at most this long, in the
            units of the ``Length`` column.

        outlet : str or list of str, optional
            Outlet(s) to connect ``node`` to. Defaults to ``outfall_node``. A shortest path (by
            length) is added to every listed outlet that ``node`` drains to.

        Returns
        ----------
        networkx.DiGraph
            The new ``subgraph``.
        """
        self._check_nodes([node])
        upstream = self.graph.reverse(copy=False)

        if max_hops is None:
            selected = set(nx.descendants(upstream, node))
        else:
            selected = set(nx.single_source_shortest_path_length(upstream, node, cutoff=max_hops))
        if max_length is not None:
            within = nx.single_source_dijkstra_path_length(
                upstream, node, cutoff=max_length, weight=self._pipe_length
            )
            selected &= set(within)
        selected.add(node)

        if max_hops is not None or max_length is not None:
            outlets = self.outfall_node if outlet is None else outlet
            outlets = [outlets] if isinstance(outlets, str) else list(outlets)
            self._check_nodes(outlets)
            reachable = [o for o in sorted(outlets) if nx.has_path(self.graph, node, o)]
            if not reachable:
                raise ValueError(f"Node {node} does not drain to any of the outlets {outlets}.")
            for o in reachable:
                selected.update(nx.dijkstra_path(self.graph, node, o, weight=self._pipe_length))

        return self.from_nodes(selected)

    def from_nodes(self, nodes):
        """Select the subnetwork made of ``nodes`` and every pipe between two of them.

        Parameters
        ----------
        nodes : iterable of str
            Node IDs.

        Returns
        ----------
        networkx.DiGraph
            The new ``subgraph``.
        """
        nodes = sorted(set(nodes))
        self._check_nodes(nodes)
        induced = self.graph.subgraph(nodes)

        subgraph = nx.DiGraph()
        subgraph.add_nodes_from(nodes)
        pipes = {}
        for up, down, oid in induced.edges(keys=True):
            pipes.setdefault((down, up), []).append(oid)
        for edge in sorted(pipes):
            subgraph.add_edge(*edge, pipes=sorted(pipes[edge]))

        self.subgraph = subgraph
        return subgraph

    def to_node_list(self, path):
        """Write the nodes of ``subgraph`` to a text file, one per line, in sorted order.

        :meth:`from_node_list` rebuilds the same subnetwork from this file.

        Parameters
        ----------
        path : str, Path
            File to write.
        """
        lines = ["# SewerNetworkx node list: one node ID per line"]
        lines += sorted(self.subgraph.nodes)
        Path(path).write_text("\n".join(lines) + "\n")

    def from_node_list(self, path):
        """Select the subnetwork whose nodes are listed in a file written by :meth:`to_node_list`.

        Blank lines and lines starting with '#' are ignored.

        Parameters
        ----------
        path : str, Path
            File to read.

        Returns
        ----------
        networkx.DiGraph
            The new ``subgraph``.
        """
        lines = (line.strip() for line in Path(path).read_text().splitlines())
        return self.from_nodes(line for line in lines if line and not line.startswith("#"))

    def _check_nodes(self, nodes):
        missing = sorted(n for n in nodes if n not in self.graph)
        if missing:
            raise ValueError(f"Nodes not in the network: {', '.join(map(str, missing))}")

    @staticmethod
    def _pipe_length(u, v, pipes):
        """Dijkstra weight: the shortest of the (possibly parallel) pipes from u to v."""
        return min(attrs["length"] for attrs in pipes.values())

    def create_sub_network(self, ignore_successors=False, ignore_predecessors=False):
        """Create subnetwork, generate node colouring, and generate node locations.

        .. deprecated::
            The selected nodes are not a well-defined function of ``depth``. Use
            :meth:`catchment`, :meth:`from_nodes` or :meth:`from_node_list` instead.

        Returns
        ----------
        ignore_successors : bool
            ignore successor branch from starting node

        ignore_predecessors : bool
            ignore predecessors from starting node
        """
        warnings.warn(
            "create_sub_network is deprecated; use catchment, from_nodes or from_node_list.",
            DeprecationWarning,
            stacklevel=2,
        )
        print("Creating subnetwork for plotting")

        G = self._initialize_digraph()

        # Generate subnetwork starting with successors of starting node
        print("Checking successor nodes")
        successors_subgraph = self._explore_nodes(G, direction="successors")

        # Generate subnetwork with predecessors of starting node
        print("Checking predecessor nodes")
        predecessors_subgraph = self._explore_nodes(G, direction="predecessors")

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

        for index, row in self.reduced_keys_data.iterrows():
            upstream_id = row["UpstreamID"]
            downstream_id = row["DownstreamID"]
            if upstream_id == self.start:
                # print('Plant node')
                coordinates = row["geometry"].coords[0]
                self.coords[downstream_id] = (
                    coordinates[0] * self.scale,
                    coordinates[1] * self.scale,
                )
            elif upstream_id in self.subgraph.nodes():
                # print('Found upstream node')
                coordinates = row["geometry"].coords[0]
                self.coords[upstream_id] = (
                    coordinates[0] * self.scale,
                    coordinates[1] * self.scale,
                )

        # Check if each node in network_nodes has a position in pos
        missing_nodes = set(self.subgraph.nodes()) - set(self.coords.keys())
        if missing_nodes:
            print(
                f"Warning: Nodes {missing_nodes} do not have positions. Attempting to find positions for missing nodes."
            )

            for index, row in self.reduced_keys_data.iterrows():
                downstream_id = row["DownstreamID"]

                if downstream_id in missing_nodes:
                    coordinates = row["geometry"].coords[1]
                    self.coords[downstream_id] = (
                        coordinates[0] * self.scale,
                        coordinates[1] * self.scale,
                    )
                    missing_nodes.remove(downstream_id)

        if missing_nodes:
            print(f"Error: Positions for nodes {missing_nodes} could not be found.")

    def explore_nodes(self, G, direction):
        """Generate sewer network based on supplied direction

        .. deprecated::
            Used by the deprecated :meth:`create_sub_network`.
        """
        warnings.warn(
            "explore_nodes is deprecated; use catchment, from_nodes or from_node_list.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self._explore_nodes(G, direction)

    def _explore_nodes(self, G, direction):
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
            for node in sorted(connected_subgraph):
                # Find the next node in the sewer network
                next_nodes = set(getattr(G, direction)(node))
                if all(
                    not any(n in node_set for n in next_nodes) for node_set in connected_subgraph
                ):
                    if len(next_nodes) > 1:
                        print(f"Branch Detected at node: {node}")
                        branch_nodes = self._explore_branch(
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
                                    branch_nodes = self._explore_branch(
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

        .. deprecated::
            Used by the deprecated :meth:`create_sub_network`.
        """
        warnings.warn(
            "explore_branch is deprecated; use catchment, from_nodes or from_node_list.",
            DeprecationWarning,
            stacklevel=2,
        )
        return SubNetwork._explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction)

    @staticmethod
    def _explore_branch(branch_nodes, new_nodes, G, direction, opposite_direction):
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
        for branch in sorted(branch_nodes):
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

        # Replace {nodes} placeholder with the built nodes content
        modified_content = content.replace("{nodes}", nodes_content).replace(
            "{outfalls}", outfall_content
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

        """

        with open(pyswmm_file) as file:
            content = file.read()

        conduits_content = ""
        xsection_content = ""

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

        # Replace variable in template file with conduit information
        modified_content = content.replace("{conduits}", conduits_content).replace(
            "{xsections}", xsection_content
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
            Path to the pyswmm model template file (default is "template.inp").

        output_folder : str
            Path to were the new .inp file will be written.

        height_correction : dict
            If heights are not correct in geojson file, this overwrites it with a user defined value.
        """

        if output_file is None:
            raise ValueError("The 'output_path' must be specified.")

        if template_path is None:
            # Get the package directory and construct the path to the template file
            package_dir = pkg_resources.resource_filename(__name__, "")
            template_path = os.path.join(package_dir, "examples", "pyswmm_template", "template.inp")

        if not os.path.exists(template_path):
            raise FileNotFoundError(
                f"The template file does not exist at the specified path: {template_path}"
            )

        # Always create a copy of the template file
        if output_folder is not None:
            output_file_path = os.path.join(output_folder, output_file)
        else:
            # If output_folder is None, create a copy in the current working directory
            output_file_path = os.path.join(os.getcwd(), output_file)

        shutil.copy(template_path, output_file_path)
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

            node, coords = list(self.coords.items())[0]
            gmap = gmplot.GoogleMapPlotter(coords[1] / scale, coords[0] / scale, zoom=18)

        for node, coords in self.coords.items():
            gmap.marker(coords[1] / 10, coords[0] / 10, "cornflowerblue", label=str(node))

        for edge in self.subgraph.edges():
            start_node_pos = self.coords[edge[0]]
            end_node_pos = self.coords[edge[1]]
            gmap.plot(
                [start_node_pos[1] / 10, end_node_pos[1] / 10],
                [start_node_pos[0] / 10, end_node_pos[0] / 10],
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
