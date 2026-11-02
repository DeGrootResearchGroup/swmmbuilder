"""Maps of a subnetwork."""

import os

import gmplot


def generate_gmap(subnetwork, file_name, origin=None):
    """Write an html file showing the subnetwork overlaid on Google Maps.

    The node positions are used as (longitude, latitude), so the network must be in a geographic
    CRS (e.g. WGS 84).

    Parameters
    ----------
    subnetwork : Subnetwork
        The subnetwork.

    file_name : str
        Name of the html file.

    origin : (float, float), optional
        (latitude, longitude) of the centre of the map. Defaults to the first node.
    """
    if os.path.splitext(file_name)[1] != ".html":
        raise ValueError("Wrong file extension used. The proper extension is: .html.")

    positions = subnetwork.node_positions()
    if origin is None:
        if not positions:
            raise ValueError("The subnetwork has no node positions to centre the map on.")
        x, y = next(iter(positions.values()))
        origin = (y, x)
    gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], zoom=18)

    for node, (x, y) in positions.items():
        gmap.marker(y, x, "cornflowerblue", label=str(node))

    for up, down in sorted(set(subnetwork.graph.edges())):
        (x0, y0), (x1, y1) = positions[up], positions[down]
        gmap.plot([y0, y1], [x0, x1], "red", edge_width=2)

    gmap.draw(file_name)
