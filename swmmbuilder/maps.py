"""Maps of a subnetwork."""

import os

import gmplot
from pyproj import CRS, Transformer

WGS84 = CRS.from_epsg(4326)


def lonlat_positions(subnetwork):
    """{node: (longitude, latitude)} of the subnetwork's nodes, in WGS 84.

    Positions in another CRS (e.g. a projected one such as :data:`TORONTO_ON`'s EPSG:2952) are
    reprojected. A network with no CRS is taken to be in longitude/latitude already.
    """
    positions = subnetwork.node_positions()
    crs = subnetwork.network.full_data.crs
    if crs is None or CRS.from_user_input(crs) == WGS84:
        return positions
    transformer = Transformer.from_crs(crs, WGS84, always_xy=True)
    return {node: transformer.transform(x, y) for node, (x, y) in positions.items()}


def generate_gmap(subnetwork, file_name, origin=None):
    """Write an html file showing the subnetwork overlaid on Google Maps.

    The node positions are reprojected to WGS 84 (see :func:`lonlat_positions`), so networks in
    a projected CRS are drawn in the right place.

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

    positions = lonlat_positions(subnetwork)
    if origin is None:
        if not positions:
            raise ValueError("The subnetwork has no node positions to centre the map on.")
        lon, lat = next(iter(positions.values()))
        origin = (lat, lon)
    gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], zoom=18)

    for node, (lon, lat) in positions.items():
        gmap.marker(lat, lon, "cornflowerblue", label=str(node))

    for up, down in sorted(set(subnetwork.graph.edges())):
        (lon0, lat0), (lon1, lat1) = positions[up], positions[down]
        gmap.plot([lat0, lat1], [lon0, lon1], "red", edge_width=2)

    gmap.draw(file_name)
