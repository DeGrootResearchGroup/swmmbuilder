"""Interactive Leaflet maps and GeoJSON export."""

from __future__ import annotations

import importlib
import os
import warnings
from dataclasses import dataclass, field, fields
from typing import TYPE_CHECKING

import geopandas as gpd
import matplotlib as mpl
import numpy as np
from matplotlib.colors import to_hex

from ._deps import _optional
from .colours import (
    COLOUR_OPTIONS,
    ColourScale,
    _aligned,
    _colour_options,
    _fitted_norm,
    _hex_colours,
)
from .crs import WGS84, _check_nodes, _with_crs, lonlat_positions

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from .._types import StrPath, Values
    from ..network import Network
    from ..subnetwork import Subnetwork


@dataclass
class PipeLayer:
    """Pipes on a :func:`to_leaflet` map, coloured by ``values`` (per pipe ID) like
    :func:`plot_pipes`, or all in ``color``. Tooltips show the pipe ID and value."""

    values: object = None
    name: str = "Pipes"
    cmap: object = "viridis"
    norm: object = "linear"
    vmin: float | None = None
    vmax: float | None = None
    color: str = "#3366cc"
    missing_color: str = "#bbbbbb"
    weight: float = 3.0
    show: bool = True
    colour: ColourScale | None = None

    def __post_init__(self):
        # A ColourScale sets cmap, norm, vmin, vmax and missing_color (see plot_pipes).
        defaults = tuple(_FIELD_DEFAULTS[type(self).__name__][n] for n in COLOUR_OPTIONS)
        given = tuple(getattr(self, n) for n in COLOUR_OPTIONS)
        for name, value in zip(
            COLOUR_OPTIONS, _colour_options(self.colour, given, defaults), strict=True
        ):
            setattr(self, name, value)


@dataclass
class GroupLayer:
    """Group polygons on a :func:`to_leaflet` map (e.g. from :func:`node_group_polygons`),
    shaded by ``values`` (per group) or all in ``color``. Tooltips show the group and value."""

    polygons: gpd.GeoDataFrame
    values: object = None
    name: str = "Groups"
    cmap: object = "viridis"
    norm: object = "linear"
    vmin: float | None = None
    vmax: float | None = None
    color: str = "#999999"
    missing_color: str = "#dddddd"
    opacity: float = 0.5
    show: bool = True
    colour: ColourScale | None = None

    def __post_init__(self):
        # A ColourScale sets cmap, norm, vmin, vmax and missing_color (see plot_pipes).
        defaults = tuple(_FIELD_DEFAULTS[type(self).__name__][n] for n in COLOUR_OPTIONS)
        given = tuple(getattr(self, n) for n in COLOUR_OPTIONS)
        for name, value in zip(
            COLOUR_OPTIONS, _colour_options(self.colour, given, defaults), strict=True
        ):
            setattr(self, name, value)


# The declared defaults of the colour fields of each layer, for ColourScale clash checks.
_FIELD_DEFAULTS = {
    cls.__name__: {f.name: f.default for f in fields(cls)} for cls in (PipeLayer, GroupLayer)
}


@dataclass
class MarkerLayer:
    """Node markers on a :func:`to_leaflet` map. Tooltips show the node ID, or ``labels[node]``."""

    node_ids: list
    name: str = "Markers"
    labels: dict = field(default_factory=dict)
    color: str = "crimson"
    radius: float = 6.0
    show: bool = True


def _colour_legend(m, cmap, norm, caption):
    branca = importlib.import_module("branca.colormap")
    stops = np.linspace(0, 1, 11)
    index = [float(norm.inverse(s)) for s in stops]
    colours = [to_hex(c) for c in mpl.colormaps.get_cmap(cmap)(stops)]
    branca.LinearColormap(
        colours, index=index, vmin=index[0], vmax=index[-1], caption=caption
    ).add_to(m)


def _value_colours(values, layer, default_colour):
    """(hex colours, fitted norm) of ``values``, or (``default_colour`` each, None); the norm
    is None too when no value can be coloured, so no legend is drawn."""
    if layer.values is None:
        return [default_colour] * len(values), None
    norm, has_data = _fitted_norm(layer.norm, layer.vmin, layer.vmax, values)
    return _hex_colours(values, layer.cmap, norm, layer.missing_color), norm if has_data else None


def _add_pipes(m, folium, net, layer):
    pipes = _with_crs(net.to_geodataframe("pipes")).to_crs(WGS84)
    data = gpd.GeoDataFrame({"geometry": pipes.geometry.to_numpy()}, crs=WGS84)
    data.insert(0, "pipe_id", pipes.index.astype(str))
    fields, aliases = ["pipe_id"], ["Pipe"]
    values = np.full(len(pipes), np.nan)
    if layer.values is not None:
        values = _aligned(layer.values, pipes.index)
        data.insert(1, "value", values)
        fields, aliases = [*fields, "value"], [*aliases, layer.name]
    colours, norm = _value_colours(values, layer, layer.color)
    data["_colour"] = colours
    group = folium.FeatureGroup(name=layer.name, show=layer.show)
    if len(data):  # folium's tooltip fails on no features
        folium.GeoJson(
            data.to_json(na="null"),
            style_function=lambda f: {"color": f["properties"]["_colour"], "weight": layer.weight},
            tooltip=folium.GeoJsonTooltip(fields=fields, aliases=aliases),
        ).add_to(group)
    group.add_to(m)
    if norm is not None:
        _colour_legend(m, layer.cmap, norm, layer.name)


def _add_groups(m, folium, layer):
    polygons = _with_crs(layer.polygons).to_crs(WGS84)
    data = gpd.GeoDataFrame({"geometry": polygons.geometry.to_numpy()}, crs=WGS84)
    data.insert(0, "group", polygons.index.astype(str))
    fields, aliases = ["group"], ["Group"]
    values = np.full(len(polygons), np.nan)
    if layer.values is not None:
        values = _aligned(layer.values, polygons.index, "group values")
        data.insert(1, "value", values)
        fields, aliases = [*fields, "value"], [*aliases, layer.name]
    colours, norm = _value_colours(values, layer, layer.color)
    data["_colour"] = colours
    group = folium.FeatureGroup(name=layer.name, show=layer.show)
    if len(data):  # folium's tooltip fails on no features
        folium.GeoJson(
            data.to_json(na="null"),
            style_function=lambda f: {
                "fillColor": f["properties"]["_colour"],
                "color": "#444444",
                "weight": 1,
                "fillOpacity": layer.opacity,
            },
            tooltip=folium.GeoJsonTooltip(fields=fields, aliases=aliases),
        ).add_to(group)
    group.add_to(m)
    if norm is not None:
        _colour_legend(m, layer.cmap, norm, layer.name)


def _add_markers(m, folium, net, layer):
    _check_nodes(net, layer.node_ids)
    nodes = _with_crs(net.to_geodataframe("nodes")).to_crs(WGS84).geometry
    unplaced = [str(n) for n in layer.node_ids if n not in nodes.index or nodes[n] is None]
    if unplaced:  # as plot_markers
        raise ValueError(f"Nodes with no position: {', '.join(unplaced)}")
    group = folium.FeatureGroup(name=layer.name, show=layer.show)
    for node in layer.node_ids:
        folium.CircleMarker(
            location=[nodes[node].y, nodes[node].x],
            radius=layer.radius,
            color=layer.color,
            fill=True,
            fill_opacity=0.9,
            tooltip=str(layer.labels.get(node, node)),
        ).add_to(group)
    group.add_to(m)


def to_leaflet(
    net: Network | Subnetwork,
    path: StrPath,
    layers: Sequence[PipeLayer | GroupLayer | MarkerLayer] | None = None,
    *,
    tiles: str | None = "OpenStreetMap",
) -> Any:
    """Write an interactive Leaflet map of ``net`` to an HTML file, with ``folium``.

    The data are embedded in the file; Leaflet and the basemap tiles are loaded from the web
    when it is opened. Needs the ``[maps]`` extra.

    Parameters
    ----------
    net : Network or Subnetwork
        The network.

    path : str or Path
        HTML file to write.

    layers : list of PipeLayer, GroupLayer and MarkerLayer, optional
        Layers, drawn in order and listed in a layer switcher. Default: the pipes.

    tiles : str
        Basemap tiles, as for :class:`folium.Map` (None for none).

    Returns
    -------
    folium.Map
        The map that was written.
    """
    folium = _optional("folium")
    layers = [PipeLayer()] if layers is None else layers
    m = folium.Map(tiles=tiles, control_scale=True)
    bounds = _with_crs(net.to_geodataframe("nodes")).to_crs(WGS84).total_bounds
    if np.isfinite(bounds).all():  # nodes with no position (or no nodes) have none
        xmin, ymin, xmax, ymax = bounds
        m.fit_bounds([[ymin, xmin], [ymax, xmax]])
    for layer in layers:
        if isinstance(layer, PipeLayer):
            _add_pipes(m, folium, net, layer)
        elif isinstance(layer, GroupLayer):
            _add_groups(m, folium, layer)
        elif isinstance(layer, MarkerLayer):
            _add_markers(m, folium, net, layer)
        else:
            raise TypeError(f"Not a map layer: {layer!r}")
    folium.LayerControl(collapsed=False).add_to(m)
    m.save(str(path))
    return m


def to_geojson(
    net: Network | Subnetwork, path: StrPath, kind: str = "pipes", values: Values | None = None
) -> None:
    """Write the pipes or nodes of ``net`` (with ``values``; see
    :meth:`~swmmbuilder.Network.to_geodataframe`) to a GeoJSON file in WGS 84."""
    data = _with_crs(net.to_geodataframe(kind, values)).to_crs(WGS84)
    data.to_file(path, driver="GeoJSON")


def generate_gmap(
    subnetwork: Subnetwork, file_name: str, origin: tuple[float, float] | None = None
) -> None:
    """Write an html file showing the subnetwork overlaid on Google Maps.

    .. deprecated::
        Use :func:`to_leaflet`, which needs no API key and colours pipes by value.

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
    warnings.warn(
        "generate_gmap is deprecated; use swmmbuilder.maps.to_leaflet instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    if os.path.splitext(file_name)[1] != ".html":
        raise ValueError("Wrong file extension used. The proper extension is: .html.")

    positions = lonlat_positions(subnetwork)
    if origin is None:
        if not positions:
            raise ValueError("The subnetwork has no node positions to centre the map on.")
        lon, lat = next(iter(positions.values()))
        origin = (lat, lon)
    gmplot = _optional("gmplot")
    gmap = gmplot.GoogleMapPlotter(origin[0], origin[1], zoom=18)

    for node, (lon, lat) in positions.items():
        gmap.marker(lat, lon, "cornflowerblue", label=str(node))

    for up, down in sorted(set(subnetwork.graph.edges())):
        (lon0, lat0), (lon1, lat1) = positions[up], positions[down]
        gmap.plot([lat0, lat1], [lon0, lon1], "red", edge_width=2)

    gmap.draw(file_name)
