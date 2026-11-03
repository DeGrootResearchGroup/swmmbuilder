"""Maps of a network or subnetwork: static plots with matplotlib and interactive Leaflet maps.

The functions here only draw: they take per-pipe, per-node or per-group values computed elsewhere
(concentrations, flows, the share of each upstream area in a pipe's load, ...), as pandas Series
or DataFrames indexed by pipe ID, node ID or group name.

Static plots take a :class:`~swmmbuilder.Network` or :class:`~swmmbuilder.Subnetwork`
(``net``) and draw in a metric CRS (see :func:`plot_crs`). Functions drawing on the same axes must
use the same ``crs``; the default is the same for a network and all its subnetworks.

Units: pipe line widths and marker sizes are in points, so they stay the same on screen when you
zoom; ribbon widths, footprints and buffers are in metres of the plot CRS, so they scale with the
map.

Each static plot styles the axes as a map (:func:`map_axes`) and returns a :class:`DrawnLayer`.
With ``legend=True`` (the default) its colour bar or legend is stacked to the right of the axes
(see :func:`add_legends`), below those of earlier plots on the same axes.

``contextily`` (web-tile basemaps) and ``folium`` (:func:`to_leaflet`) are optional: install them
with ``pip install 'swmmbuilder[maps]'``.

The code is split into modules (crs, colours, style, polygons, ribbons, static, basemap,
export); every name is re-exported here, so ``maps.X`` keeps working.
"""

import logging

from ..geometry import line_coords  # noqa: F401
from ._deps import (  # noqa: F401
    _optional,
)
from .basemap import (  # noqa: F401
    add_basemap,
)
from .colours import (  # noqa: F401
    NEUTRAL_CATEGORIES,
    NEUTRAL_COLOR,
    SCALES,
    TABLEAU10,
    _aligned,
    _colour_by,
    _fitted_norm,
    _hex_colours,
    _norm,
    _scaled,
    category_colors,
)
from .crs import (  # noqa: F401
    FOOTPRINT_EXTENT,
    FOOTPRINT_PIPES,
    WGS84,
    _check_nodes,
    _geodata,
    _in_metres,
    _network,
    _with_crs,
    footprint_distance,
    lonlat_positions,
    network_extent,
    network_footprint,
    plot_crs,
)
from .export import (  # noqa: F401
    GroupLayer,
    MarkerLayer,
    PipeLayer,
    _add_groups,
    _add_markers,
    _add_pipes,
    _colour_legend,
    _value_colours,
    generate_gmap,
    to_geojson,
    to_leaflet,
)
from .polygons import (  # noqa: F401
    _path,
    _polygons,
    node_group_polygons,
)
from .ribbons import (  # noqa: F401
    MIN_RIBBON_WIDTH,
    OTHER,
    RIBBON_WIDTH,
    _direction,
    _half_plane,
    _joined,
    _merged,
    _shares_table,
    _stems,
    mitred_strands,
    plot_ribbons,
    ribbon_strands,
)
from .static import (  # noqa: F401
    CONTEXT_COLOR,
    CONTEXT_WIDTH,
    DEFAULT_HATCHES,
    plot_groups,
    plot_markers,
    plot_pipes,
)
from .style import (  # noqa: F401
    _AXES_STATE,
    DrawnLayer,
    _cleared,
    _finish,
    _nice_length,
    _renderer,
    _state,
    add_legends,
    map_axes,
)

logger = logging.getLogger(__name__)
