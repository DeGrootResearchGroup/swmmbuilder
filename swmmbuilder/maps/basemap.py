"""Web-tile basemaps under a static plot."""

from __future__ import annotations

from typing import TYPE_CHECKING

import geopandas as gpd

from ._deps import _optional
from .crs import plot_crs

if TYPE_CHECKING:
    from typing import Any

    from matplotlib.axes import Axes

    from ..network import Network
    from ..subnetwork import Subnetwork


def add_basemap(
    ax: Axes,
    net: Network | Subnetwork,
    source=None,
    *,
    layer=None,
    attribution=None,
    color="0.75",
    linewidth=0.6,
    crs=None,
    zorder=0,
    **kwargs,
) -> Any:
    """Draw a basemap beneath what is on ``ax``: web tiles, or a vector layer for offline and
    publication use. The extent of the axes is kept.

    Call it after plotting the network, so that the tiles cover the axes' extent.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes drawn by the plots of ``net``.

    net : Network or Subnetwork
        Gives the CRS of the axes (with ``crs``; see :func:`plot_crs`).

    source : optional
        Tile provider for ``contextily`` (an ``xyzservices`` provider or a URL); default
        CartoDB Positron. The provider's attribution is written on the plot unless
        ``attribution`` is given (False for none). Needs the ``[maps]`` extra.

    layer : geopandas.GeoDataFrame, GeoSeries or path, optional
        A vector layer (roads, rivers, ...) drawn instead of tiles, in ``color`` and
        ``linewidth``. It must have a CRS.

    **kwargs
        Passed to :func:`contextily.add_basemap` or to the layer's ``plot``.
    """
    crs = plot_crs(net, crs)
    if layer is not None:
        if not isinstance(layer, gpd.GeoDataFrame | gpd.GeoSeries):
            layer = gpd.read_file(layer)
        if layer.crs is None:
            raise ValueError("The basemap layer has no CRS.")
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        layer.to_crs(crs).plot(ax=ax, color=color, linewidth=linewidth, zorder=zorder, **kwargs)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        return
    cx = _optional("contextily")
    if source is None:
        source = cx.providers.CartoDB.Positron
    cx.add_basemap(
        ax, crs=crs.to_string(), source=source, attribution=attribution, zorder=zorder, **kwargs
    )
