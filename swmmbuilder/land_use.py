"""Allocation of land-use flows (commercial, institutional, industrial, ...) to the nodes of a
network.

Design flows that are not per person are usually given per hectare of each land use, e.g. a
commercial allowance in L/ha/d. Given zoning polygons with a land use each:

1. each pipe goes to the polygon containing its midpoint or, with a distance rule, to every
   polygon within that distance of it, by the rules of :func:`~swmmbuilder.allocate_population`;
2. each polygon whose land use has a rate gets its area (optionally only the part within a given
   area, such as the network's footprint) times that rate;
3. the polygon's flow is shared among its pipes in proportion to a weight (by default the pipe
   length) and placed at their upstream nodes (or as asked);
4. optionally, every flow is scaled so that their sum is a given total, which makes the rates
   relative weights of the land uses.

Force mains take no flow unless asked. A rated polygon with no pipe contributes nothing, and is
reported, unless it is given its nearest pipe (``to_nearest``, for sites on private sewers that the
data do not show). The result's ``nodes`` is a table of nodes by land use, so that each land use can carry
its own time pattern (see :func:`~swmmbuilder.pattern_inflows`).

A pipe with a missing ``from_node`` or ``to_node`` is in no catchment, so it is left out (see :mod:`swmmbuilder.population`); the flow it would have taken is logged.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS

from ._summary import Summary
from .diagnostics import Problems
from .population import (
    PLACEMENTS,
    Placement,
    Weight,
    _assign,
    _check_polygon_values,
    _log_left_out,
    _missing_node,
    _nearest,
    _pipes_of,
    _pipes_within,
    _polygon_values,
    _quiet,
    _weights,
)
from .units import _FLOW_UNIT_LPS, _TIME_SECONDS, _VOLUME_LITRES, SI_FLOW_UNITS, _is_real

if TYPE_CHECKING:
    from .network import Network
    from .subnetwork import Subnetwork

logger = logging.getLogger(__name__)

M2_PER_HA = 1e4


@dataclass
class LandUseAllocation:
    """The result of :func:`allocate_land_use`.

    Attributes
    ----------
    nodes : pandas.DataFrame
        Flow at each node (rows, every node a share is placed at) from each land use with a rate
        (columns), in ``flow_units``.

    pipes : pandas.DataFrame
        One row per pipe and polygon it serves (one row per pipe with the midpoint rule), indexed
        by pipe_id, with ``from_node``, ``to_node``, ``length``, ``polygon`` (missing for a pipe
        that serves none), ``land_use``, ``weight``, and ``flow``.

    polygons : pandas.DataFrame
        One row per polygon whose land use has a rate, indexed by polygon ID, with ``land_use``,
        ``area_ha`` (the area counted), ``rate``, ``design_flow`` (area times rate, in
        ``flow_units``), ``pipes`` (the number of pipes serving it), and ``allocated`` (the flow
        placed on its pipes, after scaling).

    no_pipes : list
        IDs of the rated polygons with a flow but no pipe (or whose pipes all weigh 0): their
        flow is placed nowhere, and is not part of the scaling to ``total``.

    nearest : pandas.Series
        {polygon ID: distance} of the polygons that had no pipe and were given their nearest
        (with ``to_nearest``).

    scale : float
        The factor applied to every flow to reach ``total`` (1 without one).

    unallocated : float
        Design flow placed nowhere, in ``flow_units``, before scaling.

    placement : str
        The ``placement`` the shares were placed with.

    flow_units : str
        The unit of the flows.
    """

    nodes: pd.DataFrame
    pipes: pd.DataFrame
    polygons: pd.DataFrame
    no_pipes: list = field(default_factory=list)
    nearest: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    scale: float = 1.0
    unallocated: float = 0.0
    placement: Placement = "upstream"
    flow_units: str = "LPS"

    @property
    def total(self) -> float:
        """Flow placed at the nodes, all land uses."""
        return float(self.nodes.to_numpy().sum())

    def by_pipe(self) -> pd.Series:
        """Flow of each pipe, summed over the polygons it serves, indexed by pipe_id, in
        :attr:`flow_units` (see :meth:`PopulationAllocation.by_pipe
        <swmmbuilder.PopulationAllocation.by_pipe>`)."""
        return self.pipes.groupby(level=0, sort=False)["flow"].sum()

    def within(self, sub: Subnetwork) -> pd.DataFrame:
        """Flow at the nodes of ``sub`` by land use, from the shares of the subnetwork's own pipes:
        :attr:`nodes` cut to the subnetwork, of the same type (see
        :meth:`PopulationAllocation.within <swmmbuilder.PopulationAllocation.within>`)."""
        table = _pipes_within(self.pipes, sub)
        return _place_by_use(table, self.placement, list(self.nodes.columns))

    def summary(self, sub: Subnetwork | None = None) -> dict[str, float | int]:
        """Counts and totals of the allocation, to cite in a report.

        The keys follow the rules of every ``summary()`` (see the Diagnostics page of the
        documentation), and match those of :meth:`PopulationAllocation.summary
        <swmmbuilder.PopulationAllocation.summary>` where they measure the same thing.

        Parameters
        ----------
        sub : Subnetwork, optional
            A subnetwork (e.g. a catchment) of the network the flows were allocated over. The
            keys marked "within ``sub``" count only its pipes and the nodes they place flow at
            (see :meth:`within`); without it, they count every pipe and node. Scaling to a
            ``total`` is over the whole allocation, so ``total`` within ``sub`` is the
            subnetwork's share of it.

        Returns
        -------
        dict
            Plain Python numbers (``int`` for counts, ``float`` for totals, NaN where there is
            nothing to measure). For each land use with a rate, under its name as given (e.g.
            ``COMMERCIAL``):

            - ``<use>_polygons``: its polygons served by at least one pipe.
            - ``<use>_area_ha``: their area counted, in ha.
            - ``<use>_flow``: its flow placed at the nodes, in :attr:`flow_units`, within
              ``sub``.
            - ``<use>_nodes``: nodes given a flow above 0 from it, within ``sub``.

            And over all land uses:

            - ``total``: flow placed at the nodes, within ``sub``.
            - ``nodes``: nodes given a flow above 0, within ``sub``.
            - ``weight``: total weight of the pipes, within ``sub`` (each pipe counted once).
            - ``allocation_total``: flow placed at the nodes over the whole allocation
              (:attr:`total`).
            - ``polygons``: polygons with a rate.
            - ``polygons_served``: those with at least one pipe.
            - ``area_ha``: area counted of the polygons served, in ha.
            - ``nearest``: polygons given their nearest pipe (:attr:`nearest`).
            - ``nearest_median_distance``, ``nearest_max_distance``: the median and the
              farthest of their distances to it, in m (NaN if none).
            - ``no_pipes``: polygons with a flow and no pipe.
            - ``zero_weight``: polygons with a flow whose pipes all weigh 0 (counted in
              ``polygons_served``). These and ``no_pipes`` make up :attr:`no_pipes`.
            - ``unallocated_design_flow``: the design flow of both, placed nowhere and not
              scaled (:attr:`unallocated`), in :attr:`flow_units`.
            - ``scale``: the factor applied to every flow (:attr:`scale`).

            ``unallocated``, the former name of ``unallocated_design_flow``, still reads, with
            a ``DeprecationWarning``.

        Examples
        --------
        ::

            s = land.summary(net.catchment("OUT"))
            print(f"commercial: {s['COMMERCIAL_flow']:.1f} L/s at {s['COMMERCIAL_nodes']} nodes")
        """
        pipes = self.pipes
        nodes = self.nodes
        if sub is not None:
            pipes = _pipes_within(pipes, sub)
            nodes = _place_by_use(pipes, self.placement, list(self.nodes.columns))
        served = self.polygons[self.polygons["pipes"] > 0]
        lost = self.polygons.index.isin(list(self.no_pipes))
        result: dict[str, float | int] = {}
        for use in self.nodes.columns:
            mine = served["land_use"] == use
            result[f"{use}_polygons"] = int(mine.sum())
            result[f"{use}_area_ha"] = float(served.loc[mine, "area_ha"].sum())
            result[f"{use}_flow"] = float(nodes[use].sum())
            result[f"{use}_nodes"] = int((nodes[use] > 0).sum())
        nearest = self.nearest
        result.update(
            total=float(nodes.to_numpy().sum()),
            nodes=int((nodes.sum(axis=1) > 0).sum()),
            weight=float(pipes.groupby(level=0)["weight"].mean().sum()),
            allocation_total=self.total,
            polygons=len(self.polygons),
            polygons_served=len(served),
            area_ha=float(served["area_ha"].sum()),
            nearest=len(nearest),
            nearest_median_distance=float(nearest.median()) if len(nearest) else math.nan,
            nearest_max_distance=float(nearest.max()) if len(nearest) else math.nan,
            no_pipes=int((lost & (self.polygons["pipes"] == 0)).sum()),
            zero_weight=int((lost & (self.polygons["pipes"] > 0)).sum()),
            unallocated_design_flow=float(self.unallocated),
            scale=float(self.scale),
        )
        return Summary(result, {"unallocated": "unallocated_design_flow"})


def _place_by_use(table, placement, uses):
    """Nodes x land-use flows of the pipes of ``table`` (the ``pipes`` of a LandUseAllocation)."""
    table = table[table["land_use"].isin(uses) & (table["flow"] > 0)]
    if placement == "split":
        half = table.assign(flow=table["flow"] / 2)
        parts = [
            half.rename(columns={"from_node": "node"}),
            half.rename(columns={"to_node": "node"}),
        ]
        long = pd.concat([p[["node", "land_use", "flow"]] for p in parts])
    else:
        end = "from_node" if placement == "upstream" else "to_node"
        long = table.rename(columns={end: "node"})[["node", "land_use", "flow"]]
    wide = long.pivot_table(
        index="node", columns="land_use", values="flow", aggfunc="sum", fill_value=0.0
    )
    wide = wide.reindex(columns=uses, fill_value=0.0)
    wide.columns.name = None
    return wide.rename_axis("node").astype(float)


def _per_second(unit):
    """Litres per second per hectare in one ``unit`` (e.g. 'L/ha/d'); raise ValueError if the
    unit is not volume/ha/time."""
    parts = unit.split("/")
    if len(parts) != 3 or parts[1] != "ha" or parts[0] not in _VOLUME_LITRES:
        parts = None
    if parts is None or parts[2] not in _TIME_SECONDS:
        raise ValueError(
            f"unit must be one of 'L/ha/d', 'm3/ha/d', 'L/ha/s', or 'm3/ha/s', not {unit!r}."
        )
    return _VOLUME_LITRES[parts[0]] / _TIME_SECONDS[parts[2]]


def _check(
    pipes,
    polygons,
    land_use,
    rates,
    total,
    area_within,
    max_distance,
    placement,
    weight,
    max_nearest_distance=None,
):
    problems = Problems("Invalid land-use allocation")
    if max_nearest_distance is not None and not max_nearest_distance >= 0:
        problems.add(
            f"max_nearest_distance must be a distance >= 0 or None, not {max_nearest_distance!r}"
        )
    _check_polygon_values(polygons, land_use, "land_use", problems)
    if not isinstance(rates, Mapping) or not rates:
        problems.add("rates must be a non-empty mapping {land use: rate}")
    else:
        bad = [f"{u} ({r!r})" for u, r in rates.items() if not (_is_real(r) and r >= 0)]
        if bad:
            problems.add(f"rates must be non-negative numbers: {', '.join(map(str, bad))}")
    if total is not None and not (_is_real(total) and total >= 0):
        problems.add(f"total must be a non-negative number or None, not {total!r}")
    if max_distance is not None and not max_distance >= 0:
        problems.add(f"max_distance must be a distance >= 0 or None, not {max_distance!r}")
    if placement not in PLACEMENTS:
        problems.add(f"placement must be one of {', '.join(PLACEMENTS)}, not {placement!r}")
    if isinstance(weight, str) and weight not in pipes.columns:
        problems.add(f"weight: the pipes have no column {weight!r}")
    repeated = polygons.index[polygons.index.duplicated()].unique()
    if len(repeated):
        problems.add(f"polygon IDs (the index of polygons) repeated: {list(repeated)[:20]}")
    for name, layer in (("the network", pipes), ("polygons", polygons)):
        if layer.crs is None:
            problems.add(f"{name} has no CRS")
    if pipes.crs is not None:
        crs = CRS.from_user_input(pipes.crs)
        if not (crs.is_projected and math.isclose(crs.axis_info[0].unit_conversion_factor, 1.0)):
            problems.add(f"the network's CRS must be projected in metres, not {crs.name}")
    if area_within is not None and not hasattr(area_within, "area"):
        problems.add("area_within must be a shapely geometry in the network's CRS, or None")
    problems.raise_if_any()


def allocate_land_use(
    network: Network | Subnetwork,
    polygons: gpd.GeoDataFrame,
    land_use: str | pd.Series,
    rates: Mapping[Any, float],
    unit: str = "L/ha/d",
    total: float | None = None,
    area_within: Any = None,
    max_distance: float | None = None,
    to_nearest: bool = False,
    max_nearest_distance: float | None = None,
    serve_force_mains: bool = False,
    weight: Weight = None,
    placement: Placement = "upstream",
    flow_units: str = "LPS",
) -> LandUseAllocation:
    """Allocate the flow of each land use, given per hectare, to the nodes of ``network``.

    Each pipe is assigned to the polygon containing its midpoint (or to every polygon within
    ``max_distance`` of it); each polygon whose land use is in ``rates`` gets its area times the
    rate, shared among its pipes in proportion to ``weight`` and placed at their upstream nodes
    (or as ``placement`` says). See :mod:`swmmbuilder.land_use`.

    Parameters
    ----------
    network : Network or Subnetwork
        The pipes to allocate to: every pipe of a :class:`~swmmbuilder.Network` (but those with a
        missing node), or those of a :class:`~swmmbuilder.Subnetwork`. Its CRS must be projected
        in metres.

    polygons : geopandas.GeoDataFrame
        Zoning polygons, indexed by a unique ID, with a CRS (they are reprojected to the
        network's).

    land_use : str or pandas.Series
        Land use of each polygon: a column of ``polygons``, or a Series indexed by polygon ID
        (all IDs must be present and unique). Unindexed lists/arrays and missing or
        repeated IDs raise ``ValueError``.

    rates : dict
        {land use: flow per hectare, in ``unit``}, e.g. ``{'COMMERCIAL': 23000,
        'INDUSTRIAL': 20000}`` in L/ha/d. Polygons of any other land use take no flow, but their
        pipes are still assigned to them, so a pipe in a residential polygon takes nothing.

    unit : {'L/ha/d', 'm3/ha/d', 'L/ha/s', 'm3/ha/s'}
        Unit of the rates (default L/ha/d).

    total : float, optional
        Scale every flow so that the flows placed sum to this, in ``flow_units``; the rates are
        then only relative weights of the land uses. None (default) keeps the rates as given.

    area_within : shapely geometry, optional
        Count only the part of each polygon within this area, in the network's CRS, e.g. the
        network's footprint (:func:`swmmbuilder.maps.network_footprint`), so that a large polygon
        mostly away from the sewers does not get the flow of its whole area. None (default)
        counts the whole polygon.

    max_distance : float, optional
        As in :func:`~swmmbuilder.allocate_population`: None (default) for the midpoint rule, or
        a distance (m) assigning each pipe to every polygon within it.

    to_nearest : bool
        Give a polygon with a flow but no pipe by the midpoint or distance rule its nearest pipe,
        which takes the polygon's whole flow (default False: its flow is placed nowhere). Zoning
        polygons such as commercial sites are often served by private sewers that the City's
        data do not show, so their flow reaches the nearest public sewer.

    max_nearest_distance : float, optional
        With ``to_nearest``, the farthest a polygon may be from its nearest pipe, in m (default
        None: no limit). A polygon farther than this keeps no pipe, and is reported in
        ``no_pipes``.

    serve_force_mains : bool
        Let force mains take flow (default False).

    weight : str, callable or None
        Weight of each pipe in its polygon's share: a column of the pipes, a function of the pipes
        returning one value per pipe, or None (default) for the pipe length. Missing, non-numeric,
        or negative weights count as 0.

    placement : {'upstream', 'downstream', 'split'}
        Where a pipe's share is placed: its upstream node (default), its downstream node, or half
        at each.

    flow_units : {'LPS', 'CMS', 'MLD'}
        Unit of the flows returned, and of ``total`` (default LPS, the shipped template's).

    Returns
    -------
    LandUseAllocation
        The flow at each node by land use, the tables per pipe and per polygon, and the polygons
        whose flow is placed nowhere.

    Raises
    ------
    ValueError
        If an option is invalid, a polygon ID is repeated, a layer has no CRS, the network's CRS
        is not projected in metres, or ``total`` is positive but no rated polygon has a pipe;
        the message lists every problem.

    Warns
    -----
    DataWarning
        With the midpoint rule, if a pipe midpoint is in more than one polygon (it takes the
        first).

    Examples
    --------
    London's non-residential design allowances as relative weights, over the part of each zone
    within the catchment's footprint, scaled to a calibrated total of 120 L/s::

        sub = net.catchment("OUT")
        flows = allocate_land_use(
            sub, zoning, land_use="GeneralizedLandUse",
            rates={"COMMERCIAL": 23000, "INSTITUTIONAL": 23000, "INDUSTRIAL": 20000},
            total=120.0, area_within=maps.network_footprint(sub),
        )
        flows.nodes  # nodes x land uses, L/s
    """
    pipes = _pipes_of(network)
    options = dict(
        polygons=polygons,
        land_use=land_use,
        rates=rates,
        unit=unit,
        total=total,
        area_within=area_within,
        max_distance=max_distance,
        to_nearest=to_nearest,
        max_nearest_distance=max_nearest_distance,
        serve_force_mains=serve_force_mains,
        weight=weight,
        placement=placement,
        flow_units=flow_units,
    )
    missing = _missing_node(pipes)
    result = _allocate_land_use(pipes[~missing], **options)
    if missing.any():
        with _quiet(logger):
            kept = _allocate_land_use(pipes, **options)
        taken = float(kept.pipes.loc[kept.pipes.index.isin(pipes.index[missing]), "flow"].sum())
        _log_left_out(
            logger,
            pipes,
            missing,
            f"with them they would take {taken:.6g} {result.flow_units}, in no catchment",
        )
    return result


def _allocate_land_use(
    pipes,
    polygons,
    land_use,
    rates,
    unit,
    total,
    area_within,
    max_distance,
    to_nearest,
    max_nearest_distance,
    serve_force_mains,
    weight,
    placement,
    flow_units,
) -> LandUseAllocation:
    """:func:`allocate_land_use` over ``pipes`` (indexed by pipe_id)."""
    _check(
        pipes,
        polygons,
        land_use,
        rates,
        total,
        area_within,
        max_distance,
        placement,
        weight,
        max_nearest_distance,
    )
    per_second = _per_second(unit)
    flow_units = str(flow_units).upper()
    if flow_units not in _FLOW_UNIT_LPS:
        raise ValueError(
            f"flow_units must be one of {', '.join(SI_FLOW_UNITS)}, not {flow_units!r}."
        )
    to_units = per_second / _FLOW_UNIT_LPS[flow_units]
    uses = _polygon_values(polygons, land_use)
    polygons = polygons.to_crs(pipes.crs)

    if serve_force_mains or "structure_type" not in pipes.columns:
        serving = np.ones(len(pipes), dtype=bool)
    else:
        serving = (pipes["structure_type"] != "force_main").to_numpy()
    rows = _assign(pipes, polygons, max_distance, serving)

    rated = uses.isin(list(rates))
    geometry = polygons.geometry[rated]
    if area_within is not None:
        geometry = geometry.intersection(area_within)
    summary = pd.DataFrame({"land_use": uses[rated]})
    summary["area_ha"] = geometry.area / M2_PER_HA
    summary["rate"] = summary["land_use"].map(rates).astype(float)
    summary["design_flow"] = summary["area_ha"] * summary["rate"] * to_units
    nearest = pd.Series(dtype=float)
    if to_nearest:
        assigned = summary.index.isin(rows["polygon"].dropna())
        lacking = summary.index[(summary["design_flow"] > 0) & ~assigned]
        if len(lacking):
            found = _nearest(pipes, polygons.loc[lacking], serving)
            if max_nearest_distance is not None:
                too_far = found["distance"] > max_nearest_distance
                if too_far.any():
                    logger.info(
                        "%d polygons with a land-use rate are more than %g from any pipe; "
                        "not allocated",
                        int(too_far.sum()),
                        max_nearest_distance,
                    )
                found = found[~too_far]
            nearest = found["distance"]
            added = pd.DataFrame({"pipe": found["pipe"].to_numpy(), "polygon": found.index})
            # A pipe that served no polygon loses its empty row.
            empty = rows["polygon"].isna() & rows["pipe"].isin(added["pipe"])
            rows = pd.concat([rows[~empty], added], ignore_index=True)
            rows = rows.sort_values("pipe", kind="stable", ignore_index=True)
            if len(found):
                logger.info(
                    "%d polygons with a land-use rate had no pipe; given their nearest (median "
                    "%.0f, max %.0f away)",
                    len(found),
                    float(nearest.median()),
                    float(nearest.max()),
                )
    position = rows["pipe"].to_numpy()
    table = pipes[["pipe_id", "from_node", "to_node", "length"]].iloc[position]
    table = table.reset_index(drop=True)
    table["polygon"] = rows["polygon"].to_numpy()
    table["land_use"] = table["polygon"].map(uses)

    given = pd.Series(_weights(pipes, weight).to_numpy()[position])
    invalid = given.isna() | (given < 0)
    if invalid.any():
        logger.info("%d pipe weights missing or negative; counted as 0", int(invalid.sum()))
    table["weight"] = given.where(~invalid, 0.0)
    table.loc[~serving[position], "weight"] = 0.0

    in_rated = table["polygon"].isin(summary.index)
    weights = table[in_rated].groupby("polygon")["weight"].sum()
    counts = table[in_rated].groupby("polygon").size()
    summary["pipes"] = counts.reindex(summary.index, fill_value=0).astype(int)
    weights = weights.reindex(summary.index, fill_value=0.0)
    per_weight = (summary["design_flow"] / weights).where(weights > 0, 0.0)
    table["flow"] = (table["weight"] * table["polygon"].map(per_weight)).fillna(0.0)

    lost = (summary["design_flow"] > 0) & (weights == 0)
    no_pipes = list(summary.index[lost])
    unallocated = float(summary.loc[lost, "design_flow"].sum())
    placed = float(table["flow"].sum())
    scale = 1.0
    if total is not None:
        if placed == 0 and total > 0:
            raise ValueError(
                f"total is {total} {flow_units}, but no polygon with a rate has a pipe to place "
                "it on."
            )
        scale = total / placed if placed > 0 else 1.0
        table["flow"] *= scale
    allocated = table[in_rated].groupby("polygon")["flow"].sum()
    summary["allocated"] = allocated.reindex(summary.index, fill_value=0.0)

    table = table.set_index("pipe_id")
    nodes = _place_by_use(table, placement, list(rates))
    if no_pipes:
        logger.warning(
            "%d polygons with a land-use rate have no pipe; their %.6g %s is not allocated",
            len(no_pipes),
            unallocated,
            flow_units,
        )
    logger.info(
        "land-use flow: %.6g %s placed at %d nodes from %d polygons%s",
        float(nodes.to_numpy().sum()),
        flow_units,
        len(nodes),
        int((summary["allocated"] > 0).sum()),
        f" (scaled by {scale:.6g} to the total)" if total is not None else "",
    )
    return LandUseAllocation(
        nodes=nodes,
        pipes=table,
        polygons=summary,
        no_pipes=no_pipes,
        nearest=nearest,
        scale=scale,
        unallocated=unallocated,
        placement=placement,
        flow_units=flow_units,
    )


def dominant_category(
    polygons: gpd.GeoDataFrame, layer: gpd.GeoDataFrame, column: str
) -> pd.DataFrame:
    """The category of a layer (such as zoning) covering the largest area of each polygon.

    A typical use is the ``sewered`` flag of :func:`~swmmbuilder.allocate_population`, true for
    the census blocks that are mostly zoned residential. The layer is reprojected to the
    polygons' CRS and cut by them (:func:`geopandas.overlay`, keeping polygon parts only), and
    the area of each category is summed within each polygon.

    Parameters
    ----------
    polygons : geopandas.GeoDataFrame
        The polygons, indexed by a unique ID, with a CRS projected in metres.
    layer : geopandas.GeoDataFrame
        The polygons carrying the categories, with a CRS. Features with a missing category are
        ignored.
    column : str
        The column of ``layer`` holding the category, e.g. ``"GeneralizedLandUse"``.

    Returns
    -------
    pandas.DataFrame
        Indexed like ``polygons``, with ``category`` (the category covering the largest area of
        the polygon, missing where the polygon meets no feature of the layer) and ``share`` (that
        area over the polygon's area, from 0 to 1; 0 where the category is missing). A tie goes to
        the category that sorts first. Where features of the layer overlap, the area of each is
        counted, so a category's area can exceed what it covers; ``share`` is capped at 1.

    Raises
    ------
    ValueError
        If ``layer`` has no column ``column``, a polygon ID is repeated, a layer has no CRS, or
        the polygons' CRS is not projected in metres; the message lists every problem.

    Examples
    --------
    Census blocks flagged where most of the area is zoned residential::

        use = dominant_category(blocks, zoning, "GeneralizedLandUse")
        blocks["residential"] = use["category"] == "RESIDENTIAL"
        alloc = allocate_population(net, blocks, population="population", sewered="residential")
    """
    problems = Problems("Invalid dominant category")
    if column not in layer.columns:
        problems.add(f"column: the layer has no column {column!r}")
    repeated = polygons.index[polygons.index.duplicated()].unique()
    if len(repeated):
        problems.add(f"polygon IDs (the index of polygons) repeated: {list(repeated)[:20]}")
    for name, frame in (("polygons", polygons), ("the layer", layer)):
        if frame.crs is None:
            problems.add(f"{name} has no CRS")
    if polygons.crs is not None:
        crs = CRS.from_user_input(polygons.crs)
        if not (crs.is_projected and math.isclose(crs.axis_info[0].unit_conversion_factor, 1.0)):
            problems.add(f"the polygons' CRS must be projected in metres, not {crs.name}")
    problems.raise_if_any()

    result = pd.DataFrame(
        {"category": pd.Series(np.nan, index=polygons.index, dtype=object), "share": 0.0}
    )
    features = layer[layer[column].notna()]
    if polygons.empty or features.empty:
        return result
    # Positions stand in for the polygon IDs, which may clash with the layer's columns.
    left = gpd.GeoDataFrame(
        {"_polygon": np.arange(len(polygons))}, geometry=polygons.geometry.values, crs=polygons.crs
    )
    right = gpd.GeoDataFrame(
        {"_category": features[column].to_numpy()},
        geometry=features.geometry.values,
        crs=features.crs,
    ).to_crs(polygons.crs)
    parts = gpd.overlay(left, right, how="intersection", keep_geom_type=True)
    if parts.empty:
        return result
    areas = (
        parts.assign(_area=parts.area)
        .groupby(["_polygon", "_category"], sort=True)["_area"]
        .sum()
        .reset_index()
    )
    # Largest area first; a stable sort keeps the categories in sorted order within a tie.
    areas = areas.sort_values("_area", ascending=False, kind="stable")
    top = areas.drop_duplicates("_polygon").set_index("_polygon")
    position = top.index.to_numpy()
    own_area = polygons.geometry.area.to_numpy()[position]
    with np.errstate(divide="ignore", invalid="ignore"):
        share = np.minimum(top["_area"].to_numpy() / own_area, 1.0)
    result.iloc[position, 0] = top["_category"].to_numpy()
    result.iloc[position, 1] = np.nan_to_num(share, nan=0.0)
    return result
