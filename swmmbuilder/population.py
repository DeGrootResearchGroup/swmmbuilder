"""Allocation of the population of areal units (polygons) to the nodes of a network.

The population of each polygon (a census tract or dissemination area, a postal code area, a
traffic zone, or any other areal unit) is placed on the pipes in it and then on their nodes:

1. each pipe goes to the polygon containing its midpoint or, with a distance rule, to every
   polygon within that distance of it (for small areal units such as city blocks, whose edges
   follow the streets the sewers run under, so that the pipe serves the blocks on both sides);
2. each polygon's population is shared among its pipes in proportion to a weight per pipe (by
   default the pipe length; e.g. a count of service connections), with invalid or missing weights
   imputed from the polygon's median weight per metre of pipe, and a fallback weight (one
   connection, say) where that median is 0 or the polygon's pipes would weigh nothing;
3. each pipe's share is placed at its upstream node, its downstream node, or half at each.

The same polygons serve for the share and for the imputation. Force mains, which take no service
connections, serve no polygon unless asked. A polygon with no pipe contributes no one, unless it is
flagged as sewered (e.g. by its zoning), when its people go to the nearest pipe; a polygon only
partly sewered places all its people on the pipes it has. The population left out is reported.

The result's ``nodes`` ({node: population}) is the input of
:func:`~swmmbuilder.swmm.dwf.dwf_from_population`.
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import geopandas as gpd
import numpy as np
import pandas as pd

from .diagnostics import DataWarning, Problems

if TYPE_CHECKING:
    from ._types import NodeId
    from .network import Network
    from .subnetwork import Subnetwork

logger = logging.getLogger(__name__)

Weight = str | Callable[[gpd.GeoDataFrame], pd.Series] | None
Placement = Literal["upstream", "downstream", "split"]
IMPUTATIONS = ("per_length", "zero")
PLACEMENTS = ("upstream", "downstream", "split")


@dataclass
class PopulationAllocation:
    """The result of :func:`allocate_population`.

    Attributes
    ----------
    nodes : dict
        {node: population} of every node a pipe's share is placed at (0 for a node whose pipes
        carry no one), to pass to :func:`~swmmbuilder.swmm.dwf.dwf_from_population`.

    pipes : pandas.DataFrame
        One row per pipe and polygon it serves (one row per pipe with the midpoint rule), indexed
        by pipe_id, with ``from_node``, ``to_node``, ``length``, ``polygon`` (missing for a pipe
        that serves none),
        ``weight_imputed`` (bool), ``weight_fallback`` (bool: the weight is ``fallback_weight``),
        ``weight`` (after imputation), and ``population``.

    polygons : pandas.DataFrame
        One row per polygon, indexed by polygon ID, with ``population`` (as given, a missing value
        as 0), ``pipes`` (the number of pipes serving it), ``weight`` (their total weight), and
        ``allocated`` (the population placed on its pipes).

    unallocated : float
        Population placed at no node: that of the polygons of ``no_pipes``, ``no_valid_weight``,
        and ``zero_weight``, and the shares of pipes with no node at the end they are placed at.

    no_pipes : list
        IDs of the polygons with people but no pipe (and not flagged ``sewered``).

    nearest : pandas.Series
        {polygon ID: distance} of the polygons flagged ``sewered`` that had no pipe and were given
        their nearest pipe.

    no_valid_weight : list
        IDs of the polygons with people and pipes, none of whose weights is valid, left
        unallocated because there is no ``fallback_weight`` (with ``impute='per_length'``).

    zero_weight : list
        IDs of the other polygons with people and pipes, whose pipes all weigh 0, left
        unallocated because there is no ``fallback_weight``.

    fallback : list
        IDs of the polygons some of whose pipes were given ``fallback_weight``.

    placement : str
        The ``placement`` the shares were placed with.
    """

    nodes: dict[NodeId, float]
    pipes: pd.DataFrame
    polygons: pd.DataFrame
    unallocated: float
    no_pipes: list = field(default_factory=list)
    nearest: pd.Series = field(default_factory=lambda: pd.Series(dtype=float))
    no_valid_weight: list = field(default_factory=list)
    zero_weight: list = field(default_factory=list)
    fallback: list = field(default_factory=list)
    placement: Placement = "upstream"

    @property
    def total(self) -> float:
        """Population placed at the nodes."""
        return float(sum(self.nodes.values()))

    def by_pipe(self) -> pd.Series:
        """Population of each pipe, summed over the polygons it serves, indexed by pipe_id."""
        return self.pipes.groupby(level=0, sort=False)["population"].sum()

    def within(self, sub: Subnetwork) -> dict[NodeId, float]:
        """{node: population} placed by the pipes of ``sub``, e.g. a catchment of the network the
        population was allocated over.

        Only the shares of the subnetwork's own pipes count: a pipe leaving the subnetwork from
        one of its nodes (at a flow split, say) takes its people out of it.
        """
        keys = {key for _, _, key in sub.graph.edges(keys=True)}
        return _place(self.pipes[self.pipes.index.isin(keys)], self.placement)


def _place(table, placement):
    """{node: population} of the pipes of ``table`` (the ``pipes`` of a PopulationAllocation)."""
    if placement == "split":
        half = table["population"] / 2
        placed = pd.concat(
            [half.groupby(table["from_node"]).sum(), half.groupby(table["to_node"]).sum()]
        )
        placed = placed.groupby(level=0).sum()
    else:
        end = table["from_node" if placement == "upstream" else "to_node"]
        placed = table["population"].groupby(end).sum()
    return {n: float(p) for n, p in placed.items()}


def _polygon_of(points: gpd.GeoSeries, polygons: gpd.GeoDataFrame) -> pd.Series:
    """ID (index) of the polygon of ``polygons`` containing each point, missing if none. A point
    in more than one polygon takes the first of them, in the order of ``polygons``."""
    right = gpd.GeoDataFrame(
        {"_id": polygons.index, "_order": np.arange(len(polygons))},
        geometry=polygons.geometry.to_numpy(),
        crs=polygons.crs,
    )
    left = gpd.GeoDataFrame(geometry=points.to_numpy(), crs=points.crs)
    found = left.sjoin(right, how="inner", predicate="within").sort_values("_order", kind="stable")
    overlaps = int(found.index.duplicated().sum())
    if overlaps:
        warnings.warn(
            f"{overlaps} pipe midpoints are in more than one polygon; each takes the first",
            DataWarning,
            stacklevel=3,
        )
    found = found[~found.index.duplicated()]
    ids = found["_id"].reindex(range(len(points)))
    return pd.Series(ids.to_numpy(), index=points.index)


def _near(pipes, polygons, distance):
    """Rows (pipe: position in ``pipes``, polygon: ID) of each pipe within ``distance`` of each
    polygon, in the order of the pipes and then of the polygons; a pipe near no polygon (or with
    no geometry) has one row with no polygon."""
    right = gpd.GeoDataFrame(
        {"_id": polygons.index, "_order": np.arange(len(polygons))},
        geometry=polygons.geometry.buffer(distance).to_numpy(),
        crs=polygons.crs,
    )
    left = gpd.GeoDataFrame(
        {"pipe": np.arange(len(pipes))}, geometry=pipes.geometry.to_numpy(), crs=pipes.crs
    )
    found = left.sjoin(right, how="left", predicate="intersects")
    found = found.sort_values(["pipe", "_order"], kind="stable")
    return pd.DataFrame({"pipe": found["pipe"].to_numpy(), "polygon": found["_id"].to_numpy()})


def _assign(pipes, polygons, max_distance, serving):
    """Rows (pipe: position in ``pipes``, polygon: ID) assigning the ``serving`` pipes to
    polygons by the midpoint or distance rule; every other pipe, and a serving pipe with no
    polygon, has one row with no polygon. Rows are in the order of the pipes."""
    index = np.flatnonzero(serving)
    subset = pipes.iloc[index]
    if max_distance is None:
        midpoints = subset.geometry.interpolate(0.5, normalized=True)
        polygon = _polygon_of(midpoints, polygons).to_numpy()
        rows = pd.DataFrame({"pipe": np.arange(len(subset)), "polygon": polygon})
    else:
        rows = _near(subset, polygons, max_distance)
    rows["pipe"] = index[rows["pipe"].to_numpy()]
    rest = pd.DataFrame({"pipe": np.flatnonzero(~serving), "polygon": None})
    rows = pd.concat([rows, rest], ignore_index=True)
    return rows.sort_values("pipe", kind="stable", ignore_index=True)


def _nearest(pipes, polygons, serving):
    """(pipe position, distance) of the nearest ``serving`` pipe to each polygon of
    ``polygons``, as a DataFrame indexed like ``polygons`` (the first pipe on a tie)."""
    index = np.flatnonzero(serving & pipes.geometry.notna().to_numpy())
    lines = gpd.GeoDataFrame(
        {"pipe": index}, geometry=pipes.geometry.iloc[index].to_numpy(), crs=pipes.crs
    )
    found = gpd.GeoDataFrame(geometry=polygons.geometry.to_numpy(), crs=polygons.crs).sjoin_nearest(
        lines, distance_col="distance"
    )
    found = found.sort_values("pipe", kind="stable")
    found = found[~found.index.duplicated()].sort_index()
    return pd.DataFrame(
        {"pipe": found["pipe"].to_numpy(), "distance": found["distance"].to_numpy()},
        index=polygons.index[found.index],
    )


def _pipes_of(network):
    """The pipes of a Network, or of a Subnetwork, indexed by pipe_id."""
    from . import subnetwork

    data = network.pipes() if isinstance(network, subnetwork.Subnetwork) else network.full_data
    return data.set_index("pipe_id", drop=False)


def _weights(pipes, weight):
    """The weight of each pipe as given, as numbers (anything else missing)."""
    if weight is None:
        values = pipes["length"]
    elif callable(weight):
        values = pd.Series(weight(pipes), index=pipes.index)
    else:
        values = pipes[weight]
    return pd.to_numeric(values, errors="coerce").astype(float)


def _check(pipes, polygons, population, max_distance, sewered, weight, impute, placement):
    """Raise ValueError listing every invalid option."""
    problems = Problems("Invalid population allocation")
    if isinstance(sewered, str) and sewered not in polygons.columns:
        problems.add(f"sewered: polygons has no column {sewered!r}")
    if max_distance is not None and not max_distance >= 0:
        problems.add(f"max_distance must be a distance >= 0 or None, not {max_distance!r}")
    if impute not in IMPUTATIONS:
        problems.add(f"impute must be one of {', '.join(IMPUTATIONS)}, not {impute!r}")
    if placement not in PLACEMENTS:
        problems.add(f"placement must be one of {', '.join(PLACEMENTS)}, not {placement!r}")
    if isinstance(weight, str) and weight not in pipes.columns:
        problems.add(f"weight: the pipes have no column {weight!r}")
    if isinstance(population, str) and population not in polygons.columns:
        problems.add(f"population: polygons has no column {population!r}")
    repeated = polygons.index[polygons.index.duplicated()].unique()
    if len(repeated):
        problems.add(f"polygon IDs (the index of polygons) repeated: {list(repeated)[:20]}")
    for name, layer in (("the network", pipes), ("polygons", polygons)):
        if layer.crs is None:
            problems.add(f"{name} has no CRS")
    problems.raise_if_any()


def allocate_population(
    network: Network | Subnetwork,
    polygons: gpd.GeoDataFrame,
    population: str | pd.Series = "population",
    max_distance: float | None = None,
    sewered: str | pd.Series | None = None,
    serve_force_mains: bool = False,
    weight: Weight = None,
    max_weight: float | None = None,
    impute: Literal["per_length", "zero"] = "per_length",
    fallback_weight: float | None = 1.0,
    placement: Placement = "upstream",
) -> PopulationAllocation:
    """Allocate the population of each polygon to the nodes of ``network``.

    Each pipe is assigned to the polygon containing its midpoint (or to every polygon within
    ``max_distance`` of it), each polygon's population is shared among its pipes in proportion to
    ``weight``, and each pipe's share is placed at its
    upstream node (or as ``placement`` says). See :mod:`swmmbuilder.population`.

    Parameters
    ----------
    network : Network or Subnetwork
        The pipes to allocate to: every pipe of a :class:`~swmmbuilder.Network`, or those of a
        :class:`~swmmbuilder.Subnetwork`. A polygon's population is shared over the pipes given,
        so for a catchment whose boundary cuts polygons, allocate over the whole Network and take
        the catchment's nodes with :meth:`PopulationAllocation.within`; allocating over the
        catchment alone puts all the people of a cut polygon on its part inside.

    polygons : geopandas.GeoDataFrame
        The areal units, indexed by a unique ID, with a CRS (they are reprojected to the
        network's).

    population : str or pandas.Series
        Population of each polygon: a column of ``polygons``, or a Series with the same index.
        A missing or non-numeric value counts as 0.

    max_distance : float, optional
        None (default) assigns each pipe to the polygon containing its midpoint. A distance, in
        the units of the network's CRS (m for London's), assigns each pipe instead to every
        polygon within that distance of it (the pipe line meets the polygon buffered by the
        distance), so a polygon is served by the pipes in it and along its edges, and a pipe may
        serve several polygons, with its whole weight in each. Polygons with no pipe that close
        contribute no one, unless flagged ``sewered``.

    sewered : str or pandas.Series, optional
        Polygons known to be sewered (e.g. residentially zoned): a column of booleans of
        ``polygons``, or a Series with the same index. A flagged polygon left with no pipe by the
        midpoint or distance rule (a development on private sewers, say) is given its nearest
        pipe, at any distance; an unflagged one contributes no one. By default no polygon is
        flagged.

    serve_force_mains : bool
        Let force mains serve polygons (default False: they take no service connections, so they
        get weight 0, even from ``fallback_weight``, and no population).

    weight : str, callable or None
        Weight of each pipe in its polygon's share: a column of the pipes (``full_data``), e.g. a
        count of service connections such as London's ``PDC_Count``; a function of the pipes
        returning one value per pipe; or None (default) for the pipe length.

    max_weight : float, optional
        Weights above this are invalid. Missing, non-numeric, and negative weights always are.

    impute : {'per_length', 'zero'}
        What an invalid weight becomes: the pipe's length times the median weight per metre
        (weight / length) of the pipes of its polygon with valid weights (default), or 0.

    fallback_weight : float or None
        The weight of a pipe that would otherwise have none to share by (default 1, e.g. one
        connection). With ``impute='per_length'``, an invalid weight in a polygon whose median
        weight per metre is 0, or which has no valid weight, becomes ``fallback_weight``; then
        every pipe of a polygon whose pipes still weigh 0 in all gets it, so that the polygon's
        people are shared equally among its pipes. None leaves such polygons unallocated, and
        reports them in ``no_valid_weight`` and ``zero_weight``.

    placement : {'upstream', 'downstream', 'split'}
        Where a pipe's share is placed: at its upstream node (``from_node``, default), at its
        downstream node, or half at each.

    Returns
    -------
    PopulationAllocation
        The population at each node, the tables per pipe and per polygon, and the population
        left out.

    Raises
    ------
    ValueError
        If an option is invalid, a polygon ID is repeated, a population is negative, or a layer
        has no CRS; the message lists every problem.

    Warns
    -----
    DataWarning
        With the midpoint rule, if a pipe midpoint is in more than one polygon (it takes the
        first).

    Examples
    --------
    London's census dissemination areas, shared by private drain connections, with invalid
    counts imputed from each area's median connections per metre (one connection where that
    median is 0, and one per pipe in an area with no connections)::

        alloc = allocate_population(net, das, population="pop", weight="PDC_Count", max_weight=100)
        sub = net.catchment("OUT")
        dwf = dwf_from_population(alloc.within(sub), per_capita=230, unit="L/d")

    Census blocks, each served by the pipes within 15 m of it::

        alloc = allocate_population(
            net, blocks, population="pop", max_distance=15, weight="PDC_Count", max_weight=100
        )
    """
    pipes = _pipes_of(network)
    _check(pipes, polygons, population, max_distance, sewered, weight, impute, placement)
    people = polygons[population] if isinstance(population, str) else population
    people = pd.to_numeric(pd.Series(people).reindex(polygons.index), errors="coerce")
    negative = people.index[people < 0]
    if len(negative):
        raise ValueError(f"Negative population in polygons {list(negative)[:20]}")
    if people.isna().any():
        logger.info("%d polygons have no population; counted as 0", int(people.isna().sum()))
    people = people.astype(float).fillna(0.0)
    polygons = polygons.to_crs(pipes.crs)

    if serve_force_mains or "structure_type" not in pipes.columns:
        serving = np.ones(len(pipes), dtype=bool)
    else:
        serving = (pipes["structure_type"] != "force_main").to_numpy()
    rows = _assign(pipes, polygons, max_distance, serving)
    if max_distance is not None:
        served = rows.dropna().groupby("pipe").size()
        logger.info(
            "pipes within %g of a polygon: %d, of which %d serve more than one",
            max_distance,
            len(served),
            int((served > 1).sum()),
        )
    nearest = pd.Series(dtype=float)
    if sewered is not None:
        flag = polygons[sewered] if isinstance(sewered, str) else sewered
        flag = pd.Series(flag).reindex(polygons.index).fillna(False).astype(bool)
        unserved = flag & ~polygons.index.isin(rows["polygon"].dropna())
        if unserved.any():
            found = _nearest(pipes, polygons[unserved], serving)
            nearest = found["distance"]
            added = pd.DataFrame({"pipe": found["pipe"].to_numpy(), "polygon": found.index})
            # A pipe that served no polygon loses its empty row.
            empty = rows["polygon"].isna() & rows["pipe"].isin(added["pipe"])
            rows = pd.concat([rows[~empty], added], ignore_index=True)
            rows = rows.sort_values("pipe", kind="stable", ignore_index=True)
            logger.info(
                "%d polygons flagged sewered had no pipe; given their nearest (median %.0f, "
                "max %.0f away)",
                len(found),
                float(nearest.median()),
                float(nearest.max()),
            )
    # One row per pipe and polygon, positional until the end.
    position = rows["pipe"].to_numpy()
    table = pipes[["pipe_id", "from_node", "to_node", "length"]].iloc[position]
    table = table.reset_index(drop=True)
    table["polygon"] = rows["polygon"].to_numpy()

    given = pd.Series(_weights(pipes, weight).to_numpy()[position])
    invalid = given.isna() | (given < 0)
    if max_weight is not None:
        invalid |= given > max_weight
    table["weight_imputed"] = invalid
    in_polygon = table["polygon"].notna()
    fallen_back = pd.Series(False, index=table.index)
    if impute == "zero":
        table["weight"] = given.where(~invalid, 0.0)
    else:
        per_metre = (given.where(~invalid) / table["length"]).groupby(table["polygon"]).median()
        rate = table["polygon"].map(per_metre).astype(float)
        table["weight"] = given.where(~invalid, (rate * table["length"]).fillna(0.0))
        if fallback_weight is not None:
            # No rate to impute from: the median is 0, or the polygon has no valid weight.
            fallen_back = invalid & in_polygon & (rate.isna() | (rate == 0))
            table["weight"] = table["weight"].where(~fallen_back, float(fallback_weight))
    if fallback_weight is not None:
        weightless = table["polygon"].map(table.groupby("polygon")["weight"].sum()) == 0
        weightless &= in_polygon
        table["weight"] = table["weight"].where(~weightless, float(fallback_weight))
        fallen_back |= weightless
    table["weight_fallback"] = fallen_back
    # A pipe that serves no polygon by rule (a force main) takes no connections.
    table.loc[~serving[position], "weight"] = 0.0
    no_valid = table.loc[~invalid, "polygon"].dropna().unique()
    logger.info(
        "%d of %d pipe weights invalid; imputed (%s); %d pipes given the fallback weight",
        int(invalid.sum()),
        len(table),
        impute,
        int(fallen_back.sum()),
    )

    totals = table.groupby("polygon")["weight"].sum()
    per_weight = (people / totals.reindex(polygons.index)).replace([np.inf, -np.inf], np.nan)
    table["population"] = (table["weight"] * table["polygon"].map(per_weight)).fillna(0.0)

    table = table.set_index("pipe_id")
    nodes = _place(table, placement)
    placed = sum(nodes.values())
    summary = pd.DataFrame(
        {
            "population": people,
            "pipes": table.groupby("polygon").size().reindex(polygons.index, fill_value=0),
            "weight": totals.reindex(polygons.index, fill_value=0.0),
            "allocated": table.groupby("polygon")["population"]
            .sum()
            .reindex(polygons.index, fill_value=0.0),
        }
    )
    with_people = summary["population"] > 0
    no_pipes = list(summary.index[with_people & (summary["pipes"] == 0)])
    sewered = with_people & (summary["pipes"] > 0)
    no_valid_weight = (
        list(summary.index[sewered & ~summary.index.isin(no_valid)])
        if impute == "per_length" and fallback_weight is None
        else []
    )
    zero_weight = list(
        summary.index[sewered & (summary["weight"] == 0) & ~summary.index.isin(no_valid_weight)]
    )
    fallback = list(table.loc[fallen_back.to_numpy(), "polygon"].unique())
    unallocated = float(people.sum() - placed)
    if no_valid_weight:
        logger.warning(
            "%d polygons have pipes but no valid weight to impute from; not allocated: %s",
            len(no_valid_weight),
            no_valid_weight[:20],
        )
    if no_pipes or zero_weight:
        logger.warning(
            "%d polygons with no pipe and %d whose pipes all weigh 0 are not allocated",
            len(no_pipes),
            len(zero_weight),
        )
    logger.info(
        "population: %.0f of %.0f placed at %d nodes; %.0f not allocated",
        placed,
        float(people.sum()),
        sum(p > 0 for p in nodes.values()),
        unallocated,
    )
    return PopulationAllocation(
        nodes=nodes,
        pipes=table,
        polygons=summary,
        unallocated=unallocated,
        no_pipes=no_pipes,
        nearest=nearest,
        no_valid_weight=no_valid_weight,
        zero_weight=zero_weight,
        fallback=fallback,
        placement=placement,
    )
