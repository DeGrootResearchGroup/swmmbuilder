"""Allocation of the population of areal units (polygons) to the nodes of a network.

The population of each polygon (a census tract or dissemination area, a postal code area, a
traffic zone, or any other areal unit) is placed on the pipes in it and then on their nodes:

1. each pipe goes to the polygon containing its midpoint;
2. each polygon's population is shared among its pipes in proportion to a weight per pipe (by
   default the pipe length; e.g. a count of service connections), with invalid or missing weights
   imputed from a rate per metre of pipe;
3. each pipe's share is placed at its upstream node, its downstream node, or half at each.

A polygon with no pipe contributes no one, and a polygon only partly sewered places all its
people on the pipes it has. The population left out is reported.

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
RATE_GROUPS = ("network", "polygons")
RATES = ("median", "total")
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
        One row per pipe, indexed by pipe_id, with ``from_node``, ``to_node``, ``length``,
        ``polygon`` (ID of the polygon containing the midpoint, missing if none),
        ``weight_imputed`` (bool), ``rate_group`` (the group whose rate per metre imputed the
        weight), ``weight`` (after imputation), and ``population``.

    polygons : pandas.DataFrame
        One row per polygon, indexed by polygon ID, with ``population`` (as given, a missing value
        as 0), ``pipes`` (the number of pipes in it), ``weight`` (their total weight), and
        ``allocated`` (the population placed on its pipes).

    unallocated : float
        Population placed at no node: that of the polygons with no pipe or whose pipes all weigh
        0, and the shares of pipes with no node at the end they are placed at.

    no_pipes : list
        IDs of the polygons with people but no pipe.

    zero_weight : list
        IDs of the polygons with people and pipes, whose pipes all weigh 0.

    placement : str
        The ``placement`` the shares were placed with.
    """

    nodes: dict[NodeId, float]
    pipes: pd.DataFrame
    polygons: pd.DataFrame
    unallocated: float
    no_pipes: list = field(default_factory=list)
    zero_weight: list = field(default_factory=list)
    placement: Placement = "upstream"

    @property
    def total(self) -> float:
        """Population placed at the nodes."""
        return float(sum(self.nodes.values()))

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


def _check(pipes, polygons, population, weight, impute, rate_groups, rate, placement):
    """Raise ValueError listing every invalid option."""
    problems = Problems("Invalid population allocation")
    if impute not in IMPUTATIONS:
        problems.add(f"impute must be one of {', '.join(IMPUTATIONS)}, not {impute!r}")
    if isinstance(rate_groups, str) and rate_groups not in RATE_GROUPS:
        problems.add(
            f"rate_groups must be one of {', '.join(RATE_GROUPS)} or a GeoDataFrame, "
            f"not {rate_groups!r}"
        )
    if rate not in RATES:
        problems.add(f"rate must be one of {', '.join(RATES)}, not {rate!r}")
    if placement not in PLACEMENTS:
        problems.add(f"placement must be one of {', '.join(PLACEMENTS)}, not {placement!r}")
    if isinstance(weight, str) and weight not in pipes.columns:
        problems.add(f"weight: the pipes have no column {weight!r}")
    if isinstance(population, str) and population not in polygons.columns:
        problems.add(f"population: polygons has no column {population!r}")
    repeated = polygons.index[polygons.index.duplicated()].unique()
    if len(repeated):
        problems.add(f"polygon IDs (the index of polygons) repeated: {list(repeated)[:20]}")
    layers = [("the network", pipes), ("polygons", polygons)]
    if not isinstance(rate_groups, str):
        layers.append(("rate_groups", rate_groups))
    for name, layer in layers:
        if layer.crs is None:
            problems.add(f"{name} has no CRS")
    problems.raise_if_any()


def allocate_population(
    network: Network | Subnetwork,
    polygons: gpd.GeoDataFrame,
    population: str | pd.Series = "population",
    weight: Weight = None,
    max_weight: float | None = None,
    impute: Literal["per_length", "zero"] = "per_length",
    rate_groups: gpd.GeoDataFrame | Literal["network", "polygons"] = "network",
    rate: Literal["median", "total"] = "median",
    placement: Placement = "upstream",
) -> PopulationAllocation:
    """Allocate the population of each polygon to the nodes of ``network``.

    Each pipe is assigned to the polygon containing its midpoint, each polygon's population is
    shared among its pipes in proportion to ``weight``, and each pipe's share is placed at its
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

    weight : str, callable or None
        Weight of each pipe in its polygon's share: a column of the pipes (``full_data``), e.g. a
        count of service connections such as London's ``PDC_Count``; a function of the pipes
        returning one value per pipe; or None (default) for the pipe length.

    max_weight : float, optional
        Weights above this are invalid. Missing, non-numeric, and negative weights always are.

    impute : {'per_length', 'zero'}
        What an invalid weight becomes: the pipe's length times the rate per metre of its group
        (``rate_groups``, ``rate``), or 0. With ``'per_length'``, a pipe in no group, or in a
        group with no valid weight, gets 0.

    rate_groups : {'network', 'polygons'} or geopandas.GeoDataFrame
        The pipes each rate per metre is taken over: all of them (default), those of the pipe's
        own polygon, or those of the area of this GeoDataFrame containing the pipe's midpoint
        (e.g. a coarser areal unit than ``polygons``; the first area if more than one does).

    rate : {'median', 'total'}
        The rate per metre of a group: the median of weight / length over its pipes with valid
        weights (default), or their total weight over their total length.

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
        If a pipe midpoint is in more than one polygon (it takes the first).

    Examples
    --------
    London's census dissemination areas, shared by private drain connections, with invalid
    counts imputed from each forward sortation area's median connections per metre::

        alloc = allocate_population(
            net, das, population="pop", weight="PDC_Count", max_weight=100, rate_groups=fsas
        )
        sub = net.catchment("OUT")
        dwf = dwf_from_population(alloc.within(sub), per_capita=230, unit="L/d")
    """
    pipes = _pipes_of(network)
    _check(pipes, polygons, population, weight, impute, rate_groups, rate, placement)
    people = polygons[population] if isinstance(population, str) else population
    people = pd.to_numeric(pd.Series(people).reindex(polygons.index), errors="coerce")
    negative = people.index[people < 0]
    if len(negative):
        raise ValueError(f"Negative population in polygons {list(negative)[:20]}")
    if people.isna().any():
        logger.info("%d polygons have no population; counted as 0", int(people.isna().sum()))
    people = people.astype(float).fillna(0.0)
    polygons = polygons.to_crs(pipes.crs)

    midpoints = pipes.geometry.interpolate(0.5, normalized=True)
    table = pipes[["from_node", "to_node", "length"]].copy()
    table["polygon"] = _polygon_of(midpoints, polygons)

    given = _weights(pipes, weight)
    invalid = given.isna() | (given < 0)
    if max_weight is not None:
        invalid |= given > max_weight
    table["weight_imputed"] = invalid
    if impute == "zero":
        table["rate_group"] = pd.NA
        table["weight"] = given.where(~invalid, 0.0)
    else:
        if not isinstance(rate_groups, str):
            group = _polygon_of(midpoints, rate_groups.to_crs(pipes.crs))
        elif rate_groups == "polygons":
            group = table["polygon"]
        else:
            group = pd.Series("network", index=pipes.index)
        table["rate_group"] = group
        valid = given.where(~invalid)
        if rate == "median":
            per_metre = (valid / table["length"]).groupby(group).median()
        else:
            length = table["length"].where(~invalid)
            per_metre = valid.groupby(group).sum() / length.groupby(group).sum()
        imputed = group.map(per_metre).astype(float).fillna(0.0) * table["length"]
        table["weight"] = given.where(~invalid, imputed)
        no_rate = int((invalid & (imputed == 0)).sum())
        if no_rate:
            logger.info(
                "%d pipes with invalid weights have no rate per metre (no group, or no valid "
                "weight in it); weight 0",
                no_rate,
            )
    logger.info(
        "%d of %d pipe weights invalid; imputed (%s)", int(invalid.sum()), len(table), impute
    )

    totals = table.groupby("polygon")["weight"].sum()
    per_weight = (people / totals.reindex(polygons.index)).replace([np.inf, -np.inf], np.nan)
    table["population"] = (table["weight"] * table["polygon"].map(per_weight)).fillna(0.0)

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
    zero_weight = list(
        summary.index[with_people & (summary["pipes"] > 0) & (summary["weight"] == 0)]
    )
    unallocated = float(people.sum() - placed)
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
        zero_weight=zero_weight,
        placement=placement,
    )
