"""Spreading a flow over the pipes of a network, e.g. groundwater infiltration.

Infiltration is usually given as a total for a catchment, or as a rate per metre of pipe, per
millimetre-kilometre (diameter times length), or per hectare. :func:`distribute_flow` shares a
total among the pipes in proportion to a weight (by default the length), or gives each pipe a
rate times its weight, and places each pipe's flow at its upstream node (or as asked). The
result, {node: flow}, can be passed as ``dwf`` (a steady flow with no pattern) or as the
``baseline`` of :func:`~swmmbuilder.pattern_inflows`.
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Callable
from typing import TYPE_CHECKING

import geopandas as gpd
import pandas as pd

from .diagnostics import DataWarning, Problems
from .population import PLACEMENTS, Placement, _log_left_out, _missing_node, _pipes_of
from .units import _is_real

if TYPE_CHECKING:
    from ._types import NodeId
    from .network import Network
    from .subnetwork import Subnetwork

logger = logging.getLogger(__name__)

PipeWeight = str | Callable[[gpd.GeoDataFrame], pd.Series]


def pipe_weights(pipes: gpd.GeoDataFrame, weight: PipeWeight = "length") -> pd.Series:
    """The weight of each pipe: ``'length'`` (m), ``'diameter_length'`` (diameter times length,
    m^2), a column of the pipes, or a function of the pipes returning one value per pipe.
    Missing, non-numeric, or negative values are returned as missing."""
    if callable(weight):
        values = pd.Series(weight(pipes), index=pipes.index)
    elif weight == "diameter_length":
        values = pd.to_numeric(pipes["diameter"], errors="coerce") * pipes["length"]
    else:
        values = pipes[weight]
    values = pd.to_numeric(values, errors="coerce").astype(float)
    return values.where(values >= 0)


def distribute_flow(
    network: Network | Subnetwork,
    total: float | None = None,
    rate: float | None = None,
    weight: PipeWeight = "length",
    serve_force_mains: bool = False,
    placement: Placement = "upstream",
) -> dict[NodeId, float]:
    """Spread a flow over the pipes of ``network`` and place it at their nodes.

    Give either ``total``, shared among the pipes in proportion to ``weight``, or ``rate``, the
    flow of each pipe per unit of its weight.

    Parameters
    ----------
    network : Network or Subnetwork
        The pipes to spread over: every pipe of a :class:`~swmmbuilder.Network`, or those of a
        :class:`~swmmbuilder.Subnetwork`, such as a catchment. A pipe with a missing
        ``from_node`` or ``to_node`` is in no catchment, so it is left out, and the flow it would
        have taken is logged.

    total : float, optional
        The flow to share, in the unit wanted for the result (e.g. L/s, the FLOW_UNITS of the
        shipped template).

    rate : float, optional
        Flow per unit of weight, e.g. L/s per metre of pipe with ``weight='length'``.

    weight : str or callable
        ``'length'`` (default), ``'diameter_length'`` (diameter times length, both in m), a
        column of the pipes, or a function of the pipes returning one value per pipe. A pipe
        whose weight is missing, non-numeric, or negative takes no flow (with a warning).

    serve_force_mains : bool
        Let force mains take flow (default False: groundwater does not enter a pressurised
        main).

    placement : {'upstream', 'downstream', 'split'}
        Where a pipe's flow is placed: its upstream node (default), its downstream node, or half
        at each.

    Returns
    -------
    dict
        {node: flow}, in order of the node IDs as text, for every node a pipe's flow is placed at.

    Raises
    ------
    ValueError
        If neither or both of ``total`` and ``rate`` are given, either is negative or not a
        number, ``placement`` or ``weight`` is invalid, or ``total`` is positive but the pipes
        weigh nothing.

    Warns
    -----
    DataWarning
        If some pipes have no valid weight; they take no flow.

    Examples
    --------
    Infiltration of 40 L/s over the gravity pipes of a catchment, by length::

        sub = net.catchment("OUT")
        infiltration = distribute_flow(sub, total=40.0)
        inflows, series = pattern_inflows(flows.nodes, patterns, days=3, baseline=infiltration)
    """
    pipes = _pipes_of(network)
    problems = Problems("Invalid flow distribution")
    if (total is None) == (rate is None):
        problems.add("give one of total and rate")
    for name, value in (("total", total), ("rate", rate)):
        if value is not None and not (_is_real(value) and value >= 0):
            problems.add(f"{name} must be a non-negative number, not {value!r}")
    if placement not in PLACEMENTS:
        problems.add(f"placement must be one of {', '.join(PLACEMENTS)}, not {placement!r}")
    if isinstance(weight, str) and weight not in ("length", "diameter_length"):
        if weight not in pipes.columns:
            problems.add(f"weight: the pipes have no column {weight!r}")
    problems.raise_if_any()

    weights = pipe_weights(pipes, weight)
    if not serve_force_mains and "structure_type" in pipes.columns:
        weights = weights.where(pipes["structure_type"] != "force_main", 0.0)
    no_node = _missing_node(pipes)
    left_out = float(weights[no_node].fillna(0.0).sum())
    everything, pipes, weights = pipes, pipes[~no_node], weights[~no_node]
    missing = weights.isna()
    label = weight if isinstance(weight, str) else "callable"
    if missing.any():
        warnings.warn(
            f"{int(missing.sum())} pipes have no valid weight ({label}); they take no flow: "
            f"{list(weights.index[missing])[:20]}",
            DataWarning,
            stacklevel=2,
        )
        weights = weights.fillna(0.0)
    if total is not None:
        weight_sum = float(weights.sum())
        if weight_sum == 0 and total > 0:
            raise ValueError(f"total is {total}, but the pipes weigh nothing to share it by.")
        flows = weights / weight_sum * total if weight_sum > 0 else weights * 0.0
    else:
        flows = weights * rate
    if no_node.any():
        if total is not None:
            taken = total * left_out / (float(weights.sum()) + left_out) if left_out > 0 else 0.0
        else:
            taken = left_out * rate
        _log_left_out(
            logger,
            everything,
            no_node,
            f"with them they would take {taken:.6g}, in no catchment",
        )

    if placement == "split":
        half = flows / 2
        placed = (
            pd.concat(
                [
                    half.groupby(pipes["from_node"].to_numpy()).sum(),
                    half.groupby(pipes["to_node"].to_numpy()).sum(),
                ]
            )
            .groupby(level=0)
            .sum()
        )
    else:
        end = pipes["from_node" if placement == "upstream" else "to_node"].to_numpy()
        placed = flows.groupby(end).sum()
    placed = placed[placed > 0]
    logger.info(
        "distributed flow: %.6g at %d nodes over %d pipes (weight %s)",
        float(placed.sum()),
        len(placed),
        int((flows > 0).sum()),
        label,
    )
    return {node: float(placed[node]) for node in sorted(placed.index, key=str)}
