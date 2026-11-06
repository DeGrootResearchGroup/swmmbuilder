"""External inflows and time series: the ``inflows`` and ``timeseries`` options of
:func:`~swmmbuilder.swmm.build_model`, written to [INFLOWS] and [TIMESERIES]."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import TYPE_CHECKING

from ..diagnostics import Problems
from ..units import _is_finite_number

if TYPE_CHECKING:
    from .._types import Fields, NodeId


@dataclass(frozen=True)
class Inflow:
    """An external flow into a node, written to [INFLOWS] (see the ``inflows`` option of
    :func:`build_model`).

    SWMM adds it to the node's dry-weather flow: at time t the inflow is
    ``scale * series(t) + baseline * pattern(t)``, in the template's FLOW_UNITS (L/s with the
    shipped template). A node has at most one, so a node that needs a second time-varying flow
    besides its [DWF] (one pattern of each type) takes it here.

    Attributes
    ----------
    timeseries : str, optional
        Name of a time series of the ``timeseries`` option; None for a baseline only.

    scale : float
        Factor applied to the time series (default 1). With a series of unit-mean multipliers,
        the node's mean flow.

    baseline : float
        Constant flow added to the series (default 0), e.g. groundwater infiltration.

    baseline_pattern : str, optional
        Name of a pattern of the ``patterns`` option that multiplies the baseline.
    """

    timeseries: str | None = None
    scale: float = 1.0
    baseline: float = 0.0
    baseline_pattern: str | None = None

    @classmethod
    def coerce(cls, value: Inflow | Fields) -> Inflow:
        """An Inflow from an Inflow or a dict of its fields."""
        return value if isinstance(value, cls) else cls(**value)


def hourly_series(
    multipliers: Sequence[float], days: int, ramp_minutes: float = 1.0
) -> list[tuple[float, float]]:
    """(hours from the start, value) points of a time series that holds each of 24 hourly
    ``multipliers`` for the whole hour, repeated for ``days`` days.

    SWMM interpolates a time series linearly between its points, so each hour has a point at its
    start and one ``ramp_minutes`` before its end; the value changes over that last interval. A
    closing point at ``24 * days`` h (the first multiplier) keeps the end of the run inside the
    series, since SWMM takes the value beyond the last point as zero.

    Parameters
    ----------
    multipliers : sequence of float
        The 24 hourly values, from 00:00.

    days : int
        Number of days the series covers (at least the length of the run).

    ramp_minutes : float
        Length of the change from one hour's value to the next (default 1 minute), more than 0
        and less than 60.

    Returns
    -------
    list of (float, float)
        Points to pass in ``timeseries``.
    """
    values = [float(v) for v in multipliers]
    if len(values) != 24 or not all(_is_finite_number(v) for v in values):
        raise ValueError(f"hourly_series needs 24 finite multipliers, not {len(values)}.")
    if int(days) != days or days < 1:
        raise ValueError(f"days must be a positive whole number, not {days!r}.")
    if not 0 < ramp_minutes < 60:
        raise ValueError(f"ramp_minutes must be more than 0 and less than 60, not {ramp_minutes}.")
    hold = 1 - ramp_minutes / 60
    points = []
    for day in range(int(days)):
        for hour, value in enumerate(values):
            start = 24 * day + hour
            points += [(float(start), value), (start + hold, value)]
    points.append((24.0 * days, values[0]))
    return points


def resolve_timeseries(timeseries):
    """{name: ((time, value), ...)} from the ``timeseries`` option of :func:`build_model`.

    Raises
    ------
    ValueError
        If a series is empty, or a time or value is not a finite number, or the times are
        negative or not increasing (SWMM's ERROR 173); the message lists every problem.
    """
    resolved, problems = {}, []
    for name, points in (timeseries or {}).items():
        try:
            points = [(t, v) for t, v in points]
        except (TypeError, ValueError):
            problems.append(f"{name}: give a list of (time in hours, value) pairs")
            continue
        if not points:
            problems.append(f"{name}: no points")
            continue
        if not all(_is_finite_number(t) and _is_finite_number(v) for t, v in points):
            problems.append(f"{name}: times and values must be finite numbers")
            continue
        points = tuple((float(t), float(v)) for t, v in points)
        times = [t for t, _ in points]
        if times[0] < 0 or any(b <= a for a, b in pairwise(times)):
            problems.append(f"{name}: times must be at least 0 h and increasing")
            continue
        resolved[name] = points
    if problems:
        raise ValueError(Problems("Invalid timeseries", problems).message())
    return resolved


def resolve_inflows(subnetwork, inflows, timeseries, patterns, outfalls=()) -> dict[NodeId, Inflow]:
    """{node: :class:`Inflow`} from the ``inflows`` option of :func:`build_model`, checked
    against the resolved ``timeseries`` and ``patterns``.

    Raises
    ------
    ValueError
        If a node is not in the subnetwork or is an outfall, a time series or pattern is not
        defined, or the scale or baseline is negative or not finite; the message lists every
        problem.
    """
    if not inflows:
        return {}
    if not isinstance(inflows, Mapping):
        raise TypeError("inflows must be a mapping {node: Inflow or dict of its fields}.")
    inflows = subnetwork.match_nodes(inflows, "inflows")
    nodes, outfalls = set(subnetwork.nodes), set(outfalls)
    resolved, problems = {}, []
    for node, value in inflows.items():
        try:
            inflow = Inflow.coerce(value)
        except TypeError as error:
            problems.append(f"{node}: {error}")
            continue
        if node not in nodes:
            problems.append(f"{node}: not in the subnetwork")
        elif node in outfalls:
            problems.append(f"{node}: an outfall (inflows cannot enter at outfalls)")
        if inflow.timeseries is not None and inflow.timeseries not in timeseries:
            problems.append(f"{node}: time series {inflow.timeseries!r} not in timeseries")
        if inflow.baseline_pattern is not None and inflow.baseline_pattern not in patterns:
            problems.append(f"{node}: pattern {inflow.baseline_pattern!r} not in patterns")
        for field in ("scale", "baseline"):
            x = getattr(inflow, field)
            if not (_is_finite_number(x) and float(x) >= 0):
                problems.append(f"{node}: {field} must be non-negative and finite, not {x!r}")
        resolved[node] = inflow
    if problems:
        raise ValueError(Problems("Invalid inflows", problems).message())
    return resolved
