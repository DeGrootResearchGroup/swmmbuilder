"""External inflows and time series: the ``inflows`` and ``timeseries`` options of
:func:`~swmmbuilder.swmm.build_model`, written to [INFLOWS] and [TIMESERIES]."""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from itertools import pairwise
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from ..diagnostics import DesignWarning, Problems
from ..units import _is_finite_number

if TYPE_CHECKING:
    from .._types import Fields, NodeId


# Tolerance on the mean of an hourly pattern: [DWF] patterns and pattern_inflows warn beyond it,
# and require_unit_mean=True refuses a pattern beyond it.
UNIT_MEAN_TOLERANCE = 1e-3


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
        if isinstance(value, Inflow):
            return value
        return cls(**value)


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


def unit_mean_tolerance(require_unit_mean: bool | float | None) -> float | None:
    """The tolerance on a pattern's mean meant by ``require_unit_mean`` (see
    :func:`pattern_inflows`), or None for no check (False) or a warning only (None)."""
    if require_unit_mean is None:
        return None
    if isinstance(require_unit_mean, bool):
        return UNIT_MEAN_TOLERANCE if require_unit_mean else None
    if not (_is_finite_number(require_unit_mean) and float(require_unit_mean) >= 0):
        raise ValueError(
            "require_unit_mean must be True, False, or a non-negative tolerance, "
            f"not {require_unit_mean!r}."
        )
    return float(require_unit_mean)


def pattern_mean_problem(name: Any, pattern: Sequence[float], tolerance: float) -> str | None:
    """Why the mean of hourly ``pattern`` is not 1 within ``tolerance``, or None if it is (or if
    the pattern is not numbers, which is reported elsewhere)."""
    try:
        values = np.asarray(pattern, dtype=float)
    except (TypeError, ValueError):
        return None
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        return None
    mean = float(values.mean())
    if abs(mean - 1.0) > tolerance:
        return f"pattern {name} has a mean of {mean:.6g}, not 1 (tolerance {tolerance:g})"
    return None


def pattern_inflows(
    flows: Mapping[Any, Mapping[NodeId, float]] | pd.DataFrame,
    patterns: Mapping[Any, Sequence[float]],
    days: int,
    baseline: Mapping[NodeId, float] | None = None,
    prefix: str = "",
    ramp_minutes: float = 1.0,
    require_unit_mean: bool | float | None = None,
) -> tuple[dict[NodeId, Inflow], dict[str, list[tuple[float, float]]]]:
    """The ``inflows`` and ``timeseries`` options of :func:`build_model` for mean flows that each
    follow an hourly pattern, e.g. the land-use flows of :func:`~swmmbuilder.allocate_land_use`.

    Each pattern becomes a time series ``<prefix><name>`` (:func:`hourly_series` over ``days``).
    A node with flow of one kind gets an :class:`Inflow` of that series scaled by its mean flow.
    A node with flows of several kinds gets its own series ``<prefix><node>``, the mean of their
    patterns weighted by their flows, scaled by their sum. ``baseline`` adds a steady flow (e.g.
    groundwater infiltration) as the inflow's baseline, at nodes with or without a patterned
    flow.

    Parameters
    ----------
    flows : pandas.DataFrame or dict
        Mean flow at each node by kind: a DataFrame of nodes (index) by kinds (columns), such as
        :attr:`LandUseAllocation.nodes <swmmbuilder.LandUseAllocation.nodes>`, or {kind: {node:
        flow}}. Flows in the template's FLOW_UNITS.

    patterns : dict
        {kind: 24 hourly multipliers, from 00:00}, one for every kind with a flow. The
        multipliers should have a mean of 1, so that the mean flow is the one given.

    days : int
        Number of days the time series cover (at least the length of the run).

    baseline : dict, optional
        {node: steady flow}.

    prefix : str
        Prefix of the time series names (default none), e.g. 'NR_'.

    ramp_minutes : float
        See :func:`hourly_series`.

    require_unit_mean : bool, float, or None
        What to do with a pattern whose mean is not 1, which scales the mean flows: None
        (default) gives a :class:`~swmmbuilder.diagnostics.DesignWarning` when it is off by more
        than :data:`UNIT_MEAN_TOLERANCE` (0.001), as a [DWF] pattern does; True refuses it at
        that tolerance, a number at that tolerance; False accepts any mean silently.

    Returns
    -------
    (dict, dict)
        ``inflows`` ({node: Inflow}, by node in order of their IDs as text) and ``timeseries``
        ({name: points}, every pattern first, then the mixed nodes).

    Raises
    ------
    ValueError
        If a kind with a flow has no pattern, a pattern is not 24 finite numbers, a flow or
        baseline is negative or not finite, a series name would be used twice, or (with
        ``require_unit_mean``) a pattern's mean is not 1.
    """
    table = flows if isinstance(flows, pd.DataFrame) else pd.DataFrame(flows)
    table = table.fillna(0.0).astype(float)
    baseline = dict(baseline or {})
    problems = Problems("Invalid pattern inflows")
    used = [k for k in table.columns if (table[k] > 0).any()]
    missing = [str(k) for k in used if k not in patterns]
    if missing:
        problems.add(f"no pattern for {', '.join(missing)}")
    if not np.isfinite(table.to_numpy()).all() or (table.to_numpy() < 0).any():
        problems.add("flows must be non-negative and finite")
    bad = [str(n) for n, q in baseline.items() if not (_is_finite_number(q) and float(q) >= 0)]
    if bad:
        problems.add(f"baselines must be non-negative and finite: {', '.join(bad[:20])}")
    tolerance = unit_mean_tolerance(require_unit_mean)
    if tolerance is not None:
        for kind, pattern in patterns.items():
            problem = pattern_mean_problem(kind, pattern, tolerance)
            if problem:
                problems.add(problem)
    problems.raise_if_any()
    if require_unit_mean is None:
        for kind, pattern in patterns.items():
            problem = pattern_mean_problem(kind, pattern, UNIT_MEAN_TOLERANCE)
            if problem:
                warnings.warn(
                    f"{problem[0].upper()}{problem[1:]}, so the mean flows of {kind} are scaled "
                    "by it.",
                    DesignWarning,
                    stacklevel=2,
                )

    series = {f"{prefix}{k}": hourly_series(p, days, ramp_minutes) for k, p in patterns.items()}
    multipliers = {k: np.asarray(p, dtype=float) for k, p in patterns.items()}
    totals = table.sum(axis=1)
    nodes = set(totals.index[totals > 0]) | {n for n, q in baseline.items() if float(q) > 0}
    inflows = {}
    for node in sorted(nodes, key=str):
        base = float(baseline.get(node, 0.0))
        total = float(totals.get(node, 0.0))
        if total <= 0:
            inflows[node] = Inflow(baseline=base)
            continue
        q = table.loc[node]
        kinds = q[q > 0]
        if len(kinds) == 1:
            name = f"{prefix}{kinds.index[0]}"
        else:
            name = f"{prefix}{node}"
            if name in series:
                raise ValueError(
                    f"Time series name {name!r} is both a pattern and a node with mixed flows; "
                    "give another prefix."
                )
            mixed = sum(kinds[k] * multipliers[k] for k in kinds.index) / total
            series[name] = hourly_series(list(mixed), days, ramp_minutes)
        inflows[node] = Inflow(name, scale=total, baseline=base)
    return inflows, series
