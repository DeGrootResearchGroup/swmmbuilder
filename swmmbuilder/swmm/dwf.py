"""Dry-weather flow: the ``dwf``, ``patterns`` and ``dwf_patterns`` options of
:func:`~swmmbuilder.swmm.build_model`."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

from ..diagnostics import DesignWarning, Problems
from ..units import _FLOW_UNIT_LPS, _TIME_SECONDS, _VOLUME_LITRES, SI_FLOW_UNITS, _is_finite_number

if TYPE_CHECKING:
    from collections.abc import Mapping

    from .._types import NodeId

# SWMM time pattern types and their number of multipliers. A pattern given as a plain list of
# multipliers takes its type from its length; WEEKEND (24 values) must be named explicitly.
PATTERN_LENGTHS = {"MONTHLY": 12, "DAILY": 7, "HOURLY": 24, "WEEKEND": 24}


PATTERN_TYPES_BY_LENGTH = {12: "MONTHLY", 7: "DAILY", 24: "HOURLY"}


# SWMM's limit on the time patterns of a node's dry-weather flow (one of each type).
MAX_DWF_PATTERNS = 4


def resolve_patterns(patterns):
    """{name: (type, multipliers)} from the ``patterns`` option of :func:`build_model`.

    Each pattern is a list of multipliers, whose type is set by its length (24 HOURLY, 7 DAILY,
    12 MONTHLY), or a (type, multipliers) pair. A warning is given for an HOURLY or WEEKEND
    pattern whose mean is not 1, since the mean then scales the base flows.

    Raises
    ------
    ValueError
        If a pattern has an unknown type, the wrong number of multipliers for its type, or a
        multiplier that is not a positive finite number; the message lists every problem.
    """
    resolved, problems = {}, []
    for name, value in (patterns or {}).items():
        if isinstance(value, (tuple, list)) and len(value) == 2 and isinstance(value[0], str):
            kind, values = value[0].upper(), list(value[1])
            if kind not in PATTERN_LENGTHS:
                problems.append(f"{name}: type {value[0]!r} is not one of {list(PATTERN_LENGTHS)}")
                continue
            if len(values) != PATTERN_LENGTHS[kind]:
                problems.append(
                    f"{name}: a {kind} pattern needs {PATTERN_LENGTHS[kind]} multipliers, "
                    f"not {len(values)}"
                )
                continue
        else:
            values = list(value)
            kind = PATTERN_TYPES_BY_LENGTH.get(len(values))
            if kind is None:
                problems.append(
                    f"{name}: {len(values)} multipliers; give 24 (HOURLY), 7 (DAILY) or "
                    "12 (MONTHLY), or a (type, multipliers) pair"
                )
                continue
        if not all(_is_finite_number(v) and float(v) > 0 for v in values):
            problems.append(f"{name}: multipliers must be positive finite numbers")
            continue
        values = tuple(float(v) for v in values)
        mean = sum(values) / len(values)
        if kind in ("HOURLY", "WEEKEND") and abs(mean - 1) > 0.01:
            warnings.warn(
                f"Pattern {name}: the mean of the {kind} multipliers is {mean:.3f}, not 1, so the "
                "mean dry-weather flow is the base flow scaled by it.",
                DesignWarning,
                stacklevel=3,
            )
        resolved[name] = (kind, values)
    if problems:
        raise ValueError(Problems("Invalid patterns", problems).message())
    return resolved


def resolve_dwf(subnetwork, dwf, outfalls=()):
    """{node: base flow} from the ``dwf`` option of :func:`build_model`.

    Raises
    ------
    ValueError
        If a node is not in the subnetwork or is an outfall, or a flow is negative or not finite.
    """
    if not dwf:
        return {}
    dwf = subnetwork.match_nodes(dwf, "dwf")
    nodes = set(subnetwork.nodes)
    missing = sorted(str(n) for n in dwf if n not in nodes)
    if missing:
        raise ValueError(f"Dry-weather flow nodes not in the subnetwork: {', '.join(missing)}")
    at_outfalls = sorted(str(n) for n in dwf if n in set(outfalls))
    if at_outfalls:
        raise ValueError(f"Dry-weather flows cannot enter at outfalls: {', '.join(at_outfalls)}")
    bad = sorted(
        f"{n} ({q!r})" for n, q in dwf.items() if not (_is_finite_number(q) and float(q) >= 0)
    )
    if bad:
        raise ValueError(f"Dry-weather flows must be non-negative and finite: {', '.join(bad)}")
    return {node: float(flow) for node, flow in dwf.items()}


def resolve_dwf_patterns(dwf_patterns, dwf, patterns, subnetwork=None):
    """{node: [pattern names]} from the ``dwf_patterns`` option of :func:`build_model`.

    ``dwf_patterns`` is a pattern name or list of names for every node of ``dwf``, or a mapping
    {node: name or list of names}.

    Raises
    ------
    ValueError
        If a node has no dry-weather flow, a pattern is not defined, or a node has more than
        :data:`MAX_DWF_PATTERNS` patterns or two of the same type.
    """
    if not dwf_patterns:
        return {}
    if isinstance(dwf_patterns, str) or not hasattr(dwf_patterns, "items"):
        dwf_patterns = dict.fromkeys(dwf, dwf_patterns)
    elif subnetwork is not None:
        dwf_patterns = subnetwork.match_nodes(dwf_patterns, "dwf_patterns")
    unloaded = sorted(str(n) for n in dwf_patterns if n not in dwf)
    if unloaded:
        raise ValueError(f"dwf_patterns for nodes with no dry-weather flow: {', '.join(unloaded)}")
    resolved, problems = {}, []
    for node, names in dwf_patterns.items():
        names = [names] if isinstance(names, str) else list(names)
        undefined = [n for n in names if n not in patterns]
        if undefined:
            problems.append(f"{node}: patterns not defined in patterns: {', '.join(undefined)}")
            continue
        types = [patterns[n][0] for n in names]
        if len(names) > MAX_DWF_PATTERNS or len(set(types)) != len(types):
            problems.append(
                f"{node}: at most one pattern of each type ({', '.join(PATTERN_LENGTHS)}), "
                f"not {', '.join(f'{n} ({t})' for n, t in zip(names, types, strict=True))}"
            )
            continue
        resolved[node] = names
    if problems:
        raise ValueError(Problems("Invalid dwf_patterns", problems).message())
    return resolved


def dwf_from_population(
    population: Mapping[NodeId, float],
    per_capita: float,
    unit: str = "L/d",
    flow_units: str = "LPS",
) -> dict[NodeId, float]:
    """{node: base flow} in ``flow_units`` from {node: population} and a per-capita flow.

    How the population is allocated to nodes is left to the caller.

    Parameters
    ----------
    population : dict
        {node: number of people}.

    per_capita : float
        Flow per person, in ``unit``.

    unit : {'L/d', 'm3/d', 'L/s', 'm3/s'}
        Unit of ``per_capita`` (default L/d, e.g. 250 L/person/day).

    flow_units : {'LPS', 'CMS', 'MLD'}
        FLOW_UNITS of the template the flows are written with (default LPS, as in the shipped
        template).

    Returns
    -------
    dict
        {node: base flow}, to pass as ``dwf``.
    """
    volume, _, time = unit.partition("/")
    if volume not in _VOLUME_LITRES or time not in _TIME_SECONDS:
        raise ValueError(f"unit must be one of 'L/d', 'm3/d', 'L/s' or 'm3/s', not {unit!r}.")
    flow_units = flow_units.upper()
    if flow_units not in _FLOW_UNIT_LPS:
        raise ValueError(
            f"flow_units must be one of {', '.join(SI_FLOW_UNITS)}, not {flow_units!r}."
        )
    lps = per_capita * _VOLUME_LITRES[volume] / _TIME_SECONDS[time]
    scale = lps / _FLOW_UNIT_LPS[flow_units]
    return {node: people * scale for node, people in population.items()}
