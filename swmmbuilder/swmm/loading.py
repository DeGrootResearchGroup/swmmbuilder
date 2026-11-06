"""The dry-weather load of a model's nodes, from its parts: sanitary flow from population, flows
that each follow their own hourly pattern (such as land uses), and steady flows (such as
groundwater infiltration).

SWMM gives a node one [DWF] base flow with at most one hourly pattern, so only one kind of
patterned flow can go there. :class:`NodeLoading` puts the sanitary flow in [DWF] with its
diurnal pattern (:func:`~swmmbuilder.dwf_from_population`), every other patterned flow in
[INFLOWS] as time series (:func:`~swmmbuilder.pattern_inflows`), and the steady flow in the
[INFLOWS] baseline. :meth:`NodeLoading.options` gives the matching options of
:func:`~swmmbuilder.swmm.build_model`, and :meth:`NodeLoading.table` the load of each node by
part, for checking the totals against a plant's measured dry-weather flow.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from ..diagnostics import Problems
from ..units import _is_finite_number
from .dwf import dwf_from_population
from .inflows import Inflow, pattern_inflows, pattern_mean_problem, unit_mean_tolerance

if TYPE_CHECKING:
    from .._types import NodeId

# Name of the [DWF] pattern when it is given as multipliers.
DWF_PATTERN_NAME = "diurnal"


def _node_values(values: Mapping[NodeId, float] | pd.Series, what: str) -> pd.Series:
    """{node: value} as a float Series, checked to be non-negative and finite."""
    series = pd.Series(values, dtype=float)
    bad = [str(n) for n, x in series.items() if not (np.isfinite(x) and x >= 0)]
    if bad:
        raise ValueError(f"{what} must be non-negative and finite: {', '.join(bad[:20])}")
    return series


class NodeLoading:
    """The dry-weather load of each node, as the options of
    :func:`~swmmbuilder.swmm.build_model` and as a table.

    Three parts, each optional (at least one is needed):

    - ``population`` at each node, times ``per_capita``, is the sanitary flow, written to [DWF]
      (:func:`~swmmbuilder.dwf_from_population`) with the hourly ``pattern``.
    - ``patterned`` flows, one column per kind (e.g. the land uses of
      :attr:`LandUseAllocation.nodes <swmmbuilder.LandUseAllocation.nodes>`), each follow the
      hourly pattern of their kind in ``patterns`` and are written to [INFLOWS] and
      [TIMESERIES] (:func:`~swmmbuilder.pattern_inflows`).
    - ``steady`` flows (e.g. groundwater infiltration, from
      :func:`~swmmbuilder.distribute_flow`) are the [INFLOWS] baseline.

    Nodes with a population or flow of 0 are not written.

    Parameters
    ----------
    population : dict or pandas.Series, optional
        {node: number of people}, e.g. from
        :meth:`PopulationAllocation.within <swmmbuilder.PopulationAllocation.within>`.

    per_capita : float, optional
        Flow per person, in ``unit``; needed with ``population``.

    unit : {'L/d', 'm3/d', 'L/s', 'm3/s'}
        Unit of ``per_capita`` (default L/d).

    pattern : sequence of float or str, optional
        Hourly pattern of the sanitary flow: 24 multipliers from 00:00, written as the pattern
        ``'diurnal'``, or the name of a pattern in ``patterns``. None (default) writes no pattern,
        so the sanitary flow is constant.

    patterned : pandas.DataFrame or dict, optional
        Mean flow at each node by kind: a DataFrame of nodes (index) by kinds (columns), or
        {kind: {node: flow}}, in ``flow_units``.

    patterns : dict, optional
        {name: 24 hourly multipliers, from 00:00}: one for every kind of ``patterned`` with a
        flow, and the sanitary pattern if ``pattern`` names one.

    days : int, optional
        Number of days the time series of ``patterned`` cover (at least the length of the run);
        needed with ``patterned``.

    steady : dict or pandas.Series, optional
        {node: steady flow}, in ``flow_units``.

    prefix : str
        Prefix of the time series names (default none), e.g. 'NR_'; see
        :func:`~swmmbuilder.pattern_inflows`.

    ramp_minutes : float
        See :func:`~swmmbuilder.hourly_series`.

    require_unit_mean : bool or float
        Refuse an hourly pattern (the sanitary one, or one of ``patterned``) whose mean is not 1:
        True (default) for a tolerance of 0.001, a tolerance, or False to accept any mean. A
        pattern's mean scales the mean flow, so that the table would not match the model.

    flow_units : {'LPS', 'CMS', 'MLD'}
        FLOW_UNITS of the template (default LPS, the shipped template's): the unit of the
        sanitary flow, and of ``patterned`` and ``steady``.

    Raises
    ------
    ValueError
        If no part is given, ``population`` comes without ``per_capita``, ``patterned`` without
        ``days``, or ``pattern`` without ``population``; a pattern is missing or its mean is not
        1; or a population or flow is negative or not finite.

    Examples
    --------
    >>> loading = NodeLoading(
    ...     population=alloc.within(sub), per_capita=230, unit="L/d", pattern=diurnal,
    ...     patterned=land.within(sub), patterns=land_use_patterns, days=3, prefix="NR_",
    ...     steady=infiltration,
    ... )  # doctest: +SKIP
    >>> sub.to_inp("model.inp", ["OUT"], **loading.options())  # doctest: +SKIP
    """

    def __init__(
        self,
        population: Mapping[NodeId, float] | pd.Series | None = None,
        per_capita: float | None = None,
        unit: str = "L/d",
        pattern: Sequence[float] | str | None = None,
        patterned: pd.DataFrame | Mapping[Any, Mapping[NodeId, float]] | None = None,
        patterns: Mapping[Any, Sequence[float]] | None = None,
        days: int | None = None,
        steady: Mapping[NodeId, float] | pd.Series | None = None,
        prefix: str = "",
        ramp_minutes: float = 1.0,
        require_unit_mean: bool | float = True,
        flow_units: str = "LPS",
    ) -> None:
        if population is None and patterned is None and steady is None:
            raise ValueError("NodeLoading needs population, patterned, or steady flows.")
        patterns = dict(patterns or {})
        tolerance = unit_mean_tolerance(require_unit_mean)
        problems = Problems("Invalid node loading")
        if population is not None and per_capita is None:
            problems.add("population needs per_capita")
        if pattern is not None and population is None:
            problems.add("pattern is the sanitary flow's, and needs population")
        if patterned is not None and days is None:
            problems.add("patterned needs days")
        if isinstance(pattern, str) and pattern not in patterns:
            problems.add(f"pattern {pattern!r} is not in patterns")
        problems.raise_if_any()

        self.flow_units = flow_units
        self.population = (
            _node_values(population, "populations") if population is not None else None
        )
        self.steady = _node_values(steady, "steady flows") if steady is not None else None
        if patterned is None:
            self.patterned = None
        else:
            table = patterned if isinstance(patterned, pd.DataFrame) else pd.DataFrame(patterned)
            self.patterned = table.fillna(0.0).astype(float)

        # [DWF] and its pattern.
        self.dwf: dict[NodeId, float] = {}
        self.dwf_pattern: str | None = None
        self.patterns: dict[str, list[float]] = {}
        if self.population is not None:
            assert per_capita is not None  # checked above
            people = self.population[self.population > 0]
            self.dwf = dwf_from_population(
                people.to_dict(), per_capita=float(per_capita), unit=unit, flow_units=flow_units
            )
            if pattern is not None:
                name = pattern if isinstance(pattern, str) else DWF_PATTERN_NAME
                values = patterns[name] if isinstance(pattern, str) else pattern
                if not all(_is_finite_number(v) for v in values) or len(values) != 24:
                    raise ValueError(f"The sanitary pattern {name!r} needs 24 finite multipliers.")
                problem = (
                    pattern_mean_problem(name, values, tolerance) if tolerance is not None else None
                )
                if problem:
                    raise ValueError(Problems("Invalid node loading", [problem]).message())
                self.dwf_pattern = name
                self.patterns = {name: [float(v) for v in values]}

        # [INFLOWS] and [TIMESERIES].
        self.inflows: dict[NodeId, Inflow] = {}
        self.timeseries: dict[str, list[tuple[float, float]]] = {}
        if self.patterned is not None or self.steady is not None:
            flows = self.patterned if self.patterned is not None else pd.DataFrame(dtype=float)
            kinds = {k: patterns[k] for k in flows.columns if k in patterns}
            baseline = self.steady[self.steady > 0].to_dict() if self.steady is not None else None
            self.inflows, self.timeseries = pattern_inflows(
                flows,
                kinds,
                days if days is not None else 1,
                baseline=baseline,
                prefix=prefix,
                ramp_minutes=ramp_minutes,
                require_unit_mean=require_unit_mean,
            )

    def options(self) -> dict[str, Any]:
        """The options of :func:`~swmmbuilder.swmm.build_model` for this load: ``dwf``,
        ``patterns``, ``dwf_patterns``, ``inflows``, and ``timeseries``, each only when there is
        something to write, so that ``sub.to_inp(path, outfalls, **loading.options())`` works.

        Returns
        -------
        dict
            Keyword arguments for :func:`~swmmbuilder.swmm.build_model` (or
            :meth:`Subnetwork.to_inp <swmmbuilder.Subnetwork.to_inp>`).
        """
        options: dict[str, Any] = {}
        if self.dwf:
            options["dwf"] = dict(self.dwf)
            if self.dwf_pattern is not None:
                options["patterns"] = {k: list(v) for k, v in self.patterns.items()}
                options["dwf_patterns"] = self.dwf_pattern
        if self.inflows:
            options["inflows"] = dict(self.inflows)
        if self.timeseries:
            options["timeseries"] = {k: list(v) for k, v in self.timeseries.items()}
        return options

    def table(self, columns: Mapping[str, str] | None = None) -> pd.DataFrame:
        """The load of each node by part.

        Parameters
        ----------
        columns : dict, optional
            {column: new name}, to rename the columns, e.g. ``{'steady': 'infiltration_L_s'}``.

        Returns
        -------
        pandas.DataFrame
            Indexed by ``node`` (every node given a population or flow, including those of 0),
            with ``population`` (people), ``sanitary`` (the [DWF] base flow), ``steady``, and one
            column per kind of ``patterned`` (mean flows), in ``flow_units``; 0 where a node has
            none of a part.
        """
        empty = pd.Series(dtype=float)
        parts = {
            "population": self.population if self.population is not None else empty,
            "sanitary": pd.Series(self.dwf, dtype=float),
            "steady": self.steady if self.steady is not None else empty,
        }
        if self.patterned is not None:
            parts.update({k: self.patterned[k] for k in self.patterned.columns})
        table = pd.DataFrame(parts).fillna(0.0).rename_axis("node")
        return table.rename(columns=dict(columns)) if columns else table

    def summary(self) -> dict[str, float]:
        """Totals and node counts by part, for the build log or a check against plant flow.

        Returns
        -------
        dict
            ``population`` (people) and ``population_nodes``; ``sanitary``, ``steady``,
            ``patterned`` (all kinds), and ``total`` (the three) flows in ``flow_units``, each
            with the number of nodes carrying it (``sanitary_nodes``, ``steady_nodes``,
            ``patterned_nodes``, ``loaded_nodes``); and for each kind of ``patterned``,
            ``patterned_<kind>`` and ``patterned_<kind>_nodes``.
        """
        table = self.table()
        kinds = list(self.patterned.columns) if self.patterned is not None else []
        patterned = table[kinds].sum(axis=1)
        flows = table["sanitary"] + table["steady"] + patterned
        out: dict[str, float] = {}
        for name, values in [
            ("population", table["population"]),
            ("sanitary", table["sanitary"]),
            ("steady", table["steady"]),
            ("patterned", patterned),
        ]:
            out[name] = float(values.sum())
            out[f"{name}_nodes"] = int((values > 0).sum())
        out["total"] = float(flows.sum())
        out["loaded_nodes"] = int((flows > 0).sum())
        for kind in kinds:
            out[f"patterned_{kind}"] = float(table[kind].sum())
            out[f"patterned_{kind}_nodes"] = int((table[kind] > 0).sum())
        return out

    def __repr__(self) -> str:
        s = self.summary()
        return (
            f"NodeLoading(sanitary={s['sanitary']:.4g} at {s['sanitary_nodes']} nodes, "
            f"patterned={s['patterned']:.4g} at {s['patterned_nodes']}, "
            f"steady={s['steady']:.4g} at {s['steady_nodes']}, {self.flow_units})"
        )
