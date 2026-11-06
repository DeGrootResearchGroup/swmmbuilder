"""The SWMM records of a model and their rendering into an .inp template."""

from __future__ import annotations

import math
import numbers
from dataclasses import dataclass, field

import pandas as pd

# Decimals that numbers are rounded to in the .inp file: 1 micrometre for lengths in metres.
DECIMALS = 6
# Significant figures that flows are written to, whatever the FLOW_UNITS: a flow in CMS can be far
# below 10**-DECIMALS.
FLOW_DIGITS = 10


def _number(value):
    """Text of a number for the .inp file: ints as is, floats rounded to :data:`DECIMALS`
    decimals in their shortest form (so 0.30000000000000004 is written 0.3)."""
    if isinstance(value, bool):
        raise TypeError(f"Expected a number, not {value!r}.")
    if isinstance(value, numbers.Integral):
        return str(int(value))
    value = round(float(value), DECIMALS)
    return repr(value + 0.0)  # + 0.0 turns -0.0 into 0.0


def _flow(value):
    """Text of a flow for the .inp file: :data:`FLOW_DIGITS` significant figures in their
    shortest form, so that small flows in CMS keep their value."""
    if isinstance(value, bool):
        raise TypeError(f"Expected a number, not {value!r}.")
    if isinstance(value, numbers.Integral):
        return str(int(value))
    return repr(float(f"{float(value):.{FLOW_DIGITS}g}") + 0.0)


def _name(value):
    """Text of a SWMM object name; raise ValueError if SWMM could not read it as one token."""
    text = str(value)
    if not text or any(c.isspace() for c in text) or ";" in text:
        raise ValueError(
            f"SWMM name {text!r} is empty or contains whitespace or ';', which SWMM would split "
            "or read as a comment; rename it (e.g. with Network's normalize_ids)."
        )
    return text


def _line(*fields):
    """One .inp line of (value, width) fields: left-aligned in columns of at least ``width``
    characters, and always separated by a space."""
    cells = [(_number(v) if isinstance(v, numbers.Number) else str(v)).ljust(w) for v, w in fields]
    return " ".join(cells).rstrip() + "\n"


@dataclass(frozen=True)
class Junction:
    name: str
    invert: float
    max_depth: float
    init_depth: float
    surcharge_depth: float
    ponding_area: float


@dataclass(frozen=True)
class Outfall:
    name: str
    invert: float
    type: str


@dataclass(frozen=True)
class Storage:
    """A storage node of constant plan area."""

    name: str
    invert: float
    max_depth: float
    init_depth: float
    area: float


@dataclass(frozen=True)
class Conduit:
    name: str
    inlet: str
    outlet: str
    length: float
    roughness: float
    inlet_offset: float
    outlet_offset: float
    init_flow: float
    max_flow: float


@dataclass(frozen=True)
class PumpLink:
    name: str
    inlet: str
    outlet: str
    curve: str
    status: str
    startup: float
    shutoff: float


@dataclass(frozen=True)
class Curve:
    name: str
    type: str
    points: tuple


@dataclass(frozen=True)
class XSection:
    name: str
    shape: str
    geom1: float
    geom2: float
    barrels: float


@dataclass(frozen=True)
class DryWeatherFlow:
    """A node's dry-weather flow: base flow (in the FLOW_UNITS) and time pattern names."""

    node: str
    base: float
    patterns: tuple = ()


@dataclass(frozen=True)
class Pattern:
    name: str
    type: str
    multipliers: tuple


@dataclass(frozen=True)
class ExternalInflow:
    """A node's external flow ([INFLOWS]): ``scale`` x the time series ``timeseries`` (None for
    none) plus ``baseline`` x the pattern ``baseline_pattern`` (None for none), in the
    FLOW_UNITS."""

    node: str
    timeseries: str | None
    scale: float
    baseline: float
    baseline_pattern: str | None = None


@dataclass(frozen=True)
class TimeSeries:
    """A time series ([TIMESERIES]): (hours from the start of the run, value) points."""

    name: str
    points: tuple


@dataclass
class SwmmModel:
    """The SWMM objects of a subnetwork, ready to render into an .inp template.

    Attributes
    ----------
    junctions, outfalls, storage, conduits, pumps, xsections, curves, dwf, patterns : list
        :class:`Junction`, :class:`Outfall`, :class:`Storage`, :class:`Conduit`,
        :class:`PumpLink`, :class:`XSection`, :class:`Curve`, :class:`DryWeatherFlow`, and
        :class:`Pattern` records, in the order they are written.
    inflows, timeseries : list
        :class:`ExternalInflow` and :class:`TimeSeries` records, in the order they are written.
    coordinates : dict
        {node: (x, y)}, in the units of the GIS geometry.
    vertices : dict
        {link: [(x, y), ...]}: the interior vertices of each conduit or pump.
    data_gaps : pandas.DataFrame or None
        The gap report of :func:`swmmbuilder.gaps.fill_gaps`: every gap of the subnetwork's data
        and how it was filled, if it was.
    derived : list
        :class:`Derived` records: every value that comes from a design rule (e.g.
        ``wet_well_defaults``), not from the data.
    """

    junctions: list = field(default_factory=list)
    outfalls: list = field(default_factory=list)
    storage: list = field(default_factory=list)
    conduits: list = field(default_factory=list)
    pumps: list = field(default_factory=list)
    xsections: list = field(default_factory=list)
    curves: list = field(default_factory=list)
    dwf: list = field(default_factory=list)
    patterns: list = field(default_factory=list)
    inflows: list = field(default_factory=list)
    timeseries: list = field(default_factory=list)
    coordinates: dict = field(default_factory=dict)
    vertices: dict = field(default_factory=dict)
    data_gaps: pd.DataFrame | None = None
    derived: list = field(default_factory=list)

    def derived_report(self) -> pd.DataFrame:
        """DataFrame of :attr:`derived`, one row per value (element, id, attribute, value, rule,
        basis)."""
        columns = ["element", "id", "attribute", "value", "rule", "basis"]
        return pd.DataFrame([vars(d) for d in self.derived], columns=columns)

    def summary(self) -> dict[str, float | int]:
        """Counts of the model's records, gaps, and derived values, to cite in a report.

        Returns
        -------
        dict
            Plain Python numbers (``int`` for counts, ``float`` for values, NaN where there is
            nothing to give), under these keys:

            - ``junctions``, ``outfalls``, ``storage``, ``conduits``, ``pumps``: records
              written.
            - ``stations``: pump stations, the nodes that pumps draw from.
            - ``outfall_links``: outfalls added by the ``outfall_link`` option of
              :func:`~swmmbuilder.swmm.build_model`, found by their names (an outfall
              ``<node>_OUT`` reached from ``<node>`` by the conduit ``<node>_DUMMY``).
            - ``data_gaps``: gaps found in the data (rows of :attr:`data_gaps`; 0 without it).
            - ``data_gaps_filled``: those filled by a rule.
            - ``derived``: values derived by a design rule, not from the data (:attr:`derived`).
            - ``derived_<element>_<attribute>``: those of each element and attribute, by any
              rule, e.g. ``derived_junction_max_depth``, ``derived_station_invert`` (stations
              lowered by ``station_invert``), or ``derived_storage_area``.
            - ``derived_<element>_<attribute>[<rule>]``: those by each rule, with the rule as
              :meth:`derived_report` gives it, e.g.
              ``derived_junction_max_depth[max_depth='median']``.
            - ``default_depth_junctions``: junctions given the ``max_depth`` rule's default
              depth (no usable rim).
            - ``default_depth``: that depth, in m (NaN if no junction has it, or if they were
              given different depths).

            The derived keys come from the rules used, so a model with none has only
            ``derived``.

        Examples
        --------
        ::

            s = model.summary()
            print(f"{s['data_gaps']} gaps ({s['data_gaps_filled']} filled)")
        """
        gaps = self.data_gaps
        derived = self.derived_report()
        outfalls = {o.name for o in self.outfalls}
        links = sum(
            1
            for c in self.conduits
            if c.name == f"{c.inlet}_DUMMY"
            and c.outlet == f"{c.inlet}_OUT"
            and c.outlet in outfalls
        )
        result: dict[str, float | int] = {
            "junctions": len(self.junctions),
            "outfalls": len(self.outfalls),
            "storage": len(self.storage),
            "conduits": len(self.conduits),
            "pumps": len(self.pumps),
            "stations": len({p.inlet for p in self.pumps}),
            "outfall_links": links,
            "data_gaps": 0 if gaps is None else len(gaps),
            "data_gaps_filled": 0 if gaps is None else int(gaps["filled_from"].notna().sum()),
            "derived": len(derived),
        }
        for (element, attribute), group in derived.groupby(["element", "attribute"], sort=True):
            result[f"derived_{element}_{attribute}"] = len(group)
            for rule, n in group["rule"].value_counts(sort=False).sort_index().items():
                result[f"derived_{element}_{attribute}[{rule}]"] = int(n)
        default = derived[
            (derived["element"] == "junction")
            & (derived["attribute"] == "max_depth")
            & derived["rule"].astype(str).str.startswith("max_depth")
        ]
        depths = default["value"].astype(float).unique()
        result["default_depth_junctions"] = len(default)
        result["default_depth"] = float(depths[0]) if len(depths) == 1 else math.nan
        return result

    def sections(self) -> dict[str, str]:
        """{placeholder: text} for every placeholder in :data:`TEMPLATE_PLACEHOLDERS`.

        Fields are padded into columns for readability, but always separated by at least one
        space, so a long name or number cannot run into the next field. Numbers are rounded to
        :data:`DECIMALS` decimals, and flows (dry-weather, inflow, time-series, pump-curve, and
        conduit flows) to :data:`FLOW_DIGITS` significant figures.

        Raises
        ------
        ValueError
            If a name is empty or contains whitespace or ';', which SWMM would split or read as
            a comment.
        """
        return {
            "nodes": "".join(
                _line(
                    (_name(j.name), 16),
                    (j.invert, 11),
                    (j.max_depth, 11),
                    (j.init_depth, 11),
                    (j.surcharge_depth, 11),
                    (j.ponding_area, 11),
                )
                for j in self.junctions
            ),
            "outfalls": "".join(
                _line((_name(o.name), 16), (o.invert, 11), (o.type, 28), ("NO", 10))
                for o in self.outfalls
            ),
            # FUNCTIONAL area = 0 * depth^0 + area: a constant plan area.
            "storage": "".join(
                _line(
                    (_name(s.name), 16),
                    (s.invert, 11),
                    (s.max_depth, 11),
                    (s.init_depth, 11),
                    ("FUNCTIONAL", 12),
                    (0, 11),
                    (0, 11),
                    (s.area, 11),
                    (0, 11),
                    (0, 11),
                )
                for s in self.storage
            ),
            "conduits": "".join(
                _line(
                    (_name(c.name), 16),
                    (_name(c.inlet), 16),
                    (_name(c.outlet), 16),
                    (c.length, 11),
                    (c.roughness, 11),
                    (c.inlet_offset, 11),
                    (c.outlet_offset, 11),
                    (_flow(c.init_flow), 11),
                    (_flow(c.max_flow), 11),
                )
                for c in self.conduits
            ),
            "pumps": "".join(
                _line(
                    (_name(p.name), 16),
                    (_name(p.inlet), 16),
                    (_name(p.outlet), 16),
                    (p.curve if p.curve == "*" else _name(p.curve), 16),
                    (p.status, 6),
                    (p.startup, 11),
                    (p.shutoff, 11),
                )
                for p in self.pumps
            ),
            "xsections": "".join(
                _line(
                    (_name(x.name), 16),
                    (x.shape, 16),
                    (x.geom1, 16),
                    (x.geom2, 11),
                    (0.0, 11),
                    (0.0, 11),
                    (x.barrels, 11),
                )
                for x in self.xsections
            ),
            "curves": "".join(
                "".join(
                    _line(
                        (_name(c.name), 16), (c.type if i == 0 else "", 9), (x, 11), (_flow(y), 11)
                    )
                    for i, (x, y) in enumerate(c.points)
                )
                + "\n"
                for c in self.curves
            ),
            "dwf": "".join(
                _line(
                    (_name(d.node), 16),
                    ("FLOW", 16),
                    (_flow(d.base), 11),
                    *((f'"{_name(p)}"', 0) for p in d.patterns),
                )
                for d in self.dwf
            ),
            # Six multipliers a line; the type is on the first line only.
            "patterns": "".join(
                "".join(
                    _line(
                        (_name(p.name), 16),
                        (p.type if i == 0 else "", 10),
                        *((m, 8) for m in p.multipliers[i : i + 6]),
                    )
                    for i in range(0, len(p.multipliers), 6)
                )
                + "\n"
                for p in self.patterns
            ),
            # Constituent FLOW, type FLOW, units factor 1; "" for no time series.
            "inflows": "".join(
                _line(
                    (_name(i.node), 16),
                    ("FLOW", 16),
                    ('""' if i.timeseries is None else _name(i.timeseries), 16),
                    ("FLOW", 8),
                    (1.0, 8),
                    (i.scale, 8),
                    (_flow(i.baseline), 8),
                    *(() if i.baseline_pattern is None else ((_name(i.baseline_pattern), 0),)),
                )
                for i in self.inflows
            ),
            # Times in decimal hours from the start of the run.
            "timeseries": "".join(
                "".join(_line((_name(s.name), 16), (t, 11), (_flow(v), 11)) for t, v in s.points)
                + "\n"
                for s in self.timeseries
            ),
            "coordinates": "".join(
                f"{_name(node):<16} {x:<18} {y:<18}\n" for node, (x, y) in self.coordinates.items()
            ),
            "vertices": "".join(
                f"{_name(name):<16} {x:<18} {y:<18}\n"
                for name, points in self.vertices.items()
                for x, y in points
            ),
        }

    def render(self, template: str) -> str:
        """The text of ``template`` with its placeholders filled."""
        for name, text in self.sections().items():
            template = template.replace(f"{{{name}}}", text)
        return template
