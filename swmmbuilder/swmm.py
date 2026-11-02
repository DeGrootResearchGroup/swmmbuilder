"""SWMM models of a subnetwork, and their .inp files.

:func:`build_model` turns a subnetwork into a :class:`SwmmModel`: plain records of the junctions,
outfalls, storage nodes, conduits, pumps, cross-sections, curves, dry-weather flows, time patterns,
coordinates and vertices, with all the hydraulic decisions
(node inverts, conduit offsets, cross-sections, Manning's n) already made.
:meth:`SwmmModel.render` fills an .inp template with them, and :func:`write_inp` checks the
template and the topology, builds the model and writes the file.
"""

import logging
import math
import re
import warnings
from dataclasses import dataclass, field, replace
from importlib import resources
from pathlib import Path

import pandas as pd

from .gaps import _invert, fill_gaps, lowest_inverts, xsection
from .geometry import line_coords

logger = logging.getLogger(__name__)

# SWMM's SI flow units. With any of these, lengths, elevations and section sizes are in metres,
# matching the GIS data; with a US unit (CFS, GPM, MGD) SWMM would read them as feet.
SI_FLOW_UNITS = ("CMS", "LPS", "MLD")

# Placeholders that SwmmModel.render fills in the template.
TEMPLATE_PLACEHOLDERS = (
    "nodes",
    "outfalls",
    "storage",
    "conduits",
    "pumps",
    "xsections",
    "curves",
    "dwf",
    "patterns",
    "coordinates",
    "vertices",
)

# SWMM time pattern types and their number of multipliers. A pattern given as a plain list of
# multipliers takes its type from its length; WEEKEND (24 values) must be named explicitly.
PATTERN_LENGTHS = {"MONTHLY": 12, "DAILY": 7, "HOURLY": 24, "WEEKEND": 24}
PATTERN_TYPES_BY_LENGTH = {12: "MONTHLY", 7: "DAILY", 24: "HOURLY"}

# SWMM's limit on the time patterns of a node's dry-weather flow (one of each type).
MAX_DWF_PATTERNS = 4

# Litres per unit of volume, and seconds per unit of time, for dwf_from_population.
_VOLUME_LITRES = {"L": 1.0, "m3": 1000.0}
_TIME_SECONDS = {"s": 1.0, "d": 86400.0}
# Litres per second in each SI flow unit.
_FLOW_UNIT_LPS = {"LPS": 1.0, "CMS": 1000.0, "MLD": 1e6 / 86400.0}

# SWMM pump curve types, and what the x values of each are (y is flow, in the FLOW_UNITS).
PUMP_CURVE_TYPES = {
    "PUMP1": "wet well volume",
    "PUMP2": "inlet node depth",
    "PUMP3": "head difference",
    "PUMP4": "inlet node depth",
    "PUMP5": "head difference",
}


def default_template_path():
    """Path-like (Traversable) to the SI template shipped with the package."""
    return resources.files("swmmbuilder") / "templates" / "template.inp"


def read_template(path=None):
    """Text of the .inp template at ``path`` (default: :func:`default_template_path`)."""
    if path is None:
        return default_template_path().read_text()
    if not Path(path).exists():
        raise FileNotFoundError(f"The template file does not exist at the specified path: {path}")
    return Path(path).read_text()


def check_template(content):
    """Check that a .inp template is one the writer can fill.

    Parameters
    ----------
    content : str
        Text of the .inp template.

    Raises
    ----------
    ValueError
        If FLOW_UNITS is missing or not an SI unit (the GIS data are metric), or if the template
        contains a placeholder the writer does not fill (see :data:`TEMPLATE_PLACEHOLDERS`).
    """
    flow_units = template_flow_units(content)
    if flow_units is None:
        raise ValueError(
            "The template has no FLOW_UNITS option (SWMM would default to CFS, US units). "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )
    if flow_units not in SI_FLOW_UNITS:
        raise ValueError(
            f"The template sets FLOW_UNITS {flow_units}, a US unit, so SWMM would read the metric "
            "GIS lengths, elevations and sizes as feet. "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )

    unknown = sorted(set(re.findall(r"\{(\w+)\}", content)) - set(TEMPLATE_PLACEHOLDERS))
    if unknown:
        raise ValueError(
            f"The template contains placeholders that are not filled: {unknown}. "
            f"Supported placeholders: {list(TEMPLATE_PLACEHOLDERS)}."
        )


def template_flow_units(template):
    """FLOW_UNITS of the text of an .inp template, upper case; None if it is not set."""
    match = re.search(r"^\s*FLOW_UNITS\s+(\S+)", template, re.MULTILINE | re.IGNORECASE)
    return match.group(1).upper() if match else None


def template_routing(template):
    """FLOW_ROUTING of the text of an .inp template; SWMM's default, KINWAVE, if it is not set."""
    match = re.search(r"^\s*FLOW_ROUTING\s+(\w+)", template, flags=re.MULTILINE | re.IGNORECASE)
    return match.group(1).upper() if match else "KINWAVE"


def check_topology(report, routing, dead_ends):
    """Raise if a subnetwork's :class:`TopologyReport` cannot be written with ``routing``.

    Flow splits (including parallel pipes) and cycles need DYNWAVE; dead ends are refused unless
    ``dead_ends`` is 'outfall'. Every problem is listed in one ``ValueError``.
    """
    problems = []
    if routing != "DYNWAVE":
        for node, pipes in report.splits.items():
            problems.append(
                f"node {node} has more than one outlet pipe ({', '.join(map(str, pipes))})"
            )
        for cycle in report.cycles:
            problems.append(f"cycle {' -> '.join(cycle + cycle[:1])}")
    if dead_ends == "raise" and report.dead_ends:
        problems.append(
            f"dead ends with no outlet that are not outfalls: {', '.join(report.dead_ends)} "
            "(pass dead_ends='outfall' to write them as outfalls)"
        )
    if not problems:
        return
    message = "Invalid subnetwork topology:\n  " + "\n  ".join(problems)
    if routing != "DYNWAVE" and not report.is_dendritic:
        message += (
            f"\n{routing} routing requires a dendritic network; use a template with "
            "FLOW_ROUTING DYNWAVE for flow splits, parallel pipes and loops."
        )
    raise ValueError(message)


@dataclass
class Pump:
    """A pump that replaces a pipe in the .inp file (see the ``pumps`` option of
    :func:`build_model`).

    Attributes
    ----------
    curve : list of (float, float) or str
        Pump curve points (x, y), y being the flow in the template's FLOW_UNITS (L/s with the
        shipped template) and x as given by ``curve_type``; or ``'*'`` for an ideal pump, which
        pumps whatever flows into its inlet node.

    curve_type : str
        PUMP1 (x: wet well volume, m3; the inlet must be a wet well), PUMP2 (x: inlet depth, m,
        stepwise), PUMP3 (x: head difference, m), PUMP4 (x: inlet depth, m, continuous) or PUMP5
        (variable-speed PUMP3). Ignored for an ideal pump.

    status : str
        Initial status, 'ON' or 'OFF'.

    startup, shutoff : float
        Inlet node depths (m) at which the pump turns on and off; 0 for none.
    """

    curve: object
    curve_type: str = "PUMP3"
    status: str = "ON"
    startup: float = 0.0
    shutoff: float = 0.0

    @classmethod
    def coerce(cls, pump_id, value):
        """A Pump from a Pump or a dict of its fields, checked; raise ValueError if it is invalid."""
        pump = replace(value) if isinstance(value, cls) else cls(**value)
        problems = []
        if isinstance(pump.curve, str):
            if pump.curve != "*":
                problems.append(f"curve must be a list of (x, y) points or '*', not {pump.curve!r}")
        else:
            points = [tuple(pt) for pt in pump.curve]
            if not points or any(len(pt) != 2 for pt in points):
                problems.append("curve must be a non-empty list of (x, y) points")
            else:
                pump.curve = [(float(x), float(y)) for x, y in points]
            if pump.curve_type.upper() not in PUMP_CURVE_TYPES:
                problems.append(
                    f"curve_type {pump.curve_type!r} is not one of {list(PUMP_CURVE_TYPES)}"
                )
            pump.curve_type = pump.curve_type.upper()
        pump.status = pump.status.upper()
        if pump.status not in ("ON", "OFF"):
            problems.append(f"status must be 'ON' or 'OFF', not {pump.status!r}")
        if problems:
            raise ValueError(f"Pump {pump_id}: {'; '.join(problems)}.")
        return pump


@dataclass
class WetWell:
    """A node written as a [STORAGE] node of constant plan area, e.g. a pumping station wet well
    (see the ``wet_wells`` option of :func:`build_model`).

    Attributes
    ----------
    max_depth : float
        Depth (m) from the node invert to the top of the wet well.

    area : float
        Plan area (m2).

    init_depth : float
        Initial water depth (m).
    """

    max_depth: float
    area: float
    init_depth: float = 0.0

    @classmethod
    def coerce(cls, node, value):
        """A WetWell from a WetWell or a dict of its fields; raise ValueError if it is invalid."""
        well = replace(value) if isinstance(value, cls) else cls(**value)
        if not (well.max_depth > 0 and well.area > 0 and well.init_depth >= 0):
            raise ValueError(
                f"Wet well {node}: max_depth and area must be positive and init_depth not negative."
            )
        return well


def resolve_pumps(subnetwork, pumps):
    """{pipe_id: Pump} from the ``pumps`` option of :func:`build_model`.

    The keys are matched to the subnetwork pipe IDs as strings, so ``1010`` and ``'1010'`` are
    the same pipe; the result is keyed by the pipe IDs of the data.

    Raises
    ----------
    ValueError
        If a pipe is not in the subnetwork or a pump is invalid.
    """
    if not pumps:
        return {}
    ids = {str(p): p for p in subnetwork.pipes()["pipe_id"]}
    missing = sorted(str(k) for k in pumps if str(k) not in ids)
    if missing:
        raise ValueError(f"Pumped pipes not in the subnetwork: {', '.join(missing)}")
    return {ids[str(k)]: Pump.coerce(k, v) for k, v in pumps.items()}


def resolve_wet_wells(subnetwork, wet_wells, outfalls=()):
    """{node: WetWell} from the ``wet_wells`` option of :func:`build_model`.

    Raises
    ----------
    ValueError
        If a node is not in the subnetwork or is an outfall, or a wet well is invalid.
    """
    if not wet_wells:
        return {}
    missing = sorted(str(n) for n in wet_wells if n not in set(subnetwork.nodes))
    if missing:
        raise ValueError(f"Wet wells not in the subnetwork: {', '.join(missing)}")
    at_outfalls = sorted(str(n) for n in wet_wells if n in set(outfalls))
    if at_outfalls:
        raise ValueError(f"Wet wells cannot be outfalls: {', '.join(at_outfalls)}")
    return {node: WetWell.coerce(node, v) for node, v in wet_wells.items()}


def _is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def resolve_patterns(patterns):
    """{name: (type, multipliers)} from the ``patterns`` option of :func:`build_model`.

    Each pattern is a list of multipliers, whose type is set by its length (24 HOURLY, 7 DAILY,
    12 MONTHLY), or a (type, multipliers) pair. A warning is given for an HOURLY or WEEKEND
    pattern whose mean is not 1, since the mean then scales the base flows.

    Raises
    ----------
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
                stacklevel=3,
            )
        resolved[name] = (kind, values)
    if problems:
        raise ValueError("Invalid patterns:\n  " + "\n  ".join(problems))
    return resolved


def resolve_dwf(subnetwork, dwf, outfalls=()):
    """{node: base flow} from the ``dwf`` option of :func:`build_model`.

    Raises
    ----------
    ValueError
        If a node is not in the subnetwork or is an outfall, or a flow is negative or not finite.
    """
    if not dwf:
        return {}
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


def resolve_dwf_patterns(dwf_patterns, dwf, patterns):
    """{node: [pattern names]} from the ``dwf_patterns`` option of :func:`build_model`.

    ``dwf_patterns`` is a pattern name or list of names for every node of ``dwf``, or a mapping
    {node: name or list of names}.

    Raises
    ----------
    ValueError
        If a node has no dry-weather flow, a pattern is not defined, or a node has more than
        :data:`MAX_DWF_PATTERNS` patterns or two of the same type.
    """
    if not dwf_patterns:
        return {}
    if isinstance(dwf_patterns, str) or not hasattr(dwf_patterns, "items"):
        dwf_patterns = dict.fromkeys(dwf, dwf_patterns)
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
        raise ValueError("Invalid dwf_patterns:\n  " + "\n  ".join(problems))
    return resolved


def dwf_from_population(population, per_capita, unit="L/d", flow_units="LPS"):
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
    ----------
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


def node_inverts(subnetwork, height_correction=None):
    """Invert elevation of every node of a subnetwork.

    A node's invert is the lowest invert of the subnetwork pipes connected to it, including the
    drop-structure inverts of the pipes entering it, so that every pipe sits at or above its
    nodes (the pipes are then attached with offsets). A node with no subnetwork pipe takes the
    lowest invert of the network pipes connected to it.

    Parameters
    ----------
    subnetwork : Subnetwork
        The subnetwork.

    height_correction : dict, optional
        {node: invert} overriding the computed inverts.

    Returns
    ----------
    dict
        {node: invert}

    Raises
    ----------
    ValueError
        If a node has no pipe with an invert.
    """
    nodes = set(subnetwork.nodes)
    inverts = lowest_inverts(subnetwork.pipes(), nodes)
    if nodes - set(inverts):
        inverts.update(lowest_inverts(subnetwork.network.full_data, nodes - set(inverts)))
    missing = sorted(nodes - set(inverts))
    if missing:
        raise ValueError(f"No pipe invert found for nodes: {', '.join(map(str, missing))}")
    for node, invert in (height_correction or {}).items():
        if node in nodes:
            inverts[node] = float(invert)
    return inverts


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


@dataclass
class SwmmModel:
    """The SWMM objects of a subnetwork, ready to render into an .inp template.

    Attributes
    ----------
    junctions, outfalls, storage, conduits, pumps, xsections, curves, dwf, patterns : list
        :class:`Junction`, :class:`Outfall`, :class:`Storage`, :class:`Conduit`,
        :class:`PumpLink`, :class:`XSection`, :class:`Curve`, :class:`DryWeatherFlow` and
        :class:`Pattern` records, in the order they are written.
    coordinates : dict
        {node: (x, y)}, in the units of the GIS geometry.
    vertices : dict
        {link: [(x, y), ...]}: the interior vertices of each conduit or pump.
    data_gaps : pandas.DataFrame or None
        The gap report of :func:`swmmbuilder.gaps.fill_gaps`: every gap of the subnetwork's data
        and how it was filled, if it was.
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
    coordinates: dict = field(default_factory=dict)
    vertices: dict = field(default_factory=dict)
    data_gaps: pd.DataFrame | None = None

    def sections(self):
        """{placeholder: text} for every placeholder in :data:`TEMPLATE_PLACEHOLDERS`."""
        return {
            "nodes": "".join(
                f"{j.name:<16} {j.invert:<11} {j.max_depth:<11} {j.init_depth:<11}"
                f"{j.surcharge_depth:<11} {j.ponding_area:<11}\n"
                for j in self.junctions
            ),
            "outfalls": "".join(
                f"{o.name:<16} {o.invert:<11}  {o.type:<28} {'NO':<10}\n" for o in self.outfalls
            ),
            # FUNCTIONAL area = 0 * depth^0 + area: a constant plan area.
            "storage": "".join(
                f"{s.name:<16} {s.invert:<11} {s.max_depth:<11} {s.init_depth:<11} "
                f"{'FUNCTIONAL':<12} {0:<11} {0:<11} {s.area:<11} {0:<11} {0:<11}\n"
                for s in self.storage
            ),
            "conduits": "".join(
                f"{c.name:<16} {c.inlet:<16} {c.outlet:<16}"
                f"{c.length:<11} {c.roughness:<11} "
                f"{c.inlet_offset:<11} {c.outlet_offset:<11}"
                f"{c.init_flow:<11} {c.max_flow:<11}\n"
                for c in self.conduits
            ),
            "pumps": "".join(
                f"{p.name:<16} {p.inlet:<16} {p.outlet:<16} {p.curve:<16} {p.status:<6} "
                f"{p.startup:<11} {p.shutoff:<11}\n"
                for p in self.pumps
            ),
            "xsections": "".join(
                f"{x.name:<16} {x.shape:<16} {x.geom1:<16}{x.geom2:<11} {0.0:<11} {0.0:<11} "
                f"{x.barrels:<11}\n"
                for x in self.xsections
            ),
            "curves": "".join(
                "".join(
                    f"{c.name:<16} {c.type if i == 0 else '':<9} {x:<11} {y:<11}\n"
                    for i, (x, y) in enumerate(c.points)
                )
                + "\n"
                for c in self.curves
            ),
            "dwf": "".join(
                f"{d.node:<16} {'FLOW':<16} {d.base:<11} "
                + " ".join(f'"{p}"' for p in d.patterns)
                + "\n"
                for d in self.dwf
            ),
            # Six multipliers a line; the type is on the first line only.
            "patterns": "".join(
                "".join(
                    f"{p.name:<16} {p.type if i == 0 else '':<10} "
                    + " ".join(f"{m:<8}" for m in p.multipliers[i : i + 6])
                    + "\n"
                    for i in range(0, len(p.multipliers), 6)
                )
                + "\n"
                for p in self.patterns
            ),
            "coordinates": "".join(
                f"{node:<16} {x:<18} {y:<18}\n" for node, (x, y) in self.coordinates.items()
            ),
            "vertices": "".join(
                f"{name:<16} {x:<18} {y:<18}\n"
                for name, points in self.vertices.items()
                for x, y in points
            ),
        }

    def render(self, template):
        """The text of ``template`` with its placeholders filled."""
        for name, text in self.sections().items():
            template = template.replace(f"{{{name}}}", text)
        return template


def outfall_types(outfalls, subnetwork=None):
    """{node: outfall type} from a mapping, or from an iterable of nodes (all 'FREE'); None is
    the subnetwork's outfall nodes (:attr:`Subnetwork.outfalls`)."""
    if outfalls is None:
        outfalls = subnetwork.outfalls
    if isinstance(outfalls, str):
        return {outfalls: "FREE"}
    if hasattr(outfalls, "items"):
        return dict(outfalls)
    return {node: "FREE" for node in outfalls}


def build_model(
    subnetwork,
    outfalls=None,
    max_depth=0.0,
    init_depth=0.0,
    surcharge_depth=0.0,
    ponding_area=0.0,
    manning_n=0.01,
    manning_n_by_material=None,
    init_flow=0.0,
    max_flow=0.0,
    barrels=1.0,
    height_correction=None,
    pumps=None,
    wet_wells=None,
    dwf=None,
    patterns=None,
    dwf_patterns=None,
    **fills,
):
    """The :class:`SwmmModel` of a subnetwork.

    Gaps in the data are found, and filled as the fill options say, by
    :func:`swmmbuilder.gaps.fill_gaps` (see also :meth:`Subnetwork.data_gaps`). Node inverts are
    those of :func:`node_inverts`, plus any filled ones. Each conduit is attached to its nodes with
    offsets (pipe invert - node invert, for LINK_OFFSETS DEPTH), so that it keeps its GIS
    inverts; conduits that slope uphill are reported with a warning. The topology is not checked
    here (see :func:`write_inp`).

    Parameters
    ----------
    subnetwork : Subnetwork
        The subnetwork.

    outfalls : dict or iterable of str, optional
        {node: outfall type, e.g. 'FREE'}, or the outfall nodes, which are then FREE. Every other
        node is a junction. Default: the subnetwork's outfall nodes
        (:attr:`Subnetwork.outfalls`), FREE.

    max_depth, init_depth, surcharge_depth, ponding_area : float
        Junction parameters (default 0.0). A junction whose rim elevation is known (see
        :attr:`Network.nodes`) and above its invert takes rim - invert as its maximum depth
        instead of ``max_depth``.

    manning_n : float
        Manning's n of every conduit (default 0.01) or, with ``manning_n_by_material``, of those
        whose material is not in the mapping.

    manning_n_by_material : dict, optional
        material -> Manning's n, e.g. ``{'PVC': 0.011, 'CONC': 0.013, 'VIT': 0.013}``.

    init_flow, max_flow : float
        Conduit initial and maximum flow (default 0.0).

    barrels : float
        Number of barrels of each conduit (default 1.0).

    height_correction : dict, optional
        {node: invert} overriding the computed node inverts.

    pumps : dict, optional
        {pipe_id: :class:`Pump` or dict of its fields}, e.g.
        ``{1010: dict(curve=[(0, 6.0), (2, 4.0), (4, 0.0)], curve_type='PUMP3')}``. Each pipe is
        written as a pump from its upstream to its downstream node instead of a conduit, with its
        curve as ``<pipe_id>_curve``; it needs no cross-section and is not reported as adverse.
        The GIS data do not identify pumps, so they must be given here.

    wet_wells : dict, optional
        {node: :class:`WetWell` or dict of its fields}, e.g. ``{'D1': dict(max_depth=3.0,
        area=4.0)}``: storage nodes of constant plan area instead of junctions. The inlet of a
        PUMP1 pump must be a wet well.

    dwf : dict, optional
        {node: base flow}: the dry-weather flow entering each node, in the template's FLOW_UNITS
        (L/s with the shipped template), written to [DWF]. Nodes with no upstream pipes, e.g.
        pumping stations fed by systems outside the data, take point loads this way. See
        :func:`dwf_from_population` for per-capita flows.

    patterns : dict, optional
        {name: multipliers} time patterns, written to [PATTERNS]. The type is set by the number
        of multipliers (24 HOURLY, 7 DAILY, 12 MONTHLY), or given as ``(type, multipliers)``,
        e.g. ``('WEEKEND', [...24 values])``. The multipliers must be positive; an HOURLY
        pattern's mean should be 1, or it scales the base flows (a warning is given).

    dwf_patterns : str, list of str or dict, optional
        The patterns applied to the base flows: a name or list of names for every node of
        ``dwf``, or {node: name or list of names}. SWMM allows up to four per node, one of each
        type.

    **fills
        Fill rules for gaps in the data, all off by default: ``size_overrides``, ``shape_codes``,
        ``size_fill``, ``invert_fill``, ``headwater_slope`` and ``force_main_discharge``; see
        :func:`swmmbuilder.gaps.fill_gaps`.

    Raises
    ----------
    ValueError
        If an outfall is not in the subnetwork, a gap that stops the model from being written
        remains after the fills (a conduit with no cross-section, a node with no invert; every
        such gap is listed), or a pump, wet well, dry-weather flow, pattern or fill rule is
        invalid.
    """
    outfalls = outfall_types(outfalls, subnetwork)
    unknown = sorted(set(outfalls) - set(subnetwork.nodes))
    if unknown:
        raise ValueError(f"Outfalls not in the subnetwork: {', '.join(map(str, unknown))}")
    pumps = resolve_pumps(subnetwork, pumps)
    wet_wells = resolve_wet_wells(subnetwork, wet_wells, outfalls)
    dwf = resolve_dwf(subnetwork, dwf, outfalls)
    patterns = resolve_patterns(patterns)
    dwf_patterns = resolve_dwf_patterns(dwf_patterns, dwf, patterns)
    pipes = subnetwork.pipes()
    inlets = dict(zip(pipes["pipe_id"], pipes["from_node"], strict=True))
    dry = sorted(
        f"{name} (inlet {inlets[name]})"
        for name, pump in pumps.items()
        if pump.curve != "*" and pump.curve_type == "PUMP1" and inlets[name] not in wet_wells
    )
    if dry:
        raise ValueError(
            f"PUMP1 pumps need a wet well at their inlet node (pass wet_wells): {', '.join(dry)}"
        )
    filled = fill_gaps(subnetwork, pumps=pumps, height_correction=height_correction, **fills)
    remaining = filled.remaining
    if len(remaining):
        lines = [f"{row.element} {row.id}: {row.detail}" for row in remaining.itertuples()]
        raise ValueError(
            f"Cannot write the subnetwork: {len(lines)} gaps in the data remain (see "
            "Subnetwork.data_gaps() for the full report, and its fill options):\n  "
            + "\n  ".join(lines)
            + f"\nKnown shape codes in schema {subnetwork.schema.name!r}: "
            f"{list(subnetwork.schema.shapes)}."
        )
    pipes, inverts = filled.pipes, filled.inverts
    sections = {
        name: xsection(row)
        for name, row in zip(pipes["pipe_id"], pipes.to_dict("records"), strict=True)
    }
    model = SwmmModel(data_gaps=filled.report)

    for node in subnetwork.nodes:
        invert = inverts[node]
        if height_correction is not None and node in height_correction:
            logger.info("Node %s: invert set to %s by height_correction", node, invert)
        if node in outfalls:
            model.outfalls.append(Outfall(node, invert, outfalls[node]))
        elif node in wet_wells:
            well = wet_wells[node]
            model.storage.append(Storage(node, invert, well.max_depth, well.init_depth, well.area))
        else:
            depth = max_depth
            rim = subnetwork.graph.nodes[node].get("rim_elevation")
            if rim is not None and rim > invert:
                depth = round(rim - invert, 4)
            elif rim is not None:
                logger.warning(
                    "Node %s: rim elevation %s is not above the invert %s; max depth %s is used",
                    node,
                    rim,
                    invert,
                    max_depth,
                )
            model.junctions.append(
                Junction(node, invert, depth, init_depth, surcharge_depth, ponding_area)
            )

    positions = subnetwork.node_positions()
    for node in subnetwork.nodes:
        if node in positions:
            model.coordinates[node] = positions[node]
        else:
            logger.warning("Node %s has no position; it is left out of [COORDINATES]", node)

    adverse = []
    for _, row in pipes.iterrows():
        name, inlet, outlet = row["pipe_id"], row["from_node"], row["to_node"]
        interior = line_coords(row["geometry"])[1:-1]
        if interior:
            model.vertices[name] = interior

        if name in pumps:
            pump = pumps[name]
            curve = "*" if pump.curve == "*" else f"{name}_curve"
            model.pumps.append(
                PumpLink(name, inlet, outlet, curve, pump.status, pump.startup, pump.shutoff)
            )
            if curve != "*":
                model.curves.append(Curve(curve, pump.curve_type, tuple(pump.curve)))
            continue

        roughness = (manning_n_by_material or {}).get(row["material"], manning_n)

        ends = []
        for node, pipe_invert in ((inlet, row["us_invert"]), (outlet, row["ds_invert"])):
            pipe_invert = _invert(pipe_invert)
            offset = 0.0 if pipe_invert is None else round(pipe_invert - inverts[node], 4)
            if offset < 0:
                warnings.warn(
                    f"Pipe {name}: invert {pipe_invert} is below the invert of node {node} "
                    f"({inverts[node]}); its offset there is set to 0.",
                    stacklevel=2,
                )
                offset = 0.0
            ends.append((offset, inverts[node] + offset))
        (inlet_offset, inlet_elev), (outlet_offset, outlet_elev) = ends
        if outlet_elev > inlet_elev:
            adverse.append(f"{name} ({inlet} {inlet_elev:.3f} -> {outlet} {outlet_elev:.3f})")

        model.conduits.append(
            Conduit(
                name,
                inlet,
                outlet,
                row["length"],
                roughness,
                inlet_offset,
                outlet_offset,
                init_flow,
                max_flow,
            )
        )
        model.xsections.append(XSection(name, *sections[name], barrels))

    if adverse:
        warnings.warn(
            f"Conduits with an adverse slope (outlet invert above inlet invert): "
            f"{', '.join(adverse)}.",
            stacklevel=2,
        )

    model.dwf = [
        DryWeatherFlow(node, flow, tuple(dwf_patterns.get(node, ()))) for node, flow in dwf.items()
    ]
    model.patterns = [Pattern(name, kind, values) for name, (kind, values) in patterns.items()]
    return model


def write_inp(subnetwork, path, outfalls=None, template_path=None, dead_ends="raise", **options):
    """Write a subnetwork to a SWMM .inp file.

    The template is checked with :func:`check_template` and the subnetwork with
    :meth:`Subnetwork.validate`: flow splits (including parallel pipes) and cycles are refused unless the
    template uses ``FLOW_ROUTING DYNWAVE``. Nothing is written if the template, the topology or
    the model (see :func:`build_model`) is invalid.

    Parameters
    ----------
    subnetwork : Subnetwork
        The subnetwork.

    path : str, Path
        The .inp file to write.

    outfalls : dict or iterable of str, optional
        The outfalls; see :func:`build_model`.

    template_path : str, Path, optional
        The .inp template (default: :func:`default_template_path`, which uses FLOW_UNITS LPS).
        It must use an SI FLOW_UNITS (CMS, LPS or MLD) and only the placeholders in
        :data:`TEMPLATE_PLACEHOLDERS`; {storage}, {pumps} and {curves} are needed only with the
        ``wet_wells`` and ``pumps`` options, and {dwf} and {patterns} only with ``dwf`` and
        ``patterns``. The offsets are depths, so it should set
        LINK_OFFSETS DEPTH (the default).

    dead_ends : {'raise', 'outfall'}
        What to do with nodes that have no outlet pipe and are not outfalls: raise a
        ``ValueError`` (default) or write them as FREE outfalls.

    **options
        Passed to :func:`build_model`.

    Returns
    ----------
    SwmmModel
        The model that was written.
    """
    if dead_ends not in ("raise", "outfall"):
        raise ValueError("dead_ends must be 'raise' or 'outfall'.")
    template = read_template(template_path)
    check_template(template)

    outfalls = outfall_types(outfalls, subnetwork)
    report = subnetwork.validate(outfalls)
    check_topology(report, template_routing(template), dead_ends)
    if dead_ends == "outfall":
        outfalls = {**outfalls, **{node: "FREE" for node in report.dead_ends}}

    model = build_model(subnetwork, outfalls, **options)
    needed = [
        name
        for name, records in (
            ("storage", model.storage),
            ("pumps", model.pumps),
            ("curves", model.curves),
            ("dwf", model.dwf),
            ("patterns", model.patterns),
        )
        if records and f"{{{name}}}" not in template
    ]
    if needed:
        raise ValueError(
            f"The template has no {', '.join(f'{{{n}}}' for n in needed)} placeholder, so the "
            "pumps, wet wells, dry-weather flows or patterns would not be written."
        )
    if model.dwf:
        logger.info(
            "Dry-weather flow: total base flow %.6g %s at %d nodes",
            sum(d.base for d in model.dwf),
            template_flow_units(template),
            len(model.dwf),
        )
    logger.info("Writing %s", path)
    Path(path).write_text(model.render(template))
    return model
