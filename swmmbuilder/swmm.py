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
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from importlib import resources
from numbers import Real
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

# SWMM's FORCE_MAIN_EQUATION options, and what the roughness of a FORCE_MAIN cross-section is
# with each (with SI units).
FORCE_MAIN_EQUATIONS = {
    "H-W": "Hazen-Williams C",
    "D-W": "Darcy-Weisbach roughness height (mm)",
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


def template_force_main_equation(template):
    """FORCE_MAIN_EQUATION of the text of an .inp template; SWMM's default, H-W, if it is not
    set."""
    match = re.search(
        r"^\s*FORCE_MAIN_EQUATION\s+(\S+)", template, flags=re.MULTILINE | re.IGNORECASE
    )
    return match.group(1).upper() if match else "H-W"


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
        (variable-speed PUMP3). Ignored for an ideal pump. Required for a pump station's curve;
        a pumped pipe's (``pumps``) defaults to PUMP3.

    status : str
        Initial status, 'ON' or 'OFF'.

    startup, shutoff : float, optional
        Inlet node depths (m) at which the pump turns on and off; 0 for none. Not given (None),
        they are written as 0, unless a pump station's ``wet_well_defaults`` rule derives them.
    """

    curve: object
    curve_type: str | None = None
    status: str = "ON"
    startup: float | None = None
    shutoff: float | None = None

    @classmethod
    def coerce(cls, pump_id, value, curve_type=None):
        """A Pump from a Pump or a dict of its fields, checked, with ``curve_type`` if it has none;
        raise ValueError if it is invalid."""
        pump = replace(value) if isinstance(value, cls) else cls(**value)
        if pump.curve_type is None:
            pump.curve_type = curve_type
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
            if pump.curve_type is None:
                problems.append(f"a curve needs a curve_type, one of {list(PUMP_CURVE_TYPES)}")
            elif str(pump.curve_type).upper() not in PUMP_CURVE_TYPES:
                problems.append(
                    f"curve_type {pump.curve_type!r} is not one of {list(PUMP_CURVE_TYPES)}"
                )
            else:
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
    return {ids[str(k)]: Pump.coerce(k, v, curve_type="PUMP3") for k, v in pumps.items()}


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


@dataclass
class PumpStation:
    """A pumping station that discharges into one or more force mains (see the
    ``pump_stations`` option of :func:`build_model`).

    Attributes
    ----------
    pump : Pump or dict, optional
        The pump: an ideal pump (``curve='*'``) suits a variable-speed station that delivers a
        near-uniform flow, a curve with startup and shutoff depths (e.g. PUMP2 or PUMP4) a
        fixed-speed station that cycles on and off. There is no built-in default: a station
        without one takes the ``default_pump`` of :func:`build_model`, which must then be given.

    wet_well : WetWell or dict, optional
        The station node's wet well; without one the station node is a junction.

    force_mains : pipe ID or list of pipe IDs, optional
        The first pipes of the force mains to pump into (default: every force main leaving the
        station node).

    discharge : {'shared', 'separate'}, optional
        With more than one force main: one pump whose discharge node feeds them all
        (``'shared'``), or one pump per force main (``'separate'``). Required then.

    design_flow : float, optional
        Design average flow into the station, in the template's FLOW_UNITS (L/s with the shipped
        template). Used only by the ``wet_well_defaults`` rule (the wet well's plan area) and the
        ``force_main_checks`` of an ideal pump or a head-flow curve (its pumped flow); the GIS
        data do not carry it.

    suction_diameter : float, optional
        Diameter of the pump suction (mm), for the ``low_level_suction_diameters`` parameter of
        the ``wet_well_defaults`` rule.
    """

    pump: object = None
    wet_well: object = None
    force_mains: object = None
    discharge: str | None = None
    design_flow: float | None = None
    suction_diameter: float | None = None

    @classmethod
    def coerce(cls, node, value, default_pump=None):
        """A PumpStation from a PumpStation or a dict of its fields, with its pump (or
        ``default_pump``) and wet well coerced; raise ValueError if it is invalid."""
        station = replace(value) if isinstance(value, cls) else cls(**value)
        if station.pump is None:
            station.pump = default_pump
        if station.pump is None:
            raise ValueError(
                f"Pump station {node} has no pump: pass its pump, or a default_pump for the "
                "stations."
            )
        station.pump = Pump.coerce(node, station.pump)
        if station.wet_well is not None:
            station.wet_well = WetWell.coerce(node, station.wet_well)
        if station.force_mains is not None and (
            isinstance(station.force_mains, str) or not hasattr(station.force_mains, "__iter__")
        ):
            station.force_mains = [station.force_mains]
        for name in ("design_flow", "suction_diameter"):
            value = getattr(station, name)
            if value is not None:
                if not _is_finite_number(value) or not float(value) > 0:
                    raise ValueError(f"Pump station {node}: {name} must be positive.")
                setattr(station, name, float(value))
        if station.discharge not in (None, "shared", "separate"):
            raise ValueError(
                f"Pump station {node}: discharge must be 'shared' or 'separate', "
                f"not {station.discharge!r}."
            )
        return station


@dataclass(frozen=True)
class StationPump:
    """A pump of a :class:`PumpStation`, from the station node to its discharge node, and the
    force mains it feeds: {first pipe ID: pipe IDs of the chain}."""

    name: str
    station: str
    discharge: str
    force_mains: dict


def resolve_pump_stations(subnetwork, pump_stations, outfalls=(), default_pump=None):
    """({node: PumpStation}, [StationPump]) from the ``pump_stations`` option of
    :func:`build_model`.

    Raises
    ----------
    ValueError
        If a station node is not in the subnetwork, is an outfall or has no force-main outlet, a
        chosen force main does not leave it, it has several force mains and no ``discharge``,
        a pump is missing or invalid, a wet well is invalid, or a name the station needs
        (``<node>_P``, ``<node>_D``) is already a node or pipe ID. Every problem is listed.
    """
    if not pump_stations:
        return {}, []
    nodes = set(subnetwork.nodes)
    taken = {str(n) for n in subnetwork.network.graph.nodes} | {
        str(p) for p in subnetwork.network.full_data["pipe_id"]
    }
    typed = {n for n, t in subnetwork.graph.nodes(data="node_type") if t == "pump_station"}
    ends = typed | set(pump_stations)
    stations, links, problems = {}, [], []
    for node, value in pump_stations.items():
        if node not in nodes:
            problems.append(f"{node} is not in the subnetwork")
            continue
        if node in set(outfalls):
            problems.append(f"{node} is an outfall")
            continue
        station = stations[node] = PumpStation.coerce(node, value, default_pump)
        outlets = subnetwork.force_main_outlets(node)
        if not outlets:
            problems.append(f"{node} has no force-main outlet pipe")
            continue
        chosen = outlets
        if station.force_mains is not None:
            ids = {str(p): p for p in outlets}
            wrong = [str(p) for p in station.force_mains if str(p) not in ids]
            if wrong:
                problems.append(
                    f"{node}: force mains {', '.join(wrong)} do not leave it "
                    f"(its force mains: {', '.join(map(str, outlets))})"
                )
                continue
            chosen = [ids[str(p)] for p in station.force_mains]
            unused = [str(p) for p in outlets if p not in chosen]
            if unused:
                warnings.warn(
                    f"Pump station {node}: force mains {', '.join(unused)} are not pumped and "
                    "are written as ordinary conduits from the station node.",
                    stacklevel=3,
                )
        if len(chosen) > 1 and station.discharge is None:
            problems.append(
                f"{node} has {len(chosen)} force mains ({', '.join(map(str, chosen))}); pass "
                "discharge='shared' (one pump feeding all of them) or 'separate' (one pump each)"
            )
            continue
        if station.discharge == "separate":
            groups = [(f"{node}_P{i}", f"{node}_D{i}", [p]) for i, p in enumerate(chosen, 1)]
        else:
            groups = [(f"{node}_P", f"{node}_D", chosen)]
        clashes = sorted({n for name, d, _ in groups for n in (name, d) if n in taken})
        if clashes:
            problems.append(f"{node}: {', '.join(clashes)} already exist in the network")
        for name, discharge, mains in groups:
            chains = {p: subnetwork.force_main_chain(p, outfalls, True, ends)[0] for p in mains}
            links.append(StationPump(name, node, discharge, chains))
    if problems:
        raise ValueError("Invalid pump stations:\n  " + "\n  ".join(problems))
    undeclared = sorted(str(n) for n in typed - set(stations) if subnetwork.force_main_outlets(n))
    if undeclared:
        warnings.warn(
            f"Pump-station nodes with force mains not in pump_stations: {', '.join(undeclared)}; "
            "their force mains are written as ordinary conduits.",
            stacklevel=3,
        )
    return stations, links


def force_main_roughness_rule(rule):
    """A function pipe -> roughness (None if the rule has no value for the pipe) from the
    ``force_main_roughness`` option of :func:`build_model`.

    Raises
    ----------
    ValueError
        If the rule is not one of the accepted forms.
    """
    if callable(rule):
        return rule
    if isinstance(rule, Mapping):
        by, table = rule.get("by"), rule.get("table")
        unknown = sorted(set(rule) - {"by", "table", "equation"})
        if unknown or by not in ("diameter", "material") or table is None:
            raise ValueError(
                "force_main_roughness must be a number, a callable, or a dict with 'by' "
                "('diameter' or 'material'), 'table' and optionally 'equation'"
                + (f"; unknown keys {unknown}" if unknown else "")
                + "."
            )
        if by == "material":
            table = dict(table)
            return lambda pipe: table.get(pipe.get("material"))
        table = sorted((float(upper), float(value)) for upper, value in table)

        def by_diameter(pipe):
            diameter = pipe.get("diameter")
            if pd.isna(diameter):
                return None
            mm = round(1000 * diameter, 6)
            return next((value for upper, value in table if mm <= upper), None)

        return by_diameter
    value = float(rule)
    return lambda pipe: value


def force_main_roughnesses(pipes, sections, rule, by_pipe=None):
    """{pipe_id: roughness} of the force-main pipes ``pipes`` (see the ``force_main_roughness``
    and ``force_main_roughness_by_pipe`` options of :func:`build_model`).

    Raises
    ----------
    ValueError
        If a force main is not circular or has no roughness; every such pipe is listed.
    """
    if pipes.empty:
        return {}
    overrides = {str(k): float(v) for k, v in (by_pipe or {}).items()}
    rule = (lambda pipe: None) if rule is None else force_main_roughness_rule(rule)
    found, not_round, missing = {}, [], []
    for _, pipe in pipes.iterrows():
        name = pipe["pipe_id"]
        if sections[name][0] != "CIRCULAR":
            not_round.append(f"{name} ({sections[name][0]})")
        value = overrides.get(str(name))
        if value is None:
            value = rule(pipe)
        if value is None or pd.isna(value):
            missing.append(
                f"{name} (diameter {pipe.get('diameter')} m, material {pipe.get('material')!r})"
            )
        else:
            found[name] = float(value)
    problems = []
    if not_round:
        problems.append(f"not circular: {', '.join(not_round)}")
    if missing:
        problems.append(
            f"no roughness: {', '.join(missing)} (pass force_main_roughness, e.g. a "
            "Hazen-Williams C, or force_main_roughness_by_pipe)"
        )
    if problems:
        raise ValueError("Invalid force mains:\n  " + "\n  ".join(problems))
    return found


def check_force_main_equation(rule, equation):
    """Raise ValueError if a ``force_main_roughness`` rule declares an ``equation`` that is not
    the template's FORCE_MAIN_EQUATION."""
    if not isinstance(rule, Mapping) or "equation" not in rule:
        return
    declared = str(rule["equation"]).upper()
    if declared != equation:
        raise ValueError(
            f"force_main_roughness gives {FORCE_MAIN_EQUATIONS.get(declared, declared)} values "
            f"(equation {declared}), but the template sets FORCE_MAIN_EQUATION {equation}, so "
            f"SWMM would read them as {FORCE_MAIN_EQUATIONS.get(equation, equation)}."
        )


def _is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _is_real(value):
    """True for a finite real number that is not a bool or a string."""
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)


@dataclass(frozen=True)
class Derived:
    """A model value that comes from a design rule, not from the data (see
    :attr:`SwmmModel.derived`).

    Attributes
    ----------
    element : str
        'storage' or 'pump'.
    id : str
        The storage node, or the pump link (``<node>_P``, ``<node>_P1``, ...).
    attribute : str
        e.g. 'area', 'max_depth', 'startup', 'shutoff'.
    value : float
        The value written.
    rule : str
        The rule that gave it, e.g. 'wet_well_defaults'.
    basis : str
        How it was derived.
    """

    element: str
    id: str
    attribute: str
    value: float
    rule: str
    basis: str


WET_WELL_DEFAULTS = ("fill_time", "high_level_below_inlet", "low_level_above_pump")
WET_WELL_OPTIONS = ("low_level_suction_diameters",)
FORCE_MAIN_CHECKS = ("min_diameter", "velocity")


def _rule_params(name, rule, keys, required):
    """The parameters of a design rule given as a dict; raise ValueError if they are not
    ``keys``, or a required one is missing."""
    if not isinstance(rule, Mapping):
        raise ValueError(f"{name} must be a dict with keys {list(keys)}, not {rule!r}.")
    unknown = sorted(set(rule) - set(keys))
    missing = [k for k in required if k not in rule]
    if unknown or missing:
        raise ValueError(
            f"{name} takes the keys {list(keys)}"
            + (f"; unknown keys {unknown}" if unknown else "")
            + (f"; missing keys {missing}" if missing else "")
            + "."
        )
    return dict(rule)


NOT_GRAVITY = ("force_main", "overflow")


def gravity_inlets(pipes, node, force_mains=()):
    """Rows of ``pipes`` entering ``node`` that have a downstream invert and are neither in
    ``force_mains`` nor of a structure type in :data:`NOT_GRAVITY`."""
    entering = pipes[(pipes["to_node"] == node) & ~pipes["pipe_id"].isin(list(force_mains))]
    if "structure_type" in entering:
        entering = entering[~entering["structure_type"].isin(NOT_GRAVITY)]
    return entering[entering["ds_invert"].notna()]


def wet_well_params(rule):
    """(fill_time, high_level_below_inlet, low_level_above_pump, low_level_suction_diameters)
    of a ``wet_well_defaults`` rule, the last 0 if not given; raise ValueError if the rule is
    invalid."""
    params = _rule_params(
        "wet_well_defaults", rule, WET_WELL_DEFAULTS + WET_WELL_OPTIONS, WET_WELL_DEFAULTS
    )
    values = [params[k] for k in WET_WELL_DEFAULTS] + [params.get("low_level_suction_diameters", 0)]
    if not all(_is_real(v) for v in values):
        raise ValueError(f"wet_well_defaults: the parameters must be numbers, not {params}.")
    fill_time, below_inlet, above_pump, suctions = map(float, values)
    if not (fill_time > 0 and below_inlet >= 0 and above_pump >= 0 and suctions >= 0):
        raise ValueError(
            "wet_well_defaults: fill_time must be positive and the levels not negative."
        )
    return fill_time, below_inlet, above_pump, suctions


def apply_wet_well_defaults(
    subnetwork, pipes, stations, inverts, rule, flow_units, force_mains=(), pump_names=None
):
    """Fill in the stations' missing wet wells and pump on/off depths from the
    ``wet_well_defaults`` rule of :func:`build_model`, in place.

    The wet well's floor is the station node invert, taken as the pump level, and the inlet
    invert is the lowest downstream invert of the gravity pipes entering the station (see
    :func:`gravity_inlets`):

    - shutoff depth (low water level) = ``low_level_above_pump``, or, if the rule gives
      ``low_level_suction_diameters`` and the station a ``suction_diameter``, that many suction
      diameters if it is higher;
    - startup depth (high water level) = inlet invert - ``high_level_below_inlet`` - floor;
    - plan area = design flow x ``fill_time`` (minutes) / (startup - shutoff depth), so that the
      design flow fills the well from the low to the high water level in ``fill_time``;
    - max depth = rim - floor if the rim is known and above the inlet, else the crown of the
      lowest inlet sewer - floor.

    A station with a :class:`WetWell` keeps it, and a pump's own ``startup``/``shutoff`` are kept;
    the on/off depths are derived only for a pump with a curve (an ideal pump does not cycle).
    ``pump_names`` ({node: [pump link names]}, default the node) names the pumps in the
    :class:`Derived` rows.

    Returns
    ----------
    list of Derived
        Every value the rule gave.

    Raises
    ----------
    ValueError
        If the rule is invalid, or a station that needs it has no gravity inlet with an invert,
        a high water level not above its low water level, or (for its wet well) no
        ``design_flow``. Every problem is listed.
    """
    fill_time, below_inlet, above_pump, suctions = wet_well_params(rule)
    derived, problems = [], []
    for node, station in stations.items():
        pump = station.pump
        needs_well = station.wet_well is None
        needs_levels = pump.curve != "*" and (pump.startup is None or pump.shutoff is None)
        if not (needs_well or needs_levels):
            continue
        floor = inverts[node]
        inlets = gravity_inlets(pipes, node, force_mains)
        if inlets.empty:
            problems.append(f"{node} has no gravity inlet pipe with an invert")
            continue
        lowest = inlets.loc[inlets["ds_invert"].idxmin()]
        inlet = float(lowest["ds_invert"])
        low_how = f"{above_pump} m above the floor"
        low = round(above_pump, 4)
        if suctions and station.suction_diameter is not None:
            by_suction = round(suctions * station.suction_diameter / 1000, 4)
            if by_suction > low:
                low = by_suction
                low_how = (
                    f"{suctions:g} suction diameters ({station.suction_diameter:g} mm) above "
                    "the floor"
                )
        high = round(inlet - below_inlet - floor, 4)
        if high <= low:
            problems.append(
                f"{node}: the high water level ({below_inlet} m below the inlet invert "
                f"{inlet}) is not above the low water level ({low_how}, floor {floor})"
            )
            continue
        basis = f"floor {floor}, inlet invert {inlet} (pipe {lowest['pipe_id']})"
        found = []
        if needs_well:
            if station.design_flow is None:
                problems.append(f"{node} has no wet well and no design_flow to size one")
                continue
            volume = station.design_flow * _FLOW_UNIT_LPS[flow_units] / 1000 * fill_time * 60
            area = round(volume / (high - low), 4)
            rim = subnetwork.graph.nodes[node].get("rim_elevation")
            if rim is not None and rim > inlet:
                max_depth, how = round(rim - floor, 4), "rim - floor"
            else:
                size = lowest["diameter"] if lowest.get("shape") == "CIRCULAR" else None
                size = lowest["height"] if size is None or pd.isna(size) else size
                crown = inlet + (0.0 if pd.isna(size) else float(size))
                max_depth, how = round(crown - floor, 4), "inlet sewer crown - floor"
            station.wet_well = WetWell(max_depth=max_depth, area=area)
            found += [
                (
                    "storage",
                    "area",
                    area,
                    f"design flow {station.design_flow} {flow_units} for {fill_time:g} min "
                    f"between depths {low} and {high} m",
                ),
                ("storage", "max_depth", max_depth, how),
            ]
        if needs_levels and pump.startup is None:
            pump.startup = high
            found.append(("pump", "startup", high, f"{below_inlet} m below the inlet invert"))
        if needs_levels and pump.shutoff is None:
            pump.shutoff = low
            found.append(("pump", "shutoff", low, low_how))
        names = (pump_names or {}).get(node, [node])
        for element, attribute, value, how in found:
            for name in [node] if element == "storage" else names:
                derived.append(
                    Derived(element, name, attribute, value, "wet_well_defaults", f"{how}; {basis}")
                )
    if problems:
        raise ValueError("wet_well_defaults cannot be applied:\n  " + "\n  ".join(problems))
    for d in derived:
        logger.info(
            "%s %s: %s = %s from rule %s, not data (%s)",
            d.element.capitalize(),
            d.id,
            d.attribute,
            d.value,
            d.rule,
            d.basis,
        )
    return derived


HEAD_FLOW_CURVES = ("PUMP3", "PUMP5")


def pump_level_problem(name, pump, well=None):
    """Why a pump's startup and shutoff depths cannot work, or None: a startup depth not above
    the shutoff depth, or not below the max depth of its inlet's wet well ``well``."""
    startup, shutoff = pump.startup or 0.0, pump.shutoff or 0.0
    if startup > 0 and startup <= shutoff:
        return f"{name}: startup depth {startup} m is not above the shutoff depth {shutoff} m"
    if well is not None and startup >= well.max_depth:
        return (
            f"{name}: startup depth {startup} m is not below the wet well's max depth "
            f"{well.max_depth} m, so the pump never starts"
        )
    return None


def pumped_flow(station):
    """Flow (in the FLOW_UNITS) a station's pump delivers: the largest flow of its curve, or
    its ``design_flow`` (None if it has none) for an ideal pump or a head-flow curve (PUMP3,
    PUMP5), whose largest flow is at zero head, not against the force main."""
    pump = station.pump
    if pump.curve == "*" or pump.curve_type in HEAD_FLOW_CURVES:
        return station.design_flow
    return max(y for _, y in pump.curve)


def force_main_check_params(rule):
    """(min_diameter, velocity) of a ``force_main_checks`` rule, each None if not given;
    raise ValueError if the rule is invalid."""
    params = _rule_params("force_main_checks", rule, FORCE_MAIN_CHECKS, ())
    min_diameter = params.get("min_diameter")
    velocity = params.get("velocity")
    if min_diameter is not None:
        if not _is_real(min_diameter) or not min_diameter > 0:
            raise ValueError(
                f"force_main_checks: min_diameter must be a positive number (mm), "
                f"not {min_diameter!r}."
            )
        min_diameter = float(min_diameter)
    if velocity is not None:
        if (
            isinstance(velocity, (str, bytes))
            or not hasattr(velocity, "__len__")
            or len(velocity) != 2
            or not all(_is_real(v) and v >= 0 for v in velocity)
        ):
            raise ValueError(
                "force_main_checks: velocity must be (min, max), two numbers in m/s that are "
                f"not negative, not {velocity!r}."
            )
        velocity = tuple(sorted(float(v) for v in velocity))
    return min_diameter, velocity


def check_design_rule(name, rule):
    """Raise ValueError if ``rule`` is not a valid ``name`` design rule of :func:`build_model`
    (``force_main_discharge`` is checked by :func:`swmmbuilder.gaps.fill_gaps`)."""
    if name == "wet_well_defaults":
        wet_well_params(rule)
    elif name == "force_main_checks":
        force_main_check_params(rule)
    elif name == "force_main_roughness":
        force_main_roughness_rule(rule)
    elif name == "default_pump":
        Pump.coerce("default_pump", rule)


def check_force_mains(station_pumps, stations, sections, rule, flow_units):
    """Warnings (strings) for the force mains of :class:`StationPump` links that break the
    ``force_main_checks`` of :func:`build_model`: a diameter below ``min_diameter`` (mm) or a
    velocity (pumped flow / full area) outside ``velocity`` (min, max), in m/s.

    The pumped flow is :func:`pumped_flow`, summed over the pumps feeding a main shared
    downstream of a join. With ``discharge='shared'`` and several force mains the split between
    them is not known, so the velocity of those mains is not checked. Raise ValueError if the
    rule is invalid.
    """
    min_diameter, velocity = force_main_check_params(rule)
    flows = {}  # pipe -> pumped flows of the links feeding it (None: unknown)
    for link in station_pumps:
        flow = pumped_flow(stations[link.station])
        if len(link.force_mains) > 1:
            flow = None
        for chain in link.force_mains.values():
            for name in chain:
                flows.setdefault(name, []).append(flow)
    narrow, slow, fast, unchecked = [], [], [], []
    for name in sorted(flows, key=lambda p: (isinstance(p, str), p)):
        feeds = flows[name]
        diameter = sections[name][1]
        if min_diameter is not None and round(1000 * diameter, 6) < min_diameter:
            narrow.append(f"{name} ({1000 * diameter:g} mm)")
        if velocity is None:
            continue
        if None in feeds:
            unchecked.append(str(name))
            continue
        flow = sum(feeds)
        v = flow * _FLOW_UNIT_LPS[flow_units] / 1000 / (math.pi / 4 * diameter**2)
        if v < velocity[0]:
            slow.append(f"{name} ({v:.2f} m/s)")
        elif v > velocity[1]:
            fast.append(f"{name} ({v:.2f} m/s)")
    found = []
    if narrow:
        found.append(f"Force mains narrower than {min_diameter:g} mm: {', '.join(narrow)}.")
    if slow:
        found.append(f"Force mains with a velocity below {velocity[0]:g} m/s: {', '.join(slow)}.")
    if fast:
        found.append(f"Force mains with a velocity above {velocity[1]:g} m/s: {', '.join(fast)}.")
    if unchecked:
        logger.info(
            "Force-main velocity not checked (no design_flow for an ideal pump or a head-flow "
            "curve, or a shared pump feeding several force mains): %s",
            ", ".join(unchecked),
        )
    return found


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
    coordinates: dict = field(default_factory=dict)
    vertices: dict = field(default_factory=dict)
    data_gaps: pd.DataFrame | None = None
    derived: list = field(default_factory=list)

    def derived_report(self):
        """DataFrame of :attr:`derived`, one row per value (element, id, attribute, value, rule,
        basis)."""
        columns = ["element", "id", "attribute", "value", "rule", "basis"]
        return pd.DataFrame([vars(d) for d in self.derived], columns=columns)

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
    pump_stations=None,
    default_pump=None,
    force_main_roughness=None,
    force_main_roughness_by_pipe=None,
    force_main_surcharge_depth=100.0,
    dwf=None,
    patterns=None,
    dwf_patterns=None,
    wet_well_defaults=None,
    force_main_checks=None,
    flow_units="LPS",
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

    pump_stations : dict, optional
        {station node: :class:`PumpStation` or dict of its fields}, e.g. ``{'PS1':
        dict(pump=dict(curve='*'), wet_well=dict(max_depth=3.0, area=4.0))}``. For each station,
        a discharge node ``<node>_D`` is written at the station's coordinates, with the invert of
        the force main's upstream end, and a pump ``<node>_P`` from the station node to it (curve
        ``<node>_curve``). The force mains leaving the station (see
        :meth:`Subnetwork.force_main_chain`) start at ``<node>_D`` and are kept as conduits with
        their GIS lengths and a FORCE_MAIN cross-section; their interior nodes and ``<node>_D``
        get ``force_main_surcharge_depth``, so that the main can run pressurised under DYNWAVE.
        Other outlets of the station, e.g. overflows, stay conduits. With ``discharge='separate'``
        the pumps and discharge nodes are numbered: ``<node>_P1``, ``<node>_D1``, ...

    default_pump : Pump or dict, optional
        The pump of the stations in ``pump_stations`` given without one. Not applied unless
        given; London's preset, an ideal pump for its variable-speed stations (s. 4.3.5), is
        ``LONDON_ON.design['default_pump']``.

    force_main_roughness : float, dict or callable
        Roughness of the FORCE_MAIN cross-sections, read by SWMM as a Hazen-Williams C or a
        Darcy-Weisbach roughness height (mm) according to the template's FORCE_MAIN_EQUATION.
        Required with ``pump_stations``; no rule is applied unless it is given here. One of:

        - a constant;
        - ``{'by': 'diameter', 'table': [(upper diameter in mm, value), ...]}``: the value of the
          first row whose upper diameter is at or above the pipe's;
        - ``{'by': 'material', 'table': {material: value}}``;
        - a callable ``pipe -> value``, given the pipe's row of the canonical data.

        A dict may also give ``'equation': 'H-W'`` or ``'D-W'``; :func:`write_inp` then refuses a
        template with another FORCE_MAIN_EQUATION. Presets are in the schemas, e.g.
        ``LONDON_ON.design['force_main_roughness']`` (C by diameter).

    force_main_roughness_by_pipe : dict, optional
        {pipe_id: roughness}, overriding ``force_main_roughness``.

    force_main_surcharge_depth : float
        Surcharge depth (m) of the force mains' interior and discharge nodes (default 100).

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

    wet_well_defaults : dict, optional
        A design rule for stations whose wet well is not known: ``{'fill_time': minutes,
        'high_level_below_inlet': m, 'low_level_above_pump': m}``, and optionally
        ``'low_level_suction_diameters'`` (the low water level is then at least that many of a
        station's ``suction_diameter`` above the floor). A station with no wet well
        (in ``pump_stations`` or ``wet_wells``) gets one whose plan area holds its
        ``design_flow`` for ``fill_time`` between the low and high water levels, and a pump with
        a curve and no ``startup``/``shutoff`` gets them at those levels (see
        :func:`apply_wet_well_defaults`). Explicit values always win, and every derived value is
        logged and listed in :attr:`SwmmModel.derived`. Not applied unless given; London's
        preset is ``LONDON_ON.design['wet_well_defaults']``.

    force_main_checks : dict, optional
        Validation of the stations' force mains, as warnings: ``{'min_diameter': mm,
        'velocity': (min, max) m/s}`` (either key may be left out). The velocity is the pumped
        flow (the largest flow of the pump curve, or the station's ``design_flow`` for an ideal
        pump) over the full pipe area. London's preset is
        ``LONDON_ON.design['force_main_checks']``.

    flow_units : str
        The template's FLOW_UNITS, in which pump curves and ``design_flow`` are given (default
        LPS, the shipped template's; :func:`write_inp` sets it from the template).

    Raises
    ----------
    ValueError
        If an outfall is not in the subnetwork, a gap that stops the model from being written
        remains after the fills (a conduit with no cross-section, a node with no invert; every
        such gap is listed), a pump, wet well, pump station, dry-weather flow, pattern or fill
        rule is invalid, a pipe is both pumped and part of a station's force main, a force main
        has no roughness or is not circular, a pump's startup depth is not above its shutoff depth
        or not below its wet well's max depth, a PUMP1 station has no wet well, or an ideal pump
        is not the only outlet link of its inlet node (SWMM's ERROR 134).
    """
    outfalls = outfall_types(outfalls, subnetwork)
    unknown = sorted(set(outfalls) - set(subnetwork.nodes))
    if unknown:
        raise ValueError(f"Outfalls not in the subnetwork: {', '.join(map(str, unknown))}")
    pumps = resolve_pumps(subnetwork, pumps)
    wet_wells = resolve_wet_wells(subnetwork, wet_wells, outfalls)
    stations, station_pumps = resolve_pump_stations(
        subnetwork, pump_stations, outfalls, default_pump
    )
    both = sorted(str(n) for n, s in stations.items() if s.wet_well is not None and n in wet_wells)
    if both:
        raise ValueError(f"Wet wells given both in wet_wells and pump_stations: {', '.join(both)}")
    flow_units = str(flow_units).upper()
    if flow_units not in _FLOW_UNIT_LPS:
        raise ValueError(f"flow_units must be one of {list(_FLOW_UNIT_LPS)}, not {flow_units!r}.")
    wet_wells.update({n: s.wet_well for n, s in stations.items() if s.wet_well is not None})
    chains = {first: chain for link in station_pumps for first, chain in link.force_mains.items()}
    force_mains = {pipe for chain in chains.values() for pipe in chain}
    both = sorted(str(p) for p in pumps if p in force_mains)
    if both:
        raise ValueError(
            "Pipes both in pumps and in a pump station's force main: " + ", ".join(both)
        )
    curves = {f"{p}_curve" for p, pump in pumps.items() if pump.curve != "*"}
    both = sorted(
        f"{n}_curve" for n, s in stations.items() if s.pump.curve != "*" and f"{n}_curve" in curves
    )
    if both:
        raise ValueError(
            "Pump curve names of both a pumped pipe and a pump station (a pipe ID equal to a "
            f"station node): {', '.join(both)}"
        )
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
    derived = []
    if wet_well_defaults is not None:
        for node, station in stations.items():
            if station.wet_well is None and node in wet_wells:
                station.wet_well = wet_wells[node]
        derived = apply_wet_well_defaults(
            subnetwork,
            pipes,
            stations,
            inverts,
            wet_well_defaults,
            flow_units,
            force_mains,
            pump_names={
                node: [link.name for link in station_pumps if link.station == node]
                for node in stations
            },
        )
    wet_wells.update({n: s.wet_well for n, s in stations.items() if s.wet_well is not None})
    problems = [
        f"{node} has a PUMP1 pump and no wet well (pass wet_well, wet_wells or wet_well_defaults)"
        for node, station in stations.items()
        if station.pump.curve != "*"
        and station.pump.curve_type == "PUMP1"
        and node not in wet_wells
    ]
    problems += [
        problem
        for node, station in stations.items()
        if (problem := pump_level_problem(node, station.pump, wet_wells.get(node)))
    ]
    problems += [
        problem
        for name, pump in pumps.items()
        if (problem := pump_level_problem(name, pump, wet_wells.get(inlets[name])))
    ]
    if problems:
        raise ValueError("Invalid pumps:\n  " + "\n  ".join(problems))
    sections = {
        name: xsection(row)
        for name, row in zip(pipes["pipe_id"], pipes.to_dict("records"), strict=True)
    }
    roughness_of = force_main_roughnesses(
        pipes[pipes["pipe_id"].isin(list(force_mains))],
        sections,
        force_main_roughness,
        force_main_roughness_by_pipe,
    )
    outlets = dict(zip(pipes["pipe_id"], pipes["to_node"], strict=True))
    pressurised = {outlets[p] for chain in chains.values() for p in chain[:-1]}
    model = SwmmModel(data_gaps=filled.report, derived=derived)

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
            surcharge = surcharge_depth
            if node in pressurised:
                surcharge = max(surcharge_depth, force_main_surcharge_depth)
            model.junctions.append(
                Junction(node, invert, depth, init_depth, surcharge, ponding_area)
            )

    positions = subnetwork.node_positions()
    for node in subnetwork.nodes:
        if node in positions:
            model.coordinates[node] = positions[node]
        else:
            logger.warning("Node %s has no position; it is left out of [COORDINATES]", node)

    us_inverts = dict(zip(pipes["pipe_id"], pipes["us_invert"], strict=True))
    discharge_of = {}
    for link in station_pumps:
        ends = [_invert(us_inverts[p]) for p in link.force_mains]
        ends = [e for e in ends if e is not None]
        inverts[link.discharge] = min(ends) if ends else inverts[link.station]
        surcharge = max(surcharge_depth, force_main_surcharge_depth)
        model.junctions.append(
            Junction(
                link.discharge,
                inverts[link.discharge],
                max_depth,
                init_depth,
                surcharge,
                ponding_area,
            )
        )
        if link.station in positions:
            model.coordinates[link.discharge] = positions[link.station]
        discharge_of.update(dict.fromkeys(link.force_mains, link.discharge))

    adverse = []
    for _, row in pipes.iterrows():
        name, inlet, outlet = row["pipe_id"], row["from_node"], row["to_node"]
        inlet = discharge_of.get(name, inlet)
        interior = line_coords(row["geometry"])[1:-1]
        if interior:
            model.vertices[name] = interior

        if name in pumps:
            pump = pumps[name]
            curve = "*" if pump.curve == "*" else f"{name}_curve"
            model.pumps.append(
                PumpLink(
                    name,
                    inlet,
                    outlet,
                    curve,
                    pump.status,
                    pump.startup or 0.0,
                    pump.shutoff or 0.0,
                )
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
        if outlet_elev > inlet_elev and name not in force_mains:
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
        if name in force_mains:
            diameter = sections[name][1]
            model.xsections.append(
                XSection(name, "FORCE_MAIN", diameter, roughness_of[name], barrels)
            )
        else:
            model.xsections.append(XSection(name, *sections[name], barrels))

    for link in station_pumps:
        pump = stations[link.station].pump
        curve = "*" if pump.curve == "*" else f"{link.station}_curve"
        model.pumps.append(
            PumpLink(
                link.name,
                link.station,
                link.discharge,
                curve,
                pump.status,
                pump.startup or 0.0,
                pump.shutoff or 0.0,
            )
        )
        if curve != "*" and curve not in {c.name for c in model.curves}:
            model.curves.append(Curve(curve, pump.curve_type, tuple(pump.curve)))

    outlet_links = Counter(link.inlet for link in model.conduits + model.pumps)
    lone = sorted(
        f"{p.name} (inlet {p.inlet}, {outlet_links[p.inlet]} outlet links)"
        for p in model.pumps
        if p.curve == "*" and outlet_links[p.inlet] > 1
    )
    if lone:
        raise ValueError(
            "An ideal pump (curve '*') must be the only outlet link of its inlet node (SWMM "
            f"ERROR 134); use a pump curve or drop the other outlets: {', '.join(lone)}"
        )

    if force_main_checks is not None:
        for message in check_force_mains(
            station_pumps, stations, sections, force_main_checks, flow_units
        ):
            warnings.warn(message, stacklevel=2)

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
        ``wet_wells``, ``pumps`` and ``pump_stations`` options, and {dwf} and {patterns} only with
        ``dwf`` and ``patterns``. The offsets are depths, so it should set LINK_OFFSETS DEPTH
        (the default). A ``force_main_roughness`` rule that declares an ``equation`` must match
        its FORCE_MAIN_EQUATION (SWMM's default: H-W).

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
    check_force_main_equation(
        options.get("force_main_roughness"), template_force_main_equation(template)
    )

    outfalls = outfall_types(outfalls, subnetwork)
    report = subnetwork.validate(outfalls)
    check_topology(report, template_routing(template), dead_ends)
    if dead_ends == "outfall":
        outfalls = {**outfalls, **{node: "FREE" for node in report.dead_ends}}

    flow_units = template_flow_units(template)
    if str(options.get("flow_units", flow_units)).upper() != flow_units:
        raise ValueError(
            f"flow_units {options['flow_units']!r} is not the template's FLOW_UNITS {flow_units}."
        )
    model = build_model(subnetwork, outfalls, **{**options, "flow_units": flow_units})
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
            "pumps, pump stations, wet wells, dry-weather flows or patterns would not be written."
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
