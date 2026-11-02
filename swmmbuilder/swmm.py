"""SWMM models of a subnetwork, and their .inp files.

:func:`build_model` turns a subnetwork into a :class:`SwmmModel`: plain records of the junctions,
outfalls, storage nodes, conduits, pumps, cross-sections, curves, coordinates and vertices, with all the hydraulic decisions
(node inverts, conduit offsets, cross-sections, Manning's n) already made.
:meth:`SwmmModel.render` fills an .inp template with them, and :func:`write_inp` checks the
template and the topology, builds the model and writes the file.
"""

import logging
import re
import warnings
from collections import Counter
from collections.abc import Mapping
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
    "coordinates",
    "vertices",
)

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
    match = re.search(r"^\s*FLOW_UNITS\s+(\S+)", content, re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise ValueError(
            "The template has no FLOW_UNITS option (SWMM would default to CFS, US units). "
            f"Set FLOW_UNITS to one of {', '.join(SI_FLOW_UNITS)}."
        )
    flow_units = match.group(1).upper()
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


@dataclass
class PumpStation:
    """A pumping station that discharges into one or more force mains (see the
    ``pump_stations`` option of :func:`build_model`).

    Attributes
    ----------
    pump : Pump or dict
        The pump. There is no default: an ideal pump (``curve='*'``) suits a variable-speed
        station that delivers a near-uniform flow, a curve with startup and shutoff depths
        (e.g. PUMP2 or PUMP4) a fixed-speed station that cycles on and off.

    wet_well : WetWell or dict, optional
        The station node's wet well; without one the station node is a junction.

    force_mains : pipe ID or list of pipe IDs, optional
        The first pipes of the force mains to pump into (default: every force main leaving the
        station node).

    discharge : {'shared', 'separate'}, optional
        With more than one force main: one pump whose discharge node feeds them all
        (``'shared'``), or one pump per force main (``'separate'``). Required then.
    """

    pump: object
    wet_well: object = None
    force_mains: object = None
    discharge: str | None = None

    @classmethod
    def coerce(cls, node, value):
        """A PumpStation from a PumpStation or a dict of its fields, with its pump and wet well
        coerced; raise ValueError if it is invalid."""
        station = replace(value) if isinstance(value, cls) else cls(**value)
        station.pump = Pump.coerce(node, station.pump)
        if station.wet_well is not None:
            station.wet_well = WetWell.coerce(node, station.wet_well)
        if station.force_mains is not None and (
            isinstance(station.force_mains, str) or not hasattr(station.force_mains, "__iter__")
        ):
            station.force_mains = [station.force_mains]
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


def resolve_pump_stations(subnetwork, pump_stations, outfalls=()):
    """({node: PumpStation}, [StationPump]) from the ``pump_stations`` option of
    :func:`build_model`.

    Raises
    ----------
    ValueError
        If a station node is not in the subnetwork, is an outfall or has no force-main outlet, a
        chosen force main does not leave it, it has several force mains and no ``discharge``,
        a pump or wet well is invalid, a PUMP1 pump has no wet well, or a name the station needs
        (``<node>_P``, ``<node>_D``) is already a node or pipe ID. Every problem is listed.
    """
    if not pump_stations:
        return {}, []
    nodes = set(subnetwork.nodes)
    taken = {str(n) for n in subnetwork.network.graph.nodes} | {
        str(p) for p in subnetwork.network.full_data["pipe_id"]
    }
    stations, links, problems = {}, [], []
    for node, value in pump_stations.items():
        if node not in nodes:
            problems.append(f"{node} is not in the subnetwork")
            continue
        if node in set(outfalls):
            problems.append(f"{node} is an outfall")
            continue
        station = stations[node] = PumpStation.coerce(node, value)
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
        pump = station.pump
        if pump.curve != "*" and pump.curve_type == "PUMP1" and station.wet_well is None:
            problems.append(f"{node} has a PUMP1 pump and no wet well")
        if station.discharge == "separate":
            groups = [(f"{node}_P{i}", f"{node}_D{i}", [p]) for i, p in enumerate(chosen, 1)]
        else:
            groups = [(f"{node}_P", f"{node}_D", chosen)]
        clashes = sorted({n for name, d, _ in groups for n in (name, d) if n in taken})
        if clashes:
            problems.append(f"{node}: {', '.join(clashes)} already exist in the network")
        for name, discharge, mains in groups:
            chains = {p: subnetwork.force_main_chain(p, outfalls)[0] for p in mains}
            links.append(StationPump(name, node, discharge, chains))
    if problems:
        raise ValueError("Invalid pump stations:\n  " + "\n  ".join(problems))
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


@dataclass
class SwmmModel:
    """The SWMM objects of a subnetwork, ready to render into an .inp template.

    Attributes
    ----------
    junctions, outfalls, storage, conduits, pumps, xsections, curves : list
        :class:`Junction`, :class:`Outfall`, :class:`Storage`, :class:`Conduit`,
        :class:`PumpLink`, :class:`XSection` and :class:`Curve` records, in the order they are
        written.
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
    force_main_roughness=None,
    force_main_roughness_by_pipe=None,
    force_main_surcharge_depth=100.0,
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

    **fills
        Fill rules for gaps in the data, all off by default: ``size_overrides``, ``shape_codes``,
        ``size_fill``, ``invert_fill``, ``headwater_slope`` and ``force_main_discharge``; see
        :func:`swmmbuilder.gaps.fill_gaps`.

    Raises
    ----------
    ValueError
        If an outfall is not in the subnetwork, a gap that stops the model from being written
        remains after the fills (a conduit with no cross-section, a node with no invert; every
        such gap is listed), a pump, wet well, pump station or fill rule is invalid, a pipe is
        both pumped and part of a station's force main, a force main has no roughness or is not
        circular, or an ideal pump is not the only outlet link of its inlet node (SWMM's
        ERROR 134).
    """
    outfalls = outfall_types(outfalls, subnetwork)
    unknown = sorted(set(outfalls) - set(subnetwork.nodes))
    if unknown:
        raise ValueError(f"Outfalls not in the subnetwork: {', '.join(map(str, unknown))}")
    pumps = resolve_pumps(subnetwork, pumps)
    wet_wells = resolve_wet_wells(subnetwork, wet_wells, outfalls)
    stations, station_pumps = resolve_pump_stations(subnetwork, pump_stations, outfalls)
    both = sorted(str(n) for n, s in stations.items() if s.wet_well is not None and n in wet_wells)
    if both:
        raise ValueError(f"Wet wells given both in wet_wells and pump_stations: {', '.join(both)}")
    wet_wells.update({n: s.wet_well for n, s in stations.items() if s.wet_well is not None})
    chains = {first: chain for link in station_pumps for first, chain in link.force_mains.items()}
    force_mains = {pipe for chain in chains.values() for pipe in chain}
    both = sorted(str(p) for p in pumps if p in force_mains)
    if both:
        raise ValueError(
            "Pipes both in pumps and in a pump station's force main: " + ", ".join(both)
        )
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
    roughness_of = force_main_roughnesses(
        pipes[pipes["pipe_id"].isin(list(force_mains))],
        sections,
        force_main_roughness,
        force_main_roughness_by_pipe,
    )
    outlets = dict(zip(pipes["pipe_id"], pipes["to_node"], strict=True))
    pressurised = {outlets[p] for chain in chains.values() for p in chain[:-1]}
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
                pump.startup,
                pump.shutoff,
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

    if adverse:
        warnings.warn(
            f"Conduits with an adverse slope (outlet invert above inlet invert): "
            f"{', '.join(adverse)}.",
            stacklevel=2,
        )
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
        ``wet_wells`` and ``pumps`` options. The offsets are depths, so it should set
        LINK_OFFSETS DEPTH (the default). A ``force_main_roughness`` rule that declares an
        ``equation`` must match its FORCE_MAIN_EQUATION (SWMM's default: H-W).

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

    model = build_model(subnetwork, outfalls, **options)
    needed = [
        name
        for name, records in (
            ("storage", model.storage),
            ("pumps", model.pumps),
            ("curves", model.curves),
        )
        if records and f"{{{name}}}" not in template
    ]
    if needed:
        raise ValueError(
            f"The template has no {', '.join(f'{{{n}}}' for n in needed)} placeholder, so the "
            "pumps, pump stations or wet wells would not be written."
        )
    logger.info("Writing %s", path)
    Path(path).write_text(model.render(template))
    return model
