"""Design rules from a design standard (e.g. a schema's ``design`` table): force-main
roughness, wet-well defaults and force-main checks."""

import logging
import math
from collections.abc import Mapping
from dataclasses import dataclass

import pandas as pd

from .diagnostics import Problems
from .pumps import Pump, WetWell
from .units import _FLOW_UNIT_LPS, _is_real

logger = logging.getLogger("swmmbuilder.swmm")


# SWMM's FORCE_MAIN_EQUATION options, and what the roughness of a FORCE_MAIN cross-section is
# with each (with SI units).
FORCE_MAIN_EQUATIONS = {
    "H-W": "Hazen-Williams C",
    "D-W": "Darcy-Weisbach roughness height (mm)",
}


def force_main_roughness_rule(rule):
    """A function pipe -> roughness (None if the rule has no value for the pipe) from the
    ``force_main_roughness`` option of :func:`build_model`.

    Raises
    ------
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
    ------
    ValueError
        If a force main is not circular or has no roughness; every such pipe is listed.
    """
    if pipes.empty:
        return {}
    # Keyed by the pipe IDs of the data (see Subnetwork.match_pipes); text keys still match.
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
        raise ValueError(Problems("Invalid force mains", problems).message())
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


@dataclass(frozen=True)
class Derived:
    """A model value that comes from a design rule, not from the data (see
    :attr:`SwmmModel.derived`).

    Attributes
    ----------
    element : str
        'storage', 'pump', or 'junction'.
    id : str
        The storage node, the pump link (``<node>_P``, ``<node>_P1``, ...), or the junction.
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
    -------
    list of Derived
        Every value the rule gave.

    Raises
    ------
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
        raise ValueError(Problems("wet_well_defaults cannot be applied", problems).message())
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
