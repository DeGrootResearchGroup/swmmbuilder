"""Building the :class:`~swmmbuilder.swmm.SwmmModel` of a subnetwork and writing its .inp file."""

from __future__ import annotations

import logging
import warnings
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from ..design import (
    apply_wet_well_defaults,
    check_force_main_equation,
    check_force_mains,
    force_main_roughnesses,
    pump_level_problem,
)
from ..diagnostics import DataWarning, DesignWarning, Problems
from ..gaps import _invert, fill_gaps, height_corrections, lowest_inverts, xsection
from ..geometry import line_coords
from ..units import _FLOW_UNIT_LPS
from .dwf import resolve_dwf, resolve_dwf_patterns, resolve_patterns
from .pumping import resolve_pump_stations, resolve_pumps, resolve_wet_wells
from .records import (
    Conduit,
    Curve,
    DryWeatherFlow,
    Junction,
    Outfall,
    Pattern,
    PumpLink,
    Storage,
    SwmmModel,
    XSection,
)
from .template import (
    check_template,
    check_topology,
    read_template,
    template_flow_units,
    template_force_main_equation,
    template_routing,
)

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence
    from typing import Any

    from .._types import Fields, NodeId, Outfalls, PipeId, StrPath
    from ..pumps import Pump, PumpStation, WetWell
    from ..subnetwork import Subnetwork

logger = logging.getLogger(__name__)


def node_inverts(
    subnetwork: Subnetwork, height_correction: Mapping[NodeId, float] | None = None
) -> dict[NodeId, float]:
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
    -------
    dict
        {node: invert}

    Raises
    ------
    ValueError
        If a node has no pipe with an invert.
    """
    nodes = set(subnetwork.nodes)
    inverts = lowest_inverts(subnetwork.pipes(), nodes)
    if nodes - set(inverts):
        inverts.update(subnetwork.network.lowest_inverts(nodes - set(inverts)))
    missing = sorted(nodes - set(inverts))
    if missing:
        raise ValueError(f"No pipe invert found for nodes: {', '.join(map(str, missing))}")
    inverts.update(height_corrections(subnetwork, height_correction))
    return inverts


def outfall_types(
    outfalls: Outfalls | None, subnetwork: Subnetwork | None = None
) -> dict[NodeId, str]:
    """{node: outfall type} from a mapping, or from an iterable of nodes (all 'FREE'); None is
    the subnetwork's outfall nodes (:attr:`Subnetwork.outfalls`)."""
    if outfalls is None:
        if subnetwork is None:
            raise ValueError("outfall_types needs the outfalls or the subnetwork.")
        outfalls = subnetwork.outfalls
    types: dict[NodeId, str]
    if isinstance(outfalls, str) or not isinstance(outfalls, Iterable):
        types = {outfalls: "FREE"}
    elif isinstance(outfalls, Mapping):
        types = dict(outfalls)
    else:
        types = dict.fromkeys(outfalls, "FREE")
    if subnetwork is None:
        return types
    return subnetwork.match_nodes(types, "outfalls")


# ---- The phases of build_model ---------------------------------------------------------------


@dataclass(frozen=True)
class _Defaults:
    """The junction and conduit defaults of :func:`build_model`."""

    max_depth: float
    init_depth: float
    surcharge_depth: float
    ponding_area: float
    force_main_surcharge_depth: float
    manning_n: float
    manning_n_by_material: Mapping[Any, float]
    init_flow: float
    max_flow: float
    barrels: float


@dataclass(frozen=True)
class ModelInputs:
    """The options of :func:`build_model`, resolved against a subnetwork and cross-checked.

    The mappings may be filled in by later phases (``wet_wells`` and the stations' wet wells
    and pump levels, by the ``wet_well_defaults`` rule); the fields are not reassigned.
    """

    outfalls: dict
    pumps: dict
    height_correction: dict | None
    wet_wells: dict
    stations: dict
    station_pumps: list
    chains: dict
    force_mains: frozenset
    dwf: dict
    patterns: dict
    dwf_patterns: dict
    flow_units: str
    inlets: dict


def _resolve_inputs(
    subnetwork,
    outfalls,
    height_correction,
    pumps,
    wet_wells,
    pump_stations,
    default_pump,
    dwf,
    patterns,
    dwf_patterns,
    flow_units,
):
    """The :class:`ModelInputs` of the options of :func:`build_model`; raise ValueError for an
    invalid option or two options that contradict each other."""
    outfalls = outfall_types(outfalls, subnetwork)
    unknown = sorted(set(outfalls) - set(subnetwork.nodes))
    if unknown:
        raise ValueError(f"Outfalls not in the subnetwork: {', '.join(map(str, unknown))}")
    pumps = resolve_pumps(subnetwork, pumps)
    if height_correction is not None:
        height_correction = subnetwork.match_nodes(height_correction, "height_correction")
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
    force_mains = frozenset(pipe for chain in chains.values() for pipe in chain)
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
    dwf_patterns = resolve_dwf_patterns(dwf_patterns, dwf, patterns, subnetwork)
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
    return ModelInputs(
        outfalls=outfalls,
        pumps=pumps,
        height_correction=height_correction,
        wet_wells=wet_wells,
        stations=stations,
        station_pumps=station_pumps,
        chains=chains,
        force_mains=force_mains,
        dwf=dwf,
        patterns=patterns,
        dwf_patterns=dwf_patterns,
        flow_units=flow_units,
        inlets=inlets,
    )


def _refuse_remaining_gaps(subnetwork, filled):
    """Raise ValueError listing the gaps that stop the subnetwork from being written."""
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


def _apply_design_rules(subnetwork, inputs, pipes, inverts, wet_well_defaults):
    """Apply ``wet_well_defaults`` (filling ``inputs.wet_wells`` and the stations), check the
    pumps' wet wells and levels, and return the :class:`Derived` values."""
    stations, wet_wells = inputs.stations, inputs.wet_wells
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
            inputs.flow_units,
            inputs.force_mains,
            pump_names={
                node: [link.name for link in inputs.station_pumps if link.station == node]
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
        for name, pump in inputs.pumps.items()
        if (problem := pump_level_problem(name, pump, wet_wells.get(inputs.inlets[name])))
    ]
    if problems:
        raise ValueError(Problems("Invalid pumps", problems).message())
    return derived


def _add_nodes(model, subnetwork, inputs, pipes, inverts, defaults):
    """Add the junctions, outfalls, storage nodes and coordinates of ``subnetwork``, and each
    pump station's discharge node (whose invert is added to ``inverts``) to ``model``. Return
    {first pipe of a force main: discharge node}."""
    outfalls, wet_wells, height_correction = (
        inputs.outfalls,
        inputs.wet_wells,
        inputs.height_correction,
    )
    outlets = dict(zip(pipes["pipe_id"], pipes["to_node"], strict=True))
    pressurised = {outlets[p] for chain in inputs.chains.values() for p in chain[:-1]}
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
            depth = defaults.max_depth
            rim = subnetwork.graph.nodes[node].get("rim_elevation")
            if rim is not None and rim > invert:
                depth = round(rim - invert, 4)
            elif rim is not None:
                warnings.warn(
                    f"Node {node}: rim elevation {rim} is not above the invert {invert}; max "
                    f"depth {defaults.max_depth} is used.",
                    DataWarning,
                    stacklevel=3,  # the caller of build_model
                )
            surcharge = defaults.surcharge_depth
            if node in pressurised:
                surcharge = max(defaults.surcharge_depth, defaults.force_main_surcharge_depth)
            model.junctions.append(
                Junction(node, invert, depth, defaults.init_depth, surcharge, defaults.ponding_area)
            )

    positions = subnetwork.node_positions()
    for node in subnetwork.nodes:
        if node in positions:
            model.coordinates[node] = positions[node]
        else:
            logger.warning("Node %s has no position; it is left out of [COORDINATES]", node)

    us_inverts = dict(zip(pipes["pipe_id"], pipes["us_invert"], strict=True))
    discharge_of = {}
    for link in inputs.station_pumps:
        ends = [_invert(us_inverts[p]) for p in link.force_mains]
        ends = [e for e in ends if e is not None]
        inverts[link.discharge] = min(ends) if ends else inverts[link.station]
        surcharge = max(defaults.surcharge_depth, defaults.force_main_surcharge_depth)
        model.junctions.append(
            Junction(
                link.discharge,
                inverts[link.discharge],
                defaults.max_depth,
                defaults.init_depth,
                surcharge,
                defaults.ponding_area,
            )
        )
        if link.station in positions:
            model.coordinates[link.discharge] = positions[link.station]
        discharge_of.update(dict.fromkeys(link.force_mains, link.discharge))
    return discharge_of


def _add_links(model, inputs, pipes, inverts, discharge_of, sections, roughness_of, defaults):
    """Add the conduits, pumps, cross-sections, curves and vertices to ``model``. Return the
    descriptions of the conduits with an adverse slope (force mains excepted)."""
    pumps, force_mains = inputs.pumps, inputs.force_mains
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

        roughness = defaults.manning_n_by_material.get(row["material"], defaults.manning_n)

        ends = []
        for node, pipe_invert in ((inlet, row["us_invert"]), (outlet, row["ds_invert"])):
            pipe_invert = _invert(pipe_invert)
            offset = 0.0 if pipe_invert is None else round(pipe_invert - inverts[node], 4)
            if offset < 0:
                warnings.warn(
                    f"Pipe {name}: invert {pipe_invert} is below the invert of node {node} "
                    f"({inverts[node]}); its offset there is set to 0.",
                    DataWarning,
                    stacklevel=3,  # the caller of build_model
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
                defaults.init_flow,
                defaults.max_flow,
            )
        )
        if name in force_mains:
            diameter = sections[name][1]
            model.xsections.append(
                XSection(name, "FORCE_MAIN", diameter, roughness_of[name], defaults.barrels)
            )
        else:
            model.xsections.append(XSection(name, *sections[name], defaults.barrels))

    for link in inputs.station_pumps:
        pump = inputs.stations[link.station].pump
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
    return adverse


def _check_ideal_pumps(model):
    """Raise ValueError if an ideal pump is not the only outlet link of its inlet node (SWMM's
    ERROR 134)."""
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


def build_model(
    subnetwork: Subnetwork,
    outfalls: Outfalls | None = None,
    max_depth: float = 0.0,
    init_depth: float = 0.0,
    surcharge_depth: float = 0.0,
    ponding_area: float = 0.0,
    manning_n: float = 0.01,
    manning_n_by_material: Mapping[Any, float] | None = None,
    init_flow: float = 0.0,
    max_flow: float = 0.0,
    barrels: float = 1.0,
    height_correction: Mapping[NodeId, float] | None = None,
    pumps: Mapping[PipeId, Pump | Fields] | None = None,
    wet_wells: Mapping[NodeId, WetWell | Fields] | None = None,
    pump_stations: Mapping[NodeId, PumpStation | Fields] | None = None,
    default_pump: Pump | Fields | None = None,
    force_main_roughness: float | Mapping[str, Any] | Callable[..., float | None] | None = None,
    force_main_roughness_by_pipe: Mapping[PipeId, float] | None = None,
    force_main_surcharge_depth: float = 100.0,
    dwf: Mapping[NodeId, float] | None = None,
    patterns: Mapping[str, Any] | None = None,
    dwf_patterns: str | Sequence[str] | Mapping[NodeId, str | Sequence[str]] | None = None,
    wet_well_defaults: Mapping[str, float] | None = None,
    force_main_checks: Mapping[str, Any] | None = None,
    flow_units: str = "LPS",
    **fills: Any,
) -> SwmmModel:
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

        The node keys of this option and of ``height_correction``, ``wet_wells``,
        ``pump_stations``, ``dwf`` and ``dwf_patterns`` match a node ID exactly or by its text,
        so ``2151735`` and ``'2151735'`` are the same node (see :meth:`Subnetwork.match_node`).

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
    ------
    ValueError
        If an outfall is not in the subnetwork, a gap that stops the model from being written
        remains after the fills (a conduit with no cross-section, a node with no invert; every
        such gap is listed), a pump, wet well, pump station, dry-weather flow, pattern or fill
        rule is invalid, a pipe is both pumped and part of a station's force main, a force main
        has no roughness or is not circular, a pump's startup depth is not above its shutoff depth
        or not below its wet well's max depth, a PUMP1 station has no wet well, or an ideal pump
        is not the only outlet link of its inlet node (SWMM's ERROR 134).
    """
    defaults = _Defaults(
        max_depth=max_depth,
        init_depth=init_depth,
        surcharge_depth=surcharge_depth,
        ponding_area=ponding_area,
        force_main_surcharge_depth=force_main_surcharge_depth,
        manning_n=manning_n,
        manning_n_by_material=manning_n_by_material or {},
        init_flow=init_flow,
        max_flow=max_flow,
        barrels=barrels,
    )
    inputs = _resolve_inputs(
        subnetwork,
        outfalls,
        height_correction,
        pumps,
        wet_wells,
        pump_stations,
        default_pump,
        dwf,
        patterns,
        dwf_patterns,
        flow_units,
    )
    filled = fill_gaps(
        subnetwork, pumps=inputs.pumps, height_correction=inputs.height_correction, **fills
    )
    _refuse_remaining_gaps(subnetwork, filled)
    pipes, inverts = filled.pipes, filled.inverts
    derived = _apply_design_rules(subnetwork, inputs, pipes, inverts, wet_well_defaults)

    sections = {
        name: xsection(row)
        for name, row in zip(pipes["pipe_id"], pipes.to_dict("records"), strict=True)
    }
    roughness_of = force_main_roughnesses(
        pipes[pipes["pipe_id"].isin(list(inputs.force_mains))],
        sections,
        force_main_roughness,
        subnetwork.match_pipes(
            force_main_roughness_by_pipe or {}, "force_main_roughness_by_pipe", unknown="ignore"
        ),
    )
    model = SwmmModel(data_gaps=filled.report, derived=derived)
    discharge_of = _add_nodes(model, subnetwork, inputs, pipes, inverts, defaults)
    adverse = _add_links(
        model, inputs, pipes, inverts, discharge_of, sections, roughness_of, defaults
    )
    _check_ideal_pumps(model)

    # Warned here, not in the helpers, so that they point at the caller of build_model.
    if force_main_checks is not None:
        for message in check_force_mains(
            inputs.station_pumps, inputs.stations, sections, force_main_checks, inputs.flow_units
        ):
            warnings.warn(message, DesignWarning, stacklevel=2)
    if adverse:
        warnings.warn(
            f"Conduits with an adverse slope (outlet invert above inlet invert): "
            f"{', '.join(adverse)}.",
            DataWarning,
            stacklevel=2,
        )

    model.dwf = [
        DryWeatherFlow(node, flow, tuple(inputs.dwf_patterns.get(node, ())))
        for node, flow in inputs.dwf.items()
    ]
    model.patterns = [
        Pattern(name, kind, values) for name, (kind, values) in inputs.patterns.items()
    ]
    return model


def write_inp(
    subnetwork: Subnetwork,
    path: StrPath,
    outfalls: Outfalls | None = None,
    template_path: StrPath | None = None,
    dead_ends: str = "raise",
    **options: Any,
) -> SwmmModel:
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
    -------
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
