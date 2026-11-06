"""Building the :class:`~swmmbuilder.swmm.SwmmModel` of a subnetwork and writing its .inp file."""

from __future__ import annotations

import inspect
import logging
import math
import warnings
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from ..design import (
    Derived,
    apply_wet_well_defaults,
    check_force_main_equation,
    check_force_mains,
    force_main_roughnesses,
    pump_level_problem,
)
from ..diagnostics import DataWarning, DesignWarning, Problems, check_options
from ..gaps import _invert, fill_gaps, height_corrections, lowest_inverts, xsection
from ..geometry import line_coords
from ..units import _FLOW_UNIT_LPS, _is_real
from .dwf import resolve_dwf, resolve_dwf_patterns, resolve_patterns
from .inflows import resolve_inflows, resolve_timeseries
from .pumping import resolve_pump_stations, resolve_pumps, resolve_wet_wells
from .records import (
    Conduit,
    Curve,
    DryWeatherFlow,
    ExternalInflow,
    Junction,
    Outfall,
    Pattern,
    PumpLink,
    Storage,
    SwmmModel,
    TimeSeries,
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
    from .inflows import Inflow

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
    timeseries: dict
    inflows: dict
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
    timeseries,
    inflows,
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
    timeseries = resolve_timeseries(timeseries)
    inflows = resolve_inflows(subnetwork, inflows, timeseries, patterns, outfalls)
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
        timeseries=timeseries,
        inflows=inflows,
        flow_units=flow_units,
        inlets=inlets,
    )


def _lowest_inlet_inverts(subnetwork, inputs, slope):
    """({station: invert}, [Derived]) of the ``station_invert='lowest_inlet'`` rule of
    :func:`build_model`, for the stations of ``inputs`` without a ``height_correction``."""
    pipes = subnetwork.pipes()
    given = set(inputs.height_correction or {})
    stations = [n for n in inputs.stations if n not in given]
    current = lowest_inverts(pipes, set(stations))
    found, derived = {}, []
    for node in stations:
        ends = []
        arriving = pipes[pipes["to_node"] == node]
        for pipe, us, ds, length in zip(
            arriving["pipe_id"],
            arriving["us_invert"],
            arriving["ds_invert"],
            arriving["length"],
            strict=True,
        ):
            end, how = _invert(ds), "recorded downstream invert"
            if end is None and slope is not None and _invert(us) is not None:
                end, how = (
                    _invert(us) - slope * float(length),
                    f"upstream invert - {slope} x length",
                )
            if end is not None:
                ends.append((end, pipe, how))
        if not ends:
            continue
        end, pipe, how = min(ends)
        before = current.get(node)
        if before is None or end < before:
            found[node] = end
            basis = f"pipe {pipe} ({how}); invert from its pipes " + (
                "none" if before is None else f"{before:.3f}"
            )
            derived.append(
                Derived(
                    "station", node, "invert", found[node], "station_invert='lowest_inlet'", basis
                )
            )
            logger.info("Station %s: invert %s -> %.3f (%s)", node, before, end, basis)
    return found, derived


def junction_depths(subnetwork, inverts, exclude=()) -> pd.Series:
    """Rim - invert of each node of ``subnetwork`` (not in ``exclude``) whose rim elevation is
    known and above its invert in ``inverts``, rounded as written, indexed by node."""
    depths = {}
    for node in subnetwork.nodes:
        if node in exclude:
            continue
        rim = subnetwork.graph.nodes[node].get("rim_elevation")
        if rim is not None and not (isinstance(rim, float) and math.isnan(rim)):
            if rim > inverts[node]:
                depths[node] = round(rim - inverts[node], 4)
    return pd.Series(depths, dtype=float)


def _default_depth(subnetwork, inputs, inverts, rule):
    """(depth, junctions given it, basis) of a ``max_depth`` given as 'median' or a function."""
    exclude = set(inputs.outfalls) | set(inputs.wet_wells)
    depths = junction_depths(subnetwork, inverts, exclude)
    if rule == "median":
        if depths.empty:
            raise ValueError(
                "max_depth='median': no junction has a rim above its invert to take a median of."
            )
        depth = float(depths.median())
    elif callable(rule):
        depth = rule(depths)
    else:
        raise ValueError(f"max_depth must be a number, 'median', or a function, not {rule!r}.")
    if not (_is_real(depth) and depth >= 0):
        raise ValueError(f"max_depth: the rule gave {depth!r}, not a depth >= 0.")
    depth = round(float(depth), 4)
    defaulted = [n for n in subnetwork.nodes if n not in exclude and n not in depths.index]
    basis = (
        f"median rim - invert of {len(depths)} junctions"
        if rule == "median"
        else (f"function of rim - invert of {len(depths)} junctions")
    )
    logger.info(
        "max_depth %s: %.4g m, given to %d junctions with no usable rim",
        rule if rule == "median" else "function",
        depth,
        len(defaulted),
    )
    return depth, defaulted, basis


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


# The dummy conduit of the outfall_link option of build_model.
OUTFALL_LINK_DEFAULTS = {"length": 30.0, "diameter": 3.0, "slope": 0.01, "roughness": 0.01}


def outfall_inlets(model: SwmmModel) -> dict[str, list[str]]:
    """{outfall: names of the conduits and pumps that enter or leave it} of the outfalls of
    ``model`` with more than one such link (SWMM's ERROR 141 under DYNWAVE)."""
    names = {o.name for o in model.outfalls}
    links: dict[str, list[str]] = {}
    for link in model.conduits + model.pumps:
        for end in (link.outlet, link.inlet):
            if end in names:
                links.setdefault(end, []).append(link.name)
    return {node: found for node, found in links.items() if len(found) > 1}


def _add_outfall_links(model, subnetwork, option, defaults):
    """Write each outfall of ``model`` with more than one link as a junction draining to a new
    outfall through a dummy conduit (the ``outfall_link`` option of :func:`build_model`). Return
    the junctions given ``defaults.max_depth`` (no usable rim)."""
    params = dict(OUTFALL_LINK_DEFAULTS)
    if isinstance(option, Mapping):
        check_options(option, tuple(OUTFALL_LINK_DEFAULTS), "outfall_link")
        params.update(option)
    bad = [f"{k}={v!r}" for k, v in params.items() if not (_is_real(v) and v > 0)]
    if bad:
        raise ValueError(f"outfall_link values must be positive numbers: {', '.join(bad)}")
    crowded = outfall_inlets(model)
    taken = {*subnetwork.nodes, *(c.name for c in model.conduits), *(p.name for p in model.pumps)}
    clash = sorted(
        name for node in crowded for name in (f"{node}_OUT", f"{node}_DUMMY") if name in taken
    )
    if clash:
        raise ValueError(f"outfall_link: names already used in the model: {', '.join(clash)}")
    given_default = []
    for node in crowded:
        outfall = next(o for o in model.outfalls if o.name == node)
        model.outfalls.remove(outfall)
        depth = defaults.max_depth
        rim = subnetwork.graph.nodes[node].get("rim_elevation")
        if rim is not None and rim > outfall.invert:
            depth = round(rim - outfall.invert, 4)
        else:
            given_default.append(node)
        model.junctions.append(
            Junction(
                node,
                outfall.invert,
                depth,
                defaults.init_depth,
                defaults.surcharge_depth,
                defaults.ponding_area,
            )
        )
        new = f"{node}_OUT"
        drop = params["slope"] * params["length"]
        model.outfalls.append(Outfall(new, round(outfall.invert - drop, 6), outfall.type))
        link = f"{node}_DUMMY"
        model.conduits.append(
            Conduit(link, node, new, params["length"], params["roughness"], 0.0, 0.0, 0.0, 0.0)
        )
        model.xsections.append(XSection(link, "CIRCULAR", params["diameter"], 0.0, 1.0))
        if node in model.coordinates:
            x, y = model.coordinates[node]
            model.coordinates[new] = (x + 1.0, y)
        logger.info(
            "Outfall %s (%d links) written as a junction draining through %s (%g m) to %s",
            node,
            len(crowded[node]),
            link,
            params["length"],
            new,
        )
    return given_default


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
    max_depth: float | str | Callable[[pd.Series], float] = 0.0,
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
    inflows: Mapping[NodeId, Inflow | Fields] | None = None,
    timeseries: Mapping[str, Sequence[tuple[float, float]]] | None = None,
    wet_well_defaults: Mapping[str, float] | None = None,
    force_main_checks: Mapping[str, Any] | None = None,
    outfall_link: bool | Mapping[str, float] | None = None,
    station_invert: str | None = None,
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
        instead of ``max_depth``. With 0, SWMM takes the depth to the crown of the highest pipe
        at the junction, so a junction with no usable rim then overflows as soon as a pipe
        joining it runs full.

        ``max_depth`` may instead be ``'median'``: the median of rim - invert over the
        subnetwork's junctions whose rim is above their invert (the depths written), given to
        every junction with no usable rim. A function of those depths (a pandas Series indexed by
        node) returning one depth also works, e.g. ``lambda d: d.quantile(0.25)``. Each junction
        given that depth, including pump-station discharge nodes and outfalls written as
        junctions by ``outfall_link``, is listed in :attr:`SwmmModel.derived`.

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

    inflows : dict, optional
        {node: :class:`Inflow` or dict of its fields}: an external flow into each node, written to
        [INFLOWS], in the template's FLOW_UNITS: ``scale`` times a time series of
        ``timeseries`` plus a ``baseline`` (optionally times a pattern of ``patterns``), e.g.
        ``{'N1': dict(timeseries='commercial', scale=2.5, baseline=0.4)}``. SWMM adds it to the
        node's [DWF] flow, so a node can carry a second time-varying flow, or a steady one such
        as groundwater infiltration, besides its dry-weather flow.

    timeseries : dict, optional
        {name: [(hours from the start of the run, value), ...]}, written to [TIMESERIES]. SWMM
        interpolates linearly between the points and takes the value as zero after the last
        one; :func:`hourly_series` gives a series that holds hourly values for whole hours.

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

    station_invert : {'lowest_inlet'}, optional
        A design rule for the stations of ``pump_stations``: a station's wet well cannot sit
        above the pipes arriving at it. A station's invert is the lowest invert of its pipes, but
        where an arriving pipe has no downstream invert that end does not count, and the invert
        may then come from the force main alone, metres above the arriving sewer. With
        ``'lowest_inlet'``, each arriving pipe's downstream end is its recorded invert or, if
        missing, its upstream invert minus ``headwater_slope`` times its length (a pipe with
        neither, or no ``headwater_slope``, does not count), and a station whose invert is
        above the lowest of these ends, or has none, is lowered to it. ``height_correction``
        wins over the rule. Each change is listed in :attr:`SwmmModel.derived`. Not applied
        unless given.

    outfall_link : bool or dict, optional
        Under dynamic-wave routing SWMM accepts only one link into an outfall (its ERROR 141), so
        :func:`write_inp` refuses a DYNWAVE model with an outfall that has more than one. With
        this option, each such outfall is written instead as a junction (its depth set like any
        other junction's), joined by a short, wide dummy conduit ``<node>_DUMMY`` to a new
        outfall ``<node>_OUT`` of the same type, 1 m east of it and lower by the conduit's
        slope times its length. True uses :data:`OUTFALL_LINK_DEFAULTS` (a 30 m, 3 m diameter
        circular conduit at a slope of 0.01 with Manning's n 0.01; a much shorter conduit sets
        the dynamic-wave time step); a dict overrides any of ``length``, ``diameter``, ``slope``,
        and ``roughness``. Outfalls with one inlet link are left as they are.

    flow_units : str
        The template's FLOW_UNITS, in which pump curves and ``design_flow`` are given (default
        LPS, the shipped template's; :func:`write_inp` sets it from the template).

    Raises
    ------
    ValueError
        If an outfall is not in the subnetwork, a gap that stops the model from being written
        remains after the fills (a conduit with no cross-section, a node with no invert; every
        such gap is listed), a pump, wet well, pump station, dry-weather flow, pattern, inflow,
        time series, or fill rule is invalid, a pipe is both pumped and part of a station's force main, a force main
        has no roughness or is not circular, a pump's startup depth is not above its shutoff depth
        or not below its wet well's max depth, a PUMP1 station has no wet well, or an ideal pump
        is not the only outlet link of its inlet node (SWMM's ERROR 134).
    """
    check_options(fills, BUILD_OPTIONS, "build_model")
    defaults = _Defaults(
        # A rule ('median' or a function) is resolved once the inverts are known.
        max_depth=float(max_depth) if _is_real(max_depth) else 0.0,
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
        timeseries,
        inflows,
        flow_units,
    )
    lowered = []
    if station_invert is not None:
        if station_invert != "lowest_inlet":
            raise ValueError(
                f"station_invert must be 'lowest_inlet' or None, not {station_invert!r}."
            )
        found, lowered = _lowest_inlet_inverts(subnetwork, inputs, fills.get("headwater_slope"))
        if found:
            inputs = replace(
                inputs, height_correction={**found, **(inputs.height_correction or {})}
            )
    filled = fill_gaps(
        subnetwork, pumps=inputs.pumps, height_correction=inputs.height_correction, **fills
    )
    _refuse_remaining_gaps(subnetwork, filled)
    pipes, inverts = filled.pipes, filled.inverts
    derived = _apply_design_rules(subnetwork, inputs, pipes, inverts, wet_well_defaults)
    derived += lowered
    depth_rule = None
    if not _is_real(max_depth):
        depth, defaulted, basis = _default_depth(subnetwork, inputs, inverts, max_depth)
        defaults = replace(defaults, max_depth=depth)
        rule = "max_depth='median'" if max_depth == "median" else "max_depth function"
        depth_rule = (depth, rule, basis)
        derived += [
            Derived("junction", node, "max_depth", depth, rule, basis) for node in defaulted
        ]

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
    linked = []
    if outfall_link:
        linked = _add_outfall_links(model, subnetwork, outfall_link, defaults)
    if depth_rule is not None:
        # The other junctions given the rule's depth: station discharge nodes, and outfalls
        # written as junctions with no usable rim.
        depth, rule, basis = depth_rule
        model.derived.extend(
            Derived("junction", node, "max_depth", depth, rule, basis)
            for node in [link.discharge for link in inputs.station_pumps] + linked
        )

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
    model.inflows = [
        ExternalInflow(node, i.timeseries, float(i.scale), float(i.baseline), i.baseline_pattern)
        for node, i in inputs.inflows.items()
    ]
    model.timeseries = [TimeSeries(name, points) for name, points in inputs.timeseries.items()]
    return model


# The fill rules passed on to fill_gaps through **fills, and every option of build_model (so a
# misspelt option, e.g. size_overide, is named with a suggestion).
FILL_OPTIONS = tuple(
    name
    for name in inspect.signature(fill_gaps).parameters
    if name not in ("subnetwork", "pumps", "height_correction")
)
BUILD_OPTIONS = (
    *(
        name
        for name, parameter in inspect.signature(build_model).parameters.items()
        if name not in ("subnetwork", "outfalls") and parameter.kind is not parameter.VAR_KEYWORD
    ),
    *FILL_OPTIONS,
)


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
        ``wet_wells``, ``pumps`` and ``pump_stations`` options, {dwf} and {patterns} only with
        ``dwf`` and ``patterns``, and {inflows} and {timeseries} only with ``inflows`` and
        ``timeseries``. The offsets are depths, so it should set LINK_OFFSETS DEPTH
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
    check_options(options, BUILD_OPTIONS, "write_inp")  # before any other work
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
    if template_routing(template) == "DYNWAVE":
        crowded = outfall_inlets(model)
        if crowded:
            raise ValueError(
                "Under DYNWAVE routing an outfall may have only one link (SWMM ERROR 141); pass "
                "outfall_link=True to write each of these as a junction draining to a new "
                "outfall: "
                + "; ".join(f"{n} ({', '.join(map(str, links))})" for n, links in crowded.items())
            )
    needed = [
        name
        for name, records in (
            ("storage", model.storage),
            ("pumps", model.pumps),
            ("curves", model.curves),
            ("dwf", model.dwf),
            ("patterns", model.patterns),
            ("inflows", model.inflows),
            ("timeseries", model.timeseries),
        )
        if records and f"{{{name}}}" not in template
    ]
    if needed:
        raise ValueError(
            f"The template has no {', '.join(f'{{{n}}}' for n in needed)} placeholder, so the "
            "pumps, pump stations, wet wells, dry-weather flows, patterns, inflows, or time series "
            "would not be written."
        )
    if model.dwf:
        logger.info(
            "Dry-weather flow: total base flow %.6g %s at %d nodes",
            sum(d.base for d in model.dwf),
            template_flow_units(template),
            len(model.dwf),
        )
    if model.inflows:
        logger.info(
            "External inflows: total baseline %.6g %s at %d nodes, %d time series",
            sum(i.baseline for i in model.inflows),
            template_flow_units(template),
            len(model.inflows),
            len(model.timeseries),
        )
    logger.info("Writing %s", path)
    Path(path).write_text(model.render(template))
    return model
