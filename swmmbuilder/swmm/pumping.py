"""The ``pumps``, ``wet_wells`` and ``pump_stations`` options of
:func:`~swmmbuilder.swmm.build_model`, resolved against a subnetwork."""

import warnings

from ..diagnostics import OptionWarning, Problems
from ..pumps import Pump, PumpStation, StationPump, WetWell


def resolve_pumps(subnetwork, pumps):
    """{pipe_id: Pump} from the ``pumps`` option of :func:`build_model`.

    The keys are matched with :meth:`Subnetwork.match_pipe`, so ``1010`` and ``'1010'`` are the
    same pipe; the result is keyed by the pipe IDs of the data.

    Raises
    ------
    ValueError
        If a pipe is not in the subnetwork or a pump is invalid.
    """
    if not pumps:
        return {}
    return {
        pipe: Pump.coerce(pipe, value, curve_type="PUMP3")
        for pipe, value in subnetwork.match_pipes(pumps, "pumps").items()
    }


def resolve_wet_wells(subnetwork, wet_wells, outfalls=()):
    """{node: WetWell} from the ``wet_wells`` option of :func:`build_model`.

    Raises
    ------
    ValueError
        If a node is not in the subnetwork or is an outfall, or a wet well is invalid.
    """
    if not wet_wells:
        return {}
    wet_wells = subnetwork.match_nodes(wet_wells, "wet_wells")
    missing = sorted(str(n) for n in wet_wells if n not in set(subnetwork.nodes))
    if missing:
        raise ValueError(f"Wet wells not in the subnetwork: {', '.join(missing)}")
    at_outfalls = sorted(str(n) for n in wet_wells if n in set(outfalls))
    if at_outfalls:
        raise ValueError(f"Wet wells cannot be outfalls: {', '.join(at_outfalls)}")
    return {node: WetWell.coerce(node, v) for node, v in wet_wells.items()}


def resolve_pump_stations(subnetwork, pump_stations, outfalls=(), default_pump=None, excluded=()):
    """({node: PumpStation}, [StationPump]) from the ``pump_stations`` option of
    :func:`build_model`. Station nodes in ``excluded`` (left out on purpose, by
    ``pump_stations_exclude``) are not warned of as missing from ``pump_stations``.

    Raises
    ------
    ValueError
        If a station node is not in the subnetwork, is an outfall or has no force-main outlet, a
        chosen force main does not leave it, it has several force mains and no ``discharge``,
        a pump is missing or invalid, a wet well is invalid, or a name the station needs
        (``<node>_P``, ``<node>_D``) is already a node or pipe ID. Every problem is listed.
    """
    if not pump_stations:
        return {}, []
    pump_stations = subnetwork.match_nodes(pump_stations, "pump_stations")
    nodes = set(subnetwork.nodes)
    taken = {str(i) for i in subnetwork.network.ids()}
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
            matched = [subnetwork.match_pipe(p) for p in station.force_mains]
            wrong = [str(p) for p in matched if p not in outlets]
            if wrong:
                problems.append(
                    f"{node}: force mains {', '.join(wrong)} do not leave it "
                    f"(its force mains: {', '.join(map(str, outlets))})"
                )
                continue
            chosen = matched
            unused = [str(p) for p in outlets if p not in chosen]
            if unused:
                warnings.warn(
                    f"Pump station {node}: force mains {', '.join(unused)} are not pumped and "
                    "are written as ordinary conduits from the station node.",
                    OptionWarning,
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
        raise ValueError(Problems("Invalid pump stations", problems).message())
    undeclared = sorted(
        str(n) for n in typed - set(stations) - set(excluded) if subnetwork.force_main_outlets(n)
    )
    if undeclared:
        warnings.warn(
            f"Pump-station nodes with force mains not in pump_stations: {', '.join(undeclared)}; "
            "their force mains are written as ordinary conduits.",
            OptionWarning,
            stacklevel=3,
        )
    return stations, links
