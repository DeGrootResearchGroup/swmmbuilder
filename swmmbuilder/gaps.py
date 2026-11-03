"""Gaps in the pipe and node data, and the rules that may fill them.

A gap is a pipe with no usable cross-section, a pipe end with no invert, or a node with no
invert. :func:`fill_gaps` finds every gap of a subnetwork at once and applies the fill rules it is
given. Nothing is filled by default except one unambiguous inference: a round pipe with no
diameter whose height equals its width takes that value as its diameter. Every filled value is
logged and flagged in the report's ``filled_from`` column, so nothing filled looks like source
data.

The gap kinds (:data:`GAP_KINDS`):

==================  =========================================================================
``unknown_shape``   pipe whose shape code is missing or not in the schema
``no_size``         pipe whose shape is known but a size it needs is missing or not positive
``no_us_invert``    pipe with no upstream invert (written at the invert of its upstream node)
``no_ds_invert``    pipe with no downstream invert (written at the invert of its downstream node)
``no_node_invert``  node with no invert from any pipe
==================  =========================================================================

``unknown_shape``, ``no_size`` and ``no_node_invert`` stop a subnetwork from being written
(:data:`BLOCKING_GAPS`); a pipe end with no invert is attached to its node with no offset.
"""

import heapq
import logging
import warnings
from dataclasses import dataclass

import pandas as pd

from .diagnostics import OptionWarning

logger = logging.getLogger(__name__)

GAP_KINDS = ("unknown_shape", "no_size", "no_us_invert", "no_ds_invert", "no_node_invert")
BLOCKING_GAPS = frozenset({"unknown_shape", "no_size", "no_node_invert"})
REPORT_COLUMNS = ("element", "id", "gap", "detail", "filled_from", "value", "source", "order")

SIZE_FILLS = ("downstream", "upstream")
INVERT_FILLS = ("interpolate",)
SECTION_COLUMNS = ("shape", "diameter", "height", "width")


def xsection(pipe):
    """SWMM cross-section of a pipe, from its canonical shape and sizes (in m).

    CIRCULAR uses the diameter and EGG the full height; the other shapes use the full height
    and maximum width.

    Parameters
    ----------
    pipe : mapping
        A row of the canonical data (see :mod:`swmmbuilder.schema`).

    Returns
    -------
    tuple or None
        (shape, geom1, geom2), or None if the shape is unknown or missing, or a size it needs is
        missing or not positive.
    """
    shape = pipe.get("shape")
    if pd.isna(shape):
        return None
    sizes = [pipe.get(column) for column in _size_columns(shape)]
    if any(pd.isna(s) or s <= 0 for s in sizes):
        return None
    return shape, *sizes, *[0.0] * (2 - len(sizes))


def _size_columns(shape):
    """The size columns a cross-section of ``shape`` needs."""
    if shape == "CIRCULAR":
        return ("diameter",)
    if shape == "EGG":
        return ("height",)
    return ("height", "width")


def _invert(value):
    """An invert elevation, or None if it is missing."""
    return None if pd.isna(value) else float(value)


def lowest_inverts(pipes, nodes):
    """{node: lowest invert of the pipes connected to it} for the nodes in ``nodes``."""
    found = {}
    for up, down, inv_up, inv_down, drop in zip(
        pipes["from_node"],
        pipes["to_node"],
        pipes["us_invert"],
        pipes["ds_invert"],
        pipes["drop_invert"],
        strict=True,
    ):
        for node, value in ((up, inv_up), (down, inv_down), (down, drop)):
            value = _invert(value)
            if node in nodes and value is not None:
                found[node] = min(found.get(node, value), value)
    return found


def height_corrections(subnetwork, height_correction):
    """{node: invert} of the ``height_correction`` option for the nodes of ``subnetwork``.

    Keys are matched with :meth:`Subnetwork.match_node`. Nodes of the network outside the
    subnetwork are skipped, so one mapping can serve several subnetworks; keys that name no node
    of the network at all (likely typos) are skipped with an :class:`OptionWarning`.
    """
    matched = subnetwork.match_nodes(height_correction or {}, "height_correction")
    unknown = sorted(str(n) for n in matched if not subnetwork.in_network(n))
    if unknown:
        warnings.warn(
            f"height_correction: nodes not in the network, ignored: {', '.join(unknown)}",
            OptionWarning,
            stacklevel=3,
        )
    nodes = set(subnetwork.nodes)
    return {node: float(invert) for node, invert in matched.items() if node in nodes}


def force_main_discharge_rule(rule):
    """A function (pipe, receiving invert) -> downstream invert, from the ``force_main_discharge``
    option of :func:`fill_gaps`.

    Parameters
    ----------
    rule : dict or callable
        ``{'max_height_above_receiving_invert': h_max, 'height_above_receiving_invert': h}``: the
        force main discharges ``h`` (default ``h_max``) above the receiving node's invert, and
        ``0 <= h <= h_max``. A design standard's value can be taken from a schema's ``design``
        table, e.g. ``LONDON_ON.design['force_main_discharge']``. A callable is used as is.

    Raises
    ------
    ValueError
        If the dict has unknown keys, no ``max_height_above_receiving_invert``, or a height
        outside ``[0, h_max]``.
    """
    if callable(rule):
        return rule
    keys = {"max_height_above_receiving_invert", "height_above_receiving_invert"}
    unknown = sorted(set(rule) - keys)
    if unknown or "max_height_above_receiving_invert" not in rule:
        raise ValueError(
            f"force_main_discharge must be a callable or a dict with "
            f"max_height_above_receiving_invert (and optionally height_above_receiving_invert); "
            f"got keys {sorted(rule)}."
        )
    max_height = float(rule["max_height_above_receiving_invert"])
    height = float(rule.get("height_above_receiving_invert", max_height))
    if not 0 <= height <= max_height:
        raise ValueError(
            f"force_main_discharge: height {height} m must be between 0 and the maximum "
            f"height above the receiving invert, {max_height} m."
        )
    return lambda pipe, receiving_invert: receiving_invert + height


@dataclass
class FilledData:
    """The pipes and node inverts of a subnetwork after :func:`fill_gaps`.

    Attributes
    ----------
    pipes : pandas.DataFrame
        The subnetwork's pipes (see :meth:`Subnetwork.pipes`) with the filled values.

    inverts : dict
        {node: invert} of every node that has one, filled or not.

    report : pandas.DataFrame
        One row per gap, with the columns of :data:`REPORT_COLUMNS`: ``element`` ('pipe' or
        'node'), ``id``, ``gap`` (see :data:`GAP_KINDS`), ``detail``, ``filled_from`` (the rule
        that filled it, or NaN), ``value`` (the filled value: a shape, an (shape, geom1, geom2)
        cross-section or an invert), ``source`` (what the value came from) and ``order`` (the
        order in which the fills were made).
    """

    pipes: pd.DataFrame
    inverts: dict
    report: pd.DataFrame

    @property
    def remaining(self):
        """The rows of ``report`` with a gap in :data:`BLOCKING_GAPS` that was not filled."""
        report = self.report
        return report[report["gap"].isin(BLOCKING_GAPS) & report["filled_from"].isna()]


def _code_key(code):
    """A shape code as a key of ``shape_codes``: missing and blank codes are ''."""
    return "" if pd.isna(code) else str(code).strip()


def _pipe_gaps(pipe):
    """[(gap, detail)] of a pipe row."""
    gaps = []
    if pd.isna(pipe["shape"]):
        gaps.append(("unknown_shape", f"no cross-section: shape code {pipe['shape_code']!r}"))
    elif xsection(pipe) is None:
        needs = ", ".join(_size_columns(pipe["shape"]))
        gaps.append(("no_size", f"no cross-section: {pipe['shape']} needs {needs}"))
    for gap, column, node in (
        ("no_us_invert", "us_invert", "from_node"),
        ("no_ds_invert", "ds_invert", "to_node"),
    ):
        if pd.isna(pipe[column]):
            detail = f"at node {pipe[node]}"
            if pipe.get("structure_type") == "force_main":
                detail += " (force main)"
            gaps.append((gap, detail))
    return gaps


def fill_gaps(
    subnetwork,
    pumps=(),
    height_correction=None,
    size_overrides=None,
    shape_codes=None,
    size_fill=None,
    invert_fill=None,
    headwater_slope=None,
    force_main_discharge=None,
):
    """Find every gap in a subnetwork's pipes and node inverts, and fill those the rules cover.

    The fills are made in this order, each logged and recorded in the report: overrides
    (``size_overrides``, ``height_correction``); ``shape_codes``; the round-pipe inference;
    ``size_fill``; node inverts from ``invert_fill`` and ``headwater_slope``;
    ``force_main_discharge``. Every rule except the round-pipe inference is off by default.

    Parameters
    ----------
    subnetwork : Subnetwork
        The subnetwork.

    pumps : iterable
        Pipe IDs written as pumps, which need no cross-section: their shape and size gaps are
        left out.

    height_correction : dict, optional
        {node: invert} overriding the node inverts computed from the pipes.

    size_overrides : dict, optional
        {pipe_id: diameter (m) or dict of ``shape`` (a SWMM shape), ``diameter``, ``height`` and
        ``width`` (m)}, e.g. values obtained from the City. Pipe IDs match as strings.

    shape_codes : dict, optional
        {source shape code: SWMM shape} for pipes whose shape code is missing or not in the
        schema, e.g. ``{'': 'CIRCULAR'}``; missing and blank codes are ``''``.

    size_fill : {'downstream', 'upstream'}, optional
        Copy the sizes of a pipe with no size from an adjacent pipe of the same shape: the
        downstream pipe, else the largest upstream pipe ('downstream'), or the other way round
        ('upstream'). Runs of pipes with no size are filled one after another.

    invert_fill : {'interpolate'}, optional
        Fill a node with no invert by interpolating by length along the flow path between the
        nearest nodes with a known invert upstream and downstream (a constant slope).

    headwater_slope : float, optional
        Fill a node with no known invert upstream of it from the nearest known invert downstream,
        plus this slope (m/m) times the length between them.

    force_main_discharge : dict or callable, optional
        Fill the missing downstream invert of the last force-main pipe before the receiving node
        from the receiving node's invert (see :func:`force_main_discharge_rule`). The receiving
        node's invert must be known or filled.

    Returns
    -------
    FilledData

    Raises
    ------
    ValueError
        If an option is invalid or names a pipe that is not in the subnetwork.
    """
    if size_fill not in (None, *SIZE_FILLS):
        raise ValueError(f"size_fill must be one of {list(SIZE_FILLS)} or None.")
    if invert_fill not in (None, *INVERT_FILLS):
        raise ValueError(f"invert_fill must be one of {list(INVERT_FILLS)} or None.")
    if headwater_slope is not None and not headwater_slope >= 0:
        raise ValueError("headwater_slope must not be negative.")
    discharge = (
        None if force_main_discharge is None else force_main_discharge_rule(force_main_discharge)
    )

    pipes = subnetwork.pipes().copy()
    pumps = set(pumps)
    for column in SECTION_COLUMNS:
        pipes[column] = pipes[column].astype(object if column == "shape" else float)
    by_id = dict(zip(pipes["pipe_id"], pipes.index, strict=True))
    original = {
        ("pipe", p["pipe_id"], gap): detail
        for _, p in pipes.iterrows()
        for gap, detail in _pipe_gaps(p)
    }
    fills = {}

    def record(element, element_id, gap, rule, value, source=""):
        fills[(element, element_id, gap)] = (rule, value, source, len(fills) + 1)
        logger.info(
            "%s %s: %s filled by %s with %s%s",
            element.capitalize(),
            element_id,
            gap,
            rule,
            value,
            f" ({source})" if source else "",
        )

    for pipe_id, value in subnetwork.match_pipes(size_overrides or {}, "size_overrides").items():
        i = by_id[pipe_id]
        values = {"diameter": value} if not hasattr(value, "items") else dict(value)
        unknown = sorted(set(values) - set(SECTION_COLUMNS))
        if unknown:
            raise ValueError(
                f"size_overrides: pipe {pipe_id} sets {unknown}; use {list(SECTION_COLUMNS)}."
            )
        before = {gap for gap, _ in _pipe_gaps(pipes.loc[i])}
        for column, v in values.items():
            pipes.at[i, column] = v if column == "shape" else float(v)
        after = {gap for gap, _ in _pipe_gaps(pipes.loc[i])}
        for gap in sorted(before - after):
            record("pipe", pipe_id, gap, "size_overrides", xsection(pipes.loc[i]))
        if not before - after:
            logger.info("Pipe %s: cross-section set by size_overrides: %s", pipe_id, values)

    for i, pipe in pipes.iterrows():
        code = _code_key(pipe["shape_code"])
        if pd.isna(pipe["shape"]) and code in (shape_codes or {}):
            pipes.at[i, "shape"] = shape_codes[code]
            record(
                "pipe",
                pipe["pipe_id"],
                "unknown_shape",
                "shape_codes",
                shape_codes[code],
                f"shape code {pipe['shape_code']!r}",
            )

    for i, pipe in pipes.iterrows():
        if (
            pipe["shape"] == "CIRCULAR"
            and not pipe["diameter"] > 0
            and pipe["height"] > 0
            and pipe["height"] == pipe["width"]
        ):
            pipes.at[i, "diameter"] = pipe["height"]
            record(
                "pipe",
                pipe["pipe_id"],
                "no_size",
                "height_width",
                xsection(pipes.loc[i]),
                "height = width",
            )

    if size_fill is not None:
        _fill_sizes(pipes, size_fill, pumps, record)

    nodes = set(subnetwork.nodes)
    inverts = lowest_inverts(pipes, nodes)
    if nodes - set(inverts):
        inverts.update(subnetwork.network.lowest_inverts(nodes - set(inverts)))
    for node in nodes - set(inverts):
        original[("node", node, "no_node_invert")] = "no invert from any pipe"
    for node, invert in height_corrections(subnetwork, height_correction).items():
        if node not in inverts:
            record("node", node, "no_node_invert", "height_correction", invert)
        inverts[node] = invert
    unknown_nodes = sorted(nodes - set(inverts), key=str)
    if invert_fill is not None or headwater_slope is not None:
        _fill_node_inverts(
            subnetwork.graph, unknown_nodes, inverts, invert_fill, headwater_slope, record
        )

    if discharge is not None:
        filled_nodes = {key[1] for key in fills if key[0] == "node"}
        outlets_fm = {
            up
            for up, structure in zip(pipes["from_node"], pipes["structure_type"], strict=True)
            if structure == "force_main"
        }
        for i, pipe in pipes.iterrows():
            receiving = pipe["to_node"]
            if (
                pipe["structure_type"] != "force_main"
                or pd.notna(pipe["ds_invert"])
                or receiving in outlets_fm
                or receiving not in inverts
            ):
                continue
            value = round(float(discharge(pipe, inverts[receiving])), 4)
            pipes.at[i, "ds_invert"] = value
            source = f"node {receiving} invert {inverts[receiving]}"
            if receiving in filled_nodes:
                source += f", itself filled by {fills[('node', receiving, 'no_node_invert')][0]}"
            record("pipe", pipe["pipe_id"], "no_ds_invert", "force_main_discharge", value, source)

    details = {
        ("pipe", p["pipe_id"], gap): detail
        for _, p in pipes.iterrows()
        for gap, detail in _pipe_gaps(p)
    }
    details.update(original)
    for key in fills:
        details.setdefault(key, "")

    rows = []
    for key, detail in details.items():
        element, element_id, gap = key
        if element == "pipe" and element_id in pumps and gap in ("unknown_shape", "no_size"):
            continue
        rule, value, source, order = fills.get(key, (None, None, None, None))
        rows.append((element, element_id, gap, detail, rule, value, source, order))
    rank = {gap: i for i, gap in enumerate(GAP_KINDS)}
    rows.sort(key=lambda r: (r[0] != "pipe", str(r[1]), rank[r[2]]))
    report = pd.DataFrame(rows, columns=list(REPORT_COLUMNS))
    report["order"] = report["order"].astype("Int64")
    return FilledData(pipes, inverts, report)


def _fill_sizes(pipes, size_fill, pumps, record):
    """Copy sizes from adjacent pipes of the same shape to pipes with a shape but no size."""
    by_from, by_to = {}, {}
    for i, up, down in zip(pipes.index, pipes["from_node"], pipes["to_node"], strict=True):
        by_from.setdefault(up, []).append(i)
        by_to.setdefault(down, []).append(i)

    def candidates(i, indices):
        shape = pipes.at[i, "shape"]
        found = [
            j
            for j in indices
            if j != i and pipes.at[j, "shape"] == shape and xsection(pipes.loc[j]) is not None
        ]
        # The largest first; ties by pipe ID.
        return sorted(
            found, key=lambda j: (-xsection(pipes.loc[j])[1], str(pipes.at[j, "pipe_id"]))
        )

    changed = True
    while changed:
        changed = False
        for i, pipe in pipes.iterrows():
            if pipe["pipe_id"] in pumps or pd.isna(pipe["shape"]) or xsection(pipe) is not None:
                continue
            downstream = candidates(i, by_from.get(pipe["to_node"], []))
            upstream = candidates(i, by_to.get(pipe["from_node"], []))
            order = (
                [("downstream", downstream), ("upstream", upstream)]
                if size_fill == "downstream"
                else [("upstream", upstream), ("downstream", downstream)]
            )
            for direction, found in order:
                if found:
                    j = found[0]
                    for column in _size_columns(pipe["shape"]):
                        pipes.at[i, column] = pipes.at[j, column]
                    record(
                        "pipe",
                        pipe["pipe_id"],
                        "no_size",
                        "size_fill",
                        xsection(pipes.loc[i]),
                        f"{direction} pipe {pipes.at[j, 'pipe_id']}",
                    )
                    changed = True
                    break


def _nearest_known(graph, start, known, downstream):
    """(node, length) of the nearest node in ``known`` downstream (or upstream) of ``start``
    along the pipes of ``graph``, not passing through another known node; None if there is
    none."""
    heap, seen = [(0.0, str(start), start)], set()
    while heap:
        distance, _, node = heapq.heappop(heap)
        if node in seen:
            continue
        seen.add(node)
        if node != start and node in known:
            return node, distance
        edges = (
            graph.out_edges(node, data="length")
            if downstream
            else graph.in_edges(node, data="length")
        )
        for up, down, length in edges:
            other = down if downstream else up
            if other not in seen:
                heapq.heappush(heap, (distance + float(length), str(other), other))
    return None


def _fill_node_inverts(graph, unknown, inverts, invert_fill, headwater_slope, record):
    """Fill the inverts of the ``unknown`` nodes from the known ones (not from each other, so a
    run of nodes between two known inverts lies on one slope)."""
    known = {n: v for n, v in inverts.items() if n in graph}
    for node in unknown:
        down = _nearest_known(graph, node, known, downstream=True)
        up = _nearest_known(graph, node, known, downstream=False)
        if down is None:
            continue
        (down_node, down_length) = down
        if up is not None and invert_fill == "interpolate":
            up_node, up_length = up
            total = up_length + down_length
            fraction = down_length / total if total > 0 else 0.0
            value = known[down_node] + (known[up_node] - known[down_node]) * fraction
            rule = "interpolate"
            source = (
                f"nodes {up_node} ({up_length:.1f} m up) and {down_node} ({down_length:.1f} m down)"
            )
        elif up is None and headwater_slope is not None:
            value = known[down_node] + headwater_slope * down_length
            rule = "headwater_slope"
            source = f"node {down_node} ({down_length:.1f} m down) at slope {headwater_slope}"
        else:
            continue
        inverts[node] = round(value, 4)
        record("node", node, "no_node_invert", rule, inverts[node], source)
