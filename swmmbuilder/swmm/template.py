""".inp templates: reading and checking them, their options, and the topology they
accept."""

import re
from importlib import resources
from pathlib import Path

from ..units import SI_FLOW_UNITS

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
    ------
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
            problems.append(f"cycle {' -> '.join(map(str, cycle + cycle[:1]))}")
    if dead_ends == "raise" and report.dead_ends:
        problems.append(
            f"dead ends with no outlet that are not outfalls: "
            f"{', '.join(map(str, report.dead_ends))} "
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
