import contextlib
import io
import re
import shutil
from pathlib import Path

import networkx as nx
import pytest

from sewernetworkx import Network, SubNetwork

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
FIXTURE = TESTS / "data" / "fixture_network.geojson"
TEMPLATE = REPO / "examples" / "pyswmm_template" / "template.inp"
EXAMPLE_INP = REPO / "examples" / "inp_test_file" / "test.inp"

# Dendritic part of the fixture (no split, no adverse GIS pipe), as (upstream, downstream).
DENDRITIC = [
    ("A1", "A2"),
    ("A2", "J1"),
    ("B1", "B2"),
    ("B2", "J1"),
    ("C1", "J1"),
    ("J1", "J2"),
    ("J2", "OUT"),
]
OVERFLOW = [("J2", "OV1"), ("OV1", "OV2")]


def quiet(fn, *args, **kwargs):
    """Call fn with the library's print output suppressed."""
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


@pytest.fixture
def network():
    return quiet(Network, str(FIXTURE))


def subnetwork_with_edges(network, edges, outfalls=("OUT",), start="J1"):
    """A SubNetwork whose subgraph is set directly, bypassing extraction, so the .inp writer can be
    tested on a known network. SewerNetworkx stores edges as (downstream, upstream)."""
    sub = quiet(
        SubNetwork,
        {
            "network_data": network,
            "starting_node": start,
            "outfall_node": list(outfalls),
            "depth": 1,
            "outfall_node_type": ["FREE"] * len(outfalls),
            "scale": 1.0,
        },
    )
    sub.subgraph = nx.DiGraph([(down, up) for up, down in edges])
    return sub


def write_inp(sub, tmp_path, name="model.inp", **kwargs):
    quiet(
        sub.write_pyswmm,
        output_file=name,
        template_path=str(TEMPLATE),
        output_folder=str(tmp_path),
        **kwargs,
    )
    return tmp_path / name


def read_sections(path):
    """Parse an .inp into {SECTION: [token rows]}, skipping ;; comments and blank lines."""
    sections, current = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.split(";;")[0].strip()
        if not line:
            continue
        m = re.fullmatch(r"\[(\w+)\]", line)
        if m:
            current = m.group(1).upper()
            sections.setdefault(current, [])
        elif current:
            sections[current].append(line.split())
    return sections


def run_swmm(inp_path, tmp_path):
    """Run an .inp in pyswmm from a copy in tmp_path; return the ERROR lines from the report."""
    from pyswmm import Simulation

    local = Path(tmp_path) / Path(inp_path).name
    if Path(inp_path).resolve() != local.resolve():
        shutil.copy(inp_path, local)
    try:
        with Simulation(str(local)) as sim:
            for _ in sim:
                pass
    except Exception:
        pass
    rpt = local.with_suffix(".rpt")
    text = rpt.read_text() if rpt.exists() else "ERROR: no report written"
    return re.findall(r"ERROR \d+:[^\n]*", text) or (
        ["ERROR: no report written"] if not rpt.exists() else []
    )
