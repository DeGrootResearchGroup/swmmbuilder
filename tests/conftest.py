import re
import shutil
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pytest

from swmmbuilder import Network

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


@pytest.fixture
def network():
    return Network(FIXTURE)


def nodes_of(edges):
    return {n for edge in edges for n in edge}


def subnetwork(network, edges):
    """The subnetwork of the nodes of ``edges`` (upstream, downstream), checking that its pipes
    join exactly those node pairs."""
    sub = network.subnetwork(nodes_of(edges))
    assert set(sub.graph.edges()) == set(edges)
    return sub


def with_pipes(*pairs):
    """The fixture network plus a copy of pipe 1000 (A1 -> A2) for each (upstream, downstream)
    pair, with pipe IDs from 2000; e.g. to close a cycle."""
    data = gpd.read_file(FIXTURE)
    template = data[data["OBJECTID"] == 1000]
    extra = [
        template.assign(OBJECTID=2000 + i, UpstreamID=up, DownstreamID=down)
        for i, (up, down) in enumerate(pairs)
    ]
    return Network(gpd.GeoDataFrame(pd.concat([data, *extra], ignore_index=True)))


def write_inp(sub, tmp_path, name="model.inp", routing=None, outfalls=("OUT",), **kwargs):
    """Write sub to tmp_path/name from the shipped template, or from a copy of it with FLOW_ROUTING
    set to routing. The fixture's J1 -> J2 parallel pipes are a split that KINWAVE (the shipped
    template's routing) rejects, so tests that write them use routing="DYNWAVE"."""
    template = TEMPLATE
    if routing is not None:
        template = tmp_path / f"template_{routing}.inp"
        text = re.sub(
            r"(?m)^FLOW_ROUTING\s+\w+", f"FLOW_ROUTING         {routing}", TEMPLATE.read_text()
        )
        template.write_text(text)
    sub.to_inp(tmp_path / name, outfalls, template_path=template, **kwargs)
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
