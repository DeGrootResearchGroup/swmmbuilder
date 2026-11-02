"""Node and pipe coordinates (#5): [COORDINATES]/[VERTICES] in the .inp, add_coordinates and
generate_gmap."""

import json

import pytest
from conftest import DENDRITIC, FIXTURE, OVERFLOW, read_sections, subnetwork_with_edges, write_inp

from sewernetworkx import sewers


def fixture_geometry():
    """From the fixture GeoJSON: ({node: (x, y)}, {OBJECTID: [interior (x, y), ...]})."""
    nodes, vertices = {}, {}
    for f in json.loads(FIXTURE.read_text())["features"]:
        p, coords = f["properties"], [tuple(c) for c in f["geometry"]["coordinates"]]
        for node, xy in [(p["UpstreamID"], coords[0]), (p["DownstreamID"], coords[-1])]:
            assert nodes.setdefault(node, xy) == xy, node  # fixture pipes meet at their nodes
        vertices[str(p["OBJECTID"])] = coords[1:-1]
    return nodes, vertices


NODES, VERTICES = fixture_geometry()
EDGES = DENDRITIC + OVERFLOW


@pytest.fixture
def sections(network, tmp_path):
    return read_sections(write_inp(subnetwork_with_edges(network, EDGES), tmp_path))


def test_coordinates_section(sections):
    rows = {row[0]: (float(row[1]), float(row[2])) for row in sections["COORDINATES"]}
    assert set(rows) == {n for e in EDGES for n in e}
    for node, xy in rows.items():
        assert xy == pytest.approx(NODES[node], abs=1e-12), node


def test_vertices_section(sections):
    written = {}
    for row in sections.get("VERTICES", []):
        written.setdefault(row[0], []).append((float(row[1]), float(row[2])))
    conduits = [row[0] for row in sections["CONDUITS"]]
    expected = {name: VERTICES[name] for name in conduits if VERTICES[name]}
    assert expected, "the fixture should have multi-vertex pipes in this subnetwork"
    assert set(written) == set(expected)
    for name, xys in expected.items():
        flat = [v for xy in xys for v in xy]
        assert [v for xy in written[name] for v in xy] == pytest.approx(flat, abs=1e-12), name


def test_add_coordinates_uses_pipe_endpoints(network):
    sub = subnetwork_with_edges(network, EDGES)
    sub.scale = 10.0
    sub.add_coordinates()
    assert set(sub.coords) == {n for e in EDGES for n in e}
    # OV2 is only a downstream end, of the multi-vertex OV1 -> OV2 pipe: it must sit at the last
    # vertex, not the midpoint vertex.
    for node, (x, y) in sub.coords.items():
        assert (x, y) == pytest.approx((10 * NODES[node][0], 10 * NODES[node][1])), node


def test_generate_gmap_divides_by_scale(network, tmp_path, monkeypatch):
    calls = {"markers": {}, "lines": []}

    class FakePlotter:
        def __init__(self, lat, lng, zoom):
            calls["centre"] = (lat, lng)

        def marker(self, lat, lng, color, label):
            calls["markers"][label] = (lat, lng)

        def plot(self, lats, lngs, color, edge_width):
            calls["lines"].append((lats, lngs))

        def draw(self, file_name):
            pass

    monkeypatch.setattr(sewers.gmplot, "GoogleMapPlotter", FakePlotter)
    sub = subnetwork_with_edges(network, DENDRITIC)
    sub.scale = 10.0
    sub.add_coordinates()
    sub.generate_gmap(str(tmp_path / "map.html"), scale=10.0)
    for node, (lat, lng) in calls["markers"].items():
        assert (lng, lat) == pytest.approx(NODES[node]), node
    first = next(iter(sub.coords))
    assert calls["centre"][::-1] == pytest.approx(NODES[first])
    assert len(calls["lines"]) == len(DENDRITIC)
