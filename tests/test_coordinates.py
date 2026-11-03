"""Node and pipe coordinates (#5): [COORDINATES]/[VERTICES] in the .inp, node_positions and
maps.generate_gmap."""

import dataclasses
import json

import geopandas as gpd
import pytest
from conftest import DENDRITIC, FIXTURE, OVERFLOW, read_sections, subnetwork, write_inp
from pyproj import Transformer

from swmmbuilder import LONDON_ON, Network, maps


def fixture_geometry(project=True):
    """From the fixture GeoJSON: ({node: (x, y)}, {OBJECTID: [interior (x, y), ...]}), projected
    to the London schema's CRS (or in longitude/latitude, as in the file)."""
    to_utm = Transformer.from_crs(4326, LONDON_ON.crs, always_xy=True)
    nodes, vertices = {}, {}
    for f in json.loads(FIXTURE.read_text())["features"]:
        p, coords = f["properties"], [tuple(c) for c in f["geometry"]["coordinates"]]
        if project:
            coords = [to_utm.transform(x, y) for x, y in coords]
        for node, xy in [(p["UpstreamID"], coords[0]), (p["DownstreamID"], coords[-1])]:
            assert nodes.setdefault(node, xy) == xy, node  # fixture pipes meet at their nodes
        vertices[str(p["OBJECTID"])] = coords[1:-1]
    return nodes, vertices


NODES, VERTICES = fixture_geometry()
LONLAT, _ = fixture_geometry(project=False)
EDGES = DENDRITIC + OVERFLOW


@pytest.fixture
def sections(network, tmp_path):
    sub = subnetwork(network, EDGES)
    return read_sections(write_inp(sub, tmp_path, routing="DYNWAVE", outfalls=("OUT", "OV2")))


def test_coordinates_section(sections):
    rows = {row[0]: (float(row[1]), float(row[2])) for row in sections["COORDINATES"]}
    assert set(rows) == {n for e in EDGES for n in e}
    for node, xy in rows.items():
        assert xy == pytest.approx(NODES[node], abs=1e-3), node


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
        assert [v for xy in written[name] for v in xy] == pytest.approx(flat, abs=1e-3), name


def test_node_positions_use_pipe_endpoints(network):
    positions = subnetwork(network, EDGES).node_positions()
    assert set(positions) == {n for e in EDGES for n in e}
    for node, xy in positions.items():
        assert xy == pytest.approx(NODES[node]), node


def test_generate_gmap(network, tmp_path, monkeypatch):
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

    monkeypatch.setattr("gmplot.GoogleMapPlotter", FakePlotter)
    sub = subnetwork(network, DENDRITIC)
    with pytest.warns(DeprecationWarning, match="to_leaflet"):
        maps.generate_gmap(sub, str(tmp_path / "map.html"))
    for node, (lat, lng) in calls["markers"].items():
        assert (lng, lat) == pytest.approx(LONLAT[node]), node
    first = next(iter(sub.node_positions()))
    assert calls["centre"][::-1] == pytest.approx(LONLAT[first])
    assert len(calls["lines"]) == len(DENDRITIC)  # parallel J1 -> J2 pipes drawn once


def test_map_positions_without_crs_are_lonlat():
    data = gpd.read_file(FIXTURE).set_crs(None, allow_override=True)
    sub = subnetwork(Network(data, schema=dataclasses.replace(LONDON_ON, crs=None)), DENDRITIC)
    assert maps.lonlat_positions(sub) == sub.node_positions()
