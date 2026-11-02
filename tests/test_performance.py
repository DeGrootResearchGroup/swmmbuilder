"""The .inp writer must scale with the size of the subnetwork, not edges x all pipes (#7)."""

import json
import time

from conftest import quiet, read_sections, subnetwork_with_edges, write_inp

from swmmbuilder import Network

N_PIPES = 3000
# Generous: the writer takes ~0.2 s here; the old per-edge table scan took ~80 s.
TIME_LIMIT_S = 10.0


def _node(i):
    return f"N{i}"


def _depth(i):
    return (i + 1).bit_length() - 1


def binary_tree_geojson(path, n_pipes):
    """A dendritic network: node i drains to node (i - 1) // 2; N0 is the outlet."""
    features = []
    for i in range(1, n_pipes + 1):
        parent = (i - 1) // 2
        props = {
            "OBJECTID": 10_000 + i,
            "FlowType": "SAN",
            "UpstreamID": _node(i),
            "DownstreamID": _node(parent),
            "Length": 50.0,
            "PipeShape": "R" if i % 2 else "REC",
            "Height": 300,
            "Width": 450,
            "Diameter": 300,
            "Material": "PVC",
            "UpstreamInvert": 100.0 + 0.5 * _depth(i),
            "DownstreamInvert": 100.1 + 0.5 * _depth(parent),
            "Abandoned": 0,
        }
        coords = [[-81.3 + 1e-4 * i, 42.98 + 1e-4 * _depth(i)]]
        coords.append([-81.3 + 1e-4 * parent, 42.98 + 1e-4 * _depth(parent)])
        features.append(
            {
                "type": "Feature",
                "properties": props,
                "geometry": {"type": "LineString", "coordinates": coords},
            }
        )
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))
    return path


def test_writer_scales_to_thousands_of_pipes(tmp_path):
    network = quiet(Network, str(binary_tree_geojson(tmp_path / "tree.geojson", N_PIPES)))
    edges = [(_node(i), _node((i - 1) // 2)) for i in range(1, N_PIPES + 1)]
    sub = subnetwork_with_edges(network, edges, outfalls=(_node(0),), start=_node(0))

    start = time.perf_counter()
    inp = write_inp(sub, tmp_path)
    elapsed = time.perf_counter() - start

    sections = read_sections(inp)
    assert len(sections["CONDUITS"]) == N_PIPES
    assert len(sections["JUNCTIONS"]) + len(sections["OUTFALLS"]) == N_PIPES + 1
    assert elapsed < TIME_LIMIT_S, f"writing {N_PIPES} pipes took {elapsed:.1f} s"
