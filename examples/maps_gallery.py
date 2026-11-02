"""Render every static map plot with its defaults, for reviewing how the plots look.

By default it renders two Toronto test catchments, so that defaults can be checked at both
scales: the 31-pipe tests/data/toronto_catchment.geojson ("toronto_small", about 560 m across)
and the 2,692-pipe tests/data/toronto_large_catchment.geojson.gz ("toronto_large", about 9 km
across). Pass ``--pipes`` (London schema) and ``--outlet`` to render another network instead.
The values are synthetic: "flow" is the number of nodes upstream of each pipe, "concentration"
decays with the distance to the outlet, and the zones are the largest disjoint branches plus the
rest. Images are named after the network, e.g. toronto_large_ribbons.png.

    python examples/maps_gallery.py [OUT_DIR] [--pipes PATH --outlet NODE]

CI runs it and uploads OUT_DIR (default maps_gallery/) as the "maps-gallery" artifact.
"""

import argparse
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from swmmbuilder import LONDON_ON, TORONTO_ON, Network, maps

DATA = Path(__file__).resolve().parent.parent / "tests" / "data"
# (name, pipes file, outlet)
TORONTO = [
    ("toronto_small", DATA / "toronto_catchment.geojson", "MH5128918783"),
    ("toronto_large", DATA / "toronto_large_catchment.geojson.gz", "MH3948903491"),
]


def read_pipes(path):
    """The pipes in ``path``, which may be gzipped."""
    path = Path(path)
    return gpd.read_file(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path)


def synthetic_values(sub, outlet, n_zones=4):
    """(flow, concentration, zones, shares) for the pipes and nodes of ``sub``."""
    G = sub.graph
    ancestors = {n: nx.ancestors(G, n) | {n} for n in G}
    flow = pd.Series({pipe: len(ancestors[u]) for u, _, pipe in G.edges(keys=True)}, dtype=float)
    distance = nx.single_source_dijkstra_path_length(
        G.reverse(copy=False), outlet, weight=lambda u, v, d: min(e["length"] for e in d.values())
    )
    scale = max(distance.values()) / 3 or 1.0
    concentration = pd.Series(
        {pipe: 1e6 * np.exp(-distance.get(u, 0) / scale) for u, _, pipe in G.edges(keys=True)}
    )

    zones, used = {}, set()
    low, high = 0.1 * len(G), 0.45 * len(G)
    size = {n: len(a) for n, a in ancestors.items()}
    for node in sorted((n for n in G if low <= size[n] <= high), key=lambda n: -size[n]):
        if ancestors[node] & used:
            continue
        zones[f"Z{len(zones) + 1}"] = sorted(ancestors[node])
        used |= ancestors[node]
        if len(zones) == n_zones:
            break
    zones["rest"] = sorted(set(G) - used)
    zone_of = {n: z for z, nodes in zones.items() for n in nodes}
    counts = {n: pd.Series([zone_of[a] for a in ancestors[n]]).value_counts() for n in G}
    shares = (
        pd.DataFrame({pipe: counts[u] for u, _, pipe in G.edges(keys=True)})
        .T.reindex(columns=list(zones))
        .fillna(0.0)
    )
    return flow, concentration, zones, shares


def render(name, net, outlet, out_dir):
    """Write ``out_dir/{name}_{plot}.png`` for every plot of the catchment of ``outlet``."""
    sub = net.catchment(outlet)
    flow, concentration, zones, shares = synthetic_values(sub, outlet)
    colours = maps.category_colors(zones)
    polygons = maps.node_group_polygons(sub, zones)

    def figure(plot, draw):
        fig, ax = plt.subplots(figsize=(7, 6))
        draw(ax)
        ax.set_title(f"{name}: {plot.replace('_', ' ')}", fontsize=10, loc="left")
        fig.savefig(out_dir / f"{name}_{plot}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    figure("pipes", lambda ax: maps.plot_pipes(ax, sub))
    figure(
        "pipes_by_value",
        lambda ax: maps.plot_pipes(
            ax, sub, concentration, norm="log", width=flow, label="Concentration (copies/L)"
        ),
    )
    figure("ribbons", lambda ax: maps.plot_ribbons(ax, sub, shares, colours, title="Upstream zone"))
    figure("groups", lambda ax: maps.plot_groups(ax, polygons, colors=colours, label="Zone"))
    figure(
        "groups_by_value",
        lambda ax: maps.plot_groups(
            ax,
            polygons,
            {z: len(nodes) for z, nodes in zones.items()},
            hatch_groups={"Z1 + Z2": ["Z1", "Z2"]},
            label="Nodes",
        ),
    )
    figure(
        "markers",
        lambda ax: (
            maps.plot_pipes(ax, sub),
            maps.plot_markers(ax, sub, [outlet], labels={outlet: "Outlet"}, name="Outlet"),
        ),
    )
    figure(
        "concave_hulls",
        lambda ax: maps.plot_groups(
            ax, maps.node_group_polygons(sub, zones, "concave_hull"), colors=colours
        ),
    )

    def combined(ax):
        groups = maps.plot_groups(ax, polygons, colors=colours, label="Zone", legend=False)
        maps.plot_pipes(ax, sub, width=flow, context=True, legend=False)
        ribbons = maps.plot_ribbons(
            ax, sub, shares, colours, context=False, title="Upstream zone", legend=False
        )
        site = maps.plot_markers(
            ax, sub, [outlet], labels={outlet: "Outlet"}, name="Outlet", legend=False
        )
        maps.map_axes(ax, north_arrow=True)
        maps.add_legends(ax, groups, ribbons, site)

    figure("combined", combined)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_dir", nargs="?", default="maps_gallery", type=Path)
    parser.add_argument("--pipes", type=Path, help="pipes GIS file (London schema)")
    parser.add_argument("--outlet", help="outlet node of the catchment to draw")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.pipes is None:
        for name, path, outlet in TORONTO:
            render(name, Network(read_pipes(path), schema=TORONTO_ON), outlet, args.out_dir)
    else:
        if args.outlet is None:
            parser.error("--pipes needs --outlet")
        net = Network(read_pipes(args.pipes), schema=LONDON_ON)
        render(args.pipes.name.split(".")[0], net, args.outlet, args.out_dir)
    print(f"Wrote {len(list(args.out_dir.glob('*.png')))} images to {args.out_dir}/")


if __name__ == "__main__":
    main()
