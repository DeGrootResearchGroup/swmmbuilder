"""Render every static map plot with its defaults, for reviewing how the plots look.

By default the network is the Toronto test catchment (tests/data/toronto_catchment.geojson);
pass ``--pipes`` (London schema) and ``--outlet`` to render another one. The values are synthetic:
"flow" is the number of nodes upstream of each pipe, "concentration" decays with the distance to
the outlet, and the zones are the largest disjoint branches plus the rest.

    python examples/maps_gallery.py [OUT_DIR] [--pipes PATH --outlet NODE]

CI runs it and uploads OUT_DIR (default maps_gallery/) as the "maps-gallery" artifact.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

from swmmbuilder import LONDON_ON, TORONTO_ON, Network, maps

REPO = Path(__file__).resolve().parent.parent
TORONTO = REPO / "tests" / "data" / "toronto_catchment.geojson"
TORONTO_OUTLET = "MH5128918783"


def synthetic_values(sub, outlet, n_zones=4):
    """(flow, concentration, zones, shares) for the pipes and nodes of ``sub``."""
    G = sub.graph
    upstream = {n: len(nx.ancestors(G, n)) + 1 for n in G}
    flow = pd.Series({pipe: upstream[u] for u, _, pipe in G.edges(keys=True)}, dtype=float)
    distance = nx.single_source_dijkstra_path_length(
        G.reverse(copy=False), outlet, weight=lambda u, v, d: min(e["length"] for e in d.values())
    )
    scale = max(distance.values()) / 3 or 1.0
    concentration = pd.Series(
        {pipe: 1e6 * np.exp(-distance.get(u, 0) / scale) for u, _, pipe in G.edges(keys=True)}
    )

    zones, used = {}, set()
    low, high = 0.1 * len(G), 0.45 * len(G)
    for node in sorted((n for n in G if low <= upstream[n] <= high), key=lambda n: -upstream[n]):
        branch = nx.ancestors(G, node) | {node}
        if branch & used:
            continue
        zones[f"Z{len(zones) + 1}"] = sorted(branch)
        used |= branch
        if len(zones) == n_zones:
            break
    zones["rest"] = sorted(set(G) - used)
    zone_of = {n: z for z, nodes in zones.items() for n in nodes}
    shares = (
        pd.DataFrame(
            {
                pipe: pd.Series([zone_of[n] for n in nx.ancestors(G, u) | {u}]).value_counts()
                for u, _, pipe in G.edges(keys=True)
            }
        )
        .T.reindex(columns=list(zones))
        .fillna(0.0)
    )
    return flow, concentration, zones, shares


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_dir", nargs="?", default="maps_gallery", type=Path)
    parser.add_argument("--pipes", type=Path, help="pipes GIS file (London schema)")
    parser.add_argument("--outlet", help="outlet node of the catchment to draw")
    args = parser.parse_args()

    if args.pipes is None:
        net, outlet = Network(TORONTO, schema=TORONTO_ON), TORONTO_OUTLET
    else:
        net, outlet = Network(args.pipes, schema=LONDON_ON), args.outlet
    sub = net.catchment(outlet)
    flow, concentration, zones, shares = synthetic_values(sub, outlet)
    colours = maps.category_colors(zones)
    polygons = maps.node_group_polygons(sub, zones)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    def figure(name, draw):
        fig, ax = plt.subplots(figsize=(7, 6))
        draw(ax)
        ax.set_title(name.replace("_", " "), fontsize=10, loc="left")
        fig.savefig(args.out_dir / f"{name}.png", dpi=150, bbox_inches="tight")
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
        pipes = maps.plot_pipes(ax, sub, legend=False)
        ribbons = maps.plot_ribbons(
            ax, sub, shares, colours, context=False, title="Upstream zone", legend=False
        )
        site = maps.plot_markers(
            ax, sub, [outlet], labels={outlet: "Outlet"}, name="Outlet", legend=False
        )
        maps.map_axes(ax, north_arrow=True)
        maps.add_legends(ax, groups, pipes, ribbons, site)

    figure("combined", combined)
    print(f"Wrote {len(list(args.out_dir.glob('*.png')))} images to {args.out_dir}/")


if __name__ == "__main__":
    main()
