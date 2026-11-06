"""A loaded catchment model, end to end: from the GIS pipes to a SWMM .inp file with its
dry-weather flows.

The steps:

1. read the pipes into a ``Network`` and take the catchment of the outlet;
2. place the population of census-like polygons at the nodes (``allocate_population``);
3. place the flow of commercial and industrial zones at the nodes (``allocate_land_use``);
4. spread groundwater infiltration over the pipes (``distribute_flow``);
5. combine the three into the model's options and a per-node table (``NodeLoading``);
6. write the .inp file (``Subnetwork.to_inp``).

The pipes are the 52 of examples/geojson_test_data/test_data_small.geojson (London, Ontario's
schema). The polygons are synthetic: a grid of blocks over the network's extent with made-up
populations, and three zoning boxes. The rates and patterns are illustrative, not design values.

    python examples/loaded_catchment.py [OUT_DIR]

It writes OUT_DIR/loaded_catchment.inp and OUT_DIR/node_loading.csv (default OUT_DIR:
loaded_catchment/).
"""

import argparse
from pathlib import Path

import geopandas as gpd
from shapely.geometry import box

from swmmbuilder import (
    LONDON_ON,
    Network,
    NodeLoading,
    allocate_land_use,
    allocate_population,
    distribute_flow,
)

PIPES = Path(__file__).resolve().parent / "geojson_test_data" / "test_data_small.geojson"
OUTLET = "WT1420"

# Hourly multipliers from 00:00, each with a mean of 1.
DIURNAL = [0.5, 0.3, 0.2, 0.2, 0.3, 0.5, 1.0, 1.9, 1.9, 1.4, 1.2, 1.1,
           1.1, 1.0, 1.0, 1.0, 1.1, 1.3, 1.5, 1.5, 1.3, 1.1, 0.9, 0.7]  # fmt: skip
BUSINESS_HOURS = [0.4] * 7 + [1.6] * 11 + [0.6] * 6  # 07:00 to 18:00
ROUND_THE_CLOCK = [1.0] * 24

PER_CAPITA = 230  # L/person/day
LAND_USE_RATES = {"COMMERCIAL": 5_000, "INDUSTRIAL": 8_000}  # L/ha/day
LAND_USE_PATTERNS = {"COMMERCIAL": BUSINESS_HOURS, "INDUSTRIAL": ROUND_THE_CLOCK}
INFILTRATION = 0.0005  # L/s per metre of pipe
DAYS = 3  # length of the run, and of the land-use time series


def blocks(net, nx=6, ny=1, people=300):
    """A grid of nx by ny census-like blocks covering the network, ``people`` each."""
    x0, y0, x1, y1 = net.full_data.total_bounds
    x0, y0, x1, y1 = x0 - 50, y0 - 50, x1 + 50, y1 + 50
    dx, dy = (x1 - x0) / nx, (y1 - y0) / ny
    cells = [
        box(x0 + i * dx, y0 + j * dy, x0 + (i + 1) * dx, y0 + (j + 1) * dy)
        for j in range(ny)
        for i in range(nx)
    ]
    return gpd.GeoDataFrame(
        {"population": [people + 50 * k for k in range(len(cells))]},
        geometry=cells,
        crs=net.full_data.crs,
        index=[f"B{k + 1}" for k in range(len(cells))],
    )


def zones(net):
    """Three zoning boxes: commercial at the west end, industrial in the middle, and residential
    (no rate, so no flow) at the east end."""
    x0, y0, x1, y1 = net.full_data.total_bounds
    width = (x1 - x0) / 5
    return gpd.GeoDataFrame(
        {"land_use": ["COMMERCIAL", "INDUSTRIAL", "RESIDENTIAL"]},
        geometry=[
            box(x0 - 50, y0 - 50, x0 + width, y1 + 50),
            box(x0 + 2 * width, y0 - 50, x0 + 3 * width, y1 + 50),
            box(x0 + 4 * width, y0 - 50, x1 + 50, y1 + 50),
        ],
        crs=net.full_data.crs,
        index=["Z1", "Z2", "Z3"],
    )


def main(out_dir="loaded_catchment"):
    """Build the loaded model and write OUT_DIR/loaded_catchment.inp and node_loading.csv.

    Returns
    -------
    (SwmmModel, NodeLoading)
        The model written and its load.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. The network, and the catchment of its outlet.
    net = Network(PIPES, schema=LONDON_ON)
    sub = net.catchment(OUTLET)

    # 2. Population at the nodes. Allocate over the whole network, then keep the catchment's.
    people = allocate_population(net, blocks(net), population="population")

    # 3. Land-use flows (L/s) at the nodes, one column per land use.
    land = allocate_land_use(net, zones(net), land_use="land_use", rates=LAND_USE_RATES)

    # 4. Infiltration (L/s), in proportion to pipe length.
    infiltration = distribute_flow(sub, rate=INFILTRATION)

    # 5. The three parts as [DWF], [INFLOWS] series, and [INFLOWS] baselines.
    loading = NodeLoading(
        population=people.within(sub),
        per_capita=PER_CAPITA,
        unit="L/d",
        pattern=DIURNAL,
        patterned=land.within(sub),
        patterns=LAND_USE_PATTERNS,
        days=DAYS,
        prefix="LU_",
        steady=infiltration,
    )
    for name, value in loading.summary().items():
        shown = f"{value:,}" if isinstance(value, int) else f"{value:,.2f}"
        print(f"{name:>30}  {shown}")

    # 6. The model.
    model = sub.to_inp(out / "loaded_catchment.inp", outfalls=[OUTLET], **loading.options())
    loading.table().to_csv(out / "node_loading.csv", float_format="%.6f")
    print(f"wrote {out / 'loaded_catchment.inp'} and {out / 'node_loading.csv'}")
    return model, loading


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_dir", nargs="?", default="loaded_catchment")
    main(parser.parse_args().out_dir)
