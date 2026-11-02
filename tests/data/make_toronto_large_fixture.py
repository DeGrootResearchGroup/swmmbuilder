"""Make tests/data/toronto_large_catchment.geojson.gz from the City of Toronto Open Data Portal.

The fixture is the sanitary catchment draining to maintenance hole OUTLET (about 2,700 pipes,
about 9 km across): a full-size network for the map gallery (examples/maps_gallery.py), next to
the 31-pipe toronto_catchment.geojson. It is taken from the Sewer Gravity Mains dataset
(https://open.toronto.ca/dataset/sewer-gravity-mains/): the gravity mains and the trunk sewers,
which are published as a separate layer and join the mains' catchments together. Trunk sewers
are renamed to the gravity mains' columns, with 'Sewer Gravity Trunk Sewer' = 'Yes'. The output
has the properties of the dataset's GeoJSON download, MultiLineString geometry in WGS 84, and is
gzipped. The dataset is updated daily, so a rerun may differ.

Contains information licensed under the Open Government Licence - Toronto.

Run from the repository root (downloads about 60 MB unless the GeoPackages are given):

    python tests/data/make_toronto_large_fixture.py [--mains MAINS.gpkg --trunks TRUNKS.gpkg]
"""

import argparse
import gzip
import tempfile
import urllib.request
from pathlib import Path

import geopandas as gpd
import networkx as nx
import pandas as pd

DATASET = (
    "https://ckan0.cf.opendata.inter.prod-toronto.ca/dataset/71a2da32-2a26-4343-bcaf-8ca63a696060"
)
MAINS_URL = (
    f"{DATASET}/resource/2dc28550-4822-45ec-afb8-848ce13d52d4/download/"
    "Sewer%20Gravity%20Main%20-%204326.gpkg"
)
TRUNKS_URL = (
    f"{DATASET}/resource/804aef88-c695-49b8-ad8f-1621ce8e6727/download/"
    "Trunk%20Sewer%20-%204326.gpkg"
)
OUTLET = "MH3948903491"
OUT = Path(__file__).with_name("toronto_large_catchment.geojson.gz")

PIPE_ID = "Sewer Gravity Asset Identification"
UPSTREAM = "Sewer Gravity Upstream Maintenance Hole"
DOWNSTREAM = "Sewer Gravity Downstream Maintenance Hole"
FLOW_TYPE = "Sewer Gravity Flow Type"

# Trunk Sewer column -> Sewer Gravity Main column, where the names differ beyond the prefix.
TRUNK_COLUMNS = {
    "Trunk Sewer Water Type": FLOW_TYPE,
    "Trunk Sewer Maintained By": "Sewer Gravity Managed By",
    "Trunk Sewer Name": "Sewer Gravity Trunk Name",
    "Trunk Sewer": "Sewer Gravity Trunk Sewer",
}


def fetch(url, directory):
    path = Path(directory) / url.rsplit("/", 1)[-1].replace("%20", " ")
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, path)
    return path


def trunks_as_mains(trunks, columns):
    """The trunk sewers in the gravity mains' columns."""
    renamed = {
        c: TRUNK_COLUMNS.get(c, c.replace("Trunk Sewer ", "Sewer Gravity "))
        for c in trunks.columns
        if c.startswith("Trunk Sewer")
    }
    trunks = trunks.rename(columns=renamed).drop(columns="_id")
    trunks["Sewer Gravity Trunk Sewer"] = "Yes"
    return trunks.reindex(columns=[c for c in columns if c != "_id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mains", type=Path, help="Sewer Gravity Main - 4326.gpkg")
    parser.add_argument("--trunks", type=Path, help="Trunk Sewer - 4326.gpkg")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as directory:
        mains = gpd.read_file(args.mains or fetch(MAINS_URL, directory))
        trunks = gpd.read_file(args.trunks or fetch(TRUNKS_URL, directory))
    trunks = trunks_as_mains(trunks, mains.columns)
    trunks = trunks[~trunks[PIPE_ID].isin(mains[PIPE_ID])]
    pipes = pd.concat([mains, trunks], ignore_index=True)
    pipes = pipes[pipes[FLOW_TYPE] == "SAN"]

    graph = nx.DiGraph()
    graph.add_edges_from(zip(pipes[UPSTREAM], pipes[DOWNSTREAM], strict=True))
    catchment = nx.ancestors(graph, OUTLET) | {OUTLET}
    pipes = pipes[pipes[UPSTREAM].isin(catchment) & pipes[DOWNSTREAM].isin(catchment)]
    pipes = gpd.GeoDataFrame(pipes.sort_values(PIPE_ID), crs=mains.crs).reset_index(drop=True)
    for column in pipes.columns:
        if pd.api.types.is_datetime64_any_dtype(pipes[column]):
            pipes[column] = pipes[column].dt.strftime("%Y-%m-%d")

    text = pipes.to_json(na="null", drop_id=True)
    with gzip.open(OUT, "wt", encoding="utf-8") as f:
        f.write(text)
    print(f"Wrote {len(pipes)} pipes to {OUT}")


if __name__ == "__main__":
    main()
