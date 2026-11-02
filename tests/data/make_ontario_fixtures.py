"""Fetch tests/data/<city>_catchment.geojson from the open data portals of six Ontario cities.

Each fixture is the sanitary catchment draining to a maintenance hole (the city's ``outlet``), plus
the pipe leaving it (so that the catchment can discharge through an outfall with one inlet), taken
from the city's ArcGIS REST layer and written as the layer's GeoJSON query output (all fields, WGS
84), less the fields that can name a person (editors, operators). The layers are updated by the
cities, so a rerun may differ.

Contains information licensed under:
- Barrie: the City of Barrie open data terms;
- Kitchener: the Open Government Licence - The Corporation of the City of Kitchener;
- Peterborough: the City of Peterborough Open Data Terms of Use;
- Sarnia: the City of Sarnia Open Data Terms of Use;
- Waterloo: the City of Waterloo Open Data User Licence;
- Windsor: the City of Windsor Open Data Terms of Use.

Run from the repository root:  python tests/data/make_ontario_fixtures.py [city ...]
"""

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

CITIES = {
    "barrie": {
        "layer": "https://gispublic.barrie.ca/arcgis/rest/services/Open_Data/"
        "SanitaryInfrastructure/MapServer/2",
        "upstream": "FROM_ID",
        "downstream": "TO_ID",
        "outlet": "120264434",
    },
    "kitchener": {
        "layer": "https://services1.arcgis.com/qAo1OsXi67t7XgmS/arcgis/rest/services/"
        "Sanitary_Pipes/FeatureServer/0",
        "upstream": "UP_SANMANHOLEID",
        "downstream": "DN_SANMANHOLEID",
        "outlet": 2151734,
    },
    "peterborough": {
        "layer": "https://services1.arcgis.com/JOTjMnANji9YTrZ2/arcgis/rest/services/"
        "San_Sewer_Pipe/FeatureServer/0",
        "upstream": "UP_MH_ID",
        "downstream": "DOWN_MH_ID",
        "outlet": "83241",
    },
    "sarnia": {
        "layer": "https://services1.arcgis.com/ICybsLmBXrZCZV3x/arcgis/rest/services/"
        "Sanitary_Sewers_Open_Data_AGOL/FeatureServer/0",
        "upstream": "MH_Upstream",
        "downstream": "MH_Downstream",
        "outlet": "MH7923",
    },
    "waterloo": {
        "layer": "https://services.arcgis.com/ZpeBVw5o1kjit7LT/arcgis/rest/services/"
        "Sanitary_Sewer_Gravity_Mains/FeatureServer/0",
        "upstream": "FROMMH",
        "downstream": "TOMH",
        "outlet": "W2-11D",
    },
    "windsor": {
        "layer": "https://mappmycity.ca/arcgis/rest/services/OpenDataServices/Infrastructure/"
        "MapServer/14",
        "upstream": "UPSTREAM_MH",
        "downstream": "DWNSTREAM_MH",
        "outlet": "1S4102",
    },
}
# Fields that can name a person; they are left out of the fixtures.
PERSONAL = {"created_user", "last_edited_user", "Clean_Operator"}


def sql(value):
    return f"'{value}'" if isinstance(value, str) else str(value)


def pipes(layer, field, values):
    """GeoJSON features of the pipes of ``layer`` whose ``field`` is one of ``values``."""
    query = {
        "where": f"{field} IN ({', '.join(sql(v) for v in values)})",
        "outFields": "*",
        "outSR": 4326,
        "f": "geojson",
    }
    with urllib.request.urlopen(f"{layer}/query?{urllib.parse.urlencode(query)}") as response:
        return json.load(response)["features"]


def catchment(city, outlet=None):
    """Features of the pipes draining to ``outlet`` (default: the city's), plus those leaving it."""
    spec = CITIES[city]
    outlet = spec["outlet"] if outlet is None else outlet
    features = pipes(spec["layer"], spec["upstream"], [outlet])
    queue, seen = [outlet], {outlet}
    while queue:
        batch, queue = queue[:50], queue[50:]
        for feature in pipes(spec["layer"], spec["downstream"], batch):
            features.append(feature)
            node = feature["properties"][spec["upstream"]]
            if node is not None and node not in seen:
                seen.add(node)
                queue.append(node)
    return sorted(features, key=lambda f: f["properties"]["OBJECTID"])


def main(cities):
    for city in cities:
        features = catchment(city)
        for feature in features:
            feature.pop("id", None)
            for field in PERSONAL:
                feature["properties"].pop(field, None)
        out = Path(__file__).with_name(f"{city}_catchment.geojson")
        collection = {"type": "FeatureCollection", "features": features}
        out.write_text(json.dumps(collection, indent=1) + "\n")
        print(f"Wrote {len(features)} pipes to {out}")


if __name__ == "__main__":
    main(sys.argv[1:] or CITIES)
