"""Fetch tests/data/toronto_catchment.geojson from the City of Toronto Open Data Portal.

The fixture is the sanitary catchment draining to maintenance hole OUTLET, plus the pipe leaving
OUTLET (so that the catchment can discharge through an outfall with one inlet), taken from the Sewer
Gravity Mains dataset (https://open.toronto.ca/dataset/sewer-gravity-mains/) through the CKAN
datastore API, and written in the format of the dataset's GeoJSON download (same properties,
MultiLineString geometry, WGS 84). The dataset is updated daily, so a rerun may differ.

Contains information licensed under the Open Government Licence - Toronto.

Run from the repository root:  python tests/data/make_toronto_fixture.py
"""

import json
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://ckan0.cf.opendata.inter.prod-toronto.ca/api/3/action/datastore_search"
RESOURCE = "b774654a-d58d-47ec-95cb-2030c3621b1b"  # Sewer Gravity Main
OUTLET = "MH5130118747"
DOWNSTREAM = "Sewer Gravity Downstream Maintenance Hole"
UPSTREAM = "Sewer Gravity Upstream Maintenance Hole"
OUT = Path(__file__).with_name("toronto_catchment.geojson")


def sanitary_pipes(**filters):
    """Sanitary pipes matching ``filters`` ({field: value})."""
    query = {
        "resource_id": RESOURCE,
        "filters": json.dumps({**filters, "Sewer Gravity Flow Type": "SAN"}),
        "limit": 100,
    }
    with urllib.request.urlopen(f"{API}?{urllib.parse.urlencode(query)}") as response:
        return json.load(response)["result"]["records"]


def main():
    features, queue, seen = [], [OUTLET], {OUTLET}
    records = sanitary_pipes(**{UPSTREAM: OUTLET})
    while queue:
        records += sanitary_pipes(**{DOWNSTREAM: queue.pop()})
        while records:
            record = records.pop()
            geometry = json.loads(record.pop("geometry"))
            if geometry["type"] == "LineString":
                geometry = {"type": "MultiLineString", "coordinates": [geometry["coordinates"]]}
            features.append({"type": "Feature", "properties": record, "geometry": geometry})
            if record[UPSTREAM] not in seen:
                seen.add(record[UPSTREAM])
                queue.append(record[UPSTREAM])
    features.sort(key=lambda f: f["properties"]["_id"])
    collection = {
        "type": "FeatureCollection",
        "name": "Sewer Gravity Main - 4326",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
        "features": features,
    }
    OUT.write_text(json.dumps(collection, indent=1) + "\n")
    print(f"Wrote {len(features)} pipes to {OUT}")


if __name__ == "__main__":
    main()
