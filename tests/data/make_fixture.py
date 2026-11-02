"""Generate fixture_network.geojson: a small sewer network with deliberate edge cases.

Run from anywhere: ``python tests/data/make_fixture.py``. The schema matches the City of London
(Ontario) sewer GeoJSON that swmmbuilder was written for. Inverts in m, sizes in mm.

Layout (flow direction):

    A1 -> A2 -> J1 -> J2 -> OUT         main line
    B1 -> B2 -> J1                      REC section on B2-J1
    C1 -> J1                            drop structure into J1
    J1 -> J2                            EGG section, plus a parallel round pipe
    J2 -> OV1 -> OV2                    high-level overflow: flow split at J2, dead end at OV2
    D1 -> D2 -> J2                      D1-D2 is adverse in the GIS inverts
    E1 -> J2                            blank PipeShape
    X1 -> A2                            Abandoned = 1
    'Abandoned' -> A1                   abandoned by ID (the only case the library drops today)
    S1 -> S2                            storm pipe (FlowType STM)
    F1..F4, G1..G4 -> A1, B1            extra branches so that branch exploration has choices

The J2 -> OV1 overflow is listed *before* J2 -> OUT, so the "first pipe leaving the node" invert
rule picks the overflow's higher invert for J2 (the OK0451 case in the London data).
"""

import json
from pathlib import Path

X0, Y0 = -81.30, 42.98
POS = {
    "A1": (0, 4),
    "A2": (1, 4),
    "J1": (2, 3),
    "J2": (3, 2),
    "OUT": (4, 1),
    "B1": (0, 2),
    "B2": (1, 2.5),
    "C1": (2, 5),
    "OV1": (4, 3),
    "OV2": (5, 3.5),
    "D1": (2, 0),
    "D2": (2.5, 1),
    "E1": (3, 0.5),
    "X1": (0.5, 5),
    "Abandoned": (-1, 4.5),
    "S1": (5, 5),
    "S2": (6, 5),
    "F1": (-1, 6),
    "F2": (-1, 5.5),
    "F3": (-2, 6),
    "F4": (-2, 5),
    "G1": (-1, 1),
    "G2": (-1, 2.5),
    "G3": (-2, 1),
    "G4": (-2, 3),
}

# (up, down, shape, height, width, diameter, inv_up, inv_down, drop, length, material, flow, abandoned)
PIPES = [
    ("A1", "A2", "R", 200, 200, 200, 105.00, 104.50, 0.0, 80.0, "PVC", "SAN", 0),
    ("A2", "J1", "R", 250, 250, 250, 104.40, 104.00, 0.0, 60.0, "PVC", "SAN", 0),
    ("B1", "B2", "Round", 200, 200, 200, 104.80, 104.30, 0.0, 70.0, "PVC", "SAN", 0),
    ("B2", "J1", "REC", 300, 450, 450, 104.20, 103.95, 0.0, 50.0, "CONC", "SAN", 0),
    ("C1", "J1", "R", 200, 200, 200, 106.00, 105.50, 103.98, 40.0, "PVC", "SAN", 0),
    ("J2", "OV1", "R", 300, 300, 300, 104.60, 104.55, 0.0, 15.0, "CONC", "SAN", 0),
    ("J1", "J2", "EGG", 600, 400, 490, 103.90, 103.60, 0.0, 90.0, "BRCK", "SAN", 0),
    ("J1", "J2", "R", 300, 300, 300, 103.90, 103.60, 0.0, 90.0, "CONC", "SAN", 0),
    ("J2", "OUT", "R", 600, 600, 600, 103.50, 103.30, 0.0, 40.0, "CONC", "SAN", 0),
    ("OV1", "OV2", "R", 300, 300, 300, 104.50, 104.40, 0.0, 10.0, "CONC", "SAN", 0),
    ("D1", "D2", "R", 200, 200, 200, 103.80, 104.10, 0.0, 30.0, "PVC", "SAN", 0),
    ("D2", "J2", "R", 200, 200, 200, 104.00, 103.70, 0.0, 30.0, "PVC", "SAN", 0),
    ("E1", "J2", "", None, 0, 200, 104.30, 103.80, 0.0, 25.0, "", "SAN", 0),
    ("X1", "A2", "R", 200, 200, 200, 105.20, 104.60, 0.0, 30.0, "AC", "SAN", 1),
    ("Abandoned", "A1", "R", 200, 200, 200, 105.50, 105.10, 0.0, 30.0, "AC", "SAN", 0),
    ("S1", "S2", "R", 300, 300, 300, 104.00, 103.50, 0.0, 50.0, "CONC", "STM", 0),
    ("F1", "A1", "R", 200, 200, 200, 105.60, 105.10, 0.0, 40.0, "PVC", "SAN", 0),
    ("F2", "A1", "R", 200, 200, 200, 105.50, 105.10, 0.0, 40.0, "PVC", "SAN", 0),
    ("F3", "F1", "R", 200, 200, 200, 106.00, 105.70, 0.0, 40.0, "PVC", "SAN", 0),
    ("F4", "F2", "R", 200, 200, 200, 106.00, 105.60, 0.0, 40.0, "PVC", "SAN", 0),
    ("G1", "B1", "R", 200, 200, 200, 105.30, 104.90, 0.0, 40.0, "PVC", "SAN", 0),
    ("G2", "B1", "R", 200, 200, 200, 105.20, 104.90, 0.0, 40.0, "PVC", "SAN", 0),
    ("G3", "G1", "R", 200, 200, 200, 105.80, 105.40, 0.0, 40.0, "PVC", "SAN", 0),
    ("G4", "G2", "R", 200, 200, 200, 105.70, 105.30, 0.0, 40.0, "PVC", "SAN", 0),
]


def feature(i, p):
    up, down, shape, h, w, d, iu, idn, drop, length, mat, flow, aband = p
    (x1, y1), (x2, y2) = POS[up], POS[down]
    coords = [[X0 + 1e-3 * x1, Y0 + 1e-3 * y1]]
    if i % 3 == 0:  # some multi-vertex pipes
        coords.append([X0 + 1e-3 * (x1 + x2) / 2, Y0 + 1e-3 * (y1 + y2) / 2 + 2e-4])
    coords.append([X0 + 1e-3 * x2, Y0 + 1e-3 * y2])
    props = {
        "OBJECTID": 1000 + i,
        "GIS_ID": 2000 + i,
        "GIS_FeatureKey": str(2000 + i),
        "FlowType": flow,
        "StreetName": "Test St.",
        "Location": "Road",
        "UpstreamID": up,
        "DownstreamID": down,
        "Length": length,
        "PipeShape": shape,
        "Height": h,
        "Width": w,
        "Material": mat,
        "MaterialClass": "",
        "UpstreamInvert": iu,
        "DownstreamInvert": idn,
        "DropStructureInvert": drop,
        "Slope": round(100 * (iu - idn) / length, 2),
        "AverageDepth": 2.5,
        "MaximumDepth": 3.0,
        "PDC_Count": 0.0,
        "ConstructionDate": "2000-01-01T00:00:00Z",
        "Ownership": "CITY",
        "Assumed": 1,
        "Combined": 0,
        "Relief": 0,
        "Trunk": int(d >= 450),
        "LastEditDate": "2022-01-01T00:00:00Z",
        "Abandoned": aband,
        "Grouted": 0,
        "Diameter": d,
        "CCTVInspectionDate": None,
        "ShapeSTLength": length,
    }
    return {
        "type": "Feature",
        "properties": props,
        "geometry": {"type": "LineString", "coordinates": coords},
    }


if __name__ == "__main__":
    out = Path(__file__).with_name("fixture_network.geojson")
    fc = {"type": "FeatureCollection", "features": [feature(i, p) for i, p in enumerate(PIPES)]}
    out.write_text(json.dumps(fc, indent=1))
    print(f"wrote {out} ({len(PIPES)} pipes)")
