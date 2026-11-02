"""Generate fixture_network.geojson: a small sewer network with deliberate edge cases, and the
fixture_nodes_*.geojson node-layer network (see the second half of this file).

Run from anywhere: ``python tests/data/make_fixture.py``. The schema matches the City of London
(Ontario) sewer GeoJSON that SewerNetworkx was written for. Inverts in m, sizes in mm.

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


# --- Node layers (#19) -----------------------------------------------------------------------
#
# A second network in London's newer export schema (ConstructedStatus, StructureType and the
# inventory types of both end nodes; no Abandoned column), with its maintenance-hole, pumping
# station, outfall and other-node layers. Coordinates are in UTM 17N (EPSG:32617, m), so that
# the coincident nodes can be placed a fraction of a metre apart. Layout (flow direction):
#
#     M1 -> M2 -> PS1 => M3 -> M4 -> OF1    PS1 is a pumping station, PS1 => M3 a force main
#     M3 -> OF2                             overflow pipe (StructureType OVF) to an outfall
#     M8 -> DUP1   DUP2 -> M3               DUP1 (maintenance hole) and DUP2 (DATE point) are
#                                           0.2 m apart: one structure under two IDs
#     X1 -> X2 -> M4    Y1 -> M4            X2 and Y1 are 0.3 m apart in plan but 2.1 m apart
#                                           in invert: crossing pipes, not one node
#     M5 -> CP1                             capped dead end (CAP in the other-nodes layer)
#     M6 -> M7                              unexplained dead end
#     M9, M10, M11 -> M1, M1, M2            abandoned, proposed and removed pipes
#
# Maintenance holes: M1 has a rim (LidElevation), M2 only a depth, M3 an implausible rim (12 m),
# and M4 records 2 pipes out where the graph has 1.

NODES_X0, NODES_Y0 = 480000.0, 4760000.0
NODES_CRS = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::32617"}}
NODE_POS = {
    "M1": (0, 100),
    "M2": (50, 100),
    "PS1": (100, 100),
    "M3": (150, 100),
    "M4": (200, 100),
    "OF1": (250, 100),
    "OF2": (120, 140),
    "M8": (100, 200),
    "DUP1": (150, 200),
    "DUP2": (150.2, 200),
    "X1": (180, 160),
    "X2": (200, 150),
    "Y1": (200.3, 150),
    "M5": (50, 50),
    "CP1": (0, 50),
    "M6": (300, 50),
    "M7": (350, 50),
    "M9": (0, 150),
    "M10": (-50, 100),
    "M11": (50, 150),
}

# Inventory type of each node, as the pipes record it.
INVENTORY = {"PS1": "PS", "OF1": "OF", "OF2": "OF", "CP1": "CAP", "DUP2": "DATE"}

# (up, down, inv_up, inv_down, length, structure type, status)
NODE_PIPES = [
    ("M1", "M2", 249.00, 248.80, 50.0, "SL", "Built"),
    ("M2", "PS1", 248.70, 248.50, 50.0, "SL", "Built"),
    ("PS1", "M3", 248.40, 246.00, 50.0, "FM", "Built"),
    ("M3", "M4", 245.50, 245.00, 50.0, "SL", "Built"),
    ("M4", "OF1", 244.90, 240.50, 50.0, "SL", "Built"),
    ("M3", "OF2", 246.50, 246.00, 50.0, "OVF", "Built"),
    ("M8", "DUP1", 247.00, 246.00, 50.0, "SL", "Built"),
    ("DUP2", "M3", 245.95, 245.80, 100.0, "SL", "Built"),
    ("X1", "X2", 247.00, 246.50, 25.0, "SL", "Built"),
    ("X2", "M4", 246.40, 246.00, 50.0, "SL", "Built"),
    ("Y1", "M4", 248.50, 247.00, 50.0, "SL", "Built"),
    ("M5", "CP1", 250.00, 249.50, 50.0, "SL", "Built"),
    ("M6", "M7", 250.00, 249.00, 50.0, "SL", "Built"),
    ("M9", "M1", 250.00, 249.20, 50.0, "SL", "Abandoned"),
    ("M10", "M1", 250.00, 249.20, 50.0, "SL", "Proposed"),
    ("M11", "M2", 250.00, 249.00, 50.0, "SL", "Removed"),
]

# Maintenance-hole layer: node -> (LidElevation, Depth, PipesIn, PipesOut, FlowMonitor)
MAINTENANCE_HOLES = {
    "M1": (252.0, None, None, None, 0),
    "M2": (None, 3.0, None, None, 0),
    "M3": (12.0, None, None, None, 0),
    "M4": (251.0, None, 3, 2, 1),
    "M5": (252.0, None, None, None, 0),
    "M7": (252.0, None, None, None, 0),
    "M8": (250.0, None, None, None, 0),
    "DUP1": (250.0, None, None, None, 0),
}


def node_xy(node):
    x, y = NODE_POS[node]
    return [NODES_X0 + x, NODES_Y0 + y]


def node_pipe(i, p):
    up, down, inv_up, inv_down, length, structure, status = p
    props = {
        "OBJECTID": 5000 + i,
        "GIS_ID": 6000 + i,
        "FlowType": "SAN",
        "UpstreamID": up,
        "DownstreamID": down,
        "UpstreamInventoryType": INVENTORY.get(up, "MH"),
        "DownstreamInventoryType": INVENTORY.get(down, "MH"),
        "Length": length,
        "PipeShape": "R",
        "Diameter": 300,
        "Height": 300,
        "Width": 300,
        "Material": "PVC",
        "UpstreamInvert": inv_up,
        "DownstreamInvert": inv_down,
        "DropStructureInvert": 0.0,
        "ConstructedStatus": status,
        "StructureType": structure,
        "LocationArea": "TEST",
    }
    geometry = {"type": "LineString", "coordinates": [node_xy(up), node_xy(down)]}
    return {"type": "Feature", "properties": props, "geometry": geometry}


def node_point(node, props):
    props = {"GIS_FeatureKey": node, **props}
    geometry = {"type": "Point", "coordinates": node_xy(node)}
    return {"type": "Feature", "properties": props, "geometry": geometry}


def node_layers():
    """{file name suffix: features} of the pipes and the four node layers."""
    maintenance_holes = [
        node_point(
            node,
            {
                "LidElevation": lid,
                "Depth": depth,
                "DropStructures": 0,
                "Overflow": "N",
                "Weir": "N",
                "PipesIn": pipes_in,
                "PipesOut": pipes_out,
                "Location": "Road",
                "Verified": "Y",
                "FlowMonitor": monitor,
            },
        )
        for node, (lid, depth, pipes_in, pipes_out, monitor) in MAINTENANCE_HOLES.items()
    ]
    pump_stations = [node_point("PS1", {"Fullname": "Test Low Level P.S.", "StructureType": "PS"})]
    outfalls = [
        node_point("OF1", {"FlowType": "SAN", "PipeInvert": 240.0, "HeadwallElevation": 243.0}),
        node_point("OF2", {"FlowType": "STM", "PipeInvert": 245.5, "HeadwallElevation": 248.0}),
    ]
    other_nodes = [
        node_point("CP1", {"StructureType": "CAP"}),
        node_point("DUP2", {"StructureType": "DATE"}),
    ]
    return {
        "pipes": [node_pipe(i, p) for i, p in enumerate(NODE_PIPES)],
        "maintenance_holes": maintenance_holes,
        "pump_stations": pump_stations,
        "outfalls": outfalls,
        "other_nodes": other_nodes,
    }


if __name__ == "__main__":
    for name, features in node_layers().items():
        out = Path(__file__).with_name(f"fixture_nodes_{name}.geojson")
        fc = {"type": "FeatureCollection", "crs": NODES_CRS, "features": features}
        out.write_text(json.dumps(fc, indent=1))
        print(f"wrote {out} ({len(features)} features)")
