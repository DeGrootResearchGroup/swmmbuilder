# swmmbuilder

Python package that loads a city's GIS sewer network (GeoJSON or any format geopandas reads),
selects subnetworks of it, and writes them as SWMM `.inp` files for pyswmm.

```python
from swmmbuilder import Network, maps, swmm

net = Network("london.geojson")                     # schema=LONDON_ON by default
sub = net.catchment("OUT")                          # every node draining to OUT
sub = net.catchment("J1", max_hops=2, outlet="OUT") # truncated, plus the path down to OUT
sub = net.subnetwork(["J1", "J2", "OUT"])           # given nodes and the pipes between them

sub.validate(outfalls=["OUT"])                      # TopologyReport: splits, cycles, dead ends
sub.to_inp("model.inp", outfalls=["OUT"])           # or {"OUT": "FIXED"} for other types
model = swmm.build_model(sub, ["OUT"], manning_n_by_material={"PVC": 0.011})  # without a file
sub.to_node_list("nodes.txt")                       # Subnetwork.from_node_list(net, "nodes.txt")
maps.to_leaflet(sub, "map.html")                    # interactive map; see Maps below
```

A `Subnetwork`'s `graph` is a `networkx.MultiDiGraph` in the direction of flow, with one edge per
pipe keyed by pipe ID, so parallel pipes are kept.

## Maps

`sewernetworkx.maps` draws values computed elsewhere (concentrations, flows, each upstream
area's share of a pipe's load, ...) on the network: pandas Series or DataFrames indexed by pipe
ID, node ID or group name. Static plots use matplotlib and draw in a metric CRS (the network's
CRS if projected, else its UTM zone; see `maps.plot_crs`), so ribbon widths are in metres.

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
zones = maps.node_group_polygons(sub, {"north": ["A1", "A2"], "south": ["B1", "B2"]})
maps.plot_groups(ax, zones, zone_values, hatch_groups=[["north", "south"]])
maps.plot_pipes(ax, sub, concentration, norm="log", width=flow)   # missing values in grey
maps.plot_ribbons(ax, sub, shares, width=10, top_k=4)             # pipes x areas, rows sum to 1
maps.plot_markers(ax, sub, ["A1", "OUT"])
maps.add_basemap(ax, sub)                       # web tiles; or layer="roads.gpkg" for offline

maps.to_leaflet(sub, "map.html", [maps.PipeLayer(concentration, name="Concentration"),
                                  maps.GroupLayer(zones), maps.MarkerLayer(["OUT"])])
maps.to_geojson(sub, "pipes.geojson", values=concentration)
sub.to_geodataframe("pipes", values=concentration)  # or "nodes"; also on Network
```

Web-tile basemaps and `to_leaflet` need the optional `maps` extra
(`pip install "SewerNetworkx[maps]"`: contextily, xyzservices, folium). `generate_gmap`
(Google Maps) is deprecated in favour of `to_leaflet`.

## Logging

Progress and data problems (pipes dropped on loading, nodes with no position, height corrections)
are reported through `logging`, on the `swmmbuilder` logger; nothing is printed. To see them:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

## GIS schema

`Network(path, schema=...)` reads the GeoJSON and converts it to canonical columns (`pipe_id`,
`from_node`, `to_node`, `length`, `us_invert`, `ds_invert`, `drop_invert`, `shape`, `diameter`,
`height`, `width`, `flow_type`, `material`, `abandoned`, in m); everything after loading reads only
those (see `swmmbuilder/schema.py`). Two schemas ship with the package, as TOML files in
`swmmbuilder/schemas/`:

- `LONDON_ON` (the default): London, ON's GIS export.
- `TORONTO_ON`: the [Sewer Gravity Mains](https://open.toronto.ca/dataset/sewer-gravity-mains/)
  dataset of the City of Toronto Open Data Portal. It reprojects to MTM zone 10 (EPSG:2952), so the
  .inp coordinates are in metres. Toronto publishes only a diameter, so only circular pipes get a
  cross-section; the other shapes are reported on loading. Flow types other than `SAN`, `Storm`
  and `Combined` (`CSO`, `SCSO`, `FD`, `EO`) pass through unchanged, so e.g.
  `flow_types={"sanitary", "combined"}` keeps sanitary and combined sewers.

```python
from swmmbuilder import TORONTO_ON, Network

network = Network("sewer-gravity-main-4326.geojson", schema=TORONTO_ON)
```

For another city, write a schema file and pass its path as `schema=` (or build a `Schema` in
Python). The top-level keys are the parameters of `Schema`; `columns` maps canonical columns to
source columns, and `shapes` and `flow_types` map source codes:

```toml
name = "My city"
null_values = ["None"]   # values that mean "missing" in any mapped column
size_unit = "mm"         # m, cm, mm, ft or in; also length_unit and invert_unit
crs = "EPSG:32617"       # optional: reproject to this (projected) CRS on loading

[columns]
pipe_id = "PIPE"
from_node = "US_MH"
to_node = "DS_MH"
length = "LEN"
us_invert = "US_INV"
ds_invert = "DS_INV"
shape = "SHAPE"
diameter = "DIA"
height = "HT"
width = "WD"
flow_type = "SYSTEM"

[shapes]
CIRC = "CIRCULAR"
BOX = "RECT_CLOSED"

[flow_types]
Sanitary = "sanitary"
Storm = "storm"
```

```python
network = Network("my_city.geojson", schema="my_city.toml")
```

Only `pipe_id`, `from_node` and `to_node` are required; the pipes must name their upstream and
downstream nodes. Missing lengths (all of them, if `length` is not mapped) are computed from the
geometry, which must then be in a projected CRS (set `crs`). Shape codes that are not in `shapes`
give a warning on loading, and those pipes are rejected if they are written to an .inp.
`flow_types=` filters on canonical flow types (`{'sanitary'}` by default); source codes such as
`'SAN'` are translated too.

## SWMM template

`.inp` files are written by `swmmbuilder.swmm`. `build_model(subnetwork, ...)` returns a
`SwmmModel`: plain records of the junctions, outfalls, conduits, cross-sections, coordinates and
vertices, with node inverts, conduit offsets and Manning's n already decided, so the model can be
inspected or tested without a file. `write_inp(subnetwork, path, ...)` checks the template and
the topology, builds the model and renders it into the template in one write;
`Subnetwork.to_inp` calls it.

The default template is the one shipped in the package,
`swmmbuilder/templates/template.inp` (`examples/pyswmm_template/template.inp` is a copy kept in
sync by a test). It uses SI units, `FLOW_UNITS LPS`, because the GIS lengths, inverts and sizes are
in metres. A custom template (`template_path=`) must use `CMS`, `LPS` or `MLD`, and may only
contain the `{nodes}`, `{outfalls}`, `{storage}`, `{conduits}`, `{pumps}`, `{xsections}`,
`{curves}`, `{coordinates}` and `{vertices}` placeholders; anything else (including US units)
raises `ValueError`. Node `[COORDINATES]` and pipe `[VERTICES]` are written in the units of the GIS
geometry.

## Pumps

The GIS data do not identify pumping stations, so pumps are passed to `to_inp` by pipe ID.
Each pumped pipe is written to `[PUMPS]` (from its upstream to its downstream node) instead of
`[CONDUITS]`, with its curve in `[CURVES]`; it needs no cross-section and is not reported as
adverse. Wet wells can be written as `[STORAGE]` nodes of constant plan area instead of junctions:

```python
sub.to_inp(
    "model.inp",
    outfalls=["OUT"],
    pumps={1010: dict(curve=[(0, 6.0), (2, 4.0), (4, 0.0)], curve_type="PUMP3",
                      startup=0.5, shutoff=0.1)},
    wet_wells={"D1": dict(max_depth=3.0, area=4.0)},
)
```

`curve` is a list of (x, y) points, with y the flow in the template's `FLOW_UNITS` (L/s with the
shipped template) and x set by `curve_type`: `PUMP1` wet well volume (m³; the inlet must be a wet
well), `PUMP2`/`PUMP4` inlet depth (m), `PUMP3`/`PUMP5` head (m). `curve="*"` is an ideal pump.
`status` (`ON`/`OFF`), `startup` and `shutoff` (inlet depths, m) are optional. `Pump` and `WetWell`
objects can be passed instead of dicts.

## Conduit geometry

Cross-sections come from the schema's shape codes. For London (`PipeShape`, sizes in mm, written
in m): `R`/`Round` → `CIRCULAR(Diameter)`, `REC` → `RECT_CLOSED(Height, Width)`, `EGG` →
`EGG(Height)`, `ARC` → `ARCH`, `ELH`/`ELV` → `HORIZ_ELLIPSE`/`VERT_ELLIPSE` (`Height, Width`). An
unknown or missing shape or size raises `ValueError` listing the pipe IDs, before anything is
written.

A node's invert is the lowest invert of the subnetwork pipes connected to it, including the
drop-structure invert (`DropStructureInvert` in London) of pipes entering it. Each conduit's
inlet/outlet offset is its GIS pipe invert minus the node invert (`LINK_OFFSETS DEPTH`), so
conduits keep their GIS inverts; any that slope uphill are reported with a warning.
`height_correction={node: invert}` overrides a node's invert, and the offsets are measured from it.

## Tests

```bash
pip install -e ".[test]"
pytest
```

The tests use a small synthetic network (`tests/data/fixture_network.geojson`, regenerated by
`tests/data/make_fixture.py`) that contains deliberate edge cases. Known bugs are encoded as strict
`xfail` tests that name their issue; a fix should turn its tests green and remove the marker.

`tests/test_toronto.py` also loads a 31-pipe sanitary catchment from the Toronto dataset
(`tests/data/toronto_catchment.geojson`, fetched by `tests/data/make_toronto_fixture.py`; contains
information licensed under the Open Government Licence - Toronto).

On macOS, swmm-toolkit's wheel ships `libomp.dylib` and `libswmm-output.dylib` with invalid code
signatures, and `import pyswmm` is killed (exit 137). Re-sign them ad hoc:

```bash
codesign --force --sign - "$(python -c 'import swmm.toolkit, os; print(os.path.dirname(swmm.toolkit.__file__))')"/{libomp,libswmm-output}.dylib
```

## CI

`.github/workflows/ci.yml` runs ruff (`ruff check` and `ruff format --check`, rules in
`pyproject.toml`) and the test suite on Python 3.11–3.13. The `CI gate` job is the single check
to require in branch protection.
