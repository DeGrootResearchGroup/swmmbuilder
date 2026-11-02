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
maps.generate_gmap(sub, "map.html")                 # reprojected to WGS 84 for the map
```

A `Subnetwork`'s `graph` is a `networkx.MultiDiGraph` in the direction of flow, with one edge per
pipe keyed by pipe ID, so parallel pipes are kept.

## Node layers

The pipes name their end nodes, but say little about them. Node layers add what they are and what
the SWMM nodes need; each is optional:

```python
from swmmbuilder import Network, strip_leading_zeros

net = Network(
    "Sewer_Pipes.geojson",
    maintenance_holes="Maintenance_Holes.geojson",
    pump_stations="Pumping_Stations.geojson",
    outfalls="Sewer_Outfalls.geojson",
    other_nodes="Sewer_Other_Nodes.geojson",  # caps, tees, chambers, ...
    include_overflows=False,                   # dry-weather network: drop overflow pipes
    normalize_ids=strip_leading_zeros,         # 'OK0329' (older exports, node lists) -> 'OK329'
    snap_tolerance=0.5,                        # merge node IDs within 0.5 m ...
    max_invert_difference=1.0,                 # ... unless their inverts differ by more than 1 m
)
net.nodes            # GeoDataFrame: node_type, name, invert, rim_elevation, depth, ...
net.snap_report      # merged pairs: both IDs, distance, both inverts, kept ID or why not merged
sub = net.catchment("WT6011")
sub.validate()       # dead ends classified as capped / outfall / pump_station / unexplained
sub.to_inp("model.inp")  # outfalls default to the outfall nodes of the subnetwork
```

- **Node types.** `net.nodes` has one row per node of the graph, indexed by node ID, with its
  point; the same values are node attributes of the graph. A node takes the type of its layer
  (`maintenance_hole`, `pump_station`, `outfall`, or the other-nodes layer's own type such as
  `cap`, `tee` or `chamber`); a node in no layer takes the pipes' inventory type
  (`UpstreamInventoryType`/`DownstreamInventoryType` in London), if the schema maps it.
- **Pipe status and structure type.** `status={'built'}` (the default) keeps built pipes only;
  London's `ConstructedStatus` also has `abandoned`, `proposed` and `removed`. Older exports
  without it use the `Abandoned` flag, and `keep_abandoned=True` still works. `structure_type` is
  `gravity`, `overflow` or `force_main`; `include_overflows=False` drops overflow pipes, which
  carry no dry-weather flow and make spurious splits and dead ends. `Network.removed` counts
  what was dropped by flow type, status and structure type.
- **Outfalls from data.** With no `outfalls`, `validate`, `to_inp` and `swmm.build_model` use the
  outfall-type nodes of the subnetwork (`Subnetwork.outfalls`); an explicit list still works.
- **Coincident nodes.** The same structure sometimes has two IDs at one point, which disconnects
  the network there. `net.snap_nodes(tolerance=0.5, max_invert_difference=None)` (or
  `snap_tolerance=` on loading) merges node IDs whose points are within the tolerance, preferring
  the ID in the maintenance-hole layer. Pairs whose inverts differ by more than
  `max_invert_difference` (pipes crossing rather than joining) and pairs already joined by a pipe
  are not merged. The merged IDs are kept in `net.aliases`, and every method that takes node IDs
  accepts either one. Geographic coordinates are projected to the local UTM zone for this.
- **ID normalisation.** `normalize_ids` is applied to the node IDs of the pipes, the node layers,
  and those passed to `catchment`, `subnetwork` and `Subnetwork.from_node_list`.
- **Validation.** `TopologyReport` also has `dead_end_types`, `force_mains`, and
  `pipe_count_mismatches`: maintenance holes whose recorded `pipes_in`/`pipes_out` differ from
  the network's graph.
- **Maintenance-hole attributes.** `rim_elevation` is the layer's (London: `LidElevation`, or
  `HeadwallElevation` for outfalls), else `invert + depth`; rims outside the schema's
  `rim_range` (London: 150–350 m) are treated as missing. The `.inp` writer uses rim − invert
  as a junction's maximum depth where it is known, and `max_depth` elsewhere. `has_drop`,
  `overflow`, `weir`, `pipes_in`, `pipes_out`, `flow_monitor`, `access` and `verified` are passed
  through for downstream use, e.g. choosing accessible or flow-monitored sampling sites.

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
`height`, `width`, `flow_type`, `material`, `status`, `abandoned`, `structure_type`,
`us_node_type`, `ds_node_type`, in m), and the node layers to canonical node columns (`node_id`,
`node_type`, `name`, `invert`, `rim_elevation`, `depth`, `has_drop`, `overflow`, `weir`,
`pipes_in`, `pipes_out`, `flow_monitor`, `access`, `verified`); everything after loading reads
only those (see `swmmbuilder/schema.py`). These schemas ship with the package, as TOML files in
`swmmbuilder/schemas/`:

- `LONDON_ON` (the default): London, ON's GIS export.
- `TORONTO_ON`: the [Sewer Gravity Mains](https://open.toronto.ca/dataset/sewer-gravity-mains/)
  dataset of the City of Toronto Open Data Portal. It reprojects to MTM zone 10 (EPSG:2952), so the
  .inp coordinates are in metres. Toronto publishes only a diameter, so only circular pipes get a
  cross-section; the other shapes are reported on loading. Flow types other than `SAN`, `Storm`
  and `Combined` (`CSO`, `SCSO`, `FD`, `EO`) pass through unchanged, so e.g.
  `flow_types={"sanitary", "combined"}` keeps sanitary and combined sewers.
- Six more Ontario cities' open data, each reprojected to its UTM zone 17N CRS:
  - `WINDSOR_ON`: the SanitarySewer, StormSewer and CombinedSewer layers.
  - `KITCHENER_ON`: Sanitary Pipes. `CATEGORY` is the structure type.
  - `WATERLOO_ON`: Sanitary Sewer Gravity Mains. All pipes are taken as circular, and an invert
    of 111.111 means none.
  - `PETERBOROUGH_ON`: Sanitary (and Storm) Sewer Pipes. Only circular pipes get a cross-section.
  - `BARRIE_ON`: Sanitary Pipes. `TYPE` is the structure type.
  - `SARNIA_ON`: Sanitary and Storm Sewers. All pipes are taken as circular; non-circular sizes
    (`"W x H"`) are rejected if written.

  Each TOML file's header gives the dataset's URL and its quirks. The layers are ArcGIS services
  that return at most 1,000–2,000 pipes per query, so download them from the portal (GeoJSON) or
  page through the REST query.

```python
from swmmbuilder import TORONTO_ON, Network

network = Network("sewer-gravity-main-4326.geojson", schema=TORONTO_ON)
```

For another city, write a schema file and pass its path as `schema=` (or build a `Schema` in
Python). The top-level keys are the parameters of `Schema`; `columns` maps canonical columns to
source columns, and `shapes` and `flow_types` map source codes (or set `default_shape`, e.g.
`"CIRCULAR"`, where there is no shape column):

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

Pipe statuses, pipe structure types and node types are mapped from source codes by the
`[statuses]`, `[structure_types]` and `[node_types]` tables, and the node layers' columns by
`[node_columns]`, where a list names alternatives (the first one present in a layer is used);
see `swmmbuilder/schemas/london_on.toml`. The Toronto preset does not map node layers yet.

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

Pumping stations and force mains can be identified from the node layers and the pipes' structure
type (see above), but the GIS data have no pump curves, so pumps are passed to `to_inp` by pipe ID.
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
`tests/data/make_fixture.py`) that contains deliberate edge cases. The same script writes a second
network in London's newer export schema with its four node layers
(`tests/data/fixture_nodes_*.geojson`), used by `tests/test_nodes.py`. Known bugs are encoded as strict
`xfail` tests that name their issue; a fix should turn its tests green and remove the marker.

`tests/test_toronto.py` also loads a 31-pipe sanitary catchment from the Toronto dataset
(`tests/data/toronto_catchment.geojson`, fetched by `tests/data/make_toronto_fixture.py`; contains
information licensed under the Open Government Licence - Toronto). `tests/test_ontario.py` does
the same for Barrie, Kitchener, Peterborough, Sarnia, Waterloo and Windsor, with a 14–46 pipe
sanitary catchment each (`tests/data/<city>_catchment.geojson`, fetched by
`tests/data/make_ontario_fixtures.py`, whose docstring gives each city's licence).

On macOS, swmm-toolkit's wheel ships `libomp.dylib` and `libswmm-output.dylib` with invalid code
signatures, and `import pyswmm` is killed (exit 137). Re-sign them ad hoc:

```bash
codesign --force --sign - "$(python -c 'import swmm.toolkit, os; print(os.path.dirname(swmm.toolkit.__file__))')"/{libomp,libswmm-output}.dylib
```

## CI

`.github/workflows/ci.yml` runs ruff (`ruff check` and `ruff format --check`, rules in
`pyproject.toml`) and the test suite on Python 3.11–3.13. The `CI gate` job is the single check
to require in branch protection.
