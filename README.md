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

## Install

What changed between versions is in
[CHANGELOG.md](https://github.com/DeGrootResearchGroup/swmmbuilder/blob/main/CHANGELOG.md).

```bash
pip install swmmbuilder                 # from PyPI
pip install "swmmbuilder[maps]"         # + web-tile basemaps and Leaflet maps (contextily, xyzservices, folium)
pip install "swmmbuilder[simulate]"     # + pyswmm, to run the .inp files
pip install "swmmbuilder[census]"       # + pycancensus, for Canadian census polygons (CensusMapper)
pip install "swmmbuilder[elevation]"    # + rasterio, to estimate missing rims from a DEM
pip install -e .                        # from a clone of this repository (extras as above)
```

swmmbuilder only writes `.inp` files and never imports pyswmm, so pyswmm is optional. pyswmm
depends on `julian`, an sdist-only package whose `setup.py` imports `distutils`, which fails to
build on some setups (e.g. `AttributeError: install_layout`); swmmbuilder itself installs and works
without it.

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
    split_tolerance=50.0,                      # split IDs reused for structures > 50 m apart
)
net.nodes            # GeoDataFrame: node_type, name, invert, rim_elevation, depth, ...
net.snap_report      # merged pairs: both IDs, distance, both inverts, kept ID or why not merged
net.split_report     # split IDs: new ID, its pipes, distance to the other structures, attributes
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
- **Reused IDs.** One ID is sometimes used for two structures far apart (a newer subdivision
  given the IDs of an older one, or a mistyped ID), which joins the two in the graph and can make
  false cycles. `split_tolerance=` on loading groups the pipe ends that carry each ID into
  structures (ends more than the tolerance apart, with no chain of closer ends between them,
  are different structures) and gives each structure its own node, named by appending `a`, `b`,
  `c`, ... in the order of its smallest pipe ID. A pipe drawn from its downstream node to its
  upstream node has its ends matched the other way round, so it does not split its nodes. A
  node-layer point with a split ID goes to the nearest structure; a structure with no point takes
  its type and invert from its pipes. Each split is logged as a warning and listed in
  `net.split_report`, and `net.splits` maps each original ID to its new IDs. Splitting runs
  before snapping, and the tolerance must exceed `snap_tolerance`: a mistyped end that lies on
  another structure is then merged into it by the snap, and the structures of one ID are never
  merged again.
- **ID normalisation.** `normalize_ids` is applied to the node IDs of the pipes, the node layers,
  and those passed to `catchment`, `subnetwork` and `Subnetwork.from_node_list`.
- **Validation.** `TopologyReport` also has `dead_end_types`, `force_mains`, and
  `pipe_count_mismatches`: maintenance holes whose recorded `pipes_in`/`pipes_out` differ from
  the network's graph. For pumping stations it lists `stations_without_force_main`,
  `orphan_force_mains` (force mains that start at no pump station) and `force_main_ends`
  (force-main chains that join another force main or end at a dead end).
- **Maintenance-hole attributes.** `rim_elevation` is the layer's (London: `LidElevation`, or
  `HeadwallElevation` for outfalls), else `invert + depth`; rims outside the schema's
  `rim_range` (London: 150–350 m) are treated as missing. The `.inp` writer uses rim − invert
  as a junction's maximum depth where it is known, and `max_depth` elsewhere (or a rim estimated
  from a ground-elevation model: see [Rims from a DEM](#rims-from-a-dem)). `has_drop`,
  `overflow`, `weir`, `pipes_in`, `pipes_out`, `flow_monitor`, `access` and `verified` are passed
  through for downstream use, e.g. choosing accessible or flow-monitored sampling sites.

## Maps

`swmmbuilder.maps` draws values computed elsewhere (concentrations, flows, each upstream
area's share of a pipe's load, ...) on the network: pandas Series or DataFrames indexed by pipe
ID, node ID or group name. Static plots use matplotlib and draw in a metric CRS (the network's
CRS if projected, else its UTM zone; see `maps.plot_crs`). Each one styles the axes as a map
(no ticks or frame, a scale bar; `maps.map_axes`, or `style=False`) and returns a `DrawnLayer`
with its artist and its key.

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
zones = maps.node_group_polygons(sub, {"Z1": ["A1", "A2"], "Z2": ["B1", "B2"], "rest": [...]})
colours = maps.category_colors(zones.index)       # one palette for zones and ribbons
groups = maps.plot_groups(ax, zones, colors=colours, label="Zone", legend=False)
pipes = maps.plot_pipes(ax, sub, concentration, norm="log", width=flow, legend=False)
ribbons = maps.plot_ribbons(ax, sub, loads, colours, title="Upstream zone", legend=False)
site = maps.plot_markers(ax, sub, ["OUT"], labels={"OUT": "Treatment plant"}, legend=False)
maps.add_legends(ax, groups, pipes, ribbons)      # stacked to the right of the map
maps.add_basemap(ax, sub)                         # web tiles; or layer="roads.gpkg" for offline

maps.to_leaflet(sub, "map.html", [maps.PipeLayer(concentration, name="Concentration"),
                                  maps.GroupLayer(zones), maps.MarkerLayer(["OUT"])])
maps.to_geojson(sub, "pipes.geojson", values=concentration)
sub.to_geodataframe("pipes", values=concentration)  # or "nodes" (rows of net.nodes)
```

- **Pipes.** `plot_pipes` colours pipes by value (missing values in grey) and, with `width=`,
  scales line widths by the square root of the value by default (`width_scale="linear"` or
  `"log"`), drawing thick pipes last. `context=True` draws them lighter and thinner, as context
  beneath zones or ribbons.
- **Colour scales.** `cmap`, `norm`, `vmin`, `vmax` and `missing_color` can be given together
  as a `maps.ColourScale` (`colour=` on `plot_pipes` and `plot_groups`, or on `PipeLayer` and
  `GroupLayer`), e.g. so that several plots share one scale:
  `scale = maps.ColourScale(norm="log", vmin=1e-3, vmax=10)`, then
  `maps.plot_pipes(ax, sub, before, colour=scale)` and `maps.PipeLayer(after, colour=scale)`.
- **Ribbons.** `plot_ribbons` splits each pipe into strands, one per category, as wide as its
  share. By default the ribbon width is proportional to the pipe's total load, with the largest
  2.5 % of the network's extent wide (`max_width=` in metres to override); ribbons under 5 % of
  that (`min_width=`) leave their pipes as thin grey context. Pipes are joined into one ribbon
  along main stems (the widest ribbon through each junction) until the width or a share
  visibly changes (`merge_tol=`), and the joins are mitred, so a trunk reads as one band.
- **Zones.** `node_group_polygons` makes Voronoi tiles (or concave hulls) clipped to the network
  footprint, the area near the pipes, so zones cover the serviced area. Its default distance is
  twice the median pipe length, at most 5 % of the extent (`maps.footprint_distance`;
  `footprint=` in metres to override, `False` for the convex hull).
  `plot_groups` fills them lightly by category (with a legend) or by value (with a colour bar),
  and can hatch sets of zones.
- **Keys.** With `legend=True` (the default) each plot adds its colour bar or legend to a column
  to the right of the map, below those already there. To compose several layers, pass
  `legend=False` and call `maps.add_legends(ax, *layers)` once.
- **Units.** Pipe line widths and marker sizes are in points: they stay the same when you zoom.
  Ribbon widths, footprints and buffers are in metres: they scale with the map. Their defaults
  are relative to the size of the network drawn, so they suit a street and a city alike.

`python examples/maps_gallery.py` renders every static plot on two Toronto test catchments, of
31 pipes (about 560 m across) and 2,692 pipes (about 9 km across), so the defaults can be checked
at both scales; `--pipes` and `--outlet` render another network with the London schema. CI
uploads the images as the `maps-gallery` artifact. Web-tile basemaps and `to_leaflet` need the
optional `maps` extra (`pip install "swmmbuilder[maps]"`: contextily, xyzservices, folium).
`generate_gmap` (Google Maps) is deprecated in favour of `to_leaflet`.

## Diagnostics

Problems are reported in four ways, by what you can do about them (see
`swmmbuilder.diagnostics`):

| Channel | When | Examples |
|---|---|---|
| `ValueError` | the output would be wrong, or SWMM would refuse it; every problem found by one check is listed in one error | gaps that stop a pipe being written, invalid pumps, flow splits under KINWAVE |
| `warnings.warn` | the model is written as asked, but suspect data or an option you should act on affects it | `DataWarning`: adverse slopes, a pipe below its node, a rim not above the invert, unknown shape codes; `DesignWarning`: `force_main_checks`, a pattern whose mean is not 1; `OptionWarning`: force mains left unpumped, `height_correction` keys that name no node |
| `logger.warning` | input dropped or repaired on loading, or problems that do not change the hydraulics | pipes with no ID or length, nodes left out of `[COORDINATES]`, colours repeating in a plot |
| `logger.info` | an audit trail | values filled by the gap rules or derived by design rules, pipes removed by the filters |

The warning categories derive from `swmmbuilder.SwmmBuilderWarning` (a `UserWarning`), so they
can be filtered by kind, e.g. `warnings.filterwarnings("ignore", category=swmmbuilder.DataWarning)`
or turned into errors with `"error"`. The log goes to the `swmmbuilder` logger; nothing is
printed. To see it:

```python
import logging
logging.basicConfig(level=logging.INFO)
```

### Capturing the log of a build

`logging.basicConfig` does nothing once the root logger has a handler, so a second build in the
same process gets no log, and it never sees the warnings. `capture_log` collects both for one
build:

```python
with swmmbuilder.capture_log(level="INFO", path=out / "build_log.txt") as captured:
    ...  # build the model
print(captured.text)
```

It attaches a handler to the `swmmbuilder` logger only (lowering the logger's level for the block
if needed) and removes it on exit, even after an error, so each use gets its own text. Warnings of
`SwmmBuilderWarning` and its subclasses are added to the text and still shown as usual. Python
shows a given warning once per location by default, so the text holds the warnings that were
shown. Lines have the `basicConfig` format, `INFO:swmmbuilder.gaps:...`, with the category in
place of the logger for a warning, `WARNING:DataWarning:...`. With `path=`, the text is written
to that file (UTF-8) on exit.

## GIS schema

`Network(path, schema=...)` reads the GeoJSON and converts it to canonical columns (`pipe_id`,
`from_node`, `to_node`, `length`, `us_invert`, `ds_invert`, `drop_invert`, `shape`, `diameter`,
`height`, `width`, `flow_type`, `material`, `status`, `abandoned`, `structure_type`,
`us_node_type`, `ds_node_type`, in m), and the node layers to canonical node columns (`node_id`,
`node_type`, `name`, `invert`, `rim_elevation`, `depth`, `has_drop`, `overflow`, `weir`,
`pipes_in`, `pipes_out`, `flow_monitor`, `access`, `verified`); everything after loading reads
only those (see `swmmbuilder/schema.py`). These schemas ship with the package, as TOML files in
`swmmbuilder/schemas/`:

- `LONDON_ON` (the default): London, ON's GIS export. The City's GeoJSON is in WGS 84; it
  reprojects to UTM zone 17N (EPSG:2958, NAD83(CSRS)), so that missing lengths can be computed
  and the .inp coordinates are in metres.
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
crs = "EPSG:32617"       # optional: reproject to this (projected) CRS on loading (data with
                         # no CRS are taken to be in it already)

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
`{curves}`, `{dwf}`, `{patterns}`, `{coordinates}` and `{vertices}` placeholders; anything else
(including US units) raises `ValueError`. Node `[COORDINATES]` and pipe `[VERTICES]` are written in
the units of the GIS geometry.

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

An ideal pump must be the only outlet link of its inlet node (SWMM refuses it otherwise, with
ERROR 134), so a station with an overflow needs a pump curve.

### Pumping stations with force mains

`pumps=` replaces a pipe with a pump, which suits a short lift. For a pumping station with a
force main, use `pump_stations=`, keyed by the station node, so that the force main keeps its
length and volume (and so its travel time):

```python
from swmmbuilder import LONDON_ON, PumpStation

sub.to_inp(
    "model.inp",
    pump_stations={"SL6001": PumpStation(
        pump=dict(curve=[(0.0, 120.0)], curve_type="PUMP2", startup=1.5, shutoff=0.5),
        wet_well=dict(max_depth=4.0, area=12.0),
    )},
    force_main_roughness=LONDON_ON.design["force_main_roughness"],
)
```

For each station the writer adds a discharge node `<node>_D` at the station, with the invert of
the force main's upstream end, and a pump `<node>_P` from the station node (a wet well if one is
given, else a junction) to it. The force-main pipes leaving the station (`StructureType` `FM` in
London), followed downstream through every node whose only outlet is a force main
(`Subnetwork.force_main_chain`), start at `<node>_D` and stay conduits with their GIS lengths and a
`FORCE_MAIN` cross-section; their interior nodes and `<node>_D` get `force_main_surcharge_depth`
(default 100 m) so that the main can run pressurised under `DYNWAVE`. Where two stations' force
mains join, the shared main downstream is a force main of both; a chain stops at the next pump
station, whose own pump starts its force main. Other outlets of the station, such as overflows,
stay conduits (`include_overflows=False` drops them). A pump-station node with force mains that
is not in `pump_stations` gets a warning: its force mains are written as ordinary conduits.

- **Pump operation.** There is no built-in default pump. An ideal pump (`curve="*"`) approximates
  a variable-speed station that delivers a near-uniform flow; a curve with `startup`/`shutoff`
  depths (e.g. `PUMP2`) a fixed-speed station whose on/off cycling pulses the force-main flow. A
  station's curve needs its `curve_type`. `default_pump=` gives the pump of the stations given
  without one; London's preset, an ideal pump for its variable-speed stations (s. 4.3.5), is
  `LONDON_ON.design["default_pump"]`. A startup depth must be above the shutoff depth and below
  the wet well's max depth, or the writer refuses.
- **Twin force mains.** `force_mains=` picks the force mains to pump into (the others are written
  as ordinary conduits, with a warning). With more than one, `discharge="shared"` writes one pump
  feeding them all and `discharge="separate"` one pump each (`<node>_P1`, `<node>_D1`, ...);
  without it the writer refuses.
- **Roughness.** `force_main_roughness` is required and is read by SWMM according to the
  template's `FORCE_MAIN_EQUATION` (Hazen–Williams C for `H-W`, the shipped template; roughness
  height in mm for `D-W`). It is a constant, `{"by": "diameter", "table": [(upper_mm, value),
  ...]}`, `{"by": "material", "table": {material: value}}` or a callable `pipe -> value`; a dict
  may declare `"equation"`, and a template with another equation is refused (a constant or a
  callable has no declared unit, so it is not checked).
  `force_main_roughness_by_pipe={pipe_id: value}` overrides the rule. Schemas carry design-rule
  presets, never applied unless passed: `LONDON_ON.design["force_main_roughness"]` is London's
  C by diameter (Design Specifications & Requirements Manual, s. 4.8).
- **Wet-well defaults (opt-in).** When a station's wet well is not known,
  `wet_well_defaults={"fill_time": minutes, "high_level_below_inlet": m, "low_level_above_pump": m}`
  derives one from the station's `design_flow` (design average flow, in the template's
  `FLOW_UNITS`: L/s with the shipped template; the GIS data do not carry it, so it is given per
  station) and the invert of its inlet sewer (the lowest gravity pipe entering it, not a force
  main or overflow). The station
  node invert is the wet-well floor, taken as the pump level: the shutoff depth is
  `low_level_above_pump` (or, with `"low_level_suction_diameters": n` in the rule and a
  station's `suction_diameter` in mm, n suction diameters if that is higher), the startup depth
  is the inlet invert minus `high_level_below_inlet`,
  the plan area holds `design_flow` for `fill_time` between the two, and the well rises to the
  rim if known, else to the crown of the inlet sewer. On/off depths are derived only for a pump
  with a curve and no `startup`/`shutoff`. Explicit `WetWell` (or `wet_wells`) and pump values
  always win. Every derived value is logged and listed in `model.derived`
  (`model.derived_report()` as a table), so nothing derived looks like data. London's preset
  (s. 4.3.10: 30 min fill, high water 300 mm below the inlet, low water 300 mm, or two suction
  diameters, above the pump) is
  `LONDON_ON.design["wet_well_defaults"]`:

  ```python
  model = sub.to_inp(
      "model.inp",
      pump_stations={"SL6001": PumpStation(
          pump=dict(curve=[(0.0, 120.0)], curve_type="PUMP2"), design_flow=40.0,
      )},
      force_main_roughness=LONDON_ON.design["force_main_roughness"],
      wet_well_defaults=LONDON_ON.design["wet_well_defaults"],
  )
  print(model.derived_report())
  ```

- **Force-main checks (warnings).** `force_main_checks={"min_diameter": mm, "velocity": (min,
  max)}` warns of force mains narrower than `min_diameter` and of those whose velocity, the
  pumped flow over the full pipe area, is outside `velocity` (m/s). The pumped flow is the pump
  curve's largest flow, or the station's `design_flow` for an ideal pump or a head-flow curve
  (PUMP3, PUMP5, whose largest flow is at zero head; not checked without one), summed over the
  stations feeding a shared main after a join; a shared pump feeding several force mains is not
  checked, as the split is unknown.
  London's preset (100 mm, 0.6–3.0 m/s) is `LONDON_ON.design["force_main_checks"]`.
- A pipe in both `pumps` and a station's force main, a station with no force-main outlet, and a
  wet well given both in `wet_wells` and a station are errors. A schema's design-rule presets are
  checked when it is made, so a typo in a TOML preset is caught on loading.

## Dry-weather flow

`dwf={node: base_flow}` loads nodes with a dry-weather flow in the template's `FLOW_UNITS` (L/s
with the shipped template), written to `[DWF]`. This also covers point loads at nodes with no
upstream pipes, such as a pumping station fed by a system that is not in the data. Time patterns
go in `patterns={name: multipliers}` (`[PATTERNS]`); the type comes from the number of multipliers
(24 → `HOURLY`, 7 → `DAILY`, 12 → `MONTHLY`) or is given as `(type, multipliers)`, e.g.
`("WEEKEND", [...])`. `dwf_patterns` applies them: one name (or list of names) for every loaded
node, or `{node: [names]}`, up to four per node, one of each type.

```python
from swmmbuilder import dwf_from_population

dwf = dwf_from_population({"A1": 1200, "B1": 800}, per_capita=250, unit="L/d")  # -> L/s
sub.to_inp(
    "model.inp",
    outfalls=["OUT"],
    dwf={**dwf, "D1": 3.5},
    patterns={"diurnal": [0.5, 0.4, 0.35, ..., 0.8]},  # 24 hourly multipliers, mean 1
    dwf_patterns="diurnal",
)
```

Nodes must be in the subnetwork and not outfalls, flows non-negative and finite, and patterns
defined, of the right length and positive; otherwise `ValueError` is raised before anything is
written. An hourly pattern whose mean is not 1 gives a warning, since it scales the base flows.
The total base flow and the number of loaded nodes are logged, for comparison with the plant's
measured dry-weather flow. `dwf_from_population` converts a per-capita flow (`unit` `L/d`,
`m3/d`, `L/s`, or `m3/s`) to the template's `flow_units` (default `LPS`).

### Population from polygons

`allocate_population` places the population of areal units (census areas, postal code areas, or
any polygons) at nodes. Each pipe goes to the polygon containing its midpoint, each polygon's
population is shared among its pipes in proportion to a weight, and each pipe's share is placed at
its upstream node:

```python
from swmmbuilder import allocate_population

alloc = allocate_population(
    net,                 # the whole Network, not the catchment (see `within` below)
    das,                 # GeoDataFrame indexed by polygon ID
    population="pop",
    weight="PDC_Count",  # London's private drain connections per pipe; default: length
    max_weight=100,      # larger counts are invalid, as are missing and negative ones
)
sub = net.catchment("OUT")
dwf = dwf_from_population(alloc.within(sub), per_capita=230, unit="L/d")
```

- `max_distance` replaces the midpoint rule with a distance rule: each polygon is served by every
  pipe within that distance of it (in the CRS's units), with each pipe's whole weight. This suits
  small units such as census blocks, whose edges follow the streets the sewers run under: a pipe
  then serves the blocks on both sides. A polygon with no pipe that close contributes no one,
  unless it is flagged in `sewered` (a column of booleans, e.g. residential zoning): a flagged
  polygon left with no pipe, such as a development on private sewers, is given its nearest pipe,
  up to `max_nearest_distance` away (no limit by default). `dominant_category` gives such a
  flag from a zoning layer, the category covering the largest area of each polygon:

  ```python
  from swmmbuilder import dominant_category

  use = dominant_category(blocks, zoning, "GeneralizedLandUse")  # category, share
  blocks["residential"] = use["category"] == "RESIDENTIAL"
  ```

  The zoning is reprojected to the polygons' CRS (projected in metres), a tie goes to the
  category that sorts first, and a polygon meeting no zoning has a missing category (so its flag
  is False).
- Force mains take no service connections: they serve no polygon and weigh 0, unless
  `serve_force_mains=True`.
- `weight` is a pipe column, a function of the pipes, or None for the length.
- `impute` is `"per_length"` (default: an invalid weight becomes the pipe's length times the
  median weight per metre of the valid pipes of its own polygon) or `"zero"`. The polygons passed
  serve for both the share and the imputation.
- `fallback_weight` (default 1, e.g. one connection) is given to an invalid weight where the
  polygon's median is 0 or there is no valid weight, and then to every pipe of a polygon whose
  pipes still weigh nothing, so its people are shared equally; `None` leaves such polygons
  unallocated.
- `placement` is `"upstream"`, `"downstream"`, or `"split"` (half at each end).

A polygon with no pipe contributes no one, and a partly sewered polygon places all its people on
its pipes. `alloc.unallocated` is the population left out, and `alloc.no_pipes`,
`alloc.no_valid_weight`, and `alloc.zero_weight` (the last two only with no fallback) name the
polygons it comes from; `alloc.nearest` gives the flagged polygons sent to their nearest pipe,
with the distance, and `alloc.fallback` those given the fallback weight. `alloc.pipes` has
one row per pipe and polygon it serves, `alloc.by_pipe()` sums each pipe's shares, and
`alloc.polygons` has one row per polygon. `alloc.within(sub)` keeps the shares of the
subnetwork's own pipes, so people on a pipe leaving it are not counted.

### Census polygons

`swmmbuilder.census.census_polygons` gets the polygons covering a network from a source, keeps
those meeting the footprint of its gravity pipes (whole, not clipped), and returns them indexed by
ID with `population` and `parent` columns, ready for `allocate_population`. For Canada's 2021
Census, the source is CensusMapper, through pycancensus (`pip install "swmmbuilder[census]"`).
The level is `"DB"` (dissemination blocks: usually one city block, bounded by streets) or `"DA"`
(dissemination areas: groups of blocks with about 400 to 700 people).

```python
from swmmbuilder.census import CensusMapper, census_polygons

blocks = census_polygons(
    net,                                  # the whole Network
    "DB",
    CensusMapper(api_key=my_key),         # your own key; or set CENSUSMAPPER_API_KEY
    cache="data/census/blocks_2021.gpkg", # read if present: later builds are offline
    margin=200,                           # m around the pipes; at least max_nearest_distance
)
alloc = allocate_population(net, blocks, max_distance=15)
```

- **The API key is yours to provide.** swmmbuilder has none. A free key comes with a
  [CensusMapper](https://censusmapper.ca) account, and its quota is yours. It is passed in, or read
  from your `CENSUSMAPPER_API_KEY` (or `CANCENSUS_API_KEY`) environment variable, and only needed
  when there is no cache. swmmbuilder never stores, logs, or shows it.
- **Citing the data is your responsibility.** swmmbuilder ships no census data. Statistics
  Canada's data are under the
  [Statistics Canada Open Licence](https://www.statcan.gc.ca/en/reference/licence), which asks
  for a source notice, and CensusMapper asks that cancensus be cited. The cache records the
  source, dataset, regions, and fetch date.
- With no `regions=`, `CensusMapper` finds the census subdivisions (municipalities) meeting the
  network with CensusMapper's intersecting-geometries query; `regions={"CSD": ["3539036"]}` names
  them instead.
- The source is pluggable: `PolygonFile(path, id=..., population=...)` reads any local file (another
  country's census, say), and any object with a `name` and a `fetch(area, crs, level)` method
  returning such polygons will do.

## Conduit geometry

Cross-sections come from the schema's shape codes. For London (`PipeShape`, sizes in mm, written
in m): `R`/`Round` → `CIRCULAR(Diameter)`, `REC` → `RECT_CLOSED(Height, Width)`, `EGG` →
`EGG(Height)`, `ARC` → `ARCH`, `ELH`/`ELV` → `HORIZ_ELLIPSE`/`VERT_ELLIPSE` (`Height, Width`). A
round pipe with no diameter whose height equals its width takes that value as its diameter (logged
and flagged; see below). An unknown or missing shape or size raises `ValueError`, before anything
is written.

A node's invert is the lowest invert of the subnetwork pipes connected to it, including the
drop-structure invert (`DropStructureInvert` in London) of pipes entering it. Each conduit's
inlet/outlet offset is its GIS pipe invert minus the node invert (`LINK_OFFSETS DEPTH`), so
conduits keep their GIS inverts; any that slope uphill are reported with a warning.
`height_correction={node: invert}` overrides a node's invert, and the offsets are measured from it.

## Data gaps

`sub.data_gaps()` lists every gap in a subnetwork's data at once: pipes with an unknown or missing
shape code (`unknown_shape`) or a missing size (`no_size`), pipe ends with no invert
(`no_us_invert`, `no_ds_invert`; written at the node invert), and nodes with no invert from any pipe
(`no_node_invert`). `validate()` reports their counts in `data_gaps`, and `to_inp` refuses a
subnetwork with any `unknown_shape`, `no_size` or `no_node_invert` gap left, listing all of them.

Gaps can be filled by rules that the user chooses; none is applied by default. Every filled value
is logged and flagged in the report's `filled_from` column (with `value`, `source` and the `order`
of the fills), so nothing filled looks like City data. The same options go to `data_gaps` (to see
what would be filled) and `to_inp` (whose returned model keeps the report as `model.data_gaps`):

```python
from swmmbuilder import LONDON_ON

fills = dict(
    size_overrides={145041: 0.25, 207149: dict(height=0.6)},  # values from the City, in m
    shape_codes={"": "CIRCULAR"},       # missing and blank shape codes are ""
    size_fill="downstream",             # copy from the downstream pipe of the same shape,
                                        # else the largest upstream one ("upstream": reversed)
    invert_fill="interpolate",          # by length between the nearest known inverts
    headwater_slope=0.005,              # no known invert upstream: downstream invert + slope
    force_main_discharge=LONDON_ON.design["force_main_discharge"],
)
report = sub.data_gaps(**fills)
sub.to_inp("model.inp", **fills)
```

`force_main_discharge` sets the missing downstream invert of a force main's last pipe from the
receiving node's invert (known, or filled first by `invert_fill`/`headwater_slope`), plus
`height_above_receiving_invert` (default: `max_height_above_receiving_invert`); a callable
`(pipe, receiving_invert) -> invert` can be passed instead. Design values like this sit in the
schema's optional `[design]` table, apart from the column mapping, and are used only when passed:
London's (`london_on.toml`) is the Design Specifications & Requirements Manual's limit of 200 mm
above the receiving maintenance hole's flow line (Ch. 4, s. 4.8).

### Rims from a DEM

A junction's maximum depth is rim − invert, so a junction with no usable rim (none recorded, or one
at or below the invert) takes the default `max_depth`, whatever the ground above it. Every such
node is listed in `data_gaps()` as a `no_rim` gap (not blocking). Maintenance-hole lids in streets
are generally close to flush with the pavement, so a bare-earth ground-elevation model (DEM) from
airborne lidar, such as the Ontario Digital Terrain Model (Lidar-Derived), estimates the missing
rims. `Network.estimate_rims` samples a GeoTIFF (or a list of tiles, in any CRS) at the nodes,
compares it with the recorded rims, and estimates the missing ones; it does not change the network.
It needs the `elevation` extra (`pip install "swmmbuilder[elevation]"`: rasterio).

```python
est = net.estimate_rims(["tile_1.tif", "tile_2.tif"], window=1.0, offset="median")
print(est.summary())                 # recorded rim - DEM: median, IQR, 5-95 %, largest
est.largest_differences(10)          # recorded rims that are probably wrong, or raised lids
est.samples                          # DEM value, difference, estimate and status of each node
sub.to_inp("model.inp", rims=est.rims, max_depth="median")
```

- `window=0` (default) takes the cell at the node; `window=1.0` the median (or `statistic="mean"`)
  of the cells within 1 m, to soften curbs and ditches.
- Only nodes with no usable rim are estimated, as the DEM value plus `offset` (a number, or
  `"median"`: the median recorded rim − DEM). A DEM rim at or below the invert is not used (a
  `DataWarning` gives the count). Usable recorded rims are never replaced.
- The `rims` fill option (of `to_inp`, `build_model`, `data_gaps`, and `fill_gaps`) fills the
  `no_rim` gaps: rims from `est.rims` are flagged as filled from `estimate_rims`, with the DEM as
  their source, and any other `{node: rim}` mapping as filled from `rims`. A given rim at or below
  the node's invert is not used, with a `DataWarning`.
- Filled rims count as usable for `max_depth="median"` (or a function): they enter the median, and
  the junctions they fill are no longer given the default depth (nor listed in `derived`).

Caveats: lids in fields, ravines, or easements are often raised above grade, and structures under
bridges, at creek crossings, or inside buildings cannot be read from the ground surface. A
bare-earth model may smooth sharp features; the window and the comparison are the guard. The DEM
and the City's elevations must share a vertical datum (CGVD28 and CGVD2013 differ by decimetres in
Ontario): check, or convert, before comparing.

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

## API reference

The API reference is generated from the docstrings with MkDocs (`mkdocs.yml`, `docs/`):

```bash
pip install -e ".[docs]"
mkdocs serve            # http://127.0.0.1:8000, rebuilt on save; or mkdocs build --strict
```

## CI

`.github/workflows/ci.yml` runs ruff (`ruff check` and `ruff format --check`, rules in
`pyproject.toml`) and the test suite on Python 3.11–3.13. The `CI gate` job is the single check
to require in branch protection.

## Releasing

1. In a pull request, set `version` in `pyproject.toml` and turn the CHANGELOG's `[Unreleased]`
   entries into the new version's section, dated, with its compare link at the bottom.
2. After it merges, publish a GitHub release whose tag is `v` plus the version (e.g. `v0.2.0`).
   `.github/workflows/publish.yml` then builds the tag and uploads it to PyPI, after checking that
   the tag matches `pyproject.toml`; PyPI never accepts the same version twice.

The upload uses PyPI's trusted publishing, so no token is stored: PyPI trusts this repository's
`publish.yml` in the `pypi` environment. To publish a tag again by hand, run the workflow from the
Actions tab and name the tag.
