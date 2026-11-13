# Changelog

Notable changes to swmmbuilder. The format follows [Keep a Changelog](https://keepachangelog.com/),
and the project uses [semantic versioning](https://semver.org/) (before 1.0, a minor version may
change the API). Numbers in parentheses are pull requests.

0.1.0 is the first release. Its entries are the changes since the restructuring began (#24);
earlier history is in the git log.

## [Unreleased]

## [0.2.0] - 2026-11-13

### Added

- **Reused node IDs** (#97, #126, #143): `Network(split_tolerance=...)` splits a node ID used for
  structures far apart into one node per structure (`FH207a`, `FH207b`, ...), so that the network
  is not joined where the data are not. Pipe ends are grouped into structures by distance; a pipe
  drawn against its direction is matched by its other end and has its geometry reversed (listed in
  `Network.reversed_pipes`), so that node positions and .inp vertices run downstream; and each
  node-layer point goes to the nearest structure. `split_tolerance` must be larger than
  `snap_tolerance`. Every split is logged and listed in `Network.split_report` and
  `Network.splits`, and snapping never merges the structures of one ID again.
  `Subnetwork.split_report()` adds the node each new ID resolves to after snapping and whether it is
  in the subnetwork, and `Subnetwork.split_summary()` counts them.
- **Pipe corrections** (#145, #171): `Network(..., corrections=...)` applies documented
  corrections keyed by the source pipe ID (matched exactly or by text), after schema normalization
  and `normalize_ids` and before the invalid-pipe, flow-type, status, and overflow filters. `drop`
  removes pipes, `set` replaces canonical column values (e.g. `{'status': 'built'}` for a built
  pipe still marked proposed), and `swap_nodes` swaps `from_node` and `to_node`, with their node
  types, keeping the inverts and geometry as recorded. Each is logged at INFO with its reason
  (optional for `set`) and listed in `Network.corrections_report`; dropped pipes are counted in
  `Network.removed['corrections']['dropped']`. Corrections can be a mapping, a `PipeCorrections`,
  or a TOML file. Unknown pipe IDs or columns, a pipe given twice, and conflicting entries raise
  `ValueError`.
- **Implausible pipe inverts** (#116, #129): a schema's `invert_range` (opt-in, e.g.
  `invert_range=LONDON_ON.rim_range`, or the TOML key) treats pipe and drop-structure inverts
  outside it as missing, such as 25089.0 keyed for 250.89 or a drop height in the invert column.
  Each is logged, a `DataWarning` gives the count, and `Network.invert_report` lists those of the
  pipes kept, with their recorded values (inverts replaced by `corrections` are left out).
- **Rims from a ground-elevation model** (#132, #148, #168): `Network.estimate_rims(dem)` samples
  a lidar DEM (a GeoTIFF or a list of tiles, in any CRS) at the nodes, at the cell or as the
  median or mean within a window, and returns a `RimEstimate`. It compares the DEM with the
  recorded rims (`validation`: median offset and spread, and `largest_differences`) and estimates
  the rims that are missing or not above the invert, optionally after adding an offset such as
  the median difference; a DEM rim at or below the invert is not used. With `flag=...` it flags
  the recorded rims that deviate from the DEM by more than a threshold, measured from the median
  recorded rim − DEM so that a datum shift does not count, and `replace_above` and
  `replace_below` replace by the DEM estimate those recorded more than that far above or below
  the ground (separate thresholds, since raised lids are commoner than sunk ones). It does not
  change the network: the new `rims` fill option (`to_inp(rims=estimate.rims)`) fills them,
  listed in `data_gaps()` as filled from `estimate_rims`, with flagged rims as the new
  `implausible_rim` gap. A given rim at or below the invert is not used and gives a
  `DataWarning`. Filled rims count for `max_depth='median'` and `wet_well_defaults`.
  `swmmbuilder.elevation.sample_dem` reads rasters masked (nodata, mask band, mask file, or alpha
  band). `sample_dem` and `estimate_rims` take `z_factor` (metres per elevation unit); by default
  a DEM whose vertical CRS or band units are in feet is converted to metres, and units it does not
  know give a `DataWarning` and are taken as metres. Reading rasters needs the new `elevation`
  extra (rasterio).
- **Population allocation** (#96, #127, #149, #169): `allocate_population` places the population
  of polygons (census areas, or any other areal units) at network nodes. Each pipe goes to the
  polygon containing its midpoint (the first, for a midpoint on an edge two polygons share) or,
  with `max_distance`, to every polygon within that distance (for census blocks). A polygon
  flagged `sewered` (booleans, or 0 and 1) with no pipe goes to its nearest pipe, and force mains
  take no connections. A polygon's people are shared among its pipes by a weight
  (`swmmbuilder.population.Weight`: the length, `"diameter_length"`, a column such as a count of
  service connections, or a function returning one value per pipe). Invalid weights are imputed
  from the polygon's median weight per metre, with a `DataWarning`, and a fallback weight (by
  default one) is used where there is none to share by. Each share is placed at the pipe's
  upstream node, downstream node, or both. Given a `Network`, pipes with a missing `from_node` or
  `to_node` are left out, and a warning is logged listing them with the people they would have
  taken. `population` and `sewered` must be columns, or Series indexed by every polygon ID; a
  negative population and a network CRS not projected in metres are refused; and every invalid
  option is reported in one `ValueError`. The result reports the population of polygons with no
  pipe, and its `nodes` feed `dwf_from_population`.
- **Land-use flows** (#109, #127, #149, #169, #170): `allocate_land_use` places flows given per
  hectare of each land use (e.g. commercial, institutional, and industrial design allowances) on
  the nodes of a network. Pipes go to zoning polygons by the rules of `allocate_population`; each
  rated polygon's area, optionally clipped to a given area such as the network's footprint, times
  its rate is shared among its pipes by a weight as in `allocate_population` (by default the
  length, with invalid weights imputed and a fallback weight); and the flows can be scaled to a
  total, so that the rates are only relative weights. With `to_nearest`, a zone with no pipe of
  its own (a site on private sewers, say) sends its flow to its nearest pipe. A key of `rates`
  that is the land use of no polygon gives a `DataWarning`, and `land_use` must be a column or a
  Series indexed by every polygon ID. The result's `nodes` is a table of nodes by land use, and
  `by_pipe()` gives each pipe's share. `pattern_inflows` turns the table into `inflows` and
  `timeseries` with each land use's hourly pattern (a node with several land uses gets one
  flow-weighted series), plus an optional steady baseline.
- **Dominant category** (#142): `dominant_category(polygons, layer, column)` gives, for each
  polygon, the category of a layer (such as zoning) covering its largest area, and that area's
  share of the polygon. It is the usual source of `allocate_population`'s `sewered` flag, e.g.
  census blocks that are mostly zoned residential. A tie goes to the category that sorts first,
  and a polygon meeting no feature of the layer has a missing category. The polygons' CRS must be
  projected in metres.
- **Census polygons** (#130, #168): `swmmbuilder.census.census_polygons` gets the polygons
  covering a network (or any area) from a pluggable source, keeps those meeting the footprint of
  its gravity pipes, and returns them with `population` and `parent` columns for
  `allocate_population`. A GeoPackage cache, with the source, its parameters, and the fetch date,
  makes later builds offline. Sources: `CensusMapper` (Canadian census dissemination blocks or
  areas, through pycancensus, in the new `census` extra) and `PolygonFile` (any local file). Each
  has a `provenance()` method (for `CensusMapper`, the dataset, resolution, and regions or region
  level; for `PolygonFile`, the file name and columns). A cache read with a source is refused if
  it was fetched with other parameters, naming the one that differs, and gives a `DataWarning` if
  it does not record one. The CensusMapper API key is the user's to provide, and citing the data
  is the user's responsibility.
- **Spreading a flow over the pipes** (#110, #127, #169): `distribute_flow` shares a total flow
  (e.g. groundwater infiltration) among the pipes of a network or catchment by a weight as in
  `allocate_population` (by default the length), or gives each pipe a rate per unit of weight,
  and places it at the upstream nodes (or downstream, or half at each). Force mains take none
  unless asked. `max_weight` marks larger weights invalid, and a pipe with no valid weight takes
  no flow, with a `DataWarning`. `pipe_weights` gives the weight of each pipe. Given a `Network`,
  pipes with a missing node are left out.
- **External inflows and time series** (#108): the `inflows` and `timeseries` options of
  `to_inp` write [INFLOWS] and [TIMESERIES]. An `Inflow` is a time series times a scale plus a
  baseline, optionally times a pattern, so a node can carry a second time-varying flow (e.g.
  commercial flow with its own hourly pattern) or a steady one (e.g. infiltration) besides its
  [DWF] flow. `hourly_series` gives a series that holds each hourly value for the whole hour. The
  shipped template has the `{inflows}` and `{timeseries}` placeholders; a template needs them
  only when these options are used.
- **Node loading** (#146, #166, #171): `NodeLoading` combines the three parts of a dry-weather
  load: sanitary flow from population (`[DWF]` with one diurnal pattern), flows that each follow
  their own hourly pattern such as land uses (`[INFLOWS]` time series), and steady flows such as
  infiltration (the `[INFLOWS]` baseline). `options()` gives the matching `build_model` options,
  always with `flow_units`, so `sub.to_inp(path, outfalls, **loading.options())` writes them and
  refuses a template whose `FLOW_UNITS` differ. `flow_units` must be `LPS`, `CMS`, or `MLD` (any
  case). `table()` and `summary()` give the load of each node by part and the totals. Zero loads
  are not written, and a patterned kind cannot be named after a part (`population`, `sanitary`,
  `steady`, `patterned`). Its hourly patterns must have a mean within 0.001 of 1 unless
  `require_unit_mean=False` (or another tolerance is given). `pattern_inflows` takes
  `require_unit_mean`: by default a pattern whose mean is off 1 by more than 0.001 gives a
  `DesignWarning`, `True` or a tolerance refuses it, and `False` accepts it silently. The new
  `examples/loaded_catchment.py`, and a README section, build a loaded catchment model from the
  example pipes, from `Network` to `.inp`.
- **Design rules from the schema** (#150, #171): `build_model`, `write_inp`,
  `Subnetwork.to_inp`, and `Subnetwork.data_gaps` take `design="schema"`, which applies the
  schema's `force_main_discharge`, `force_main_roughness`, `force_main_checks`, and `default_pump`
  (those its `[design]` table holds), or a list of rule names, which applies just those
  (`wet_well_defaults` only by name, since it needs station data). An option passed explicitly
  wins, and the rules applied are logged. The default, None, applies none.
- **Finding pump stations** (#115): `Subnetwork.find_pump_stations()` returns a
  `PumpStationReport` listing the station nodes with a force main leaving them (to pass to
  `pump_stations`), those with none, and, as likely lift stations missing from the station layer,
  force mains that climb (by at least `min_rise`) from a node that is not a station.
- **Pump every station found** (#151, #171): `pump_stations` also takes `"found"` or a
  `PumpStationReport`, pumping every station of `Subnetwork.find_pump_stations()` with
  `default_pump` (given, or the schema's with `design="schema"`). `pump_stations_exclude`
  ({node: reason}) leaves out stations, stations with no force main, and unlisted climbing force
  mains; `pump_stations_extra` adds stations (a reason for the default pump, or a `PumpStation`
  used as given). A climbing force main left unpumped gives a `DataWarning`
  (`pump_stations_unlisted="raise"` or `"ignore"` to change that). The choices are logged and
  kept in `SwmmModel.station_choices`. `Subnetwork.data_gaps()` takes the same station options,
  so that its report is the one `to_inp` gives. An ideal pump must be its node's only outlet
  (SWMM's ERROR 134), so with an ideal default pump, such as London's, a station with an overflow
  pipe raises `ValueError`; give it a pump with a curve in `pump_stations_extra`, or exclude it.
- **Wet wells below the arriving sewers** (#114, #171): the opt-in design rule
  `station_invert='lowest_inlet'` lowers a pump station's invert to the lowest end of the pipes
  arriving at it, extrapolating a missing downstream invert at `headwater_slope`, so that a
  station whose invert would otherwise come from its force main alone does not sit above its
  inlet sewer. `wet_well_defaults` then takes the same estimated inlet end
  (`swmmbuilder.design.inlet_end`). Each change is listed in `SwmmModel.derived`, and the
  station's invert is reported in `data_gaps` as filled by this rule (`fill_gaps` takes
  `correction_rules`; `Subnetwork.data_gaps()` takes `station_invert`). `height_correction` wins
  over the rule.
- **Default junction depth from the recorded depths** (#113, #129): `max_depth='median'` gives
  every junction with no usable rim the median rim − invert of the junctions that have one (the
  depths written); a function of those depths also works. Each junction given it, including
  pump-station discharge nodes and outfalls written as junctions by `outfall_link`, is listed in
  `SwmmModel.derived`.
- **Linking outfalls under DYNWAVE** (#111, #171): with `outfall_link=True` (or a dict of
  `length`, `diameter`, `slope`, and `roughness`; an empty dict takes the defaults in
  `swmm.OUTFALL_LINK_DEFAULTS`), each outfall with more than one link in a dynamic-wave model is
  written as a junction draining through a dummy conduit `<node>_DUMMY` to a new outfall
  `<node>_OUT`. See Changed for the refusal without it.
- **Upstream totals** (#118, #124, #129): `Subnetwork.accumulate(values, split='equal')` gives the
  total of a node value (population, flow, load) upstream of every node and pipe, as an
  `UpstreamTotals` (`nodes`, `pipes`, `fractions`). A node's total is shared equally among its
  outlets (`'equal'`, conserving), counted once per upstream node (`'full'`), or split by given
  fractions `{pipe: fraction}`, e.g. from simulated flows; outlets not given share the rest
  equally. Fractions above 1 at a node raise `ValueError`, and where every outlet is given and
  they add up to less than 1, a `DataWarning` says the rest leaves the totals. Loops are solved
  exactly; with shared splits a loop is refused (`ValueError`) only if some of its nodes cannot
  reach an exit. A DataFrame of values is accumulated column by column.
- **Pipes crossing a subnetwork's boundary** (#164): `Subnetwork.boundary_pipes(values=None,
  ends=True)` lists the pipes with one end in the subnetwork, which it leaves out, classified as
  `force_main`, `split`, `summit`, or `inflow`, with the number of subnetwork nodes upstream (and
  the sum of `values` over them), where a leaving pipe's flow ends, and how far that is from the
  subnetwork's outlet. `TopologyReport.boundary` holds it. `Network.ends()` lists every node with
  no outlet pipe, with its node type and the number of nodes draining to it. See Changed for the
  new warnings of `Network.catchment`.
- **Summaries** (#144, #170, #171): `summary()` on `PopulationAllocation` and `LandUseAllocation`
  (each with an optional `sub`, to count one catchment), `NodeLoading`, `UpstreamTotals`,
  `PumpStationReport`, `RimEstimate`, `TopologyReport`, and `SwmmModel` returns a flat dict of the
  counts and totals a report cites, as plain numbers (NaN where undefined). The keys follow one
  pattern: a bare plural noun counts things, `total` is the amount placed at the nodes
  (`allocation_total` over the whole allocation), and `<part>_<measure>` measures one part, part
  first (`COMMERCIAL_flow`). The Diagnostics page of the docs lists every key. The population
  summary counts a pipe serving several polygons once, with the mean of its imputed weights; the
  land-use summary counts polygons whose pipes all weigh 0 under `zero_weight`; the model summary
  counts derived values by element, attribute, and rule, and pump stations (the new
  `SwmmModel.stations`); and the model and topology summaries count `data_gaps_blocking` and
  `data_gaps_no_rim`, and loops under `cycles`. `RimEstimate.report()` (also `str()`) gives its
  text. `TopologyReport` gains `nodes` and `pipes`, and `Network.pipes()` gives the pipes of the
  graph, as `Subnetwork.pipes()` does. The Network page of the docs compares `Network` and
  `Subnetwork`.
- **Capturing the log of a build** (#141): `swmmbuilder.capture_log(level="INFO", path=None)` is
  a context manager that collects the `swmmbuilder` log and the `SwmmBuilderWarning` warnings
  shown in its block as text (a `CapturedLog`, `captured.text`), and optionally writes it to a
  file. It attaches its own handler to the `swmmbuilder` logger and removes it on exit, so
  repeated builds in one process each get their log, unlike `logging.basicConfig`. Warnings are
  still shown as usual.
- **Footprint over some of the pipes** (#117): `maps.network_footprint`, `footprint_distance`,
  and `network_extent` take `pipes`: `'gravity'` (every pipe but the force mains, for a service
  area), a list of pipe IDs, or a function of the pipe rows. The default keeps every pipe.
- Releases are published to PyPI from GitHub releases (`publish.yml`, trusted publishing), and the
  README installs from PyPI (#95).

### Changed

- **Status and structure types for more Ontario presets** (#167): Barrie, Kitchener, Sarnia,
  Waterloo, and Windsor map their status columns, and Toronto and Waterloo their structure types.
  Windsor's `TEMP. OUT-NEW CONSTRUCTION` and `PRIVATE` pipes are now dropped by the default status
  filter, and Toronto's overflow pipes (`OVF`) by `include_overflows=False`. Rarer codes
  (Toronto's `SP`, `PF`, ..., Waterloo's `WWTP Outfall`) pass through unmapped. Each TOML header
  says which columns a dataset lacks.
- **One cycle per loop** (#165): `TopologyReport.cycles` holds one cycle per loop (the shortest
  through the loop's smallest node) instead of every simple cycle. The new `TopologyReport.loops`
  lists the nodes of each loop, a set of nodes that can all reach one another. A network with no
  cycles reports the same as before.
- **Data-gaps reports list nodes with no usable rim** (#132): `data_gaps()` and `model.data_gaps`
  have a `no_rim` row (not blocking), and `TopologyReport.data_gaps` a `no_rim` count, for every
  node with no rim recorded, or one at or below the invert, so existing reports gain rows. Unless
  the `rims` option fills it, such a junction takes the default `max_depth`, as before.
- **Outfalls with more than one link under DYNWAVE** (#111): `write_inp` now refuses, before
  writing, a dynamic-wave model with an outfall that has more than one link (SWMM's ERROR 141),
  naming the outfalls and their links. `outfall_link` (see Added) writes such outfalls so that
  SWMM can run them.
- **Templates and routing** (#166): `check_template`, and so `write_inp`, refuses a template
  whose `LINK_OFFSETS` is not `DEPTH`, since offsets are written as depths. Force mains written
  under routing other than `DYNWAVE` give an `OptionWarning`, since SWMM then treats them as open
  channels.
- **[DWF] pattern mean** (#171): an HOURLY or WEEKEND [DWF] pattern gives a `DesignWarning` when
  its mean is off 1 by more than 0.001 (it was 0.01), the tolerance of `pattern_inflows` and
  `NodeLoading`.
- **Catchment warnings** (#164): `Network.catchment` gives a `DataWarning` about splits and force
  mains leaving the catchment, and about other networks ending within `near` (a new parameter,
  default 100 m) of its outlet, a likely second inlet.

### Fixed

- **A missing invert no longer blocks a node merge** (#173): `snap_nodes` with
  `max_invert_difference` compares only known inverts, for each pair and for each merged group.
  A missing invert made the comparison fail, so two IDs at one point were never merged when one
  had no invert, and the report gave the reason as an invert difference. Those pairs are now
  merged on distance and the pipe rule, which changes the network for data with missing inverts.
  The report gains `invert_known` (both inverts present).
- **Topology checks on looped areas** (#165): `Subnetwork.validate()` finds cycles as strongly
  connected components rather than listing every simple cycle, whose number grows exponentially
  with the size of a looped area (a 7 × 7 grid of pipes drawn both ways did not finish). It now
  takes linear time.
- **Pipes with no ID** (#125): two or more pipes with a missing ID (or one of the schema's
  `null_values`) no longer stop the data from loading as duplicates; they are dropped and counted
  in `Network.removed['invalid']['no pipe_id']`, as one such pipe already was.
- **Disjoint MultiLineStrings** (#126): `line_coords` keeps the parts of a pipe whose parts do not
  join in their order, as documented; `shapely.line_merge` could reorder them. A MultiLineString
  is merged only if it merges into one line.
- **Small flows in the .inp file** (#129): dry-weather, pump-curve, and conduit flows, and the new
  inflows and time series, are written to 10 significant figures (`swmm.FLOW_DIGITS`) instead of
  6 decimals, so a flow below a micro-unit (common in CMS) is no longer rounded off or written
  as 0.
- **`Subnetwork.in_network`** builds its index of node IDs as text once, not on every call, which
  made `height_correction` slow on large networks (#129).

## [0.1.0] - 2026-11-03

### Added

- **Pumping.** Pumps and wet wells in the .inp writer (#30); pump stations whose force mains are
  kept as conduits (#41); opt-in `wet_well_defaults` and `force_main_checks` design rules (#42).
- **Dry-weather flow** loading: `dwf`, `patterns`, `dwf_patterns` and `dwf_from_population` (#40).
- **Data gaps**: `Subnetwork.data_gaps()` reports missing sizes, shape codes and inverts, and opt-in
  fill rules fill them (#39).
- **Node layers** (maintenance holes, pump stations, outfalls, other nodes), pipe status and
  structure type, and snapping of coincident node IDs (#37).
- **Schemas** for Barrie, Kitchener, Peterborough, Sarnia, Waterloo and Windsor (#38).
- **Maps**: pipes coloured and sized by value, ribbons, group polygons, basemaps and Leaflet export
  (#36); `maps.ColourScale` to share one colour scale between plots and Leaflet layers (#87).
- **Warning categories** `DataWarning`, `DesignWarning` and `OptionWarning` (all
  `SwmmBuilderWarning`, a `UserWarning`), and a written rule for raise / warn / log (#89).
- **Type hints** on the public API, with `NodeId` / `PipeId` aliases and a `py.typed` marker, checked
  by mypy in CI (#90).
- `Subnetwork.match_pipe` / `match_pipes`, `Network.check_nodes`, `Network.ids` and
  `Network.lowest_inverts` (#88); `Network.version` (#83).
- A misspelt option of `to_inp`, `write_inp`, `build_model` or `Subnetwork.data_gaps` raises a
  `TypeError` naming the closest options (e.g. `size_overide` -> `size_overrides`) (#91).
- The `simulate` extra (pyswmm) (#66).
- An **API reference** built with MkDocs from the docstrings (`docs` extra), built strictly in CI
  (#92).

### Changed

- The package is `swmmbuilder` (formerly SewerNetworkx) (#35). `Network.catchment` returns a
  `Subnetwork` (#27); models are built as records and rendered once (#25); progress is logged, not
  printed (#24).
- **Package layout.** `swmmbuilder.swmm` and `swmmbuilder.maps` are packages, with the design rules,
  pumps and units in their own modules (#84, #85), and `build_model` runs in phases (#86). Every
  former name is re-exported, so imports keep working.
- **pyswmm and gmplot are optional**: install `swmmbuilder[simulate]` to run models and
  `swmmbuilder[maps]` for `generate_gmap` (#66).
- **.inp output**: fields are always separated, numbers are rounded to 6 decimals in their shortest
  form, and names that SWMM cannot read (empty, whitespace or `;`) raise `ValueError` (#67).
- **Node and pipe IDs** given to options match the data's IDs exactly or by their text (`2151735`
  and `'2151735'`), and through `normalize_ids` and snap aliases (#68, #71, #88). An unknown pipe
  in `pumps` is reported as `pumps: pipes not in the subnetwork: ...`, and two keys naming one ID
  with different values raise `ValueError` (#88).
- **Loading**: a pipe with no geometry but a length is kept; pipes with no ID or no length are
  dropped and counted in `Network.removed['invalid']` (#70).
- **Snapping** applies its pipe and invert rules to the whole merged group, so groups of three or
  more nodes can come out smaller (#69). A `Subnetwork` made before a snap that merges nodes raises
  `RuntimeError` when used (#83).
- **Maps** always draw in a CRS projected in metres (#72); Leaflet markers for nodes with no position
  raise like `plot_markers` (#74).
- **Diagnostics**: a rim not above the invert, `height_correction` keys naming no node and plot
  values with no known ID are now warnings, not log messages (#89).
- The schema presets other than `LONDON_ON` are read on first use, not on import (#91).

### Fixed

- Long names or numbers running into the next .inp field (#67).
- Numeric node IDs (Kitchener) breaking `catchment`, node lists and error messages (#68).
- `snap_nodes` turning a pipe into a self-loop through a chain of merges (#69).
- Crashes on missing geometry, NaN pipe IDs, numeric node-layer columns and zero-length invert
  interpolation; multi-part pipes drawn with a bridging segment (#70).
- Maps: scale bar and distances in feet labelled as metres (#72); log norms, all-NaN, all-equal and
  all-zero values (#73); reusing an axes after `ax.clear()`, list colours shifted by empty groups,
  quadratic ribbons, `to_leaflet` with no pipes, and `to_geodataframe('nodes')` returning the
  network's own table (#74).
- Pump-station review findings (#56).
- Regular (non-editable) installs left out the subpackages (#84).

### Deprecated

- `maps.generate_gmap`; use `maps.to_leaflet` (#36).

[Unreleased]: https://github.com/DeGrootResearchGroup/swmmbuilder/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/DeGrootResearchGroup/swmmbuilder/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/DeGrootResearchGroup/swmmbuilder/releases/tag/v0.1.0
