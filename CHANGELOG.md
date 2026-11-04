# Changelog

Notable changes to swmmbuilder. The format follows [Keep a Changelog](https://keepachangelog.com/),
and the project uses [semantic versioning](https://semver.org/) (before 1.0, a minor version may
change the API). Numbers in parentheses are pull requests.

0.1.0 is the first release. Its entries are the changes since the restructuring began (#24);
earlier history is in the git log.

## [Unreleased]

### Added

- Releases are published to PyPI from GitHub releases (`publish.yml`, trusted publishing), and the
  README installs from PyPI.
- **Population allocation**: `allocate_population` places the population of polygons (census
  areas, or any other areal units) at network nodes. Each pipe goes to the polygon containing its
  midpoint or, with `max_distance`, to every polygon within that distance (for census blocks). A
  polygon's people are shared among its pipes by a weight (a column such as a count of service
  connections, a function, or the length). Invalid weights are imputed from the polygon's median
  weight per metre, with a fallback weight (by default one) where there is none to share by. Each
  share is placed at the pipe's upstream node, downstream node, or both. The result reports the
  population of polygons with no pipe, and its `nodes` feed `dwf_from_population`.

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

[Unreleased]: https://github.com/DeGrootResearchGroup/swmmbuilder/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/DeGrootResearchGroup/swmmbuilder/releases/tag/v0.1.0
