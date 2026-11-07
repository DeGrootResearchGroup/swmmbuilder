# Diagnostics

## Summaries

Every main result has a `summary()`: a flat dict of plain Python numbers (`int` for counts,
`float` for amounts, NaN where there is nothing to measure), to cite in a report or a paper.

### Key rules

The keys of every `summary()` follow the same rules.

- **A bare plural noun counts things** over the whole result: `nodes`, `pipes`, `polygons`,
  `stations`. `nodes` counts the nodes given an amount above 0.
- **`total` is the amount placed at the nodes**, in the result's own quantity: people for a
  population, flow (in the result's `flow_units`) otherwise. With a subnetwork `sub`, `total`
  and `nodes` are within it, and `allocation_total` is the whole allocation's.
- **`<part>_<measure>` measures one part**, part first: the part is a land use, a patterned kind,
  a part of a load (`sanitary`, `patterned`, `steady`), or a group of polygons (`nearest`,
  `no_pipes`); the measure is `nodes`, `polygons`, `area_ha`, `flow`, `population`,
  `design_flow`, or a statistic (`max`, `median_distance`, `max_distance`). So
  `COMMERCIAL_flow` and `COMMERCIAL_nodes` in both the land-use allocation and the node loading.
- **Units** are those of the result: areas in ha (`_ha`), distances and elevations in m, flows
  in `flow_units`.

A key renamed to fit these rules before the first release still reads, with a
`DeprecationWarning` naming its new key, but is not listed.

### Every key

[`PopulationAllocation.summary(sub=None)`][swmmbuilder.population.PopulationAllocation.summary]
(keys marked * are within `sub` when it is given):

| Key | Meaning |
| --- | --- |
| `total`* | people placed at the nodes (formerly `population`) |
| `nodes`* | nodes given people |
| `weight`* | total weight of the pipes, each pipe once |
| `allocation_total` | people placed over the whole allocation (formerly `population_total`) |
| `unallocated` | people placed at no node |
| `polygons`, `polygons_served` | polygons given, and those with a pipe |
| `nearest`, `nearest_population`, `nearest_max_distance` | polygons given their nearest pipe, their people, and the farthest distance (m) |
| `no_pipes`, `no_pipes_population` | polygons with people and no pipe, and their people |
| `no_valid_weight`, `zero_weight` | polygons served but left unallocated, with no valid weight or pipes that all weigh 0 |
| `fallback` | polygons some of whose pipes took the fallback weight |

[`LandUseAllocation.summary(sub=None)`][swmmbuilder.land_use.LandUseAllocation.summary]
(keys marked * are within `sub` when it is given):

| Key | Meaning |
| --- | --- |
| `<use>_polygons`, `<use>_area_ha` | the land use's polygons served, and their area counted (ha) |
| `<use>_flow`*, `<use>_nodes`* | its flow placed, and the nodes given it |
| `total`*, `nodes`*, `weight`* | flow placed at the nodes, nodes given flow, and the weight of the pipes |
| `allocation_total` | flow placed over the whole allocation |
| `polygons`, `polygons_served`, `area_ha` | polygons with a rate, those with a pipe, and their area (ha) |
| `nearest`, `nearest_median_distance`, `nearest_max_distance` | polygons given their nearest pipe, and their distances (m) |
| `no_pipes`, `zero_weight` | polygons with a flow and no pipe, or with pipes that all weigh 0 (counted as served) |
| `fallback` | polygons some of whose pipes took the fallback weight |
| `unallocated_design_flow` | the design flow of both, placed nowhere and not scaled (formerly `unallocated`) |
| `scale` | the factor applied to every flow to reach a given total |

[`NodeLoading.summary()`][swmmbuilder.swmm.loading.NodeLoading.summary]:

| Key | Meaning |
| --- | --- |
| `total`, `nodes` | flow of the three parts, and the nodes given any (formerly `loaded_nodes`) |
| `population`, `population_nodes` | people, and the nodes with any |
| `sanitary_flow`, `sanitary_nodes` | the [DWF] base flow (formerly `sanitary`) |
| `patterned_flow`, `patterned_nodes` | the patterned flows, all kinds (formerly `patterned`) |
| `steady_flow`, `steady_nodes` | the steady flows (formerly `steady`) |
| `<kind>_flow`, `<kind>_nodes` | each patterned kind (formerly `patterned_<kind>` and `patterned_<kind>_nodes`) |

[`UpstreamTotals.summary()`][swmmbuilder.subnetwork.UpstreamTotals.summary]:

| Key | Meaning |
| --- | --- |
| `nodes`, `pipes` | nodes and pipes with a total |
| `shared_pipes` | pipes carrying less than the whole total of their upstream node |
| `max` (or `<column>_max`) | the largest total at a node (per column, for values given as a DataFrame) |

[`PumpStationReport.summary()`][swmmbuilder.subnetwork.PumpStationReport.summary]:

| Key | Meaning |
| --- | --- |
| `stations` | stations with a force main |
| `no_force_main` | stations with none |
| `unlisted` | likely lift stations missing from the station layer, one per force main |

[`RimEstimate.summary()`][swmmbuilder.elevation.RimEstimate.summary] (`report()` and `str()`
give the same as text):

| Key | Meaning |
| --- | --- |
| `nodes` | nodes sampled |
| `compared`, `implausible`, `replaced` | recorded rims compared with the DEM, flagged as implausible, and replaced |
| `estimated`, `no_dem_value`, `dem_not_above_invert` | nodes with no usable rim given an estimate, nodes with no DEM value, and estimates not above the invert |
| `offset` | added to the DEM values (m) |
| `difference_<statistic>` | recorded rim − DEM (m): `median`, `mean`, `std`, `iqr`, `p05`, `p95`, `max_abs` |

[`TopologyReport.summary()`][swmmbuilder.subnetwork.TopologyReport.summary]:

| Key | Meaning |
| --- | --- |
| `nodes`, `pipes` | nodes and pipes of the subnetwork |
| `splits`, `cycles`, `dead_ends`, `parallel_pipes` | nodes with more than one outlet, directed cycles, dead ends, and node pairs joined by more than one pipe |
| `force_mains`, `stations_without_force_main`, `orphan_force_mains`, `force_main_ends` | force mains, and their problems |
| `pipe_count_mismatches` | nodes whose recorded pipe counts differ from the graph |
| `data_gaps` | gaps in the data, all kinds |
| `boundary_outflows`, `boundary_inflows`, `boundary_splits` | pipes leaving and entering the subnetwork, and the leaving ones whose flow the model sends elsewhere |

[`SwmmModel.summary()`][swmmbuilder.swmm.records.SwmmModel.summary]:

| Key | Meaning |
| --- | --- |
| `junctions`, `outfalls`, `storage`, `conduits`, `pumps` | records written |
| `stations`, `stations_extra`, `stations_excluded` | pump stations, and those added and left out |
| `outfall_links` | outfalls added by `outfall_link` |
| `data_gaps`, `data_gaps_filled` | gaps found, and those filled |
| `derived`, `derived_<element>_<attribute>`, `derived_<element>_<attribute>[<rule>]` | values derived by a design rule, in all, by element and attribute, and by rule |
| `default_depth_junctions`, `default_depth` | junctions given the default depth, and that depth (m) |

`Subnetwork.split_summary()` counts the reused IDs split (`ids_split`, `nodes_made`,
`ids_in_subnetwork`).

## Warnings, errors, and the log

The warning categories, and `capture_log`, which collects the log and the warnings of one build
as text (also importable as `swmmbuilder.capture_log`).

::: swmmbuilder.diagnostics
