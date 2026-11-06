# Diagnostics

## Summaries

Each main result reduces to a flat dict of counts and totals, under the keys listed in its
docstring:
[`PopulationAllocation.summary`][swmmbuilder.population.PopulationAllocation.summary],
[`LandUseAllocation.summary`][swmmbuilder.land_use.LandUseAllocation.summary],
[`TopologyReport.summary`][swmmbuilder.subnetwork.TopologyReport.summary], and
[`SwmmModel.summary`][swmmbuilder.swmm.records.SwmmModel.summary].

## Warnings, errors, and the log

The warning categories, and `capture_log`, which collects the log and the warnings of one build
as text (also importable as `swmmbuilder.capture_log`).

::: swmmbuilder.diagnostics
