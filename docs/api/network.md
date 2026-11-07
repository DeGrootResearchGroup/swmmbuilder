# Network and Subnetwork

A `Network` is the whole city's data; a `Subnetwork` (usually a catchment) is a set of its nodes
and the pipes between them. Most functions take either. Where they differ:

| | `Network` | `Subnetwork` |
| --- | --- | --- |
| Nodes | `nodes`: a GeoDataFrame of node data, indexed by ID | `nodes`: a list of node IDs (their data are `sub.network.nodes.loc[sub.nodes]`) |
| Pipes | `pipes()`, or `full_data` with the pipes left out of the graph (a missing node) | `pipes()` |
| Split IDs | `split_report`, an attribute (None without `split_tolerance`) | `split_report()`, a method adding `node` and `in_subnetwork`, and `split_summary()` |
| Rims from a DEM | `estimate_rims(dem, nodes=None)` | none: `net.estimate_rims(dem, nodes=sub.nodes)`, then `sub.to_inp(..., rims=est.rims)` |


::: swmmbuilder.network.Network

::: swmmbuilder.network.strip_leading_zeros

::: swmmbuilder.corrections.PipeCorrections

::: swmmbuilder.subnetwork.Subnetwork

::: swmmbuilder.subnetwork.TopologyReport

::: swmmbuilder.subnetwork.PumpStationReport

::: swmmbuilder.subnetwork.UpstreamTotals
