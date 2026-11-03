# swmmbuilder

swmmbuilder loads a city's GIS sewer network (GeoJSON, or any format geopandas reads), selects
subnetworks of it, and writes them as SWMM `.inp` files.

```python
from swmmbuilder import Network

net = Network("london.geojson")                # schema=LONDON_ON by default
sub = net.catchment("OUT")                     # every node draining to OUT
sub.to_inp("model.inp", outfalls=["OUT"])
```

The [README](https://github.com/DeGrootResearchGroup/swmmbuilder#readme) is the user guide: node
layers, schemas, pumping stations, dry-weather flow, data gaps, maps and diagnostics. This site is
the API reference, generated from the docstrings:

- [Network and Subnetwork](api/network.md): loading, selecting and validating subnetworks;
- [Schemas](api/schema.md): how a city's GIS columns map to the canonical ones, and the presets;
- [SWMM model](api/swmm.md): `build_model`, `write_inp` and the model records;
- [Pumps and design rules](api/pumps.md);
- [Data gaps](api/gaps.md): reporting and filling missing sizes and inverts;
- [Diagnostics](api/diagnostics.md): how problems are reported, and the warning categories;
- [Maps](api/maps.md): static plots and Leaflet maps.
