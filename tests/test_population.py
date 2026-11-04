"""Allocation of polygon populations to network nodes (swmmbuilder.population)."""

import math

import geopandas as gpd
import pandas as pd
import pytest
from conftest import FIXTURE
from shapely.geometry import LineString, box

from swmmbuilder import DataWarning, Network, allocate_population, dwf_from_population

CRS = "EPSG:2958"
X0, Y0 = 480_000.0, 4_760_000.0

# A straight line of pipes along y = 50 (x in m from X0), as (pipe, up, down, x_up, x_down,
# PDC_Count). Their lengths are 50, 30, 60, 40 and 60 m; midpoints at x = 35, 75, 120, 170
# and 220.
PIPES = [
    (1, "N1", "N2", 10, 60, 2),
    (2, "N2", "N3", 60, 90, 6),
    (3, "N3", "N4", 90, 150, None),  # missing: imputed
    (4, "N4", "N5", 150, 190, 2),
    (5, "N5", "N6", 190, 250, 999),  # implausible: imputed; midpoint in no polygon
]
# P holds pipes 1 and 2, Q pipes 3 and 4, R none.
POLYGONS = {"P": (0, 100, 80.0), "Q": (100, 200, 100.0), "R": (300, 400, 30.0)}


def make_network(pipes=PIPES):
    template = gpd.read_file(FIXTURE).iloc[[0]].drop(columns="geometry")
    rows = []
    for pipe, up, down, x_up, x_down, pdc in pipes:
        row = template.iloc[0].to_dict()
        row.update(
            OBJECTID=pipe,
            UpstreamID=up,
            DownstreamID=down,
            Length=float(x_down - x_up),
            PDC_Count=pdc,
            geometry=LineString([(X0 + x_up, Y0 + 50), (X0 + x_down, Y0 + 50)]),
        )
        rows.append(row)
    return Network(gpd.GeoDataFrame(rows, crs=CRS))


def squares(spec=POLYGONS):
    """Polygons 0-100 m in y, x from x0 to x1, indexed by name, with a population column."""
    return gpd.GeoDataFrame(
        {"population": [pop for _, _, pop in spec.values()]},
        geometry=[box(X0 + x0, Y0, X0 + x1, Y0 + 100) for x0, x1, _ in spec.values()],
        index=pd.Index(list(spec), name="id"),
        crs=CRS,
    )


@pytest.fixture
def net():
    return make_network()


def allocate(net, polygons=None, **options):
    options = {"weight": "PDC_Count", "max_weight": 100, **options}
    return allocate_population(net, squares() if polygons is None else polygons, **options)


def approx(d):
    return pytest.approx(d, abs=1e-9)


def test_shared_by_weight_at_upstream_nodes(net):
    # Q's only valid weight per metre is pipe 4's, 2/40, so pipe 3 weighs 0.05 x 60 = 3. Pipe 5
    # is in no polygon, so has no rate: weight 0. P (weights 2 and 6) and Q (3 and 2).
    alloc = allocate(net)
    assert alloc.nodes == approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})
    assert alloc.pipes["polygon"].tolist()[:4] == ["P", "P", "Q", "Q"]
    assert pd.isna(alloc.pipes.loc[5, "polygon"])
    assert alloc.pipes["weight"].tolist() == approx([2, 6, 3, 2, 0])
    assert alloc.pipes["weight_imputed"].tolist() == [False, False, True, False, True]
    assert alloc.total == pytest.approx(180)
    assert alloc.unallocated == pytest.approx(30)
    assert alloc.no_pipes == ["R"]
    assert alloc.zero_weight == []
    assert alloc.polygons.loc["Q"].to_dict() == approx(
        {"population": 100, "pipes": 2, "weight": 5, "allocated": 100}
    )
    assert alloc.polygons.loc["R", "allocated"] == 0


def test_feeds_dwf_from_population(net):
    dwf = dwf_from_population(allocate(net).nodes, per_capita=86_400, unit="L/d")
    assert dwf == approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})  # L/s


@pytest.mark.parametrize(
    "placement, expected",
    [
        ("downstream", {"N2": 20, "N3": 60, "N4": 60, "N5": 40, "N6": 0}),
        ("split", {"N1": 10, "N2": 40, "N3": 60, "N4": 50, "N5": 20, "N6": 0}),
    ],
)
def test_placement(net, placement, expected):
    assert allocate(net, placement=placement).nodes == approx(expected)


def test_length_weights_by_default(net):
    alloc = allocate_population(net, squares())
    assert alloc.nodes == approx({"N1": 50, "N2": 30, "N3": 60, "N4": 40, "N5": 0})
    assert not alloc.pipes["weight_imputed"].any()


def test_callable_weight(net):
    alloc = allocate(net, weight=lambda pipes: 10 * pipes["PDC_Count"], max_weight=1000)
    assert alloc.nodes == approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})


def test_imputed_from_the_pipes_own_polygon():
    # P: pipe 1 has 3 / 30 per metre, so pipe 2 (20 m, missing) weighs 2; the network median
    # (of 0.1 and Q's 0.05) would give it 1.5.
    pipes = [
        (1, "N1", "N2", 10, 40, 3),
        (2, "N2", "N3", 40, 60, None),
        (3, "N4", "N5", 110, 150, 2),
    ]
    alloc = allocate(make_network(pipes))
    assert alloc.pipes["weight"].tolist() == approx([3, 2, 2])
    assert alloc.nodes == approx({"N1": 48, "N2": 32, "N4": 100})


def test_polygon_with_no_valid_weight_is_reported():
    # Q's only pipe has no count: nothing to impute from, so Q's people are not allocated.
    pipes = [(1, "N1", "N2", 10, 60, 2), (3, "N3", "N4", 110, 150, None)]
    alloc = allocate(make_network(pipes))
    assert alloc.no_valid_weight == ["Q"]
    assert alloc.zero_weight == []
    assert alloc.pipes.loc[3, "weight"] == 0
    assert alloc.nodes == approx({"N1": 80, "N3": 0})
    assert alloc.unallocated == pytest.approx(130)


def test_impute_zero(net):
    alloc = allocate(net, impute="zero")
    assert alloc.pipes["weight"].tolist() == approx([2, 6, 0, 2, 0])
    assert alloc.nodes == approx({"N1": 20, "N2": 60, "N3": 0, "N4": 100, "N5": 0})


def test_partly_sewered_and_zero_weight_polygons():
    # S holds only pipe 4, whose weight is 0; Q is only partly sewered (its east half has no pipe).
    pipes = [(1, "N1", "N2", 110, 140, 3), (4, "N4", "N5", 210, 240, 0)]
    alloc = allocate(make_network(pipes), squares({"Q": (100, 200, 50.0), "S": (200, 300, 7.0)}))
    assert alloc.nodes == approx({"N1": 50, "N4": 0})
    assert alloc.zero_weight == ["S"]
    assert alloc.unallocated == pytest.approx(7)


def test_missing_population_counts_as_zero(net):
    polygons = squares()
    polygons.loc["P", "population"] = None
    alloc = allocate(net, polygons, population=polygons["population"].astype(object))
    assert alloc.nodes["N1"] == 0 and alloc.polygons.loc["P", "population"] == 0
    assert alloc.no_pipes == ["R"]


def test_within_counts_only_the_subnetworks_pipes(net):
    alloc = allocate(net)
    sub = net.catchment("N4")  # pipes 1-3
    assert alloc.within(sub) == approx({"N1": 20, "N2": 60, "N3": 60})
    # Allocating over the subnetwork alone (here by length) puts all of Q's people on pipe 3.
    alone = allocate_population(net.catchment("N4"), squares())
    assert alone.nodes == approx({"N1": 50, "N2": 30, "N3": 100})


def test_overlapping_polygons_take_the_first(net):
    polygons = squares({"P": (0, 100, 80.0), "Big": (0, 200, 10.0)})
    with pytest.warns(DataWarning, match="2 pipe midpoints are in more than one polygon"):
        alloc = allocate(net, polygons)
    assert alloc.pipes["polygon"].tolist()[:4] == ["P", "P", "Big", "Big"]


def test_reprojects_polygons(net):
    alloc = allocate(net, squares().to_crs("EPSG:4326"))
    assert alloc.nodes == pytest.approx({"N1": 20, "N2": 60, "N3": 60, "N4": 40, "N5": 0})


def test_invalid_options_listed_together(net):
    with pytest.raises(ValueError) as error:
        allocate(net, weight="Connections", placement="middle", impute="mean")
    message = str(error.value)
    for name in ("weight", "placement", "impute"):
        assert name in message


def test_repeated_ids_and_missing_crs(net):
    polygons = squares()
    polygons.index = ["P", "P", "R"]
    polygons = polygons.set_crs(None, allow_override=True)
    with pytest.raises(ValueError, match="repeated") as error:
        allocate(net, polygons)
    assert "polygons has no CRS" in str(error.value)


def test_negative_population(net):
    polygons = squares()
    polygons.loc["Q", "population"] = -1
    with pytest.raises(ValueError, match=r"Negative population in polygons \['Q'\]"):
        allocate(net, polygons)


def test_shares_sum_to_each_polygon(net):
    alloc = allocate(net, placement="split")
    sewered = alloc.polygons[alloc.polygons["pipes"] > 0]
    assert sewered["allocated"].tolist() == approx(sewered["population"].tolist())
    assert math.isclose(alloc.total + alloc.unallocated, squares()["population"].sum())
