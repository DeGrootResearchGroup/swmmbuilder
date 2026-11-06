"""network_footprint, footprint_distance, and network_extent over some of the pipes (#106)."""

import pytest
from shapely.geometry import Point
from test_pump_stations import make_network

from swmmbuilder.maps import footprint_distance, network_extent, network_footprint

# Force mains 2 (PS -> F1, x 50-350) and 3 (F1 -> MH, x 350-550) along y = 0; gravity pipes 1
# (U1 -> PS, x 0-50), 4 (MH -> OUT, x 550-600), and 5 (PS -> OVO, up to y = 50).
X0, Y0 = 480000, 4760000


@pytest.fixture(scope="module")
def sub():
    return make_network().subnetwork(["U1", "PS", "F1", "MH", "OUT", "OVO"])


def test_gravity_leaves_out_the_force_mains(sub):
    every = network_footprint(sub, distance=10)
    gravity = network_footprint(sub, distance=10, pipes="gravity")
    middle = Point(X0 + 300, Y0)  # on force main 2 only
    assert every.contains(middle) and not gravity.contains(middle)
    assert gravity.area < every.area


def test_distance_and_extent_follow_the_pipes(sub):
    assert network_extent(sub) == pytest.approx(600)
    assert network_extent(sub, pipes="gravity") == pytest.approx(600)  # 1 and 4 are 600 m apart
    assert network_extent(sub, pipes=[1, 5]) == pytest.approx(50)
    # Gravity pipes are all 50 m long: twice the median is 100 m, capped at 5 % of 600 m.
    assert footprint_distance(sub, pipes="gravity") == pytest.approx(30)
    assert footprint_distance(sub, pipes=[1, 5]) == pytest.approx(2.5)


def test_ids_and_a_function(sub):
    by_ids = network_footprint(sub, distance=5, pipes=["1", 4])
    by_function = network_footprint(sub, distance=5, pipes=lambda p: p.index.isin([1, 4]))
    assert by_ids.equals(by_function)
    assert not by_ids.contains(Point(X0 + 50, Y0 + 40))  # on pipe 5


def test_invalid(sub):
    with pytest.raises(ValueError, match="pipes must be 'gravity'"):
        network_footprint(sub, pipes="sanitary")
