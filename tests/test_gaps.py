"""Reporting, and optionally filling, gaps in pipe sizes, shape codes and inverts (#32)."""

import geopandas as gpd
import pandas as pd
import pytest
from conftest import DENDRITIC, FIXTURE, read_sections, run_swmm, write_inp
from shapely.geometry import LineString

from swmmbuilder import LONDON_ON, Network, Schema
from swmmbuilder.gaps import force_main_discharge_rule

X0, Y0 = -81.30, 42.98
# Fixture positions (see tests/data/make_fixture.py) and those of the nodes added here.
POS = {
    "A1": (0, 4),
    "A2": (1, 4),
    "J1": (2, 3),
    "H0": (-3, 7),
    "H1": (-2.5, 7),
    "H2": (-2, 6.5),
    "H3": (-1.5, 6),
    "H4": (-1, 5),
    "Q1": (1, 5),
    "P1": (3, 4),
    "P2": (2.5, 3.5),
}

# Pipes added to the fixture: (pipe_id, upstream, downstream, us_invert, ds_invert, length,
# structure type). Zero inverts are missing (as in the London data).
ADDED = [
    # H1 -> H2 -> H3 -> H4: three pipes with zero inverts between the known H1 and H4.
    (3000, "H0", "H1", 106.50, 106.20, 40.0, "SL"),
    (3001, "H1", "H2", 0.0, 0.0, 30.0, "SL"),
    (3002, "H2", "H3", 0.0, 0.0, 30.0, "SL"),
    (3003, "H3", "H4", 0.0, 0.0, 40.0, "SL"),
    (3004, "H4", "A1", 105.40, 105.10, 40.0, "SL"),
    # A headwater with no invert.
    (3005, "Q1", "A2", 0.0, 0.0, 50.0, "SL"),
    # A two-pipe force main from P1, with no inverts past P1, into the maintenance hole J1.
    (3006, "P1", "P2", 101.00, 0.0, 100.0, "FM"),
    (3007, "P2", "J1", 0.0, 0.0, 100.0, "FM"),
]
GAPS = DENDRITIC + [(up, down) for _, up, down, *_ in ADDED]

# A1 -> A2 (1000) is round with height = width = 200 mm but no diameter; A2 -> J1 (1001) has
# no size at all; E1 -> J2 (1012) has no shape code but a 200 mm diameter.
NO_DIAMETER, NO_SIZE, NO_SHAPE = 1000, 1001, 1012

FILLS = dict(
    shape_codes={"": "CIRCULAR"},
    size_fill="downstream",
    invert_fill="interpolate",
    headwater_slope=0.01,
    force_main_discharge=LONDON_ON.design["force_main_discharge"],
)


def point(node):
    x, y = POS[node]
    return X0 + 1e-3 * x, Y0 + 1e-3 * y


@pytest.fixture
def network():
    data = gpd.read_file(FIXTURE)
    data.loc[data["OBJECTID"] == NO_DIAMETER, "Diameter"] = None
    data.loc[data["OBJECTID"] == NO_SIZE, ["Diameter", "Height", "Width"]] = None
    template = data[data["OBJECTID"] == 1000].iloc[0]
    added = [
        template.to_dict()
        | dict(
            OBJECTID=pipe_id,
            UpstreamID=up,
            DownstreamID=down,
            UpstreamInvert=us,
            DownstreamInvert=ds,
            Length=length,
            Diameter=200,
            StructureType=structure,
            geometry=LineString([point(up), point(down)]),
        )
        for pipe_id, up, down, us, ds, length, structure in ADDED
    ]
    data = pd.concat([data, gpd.GeoDataFrame(added, crs=data.crs)], ignore_index=True)
    return Network(gpd.GeoDataFrame(data, crs="EPSG:4326"))


@pytest.fixture
def sub(network):
    return network.subnetwork({n for edge in GAPS + [("E1", "J2")] for n in edge})


def gaps_of(report):
    """{(element, id, gap): filled_from or None}"""
    return {
        (r.element, r.id, r.gap): None if pd.isna(r.filled_from) else r.filled_from
        for r in report.itertuples()
    }


def test_report_lists_every_gap_and_fills_nothing_by_default(sub):
    gaps = gaps_of(sub.data_gaps())
    blocking = {key: rule for key, rule in gaps.items() if key[2] != "no_us_invert"}
    assert {k: v for k, v in blocking.items() if k[2] != "no_ds_invert"} == {
        ("pipe", NO_DIAMETER, "no_size"): "height_width",  # the one unambiguous inference
        ("pipe", NO_SIZE, "no_size"): None,
        ("pipe", NO_SHAPE, "unknown_shape"): None,
        ("node", "H2", "no_node_invert"): None,
        ("node", "H3", "no_node_invert"): None,
        ("node", "Q1", "no_node_invert"): None,
        ("node", "P2", "no_node_invert"): None,
    }
    ends = {k for k, v in gaps.items() if k[2] in ("no_us_invert", "no_ds_invert") and v is None}
    assert ("pipe", 3007, "no_ds_invert") in ends and ("pipe", 3001, "no_us_invert") in ends
    assert ("pipe", 3006, "no_us_invert") not in ends


def test_report_details(sub):
    report = sub.data_gaps().set_index(["element", "id", "gap"])
    assert report.loc[("pipe", NO_SHAPE, "unknown_shape"), "detail"] == (
        "no cross-section: shape code ''"
    )
    assert "CIRCULAR needs diameter" in report.loc[("pipe", NO_SIZE, "no_size"), "detail"]
    assert report.loc[("pipe", 3007, "no_ds_invert"), "detail"] == "at node J1 (force main)"
    assert report.loc[("pipe", NO_DIAMETER, "no_size"), "value"] == ("CIRCULAR", 0.2, 0.0)


def test_validate_counts_gaps(sub):
    counts = sub.validate(["OUT"]).data_gaps
    assert counts["no_size"] == 2 and counts["unknown_shape"] == 1
    assert counts["no_node_invert"] == 4


def test_to_inp_lists_every_remaining_gap(sub, tmp_path):
    with pytest.raises(ValueError, match="6 gaps in the data remain") as error:
        write_inp(sub, tmp_path, routing="DYNWAVE")
    message = str(error.value)
    assert "Subnetwork.data_gaps()" in message
    for item in (f"pipe {NO_SIZE}", f"pipe {NO_SHAPE}", "node H2", "node H3", "node Q1", "P2"):
        assert item in message
    assert f"pipe {NO_DIAMETER}" not in message
    assert not (tmp_path / "model.inp").exists()


def test_shape_codes_fill(sub):
    report = sub.data_gaps(shape_codes={"": "CIRCULAR"}).set_index(["element", "id", "gap"])
    row = report.loc[("pipe", NO_SHAPE, "unknown_shape")]
    assert row["filled_from"] == "shape_codes" and row["value"] == "CIRCULAR"


@pytest.mark.parametrize(
    "direction, diameter, source",
    [("downstream", 0.3, "downstream pipe 1007"), ("upstream", 0.2, "upstream pipe 1000")],
)
def test_size_fill(sub, direction, diameter, source):
    """The J1 -> J2 pipes downstream are EGG and a 300 mm round pipe; upstream is A1 -> A2, whose
    200 mm diameter is itself inferred from its height and width."""
    report = sub.data_gaps(size_fill=direction).set_index(["element", "id", "gap"])
    row = report.loc[("pipe", NO_SIZE, "no_size")]
    assert row["filled_from"] == "size_fill"
    assert row["value"] == ("CIRCULAR", pytest.approx(diameter), 0.0)
    assert row["source"] == source


def test_size_overrides(sub):
    report = sub.data_gaps(size_overrides={str(NO_SIZE): 0.25, NO_SHAPE: {"shape": "CIRCULAR"}})
    report = report.set_index(["element", "id", "gap"])
    assert report.loc[("pipe", NO_SIZE, "no_size"), "filled_from"] == "size_overrides"
    assert report.loc[("pipe", NO_SIZE, "no_size"), "value"] == ("CIRCULAR", 0.25, 0.0)
    assert report.loc[("pipe", NO_SHAPE, "unknown_shape"), "filled_from"] == "size_overrides"


def test_interpolated_inverts_lie_on_one_slope(sub):
    """H1 (106.20) to H4 (105.40) over 100 m: a slope of 0.008."""
    report = sub.data_gaps(invert_fill="interpolate").set_index(["element", "id", "gap"])
    h2 = report.loc[("node", "H2", "no_node_invert")]
    h3 = report.loc[("node", "H3", "no_node_invert")]
    assert (h2["filled_from"], h3["filled_from"]) == ("interpolate", "interpolate")
    assert h2["value"] == pytest.approx(105.40 + 0.008 * 70)
    assert h3["value"] == pytest.approx(105.40 + 0.008 * 40)
    assert "H1" in h2["source"] and "H4" in h2["source"]
    # A headwater has nothing to interpolate from.
    assert pd.isna(report.loc[("node", "Q1", "no_node_invert"), "filled_from"])


def test_headwater_slope(sub):
    """Q1 is 50 m upstream of A2, whose invert is 104.40."""
    report = sub.data_gaps(headwater_slope=0.01).set_index(["element", "id", "gap"])
    row = report.loc[("node", "Q1", "no_node_invert")]
    assert row["filled_from"] == "headwater_slope"
    assert row["value"] == pytest.approx(104.40 + 0.5)
    # Nodes with a known invert upstream are not headwaters.
    assert pd.isna(report.loc[("node", "H2", "no_node_invert"), "filled_from"])


def test_force_main_discharge(sub):
    """The force main ends at J1 (invert 103.90): it discharges 0.2 m above it, after P2 is
    interpolated, and only the last force-main pipe is filled."""
    report = sub.data_gaps(
        invert_fill="interpolate", force_main_discharge=FILLS["force_main_discharge"]
    ).set_index(["element", "id", "gap"])
    row = report.loc[("pipe", 3007, "no_ds_invert")]
    assert row["filled_from"] == "force_main_discharge"
    assert row["value"] == pytest.approx(104.10)
    assert row["source"] == "node J1 invert 103.9"
    assert pd.isna(report.loc[("pipe", 3006, "no_ds_invert"), "filled_from"])
    p2 = report.loc[("node", "P2", "no_node_invert")]
    assert p2["value"] == pytest.approx((101.00 + 103.90) / 2)
    assert p2["order"] < row["order"]


def test_force_main_discharge_rule():
    assert force_main_discharge_rule({"max_height_above_receiving_invert": 0.2})(None, 100) == (
        pytest.approx(100.2)
    )
    rule = {"max_height_above_receiving_invert": 0.2, "height_above_receiving_invert": 0.1}
    assert force_main_discharge_rule(rule)(None, 100) == pytest.approx(100.1)
    with pytest.raises(ValueError, match="between 0"):
        force_main_discharge_rule({**rule, "height_above_receiving_invert": 0.3})
    with pytest.raises(ValueError, match="max_height"):
        force_main_discharge_rule({"height": 0.1})


def test_london_design_table():
    assert LONDON_ON.design["force_main_discharge"] == {"max_height_above_receiving_invert": 0.2}
    with pytest.raises(ValueError, match="unknown design rules"):
        Schema(name="x", columns=LONDON_ON.columns, design={"min_slope": 0.01})


@pytest.mark.parametrize(
    "option, match",
    [
        (dict(size_fill="largest"), "size_fill"),
        (dict(invert_fill="nearest"), "invert_fill"),
        (dict(headwater_slope=-0.01), "headwater_slope"),
        (dict(size_overrides={9999: 0.2}), "not in the subnetwork: 9999"),
        (dict(size_overrides={NO_SIZE: {"radius": 0.1}}), "radius"),
    ],
)
def test_invalid_fill_options(sub, option, match):
    with pytest.raises(ValueError, match=match):
        sub.data_gaps(**option)


def test_filled_model_is_flagged_and_runs_in_swmm(sub, tmp_path):
    path = tmp_path / "model.inp"
    with pytest.warns(UserWarning, match="adverse"):  # the force main runs uphill
        model = sub.to_inp(path, ["OUT"], template_path=_dynwave(tmp_path), **FILLS)
    filled = model.data_gaps.dropna(subset=["filled_from"])
    assert set(filled["filled_from"]) == {
        "height_width",
        "shape_codes",
        "size_fill",
        "interpolate",
        "headwater_slope",
        "force_main_discharge",
    }
    assert list(filled.sort_values("order")["order"]) == list(range(1, len(filled) + 1))

    s = read_sections(path)
    inverts = {row[0]: float(row[1]) for row in s["JUNCTIONS"]}
    assert inverts["H2"] == pytest.approx(105.96) and inverts["Q1"] == pytest.approx(104.90)
    conduits = {row[0]: row for row in s["CONDUITS"]}
    assert float(conduits["3007"][6]) == pytest.approx(0.2)  # outlet offset above J1
    xsections = {row[0]: row for row in s["XSECTIONS"]}
    assert xsections[str(NO_SHAPE)][1] == "CIRCULAR"
    assert float(xsections[str(NO_SIZE)][2]) == pytest.approx(0.3)
    assert run_swmm(path, tmp_path) == []


def _dynwave(tmp_path):
    from conftest import TEMPLATE

    path = tmp_path / "dynwave.inp"
    path.write_text(TEMPLATE.read_text().replace("KINWAVE", "DYNWAVE"))
    return path
