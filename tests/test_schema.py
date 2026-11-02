import math
import warnings

import geopandas as gpd
import pytest
from conftest import DENDRITIC, FIXTURE, quiet, read_sections, subnetwork_with_edges, write_inp
from shapely.geometry import LineString

from sewernetworkx import LONDON_ON, Network, Schema
from sewernetworkx.schema import CANONICAL_COLUMNS
from sewernetworkx.sewers import xsection


def london_rows(*rows):
    """GIS rows in London's schema, as a GeoDataFrame; each row overrides a default pipe."""
    default = {"UpstreamID": "U", "DownstreamID": "D", "Length": 10.0}
    records = [{"OBJECTID": i, **default, **row} for i, row in enumerate(rows)]
    geometry = [LineString([(0, 0), (1, 0)])] * len(records)
    return gpd.GeoDataFrame(records, geometry=geometry)


def london_pipe(**props):
    """One London GIS row, normalized."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return LONDON_ON.normalize(london_rows(props)).iloc[0]


@pytest.mark.parametrize(
    "code, expected",
    [
        ("R", ("CIRCULAR", 0.25, 0.0)),
        ("Round", ("CIRCULAR", 0.25, 0.0)),
        ("REC", ("RECT_CLOSED", 0.3, 0.45)),
        ("EGG", ("EGG", 0.3, 0.0)),
        ("ARC", ("ARCH", 0.3, 0.45)),
        ("ELH", ("HORIZ_ELLIPSE", 0.3, 0.45)),
        ("ELV", ("VERT_ELLIPSE", 0.3, 0.45)),
        ("S", None),
        ("", None),
        (None, None),
    ],
)
def test_london_shape_map(code, expected):
    pipe = london_pipe(PipeShape=code, Height=300, Width=450, Diameter=250)
    assert xsection(pipe) == (pytest.approx(expected) if expected else None)


@pytest.mark.parametrize("size", [None, float("nan"), 0])
def test_xsection_missing_size(size):
    assert xsection(london_pipe(PipeShape="R", Height=300, Width=300, Diameter=size)) is None


def test_london_normalize():
    pipe = london_pipe(
        UpstreamID="Abandoned",
        UpstreamInvert=0,
        DownstreamInvert=101.5,
        FlowType="STM",
        Abandoned=1,
        StreetName="Main St",
    )
    assert math.isnan(pipe["from_node"]) and pipe["to_node"] == "D"
    assert math.isnan(pipe["us_invert"]) and pipe["ds_invert"] == 101.5
    assert pipe["flow_type"] == "storm" and bool(pipe["abandoned"])
    assert pipe["StreetName"] == "Main St"
    assert "OBJECTID" not in pipe


def test_normalize_column_order():
    data = LONDON_ON.normalize(london_rows({"StreetName": "Main St"}))
    assert list(data.columns) == [*CANONICAL_COLUMNS, "StreetName", "geometry"]


def test_unmapped_flow_type_passes_through():
    assert london_pipe(FlowType="XYZ")["flow_type"] == "XYZ"


def test_unknown_shape_code_warns():
    with pytest.warns(UserWarning, match=r"\{'S': 2\}"):
        LONDON_ON.normalize(london_rows({"PipeShape": "S"}, {"PipeShape": "S"}, {"PipeShape": ""}))


def test_missing_required_column():
    with pytest.raises(ValueError, match=r"missing required columns \['Length'\]"):
        LONDON_ON.normalize(london_rows({}).drop(columns="Length"))


def test_duplicate_pipe_ids():
    data = london_rows({}, {})
    data["OBJECTID"] = 7
    with pytest.raises(ValueError, match=r"duplicate OBJECTID values \[7\]"):
        LONDON_ON.normalize(data)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({"columns": {"pipe_id": "id", "from_node": "a", "to_node": "b"}}, "required"),
        ({"columns": {**LONDON_ON.columns, "colour": "Colour"}}, "colour"),
        ({"columns": {**LONDON_ON.columns, "abandoned": "Abandoned"}}, "abandoned"),
        ({"columns": LONDON_ON.columns, "size_unit": "furlong"}, "furlong"),
    ],
)
def test_invalid_schema(kwargs, match):
    with pytest.raises(ValueError, match=match):
        Schema(name="bad", **kwargs)


# The fixture in another city's schema: other column names and codes, sizes in m, lengths in cm.
OTHER = Schema(
    name="Other",
    columns={
        "pipe_id": "PIPE",
        "from_node": "US_MH",
        "to_node": "DS_MH",
        "length": "LEN_CM",
        "us_invert": "US_INV",
        "ds_invert": "DS_INV",
        "drop_invert": "DROP_INV",
        "shape": "SHAPE",
        "diameter": "DIA_M",
        "height": "HT_M",
        "width": "WD_M",
        "flow_type": "SYSTEM",
        "material": "MAT",
    },
    shapes={"CIRC": "CIRCULAR", "BOX": "RECT_CLOSED", "OVOID": "EGG"},
    flow_types={"Sanitary": "sanitary", "Storm": "storm"},
    length_unit="cm",
)


def other_city(path):
    data = gpd.read_file(FIXTURE)
    data = data[data["Abandoned"] != 1]
    data = data[~data[["UpstreamID", "DownstreamID"]].isin(["Abandoned"]).any(axis=1)]
    data = data.rename(columns={LONDON_ON.columns[k]: v for k, v in OTHER.columns.items()})
    data["LEN_CM"] = data["LEN_CM"] * 100
    for column in ("DIA_M", "HT_M", "WD_M"):
        data[column] = data[column] / 1000
    for column in ("US_INV", "DS_INV", "DROP_INV"):
        data[column] = data[column].where(data[column] > 0)
    data["SHAPE"] = data["SHAPE"].map({"R": "CIRC", "Round": "CIRC", "REC": "BOX", "EGG": "OVOID"})
    data["SYSTEM"] = data["SYSTEM"].map({"SAN": "Sanitary", "STM": "Storm"})
    data.to_file(path, driver="GeoJSON")
    return path


def numbers(sections):
    """{SECTION: rows}, with numeric tokens as floats."""

    def value(token):
        try:
            return pytest.approx(float(token))
        except ValueError:
            return token

    return {k: [[value(t) for t in row] for row in rows] for k, rows in sections.items()}


def test_other_schema_writes_same_inp(network, tmp_path):
    other = quiet(Network, str(other_city(tmp_path / "other.geojson")), schema=OTHER)
    assert set(other.graph.edges) == set(network.graph.edges)
    london = read_sections(
        write_inp(subnetwork_with_edges(network, DENDRITIC), tmp_path, "a.inp", "DYNWAVE")
    )
    path = write_inp(subnetwork_with_edges(other, DENDRITIC), tmp_path, "b.inp", "DYNWAVE")
    assert numbers(read_sections(path)) == numbers(london)
