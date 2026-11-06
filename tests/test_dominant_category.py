"""The category of a layer covering most of each polygon (swmmbuilder.dominant_category, #139)."""

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

from swmmbuilder import dominant_category

CRS = "EPSG:2958"
X0, Y0 = 480000.0, 4760000.0


def frame(boxes, crs=CRS, index=None, **columns):
    geometry = [box(X0 + x0, Y0 + y0, X0 + x1, Y0 + y1) for x0, y0, x1, y1 in boxes]
    return gpd.GeoDataFrame(columns, geometry=geometry, index=index, crs=crs)


# Blocks of 100 x 100 m: a is 70 % residential, b is split evenly between two uses (a tie), c is
# all commercial over two features, and d meets no zoning.
BLOCKS = frame(
    [(0, 0, 100, 100), (100, 0, 200, 100), (200, 0, 300, 100), (1000, 0, 1100, 100)],
    index=pd.Index(["a", "b", "c", "d"], name="id"),
    population=[10, 20, 30, 40],
)
ZONING = frame(
    [
        (0, 0, 70, 100),
        (70, 0, 150, 100),
        (150, 0, 200, 100),
        (200, 0, 250, 100),
        (250, 0, 300, 100),
    ],
    GeneralizedLandUse=["RESIDENTIAL", "COMMERCIAL", "RESIDENTIAL", "COMMERCIAL", "COMMERCIAL"],
)


def by_pivot(blocks, zoning):
    """The paper's build script: the residential flag from a pivot table and idxmax."""
    zoning = zoning.to_crs(blocks.crs)
    parts = gpd.overlay(
        blocks[["geometry"]].reset_index(),
        zoning[["GeneralizedLandUse", "geometry"]],
        how="intersection",
        keep_geom_type=True,
    )
    area = parts.assign(a=parts.area).pivot_table(
        index="id", columns="GeneralizedLandUse", values="a", aggfunc="sum"
    )
    return (area.idxmax(axis=1) == "RESIDENTIAL").reindex(blocks.index, fill_value=False)


def test_category_and_share():
    use = dominant_category(BLOCKS, ZONING, "GeneralizedLandUse")
    assert list(use.index) == ["a", "b", "c", "d"]
    assert list(use.columns) == ["category", "share"]
    assert use.loc["a", "category"] == "RESIDENTIAL"
    assert use.loc["b", "category"] == "COMMERCIAL"  # a tie goes to the first in sorted order
    assert use.loc["c", "category"] == "COMMERCIAL"
    assert pd.isna(use.loc["d", "category"])
    np.testing.assert_allclose(use["share"], [0.7, 0.5, 1.0, 0.0])


def test_flag_matches_pivot_table():
    flag = dominant_category(BLOCKS, ZONING, "GeneralizedLandUse")["category"] == "RESIDENTIAL"
    pd.testing.assert_series_equal(flag, by_pivot(BLOCKS, ZONING), check_names=False)
    assert list(flag) == [True, False, False, False]


def test_layer_is_reprojected():
    # Block b's exact tie does not survive the round trip, so it is left out.
    blocks = BLOCKS.loc[["a", "c", "d"]]
    use = dominant_category(blocks, ZONING.to_crs("EPSG:4326"), "GeneralizedLandUse")
    assert list(use["category"].fillna("")) == ["RESIDENTIAL", "COMMERCIAL", ""]
    np.testing.assert_allclose(use["share"], [0.7, 1.0, 0.0], atol=1e-3)


def test_overlapping_features_cap_the_share():
    zoning = frame([(0, 0, 100, 100), (0, 0, 100, 100)], GeneralizedLandUse=["PARK", "PARK"])
    use = dominant_category(BLOCKS.loc[["a"]], zoning, "GeneralizedLandUse")
    assert use.loc["a", "category"] == "PARK"
    assert use.loc["a", "share"] == 1.0


def test_missing_categories_are_ignored():
    zoning = ZONING.copy()
    zoning.loc[0, "GeneralizedLandUse"] = None
    use = dominant_category(BLOCKS, zoning, "GeneralizedLandUse")
    assert use.loc["a", "category"] == "COMMERCIAL"
    assert use.loc["a", "share"] == pytest.approx(0.3)


def test_no_overlap_gives_missing_categories():
    use = dominant_category(BLOCKS.loc[["d"]], ZONING, "GeneralizedLandUse")
    assert pd.isna(use.loc["d", "category"])
    assert use.loc["d", "share"] == 0.0


def test_invalid_input_lists_every_problem():
    repeated = BLOCKS.set_axis(["a", "a", "b", "c"]).to_crs("EPSG:4326")
    with pytest.raises(ValueError) as error:
        dominant_category(repeated, ZONING.set_crs(None, allow_override=True), "use")
    message = str(error.value)
    assert "no column 'use'" in message
    assert "repeated: ['a']" in message
    assert "the layer has no CRS" in message
    assert "projected in metres" in message
