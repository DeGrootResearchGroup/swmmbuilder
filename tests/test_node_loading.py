"""NodeLoading: population, patterned, and steady flows as build_model options and a table
(#138), and require_unit_mean in pattern_inflows."""

import importlib.util

import numpy as np
import pandas as pd
import pytest
from conftest import DENDRITIC, REPO, TEMPLATE, read_sections, subnetwork, write_inp

from swmmbuilder import (
    DesignWarning,
    Inflow,
    NodeLoading,
    dwf_from_population,
    pattern_inflows,
)

FLAT = [1.0] * 24
DAY = [0.5] * 12 + [1.5] * 12
DIURNAL = [0.6] * 6 + [1.4] * 12 + [0.6] * 6
USES = ["COMMERCIAL", "INSTITUTIONAL", "INDUSTRIAL"]
PATTERNS = {"COMMERCIAL": DAY, "INSTITUTIONAL": FLAT, "INDUSTRIAL": DAY[::-1]}


@pytest.fixture
def parts():
    """(population, infiltration, non-residential flows) in the shapes of a build script: Series
    with zeros, and a DataFrame of nodes by land use with a mixed node and an empty one."""
    population = pd.Series({"A1": 120.0, "A2": 0.0, "B1": 80.0, "C1": 40.5, "J1": 0.0})
    infiltration = pd.Series({"A1": 0.02, "A2": 0.0, "B2": 0.015, "J1": 0.01, "J2": 0.0})
    nonres = pd.DataFrame(
        {
            "COMMERCIAL": {"A2": 0.3, "B2": 0.0, "C1": 0.2, "J2": 0.0},
            "INSTITUTIONAL": {"A2": 0.0, "B2": 0.4, "C1": 0.1, "J2": 0.0},
            "INDUSTRIAL": {"A2": 0.0, "B2": 0.0, "C1": 0.0, "J2": 0.0},
        }
    )
    return population.sort_index(), infiltration, nonres


def by_hand(population, infiltration, nonres, n_days=2):
    """The options and node table as a build script assembles them without NodeLoading."""
    dwf = dwf_from_population(population[population > 0].to_dict(), per_capita=230, unit="L/d")
    inflows, series = pattern_inflows(
        nonres, PATTERNS, n_days, prefix="NR_", baseline=infiltration[infiltration > 0].to_dict()
    )
    options = dict(
        dwf=dwf,
        patterns={"diurnal": DIURNAL},
        dwf_patterns="diurnal",
        inflows=inflows,
        timeseries=series,
    )
    table = (
        pd.DataFrame(
            {
                "population": population,
                "sanitary_L_s": pd.Series(dwf),
                "infiltration_L_s": infiltration,
                **{f"{u.lower()}_L_s": nonres[u] for u in USES},
            }
        )
        .fillna(0.0)
        .rename_axis("node")
    )
    return options, table


def loading_of(population, infiltration, nonres, n_days=2):
    return NodeLoading(
        population=population,
        per_capita=230,
        unit="L/d",
        pattern=DIURNAL,
        patterned=nonres,
        patterns=PATTERNS,
        days=n_days,
        prefix="NR_",
        steady=infiltration,
    )


RENAME = {
    "sanitary": "sanitary_L_s",
    "steady": "infiltration_L_s",
    **{u: f"{u.lower()}_L_s" for u in USES},
}


def test_matches_the_assembly_by_hand(parts):
    expected, expected_table = by_hand(*parts)
    loading = loading_of(*parts)
    options = loading.options()
    assert options.pop("flow_units") == "LPS"
    assert list(options) == list(expected)
    assert options == expected
    assert all(isinstance(i, Inflow) for i in options["inflows"].values())
    assert list(options["timeseries"]) == list(expected["timeseries"])
    assert "NR_C1" in options["timeseries"]  # the mixed node has its own series
    assert "A2" not in options["dwf"] and "A2" in options["inflows"]  # zero loads not written
    assert "J2" not in options["inflows"]
    pd.testing.assert_frame_equal(loading.table(RENAME), expected_table)


def test_writes_the_same_inp_file(parts, network, tmp_path):
    sub = subnetwork(network, DENDRITIC)
    expected, _ = by_hand(*parts)
    a = write_inp(sub, tmp_path, "by_hand.inp", routing="DYNWAVE", **expected)
    b = write_inp(sub, tmp_path, "loading.inp", routing="DYNWAVE", **loading_of(*parts).options())
    assert a.read_text() == b.read_text()
    sections = read_sections(b)
    assert len(sections["DWF"]) == 3 and len(sections["INFLOWS"]) == 5


def test_table_and_summary(parts):
    loading = loading_of(*parts)
    table = loading.table()
    assert table.index.name == "node"
    assert list(table.columns) == ["population", "sanitary", "steady", *USES]
    assert table.loc["J2"].sum() == 0  # given, with nothing
    s = loading.summary()
    sanitary = 240.5 * 230 / 86400
    assert s["population"] == pytest.approx(240.5) and s["population_nodes"] == 3
    assert s["sanitary_flow"] == pytest.approx(sanitary) and s["sanitary_nodes"] == 3
    assert s["steady_flow"] == pytest.approx(0.045) and s["steady_nodes"] == 3
    assert s["patterned_flow"] == pytest.approx(1.0) and s["patterned_nodes"] == 3
    assert s["COMMERCIAL_flow"] == pytest.approx(0.5)
    assert s["COMMERCIAL_nodes"] == 2
    assert s["INDUSTRIAL_flow"] == 0 and s["INDUSTRIAL_nodes"] == 0
    assert s["total"] == pytest.approx(sanitary + 1.045)
    assert s["nodes"] == 6  # A1, A2, B1, B2, C1, J1
    assert list(s)[:4] == ["total", "nodes", "population", "population_nodes"]
    assert all(type(v) in (int, float) for v in s.values())
    assert "sanitary=" in repr(loading)


def test_summary_former_keys_still_read(parts):
    s = loading_of(*parts).summary()
    for old, new in [
        ("sanitary", "sanitary_flow"),
        ("steady", "steady_flow"),
        ("patterned", "patterned_flow"),
        ("loaded_nodes", "nodes"),
        ("patterned_COMMERCIAL", "COMMERCIAL_flow"),
        ("patterned_COMMERCIAL_nodes", "COMMERCIAL_nodes"),
    ]:
        with pytest.warns(DeprecationWarning, match=f"'{old}' is now '{new}'"):
            assert s[old] == s[new]
        assert old not in s


@pytest.mark.parametrize("kind", ["sanitary", "steady", "population", "patterned"])
def test_a_patterned_kind_may_not_take_the_name_of_a_part(kind):
    with pytest.raises(ValueError, match=f"patterned kind '{kind}' is the name of a part"):
        NodeLoading(patterned={kind: {"A1": 1.0}}, patterns={kind: [1.0] * 24}, days=1)


def test_only_the_options_needed():
    assert set(NodeLoading(steady={"A1": 0.1, "B1": 0.0}).options()) == {"inflows", "flow_units"}
    assert NodeLoading(steady={"A1": 0.1}).options()["inflows"] == {"A1": Inflow(baseline=0.1)}
    options = NodeLoading(population={"A1": 100}, per_capita=0.001, unit="m3/s").options()
    assert options == {"dwf": {"A1": pytest.approx(0.1 * 1000)}, "flow_units": "LPS"}
    patterned = NodeLoading(patterned={"COM": {"A1": 0.2}}, patterns={"COM": FLAT}, days=1)
    assert set(patterned.options()) == {"inflows", "timeseries", "flow_units"}
    assert NodeLoading(steady={"A1": 0.0}).options() == {"flow_units": "LPS"}


def test_sanitary_pattern_by_name():
    loading = NodeLoading(
        population={"A1": 100},
        per_capita=250,
        pattern="dry_day",
        patterns={"dry_day": DIURNAL, "COM": FLAT},
        patterned={"COM": {"B1": 0.2}},
        days=1,
    )
    options = loading.options()
    assert options["patterns"] == {"dry_day": DIURNAL} and options["dwf_patterns"] == "dry_day"
    assert list(options["timeseries"]) == ["COM"]


def test_flow_units():
    loading = NodeLoading(population={"A1": 1000}, per_capita=86.4, flow_units="cms")
    assert loading.options()["dwf"]["A1"] == pytest.approx(0.001)
    assert loading.options()["flow_units"] == "CMS" == loading.flow_units


@pytest.mark.parametrize("kwargs", [{"steady": {"A1": 1.0}}, {"population": {"A1": 10}}])
def test_invalid_flow_units_refused(kwargs):
    # Refused without a population too, which used to be the only part that checked it (#158).
    with pytest.raises(ValueError, match="flow_units must be one of CMS, LPS, MLD, not 'CFS'"):
        NodeLoading(**kwargs, per_capita=200, flow_units="cfs")


def units_template(tmp_path, units):
    path = tmp_path / f"template_{units}.inp"
    path.write_text(
        TEMPLATE.read_text().replace("FLOW_UNITS           LPS", f"FLOW_UNITS           {units}")
    )
    return path


@pytest.mark.parametrize(
    "loading_units, template_units", [("LPS", "CMS"), ("CMS", "LPS"), ("MLD", "LPS")]
)
def test_load_in_other_units_than_the_template_refused(
    network, tmp_path, loading_units, template_units
):
    # An L/s load in a CMS template would be 1000 times too large (#158).
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1")])
    loading = NodeLoading(steady={"A1": 1.0}, flow_units=loading_units)
    with pytest.raises(ValueError, match=f"not the template's FLOW_UNITS {template_units}"):
        sub.to_inp(
            tmp_path / "m.inp",
            ["J1"],
            template_path=units_template(tmp_path, template_units),
            **loading.options(),
        )
    assert not (tmp_path / "m.inp").exists()


def test_load_in_the_template_units_written(network, tmp_path):
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1")])
    loading = NodeLoading(population={"A1": 1000}, per_capita=86.4, flow_units="CMS")
    path = tmp_path / "m.inp"
    sub.to_inp(path, ["J1"], template_path=units_template(tmp_path, "CMS"), **loading.options())
    (dwf,) = read_sections(path)["DWF"]
    assert dwf[0] == "A1" and float(dwf[2]) == pytest.approx(0.001)


@pytest.mark.parametrize(
    "kwargs, match",
    [
        ({}, "needs population, patterned, or steady"),
        ({"population": {"A1": 10}}, "population needs per_capita"),
        ({"steady": {"A1": 1.0}, "pattern": FLAT}, "needs population"),
        ({"patterned": {"COM": {"A1": 1.0}}, "patterns": {"COM": FLAT}}, "patterned needs days"),
        ({"population": {"A1": 10}, "per_capita": 200, "pattern": "x"}, "'x' is not in patterns"),
        ({"steady": {"A1": -1.0}}, "steady flows must be non-negative and finite: A1"),
        ({"population": {"A1": np.nan}, "per_capita": 200}, "populations must be non-negative"),
        ({"population": {"A1": 10}, "per_capita": 200, "pattern": [1.0] * 12}, "24 finite"),
        ({"population": {"A1": 10}, "per_capita": 200, "pattern": DAY[:23] + [2.0]}, "mean"),
        (
            {"patterned": {"COM": {"A1": 1.0}}, "patterns": {"COM": [1.1] * 24}, "days": 1},
            "pattern COM has a mean of 1.1",
        ),
        ({"patterned": {"COM": {"A1": 1.0}}, "days": 1}, "no pattern for COM"),
        ({"steady": {"A1": 1.0}, "require_unit_mean": -1}, "require_unit_mean"),
    ],
)
def test_refused(kwargs, match):
    with pytest.raises(ValueError, match=match):
        NodeLoading(**kwargs)


def test_require_unit_mean_can_be_relaxed():
    patterned = {"COM": {"A1": 1.0}}
    kwargs = dict(patterned=patterned, patterns={"COM": [1.1] * 24}, days=1)
    assert NodeLoading(**kwargs, require_unit_mean=False).options()["inflows"]["A1"].scale == 1
    NodeLoading(**kwargs, require_unit_mean=0.2)
    with pytest.raises(ValueError, match=r"tolerance 0\.05"):
        NodeLoading(**kwargs, require_unit_mean=0.05)


# ---- pattern_inflows(require_unit_mean=...) -----------------------------------------------------


def test_pattern_inflows_warns_of_a_mean_by_default(recwarn):
    # #160: warned at the tolerance of [DWF] patterns and NodeLoading, 0.001.
    with pytest.warns(DesignWarning, match=r"Pattern COM has a mean of 1\.002, not 1"):
        inflows, _ = pattern_inflows({"COM": {"A": 1.0}}, {"COM": [1.002] * 24}, days=1)
    assert inflows["A"].scale == 1.0
    recwarn.clear()
    pattern_inflows({"COM": {"A": 1.0}}, {"COM": [2.0] * 24}, days=1, require_unit_mean=False)
    pattern_inflows({"COM": {"A": 1.0}}, {"COM": [1.0005] * 24}, days=1)
    assert not [w for w in recwarn if issubclass(w.category, DesignWarning)]


@pytest.mark.parametrize("require", [True, 1e-3, 0.0])
def test_pattern_inflows_require_unit_mean_names_the_pattern(require):
    patterns = {"COM": FLAT, "IND": [1.01] * 24}
    with pytest.raises(ValueError, match=r"pattern IND has a mean of 1\.01, not 1"):
        pattern_inflows({"COM": {"A": 1.0}}, patterns, days=1, require_unit_mean=require)
    pattern_inflows({"COM": {"A": 1.0}}, {"COM": FLAT}, days=1, require_unit_mean=require)


def test_pattern_inflows_tolerance():
    patterns = {"COM": [1.01] * 24}
    pattern_inflows({"COM": {"A": 1.0}}, patterns, days=1, require_unit_mean=0.02)
    with pytest.raises(ValueError, match="require_unit_mean must be"):
        pattern_inflows({"COM": {"A": 1.0}}, patterns, days=1, require_unit_mean="yes")


# ---- the end-to-end example ----------------------------------------------------------------------


def test_loaded_catchment_example(tmp_path, capsys):
    path = REPO / "examples" / "loaded_catchment.py"
    spec = importlib.util.spec_from_file_location("loaded_catchment", path)
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    _, loading = example.main(tmp_path)
    sections = read_sections(tmp_path / "loaded_catchment.inp")
    summary = loading.summary()
    assert len(sections["DWF"]) == summary["sanitary_nodes"] > 0
    assert len(sections["INFLOWS"]) == len(loading.options()["inflows"]) > 0
    assert {row[0] for row in sections["PATTERNS"]} == {"diurnal"}
    assert {"LU_COMMERCIAL", "LU_INDUSTRIAL"} <= {row[0] for row in sections["TIMESERIES"]}
    assert summary["population"] == pytest.approx(sum(300 + 50 * k for k in range(6)))
    assert summary["COMMERCIAL_flow"] > 0 and summary["INDUSTRIAL_flow"] > 0
    table = pd.read_csv(tmp_path / "node_loading.csv", index_col="node")
    assert table["sanitary"].sum() == pytest.approx(summary["sanitary_flow"], rel=1e-5)
    assert "wrote" in capsys.readouterr().out
