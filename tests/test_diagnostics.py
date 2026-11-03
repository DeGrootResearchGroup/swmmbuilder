"""How problems are reported (#79): the Problems collector, the warning categories, and the
channel each kind of problem goes to (see swmmbuilder.diagnostics)."""

import ast
import warnings

import pandas as pd
import pytest
from conftest import DENDRITIC, REPO, subnetwork
from test_nodes import load

import swmmbuilder
from swmmbuilder import maps
from swmmbuilder.diagnostics import Problems
from swmmbuilder.swmm import build_model


def test_problems_collects_and_raises_one_error():
    problems = Problems("Invalid things")
    problems.raise_if_any()  # nothing added: no error
    assert not problems and len(problems) == 0
    problems.add("a is wrong")
    problems.extend(["b is wrong", "c is wrong"])
    assert len(problems) == 3
    with pytest.raises(ValueError) as error:
        problems.raise_if_any()
    assert str(error.value) == "Invalid things:\n  a is wrong\n  b is wrong\n  c is wrong"
    assert Problems("T", ["x"]).message() == "T:\n  x"


@pytest.mark.parametrize("name", ["DataWarning", "DesignWarning", "OptionWarning"])
def test_warning_categories(name):
    category = getattr(swmmbuilder, name)
    assert issubclass(category, swmmbuilder.SwmmBuilderWarning)
    assert issubclass(category, UserWarning)  # pytest.warns(UserWarning) keeps working


def test_every_warning_has_a_category():
    # warnings.warn(message) alone is a bare UserWarning, which users cannot filter by kind.
    bare = []
    for path in sorted((REPO / "swmmbuilder").rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr == "warn"
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "warnings"
                and len(node.args) < 2
                and not any(k.arg == "category" for k in node.keywords)
            ):
                bare.append(f"{path.relative_to(REPO)}:{node.lineno}")
    assert not bare


def test_rim_below_invert_is_a_data_warning():
    # Raising M2's invert above its rim (247.9): the junction falls back to max_depth.
    sub = load(include_overflows=False, snap_tolerance=0.5, max_invert_difference=1.0).catchment(
        "OF1"
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = build_model(sub, max_depth=1.5, height_correction={"M2": 400.0})
    assert {j.name: j.max_depth for j in model.junctions}["M2"] == 1.5
    messages = [str(w.message) for w in caught]
    assert any(m.startswith("Node M2: rim elevation") for m in messages)
    # The rim, and the pipes now below M2 and sloping up to it: all suspect data.
    assert {type(w.message) for w in caught} == {swmmbuilder.DataWarning}


def test_plot_values_with_unknown_ids_are_an_option_warning(network):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    try:
        with pytest.warns(swmmbuilder.OptionWarning, match="None of the IDs of the values"):
            maps.plot_pipes(ax, subnetwork(network, DENDRITIC), pd.Series({"nope": 1.0}))
    finally:
        plt.close(fig)


def test_adverse_slope_is_a_data_warning(network):
    # Pipe 1000 (A1 -> A2) made to slope uphill by raising A2 above A1's pipe invert.
    sub = subnetwork(network, DENDRITIC)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        build_model(sub, ["OUT"], height_correction={"A2": 200.0})
    kinds = {type(w.message).__name__ for w in caught if "adverse slope" in str(w.message)}
    assert kinds == {"DataWarning"}
