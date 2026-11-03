"""The package layout (#76): swmmbuilder.swmm keeps every name of the former single module, and
the column-mapping layer does not depend on the .inp writer."""

import ast

import pytest
from conftest import REPO

import swmmbuilder.swmm

# Every top-level name of swmmbuilder/swmm.py before it was split, and the names it imported.
FORMER_NAMES = """
Conduit Curve DECIMALS Derived DryWeatherFlow FORCE_MAIN_CHECKS FORCE_MAIN_EQUATIONS
HEAD_FLOW_CURVES Junction MAX_DWF_PATTERNS NOT_GRAVITY Outfall PATTERN_LENGTHS
PATTERN_TYPES_BY_LENGTH PUMP_CURVE_TYPES Pattern Pump PumpLink PumpStation SI_FLOW_UNITS
StationPump Storage SwmmModel TEMPLATE_PLACEHOLDERS WET_WELL_DEFAULTS WET_WELL_OPTIONS WetWell
XSection _FLOW_UNIT_LPS _TIME_SECONDS _VOLUME_LITRES _is_finite_number _is_real _line _name
_number _rule_params apply_wet_well_defaults build_model check_design_rule
check_force_main_equation check_force_mains check_template check_topology default_template_path
dwf_from_population force_main_check_params force_main_roughness_rule force_main_roughnesses
gravity_inlets logger node_inverts outfall_types pump_level_problem pumped_flow read_template
resolve_dwf resolve_dwf_patterns resolve_patterns resolve_pump_stations resolve_pumps
resolve_wet_wells template_flow_units template_force_main_equation template_routing
wet_well_params write_inp
_invert fill_gaps height_corrections lowest_inverts xsection line_coords
""".split()


@pytest.mark.parametrize("name", FORMER_NAMES)
def test_former_swmm_names_are_importable(name):
    assert hasattr(swmmbuilder.swmm, name)


def _package_imports(module):
    """The swmmbuilder modules that ``module`` (a path under swmmbuilder/) imports."""
    tree = ast.parse((REPO / "swmmbuilder" / module).read_text())
    return {
        node.module or ""
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
        if node.level > 0
    }


@pytest.mark.parametrize("module", ["schema.py", "design.py", "pumps.py", "units.py"])
def test_data_and_design_layers_do_not_import_the_writer(module):
    assert not any(m == "swmm" or m.startswith("swmm.") for m in _package_imports(module))


def test_every_subpackage_is_packaged():
    # The tests run on an editable install, which finds every module; a regular install takes
    # only the packages pyproject.toml selects (#84: swmmbuilder.swmm was left out). Checked
    # with the standard library only: setuptools is not installed in every test environment.
    import tomllib
    from fnmatch import fnmatchcase

    config = tomllib.loads((REPO / "pyproject.toml").read_text())["tool"]["setuptools"]
    on_disk = {
        ".".join(p.relative_to(REPO).parts[:-1])
        for p in (REPO / "swmmbuilder").rglob("__init__.py")
    }
    if isinstance(config.get("packages"), list):
        missing = on_disk - set(config["packages"])
    else:
        find = config["packages"]["find"]
        include, exclude = find.get("include", ["*"]), find.get("exclude", [])
        missing = {
            package
            for package in on_disk
            if not any(fnmatchcase(package, pattern) for pattern in include)
            or any(fnmatchcase(package, pattern) for pattern in exclude)
        }
    assert not missing, f"Subpackages a regular install leaves out: {sorted(missing)}"
