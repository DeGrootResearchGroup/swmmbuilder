import subprocess
import sys


def test_import_needs_no_optional_dependency():
    # pyswmm (only for running .inp files), gmplot (deprecated generate_gmap) and pycancensus
    # (census polygons) are optional extras, so importing the package must not import them.
    optional = ("pyswmm", "gmplot", "folium", "contextily", "pycancensus")
    code = (
        "import sys, swmmbuilder, swmmbuilder.maps, swmmbuilder.census; "
        f"print(sorted(m for m in {optional!r} if m in sys.modules))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "[]"
