import subprocess
import sys


def test_import_needs_no_optional_dependency():
    # pyswmm (only for running .inp files) and gmplot (deprecated generate_gmap) are optional
    # extras, so importing the package must not import them.
    code = (
        "import sys, swmmbuilder, swmmbuilder.maps; "
        "print(sorted(m for m in ('pyswmm', 'gmplot', 'folium', 'contextily') if m in sys.modules))"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert result.stdout.strip() == "[]"
