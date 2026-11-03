"""SWMM's SI flow units and conversions, and checks on numbers."""

import math
from numbers import Real

# SWMM's SI flow units. With any of these, lengths, elevations and section sizes are in metres,
# matching the GIS data; with a US unit (CFS, GPM, MGD) SWMM would read them as feet.
SI_FLOW_UNITS = ("CMS", "LPS", "MLD")


# Litres per unit of volume, and seconds per unit of time, for dwf_from_population.
_VOLUME_LITRES = {"L": 1.0, "m3": 1000.0}


_TIME_SECONDS = {"s": 1.0, "d": 86400.0}


# Litres per second in each SI flow unit.
_FLOW_UNIT_LPS = {"LPS": 1.0, "CMS": 1000.0, "MLD": 1e6 / 86400.0}


def _is_finite_number(value):
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _is_real(value):
    """True for a finite real number that is not a bool or a string."""
    return isinstance(value, Real) and not isinstance(value, bool) and math.isfinite(value)
