from ._types import NodeId, PipeId
from .diagnostics import DataWarning, DesignWarning, OptionWarning, SwmmBuilderWarning
from .distribute import distribute_flow
from .land_use import LandUseAllocation, allocate_land_use, dominant_category
from .network import Network, strip_leading_zeros
from .population import PopulationAllocation, allocate_population
from .schema import LAZY_PRESETS, LONDON_ON, Schema
from .subnetwork import PumpStationReport, Subnetwork, TopologyReport, UpstreamTotals
from .swmm import (
    Derived,
    Inflow,
    Pump,
    PumpStation,
    SwmmModel,
    WetWell,
    dwf_from_population,
    hourly_series,
    pattern_inflows,
)

__all__ = [
    "BARRIE_ON",
    "KITCHENER_ON",
    "LONDON_ON",
    "PETERBOROUGH_ON",
    "SARNIA_ON",
    "TORONTO_ON",
    "WATERLOO_ON",
    "WINDSOR_ON",
    "DataWarning",
    "Derived",
    "DesignWarning",
    "Inflow",
    "LandUseAllocation",
    "Network",
    "NodeId",
    "OptionWarning",
    "PipeId",
    "PopulationAllocation",
    "Pump",
    "PumpStation",
    "PumpStationReport",
    "Schema",
    "Subnetwork",
    "SwmmBuilderWarning",
    "SwmmModel",
    "TopologyReport",
    "UpstreamTotals",
    "WetWell",
    "allocate_land_use",
    "allocate_population",
    "distribute_flow",
    "dominant_category",
    "dwf_from_population",
    "hourly_series",
    "pattern_inflows",
    "strip_leading_zeros",
]


def __getattr__(name: str) -> Schema:
    """The schema presets other than LONDON_ON, read on first use (see swmmbuilder.schema)."""
    if name in LAZY_PRESETS:
        from . import schema

        value = globals()[name] = getattr(schema, name)
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted([*globals(), *LAZY_PRESETS])
