"""SWMM models of a subnetwork, and their .inp files.

:func:`build_model` turns a subnetwork into a :class:`SwmmModel`: plain records of the junctions,
outfalls, storage nodes, conduits, pumps, cross-sections, curves, dry-weather flows, time patterns,
external inflows, time series, coordinates and vertices, with all the hydraulic decisions
(node inverts, conduit offsets, cross-sections, Manning's n) already made.
:meth:`SwmmModel.render` fills an .inp template with them, and :func:`write_inp` checks the
template and the topology, builds the model and writes the file.

The code is split into modules (template, pumping, dwf, inflows, records, build) and the
shared :mod:`swmmbuilder.units`, :mod:`swmmbuilder.pumps` and :mod:`swmmbuilder.design`; every
name is re-exported here, so ``from swmmbuilder.swmm import X`` keeps working.
"""

# Imported by the old single module, and so reachable as swmmbuilder.swmm.<name>.
from ..design import (  # noqa: F401
    FORCE_MAIN_CHECKS,
    FORCE_MAIN_EQUATIONS,
    HEAD_FLOW_CURVES,
    NOT_GRAVITY,
    WET_WELL_DEFAULTS,
    WET_WELL_OPTIONS,
    Derived,
    _rule_params,
    apply_wet_well_defaults,
    check_design_rule,
    check_force_main_equation,
    check_force_mains,
    force_main_check_params,
    force_main_roughness_rule,
    force_main_roughnesses,
    gravity_inlets,
    pump_level_problem,
    pumped_flow,
    wet_well_params,
)
from ..gaps import _invert, fill_gaps, height_corrections, lowest_inverts, xsection  # noqa: F401
from ..geometry import line_coords  # noqa: F401
from ..pumps import (  # noqa: F401
    PUMP_CURVE_TYPES,
    Pump,
    PumpStation,
    StationPump,
    WetWell,
)
from ..units import (  # noqa: F401
    _FLOW_UNIT_LPS,
    _TIME_SECONDS,
    _VOLUME_LITRES,
    SI_FLOW_UNITS,
    _is_finite_number,
    _is_real,
)
from .build import (  # noqa: F401
    BUILD_OPTIONS,
    FILL_OPTIONS,
    build_model,
    logger,
    node_inverts,
    outfall_types,
    write_inp,
)
from .dwf import (  # noqa: F401
    MAX_DWF_PATTERNS,
    PATTERN_LENGTHS,
    PATTERN_TYPES_BY_LENGTH,
    dwf_from_population,
    resolve_dwf,
    resolve_dwf_patterns,
    resolve_patterns,
)
from .inflows import (  # noqa: F401
    Inflow,
    hourly_series,
    pattern_inflows,
    resolve_inflows,
    resolve_timeseries,
)
from .pumping import (  # noqa: F401
    resolve_pump_stations,
    resolve_pumps,
    resolve_wet_wells,
)
from .records import (  # noqa: F401
    DECIMALS,
    Conduit,
    Curve,
    DryWeatherFlow,
    ExternalInflow,
    Junction,
    Outfall,
    Pattern,
    PumpLink,
    Storage,
    SwmmModel,
    TimeSeries,
    XSection,
    _line,
    _name,
    _number,
)
from .template import (  # noqa: F401
    TEMPLATE_PLACEHOLDERS,
    check_template,
    check_topology,
    default_template_path,
    read_template,
    template_flow_units,
    template_force_main_equation,
    template_routing,
)
