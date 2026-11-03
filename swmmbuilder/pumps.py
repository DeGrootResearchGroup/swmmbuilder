"""Pumps, wet wells and pumping stations, as given to the .inp writer."""

from dataclasses import dataclass, replace

from .units import _is_finite_number

# SWMM pump curve types, and what the x values of each are (y is flow, in the FLOW_UNITS).
PUMP_CURVE_TYPES = {
    "PUMP1": "wet well volume",
    "PUMP2": "inlet node depth",
    "PUMP3": "head difference",
    "PUMP4": "inlet node depth",
    "PUMP5": "head difference",
}


@dataclass
class Pump:
    """A pump that replaces a pipe in the .inp file (see the ``pumps`` option of
    :func:`build_model`).

    Attributes
    ----------
    curve : list of (float, float) or str
        Pump curve points (x, y), y being the flow in the template's FLOW_UNITS (L/s with the
        shipped template) and x as given by ``curve_type``; or ``'*'`` for an ideal pump, which
        pumps whatever flows into its inlet node.

    curve_type : str
        PUMP1 (x: wet well volume, m3; the inlet must be a wet well), PUMP2 (x: inlet depth, m,
        stepwise), PUMP3 (x: head difference, m), PUMP4 (x: inlet depth, m, continuous) or PUMP5
        (variable-speed PUMP3). Ignored for an ideal pump. Required for a pump station's curve;
        a pumped pipe's (``pumps``) defaults to PUMP3.

    status : str
        Initial status, 'ON' or 'OFF'.

    startup, shutoff : float, optional
        Inlet node depths (m) at which the pump turns on and off; 0 for none. Not given (None),
        they are written as 0, unless a pump station's ``wet_well_defaults`` rule derives them.
    """

    curve: object
    curve_type: str | None = None
    status: str = "ON"
    startup: float | None = None
    shutoff: float | None = None

    @classmethod
    def coerce(cls, pump_id, value, curve_type=None):
        """A Pump from a Pump or a dict of its fields, checked, with ``curve_type`` if it has none;
        raise ValueError if it is invalid."""
        pump = replace(value) if isinstance(value, cls) else cls(**value)
        if pump.curve_type is None:
            pump.curve_type = curve_type
        problems = []
        if isinstance(pump.curve, str):
            if pump.curve != "*":
                problems.append(f"curve must be a list of (x, y) points or '*', not {pump.curve!r}")
        else:
            points = [tuple(pt) for pt in pump.curve]
            if not points or any(len(pt) != 2 for pt in points):
                problems.append("curve must be a non-empty list of (x, y) points")
            else:
                pump.curve = [(float(x), float(y)) for x, y in points]
            if pump.curve_type is None:
                problems.append(f"a curve needs a curve_type, one of {list(PUMP_CURVE_TYPES)}")
            elif str(pump.curve_type).upper() not in PUMP_CURVE_TYPES:
                problems.append(
                    f"curve_type {pump.curve_type!r} is not one of {list(PUMP_CURVE_TYPES)}"
                )
            else:
                pump.curve_type = pump.curve_type.upper()
        pump.status = pump.status.upper()
        if pump.status not in ("ON", "OFF"):
            problems.append(f"status must be 'ON' or 'OFF', not {pump.status!r}")
        if problems:
            raise ValueError(f"Pump {pump_id}: {'; '.join(problems)}.")
        return pump


@dataclass
class WetWell:
    """A node written as a [STORAGE] node of constant plan area, e.g. a pumping station wet well
    (see the ``wet_wells`` option of :func:`build_model`).

    Attributes
    ----------
    max_depth : float
        Depth (m) from the node invert to the top of the wet well.

    area : float
        Plan area (m2).

    init_depth : float
        Initial water depth (m).
    """

    max_depth: float
    area: float
    init_depth: float = 0.0

    @classmethod
    def coerce(cls, node, value):
        """A WetWell from a WetWell or a dict of its fields; raise ValueError if it is invalid."""
        well = replace(value) if isinstance(value, cls) else cls(**value)
        if not (well.max_depth > 0 and well.area > 0 and well.init_depth >= 0):
            raise ValueError(
                f"Wet well {node}: max_depth and area must be positive and init_depth not negative."
            )
        return well


@dataclass
class PumpStation:
    """A pumping station that discharges into one or more force mains (see the
    ``pump_stations`` option of :func:`build_model`).

    Attributes
    ----------
    pump : Pump or dict, optional
        The pump: an ideal pump (``curve='*'``) suits a variable-speed station that delivers a
        near-uniform flow, a curve with startup and shutoff depths (e.g. PUMP2 or PUMP4) a
        fixed-speed station that cycles on and off. There is no built-in default: a station
        without one takes the ``default_pump`` of :func:`build_model`, which must then be given.

    wet_well : WetWell or dict, optional
        The station node's wet well; without one the station node is a junction.

    force_mains : pipe ID or list of pipe IDs, optional
        The first pipes of the force mains to pump into (default: every force main leaving the
        station node).

    discharge : {'shared', 'separate'}, optional
        With more than one force main: one pump whose discharge node feeds them all
        (``'shared'``), or one pump per force main (``'separate'``). Required then.

    design_flow : float, optional
        Design average flow into the station, in the template's FLOW_UNITS (L/s with the shipped
        template). Used only by the ``wet_well_defaults`` rule (the wet well's plan area) and the
        ``force_main_checks`` of an ideal pump or a head-flow curve (its pumped flow); the GIS
        data do not carry it.

    suction_diameter : float, optional
        Diameter of the pump suction (mm), for the ``low_level_suction_diameters`` parameter of
        the ``wet_well_defaults`` rule.
    """

    pump: object = None
    wet_well: object = None
    force_mains: object = None
    discharge: str | None = None
    design_flow: float | None = None
    suction_diameter: float | None = None

    @classmethod
    def coerce(cls, node, value, default_pump=None):
        """A PumpStation from a PumpStation or a dict of its fields, with its pump (or
        ``default_pump``) and wet well coerced; raise ValueError if it is invalid."""
        station = replace(value) if isinstance(value, cls) else cls(**value)
        if station.pump is None:
            station.pump = default_pump
        if station.pump is None:
            raise ValueError(
                f"Pump station {node} has no pump: pass its pump, or a default_pump for the "
                "stations."
            )
        station.pump = Pump.coerce(node, station.pump)
        if station.wet_well is not None:
            station.wet_well = WetWell.coerce(node, station.wet_well)
        if station.force_mains is not None and (
            isinstance(station.force_mains, str) or not hasattr(station.force_mains, "__iter__")
        ):
            station.force_mains = [station.force_mains]
        for name in ("design_flow", "suction_diameter"):
            value = getattr(station, name)
            if value is not None:
                if not _is_finite_number(value) or not float(value) > 0:
                    raise ValueError(f"Pump station {node}: {name} must be positive.")
                setattr(station, name, float(value))
        if station.discharge not in (None, "shared", "separate"):
            raise ValueError(
                f"Pump station {node}: discharge must be 'shared' or 'separate', "
                f"not {station.discharge!r}."
            )
        return station


@dataclass(frozen=True)
class StationPump:
    """A pump of a :class:`PumpStation`, from the station node to its discharge node, and the
    force mains it feeds: {first pipe ID: pipe IDs of the chain}."""

    name: str
    station: str
    discharge: str
    force_mains: dict
