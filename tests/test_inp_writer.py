import time

import pytest
from conftest import (
    DENDRITIC,
    OVERFLOW,
    read_sections,
    run_swmm,
    subnetwork,
    with_pipes,
    write_inp,
)

# GIS values for the fixture pipes used below (see tests/data/make_fixture.py).
GIS = {
    ("A1", "A2"): dict(shape="R", d=0.200, inv=(105.00, 104.50)),
    ("A2", "J1"): dict(shape="R", d=0.250, inv=(104.40, 104.00)),
    ("B1", "B2"): dict(shape="Round", d=0.200, inv=(104.80, 104.30)),
    ("B2", "J1"): dict(shape="REC", h=0.300, w=0.450, inv=(104.20, 103.95)),
    ("C1", "J1"): dict(shape="R", d=0.200, inv=(106.00, 105.50), drop=103.98),
    ("J2", "OUT"): dict(shape="R", d=0.600, inv=(103.50, 103.30)),
}


def conduits_by_pair(sections):
    """{(inlet, outlet): [(conduit row, xsection row), ...]}"""
    xs = {row[0]: row for row in sections.get("XSECTIONS", [])}
    out = {}
    for row in sections.get("CONDUITS", []):
        out.setdefault((row[1], row[2]), []).append((row, xs.get(row[0])))
    return out


def node_inverts(sections):
    rows = sections.get("JUNCTIONS", []) + sections.get("OUTFALLS", [])
    return {row[0]: float(row[1]) for row in rows}


@pytest.fixture
def dendritic_inp(network, tmp_path):
    return write_inp(subnetwork(network, DENDRITIC), tmp_path, routing="DYNWAVE")


def test_every_node_and_pipe_is_written(dendritic_inp):
    s = read_sections(dendritic_inp)
    nodes = {n for e in DENDRITIC for n in e}
    assert set(node_inverts(s)) == nodes
    assert set(conduits_by_pair(s)) == set(DENDRITIC)


def test_template_placeholders_filled(dendritic_inp):
    assert "{" not in dendritic_inp.read_text()


def test_flow_units_are_si(dendritic_inp):
    opts = {row[0]: row[1] for row in read_sections(dendritic_inp)["OPTIONS"] if len(row) > 1}
    assert opts["FLOW_UNITS"] in {"CMS", "LPS", "MLD"}


def test_round_pipes_are_circular(dendritic_inp):
    pairs = conduits_by_pair(read_sections(dendritic_inp))
    for pair in [("A1", "A2"), ("A2", "J1"), ("B1", "B2"), ("J2", "OUT")]:
        ((_, xs),) = pairs[pair]
        assert xs[1] == "CIRCULAR" and float(xs[2]) == pytest.approx(GIS[pair]["d"]), pair


def test_rectangular_and_egg_sections(dendritic_inp):
    pairs = conduits_by_pair(read_sections(dendritic_inp))
    ((_, rec),) = pairs[("B2", "J1")]
    assert (
        rec[1] == "RECT_CLOSED"
        and float(rec[2]) == pytest.approx(0.300)
        and float(rec[3]) == pytest.approx(0.450)
    )
    shapes = sorted(xs[1] for _, xs in pairs[("J1", "J2")])  # EGG pipe and parallel round pipe
    assert shapes == ["CIRCULAR", "EGG"]


def test_missing_shape_raises(network, tmp_path):
    sub = subnetwork(network, DENDRITIC + [("E1", "J2")])
    with pytest.raises(ValueError, match="1012"):
        write_inp(sub, tmp_path, routing="DYNWAVE")


def test_conduit_inverts_reproduce_gis(dendritic_inp):
    s = read_sections(dendritic_inp)
    inv = node_inverts(s)
    pairs = conduits_by_pair(s)
    for pair, g in GIS.items():
        ((row, _),) = pairs[pair]
        up_inv = inv[pair[0]] + float(row[5])
        down_inv = inv[pair[1]] + float(row[6])
        assert (up_inv, down_inv) == pytest.approx(g["inv"], abs=1e-3), pair
    assert inv["J1"] == pytest.approx(103.90)  # lowest pipe invert at J1 (outgoing EGG/round pipes)


def test_no_adverse_conduits_in_dendritic_network(dendritic_inp):
    s = read_sections(dendritic_inp)
    inv = node_inverts(s)
    for (a, b), rows in conduits_by_pair(s).items():
        for row, _ in rows:
            assert inv[a] + float(row[5]) >= inv[b] + float(row[6]), (a, b)


def test_written_inp_runs_in_swmm(dendritic_inp, tmp_path):
    assert run_swmm(dendritic_inp, tmp_path) == []


def test_validate_reports_topology(network):
    sub = subnetwork(network, DENDRITIC + OVERFLOW + [("D2", "J2")])
    report = sub.validate(["OUT"])
    assert report.splits == {"J1": [1006, 1007], "J2": [1005, 1008]}
    assert report.parallel_pipes == {("J1", "J2"): [1006, 1007]}
    assert report.dead_ends == ["OV2"]  # D2 has an outlet; OUT is an outfall
    assert report.cycles == []
    assert not report.is_dendritic


def test_validate_finds_cycles():
    sub = with_pipes(("J1", "A1")).subnetwork(["A1", "A2", "J1"])
    assert sub.validate().cycles == [["A1", "A2", "J1"]]


def test_validate_reports_each_loop_once():
    # Two loops, A1 <-> A2 and B1 -> B2 -> J1 -> B1, report one cycle each.
    sub = with_pipes(("A2", "A1"), ("J1", "B1")).subnetwork(["A1", "A2", "B1", "B2", "J1"])
    report = sub.validate()
    assert report.loops == [["A1", "A2"], ["B1", "B2", "J1"]]
    assert report.cycles == [["A1", "A2"], ["B1", "B2", "J1"]]
    assert report.summary()["cycles"] == 2
    assert not report.is_dendritic


def test_validate_reports_a_pipe_to_itself_as_a_loop():
    report = with_pipes(("A1", "A1")).subnetwork(["A1", "A2"]).validate()
    assert report.loops == [["A1"]] and report.cycles == [["A1"]]


def test_validate_handles_a_grid_of_two_way_pipes():
    # A 7 x 7 grid with every pipe drawn both ways holds far too many simple cycles to list
    # (a 6 x 6 one has about 2.4 million); it is one loop, shown by its shortest cycle.
    n = 7
    pairs = []
    for i in range(n):
        for j in range(n):
            here = f"G{i}{j}"
            for there in ([f"G{i + 1}{j}"] if i + 1 < n else []) + (
                [f"G{i}{j + 1}"] if j + 1 < n else []
            ):
                pairs += [(here, there), (there, here)]
    sub = with_pipes(*pairs).subnetwork({node for pair in pairs for node in pair})
    start = time.perf_counter()
    report = sub.validate()
    assert time.perf_counter() - start < 5
    assert report.loops == [sorted(sub.graph.nodes)]
    assert report.cycles == [["G00", "G01"]]
    assert report.summary()["cycles"] == 1


def test_dendritic_network_validates_clean(network):
    edges = [("A1", "A2"), ("A2", "J1"), ("B2", "J1"), ("C1", "J1")]
    report = subnetwork(network, edges).validate(["J1"])
    assert report.is_dendritic and not report.dead_ends and not report.parallel_pipes


def test_split_rejected_under_kinematic_wave(network, tmp_path):
    sub = subnetwork(network, DENDRITIC + OVERFLOW)
    with pytest.raises(ValueError, match=r"(?s)node J2 .*DYNWAVE"):
        write_inp(sub, tmp_path, outfalls=("OUT", "OV2"))
    assert not (tmp_path / "model.inp").exists()


def test_parallel_pipes_rejected_under_kinematic_wave(network, tmp_path):
    with pytest.raises(ValueError, match=r"node J1 has more than one outlet pipe \(1006, 1007\)"):
        write_inp(subnetwork(network, DENDRITIC), tmp_path)


def test_cycle_rejected_under_kinematic_wave(tmp_path):
    sub = with_pipes(("A2", "A1")).subnetwork(["A1", "A2"])
    with pytest.raises(ValueError, match="cycle A1 -> A2 -> A1"):
        write_inp(sub, tmp_path, outfalls=("A1",))


def test_multiple_outfalls_written_once(network, tmp_path):
    sub = subnetwork(network, DENDRITIC + OVERFLOW)
    outfalls = {"OUT": "FREE", "OV2": "NORMAL"}
    rows = read_sections(write_inp(sub, tmp_path, routing="DYNWAVE", outfalls=outfalls))["OUTFALLS"]
    # Invert = lowest downstream invert of the pipes entering the outfall.
    assert sorted((r[0], float(r[1]), r[2]) for r in rows) == [
        ("OUT", 103.30, "FREE"),
        ("OV2", 104.40, "NORMAL"),
    ]


def test_dead_end_rejected(network, tmp_path):
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1")])  # J1 has no outlet; OUT absent
    with pytest.raises(ValueError, match=r"dead ends .*: J1"):
        write_inp(sub, tmp_path, outfalls=())


def test_dead_end_written_as_outfall(network, tmp_path):
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1"), ("B2", "J1")])
    s = read_sections(write_inp(sub, tmp_path, outfalls=(), dead_ends="outfall"))
    assert [(r[0], float(r[1]), r[2]) for r in s["OUTFALLS"]] == [("J1", 103.95, "FREE")]
    assert {r[0] for r in s["JUNCTIONS"]} == {"A1", "A2", "B2"}


def test_manning_n_by_material(network, tmp_path):
    sub = subnetwork(network, [("A1", "A2"), ("A2", "J1"), ("B2", "J1")])
    path = write_inp(
        sub, tmp_path, outfalls=("J1",), manning_n=0.02, manning_n_by_material={"PVC": 0.011}
    )
    n = {r[0]: float(r[4]) for r in read_sections(path)["CONDUITS"]}
    # A1-A2 and A2-J1 are PVC; B2-J1 is CONC, which is not in the mapping, so it falls back to
    # manning_n. (E1-J2, the fixture pipe with no Material, also has no PipeShape and is
    # rejected, #4.)
    assert n == {"1000": 0.011, "1001": 0.011, "1003": 0.02}


def test_adverse_gis_pipe_is_reported(network, tmp_path):
    """D1-D2 is adverse in the GIS inverts: it is written as is, and reported."""
    sub = subnetwork(network, DENDRITIC + [("D1", "D2"), ("D2", "J2")])
    with pytest.warns(UserWarning, match=r"adverse slope.*1010 \(D1") as record:
        path = write_inp(sub, tmp_path, routing="DYNWAVE")
    assert len([w for w in record if "adverse" in str(w.message)]) == 1
    s = read_sections(path)
    inv = node_inverts(s)
    ((row, _),) = conduits_by_pair(s)[("D1", "D2")]
    assert inv["D1"] + float(row[5]) == pytest.approx(103.80)
    assert inv["D2"] + float(row[6]) == pytest.approx(104.10)


def test_height_correction_overrides_node_invert(network, tmp_path):
    path = write_inp(
        subnetwork(network, DENDRITIC),
        tmp_path,
        routing="DYNWAVE",
        height_correction={"J2": 103.40},
    )
    s = read_sections(path)
    assert node_inverts(s)["J2"] == pytest.approx(103.40)
    pairs = conduits_by_pair(s)
    ((row, _),) = pairs[("J2", "OUT")]
    assert float(row[5]) == pytest.approx(0.10)  # pipe keeps its GIS invert of 103.50


def test_height_correction_above_pipe_invert_warns(network, tmp_path):
    with pytest.warns(UserWarning, match="below the invert of node J2"):
        write_inp(
            subnetwork(network, DENDRITIC),
            tmp_path,
            routing="DYNWAVE",
            height_correction={"J2": 103.55},
        )


# ---- Field separation (#57) ----------------------------------------------------------------

# Columns of each section as written: a field that ran into the next one would change the count.
FIELD_COUNTS = {"JUNCTIONS": 6, "OUTFALLS": 4, "CONDUITS": 9, "XSECTIONS": 7, "COORDINATES": 3}


def test_long_names_and_numbers_stay_separate_fields(tmp_path):
    # Node IDs of 20+ characters (wider than the 16-character name columns) and sizes with a
    # long float repr (wider than the 16-character geom1 column).
    from conftest import FIXTURE

    from swmmbuilder import Network

    longer = lambda n: n if not isinstance(n, str) else f"MAINTENANCE_HOLE_{n}_X"  # noqa: E731
    network = Network(FIXTURE, normalize_ids=longer)
    sub = network.subnetwork([n for e in DENDRITIC for n in e])  # normalized on lookup
    pipe = sub.pipes()["pipe_id"].iloc[0]
    path = write_inp(
        sub,
        tmp_path,
        routing="DYNWAVE",
        outfalls=("OUT",),  # normalized like the other node IDs
        size_overrides={pipe: 0.1 + 0.2},  # 0.30000000000000004
    )
    sections = read_sections(path)
    for section, count in FIELD_COUNTS.items():
        assert sections[section], section
        assert {len(row) for row in sections[section]} == {count}, section
    xs = {row[0]: row for row in sections["XSECTIONS"]}
    assert xs[str(pipe)][2:4] == ["0.3", "0.0"]
    assert all(len(row[1]) > 16 for row in sections["CONDUITS"])
    assert run_swmm(path, tmp_path) == []


def test_numbers_are_rounded_and_shortest():
    from swmmbuilder.swmm import _line

    assert _line(("P1", 4), (0.1 + 0.2, 0), (12 / (1 / 0.0254), 0), (-0.0, 0), (3, 0)) == (
        "P1   0.3 0.3048 0.0 3\n"
    )
    assert _line(("A", 0), (245.123456789, 0)) == "A 245.123457\n"


@pytest.mark.parametrize("name", ["MH 12", "MH\t12", "MH;12", ""])
def test_names_swmm_cannot_read_are_refused(name):
    from swmmbuilder import swmm

    model = swmm.SwmmModel(junctions=[swmm.Junction(name, 100.0, 0.0, 0.0, 0.0, 0.0)])
    with pytest.raises(ValueError, match="whitespace or ';'"):
        model.sections()


def test_flows_keep_their_significant_figures():
    """A small flow in CMS is not rounded to DECIMALS decimals (#123)."""
    from swmmbuilder.swmm.records import _flow, _line

    assert _line(("A", 0), (3.4e-6, 0)) == "A 3e-06\n"  # a length: 6 decimals
    assert _flow(3.4e-6) == "3.4e-06" and _flow(4e-7) == "4e-07"
    assert _flow(0.1 + 0.2) == "0.3" and _flow(12) == "12" and _flow(-0.0) == "0.0"
    with pytest.raises(TypeError):
        _flow(True)


def test_small_flows_are_written_in_full(network, tmp_path):
    from swmmbuilder import swmm

    model = swmm.build_model(subnetwork(network, DENDRITIC), ["OUT"], dwf={"A1": 3.4e-6})
    assert model.sections()["dwf"].split()[2] == "3.4e-06"
