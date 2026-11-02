import pytest
from conftest import DENDRITIC, OVERFLOW, read_sections, run_swmm, subnetwork_with_edges, write_inp

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
    return write_inp(subnetwork_with_edges(network, DENDRITIC), tmp_path, routing="DYNWAVE")


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
    sub = subnetwork_with_edges(network, DENDRITIC + [("E1", "J2")])
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
    sub = subnetwork_with_edges(network, DENDRITIC + OVERFLOW + [("D2", "J2")], outfalls=("OUT",))
    report = sub.validate()
    assert report.splits == {"J1": [1006, 1007], "J2": [1005, 1008]}
    assert report.parallel_pipes == {("J1", "J2"): [1006, 1007]}
    assert report.dead_ends == ["OV2"]  # D2 has an outlet; OUT is an outfall
    assert report.cycles == []
    assert not report.is_dendritic


def test_validate_finds_cycles(network):
    sub = subnetwork_with_edges(network, [("A1", "A2"), ("A2", "J1"), ("J1", "A1"), ("A2", "OUT")])
    assert sub.validate().cycles == [["A1", "A2", "J1"]]


def test_dendritic_network_validates_clean(network):
    edges = [("A1", "A2"), ("A2", "J1"), ("B2", "J1"), ("C1", "J1")]
    report = subnetwork_with_edges(network, edges, outfalls=("J1",)).validate()
    assert report.is_dendritic and not report.dead_ends and not report.parallel_pipes


def test_split_rejected_under_kinematic_wave(network, tmp_path):
    sub = subnetwork_with_edges(network, DENDRITIC + OVERFLOW, outfalls=("OUT", "OV2"))
    with pytest.raises(ValueError, match=r"(?s)node J2 .*DYNWAVE"):
        write_inp(sub, tmp_path)
    assert not (tmp_path / "model.inp").exists()


def test_parallel_pipes_rejected_under_kinematic_wave(network, tmp_path):
    with pytest.raises(ValueError, match=r"node J1 has more than one outlet pipe \(1006, 1007\)"):
        write_inp(subnetwork_with_edges(network, DENDRITIC), tmp_path)


def test_cycle_rejected_under_kinematic_wave(network, tmp_path):
    sub = subnetwork_with_edges(network, [("A1", "A2"), ("A2", "A1")], outfalls=("A1",))
    with pytest.raises(ValueError, match="cycle A1 -> A2 -> A1"):
        write_inp(sub, tmp_path)


def test_multiple_outfalls_written_once(network, tmp_path):
    sub = subnetwork_with_edges(network, DENDRITIC + OVERFLOW, outfalls=("OUT", "OV2"))
    sub.outfall_node_type = ["FREE", "NORMAL"]
    rows = read_sections(write_inp(sub, tmp_path, routing="DYNWAVE"))["OUTFALLS"]
    # Invert = lowest downstream invert of the pipes entering the outfall.
    assert sorted((r[0], float(r[1]), r[2]) for r in rows) == [
        ("OUT", 103.30, "FREE"),
        ("OV2", 104.40, "NORMAL"),
    ]


def test_dead_end_rejected(network, tmp_path):
    sub = subnetwork_with_edges(
        network, [("A1", "A2"), ("A2", "J1")]
    )  # J1 has no outlet; OUT absent
    with pytest.raises(ValueError, match=r"dead ends .*: J1"):
        write_inp(sub, tmp_path)


def test_dead_end_written_as_outfall(network, tmp_path):
    sub = subnetwork_with_edges(network, [("A1", "A2"), ("A2", "J1"), ("B2", "J1")])
    s = read_sections(write_inp(sub, tmp_path, dead_ends="outfall"))
    assert [(r[0], float(r[1]), r[2]) for r in s["OUTFALLS"]] == [("J1", 103.95, "FREE")]
    assert {r[0] for r in s["JUNCTIONS"]} == {"A1", "A2", "B2"}


def test_manning_n_by_material(network, tmp_path):
    sub = subnetwork_with_edges(network, [("A1", "A2"), ("B2", "J1")], outfalls=())
    path = write_inp(
        sub, tmp_path, dead_ends="outfall", manning_n=0.02, manning_n_by_material={"PVC": 0.011}
    )
    n = {r[0]: float(r[4]) for r in read_sections(path)["CONDUITS"]}
    # A1-A2 is PVC; B2-J1 is CONC, which is not in the mapping, so it falls back to manning_n.
    # (E1-J2, the fixture pipe with no Material, also has no PipeShape and is rejected, #4.)
    assert n == {"1000": 0.011, "1003": 0.02}


def test_adverse_gis_pipe_is_reported(network, tmp_path):
    """D1-D2 is adverse in the GIS inverts: it is written as is, and reported."""
    sub = subnetwork_with_edges(network, DENDRITIC + [("D1", "D2"), ("D2", "J2")])
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
        subnetwork_with_edges(network, DENDRITIC),
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
            subnetwork_with_edges(network, DENDRITIC),
            tmp_path,
            routing="DYNWAVE",
            height_correction={"J2": 103.55},
        )


def test_offset_arguments_are_deprecated(network, tmp_path):
    with pytest.warns(DeprecationWarning, match="inlet_offset"):
        write_inp(
            subnetwork_with_edges(network, DENDRITIC), tmp_path, routing="DYNWAVE", inlet_offset=0.1
        )
