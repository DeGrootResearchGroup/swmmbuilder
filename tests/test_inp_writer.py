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
    return write_inp(subnetwork_with_edges(network, DENDRITIC), tmp_path)


def test_every_node_and_pipe_is_written(dendritic_inp):
    s = read_sections(dendritic_inp)
    nodes = {n for e in DENDRITIC for n in e}
    assert set(node_inverts(s)) == nodes
    assert set(conduits_by_pair(s)) == set(DENDRITIC)


@pytest.mark.xfail(strict=True, reason="#5: the {pumps} template placeholder is never filled")
def test_template_placeholders_filled(dendritic_inp):
    assert "{" not in dendritic_inp.read_text()


@pytest.mark.xfail(
    strict=True, reason="#5: template uses FLOW_UNITS MGD (US units) while geometry is metric"
)
def test_flow_units_are_si(dendritic_inp):
    opts = {row[0]: row[1] for row in read_sections(dendritic_inp)["OPTIONS"] if len(row) > 1}
    assert opts["FLOW_UNITS"] in {"CMS", "LPS", "MLD"}


@pytest.mark.xfail(strict=True, reason="#4: PipeShape 'R' (round) is written as RECT_CLOSED")
def test_round_pipes_are_circular(dendritic_inp):
    pairs = conduits_by_pair(read_sections(dendritic_inp))
    for pair in [("A1", "A2"), ("A2", "J1"), ("B1", "B2"), ("J2", "OUT")]:
        ((_, xs),) = pairs[pair]
        assert xs[1] == "CIRCULAR" and float(xs[2]) == pytest.approx(GIS[pair]["d"]), pair


@pytest.mark.xfail(strict=True, reason="#4: only 'R' and 'Round' shapes are handled")
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


@pytest.mark.xfail(
    strict=True, reason="#4: missing PipeShape silently reuses the previous pipe's section"
)
def test_missing_shape_raises(network, tmp_path):
    sub = subnetwork_with_edges(network, DENDRITIC + [("E1", "J2")])
    with pytest.raises(ValueError, match="1012"):
        write_inp(sub, tmp_path)


@pytest.mark.xfail(
    strict=True, reason="#4: node invert = first outgoing pipe; offsets 0; drop structures ignored"
)
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


@pytest.mark.xfail(
    strict=True, reason="#4: J2 takes the overflow pipe's invert, making J1->J2 adverse"
)
def test_no_adverse_conduits_in_dendritic_network(dendritic_inp):
    s = read_sections(dendritic_inp)
    inv = node_inverts(s)
    for (a, b), rows in conduits_by_pair(s).items():
        for row, _ in rows:
            assert inv[a] + float(row[5]) >= inv[b] + float(row[6]), (a, b)


@pytest.mark.xfail(
    strict=True,
    reason="#4, #5: written .inp fails in SWMM (ERROR 203 {pumps}; ERROR 115 adverse slope)",
)
def test_written_inp_runs_in_swmm(dendritic_inp, tmp_path):
    assert run_swmm(dendritic_inp, tmp_path) == []


@pytest.mark.xfail(
    strict=True, reason="#3: flow splits are written without validation (KINWAVE ERROR 133)"
)
def test_split_rejected_under_kinematic_wave(network, tmp_path):
    sub = subnetwork_with_edges(network, DENDRITIC + OVERFLOW, outfalls=("OUT", "OV2"))
    with pytest.raises(ValueError, match="J2"):
        write_inp(sub, tmp_path)


@pytest.mark.xfail(strict=True, reason="#3: each outfall node is written once per outfall")
def test_multiple_outfalls_written_once(network, tmp_path):
    sub = subnetwork_with_edges(network, DENDRITIC + OVERFLOW, outfalls=("OUT", "OV2"))
    try:
        path = write_inp(sub, tmp_path)
    except (
        ValueError
    ):  # once #3 rejects the split under KINWAVE, the outfall rows are checked elsewhere
        pytest.skip("split rejected before outfalls are written")
    names = [row[0] for row in read_sections(path)["OUTFALLS"]]
    assert sorted(names) == ["OUT", "OV2"]


@pytest.mark.xfail(strict=True, reason="#3: dead ends are written as junctions with no outlet")
def test_dead_end_rejected(network, tmp_path):
    sub = subnetwork_with_edges(
        network, [("A1", "A2"), ("A2", "J1")]
    )  # J1 has no outlet; OUT absent
    with pytest.raises(ValueError, match="J1"):
        write_inp(sub, tmp_path)
