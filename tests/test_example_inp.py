from conftest import EXAMPLE_INP, run_swmm


def test_example_inp_runs_in_swmm(tmp_path):
    assert run_swmm(EXAMPLE_INP, tmp_path) == []
