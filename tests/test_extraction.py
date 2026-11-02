import os
import subprocess
import sys

import pytest

from conftest import FIXTURE, REPO

PROBE = """
import contextlib, io
from sewernetworkx import Network, SubNetwork
with contextlib.redirect_stdout(io.StringIO()):
    n = Network({fixture!r})
    s = SubNetwork({{"network_data": n, "starting_node": "J1", "outfall_node": ["OUT"], "depth": 2,
                    "outfall_node_type": ["FREE"], "scale": 1.0}})
    s.create_sub_network()
print(",".join(sorted(s.subgraph.nodes)))
"""


def extract_with_seed(seed):
    env = {**os.environ, "PYTHONHASHSEED": str(seed), "PYTHONPATH": str(REPO)}
    out = subprocess.run([sys.executable, "-W", "ignore", "-c", PROBE.format(fixture=str(FIXTURE))],
                         env=env, capture_output=True, text=True, check=True)
    return out.stdout.strip().splitlines()[-1]


@pytest.mark.xfail(strict=True, reason="#2: create_sub_network depends on set iteration order (hash seed)")
def test_extraction_is_deterministic_across_hash_seeds():
    results = {extract_with_seed(seed) for seed in range(8)}
    assert len(results) == 1, f"{len(results)} different node sets: {sorted(results)}"
