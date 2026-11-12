#!/usr/bin/env python3
"""Check HGR early-terminated SciFort objectives against independent SciPy.

Stopping counts and sampling order are intentionally not compared with SciPy:
SciFort's count is a sampled graph-minimum-birth heuristic, not scipy.shgo's
minhgrd option. This verifies feasible results and objective quality only.
"""
from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy.optimize import shgo


def objective(case_id: int, x: np.ndarray) -> float:
    if case_id == 1:
        return float((x[0] ** 2 - 0.25) ** 2)
    if case_id == 2:
        return float((x[0] - 0.3) ** 2)
    if case_id == 3:
        return float(np.dot(x, x))
    if case_id == 4:
        return float(min((x[0] - 0.5) ** 2 + (x[1] + 0.5) ** 2,
                         (x[0] + 0.5) ** 2 + (x[1] - 0.5) ** 2))
    raise ValueError(case_id)


def main(driver: str) -> None:
    lines = subprocess.check_output([driver], text=True).splitlines()
    assert len(lines) == 4, lines
    for case_id, line in enumerate(lines, 1):
        fields = line.split()
        cid, dim, status, nit, nfev = map(int, fields[:5])
        val = float(fields[5])
        hgr, streak, nmesh = map(int, fields[6:9])
        x = np.array([float(v) for v in fields[9:]])
        assert cid == case_id and len(x) == dim
        assert status == 5, (case_id, status)
        assert nit < 5 and nfev < (33 ** dim), (case_id, nit, nfev)
        assert nmesh == nfev and hgr >= 1 and streak >= 1
        assert np.all(np.abs(x) <= 1.0 + 1e-12)
        assert abs(objective(case_id, x) - val) < 1e-12, (case_id, val, x)
        bounds = [(-1.0, 1.0)] * dim
        scipy = shgo(lambda z: objective(case_id, z), bounds,
                     iters=2, sampling_method="simplicial")
        assert scipy.success and np.isfinite(scipy.fun), (case_id, scipy.message)
        tolerance = 3e-3 if case_id == 2 else 1e-10
        assert val <= float(scipy.fun) + tolerance, (
            case_id, val, scipy.fun, tolerance
        )
        print(f"case {case_id}: HGR={val:.7g}, SciPy={scipy.fun:.7g}, "
              f"HGR nfev={nfev}, full-grid nfev={33**dim}")
    print("4 HGR-terminated SciFort/SciPy objective comparisons passed")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_optimize_shgo_hgr_scipy.py DRIVER")
    main(sys.argv[1])
