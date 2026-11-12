#!/usr/bin/env python3
"""Check SHGO feasibility restoration against independent SciPy SLSQP solutions.

Restoration is not local minimization. Accordingly, we compare feasible points,
true objective values, and lower bounds from numerical constrained minimization,
not claimed equivalence of final objective values.
"""
import subprocess
import sys

import numpy as np
from scipy.optimize import minimize


def run_native(exe, i):
    out = subprocess.run([exe, str(i)], check=True, capture_output=True, text=True)
    lines = out.stdout.strip().splitlines()
    fields = lines[-2].split()
    assert int(fields[0]) == i
    status = int(fields[1])
    success = fields[2] == "T"
    fun, cv = map(float, fields[3:5])
    nfev, ncev, ncjev, nrestore, nrestored = map(int, fields[5:10])
    x = np.fromstring(lines[-1], sep=" ")
    return status, success, fun, cv, nfev, ncev, ncjev, nrestore, nrestored, x


def main(exe):
    obj = lambda x: np.sum((x - 0.25)**2)
    funcs = [lambda x: np.sum(x), lambda x: np.sum(x),
             lambda x: x @ x, lambda x: np.prod(x), lambda x: x @ x]
    targets = [0.7, 0.7, 0.13, 0.11, 0.13]
    compared = 0
    for i in range(1, 6):
        status, success, f, cv, nf, nc, nj, nr, nok, x = run_native(exe, i)
        assert status == 0 and success and np.all(np.isfinite(x)), (i, status, x)
        assert x.shape == (2,) and np.all((0.0 <= x) & (x <= 1.0))
        assert abs(funcs[i-1](x) - targets[i-1]) <= 1.01e-8, (i, x)
        assert abs(f - obj(x)) < 1e-12 and cv <= 1.01e-8
        assert 0 < nf <= 200 and nc > 0 and nr > 0 and nok > 0
        assert (nj == 0) if i in (1, 5) else (nj > 0)
        ref = minimize(obj, [0.4, 0.4], method="SLSQP", bounds=[(0, 1)]*2,
                       constraints=[{"type": "eq", "fun": lambda z, fn=funcs[i-1], t=targets[i-1]: fn(z)-t}],
                       options={"ftol": 1e-12, "maxiter": 200})
        assert ref.success and abs(funcs[i-1](ref.x)-targets[i-1]) < 1e-8, (i, ref.message)
        assert f >= ref.fun - 3e-6, (i, f, ref.fun)
        assert f < ref.fun + 0.55, (i, f, ref.fun)
        compared += 1
    status, success, f, cv, nf, nc, nj, nr, nok, x = run_native(exe, 6)
    assert status == 6 and not success and nf == 0 and nok == 0 and cv >= 1.0
    print(f"{compared} SciPy SLSQP feasibility/objective comparisons passed; "
          "impossible equality correctly rejected")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_optimize_shgo_restoration_scipy.py NATIVE_EXE")
    main(sys.argv[1])
