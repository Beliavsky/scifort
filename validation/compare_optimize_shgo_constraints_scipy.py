#!/usr/bin/env python3
"""Compare bounded constrained SHGO objective/feasibility against SciPy.

The Fortran driver supplies exact known-solution polishers, *not* a mock of
SciFort's public SLSQP adapter, which needs full fpm integration testing.
"""
import subprocess
import sys

import numpy as np
from scipy.optimize import shgo, LinearConstraint, NonlinearConstraint


def run_case(number):
    proc = subprocess.run([sys.argv[1], str(number)], text=True,
                          capture_output=True, check=True)
    lines = proc.stdout.strip().splitlines()
    values = lines[-2].split()
    index, status = map(int, values[:2])
    f, cv = map(float, values[2:4])
    nfev, ncev, nlocal = map(int, values[4:7])
    x = np.fromstring(lines[-1], sep=" ")
    assert index == number and status == 0
    assert np.isfinite(f) and np.isfinite(cv) and cv <= 1.01e-8
    assert 0 < nfev <= (300 if number == 6 else 200)
    assert ncev > 0 and nlocal > 0
    return f, x


def main():
    cases = []
    sq = [(0.0, 1.0)] * 2
    cases.append((lambda x: (x[0] - 0.2) ** 2 + (x[1] - 0.1) ** 2,
                  sq, [LinearConstraint([[1, 1]], 0.9, np.inf)],
                  0.18))
    cases.append((lambda x: (x[0] - 0.55) ** 2 + (x[1] - 0.15) ** 2,
                  sq, [LinearConstraint([[1, 1]], 0.7, 0.7)], 0.0))
    cases.append((lambda x: np.sum((x - 0.8) ** 2),
                  sq, [NonlinearConstraint(lambda x: np.dot(x, x), -np.inf, 0.4)],
                  2 * (0.8 - np.sqrt(0.2)) ** 2))
    cases.append((lambda x: (x[0] - 0.3) ** 2 + (x[1] - 0.6) ** 2,
                  sq, [NonlinearConstraint(lambda x: np.dot(x, x), 0.5, 0.5)],
                  (np.sqrt(0.5) - np.sqrt(0.45)) ** 2))
    cases.append((lambda x: (x[0] - 0.55) ** 2 + (x[1] - 0.15) ** 2,
                  sq, [LinearConstraint([[1, 1]], 0.7, 0.7),
                       NonlinearConstraint(lambda x: np.dot(x, x), -np.inf, 0.5)], 0.0))
    cases.append((lambda x: np.sum((x - [0.4, 0.5, 0.6]) ** 2),
                  [(0.0, 1.0)] * 3,
                  [LinearConstraint([[1, 1, 1]], 1.5, 1.5)], 0.0))
    compared = 0
    for i, (fun, bounds, con, reference) in enumerate(cases, start=1):
        f, x = run_case(i)
        assert abs(f - fun(x)) < 2e-11, (i, f, fun(x))
        assert abs(f - reference) < 2e-8, (i, f, reference)
        scipy_res = shgo(fun, bounds, iters=2, constraints=con,
                         minimizer_kwargs={"method": "SLSQP",
                                           "options": {"ftol": 1e-10}},
                         sampling_method="simplicial")
        if scipy_res.success:
            assert abs(f - scipy_res.fun) < 2e-5, (i, f, scipy_res.fun)
            compared += 1
        else:
            print(f"case {i}: SciPy returned {scipy_res.message!r}; "
                  "known analytic optimum and feasibility verified")
    print(f"6 independent constrained analytic/reference checks passed; "
          f"{compared} successful SciPy SHGO comparisons")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: compare_optimize_shgo_constraints_scipy.py DRIVER")
    main()
