#!/usr/bin/env python3
"""Compare multi-block constrained SHGO with independently checked SciPy solutions.

Cases 1, 3, 4, and 5 use scipy.optimize.shgo. Case 2 uses SLSQP
because the installed SciPy SHGO can signal success for an infeasible
intersection of a linear and a curved equality; each comparison verifies
feasibility explicitly rather than accepting scipy.success on its own.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
from scipy.optimize import NonlinearConstraint, minimize, shgo


def constraint_violation(x: np.ndarray, objects: list[NonlinearConstraint]) -> float:
    worst = 0.0
    for con in objects:
        v = np.asarray(con.fun(x), dtype=float)
        lb = np.broadcast_to(con.lb, v.shape)
        ub = np.broadcast_to(con.ub, v.shape)
        worst = max(worst, float(np.max(np.maximum(lb - v, 0.0))),
                    float(np.max(np.maximum(v - ub, 0.0))))
    return worst


def main(path: str) -> None:
    f = lambda x: float(np.sum((x - 0.25) ** 2))
    csum = NonlinearConstraint(lambda x: np.array([np.sum(x)]), [0.7], [0.7],
                               jac=lambda x: np.ones((1, len(x))))
    cdiff = NonlinearConstraint(lambda x: np.array([x[0] - x[1]]), [0.1], [0.1])
    cradius = NonlinearConstraint(lambda x: np.array([np.sum(x * x)]), [0.25], [0.25])
    cfirst2 = NonlinearConstraint(lambda x: np.array([x[0]]), [0.4], [0.4],
                                  jac=lambda x: np.array([[1.0, 0.0]]))
    cfirst3 = NonlinearConstraint(lambda x: np.array([x[0]]), [0.2], [0.2],
                                  jac=lambda x: np.array([[1.0, 0.0, 0.0]]))
    cdiff23 = NonlinearConstraint(lambda x: np.array([x[1] - x[2]]), [0.1], [0.1])
    cvector = NonlinearConstraint(lambda x: np.array([x[0]*x[1], np.sum(x*x)]),
                                  [0.05, -np.inf], [np.inf, 0.35])
    cases = [
        ([csum, cdiff], 2, "shgo"),
        ([csum, cradius], 2, "slsqp"),
        ([csum, cfirst2], 2, "shgo"),
        ([csum, cfirst3, cdiff23], 3, "shgo"),
        ([csum, cdiff, cvector], 2, "shgo"),
    ]
    for case_num, (constraints, nvar, comparator) in enumerate(cases, 1):
        out = subprocess.check_output([path], input=f"{case_num}\n", text=True)
        fields = out.split()
        assert int(fields[0]) == case_num, (case_num, out)
        success = fields[1] == "T"
        objective = float(fields[2])
        maxcv = float(fields[3])
        nfev, ncjev = (int(fields[4]), int(fields[5]))
        x = np.array([float(z) for z in fields[6:]])
        assert success and x.size == nvar and np.all(np.isfinite(x)), (case_num, out)
        assert np.all(x >= -1e-10) and np.all(x <= 1+1e-10), (case_num, out)
        assert maxcv <= 1.1e-8, (case_num, out)
        assert constraint_violation(x, constraints) <= 1.1e-8, (case_num, out)
        assert abs(f(x) - objective) <= 1e-12, (case_num, out)
        assert 1 <= nfev <= 150 and ncjev > 0, (case_num, out)
        if comparator == "shgo":
            # In this specific 3D case Sobol local search recovers the feasible
            # solution while SciPy's simplicial path can return a false success.
            sampling = "sobol" if case_num == 4 else "simplicial"
            scipy_result = shgo(f, [(0.0, 1.0)] * nvar,
                                constraints=constraints, n=64, iters=1,
                                sampling_method=sampling)
        else:
            scipy_result = minimize(f, [0.5, 0.2], bounds=[(0.0, 1.0)] * nvar,
                                    constraints=constraints, method="SLSQP",
                                    options={"ftol": 1e-11, "maxiter": 300})
        assert scipy_result.success, (case_num, scipy_result.message)
        y = np.asarray(scipy_result.x)
        assert constraint_violation(y, constraints) <= 1e-6, (case_num, y)
        assert np.isclose(scipy_result.fun, objective, atol=5e-6, rtol=0), (
            case_num, scipy_result.fun, objective)
        print(f"SHGO multi-constraint SciPy comparison {case_num}: "
              f"ok f={objective:.10f} cv={maxcv:.2e} nfev={nfev}")
    print("SHGO multi-constraint SciPy comparisons passed: 5")


if __name__ == "__main__":
    main(sys.argv[1])
