#!/usr/bin/env python3
"""Check SciPy 1.17 trust-constr callback timing and stop semantics.

Run from the repository root with SciPy installed.  Unlike finite-difference
trajectory regressions, this checks callback observable behavior, which does
not depend on implementation-specific SQP trial steps.
"""

import numpy as np
from scipy import __version__ as scipy_version
from scipy.optimize import Bounds, NonlinearConstraint, minimize


def objective(x):
    return np.sum((x - np.array([1.0, 2.0])) ** 2)


def gradient(x):
    return 2.0 * (x - np.array([1.0, 2.0]))


def hessian(x):
    return 2.0 * np.eye(2)


def check(label, x0, constraints=(), bounds=None):
    visited = []

    def callback(x, state):
        visited.append((state.nit, x.copy()))
        return state.nit >= 2

    result = minimize(
        objective,
        np.asarray(x0, dtype=float),
        method="trust-constr",
        jac=gradient,
        hess=hessian,
        constraints=constraints,
        bounds=bounds,
        callback=callback,
    )
    assert result.status == 3 and not result.success, (label, result.message)
    assert result.nit == 2
    assert [iteration for iteration, _ in visited] == [1, 2]
    assert all(point.shape == (2,) for _, point in visited)
    np.testing.assert_allclose(result.x, visited[-1][1])
    print(f"{label}: status={result.status} nit={result.nit} callback events={len(visited)}")


def main():
    print("SciPy", scipy_version)
    eq = NonlinearConstraint(lambda x: [np.sum(x)], 1.0, 1.0, jac=lambda x: [[1.0, 1.0]])
    bounds = Bounds([0.0, 0.0], [np.inf, np.inf])
    check("unconstrained", [3.0, 4.0])
    check("equality", [0.0, 0.0], constraints=(eq,))
    check("variable bounds", [0.5, 0.4], bounds=bounds)
    check("equality with bounds", [0.5, 0.4], constraints=(eq,), bounds=bounds)


if __name__ == "__main__":
    main()
