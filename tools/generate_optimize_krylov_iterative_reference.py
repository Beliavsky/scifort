#!/usr/bin/env python3
"""Reproduce SciPy 1.17.0 matrix-free inner Krylov work-count references."""

import numpy as np
from scipy import __version__ as scipy_version
from scipy.optimize import root
from scipy.sparse.linalg import LinearOperator, bicgstab, cgs, minres, tfqmr


LINEAR = np.array(
    [[4.0, 1.0, 0.0, 0.0],
     [0.0, 3.0, 1.0, 0.0],
     [0.0, 0.0, 2.0, 1.0],
     [1.0, 0.0, 0.0, 5.0]]
)
SYMMETRIC = np.array(
    [[4.0, 1.0, 0.0, 0.0],
     [1.0, -3.0, 1.0, 0.0],
     [0.0, 1.0, 2.0, 1.0],
     [0.0, 0.0, 1.0, -4.0]]
)
TARGET = np.array([1.0, -2.0, 3.0, -1.0])
METHODS = (
    ("bicgstab", bicgstab, 7, 9, 2, 34, 9),
    ("cgs", cgs, 8, 10, 2, 30, 8),
    ("minres", minres, 4, 6, 2, 27, 9),
    ("tfqmr", tfqmr, 7, 9, 2, 37, 9),
)


def solve_reference(method, residual, x0, fatol):
    """Run a SciPy nonlinear root with an explicit inner method."""
    return root(
        residual,
        x0,
        method="krylov",
        options={
            "fatol": fatol,
            "line_search": None,
            "jac_options": {
                "method": method,
                "inner_maxiter": 30,
                "rdiff": 1e-6,
            },
        },
    )


def check_result(label, result, target, nit, nfev):
    """Check work counts and numerical root against the pinned reference."""
    np.testing.assert_allclose(result.x, target, atol=1e-7, rtol=0)
    assert result.success, (label, result.message)
    assert (result.nit, result.nfev) == (nit, nfev), (
        label, result.nit, result.nfev
    )
    print(f"{label:20s}: nit={result.nit:2d}, nfev={result.nfev:2d}")


def main():
    """Generate SciPy linear, nonlinear, and direct matrix-free references."""
    print(f"SciPy {scipy_version}")
    for method, linear_solver, direct_nmatvec, linear_nfev, linear_nit, curved_nfev, curved_nit in METHODS:
        a = SYMMETRIC if method == "minres" else LINEAR
        b = a @ TARGET
        calls = [0]

        def matrix_product(v):
            calls[0] += 1
            return a @ v

        op = LinearOperator((4, 4), matvec=matrix_product, dtype=float)
        solution, info = linear_solver(op, b, rtol=1e-10, maxiter=60)
        assert info == 0, (method, info)
        np.testing.assert_allclose(solution, TARGET, rtol=0, atol=1e-8)
        assert calls[0] == direct_nmatvec, (method, calls[0], direct_nmatvec)
        print(f"{method:8s} direct matvec calls={calls[0]}")

        result = solve_reference(method, lambda x: a @ x - b, np.zeros(4), 1e-7)
        check_result(method + " linear", result, TARGET, linear_nit, linear_nfev)

        result = solve_reference(
            method,
            lambda x: x * x - [4.0, 9.0, 16.0, 25.0],
            np.ones(4),
            1e-8,
        )
        check_result(
            method + " nonlinear", result, [2.0, 3.0, 4.0, 5.0],
            curved_nit, curved_nfev,
        )


if __name__ == "__main__":
    main()
