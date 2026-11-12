#!/usr/bin/env python3
"""Reproduce the SciPy 1.17.0 inner_M / adaptive Krylov reference cases."""

import numpy as np
from scipy import __version__ as scipy_version
from scipy.optimize import root
from scipy.sparse.linalg import LinearOperator


def check(label, result, expected_x, expected_nit, expected_nfev):
    """Check and display reproducible nonlinear iteration work counts."""
    assert result.success, (label, result.message)
    np.testing.assert_allclose(result.x, expected_x, rtol=0, atol=1e-7)
    assert (result.nit, result.nfev) == (expected_nit, expected_nfev), (
        label, result.nit, result.nfev
    )
    print(f"{label}: nit={result.nit}, nfev={result.nfev}, x={result.x}")


class DiagonalInverse(LinearOperator):
    """Preconditioner for a diagonally scaled linear system."""

    def __init__(self, diagonal):
        self.diagonal = np.asarray(diagonal, dtype=np.float64)
        super().__init__(dtype=np.float64, shape=(len(diagonal), len(diagonal)))

    def _matvec(self, vector):
        return vector / self.diagonal


class AdaptiveInverse(LinearOperator):
    """Update the approximate inverse diagonal after accepted Newton steps."""

    def __init__(self, x):
        self.x = np.asarray(x, dtype=np.float64)
        self.updates = 0
        super().__init__(dtype=np.float64, shape=(len(x), len(x)))

    def _matvec(self, vector):
        return vector / (2 * self.x)

    def update(self, x, f):
        self.x = x.copy()
        self.updates += 1


def main():
    """Generate the matching root/GMRES/LGMRES reference results."""
    print(f"SciPy {scipy_version}")
    diagonal = np.array([1.0, 10.0, 100.0, 1000.0])
    rhs = np.array([2.0, -3.0, 4.0, -5.0])
    for method, expected_nfev in (("gmres", 4), ("lgmres", 3)):
        result = root(
            lambda x: diagonal * x - rhs,
            np.zeros(4),
            method="krylov",
            options={
                "fatol": 1e-7,
                "line_search": None,
                "jac_options": {
                    "method": method,
                    "inner_maxiter": 1,
                    "rdiff": 1e-6,
                    "inner_M": DiagonalInverse(diagonal),
                },
            },
        )
        check(method, result, rhs / diagonal, 2, expected_nfev)

    initial_x = np.ones(2)
    preconditioner = AdaptiveInverse(initial_x)
    result = root(
        lambda x: x * x - [4.0, 9.0],
        initial_x,
        method="krylov",
        options={
            "fatol": 1e-8,
            "line_search": None,
            "jac_options": {
                "rdiff": 1e-6,
                "method": "lgmres",
                "inner_maxiter": 2,
                "outer_k": 1,
                "inner_M": preconditioner,
            },
        },
    )
    check("adaptive lgmres", result, [2.0, 3.0], 7, 16)
    assert preconditioner.updates == 6
    print(f"adaptive preconditioner update calls={preconditioner.updates}")


if __name__ == "__main__":
    main()
