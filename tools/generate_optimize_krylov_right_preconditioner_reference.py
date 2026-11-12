#!/usr/bin/env python3
"""Reproduce TFQMR right-preconditioning reference solves using SciPy 1.17.0.

SciPy's tfqmr() has a left-preconditioning callback, not a right-side switch.
To check the SciFort extension without depending on its implementation,
solve A*M*y=b with SciPy's unpreconditioned TFQMR and then set x=M*y.
The matrix-free operator counts below exclude SciFort's optional final check
of the original residual A*x-b (one additional operator application).
"""

import numpy as np
import scipy
from scipy.sparse.linalg import LinearOperator, tfqmr


def make_problem(n):
    """Create a reproducible nonsymmetric matrix and known solution."""
    a = np.zeros((n, n))
    target = np.sin(0.3 * np.arange(1, n + 1, dtype=float))
    for j in range(n):
        a[j, j] = 3.0 + 0.1 * (j + 1)
    for j in range(n - 1):
        a[j, j + 1] = 0.22
        a[j + 1, j] = 0.37
    for j in range(n - 2):
        a[j, j + 2] = 0.10
    return a, a @ target, target


def make_preconditioner(a, kind):
    """Return a diagonal or upper-bidiagonal inverse approximation."""
    diagonal = np.diag(a)

    def apply(v):
        product = v / diagonal
        if kind == "triangular":
            product[:-1] -= 0.06 * v[1:] / diagonal[1:]
        return product

    return apply


def run_reference(n, kind, expected_calls):
    """Check transformed-system accuracy and deterministic operator count."""
    a, b, target = make_problem(n)
    apply_m = make_preconditioner(a, kind)
    calls = 0

    def transformed_matvec(y):
        nonlocal calls
        calls += 1
        return a @ apply_m(y)

    transformed_a = LinearOperator((n, n), matvec=transformed_matvec, dtype=float)
    y, info = tfqmr(transformed_a, b, rtol=1e-10, atol=0.0, maxiter=120)
    solution = apply_m(y)
    relres = np.linalg.norm(b - a @ solution) / np.linalg.norm(b)
    assert info == 0, (n, kind, info)
    assert calls == expected_calls, (n, kind, calls, expected_calls)
    assert relres < 1e-10, (n, kind, relres)
    np.testing.assert_allclose(solution, target, atol=1e-8, rtol=0)
    print(f"n={n:2d} {kind:10s} right: matvecs={calls:2d}, true residual={relres:.3e}")

    # SciPy's direct M= path is the backward-compatible SciFort default.
    if kind == "diagonal":
        direct_a = LinearOperator((n, n), matvec=lambda v: a @ v, dtype=float)
        m_op = LinearOperator((n, n), matvec=apply_m, dtype=float)
        legacy_x, legacy_info = tfqmr(
            direct_a, b, M=m_op, rtol=1e-10, atol=0.0, maxiter=120
        )
        legacy_relres = np.linalg.norm(b - a @ legacy_x) / np.linalg.norm(b)
        assert legacy_info == 0, (n, legacy_info)
        assert legacy_relres > 1e-5, (n, legacy_relres)
        print(f"n={n:2d} diagonal   SciPy M= true residual={legacy_relres:.3e}")


def main():
    """Reproduce parity cases used by the Fortran test suite."""
    print(f"SciPy {scipy.__version__}")
    run_reference(16, "diagonal", 10)
    run_reference(16, "triangular", 11)
    run_reference(64, "diagonal", 10)
    run_reference(64, "triangular", 10)


if __name__ == "__main__":
    main()
