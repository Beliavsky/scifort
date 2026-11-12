#!/usr/bin/env python3
"""Reproduce the SciPy TFQMR preconditioned-residual caveat.

SciFort's opt-in verify_residual option is not part of SciPy: the checks here
pin the unchanged SciPy-compatible default on the 16-dimensional test system.
"""

import numpy as np
from scipy import __version__ as scipy_version
from scipy.sparse.linalg import LinearOperator, tfqmr


def main():
    """Compare approximate-diagonal and identity TFQMR inner solves."""
    n = 16
    a = np.zeros((n, n))
    for j in range(n):
        a[j, j] = 3.0 + 0.1 * (j + 1)
        if j + 1 < n:
            a[j, j + 1] = 0.22
            a[j + 1, j] = 0.37
        if j + 2 < n:
            a[j, j + 2] = 0.10
    target = np.sin(np.arange(1, n + 1, dtype=float) * 0.3)
    rhs = a @ target
    operator = LinearOperator((n, n), matvec=lambda x: a @ x, dtype=float)
    inverse_diag = LinearOperator(
        (n, n), matvec=lambda x: x / np.diag(a), dtype=float
    )

    print(f"SciPy {scipy_version}")
    for name, preconditioner in (("identity", None), ("diagonal", inverse_diag)):
        sol, info = tfqmr(
            operator, rhs, rtol=1.0e-10, maxiter=60, M=preconditioner
        )
        true_relative_residual = np.linalg.norm(rhs - a @ sol) / np.linalg.norm(rhs)
        print(f"{name:9s}: info={info}, true relative residual={true_relative_residual:.8e}")
        assert info == 0, (name, info)
        if preconditioner is None:
            assert true_relative_residual < 1.0e-10
        else:
            assert true_relative_residual > 1.0e-5


if __name__ == "__main__":
    main()
