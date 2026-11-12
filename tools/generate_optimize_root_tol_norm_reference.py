#!/usr/bin/env python3
"""Print reproducible SciPy reference work counts for custom root norms.

Tested against SciPy 1.17.0; no generated Fortran source is needed.
"""

import warnings

import numpy as np
import scipy
from scipy.optimize import root


def affine(x):
    """Two-component affine residual with different component magnitudes."""
    return np.asarray(x) - [0.2, 0.02]


def weighted_norm(v):
    """Anisotropic maximum norm of a two-component vector."""
    return np.max(np.abs(v) * [0.1, 1.0])


def amplified_norm(v):
    """Large norm to test convergence of small Krylov residuals."""
    return 1.0e9 * np.max(np.abs(v))


def compare(label, method, options):
    """Print status, iterations, work counts, and terminal residuals."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        r = root(affine, [0.0, 0.0], method=method, options=options)
    print(
        f"{label:27s} method={method:15s} status={r.status} "
        f"nit={r.nit:2d} nfev={r.nfev:3d} residual={r.fun}"
    )


def main():
    print(f"SciPy {scipy.__version__}")
    compare("broyden1 default", "broyden1", {"fatol": 0.08})
    for method in (
        "broyden1", "broyden2", "anderson", "linearmixing",
        "diagbroyden", "excitingmixing",
    ):
        compare(
            "weighted tolerance", method,
            {"fatol": 0.08, "tol_norm": weighted_norm},
        )
    compare("default step norm", "broyden1", {"fatol": 1.0, "xatol": 0.08})
    compare(
        "weighted step norm", "broyden1",
        {"fatol": 1.0, "xatol": 0.08, "tol_norm": weighted_norm},
    )
    compare("krylov default", "krylov", {"fatol": 0.001})
    compare(
        "krylov amplified", "krylov",
        {"fatol": 0.001, "tol_norm": amplified_norm},
    )


if __name__ == "__main__":
    main()
