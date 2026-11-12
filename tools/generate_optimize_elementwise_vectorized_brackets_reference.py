#!/usr/bin/env python3
"""Generate rank-one bracketing reference results with SciPy 1.17."""
from pathlib import Path

import numpy as np
from scipy.optimize._elementwise import bracket_minimum, bracket_root

xl0 = np.array([0.0, 0.0, 2.0, -3.0, 6.0, 0.0, 1.0])
xr0 = xl0 + 0.5
xmin = np.full(7, -10.0)
xmax_root = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 0.75, 10.0])
factors = np.array([2.0, 1.5, 2.0, 2.0, 2.0, 2.0, 2.0])
root = bracket_root(
    lambda x: x - 1.0,
    xl0, xr0=xr0, xmin=xmin, xmax=xmax_root,
    factor=factors, maxiter=30,
)

xm0 = np.array([0.0, 1.5, 5.0, -4.0, 9.5, 7.0, 1.0])
xmax_min = np.array([10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 1.1])
minimum = bracket_minimum(
    lambda x: (x - 1.5)**2 + 0.02 * (x - 1.5)**4,
    xm0, xmin=xmin, xmax=xmax_min, maxiter=25,
)


def line(name: str, values: np.ndarray) -> str:
    numbers = ", ".join(str(int(v)) for v in values)
    return f"    integer, parameter :: {name}(7) = [{numbers}]"


text = "\n".join([
    "! SPDX-License-Identifier: MIT",
    "! Regenerate: python tools/generate_optimize_elementwise_vectorized_brackets_reference.py",
    "! SciPy 1.17.0 scipy.optimize._elementwise bracket_root/bracket_minimum",
    "module optimize_elementwise_vectorized_brackets_reference",
    "    implicit none",
    "    private",
    "    public :: vector_bracket_root_status, vector_bracket_root_nit, vector_bracket_root_nfev",
    "    public :: vector_bracket_min_status, vector_bracket_min_nit, vector_bracket_min_nfev",
    line("vector_bracket_root_status", root.status),
    line("vector_bracket_root_nit", root.nit),
    line("vector_bracket_root_nfev", root.nfev),
    line("vector_bracket_min_status", minimum.status),
    line("vector_bracket_min_nit", minimum.nit),
    line("vector_bracket_min_nfev", minimum.nfev),
    "end module optimize_elementwise_vectorized_brackets_reference",
    "",
])
path = (Path(__file__).resolve().parents[1] / "test" /
        "optimize_elementwise_vectorized_brackets_reference.f90")
path.write_text(text)
print(f"Wrote {path}")
