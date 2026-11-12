#!/usr/bin/env python3
"""Generate signed constraint-Jacobian finite-difference references (SciPy 1.17)."""

from pathlib import Path
import numpy as np
from scipy.optimize._numdiff import approx_derivative
from scipy.sparse import eye

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "test" / "optimize_trust_constr_constraint_rel_step_reference.f90"


def source_array(v):
    return ", ".join(f"{float(y):.17e}_dp" for y in v)


def case(label, x0, relative_step, scheme, colored=False, bounds=None):
    probes = []

    def fun(x):
        probes.append(x.copy())
        return x * x

    kwargs = {}
    if colored:
        kwargs["sparsity"] = eye(2, format="csr")
    if bounds is not None:
        kwargs["bounds"] = bounds
    jac = approx_derivative(fun, np.asarray(x0, dtype=float), method=scheme,
                            rel_step=relative_step, **kwargs)
    if hasattr(jac, "toarray"):
        jac = jac.toarray()
    print(label, "diag(J)=", np.diag(jac), "calls=", len(probes))
    return np.diag(jac), len(probes)


cases = [
    ("dense_two", [2, -3], [-0.1, 0.2], "2-point", False, None),
    ("sparse_two", [2, -3], [-0.1, 0.2], "2-point", True, None),
    ("dense_three", [2, -3], [-0.1, 0.2], "3-point", False, None),
    ("sparse_three", [2, -3], [-0.1, 0.2], "3-point", True, None),
    ("bounded_two", [2, -3], [-0.1, 0.2], "2-point", False,
     ([2, -np.inf], [np.inf, -3])),
    ("bounded_sparse_two", [2, -3], [-0.1, 0.2], "2-point", True,
     ([2, -np.inf], [np.inf, -3])),
    ("bounded_sparse_three", [2, -3], [-0.1, 0.2], "3-point", True,
     ([2, -np.inf], [np.inf, -3])),
    ("zero_step_two", [0, -3], [0, 0.2], "2-point", False, None),
    ("small_two", [0.5, -0.5], [0.2, 0.2], "2-point", False, None),
    ("reverse_two", [2, -3], [0.1, -0.2], "2-point", False, None),
]

lines = [
    "! Generated from SciPy 1.17 _numdiff.approx_derivative by",
    "! tools/generate_optimize_trust_constr_constraint_rel_step_reference.py",
    "module test_optimize_trust_constr_constraint_rel_step_reference",
    "    use scifort_kinds, only : dp",
    "    implicit none",
]
for label, *parameters in cases:
    jac, nfev = case(label, *parameters)
    lines.append(f"    real(dp), parameter :: {label}_jac(2) = [{source_array(jac)}]")
    lines.append(f"    integer, parameter :: {label}_nfev = {nfev}")
lines.extend(["end module test_optimize_trust_constr_constraint_rel_step_reference", ""])
OUT.write_text("\n".join(lines), encoding="ascii")
print("Wrote", OUT)
