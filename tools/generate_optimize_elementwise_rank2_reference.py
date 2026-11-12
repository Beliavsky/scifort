"""Generate SciPy 1.17 references for rank-two elementwise Fortran tests."""

from pathlib import Path

import numpy as np
from scipy.optimize import elementwise


def fortran_array(name: str, array: np.ndarray) -> str:
    """Return Fortran rank-two parameter data in column-major order."""
    array = np.asarray(array)
    shape = array.shape
    assert len(shape) == 2
    flat = array.ravel(order="F")
    if array.dtype.kind in "iu":
        values = ", ".join(str(int(value)) for value in flat)
        prefix = "integer"
    else:
        values = ", ".join(f"{float(value):.16e}_dp" for value in flat)
        prefix = "real(dp)"
    tokens = values.split(", ")
    lines = [", ".join(tokens[i:i+2]) for i in range(0, len(tokens), 2)]
    contents = ", &\n            ".join(lines)
    return (
        f"    {prefix}, parameter :: {name}({shape[0]},{shape[1]}) = &\n"
        f"        reshape([ &\n            {contents}], [{shape[0]}, {shape[1]}])\n"
    )


def generate() -> None:
    """Write root and minimum bracketing/refinement array references."""
    root = lambda x: x * x - 4
    minimum = lambda x: (x - 1.5) ** 2
    starts = np.array([[-4.0, 0.0, 3.0], [-4.0, 0.0, 3.0]])
    right = np.array([[-3.0, 1.0, 3.5], [-3.0, 1.0, 3.5]])
    lower = np.array([[-10.0, 0.0, 1.0], [-10.0, 0.0, 1.0]])
    upper = np.array([[-1.0, 10.0, 4.0], [-1.0, 10.0, 4.0]])
    br = elementwise.bracket_root(
        root, starts, xr0=right, xmin=lower, xmax=upper, factor=2.0, maxiter=15
    )
    root_left = np.array([[-4.0, 1.0, 3.0], [4.0, -4.0, 1.0]])
    root_right = np.array([[-1.0, 4.0, 4.0], [1.0, -1.0, 4.0]])
    rr = elementwise.find_root(root, (root_left, root_right), tolerances={"xatol": 1e-12})
    mid = np.array([[-3.0, 0.5, 1.5], [-3.0, 0.5, 1.5]])
    bm = elementwise.bracket_minimum(
        minimum, mid, xl0=mid - 0.3, xr0=mid + 0.3,
        xmin=-10.0, xmax=10.0
    )
    min_left = np.array([[0.0, -3.0, 3.0], [3.0, 0.0, -3.0]])
    min_mid = np.array([[1.0, 0.0, 4.0], [1.0, 1.0, 0.0]])
    min_right = np.array([[3.0, 4.0, 5.0], [0.0, 3.0, 4.0]])
    rm = elementwise.find_minimum(
        minimum, (min_left, min_mid, min_right), tolerances={"xatol": 1e-9}
    )
    parts = [
        "! Generated from scipy.optimize.elementwise, SciPy 1.17.0; do not edit.\n",
        "! SPDX-License-Identifier: MIT\n",
        "module optimize_elementwise_rank2_reference\n",
        "    use scifort_kinds, only : dp\n",
        "    implicit none\n",
    ]
    for prefix, result, points in [
        ("br_root", br, {"xl": br.bracket[0], "xr": br.bracket[1]}),
        ("find_root", rr, {"x": np.where(rr.success, rr.x, 0.0)}),
        ("br_min", bm, {"xm": bm.bracket[1]}),
        ("find_min", rm, {"x": np.where(rm.success, rm.x, 0.0)}),
    ]:
        for key in ("status", "nit", "nfev"):
            parts.append(fortran_array(f"{prefix}_{key}_ref", getattr(result, key)))
        for key, value in points.items():
            parts.append(fortran_array(f"{prefix}_{key}_ref", value))
    parts.append("end module optimize_elementwise_rank2_reference\n")
    output = Path(__file__).resolve().parents[1] / "test" / "optimize_elementwise_rank2_reference.f90"
    output.write_text("".join(parts), encoding="ascii")
    print(f"Generated {output} from SciPy 1.17 rank-two reference cases")


if __name__ == "__main__":
    generate()
