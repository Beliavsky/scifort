"""Generate SciPy 1.17 rank-one elementwise references for Fortran regressions."""

from pathlib import Path

import numpy as np
from scipy.optimize import elementwise


def fortran_integer(name: str, arr: np.ndarray) -> str:
    values = ", ".join(str(int(x)) for x in np.ravel(arr))
    return f"    integer, parameter :: {name}({np.size(arr)}) = [{values}]\n"


def fortran_real(name: str, arr: np.ndarray) -> str:
    values = ", &\n        ".join(f"{float(x):.16e}_dp" for x in np.ravel(arr))
    return f"    real(dp), parameter :: {name}({np.size(arr)}) = [ &\n        {values}]\n"


def generate() -> None:
    root = lambda x: x * x - 4
    minimum = lambda x: (x - 1.5) ** 2

    br = elementwise.bracket_root(
        root,
        np.array([-4.0, 0.0, 3.0]),
        xr0=np.array([-3.0, 1.0, 3.5]),
        xmin=np.array([-10.0, 0.0, 1.0]),
        xmax=np.array([-1.0, 10.0, 4.0]),
        factor=2.0,
        maxiter=15,
    )
    rr = elementwise.find_root(
        root,
        (np.array([-4.0, 1.0, 3.0, 4.0]), np.array([-1.0, 4.0, 4.0, 1.0])),
        tolerances={"xatol": 1e-12},
    )
    bm = elementwise.bracket_minimum(
        minimum,
        np.array([-3.0, 0.5, 1.5]),
        xl0=np.array([-3.3, 0.2, 1.2]),
        xr0=np.array([-2.7, 0.8, 1.8]),
        xmin=np.array([-10.0, -5.0, -5.0]),
        xmax=np.array([10.0, 10.0, 10.0]),
    )
    rm = elementwise.find_minimum(
        minimum,
        (
            np.array([0.0, -3.0, 3.0, 3.0]),
            np.array([1.0, 0.0, 4.0, 1.0]),
            np.array([3.0, 4.0, 5.0, 0.0]),
        ),
        tolerances={"xatol": 1e-9},
    )
    text = (
        "! Generated from SciPy 1.17.0; do not edit by hand.\n"
        "! SPDX-License-Identifier: MIT\n"
        "module test_optimize_elementwise_arrays_reference\n"
        "    use scifort_kinds, only : dp\n"
        "    implicit none\n"
        "    private\n"
        "    public :: br_root_status_ref, br_root_nit_ref, br_root_nfev_ref, br_root_xl_ref, br_root_xr_ref\n"
        "    public :: root_status_ref, root_nit_ref, root_nfev_ref, root_x_valid_ref\n"
        "    public :: br_min_status_ref, br_min_nit_ref, br_min_nfev_ref, br_min_xm_ref\n"
        "    public :: min_status_ref, min_nit_ref, min_nfev_ref, min_x_valid_ref\n"
        + fortran_integer("br_root_status_ref", br.status)
        + fortran_integer("br_root_nit_ref", br.nit)
        + fortran_integer("br_root_nfev_ref", br.nfev)
        + fortran_real("br_root_xl_ref", br.bracket[0])
        + fortran_real("br_root_xr_ref", br.bracket[1])
        + fortran_integer("root_status_ref", rr.status)
        + fortran_integer("root_nit_ref", rr.nit)
        + fortran_integer("root_nfev_ref", rr.nfev)
        + fortran_real("root_x_valid_ref", rr.x[np.asarray(rr.success)])
        + fortran_integer("br_min_status_ref", bm.status)
        + fortran_integer("br_min_nit_ref", bm.nit)
        + fortran_integer("br_min_nfev_ref", bm.nfev)
        + fortran_real("br_min_xm_ref", bm.bracket[1])
        + fortran_integer("min_status_ref", rm.status)
        + fortran_integer("min_nit_ref", rm.nit)
        + fortran_integer("min_nfev_ref", rm.nfev)
        + fortran_real("min_x_valid_ref", rm.x[np.asarray(rm.success)])
        + "end module test_optimize_elementwise_arrays_reference\n"
    )
    out = Path(__file__).resolve().parents[1] / "test" / "optimize_elementwise_arrays_reference.f90"
    out.write_text(text, encoding="ascii")
    print(f"Generated {out.name} from SciPy rank-one cases")


if __name__ == "__main__":
    generate()
