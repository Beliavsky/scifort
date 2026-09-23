#!/usr/bin/env python3
"""Generate SciPy 1.17.0 references for four closed-form continuous families."""
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import bradford, dweibull, fisk, truncexpon

if scipy.__version__ != "1.17.0":
    raise SystemExit(f"expected SciPy 1.17.0, got {scipy.__version__}")

families = {
    "br": (bradford, 1.7, np.array([0.05, 0.2, 0.5, 0.8, 0.95])),
    "te": (truncexpon, 1.7, np.array([0.05, 0.2, 0.5, 1.0, 1.6])),
    "fi": (fisk, 1.7, np.array([0.1, 0.3, 1.0, 2.5, 10.0])),
    "dw": (dweibull, 1.7, np.array([-3.0, -1.0, -0.2, 0.7, 2.5])),
}
probs = np.array([0.1, 0.5, 0.9])


def fmt(value):
    return f"{float(value):.17e}_dp"


def arr(name, values):
    body = ", &\n        ".join(fmt(v) for v in values)
    return f"    real(dp), parameter :: {name}({len(values)}) = [ &\n        {body} ]\n"


out = [
    "! SPDX-License-Identifier: MIT\n",
    "! Copyright (c) 2026 SciFort contributors\n\n",
    "module test_closed_form_reference\n",
    "    use scifort_kinds, only : dp\n",
    "    implicit none\n",
    "    private\n\n",
    "    public :: reference_shape, reference_probs\n",
    "    public :: br_x, br_pdf, br_cdf, br_sf, br_logcdf, br_logsf, br_ppf, br_isf\n",
    "    public :: te_x, te_pdf, te_cdf, te_sf, te_logcdf, te_logsf, te_ppf, te_isf\n",
    "    public :: fi_x, fi_pdf, fi_cdf, fi_sf, fi_logcdf, fi_logsf, fi_ppf, fi_isf\n",
    "    public :: dw_x, dw_pdf, dw_cdf, dw_sf, dw_logcdf, dw_logsf, dw_ppf, dw_isf\n\n",
    "    real(dp), parameter :: reference_shape = 1.7_dp\n",
    arr("reference_probs", probs),
]

for prefix, (dist, shape, xs) in families.items():
    out.append(arr(f"{prefix}_x", xs))
    values = {
        "pdf": dist.pdf(xs, shape),
        "cdf": dist.cdf(xs, shape),
        "sf": dist.sf(xs, shape),
        "logcdf": dist.logcdf(xs, shape),
        "logsf": dist.logsf(xs, shape),
        "ppf": dist.ppf(probs, shape),
        "isf": dist.isf(probs, shape),
    }
    for suffix, vals in values.items():
        out.append(arr(f"{prefix}_{suffix}", vals))
    out.append("\n")

out.append("end module test_closed_form_reference\n")
Path("test/closed_form_reference.f90").write_text("".join(out), encoding="ascii")
print("generated with scipy", scipy.__version__)
