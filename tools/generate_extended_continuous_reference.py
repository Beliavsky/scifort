#!/usr/bin/env python3
"""Generate SciPy 1.17.0 references for Milestone 7 continuous families."""
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import genpareto, gumbel_l, gumbel_r, powerlaw, triang

if scipy.__version__ != "1.17.0":
    raise SystemExit(f"expected SciPy 1.17.0, got {scipy.__version__}")

families = {
    "gr": (gumbel_r, np.array([-12.0, -2.0, 0.0, 1.5, 10.0, 40.0]), None),
    "gl": (gumbel_l, np.array([-40.0, -10.0, -1.5, 0.0, 2.0, 12.0]), None),
    "pw": (powerlaw, np.array([1e-12, 0.01, 0.2, 0.7, 0.999999]),
           np.array([0.3, 0.7, 1.0, 2.5, 5.0])),
    "tr": (triang, np.array([0.1, 0.05, 0.5, 0.9, 0.7]),
           np.array([0.0, 0.2, 0.5, 0.8, 1.0])),
    "gp": (genpareto, np.array([0.2, 1.0, 1.3, 2.0, 5.0]),
           np.array([-2.0, -0.5, 0.0, 0.2, 2.0])),
}

def fmt(v):
    if np.isposinf(v): return "huge(1.0_dp)"
    if np.isneginf(v): return "-huge(1.0_dp)"
    return f"{float(v):.17e}_dp".replace("e+", "e+").replace("e-", "e-")

def arr(name, vals):
    body = ", &\n        ".join(fmt(v) for v in vals)
    return f"    real(dp), parameter :: {name}({len(vals)}) = [ &\n        {body} ]\n"

out = [
"! SPDX-License-Identifier: MIT\n",
"! Copyright (c) 2026 SciFort contributors\n\n",
"module test_extended_continuous_reference\n",
"    use scifort_kinds, only : dp\n",
"    implicit none\n",
"    private\n\n",
]
for prefix, (dist, xs, shapes) in families.items():
    out.append(f"    public :: {prefix}_x, {prefix}_pdf, {prefix}_cdf, {prefix}_sf, &\n")
    out.append(f"        {prefix}_logcdf, {prefix}_logsf")
    if shapes is not None:
        out.append(f", {prefix}_shape")
    out.append("\n")
    out.append(arr(f"{prefix}_x", xs))
    if shapes is None:
        args = ()
        pdf = dist.pdf(xs)
        cdf = dist.cdf(xs)
        sf = dist.sf(xs)
        lc = dist.logcdf(xs)
        ls = dist.logsf(xs)
    else:
        out.append(arr(f"{prefix}_shape", shapes))
        pdf = np.array([dist.pdf(x, s) for x,s in zip(xs, shapes)])
        cdf = np.array([dist.cdf(x, s) for x,s in zip(xs, shapes)])
        sf = np.array([dist.sf(x, s) for x,s in zip(xs, shapes)])
        lc = np.array([dist.logcdf(x, s) for x,s in zip(xs, shapes)])
        ls = np.array([dist.logsf(x, s) for x,s in zip(xs, shapes)])
    for suffix, vals in [("pdf",pdf),("cdf",cdf),("sf",sf),("logcdf",lc),("logsf",ls)]:
        out.append(arr(f"{prefix}_{suffix}", vals))
    out.append("\n")
out.append("end module test_extended_continuous_reference\n")
Path("test/extended_continuous_reference.f90").write_text("".join(out), encoding="ascii")
print("generated with scipy", scipy.__version__)
