#!/usr/bin/env python3
"""Generate SciPy 1.17.0 and high-precision tail references for four discrete laws."""
from pathlib import Path
from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb
import numpy as np
import scipy
from scipy.stats import betabinom, hypergeom, nhypergeom, boltzmann

if scipy.__version__ != "1.17.0":
    raise SystemExit(f"expected SciPy 1.17.0, found {scipy.__version__}")

PROBS = np.array([1e-6, 0.05, 0.5, 0.95, 1-1e-6], dtype=float)
FAMILIES = [
    ("bb", betabinom, np.array([0., 2., 5., 8., 11.]), (12, 2.5, 4.0)),
    ("hg", hypergeom, np.array([0., 2., 4., 7., 10.]), (40, 15, 12)),
    ("nh", nhypergeom, np.array([0., 2., 4., 7., 12.]), (40, 15, 8)),
    ("bo", boltzmann, np.array([0., 1., 3., 7., 10.]), (0.35, 12)),
]

def fval(x: float) -> str:
    return f"{float(x):.17e}_dp"

def emit_array(name: str, values: np.ndarray) -> list[str]:
    values = np.asarray(values)
    lines = [f"    real(dp), parameter, public :: {name}({values.size}) = [ &"]
    for i, value in enumerate(values):
        suffix = ", &" if i + 1 < values.size else "]"
        lines.append(f"        {fval(value)}{suffix}")
    return lines

lines = [
    "! SPDX-License-Identifier: MIT",
    "! Generated from SciPy 1.17.0, with exact/high-precision formulas for two stress tails.",
    "module test_betabinom_hypergeom_nhypergeom_boltzmann_reference",
    "    use scifort_kinds, only : dp",
    "    implicit none",
    "    private",
]
lines += emit_array("reference_probs", PROBS)
for prefix, dist, x, args in FAMILIES:
    lines += emit_array(f"{prefix}_x", x)
    for method in ("pmf", "logpmf", "cdf", "sf", "logcdf", "logsf"):
        lines += emit_array(f"{prefix}_{method}", getattr(dist, method)(x, *args))
    lines += emit_array(f"{prefix}_ppf", dist.ppf(PROBS, *args))
    lines += emit_array(f"{prefix}_isf", dist.isf(PROBS, *args))

lines += emit_array("bb_stress_sf", betabinom.sf(np.array([150.0]), 200, 0.2, 3.5))
lines += emit_array("bb_stress_isf", betabinom.isf(np.array([1e-10]), 200, 0.2, 3.5))
lines += emit_array("hg_stress_sf", hypergeom.sf(np.array([140.0]), 10000, 2500, 400))
lines += emit_array("hg_stress_ppf", hypergeom.ppf(np.array([1-1e-10]), 10000, 2500, 400))
# The two following SF stress values use direct high-precision defining formulas.
# SciPy 1.17.0 loses several digits in these regimes.
nh_num = Fraction(0, 1)
nh_den = comb(500, 120)
for j in range(61, 121):
    nh_num += Fraction(comb(j + 99, j) * comb(400 - j, 120 - j), nh_den)
lines += emit_array("nh_stress_sf", np.array([float(nh_num)]))
lines += emit_array("nh_stress_isf", nhypergeom.isf(np.array([1e-10]), 500, 120, 100))
lines += emit_array("bo_stress_ppf", boltzmann.ppf(np.array([0.999999]), 1e-8, 100000))
getcontext().prec = 80
lam = Decimal("1e-8")
bo_sf = (-lam * Decimal(99991)).exp() * (Decimal(1) - (-lam * Decimal(9)).exp()) / \
    (Decimal(1) - (-lam * Decimal(100000)).exp())
lines += emit_array("bo_stress_sf", np.array([float(bo_sf)]))
lines.append("end module test_betabinom_hypergeom_nhypergeom_boltzmann_reference")
lines.append("")

out = Path(__file__).resolve().parents[1] / "test" / "betabinom_hypergeom_nhypergeom_boltzmann_reference.f90"
out.write_text("\n".join(lines), encoding="ascii")
