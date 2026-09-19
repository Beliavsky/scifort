# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Generate test/continuous_reference.f90 from high-precision mpmath values.

Usage:
    python tools/generate_continuous_reference.py > test/continuous_reference.f90

Generates reference values and tolerances for lognormal, weibull, pareto,
and rayleigh distributions.
"""

import sys
import mpmath as mp

mp.mp.dps = 50
EPS = mp.mpf(2) ** -52
EPS_FACTOR = 20 * EPS

def round_bin64(val_str):
    return mp.mpf(float(val_str))

def fstr(val):
    s = mp.nstr(val, 19, min_fixed=None, max_fixed=None)
    if "e" in s:
        mant, exp = s.split("e")
        if "." not in mant:
            mant += ".0"
        return f"{mant}e{exp}_dp"
    if "." not in s:
        s += ".0"
    return f"{s}_dp"

LOGNORMAL_CASES = [
    ("0.25", "0.5"),
    ("0.25", "1.0"),
    ("0.25", "2.0"),
    ("0.5", "0.1"),
    ("0.5", "1.0"),
    ("0.5", "2.5"),
    ("1.0", "0.01"),
    ("1.0", "0.5"),
    ("1.0", "1.0"),
    ("1.0", "2.0"),
    ("1.0", "10.0"),
    ("2.0", "0.001"),
    ("2.0", "1.0"),
    ("2.0", "50.0"),
    ("3.0", "1.0"),
]

WEIBULL_CASES = [
    ("0.5", "0.01"),
    ("0.5", "0.5"),
    ("0.5", "1.0"),
    ("0.5", "4.0"),
    ("1.0", "0.5"),
    ("1.0", "1.0"),
    ("1.0", "2.0"),
    ("1.5", "0.2"),
    ("1.5", "1.0"),
    ("1.5", "3.0"),
    ("2.0", "0.1"),
    ("2.0", "1.0"),
    ("2.0", "2.5"),
    ("5.0", "0.8"),
    ("5.0", "1.0"),
    ("5.0", "1.2"),
]

PARETO_CASES = [
    ("0.5", "1.05"),
    ("0.5", "2.0"),
    ("0.5", "10.0"),
    ("1.0", "1.1"),
    ("1.0", "2.0"),
    ("1.0", "5.0"),
    ("2.0", "1.01"),
    ("2.0", "1.5"),
    ("2.0", "3.0"),
    ("2.0", "10.0"),
    ("3.5", "1.1"),
    ("3.5", "2.0"),
    ("5.0", "1.05"),
    ("5.0", "1.5"),
    ("10.0", "1.2"),
]

RAYLEIGH_CASES = [
    "0.01",
    "0.2",
    "0.5",
    "1.0",
    "1.4142135623730951",
    "2.0",
    "3.0",
    "5.0",
]

def generate():
    lines = []
    lines.append("! SPDX-License-Identifier: MIT")
    lines.append("! Copyright (c) 2026 SciFort contributors")
    lines.append("!")
    lines.append("! Reference values for lognormal, weibull, pareto, and rayleigh distributions.")
    lines.append(f"! Generated with mpmath {mp.__version__} at {mp.mp.dps} digits.")
    lines.append("")
    lines.append("module test_continuous_reference")
    lines.append("    use scifort_kinds, only : dp")
    lines.append("    implicit none")
    lines.append("    private")
    lines.append("")


    # LOGNORMAL
    n_ln = len(LOGNORMAL_CASES)
    lines.append(f"    integer, parameter :: n_ln = {n_ln}")
    lines.append(f"    real(dp), dimension(n_ln), parameter, public :: ln_s = [&")
    for i, (s_str, x_str) in enumerate(LOGNORMAL_CASES):
        comma = "," if i < n_ln - 1 else ""
        lines.append(f"        {fstr(round_bin64(s_str))}{comma} &")
    lines.append("    ]")

    lines.append(f"    real(dp), dimension(n_ln), parameter, public :: ln_x = [&")
    for i, (s_str, x_str) in enumerate(LOGNORMAL_CASES):
        comma = "," if i < n_ln - 1 else ""
        lines.append(f"        {fstr(round_bin64(x_str))}{comma} &")
    lines.append("    ]")

    ln_pdf_vals, ln_cdf_vals, ln_sf_vals, ln_logcdf_vals, ln_logsf_vals, ln_tol_vals = [], [], [], [], [], []
    for s_str, x_str in LOGNORMAL_CASES:
        s = round_bin64(s_str)
        x = round_bin64(x_str)
        w = mp.log(x) / s
        pdf = (1 / (s * x * mp.sqrt(2 * mp.pi))) * mp.exp(-0.5 * w * w)
        cdf = 0.5 * mp.erfc(-w / mp.sqrt(2))
        sf = 0.5 * mp.erfc(w / mp.sqrt(2))
        logcdf = mp.log(cdf)
        logsf = mp.log(sf)
        tol = EPS_FACTOR * (1 + abs(w))
        ln_pdf_vals.append(pdf)
        ln_cdf_vals.append(cdf)
        ln_sf_vals.append(sf)
        ln_logcdf_vals.append(logcdf)
        ln_logsf_vals.append(logsf)
        ln_tol_vals.append(tol)

    for name, arr in [("ln_pdf", ln_pdf_vals), ("ln_cdf", ln_cdf_vals), ("ln_sf", ln_sf_vals),
                      ("ln_logcdf", ln_logcdf_vals), ("ln_logsf", ln_logsf_vals), ("ln_tol", ln_tol_vals)]:
        lines.append(f"    real(dp), dimension(n_ln), parameter, public :: {name} = [&")
        for i, val in enumerate(arr):
            comma = "," if i < n_ln - 1 else ""
            lines.append(f"        {fstr(val)}{comma} &")
        lines.append("    ]")
    lines.append("")

    # WEIBULL
    n_wb = len(WEIBULL_CASES)
    lines.append(f"    integer, parameter :: n_wb = {n_wb}")
    lines.append(f"    real(dp), dimension(n_wb), parameter, public :: wb_c = [&")
    for i, (c_str, x_str) in enumerate(WEIBULL_CASES):
        comma = "," if i < n_wb - 1 else ""
        lines.append(f"        {fstr(round_bin64(c_str))}{comma} &")
    lines.append("    ]")

    lines.append(f"    real(dp), dimension(n_wb), parameter, public :: wb_x = [&")
    for i, (c_str, x_str) in enumerate(WEIBULL_CASES):
        comma = "," if i < n_wb - 1 else ""
        lines.append(f"        {fstr(round_bin64(x_str))}{comma} &")
    lines.append("    ]")

    wb_pdf_vals, wb_cdf_vals, wb_sf_vals, wb_logcdf_vals, wb_logsf_vals, wb_tol_vals = [], [], [], [], [], []
    for c_str, x_str in WEIBULL_CASES:
        c = round_bin64(c_str)
        x = round_bin64(x_str)
        t = x ** c
        pdf = c * (x ** (c - 1)) * mp.exp(-t)
        sf = mp.exp(-t)
        cdf = -mp.expm1(-t)
        logsf = -t
        logcdf = mp.log(cdf)
        tol = EPS_FACTOR * (1 + t)
        wb_pdf_vals.append(pdf)
        wb_cdf_vals.append(cdf)
        wb_sf_vals.append(sf)
        wb_logcdf_vals.append(logcdf)
        wb_logsf_vals.append(logsf)
        wb_tol_vals.append(tol)

    for name, arr in [("wb_pdf", wb_pdf_vals), ("wb_cdf", wb_cdf_vals), ("wb_sf", wb_sf_vals),
                      ("wb_logcdf", wb_logcdf_vals), ("wb_logsf", wb_logsf_vals), ("wb_tol", wb_tol_vals)]:
        lines.append(f"    real(dp), dimension(n_wb), parameter, public :: {name} = [&")
        for i, val in enumerate(arr):
            comma = "," if i < n_wb - 1 else ""
            lines.append(f"        {fstr(val)}{comma} &")
        lines.append("    ]")
    lines.append("")

    # PARETO
    n_pa = len(PARETO_CASES)
    lines.append(f"    integer, parameter :: n_pa = {n_pa}")
    lines.append(f"    real(dp), dimension(n_pa), parameter, public :: pa_b = [&")
    for i, (b_str, x_str) in enumerate(PARETO_CASES):
        comma = "," if i < n_pa - 1 else ""
        lines.append(f"        {fstr(round_bin64(b_str))}{comma} &")
    lines.append("    ]")

    lines.append(f"    real(dp), dimension(n_pa), parameter, public :: pa_x = [&")
    for i, (b_str, x_str) in enumerate(PARETO_CASES):
        comma = "," if i < n_pa - 1 else ""
        lines.append(f"        {fstr(round_bin64(x_str))}{comma} &")
    lines.append("    ]")

    pa_pdf_vals, pa_cdf_vals, pa_sf_vals, pa_logcdf_vals, pa_logsf_vals, pa_tol_vals = [], [], [], [], [], []
    for b_str, x_str in PARETO_CASES:
        b = round_bin64(b_str)
        x = round_bin64(x_str)
        pdf = b * (x ** (-b - 1))
        sf = x ** (-b)
        cdf = 1 - sf
        logsf = -b * mp.log(x)
        logcdf = mp.log(cdf)
        tol = EPS_FACTOR * (1 + b * mp.log(x))
        pa_pdf_vals.append(pdf)
        pa_cdf_vals.append(cdf)
        pa_sf_vals.append(sf)
        pa_logcdf_vals.append(logcdf)
        pa_logsf_vals.append(logsf)
        pa_tol_vals.append(tol)

    for name, arr in [("pa_pdf", pa_pdf_vals), ("pa_cdf", pa_cdf_vals), ("pa_sf", pa_sf_vals),
                      ("pa_logcdf", pa_logcdf_vals), ("pa_logsf", pa_logsf_vals), ("pa_tol", pa_tol_vals)]:
        lines.append(f"    real(dp), dimension(n_pa), parameter, public :: {name} = [&")
        for i, val in enumerate(arr):
            comma = "," if i < n_pa - 1 else ""
            lines.append(f"        {fstr(val)}{comma} &")
        lines.append("    ]")
    lines.append("")

    # RAYLEIGH
    n_ry = len(RAYLEIGH_CASES)
    lines.append(f"    integer, parameter :: n_ry = {n_ry}")
    lines.append(f"    real(dp), dimension(n_ry), parameter, public :: ry_x = [&")
    for i, x_str in enumerate(RAYLEIGH_CASES):
        comma = "," if i < n_ry - 1 else ""
        lines.append(f"        {fstr(round_bin64(x_str))}{comma} &")
    lines.append("    ]")

    ry_pdf_vals, ry_cdf_vals, ry_sf_vals, ry_logcdf_vals, ry_logsf_vals, ry_tol_vals = [], [], [], [], [], []
    for x_str in RAYLEIGH_CASES:
        x = round_bin64(x_str)
        t = 0.5 * x * x
        pdf = x * mp.exp(-t)
        sf = mp.exp(-t)
        cdf = -mp.expm1(-t)
        logsf = -t
        logcdf = mp.log(cdf)
        tol = EPS_FACTOR * (1 + t)
        ry_pdf_vals.append(pdf)
        ry_cdf_vals.append(cdf)
        ry_sf_vals.append(sf)
        ry_logcdf_vals.append(logcdf)
        ry_logsf_vals.append(logsf)
        ry_tol_vals.append(tol)

    for name, arr in [("ry_pdf", ry_pdf_vals), ("ry_cdf", ry_cdf_vals), ("ry_sf", ry_sf_vals),
                      ("ry_logcdf", ry_logcdf_vals), ("ry_logsf", ry_logsf_vals), ("ry_tol", ry_tol_vals)]:
        lines.append(f"    real(dp), dimension(n_ry), parameter, public :: {name} = [&")
        for i, val in enumerate(arr):
            comma = "," if i < n_ry - 1 else ""
            lines.append(f"        {fstr(val)}{comma} &")
        lines.append("    ]")

    lines.append("end module test_continuous_reference")
    return "\n".join(lines) + "\n"

if __name__ == "__main__":
    sys.stdout.write(generate())
