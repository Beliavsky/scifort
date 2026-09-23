# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Generate references for geometric and negative-binomial tests.

Usage:
    python tools/generate_discrete_extended_reference.py > \
        test/discrete_extended_reference.f90

The formulas are evaluated with mpmath after decimal inputs have first been
rounded to binary64. Shape values n are chosen to be exactly representable in
binary64; tolerances account for binary64 rounding of p and for conditioning
of the requested probability.
"""

import sys

import mpmath as mp

mp.mp.dps = 70
EPS = mp.mpf(2) ** -52
EPS_FACTOR = 40 * EPS

GEOMETRIC_CASES = [
    ("1", "0.25"),
    ("4", "0.25"),
    ("30", "0.25"),
    ("1", "1e-8"),
    ("100000", "1e-5"),
    ("100", "0.999"),
]

NBINOM_CASES = [
    ("0", "0.5", "0.25"),
    ("3", "0.5", "0.25"),
    ("20", "0.5", "0.25"),
    ("0", "2.5", "0.3"),
    ("5", "2.5", "0.3"),
    ("50", "2.5", "0.3"),
    ("100", "50", "0.5"),
    ("1000", "1000", "0.5"),
    ("100000", "2.5", "1e-5"),
]

GEOMETRIC_QUANTILES = [
    ("0.01", "0.25"),
    ("0.5", "0.25"),
    ("0.99", "0.25"),
    ("1e-12", "1e-5"),
    ("0.999999", "1e-5"),
]

NBINOM_QUANTILES = [
    ("0.01", "0.5", "0.25"),
    ("0.5", "0.5", "0.25"),
    ("0.99", "2.5", "0.3"),
    ("1e-10", "50", "0.5"),
    ("0.999", "1000", "0.5"),
]


def mpf(text):
    return mp.mpf(float(text))


def geom_pmf(k, p):
    return p * (1 - p) ** (k - 1)


def geom_cdf(k, p):
    return 1 - (1 - p) ** k


def geom_sf(k, p):
    return (1 - p) ** k


def nbinom_pmf(k, n, p):
    coeff = mp.gamma(k + n) / (mp.gamma(n) * mp.gamma(k + 1))
    return coeff * p**n * (1 - p) ** k


def nbinom_cdf(k, n, p):
    return mp.betainc(n, k + 1, 0, p, regularized=True)


def nbinom_sf(k, n, p):
    return mp.betainc(k + 1, n, 0, 1 - p, regularized=True)


def tolerance(value, derivative_term):
    if value == 0:
        return mp.mpf(0)
    sensitivity = abs(derivative_term) / value
    return EPS_FACTOR * (1 + abs(mp.log(value)) + sensitivity) * value


def fmt(value):
    value = mp.mpf(value)
    if value == 0:
        return "0.0_dp"
    text = mp.nstr(value, 20, min_fixed=1, max_fixed=0)
    mantissa, _, exponent = text.partition("e")
    if "." not in mantissa:
        mantissa += ".0"
    exponent = int(exponent) if exponent else 0
    return f"{mantissa}e{exponent}_dp"


def emit_array(name, values):
    lines = [f"    real(dp), parameter, public :: {name}({len(values)}) = [ &"]
    for i, value in enumerate(values):
        sep = ", &" if i + 1 < len(values) else "]"
        lines.append(f"        {fmt(value)}{sep}")
    return "\n".join(lines)


def table(columns, rows):
    cols = {name: [] for name in columns}
    for row in rows:
        for name, value in zip(columns, row):
            cols[name].append(value)
    return [emit_array(name, values) for name, values in cols.items()]


def integer_quantile(cdf, target):
    lo = -1
    hi = 1
    while cdf(hi) < target:
        lo = hi
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if cdf(mid) >= target:
            hi = mid
        else:
            lo = mid
    return mp.mpf(hi)


def integer_isf(sf, target):
    lo = -1
    hi = 1
    while sf(hi) > target:
        lo = hi
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if sf(mid) <= target:
            hi = mid
        else:
            lo = mid
    return mp.mpf(hi)


def main(out):
    arrays = []

    rows = []
    for k_text, p_text in GEOMETRIC_CASES:
        k, p = mpf(k_text), mpf(p_text)
        pmf = geom_pmf(k, p)
        cdf = geom_cdf(k, p)
        sf = geom_sf(k, p)
        q = 1 - p
        d_pmf = pmf * (1 / p - (k - 1) / q) if q != 0 else mp.mpf(0)
        d_tail = k * q ** (k - 1) if q != 0 else mp.mpf(0)
        rows.append([
            k, p, pmf, cdf, sf,
            tolerance(pmf, p * d_pmf),
            tolerance(cdf, p * d_tail),
            tolerance(sf, p * d_tail),
        ])
    arrays += table([
        "ge_k", "ge_p", "ge_pmf", "ge_cdf", "ge_sf",
        "ge_tol_pmf", "ge_tol_cdf", "ge_tol_sf",
    ], rows)

    rows = []
    for k_text, n_text, p_text in NBINOM_CASES:
        k, n, p = mpf(k_text), mpf(n_text), mpf(p_text)
        pmf = nbinom_pmf(k, n, p)
        cdf = nbinom_cdf(k, n, p)
        sf = nbinom_sf(k, n, p)
        beta_density = p ** (n - 1) * (1 - p) ** k / mp.beta(n, k + 1)
        d_pmf = pmf * (n / p - k / (1 - p))
        rows.append([
            k, n, p, pmf, cdf, sf,
            tolerance(pmf, p * d_pmf),
            tolerance(cdf, p * beta_density),
            tolerance(sf, p * beta_density),
        ])
    arrays += table([
        "nb_k", "nb_n", "nb_p", "nb_pmf", "nb_cdf", "nb_sf",
        "nb_tol_pmf", "nb_tol_cdf", "nb_tol_sf",
    ], rows)

    rows = []
    for probability_text, p_text in GEOMETRIC_QUANTILES:
        probability, p = mpf(probability_text), mpf(p_text)
        ppf = integer_quantile(lambda k: geom_cdf(k, p), probability)
        isf = integer_isf(lambda k: geom_sf(k, p), probability)
        rows.append([probability, p, ppf, isf])
    arrays += table(["geq_prob", "geq_p", "geq_ppf", "geq_isf"], rows)

    rows = []
    for probability_text, n_text, p_text in NBINOM_QUANTILES:
        probability, n, p = mpf(probability_text), mpf(n_text), mpf(p_text)
        ppf = integer_quantile(lambda k: nbinom_cdf(k, n, p), probability)
        isf = integer_isf(lambda k: nbinom_sf(k, n, p), probability)
        rows.append([probability, n, p, ppf, isf])
    arrays += table(["nbq_prob", "nbq_n", "nbq_p", "nbq_ppf", "nbq_isf"], rows)

    out.write("! SPDX-License-Identifier: MIT\n")
    out.write("! Copyright (c) 2026 SciFort contributors\n\n")
    out.write("! Generated by tools/generate_discrete_extended_reference.py.\n")
    out.write(f"! mpmath {mp.__version__}, {mp.mp.dps} decimal digits.\n\n")
    out.write("module test_discrete_extended_reference\n")
    out.write("    use scifort_kinds, only : dp\n")
    out.write("    implicit none\n")
    out.write("    private\n\n")
    out.write("\n\n".join(arrays))
    out.write("\n\nend module test_discrete_extended_reference\n")


if __name__ == "__main__":
    main(sys.stdout)
