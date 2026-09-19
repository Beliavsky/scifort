# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Generate test/discrete_reference.f90 from high-precision mpmath values.

Usage:
    python tools/generate_discrete_reference.py > test/discrete_reference.f90

Requires mpmath. Inputs are rounded to binary64 before the references are
computed. The counts and the number of trials are exact integers, so the only
inexact inputs are mu and p; tolerances therefore combine

- the relative error of the underlying special function,
  EPS_FACTOR * (1 + |log V|) for a tail V, and
- the sensitivity of V to the rounding of the continuous parameter,
  EPS_FACTOR * |parameter dV/dparameter| / V.

For the Poisson distribution d/dmu of P(X <= k) is -pmf(k); for the binomial
distribution d/dp of P(X <= k) is -n pmf(k; n - 1, p). A mass function value
m gets the relative tolerance EPS_FACTOR * (1 + |log m| + sensitivity).
"""

import os
import sys

import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_beta_reference import continued_fraction_tail  # noqa: E402

mp.mp.dps = 50
EPS = mp.mpf(2) ** -52
EPS_FACTOR = 20 * EPS

# (k, mu)
POISSON_CASES = [
    ("0", "2.5"),
    ("3", "2.5"),
    ("10", "2.5"),
    ("1", "1e-8"),
    ("0", "1e-8"),
    ("50", "1"),
    ("100", "120"),
    ("200", "100"),
    ("1000", "1000"),
    ("999000", "1000000"),
    ("7", "0.001"),
]

# (k, n, p)
BINOMIAL_CASES = [
    ("0", "10", "0.25"),
    ("3", "10", "0.25"),
    ("10", "10", "0.25"),
    ("2", "5", "0.5"),
    ("500", "1000", "0.5"),
    ("900", "1000", "0.5"),
    ("1", "1000000", "1e-6"),
    ("30", "40", "0.9"),
    ("0", "3", "1e-12"),
    ("100000", "1000000", "0.1"),
]

# (p, mu) for the Poisson quantiles
POISSON_QUANTILE_CASES = [
    ("0.5", "2.5"),
    ("0.01", "10"),
    ("0.99", "10"),
    ("1e-12", "100"),
    ("0.975", "1000"),
    ("0.3", "0.5"),
]

# (p, n, prob) for the binomial quantiles
BINOMIAL_QUANTILE_CASES = [
    ("0.5", "10", "0.25"),
    ("0.05", "100", "0.3"),
    ("0.95", "100", "0.3"),
    ("1e-9", "1000", "0.5"),
    ("0.999", "50", "0.02"),
]


def mpf(text):
    return mp.mpf(float(text))


def poisson_pmf(k, mu):
    if mu == 0:
        return mp.mpf(1) if k == 0 else mp.mpf(0)
    return mp.exp(k * mp.log(mu) - mu - mp.loggamma(k + 1))


def poisson_cdf(k, mu):
    """P(X <= k) = Q(k + 1, mu) (DLMF 8.4.8)."""
    if mu == 0:
        return mp.mpf(1)
    return mp.gammainc(k + 1, mu, mp.inf, regularized=True)


def poisson_sf(k, mu):
    if mu == 0:
        return mp.mpf(0)
    return mp.gammainc(k + 1, 0, mu, regularized=True)


def binomial_pmf(k, n, p):
    if p == 0:
        return mp.mpf(1) if k == 0 else mp.mpf(0)
    if p == 1:
        return mp.mpf(1) if k == n else mp.mpf(0)
    return mp.binomial(n, k) * p ** k * (1 - p) ** (n - k)


def regularized_beta(a, b, x):
    """I_x(a, b), falling back to the continued fraction of DLMF 8.17.22."""
    try:
        return mp.betainc(a, b, 0, x, regularized=True)
    except mp.libmp.libhyper.NoConvergence:
        return continued_fraction_tail(a, b, x)


def binomial_sf(k, n, p):
    """P(X > k) = I_p(k + 1, n - k) (DLMF 8.17.5)."""
    if k >= n or p == 0:
        return mp.mpf(0)
    if p == 1:
        return mp.mpf(1)
    return regularized_beta(k + 1, n - k, p)


def binomial_cdf(k, n, p):
    if k >= n or p == 0:
        return mp.mpf(1)
    if p == 1:
        return mp.mpf(0)
    return regularized_beta(n - k, k + 1, 1 - p)


def check_sum(cdf, sf, pmf_terms, label):
    """Confirm the closed forms against a direct sum for small cases."""
    if abs(cdf + sf - 1) > mp.mpf(10) ** -40:
        raise RuntimeError(f"tails do not sum to one for {label}")
    if pmf_terms is not None and abs(sum(pmf_terms) - cdf) > mp.mpf(10) ** -35:
        raise RuntimeError(f"cdf does not match the direct sum for {label}")


def tolerance(value, derivative_term):
    """Absolute tolerance for value, given |parameter dV/dparameter|.

    An exact zero (an endpoint of the support) is required exactly.
    """
    if value == 0:
        return mp.mpf(0)
    sensitivity = derivative_term / value
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
    return [emit_array(k, v) for k, v in cols.items()]


def main(out):
    arrays = []

    rows = []
    for k_text, mu_text in POISSON_CASES:
        k, mu = mpf(k_text), mpf(mu_text)
        with mp.workdps(60):
            pmf = poisson_pmf(k, mu)
            cdf = poisson_cdf(k, mu)
            sf = poisson_sf(k, mu)
            terms = [poisson_pmf(mp.mpf(j), mu) for j in range(int(k) + 1)] \
                if k <= 300 else None
            check_sum(cdf, sf, terms, f"poisson({k}, {mu})")
            # d/dmu of the mass function is pmf(k - 1) - pmf(k).
            d_pmf = (poisson_pmf(k - 1, mu) if k >= 1 else mp.mpf(0)) - pmf
            rows.append([k, mu, pmf, cdf, sf,
                         tolerance(pmf, mu * abs(d_pmf)),
                         tolerance(cdf, mu * pmf),
                         tolerance(sf, mu * pmf)])
    arrays += table(["po_k", "po_mu", "po_pmf", "po_cdf", "po_sf", "po_tol_pmf",
                     "po_tol_cdf", "po_tol_sf"], rows)

    rows = []
    for k_text, n_text, p_text in BINOMIAL_CASES:
        k, n, p = mpf(k_text), mpf(n_text), mpf(p_text)
        with mp.workdps(60):
            pmf = binomial_pmf(k, n, p)
            cdf = binomial_cdf(k, n, p)
            sf = binomial_sf(k, n, p)
            terms = [binomial_pmf(mp.mpf(j), n, p) for j in range(int(k) + 1)] \
                if k <= 300 else None
            check_sum(cdf, sf, terms, f"binomial({k}, {n}, {p})")
            # d/dp of the mass function is n (pmf(k - 1; n - 1) - pmf(k; n - 1)).
            lower = binomial_pmf(k - 1, n - 1, p) if k >= 1 else mp.mpf(0)
            upper = binomial_pmf(k, n - 1, p) if k <= n - 1 else mp.mpf(0)
            d_pmf = n * (lower - upper)
            d_cdf = n * binomial_pmf(k, n - 1, p) if k <= n - 1 else mp.mpf(0)
            rows.append([k, n, p, pmf, cdf, sf,
                         tolerance(pmf, p * abs(d_pmf)),
                         tolerance(cdf, p * d_cdf),
                         tolerance(sf, p * d_cdf)])
    arrays += table(["bi_k", "bi_n", "bi_p", "bi_pmf", "bi_cdf", "bi_sf",
                     "bi_tol_pmf", "bi_tol_cdf", "bi_tol_sf"], rows)

    rows = []
    for p_text, mu_text in POISSON_QUANTILE_CASES:
        p, mu = mpf(p_text), mpf(mu_text)
        with mp.workdps(60):
            k = mp.mpf(0)
            while poisson_cdf(k, mu) < p:
                k += 1
            upper = mp.mpf(0)
            while poisson_sf(upper, mu) > p:
                upper += 1
            rows.append([p, mu, k, upper])
    arrays += table(["poq_p", "poq_mu", "poq_ppf", "poq_isf"], rows)

    rows = []
    for p_text, n_text, prob_text in BINOMIAL_QUANTILE_CASES:
        p, n, prob = mpf(p_text), mpf(n_text), mpf(prob_text)
        with mp.workdps(60):
            k = mp.mpf(0)
            while k < n and binomial_cdf(k, n, prob) < p:
                k += 1
            upper = mp.mpf(0)
            while upper < n and binomial_sf(upper, n, prob) > p:
                upper += 1
            rows.append([p, n, prob, k, upper])
    arrays += table(["biq_p", "biq_n", "biq_prob", "biq_ppf", "biq_isf"], rows)

    out.write("! SPDX-License-Identifier: MIT\n")
    out.write("! Copyright (c) 2026 SciFort contributors\n\n")
    out.write("! Generated by tools/generate_discrete_reference.py. Do not edit by hand.\n")
    out.write(f"! mpmath {mp.__version__}, at least {mp.mp.dps} decimal digits.\n\n")
    out.write("module test_discrete_reference\n")
    out.write("    use scifort_kinds, only : dp\n")
    out.write("    implicit none\n")
    out.write("    private\n\n")
    out.write("\n\n".join(arrays))
    out.write("\n\nend module test_discrete_reference\n")


if __name__ == "__main__":
    main(sys.stdout)
