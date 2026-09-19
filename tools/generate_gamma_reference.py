# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Generate test/gamma_reference.f90 from high-precision mpmath values.

Usage:
    python tools/generate_gamma_reference.py > test/gamma_reference.f90

Requires mpmath. The generated file records the mpmath version and working
precision. Each reference value carries a tolerance derived from the
conditioning of the problem rather than one blanket tolerance:

- a probability V (P or Q) with x f(x) = x**a exp(-x) / Gamma(a) gets
  relative tolerance r(V) = EPS_FACTOR * (1 + |log V| + x f(x) / V), where
  |log V| accounts for rounding in the exponent of tiny values and
  x f(x) / V is the condition number with respect to x. Because one tail
  may be computed as the complement of the other, the absolute tolerance is
  min(r(P) P, r(Q) Q), floored at one unit of V;
- log(V) gets the absolute tolerance of V divided by V;
- an inverse x solving V(a, x) = v gets relative tolerance
  EPS_FACTOR * (1 + (1 + |log v|) v / (x f(x)));
- inputs are rounded to binary64 before the references are computed, so
  the tests compare against the exact function of the arguments they pass;
- a log-density value L at standardized z gets absolute tolerance
  EPS_FACTOR * (1 + |L| + |a - 1 - z|).
"""

import sys

import mpmath as mp

mp.mp.dps = 50
EPS = mp.mpf(2) ** -52
EPS_FACTOR = 20 * EPS

INCGAMMA_CASES = [
    ("1e-10", "1e-5"),
    ("1e-10", "1.4"),
    ("1e-10", "3"),
    ("0.01", "0.001"),
    ("0.1", "1e-300"),
    ("0.5", "0.5"),
    ("0.5", "1.5"),
    ("0.5", "10"),
    ("0.9", "1.4"),
    ("1", "1e-20"),
    ("1", "1"),
    ("1", "50"),
    ("2.5", "0.1"),
    ("3", "700"),
    ("5", "5"),
    ("10", "3"),
    ("10", "10"),
    ("10", "25"),
    ("64.5", "0.005"),
    ("100", "90"),
    ("100", "100"),
    ("100", "130"),
    ("1000", "900"),
    ("10000", "10000"),
    ("10000", "10500"),
    ("100000", "100300"),
    ("1000000", "1000000"),
]

INVERSE_CASES = [
    ("0.01", "1e-3"),
    ("0.1", "1e-10"),
    ("0.5", "0.25"),
    ("0.5", "0.999"),
    ("1", "0.5"),
    ("1", "1e-300"),
    ("2.5", "0.01"),
    ("3", "0.999999999999"),
    ("10", "0.5"),
    ("10", "1e-100"),
    ("100", "0.05"),
    ("1e4", "0.975"),
    ("1e5", "1e-20"),
]

# (x, a, loc, scale)
GAMMA_DENSITY_CASES = [
    ("2", "3", "0", "1"),
    ("0.5", "0.5", "0", "1"),
    ("1e-3", "0.5", "0", "1"),
    ("250", "200", "0", "1"),
    ("7", "2.5", "1", "2"),
    ("100100", "100000", "0", "1"),
    ("30", "1", "0", "1"),
]

# (x, df)
CHI2_CASES = [
    ("3.84", "1"),
    ("0.001", "1"),
    ("18.3", "10"),
    ("1", "2"),
    ("150", "100"),
    ("400", "3"),
]


def mpf(text):
    """Return the binary64 value that the Fortran tests receive for text."""
    return mp.mpf(float(text))


def lower(a, x):
    try:
        return mp.gammainc(a, 0, x, regularized=True)
    except mp.libmp.libhyper.NoConvergence:
        # DLMF 8.7.1 summed by mpmath with a larger term limit.
        return (mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1))
                * mp.hyp1f1(1, a + 1, x, maxterms=10**8))


def upper(a, x):
    try:
        return mp.gammainc(a, x, mp.inf, regularized=True)
    except mp.libmp.libhyper.NoConvergence:
        with mp.workdps(mp.mp.dps + 400):
            return +(1 - lower(a, x))


def pq(a, x):
    """Return P, Q, log P, log Q, each computed without complementing."""
    p = lower(a, x)
    q = upper(a, x)
    if abs(p + q - 1) > mp.mpf(10) ** -40:
        raise RuntimeError(f"P + Q != 1 for a={a}, x={x}")
    return p, q, mp.log(p), mp.log(q)


def x_density(a, x):
    return mp.exp(a * mp.log(x) - x - mp.loggamma(a))


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


def prob_tolerances(a, x, p, q):
    xf = x_density(a, x)
    rp = EPS_FACTOR * (1 + abs(mp.log(p)) + xf / p)
    rq = EPS_FACTOR * (1 + abs(mp.log(q)) + xf / q)
    abs_tol = min(rp * p, rq * q)
    tol_p = max(abs_tol, EPS * p)
    tol_q = max(abs_tol, EPS * q)
    return tol_p, tol_q, tol_p / p, tol_q / q


def invert(h):
    """Root in x = exp(t) of h(t), which is increasing in t.

    Bisection on t in [-745, 745] followed by secant refinement.
    """
    lo, hi = mp.mpf(-745), mp.mpf(745)
    if not (h(lo) < 0 < h(hi)):
        raise RuntimeError("root is not bracketed")
    for _ in range(64):
        mid = (lo + hi) / 2
        if h(mid) < 0:
            lo = mid
        else:
            hi = mid
    t = mp.findroot(h, (lo, hi), solver="secant", tol=mp.mpf(10) ** -70)
    if not (lo - (hi - lo) <= t <= hi + (hi - lo)):
        raise RuntimeError("secant refinement left the bracket")
    return mp.exp(t)


def standard_logpdf(a, z):
    return (a - 1) * mp.log(z) - z - mp.loggamma(a)


def main(out):
    arrays = []

    cols = {k: [] for k in ["ig_a", "ig_x", "ig_p", "ig_q", "ig_logp", "ig_logq",
                            "ig_tol_p", "ig_tol_q", "ig_tol_logp", "ig_tol_logq"]}
    for a_text, x_text in INCGAMMA_CASES:
        a, x = mpf(a_text), mpf(x_text)
        p, q, logp, logq = pq(a, x)
        tol_p, tol_q, tol_logp, tol_logq = prob_tolerances(a, x, p, q)
        for key, value in zip(cols, [a, x, p, q, logp, logq, tol_p, tol_q,
                                     tol_logp, tol_logq]):
            cols[key].append(value)
    arrays += [emit_array(k, v) for k, v in cols.items()]

    cols = {k: [] for k in ["inv_a", "inv_p", "inv_x", "inv_tol",
                            "cinv_x", "cinv_tol"]}
    for a_text, p_text in INVERSE_CASES:
        a, p = mpf(a_text), mpf(p_text)
        with mp.workdps(80):
            x_lower = invert(lambda t: mp.log(lower(a, mp.exp(t))) - mp.log(p))
            x_upper = invert(lambda t: mp.log(p) - mp.log(upper(a, mp.exp(t))))
        for key, x in [("inv", x_lower), ("cinv", x_upper)]:
            v = min(p, 1 - p)
            tol = EPS_FACTOR * (1 + (1 + abs(mp.log(v))) * v / x_density(a, x))
            cols[f"{key}_x"].append(x)
            cols[f"{key}_tol"].append(tol)
        cols["inv_a"].append(a)
        cols["inv_p"].append(p)
    arrays += [emit_array(k, v) for k, v in cols.items()]

    cols = {k: [] for k in ["gd_x", "gd_a", "gd_loc", "gd_scale", "gd_logpdf",
                            "gd_tol_logpdf", "gd_cdf", "gd_sf", "gd_tol_cdf",
                            "gd_tol_sf"]}
    for x_text, a_text, loc_text, scale_text in GAMMA_DENSITY_CASES:
        x, a, loc, scale = map(mpf, (x_text, a_text, loc_text, scale_text))
        z = (x - loc) / scale
        logpdf = standard_logpdf(a, z) - mp.log(scale)
        p, q, _, _ = pq(a, z)
        tol_p, tol_q, _, _ = prob_tolerances(a, z, p, q)
        values = [x, a, loc, scale, logpdf,
                  EPS_FACTOR * (1 + abs(logpdf) + abs(a - 1 - z)), p, q, tol_p, tol_q]
        for key, value in zip(cols, values):
            cols[key].append(value)
    arrays += [emit_array(k, v) for k, v in cols.items()]

    cols = {k: [] for k in ["c2_x", "c2_df", "c2_logpdf", "c2_tol_logpdf",
                            "c2_cdf", "c2_sf", "c2_tol_cdf", "c2_tol_sf"]}
    for x_text, df_text in CHI2_CASES:
        x, df = mpf(x_text), mpf(df_text)
        a, z = df / 2, x / 2
        logpdf = standard_logpdf(a, z) - mp.log(2)
        p, q, _, _ = pq(a, z)
        tol_p, tol_q, _, _ = prob_tolerances(a, z, p, q)
        values = [x, df, logpdf, EPS_FACTOR * (1 + abs(logpdf) + abs(a - 1 - z)),
                  p, q, tol_p, tol_q]
        for key, value in zip(cols, values):
            cols[key].append(value)
    arrays += [emit_array(k, v) for k, v in cols.items()]

    out.write("! SPDX-License-Identifier: MIT\n")
    out.write("! Copyright (c) 2026 SciFort contributors\n\n")
    out.write("! Generated by tools/generate_gamma_reference.py. Do not edit by hand.\n")
    out.write(f"! mpmath {mp.__version__}, {mp.mp.dps} decimal digits.\n\n")
    out.write("module test_gamma_reference\n")
    out.write("    use scifort_kinds, only : dp\n")
    out.write("    implicit none\n")
    out.write("    private\n\n")
    out.write("\n\n".join(arrays))
    out.write("\n\nend module test_gamma_reference\n")


if __name__ == "__main__":
    main(sys.stdout)
