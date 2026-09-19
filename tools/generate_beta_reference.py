# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Generate test/beta_reference.f90 from high-precision mpmath values.

Usage:
    python tools/generate_beta_reference.py > test/beta_reference.f90

Requires mpmath. Inputs are rounded to binary64 before the references are
computed, so the tests compare against the exact function of the arguments
they pass. Tolerances follow the conditioning of each problem:

- a tail V of I_x(a, b) with density f gets relative tolerance
  EPS_FACTOR * (1 + |log V| + z f / V), z = min(x, 1 - x), because the
  implementation consumes x when x <= 1/2 and the exact 1 - x otherwise;
  the absolute tolerance is min(r(P) P, r(Q) Q), floored at one unit of V,
  because either tail may be formed as the complement of the other;
- log(V) gets the absolute tolerance of V divided by V;
- an inverse gets an absolute tolerance EPS_FACTOR * z * (1 + (1 + |log v|)
  v / (z f)) with z = min(x, 1 - x) and v = min(p, 1 - p);
- a log density L gets absolute tolerance EPS_FACTOR * (1 + |L| + s), where
  s is the relative sensitivity of the density to its argument.

Each tail is computed independently (never as a complement) at increasing
precision until two precisions agree to 25 digits.
"""

import sys

import mpmath as mp

mp.mp.dps = 50
EPS = mp.mpf(2) ** -52
EPS_FACTOR = 20 * EPS

# (a, b, x)
BETAINC_CASES = [
    ("0.5", "0.5", "0.3"),
    ("1", "1", "0.25"),
    ("2", "3", "0.4"),
    ("0.01", "1", "0.3"),
    ("1e-5", "1", "0.999"),
    ("6.77", "1.08e-5", "0.9999984"),
    ("3.9e-6", "3.4", "1.16e-6"),
    ("10", "10", "0.5"),
    ("50", "7", "0.9"),
    ("100", "0.5", "0.99"),
    ("1000", "1000", "0.49"),
    ("10000", "2", "0.9995"),
    ("100000", "100000", "0.501"),
    ("5", "400", "0.02"),
    ("0.5", "20000", "1e-4"),
    ("150000", "0.5", "0.999982"),
    ("30", "30", "0.05"),
    ("2", "5", "1e-100"),
    ("18.6", "2324", "0.0087"),
]

# (a, b, p)
INVERSE_CASES = [
    ("0.5", "0.5", "0.3"),
    ("2", "3", "1e-10"),
    ("2", "3", "0.999999"),
    ("0.1", "0.2", "0.5"),
    ("50", "7", "0.05"),
    ("1000", "0.5", "0.975"),
    ("0.05", "5", "0.9"),
    ("300", "400", "1e-200"),
]

# (x, df)
T_CASES = [
    ("0.5", "1"),
    ("-3", "2.5"),
    ("2.2", "30"),
    ("-40", "4"),
    ("1e-8", "10"),
    ("1e30", "3"),
    ("-2.5", "1e6"),
    ("4", "0.3"),
]

# (x, dfn, dfd)
F_CASES = [
    ("1", "1", "1"),
    ("3.5", "5", "12"),
    ("0.01", "3", "40"),
    ("250", "2", "7"),
    ("1.1", "1000", "900"),
    ("0.9", "4", "100000"),
    ("7200", "800000", "0.25"),
]

# (x, a, b, loc, scale)
BETA_DIST_CASES = [
    ("0.3", "2", "5", "0", "1"),
    ("0.999", "0.5", "0.5", "0", "1"),
    ("3.5", "1.5", "4", "2", "3"),
    ("0.145", "1186", "7778", "0", "1"),
    ("1e-6", "0.2", "3", "0", "1"),
]


def mpf(text):
    return mp.mpf(float(text))


def continued_fraction_tail(a, b, x):
    """I_x(a, b) from DLMF 8.17.22 evaluated in the current precision.

    Used when mpmath's hypergeometric evaluation does not converge. The
    iteration stops only after 50 consecutive steps with |delta - 1| below
    the working tolerance, so a transient delta = 1 cannot end it early.
    """
    tiny = mp.mpf(10) ** (-3 * mp.mp.dps)
    tolerance = mp.mpf(10) ** (-mp.mp.dps + 10)
    f = mp.mpf(1)
    c = f
    d = mp.mpf(0)
    run = 0
    n = 0
    while run < 50:
        n += 1
        if n % 2 == 0:
            m = n // 2
            num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            m = (n - 1) // 2
            num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1 + num * d
        d = d if d != 0 else tiny
        c = 1 + num / c
        c = c if c != 0 else tiny
        d = 1 / d
        delta = c * d
        f *= delta
        run = run + 1 if abs(delta - 1) < tolerance else 0
        if n > 10 ** 7:
            raise RuntimeError("continued fraction did not converge")
    log_kernel = a * mp.log(x) + b * mp.log1p(-x) - mp.log(mp.beta(a, b))
    return mp.exp(log_kernel) / (a * f)


def lower_raw(a, b, x):
    try:
        return mp.betainc(a, b, 0, x, regularized=True)
    except mp.libmp.libhyper.NoConvergence:
        return continued_fraction_tail(a, b, x)


def tail(a, b, x):
    """I_x(a, b) computed directly, with precision escalation."""
    previous = None
    for dps in (50, 100, 200, 400):
        with mp.workdps(dps):
            value = lower_raw(mp.mpf(a), mp.mpf(b), mp.mpf(x))
        if previous is not None and value > 0 and abs(value / previous - 1) < mp.mpf(10) ** -25:
            return +value
        previous = value
    raise RuntimeError(f"no reference for I_{x}({a}, {b})")


def tails(a, b, x):
    """P = I_x(a, b) and Q = I_(1-x)(b, a), each computed directly."""
    y = 1 - x
    p = tail(a, b, x)
    q = tail(b, a, y)
    if abs(p + q - 1) > mp.mpf(10) ** -30:
        raise RuntimeError("P + Q != 1")
    return p, q


def density(a, b, x):
    return mp.exp((a - 1) * mp.log(x) + (b - 1) * mp.log1p(-x) - mp.log(mp.beta(a, b)))


def tolerances(p, q, z_f):
    rp = EPS_FACTOR * (1 + abs(mp.log(p)) + z_f / p)
    rq = EPS_FACTOR * (1 + abs(mp.log(q)) + z_f / q)
    abs_tol = min(rp * p, rq * q)
    tol_p = max(abs_tol, EPS * p)
    tol_q = max(abs_tol, EPS * q)
    return tol_p, tol_q, tol_p / p, tol_q / q


def invert(a, b, target, complement):
    """Root x of I_x(a, b) = target, complement = 1 - target exactly.

    The root is found in log(x) when it is at most 1/2 and in log(1 - x)
    otherwise, with the residual taken in the tail that is represented
    exactly.
    """
    if target <= complement:
        lower_side = tail(a, b, mp.mpf("0.5")) >= target
    else:
        lower_side = tail(b, a, mp.mpf("0.5")) <= complement
    if lower_side:
        def h(t):
            return mp.log(tail(a, b, mp.exp(t))) - mp.log(target)
    else:
        def h(t):
            return mp.log(complement) - mp.log(tail(b, a, mp.exp(t)))
    lo, hi = mp.mpf(-745), mp.log(mp.mpf("0.5"))
    if not lower_side:
        # h increases with x, so decreases with t = log(1 - x).
        def g(t):
            return -h(t)
    else:
        g = h
    if not (g(lo) < 0 < g(hi)):
        raise RuntimeError("root not bracketed")
    for _ in range(60):
        mid = (lo + hi) / 2
        if g(mid) < 0:
            lo = mid
        else:
            hi = mid
    t = mp.findroot(g, (lo, hi), solver="secant", tol=mp.mpf(10) ** -60)
    z = mp.exp(t)
    return (z, 1 - z) if lower_side else (1 - z, z)


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
    for a_text, b_text, x_text in BETAINC_CASES:
        a, b, x = mpf(a_text), mpf(b_text), mpf(x_text)
        with mp.workdps(60):
            p, q = tails(a, b, x)
            z_f = min(x, 1 - x) * density(a, b, x)
            tol_p, tol_q, tol_logp, tol_logq = tolerances(p, q, z_f)
            rows.append([a, b, x, p, q, mp.log(p), mp.log(q), tol_p, tol_q,
                         tol_logp, tol_logq])
    arrays += table(["ib_a", "ib_b", "ib_x", "ib_p", "ib_q", "ib_logp", "ib_logq",
                     "ib_tol_p", "ib_tol_q", "ib_tol_logp", "ib_tol_logq"], rows)

    rows = []
    for a_text, b_text, p_text in INVERSE_CASES:
        a, b, p = mpf(a_text), mpf(b_text), mpf(p_text)
        with mp.workdps(60):
            x, y = invert(a, b, p, 1 - p)
            z = min(x, y)
            v = min(p, 1 - p)
            z_f = z * density(a, b, x)
            tol = EPS_FACTOR * z * (1 + (1 + abs(mp.log(v))) * v / z_f)
            # Complement inverse: 1 - I_x = p has the root x' with I_x' = 1 - p.
            xc, yc = invert(a, b, 1 - p, p)
            zc = min(xc, yc)
            z_fc = zc * density(a, b, xc)
            tolc = EPS_FACTOR * zc * (1 + (1 + abs(mp.log(v))) * v / z_fc)
            rows.append([a, b, p, x, tol, xc, tolc])
    arrays += table(["ibinv_a", "ibinv_b", "ibinv_p", "ibinv_x", "ibinv_tol",
                     "ibcinv_x", "ibcinv_tol"], rows)

    rows = []
    for x_text, df_text in T_CASES:
        t, df = mpf(x_text), mpf(df_text)
        with mp.workdps(60):
            a = df / 2
            xb = df / (df + t * t)
            far = tail(a, mp.mpf("0.5"), xb) / 2
            near = 1 - far
            cdf, sf = (far, near) if t < 0 else (near, far)
            logpdf = (mp.loggamma((df + 1) / 2) - mp.loggamma(df / 2) - mp.log(df * mp.pi) / 2
                      - (df + 1) / 2 * mp.log1p(t * t / df))
            z_f = abs(t) * mp.exp(logpdf)
            tol_cdf, tol_sf, _, _ = tolerances(cdf, sf, z_f)
            sensitivity = (df + 1) * t * t / (df + t * t)
            rows.append([t, df, logpdf, EPS_FACTOR * (1 + abs(logpdf) + sensitivity),
                         cdf, sf, tol_cdf, tol_sf])
    arrays += table(["st_x", "st_df", "st_logpdf", "st_tol_logpdf", "st_cdf", "st_sf",
                     "st_tol_cdf", "st_tol_sf"], rows)

    rows = []
    for x_text, dfn_text, dfd_text in F_CASES:
        z, dfn, dfd = mpf(x_text), mpf(dfn_text), mpf(dfd_text)
        with mp.workdps(60):
            a, b = dfn / 2, dfd / 2
            xb = dfn * z / (dfn * z + dfd)
            yb = dfd / (dfn * z + dfd)
            p = tail(a, b, xb)
            q = tail(b, a, yb)
            logpdf = a * mp.log(xb) + b * mp.log(yb) - mp.log(mp.beta(a, b)) - mp.log(z)
            z_f = z * mp.exp(logpdf)
            tol_p, tol_q, _, _ = tolerances(p, q, z_f)
            sensitivity = abs((a - 1) - (a + b) * xb) + 1
            rows.append([z, dfn, dfd, logpdf, EPS_FACTOR * (1 + abs(logpdf) + sensitivity),
                         p, q, tol_p, tol_q])
    arrays += table(["fd_x", "fd_dfn", "fd_dfd", "fd_logpdf", "fd_tol_logpdf", "fd_cdf",
                     "fd_sf", "fd_tol_cdf", "fd_tol_sf"], rows)

    rows = []
    for x_text, a_text, b_text, loc_text, scale_text in BETA_DIST_CASES:
        xv, a, b, loc, scale = map(mpf, (x_text, a_text, b_text, loc_text, scale_text))
        with mp.workdps(60):
            z = (xv - loc) / scale
            p, q = tails(a, b, z)
            logpdf = mp.log(density(a, b, z)) - mp.log(scale)
            z_f = min(z, 1 - z) * density(a, b, z)
            tol_p, tol_q, _, _ = tolerances(p, q, z_f)
            sensitivity = abs(a - 1) + abs(b - 1) * z / (1 - z)
            rows.append([xv, a, b, loc, scale, logpdf,
                         EPS_FACTOR * (1 + abs(logpdf) + sensitivity), p, q, tol_p, tol_q])
    arrays += table(["bd_x", "bd_a", "bd_b", "bd_loc", "bd_scale", "bd_logpdf",
                     "bd_tol_logpdf", "bd_cdf", "bd_sf", "bd_tol_cdf", "bd_tol_sf"], rows)

    out.write("! SPDX-License-Identifier: MIT\n")
    out.write("! Copyright (c) 2026 SciFort contributors\n\n")
    out.write("! Generated by tools/generate_beta_reference.py. Do not edit by hand.\n")
    out.write(f"! mpmath {mp.__version__}, at least {mp.mp.dps} decimal digits.\n\n")
    out.write("module test_beta_reference\n")
    out.write("    use scifort_kinds, only : dp\n")
    out.write("    implicit none\n")
    out.write("    private\n\n")
    out.write("\n\n".join(arrays))
    out.write("\n\nend module test_beta_reference\n")


if __name__ == "__main__":
    main(sys.stdout)
