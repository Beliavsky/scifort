# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Print SciPy reference values for Milestone 6 special-function tests."""

import numpy as np
import scipy
import mpmath as mp
from scipy import special


def main():
    print(f"numpy {np.__version__}")
    print(f"scipy {scipy.__version__}")

    for x in (-5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0):
        print("erf", x, special.erf(x), special.erfc(x))
    for y in (-0.999999999999, -0.9, -0.1, 0.0, 0.1, 0.9, 0.999999999999):
        print("erfinv", y, special.erfinv(y))
    for y in (1e-200, 1e-12, 0.1, 1.0, 1.9, 2.0 - 1e-12):
        print("erfcinv", y, special.erfcinv(y))

    for x in (-10.0, -2.0, -0.5, 0.0, 0.5, 2.0, 6.0, 10.0, 20.0):
        print("dawson_error", x, special.dawsn(x), special.erfi(x), special.erfcx(x))
    for x in (-6.0, -4.0, -2.0, -0.5, 0.0, 0.5, 2.0, 4.0, 6.0, 20.0):
        print("fresnel", x, *special.fresnel(x))
    for x in (0.0, 0.5, 1.0, 3.0, 10.0, 1000.0):
        print("modfresnel_scipy", x, *special.modfresnelp(x), *special.modfresnelm(x))
    mp.mp.dps = 80
    for x in (-10.0, -3.0, -1.0, -0.5, 0.0, 0.5, 1.0, 3.0, 10.0,
              1000.0, 100000000.25):
        numerator, denominator = float(x).as_integer_ratio()
        xx = mp.mpf(numerator) / denominator
        fp = (mp.sqrt(mp.pi) / 2 * mp.exp(1j * mp.pi / 4) *
              mp.erfc(mp.exp(-1j * mp.pi / 4) * xx))
        kp = fp / mp.sqrt(mp.pi) * mp.exp(-1j * (xx * xx + mp.pi / 4))
        print("modfresnel_mpmath", x, fp, kp, mp.conj(fp), mp.conj(kp))
    for z in (0.0 + 0.0j, 0.25 + 0.5j, -1.25 + 0.75j,
              2.0 + 0.1j, 3.0 + 2.0j, 8.0 + 3.0j,
              2.0 - 0.5j, -3.0 - 1.25j):
        print("faddeeva", z, special.wofz(z), special.erfcx(z))
    mp.mp.dps = 80
    for x in (1000000.125, 100000000.25, 134262676.7054991,
              4503599627370497.0, 9007199254740992.0):
        numerator, denominator = float(x).as_integer_ratio()
        xx = mp.mpf(numerator) / denominator
        print("fresnel_mpmath", x, mp.fresnels(xx), mp.fresnelc(xx))

    for x in (-4.75, -2.5, -0.5, 0.5, 1.0, 3.5, 20.0):
        print("gammaln", x, special.gammaln(x))
    for x in (-170.75, -170.5, -10.25, -4.75, -2.5, -0.5,
              0.1, 0.5, 1.0, 5.0, 20.0, 171.0, 172.0):
        print("gamma_family", x, special.gamma(x), special.rgamma(x),
              special.gammasgn(x), special.loggamma(x))
    for z in (1.5 + 2.0j, 2.5 + 1.0j, -0.5 + 0.0j,
              -3.2 + 0.0j, -3.2 + 0.7j, 0.25 - 4.0j,
              10.0 + 15.0j, -100.0 + 50.0j, 100.0 + 200.0j):
        print("gamma_family_complex", z, special.gamma(z), special.rgamma(z),
              special.loggamma(z))
    for z in (3.0 + 4.0j, -2.75 + 0.625j, 0.5 + 0.25j,
              -3.0 + 1.0e-8j, -3.0 + 1.0e-300j, 100.0 + 200.0j):
        print("digamma_complex", z, special.digamma(z))
    for a, d in ((1.0, 1), (2.5, 2), (3.5, 3), (10.0, 5), (100.0, 8)):
        print("multigammaln", a, d, special.multigammaln(a, d))
    mp.mp.dps = 80
    for x in (-170.75, -171.5, -172.5, -177.5, -178.5, 172.0, 179.0):
        xx = mp.mpf(str(x))
        print("gamma_family_mpmath", x, mp.gamma(xx), 1 / mp.gamma(xx))
    for a, b in ((0.1, 0.2), (2.0, 3.0), (100.0, 200.0), (1e-3, 2e3)):
        print("betaln", a, b, special.betaln(a, b))
    for n, x in ((1, 1.0), (2, 1.0), (3, 0.5), (1, -0.5),
                 (4, -1.5), (8, -20.5), (10, 3.25)):
        print("polygamma", n, x, special.polygamma(n, x))
    for a, m in ((0.5, 0.5), (-0.5, 1.5), (-2.0, 3.0), (-3.0, -2.0),
                 (10.25, -3.5), (100.5, 20.25)):
        print("poch", a, m, special.poch(a, m))
    for a, b in ((-0.5, 1.5), (-0.5, 0.5), (-2.0, 2.0),
                 (-4.0, 3.0), (-2.25, 3.5)):
        print("beta_extended", a, b, special.beta(a, b), special.betaln(a, b))
    for n, k in ((5.0, 2.0), (5.5, 6.0), (-2.5, 3.0),
                 (10.25, 3.5), (100.5, 20.25)):
        print("binom", n, k, special.binom(n, k))
    for n in (-3.0, -0.5, 0.0, 3.5, 10.0, 170.5, 171.0):
        print("factorial", n, special.factorial(n))
    for n in (-5, 0, 1, 7, 8, 300, 301):
        print("factorial2", n, special.factorial2(n))
    for n, k in ((0, 3), (5, 3), (17, 4), (170, 1), (171, 1),
                 (1000, 999), (1000, 1001)):
        print("factorialk", n, k, special.factorialk(n, k))
    for n, k in ((10.0, 3.0), (10.0, 3.5), (3.0, 4.0),
                 (2.2, 2.1), (-3.0, 2.0)):
        print("comb_perm", n, k, special.comb(n, k), special.perm(n, k))
    for n, k in ((10.0, 3.0), (0.0, 0.0), (1.0, 2.0),
                 (2.2, 3.1), (-3.0, 2.0)):
        print("comb_repetition", n, k, special.comb(n, k, repetition=True))

    for n in (0, 1, 4, 20, 50, 100):
        print("bernoulli_sequence", n, special.bernoulli(n).tolist())
    for n in (0, 2, 10, 22, 50, 100):
        print("euler_sequence", n, special.euler(n).tolist())
    for n, k in ((0, 0), (10, 3), (50, 25), (100, 50)):
        print("stirling2_scipy", n, k, special.stirling2(n, k))

    for a, b, x in ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0),
                     (1.0, 1.0, 1.0), (1.5, 1.1, 2.5),
                     (0.5, 0.25, 10.0), (0.1, 2.0, 100.0)):
        print("wright_bessel", a, b, x, special.wright_bessel(a, b, x))
        print("log_wright_bessel", a, b, x, special.log_wright_bessel(a, b, x))

    for x in (0.1, 1.0, 5.0, 10.0, 25.0, 100.0):
        print("kelvin_scipy", x, special.ber(x), special.bei(x), special.ker(x),
              special.kei(x), special.berp(x), special.beip(x), special.kerp(x),
              special.keip(x))
    print("kelvin_zeros_scipy", *(arr.tolist() for arr in special.kelvin_zeros(8)))

    for m, q in ((0, 1.0), (1, 1.0), (2, 1.0), (3, 10.0), (8, 100.0)):
        print("mathieu_characteristic", m, q, special.mathieu_a(m, q),
              special.mathieu_b(m, q) if m > 0 else None)
    for m, q, x in ((1, 1.0, 30.0), (2, 10.0, 45.0), (8, 100.0, 12.3),
                    (3, -10.0, 30.0)):
        print("mathieu_angular", m, q, x, special.mathieu_cem(m, q, x),
              special.mathieu_sem(m, q, x))
    print("mathieu_even_coef", special.mathieu_even_coef(0, 1.0).tolist())
    print("mathieu_odd_coef", special.mathieu_odd_coef(1, 1.0).tolist())
    for m, q, x in ((0, 0.1, 0.5), (1, 1.0, 0.5), (2, 10.0, 1.0),
                    (3, 1.0, 0.0)):
        print("mathieu_radial_even", m, q, x, special.mathieu_modcem1(m, q, x),
              special.mathieu_modcem2(m, q, x))
        if m > 0:
            print("mathieu_radial_odd", m, q, x, special.mathieu_modsem1(m, q, x),
                  special.mathieu_modsem2(m, q, x))

    for v, x in ((0.3, -3.0), (0.3, 2.0), (2.3, 2.0), (5.0, 0.0),
                 (10.2, 10.0)):
        vl, dl = special.lmbda(v, x)
        print("lmbda", v, x, vl.tolist(), dl.tolist())
    for a, lmb, nu in ((1.0, 1.0, 1.0), (1.0, -0.5, 1.0),
                       (1.0, 1.0, -0.5), (10.0, 1.0, 2.0),
                       (100.0, 1.0, 1.0), (100.0, 0.5, 2.3)):
        print("besselpoly_scipy", a, lmb, nu, special.besselpoly(a, lmb, nu))
        mp_a, mp_lmb, mp_nu = map(mp.mpf, (a, lmb, nu))
        reference = mp.quad(lambda t: t**mp_lmb * mp.besselj(mp_nu, 2*mp_a*t), [0, 1])
        print("besselpoly_mpmath", a, lmb, nu, reference)
    for nt in (1, 5, 10, 76, 100):
        print("erf_zeros", nt, special.erf_zeros(nt)[-1])

    for x in (0.0, 1e-12, 0.1, 0.5, 1.0, 2.0, 10.0, 1e6, 1e100):
        print("spence", x, special.spence(x))
    for z in (1 - 1j, 1.5 + 2j, 3j, -10 - 5j,
              -2 + 1e-8j, -2 - 1e-8j, 1e100 + 1e100j):
        print("spence_complex", z, special.spence(z))
    for signed_zero in (0.0, -0.0):
        z = complex(-2.0, signed_zero)
        print("spence_cut", repr(z), special.spence(z))
    for s_value in (-100.5, -19.5, -5.0, -3.0, -2.0, -1.0, -0.5,
                    0.0, 0.1, 0.5, 0.9, 2.0, 10.0, 100.0):
        print("zeta_riemann", s_value, special.zeta(s_value), special.zetac(s_value))
    for s_value, q_value in ((2.0, 0.5), (3.0, 0.5), (2.0, 2.0),
                             (10.0, 2.0), (2.0, -0.5), (3.0, -0.5),
                             (4.0, -3.5), (2.0, -0.1)):
        print("zeta_hurwitz", s_value, q_value, special.zeta(s_value, q_value))
    mp.mp.dps = 400
    print("zetac_mpmath", 1000.0, mp.zeta(1000) - 1)
    mp.mp.dps = 80
    for s_value in (-100.5, -19.5, -5.0, -1.0, -0.5, -1e-6, 0.5):
        numerator, denominator = float(s_value).as_integer_ratio()
        ss = mp.mpf(numerator) / denominator
        print("zeta_mpmath", s_value, mp.zeta(ss))
    for x in (0.5, 2.0, 10.0, 1e6):
        xx = mp.mpf(x)
        print("spence_mpmath", x, mp.polylog(2, 1 - xx))
    for z in (-2 + 1e-16j, -2 - 1e-16j, 0.5 + 0.8660254037844386j,
              1e100 + 1e100j):
        zz = mp.mpc(z.real, z.imag)
        print("spence_complex_mpmath", z, mp.polylog(2, 1 - zz))

    for x in (-4.75, -2.5, -1.5, -0.5, 0.1, 0.5, 1.0, 5.0, 50.0):
        print("digamma", x, special.digamma(x))
    for x in (-30.0, -10.0, -2.0, 0.0, 2.0, 10.0):
        print("ndtr", x, special.ndtr(x), special.log_ndtr(x))
    for p in (1e-200, 1e-20, 0.001, 0.5, 0.999, 1.0 - 1e-12):
        print("ndtri", p, special.ndtri(p))

    for x in (-30.0, -2.0, 0.0, 2.0, 30.0):
        print("logistic", x, special.expit(x), special.log_expit(x))
    for p in (1e-200, 1e-12, 0.1, 0.5, 0.9, 1.0 - 1e-12):
        print("logit", p, special.logit(p))

    for x, y in ((0.0, 0.0), (0.1, 0.2), (2.0, 3.0), (10.0, 0.5)):
        print("entropy", x, y, special.entr(x), special.rel_entr(x, y), special.kl_div(x, y))

    for x in (-100.0, -1.0, -1e-10, 0.0, 1e-10, 1.0, 100.0):
        print("stable", x, special.exprel(x), special.cosm1(x))

    for x in (1e-300, 1e-12, 0.1, 0.5, 1.0, 10.0, 64.0, 100.0, 500.0, 700.0):
        print("exp_integrals", x, special.exp1(x), special.expi(x))
    for z in (0.5 + 0.25j, -1.0 + 0.5j, 3.0 - 2.0j, -10.0 + 0.1j,
              10.0 + 10.0j, -20.0 + 1e-8j, 10.0 + 1e-12j):
        print("complex_exp_integrals", z, special.exp1(z), special.expi(z))
    for signed_zero in (0.0, -0.0):
        z = complex(-1.0, signed_zero)
        print("complex_exp_cut", repr(z), special.exp1(z), special.expi(z))

    for n, x in ((0, 0.5), (1, 0.5), (2, 0.5), (10, 20.0),
                 (20, 1.0), (100, 0.1), (5000, 0.5)):
        print("expn", n, x, special.expn(n, x))
    for x in (-100.0, -16.0, -10.0, -1.0, -1e-12, 1e-12,
              1.0, 10.0, 16.0, 100.0, 1e6):
        print("sici", x, *special.sici(x))
    for x in (-100.0, -10.0, -1.0, -1e-12, 1e-12, 1.0, 10.0, 100.0, 700.0):
        print("shichi", x, *special.shichi(x))
    for z in (2.5 + 3j, -2.0 + 3j, -2.0 - 3j, 1e-12 + 2e-12j,
              9.0 + 2j, 20.0 + 15j, 2j, -2j):
        print("complex_sici", z, *special.sici(z))
        print("complex_shichi", z, *special.shichi(z))
    for signed_zero in (0.0, -0.0):
        z = complex(-2.0, signed_zero)
        print("complex_sici_cut", repr(z), *special.sici(z))
        print("complex_shichi_cut", repr(z), *special.shichi(z))
    mp.mp.dps = 80
    for x in (710.0, 715.0, 716.0):
        print("expi_mpmath", x, mp.ei(x))

    for x, lam in ((0.25, -2.0), (0.25, 0.0), (0.25, 1e-10), (2.0, 0.5), (10.0, 2.0)):
        y = special.boxcox(x, lam)
        print("boxcox", x, lam, y, special.inv_boxcox(y, lam))

    for v, z in ((2.5, 3.0), (3.0, -100.0), (3.0, 10000.0),
                 (-0.5, -100.0), (-0.5, 100.0),
                 (2.5, 3.0 + 4.0j), (3.0, -100.0 + 20.0j),
                 (-0.5, -100.0 + 20.0j), (10.0, -250.0 + 15.0j)):
        print("hyp0f1", v, z, special.hyp0f1(v, z))

    for a, b, z in ((2.5, 3.5, 4.0), (2.3, 4.7, -100.0),
                      (1.0, 2.0, 700.0), (-0.5, 2.5, -100.0),
                      (-3.0, 2.5, 4.0), (-2.0, -3.0, 2.0),
                      (2.5, 3.5, 3.0 + 4.0j),
                      (2.3, 4.7, -50.0 + 30.0j),
                      (2.3, 4.7, 100.0 + 50.0j),
                      (2.3, 4.7, -100.0 + 50.0j),
                      (-1.5, 3.2, 20.0 + 10.0j)):
        print("hyp1f1", a, b, z, special.hyp1f1(a, b, z))

    for x in (-100.0, -20.0, -10.0, -2.0, -1.00000001, -1.0,
              -0.99999999, -0.1, 0.0, 0.1, 0.99999999, 1.0,
              1.00000001, 2.0, 10.0, 100.0, 1000.0):
        print("airy", x, *special.airy(x))
        print("airye", x, *special.airye(x))

    for z in (0.5 + 0.25j, 1.0 + 2.0j, 3.0 - 2.0j,
              -4.0 + 1.0j, 5.0 + 5.0j, 10.0 + 3.0j):
        print("airy_complex", z, *special.airy(z))

    for x in (-20.0, -10.0, -5.0, -1.0, 0.0, 1.0, 5.0, 8.0, 10.0, 20.0, 100.0):
        print("itairy_scipy", x, *special.itairy(x))

    mp.mp.dps = 80
    for x in (1.0, 5.0, 8.0, 10.0, 12.0, 20.0):
        xx = mp.mpf(x)
        print("itairy_mpmath", x,
              mp.quad(lambda t: mp.airyai(t), [0, xx]),
              mp.quad(lambda t: mp.airybi(t), [0, xx]),
              mp.quad(lambda t: mp.airyai(-t), [0, xx]),
              mp.quad(lambda t: mp.airybi(-t), [0, xx]))

    print("ai_zeros", *(arr.tolist() for arr in special.ai_zeros(5)))
    print("bi_zeros", *(arr.tolist() for arr in special.bi_zeros(5)))

    for x, y in ((0.0, 1.0), (1.0, 2.0), (2.0, -1.0)):
        print("elliprc", x, y, special.elliprc(x, y))
    for x, y, z in ((0.0, 1.0, 2.0), (1.0, 2.0, 3.0), (0.0, 2.0, 1.0)):
        print("elliprf_rd_rg", x, y, z, special.elliprf(x, y, z),
              special.elliprd(x, y, z), special.elliprg(x, y, z))
    for x, y, z, p in ((0.0, 1.0, 2.0, 3.0), (1.0, 2.0, 3.0, 4.0),
                       (1.0, 2.0, 3.0, -1.0)):
        print("elliprj", x, y, z, p, special.elliprj(x, y, z, p))

    for x, y in ((1.2 + 3.4j, 5.0 + 0.0j),
                 (1.0 + 2.0j, -2.0 + 0.0j)):
        print("elliprc_complex", x, y, special.elliprc(x, y))
    for x, y, z in ((1.2 + 3.4j, 5.0 + 0.0j, 6.0 + 0.0j),
                    (-1.0 + 2.0j, -1.0 - 2.0j, 3.0 + 0.0j),
                    (0.0 + 0.0j, 1.0 + 2.0j, 3.0 + 0.0j)):
        print("elliprf_rd_rg_complex", x, y, z, special.elliprf(x, y, z),
              special.elliprd(x, y, z), special.elliprg(x, y, z))
    for x, y, z, p in ((1.2 + 3.4j, 5.0 + 0.0j, 6.0 + 0.0j, 7.0 + 0.0j),
                       (-1.0 + 2.0j, -1.0 - 2.0j, 3.0 + 0.0j, 4.0 + 1.0j),
                       (1.0 + 0.0j, 2.0 + 0.0j, 3.0 + 0.0j, -0.5 + 0.25j)):
        print("elliprj_complex", x, y, z, p, special.elliprj(x, y, z, p))

    mp.mp.dps = 80
    for scale in (mp.mpf("1e-200"), mp.mpf("1e200")):
        print("carlson_mpmath", scale,
              mp.elliprc(scale, 2 * scale),
              mp.elliprf(scale, 2 * scale, 3 * scale),
              mp.elliprd(scale, 2 * scale, 3 * scale),
              mp.elliprj(scale, 2 * scale, 3 * scale, 4 * scale),
              mp.elliprg(scale, 2 * scale, 3 * scale))

    for m in (-1e200, -10.0, -2.0, 0.0, 0.5, 0.9, 1.0 - 1e-12, 1.0):
        print("legendre_complete", m, special.ellipk(m), special.ellipe(m))
    for p in (5e-324, 1e-300, 1e-100, 1e-12, 0.1, 1.0, 1e100):
        print("ellipkm1", p, special.ellipkm1(p))
    for phi, m in ((0.1, 0.5), (1.0, 0.5), (2.0, 0.5), (-2.0, 0.5),
                   (10.0, -2.0), (1000.0, 0.9)):
        print("legendre_incomplete", phi, m, special.ellipkinc(phi, m),
              special.ellipeinc(phi, m))

    mp.mp.dps = 100
    for p_float in (5e-324, 1e-300, 1e-100):
        numerator, denominator = p_float.as_integer_ratio()
        p_value = mp.mpf(numerator) / denominator
        print("ellipkm1_mpmath", p_float, mp.elliprf(0, p_value, 1))
    for m_text in ("-1e100", "-1e200"):
        m_value = mp.mpf(m_text)
        print("ellipe_mpmath", m_text, mp.ellipe(m_value))
    for u, m in ((0.5, 0.0), (0.5, 0.5), (2.0, 0.5), (-2.0, 0.5),
                 (10.0, 0.5), (1.0, 0.2), (3.0, 0.9),
                 (10.0, 1e-12), (10.0, 1.0 - 1e-12)):
        print("ellipj", u, m, *special.ellipj(u, m))

    mp.mp.dps = 100
    for u_float, m_float in ((100.0, 1.0 - 1e-12),
                             (500.0, 1.0 - 1e-12),
                             (1000.0, 1.0 - 1e-12),
                             (1e6, 0.5), (1e6, 0.9)):
        un, ud = u_float.as_integer_ratio()
        mn, md = m_float.as_integer_ratio()
        uu = mp.mpf(un) / ud
        mm = mp.mpf(mn) / md
        print("ellipj_mpmath", u_float, m_float,
              mp.ellipfun("sn", uu, mm),
              mp.ellipfun("cn", uu, mm),
              mp.ellipfun("dn", uu, mm))

    for k in range(1, 6):
        print("airy_zeros_mpmath", k, mp.airyaizero(k), mp.airyaizero(k, 1),
              mp.airybizero(k), mp.airybizero(k, 1))

    values = np.array([1000.0, 1001.0, 999.0, -1000.0])
    print("logsumexp", special.logsumexp(values))
    print("softmax", special.softmax(values))
    print("log_softmax", special.log_softmax(values))

    for n in (0, 2, 10, 100, 1000):
        print("bessel_zeros", n,
              special.jn_zeros(n, 5).tolist(),
              special.jnp_zeros(n, 5).tolist(),
              special.yn_zeros(n, 5).tolist(),
              special.ynp_zeros(n, 5).tolist())
    zeros_1200 = special.jnyn_zeros(0, 1200)
    print("bessel_zeros_1200_last", *(float(z[-1]) for z in zeros_1200))
    jnjnp_1200 = special.jnjnp_zeros(1200)
    print("jnjnp_zeros_20", *(arr[:20].tolist() for arr in jnjnp_1200))
    print("jnjnp_zeros_1200_last", *(arr[-1].item() for arr in jnjnp_1200))
    print("jnjnp_zeros_1200_max_order", int(np.max(jnjnp_1200[1])))
    for name in ("y0_zeros", "y1_zeros", "y1p_zeros"):
        roots, values = getattr(special, name)(4)
        print(name, roots.tolist(), values.tolist())
    y0_ext = special.yn_zeros(0, 141)
    y1_ext = special.yn_zeros(1, 141)
    y1p_ext = special.ynp_zeros(1, 141)
    print("y0_zeros_real_141", y0_ext[-1], special.yvp(0, y0_ext[-1]))
    print("y1_zeros_real_141", y1_ext[-1], special.yvp(1, y1_ext[-1]))
    print("y1p_zeros_real_141", y1p_ext[-1], special.yv(1, y1p_ext[-1]))

    for x in (-10.0, -1.0, 0.0, 0.1, 1.0, 10.0, 100.0, 700.0):
        print("bessel_j", x, special.j0(x), special.j1(x))
    for x in (0.1, 1.0, 10.0, 100.0, 700.0):
        print("bessel_y", x, special.y0(x), special.y1(x))
        print("bessel_k", x, special.k0(x), special.k1(x))
        print("bessel_ke", x, special.k0e(x), special.k1e(x))
    for x in (-10.0, -1.0, 0.0, 0.1, 1.0, 10.0, 100.0, 700.0):
        print("bessel_i", x, special.i0(x), special.i1(x))

    for n, x in ((-3, 0.1), (-2, 10.0), (4, 1.0), (20, 100.0)):
        print("bessel_yn_kn", n, x, special.yn(n, x), special.kn(n, x))
    for v, x in ((2.5, 1.2), (-0.5, 1.0), (-1.5, 1.0), (-1.5, 10.0),
                 (2.0, -1.2), (-3.0, -1.2), (20.0, 100.0)):
        print("bessel_iv_ive", v, x, special.iv(v, x), special.ive(v, x))
    for v, x in ((2.5, 1.2), (-2.5, 1.2), (5.5, 0.1), (20.0, 100.0),
                 (2.5, 700.0)):
        print("bessel_kv_kve", v, x, special.kv(v, x), special.kve(v, x))

    for v, x in ((0.5, 1.0), (20.75, 5.0), (100.5, 1000.0),
                 (-1.5, 1.0), (-0.999999, 1e-6), (-1.000001, 1e-6),
                 (3.0, -2.0), (-3.0, -2.0)):
        print("bessel_jv_jve", v, x, special.jv(v, x), special.jve(v, x))
    for v, x in ((0.25, 0.1), (5.5, 10.0), (100.5, 1000.0),
                 (-1.5, 1.0), (-0.999999, 1e-6), (-100.5, 0.1)):
        print("bessel_yv_yve", v, x, special.yv(v, x), special.yve(v, x))

    for v, x, n in ((0.0, 1.0, 1), (0.0, 1.0, 2), (2.5, 3.0, 4),
                    (-1.5, 2.0, 2), (0.0, 1.0, 20)):
        print("bessel_derivatives", v, x, n, special.jvp(v, x, n),
              special.yvp(v, x, n), special.ivp(v, x, n), special.kvp(v, x, n))

    for n, x in ((0, 0.1), (1, 0.1), (2, 1.0), (5, 1.0),
                 (5, 10.0), (20, 5.0), (20, 100.0), (100, 100.0),
                 (5, -1.0), (5, -10.0)):
        print("spherical_jn", n, x, special.spherical_jn(n, x),
              special.spherical_jn(n, x, derivative=True))
        print("spherical_yn", n, x, special.spherical_yn(n, x),
              special.spherical_yn(n, x, derivative=True))
        print("spherical_in", n, x, special.spherical_in(n, x),
              special.spherical_in(n, x, derivative=True))
        print("spherical_kn", n, x, special.spherical_kn(n, x),
              special.spherical_kn(n, x, derivative=True))

    mp.mp.dps = 80
    for x in (710.0, 713.0):
        print("bessel_i_mpmath", x, mp.besseli(0, x), mp.besseli(1, x))


if __name__ == "__main__":
    main()
