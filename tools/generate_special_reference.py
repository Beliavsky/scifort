# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Print SciPy reference values for Milestone 6 special-function tests."""

import numpy as np
import scipy
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

    for x in (-4.75, -2.5, -0.5, 0.5, 1.0, 3.5, 20.0):
        print("gammaln", x, special.gammaln(x))
    for a, b in ((0.1, 0.2), (2.0, 3.0), (100.0, 200.0), (1e-3, 2e3)):
        print("betaln", a, b, special.betaln(a, b))
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

    for x, lam in ((0.25, -2.0), (0.25, 0.0), (0.25, 1e-10), (2.0, 0.5), (10.0, 2.0)):
        y = special.boxcox(x, lam)
        print("boxcox", x, lam, y, special.inv_boxcox(y, lam))

    values = np.array([1000.0, 1001.0, 999.0, -1000.0])
    print("logsumexp", special.logsumexp(values))
    print("softmax", special.softmax(values))
    print("log_softmax", special.log_softmax(values))


if __name__ == "__main__":
    main()
