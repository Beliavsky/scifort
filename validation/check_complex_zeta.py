#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare native complex zeta with 90-digit mpmath and SciPy references."""
import argparse
from pathlib import Path
import random
import re
import subprocess

import mpmath as mp
import numpy as np
import scipy
from scipy.special import zeta


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    args = parser.parse_args()
    mp.mp.dps = 90
    rng = random.Random(131126)
    points = [complex(x, y) for x in (-100, -50, -20, -4, -2, -0.251,
              -0.25, -0.249, -0.01, 0, 0.5, 1, 2, 10, 32, 100)
              for y in (1e-8, 0.1, 1, 10, 25, 100)]
    points += [complex(rng.uniform(-30, 30), rng.uniform(-200, 200)) for _ in range(120)]
    points += [complex(x, y) for x in (-10, 0.5, 2)
               for y in (1000, 10000, 100000, -100000)]
    points += [complex(1 + d, d) for d in (1e-4, 1e-8, 1e-12, 1e-100, 1e-300)]
    points += [complex(x, y) for x in (-20, -4, -2, 0)
               for y in (1e-100, -1e-100)]
    points += [complex(0.5, float(mp.im(mp.zetazero(n)))) for n in range(1, 6)]
    points += [complex(10000, 1), complex(0.5, -10000)]
    # Independent verification of every adapted Bernoulli coefficient.
    source = Path(__file__).resolve().parents[1] / "src/special/complex_zeta.f90"
    block = source.read_text().split("bernoulli_coeff(25) = [", 1)[1].split("]", 1)[0]
    coeff = re.findall(r"([-+]?\d+\.\d+e[-+]\d+)_dp", block)
    assert len(coeff) == 25
    for k, text in enumerate(coeff, 1):
        expected = mp.bernoulli(2*k)/mp.factorial(2*k)
        assert abs(mp.mpf(text)/expected - 1) < mp.mpf("1e-35")
        assert float(text) == float(expected)
    payload = str(len(points)) + "\n" + "\n".join(
        f"{s.real:.17g} {s.imag:.17g}" for s in points) + "\n"
    result = subprocess.run([str(args.executable.resolve())], input=payload,
                            text=True, capture_output=True, check=True)
    values = [complex(*map(float, line.split())) for line in result.stdout.splitlines()]
    assert len(values) == len(points)
    worst = (0.0, None)
    scipy_worst = 0.0
    scipy_outliers = 0
    scipy_outlier_points = []
    ordinary_worst = 0.0
    for s, got in zip(points, values):
        ref = complex(mp.zeta(mp.mpc(s.real, s.imag)))
        # Absolute control near zeros; relative control for large magnitudes.
        bound = 5e-13 if abs(s.imag) <= 200 else 3e-9
        error = abs(got-ref)/(1+abs(ref))
        assert np.isfinite(got), (s, got)
        assert error <= bound, (s, got, ref, error, bound)
        if abs(s.imag) <= 200:
            ordinary_worst = max(ordinary_worst, error)
        if error > worst[0]:
            worst = (error, s)
        other = complex(zeta(s))
        if np.isfinite(other):
            delta = abs(got-other)/(1+abs(ref))
            if abs(other-ref)/(1+abs(ref)) <= bound:
                assert delta <= 2*bound, (s, got, other, delta)
            else:
                # mpmath remains the oracle when SciPy loses accuracy, notably
                # near zero. Never relax the native-versus-mpmath assertion.
                scipy_outliers += 1
                scipy_outlier_points.append(s)
            scipy_worst = max(scipy_worst, delta)
    print(f"PASS: {len(points)} points, 25 coefficients; mpmath {mp.__version__}, "
          f"SciPy {scipy.__version__}")
    print(f"Maximum |error|/(1+|reference|): {worst[0]:.6g} at {worst[1]}")
    print(f"Maximum normalized SciPy difference: {scipy_worst:.6g}")
    print(f"Maximum normalized error for |Im(s)| <= 200: {ordinary_worst:.6g}")
    print(f"SciPy points outside the mpmath bound: {scipy_outliers}")
    print(f"SciPy outlier inputs: {scipy_outlier_points}")


if __name__ == "__main__":
    main()
