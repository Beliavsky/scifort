# SPDX-License-Identifier: MIT
# Copyright (c) 2026 SciFort contributors
"""Compare every native Lebedev rule with SciPy and independent sphere moments."""
import argparse
from fractions import Fraction
from functools import lru_cache
from math import factorial
from pathlib import Path
import subprocess

import numpy as np
import scipy
from scipy.integrate import lebedev_rule
from scipy.special import eval_legendre

ORDERS = (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35,
          41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131)


@lru_cache(None)
def moment(a, b, c):
    """Exact rational factor multiplying the unit sphere's surface area."""
    s = a+b+c
    ratio = Fraction(factorial(2*a)*factorial(2*b)*factorial(2*c)*factorial(s),
                     factorial(a)*factorial(b)*factorial(c)*factorial(2*s+1))
    return 4*np.pi*float(ratio)


def symmetry(x, w):
    """Inversion, individual reflections, and generators of coordinate permutations."""
    rounded = np.round(x.T, 14)
    locations = {tuple(p): i for i, p in enumerate(rounded)}
    assert len(locations) == len(w), 'Duplicate nodes'
    checks = 0
    for i, point in enumerate(rounded):
        transforms = [-point, point[[1, 0, 2]], point[[1, 2, 0]]]
        for axis in range(3):
            reflected = point.copy()
            reflected[axis] *= -1
            transforms.append(reflected)
        for transformed in transforms:
            j = locations[tuple(transformed)]
            assert w[i] == w[j], (i, j, w[i], w[j])
            checks += 1
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('executable', type=Path)
    args = parser.parse_args()
    result = subprocess.run([str(args.executable.resolve())], capture_output=True,
                            text=True, check=True)
    lines = result.stdout.splitlines()
    offset = 0
    coordinates = weights = moments = symmetries = 0
    worst = 0.
    eps = np.finfo(float).eps
    for order in ORDERS:
        actual_order, m = map(int, lines[offset].split())
        assert actual_order == order
        offset += 1
        data = np.loadtxt(lines[offset:offset+m]).T
        offset += m
        x, w = data[:3], data[3]
        expected_x, expected_w = lebedev_rule(order)
        np.testing.assert_allclose(x, expected_x, atol=4*eps, rtol=0)
        np.testing.assert_allclose(w, expected_w, atol=0, rtol=8*eps)
        np.testing.assert_allclose(np.sum(x*x, axis=0), 1, atol=4*eps, rtol=0)
        np.testing.assert_allclose(np.sum(w), 4*np.pi, atol=5e-13, rtol=0)
        coordinates += x.size
        weights += w.size
        symmetries += symmetry(x, w)
        powers = np.ones((order//2+1, 3, m))
        for k in range(1, len(powers)):
            powers[k] = powers[k-1]*x*x
        for a in range(order//2+1):
            for b in range(a+1):
                for c in range(b+1):
                    if a+b+c > order//2:
                        continue
                    terms = w*powers[a, 0]*powers[b, 1]*powers[c, 2]
                    expected = moment(a, b, c)
                    error = abs(np.sum(terms)-expected)
                    bound = 3e-13*expected + 16*eps*np.sum(abs(terms))
                    assert error <= bound, (order, a, b, c, error, expected)
                    worst = max(worst, error/expected)
                    moments += 1
        # Nonconstant zonal spherical harmonics have zero surface integral.
        for degree in sorted({1, 2, order-1, order}):
            integral = np.dot(w, eval_legendre(degree, x[2]))
            assert abs(integral) <= 1e-13, (order, degree, integral)
    assert offset == len(lines), 'Unexpected trailing output'
    print(f'SciPy {scipy.__version__}: all {len(ORDERS)} rules passed; '
          f'{coordinates} coordinates, {weights} weights, '
          f'{symmetries} symmetry checks, {moments} polynomial moments.')
    print(f'Maximum polynomial relative error: {worst:.4g}')


if __name__ == '__main__':
    main()
