#!/usr/bin/env python3
"""Independent Hermite/BPoly `orders` validation against SciPy 1.17.

SciPy 1.17's `BPoly.from_derivatives` contains an asymmetric-constraint
selection defect (`len(y2)` instead of `len(y1)`). We validate *every* case
against an independent Hermite linear system, and also compare directly
against SciPy's public API wherever SciPy agrees with that reference.

Compiles fresh checked and optimized GNU Fortran binaries, without FPM.
"""
from __future__ import annotations
import argparse
import subprocess
import tempfile
from pathlib import Path
import numpy as np
from scipy.interpolate import PPoly, BPoly
from scipy.special import factorial

ROOT = Path(__file__).resolve().parents[2]
SOURCES = [
    'src/kinds.f90', 'src/interpolate_piecewise.f90',
    'src/interpolate_bspline.f90', 'src/interpolate_regular_grid.f90',
    'src/interpolate_barycentric.f90', 'src/interpolate_global.f90', 'src/interpolate_kdtree.f90', 'src/interpolate_rbf.f90', 'src/interpolate_legacy.f90', 'src/interpolate_scattered.f90', 'src/interpolate_clough_tocher.f90', 'src/interpolate_fitpack_bivariate.f90', 'src/interpolate_fitpack_compat.f90', 'src/interpolate.f90',
]


def fstr(a):
    return ' '.join(format(float(t), '.17g') for t in np.ravel(a, order='F'))


def samples(seed, number):
    rng = np.random.default_rng(seed)
    for i in range(number):
        n = 2 + i % 6
        nch = 1 + (i // 2) % 4
        counts = rng.integers(1, 5, size=n, dtype=int)
        if i % 6 == 0:
            counts[:] = 3  # balanced cases are unaffected by SciPy bug
        maxd = int(max(counts))
        x = np.r_[0.0, np.cumsum(rng.uniform(0.5, 1.5, n - 1))]
        x += rng.uniform(-3, 3)
        y = rng.normal(size=(n, maxd, nch))
        for j, count in enumerate(counts):
            y[j, count:, :] = 0.0
        maxdegrees = counts[:-1] + counts[1:] - 1
        scalar = bool(i % 7 == 0)
        if scalar:
            order = int(rng.integers(1, min(maxdegrees) + 1))
            orders = np.full(n - 1, order, dtype=int)
        else:
            orders = np.array([rng.integers(1, m + 1) for m in maxdegrees], dtype=int)
            order = 0
        q = np.r_[rng.uniform(x[0] - 0.3, x[-1] + 0.3, 9), x]
        yield x, y, counts, orders, scalar, order, q


def left_right(nl, nr, deg):
    total = deg + 1
    left = min(total // 2, nl)
    right = min(total - left, nr)
    left = min(total - right, nl)
    assert left + right == total
    return left, right


def independent_bernstein(x, y, counts, orders):
    """Solve endpoint derivative conditions in t=(x-x[i])/dx.

    Construct a local monomial Hermite polynomial via a dense linear solve,
    then use SciPy's independent power-to-Bernstein conversion.
    """
    highest = max(max(counts[:-1] + counts[1:] - 1), max(orders))
    pieces = []
    for i, degree in enumerate(orders):
        degree = int(degree)
        nl, nr = left_right(counts[i], counts[i + 1], degree)
        dx = x[i + 1] - x[i]
        mat = np.zeros((degree + 1, degree + 1))
        rhs = np.empty((degree + 1, y.shape[2]))
        row = 0
        for j in range(nl):
            mat[row, j] = factorial(j, exact=True)
            rhs[row] = y[i, j] * dx ** j
            row += 1
        for j in range(nr):
            for p in range(j, degree + 1):
                mat[row, p] = factorial(p, exact=True) / factorial(p - j, exact=True)
            rhs[row] = y[i + 1, j] * dx ** j
            row += 1
        power = np.linalg.solve(mat, rhs)
        polys = power[::-1, :] / dx ** np.arange(degree, -1, -1)[:, None]
        basis = BPoly.from_power_basis(PPoly(polys[:, None, :], x[i:i + 2])).c[:, 0, :]
        if len(basis) < highest + 1:
            basis = BPoly._raise_degree(basis, highest + 1 - len(basis))
        pieces.append(basis)
    return BPoly(np.stack(pieces, axis=1), x)


def serialize(cases):
    data = [str(len(cases))]
    for x, y, counts, orders, scalar, order, q in cases:
        data.extend([
            f'{len(x)} {y.shape[1]} {y.shape[2]} {len(q)} {int(scalar)} {order}',
            fstr(x), fstr(counts), fstr(orders), fstr(y), fstr(q),
        ])
    return '\n'.join(data) + '\n'


def build(path, flags):
    cmd = [
        'gfortran', '-std=f2018', '-ffree-line-length-none', *flags,
        '-J' + str(path), '-I' + str(path),
        *(str(ROOT / file) for file in SOURCES),
        str(ROOT / 'tools/interpolate_reference/bernstein_orders_driver.f90'),
        '-o', str(path / 'bernstein_orders'),
    ]
    cp = subprocess.run(cmd, cwd=path, text=True, capture_output=True)
    if cp.returncode:
        raise RuntimeError(cp.stderr[-6000:])
    return path / 'bernstein_orders'


def test_binary(executable, cases):
    cp = subprocess.run([str(executable)], input=serialize(cases), text=True,
                        capture_output=True)
    if cp.returncode:
        raise RuntimeError(cp.stderr[-5000:])
    output = iter(cp.stdout.splitlines())
    err = 0.0
    nvals = 0
    direct = 0
    upstream_mismatch = 0
    for idx, (x, y, counts, orders, scalar, order, q) in enumerate(cases, 1):
        ref = independent_bernstein(x, y, counts, orders)
        given = next(output).split()
        assert given[:2] == ['CASE', str(idx)], (given, idx)
        shape = (int(given[2]), int(given[3]), y.shape[2])
        got_c = np.fromstring(next(output), sep=' ').reshape(shape, order='F')
        np.testing.assert_allclose(got_c, ref.c, atol=5e-9, rtol=1e-9,
                                   err_msg=f'local Hermite coefficients case {idx}')
        nvals += got_c.size
        err = max(err, float(np.max(np.abs(got_c - ref.c))))
        for der in range(3):
            for t in q:
                got = np.fromstring(next(output), sep=' ')
                expected = np.asarray(ref(t, nu=der)).ravel()
                np.testing.assert_allclose(got, expected, atol=2e-8, rtol=2e-8,
                                           err_msg=f'case {idx}, derivative={der}, x={t}')
                err = max(err, float(np.max(np.abs(got - expected))))
                nvals += got.size
        yi = [y[i, :counts[i], :].copy() for i in range(len(x))]
        try:
            scipy_ref = BPoly.from_derivatives(x, yi, orders=order if scalar else list(orders))
        except (ValueError, IndexError):
            upstream_mismatch += 1
            continue
        if scipy_ref.c.shape != ref.c.shape or not np.allclose(
            scipy_ref.c, ref.c, rtol=5e-9, atol=5e-9
        ):
            upstream_mismatch += 1
            continue
        np.testing.assert_allclose(got_c, scipy_ref.c, atol=5e-9, rtol=1e-9)
        direct += 1
    assert next(output, None) is None, 'unexpected additional Fortran output'
    return nvals, err, direct, upstream_mismatch


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=20261031)
    parser.add_argument('--cases', type=int, default=320)
    args = parser.parse_args()
    cases = list(samples(args.seed, args.cases))
    for label, flags in [('checked', ['-O0', '-g', '-fcheck=all', '-fbacktrace']),
                         ('optimized', ['-O3'])]:
        with tempfile.TemporaryDirectory() as folder:
            executable = build(Path(folder), flags)
            nvals, err, direct, upstream = test_binary(executable, cases)
        print(f'{label}: PASS {len(cases)} cases, {nvals} values; '
              f'{direct} direct SciPy matches, {upstream} upstream mismatch/exception; '
              f'max_abs_error={err:.3g}')


if __name__ == '__main__':
    main()
