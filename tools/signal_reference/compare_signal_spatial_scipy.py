#!/usr/bin/env python3
"""Compare SciFort order_filter and sepfir2d with SciPy 1.17.

Build spatial_reference_driver.f90 with kinds.f90 and signal_spatial.f90.
Usage: python tools/signal_reference/compare_signal_spatial_scipy.py ./spatial_driver
"""
from __future__ import annotations

import subprocess
import sys
import numpy as np
import scipy
from scipy import signal, ndimage


def payload(values):
    return ' '.join(format(float(v), '.17g') for v in np.ravel(values, order='F'))


def main():
    if len(sys.argv) != 2:
        raise SystemExit('usage: compare_signal_spatial_scipy.py FORTRAN_DRIVER')
    rng = np.random.default_rng(20261019)
    cases = []
    scipy_order_agreement = 0
    scipy_order_disagreement = 0
    # All six (real/integer input, integer/real/logical footprint) combinations.
    for ndim in (1, 2, 3):
        for xtype in (0, 1):
            for dtype in (0, 1, 2):
                for i in range(28):
                    shape = tuple(int(q) for q in rng.integers(1, 8, size=ndim))
                    footprint = tuple(int(rng.choice([q for q in [1, 3, 5, 7] if q <= n])) for n in shape)
                    mask = rng.choice([False, False, True], size=footprint)
                    if not np.any(mask):
                        mask.flat[0] = True
                    # Non-binary values confirm that footprint inclusion means nonzero.
                    if dtype == 0:
                        dom = mask.astype(np.int32)*rng.choice([-4, 2, 7], size=footprint)
                    elif dtype == 1:
                        dom = mask.astype(float)*rng.uniform(0.2, 2.0, size=footprint)
                    else:
                        dom = mask
                    rank = int(rng.integers(0, np.count_nonzero(dom)))
                    if xtype:
                        x = rng.integers(-100, 100, size=shape, dtype=np.int32)
                    else:
                        x = rng.normal(size=shape)
                    # SciPy 1.17's median-filter shortcut can mishandle
                    # sparse footprints. Validate the documented operation
                    # against ndimage.generic_filter in those cases.
                    expected = ndimage.generic_filter(x.astype(float),
                        lambda win: np.partition(win, rank)[rank],
                        footprint=np.asarray(dom) != 0, mode='constant')
                    scipy_result = signal.order_filter(x, dom, rank)
                    if np.array_equal(scipy_result, expected):
                        scipy_order_agreement += 1
                    else:
                        scipy_order_disagreement += 1
                    cases.append((0, xtype, dtype, ndim, rank, shape, footprint,
                                  x, dom, expected))
    # 2D separable filtering: real, complex and integer inputs; odd asymmetric taps.
    # SciPy rejects complex taps with real/integer images. For these extensions,
    # explicitly promote the reference image to complex before calling SciPy.
    for xtype in (0, 1, 2):
        for row_complex in (0, 1):
            for col_complex in (0, 1):
                for i in range(32):
                    shape = tuple(int(v) for v in rng.integers(1, 9, size=2))
                    lengths = (int(rng.choice([v for v in [1, 3, 5, 7] if v <= shape[1]])),
                               int(rng.choice([v for v in [1, 3, 5, 7] if v <= shape[0]])))
                    if xtype == 2:
                        x = rng.integers(-7, 9, size=shape, dtype=np.int32)
                    else:
                        x = rng.normal(size=shape)
                        if xtype == 1:
                            x = x + 1j*rng.normal(size=shape)
                    row = rng.normal(size=lengths[0])
                    col = rng.normal(size=lengths[1])
                    if row_complex:
                        row = row + 1j*rng.normal(size=lengths[0])
                    if col_complex:
                        col = col + 1j*rng.normal(size=lengths[1])
                    ref_x = x.astype(np.complex128) if (row_complex or col_complex) and xtype != 1 else x
                    expected = signal.sepfir2d(ref_x, row, col)
                    cases.append((1, xtype, row_complex + 2*col_complex, 2, 0,
                                  shape, lengths, x, (row, col), expected))

    content = [str(len(cases))]
    for method, xtype, dtype, ndim, rank, shape, footprint, x, domain, _ in cases:
        content.append(f'{method} {xtype} {dtype} {ndim} {rank}')
        content.append(' '.join(map(str, shape + (1,)*(3-ndim))))
        content.append(' '.join(map(str, footprint + (1,)*(3-ndim))))
        content.append(payload(np.real(x)))
        content.append(payload(np.imag(x)))
        if method == 0:
            content.append(payload(domain))
        else:
            row, col = domain
            content.extend((payload(np.real(row)), payload(np.imag(row)),
                            payload(np.real(col)), payload(np.imag(col))))
    run = subprocess.run([sys.argv[1]], input='\n'.join(content)+'\n',
                         text=True, capture_output=True, check=True)
    lines = run.stdout.splitlines()
    if len(lines) != len(cases):
        raise AssertionError(f'{len(lines)} lines for {len(cases)} cases; {run.stderr[-400:]}')
    max_error = 0.0
    for i, (line, case) in enumerate(zip(lines, cases), 1):
        method, xtype, dtype, ndim, rank, shape, footprint, x, domain, expected = case
        values = np.asarray([float(v.replace('D','E')) for v in line.split()])
        if int(values[0]) != 0:
            raise AssertionError(f'case {i}: Fortran status {values[0]}')
        actual = (values[1::2] + 1j*values[2::2]).reshape(shape, order='F')
        try:
            np.testing.assert_allclose(actual, expected, atol=2e-11, rtol=5e-12)
        except AssertionError as exc:
            raise AssertionError(f'case {i}, method={method}, xtype={xtype}, '
                                 f'dtype={dtype}, shape={shape}, footprint={footprint}, '
                                 f'rank={rank}: {exc}') from exc
        if actual.size:
            max_error = max(max_error,float(np.max(abs(actual-expected))))
    print(f'SciPy {scipy.__version__}: {len(cases)}/{len(cases)} spatial reference cases pass; '
          f'max abs error {max_error:.3e}; '
          f'order_filter agrees with SciPy directly in {scipy_order_agreement}/504 cases, '
          f'{scipy_order_disagreement} sparse-footprint cases checked against generic_filter')


if __name__ == '__main__':
    main()
