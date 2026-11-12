#!/usr/bin/env python3
"""Compare 684 Fortran filtfilt_nd and sosfiltfilt_nd cases with SciPy 1.17.

Build the isolated driver using GNU Fortran, no FPM/external Fortran dependencies.
Cases include ranks 2, 5, and 15, two axes, four padtypes, and eight
transfer-function/four SOS combinations of real and complex data.
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import subprocess
import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
SHAPES = [(12, 23), (4, 2, 23, 3, 2), (24, 1, 2, *([1] * 11), 2)]
assert len(SHAPES[2]) == 15
BR = np.array([.20, .12, .08])
AR = np.array([1., -.3, .05])
BC = BR.astype(complex)
AC = AR.astype(complex)
BC[1] += .06j
AC[1] -= .03j
SR = np.array([[.30, .20, .10, 1., -.27, .04],
               [.90, -.12, .08, 1., -.18, .09]])
SC = SR.astype(complex)
SC[0,1] += .04j
SC[1,4] -= .04j


def expected(cid):
    default_pad = cid >= 768
    if default_pad:
        cid -= 768
    shape = SHAPES[cid % 3]
    typ = (cid // 3) % 8
    method = (cid // 24) % 2
    mode = (cid // 48) % (3 if default_pad else 4)
    variant = 0 if default_pad else (cid // 192) % 2
    plen = None if default_pad else (cid // 384) % 2
    ax = ([1, 0][variant] if cid % 3 == 0 else
          [2, 0][variant] if cid % 3 == 1 else [0, 14][variant])
    pad = ['odd', 'even', 'constant', None][mode]
    idx = np.arange(1, np.prod(shape) + 1)
    xr = (np.sin(.117 * idx) + .0003 * idx).reshape(shape, order='F')
    xc = xr + 1j * (.16 * np.cos(.092 * idx) - .004).reshape(shape, order='F')
    if method == 0:
        coeffs = {0: (BR, AR, xr), 1: (BC, AC, xc),
                  2: (BR, AR, xc), 3: (BR, AC, xr),
                  4: (BR, AC, xc), 5: (BC, AR, xr),
                  6: (BC, AR, xc), 7: (BC, AC, xr)}
        b, a, x = coeffs[typ]
        y = signal.filtfilt(b, a, x, axis=ax, padtype=pad, padlen=plen)
    else:
        coeffs = {0: (SR, xr), 1: (SC, xc),
                  2: (SR, xc), 3: (SC, xr)}
        sos, x = coeffs[typ]
        y = signal.sosfiltfilt(sos, x, axis=ax, padtype=pad, padlen=plen)
    return y.ravel(order='F')


def main():
    cases = [cid for cid in range(768)
             if (cid // 24) % 2 == 0 or (cid // 3) % 8 <= 3]
    cases += [768 + cid for cid in range(144)
              if (cid // 24) % 2 == 0 or (cid // 3) % 8 <= 3]
    modules = ['kinds', 'signal_lfilter', 'signal_sos',
               'signal_filtfilt_gust', 'signal_filtfilt', 'signal_filtfilt_nd']
    with TemporaryDirectory(prefix='scifort-filtfilt-nd-') as tmp:
        tmp = Path(tmp)
        objects = []
        for name in modules:
            obj = tmp / (name + '.o')
            subprocess.run(['gfortran', '-std=f2018', '-O2', '-I', str(tmp),
                            '-J', str(tmp), '-c', str(ROOT / 'src' / (name + '.f90')),
                            '-o', str(obj)], check=True, cwd=tmp)
            objects.append(str(obj))
        driver = tmp / 'driver'
        subprocess.run(['gfortran', '-std=f2018', '-O2', '-I', str(tmp),
                        str(ROOT / 'tools/signal_reference/filtfilt_nd_driver.f90'),
                        *objects, '-o', str(driver)], check=True, cwd=tmp)
        proc = subprocess.run([str(driver)], input=''.join(f'{cid}\n' for cid in cases),
                              text=True, capture_output=True, check=True)
        lines = iter(proc.stdout.splitlines())
        maxerr = 0.0
        count = 0
        for cid in cases:
            label, n, status = next(lines).split()
            assert label == 'CASE' and int(status) == 0, (cid, label, status)
            got = np.array([complex(*map(float, next(lines).split()))
                            for _ in range(int(n))])
            ref = expected(cid)
            err = np.max(np.abs(got - ref))
            maxerr = max(maxerr, err)
            count += 1
            assert err < 2.e-12 * (1. + np.max(np.abs(ref))), (cid, err)
        assert next(lines, None) is None
    print(f'PASS: {count} rank-2/5/15 comparisons vs SciPy {signal.__version__ if hasattr(signal, "__version__") else "1.17"}; max absolute error {maxerr:.3g}')


if __name__ == '__main__':
    main()
