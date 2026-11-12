#!/usr/bin/env python3
"""Compare affine/displacement analytic Hessians with SciPy numerical derivatives.

SciPy does not expose interpolation Hessians, so the independently computed
map_coordinates values are differentiated with central differences at smooth
non-knot coordinates. The Fortran program produces its own analytic results.
"""
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
MODES = ['mirror', 'grid-wrap', 'nearest', 'grid-constant']
A = np.array([[0.32, 0.12], [0.09, 0.26]])
SOURCE = np.fromfunction(lambda i, j: np.sin(0.19 * (i + 1) + 0.29 * (j + 1)) +
                         0.015 * (i + 1) * (j + 1), (8, 7), dtype=float)
EPS = 2.0e-3


def evaluate(coords, degree, mode):
    return ndimage.map_coordinates(SOURCE, np.asarray(coords, dtype=float) - 1,
                                   order=degree, mode=mode, cval=0.37,
                                   prefilter=True)


def derivative(coords, degree, mode):
    unit_x = np.array([[EPS], [0.]])
    unit_y = np.array([[0.], [EPS]])
    f = evaluate(coords, degree, mode)
    fxp = evaluate(coords + unit_x, degree, mode)
    fxm = evaluate(coords - unit_x, degree, mode)
    fyp = evaluate(coords + unit_y, degree, mode)
    fym = evaluate(coords - unit_y, degree, mode)
    grad_x = (fxp - fxm) / (2 * EPS)
    hxx = (fxp - 2 * f + fxm) / EPS**2
    hxy = (evaluate(coords+unit_x+unit_y,degree,mode) -
           evaluate(coords+unit_x-unit_y,degree,mode) -
           evaluate(coords-unit_x+unit_y,degree,mode) +
           evaluate(coords-unit_x-unit_y,degree,mode)) / (4 * EPS**2)
    return f, grad_x, hxx, hxy


def main():
    with tempfile.TemporaryDirectory(prefix='scifort-spline-parameter-hessians-') as tmp:
        tmp = Path(tmp)
        objs = []
        for path in ['src/kinds.f90', 'src/signal_spline.f90', 'src/signal_spline_nd.f90']:
            obj = tmp / (Path(path).stem + '.o')
            subprocess.run(['gfortran', '-std=f2018', '-O0', '-Werror=implicit-interface',
                            '-fcheck=all', '-J', str(tmp), '-I', str(tmp), '-c',
                            str(ROOT/path), '-o', str(obj)], cwd=tmp, check=True,
                           capture_output=True)
            objs.append(str(obj))
        exe = tmp / 'parameter_hessian_reference'
        subprocess.run(['gfortran', *objs, str(ROOT/'tools/signal_reference'/
                       'spline_parameter_hessian_reference.f90'), '-o', str(exe)],
                       cwd=tmp, check=True, capture_output=True)
        records = np.loadtxt(subprocess.check_output([str(exe)], cwd=tmp).splitlines())
    count = 0
    max_value_error = 0.0
    max_hessian_error = 0.0
    for mode_idx, mode in enumerate(MODES, start=1):
        for degree in (2, 3, 5):
            part = records[(records[:, 0] == mode_idx) & (records[:, 1] == degree)]
            assert part.shape == (20, 18), part.shape
            affine = part[:, 3:5].T
            displaced = part[:, 5:7].T
            f, gx, hxx, hxy = derivative(affine, degree, mode)
            fw, gwx, hwxx, hwxy = derivative(displaced, degree, mode)
            # Columns: m, degree, k, ax, ay, wx, wy, value, dA11,
            # db1, d2/db1db1, d2/db1db2, d2/dA11db2, d2/dA11dA22,
            # warped, warped gradient 1, warped Hessian 11, warped Hessian 12
            value_err = max(np.max(abs(part[:, 7] - f)),
                            np.max(abs(part[:, 14] - fw)))
            grad_err = max(np.max(abs(part[:, 9] - gx)),
                           np.max(abs(part[:, 8] - gx * ((part[:, 2] - 1) % 5))),
                           np.max(abs(part[:, 15] - gwx)))
            hess_err = max(np.max(abs(part[:, 10] - hxx)),
                           np.max(abs(part[:, 11] - hxy)),
                           np.max(abs(part[:, 12] - hxy * ((part[:, 2]-1) % 5))),
                           np.max(abs(part[:, 13] - hxy * ((part[:, 2]-1) % 5) *
                                      ((part[:, 2]-1) // 5))),
                           np.max(abs(part[:, 16] - hwxx)),
                           np.max(abs(part[:, 17] - hwxy)))
            max_value_error = max(max_value_error, value_err)
            max_hessian_error = max(max_hessian_error, hess_err)
            assert value_err < 1e-4, (mode, degree, 'value', value_err)
            assert grad_err < 1e-3, (mode, degree, 'gradient', grad_err)
            assert hess_err < 3e-3, (mode, degree, 'hessian', hess_err)
            count += 20 * 11
    print(f'PASS: {count} spline affine/displacement SciPy comparisons; '
          f'max value error {max_value_error:.3g}, Hessian error {max_hessian_error:.3g}')


if __name__ == '__main__':
    main()
