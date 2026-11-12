#!/usr/bin/env python3
"""Validate rank-4 real/complex Savitzky-Golay filtering against SciPy.

SciFort's complex-input extension is checked against separate SciPy filters
for the real and imaginary components, as SciPy 1.17 discards imaginary data.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[1]
DRIVER = r'''
program savgol_rank4_driver
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_filter
    implicit none
    integer :: n1, n2, n3, n4, axis, window, degree, deriv, kind, status
    integer :: i, j, k, l
    real(dp) :: spacing, cv_r, cv_i, xr, xi
    real(dp), allocatable :: x(:,:,:,:), y(:,:,:,:)
    complex(dp), allocatable :: z(:,:,:,:), w(:,:,:,:)
    character(len=12) :: mode

    read(*,*) n1,n2,n3,n4,axis,window,degree,deriv,kind
    read(*,*) spacing, mode, cv_r, cv_i
    allocate(z(n1,n2,n3,n4))
    do l=1,n4
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    read(*,*) xr,xi
                    z(i,j,k,l)=cmplx(xr,xi,dp)
                end do
            end do
        end do
    end do
    if (kind==0) then
        x=real(z,dp)
        y=savgol_filter(x,window,degree,deriv=deriv,delta=spacing, &
            mode=trim(mode),cval=cv_r,axis=axis,status=status)
        w=cmplx(y,0.0_dp,dp)
    else
        w=savgol_filter(z,window,degree,deriv=deriv,delta=spacing, &
            mode=trim(mode),cval=cmplx(cv_r,cv_i,dp),axis=axis,status=status)
    end if
    print *,status
    do l=1,n4
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    print '(2(es25.16,1x))',real(w(i,j,k,l),dp),aimag(w(i,j,k,l))
                end do
            end do
        end do
    end do
end program savgol_rank4_driver
'''


def main() -> None:
    rng = np.random.default_rng(20261019)
    max_error = 0.0
    max_scaled = 0.0
    nc = 400
    with tempfile.TemporaryDirectory(prefix='scifort-savgol-rank4-') as directory:
        wd = Path(directory)
        (wd / 'driver.f90').write_text(DRIVER)
        exe = wd / 'driver'
        cmd = [os.environ.get('FC', 'gfortran'), '-std=f2018', '-O2',
               '-Werror=implicit-interface', str(ROOT / 'src/kinds.f90'),
               str(ROOT / 'src/signal_savgol.f90'), str(wd / 'driver.f90'),
               '-o', str(exe)]
        subprocess.run(cmd, cwd=wd, check=True)
        for case in range(nc):
            shape = tuple(int(n) for n in rng.integers(5, 10, size=4))
            axis = case % 4
            mode = ('interp', 'mirror', 'nearest', 'wrap', 'constant')[case % 5]
            window = int(rng.integers(2, min(shape[axis], 9) + 1))
            degree = int(rng.integers(0, window))
            deriv = int(rng.integers(0, degree + 3))
            delta = float(rng.uniform(0.3, 2.0))
            kind = (case // 4) % 2
            x = rng.normal(size=shape)
            if kind:
                x = x + 1j * rng.normal(size=shape)
            cval = complex(rng.normal(), rng.normal()) if kind else complex(float(rng.normal()))
            # Use normal Fortran list-directed integers without tuple punctuation.
            inp = ' '.join(str(n) for n in (*shape, axis+1, window, degree, deriv, kind)) + '\n'
            inp += f'{delta:.17g} {mode} {cval.real:.17g} {cval.imag:.17g}\n'
            inp += ''.join(f'{z.real:.17g} {z.imag:.17g}\n'
                           for z in x.ravel(order='F').astype(complex))
            run = subprocess.run([str(exe)], input=inp, text=True,
                                 capture_output=True, check=True)
            lines = run.stdout.strip().splitlines()
            if int(lines[0]) != 0:
                raise AssertionError(f'Unexpected nonzero Fortran status in case {case}: {lines[0]}')
            nums = np.fromstring(' '.join(lines[1:]), sep=' ').reshape(-1, 2)
            got = (nums[:, 0] + 1j * nums[:, 1]).reshape(shape, order='F')
            if kind:
                expected = savgol_filter(x.real, window, degree, deriv=deriv,
                                         delta=delta, axis=axis, mode=mode, cval=cval.real)
                expected = expected + 1j * savgol_filter(x.imag, window, degree,
                    deriv=deriv, delta=delta, axis=axis, mode=mode, cval=cval.imag)
            else:
                expected = savgol_filter(x, window, degree, deriv=deriv,
                                         delta=delta, axis=axis, mode=mode, cval=cval.real)
            difference = np.abs(got - expected)
            max_error = max(max_error, float(difference.max()))
            max_scaled = max(max_scaled,
                             float((difference / (1.0 + np.abs(expected))).max()))
            np.testing.assert_allclose(got, expected, atol=3e-10, rtol=3e-10,
                err_msg=f'case={case}, shape={shape}, axis={axis}, mode={mode}')
    print(f'PASS: {nc} rank-4 SciPy comparisons, '
          f'max abs={max_error:.3e}, max scaled={max_scaled:.3e}')


if __name__ == '__main__':
    main()
