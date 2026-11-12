#!/usr/bin/env python3
"""Check assumed-rank Savitzky-Golay filtering at ranks 5 and 15 against SciPy."""
from __future__ import annotations

from pathlib import Path
import os
import subprocess
import tempfile

import numpy as np
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[1]
DRIVER = r'''
program savgol_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_filter_nd
    implicit none
    real(dp) :: x5(5,4,3,2,6), y5(5,4,3,2,6)
    real(dp) :: x15(5,2,1,1,3,1,1,1,1,1,1,1,1,1,7)
    real(dp) :: y15(5,2,1,1,3,1,1,1,1,1,1,1,1,1,7)
    complex(dp) :: c5(5,4,3,2,6), d5(5,4,3,2,6)
    complex(dp) :: c15(5,2,1,1,3,1,1,1,1,1,1,1,1,1,7)
    complex(dp) :: d15(5,2,1,1,3,1,1,1,1,1,1,1,1,1,7)
    complex(dp), allocatable :: data(:), output(:)
    real(dp) :: delta, cv_r, cv_i, xr, xi
    integer :: ndim, axis, window, degree, deriv, complex_flag, code, i, n
    character(len=12) :: mode

    read(*,*) ndim, axis, window, degree, deriv, complex_flag
    read(*,*) delta, mode, cv_r, cv_i
    if (ndim == 5) then
        n = size(x5)
    else
        n = size(x15)
    end if
    allocate(data(n), output(n))
    do i = 1, n
        read(*,*) xr, xi
        data(i) = cmplx(xr,xi,dp)
    end do
    if (ndim == 5) then
        if (complex_flag == 0) then
            x5 = reshape(real(data,dp),shape(x5))
            call savgol_filter_nd(x5,y5,window,degree,deriv=deriv,delta=delta, &
                mode=trim(mode),cval=cv_r,axis=axis,status=code)
            output = cmplx(reshape(y5,[n]),0.0_dp,dp)
        else
            c5 = reshape(data,shape(c5))
            call savgol_filter_nd(c5,d5,window,degree,deriv=deriv,delta=delta, &
                mode=trim(mode),cval=cmplx(cv_r,cv_i,dp),axis=axis,status=code)
            output = reshape(d5,[n])
        end if
    else
        if (complex_flag == 0) then
            x15 = reshape(real(data,dp),shape(x15))
            call savgol_filter_nd(x15,y15,window,degree,deriv=deriv,delta=delta, &
                mode=trim(mode),cval=cv_r,axis=axis,status=code)
            output = cmplx(reshape(y15,[n]),0.0_dp,dp)
        else
            c15 = reshape(data,shape(c15))
            call savgol_filter_nd(c15,d15,window,degree,deriv=deriv,delta=delta, &
                mode=trim(mode),cval=cmplx(cv_r,cv_i,dp),axis=axis,status=code)
            output = reshape(d15,[n])
        end if
    end if
    print *, code
    do i = 1, n
        print '(2(es25.16,1x))',real(output(i),dp),aimag(output(i))
    end do
end program savgol_nd_driver
'''


def main() -> None:
    rng = np.random.default_rng(20261018)
    maximum = 0.0
    total = 0
    shape5 = (5,4,3,2,6)
    shape15 = (5,2,1,1,3,1,1,1,1,1,1,1,1,1,7)
    with tempfile.TemporaryDirectory(prefix='scifort-savgol-nd-') as directory:
        work = Path(directory)
        (work / 'driver.f90').write_text(DRIVER)
        exe = work / 'driver'
        subprocess.run([os.environ.get('FC','gfortran'), '-std=f2018', '-O2',
                        '-Werror=implicit-interface', str(ROOT/'src/kinds.f90'),
                        str(ROOT/'src/signal_savgol.f90'), str(work/'driver.f90'),
                        '-o',str(exe)], cwd=work, check=True)
        for ndim, cases in ((5,180),(15,120)):
            shape = shape5 if ndim == 5 else shape15
            for test in range(cases):
                axis = (test % ndim)
                modes = ('interp','mirror','nearest','wrap','constant')
                mode = modes[test % len(modes)]
                window = int(rng.integers(1,min(shape[axis],5)+1))
                degree = int(rng.integers(0,window))
                deriv = int(rng.integers(0,degree+3))
                spacing = float(rng.uniform(0.4,2.0))
                complex_flag = (test // ndim) % 2
                x = rng.normal(size=shape)
                if complex_flag:
                    x = x + 1j * rng.normal(size=shape)
                cval = complex(rng.normal(),rng.normal())
                args = [ndim,axis+1,window,degree,deriv,complex_flag]
                inp = ' '.join(map(str,args)) + '\n'
                inp += f'{spacing:.17g} {mode} {cval.real:.17g} {cval.imag:.17g}\n'
                inp += ''.join(f'{z.real:.17g} {z.imag:.17g}\n'
                               for z in x.ravel(order='F').astype(complex))
                run = subprocess.run([str(exe)], input=inp, text=True,
                                     capture_output=True, check=True)
                lines = run.stdout.strip().splitlines()
                if int(lines[0]) != 0:
                    raise AssertionError(f'Nonzero status in case {ndim}:{test}: {lines[0]}, args={args}, mode={mode}, delta={spacing}, cval={cval}')
                values = np.fromstring(' '.join(lines[1:]),sep=' ').reshape(-1,2)
                got = (values[:,0] + 1j*values[:,1]).reshape(shape,order='F')
                expected = savgol_filter(x.real,window,degree,deriv=deriv,
                                         delta=spacing,axis=axis,mode=mode,cval=cval.real)
                if complex_flag:
                    expected = expected + 1j*savgol_filter(x.imag,window,degree,
                        deriv=deriv,delta=spacing,axis=axis,mode=mode,cval=cval.imag)
                diff = np.max(np.abs(got-expected))
                maximum = max(maximum,float(diff))
                np.testing.assert_allclose(got,expected,atol=1e-10,rtol=1e-10,
                    err_msg=f'ndim={ndim}, test={test}, axis={axis}, mode={mode}, delta={spacing}, cval={cval}')
                total += 1
    print(f'PASS: {total} SciPy comparisons (rank 5 and 15), max abs={maximum:.3e}')


if __name__ == '__main__':
    main()
