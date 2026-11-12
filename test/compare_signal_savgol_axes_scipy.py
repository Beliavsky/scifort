#!/usr/bin/env python3
"""Compare rank-2/rank-3 Savitzky-Golay filtering with SciPy 1.17.

Uses a tiny Fortran driver compiled against the existing SciFort modules.
Complex SciFort filtering is compared to separate SciPy real/imaginary filtering;
SciPy itself casts complex signals to real in this version.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
FORTRAN = os.environ.get('FC', 'gfortran')
DRIVER = r'''
program test_driver
    use scifort_kinds, only : dp
    use scifort_signal_savgol, only : savgol_filter
    implicit none
    integer :: ndim,n1,n2,n3,along,w,p,d,kind,st,i,j,k,n
    real(dp) :: delta,cr,ci,vr,vi
    real(dp), allocatable :: x1(:),x2(:,:),x3(:,:,:),y1(:),y2(:,:),y3(:,:,:)
    complex(dp), allocatable :: z1(:),z2(:,:),z3(:,:,:), &
        v1(:),v2(:,:),v3(:,:,:)
    character(len=12) :: mode
    read(*,*) ndim,n1,n2,n3,along,w,p,d,kind
    read(*,*) delta,mode,cr,ci
    select case (ndim)
    case (1)
        allocate(z1(n1))
        do i=1,n1
            read(*,*) vr,vi
            z1(i)=cmplx(vr,vi,dp)
        end do
        if (kind==0) then
            x1=real(z1,dp)
            y1=savgol_filter(x1,w,p,deriv=d,delta=delta,mode=trim(mode),cval=cr,status=st)
            v1=cmplx(y1,0.0_dp,dp)
        else
            v1=savgol_filter(z1,w,p,deriv=d,delta=delta,mode=trim(mode), &
                cval=cmplx(cr,ci,dp),status=st)
        end if
        print *,st
        do i=1,n1
            print '(2(es25.16,1x))',real(v1(i),dp),aimag(v1(i))
        end do
    case (2)
        allocate(z2(n1,n2))
        do j=1,n2
            do i=1,n1
                read(*,*) vr,vi
                z2(i,j)=cmplx(vr,vi,dp)
            end do
        end do
        if (kind==0) then
            x2=real(z2,dp)
            y2=savgol_filter(x2,w,p,deriv=d,delta=delta,mode=trim(mode), &
                cval=cr,status=st,axis=along)
            v2=cmplx(y2,0.0_dp,dp)
        else
            v2=savgol_filter(z2,w,p,deriv=d,delta=delta,mode=trim(mode), &
                cval=cmplx(cr,ci,dp),status=st,axis=along)
        end if
        print *,st
        do j=1,n2
            do i=1,n1
                print '(2(es25.16,1x))',real(v2(i,j),dp),aimag(v2(i,j))
            end do
        end do
    case (3)
        allocate(z3(n1,n2,n3))
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    read(*,*) vr,vi
                    z3(i,j,k)=cmplx(vr,vi,dp)
                end do
            end do
        end do
        if (kind==0) then
            x3=real(z3,dp)
            y3=savgol_filter(x3,w,p,deriv=d,delta=delta,mode=trim(mode), &
                cval=cr,status=st,axis=along)
            v3=cmplx(y3,0.0_dp,dp)
        else
            v3=savgol_filter(z3,w,p,deriv=d,delta=delta,mode=trim(mode), &
                cval=cmplx(cr,ci,dp),status=st,axis=along)
        end if
        print *,st
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    print '(2(es25.16,1x))',real(v3(i,j,k),dp),aimag(v3(i,j,k))
                end do
            end do
        end do
    end select
end program test_driver
'''


def main() -> None:
    rng = np.random.default_rng(20261018)
    errors = []
    checks = 0
    with tempfile.TemporaryDirectory(prefix='scifort_savgol_') as tmp:
        tmpdir = Path(tmp)
        driver = tmpdir / 'driver.f90'
        driver.write_text(DRIVER)
        exe = tmpdir / 'driver'
        flags = ['-std=f2018', '-O2', '-Werror=implicit-interface']
        subprocess.run([FORTRAN, *flags, str(ROOT / 'src/kinds.f90'),
                        str(ROOT / 'src/signal_savgol.f90'), str(driver),
                        '-o', str(exe)], cwd=tmpdir, check=True)
        for case in range(600):
            ndim = [1, 2, 3][case % 3]
            shape = tuple(int(v) for v in rng.integers(5, 18, size=ndim))
            axis = int(rng.integers(0, ndim))
            window = int(rng.integers(2, min(shape[axis], 9) + 1))
            degree = int(rng.integers(0, window))
            derivative = int(rng.integers(0, degree+3))
            delta = float(rng.uniform(0.2, 2.0))
            mode = ['interp','mirror','nearest','wrap','constant'][case % 5]
            value = complex(rng.normal(), rng.normal())
            kind = (case // 3) % 2
            x = rng.normal(size=shape)
            if kind:
                x = x + 1j*rng.normal(size=shape)
            cval = complex(rng.normal(), rng.normal()) if kind else complex(float(rng.normal()),0.0)
            settings = f'{ndim} {shape[0]} {shape[1] if ndim>=2 else 1} ' \
                f'{shape[2] if ndim>=3 else 1} {axis+1} {window} {degree} {derivative} {kind}\n'
            settings += f'{delta:.17g} {mode} {cval.real:.17g} {cval.imag:.17g}\n'
            samples = ''.join(f'{z.real:.17g} {z.imag:.17g}\n' for z in x.ravel(order='F').astype(complex))
            completed = subprocess.run([str(exe)], input=settings+samples, capture_output=True,
                                       text=True, check=True)
            lines = completed.stdout.strip().splitlines()
            status = int(lines[0])
            if status != 0:
                raise AssertionError(f'Unexpected Fortran status {status}: {case}')
            nums = np.array([list(map(float,line.split())) for line in lines[1:]])
            got = (nums[:,0]+1j*nums[:,1]).reshape(shape,order='F')
            if kind:
                expected = signal.savgol_filter(x.real,window,degree,deriv=derivative,
                                                delta=delta,axis=axis,mode=mode,cval=cval.real)
                expected = expected + 1j*signal.savgol_filter(x.imag,window,degree,
                                                deriv=derivative,delta=delta,axis=axis,
                                                mode=mode,cval=cval.imag)
            else:
                expected = signal.savgol_filter(x,window,degree,deriv=derivative,
                                                delta=delta,axis=axis,mode=mode,cval=cval.real)
            error = float(np.max(np.abs(got-expected)))
            errors.append(error)
            np.testing.assert_allclose(got, expected, atol=2e-10, rtol=3e-10,
                                       err_msg=f'case {case}, shape {shape}, mode {mode}, axis {axis}')
            checks += 1
    print(f'Passed {checks} SciPy Savitzky-Golay comparisons; max absolute error {max(errors):.3e}')

if __name__ == '__main__':
    main()
