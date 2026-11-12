#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare rank-2/5/15 SOS filtering and final states against SciPy 1.17."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import sosfilt

ROOT = Path(__file__).resolve().parents[1]
SHAPES = {
    2: (8, 9),
    5: (4, 3, 5, 2, 6),
    15: (3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8),
}


def declarations(rank: int, complex_data: bool) -> str:
    suffix = 'c' if complex_data else 'r'
    dtype = 'complex(dp)' if complex_data else 'real(dp)'
    rstate = 1 if rank == 15 else rank + 1
    return '\n'.join(f'    {dtype}, allocatable :: {name}{rank}{suffix}({dims})'
                     for name, dims in (
                         ('x', ','.join(':' for _ in range(rank))),
                         ('y', ','.join(':' for _ in range(rank))),
                         ('zi', ','.join(':' for _ in range(rstate))),
                         ('zf', ','.join(':' for _ in range(rstate))),
                     ))


def branch(rank: int, c_sos: bool, c_x: bool) -> str:
    dtype_x = 'c' if c_x else 'r'
    dtype_y = 'c' if (c_sos or c_x) else 'r'
    x = f'x{rank}{dtype_x}'
    y = f'y{rank}{dtype_y}'
    zi = f'zi{rank}{dtype_y}'
    zf = f'zf{rank}{dtype_y}'
    input_shape = ','.join(f'dims({i})' for i in range(1, rank + 1))
    state_shape = 'nstates' if rank == 15 else ','.join(f'sh({i})' for i in range(1, rank + 2))
    coeffs = 'sosc' if c_sos else 'sosr'
    cast_x = 'xdata' if c_x else 'real(xdata,dp)'
    cast_zi = 'stdata' if dtype_y == 'c' else 'real(stdata,dp)'
    lines = [
        f'allocate({x}({input_shape}),{y}({input_shape}))',
        f'allocate({zi}({state_shape}),{zf}({state_shape}))',
        f'{x} = reshape({cast_x},dims(:{rank}))',
        f'{zi} = reshape({cast_zi},shape({zi}))',
        'if (with_zi == 1) then',
        f'    call sosfilt_nd({coeffs},{x},{y},zi={zi},zf={zf},axis=axis,status=status)',
        'else',
        f'    call sosfilt_nd({coeffs},{x},{y},zf={zf},axis=axis,status=status)',
        'end if',
        f'ydata = reshape({y},[size({y})])' if dtype_y == 'c'
        else f'ydata = cmplx(reshape({y},[size({y})]),0.0_dp,dp)',
        f'stdata = reshape({zf},[size({zf})])' if dtype_y == 'c'
        else f'stdata = cmplx(reshape({zf},[size({zf})]),0.0_dp,dp)',
    ]
    return '\n'.join('        '+line for line in lines)


def driver_source() -> str:
    decls = '\n'.join(declarations(rank, iscomplex)
                      for rank in SHAPES for iscomplex in (False, True))
    branches = []
    for rank in SHAPES:
        for combination in range(4):
            branches.append(f'    case ({rank*10+combination})\n' +
                            branch(rank, bool(combination & 1), bool(combination & 2)))
    cases = '\n'.join(branches)
    return f'''! Original SciFort numerical reference driver, no SciPy source copied.
program sos_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_sos, only : sosfilt_nd
    implicit none
{decls}
    complex(dp), allocatable :: sosc(:,:), xdata(:), ydata(:), stdata(:)
    real(dp), allocatable :: sosr(:,:)
    integer :: r, axis, combination, with_zi, nsec, dims(15), sh(16)
    integer :: nstates, status, i, j, k
    real(dp) :: re, im

    read(*,*) r, axis, combination, with_zi, nsec
    read(*,*) dims
    sh = [nsec,dims]
    sh(axis+1) = 2
    nstates = product(sh(:r+1))
    allocate(sosc(nsec,6),sosr(nsec,6))
    allocate(xdata(product(dims(:r))),ydata(product(dims(:r))),stdata(nstates))
    do j = 1, 6
        do i = 1, nsec
            read(*,*) re, im
            sosc(i,j) = cmplx(re,im,dp)
        end do
    end do
    sosr = real(sosc,dp)
    do i = 1, size(xdata)
        read(*,*) re, im
        xdata(i) = cmplx(re,im,dp)
    end do
    do i = 1, nstates
        read(*,*) re, im
        stdata(i) = cmplx(re,im,dp)
    end do
    select case (r*10+combination)
{cases}
    case default
        error stop 'unsupported rank or type combination'
    end select
    print '(i0)', status
    do i = 1, size(ydata)
        print '(2(es26.17,1x))', real(ydata(i),dp), aimag(ydata(i))
    end do
    do i = 1, size(stdata)
        print '(2(es26.17,1x))', real(stdata(i),dp), aimag(stdata(i))
    end do
end program sos_nd_driver
'''


def main() -> None:
    rng = np.random.default_rng(20261019)
    tested = 0
    worst = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort-sos-nd-') as tmp:
        wd = Path(tmp)
        driver = wd/'driver.f90'
        driver.write_text(driver_source())
        exe = wd/'driver'
        subprocess.run([os.environ.get('FC','gfortran'), '-std=f2018', '-O2',
                        '-fcheck=all', '-Werror=implicit-interface',
                        '-ffree-line-length-none', str(ROOT/'src/kinds.f90'),
                        str(ROOT/'src/signal_lfilter.f90'),
                        str(ROOT/'src/signal_sos.f90'), str(driver), '-o',str(exe)],
                       cwd=wd,check=True)
        for rank, trials in ((2,100),(5,120),(15,100)):
            dims = SHAPES[rank]
            for case in range(trials):
                axis = case % rank
                combination = (case//rank)%4
                with_zi = (case%3 != 0)
                nsec = int(rng.integers(1,5))
                sos = np.zeros((nsec,6))
                sos[:,0] = rng.uniform(0.3,1.0,size=nsec)
                sos[:,1:3] = rng.normal(size=(nsec,2))*0.11
                sos[:,3] = 1.0
                sos[:,4:6] = rng.normal(size=(nsec,2))*0.12
                if combination & 1:
                    sos = sos.astype(np.complex128)
                    sos[:,0:3] += 1j*rng.normal(size=(nsec,3))*0.04
                    sos[:,4:6] += 1j*rng.normal(size=(nsec,2))*0.03
                x = rng.normal(size=dims)
                if combination & 2:
                    x = x + 1j*rng.normal(size=dims)*0.4
                sh = (nsec, *dims[:axis], 2, *dims[axis+1:])
                zi = rng.normal(size=sh)*0.09
                if combination:
                    zi = zi + 1j*rng.normal(size=sh)*0.08
                actual_zi = zi if with_zi else np.zeros_like(zi)
                expected, expected_zf = sosfilt(sos,x,zi=actual_zi,axis=axis)
                head = [f'{rank} {axis+1} {combination} {int(with_zi)} {nsec}',
                        ' '.join(map(str,dims+(1,)*(15-rank)))]
                for arr in (sos.ravel(order='F'),x.ravel(order='F'),zi.ravel(order='F')):
                    head.extend(f'{v.real:.17g} {v.imag:.17g}' for v in arr.astype(complex))
                run = subprocess.run([str(exe)],cwd=wd,input='\n'.join(head)+'\n',
                                     text=True,capture_output=True,check=True)
                lines = run.stdout.strip().splitlines()
                if int(lines[0]) != 0:
                    raise AssertionError(f'failed status for rank {rank}, case {case}')
                values = np.fromstring(' '.join(lines[1:]),sep=' ').reshape(-1,2)
                values = values[:,0]+1j*values[:,1]
                y = values[:x.size].reshape(dims,order='F')
                zf = values[x.size:].reshape(sh,order='F')
                err = max(float(np.max(np.abs(y-expected))),
                          float(np.max(np.abs(zf-expected_zf))))
                worst = max(worst,err)
                np.testing.assert_allclose(y,expected,rtol=2.e-13,atol=2.e-13,
                    err_msg=f'output rank={rank} case={case}')
                np.testing.assert_allclose(zf,expected_zf,rtol=2.e-13,atol=2.e-13,
                    err_msg=f'state rank={rank} case={case}')
                tested += 1
    print(f'PASS: {tested} SciPy SOS ND cases, max abs={worst:.3e}')


if __name__ == '__main__':
    main()
