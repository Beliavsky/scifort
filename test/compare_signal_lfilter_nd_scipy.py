#!/usr/bin/env python3
"""Compare rank-5/15 assumed-rank lfilter outputs and final states with SciPy."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.signal import lfilter

ROOT = Path(__file__).resolve().parents[1]
SHAPES = {
    5: (4, 3, 2, 2, 6),
    15: (5, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7),
}


def fortran_decl(name: str, rank: int, kind: str) -> str:
    return f'    {kind}, allocatable :: {name}({",".join(":" for _ in range(rank))})'


def branch(rank: int, b_complex: bool, a_complex: bool, x_complex: bool) -> str:
    bsuff = 'b' if b_complex else 'real(b,dp)'
    asuff = 'a' if a_complex else 'real(a,dp)'
    suffix = f'{rank}{"c" if x_complex else "r"}'
    out_suffix = f'{rank}{"c" if (b_complex or a_complex or x_complex) else "r"}'
    dims = ','.join(f'dims({i})' for i in range(1, rank + 1))
    state_dims = ','.join(f'shape_state({i})' for i in range(1, rank + 1))
    complex_out = b_complex or a_complex or x_complex
    lines = [
        f'allocate(x{suffix}({dims}), y{out_suffix}({dims}))',
        f'allocate(zi{out_suffix}({state_dims}), zf{out_suffix}({state_dims}))',
        f'x{suffix} = reshape({"xdata" if x_complex else "real(xdata,dp)"},dims(:{rank}))',
        f'zi{out_suffix} = reshape({"zdata" if complex_out else "real(zdata,dp)"},shape_state(:{rank}))',
        f'if (with_zi == 1) then',
        f'    call lfilter_nd({bsuff},{asuff},x{suffix},y{out_suffix}, &',
        f'        zi=zi{out_suffix},zf=zf{out_suffix},axis=axis,status=status)',
        'else',
        f'    call lfilter_nd({bsuff},{asuff},x{suffix},y{out_suffix}, &',
        f'        zf=zf{out_suffix},axis=axis,status=status)',
        'end if',
        f'ydata = reshape(y{out_suffix},[size(y{out_suffix})])' if complex_out
        else f'ydata = cmplx(reshape(y{out_suffix},[size(y{out_suffix})]),0.0_dp,dp)',
        f'zdata = reshape(zf{out_suffix},[size(zf{out_suffix})])' if complex_out
        else f'zdata = cmplx(reshape(zf{out_suffix},[size(zf{out_suffix})]),0.0_dp,dp)',
    ]
    return '\n'.join('        '+line for line in lines)


def driver_source() -> str:
    declarations = []
    for rank in (5, 15):
        for suffix, kind in ((f'{rank}r','real(dp)'),(f'{rank}c','complex(dp)')):
            for prefix in ('x','y','zi','zf'):
                declarations.append(fortran_decl(prefix + suffix,rank,kind))
    case_parts = []
    for rank in (5,15):
        for bits in range(8):
            lines = branch(rank,bool(bits & 1),bool(bits & 2),bool(bits & 4))
            case_parts.append(f'    case ({rank*10+bits})\n' + lines)
    cases='\n'.join(case_parts)
    return f'''program lfilter_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter_nd
    implicit none
    {chr(10).join(declarations)}
    complex(dp), allocatable :: xdata(:), ydata(:), zdata(:), b(:), a(:)
    integer :: rank_value, axis, complex_flag, with_zi, nb, na, n, ns, code, status, i
    integer :: dims(15), shape_state(15)
    real(dp) :: re, im

    read(*,*) rank_value, axis, complex_flag, with_zi, nb, na
    read(*,*) dims
    n = product(dims(:rank_value))
    ns = max(nb,na)-1
    shape_state = dims
    shape_state(axis) = ns
    allocate(b(nb),a(na),xdata(n),zdata(product(shape_state(:rank_value))))
    allocate(ydata(n))
    do i = 1, nb
        read(*,*) re, im
        b(i) = cmplx(re,im,dp)
    end do
    do i = 1, na
        read(*,*) re, im
        a(i) = cmplx(re,im,dp)
    end do
    do i = 1, n
        read(*,*) re, im
        xdata(i) = cmplx(re,im,dp)
    end do
    do i = 1, size(zdata)
        read(*,*) re, im
        zdata(i) = cmplx(re,im,dp)
    end do
    select case (rank_value*10+complex_flag)
{cases}
    case default
        error stop 'unsupported rank'
    end select
    print '(i0)', status
    do i = 1, size(ydata)
        print '(2(es26.17,1x))', real(ydata(i),dp), aimag(ydata(i))
    end do
    do i = 1, size(zdata)
        print '(2(es26.17,1x))', real(zdata(i),dp), aimag(zdata(i))
    end do
end program lfilter_nd_driver
'''


def main() -> None:
    rng = np.random.default_rng(20261018)
    maximum = 0.0
    count = 0
    with tempfile.TemporaryDirectory(prefix='scifort-lfilter-nd-') as tmp:
        work = Path(tmp)
        (work/'driver.f90').write_text(driver_source())
        exe = work/'driver'
        subprocess.run([os.environ.get('FC','gfortran'),'-std=f2018','-O2',
            '-Werror=implicit-interface','-fcheck=all','-ffree-line-length-none',
            str(ROOT/'src/kinds.f90'),str(ROOT/'src/signal_lfilter.f90'),
            str(work/'driver.f90'),
            '-o',str(exe)],cwd=work,check=True)
        for rank, samples in ((5,180),(15,120)):
            shape = SHAPES[rank]
            for case in range(samples):
                axis = case % rank
                cflag = (case // rank) % 8
                with_zi = case % 3 != 0
                nb = int(rng.integers(1,5))
                na = int(rng.integers(1,5))
                nstate = max(nb,na)-1
                stshape = list(shape)
                stshape[axis] = nstate
                b = rng.normal(size=nb) * 0.15
                a = rng.normal(size=na) * 0.12
                a[0] = float(rng.uniform(0.8,2.0))
                x = rng.normal(size=shape)
                zi = rng.normal(size=stshape)*0.1
                if cflag & 1:
                    b = b + 1j*rng.normal(size=nb)*0.05
                if cflag & 2:
                    a = a + 1j*rng.normal(size=na)*0.05
                if cflag & 4:
                    x = x + 1j*rng.normal(size=shape)
                if cflag:
                    zi = zi + 1j*rng.normal(size=stshape)*0.1
                # The reference always supplies the same initial states; when
                # SciFort starts from zero, explicitly pass zero states to SciPy.
                actual_zi = zi if with_zi else np.zeros_like(zi)
                expected, expected_zf = lfilter(b,a,x,axis=axis,zi=actual_zi)
                parts = [f'{rank} {axis+1} {cflag} {int(with_zi)} {nb} {na}',
                         ' '.join(map(str,shape+ (1,)*(15-rank)))]
                for arr in (b,a,x.ravel(order='F'),zi.ravel(order='F')):
                    parts.extend(f'{z.real:.17g} {z.imag:.17g}' for z in arr.astype(complex))
                run = subprocess.run([str(exe)],cwd=work,input='\n'.join(parts)+'\n',
                    text=True,capture_output=True,check=True)
                lines=run.stdout.strip().splitlines()
                if int(lines[0]) != 0:
                    raise AssertionError(f'nonzero status in rank {rank} case {case}')
                vals=np.fromstring(' '.join(lines[1:]),sep=' ').reshape(-1,2)
                values=vals[:,0]+1j*vals[:,1]
                y=values[:x.size].reshape(shape,order='F')
                zf=values[x.size:].reshape(stshape,order='F')
                error=max(float(np.max(np.abs(y-expected))),
                          float(np.max(np.abs(zf-expected_zf))) if zf.size else 0.)
                maximum=max(maximum,error)
                np.testing.assert_allclose(y,expected,rtol=1e-12,atol=2e-12,
                    err_msg=f'output rank={rank}, case={case}')
                np.testing.assert_allclose(zf,expected_zf,rtol=1e-12,atol=2e-12,
                    err_msg=f'state rank={rank}, case={case}')
                count+=1
    print(f'PASS: {count} SciPy rank-5/rank-15 lfilter cases, max abs={maximum:.3e}')


if __name__=='__main__':
    main()
