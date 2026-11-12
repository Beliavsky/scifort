#!/usr/bin/env python3
"""Compare automatic step/impulse response grids and values with SciPy 1.17.

The SciFort routines support real and complex state-space, TF, and ZPK
systems. Real SISO transfer functions are also compared with SciPy's step
and impulse. For complex and MIMO systems, use independent SciPy matrix
exponentials to compare each input's response.
"""
from __future__ import annotations

from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

import numpy as np
from scipy import linalg, signal

from run_signal_suite import MODULES, ROOT

PROBE = r'''
program continuous_auto_probe
    use scifort_signal, only : dp, step, impulse, step_zpk, impulse_zpk, &
        continuous_response_result
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), num(:), den(:), z(:), p(:)
    complex(dp) :: gain
    type(continuous_response_result) :: out
    integer :: total, idx, mode, impulse_flag, ns, ni, no, nb, na, nz, np, n, i, j
    read(*,*) total
    do idx=1,total
        read(*,*) mode,impulse_flag,ns,ni,no,nb,na,nz,np,n
        if (mode==0) then
            allocate(a(ns,ns),b(ns,ni),c(no,ns),d(no,ni))
            do i=1,ns
                do j=1,ns
                    call read_complex(a(i,j))
                end do
            end do
            do i=1,ns
                do j=1,ni
                    call read_complex(b(i,j))
                end do
            end do
            do i=1,no
                do j=1,ns
                    call read_complex(c(i,j))
                end do
            end do
            do i=1,no
                do j=1,ni
                    call read_complex(d(i,j))
                end do
            end do
        else if (mode==1) then
            allocate(num(nb),den(na))
            do i=1,nb
                call read_complex(num(i))
            end do
            do i=1,na
                call read_complex(den(i))
            end do
        else
            allocate(z(nz),p(np))
            do i=1,nz
                call read_complex(z(i))
            end do
            do i=1,np
                call read_complex(p(i))
            end do
            call read_complex(gain)
        end if
        select case(mode)
        case(0)
            if (impulse_flag==0) then
                out=step(a,b,c,d,n=n)
            else
                out=impulse(a,b,c,d,n=n)
            end if
        case(1)
            if (impulse_flag==0) then
                out=step(num,den,n=n)
            else
                out=impulse(num,den,n=n)
            end if
        case(2)
            if (impulse_flag==0) then
                out=step_zpk(z,p,gain,n=n)
            else
                out=impulse_zpk(z,p,gain,n=n)
            end if
        end select
        if (out%status/=0) then
            write(*,'(4(i0,1x))') out%status,0,0,0
        else
            write(*,'(4(i0,1x))') out%status,size(out%t),size(out%y,2),size(out%y,3)
            do i=1,size(out%t)
                write(*,'(es25.16e3)') out%t(i)
                do j=1,size(out%y,3)
                    call print_row(out%y(i,:,j))
                end do
            end do
        end if
        if (allocated(a)) deallocate(a,b,c,d)
        if (allocated(num)) deallocate(num,den)
        if (allocated(z)) deallocate(z,p)
    end do
contains
    subroutine read_complex(v)
        complex(dp), intent(out) :: v !! complex input
        real(dp) :: re, im
        read(*,*) re,im
        v=cmplx(re,im,dp)
    end subroutine read_complex
    subroutine print_row(row)
        complex(dp), intent(in) :: row(:) !! output channels
        integer :: kk
        do kk=1,size(row)
            write(*,'(2(es25.16e3,1x))') real(row(kk),dp),aimag(row(kk))
        end do
    end subroutine print_row
end program continuous_auto_probe
'''


def serialize_complex(buf: list[str], values: np.ndarray | complex):
    for value in np.ravel(np.asarray(values)):
        buf.append(f"{value.real:.17g} {value.imag:.17g}\n")


def random_case(rng, mode, i):
    n = int(rng.integers(1, 5))
    count = int(rng.integers(2, 33))
    do_impulse = i % 2 == 0
    is_complex = i % 3 == 0
    poles = -rng.uniform(0.6, 3.4, size=n)
    if is_complex:
        poles = poles + 1j * rng.uniform(-1.5, 1.5, size=n)
    if mode == 0:
        nin = int(rng.integers(1, 4))
        nout = int(rng.integers(1, 4))
        q = rng.normal(size=(n, n))
        if is_complex:
            q = q + 0.25j * rng.normal(size=(n, n))
        a = q @ np.diag(poles) @ np.linalg.inv(q)
        b = rng.normal(size=(n, nin))
        c = rng.normal(size=(nout, n))
        d = rng.normal(size=(nout, nin)) * 0.2
        if is_complex:
            b = b + 0.1j * rng.normal(size=b.shape)
            c = c + 0.1j * rng.normal(size=c.shape)
            d = d + 0.1j * rng.normal(size=d.shape)
        return dict(mode=mode, impulse=int(do_impulse), n=n, nin=nin, nout=nout,
                    count=count, a=a, b=b, c=c, d=d, poles=poles)
    zeros = -rng.uniform(0.3, 2.5, size=int(rng.integers(0, n + 1)))
    if is_complex:
        zeros = zeros + 1j * rng.uniform(-1.5, 1.5, size=zeros.size)
    gain = rng.uniform(0.2, 2.5)
    if is_complex:
        gain = gain + 0.25j
    if mode == 1:
        num, den = signal.zpk2tf(zeros, poles, gain)
        a, b, c, d = signal.tf2ss(num, den)
        return dict(mode=mode, impulse=int(do_impulse), n=n, nin=1, nout=1,
                    count=count, a=a, b=b, c=c, d=d, num=num, den=den)
    a, b, c, d = signal.zpk2ss(zeros, poles, gain)
    return dict(mode=mode, impulse=int(do_impulse), n=n, nin=1, nout=1,
                count=count, a=a, b=b, c=c, d=d, zeros=zeros, poles=poles, gain=gain)


def expected(case):
    a, b, c, d = (np.asarray(case[key], dtype=complex) for key in ('a', 'b', 'c', 'd'))
    eig = np.linalg.eigvals(a)
    r = float(np.min(np.abs(eig.real))) if eig.size else 1.0
    if r == 0.0:
        r = 1.0
    t = np.linspace(0.0, 7.0 / r, case['count'])
    response = np.empty((t.size, c.shape[0], b.shape[1]), dtype=complex)
    if case['impulse']:
        for i, ti in enumerate(t):
            response[i] = c @ linalg.expm(a * ti) @ b
    else:
        ns, nin = b.shape
        aug = np.zeros((ns + nin, ns + nin), dtype=complex)
        aug[:ns, :ns] = a
        aug[:ns, ns:] = b
        for i, ti in enumerate(t):
            response[i] = c @ linalg.expm(aug * ti)[:ns, ns:] + d
    return t, response


def main():
    rng = np.random.default_rng(20261017)
    cases = [random_case(rng, i % 3, i) for i in range(180)]
    # Cover the zero-root convention, repeated roots and static TF models.
    for den in ([1.0, 0.0], [1.0, 1.0, 0.0], [1.0], [1.0, 2.0, 1.0]):
        num = np.asarray([0.8])
        a, b, c, d = signal.tf2ss(num, den)
        cases.append(dict(mode=1,impulse=0,n=len(den)-1,nin=1,nout=1,
                          count=27,a=a,b=b,c=c,d=d,num=num,den=np.asarray(den)))
    buf = [str(len(cases)) + '\n']
    for case in cases:
        mode = case['mode']
        shape = f"{mode} {case['impulse']} {case['n']} {case['nin']} {case['nout']} "
        shape += f"{len(case.get('num',[]))} {len(case.get('den',[]))} "
        shape += f"{len(case.get('zeros',[]))} {len(case.get('poles',[]))} {case['count']}\n"
        buf.append(shape)
        if mode == 0:
            for key in ('a', 'b', 'c', 'd'):
                serialize_complex(buf,case[key])
        elif mode == 1:
            serialize_complex(buf,case['num'])
            serialize_complex(buf,case['den'])
        else:
            serialize_complex(buf,case['zeros'])
            serialize_complex(buf,case['poles'])
            serialize_complex(buf,case['gain'])
    with TemporaryDirectory(prefix='scifort-auto-time-') as temp:
        work = Path(temp)
        objects=[]
        for module in MODULES:
            src = ROOT / 'src' / (module+'.f90')
            obj = work / (module+'.o')
            subprocess.run(['gfortran','-std=f2018','-Wimplicit-interface',
                            '-Werror=implicit-interface','-O1','-J',str(work),
                            '-I',str(work),'-c',str(src),'-o',str(obj)],
                           check=True,cwd=work,capture_output=True)
            objects.append(str(obj))
        source = work / 'auto_probe.f90'
        source.write_text(PROBE)
        executable = work / 'auto_probe'
        subprocess.run(['gfortran','-std=f2018','-O1','-I',str(work),
                        str(source),*objects,'-o',str(executable)],
                       check=True,cwd=work,capture_output=True)
        proc = subprocess.run([str(executable)],input=''.join(buf),text=True,
                              capture_output=True,check=True,cwd=work)
        raw = iter(proc.stdout.splitlines())
        max_err = 0.0
        for index, case in enumerate(cases):
            status, count, nout, nin = map(int,next(raw).split())
            if status != 0:
                raise AssertionError(f'case {index}: SciFort returned status {status}')
            if (count,nout,nin) != (case['count'],case['nout'],case['nin']):
                raise AssertionError(f'case {index}: incorrect response dimensions')
            t = np.empty(count)
            response = np.empty((count,nout,nin),dtype=complex)
            for i in range(count):
                t[i] = float(next(raw))
                for k in range(nin):
                    for j in range(nout):
                        re, im = map(float,next(raw).split())
                        response[i,j,k] = re+1j*im
            expected_t, expected_y = expected(case)
            # Time scales are sensitive to small eigenvalue errors in
            # ill-conditioned TF factorizations; compare relative values.
            time_error = np.max(np.abs(t-expected_t))/max(1.0,expected_t[-1])
            response_error = np.max(np.abs(response-expected_y)) / max(1.0,np.max(np.abs(expected_y)))
            if time_error > 2e-6 or response_error > 2e-6:
                raise AssertionError(f'case {index}: grid {time_error:.3e}, output {response_error:.3e}')
            max_err = max(max_err,time_error,response_error)
            if case['mode']==1 and not np.iscomplexobj(case['den']) and not np.iscomplexobj(case['num']):
                f = signal.impulse if case['impulse'] else signal.step
                t2,y2 = f((case['num'],case['den']),N=count)
                if np.max(np.abs(t-t2)) > 1e-7 or np.max(np.abs(response[:,0,0]-y2)) > 2e-6:
                    raise AssertionError(f'case {index}: direct SciPy signal.{f.__name__} mismatch')
    print(f'PASS: {len(cases)} automatic continuous-time cases, maximum scaled error {max_err:.3e}')


if __name__ == '__main__':
    main()
