#!/usr/bin/env python3
"""Randomized SciPy 1.17 comparison for SciFort TF/ZPK/digital SOS conversions.

Run from any directory with NumPy, SciPy and gfortran. Does not require FPM.
Checks roots/gain, polynomial coefficients, and digital frequency responses;
SOS row ordering may differ from SciPy's pole-zero pairing.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
import sys
import numpy as np
from scipy import signal
from scipy.optimize import linear_sum_assignment
from run_signal_suite import ROOT, MODULES


def generate_cases():
    rng = np.random.default_rng(173417)
    cases = []
    for n in range(10):
        for i in range(84):
            m = (i * 7 + 3) % 10
            if n == 0 and m == 0:
                z, p = np.array([],complex), np.array([],complex)
            elif i % 6 == 0:
                # independent complex coefficients for the transfer function
                z = rng.uniform(-.76,.76,m) + 1j*rng.uniform(-.6,.6,m)
                p = rng.uniform(-.67,.67,n) + 1j*rng.uniform(-.4,.4,n)
            else:
                def symmetric_roots(degree):
                    roots = []
                    for _ in range(degree // 2):
                        if i % 4 == 0:
                            roots.extend(rng.uniform(-.8,.8,2))
                        else:
                            x, y = rng.uniform(-.8,.8),rng.uniform(.03,.7)
                            roots.extend([x+1j*y,x-1j*y])
                    if degree % 2:
                        roots.append(rng.uniform(-.8,.8))
                    return np.array(roots,dtype=complex)
                z, p = symmetric_roots(m), symmetric_roots(n)
            k = complex(rng.uniform(.05,2.) * (-1 if i % 7 == 0 else 1))
            if i % 6 == 0:
                k += 1j*rng.uniform(-.3,.3)
            cases.append((z,p,k))
    return cases


def print_complex_lines(values):
    return [f'{v.real:.17e} {v.imag:.17e}' for v in values]


def fortran_input(cases):
    lines = [str(len(cases))]
    for z,p,k in cases:
        lines.append(f'{len(z)} {len(p)}')
        lines.extend(print_complex_lines([k]))
        lines.extend(print_complex_lines(z))
        lines.extend(print_complex_lines(p))
    return '\n'.join(lines)+'\n'


def probe_source():
    return r'''program check_zpk
    use scifort_signal, only : dp, zpk2tf, tf2zpk, zpk2sos, &
        filter_zpk_result, filter_tf_result, filter_sos_result
    implicit none
    complex(dp), allocatable :: z(:), p(:)
    complex(dp) :: k
    type(filter_tf_result) :: tf
    type(filter_zpk_result) :: zp
    type(filter_sos_result) :: ss
    integer :: count, i, j, nz, np
    real(dp) :: x, y
    read (*,*) count
    do i = 1, count
        read (*,*) nz, np
        allocate(z(nz),p(np))
        read (*,*) x, y
        k = cmplx(x,y,dp)
        do j = 1, nz
            read (*,*) x, y
            z(j) = cmplx(x,y,dp)
        end do
        do j = 1, np
            read (*,*) x, y
            p(j) = cmplx(x,y,dp)
        end do
        tf = zpk2tf(z,p,k)
        write (*,'(3(i0,1x))') tf%status, size(tf%b), size(tf%a)
        call print_vector(tf%b)
        call print_vector(tf%a)
        zp = tf2zpk(tf%b,tf%a)
        write (*,'(3(i0,1x))') zp%status, size(zp%zeros), size(zp%poles)
        write (*,'(2(es25.16e3,1x))') real(zp%gain,dp), aimag(zp%gain)
        call print_vector(zp%zeros)
        call print_vector(zp%poles)
        ss = zpk2sos(z,p,k)
        write (*,'(2(i0,1x))') ss%status, size(ss%sos,1)
        if (ss%status == 0) then
            do j = 1, size(ss%sos,1)
                call print_vector(ss%sos(j,:))
            end do
        end if
        deallocate(z,p)
    end do
contains
    subroutine print_vector(v)
        complex(dp), intent(in) :: v(:) !! values written as real/imag pairs
        integer :: j
        do j = 1, size(v)
            write (*,'(2(es25.16e3,1x))') real(v(j),dp), aimag(v(j))
        end do
    end subroutine print_vector
end program check_zpk
'''


def read_results(text, cases):
    lines = text.splitlines()
    pos = 0
    results = []

    def vector(n):
        nonlocal pos
        values = np.empty(n,dtype=complex)
        for i in range(n):
            r,im = map(float,lines[pos].split())
            values[i] = complex(r,im)
            pos += 1
        return values

    for z,p,k in cases:
        st, nb, na = map(int,lines[pos].split()); pos += 1
        b = vector(nb)
        a = vector(na)
        root_st, nz, npoles = map(int,lines[pos].split()); pos += 1
        gain = vector(1)[0]
        zz = vector(nz)
        pp = vector(npoles)
        sos_st, sections = map(int,lines[pos].split()); pos += 1
        sos = np.array([vector(6) for _ in range(sections)]) if sos_st == 0 else np.empty((0,6),complex)
        results.append((st,b,a,root_st,gain,zz,pp,sos_st,sos))
    if pos != len(lines):
        raise AssertionError(f'unparsed output lines: {len(lines)-pos}')
    return results


def compare_roots(actual, expected, label, i):
    if len(actual) != len(expected):
        raise AssertionError(f'case {i}: {label} count {len(actual)} != {len(expected)}')
    if len(actual):
        costs = abs(actual[:,None] - expected[None,:])
        ii,jj = linear_sum_assignment(costs)
        if max(costs[ii,jj]) > 1.5e-5:
            raise AssertionError(f'case {i}: {label} roots differ, max {max(costs[ii,jj])}')


def compare_case(case, output, idx):
    z,p,k = case
    status,b,a,root_st,gain,zz,pp,sos_st,sos = output
    if status or root_st or sos_st:
        raise AssertionError(f'case {idx}: statuses = {status}, {root_st}, {sos_st}')
    expected_b,expected_a = signal.zpk2tf(z,p,k)
    if not np.allclose(b,expected_b,rtol=2e-11,atol=3e-12):
        raise AssertionError(f'case {idx}: polynomial b differs')
    if not np.allclose(a,expected_a,rtol=2e-11,atol=3e-12):
        raise AssertionError(f'case {idx}: polynomial a differs')
    ref_z,ref_p,ref_k = signal.tf2zpk(expected_b,expected_a)
    if abs(ref_k-gain) > 1e-8*max(1.,abs(ref_k)):
        raise AssertionError(f'case {idx}: gain {gain} != {ref_k}')
    compare_roots(zz,ref_z,'zeros',idx)
    compare_roots(pp,ref_p,'poles',idx)
    if len(z)==len(p):
        w = np.linspace(0.025,np.pi-0.025,21)
        ours = signal.sosfreqz(sos,worN=w)[1]
        zzw = np.exp(1j*w)
        reference = k*np.prod(zzw[:,None]-z[None,:],axis=1) / np.prod(zzw[:,None]-p[None,:],axis=1)
        tol = 5e-10*(1.+np.maximum(np.abs(ours),np.abs(reference)))
        if np.any(np.abs(ours-reference)>tol):
            raise AssertionError(f'case {idx}: SOS response differs')
    if abs(k.imag) < 1e-14:
        ref_sos = signal.zpk2sos(z,p,k.real)
        w = np.linspace(0.035,2.96,21)
        v1 = signal.sosfreqz(sos,worN=w)[1]
        v2 = signal.sosfreqz(ref_sos,worN=w)[1]
        if not np.allclose(v1,v2,rtol=2e-10,atol=2e-10):
            raise AssertionError(f'case {idx}: SOS response not SciPy-compatible')
    return True


def main():
    cases = generate_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-zpk-') as d:
        work = Path(d)
        module_objects=[]
        flags=['-std=f2018','-O0','-Wimplicit-interface','-Werror=implicit-interface','-fcheck=all']
        for name in MODULES:
            source=ROOT/'src'/f'{name}.f90'
            obj=work/f'{name}.o'
            subprocess.run(['gfortran',*flags,'-J',str(work),'-I',str(work),'-c',str(source),'-o',str(obj)],cwd=work,check=True)
            module_objects.append(str(obj))
        source=work/'probe.f90'
        source.write_text(probe_source())
        exe=work/'probe'
        subprocess.run(['gfortran',*flags,'-J',str(work),'-I',str(work),str(source),*module_objects,'-o',str(exe)],cwd=work,check=True)
        result = subprocess.run([str(exe)],input=fortran_input(cases),text=True,capture_output=True,cwd=work,check=True)
        outputs=read_results(result.stdout,cases)
        for i,(c,o) in enumerate(zip(cases,outputs)):
            compare_case(c,o,i)
    print(f'PASS: {len(cases)} randomized TF/ZPK/SOS cases vs SciPy 1.17')

if __name__=='__main__':
    main()
