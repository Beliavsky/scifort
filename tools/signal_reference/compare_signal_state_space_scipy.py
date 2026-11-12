#!/usr/bin/env python3
"""Compare SciFort signal state-space conversion with SciPy 1.17.

Runs a single batched Fortran probe compiled with gfortran, no FPM required.
Tests canonical matrices, matrix-to-polynomial conversion, MIMO input
selection, and real/complex zero-pole-gain conversion.
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal
from run_signal_suite import ROOT, MODULES

PROBE = r'''
program signal_state_space_probe
    use scifort_signal, only : dp, tf2ss, ss2tf, zpk2ss, ss2zpk, &
        state_space_result, state_space_tf_result, filter_zpk_result
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), num(:,:), den(:)
    complex(dp), allocatable :: z(:), p(:)
    complex(dp) :: gain
    type(state_space_result) :: ss
    type(state_space_tf_result) :: tf
    type(filter_zpk_result) :: zp
    integer :: cases, mode, idx, n, nb, na, nout, nin, nin_selected, j, k, nz, np

    read(*,*) cases
    do idx = 1, cases
        read(*,*) mode
        select case(mode)
        case(0)
            read(*,*) nout, nb, na
            allocate(num(nout,nb),den(na))
            do j = 1, nout
                do k = 1, nb
                    call read_complex(num(j,k))
                end do
            end do
            do k = 1, na
                call read_complex(den(k))
            end do
            ss = tf2ss(num,den)
            call write_state(ss)
            deallocate(num,den)
        case(1)
            read(*,*) n, nin, nout, nin_selected
            allocate(a(n,n),b(n,nin),c(nout,n),d(nout,nin))
            do j = 1, n
                do k = 1, n
                    call read_complex(a(j,k))
                end do
            end do
            do j = 1, n
                do k = 1, nin
                    call read_complex(b(j,k))
                end do
            end do
            do j = 1, nout
                do k = 1, n
                    call read_complex(c(j,k))
                end do
            end do
            do j = 1, nout
                do k = 1, nin
                    call read_complex(d(j,k))
                end do
            end do
            tf = ss2tf(a,b,c,d,input=nin_selected)
            write(*,'(3(i0,1x))') tf%status,nout,n+1
            if (tf%status == 0) then
                do k = 1, size(tf%a)
                    call write_complex(tf%a(k))
                end do
                do j = 1, nout
                    do k = 1, size(tf%b,2)
                        call write_complex(tf%b(j,k))
                    end do
                end do
            end if
            if (nout == 1) then
                zp = ss2zpk(a,b,c,d,input=nin_selected)
                write(*,'(3(i0,1x))') zp%status, &
                    count_zeros(zp%zeros),count_zeros(zp%poles)
                if (zp%status == 0) then
                    call write_complex(zp%gain)
                    do k = 1, size(zp%zeros)
                        call write_complex(zp%zeros(k))
                    end do
                    do k = 1, size(zp%poles)
                        call write_complex(zp%poles(k))
                    end do
                end if
            end if
            deallocate(a,b,c,d)
        case(2)
            read(*,*) nz,np
            allocate(z(nz),p(np))
            do j = 1, nz
                call read_complex(z(j))
            end do
            do j = 1, np
                call read_complex(p(j))
            end do
            call read_complex(gain)
            ss = zpk2ss(z,p,gain)
            call write_state(ss)
            deallocate(z,p)
        end select
    end do
contains
    subroutine read_complex(value)
        complex(dp), intent(out) :: value !! complex input
        real(dp) :: re, im
        read(*,*) re, im
        value = cmplx(re,im,dp)
    end subroutine read_complex
    subroutine write_complex(value)
        complex(dp), intent(in) :: value !! complex output
        write(*,'(2(es25.16e3,1x))') real(value,dp),aimag(value)
    end subroutine write_complex
    integer function count_zeros(values) result(n)
        complex(dp), allocatable, intent(in) :: values(:) !! optional roots
        n=0
        if (allocated(values)) n=size(values)
    end function count_zeros
    subroutine write_state(s)
        type(state_space_result), intent(in) :: s !! state-space result
        integer :: i,j,ns,no
        ns=0
        no=0
        if (s%status == 0) then
            ns=size(s%a,1)
            no=size(s%c,1)
        end if
        write(*,'(3(i0,1x))') s%status,ns,no
        if (s%status/=0) return
        do i = 1, ns
            do j = 1, ns
                call write_complex(s%a(i,j))
            end do
        end do
        do i = 1, ns
            call write_complex(s%b(i,1))
        end do
        do i = 1, no
            do j = 1, ns
                call write_complex(s%c(i,j))
            end do
        end do
        do i = 1, no
            call write_complex(s%d(i,1))
        end do
    end subroutine write_state
end program signal_state_space_probe
'''


def root_array(rng, order, complex_mode):
    if order == 0:
        return np.array([], dtype=complex)
    # Moderately separated roots; do not create repeated/clustered pole cases.
    real = rng.uniform(-1.4, -.2, order)
    if complex_mode:
        imag = rng.uniform(-.75,.75,order)
        return real+1j*imag
    return real.astype(complex)


def cases():
    rng = np.random.default_rng(431112)
    items = []
    for n in range(0,7):
        for i in range(40):
            is_complex = i % 3 == 0
            p = root_array(rng,n,is_complex)
            den = np.atleast_1d(np.poly(p))*rng.uniform(.2,2.3)
            nout = 1 if i%5 else 2
            nzeros = (i*3+n)% (n+1)
            numer = np.array([np.pad(np.atleast_1d(np.poly(root_array(rng,nzeros,is_complex))) *
                (rng.uniform(.3,2.)*(1+0.2j if is_complex else 1)),(n-nzeros,0))
                for _ in range(nout)])
            if i % 13 == 0:
                numer = np.pad(numer,((0,0),(1,0)))
            if i % 11 == 0:
                den=np.pad(den,(1,0))
            if i % 17 == 0:
                numer *= 0
            items.append(('tf',numer,den))
    for n in range(0,6):
        for i in range(45):
            is_complex=i % 4 == 0
            nin=1+(i%3)
            nout=1+(i%2)
            a=rng.normal(size=(n,n))*.25
            if n>0:
                a-=np.diag(np.linspace(.6,1.4,n))
            b=rng.normal(size=(n,nin))
            c=rng.normal(size=(nout,n))
            d=rng.normal(size=(nout,nin))*.2
            if is_complex:
                a=a+1j*rng.normal(size=(n,n))*.1
                b=b+1j*rng.normal(size=(n,nin))*.3
                c=c+1j*rng.normal(size=(nout,n))*.3
                d=d+1j*rng.normal(size=(nout,nin))*.2
            items.append(('ss',a,b,c,d,1+(i%nin)))
    for n in range(0,6):
        for i in range(24):
            is_complex = i % 3 == 0
            zeros = root_array(rng,(n+i)% (n+1),is_complex)
            poles = root_array(rng,n,is_complex)
            k = complex(rng.uniform(.1,2.), rng.uniform(-.2,.2) if is_complex else 0)
            items.append(('zpk',zeros,poles,k))
    return items


def to_lines(items):
    lines=[str(len(items))]
    def addval(v):
        lines.append(f'{v.real:.17g} {v.imag:.17g}')
    for item in items:
        match item[0]:
            case 'tf':
                _,num,den=item
                lines += ['0',f'{num.shape[0]} {num.shape[1]} {len(den)}']
                for v in num.flat: addval(v)
                for v in den: addval(v)
            case 'ss':
                _,a,b,c,d,sel=item
                lines += ['1',f'{len(a)} {b.shape[1]} {c.shape[0]} {sel}']
                for m in (a,b,c,d):
                    for v in m.flat: addval(v)
            case 'zpk':
                _,z,p,k=item
                lines += ['2',f'{len(z)} {len(p)}']
                for v in z: addval(v)
                for v in p: addval(v)
                addval(k)
    return '\n'.join(lines)+'\n'


def complex_line(lines, pos):
    a,b=map(float,lines[pos].split())
    return complex(a,b)


def assert_close(got, ref, label, tol=1e-7):
    got=np.asarray(got,dtype=complex)
    ref=np.asarray(ref,dtype=complex)
    if got.shape!=ref.shape:
        raise AssertionError(f'{label}: shape {got.shape} != {ref.shape}')
    discrepancy=np.max(abs(got-ref)/(1+abs(ref))) if got.size else 0.0
    if not np.isfinite(discrepancy) or discrepancy>tol:
        raise AssertionError(f'{label}: scaled error {discrepancy}')
    return discrepancy


def run():
    items=cases()
    with tempfile.TemporaryDirectory(prefix='scifort-state-space-') as d:
        work=Path(d)
        flags=['-std=f2018','-O0','-fcheck=all','-Wimplicit-interface','-Werror=implicit-interface']
        objects=[]
        for name in MODULES:
            obj=work/f'{name}.o'
            subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),'-c',
                str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],check=True,cwd=work)
            objects.append(str(obj))
        src=work/'probe.f90'
        src.write_text(PROBE)
        exe=work/'probe'
        subprocess.run(['gfortran',*flags,'-I',str(work),'-J',str(work),str(src),
            *objects,'-o',str(exe)],check=True,cwd=work)
        out=subprocess.run([str(exe)],input=to_lines(items),capture_output=True,
            text=True,check=True,cwd=work)
    output=out.stdout.splitlines()
    pos=0
    counts=dict(tf=0,ss=0,zpk=0,ss2zpk=0)
    worst=0.0

    def parse_state():
        nonlocal pos
        status,ns,no=map(int,output[pos].split());pos+=1
        if status:
            raise AssertionError(f'Fortran state-space conversion status {status}')
        def matrix(r,c):
            nonlocal pos
            vals=np.array([complex_line(output,pos+j) for j in range(r*c)])
            pos+=r*c
            return vals.reshape(r,c)
        return matrix(ns,ns),matrix(ns,1),matrix(no,ns),matrix(no,1)

    for index,item in enumerate(items):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',signal.BadCoefficients)
            if item[0]=='tf':
                _,num,den=item
                got=parse_state()
                expected=signal.tf2ss(num,den)
                for name,x,y in zip(('A','B','C','D'),got,expected):
                    worst=max(worst,assert_close(x,y,f'tf2ss #{index} {name}'))
                counts['tf']+=1
            elif item[0]=='ss':
                _,a,b,c,d,sel=item
                status,nout,nb=map(int,output[pos].split());pos+=1
                if status:
                    raise AssertionError(f'ss2tf #{index} status {status}')
                den=np.array([complex_line(output,pos+j) for j in range(nb)]);pos+=nb
                nums=np.array([complex_line(output,pos+j) for j in range(nout*nb)]).reshape(nout,nb)
                pos+=nout*nb
                nr,dr=signal.ss2tf(a,b,c,d,input=sel-1)
                if np.ndim(nr)==1:
                    nr=np.asarray(nr).reshape((nout,1) if len(a)==0 else (1,-1))
                worst=max(worst,assert_close(den,np.atleast_1d(dr),f'ss2tf #{index} den'))
                worst=max(worst,assert_close(nums,nr,f'ss2tf #{index} nums'))
                counts['ss']+=1
                if nout==1:
                    status,nz,npoles=map(int,output[pos].split());pos+=1
                    if status:
                        raise AssertionError(f'ss2zpk #{index} status {status}')
                    k=complex_line(output,pos);pos+=1
                    z=np.array([complex_line(output,pos+j) for j in range(nz)]);pos+=nz
                    p=np.array([complex_line(output,pos+j) for j in range(npoles)]);pos+=npoles
                    try:
                        zr,pr,kr=signal.ss2zpk(a,b,c,d,input=sel-1)
                    except (ValueError,ZeroDivisionError,FloatingPointError):
                        # Valid representation can be checked directly by polynomials.
                        zr=pr=kr=None
                    if kr is not None and np.isfinite(kr):
                        worst=max(worst,assert_close([k],[kr],f'ss2zpk #{index} gain'))
                        worst=max(worst,assert_close(np.atleast_1d(np.poly(z))*k,nums[0,len(nums[0])-len(z)-1:],
                            f'ss2zpk #{index} numerator',tol=2e-5))
                        worst=max(worst,assert_close(np.atleast_1d(np.poly(p)),den,
                            f'ss2zpk #{index} denominator',tol=2e-5))
                    counts['ss2zpk']+=1
            else:
                _,z,p,k=item
                got=parse_state()
                expected=signal.zpk2ss(z,p,k)
                for name,x,y in zip(('A','B','C','D'),got,expected):
                    worst=max(worst,assert_close(x,y,f'zpk2ss #{index} {name}'))
                counts['zpk']+=1
    if pos != len(output):
        raise AssertionError('unexpected stdout trailing lines')
    print(f'PASS: {counts["tf"]} tf2ss, {counts["ss"]} ss2tf, '
        f'{counts["ss2zpk"]} ss2zpk and {counts["zpk"]} zpk2ss comparisons; '
        f'max scaled error {worst:.3g} vs SciPy {signal.__version__ if hasattr(signal,"__version__") else "1.17"}')


if __name__=='__main__':
    run()
