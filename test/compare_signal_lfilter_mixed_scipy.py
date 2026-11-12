#!/usr/bin/env python3
"""Seeded cross-check of mixed real/complex lfilter inputs against SciPy 1.17.

Compiles a temporary Fortran driver; verifies output samples and final filter
states for 1D, rank-two, and rank-three arrays with all six mixed type triples.
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
COMBOS = ('rrc','rcr','rcc','crr','crc','ccr')


def make_driver() -> str:
    header = '''program mixed_driver
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter
    implicit none
    integer :: rank,n1,n2,n3,ax,nb,na,ns,code,i,j,k,itype
    integer :: nd(3), sd(3)
    real(dp) :: re,im
    complex(dp), allocatable :: bc(:),ac(:),xc1(:),xc2(:,:),xc3(:,:,:)
    complex(dp), allocatable :: zi1(:),zi2(:,:),zi3(:,:,:)
    complex(dp), allocatable :: zf1(:),zf2(:,:),zf3(:,:,:)
    complex(dp), allocatable :: y1(:),y2(:,:),y3(:,:,:)
    real(dp), allocatable :: br(:),ar(:),xr1(:),xr2(:,:),xr3(:,:,:)
    read(*,*) rank,n1,n2,n3,ax,nb,na,itype
    ns=max(nb,na)-1
    allocate(bc(nb),ac(na))
    do i=1,nb
        read(*,*) re,im
        bc(i)=cmplx(re,im,dp)
    end do
    do i=1,na
        read(*,*) re,im
        ac(i)=cmplx(re,im,dp)
    end do
    br=real(bc,dp)
    ar=real(ac,dp)
    nd=[n1,n2,n3]
    sd=nd
    sd(ax)=ns
    select case(rank)
    case(1)
        allocate(xc1(n1),zi1(ns))
        do i=1,n1
            read(*,*) re,im
            xc1(i)=cmplx(re,im,dp)
        end do
        do i=1,ns
            read(*,*) re,im
            zi1(i)=cmplx(re,im,dp)
        end do
        xr1=real(xc1,dp)
'''
    for n in (1,2,3):
        if n==1:
            xcomp='xc1'; xreal='xr1'; ziprop='zi1'; zfprop='zf1'
        elif n==2:
            header+='''    case(2)
        allocate(xc2(n1,n2),zi2(sd(1),sd(2)))
        do j=1,n2
            do i=1,n1
                read(*,*) re,im
                xc2(i,j)=cmplx(re,im,dp)
            end do
        end do
        do j=1,sd(2)
            do i=1,sd(1)
                read(*,*) re,im
                zi2(i,j)=cmplx(re,im,dp)
            end do
        end do
        xr2=real(xc2,dp)
'''
            xcomp='xc2'; xreal='xr2'; ziprop='zi2'; zfprop='zf2'
        else:
            header+='''    case(3)
        allocate(xc3(n1,n2,n3),zi3(sd(1),sd(2),sd(3)))
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    read(*,*) re,im
                    xc3(i,j,k)=cmplx(re,im,dp)
                end do
            end do
        end do
        do k=1,sd(3)
            do j=1,sd(2)
                do i=1,sd(1)
                    read(*,*) re,im
                    zi3(i,j,k)=cmplx(re,im,dp)
                end do
            end do
        end do
        xr3=real(xc3,dp)
'''
            xcomp='xc3'; xreal='xr3'; ziprop='zi3'; zfprop='zf3'
        header += '        select case(itype)\n'
        for num,kind in enumerate(COMBOS,1):
            b='br' if kind[0]=='r' else 'bc'
            a='ar' if kind[1]=='r' else 'ac'
            x=xreal if kind[2]=='r' else xcomp
            axis=',axis=ax' if n>1 else ''
            header+=f'        case({num})\n            y{n}=lfilter({b},{a},{x},zi={ziprop},zf={zfprop},status=code{axis})\n'
        header += '''        end select
        print *,code
'''
        if n==1:
            header += '''        do i=1,n1
            write(*,'(2(es25.16,1x))') real(y1(i),dp),aimag(y1(i))
        end do
        do i=1,ns
            write(*,'(2(es25.16,1x))') real(zf1(i),dp),aimag(zf1(i))
        end do
'''
        elif n==2:
            header += '''        do j=1,n2
            do i=1,n1
                write(*,'(2(es25.16,1x))') real(y2(i,j),dp),aimag(y2(i,j))
            end do
        end do
        do j=1,sd(2)
            do i=1,sd(1)
                write(*,'(2(es25.16,1x))') real(zf2(i,j),dp),aimag(zf2(i,j))
            end do
        end do
'''
        else:
            header += '''        do k=1,n3
            do j=1,n2
                do i=1,n1
                    write(*,'(2(es25.16,1x))') real(y3(i,j,k),dp),aimag(y3(i,j,k))
                end do
            end do
        end do
        do k=1,sd(3)
            do j=1,sd(2)
                do i=1,sd(1)
                    write(*,'(2(es25.16,1x))') real(zf3(i,j,k),dp),aimag(zf3(i,j,k))
                end do
            end do
        end do
'''
    return header+'    end select\nend program mixed_driver\n'


def main() -> None:
    rng = np.random.default_rng(20261018)
    checked = 0
    max_error = 0.0
    with tempfile.TemporaryDirectory(prefix='scifort_mixed_lfilter_') as tmp:
        work=Path(tmp)
        driver=work/'driver.f90'
        driver.write_text(make_driver())
        exe=work/'driver'
        subprocess.run([FORTRAN,'-std=f2018','-O2','-Werror=implicit-interface',
                        '-ffree-line-length-none',str(ROOT/'src/kinds.f90'),
                        str(ROOT/'src/signal_lfilter.f90'),str(driver),'-o',str(exe)],
                       cwd=work,check=True)
        for case in range(540):
            rank=case%3+1
            combo=COMBOS[(case//3)%len(COMBOS)]
            shape=tuple(int(v) for v in rng.integers(5,11,size=rank))
            axis=case//18 % rank
            na=int(rng.integers(1,5))
            nb=int(rng.integers(1,5))
            ns=max(nb,na)-1
            b=rng.normal(size=nb)*0.5
            a=np.r_[1.0, rng.normal(size=na-1)*0.12]
            if combo[0]=='c': b=b+1j*rng.normal(size=nb)*0.15
            if combo[1]=='c': a=a+1j*np.r_[0.0, rng.normal(size=na-1)*0.1]
            x=rng.normal(size=shape)
            if combo[2]=='c': x=x+1j*rng.normal(size=shape)
            zshape=list(shape)
            zshape[axis]=ns
            zi=rng.normal(size=zshape)+1j*rng.normal(size=zshape)
            yref,zref=signal.lfilter(b,a,x,axis=axis,zi=zi)
            desc=f'{rank} {shape[0]} {shape[1] if rank>1 else 1} {shape[2] if rank>2 else 1} '
            desc+=f'{axis+1} {nb} {na} {COMBOS.index(combo)+1}\n'
            arrays=(b.ravel(order='F'),a.ravel(order='F'),x.ravel(order='F'),zi.ravel(order='F'))
            desc+=''.join(f'{z.real:.17g} {z.imag:.17g}\n' for arr in arrays for z in arr.astype(complex))
            proc=subprocess.run([str(exe)],input=desc,capture_output=True,text=True,check=True)
            lines=proc.stdout.strip().splitlines()
            st=int(lines[0])
            if st!=0: raise AssertionError(f'status {st} for {combo}, rank {rank}')
            values=np.array([list(map(float,line.split())) for line in lines[1:]])
            ny=x.size
            got=(values[:ny,0]+1j*values[:ny,1]).reshape(shape,order='F')
            zgot=(values[ny:,0]+1j*values[ny:,1]).reshape(zshape,order='F')
            np.testing.assert_allclose(got,yref,rtol=2e-12,atol=5e-13,
                err_msg=f'{combo} rank {rank}, axis {axis}, case {case}')
            np.testing.assert_allclose(zgot,zref,rtol=2e-12,atol=5e-13,
                err_msg=f'{combo} final state rank {rank}, case {case}')
            max_error=max(max_error,float(np.max(np.abs(got-yref))),
                          float(np.max(np.abs(zgot-zref))) if zgot.size else 0.0)
            checked+=1
    print(f'Passed {checked} SciPy mixed-coefficient lfilter cases; maximum absolute error {max_error:.3e}')

if __name__=='__main__':
    main()
