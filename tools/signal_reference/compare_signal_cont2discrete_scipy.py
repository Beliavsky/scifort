#!/usr/bin/env python3
"""Compare SciFort continuous-to-discrete conversions with SciPy 1.17.

Compiles one batched Fortran probe (no FPM needed). Tests real/complex,
SISO/MIMO, ZOH, FOH, impulse, GBT, bilinear, Euler, backward difference.
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy.signal import BadCoefficients, cont2discrete
from run_signal_suite import MODULES, ROOT

PROBE = r'''
program test_cont2discrete_probe
    use scifort_signal, only : dp, cont2discrete, cont2discrete_zpk, &
        discrete_state_space_result, discrete_tf_result, discrete_zpk_result
    implicit none
    type(discrete_state_space_result) :: ss
    type(discrete_tf_result) :: tf
    type(discrete_zpk_result) :: zp
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), num(:), den(:), z(:), p(:)
    complex(dp) :: gain
    integer :: ncase, kcase, mode, n, ni, no, j, k, nb, na, nz, np
    character(24) :: method
    real(dp) :: dt, alpha
    read(*,*) ncase
    do kcase = 1, ncase
        read(*,*) mode, method, dt, alpha
        select case(mode)
        case(0)
            read(*,*) n, ni, no
            allocate(a(n,n),b(n,ni),c(no,n),d(no,ni))
            do j=1,n
                do k=1,n
                    call read_value(a(j,k))
                end do
            end do
            do j=1,n
                do k=1,ni
                    call read_value(b(j,k))
                end do
            end do
            do j=1,no
                do k=1,n
                    call read_value(c(j,k))
                end do
            end do
            do j=1,no
                do k=1,ni
                    call read_value(d(j,k))
                end do
            end do
            ss=cont2discrete(a,b,c,d,dt,method=trim(method),alpha=alpha)
            write(*,'(4(i0,1x))') ss%status,n,ni,no
            if (ss%status==0) then
                do j=1,n
                    do k=1,n
                        call write_value(ss%a(j,k))
                    end do
                end do
                do j=1,n
                    do k=1,ni
                        call write_value(ss%b(j,k))
                    end do
                end do
                do j=1,no
                    do k=1,n
                        call write_value(ss%c(j,k))
                    end do
                end do
                do j=1,no
                    do k=1,ni
                        call write_value(ss%d(j,k))
                    end do
                end do
            end if
            deallocate(a,b,c,d)
        case(1)
            read(*,*) nb,na
            allocate(num(nb),den(na))
            do j=1,nb
                call read_value(num(j))
            end do
            do j=1,na
                call read_value(den(j))
            end do
            tf=cont2discrete(num,den,dt,method=trim(method),alpha=alpha)
            n=0
            if (tf%status==0) n=size(tf%a)
            write(*,'(2(i0,1x))') tf%status,n
            if (tf%status==0) then
                do j=1,n
                    call write_value(tf%a(j))
                end do
                do j=1,n
                    call write_value(tf%b(1,j))
                end do
            end if
            deallocate(num,den)
        case(3)
            read(*,*) no,nb,na
            allocate(b(no,nb),den(na))
            do j=1,no
                do k=1,nb
                    call read_value(b(j,k))
                end do
            end do
            do j=1,na
                call read_value(den(j))
            end do
            tf=cont2discrete(b,den,dt,method=trim(method),alpha=alpha)
            n=0
            if (tf%status==0) n=size(tf%a)
            write(*,'(3(i0,1x))') tf%status,no,n
            if (tf%status==0) then
                do j=1,n
                    call write_value(tf%a(j))
                end do
                do j=1,no
                    do k=1,n
                        call write_value(tf%b(j,k))
                    end do
                end do
            end if
            deallocate(b,den)
        case(2)
            read(*,*) nz,np
            allocate(z(nz),p(np))
            do j=1,nz
                call read_value(z(j))
            end do
            do j=1,np
                call read_value(p(j))
            end do
            call read_value(gain)
            zp=cont2discrete_zpk(z,p,gain,dt,method=trim(method),alpha=alpha)
            n=0
            ni=0
            if (zp%status==0) then
                n=size(zp%zeros)
                ni=size(zp%poles)
            end if
            write(*,'(3(i0,1x))') zp%status,n,ni
            if (zp%status==0) then
                call write_value(zp%gain)
                do j=1,n
                    call write_value(zp%zeros(j))
                end do
                do j=1,ni
                    call write_value(zp%poles(j))
                end do
            end if
            deallocate(z,p)
        end select
    end do
contains
    subroutine read_value(x)
        complex(dp), intent(out) :: x !! complex input
        real(dp) :: re,im
        read(*,*) re,im
        x=cmplx(re,im,dp)
    end subroutine read_value
    subroutine write_value(x)
        complex(dp), intent(in) :: x !! complex output
        write(*,'(2(es25.16e3,1x))') real(x,dp),aimag(x)
    end subroutine write_value
end program test_cont2discrete_probe
'''

METHODS = ['zoh', 'foh', 'gbt', 'bilinear', 'euler', 'backward_diff', 'impulse']


def generate():
    rng = np.random.default_rng(20261017)
    items = []
    for n in range(1, 6):
        for trial in range(18):
            method = METHODS[(trial+n) % len(METHODS)]
            ni = 1 + trial % 3
            no = 1 + (trial // 3) % 3
            scale = 0.25 if trial % 2 else 0.9
            a = rng.normal(size=(n,n)) * scale
            a -= np.diag(rng.uniform(.1,2,n))
            b = rng.normal(size=(n,ni))
            c = rng.normal(size=(no,n))
            d = np.zeros((no,ni)) if method == 'impulse' else rng.normal(size=(no,ni))*.4
            if trial % 4 == 1:
                a = a + 1j*rng.normal(size=(n,n))*.25
                b = b + 1j*rng.normal(size=(n,ni))*.4
                c = c + 1j*rng.normal(size=(no,n))*.3
                d = d + 1j*rng.normal(size=(no,ni))*.1
            if method == 'impulse':
                d = np.zeros((no,ni), dtype=complex if trial % 4 == 1 else float)
            dt = float(rng.uniform(.005,.75))
            alpha = float(rng.uniform(0,1))
            items.append(('ss',a,b,c,d,dt,method,alpha))
    # Stress the exponential scaling/squaring branch and nonnormal matrices.
    for j in range(30):
        n = 1 + j % 5
        a = rng.normal(size=(n,n))*1.9 - np.diag(rng.uniform(1,5,n))
        b = rng.normal(size=(n,2))
        c = rng.normal(size=(2,n))
        d = rng.normal(size=(2,2))*.1
        method = ['zoh','foh','bilinear'][j%3]
        items.append(('ss',a,b,c,d,float(rng.uniform(.4,2)),method,0.5))
    # TF systems: well-separated stable poles avoid sensitivity of
    # characteristic polynomial reconstruction at high order.
    for n in range(1,6):
        for j in range(26):
            method = METHODS[(j+2*n)%len(METHODS)]
            complex_mode = j%4 == 0
            poles = -rng.uniform(.2,2.0,n)
            if complex_mode:
                poles = poles+1j*rng.uniform(-.3,.3,n)
            den = np.poly(poles)
            nzero = min(n, j % (n+1))
            if method == 'impulse':
                nzero = min(nzero,n-1)
            zeros = -rng.uniform(.1,1.7,nzero)
            if complex_mode:
                zeros = zeros+1j*rng.uniform(-.4,.4,nzero)
            num = np.atleast_1d(np.poly(zeros)) * (rng.uniform(.2,2)*(1+.2j if complex_mode else 1))
            dt = float(rng.uniform(.01,.8))
            items.append(('tf',num,den,dt,method,float(rng.uniform(0,1))))
    for n in range(1,5):
        for j in range(15):
            method=METHODS[(2*j+n)%len(METHODS)]
            den=np.poly(-rng.uniform(.3,1.8,n))
            num=np.vstack([np.pad(np.atleast_1d(np.poly(-rng.uniform(.1,1.6,jj%n))),
                                  (n-(jj%n),0)) * rng.uniform(.5,2.0)
                           for jj in [j,j+1]])
            if method=='impulse':
                num[:,0]=0.0
            dt=float(rng.uniform(.03,.6))
            items.append(('simo',num,den,dt,method,float(rng.uniform(0,1))))
    for n in range(1,5):
        for j in range(20):
            method=METHODS[(j+n)%len(METHODS)]
            roots = -rng.uniform(.3,2.0,n)
            numz = min(j%(n+1), n-(method=='impulse'))
            zeros = -rng.uniform(.1,1.8,numz)
            if j%5 == 0:
                roots = roots+1j*rng.uniform(-.3,.3,n)
                zeros = zeros+1j*rng.uniform(-.3,.3,numz)
            gain = complex(rng.uniform(.2,1.6),0.1 if j%5 == 0 else 0)
            items.append(('zpk',zeros,roots,gain,float(rng.uniform(.02,.7)),method,
                          float(rng.uniform(0,1))))
    return items


def serialize(items):
    lines=[str(len(items))]
    def val(x):
        x=complex(x)
        lines.append(f'{x.real:.17g} {x.imag:.17g}')
    def mat(m):
        for x in np.asarray(m).flat:
            val(x)
    for item in items:
        if item[0]=='ss':
            _,a,b,c,d,dt,method,alpha=item
            lines.append(f'0 {method} {dt:.17g} {alpha:.17g}')
            lines.append(f'{len(a)} {b.shape[1]} {c.shape[0]}')
            for m in [a,b,c,d]: mat(m)
        elif item[0]=='tf':
            _,num,den,dt,method,alpha=item
            lines.append(f'1 {method} {dt:.17g} {alpha:.17g}')
            lines.append(f'{len(num)} {len(den)}')
            mat(num); mat(den)
        elif item[0]=='simo':
            _,num,den,dt,method,alpha=item
            lines.append(f'3 {method} {dt:.17g} {alpha:.17g}')
            lines.append(f'{num.shape[0]} {num.shape[1]} {len(den)}')
            mat(num);mat(den)
        else:
            _,z,p,k,dt,method,alpha=item
            lines.append(f'2 {method} {dt:.17g} {alpha:.17g}')
            lines.append(f'{len(z)} {len(p)}')
            mat(z);mat(p);val(k)
    return '\n'.join(lines)+'\n'


def read_array(lines,ptr,shape):
    n=int(np.prod(shape))
    data=[complex(*map(float,lines[ptr+i].split())) for i in range(n)]
    return np.asarray(data).reshape(shape),ptr+n


def max_scaled(got,expected):
    assert got.shape == expected.shape,(got.shape,expected.shape)
    if not got.size: return 0.
    return float(np.max(abs(got-expected)/(1+abs(expected))))


def scipy_ref(item):
    if item[0]=='ss':
        _,a,b,c,d,dt,method,alpha=item
        return cont2discrete((a,b,c,d),dt,method=method,**({'alpha':alpha} if method=='gbt' else {}))
    if item[0]=='tf':
        _,num,den,dt,method,alpha=item
        return cont2discrete((num,den),dt,method=method,**({'alpha':alpha} if method=='gbt' else {}))
    if item[0]=='simo':
        _,num,den,dt,method,alpha=item
        return cont2discrete((num,den),dt,method=method,**({'alpha':alpha} if method=='gbt' else {}))
    _,z,p,k,dt,method,alpha=item
    return cont2discrete((z,p,k),dt,method=method,**({'alpha':alpha} if method=='gbt' else {}))


def main():
    items=generate()
    with tempfile.TemporaryDirectory(prefix='scifort-cont2discrete-') as temp:
        w=Path(temp)
        flags=['-std=f2018','-O2','-Wimplicit-interface','-Werror=implicit-interface']
        objects=[]
        for name in MODULES:
            obj=w/f'{name}.o'
            subprocess.run(['gfortran',*flags,'-I',str(w),'-J',str(w),'-c',
                            str(ROOT/'src'/f'{name}.f90'),'-o',str(obj)],check=True,cwd=w)
            objects.append(str(obj))
        probe=w/'compare.f90'
        probe.write_text(PROBE,encoding='ascii')
        exe=w/'compare'
        subprocess.run(['gfortran',*flags,'-I',str(w),'-J',str(w),str(probe),*objects,
                        '-o',str(exe)],check=True,cwd=w)
        output=subprocess.run([str(exe)],input=serialize(items),capture_output=True,
                              text=True,cwd=w,check=True).stdout.splitlines()
    ptr=0
    biggest=0.
    counts={'ss':0,'tf':0,'simo':0,'zpk':0}
    for i,item in enumerate(items):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BadCoefficients)
            ref=scipy_ref(item)
        desc=f'case={i} representation={item[0]} method={item[-2]}'
        head=list(map(int,output[ptr].split()))
        ptr+=1
        if head[0]!=0:
            raise AssertionError(f'{desc}: Fortran status {head[0]}')
        if item[0]=='ss':
            _,n,ni,no=head
            fields=[]
            for shape in [(n,n),(n,ni),(no,n),(no,ni)]:
                field,ptr=read_array(output,ptr,shape)
                fields.append(field)
            error=max(max_scaled(x,np.asarray(y)) for x,y in zip(fields,ref[:4]))
        elif item[0]=='tf':
            n=head[1]
            den,ptr=read_array(output,ptr,(n,))
            num,ptr=read_array(output,ptr,(n,))
            ref_num=np.asarray(ref[0]).reshape(-1)
            ref_den=np.asarray(ref[1]).reshape(-1)
            # SciPy can trim leading numerator/denominator zero coefficients.
            # Align to common padded form for accurate comparison.
            nr=max(n,len(ref_den),len(ref_num))
            error=max(max_scaled(np.pad(num,(nr-n,0)),np.pad(ref_num,(nr-len(ref_num),0))),
                      max_scaled(np.pad(den,(nr-n,0)),np.pad(ref_den,(nr-len(ref_den),0))))
        elif item[0]=='simo':
            _,nout,n=head
            den,ptr=read_array(output,ptr,(n,))
            num,ptr=read_array(output,ptr,(nout,n))
            ref_num=np.asarray(ref[0]).reshape((nout,-1))
            ref_den=np.asarray(ref[1]).reshape(-1)
            nr=max(n,len(ref_den),ref_num.shape[1])
            error=max(max_scaled(np.pad(num,((0,0),(nr-n,0))),
                                 np.pad(ref_num,((0,0),(nr-ref_num.shape[1],0)))),
                      max_scaled(np.pad(den,(nr-n,0)),np.pad(ref_den,(nr-len(ref_den),0))))
        else:
            _,nz,npoles=head
            gain,ptr=read_array(output,ptr,(1,))
            zeros,ptr=read_array(output,ptr,(nz,))
            poles,ptr=read_array(output,ptr,(npoles,))
            rz,rp,rk,_=ref
            # Use polynomial evaluations rather than comparing root order.
            # The root solver may order conjugates differently from SciPy.
            points=1.1*np.exp(1j*np.linspace(-2.5,2.5,21))
            gv=gain[0]*np.prod(points[:,None]-zeros,axis=1)/np.prod(points[:,None]-poles,axis=1)
            rv=rk*np.prod(points[:,None]-rz,axis=1)/np.prod(points[:,None]-rp,axis=1)
            error=max_scaled(gv,rv)
        if error>2e-8:
            raise AssertionError(f'{desc}: scaled discrepancy {error:.3g}')
        biggest=max(biggest,error)
        counts[item[0]]+=1
    assert ptr==len(output),(ptr,len(output))
    print(f'PASS: {len(items)} SciPy 1.17 comparisons, {counts}, max scaled error {biggest:.3g}')


if __name__=='__main__':
    main()
