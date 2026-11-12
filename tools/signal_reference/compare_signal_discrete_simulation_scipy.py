#!/usr/bin/env python3
"""Compare SciFort discrete-time LTI simulation against SciPy 1.17.

The SciPy 1.17 dlsim implementation casts complex state/output to real.
For complex systems we instead use the exact same recurrence with SciPy's
state-space realizations as a reference (and check real systems via SciPy).
"""
from __future__ import annotations

import subprocess
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy import signal
from run_signal_suite import MODULES, ROOT

PROBE = r'''
program discrete_simulation_probe
    use scifort_signal, only : dp, dlsim, dstep, dimpulse, dlsim_zpk, dstep_zpk, dimpulse_zpk, &
        discrete_simulation_result, discrete_response_result
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), u(:,:), x0(:), &
        num(:), den(:), z(:), p(:), u1(:)
    complex(dp) :: gain
    type(discrete_simulation_result) :: sim
    type(discrete_response_result) :: resp
    real(dp), allocatable :: times(:)
    real(dp) :: dt
    integer :: total, idx, mode, nt, has_time, ns, ni, no, nb, na, nz, np, i, j

    read(*,*) total
    do idx = 1, total
        read(*,*) mode, nt, has_time, dt, ns, ni, no, nb, na, nz, np
        if (has_time == 1) then
            allocate(times(nt))
            read(*,*) times
        end if
        if (mode <= 2) then
            allocate(a(ns,ns),b(ns,ni),c(no,ns),d(no,ni),x0(ns))
            do i = 1,ns
                do j = 1,ns
                    call read_value(a(i,j))
                end do
            end do
            do i = 1,ns
                do j = 1,ni
                    call read_value(b(i,j))
                end do
            end do
            do i = 1,no
                do j = 1,ns
                    call read_value(c(i,j))
                end do
            end do
            do i = 1,no
                do j = 1,ni
                    call read_value(d(i,j))
                end do
            end do
            do i = 1,ns
                call read_value(x0(i))
            end do
        else if (mode <= 5) then
            allocate(num(nb),den(na))
            do i = 1,nb
                call read_value(num(i))
            end do
            do i = 1,na
                call read_value(den(i))
            end do
        else
            allocate(z(nz),p(np))
            do i = 1,nz
                call read_value(z(i))
            end do
            do i = 1,np
                call read_value(p(i))
            end do
            call read_value(gain)
        end if
        if (mod(mode,3) == 0) then
            allocate(u(nt,ni))
            do i = 1,nt
                do j = 1,ni
                    call read_value(u(i,j))
                end do
            end do
            if (mode >= 3) then
                allocate(u1(nt))
                u1 = u(:,1)
            end if
        end if
        select case (mode)
        case (0)
            if (has_time==1) then
                sim = dlsim(a,b,c,d,u,dt,t=times,x0=x0)
            else
                sim = dlsim(a,b,c,d,u,dt,x0=x0)
            end if
        case (1)
            if (has_time==1) then
                resp = dstep(a,b,c,d,dt,t=times,x0=x0)
            else
                resp = dstep(a,b,c,d,dt,n=nt,x0=x0)
            end if
        case (2)
            if (has_time==1) then
                resp = dimpulse(a,b,c,d,dt,t=times,x0=x0)
            else
                resp = dimpulse(a,b,c,d,dt,n=nt,x0=x0)
            end if
        case (3)
            if (has_time==1) then
                sim = dlsim(num,den,u1,dt,t=times)
            else
                sim = dlsim(num,den,u1,dt)
            end if
        case (4)
            if (has_time==1) then
                resp = dstep(num,den,dt,t=times)
            else
                resp = dstep(num,den,dt,n=nt)
            end if
        case (5)
            if (has_time==1) then
                resp = dimpulse(num,den,dt,t=times)
            else
                resp = dimpulse(num,den,dt,n=nt)
            end if
        case (6)
            if (has_time==1) then
                sim = dlsim_zpk(z,p,gain,u1,dt,t=times)
            else
                sim = dlsim_zpk(z,p,gain,u1,dt)
            end if
        case (7)
            if (has_time==1) then
                resp = dstep_zpk(z,p,gain,dt,t=times)
            else
                resp = dstep_zpk(z,p,gain,dt,n=nt)
            end if
        case (8)
            if (has_time==1) then
                resp = dimpulse_zpk(z,p,gain,dt,t=times)
            else
                resp = dimpulse_zpk(z,p,gain,dt,n=nt)
            end if
        end select
        if (mod(mode,3) == 0) then
            call write_sim(sim)
        else
            call write_response(resp)
        end if
        if (allocated(a)) deallocate(a,b,c,d,x0)
        if (allocated(num)) deallocate(num,den)
        if (allocated(z)) deallocate(z,p)
        if (allocated(times)) deallocate(times)
        if (allocated(u)) deallocate(u)
        if (allocated(u1)) deallocate(u1)
    end do
contains
    subroutine read_value(v)
        complex(dp), intent(out) :: v !! complex input value
        real(dp) :: re, im
        read(*,*) re,im
        v=cmplx(re,im,dp)
    end subroutine read_value
    subroutine write_value(v)
        complex(dp), intent(in) :: v !! complex output value
        write(*,'(2(es25.16e3,1x))') real(v,dp),aimag(v)
    end subroutine write_value
    subroutine write_sim(s)
        type(discrete_simulation_result), intent(in) :: s !! simulation result
        integer :: k,l
        if (s%status/=0) then
            write(*,'(4(i0,1x))') s%status,0,0,0
            return
        end if
        write(*,'(4(i0,1x))') s%status,size(s%t),size(s%y,2),size(s%x,2)
        do k=1,size(s%t)
            write(*,'(es25.16e3)') s%t(k)
            do l=1,size(s%y,2)
                call write_value(s%y(k,l))
            end do
            do l=1,size(s%x,2)
                call write_value(s%x(k,l))
            end do
        end do
    end subroutine write_sim
    subroutine write_response(s)
        type(discrete_response_result), intent(in) :: s !! response result
        integer :: k,l,m
        if (s%status/=0) then
            write(*,'(4(i0,1x))') s%status,0,0,0
            return
        end if
        write(*,'(4(i0,1x))') s%status,size(s%t),size(s%y,2),size(s%y,3)
        do k=1,size(s%t)
            write(*,'(es25.16e3)') s%t(k)
            do m=1,size(s%y,3)
                do l=1,size(s%y,2)
                    call write_value(s%y(k,l,m))
                end do
            end do
        end do
    end subroutine write_response
end program discrete_simulation_probe
'''


def build_cases():
    rng = np.random.default_rng(240682)
    cases = []
    # Moderate stable matrices, both real and complex, including MIMO.
    for mode in range(3):
        for k in range(72):
            n = 1 + k % 4
            ni = 1 + k % 3
            no = 1 + k % 2
            cmplx = k % 4 == 1
            A = rng.normal(size=(n,n))*.13 + np.eye(n)*.4
            B = rng.normal(size=(n,ni))*.5
            C = rng.normal(size=(no,n))*.6
            D = rng.normal(size=(no,ni))*.2
            x0 = rng.normal(size=n)*.4
            nsample = 3+k%13
            u = rng.normal(size=(nsample,ni))
            if cmplx:
                A = A+1j*rng.normal(size=(n,n))*.2
                B = B+1j*rng.normal(size=(n,ni))*.3
                C = C+1j*rng.normal(size=(no,n))*.3
                D = D+1j*rng.normal(size=(no,ni))*.2
                x0 = x0+1j*rng.normal(size=n)*.2
                u = u+1j*rng.normal(size=(nsample,ni))*.6
            dt = .05+.1*(k%5)
            ti = None
            if k % 4 == 3:
                ti = np.cumsum(rng.uniform(.6,1.5,size=nsample))*dt
                ti -= ti[0]
            cases.append((mode,dt,ti,(A,B,C,D,x0,u,nsample)))
    for mode in range(3,6):
        for k in range(54):
            n = k%5
            roots = rng.uniform(-.6,.6,size=n)
            z = rng.uniform(-.7,.7,size=k%(n+1))
            cmplx = k%3==0
            if cmplx:
                roots = roots+1j*rng.normal(size=n)*.16
                z = z+1j*rng.normal(size=len(z))*.15
            num=np.atleast_1d(np.poly(z))*(1+.2j if cmplx else 1.)
            den=np.atleast_1d(np.poly(roots))
            nsample = 3+k%15
            u=rng.normal(size=nsample)
            if cmplx: u=u+1j*rng.normal(size=nsample)*.5
            dt=.1+.07*(k%5)
            ti=None
            if k%4==2:
                ti = np.cumsum(rng.uniform(.6,1.5,size=nsample))*dt
                ti -= ti[0]
            cases.append((mode,dt,ti,(num,den,u,nsample)))
    for mode in range(6,9):
        for k in range(40):
            n = k%5
            p = rng.uniform(-.6,.6,size=n)
            z = rng.uniform(-.5,.5,size=k%(n+1))
            cmplx=k%3==1
            if cmplx:
                p=p+1j*rng.normal(size=len(p))*.1
                z=z+1j*rng.normal(size=len(z))*.1
            gain=complex(rng.uniform(.2,1.5),rng.uniform(-.2,.2) if cmplx else 0.)
            nsample=3+k%13
            u=rng.normal(size=nsample)
            if cmplx: u=u+1j*rng.normal(size=nsample)*.2
            dt=.05+.08*(k%4)
            ti = None
            if k%4==3:
                ti=np.cumsum(rng.uniform(.5,1.4,size=nsample))*dt
                ti-=ti[0]
            cases.append((mode,dt,ti,(z,p,gain,u,nsample)))
    return cases


def reference_sim(A,B,C,D,u,dt,ti,x0):
    if ti is None:
        tout=np.arange(len(u))*dt
        ui=u
    else:
        stop=ti[-1]
        nt=int(np.floor(stop/dt))+1
        tout=np.linspace(0,stop,nt)
        # SciPy uses a linear B-spline which extrapolates linearly.
        from scipy.interpolate import make_interp_spline
        ui=make_interp_spline(ti,u,k=1)(tout)
    x=np.zeros((len(tout),len(A)),dtype=complex)
    y=np.zeros((len(tout),C.shape[0]),dtype=complex)
    x[0]=x0
    for i in range(len(tout)):
        y[i]=C@x[i]+D@ui[i]
        if i+1<len(tout): x[i+1]=A@x[i]+B@ui[i]
    return tout,y,x


def reference(case):
    mode,dt,ti,data=case
    if mode<=2:
        A,B,C,D,x0,u,n=data
        original=(A,B,C,D,dt)
        is_complex=any(np.iscomplexobj(item) for item in (A,B,C,D,x0,u))
    elif mode<=5:
        num,den,u,n=data
        A,B,C,D=signal.tf2ss(num,den)
        x0=np.zeros(A.shape[0])
        original=(num,den,dt)
        is_complex=any(np.iscomplexobj(item) for item in (num,den,u))
    else:
        z,p,gain,u,n=data
        A,B,C,D=signal.zpk2ss(z,p,gain)
        x0=np.zeros(A.shape[0])
        original=(z,p,gain,dt)
        is_complex=any(np.iscomplexobj(item) for item in (z,p)) or gain.imag!=0 or np.iscomplexobj(u)
    if mode%3==0:
        ui=np.atleast_2d(u).T if u.ndim==1 else u
        if not is_complex:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                if mode==0:
                    return signal.dlsim(original,ui,t=ti,x0=x0)
                else:
                    tout,y=signal.dlsim(original,ui,t=ti)
                    _,_,x=reference_sim(A,B,C,D,ui,dt,ti,x0)
                    return tout,y,x
        return reference_sim(A,B,C,D,ui,dt,ti,x0)
    noinputs=B.shape[1]
    if not is_complex:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if mode%3==1:
                result = signal.dstep(original,t=ti,n=n,x0=x0)
            else:
                result = signal.dimpulse(original,t=ti,n=n,x0=x0)
            # On some dt values SciPy's endpoint=False floating grid rounds
            # (n-1)*dt slightly low, dropping the final sample via floor.
            # Retain precisely n samples as the documented API requests.
            if ti is not None or len(result[0]) == n:
                return result
    ys=[]
    for channel in range(noinputs):
        ui=np.zeros((n,noinputs),dtype=complex)
        if mode%3==1: ui[:,channel]=1
        else: ui[0,channel]=1
        out=reference_sim(A,B,C,D,ui,dt,ti,x0)
        ys.append(out[1])
    return out[0],tuple(ys)


def to_stdin(cases):
    lines=[str(len(cases))]
    def val(v):
        v=complex(v)
        lines.append(f'{v.real:.17g} {v.imag:.17g}')
    for mode,dt,ti,data in cases:
        if mode<=2:
            A,B,C,D,x0,u,n=data
            ns=A.shape[0]
            ni=B.shape[1]
            no=C.shape[0]
            nb=na=nz=npoles=0
        elif mode<=5:
            num,den,u,n=data
            ns=ni=no=1
            nb=len(num)
            na=len(den)
            nz=npoles=0
        else:
            z,p,k,u,n=data
            ns=ni=no=1
            nz=len(z)
            npoles=len(p)
            nb=na=0
        lines.append(f'{mode} {n} {int(ti is not None)} {dt:.17g} {ns} {ni} {no} {nb} {na} {nz} {npoles}')
        if ti is not None: lines.append(' '.join(f'{t:.17g}' for t in ti))
        if mode<=2:
            for mat in (A,B,C,D):
                for v in mat.flat:val(v)
            for v in x0:val(v)
        elif mode<=5:
            for v in num:val(v)
            for v in den:val(v)
        else:
            for v in z:val(v)
            for v in p:val(v)
            val(k)
        if mode%3==0:
            for v in u.flat:val(v)
    return '\n'.join(lines)+'\n'


def parse_stdout(output,cases):
    values=iter(output.splitlines())
    result=[]
    for case in cases:
        mode=case[0]
        status,nt,no,nlast=map(int,next(values).split())
        if status!=0:
            raise AssertionError(f'{mode}: nonzero Fortran status {status}')
        ts=np.empty(nt)
        if mode%3==0:
            y=np.empty((nt,no),complex)
            x=np.empty((nt,nlast),complex)
        else:
            y=np.empty((nt,no,nlast),complex)
        def readvalue():
            re,im=map(float,next(values).split())
            return complex(re,im)
        for i in range(nt):
            ts[i]=float(next(values))
            if mode%3==0:
                for j in range(no):y[i,j]=readvalue()
                for j in range(nlast):x[i,j]=readvalue()
            else:
                for k in range(nlast):
                    for j in range(no):y[i,j,k]=readvalue()
        if mode%3==0:result.append((ts,y,x))
        else:result.append((ts,y))
    return result


def main():
    cases=build_cases()
    with tempfile.TemporaryDirectory(prefix='scifort-discrete-sim-') as directory:
        build=Path(directory)
        objects=[]
        flags=['-std=f2018','-O0','-fcheck=all','-Wimplicit-interface','-Werror=implicit-interface']
        for module in MODULES:
            obj=build/(module+'.o')
            subprocess.run(['gfortran',*flags,'-J',str(build),'-I',str(build),
                '-c',str(ROOT/'src'/f'{module}.f90'),'-o',str(obj)],check=True,cwd=build)
            objects.append(str(obj))
        src=build/'probe.f90'
        src.write_text(PROBE)
        exe=build/'probe.exe'
        subprocess.run(['gfortran',*flags,'-I',str(build),str(src),*objects,'-o',str(exe)],check=True,cwd=build)
        proc=subprocess.run([str(exe)],input=to_stdin(cases),capture_output=True,text=True,cwd=build)
        if proc.returncode:
            raise AssertionError(f'Fortran probe failed: {proc.stderr[-2000:]}')
        results=parse_stdout(proc.stdout,cases)
        max_scaled=0.
        max_absolute=0.
        for i,(case,actual) in enumerate(zip(cases,results)):
            mode=case[0]
            ref=reference(case)
            if mode%3==0:
                ta,ya,xa=actual
                tr,yr,xr=ref
                arrays=[(ta,tr),(ya,yr),(xa,xr)]
            else:
                ta,ya=actual
                tr,yr=ref
                arrays=[(ta,tr),(ya,np.stack(yr,axis=-1))]
            for got,expect in arrays:
                if got.shape!=expect.shape:
                    raise AssertionError(f'case {i} mode {mode}: shape {got.shape} vs {expect.shape}')
                if not got.size:continue
                absolute=float(np.max(np.abs(got-expect)))
                scaled=float(np.max(np.abs(got-expect)/(1+np.abs(expect))))
                max_absolute=max(max_absolute,absolute)
                max_scaled=max(max_scaled,scaled)
                if not np.allclose(got,expect,rtol=3.e-10,atol=1.e-10):
                    raise AssertionError(f'case {i} mode {mode}: max error {absolute:.3g}')
    print(f'PASS: {len(cases)} discrete simulation cases against SciPy 1.17 (real),')
    print(f'      SciPy state-space recurrence (complex); max absolute={max_absolute:.3e},')
    print(f'      max scaled={max_scaled:.3e}')


if __name__=='__main__':
    main()
