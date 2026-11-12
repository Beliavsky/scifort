#!/usr/bin/env python3
"""Randomized SciFort continuous-time LTI comparisons with SciPy 1.17.

Compile one Fortran probe, then compare state-space, TF and ZPK lsim,
step and impulse outputs. For MIMO response curves SciPy currently does not
support independent input experiments, so use scipy.linalg.expm directly.
"""
from __future__ import annotations

from pathlib import Path
import warnings
from tempfile import TemporaryDirectory
import subprocess
import numpy as np
from scipy import linalg, signal

from run_signal_suite import MODULES, ROOT

PROBE = r'''
program continuous_lti_probe
    use scifort_signal, only : dp, lsim, step, impulse, lsim_zpk, step_zpk, impulse_zpk, &
        continuous_simulation_result, continuous_response_result
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:), u(:,:), x0(:), num(:), den(:), &
        z(:), p(:)
    complex(dp) :: gain
    real(dp), allocatable :: t(:)
    type(continuous_simulation_result) :: sim
    type(continuous_response_result) :: resp
    integer :: idx, total, mode, ns, ni, no, nt, nb, na, nz, np, interp, i, j

    read(*,*) total
    do idx=1,total
        read(*,*) mode, ns, ni, no, nt, nb, na, nz, np, interp
        allocate(t(nt),x0(ns),u(nt,ni))
        read(*,*) t
        if (mode <= 2) then
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
        else if (mode <= 5) then
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
        do i=1,ns
            call read_complex(x0(i))
        end do
        if (mod(mode,3)==0) then
            do i=1,nt
                do j=1,ni
                    call read_complex(u(i,j))
                end do
            end do
        end if
        select case(mode)
        case(0)
            sim=lsim(a,b,c,d,u,t,x0,interp==1)
        case(1)
            resp=step(a,b,c,d,t,x0)
        case(2)
            resp=impulse(a,b,c,d,t,x0)
        case(3)
            sim=lsim(num,den,u(:,1),t,x0,interp==1)
        case(4)
            resp=step(num,den,t,x0)
        case(5)
            resp=impulse(num,den,t,x0)
        case(6)
            sim=lsim_zpk(z,p,gain,u(:,1),t,x0,interp==1)
        case(7)
            resp=step_zpk(z,p,gain,t,x0)
        case(8)
            resp=impulse_zpk(z,p,gain,t,x0)
        end select
        if (mod(mode,3)==0) then
            call write_sim(sim)
        else
            call write_response(resp)
        end if
        if (allocated(a)) deallocate(a,b,c,d)
        if (allocated(num)) deallocate(num,den)
        if (allocated(z)) deallocate(z,p)
        deallocate(t,x0,u)
    end do
contains
    subroutine read_complex(v)
        complex(dp), intent(out) :: v !! complex input
        real(dp) :: re, im
        read(*,*) re,im
        v=cmplx(re,im,dp)
    end subroutine read_complex
    subroutine write_complex(v)
        complex(dp), intent(in) :: v !! complex output
        write(*,'(2(es25.16e3,1x))') real(v,dp),aimag(v)
    end subroutine write_complex
    subroutine write_sim(r)
        type(continuous_simulation_result), intent(in) :: r !! simulation results
        integer :: ii,jj
        if (r%status /= 0) then
            write(*,'(4(i0,1x))') r%status,0,0,0
            return
        end if
        write(*,'(4(i0,1x))') r%status,size(r%t),size(r%y,2),size(r%x,2)
        do ii=1,size(r%t)
            write(*,'(es25.16e3)') r%t(ii)
            do jj=1,size(r%y,2)
                call write_complex(r%y(ii,jj))
            end do
            do jj=1,size(r%x,2)
                call write_complex(r%x(ii,jj))
            end do
        end do
    end subroutine write_sim
    subroutine write_response(r)
        type(continuous_response_result), intent(in) :: r !! response results
        integer :: ii,jj,kk
        if (r%status /= 0) then
            write(*,'(4(i0,1x))') r%status,0,0,0
            return
        end if
        write(*,'(4(i0,1x))') r%status,size(r%t),size(r%y,2),size(r%y,3)
        do ii=1,size(r%t)
            write(*,'(es25.16e3)') r%t(ii)
            do kk=1,size(r%y,3)
                do jj=1,size(r%y,2)
                    call write_complex(r%y(ii,jj,kk))
                end do
            end do
        end do
    end subroutine write_response
end program continuous_lti_probe
'''


def pack(case):
    mode, t, x0 = case['mode'], case['t'], case['x0']
    if mode <= 2:
        a, b, c, d = (case[k] for k in ('a', 'b', 'c', 'd'))
        ns, ni, no = a.shape[0], b.shape[1], c.shape[0]
        params = [a, b, c, d]
        nb = na = nz = npoles = 0
    elif mode <= 5:
        num, den = case['num'], case['den']
        ns, ni, no = len(den) - 1, 1, 1
        nb, na, nz, npoles = len(num), len(den), 0, 0
        params = [num, den]
    else:
        z, p, gain = case['z'], case['p'], case['gain']
        ns, ni, no = len(p), 1, 1
        nb = na = 0
        nz, npoles = len(z), len(p)
        params = [z, p, np.array([gain])]
    lines = [f'{mode} {ns} {ni} {no} {len(t)} {nb} {na} {nz} {npoles} {int(case.get("interp",True))}',
             ' '.join(format(v,'.17g') for v in t)]
    for values in [*params, x0, case.get('u', np.empty((0, ni)))]:
        for v in np.asarray(values).ravel():
            lines.append(f'{v.real:.17g} {v.imag:.17g}')
    return '\n'.join(lines)


def parse(lines, case):
    status, nt, no, nx = map(int, next(lines).split())
    assert status == 0, f'Fortran failed with status {status} on case {case["mode"]}'
    assert nt == len(case['t'])
    y = np.empty((nt, no, nx if case['mode']%3 else 1), complex)
    x = np.empty((nt, nx if case['mode']%3 == 0 else 0), complex)
    for k in range(nt):
        time = float(next(lines))
        assert abs(time-case['t'][k]) <= 1e-11*max(1,abs(time))
        if case['mode'] % 3:
            for j in range(nx):
                for i in range(no):
                    y[k,i,j] = complex(*map(float, next(lines).split()))
        else:
            for i in range(no):
                y[k,i,0] = complex(*map(float, next(lines).split()))
            for i in range(nx):
                x[k,i] = complex(*map(float, next(lines).split()))
    return y, x


def scipy_expm_lsim(a,b,c,d,u,t,x0,interp):
    """Independent exact sampled-input reference using scipy.linalg.expm."""
    n,ni=b.shape
    out_x=np.empty((len(t),n),complex)
    out_x[0]=linalg.expm(a*t[0])@x0
    if len(t)>1:
        dt=t[1]-t[0]
        extra=ni if interp else 0
        big=np.zeros((n+ni+extra,n+ni+extra),complex)
        big[:n,:n]=a*dt
        big[:n,n:n+ni]=b*dt
        if interp:big[n:n+ni,n+ni:]=np.eye(ni)
        em=linalg.expm(big)
        ad=em[:n,:n]
        rb=em[:n,n+ni:] if interp else np.zeros((n,ni),complex)
        lb=em[:n,n:n+ni]-rb
        for j in range(1,len(t)):
            out_x[j]=ad@out_x[j-1]+lb@u[j-1]+rb@u[j]
    out_y=out_x@c.T+u@d.T
    return out_y.reshape(len(t),c.shape[0],1),out_x

def expected(case):
    mode, t, x0 = case['mode'], case['t'], case['x0']
    if mode <= 2:
        a,b,c,d = (case[k] for k in ('a','b','c','d'))
        if mode == 0:
            if np.any(np.imag(case['u'])) or np.any(np.imag(x0)):
                return scipy_expm_lsim(a,b,c,d,case['u'],t,x0,case['interp'])
            tr,y,x=signal.lsim((a,b,c,d),U=case['u'],T=t, X0=x0,interp=case['interp'])
            return np.asarray(y).reshape(len(t),c.shape[0],1), np.asarray(x).reshape(len(t),len(x0))
        outputs = np.empty((len(t), c.shape[0], b.shape[1]),complex)
        ns,ni = b.shape
        if mode == 1:
            big=np.zeros((ns+ni,ns+ni),complex)
            big[:ns,:ns]=a
            big[:ns,ns:]=b
            for k,tk in enumerate(t):
                exp_big=linalg.expm(big*tk)
                states=exp_big[:ns,:ns]@x0[:,None] + exp_big[:ns,ns:]
                outputs[k,:,:] = c@states+d
        else:
            for k,tk in enumerate(t):
                states=linalg.expm(a*tk)@(x0[:,None]+b)
                outputs[k,:,:]=c@states
        return outputs, None
    if mode <= 5:
        system = (case['num'],case['den'])
    else:
        system = (case['z'],case['p'],case['gain'])
    if mode % 3 == 0:
        if mode <= 5:
            a,b,c,d=signal.tf2ss(*system)
        else:
            a,b,c,d=signal.zpk2ss(*system)
        if np.any(np.imag(case['u'])) or np.any(np.imag(x0)):
            ey,_=scipy_expm_lsim(a,b,c,d,case['u'],t,x0,case['interp'])
            return ey,None
        tr,y,x=signal.lsim(system,U=case['u'][:,0],T=t,X0=x0,interp=case['interp'])
        return np.asarray(y).reshape(len(t),1,1), None
    if mode % 3 == 1:
        tr,y=signal.step(system,T=t,X0=x0)
    else:
        tr,y=signal.impulse(system,T=t)
    return np.asarray(y).reshape(len(t),1,1),None


def main():
    rng=np.random.default_rng(20261017)
    cases=[]
    for kind in range(9):
        for k in range(58):
            nt=int(rng.integers(3,21))
            if kind <=2:
                n=int(rng.integers(1,5))
                ni=int(rng.integers(1,4))
                no=int(rng.integers(1,4))
                cx=k%4==0
                factor=1j if cx else 0.
                a=-(0.7+0.4*rng.random())*np.eye(n)+0.12*rng.standard_normal((n,n))
                a=a+factor*0.15*rng.standard_normal((n,n))
                b=rng.standard_normal((n,ni))+factor*rng.standard_normal((n,ni))
                c=rng.standard_normal((no,n))+factor*rng.standard_normal((no,n))
                d=0.2*rng.standard_normal((no,ni))+factor*0.2*rng.standard_normal((no,ni))
                x0=rng.standard_normal(n)+factor*0.3*rng.standard_normal(n)
                payload=dict(a=a,b=b,c=c,d=d,x0=x0)
            else:
                n=int(rng.integers(1,5))
                poles=-rng.uniform(0.5,2.,n)
                zeros=-rng.uniform(0.2,2.,int(rng.integers(0,n+1)))
                gain=rng.uniform(0.1,2.)
                if k%5==0:
                    poles=poles+0.15j
                    zeros=zeros+0.2j
                    gain=gain+0.12j
                if kind <=5:
                    payload=dict(num=np.atleast_1d(np.poly(zeros)*gain),den=np.asarray(np.poly(poles)))
                else:
                    payload=dict(z=np.asarray(zeros),p=np.asarray(poles),gain=gain)
                payload['x0']=np.zeros(n,complex)
            if kind%3==0:
                t0 = 0. if k%3 else 0.12
                t=t0+np.arange(nt)*(0.03+rng.uniform(0.01,0.12))
                u=rng.standard_normal((nt, payload['b'].shape[1] if kind<=2 else 1))
                if k%4==0:u=u+1j*rng.standard_normal(u.shape)
                payload.update(u=u,interp=bool(k%2))
            else:
                t=np.cumsum(rng.uniform(0.03,0.2,nt))
                if k%3==0: t[0]=0.
                if kind>2: t=np.linspace(0.,1.5,nt)
                payload['u']=np.empty((0,payload['b'].shape[1] if kind<=2 else 1))
            payload.update(mode=kind,t=t)
            cases.append(payload)
    with TemporaryDirectory(prefix='scifort-continuous-check-') as d:
        temp=Path(d)
        sources=[]
        for name in MODULES:
            f=ROOT/'src'/f'{name}.f90'
            obj=temp/f'{name}.o'
            subprocess.run(['gfortran','-std=f2018','-O2','-Wimplicit-interface','-Werror=implicit-interface',
                            '-I',str(temp),'-J',str(temp),'-c',str(f),'-o',str(obj)],cwd=temp,check=True)
            sources.append(str(obj))
        probe=temp/'probe.f90'
        probe.write_text(PROBE)
        exe=temp/'probe'
        subprocess.run(['gfortran','-std=f2018','-O2','-I',str(temp),'-J',str(temp),str(probe),
                        *sources,'-o',str(exe)],cwd=temp,check=True)
        command_input=str(len(cases))+'\n'+'\n'.join(pack(c) for c in cases)+'\n'
        result=subprocess.run([str(exe)],input=command_input,text=True,capture_output=True,check=True)
        lines=iter(result.stdout.splitlines())
        worst=0.
        for index,case in enumerate(cases):
            y,x=parse(lines,case)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=np.exceptions.ComplexWarning)
                ey,ex=expected(case)
            error=np.max(np.abs(y-ey)/(1.+np.abs(ey)))
            if ex is not None:
                error=max(error,np.max(np.abs(x-ex)/(1.+np.abs(ex))))
            worst=max(worst,error)
            if not np.isfinite(error) or error>2e-8:
                raise AssertionError(f'case {index} mode {case["mode"]}: relative scaled error {error}')
        print(f'PASS: {len(cases)} SciPy continuous LTI comparisons, worst scaled error={worst:.3e}')

if __name__=='__main__':
    main()
