! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_continuation
    use scifort_kinds, only: dp
    use scifort_interpolate, only: splrep,splprep,splev,fitpack_tck_type, &
        fitpack_parametric_tck_type,fitpack_continuation_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    integer, parameter :: m=61
    type(fitpack_tck_type) :: fit,again
    type(fitpack_parametric_tck_type) :: curve
    type(fitpack_continuation_type) :: state,copy,other
    real(dp) :: u(m),y(m),w(m),x(m,3),bad(m),target,energy,fp,nan
    real(dp), allocatable :: v(:,:),saved(:)
    real(dp), parameter :: factors(6)=[10.0_dp,0.001_dp,0.002_dp,10.0_dp,0.0_dp,0.0015_dp]
    integer :: i,j,k,mode,code,checks,task,used
    logical :: periodic
    checks=0
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]**1.1_dp
    y=cos(2*acos(-1.0_dp)*u)+0.08_dp*sin(14*acos(-1.0_dp)*u)
    y(m)=y(1)
    x(:,1)=y
    x(:,2)=2*y
    x(:,3)=-0.5_dp*y
    w=0.8_dp+u
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    call require(.not.state%valid(),'default state invalid')
    fit=splrep(u,y,s=0.1_dp,task=1,state=state,status=code)
    call require(code==1 .and. .not.fit%valid(),'continuation needs initialized state')
    do k=1,5
        do mode=0,1
            periodic=mode==1
            used=m-mode
            energy=sum((w(:used)*y(:used))**2)
            call state%reset()
            call other%reset()
            do j=1,size(factors)
                task=1
                if (j==1) task=0
                target=factors(j)*energy
                fit=splrep(u,y,w=w,k=k,s=target,per=periodic,task=task,state=state,status=code)
                if (code/=0) write(*,*) 'Failed scalar degree, periodic mode, step, status:',k,mode,j,code
                call require(code==0 .and. fit%valid() .and. state%valid(),'scalar sequence success')
                fp=sum((w(:used)*(splev(u(:used),fit)-y(:used)))**2)
                call require(fp<=1.0011_dp*target+1e-22_dp,'scalar sequence residual')
                curve=splprep(x,u=u,w=w,k=k,s=5.25_dp*target,per=periodic,task=task, &
                    state=other,status=code)
                if (code/=0) write(*,*) 'Failed parametric degree, periodic mode, step, status:',k,mode,j,code
                call require(code==0 .and. curve%valid() .and. other%valid(),'parametric sequence success')
                v=splev(u(:used),curve)
                fp=sum((spread(w(:used),2,3)*(v-x(:used,:)))**2)
                call require(fp<=1.0011_dp*5.25_dp*target+1e-21_dp,'shared sequence residual')
                call require(all(curve%u==u),'parameter sites retained')
                do i=0,k-1
                    if (.not.periodic) exit
                    call require(abs(splev(0.0_dp,fit,der=i)-splev(1.0_dp,fit,der=i)) < &
                        1e-8_dp*max(1.0_dp,maxval(abs(splev(u,fit,der=i)))),'continued periodic seam')
                end do
            end do
            copy=state
            fit=splrep(u,y,w=w,k=k,s=0.001_dp*energy,per=periodic,task=1,state=state,status=code)
            again=splrep(u,y,w=w,k=k,s=0.001_dp*energy,per=periodic,task=1,state=copy,status=code)
            call require(code==0 .and. again%valid(),'copied state resumes independently')
            call require(all(fit%t==again%t) .and. all(fit%c==again%c),'copied state deterministic')
            saved=splev(u,fit)
            bad=y
            bad(10)=bad(10)+0.1_dp
            again=splrep(u,bad,w=w,k=k,s=target,per=periodic,task=1,state=state,status=code)
            call require(code==1 .and. .not.again%valid() .and. state%valid(),'changed data rejected transactionally')
            bad=w
            bad(10)=2*bad(10)
            again=splrep(u,y,w=bad,k=k,s=target,per=periodic,task=1,state=state,status=code)
            call require(code==1 .and. .not.again%valid(),'changed weights rejected')
            again=splrep(u,y,w=w,k=modulo(k,5)+1,s=target,per=periodic, &
                task=1,state=state,status=code)
            call require(code==1 .and. .not.again%valid(),'changed degree rejected')
            again=splrep(u,y,w=w,k=k,s=target,per=.not.periodic,task=1,state=state,status=code)
            call require(code==1 .and. .not.again%valid(),'changed boundary mode rejected')
            bad=u
            bad(10)=bad(10)+1e-4_dp
            again=splrep(bad,y,w=w,k=k,s=target,per=periodic,task=1,state=state,status=code)
            call require(code==1 .and. .not.again%valid(),'changed sites rejected')
            curve=splprep(x,u=u,w=w,k=k,s=target,per=periodic,task=1,state=state,status=code)
            call require(code==1 .and. .not.curve%valid(),'scalar state cannot resume parametric fitter')
            again=splrep(u,y,w=w,k=k,s=tiny(1.0_dp),per=periodic,task=1,state=state,status=code)
            call require(code==3 .and. .not.again%valid() .and. state%valid(), &
                'numerical range failure preserves state')
            again=splrep(u,y,w=w,k=k,s=nan,per=periodic,task=1,state=state,status=code)
            call require(code==1 .and. .not.again%valid(),'invalid target preserves state')
            again=splrep(u,y,w=w,k=k,s=0.001_dp*energy,per=periodic,task=1,state=state,status=code)
            call require(code==0 .and. maxval(abs(splev(u,again)-saved))<1e-10_dp,'state usable after failure')
            if (periodic) then
                bad=y
                bad(m)=nan
                again=splrep(u,bad,w=w,k=k,s=target,per=.true.,task=1,state=state,status=code)
                call require(code==0 .and. again%valid(),'ignored closing data may change')
                bad=w
                bad(m)=nan
                again=splrep(u,y,w=bad,k=k,s=target,per=.true.,task=1,state=state,status=code)
                call require(code==0 .and. again%valid(),'ignored closing weight may change')
            end if
        end do
    end do
    fit=splrep(u,y,task=1,status=code)
    call require(code==1 .and. .not.fit%valid(),'task one without state rejected')
    fit=splrep(u,y,s=0.1_dp,task=-1,t=[0.5_dp],state=state,status=code)
    call require(code==1 .and. .not.fit%valid(),'fixed knots cannot initialize continuation')
    call state%reset()
    call state%reset()
    call require(.not.state%valid(),'reset is idempotent')
    curve=splprep(x,w=w,k=3,s=0.05_dp,state=state,status=code)
    call require(code==0 .and. curve%valid(),'initial automatic parameter state')
    copy=state
    curve=splprep(x,w=w,k=3,s=0.02_dp,task=1,state=state,status=code)
    call require(code==0 .and. curve%valid(),'continued automatic parameters')
    fit=splrep(u,y,k=3,s=0.05_dp,bbox=[-0.1_dp,1.1_dp],state=state,status=code)
    call require(code==0 .and. fit%valid(),'fresh state replaces prior fitting family')
    again=splrep(u,y,k=3,s=0.02_dp,task=1,state=state,status=code)
    call require(code==1 .and. .not.again%valid(),'changed enclosing bounds rejected')
    again=splrep(u,y,s=-1.0_dp,state=state,status=code)
    call require(code==1 .and. state%valid(),'failed fresh initialization preserves prior state')
    again=splrep(u,y,k=3,s=0.02_dp,bbox=[-0.1_dp,1.1_dp],task=1,state=state,status=code)
    call require(code==0 .and. again%valid(),'enlarged-bound state survives rejected calls')
    write(*,'(a,i0,a)') 'Explicit spline continuation: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Whether the contract holds.
        character(*), intent(in) :: label !! Failure description.
        checks=checks+1
        if (.not.ok) then
            write(*,'(a)') 'FAIL: '//label
            stop 1
        end if
    end subroutine require
end program test_interpolate_continuation
