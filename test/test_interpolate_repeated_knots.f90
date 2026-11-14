! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_repeated_knots
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep,make_lsq_spline, &
        bspline_interpolator_type,bspline_array_interpolator_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
    implicit none
    integer, parameter :: m=81
    type(bspline_interpolator_type) :: fit,least
    type(bspline_array_interpolator_type) :: curve
    real(dp) :: u(m),y(m),w(m),x(m,3),fp,target,step(m),nan
    real(dp), allocatable :: t(:),bad(:),p(:),v(:,:)
    integer :: i,k,r,code,checks,form,power
    checks=0
    u=[(real(i-1,dp)/real(m-1,dp),i=1,m)]
    w=0.8_dp+u
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    do k=1,5
        do r=2,k+1
            do form=1,2
                t=[(0.0_dp,i=1,k+1),(0.5_dp,i=1,r),(1.0_dp,i=1,k+1)]
                if (form==2) then
                    t(:k)=[(-0.1_dp*real(k+1-i,dp),i=1,k)]
                    t(size(t)-k+1:)=[(1.0_dp+0.1_dp*real(i,dp),i=1,k)]
                end if
                ! A lower-degree corner has no kth-derivative jump.
                power=k-r+1
                if (power==0) then
                    y=merge(1.0_dp,0.0_dp,u>=0.5_dp)
                else
                    y=max(u-0.5_dp,0.0_dp)**power
                end if
                fit=make_splrep(u,y,k=k,w=w,s=1e-6_dp,t=t,status=code)
                call require(code==0 .and. fit%valid(),'repeated-knot null-space scalar fit')
                call require(all(fit%knots()==t),'supplied repeated knots retained')
                call require(maxval(abs(fit%evaluate(u)-y))<2e-8_dp,'zero-penalty corner retained')
                call require(maxval(abs(fit%derivative([0.1_dp,0.9_dp],k)))<2e-7_dp, &
                    'zero kth derivatives on either piece')
                x(:,1)=y
                x(:,2)=2*y
                x(:,3)=-0.5_dp*y
                curve=make_splprep(x,u=u,k=k,w=w,s=1e-6_dp,t=t,parameters=p,status=code)
                call require(code==0 .and. curve%valid(),'shared repeated-knot null-space fit')
                call require(all(p==u) .and. all(curve%knots()==t),'shared knots and sites retained')
                v=curve%evaluate(u)
                call require(maxval(abs(v-x))<4e-8_dp,'shared zero-penalty reproduction')
                y=sin(8*u)+0.01_dp*cos(31*u)
                least=make_lsq_spline(u,y,t,k=k,w=w,status=code)
                call require(code==0,'repeated-knot unregularized LSQ')
                fp=sum((w*(least%evaluate(u)-y))**2)
                target=fp+0.01_dp*sum((w*y)**2)
                fit=make_splrep(u,y,k=k,w=w,s=target,t=t,status=code)
                call require(code==0 .and. fit%valid(),'feasible repeated-knot smoothing')
                call require(sum((w*(fit%evaluate(u)-y))**2)<=1.0011_dp*target, &
                    'repeated-knot residual bound')
            end do
        end do
        ! Discontinuous interior pieces may still have matching periodic endpoints.
        t=[(0.0_dp,i=1,k+1),(0.25_dp,i=1,k+1),(0.75_dp,i=1,k+1), &
            (1.0_dp,i=1,k+1)]
        step=merge(1.0_dp,0.0_dp,u>=0.25_dp .and. u<0.75_dp)
        y=step
        fit=make_splrep(u,y,k=k,w=w,s=1e-6_dp,t=t,bc_type='periodic',status=code)
        call require(code==0 .and. fit%valid(),'periodic discontinuous zero-penalty fit')
        call require(maxval(abs(fit%evaluate(u)-y))<2e-8_dp,'periodic step and right-continuity')
        do r=0,k-1
            call require(abs(fit%derivative(0.0_dp,r)-fit%derivative(1.0_dp,r)) < &
                1e-7_dp,'periodic endpoint constraints with repeated interiors')
        end do
        x(:,1)=y
        x(:,2)=2*y
        x(:,3)=-0.5_dp*y
        curve=make_splprep(x,u=u,k=k,w=w,s=1e-6_dp,t=t,bc_type='periodic', &
            status=code,parameters=p)
        call require(code==0 .and. curve%valid(),'periodic shared discontinuous fit')
        v=curve%evaluate(u)
        call require(maxval(abs(v-x))<4e-8_dp,'periodic shared step recovery')
        y(10)=y(10)+1.0_dp
        fit=make_splrep(u,y,k=k,s=1e-14_dp,t=t,status=code)
        call require(code==3 .and. .not.fit%valid(),'infeasible repeated-knot target')
        x(:,1)=y
        curve=make_splprep(x,u=u,k=k,s=1e-14_dp,t=t,status=code,parameters=p)
        call require(code==3 .and. .not.curve%valid() .and. .not.allocated(p), &
            'infeasible shared target cleans parameters')
        y=step
        bad=[(0.0_dp,i=1,k+1),(0.5_dp,i=1,k+2),(1.0_dp,i=1,k+1)]
        fit=make_splrep(u,y,k=k,s=1.0_dp,t=bad,status=code)
        call require(code==1 .and. .not.fit%valid(),'multiplicity exceeding k+1 rejected')
        curve=make_splprep(x,u=u,k=k,s=1.0_dp,t=bad,status=code,parameters=p)
        call require(code==1 .and. .not.curve%valid() .and. .not.allocated(p), &
            'invalid shared multiplicity cleans parameters')
        bad=t
        bad(k+2)=nan
        fit=make_splrep(u,y,k=k,s=1.0_dp,t=bad,status=code)
        call require(code==1 .and. .not.fit%valid(),'nonfinite repeated knots rejected')
        bad=[(0.0_dp,i=1,k+2),0.5_dp,(1.0_dp,i=1,k+1)]
        fit=make_splrep(u,y,k=k,s=1.0_dp,t=bad,status=code)
        call require(code==1 .and. .not.fit%valid(),'extra base endpoint knot rejected')
    end do
    write(*,'(a,i0,a)') 'Repeated-knot smoothing: ',checks,' checks passed'
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
end program test_interpolate_repeated_knots
