! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_periodic_fixed
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep, &
        bspline_interpolator_type,bspline_array_interpolator_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: m=81
    real(dp) :: u(m),y(m),x(m,12),w(m),fp,constant_fp,mean,delta,nan,inf
    real(dp), allocatable :: t(:),bad(:),v(:,:),p(:)
    type(bspline_interpolator_type) :: fit
    type(bspline_array_interpolator_type) :: curve
    integer :: i,j,k,ni,code,checks,d
    checks=0
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    do i=1,m
        u(i)=real(i-1,dp)/real(m-1,dp)
        y(i)=cos(2*acos(-1.0_dp)*u(i))+0.15_dp*sin(6*acos(-1.0_dp)*u(i))
        w(i)=0.7_dp+u(i)
    end do
    y(m)=y(1)
    do j=1,12
        x(:,j)=y*real(j,dp)/12
    end do
    mean=sum(w(:m-1)**2*y(:m-1))/sum(w(:m-1)**2)
    constant_fp=sum((w(:m-1)*(mean-y(:m-1)))**2)
    do k=1,5
        do ni=1,12
            if (allocated(t)) deallocate(t)
            allocate(t(ni+2*k+1))
            do i=1,size(t)
                t(i)=real(i-k-1,dp)/real(ni,dp)
            end do
            fit=make_splrep(u,y,w=w,k=k,s=2*constant_fp,t=t, &
                bc_type='periodic',status=code)
            call require(code==0 .and. fit%valid(),'constant scalar fixed knots')
            call require(all(fit%knots()==t),'retain scalar full knots')
            call require(maxval(abs(fit%evaluate(u)-mean))<2e-13_dp,'weighted constant')
            curve=make_splprep(x,u=u,w=w,k=k,s=20*constant_fp,t=t, &
                bc_type='periodic',status=code,parameters=p)
            call require(code==0 .and. curve%valid(),'constant twelve-channel fixed fit')
            call require(all(curve%knots()==t) .and. all(p==u),'retain curve knots and sites')
            v=curve%evaluate(u)
            do j=1,12
                call require(maxval(abs(v(:,j)-mean*real(j,dp)/12))<3e-13_dp, &
                    'weighted constant coordinate')
            end do
            if (ni<4) then
                fit=make_splrep(u,y,w=w,k=k,s=1e-8_dp,t=t, &
                    bc_type='periodic',status=code)
                call require(code==3 .and. .not.fit%valid(),'infeasible fixed residual')
            end if
            if (ni<4) cycle
            delta=0.7_dp*constant_fp
            fit=make_splrep(u,y,w=w,k=k,s=delta,t=t,bc_type='periodic',status=code)
            call require(code==0 .and. fit%valid(),'active scalar supplied knots')
            call require(all(fit%knots()==t),'active scalar knots unchanged')
            fp=sum((w(:m-1)*(fit%evaluate(u(:m-1))-y(:m-1)))**2)
            call require(abs(fp-delta)<0.0011_dp*delta,'scalar fixed residual target')
            do d=0,k-1
                call require(abs(fit%derivative(0.0_dp,d)-fit%derivative(1.0_dp,d)) < &
                    1e-8_dp*max(1.0_dp,abs(fit%derivative(0.0_dp,d))), &
                    'unwrapped fixed scalar seam')
            end do
            curve=make_splprep(x,u=u,w=w,k=k,s=delta,t=t, &
                bc_type='periodic',status=code)
            call require(code==0 .and. curve%valid(),'active shared fixed residual')
            v=curve%evaluate(u)
            fp=sum((spread(w(:m-1),2,12)*(v(:m-1,:)-x(:m-1,:)))**2)
            call require(abs(fp-delta)<0.0011_dp*delta,'shared fixed residual target')
            call require(all(curve%knots()==t),'active curve knots unchanged')
        end do
        fit=make_splrep(u,y,k=k,s=1.0_dp,t=t,nest=m+2*k, &
            bc_type='periodic',status=code)
        call require(code==1 .and. .not.fit%valid(),'t and nest mutually exclusive')
        fit=make_splrep(u,y,k=k,t=t,bc_type='periodic',status=code)
        call require(code==1 .and. .not.fit%valid(),'s zero rejects knots')
        bad=t
        bad(1)=bad(2)+0.01_dp
        call invalid(bad,k,'unordered exterior')
        bad=t
        bad(k+1)=0.001_dp
        call invalid(bad,k,'incorrect endpoint')
        bad=t
        bad(k+2)=bad(k+1)
        call invalid(bad,k,'repeated knot')
        bad=t
        bad(1)=nan
        call invalid(bad,k,'nan knot')
        bad(1)=inf
        call invalid(bad,k,'infinite knot')
        call invalid(t(:2*k+1),k,'too few knots')
    end do
    ! Linear cyclic basis with no observations in two narrow knot spans.
    t=[-0.998_dp,0.0_dp,0.001_dp,0.002_dp,1.0_dp,1.001_dp]
    fit=make_splrep(u,y,k=1,s=1.0_dp,t=t,bc_type='periodic',status=code)
    call require(code==3 .and. .not.fit%valid(),'unattainable rank-deficient scalar target')
    curve=make_splprep(x,u=u,k=1,s=1.0_dp,t=t,bc_type='periodic',status=code)
    call require(code==3 .and. .not.curve%valid(),'unattainable rank-deficient shared target')
    t=[(real(j-4,dp)/5,j=1,12)]
    fit=make_splrep(u,y,w=w,k=3,s=0.7_dp*constant_fp,t=t, &
        bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'base fixed fit for scaling')
    curve=make_splprep(x*1e150_dp,u=u,w=w*1e-150_dp,k=3,s=constant_fp,t=t, &
        bc_type='periodic',status=code)
    call require(code==0 .and. curve%valid(),'fixed common coordinate weight scaling')
    fit=make_splrep(u*1e100_dp,y,w=w,k=3,s=0.7_dp*constant_fp,t=t*1e100_dp, &
        bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'large physical fixed knot scale')
    call require(all(fit%knots()==t*1e100_dp),'large supplied knots retained')
    fit=make_splrep(u*1e-100_dp,y,w=w,k=3,s=0.7_dp*constant_fp,t=t*1e-100_dp, &
        bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'small physical fixed knot scale')
    call require(all(fit%knots()==t*1e-100_dp),'small supplied knots retained')
    print '(a,i0,a)','Periodic supplied knots: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Required assertion.
        character(*), intent(in) :: label !! Failure description.
        checks=checks+1
        if (.not.ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine require
    subroutine invalid(knots,degree,label)
        implicit none
        real(dp), intent(in) :: knots(:) !! Invalid full periodic knot vector.
        integer, intent(in) :: degree !! Degree from one through five.
        character(*), intent(in) :: label !! Failure description.
        fit=make_splrep(u,y,k=degree,s=1.0_dp,t=knots,bc_type='periodic',status=code)
        call require(code==1 .and. .not.fit%valid(),label)
        curve=make_splprep(x,u=u,k=degree,s=1.0_dp,t=knots, &
            bc_type='periodic',status=code,parameters=p)
        call require(code==1 .and. .not.curve%valid(),label//' curve')
        call require(.not.allocated(p),'failure leaves no parameters')
    end subroutine invalid
end program test_interpolate_periodic_fixed
