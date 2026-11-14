! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_make_periodic
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_splrep,make_splprep,make_interp_spline,bspline, &
        bspline_interpolator_type,bspline_array_interpolator_type,splrep,splprep,splev, &
        fitpack_tck_type,fitpack_parametric_tck_type
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: m=61
    real(dp) :: u(m),x(m,3),w(m),y(m),bad(m),badx(m,3),fp,nan,inf
    real(dp), allocatable :: p(:),t(:),v(:,:),expected(:,:),c(:,:)
    type(bspline_interpolator_type) :: fit,reference,raw
    type(bspline_array_interpolator_type) :: curve,raw_curve
    type(fitpack_tck_type) :: legacy
    type(fitpack_parametric_tck_type) :: path
    integer :: i,k,d,code,checks
    checks=0
    do i=1,m
        u(i)=real(i-1,dp)/real(m-1,dp)
        x(i,:)=[cos(2*acos(-1.0_dp)*u(i))+0.06_dp*sin(26*acos(-1.0_dp)*u(i)), &
            2*sin(2*acos(-1.0_dp)*u(i)),0.2_dp*cos(6*acos(-1.0_dp)*u(i))]
        w(i)=0.5_dp+u(i)
    end do
    x(m,:)=x(1,:)
    y=x(:,1)
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    do k=1,5
        fit=make_splrep(u,y,k=k,bc_type='periodic',status=code)
        call require(code==0 .and. fit%valid(),'scalar periodic interpolation')
        call require(maxval(abs(fit%evaluate(u)-y))<2e-12_dp,'scalar interpolated values')
        reference=make_interp_spline(u,y,k=k,bc_type='periodic')
        call require(maxval(abs(fit%evaluate(u)-reference%evaluate(u)))<2e-12_dp, &
            'scalar native periodic interpolator agreement')
        raw=bspline(fit%knots(),fit%coefficients(),k)
        call scalar_seam(raw,k)
        if (k>1) then
            call require(abs(fit%evaluate(1.231_dp)-fit%evaluate(0.231_dp))<2e-12_dp, &
                'interpolation wraps queries')
        end if
        fit=make_splrep(u,y,w=w,k=k,s=0.03_dp,bc_type='periodic',status=code)
        call require(code==0 .and. fit%valid(),'scalar periodic smoothing')
        fp=sum((w(:m-1)*(fit%evaluate(u(:m-1))-y(:m-1)))**2)
        call require(abs(fp-0.03_dp)<0.000033_dp,'scalar residual target')
        call scalar_seam(fit,k)
        legacy=splrep(u,y,w=w,k=k,s=0.03_dp,per=.true.)
        call require(maxval(abs(fit%evaluate(u)-splev(u,legacy)))<2e-12_dp, &
            'scalar legacy automatic agreement')
        raw=bspline(fit%knots(),fit%coefficients(),k)
        call require(abs(fit%evaluate(1.231_dp)-raw%evaluate(1.231_dp))<2e-12_dp, &
            'smoothing retains polynomial extrapolation')
        fit=make_splrep(u,y,w=w,k=k,s=100.0_dp,nest=2*(k+1), &
            bc_type='periodic',status=code)
        call require(code==0 .and. fit%valid(),'scalar minimum capacity constant fit')
        call require(abs(fit%evaluate(0.231_dp)-sum(w(:m-1)**2*y(:m-1))/ &
            sum(w(:m-1)**2))<2e-13_dp,'scalar constant limit')
        fit=make_splrep(u,y,k=k,s=1e-9_dp,nest=2*(k+1),bc_type='periodic',status=code)
        call require(code==3 .and. .not.fit%valid(),'scalar insufficient capacity')
        curve=make_splprep(x,u=u,k=k,bc_type='periodic',status=code,parameters=p)
        call require(code==0 .and. curve%valid(),'parametric periodic interpolation')
        call require(all(p==u),'supplied periodic parameters retained')
        call require(maxval(abs(curve%evaluate(u)-x))<3e-12_dp,'closed interpolation values')
        raw_curve=bspline(curve%knots(),curve%coefficients(),k)
        do d=0,k-1
            v=raw_curve%derivative([0.0_dp,1.0_dp],d)
            call require(maxval(abs(v(1,:)-v(2,:)))<2e-8_dp*max(1.0_dp,maxval(abs(v))), &
                'unwrapped curve seam')
        end do
        curve=make_splprep(x,k=k,bc_type='periodic',status=code,parameters=p)
        call require(code==0 .and. curve%valid(),'automatic periodic interpolation sites')
        call require(maxval(abs(curve%evaluate(p)-x))<3e-12_dp,'chord periodic interpolation')
        curve=make_splprep(x,u=u,w=w,k=k,s=0.05_dp,bc_type='periodic', &
            status=code,parameters=p)
        call require(code==0 .and. curve%valid(),'parametric periodic smoothing')
        v=curve%evaluate(u)
        fp=sum((spread(w(:m-1),2,3)*(v(:m-1,:)-x(:m-1,:)))**2)
        call require(abs(fp-0.05_dp)<0.000055_dp,'shared periodic residual')
        path=splprep(x,u=u,w=w,k=k,s=0.05_dp,per=.true.)
        expected=splev(u,path)
        call require(maxval(abs(v-expected))<3e-12_dp,'parametric legacy agreement')
        t=curve%knots()
        c=curve%coefficients()
        call require(size(c,1)==size(t)-k-1 .and. size(c,2)==3,'unpadded periodic coefficients')
        curve=make_splprep(x,u=u,k=k,s=1e-9_dp,nest=2*(k+1),bc_type='periodic', &
            status=code,parameters=p)
        call require(code==3 .and. .not.curve%valid(),'parametric insufficient capacity')
        call require(.not.allocated(p),'no periodic parameters on failure')
        curve=make_splprep(x,u=u,w=w,k=k,s=1000.0_dp,nest=2*(k+1), &
            bc_type='periodic',status=code)
        call require(code==0 .and. curve%valid(),'constant-vector minimum capacity')
        v=curve%evaluate([0.0_dp,0.4_dp,1.0_dp])
        do i=1,3
            call require(maxval(abs(v(:,i)-sum(w(:m-1)**2*x(:m-1,i))/ &
                sum(w(:m-1)**2)))<3e-13_dp,'constant-vector weighted mean')
        end do
        curve=make_splprep(x,u=u,w=w,k=k,s=0.05_dp,nest=m+2*k, &
            bc_type='periodic',status=code)
        call require(code==0 .and. curve%valid(),'full periodic capacity')
        curve=make_splprep(x,u=u,s=0.05_dp,k=k,t=t,bc_type='periodic',status=code)
        call require(code==2 .and. .not.curve%valid(),'periodic supplied knots unsupported')
    end do
    bad=y
    bad(m)=bad(1)+0.01_dp
    fit=make_splrep(u,bad,s=0.03_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.fit%valid(),'scalar endpoint mismatch')
    badx=x
    badx(m,2)=0.01_dp
    curve=make_splprep(badx,s=0.03_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.curve%valid(),'coordinate endpoint mismatch')
    bad(m)=nan
    fit=make_splrep(u,bad,s=0.03_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.fit%valid(),'nonfinite closing ordinate rejected')
    badx(m,2)=inf
    curve=make_splprep(badx,s=0.03_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.curve%valid(),'nonfinite closing coordinate rejected')
    bad=w
    bad(m)=nan
    fit=make_splrep(u,y,w=bad,s=0.03_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.fit%valid(),'invalid closing weight rejected')
    curve=make_splprep(x,w=bad,s=0.05_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.curve%valid(),'invalid closing curve weight rejected')
    fit=make_splrep(u,y,w=w,bc_type='periodic',status=code)
    call require(code==1 .and. .not.fit%valid(),'interpolation weights rejected')
    fit=make_splrep(u,y,s=0.03_dp,xb=-0.1_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.fit%valid(),'periodic scalar bounds')
    curve=make_splprep(x,u=u,s=0.05_dp,ue=1.1_dp,bc_type='periodic',status=code)
    call require(code==1 .and. .not.curve%valid(),'periodic parameter bounds')
    fit=make_splrep(u,y,s=0.03_dp,nest=7,bc_type='periodic',status=code)
    call require(code==1 .and. .not.fit%valid(),'small periodic capacity')
    curve=make_splprep(x,s=0.05_dp,nest=m+7,bc_type='periodic',status=code)
    call require(code==1 .and. .not.curve%valid(),'excess periodic capacity')
    fit=make_splrep(u,y,s=0.03_dp,t=t,bc_type='periodic',status=code)
    call require(code==2 .and. .not.fit%valid(),'scalar periodic fixed knots unsupported')
    fit=make_splrep(u*1e-310_dp,y,bc_type='periodic',status=code)
    call require(code==3 .and. .not.fit%valid(),'periodic subnormal span')
    curve=make_splprep(x,u=u*1e-310_dp,bc_type='periodic',status=code)
    call require(code==3 .and. .not.curve%valid(),'parametric subnormal span')
    curve=make_splprep(x*1e150_dp,u=u,w=w*1e-150_dp,s=0.05_dp, &
        bc_type='periodic',status=code)
    call require(code==0 .and. curve%valid(),'periodic common scaling')
    fit=make_splrep(u,y,w=w,s=0.03_dp,bc_type='periodic',status=code)
    bad=w
    bad(m)=1e200_dp
    reference=make_splrep(u,y,w=bad,s=0.03_dp,bc_type='periodic',status=code)
    call require(code==0 .and. reference%valid(),'finite closing weight accepted')
    call require(maxval(abs(fit%evaluate(u)-reference%evaluate(u)))<2e-13_dp, &
        'closing weight excluded from smoothing target')
    fit=make_splrep(u*1e-100_dp,y,bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'normalized tiny-interval periodic interpolation')
    call require(maxval(abs(fit%evaluate(u*1e-100_dp)-y))<3e-12_dp, &
        'tiny-interval interpolation values')
    fit=make_splrep(u*1e100_dp,y,bc_type='periodic',status=code)
    call require(code==0 .and. fit%valid(),'normalized huge-interval periodic interpolation')
    block
        real(dp) :: high(m,12)
        do i=1,12
            high(:,i)=x(:,mod(i,3)+1)*real(i,dp)
        end do
        curve=make_splprep(high,u=u,s=3.0_dp,bc_type='periodic',status=code)
        call require(code==0 .and. curve%channels()==12,'twelve periodic channels')
        curve=make_splprep(high,u=u,bc_type='periodic',status=code)
        call require(code==0 .and. curve%channels()==12,'twelve interpolation channels')
    end block
    print '(a,i0,a)','Periodic native constructors: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Required assertion.
        character(*), intent(in) :: label !! Failure description.
        checks=checks+1
        if (.not.ok) then
            print '(a)','FAIL: '//label
            error stop 1
        end if
    end subroutine require
    subroutine scalar_seam(spline,degree)
        implicit none
        type(bspline_interpolator_type), intent(in) :: spline !! Unwrapped periodic fit.
        integer, intent(in) :: degree !! Spline degree.
        real(dp) :: a,b
        integer :: order
        do order=0,degree-1
            a=spline%derivative(0.0_dp,order)
            b=spline%derivative(1.0_dp,order)
            call require(abs(a-b)<2e-8_dp*max(1.0_dp,abs(a),abs(b)),'scalar seam')
        end do
    end subroutine scalar_seam
end program test_interpolate_make_periodic
