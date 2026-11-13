! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_smooth_sphere
    use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan, &
        ieee_positive_inf,ieee_is_finite,ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_constants, only : pi => scifort_pi
    use scifort_interpolate, only : smooth_sphere_bivariate_spline, &
        lsq_sphere_bivariate_spline,sphere_bivariate_spline_type, &
        sphere_spline_invalid_input,sphere_spline_numerical_failure,sphere_spline_size_limit
    implicit none
    type(sphere_bivariate_spline_type) :: fit,other,bad,lsq
    real(dp) :: t(120),p(120),r(120),w(120),cubic(120),targets(3),fp,base,scale,nan,inf
    real(dp), allocatable :: tt(:),tp(:),values(:),grid(:,:)
    integer :: i,k,status,checks
    checks = 0
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    do i = 1, size(t)
        t(i) = pi*modulo(real(i,dp)*0.618033988749895_dp,1.0_dp)
        p(i) = 2*pi*modulo(real(i,dp)*0.414213562373095_dp,1.0_dp)
        cubic(i) = 1+2*(t(i)/pi)**2*(3-2*t(i)/pi)
        r(i) = cos(t(i))+0.7_dp*sin(t(i))*cos(p(i)) &
            +0.04_dp*sin(17*t(i)+5*p(i))
        w(i) = 1+0.4_dp*sin(real(i,dp))
    end do
    fit = smooth_sphere_bivariate_spline(t,p,cubic,s=1e-12_dp,status=status)
    call require(status == 0 .and. fit%valid(),'exact cubic fit')
    call near(fit%value(0.8_dp,0.3_dp),1+2*(0.8_dp/pi)**2*(3-1.6_dp/pi), &
        3e-12_dp,'exact cubic value')
    call near(fit%value(0.0_dp,0.3_dp),1.0_dp,3e-12_dp,'north pole')
    call near(fit%value(pi,0.3_dp),3.0_dp,3e-12_dp,'south pole')
    call require(fit%rank() == -1 .and. fit%degrees_of_freedom() == 2,'polynomial metadata')
    call near(fit%value(0.0_dp,0.3_dp,dtheta=1),0.0_dp,3e-12_dp,'C1 pole derivative')
    targets = [0.1_dp,1.0_dp,1000.0_dp]
    do k = 1, size(targets)
        fit = smooth_sphere_bivariate_spline(t,p,r,w=w,s=targets(k),status=status)
        call require(status == 0 .and. fit%valid(),'weighted smoothing')
        values = fit%values(t,p,status=status)
        call require(status == 0,'paired evaluation')
        fp = sum((w*(r-values))**2)
        call near(fit%residual(),fp,3e-11_dp*(1+fp),'weighted physical residual')
        call require(fp <= 1.0011_dp*targets(k),'residual budget')
        if (k < 3) call near(fp,targets(k),0.0011_dp*targets(k),'active target')
        call near(fit%value(0.8_dp,0.3_dp+4*pi),fit%value(0.8_dp,0.3_dp), &
            3e-12_dp,'periodic query')
        call require(ieee_is_finite(fit%value(0.8_dp,0.3_dp,dtheta=2,dphi=1)), &
            'finite derivative')
        grid = fit%grid([0.5_dp,1.5_dp],[0.2_dp,0.7_dp,1.3_dp])
        call near(grid(2,3),fit%value(1.5_dp,1.3_dp),0.0_dp,'grid agreement')
        call fit%get_knots(tt,tp)
        if (k < 3) then
            call require(fit%degrees_of_freedom() == 6+(size(tt)-8)*(size(tp)-7), &
                'adaptive model dimension')
        end if
    end do
    fit = smooth_sphere_bivariate_spline(t,p,r,w=w,s=0.1_dp,status=status)
    base = fit%value(0.8_dp,0.3_dp)
    fp = fit%residual()
    do k = 1, 2
        scale = 1e150_dp
        if (k == 2) scale = 1e-150_dp
        other = smooth_sphere_bivariate_spline(t,p,r*scale,w=w, &
            s=0.1_dp*scale*scale,status=status)
        call require(status == 0 .and. other%valid(),'scaled observations')
        call near(other%value(0.8_dp,0.3_dp)/scale,base,3e-9_dp,'scaled value')
        call near((other%residual()/scale)/scale,fp,3e-9_dp,'scaled residual')
        other = smooth_sphere_bivariate_spline(t,p,r,w=w*scale, &
            s=0.1_dp*scale*scale,status=status)
        call require(status == 0 .and. other%valid(),'scaled residual weights')
        call near(other%value(0.8_dp,0.3_dp),base,3e-9_dp,'weight scale invariance')
    end do
    other = smooth_sphere_bivariate_spline(t,p,r,s=huge(1.0_dp),status=status)
    call require(status == 0 .and. other%valid(),'huge smoothing target')
    call require(other%degrees_of_freedom() == 2,'huge-target polynomial limit')
    other = smooth_sphere_bivariate_spline(t,p,0*r,status=status)
    call require((status == 0 .or. status == 2) .and. other%valid(),'zero interpolation')
    call near(other%value(0.8_dp,0.3_dp),0.0_dp,0.0_dp,'zero value')

    ! Compare the banded minimum-norm path with the existing dense SVD.
    block
        real(dp) :: rt(40),rp(40),rr(40),qt(3),qp(3)
        real(dp), allocatable :: expected(:),got(:)
        rt = pi/2
        rp = [(2*pi*real(i-1,dp)/40.0_dp,i=1,40)]
        rr = 2+sin(rp)
        lsq = lsq_sphere_bivariate_spline(rt,rp,rr,[pi/2], &
            pi*[0.5_dp,1.0_dp,1.5_dp],status=status)
        call require(status == 2 .and. lsq%rank() == 4,'rank-deficient SVD reference')
        other = smooth_sphere_bivariate_spline(rt,rp,rr, &
            s=lsq%residual()*1.00001_dp,status=status)
        call require(status == 2 .and. other%valid(),'banded rank warning retains fit')
        call require(other%rank() == 4 .and. other%degrees_of_freedom() == 10, &
            'banded rank metadata')
        qt = [0.3_dp,pi/2,2.6_dp]
        qp = [0.2_dp,0.7_dp,4.1_dp]
        expected = lsq%values(qt,qp)
        got = other%values(qt,qp)
        do i = 1, 3
            call near(got(i),expected(i),3e-10_dp,'banded/SVD minimum-norm agreement')
        end do
        call near(other%residual(),lsq%residual(),3e-11_dp,'rank-reduced residual')
    end block

    bad = smooth_sphere_bivariate_spline(t(:1),p(:1),r(:1),status=status)
    call rejected(sphere_spline_invalid_input,'minimum observation count')
    bad = smooth_sphere_bivariate_spline(t,p(:119),r,status=status)
    call rejected(sphere_spline_invalid_input,'paired shapes')
    bad = smooth_sphere_bivariate_spline(t,p,r,w=w(:119),status=status)
    call rejected(sphere_spline_invalid_input,'weight shape')
    bad = smooth_sphere_bivariate_spline(t,p,r,s=-1.0_dp,status=status)
    call rejected(sphere_spline_invalid_input,'negative target')
    bad = smooth_sphere_bivariate_spline(t,p,r,s=nan,status=status)
    call rejected(sphere_spline_invalid_input,'NaN target')
    bad = smooth_sphere_bivariate_spline(t,p,r,s=inf,status=status)
    call rejected(sphere_spline_invalid_input,'infinite target')
    bad = smooth_sphere_bivariate_spline(t,p,r,eps=0.0_dp,status=status)
    call rejected(sphere_spline_invalid_input,'rank threshold endpoint')
    bad = smooth_sphere_bivariate_spline(t,p,r,eps=nan,status=status)
    call rejected(sphere_spline_invalid_input,'NaN rank threshold')
    bad = smooth_sphere_bivariate_spline(t,p,r,ntest=7,status=status)
    call rejected(sphere_spline_invalid_input,'invalid knot capacity')
    bad = smooth_sphere_bivariate_spline(t,p,r,s=0.1_dp,ntest=8,npest=8,status=status)
    call rejected(sphere_spline_size_limit,'insufficient knot capacity')
    bad = smooth_sphere_bivariate_spline(t,p,r,npest=huge(0),status=status)
    call rejected(sphere_spline_size_limit,'workspace integer guard')
    bad = smooth_sphere_bivariate_spline(t,p,r,ntest=100,npest=100,status=status)
    call rejected(sphere_spline_size_limit,'workspace budget')
    w(1) = 0
    bad = smooth_sphere_bivariate_spline(t,p,r,w=w,status=status)
    call rejected(sphere_spline_invalid_input,'zero weight')
    w(1) = inf
    bad = smooth_sphere_bivariate_spline(t,p,r,w=w,status=status)
    call rejected(sphere_spline_invalid_input,'infinite weight')
    w = 1e300_dp
    w(1) = 1e-300_dp
    bad = smooth_sphere_bivariate_spline(t,p,r,w=w,status=status)
    call rejected(sphere_spline_numerical_failure,'relative weight underflow')
    bad = smooth_sphere_bivariate_spline(t,p,r*1e200_dp,s=1e-200_dp,status=status)
    call rejected(sphere_spline_numerical_failure,'relative target underflow')
    bad = smooth_sphere_bivariate_spline([0.0_dp,0.0_dp],[0.0_dp,0.0_dp], &
        [1.0_dp,3.0_dp],s=0.01_dp,ntest=12,npest=12,status=status)
    call rejected(sphere_spline_numerical_failure,'conflicting repeated observations')
    r(1) = nan
    bad = smooth_sphere_bivariate_spline(t,p,r,status=status)
    call rejected(sphere_spline_invalid_input,'NaN observation')
    r(1) = 1
    t(1) = -1
    bad = smooth_sphere_bivariate_spline(t,p,r,status=status)
    call rejected(sphere_spline_invalid_input,'colatitude outside support')
    t(1) = 0
    p(1) = 2*pi+0.1_dp
    bad = smooth_sphere_bivariate_spline(t,p,r,status=status)
    call rejected(sphere_spline_invalid_input,'training longitude outside support')
    call require(ieee_is_nan(bad%residual()),'invalid fit residual')
    print '(a,i0,a)', 'Scattered spherical smoothing: ',checks,' checks passed.'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
    subroutine near(got,ref,tol,label)
        implicit none
        real(dp), intent(in) :: got !! Computed result.
        real(dp), intent(in) :: ref !! Analytic or independent solver reference.
        real(dp), intent(in) :: tol !! Absolute error bound.
        character(len=*), intent(in) :: label !! Failure description.
        call require(ieee_is_finite(got),label)
        call require(abs(got-ref) <= tol,label)
    end subroutine near
    subroutine rejected(code,label)
        implicit none
        integer, intent(in) :: code !! Expected constructor status.
        character(len=*), intent(in) :: label !! Failure description.
        call require(status == code,label)
        call require(.not. bad%valid(),label)
    end subroutine rejected
end program test_interpolate_smooth_sphere
