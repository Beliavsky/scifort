! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_sphere
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, &
        ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate, only : lsq_sphere_bivariate_spline, sphere_bivariate_spline_type, &
        sphere_spline_success, sphere_spline_invalid_input, sphere_spline_rank_deficient, &
        sphere_spline_numerical_failure, sphere_spline_size_limit
    implicit none
    type(sphere_bivariate_spline_type) :: fit, other, bad, invalid
    real(dp) :: theta(108), phi(108), r(108), w(108), tt(1), tp(3), pi, u, q, ref, nan, inf
    real(dp), allocatable :: values(:), grid(:,:), tx(:), ty(:), c(:)
    integer :: i, j, k, status, checks
    checks = 0
    pi = scifort_pi
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    tt = [pi/2.0_dp]
    tp = [pi/2.0_dp,pi,1.5_dp*pi]
    k = 0
    do j = 1, 12
        do i = 1, 9
            k = k+1
            theta(k) = real(i-1,dp)*pi/8.0_dp
            phi(k) = real(j-1,dp)*2.0_dp*pi/11.0_dp
            u = theta(k)/pi
            r(k) = 1.0_dp+2.0_dp*(3.0_dp*u*u-2.0_dp*u*u*u)
            w(k) = 0.5_dp+real(mod(k,7),dp)/5.0_dp
        end do
    end do
    fit = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w,status=status)
    call require(status == sphere_spline_success .and. fit%valid(),'valid full-rank fit')
    call require(fit%rank() == 10 .and. fit%degrees_of_freedom() == 10,'dimension and rank')
    call require(fit%residual() < 1e-23_dp,'analytic residual')
    values = fit%values(theta,phi,status=status)
    call require(status == 0 .and. maxval(abs(values-r)) < 2e-12_dp,'analytic paired values')
    q = 0.73_dp
    u = q/pi
    call near(fit%value(q,0.8_dp),1.0_dp+2.0_dp*(3*u*u-2*u*u*u),2e-12_dp,'cubic value')
    call near(fit%value(q,0.8_dp,dtheta=1),12*u*(1-u)/pi,2e-12_dp,'theta derivative')
    call near(fit%value(q,0.8_dp,dtheta=2),(12-24*u)/(pi*pi),3e-12_dp,'theta second')
    call near(fit%value(q,0.8_dp,dphi=1),0.0_dp,2e-12_dp,'longitude derivative')
    call near(fit%value(q,0.8_dp,dtheta=1,dphi=1),0.0_dp,2e-12_dp,'mixed derivative')
    do i = -3, 3
        call near(fit%value(q,0.8_dp+real(i,dp)*2*pi),fit%value(q,0.8_dp), &
            3e-13_dp,'periodic wrapping')
        call near(fit%value(0.0_dp,real(i,dp)),1.0_dp,2e-12_dp,'north pole')
        call near(fit%value(pi,real(i,dp)),3.0_dp,2e-12_dp,'south pole')
    end do
    grid = fit%grid([pi,q,0.0_dp],[1.0_dp,0.0_dp,2*pi],status=status)
    call require(status == 0 .and. all(shape(grid) == [3,3]),'grid shape and status')
    do i = 1, 3
        call near(grid(1,i),3.0_dp,2e-12_dp,'grid south pole')
    end do
    call fit%get_knots(tx,ty)
    c = fit%coefficients()
    call require(size(tx) == 9 .and. size(ty) == 11 .and. size(c) == 35,'TCK dimensions')
    call near(sum(c(:7))/7.0_dp,1.0_dp,2e-12_dp,'north coefficient row')
    other = lsq_sphere_bivariate_spline(theta(108:1:-1),phi(108:1:-1),r(108:1:-1),tt,tp, &
        w(108:1:-1),status=status)
    call require(status == 0,'unsorted reversed data')
    call near(other%value(q,0.4_dp),fit%value(q,0.4_dp),3e-12_dp,'order invariance')
    other = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,status=status)
    call require(status == 0,'default weights')
    call near(other%value(q,0.4_dp),fit%value(q,0.4_dp),3e-12_dp,'unit weights exact model')
    other = lsq_sphere_bivariate_spline(theta,phi,r*1e200_dp,tt,tp,w*1e200_dp,status=status)
    call require(status == 0 .and. other%valid(),'large finite data and weights')
    call near(other%value(q,0.4_dp)/1e200_dp,fit%value(q,0.4_dp),3e-12_dp,'large scaling')
    call require(.not. ieee_is_nan(other%residual()),'residual overflow is deliberate')
    other = lsq_sphere_bivariate_spline(theta,phi,r*1e-200_dp,tt,tp,w*1e-200_dp,status=status)
    call require(status == 0,'small finite data and weights')
    call near(other%value(q,0.4_dp)/1e-200_dp,fit%value(q,0.4_dp),3e-12_dp,'small scaling')
    other = lsq_sphere_bivariate_spline(theta,phi,0*r,tt,tp,status=status)
    call require(status == 0 .and. other%residual() == 0.0_dp,'zero observations')
    call near(other%value(q,0.4_dp),0.0_dp,0.0_dp,'zero fit')
    other = lsq_sphere_bivariate_spline([0.0_dp,0.0_dp],[0.0_dp,1.0_dp], &
        [2.0_dp,4.0_dp],tt,tp,w=[1.0_dp,2.0_dp],status=status)
    call require(status == sphere_spline_rank_deficient .and. other%valid(),'rank warning valid')
    call require(other%rank() == 1,'rank one')
    call near(other%value(0.0_dp,2.0_dp),3.6_dp,2e-13_dp,'weighted duplicate pole mean')
    call near(other%residual(),3.2_dp,3e-13_dp,'weighted residual convention')

    bad = lsq_sphere_bivariate_spline(theta(:1),phi(:1),r(:1),tt,tp,status=status)
    call rejected(sphere_spline_invalid_input,'too few points')
    bad = lsq_sphere_bivariate_spline(theta,phi(:2),r,tt,tp,status=status)
    call rejected(sphere_spline_invalid_input,'mismatched coordinates')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w=w(:2),status=status)
    call rejected(sphere_spline_invalid_input,'mismatched weights')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w=0*w,status=status)
    call rejected(sphere_spline_invalid_input,'zero weights')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,eps=nan,status=status)
    call rejected(sphere_spline_invalid_input,'NaN epsilon')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,eps=0.0_dp,status=status)
    call rejected(sphere_spline_invalid_input,'invalid epsilon')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp(:2),status=status)
    call rejected(sphere_spline_invalid_input,'longitude knot minimum')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp(3:1:-1),status=status)
    call rejected(sphere_spline_invalid_input,'unordered knots')
    bad = lsq_sphere_bivariate_spline(theta,phi,r,[0.0_dp],tp,status=status)
    call rejected(sphere_spline_invalid_input,'endpoint interior knot')
    theta(1) = nan
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,status=status)
    call rejected(sphere_spline_invalid_input,'NaN sample')
    theta(1) = -1.0_dp
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,status=status)
    call rejected(sphere_spline_invalid_input,'out-of-domain sample')
    theta(1) = 0.0_dp
    bad = lsq_sphere_bivariate_spline(theta,phi,r, &
        [(real(i,dp)*pi/101.0_dp,i=1,100)],tp,status=status)
    call rejected(sphere_spline_size_limit,'basis dimension limit')
    w(1) = tiny(1.0_dp)
    w(2) = huge(1.0_dp)
    bad = lsq_sphere_bivariate_spline(theta,phi,r,tt,tp,w,status=status)
    call rejected(sphere_spline_numerical_failure,'unrepresentable relative weights')
    bad = lsq_sphere_bivariate_spline([pi/3.0_dp,pi/3.0_dp+0.001_dp], &
        [0.0_dp,0.0_dp],[-1e308_dp,1e308_dp],tt,tp,status=status)
    call rejected(sphere_spline_numerical_failure,'unrepresentable physical coefficients')

    ref = invalid%value(q,0.0_dp,status=status)
    call require(status == 1 .and. ieee_is_nan(ref),'invalid object query')
    ref = fit%value(-1.0_dp,0.0_dp,status=status)
    call require(status == 1 .and. ieee_is_nan(ref),'invalid latitude query')
    ref = fit%value(q,inf,status=status)
    call require(status == 1 .and. ieee_is_nan(ref),'infinite longitude query')
    ref = fit%value(q,0.0_dp,dtheta=3,status=status)
    call require(status == 1 .and. ieee_is_nan(ref),'invalid derivative query')
    values = fit%values([q,q],[0.0_dp],status=status)
    call require(status == 1 .and. all(ieee_is_nan(values)),'mismatched paired queries')
    values = fit%values([q,-1.0_dp],[0.0_dp,0.0_dp],status=status)
    call require(status == 1 .and. ieee_is_finite(values(1)) .and. &
        ieee_is_nan(values(2)),'per-query invalid handling')
    values = fit%values([real(dp)::],[real(dp)::],dphi=-1,status=status)
    call require(status == 1 .and. size(values) == 0,'empty query invalid option')
    values = fit%values([real(dp)::],[real(dp)::],status=status)
    call require(status == 0 .and. size(values) == 0,'empty query')
    call invalid%get_knots(tx,ty)
    c = invalid%coefficients()
    call require(size(tx) == 0 .and. size(ty) == 0 .and. size(c) == 0,'invalid metadata')
    call require(ieee_is_nan(invalid%residual()),'invalid residual')
    print '(a,i0,a)', 'Spherical spline: ',checks,' checks passed.'

contains
    subroutine near(got, expected, tolerance, label)
        implicit none
        real(dp), intent(in) :: got !! Computed value.
        real(dp), intent(in) :: expected !! Analytic or invariance reference.
        real(dp), intent(in) :: tolerance !! Absolute error bound.
        character(len=*), intent(in) :: label !! Failure description.
        call require(ieee_is_finite(got),label)
        call require(abs(got-expected) <= tolerance,label)
    end subroutine near

    subroutine rejected(expected_status, label)
        implicit none
        integer, intent(in) :: expected_status !! Required construction failure code.
        character(len=*), intent(in) :: label !! Failure description.
        call require(status == expected_status,label)
        call require(.not. bad%valid(),label)
    end subroutine rejected

    subroutine require(ok, label)
        implicit none
        logical, intent(in) :: ok !! Assertion that must hold.
        character(len=*), intent(in) :: label !! Failure description.
        if (.not. ok) then
            print '(a,a)', 'FAIL: ',label
            error stop 1
        end if
        checks = checks+1
    end subroutine require
end program test_interpolate_sphere
