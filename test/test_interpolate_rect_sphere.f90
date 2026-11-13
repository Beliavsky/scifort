! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_rect_sphere
    use, intrinsic :: ieee_arithmetic, only : ieee_value,ieee_quiet_nan, &
        ieee_positive_inf,ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate, only : rect_sphere_bivariate_spline,sphere_bivariate_spline_type, &
        sphere_spline_success,sphere_spline_invalid_input,sphere_spline_numerical_failure, &
        sphere_spline_size_limit
    implicit none
    type(sphere_bivariate_spline_type) :: fit,other,bad
    real(dp) :: u(6),v(8),r(6,8),pi,t,q,expected,nan,inf
    real(dp), allocatable :: values(:,:),tx(:),ty(:),c(:)
    integer :: i,j,status,checks
    checks = 0
    pi = scifort_pi
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    u = [(real(i,dp)*pi/7.0_dp,i=1,6)]
    v = [(-0.7_dp+real(i-1,dp)*2*pi/8.0_dp,i=1,8)]
    do j = 1, 8
        do i = 1, 6
            t = u(i)/pi
            r(i,j) = 1.0_dp+2.0_dp*(3*t*t-2*t*t*t)
        end do
    end do
    fit = rect_sphere_bivariate_spline(u,v,r,pole_continuity=[.true.,.true.], &
        pole_values=[1.0_dp,3.0_dp],pole_exact=[.true.,.true.], &
        pole_flat=[.true.,.true.],status=status)
    call require(status == sphere_spline_success .and. fit%valid(),'interpolant valid')
    call require(fit%rank() == -1,'banded rank diagnostic unavailable')
    call require(fit%residual() == 0.0_dp,'interpolation residual convention')
    values = fit%grid(u,v,status=status)
    call require(status == 0 .and. maxval(abs(values-r)) < 3e-13_dp,'grid interpolation')
    q = 0.91_dp
    t = q/pi
    expected = 1.0_dp+2.0_dp*(3*t*t-2*t*t*t)
    do i = -2, 2
        call near(fit%value(q,0.3_dp+real(i,dp)*2*pi),expected,3e-13_dp,'periodic value')
    end do
    call near(fit%value(q,0.3_dp,dtheta=1),12*t*(1-t)/pi,3e-13_dp,'theta derivative')
    call near(fit%value(q,0.3_dp,dtheta=2),(12-24*t)/(pi*pi),1e-12_dp,'second derivative')
    call near(fit%value(q,0.3_dp,dphi=1),0.0_dp,5e-13_dp,'longitude derivative')
    call near(fit%value(0.0_dp,0.1_dp),1.0_dp,3e-13_dp,'north value')
    call near(fit%value(pi,2.0_dp),3.0_dp,3e-13_dp,'south value')
    call near(fit%value(0.0_dp,2.0_dp,dtheta=1),0.0_dp,3e-13_dp,'north flat')
    call near(fit%value(pi,2.0_dp,dtheta=1),0.0_dp,3e-13_dp,'south flat')
    call fit%get_knots(tx,ty)
    c = fit%coefficients()
    call require(size(tx) == 14 .and. size(ty) == 15 .and. size(c) == 110,'knots/coefficients')
    call near(ty(4),v(1),0.0_dp,'longitude origin')
    call near(ty(size(ty)-3),v(1)+2*pi,0.0_dp,'longitude end')
    other = rect_sphere_bivariate_spline(u,v,r,pole_continuity=[.true.,.true.], &
        pole_values=[nan,3.0_dp],pole_flat=[.true.,.true.],status=status)
    call require(status == 0 .and. other%valid(),'one supplied pole')
    call near(other%value(q,0.3_dp),expected,3e-12_dp,'estimated opposite pole')
    other = rect_sphere_bivariate_spline(u,v,r*1e200_dp,pole_continuity=[.true.,.true.], &
        pole_values=[1e200_dp,3e200_dp],pole_flat=[.true.,.true.],status=status)
    call require(status == 0,'large data scaling')
    call near(other%value(q,0.3_dp)/1e200_dp,expected,3e-12_dp,'large values')
    other = rect_sphere_bivariate_spline(u,v,r*1e-200_dp,pole_continuity=[.true.,.true.], &
        pole_values=[1e-200_dp,3e-200_dp],pole_flat=[.true.,.true.],status=status)
    call require(status == 0,'small data scaling')
    call near(other%value(q,0.3_dp)/1e-200_dp,expected,3e-12_dp,'small values')
    other = rect_sphere_bivariate_spline(u,v,0*r,status=status)
    call require(status == 0,'zero grid default options')
    call near(other%value(q,0.3_dp),0.0_dp,0.0_dp,'zero interpolant')
    other = rect_sphere_bivariate_spline(u,v,0*r+2.5_dp)
    call near(other%value(q,0.3_dp),2.5_dp,3e-13_dp,'constant default options')
    call near(other%value(pi,0.3_dp),2.5_dp,3e-13_dp,'constant estimated pole')

    other = rect_sphere_bivariate_spline(u,v,r,s=1.0_dp,status=status)
    call require(status == 0,'positive smoothing accepted')
    call require(other%valid(),'positive smoothing valid')
    bad = rect_sphere_bivariate_spline(u,v,r,s=-1.0_dp,status=status)
    call rejected(sphere_spline_invalid_input,'negative smoothing')
    bad = rect_sphere_bivariate_spline(u,v,r,s=nan,status=status)
    call rejected(sphere_spline_invalid_input,'NaN smoothing')
    bad = rect_sphere_bivariate_spline(u(:3),v,r(:3,:),status=status)
    call rejected(sphere_spline_invalid_input,'minimum latitude count')
    bad = rect_sphere_bivariate_spline(u,v(:3),r(:,:3),status=status)
    call rejected(sphere_spline_invalid_input,'minimum longitude count')
    bad = rect_sphere_bivariate_spline(u,v,r(:5,:),status=status)
    call rejected(sphere_spline_invalid_input,'grid shape')
    bad = rect_sphere_bivariate_spline(u(6:1:-1),v,r,status=status)
    call rejected(sphere_spline_invalid_input,'descending grid')
    bad = rect_sphere_bivariate_spline(u,v,r,pole_flat=[.true.,.false.],status=status)
    call rejected(sphere_spline_invalid_input,'flat requires continuity')
    bad = rect_sphere_bivariate_spline(u,v,r,pole_continuity=[.true.],status=status)
    call rejected(sphere_spline_invalid_input,'option pair shape')
    bad = rect_sphere_bivariate_spline(u,v,r,pole_values=[inf,0.0_dp],status=status)
    call rejected(sphere_spline_invalid_input,'infinite pole')
    r(1,1) = nan
    bad = rect_sphere_bivariate_spline(u,v,r,status=status)
    call rejected(sphere_spline_invalid_input,'NaN grid')
    r(1,1) = 1.0_dp
    u(1) = 0.0_dp
    bad = rect_sphere_bivariate_spline(u,v,r,status=status)
    call rejected(sphere_spline_invalid_input,'latitude endpoint')
    u(1) = 1e-20_dp
    bad = rect_sphere_bivariate_spline(u,v,r,status=status)
    call rejected(sphere_spline_numerical_failure,'numerical separation guard')
    u(1) = pi/7
    v(8) = v(1)+2*pi
    bad = rect_sphere_bivariate_spline(u,v,r,status=status)
    call rejected(sphere_spline_invalid_input,'duplicated seam endpoint')
    block
        real(dp), allocatable :: bigu(:),bigv(:),bigr(:,:)
        allocate(bigu(800),bigv(800),bigr(800,800))
        bigu = [(real(i,dp)*pi/801.0_dp,i=1,800)]
        bigv = [(real(i-1,dp)*2*pi/800.0_dp,i=1,800)]
        bigr = 1.0_dp
        bad = rect_sphere_bivariate_spline(bigu,bigv,bigr,status=status)
        call rejected(sphere_spline_size_limit,'workspace limit')
    end block
    print '(a,i0,a)', 'Rectangular spherical spline: ',checks,' checks passed.'

contains
    subroutine near(got,ref,tol,label)
        implicit none
        real(dp), intent(in) :: got !! Computed result.
        real(dp), intent(in) :: ref !! Analytic reference.
        real(dp), intent(in) :: tol !! Absolute error bound.
        character(len=*), intent(in) :: label !! Failure description.
        call require(.not. ieee_is_nan(got),label)
        call require(abs(got-ref) <= tol,label)
    end subroutine near
    subroutine rejected(code,label)
        implicit none
        integer, intent(in) :: code !! Expected construction status.
        character(len=*), intent(in) :: label !! Failure description.
        call require(status == code,label)
        call require(.not. bad%valid(),label)
    end subroutine rejected
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
end program test_interpolate_rect_sphere
