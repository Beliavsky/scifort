! SPDX-License-Identifier: MIT
program test_fitpack_periodic_parametric
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use scifort_kinds, only: dp
    use scifort_interpolate_fitpack_compat, only: fitpack_tck_type, &
        fitpack_parametric_tck_type, splev
    use scifort_optional_fitpack, only: fitpack_smooth_success, &
        fitpack_smooth_invalid_input, periodic_univariate_spline, parametric_smoothing_spline
    implicit none
    integer, parameter :: n=51
    integer :: i, status
    real(dp) :: x(n), y(n), pts(n,2), r, pi, d, w(n)
    real(dp), allocatable :: values(:,:), fitted(:,:)
    type(fitpack_tck_type) :: t
    type(fitpack_parametric_tck_type) :: p
    pi=acos(-1.0_dp)
    do i=1,n
        x(i)=real(i-1,dp)/real(n-1,dp)
        y(i)=sin(2.0_dp*pi*x(i)) + 0.02_dp*sin(9.0_dp*pi*x(i))
        pts(i,1)=cos(2.0_dp*pi*x(i))
        pts(i,2)=sin(2.0_dp*pi*x(i))
        w(i)=1.0_dp+real(i,dp)/real(n,dp)
    end do
    y(n) = y(n) + 0.27_dp ! Periodic FITPACK ignores the final ordinate.
    t=periodic_univariate_spline(x,y,k=3,s=0.0_dp,status=status,residual=r)
    if (status/=fitpack_smooth_success .or. .not.t%valid()) error stop 'periodic s=0 status'
    if (abs(splev(x(1),t)-splev(x(n),t))>1.0e-10_dp) error stop 'periodic endpoint value'
    if (abs(splev(x(1),t,der=1)-splev(x(n),t,der=1))>1.0e-7_dp) &
        error stop 'periodic derivative endpoint'
    if (maxval(abs(splev(x(1:n-1),t)-y(1:n-1)))>1.0e-7_dp) error stop 'periodic interpolation'
    if (abs(r)>1.0e-8_dp) error stop 'periodic interpolation residual'
    t=periodic_univariate_spline(x,y,w=w,k=3,s=0.02_dp,status=status,residual=r)
    if (status/=fitpack_smooth_success .or. .not.t%valid()) error stop 'periodic smooth status'
    if (r>0.0201_dp) error stop 'periodic smoothing residual'
    p=parametric_smoothing_spline(pts,s=0.0_dp,k=3,status=status,residual=r)
    if (status/=fitpack_smooth_success .or. .not.p%valid()) error stop 'parametric s=0 status'
    values=splev(p%u,p)
    if (maxval(abs(values-pts))>1.0e-7_dp) error stop 'parametric s=0 interpolation'
    if (r>1.0e-8_dp) error stop 'parametric interpolation residual'
    p=parametric_smoothing_spline(pts,u=x,w=w,s=0.08_dp,k=3,status=status,residual=r)
    if (status/=fitpack_smooth_success .or. .not.p%valid()) error stop 'parametric smoothing status'
    fitted=splev(x,p)
    d=sum((spread(w,2,2)*(pts-fitted))**2)
    if (abs(d-r)>1.0e-7_dp) error stop 'parametric weighted residual'
    if (r>0.0801_dp) error stop 'parametric target'
    pts(n,1) = pts(n,1) + 0.17_dp ! Last closed point is overridden.
    p=parametric_smoothing_spline(pts,u=x,s=0.01_dp,k=3,per=.true.,status=status,residual=r)
    if (status/=fitpack_smooth_success .or. .not.p%valid()) error stop 'closed status'
    if (maxval(abs(splev(x(1),p)-splev(x(n),p)))>1.0e-8_dp) error stop 'closed endpoints'
    if (maxval(abs(splev(x(1),p,der=1)-splev(x(n),p,der=1)))>1.0e-7_dp) &
        error stop 'closed derivative endpoint'
    if (.not.all(ieee_is_finite(splev(x,t)))) error stop 'nonfinite periodic values'
    t=periodic_univariate_spline(x,y,k=6,status=status)
    if (status/=fitpack_smooth_invalid_input .or. t%valid()) error stop 'invalid spline degree'
    p=parametric_smoothing_spline(pts,u=x(1:n-1),status=status)
    if (status/=fitpack_smooth_invalid_input .or. p%valid()) error stop 'invalid parameter size'
    print *, 'Optional FITPACK periodic and parametric PASS'
end program test_fitpack_periodic_parametric
