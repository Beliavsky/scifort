! SPDX-License-Identifier: MIT
program test_fitpack_smoothing
    use scifort_kinds, only : dp
    use scifort_optional_fitpack, only : univariate_spline, fitpack_smooth_success
    use scifort_interpolate_fitpack_compat, only : fitpack_tck_type,splev,splint
    implicit none
    integer, parameter :: n=60
    real(dp) :: x(n),y(n),w(n),residual
    type(fitpack_tck_type) :: fit
    integer :: i,status
    do i=1,n
        x(i)=real(i-1,dp)/real(n-1,dp)*6.0_dp
        y(i)=sin(x(i))+0.05_dp*cos(11.0_dp*x(i))
        w(i)=1.0_dp+0.02_dp*real(i,dp)
    end do
    fit=univariate_spline(x,y,s=1.0_dp,w=w,k=3,status=status,residual=residual)
    if (status /= fitpack_smooth_success .or. .not. fit%valid()) error stop 'fitpack smoothing'
    if (abs(sum((w*(y-splev(x,fit)))**2)-residual) > 1.0e-7_dp) &
        error stop 'FITPACK weighted residual'
    if (residual > 1.01_dp) error stop 'residual target not reached'
    if (abs(splint(0.0_dp,6.0_dp,fit)) > 10.0_dp) error stop 'bad spline integral'
    fit=univariate_spline(x,y,s=0.0_dp,status=status)
    if (status /= fitpack_smooth_success) error stop 's=0 FITPACK fit'
    if (maxval(abs(splev(x,fit)-y)) > 1.0e-7_dp) error stop 's=0 exact interpolation'
    print *, 'Optional FITPACK adaptive smoothing PASS'
end program test_fitpack_smoothing
