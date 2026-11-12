! SPDX-License-Identifier: MIT
! Smoke test re-export from the public aggregate interpolation module.
program test_interpolate_make_lsq_public
    use scifort_kinds, only: dp
    use scifort_interpolate, only: make_lsq_spline, bspline_interpolator_type, &
        bspline_array_interpolator_type, bspline_status_success
    implicit none
    real(dp) :: x(12), y(12), yy(12,2), t(10)
    type(bspline_interpolator_type) :: single
    type(bspline_array_interpolator_type) :: multi
    real(dp), allocatable :: val(:)
    integer :: i, code

    do i=1,12
        x(i) = real(i-1,dp)/11.0_dp
    end do
    y = 1.0_dp + 2.0_dp*x - 0.5_dp*x*x
    yy(:,1) = y
    yy(:,2) = 2.0_dp*y
    t = [spread(x(1),1,4),0.25_dp,0.75_dp,spread(x(12),1,4)]
    single = make_lsq_spline(x,y,t,status=code)
    if (code /= bspline_status_success) error stop 'public make_lsq_spline scalar'
    multi = make_lsq_spline(x,yy,t,status=code)
    if (code /= bspline_status_success) error stop 'public make_lsq_spline channels'
    val = multi%evaluate(0.37_dp)
    if (abs(single%evaluate(0.37_dp)-(1.0_dp+2.0_dp*0.37_dp-0.5_dp*0.37_dp**2)) > 1.0e-11_dp) &
        error stop 'public make_lsq_spline polynomial fit'
    if (maxval(abs(val-[single%evaluate(0.37_dp),2.0_dp*single%evaluate(0.37_dp)])) > 2.0e-11_dp) &
        error stop 'public make_lsq_spline channels'
    print '(a)', 'PASS: public make_lsq_spline scalar and channels'
end program test_interpolate_make_lsq_public
