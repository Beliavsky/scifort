! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_interpolate_bspline
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use scifort_kinds, only: dp
    use scifort_interpolate, only: bspline_interpolator_type, bspline_array_interpolator_type, &
        bspline, make_interp_spline, bspline_design_matrix, bspline_status_success, &
        bspline_status_invalid_input
    implicit none
    type(bspline_interpolator_type) :: spl, deriv, anti, inserted, bad
    type(bspline_array_interpolator_type) :: channels, inserted_channels
    real(dp) :: x(7), y(7), yy(7,2), q(5), t(11), coef(7), values(5)
    real(dp), allocatable :: grid(:,:), cc(:,:), knots(:)
    integer :: i, code

    x = [0.0_dp,0.2_dp,0.7_dp,1.1_dp,1.8_dp,2.6_dp,3.0_dp]
    y = 1.0_dp-2.0_dp*x+0.5_dp*x**2-0.2_dp*x**3
    q = [0.0_dp,0.31_dp,1.1_dp,2.7_dp,3.0_dp]

    spl = make_interp_spline(x,y,k=3,status=code)
    if (code /= bspline_status_success .or. .not. spl%valid()) error stop 'bspline: create cubic'
    if (maxval(abs(spl%evaluate(x)-y)) > 2.0e-12_dp) error stop 'bspline: interpolate samples'
    if (maxval(abs(spl%evaluate(q)-(1.0_dp-2.0_dp*q+0.5_dp*q**2-0.2_dp*q**3))) > 2.0e-12_dp) &
        error stop 'bspline: exact polynomial'
    if (maxval(abs(spl%derivative(q)-( -2.0_dp+q-0.6_dp*q**2))) > 2.0e-12_dp) &
        error stop 'bspline: cubic first derivative'
    if (maxval(abs(spl%derivative(q,nu=2)-(1.0_dp-1.2_dp*q))) > 2.0e-12_dp) &
        error stop 'bspline: cubic second derivative'
    if (abs(spl%integrate(0.0_dp,3.0_dp)-(-5.55_dp)) > 3.0e-12_dp) &
        error stop 'bspline: exact integral'
    if (abs(spl%integrate(3.0_dp,0.0_dp)-5.55_dp) > 3.0e-12_dp) &
        error stop 'bspline: reversed exact integral'
    deriv = spl%derivative_spline(status=code)
    if (code /= 0 .or. deriv%degree() /= 2) error stop 'bspline: derivative object'
    if (maxval(abs(deriv%evaluate(q)-spl%derivative(q))) > 2.0e-12_dp) &
        error stop 'bspline: derivative object evaluation'
    anti = spl%antiderivative(status=code)
    if (code /= 0 .or. anti%degree() /= 4) error stop 'bspline: antiderivative object'
    if (maxval(abs(anti%derivative(q)-spl%evaluate(q))) > 2.0e-12_dp) &
        error stop 'bspline: antiderivative differential identity'
    inserted = spl%insert_knot(1.4_dp,m=2,status=code)
    if (code /= 0 .or. size(inserted%knots()) /= size(spl%knots())+2) &
        error stop 'bspline: insert knots'
    if (maxval(abs(inserted%evaluate(q)-spl%evaluate(q))) > 2.0e-12_dp) &
        error stop 'bspline: shape-preserving insertion'
    grid = bspline_design_matrix(q,spl%knots(),3,status=code)
    if (code /= 0 .or. any(shape(grid) /= [size(q),size(spl%coefficients())])) &
        error stop 'bspline: design matrix dimensions'
    values = matmul(grid,spl%coefficients())
    if (maxval(abs(values-spl%evaluate(q))) > 5.0e-13_dp) &
        error stop 'bspline: design matrix identity'

    spl = make_interp_spline(x,y,bc_type='natural',status=code)
    if (code /= 0) error stop 'bspline: natural construction'
    if (abs(spl%derivative(x(1),2)) > 5.0e-12_dp .or. &
        abs(spl%derivative(x(7),2)) > 5.0e-12_dp) error stop 'bspline: natural boundaries'
    spl = make_interp_spline(x,y,bc_type='clamped',status=code)
    if (code /= 0) error stop 'bspline: clamped construction'
    if (abs(spl%derivative(x(1))) > 5.0e-12_dp .or. &
        abs(spl%derivative(x(7))) > 5.0e-12_dp) error stop 'bspline: clamped boundaries'

    yy(:,1) = y
    yy(:,2) = 2.0_dp*y+1.0_dp
    channels = make_interp_spline(x,yy,k=5,status=code)
    if (code /= 0 .or. channels%channels() /= 2) error stop 'bspline: array constructor'
    cc = channels%evaluate(q)
    if (maxval(abs(cc(:,2)-2.0_dp*cc(:,1)-1.0_dp)) > 5.0e-12_dp) &
        error stop 'bspline: multichannel evaluation'
    inserted_channels = channels%insert_knot(1.21_dp,status=code)
    if (code /= 0 .or. maxval(abs(inserted_channels%evaluate(q)-cc)) > 5.0e-12_dp) &
        error stop 'bspline: multichannel knot insertion'

    ! Direct construction with repeated interior knots.
    t = [0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,2.0_dp,3.0_dp,3.0_dp,3.0_dp,3.0_dp]
    coef = [1.0_dp,1.0_dp,2.0_dp,0.0_dp,-1.0_dp,-2.0_dp,1.0_dp]
    spl = bspline(t,coef,3,status=code,extrapolate=.false.)
    if (code /= 0 .or. .not. ieee_is_nan(spl%evaluate(-1.0_dp))) &
        error stop 'bspline: out-of-bounds NaN'
    if (.not. ieee_is_nan(spl%evaluate(4.0_dp))) error stop 'bspline: outside right'
    if (abs(spl%integrate(-2.0_dp,-1.0_dp)) > 1.0e-12_dp) &
        error stop 'bspline: clipped integral'
    bad = bspline([0.0_dp,0.0_dp,2.0_dp,1.0_dp],[1.0_dp,2.0_dp],1,status=code)
    if (code /= bspline_status_invalid_input .or. bad%valid()) &
        error stop 'bspline: unsorted knots should be rejected'
    bad = make_interp_spline(x,y,k=5,bc_type='natural',status=code)
    if (code /= bspline_status_invalid_input .or. bad%valid()) &
        error stop 'bspline: unsupported boundary condition combination'

    ! k=0 is a step function; k=1 is piecewise linear interpolation.
    spl = make_interp_spline(x,y,k=0,status=code)
    if (code /= 0 .or. spl%evaluate(x(2)+0.01_dp) /= y(2)) &
        error stop 'bspline: degree zero step'
    spl = make_interp_spline(x,y,k=1,status=code)
    if (code /= 0 .or. abs(spl%evaluate(0.45_dp)- &
        (y(2)+(0.45_dp-x(2))*(y(3)-y(2))/(x(3)-x(2)))) > 2.0e-12_dp) &
        error stop 'bspline: linear interpolation'
    print '(a)', 'B-spline native regression passed'
end program test_interpolate_bspline
