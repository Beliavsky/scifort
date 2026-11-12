! SPDX-License-Identifier: MIT
program test_interpolate_cubic_spline
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : piecewise_cubic_interpolator, &
        piecewise_cubic_array_interpolator, cubic_spline, &
        interpolate_status_success, interpolate_status_invalid_input
    implicit none
    type(piecewise_cubic_interpolator) :: poly, periodic, invalid
    type(piecewise_cubic_array_interpolator) :: matrix_poly
    real(dp) :: x(5), y(5), vals(5,3), q(7), expected(7), pi, h, a, b
    real(dp), allocatable :: c(:,:,:), d(:,:), z(:)
    integer :: i, status, mode

    pi = acos(-1.0_dp)
    x = [0.0_dp, 0.75_dp, 1.5_dp, 3.0_dp, 5.0_dp]
    y = 1.2_dp + 0.3_dp*x -0.4_dp*x*x + 0.2_dp*x**3
    q = [-0.25_dp,0.0_dp,0.5_dp,1.5_dp,2.7_dp,5.0_dp,6.0_dp]
    expected = 1.2_dp+0.3_dp*q-0.4_dp*q*q+0.2_dp*q**3
    poly = cubic_spline(x,y,status=status)
    if (status /= interpolate_status_success .or. .not. poly%valid()) &
        error stop 'not-a-knot spline fit failed'
    if (maxval(abs(poly%evaluate(q)-expected)) > 3.0e-12_dp) &
        error stop 'not-a-knot cubic exactness'
    if (maxval(abs(poly%derivative(q)- (0.3_dp-0.8_dp*q+0.6_dp*q*q))) &
        > 5.0e-12_dp) error stop 'not-a-knot first derivative'
    if (maxval(abs(poly%derivative(q,2)-(-0.8_dp+1.2_dp*q))) &
        > 5.0e-12_dp) error stop 'not-a-knot second derivative'
    if (maxval(abs(poly%derivative(q,3)-1.2_dp)) > 2.0e-12_dp) &
        error stop 'not-a-knot third derivative'
    if (abs(poly%integrate(0.0_dp,5.0_dp)- &
        (1.2_dp*5.0_dp+0.15_dp*25.0_dp-0.4_dp/3.0_dp*125.0_dp+0.05_dp*625.0_dp)) &
        > 5.0e-12_dp) error stop 'not-a-knot exact integral'

    poly = cubic_spline([0.0_dp,1.0_dp], [2.0_dp,3.0_dp],status=status)
    if (status /= 0 .or. abs(poly%evaluate(2.0_dp)-4.0_dp) > 1.0e-14_dp) &
        error stop 'two-knot not-a-knot reduces to line'
    poly = cubic_spline([0.0_dp,0.75_dp,2.0_dp], &
        [1.0_dp,1.0_dp+0.75_dp**2,5.0_dp],status=status)
    if (status /= 0 .or. abs(poly%evaluate(1.5_dp)-3.25_dp) > 1.0e-13_dp) &
        error stop 'three-knot not-a-knot quadratic'

    ! Natural endpoints have vanishing second derivatives.
    poly = cubic_spline(x,sin(x),'natural',status=status)
    if (status /= 0 .or. abs(poly%derivative(x(1),2)) > 1.0e-12_dp &
        .or. abs(poly%derivative(x(5),2)) > 1.0e-12_dp) &
        error stop 'natural boundary conditions'

    poly = cubic_spline(x,sin(x),'clamped',status=status)
    if (status /= 0 .or. abs(poly%derivative(x(1),1)) > 1.0e-12_dp &
        .or. abs(poly%derivative(x(5),1)) > 1.0e-12_dp) &
        error stop 'clamped zero slopes'
    poly = cubic_spline(x,sin(x),bc_left='first',bc_right='second', &
        left_value=1.5_dp,right_value=-0.7_dp,status=status)
    if (status /= 0 .or. abs(poly%derivative(x(1),1)-1.5_dp) > 1.0e-12_dp &
        .or. abs(poly%derivative(x(5),2)+0.7_dp) > 1.0e-12_dp) &
        error stop 'asymmetric derivative boundary conditions'

    ! Periodic extrapolation and exact integration over full periods.
    x = [0.0_dp,0.5_dp,1.75_dp,4.0_dp,2.0_dp*pi]
    y = sin(x)
    y(5) = y(1)
    periodic = cubic_spline(x,y,'periodic',status=status)
    if (status /= 0 .or. .not. periodic%valid()) &
        error stop 'periodic fit failed'
    h = x(5)-x(1)
    do i = -5,5
        a = 1.2_dp+real(i,dp)*h
        if (abs(periodic%evaluate(a)-periodic%evaluate(1.2_dp)) > 5.0e-13_dp) &
            error stop 'periodic extension'
        if (abs(periodic%derivative(a,2)-periodic%derivative(1.2_dp,2)) > 5.0e-13_dp) &
            error stop 'periodic second derivative extension'
    end do
    do mode = -3,3
        a = -1.7_dp+real(mode,dp)*h
        b = a+h
        if (abs(periodic%integrate(a,b)-periodic%integrate(0.0_dp,h)) &
            > 1.0e-12_dp) error stop 'periodic full-period integral'
    end do
    if (abs(periodic%derivative(x(1))-periodic%derivative(x(5))) > 2.0e-12_dp) &
        error stop 'periodic slope matching'
    if (abs(periodic%derivative(x(1),2)-periodic%derivative(x(5),2)) > 2.0e-12_dp) &
        error stop 'periodic curvature matching'
    poly = cubic_spline(x,y,'periodic',status=status,extrapolate=.false.)
    if (status /= 0 .or. .not. ieee_is_nan(poly%evaluate(-1.0_dp))) &
        error stop 'periodic strict extrapolation override'
    if (ieee_is_nan(periodic%evaluate(-1.0_dp)) .or. &
        ieee_is_nan(periodic%integrate(-10.0_dp,10.0_dp))) &
        error stop 'periodic extension integral'
    if (abs(periodic%integrate(-1.0_dp,4.0_dp)+periodic%integrate(4.0_dp,-1.0_dp)) &
        > 2.0e-12_dp) error stop 'periodic reversed integral'

    ! Matrix ordinates: y(:,k) is the kth independent output channel.
    vals(:,1) = 2.0_dp + 0.5_dp*x
    vals(:,2) = x*x
    vals(:,3) = 5.0_dp-3.0_dp*x+0.5_dp*x**3
    matrix_poly = cubic_spline(x,vals,status=status)
    if (status /= 0 .or. .not. matrix_poly%valid()) &
        error stop 'matrix spline fit'
    c = matrix_poly%coefficients()
    if (any(shape(c) /= [4,4,3])) error stop 'matrix coefficients shape'
    d = matrix_poly%evaluate([0.25_dp,2.5_dp])
    if (any(shape(d) /= [2,3])) error stop 'matrix evaluate array shape'
    if (abs(d(1,1)-2.125_dp) > 5.0e-12_dp .or. &
        abs(d(2,2)-6.25_dp) > 5.0e-12_dp .or. &
        abs(d(1,3)-(5.0_dp-0.75_dp+0.5_dp*(0.25_dp**3))) > 5.0e-12_dp) &
        error stop 'matrix analytic values'
    z = matrix_poly%evaluate(1.5_dp)
    if (size(z) /= 3 .or. abs(z(2)-2.25_dp) > 1.0e-12_dp) &
        error stop 'matrix scalar evaluate'
    z = matrix_poly%derivative(1.5_dp)
    if (abs(z(3)-(-3.0_dp+1.5_dp*1.5_dp*1.5_dp)) > 5.0e-12_dp) &
        error stop 'matrix first derivative'
    d = matrix_poly%derivative([1.5_dp,2.5_dp],order=2)
    if (abs(d(2,2)-2.0_dp) > 5.0e-12_dp) &
        error stop 'matrix second derivative'
    z = matrix_poly%integrate(0.0_dp,h)
    if (abs(z(1)-(2.0_dp*h+0.25_dp*h*h)) > 1.0e-11_dp) &
        error stop 'matrix integral'

    invalid = cubic_spline([0.0_dp,1.0_dp],[0.0_dp,3.0_dp], &
        bc_type='unknown',status=status)
    if (status /= interpolate_status_invalid_input .or. invalid%valid()) &
        error stop 'reject unknown boundary conditions'
    invalid = cubic_spline([0.0_dp,1.0_dp],[0.0_dp,3.0_dp], &
        bc_type='periodic',status=status)
    if (status /= interpolate_status_invalid_input .or. invalid%valid()) &
        error stop 'reject mismatched periodic values'
    invalid = cubic_spline([0.0_dp,1.0_dp],[0.0_dp,3.0_dp], &
        bc_left='periodic',bc_right='natural',status=status)
    if (status /= interpolate_status_invalid_input) &
        error stop 'reject mixed periodic endpoints'
    invalid = cubic_spline([0.0_dp,1.0_dp],[0.0_dp,3.0_dp], &
        bc_type='natural',left_value=ieee_value(0.0_dp,ieee_quiet_nan),status=status)
    if (status /= interpolate_status_invalid_input .or. invalid%valid()) &
        error stop 'reject nonfinite boundary parameter'
    matrix_poly = cubic_spline(x, vals(:,1:2),left_values=[1.0_dp],status=status)
    if (status /= interpolate_status_invalid_input .or. matrix_poly%valid()) &
        error stop 'reject wrong number of boundary values'
    matrix_poly = cubic_spline(x, vals(:,1:0),status=status)
    if (status /= interpolate_status_invalid_input .or. matrix_poly%valid()) &
        error stop 'reject empty output channel dimension'

    print '(a)', 'PASS CubicSpline native boundary, periodic, matrix and invalid-input tests'
end program test_interpolate_cubic_spline
