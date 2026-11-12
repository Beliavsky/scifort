! SPDX-License-Identifier: MIT
program test_interpolate_bivariate_calculus
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : fitpack_bivariate_tck_type, rect_bivariate_spline, &
        bispline_partial_derivative, bispline_insert_knot, bisplev, &
        bispline_success, bispline_invalid_input
    implicit none
    type(fitpack_bivariate_tck_type) :: original, derivative, in_x, in_y, both, bad, copy, padded
    real(dp) :: x(9), y(11), z(9,11), qx, qy, exact, value, nanval
    integer :: i, j, status, ix, iy

    nanval = ieee_value(0.0_dp, ieee_quiet_nan)
    x = [-1.3_dp,-1.0_dp,-0.7_dp,-0.1_dp,0.25_dp,0.7_dp,1.0_dp,1.65_dp,2.0_dp]
    y = [-2.0_dp,-1.4_dp,-0.9_dp,-0.4_dp,0.1_dp,0.3_dp,0.7_dp,1.1_dp,1.8_dp,2.2_dp,2.9_dp]
    do j = 1, size(y)
        do i = 1, size(x)
            z(i,j) = 1.0_dp + 2.0_dp*x(i) - 3.0_dp*y(j) + &
                x(i)*y(j)**2 + 0.4_dp*x(i)**2*y(j)
        end do
    end do
    original = rect_bivariate_spline(x,y,z,status=status)
    if (status /= bispline_success) error stop 'original fit'

    derivative = original%partial_derivative(1,1,status)
    if (status /= bispline_success .or. .not. derivative%valid()) &
        error stop 'mixed partial derivative spline'
    if (derivative%kx /= 2 .or. derivative%ky /= 2) error stop 'derivative degree'
    if (size(derivative%tx) /= size(original%tx)-2) error stop 'derivative x knots'
    if (size(derivative%ty) /= size(original%ty)-2) error stop 'derivative y knots'
    if (any(derivative%tx /= original%tx(2:size(original%tx)-1))) error stop 'x trim'
    if (any(derivative%ty /= original%ty(2:size(original%ty)-1))) error stop 'y trim'

    do iy = 0, 7
        do ix = 0, 8
            qx = x(1) + (x(size(x))-x(1))*real(ix,dp)/8.0_dp
            qy = y(1) + (y(size(y))-y(1))*real(iy,dp)/7.0_dp
            exact = 2.0_dp*qy+0.8_dp*qx
            value = bisplev(qx,qy,derivative)
            if (abs(value-exact) > 5.0e-10_dp) &
                error stop 'derivative spline analytic values'
            if (abs(value-bisplev(qx,qy,original,dx=1,dy=1)) > 5.0e-10_dp) &
                error stop 'direct versus transformed mixed derivative'
        end do
    end do
    value = derivative%integral(-0.8_dp,1.3_dp,-0.5_dp,1.5_dp,status=status)
    exact = (1.3_dp+0.8_dp)*(1.5_dp**2-(-0.5_dp)**2) &
        + 0.4_dp*(1.3_dp**2-(-0.8_dp)**2)*(1.5_dp+0.5_dp)
    if (status /= bispline_success .or. abs(value-exact) > 8.0e-10_dp) &
        error stop 'integrated derivative analytic formula'
    derivative = bispline_partial_derivative(original,2,1,status)
    if (status /= bispline_success) error stop 'functional partial derivative'
    if (abs(bisplev(0.25_dp,1.7_dp,derivative)-0.8_dp) > 2.0e-9_dp) &
        error stop 'higher order derivative'

    in_x = original%insert_knot(1,0.48_dp,m=2,status=status)
    if (status /= bispline_success .or. .not. in_x%valid()) error stop 'insert x'
    if (size(in_x%tx) /= size(original%tx)+2) error stop 'x knot length'
    if (size(in_x%c) /= size(original%c)+2*(size(original%ty)-original%ky-1)) &
        error stop 'x coefficient array length'
    in_y = bispline_insert_knot(in_x,2,0.12_dp,m=2,status=status)
    if (status /= bispline_success .or. .not. in_y%valid()) error stop 'insert y'
    if (size(in_y%ty) /= size(original%ty)+2) error stop 'y knot length'
    both = in_y
    do iy = 0, 8
        do ix = 0, 8
            qx = x(1) + (x(size(x))-x(1))*real(ix,dp)/8.0_dp
            qy = y(1) + (y(size(y))-y(1))*real(iy,dp)/8.0_dp
            if (abs(bisplev(qx,qy,both)-bisplev(qx,qy,original)) > 2.0e-11_dp) &
                error stop 'knot insertion changed surface values'
            if (abs(bisplev(qx,qy,both,dx=1,dy=1)- &
                    bisplev(qx,qy,original,dx=1,dy=1)) > 2.0e-9_dp) &
                error stop 'knot insertion changed mixed derivatives'
        end do
    end do
    if (abs(both%integral(-1.0_dp,1.5_dp,-1.0_dp,2.0_dp) - &
            original%integral(-1.0_dp,1.5_dp,-1.0_dp,2.0_dp)) > 5.0e-10_dp) &
        error stop 'knot insertion changed integral'
    derivative = both%partial_derivative(1,1,status)
    if (status /= bispline_success) error stop 'derivative after refinement'
    if (abs(bisplev(0.7_dp,0.5_dp,derivative)- &
            bisplev(0.7_dp,0.5_dp,original,dx=1,dy=1)) > 2.0e-10_dp) &
        error stop 'post-insertion partial differentiation'
    copy = original%insert_knot(2,0.1_dp,m=0,status=status)
    if (status /= bispline_success) error stop 'zero knot insert'
    if (any(copy%c /= original%c) .or. any(copy%ty /= original%ty)) &
        error stop 'zero knot insert changed TCK'

    ! Upstream FITPACK may pad coefficient vectors past the active nx*ny
    ! tensor region. Both transformations must ignore trailing coefficients.
    padded = original
    padded%c = [original%c, 12345.0_dp, -999.0_dp]
    derivative = padded%partial_derivative(1,1,status)
    if (status /= bispline_success) error stop 'padded TCK derivative status'
    if (abs(bisplev(0.37_dp,0.81_dp,derivative)- &
            bisplev(0.37_dp,0.81_dp,original,dx=1,dy=1)) > 1.0e-10_dp) &
        error stop 'padded TCK derivative value'
    in_x = padded%insert_knot(1,padded%tx(5),m=2,status=status)
    if (status /= bispline_success .or. .not. in_x%valid()) &
        error stop 'insert at existing interior knot of padded TCK'
    if (abs(bisplev(0.37_dp,0.81_dp,in_x)-bisplev(0.37_dp,0.81_dp,original)) > 1.0e-11_dp) &
        error stop 'padded TCK insertion value'
    ! An internal knot of multiplicity k+1 has no valid differentiable
    ! TCK under the usual univariate/FITPACK derivative transformation.
    in_x = original%insert_knot(1,0.48_dp,m=4,status=status)
    if (status /= bispline_success) error stop 'full insertion multiplicity'
    bad = in_x%partial_derivative(1,0,status)
    if (status /= bispline_invalid_input .or. bad%valid()) &
        error stop 'singular derivative at fully repeated interior knot'

    bad = original%partial_derivative(3,0,status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'derivative order >= degree'
    bad = original%partial_derivative(-1,0,status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'negative derivative'
    bad = original%partial_derivative(0,4,status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'large y derivative'
    bad = original%insert_knot(3,0.4_dp,status=status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'bad refinement axis'
    bad = original%insert_knot(1,0.4_dp,m=-1,status=status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'negative repetitions'
    bad = original%insert_knot(1,nanval,status=status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'nonfinite insertion site'
    bad = original%insert_knot(1,x(1),status=status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'endpoint already full'
    bad = original%insert_knot(2,y(size(y))+1.0_dp,status=status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'insertion out of bounds'
    bad = original%insert_knot(1,0.48_dp,m=5,status=status)
    if (status /= bispline_invalid_input .or. bad%valid()) error stop 'insertion exceeds multiplicity'
    bad = original%insert_knot(1,x(1),m=0,status=status)
    if (status /= bispline_success .or. .not. bad%valid()) error stop 'zero insert endpoint valid'

    print '(a)', 'PASS: bivariate spline TCK differentiation, knot insertion, and integration'
end program test_interpolate_bivariate_calculus
