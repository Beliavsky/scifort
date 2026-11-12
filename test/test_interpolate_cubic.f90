! SPDX-License-Identifier: MIT
program test_interpolate_cubic
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : piecewise_cubic_interpolator, &
        pchip_interpolator, akima1d_interpolator, cubic_hermite_spline, pchip_interpolate, &
        interpolate_status_success, interpolate_status_invalid_input
    implicit none
    type(piecewise_cubic_interpolator) :: poly, akima, makima
    real(dp) :: x(4), y(4), q(6), expect(6), c(4,3)
    integer :: status, i

    ! Both schemes reproduce straight lines exactly, including derivatives.
    x = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
    y = 1.5_dp + 2.25_dp*x
    q = [-0.5_dp, 0.0_dp, 0.5_dp, 1.5_dp, 3.5_dp, 5.0_dp]
    expect = 1.5_dp + 2.25_dp*q
    poly = pchip_interpolator(x, y, status)
    if (status /= interpolate_status_success .or. .not. poly%valid()) &
        error stop 'pchip fit failed'
    if (maxval(abs(poly%evaluate(q)-expect)) > 1.0e-12_dp) &
        error stop 'pchip line interpolation'
    if (maxval(abs(poly%derivative(q)-2.25_dp)) > 1.0e-12_dp) &
        error stop 'pchip line first derivative'
    if (maxval(abs(poly%derivative(q, 2))) > 1.0e-12_dp) &
        error stop 'pchip line second derivative'
    if (maxval(abs(poly%derivative(q, 3))) > 1.0e-12_dp) &
        error stop 'pchip line third derivative'
    if (maxval(abs(poly%derivative(q, 4))) > 1.0e-12_dp) &
        error stop 'pchip line fourth derivative'
    if (abs(poly%integrate(-0.5_dp, 5.0_dp) - &
        (1.5_dp*5.5_dp + 1.125_dp*(25.0_dp-0.25_dp))) > 1.0e-12_dp) &
        error stop 'pchip exact line integral'
    if (abs(poly%integrate(5.0_dp,-0.5_dp)+poly%integrate(-0.5_dp,5.0_dp)) &
        > 1.0e-12_dp) error stop 'reverse integration sign'
    if (maxval(abs(poly%knots()-x)) > 0.0_dp) error stop 'knots accessor'
    c = poly%coefficients()
    if (maxval(abs(c(1,:)-y(:3))) > 0.0_dp) error stop 'coefficient accessor'
    if (maxval(abs(c(3:4,:))) > 1.0e-12_dp) error stop 'linear coefficients'
    if (abs(pchip_interpolate(x,y,1.5_dp)-expect(4)) > 1.0e-12_dp) &
        error stop 'scalar convenience'
    if (maxval(abs(pchip_interpolate(x,y,q)-expect)) > 1.0e-12_dp) &
        error stop 'array convenience'
    if (maxval(abs(pchip_interpolate(x,y,q,order=1)-2.25_dp)) > 1.0e-12_dp) &
        error stop 'array derivative convenience'

    akima = akima1d_interpolator(x,y,status=status)
    makima = akima1d_interpolator(x,y,'makima',status)
    if (.not. akima%valid() .or. .not. makima%valid() .or. status /= 0) &
        error stop 'akima fit failed'
    if (maxval(abs(akima%evaluate(q(2:5))-expect(2:5))) > 1.0e-12_dp) &
        error stop 'akima line interpolation'
    if (maxval(abs(makima%evaluate(q(2:5))-expect(2:5))) > 1.0e-12_dp) &
        error stop 'makima line interpolation'
    if (.not. ieee_is_nan(akima%evaluate(q(1)))) &
        error stop 'akima default left boundary policy'
    if (.not. ieee_is_nan(makima%evaluate(q(6)))) &
        error stop 'makima default right boundary policy'
    if (abs(akima%evaluate(q(1),extrapolate=.true.) - expect(1)) > 1.0e-12_dp) &
        error stop 'per-query extrapolation override'
    if (.not. ieee_is_nan(akima%integrate(-0.5_dp,1.0_dp))) &
        error stop 'akima integration outside support'
    if (abs(akima%integrate(-0.5_dp,1.0_dp,extrapolate=.true.) - &
        poly%integrate(-0.5_dp,1.0_dp)) > 1.0e-12_dp) &
        error stop 'akima integral override'

    ! Exact cubic Hermite spline from known values and first derivatives.
    y = x**3 - 2.0_dp*x + 5.0_dp
    poly = cubic_hermite_spline(x, y, 3.0_dp*x*x - 2.0_dp, status)
    if (status /= 0 .or. .not. poly%valid()) error stop 'Hermite fit'
    if (maxval(abs(poly%evaluate(q)-(q**3-2.0_dp*q+5.0_dp))) &
        > 1.0e-11_dp) error stop 'Hermite exact cubic values'
    if (maxval(abs(poly%derivative(q)-(3.0_dp*q*q-2.0_dp))) &
        > 1.0e-11_dp) error stop 'Hermite cubic derivative'
    if (abs(poly%integrate(-0.5_dp,5.0_dp) - &
        ((5.0_dp**4-(-0.5_dp)**4)/4.0_dp - &
        (5.0_dp**2-(-0.5_dp)**2) + 5.0_dp*5.5_dp)) &
        > 1.0e-10_dp) error stop 'Hermite polynomial integral'
    poly = cubic_hermite_spline(x, y, [1.0_dp], status)
    if (status /= interpolate_status_invalid_input .or. poly%valid()) &
        error stop 'reject invalid Hermite slopes'
    poly = cubic_hermite_spline([0.0_dp,1.0_dp], [0.0_dp,1.0_dp], &
        [ieee_value(0.0_dp,ieee_quiet_nan),1.0_dp], status)
    if (status /= interpolate_status_invalid_input .or. poly%valid()) &
        error stop 'reject NaN Hermite slopes'

    ! PCHIP must interpolate monotonic data without overshooting knot values.
    y = [0.0_dp, 0.1_dp, 3.0_dp, 3.0_dp]
    poly = pchip_interpolator(x,y,status)
    if (status /= 0) error stop 'monotonic fit failed'
    do i = 0, 400
        if (poly%evaluate(real(i,dp)*0.01_dp) < -1.0e-12_dp) &
            error stop 'monotone pchip undershoot'
        if (poly%evaluate(real(i,dp)*0.01_dp) > 3.0_dp+1.0e-12_dp) &
            error stop 'monotone pchip overshoot'
    end do
    if (abs(poly%derivative(4.0_dp)) > 1.0e-12_dp) &
        error stop 'flat last pchip interval'

    poly = pchip_interpolator([1.0_dp], [3.0_dp], status)
    if (status /= interpolate_status_invalid_input .or. poly%valid()) &
        error stop 'reject one knot'
    if (.not. ieee_is_nan(poly%evaluate(1.0_dp))) &
        error stop 'invalid object must return NaN'
    poly = pchip_interpolator([0.0_dp, 1.0_dp, 1.0_dp], [0.0_dp,2.0_dp,3.0_dp],status)
    if (status /= interpolate_status_invalid_input .or. poly%valid()) &
        error stop 'reject repeated knot'
    poly = pchip_interpolator([0.0_dp, -1.0_dp], [1.0_dp,3.0_dp],status)
    if (status /= interpolate_status_invalid_input) error stop 'reject descending knots'
    poly = pchip_interpolator([0.0_dp,1.0_dp], [3.0_dp],status)
    if (status /= interpolate_status_invalid_input) error stop 'reject sample size mismatch'
    poly = pchip_interpolator([0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)], &
        [3.0_dp,4.0_dp],status)
    if (status /= interpolate_status_invalid_input) error stop 'reject NaN knots'
    akima = akima1d_interpolator([0.0_dp,1.0_dp], [1.0_dp,2.0_dp], &
        method='bad',status=status)
    if (status /= interpolate_status_invalid_input .or. akima%valid()) &
        error stop 'reject invalid Akima method'
    poly = pchip_interpolator(x,y,status,extrapolate=.false.)
    if (.not. ieee_is_nan(poly%evaluate(-1.0_dp))) &
        error stop 'pchip extrapolation policy'
    if (.not. ieee_is_nan(poly%integrate(0.0_dp, 5.0_dp))) &
        error stop 'pchip integration policy'
    if (.not. ieee_is_nan(poly%derivative(1.0_dp,-1))) &
        error stop 'negative derivative order'

    print '(a)', 'PASS PCHIP/Akima/Makima/CubicHermite native tests'
end program test_interpolate_cubic
