! SPDX-License-Identifier: MIT
! SciPy-compatible Gaussian B-spline approximation and Kaiser FIR helpers.
program test_signal_kaiser_spline
    use scifort_signal, only : dp, gauss_spline, kaiser_beta, kaiser_atten, &
        kaiserord, kaiser_order_result, get_window, firwin, &
        design_status_success, design_status_invalid_input
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: win(:), taps(:)
    real(dp) :: x(3), y(2,3)
    type(kaiser_order_result) :: result
    integer :: status

    x = [-1.0_dp, 0.0_dp, 1.0_dp]
    if (maxval(abs(gauss_spline(x,3)- &
        [0.15418032980376928_dp, 0.6909882989426710_dp, 0.15418032980376928_dp])) > 1e-13_dp) &
        error stop 'gauss_spline example'
    y = reshape([0.0_dp, 0.5_dp, -0.5_dp, 2.0_dp, 1.0_dp, -1.0_dp], shape(y))
    if (maxval(abs(gauss_spline(y,7) - &
        exp(-y*y/(2.0_dp*(8.0_dp/12.0_dp)))/sqrt(2.0_dp*pi*(8.0_dp/12.0_dp)))) > 1e-14_dp) &
        error stop 'rank-2 elemental spline'
    if (.not.ieee_is_nan(gauss_spline(1.0_dp,-1))) error stop 'invalid spline order'
    if (gauss_spline(huge(0.0_dp),2) /= 0.0_dp) error stop 'distant spline knot'

    if (abs(kaiser_beta(65.0_dp)-6.20426_dp) > 1e-13_dp) error stop 'kaiser beta example'
    if (kaiser_beta(21.0_dp) /= 0.0_dp) error stop 'kaiser beta boundary'
    if (abs(kaiser_beta(50.0_dp)-4.533514120981248_dp) > 1e-13_dp) error stop 'kaiser beta 50'
    if (abs(kaiser_atten(211,0.0375_dp)-64.48099630593983_dp) > 1e-12_dp) &
        error stop 'kaiser attenuation example'

    result = kaiserord(65.0_dp, 24.0_dp/500.0_dp)
    if (result%status /= design_status_success) error stop 'kaiserord status'
    if (result%numtaps /= 167) error stop 'kaiserord tap count'
    if (abs(result%beta-6.20426_dp) > 1e-13_dp) error stop 'kaiserord beta'
    result = kaiserord(7.0_dp,0.05_dp)
    if (result%status /= design_status_invalid_input) error stop 'bad ripple'
    result = kaiserord(65.0_dp,0.0_dp)
    if (result%status /= design_status_invalid_input) error stop 'zero width'
    result = kaiserord(65.0_dp,-0.01_dp)
    if (result%status /= design_status_invalid_input) error stop 'negative width'
    result = kaiserord(100.0_dp,1.0e-30_dp)
    if (result%status /= design_status_invalid_input) error stop 'tap overflow'
    result = kaiserord(-65.0_dp,0.048_dp)
    if (result%numtaps /= 167) error stop 'negative ripple magnitude'

    win = get_window('kaiser', 31, fftbins=.false., parameter=6.20426_dp, status=status)
    if (status /= 0) error stop 'Kaiser window with parameter'
    taps = firwin(31, [0.25_dp], window='kaiser', window_parameter=6.20426_dp, status=status)
    if (status /= design_status_success .or. size(taps) /= 31) error stop 'Kaiser FIR design'
    if (abs(sum(taps)-1.0_dp) > 2e-14_dp) error stop 'Kaiser FIR normalization'
    if (maxval(abs(taps-taps(31:1:-1))) > 1e-14_dp) error stop 'Kaiser FIR symmetry'
    taps = firwin(31, [0.25_dp], window='hann', window_parameter=0.9_dp, status=status)
    if (status /= design_status_invalid_input) error stop 'unexpected Hann shape parameter'

    print '(a)', 'test_signal_kaiser_spline passed'
end program test_signal_kaiser_spline
