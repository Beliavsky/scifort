! SPDX-License-Identifier: MIT
! Native regression for automatic Whittaker-Henderson REML selection.
program test_signal_whittaker_reml
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_whittaker, only : whittaker_henderson, whittaker_result, &
        whittaker_status_success, whittaker_status_invalid_input, whittaker_status_singular
    implicit none
    real(dp) :: y(35), w(35), nan
    type(whittaker_result) :: a, b, c
    integer :: j

    do j=1,size(y)
        y(j) = 0.03_dp*(j-1)+sin(0.12_dp*(j-1))+0.23_dp*cos(0.51_dp*(j-1)) + &
            0.2_dp*sin(2.1_dp*(j-1))
    end do
    a = whittaker_henderson(y)
    if (a%status /= whittaker_status_success) error stop 'REML fit'
    if (a%lamb <= 0.0_dp .or. a%lamb > 1.0e10_dp) error stop 'REML lambda range'
    b = whittaker_henderson(y,a%lamb)
    if (b%status /= whittaker_status_success) error stop 'fixed lambda fit'
    if (maxval(abs(a%x-b%x)) > 1.0e-11_dp) error stop 'REML fixed penalty identity'
    w = 1.0_dp
    w(8:12) = 0.0_dp
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    y(8:12) = nan
    c = whittaker_henderson(y,order=2,weights=w)
    if (c%status /= whittaker_status_success) error stop 'REML masked NaN'
    if (c%lamb <= 0.0_dp) error stop 'REML masked penalty'
    if (any(.not. ieee_is_finite(c%x))) error stop 'REML invalid smooth'
    c = whittaker_henderson(y,order=2)
    if (c%status /= whittaker_status_invalid_input) error stop 'REML unmasked NaN'
    c = whittaker_henderson(y,order=2,weights=0.0_dp*w)
    if (c%status /= whittaker_status_singular) error stop 'REML singular weights'
    print *, 'test_signal_whittaker_reml passed'
end program test_signal_whittaker_reml
