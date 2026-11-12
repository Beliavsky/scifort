! SPDX-License-Identifier: MIT
program test_signal_band_stop_obj
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_band_stop_obj, only : band_stop_obj, &
        band_stop_success, band_stop_invalid_input
    implicit none
    integer :: status
    real(dp) :: f, expected, passb(2), stopb(2), nat, gpass, gstop

    passb = [0.5_dp,4.0_dp]
    stopb = [1.5_dp,2.5_dp]
    gpass = 1.0_dp
    gstop = 40.0_dp
    nat = minval(abs(stopb*(passb(1)-passb(2))/(stopb**2-passb(1)*passb(2))))
    expected = log((10.0_dp**(gstop/10.0_dp)-1.0_dp)/ &
                   (10.0_dp**(gpass/10.0_dp)-1.0_dp))/(2.0_dp*log(nat))
    f = band_stop_obj(passb(1),0,passb,stopb,gpass,gstop,'butter',status)
    if (status /= band_stop_success) error stop 'butter status'
    if (abs(f-expected) > 1.0e-13_dp) error stop 'butter value'
    expected = acosh(sqrt((10.0_dp**(gstop/10.0_dp)-1.0_dp)/ &
                          (10.0_dp**(gpass/10.0_dp)-1.0_dp)))/acosh(nat)
    f = band_stop_obj(passb(2),1,passb,stopb,gpass,gstop,'cheby',status)
    if (status /= band_stop_success) error stop 'cheby status'
    if (abs(f-expected) > 1.0e-13_dp) error stop 'cheby value'
    f = band_stop_obj(passb(1),0,passb,stopb,3.0_dp,3.0_dp,'ellip',status)
    if (status /= band_stop_success .or. abs(f) > 1.0e-14_dp) &
        error stop 'equal ripple/attenuation should yield zero order'
    f = band_stop_obj(passb(1),0,passb,stopb,gpass,gstop,'ellip',status)
    if (status /= band_stop_success .or. .not. (f>0.0_dp .and. f<6.0_dp)) &
        error stop 'ellip sanity'
    f = band_stop_obj(passb(1),3,passb,stopb,gpass,gstop,'butter',status)
    if (status /= band_stop_invalid_input .or. .not. ieee_is_nan(f)) &
        error stop 'invalid index'
    f = band_stop_obj(passb(1),0,passb,stopb,0.0_dp,gstop,'butter',status)
    if (status /= band_stop_invalid_input .or. .not. ieee_is_nan(f)) &
        error stop 'invalid gpass'
    f = band_stop_obj(passb(1),0,passb,stopb,gstop,gpass,'butter',status)
    if (status /= band_stop_invalid_input .or. .not. ieee_is_nan(f)) &
        error stop 'invalid gstop'
    f = band_stop_obj(passb(1),0,passb,stopb,gpass,gstop,'other',status)
    if (status /= band_stop_invalid_input .or. .not. ieee_is_nan(f)) &
        error stop 'invalid filter'
    print *, 'PASS: native band_stop_obj'
end program test_signal_band_stop_obj
