! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_noise_floor
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input, &
        spectral_result
    use scifort_signal_multitaper, only : multitaper_noise_floor, multitaper_noise_floor_result, &
        multitaper_psd
    implicit none
    real(dp) :: x(96), nan_input(96), t, pi
    complex(dp) :: z(96)
    type(multitaper_noise_floor_result) :: floor
    type(spectral_result) :: raw
    integer :: i

    pi = acos(-1.0_dp)
    do i = 1, size(x)
        t = real(i-1,dp)
        x(i) = 1.2_dp*sin(2.0_dp*pi*17.0_dp*t/96.0_dp) &
            +0.18_dp*cos(2.0_dp*pi*4.0_dp*t/96.0_dp) + 0.001_dp*t
        z(i) = cmplx(x(i),0.3_dp*sin(2.0_dp*pi*11.0_dp*t/96.0_dp),dp)
    end do
    raw = multitaper_psd(x,3.5_dp,kmax=5,fs=32.0_dp,nfft=112, &
        detrend_mode='linear',weighting='eigen')
    floor = multitaper_noise_floor(x,3.5_dp,5,kmax=5,fs=32.0_dp,nfft=112, &
        detrend_mode='linear',weighting='eigen')
    if (floor%status /= spectral_status_success) error stop 'real floor status'
    if (size(floor%power) /= 57) error stop 'real floor shape'
    if (maxval(abs(floor%frequencies-raw%frequencies)) > 1.0e-13_dp) error stop 'real frequency'
    if (maxval(abs(floor%power-raw%power)) > 1.0e-12_dp) error stop 'real power'
    if (any(floor%background < 0.0_dp)) error stop 'negative background'
    if (any(floor%whitened < 0.0_dp)) error stop 'negative whitened'
    do i = 1, size(floor%power)
        if (floor%background(i) > 0.0_dp) then
            if (abs(floor%whitened(i)*floor%background(i)-floor%power(i)) &
                > 1.0e-10_dp) error stop 'real ratio'
        end if
    end do

    floor = multitaper_noise_floor(z,2.5_dp,7,kmax=4,weighting='adaptive', &
        detrend_mode='constant',nfft=127)
    raw = multitaper_psd(z,2.5_dp,kmax=4,weighting='adaptive', &
        detrend_mode='constant',nfft=127)
    if (floor%status /= spectral_status_success) error stop 'complex floor status'
    if (size(floor%power) /= 127) error stop 'complex floor shape'
    if (maxval(abs(floor%power-raw%power)) > 1.0e-12_dp) error stop 'complex adaptive power'
    if (any(floor%background < 0.0_dp)) error stop 'complex background'

    ! Check a strided input and the zero-signal convention.
    floor = multitaper_noise_floor(x(1::2),2.0_dp,3,nfft=64,detrend_mode='none')
    if (floor%status /= spectral_status_success) error stop 'strided input'
    floor = multitaper_noise_floor(0.0_dp*x,2.5_dp,5)
    if (floor%status /= spectral_status_success) error stop 'zero spectrum status'
    if (any(floor%background /= 0.0_dp)) error stop 'zero background'
    if (any(floor%whitened /= 0.0_dp)) error stop 'zero normalized power'

    floor = multitaper_noise_floor(x,2.5_dp,4)
    if (floor%status /= spectral_status_invalid_input) error stop 'even median accepted'
    if (size(floor%power) /= 0) error stop 'even median output'
    floor = multitaper_noise_floor(x,2.5_dp,1)
    if (floor%status /= spectral_status_invalid_input) error stop 'short median accepted'
    floor = multitaper_noise_floor(x,2.5_dp,99)
    if (floor%status /= spectral_status_invalid_input) error stop 'wide median accepted'
    nan_input = x
    nan_input(3) = ieee_value(0.0_dp,ieee_quiet_nan)
    floor = multitaper_noise_floor(nan_input,2.5_dp,5)
    if (floor%status /= spectral_status_invalid_input) error stop 'nan signal accepted'
    floor = multitaper_noise_floor(x,2.5_dp,5,weighting='invalid')
    if (floor%status /= spectral_status_invalid_input) error stop 'invalid weighting'

    print '(a)', 'test_signal_multitaper_noise_floor passed'
end program test_signal_multitaper_noise_floor
