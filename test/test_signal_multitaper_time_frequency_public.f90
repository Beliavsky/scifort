! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_time_frequency_public
    use scifort_signal, only : dp, multitaper_spectrogram, multitaper_harmonic_ftest, &
        multitaper_spectrogram_result, multitaper_spectrogram_channels_result, &
        multitaper_harmonic_result, spectral_status_success
    implicit none
    type(multitaper_spectrogram_result) :: s
    type(multitaper_spectrogram_channels_result) :: sc
    type(multitaper_harmonic_result) :: f
    real(dp) :: x(48), y(48,2)
    integer :: i
    do i = 1, 48
        x(i) = cos(acos(-1.0_dp)*real(i,dp)/6.0_dp)
    end do
    y(:,1) = x
    y(:,2) = 2.0_dp*x
    s = multitaper_spectrogram(x,2.5_dp,24,kmax=4)
    sc = multitaper_spectrogram(y,2.5_dp,24,kmax=4)
    f = multitaper_harmonic_ftest(x,2.5_dp,kmax=4)
    if (s%status /= spectral_status_success .or. sc%status /= spectral_status_success .or. &
        f%status /= spectral_status_success) error stop 'multitaper public export'
    if (size(sc%power,3) /= 2) error stop 'multitaper public channel shape'
    print '(a)', 'Multitaper time-frequency public API passed'
end program test_signal_multitaper_time_frequency_public
