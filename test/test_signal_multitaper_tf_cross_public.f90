! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_tf_cross_public
    use scifort_signal, only : dp, spectral_status_success, &
        multitaper_cross_spectrogram, multitaper_tf_coherence, multitaper_spectral_matrix_spectrogram, &
        multitaper_cross_spectrogram_result, multitaper_tf_coherence_result, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    real(dp) :: x(48), y(48,2)
    type(multitaper_cross_spectrogram_result) :: s
    type(multitaper_tf_coherence_result) :: c
    type(multitaper_spectral_matrix_spectrogram_result) :: m
    integer :: i
    do i = 1, 48
        x(i) = sin(2.0_dp*acos(-1.0_dp)*real(i,dp)/12.0_dp)
    end do
    y(:,1) = x
    y(:,2) = 3.0_dp*x
    s = multitaper_cross_spectrogram(y(:,1),y(:,2),2.5_dp,24,kmax=4)
    c = multitaper_tf_coherence(y(:,1),y(:,2),2.5_dp,24,kmax=4)
    m = multitaper_spectral_matrix_spectrogram(y,2.5_dp,24,kmax=4)
    if (s%status /= spectral_status_success .or. c%status /= spectral_status_success .or. &
        m%status /= spectral_status_success) error stop 'public multitaper time-frequency exports'
    if (size(m%power,4) /= 2 .or. maxval(abs(s%power-m%power(:,:,1,2))) > 1.0e-12_dp) &
        error stop 'public matrix/cross result'
    print '(a)', 'Multitaper cross time-frequency public API passed'
end program test_signal_multitaper_tf_cross_public
