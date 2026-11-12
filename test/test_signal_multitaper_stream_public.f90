! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_stream_public
    use scifort_signal, only : dp, multitaper_stream, multitaper_spectral_matrix_spectrogram_result, &
        multitaper_matrix_result, spectral_status_success
    implicit none
    type(multitaper_stream) :: stream
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(multitaper_matrix_result) :: running
    real(dp) :: x(50,2)
    integer :: status, i
    do i = 1, size(x,1)
        x(i,1) = sin(real(i,dp)*0.1_dp)
        x(i,2) = cos(real(i,dp)*0.2_dp)
    end do
    call stream%initialize(2.5_dp,16,2,status,noverlap=8,kmax=4)
    if (status /= spectral_status_success) error stop 'stream public initialization'
    spectra = stream%push(x)
    if (spectra%status /= spectral_status_success) error stop 'stream public push'
    if (size(spectra%times) /= 5) error stop 'stream public number of windows'
    running = stream%mean_matrix()
    if (running%status /= spectral_status_success) error stop 'stream public running mean'
    if (maxval(abs(running%power-sum(spectra%power,dim=2)/5.0_dp)) > 1.0e-12_dp) &
        error stop 'stream public mean is inconsistent'
    print '(a)', 'Streaming multitaper public interface passed'
end program test_signal_multitaper_stream_public
