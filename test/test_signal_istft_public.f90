! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_istft_public
    use scifort_signal, only : dp, stft, stft_result, istft, check_cola, check_nola, &
        istft_status_success
    implicit none
    real(dp) :: samples(32)
    real(dp), allocatable :: reconstructed(:)
    type(stft_result) :: analysis
    integer :: j, status

    samples = [(cos(real(j,dp)/6.0_dp), j=1,32)]
    analysis = stft(samples, nperseg=8, noverlap=4)
    if (analysis%status /= 0) error stop 'public STFT interface'
    call istft(analysis%values, reconstructed, status, nperseg=8, noverlap=4)
    if (status /= istft_status_success) error stop 'public ISTFT interface'
    if (size(reconstructed) /= size(samples)) error stop 'public reconstruction shape'
    if (maxval(abs(reconstructed - samples)) > 2.0e-12_dp) error stop 'public reconstruction'
    if (.not. check_cola('hann', 8, 4)) error stop 'public COLA'
    if (.not. check_nola('hann', 8, 4)) error stop 'public NOLA'
    print '(a)', 'test_signal_istft_public passed'
end program test_signal_istft_public
