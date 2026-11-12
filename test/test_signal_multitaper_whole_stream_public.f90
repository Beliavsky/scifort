! SPDX-License-Identifier: MIT
program test_signal_multitaper_whole_stream_public
    use scifort_signal, only : dp, multitaper_stream, multitaper_psd, &
        multitaper_spectral_matrix, multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_matrix_result, &
        spectral_result, spectral_status_success, spectral_status_invalid_input
    implicit none
    integer, parameter :: n = 64
    real(dp) :: x(n,2)
    type(multitaper_stream) :: stream
    type(multitaper_spectral_matrix_spectrogram_result) :: expected, actual
    type(multitaper_matrix_result) :: matrix
    type(spectral_result) :: psd
    integer :: i, status

    do i = 1,n
        x(i,1) = real(i,dp)*0.05_dp + sin(0.23_dp*real(i,dp))
        x(i,2) = real(i,dp)*0.13_dp + cos(0.3_dp*real(i,dp))
    end do
    psd = multitaper_psd(x(:,1),2.5_dp,kmax=4,detrend_mode='linear')
    matrix = multitaper_spectral_matrix(x,2.5_dp,kmax=4,detrend_mode='linear')
    if (psd%status /= spectral_status_success .or. matrix%status /= spectral_status_success) &
        error stop 'public whole-record multitaper status'
    if (maxval(abs(psd%power-real(matrix%power(:,1,1),dp))) > 1.0e-10_dp) &
        error stop 'public whole-record PSD/matrix mismatch'
    expected=multitaper_spectral_matrix_spectrogram(x,2.5_dp,16,kmax=4, &
        noverlap=8,nfft=22,detrend_mode='linear')
    call stream%initialize(2.5_dp,16,2,status,noverlap=8,kmax=4,nfft=22, &
        detrend_mode='linear')
    if (status /= spectral_status_success) error stop 'public stream init'
    actual=stream%push(x)
    if (actual%status /= spectral_status_success .or. &
        maxval(abs(actual%power-expected%power)) > 1.0e-10_dp) &
        error stop 'public streaming multitaper mismatch'
    call stream%initialize(2.5_dp,16,2,status,detrend_mode='quadratic')
    if (status /= spectral_status_invalid_input) error stop 'public invalid detrend mode'
    print '(a)', 'Public multitaper linear whole-record/stream API passed'
end program test_signal_multitaper_whole_stream_public
