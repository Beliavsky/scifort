! SPDX-License-Identifier: MIT
! Public generic interface smoke test (also compiles with -fsyntax-only).
program test_signal_mixed_cross_public
    use scifort_signal, only : dp, csd, coherence, short_time_fft, &
        cross_spectral_result, spectral_result, short_fft_success, spectral_status_success
    implicit none
    real(dp) :: x(32)
    complex(dp) :: z(32)
    complex(dp), allocatable :: sxy(:,:)
    type(cross_spectral_result) :: cx, cy
    type(spectral_result) :: coh
    type(short_time_fft) :: plan
    integer :: i, status

    do i=1,32
        x(i) = cos(0.16_dp*i) + 0.2_dp*sin(0.37_dp*i)
        z(i) = cmplx(x(i),sin(0.4_dp*i),dp)
    end do
    cx = csd(x,z,nperseg=12,noverlap=4,nfft=16,window='hann')
    cy = csd(z,x,nperseg=12,noverlap=4,nfft=16,window='hann')
    coh = coherence(x,z,nperseg=12,noverlap=4)
    if (cx%status /= spectral_status_success .or. cy%status /= spectral_status_success) &
        error stop 'public CSD generic'
    if (coh%status /= spectral_status_success) error stop 'public coherence generic'
    call plan%initialize_window('hann',8,3,12.0_dp,mfft=12,fft_mode='twosided')
    if (plan%status /= short_fft_success) error stop 'public ShortTimeFFT initialization'
    call plan%cross_spectrogram(x,z,sxy,status)
    if (status /= short_fft_success) error stop 'public mixed cross spectrogram'
    print '(a)', 'test_signal_mixed_cross_public passed'
end program test_signal_mixed_cross_public
