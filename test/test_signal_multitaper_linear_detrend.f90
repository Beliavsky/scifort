! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_linear_detrend
    use scifort_kinds, only : dp
    use scifort_signal, only : multitaper_spectrogram, multitaper_spectral_matrix_spectrogram, &
        multitaper_cross_spectrogram, multitaper_tf_coherence, &
        multitaper_spectrogram_result, multitaper_spectrogram_channels_result, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_cross_spectrogram_result, &
        multitaper_tf_coherence_result, spectral_status_success, spectral_status_invalid_input
    implicit none
    integer, parameter :: n = 96, nperseg = 32
    real(dp) :: x(n), y(n), trend(n), noise(n), signals(n,2), pi
    complex(dp) :: z(n), zsignals(n,2)
    type(multitaper_spectrogram_result) :: p, q, r, bad
    type(multitaper_spectrogram_channels_result) :: channels, transposed, complex_channels
    type(multitaper_spectral_matrix_spectrogram_result) :: matrix, complex_matrix
    type(multitaper_cross_spectrogram_result) :: cross, mixed
    type(multitaper_tf_coherence_result) :: coherence
    integer :: i

    pi = acos(-1.0_dp)
    do i = 1,n
        trend(i) = 8.0_dp+0.12_dp*real(i,dp)
        noise(i) = sin(2.0_dp*pi*real(i,dp)/9.0_dp)+0.31_dp*cos(2.0_dp*pi*real(i,dp)/13.0_dp)
        x(i) = noise(i)+trend(i)
        y(i) = 2.0_dp*cos(2.0_dp*pi*real(i,dp)/15.0_dp)-0.6_dp*trend(i)
        z(i) = cmplx(x(i),y(i),dp)
    end do
    signals(:,1) = x
    signals(:,2) = y
    zsignals(:,1) = z
    zsignals(:,2) = cmplx(y,x,dp)

    p = multitaper_spectrogram(trend,3.0_dp,nperseg, &
        noverlap=16,kmax=4,fs=80.0_dp,nfft=48,detrend_mode='linear')
    if (p%status /= spectral_status_success) error stop 'real linear detrend status'
    if (maxval(abs(p%power)) > 1.0e-22_dp) error stop 'linear trend not removed'
    q = multitaper_spectrogram(trend,3.0_dp,nperseg, &
        noverlap=16,kmax=4,fs=80.0_dp,nfft=48,detrend_mode='constant')
    if (q%status /= spectral_status_success) error stop 'constant detrend status'
    if (maxval(q%power) < 1.0e-3_dp) error stop 'constant detrend wrongly removed slope'

    p = multitaper_spectrogram(x,3.0_dp,nperseg,noverlap=16,kmax=4, &
        detrend_mode='linear',weighting='eigen',nfft=48)
    q = multitaper_spectrogram(noise,3.0_dp,nperseg,noverlap=16,kmax=4, &
        detrend_mode='linear',weighting='eigen',nfft=48)
    if (p%status /= spectral_status_success .or. q%status /= spectral_status_success) &
        error stop 'comparison status'
    if (maxval(abs(p%power-q%power)) > 2.0e-11_dp) error stop 'linear detrend shift invariance'
    r = multitaper_spectrogram(x,3.0_dp,nperseg,noverlap=16,kmax=4, &
        detrend=.false.,detrend_mode='linear',weighting='eigen',nfft=48)
    if (maxval(abs(r%power-p%power)) > 1.0e-13_dp) error stop 'detrend mode precedence'
    r = multitaper_spectrogram(x,3.0_dp,nperseg,noverlap=16,kmax=4, &
        detrend=.true.,detrend_mode='constant',weighting='eigen',nfft=48)
    q = multitaper_spectrogram(x,3.0_dp,nperseg,noverlap=16,kmax=4, &
        detrend=.true.,weighting='eigen',nfft=48)
    if (maxval(abs(r%power-q%power)) > 1.0e-13_dp) error stop 'default constant regression'
    bad = multitaper_spectrogram(x,3.0_dp,nperseg,detrend_mode='quadratic')
    if (bad%status /= spectral_status_invalid_input) error stop 'unknown detrend mode accepted'
    if (size(bad%power) /= 0) error stop 'invalid mode array shape'

    channels = multitaper_spectrogram(signals,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=48,weighting='equal')
    transposed = multitaper_spectrogram(transpose(signals),3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=48,weighting='equal',axis=2)
    if (channels%status /= spectral_status_success .or. &
        transposed%status /= spectral_status_success) error stop 'channel detrend status'
    if (maxval(abs(channels%power-transposed%power)) > 1.0e-12_dp) &
        error stop 'channel axis parity'

    matrix = multitaper_spectral_matrix_spectrogram(signals,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=48,weighting='equal')
    if (matrix%status /= spectral_status_success) error stop 'matrix detrend status'
    if (maxval(abs(real(matrix%power(:,:,1,1),dp)-channels%power(:,:,1))) > 1.0e-12_dp) &
        error stop 'matrix PSD diagonal mismatch'
    if (maxval(abs(matrix%power(:,:,1,2)-conjg(matrix%power(:,:,2,1)))) > 1.0e-12_dp) &
        error stop 'matrix not Hermitian'
    cross = multitaper_cross_spectrogram(x,y,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=48,weighting='equal')
    if (cross%status /= spectral_status_success) error stop 'cross detrend status'
    if (maxval(abs(cross%power-matrix%power(:,:,1,2))) > 1.0e-12_dp) &
        error stop 'cross matrix mismatch'
    coherence = multitaper_tf_coherence(x,y,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=48,weighting='equal')
    if (coherence%status /= spectral_status_success) error stop 'coherence detrend status'
    if (any(coherence%coherence < 0.0_dp) .or. any(coherence%coherence > 1.0_dp)) &
        error stop 'coherence range'

    p = multitaper_spectrogram(z,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=49,weighting='equal')
    complex_channels = multitaper_spectrogram(zsignals,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=49,weighting='equal')
    if (maxval(abs(p%power-complex_channels%power(:,:,1))) > 1.0e-12_dp) &
        error stop 'complex channels parity'
    complex_matrix = multitaper_spectral_matrix_spectrogram(zsignals,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=49,weighting='equal')
    if (complex_matrix%status /= spectral_status_success) error stop 'complex spectral matrix status'
    if (maxval(abs(real(complex_matrix%power(:,:,1,1),dp)-p%power)) > 1.0e-12_dp) &
        error stop 'complex diagonal PSD'
    mixed = multitaper_cross_spectrogram(x,z,3.0_dp,nperseg, &
        detrend_mode='linear',kmax=4,nfft=49,weighting='equal')
    if (mixed%status /= spectral_status_success) error stop 'mixed cross detrend status'

    print '(a)', 'Multitaper linear detrending passed'
end program test_signal_multitaper_linear_detrend
