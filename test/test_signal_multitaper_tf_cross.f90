! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_tf_cross
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result, multitaper_cross_spectrogram, &
        multitaper_cross_spectrogram_result, multitaper_tf_coherence, &
        multitaper_tf_coherence_result, multitaper_spectrogram, &
        multitaper_spectrogram_channels_result
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: data(96,3), x(96), y(96), zeros(96)
    complex(dp) :: cmplx_data(96,3), xc(96), yc(96)
    type(multitaper_spectral_matrix_spectrogram_result) :: mat, other
    type(multitaper_cross_spectrogram_result) :: cross, mixed
    type(multitaper_tf_coherence_result) :: coh
    type(multitaper_spectrogram_channels_result) :: spec
    integer :: i, c, d

    do i = 1, 96
        x(i) = sin(2.0_dp*pi*real(i,dp)/13.0_dp)+0.3_dp*cos(2.0_dp*pi*real(i,dp)/7.0_dp)
        y(i) = 0.75_dp*x(i)+0.5_dp*sin(2.0_dp*pi*real(i,dp)/9.0_dp)
        xc(i) = cmplx(x(i),0.2_dp*y(i),dp)
        yc(i) = cmplx(y(i),-0.1_dp*x(i),dp)
    end do
    zeros = 0.0_dp
    data(:,1) = x
    data(:,2) = y
    data(:,3) = 2.0_dp*x
    cmplx_data(:,1) = xc
    cmplx_data(:,2) = yc
    cmplx_data(:,3) = 2.0_dp*xc

    mat = multitaper_spectral_matrix_spectrogram(data,2.5_dp,24,noverlap=12,kmax=4, &
        fs=5.0_dp,nfft=31,weighting='adaptive')
    spec = multitaper_spectrogram(data,2.5_dp,24,noverlap=12,kmax=4, &
        fs=5.0_dp,nfft=31,weighting='adaptive')
    if (mat%status /= spectral_status_success) error stop 'real matrix status'
    if (size(mat%power,1) /= 16 .or. size(mat%power,2) /= 7) error stop 'real matrix shape'
    if (size(mat%power,3) /= 3 .or. size(mat%power,4) /= 3) error stop 'real matrix channels'
    if (abs(mat%times(1)-12.0_dp/5.0_dp) > 1.0e-14_dp) error stop 'real matrix time'
    if (maxval(abs(mat%frequencies-spec%frequencies)) > 1.0e-14_dp) error stop 'real matrix freq'
    do c = 1, 3
        if (maxval(abs(real(mat%power(:,:,c,c))-spec%power(:,:,c))) > 1.0e-11_dp) &
            error stop 'PSD/spectrogram mismatch'
        do d = 1, 3
            if (maxval(abs(mat%power(:,:,c,d)-conjg(mat%power(:,:,d,c)))) > 1.0e-12_dp) &
                error stop 'real matrix not Hermitian'
        end do
    end do
    other = multitaper_spectral_matrix_spectrogram(transpose(data),2.5_dp,24,noverlap=12,kmax=4, &
        fs=5.0_dp,nfft=31,weighting='adaptive',axis=2)
    if (other%status /= spectral_status_success) error stop 'transposed matrix status'
    if (maxval(abs(other%power-mat%power)) > 1.0e-12_dp) error stop 'axis=2 mismatch'

    cross = multitaper_cross_spectrogram(x,y,2.5_dp,24,noverlap=12,kmax=4, &
        fs=5.0_dp,nfft=31,weighting='adaptive')
    coh = multitaper_tf_coherence(x,y,2.5_dp,24,noverlap=12,kmax=4, &
        fs=5.0_dp,nfft=31,weighting='adaptive')
    if (cross%status /= spectral_status_success .or. coh%status /= spectral_status_success) &
        error stop 'pair status'
    if (maxval(abs(cross%power-mat%power(:,:,1,2))) > 1.0e-12_dp) error stop 'cross/matrix mismatch'
    if (maxval(abs(cross%auto_x-real(mat%power(:,:,1,1)))) > 1.0e-12_dp) error stop 'Pxx mismatch'
    if (maxval(abs(cross%auto_y-real(mat%power(:,:,2,2)))) > 1.0e-12_dp) error stop 'Pyy mismatch'
    if (any(coh%coherence < -1.0e-14_dp) .or. any(coh%coherence > 1.0_dp+1.0e-14_dp)) &
        error stop 'coherence out of bounds'

    coh = multitaper_tf_coherence(x,2.0_dp*x,2.5_dp,24,kmax=4,weighting='adaptive')
    if (coh%status /= spectral_status_success) error stop 'perfect coherence status'
    if (maxval(abs(coh%coherence-1.0_dp)) > 1.0e-10_dp) error stop 'perfect coherence'
    coh = multitaper_tf_coherence(zeros,y,2.5_dp,24,kmax=4)
    if (coh%status /= spectral_status_success .or. any(coh%coherence /= 0.0_dp)) &
        error stop 'zero-power coherence'

    mat = multitaper_spectral_matrix_spectrogram(cmplx_data,2.5_dp,24,noverlap=14, &
        kmax=4,weighting='eigen',nfft=32)
    if (mat%status /= spectral_status_success) error stop 'complex matrix status'
    if (size(mat%frequencies) /= 32 .or. mat%frequencies(32) >= 0.0_dp) &
        error stop 'complex frequencies'
    do c = 1, 3
        do d = 1, 3
            if (maxval(abs(mat%power(:,:,c,d)-conjg(mat%power(:,:,d,c)))) > 1.0e-12_dp) &
                error stop 'complex matrix not Hermitian'
        end do
    end do
    cross = multitaper_cross_spectrogram(xc,yc,2.5_dp,24,noverlap=14, &
        kmax=4,weighting='eigen',nfft=32)
    if (maxval(abs(cross%power-mat%power(:,:,1,2))) > 1.0e-12_dp) error stop 'complex cross/matrix mismatch'
    mixed = multitaper_cross_spectrogram(x,yc,2.5_dp,24,nfft=32)
    if (mixed%status /= spectral_status_success .or. size(mixed%frequencies) /= 32) &
        error stop 'mixed spectrum must be two-sided'
    mixed = multitaper_cross_spectrogram(xc,y,2.5_dp,24,nfft=32)
    if (mixed%status /= spectral_status_success) error stop 'complex/real status'
    coh = multitaper_tf_coherence(xc,y,2.5_dp,24,nfft=32)
    if (coh%status /= spectral_status_success) error stop 'complex/real coherence status'
    coh = multitaper_tf_coherence(x,yc,2.5_dp,24,nfft=32)
    if (coh%status /= spectral_status_success) error stop 'real/complex coherence status'

    cross = multitaper_cross_spectrogram(x,y(1:95),2.5_dp,24)
    if (cross%status /= spectral_status_invalid_input) error stop 'different input lengths accepted'
    cross = multitaper_cross_spectrogram(x,y,2.5_dp,24,noverlap=24)
    if (cross%status /= spectral_status_invalid_input) error stop 'invalid overlap accepted'
    coh = multitaper_tf_coherence(x,y,2.5_dp,24,weighting='bogus')
    if (coh%status /= spectral_status_invalid_input) error stop 'invalid weighting accepted'
    mat = multitaper_spectral_matrix_spectrogram(data,2.5_dp,24,axis=3)
    if (mat%status /= spectral_status_invalid_input) error stop 'invalid axis accepted'
    print '(a)', 'Multitaper cross-spectrogram and time-frequency coherence passed'
end program test_signal_multitaper_tf_cross
