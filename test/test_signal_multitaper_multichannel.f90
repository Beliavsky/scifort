! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_multichannel
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_psd, &
        multitaper_csd, multitaper_coherence, multitaper_spectral_matrix, &
        multitaper_channels_result, multitaper_matrix_result
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        spectral_status_success, spectral_status_invalid_input
    implicit none

    real(dp) :: x(64,3), p, pi
    complex(dp) :: z(64,3)
    type(spectral_result) :: real_psd, complex_psd, coher
    type(cross_spectral_result) :: cross, mixed, inverse_cross
    type(multitaper_matrix_result) :: matrix, complex_matrix, invalid
    type(multitaper_channels_result) :: channels, rotated, complex_channels
    integer :: i, j, k

    pi = acos(-1.0_dp)
    do i = 1, 64
        p = real(i-1,dp)/64.0_dp
        x(i,1) = sin(2.0_dp*pi*8.0_dp*p) + 0.1_dp*cos(2.0_dp*pi*15.0_dp*p)
        x(i,2) = 2.0_dp*x(i,1)
        x(i,3) = cos(2.0_dp*pi*11.0_dp*p)
    end do
    z = cmplx(x,0.0_dp,dp)
    z(:,1) = cmplx(x(:,1),x(:,3),dp)
    z(:,2) = (2.0_dp,-1.0_dp)*z(:,1)

    matrix = multitaper_spectral_matrix(x,3.5_dp,kmax=5,fs=64.0_dp,nfft=128, &
        detrend=.false.,weighting='eigen')
    if (matrix%status /= spectral_status_success) error stop 'real spectral matrix status'
    if (any(shape(matrix%power) /= [65,3,3])) error stop 'real spectral matrix shape'
    if (maxval(abs(aimag(matrix%power(:,1,1)))) > 1.0e-12_dp) error stop 'auto spectrum imaginary'
    if (maxval(abs(matrix%power(:,1,2)-conjg(matrix%power(:,2,1)))) > 1.0e-12_dp) &
        error stop 'Hermitian spectral matrix'
    if (maxval(abs(matrix%power(:,2,2)-4.0_dp*matrix%power(:,1,1))) > 1.0e-11_dp) &
        error stop 'real scale invariance'

    channels = multitaper_psd(x,3.5_dp,kmax=5,fs=64.0_dp,nfft=128, &
        detrend=.false.,weighting='eigen')
    if (channels%status /= spectral_status_success) error stop 'multichannel PSD status'
    if (any(shape(channels%power) /= [65,3])) error stop 'multichannel PSD shape'
    do j = 1, 3
        if (maxval(abs(channels%power(:,j)-real(matrix%power(:,j,j),dp))) > 1.0e-12_dp) &
            error stop 'PSD diagonal'
    end do
    rotated = multitaper_psd(transpose(x),3.5_dp,kmax=5,fs=64.0_dp, &
        nfft=128,detrend=.false.,weighting='eigen',axis=2)
    if (maxval(abs(channels%power-rotated%power)) > 1.0e-12_dp) error stop 'axis 2 PSD'
    real_psd = multitaper_psd(x(:,1),3.5_dp,kmax=5,fs=64.0_dp, &
        nfft=128,detrend=.false.,weighting='eigen')
    if (maxval(abs(real_psd%power-channels%power(:,1))) > 1.0e-12_dp) &
        error stop 'legacy 1-D PSD compatibility'
    cross = multitaper_csd(x(:,1),x(:,2),3.5_dp,kmax=5,fs=64.0_dp, &
        nfft=128,detrend=.false.,weighting='eigen')
    if (cross%status /= spectral_status_success) error stop 'real CSD status'
    if (maxval(abs(cross%power-matrix%power(:,1,2))) > 1.0e-12_dp) error stop 'CSD matrix component'
    coher = multitaper_coherence(x(:,1),x(:,2),3.5_dp,kmax=5, &
        nfft=128,detrend=.false.)
    if (coher%status /= spectral_status_success) error stop 'coherence status'
    if (maxval(abs(coher%power-1.0_dp)) > 1.0e-10_dp) error stop 'identical-signal coherence'

    complex_matrix = multitaper_spectral_matrix(z,3.5_dp,kmax=5,fs=64.0_dp, &
        nfft=127,weighting='equal')
    if (complex_matrix%status /= spectral_status_success) error stop 'complex spectral matrix status'
    if (any(shape(complex_matrix%power) /= [127,3,3])) error stop 'complex spectral matrix shape'
    if (complex_matrix%frequencies(65) >= 0.0_dp) error stop 'complex two-sided fftfreq'
    if (maxval(abs(complex_matrix%power(:,2,1)-conjg(complex_matrix%power(:,1,2)))) > 1.0e-11_dp) &
        error stop 'complex spectral Hermitian'
    complex_channels = multitaper_psd(z,3.5_dp,kmax=5,fs=64.0_dp,nfft=127)
    if (complex_channels%status /= spectral_status_success) error stop 'complex channels'
    if (maxval(abs(complex_channels%power(:,1)-real(complex_matrix%power(:,1,1),dp))) > 1.0e-11_dp) &
        error stop 'complex diagonal PSD'
    complex_psd = multitaper_psd(z(:,1),3.5_dp,kmax=5,fs=64.0_dp,nfft=127)
    if (maxval(abs(complex_psd%power-complex_channels%power(:,1))) > 1.0e-11_dp) &
        error stop 'complex 1-D PSD'
    cross = multitaper_csd(z(:,1),z(:,2),3.5_dp,kmax=5,fs=64.0_dp,nfft=127)
    inverse_cross = multitaper_csd(z(:,2),z(:,1),3.5_dp,kmax=5,fs=64.0_dp,nfft=127)
    if (maxval(abs(cross%power-conjg(inverse_cross%power))) > 1.0e-11_dp) &
        error stop 'reverse CSD conjugation'
    if (maxval(abs(cross%power-complex_matrix%power(:,1,2))) > 1.0e-11_dp) &
        error stop 'complex CSD matrix component'
    mixed = multitaper_csd(x(:,1),z(:,1),3.5_dp,kmax=5,nfft=127)
    inverse_cross = multitaper_csd(z(:,1),x(:,1),3.5_dp,kmax=5,nfft=127)
    if (mixed%status /= spectral_status_success .or. size(mixed%power) /= 127) &
        error stop 'mixed-type CSD status'
    if (maxval(abs(mixed%power-conjg(inverse_cross%power))) > 1.0e-11_dp) &
        error stop 'mixed-type CSD conjugation'

    coher = multitaper_coherence(0.0_dp*x(:,1),0.0_dp*x(:,2),3.5_dp,kmax=3)
    if (coher%status /= spectral_status_success) error stop 'zero coherence status'
    if (any(coher%power /= 0.0_dp)) error stop 'zero coherence output'
    channels = multitaper_psd(x(:,1:3:2),3.5_dp,kmax=3)
    if (channels%status /= spectral_status_success .or. size(channels%power,2) /= 2) &
        error stop 'noncontiguous channels'
    complex_channels = multitaper_psd(z(1:63:2,1:3:2),2.5_dp,kmax=3,axis=1)
    if (complex_channels%status /= spectral_status_success .or. size(complex_channels%power,2) /= 2) &
        error stop 'strided complex samples'
    invalid = multitaper_spectral_matrix(x,3.5_dp,axis=3)
    if (invalid%status /= spectral_status_invalid_input) error stop 'invalid axis'
    if (size(invalid%power) /= 0) error stop 'invalid axis result shape'
    invalid = multitaper_spectral_matrix(x,3.5_dp,kmax=0)
    if (invalid%status /= spectral_status_invalid_input) error stop 'invalid taper count'
    invalid = multitaper_spectral_matrix(x,3.5_dp,nfft=16)
    if (invalid%status /= spectral_status_invalid_input) error stop 'invalid nfft'
    cross = multitaper_csd(x(:,1),x(:32,2),3.5_dp)
    if (cross%status /= spectral_status_invalid_input) error stop 'mismatched signal lengths'
    coher = multitaper_coherence(x(:,1),x(:32,2),3.5_dp)
    if (coher%status /= spectral_status_invalid_input) error stop 'mismatched coherence lengths'
    print '(a)', 'Multitaper multichannel, complex and CSD regression passed'
end program test_signal_multitaper_multichannel
