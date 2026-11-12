! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_multichannel_public
    use scifort_signal
    implicit none
    real(dp) :: samples(32,2)
    complex(dp) :: z(32,2)
    type(multitaper_matrix_result) :: matrix
    type(multitaper_channels_result) :: channels
    type(cross_spectral_result) :: cross
    type(spectral_result) :: coh
    type(multitaper_confidence_result) :: uncertainty
    integer :: i

    do i = 1, 32
        samples(i,1) = sin(2.0_dp*acos(-1.0_dp)*real(i-1,dp)*4.0_dp/32.0_dp)
        samples(i,2) = 0.5_dp*samples(i,1)
    end do
    z = cmplx(samples,0.0_dp,dp)
    z(:,1) = cmplx(samples(:,1),samples(:,2),dp)
    matrix = multitaper_spectral_matrix(samples,2.5_dp,kmax=3)
    if (matrix%status /= spectral_status_success) error stop 'public real matrix'
    if (any(shape(matrix%power) /= [17,2,2])) error stop 'public real shape'
    channels = multitaper_psd(z,2.5_dp,kmax=3)
    if (channels%status /= spectral_status_success) error stop 'public complex channels'
    if (any(shape(channels%power) /= [32,2])) error stop 'public complex shape'
    cross = multitaper_csd(samples(:,1),samples(:,2),2.5_dp,kmax=3)
    if (cross%status /= spectral_status_success) error stop 'public CSD'
    coh = multitaper_coherence(samples(:,1),samples(:,2),2.5_dp,kmax=3)
    if (coh%status /= spectral_status_success) error stop 'public coherence'
    if (maxval(abs(coh%power-1.0_dp)) > 1.0e-10_dp) error stop 'public coherence match'
    uncertainty = multitaper_psd_ci(samples(:,1),2.5_dp,kmax=3,weighting='adaptive')
    if (uncertainty%status /= spectral_status_success) error stop 'public adaptive CI'
    if (size(uncertainty%upper) /= 17) error stop 'public adaptive CI shape'
    print '(a)', 'Public multichannel multitaper interfaces passed'
end program test_signal_multitaper_multichannel_public
