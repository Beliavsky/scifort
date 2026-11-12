! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_detrend_harmonic_public
    use scifort_signal, only : dp, detrend, multitaper_harmonic_scan, &
        multitaper_harmonic_scan_result, multitaper_noise_floor, multitaper_noise_floor_result, &
        detrend_status_success, spectral_status_success
    implicit none
    real(dp) :: x(32), channels(32,2)
    real(dp), allocatable :: residual(:,:)
    type(multitaper_harmonic_scan_result) :: detection
    type(multitaper_noise_floor_result) :: floor
    integer :: i, status

    do i=1,32
        x(i) = 2.0_dp+0.5_dp*real(i,dp)
    end do
    channels(:,1) = x
    channels(:,2) = 3.0_dp*x
    residual = detrend(channels,axis=1,status=status)
    if (status /= detrend_status_success) error stop 'public detrend status'
    if (maxval(abs(residual)) > 1.0e-11_dp) error stop 'public detrend values'
    detection = multitaper_harmonic_scan(x,2.5_dp,kmax=4,method='bh')
    if (detection%status /= spectral_status_success) error stop 'public harmonic status'
    if (size(detection%adjusted_pvalue) /= 17) error stop 'public harmonic length'
    if (detection%num_tests /= 15) error stop 'public harmonic test count'
    detection = multitaper_harmonic_scan(x,2.5_dp,kmax=4,method='holm')
    if (detection%status /= spectral_status_success) error stop 'public Holm status'
    floor = multitaper_noise_floor(x,2.5_dp,5,kmax=4)
    if (floor%status /= spectral_status_success) error stop 'public floor status'
    if (size(floor%background) /= 17) error stop 'public floor shape'
    print '(a)', 'test_signal_detrend_harmonic_public passed'
end program test_signal_detrend_harmonic_public
