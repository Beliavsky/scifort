! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_adaptive
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_psd, multitaper_psd_ci, &
        multitaper_spectral_matrix, multitaper_csd, multitaper_coherence, &
        multitaper_confidence_result, multitaper_channels_result, multitaper_matrix_result
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        spectral_status_success, spectral_status_invalid_input
    implicit none

    real(dp) :: x(80), y(80), data(80,3), pi, phase
    complex(dp) :: z(80), matrix_data(80,3)
    integer :: i
    type(multitaper_confidence_result) :: ci, ci2, empty, zeros, complex_ci
    type(spectral_result) :: spectrum, psd, coherence
    type(multitaper_channels_result) :: channels, turned
    type(multitaper_matrix_result) :: matrix, matrix_z
    type(cross_spectral_result) :: cross

    pi = acos(-1.0_dp)
    do i = 1, 80
        phase = 2.0_dp*pi*real(i-1,dp)/80.0_dp
        x(i) = sin(7.0_dp*phase)+0.33_dp*sin(19.0_dp*phase)+0.11_dp*cos(37.0_dp*phase)
        y(i) = 0.2_dp*sin(8.0_dp*phase)+0.7_dp*x(i)
    end do
    z = cmplx(x,y,dp)
    data(:,1) = x
    data(:,2) = y
    data(:,3) = 3.0_dp*x
    matrix_data(:,1) = z
    matrix_data(:,2) = (0.5_dp,-2.0_dp)*z
    matrix_data(:,3) = cmplx(y,x,dp)

    spectrum = multitaper_psd(x,3.5_dp,kmax=5,fs=80.0_dp,nfft=128, &
        weighting='adaptive',detrend=.false.)
    ci = multitaper_psd_ci(x,3.5_dp,kmax=5,fs=80.0_dp,nfft=128, &
        weighting='adaptive',detrend=.false.)
    if (spectrum%status /= spectral_status_success) error stop 'adaptive 1-D status'
    if (ci%status /= spectral_status_success) error stop 'jackknife 1-D status'
    if (size(ci%power) /= 65) error stop 'jackknife shape'
    if (maxval(abs(ci%power-spectrum%power)) > 5.0e-12_dp) error stop 'jackknife PSD parity'
    if (any(ci%lower < 0.0_dp) .or. any(ci%upper < ci%lower)) error stop 'CI ordering'
    if (any(ci%power < ci%lower) .or. any(ci%power > ci%upper)) error stop 'CI coverage'
    if (any(ci%stderr_log < 0.0_dp)) error stop 'CI stderr'
    if (any(ci%equivalent_dof < 2.0_dp-1.0e-12_dp) .or. &
        any(ci%equivalent_dof > 10.0_dp+1.0e-12_dp)) error stop 'equivalent DOF'
    ci2 = multitaper_psd_ci(x,3.5_dp,kmax=5,fs=80.0_dp,nfft=128, &
        weighting='adaptive',detrend=.false.,confidence=0.90_dp)
    if (any(ci2%lower < ci%lower-1.0e-12_dp) .or. &
        any(ci2%upper > ci%upper+1.0e-12_dp)) error stop 'CI monotonicity'

    channels = multitaper_psd(data,3.5_dp,kmax=5,fs=80.0_dp, &
        weighting='adaptive',nfft=128,detrend=.false.)
    turned = multitaper_psd(transpose(data),3.5_dp,kmax=5,fs=80.0_dp, &
        weighting='adaptive',nfft=128,detrend=.false.,axis=2)
    matrix = multitaper_spectral_matrix(data,3.5_dp,kmax=5,fs=80.0_dp, &
        weighting='adaptive',nfft=128,detrend=.false.)
    if (channels%status /= spectral_status_success .or. matrix%status /= spectral_status_success) &
        error stop 'adaptive channels status'
    if (maxval(abs(turned%power-channels%power)) > 1.0e-11_dp) error stop 'adaptive axis parity'
    if (maxval(abs(real(matrix%power(:,1,1),dp)-spectrum%power)) > 1.0e-11_dp) &
        error stop 'adaptive matrix 1-D parity'
    if (maxval(abs(real(matrix%power(:,3,3),dp)-9.0_dp*spectrum%power)) > 1.0e-10_dp) &
        error stop 'adaptive amplitude scaling'
    if (maxval(abs(matrix%power(:,1,2)-conjg(matrix%power(:,2,1)))) > 1.0e-10_dp) &
        error stop 'Hermitian adaptive matrix'
    if (maxval(abs(real(matrix%power(:,1,1),dp)-channels%power(:,1))) > 1.0e-11_dp) &
        error stop 'adaptive diagonal parity'
    cross = multitaper_csd(x,y,3.5_dp,kmax=5,nfft=128,detrend=.false.,weighting='adaptive')
    if (cross%status /= spectral_status_success) error stop 'adaptive CSD status'
    coherence = multitaper_coherence(x,3.0_dp*x,3.5_dp,kmax=5,weighting='adaptive')
    if (coherence%status /= spectral_status_success) error stop 'adaptive coherence status'
    if (maxval(abs(coherence%power-1.0_dp)) > 1.0e-11_dp) error stop 'adaptive coherence scaling'

    matrix_z = multitaper_spectral_matrix(matrix_data,3.5_dp,kmax=5, &
        weighting='adaptive',nfft=127,detrend=.false.)
    if (matrix_z%status /= spectral_status_success) error stop 'complex spectral matrix status'
    if (maxval(abs(matrix_z%power(:,1,2)-conjg(matrix_z%power(:,2,1)))) > 1.0e-10_dp) &
        error stop 'complex Hermitian adaptive matrix'
    complex_ci = multitaper_psd_ci(z,3.5_dp,kmax=5,nfft=127, &
        weighting='adaptive',detrend=.false.)
    psd = multitaper_psd(z,3.5_dp,kmax=5,nfft=127, &
        weighting='adaptive',detrend=.false.)
    if (complex_ci%status /= spectral_status_success) error stop 'complex CI status'
    if (maxval(abs(complex_ci%power-psd%power)) > 1.0e-11_dp) error stop 'complex CI PSD'
    if (maxval(abs(psd%power-real(matrix_z%power(:,1,1),dp))) > 1.0e-11_dp) &
        error stop 'complex PSD matrix diagonal'

    zeros = multitaper_psd_ci(0.0_dp*x,3.5_dp,kmax=5,weighting='adaptive')
    if (zeros%status /= spectral_status_success) error stop 'zero input CI status'
    if (any(zeros%power /= 0.0_dp) .or. any(zeros%upper /= 0.0_dp)) error stop 'zero input CI'
    empty = multitaper_psd_ci(x,3.5_dp,kmax=1)
    if (empty%status /= spectral_status_invalid_input) error stop 'jackknife one taper rejected'
    empty = multitaper_psd_ci(x,3.5_dp,kmax=5,confidence=1.0_dp)
    if (empty%status /= spectral_status_invalid_input) error stop 'invalid confidence'
    empty = multitaper_psd_ci(x,3.5_dp,kmax=5,weighting='bogus')
    if (empty%status /= spectral_status_invalid_input) error stop 'invalid weighting'
    ci2 = multitaper_psd_ci(x,3.5_dp,kmax=2,weighting='adaptive')
    if (ci2%status /= spectral_status_success) error stop 'two-taper fallback'
    print '(a)', 'Adaptive multitaper and jackknife regression passed'
end program test_signal_multitaper_adaptive
