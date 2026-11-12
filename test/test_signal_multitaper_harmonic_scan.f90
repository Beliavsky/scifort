! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_multitaper_harmonic_scan
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_status_success, spectral_status_invalid_input
    use scifort_signal_multitaper, only : multitaper_harmonic_ftest, multitaper_harmonic_scan, &
        multitaper_harmonic_result, multitaper_harmonic_scan_result
    implicit none
    real(dp) :: x(64), pi
    complex(dp) :: z(64)
    type(multitaper_harmonic_result) :: raw
    type(multitaper_harmonic_scan_result) :: scan
    integer :: i,m
    real(dp) :: expect, c

    pi = acos(-1.0_dp)
    do i=1,64
        x(i) = sin(2.0_dp*pi*real(7*(i-1),dp)/64.0_dp) &
            +0.24_dp*cos(2.0_dp*pi*real(13*(i-1),dp)/64.0_dp)
        z(i) = cmplx(x(i), 0.17_dp*sin(2.0_dp*pi*real(11*(i-1),dp)/64.0_dp),dp)
    end do
    raw = multitaper_harmonic_ftest(x,2.7_dp,kmax=4,nfft=96)
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,nfft=96,method='bonferroni')
    if (scan%status /= spectral_status_success) error stop 'bonf status'
    if (scan%num_tests /= 47) error stop 'bonf tested bins'
    if (scan%adjusted_pvalue(1) /= 1.0_dp) error stop 'DC excluded'
    if (scan%adjusted_pvalue(49) /= 1.0_dp) error stop 'Nyquist excluded'
    if (maxval(abs(scan%pvalue-raw%pvalue)) > 1.0e-14_dp) error stop 'raw changed'
    do i=2,48
        expect = min(1.0_dp,real(scan%num_tests,dp)*raw%pvalue(i))
        if (abs(scan%adjusted_pvalue(i)-expect) > 1.0e-13_dp) error stop 'bonf discrepancy'
    end do
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,nfft=96,method='none',test_bins=[2,4,10])
    if (scan%num_tests /= 3) error stop 'selected count'
    if (scan%adjusted_pvalue(3) /= 1.0_dp) error stop 'selected excluded'
    if (abs(scan%adjusted_pvalue(2)-scan%pvalue(2)) > 1.0e-14_dp) error stop 'none method'
    scan = multitaper_harmonic_scan(z,2.7_dp,kmax=4,nfft=96,method='bh')
    if (scan%status /= spectral_status_success) error stop 'complex BH status'
    if (scan%num_tests /= 95) error stop 'complex BH count'
    if (any(scan%adjusted_pvalue < 0.0_dp .or. scan%adjusted_pvalue > 1.0_dp)) &
        error stop 'BH range'
    scan = multitaper_harmonic_scan(z,2.7_dp,kmax=4,nfft=96,method='by')
    if (scan%status /= spectral_status_success) error stop 'BY status'
    if (scan%num_tests /= 95) error stop 'BY count'
    scan = multitaper_harmonic_scan(z,2.7_dp,kmax=4,nfft=96,method='holm')
    if (scan%status /= spectral_status_success) error stop 'Holm status'
    if (scan%num_tests /= 95) error stop 'Holm count'
    if (any(scan%adjusted_pvalue < 0.0_dp .or. scan%adjusted_pvalue > 1.0_dp)) &
        error stop 'Holm range'
    if (any(scan%significant .and. scan%adjusted_pvalue > 0.05_dp)) error stop 'Holm threshold'
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,nfft=96,method='holm',test_bins=[2,4,10])
    if (scan%status /= spectral_status_success) error stop 'selected Holm status'
    if (scan%num_tests /= 3) error stop 'selected Holm count'
    if (scan%adjusted_pvalue(3) /= 1.0_dp) error stop 'selected Holm excluded'
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,method='bh',test_bins=[2,2])
    if (scan%status /= spectral_status_invalid_input) error stop 'duplicate bins'
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,method='bh',test_bins=[1000])
    if (scan%status /= spectral_status_invalid_input) error stop 'invalid bin'
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,alpha=0.0_dp)
    if (scan%status /= spectral_status_invalid_input) error stop 'invalid alpha'
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,method='oops')
    if (scan%status /= spectral_status_invalid_input) error stop 'invalid method'
    scan = multitaper_harmonic_scan(x,2.7_dp,kmax=4,test_bins=[integer :: ])
    if (scan%status /= spectral_status_success) error stop 'empty selection'
    if (scan%num_tests /= 0) error stop 'empty selection count'
    print '(a)', 'test_signal_multitaper_harmonic_scan passed'
end program test_signal_multitaper_harmonic_scan
