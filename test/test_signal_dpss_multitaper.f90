! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_dpss_multitaper
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result, &
        dpss_status_success, dpss_status_invalid_input
    use scifort_signal_multitaper, only : multitaper_psd
    use scifort_signal_spectral, only : spectral_result, spectral_status_success, &
        spectral_status_invalid_input
    implicit none
    type(dpss_result) :: result, periodic, single, invalid
    type(spectral_result) :: psd, eig, zero, bad
    real(dp) :: x(64), pi, delta, norm_diff
    integer :: i, j, k

    pi = acos(-1.0_dp)
    result = dpss(64,3.5_dp,kmax=6)
    if (result%status /= dpss_status_success) error stop 'dpss status'
    if (any(shape(result%tapers) /= [64,6])) error stop 'dpss shape'
    if (size(result%ratios) /= 6) error stop 'dpss concentration shape'
    do j = 1, 6
        if (abs(sum(result%tapers(:,j)**2)-1.0_dp) > 3.0e-12_dp) error stop 'orthonormal norm'
        do k = 1, j-1
            if (abs(dot_product(result%tapers(:,j),result%tapers(:,k))) > 3.0e-12_dp) &
                error stop 'orthogonality'
        end do
        if (result%ratios(j) < 0.0_dp .or. result%ratios(j) > 1.0_dp) &
            error stop 'concentration bounds'
    end do
    if (any(result%ratios(1:5) < result%ratios(2:6))) &
        error stop 'ordered spectral concentrations'
    do j = 1, 6
        if (mod(j,2) == 1) then
            if (maxval(abs(result%tapers(:,j)-result%tapers(64:1:-1,j))) > 1.0e-10_dp) &
                error stop 'even-parity taper'
        else
            if (maxval(abs(result%tapers(:,j)+result%tapers(64:1:-1,j))) > 1.0e-10_dp) &
                error stop 'odd-parity taper'
        end if
    end do
    single = dpss(64,3.5_dp)
    if (single%status /= dpss_status_success) error stop 'default dpss status'
    if (maxval(single%tapers(:,1)) > 1.0_dp + 1.0e-12_dp) error stop 'default normalization'
    if (abs(maxval(single%tapers(:,1))-64.0_dp**2/(64.0_dp**2+3.5_dp)) > 1.0e-12_dp) &
        error stop 'even-length approximate correction'
    periodic = dpss(64,3.5_dp,kmax=6,sym=.false.)
    if (periodic%status /= dpss_status_success) error stop 'periodic status'
    norm_diff = maxval(abs(periodic%tapers(:,1)-result%tapers(:,1)))
    if (norm_diff < 1.0e-5_dp) error stop 'periodic must differ from symmetric'

    invalid = dpss(16,-1.0_dp,kmax=2)
    if (invalid%status /= dpss_status_invalid_input) error stop 'invalid NW'
    invalid = dpss(16,8.0_dp,kmax=2)
    if (invalid%status /= dpss_status_invalid_input) error stop 'Nyquist NW'
    invalid = dpss(16,2.5_dp,kmax=17)
    if (invalid%status /= dpss_status_invalid_input) error stop 'invalid k'
    invalid = dpss(16,2.5_dp,kmax=3,norm='invalid')
    if (invalid%status /= dpss_status_invalid_input) error stop 'invalid normalization'
    invalid = dpss(1,0.3_dp)
    if (invalid%status /= dpss_status_success) error stop 'length 1'
    if (invalid%tapers(1,1) /= 1.0_dp) error stop 'length 1 value'
    invalid = dpss(0,0.3_dp)
    if (invalid%status /= dpss_status_success) error stop 'length 0'

    x = [(cos(2.0_dp*pi*real(i-1,dp)*7.0_dp/64.0_dp),i=1,64)]
    psd = multitaper_psd(x,3.5_dp,kmax=5,fs=64.0_dp,nfft=128,detrend=.false.)
    if (psd%status /= spectral_status_success) error stop 'multitaper status'
    if (size(psd%power) /= 65 .or. size(psd%frequencies) /= 65) error stop 'multitaper size'
    if (abs(psd%frequencies(15)-7.0_dp) > 1.0e-12_dp) error stop 'frequency axis'
    if (maxloc(psd%power,dim=1) /= 15) error stop 'sinusoid peak'
    if (any(psd%power < -1.0e-15_dp)) error stop 'negative PSD'
    eig = multitaper_psd(x,3.5_dp,kmax=5,fs=64.0_dp,nfft=128,weighting='eigen')
    if (eig%status /= spectral_status_success) error stop 'eigenweight status'
    if (maxval(abs(eig%power-psd%power)) < 1.0e-9_dp) &
        error stop 'eigen weighting must affect estimate'
    zero = multitaper_psd(spread(3.0_dp,1,64),2.5_dp,kmax=3)
    if (zero%status /= spectral_status_success .or. maxval(abs(zero%power)) > 1.0e-20_dp) &
        error stop 'constant detrending'
    zero = multitaper_psd(spread(1.0_dp,1,64),2.5_dp,kmax=3,detrend=.false.)
    delta = zero%frequencies(2)-zero%frequencies(1)
    if (abs(sum(zero%power)*delta-1.0_dp) > 1.0e-10_dp) error stop 'Parseval constant'
    bad = multitaper_psd(x,3.5_dp,fs=-1.0_dp)
    if (bad%status /= spectral_status_invalid_input) error stop 'invalid fs'
    bad = multitaper_psd(x,3.5_dp,nfft=32)
    if (bad%status /= spectral_status_invalid_input) error stop 'invalid FFT length'
    bad = multitaper_psd(x,3.5_dp,kmax=0)
    if (bad%status /= spectral_status_invalid_input) error stop 'invalid taper count'
    bad = multitaper_psd(x,3.5_dp,weighting='bogus')
    if (bad%status /= spectral_status_invalid_input) error stop 'invalid weighting'
    print '(a)', 'DPSS and multitaper regression passed'
end program test_signal_dpss_multitaper
