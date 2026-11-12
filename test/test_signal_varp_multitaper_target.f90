! SPDX-License-Identifier: MIT
program test_signal_varp_multitaper_target
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result, dpss_status_success
    use scifort_signal_multitaper, only : multitaper_matrix_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_matrix, &
        canonical_coherence_result
    use scifort_signal_varp_taper_spectrum, only : varp_taper_spectrum_result, &
        varp_expected_taper_spectrum, varp_expected_band_taper_spectrum, &
        varp_spectrum_success, varp_spectrum_invalid_input, varp_spectrum_nonstationary
    use scifort_signal_varp_multitaper_target, only : varp_expected_multitaper_spectrum, &
        varp_expected_band_multitaper_spectrum, &
        varp_expected_multitaper_canonical_coherence, &
        varp_expected_band_multitaper_canonical_coherence
    implicit none

    type(varp_taper_spectrum_result) :: out, one, integrated, direct, bad
    type(canonical_coherence_result) :: mode, reference
    type(multitaper_matrix_result) :: matrix
    type(dpss_result) :: bank
    real(dp) :: a(5,5,3), sigma(5,5), taper(13,3), frequency(5), bands(2,4)
    real(dp) :: weights(3), expected_variance, nan_value
    complex(dp), allocatable :: expected(:,:,:)
    integer :: i, j, k, q

    a = 0.0_dp
    sigma = 0.0_dp
    do i = 1, 5
        a(i,i,1) = 0.27_dp
        a(i,i,2) = -0.13_dp
        a(i,i,3) = 0.045_dp
        sigma(i,i) = 1.0_dp+0.1_dp*i
    end do
    a(1,4,1) = 0.14_dp
    a(2,3,2) = -0.21_dp
    a(4,1,3) = 0.09_dp
    sigma(1,3) = 0.48_dp
    sigma(3,1) = sigma(1,3)
    sigma(2,4) = -0.27_dp
    sigma(4,2) = sigma(2,4)
    sigma(5,1) = 0.22_dp
    sigma(1,5) = sigma(5,1)
    frequency = [-0.37_dp,0.0_dp,0.21_dp,0.49_dp,1.21_dp]
    bands(:,1) = [-0.5_dp,0.0_dp]
    bands(:,2) = [0.0_dp,0.5_dp]
    bands(:,3) = [-0.5_dp,0.5_dp]
    bands(:,4) = [0.18_dp,0.18001_dp]
    weights = [1.0_dp,3.0_dp,2.0_dp]
    do j = 1, 3
        do i = 1, 13
            taper(i,j) = (0.6_dp-0.3_dp*cos(2.0_dp*acos(-1.0_dp)* &
                real(i-1,dp)/13.0_dp))*real(j,dp)+ &
                0.03_dp*sin(real(i*j,dp))
        end do
    end do

    ! A one-column bank must agree with the existing finite-taper target.
    out = varp_expected_multitaper_spectrum(a,sigma,taper(:,1:1),frequency,normalized=.true.)
    one = varp_expected_taper_spectrum(a,sigma,taper(:,1),frequency,normalized=.true.)
    call check(out%status == varp_spectrum_success,'one-taper spectrum accepted')
    call check(maxval(abs(out%power-one%power))<1.0e-13_dp,'one-taper exact agreement')
    integrated = varp_expected_band_multitaper_spectrum(a,sigma,taper(:,1:1), &
        bands,average=.true.,normalized=.true.)
    direct = varp_expected_band_taper_spectrum(a,sigma,taper(:,1),bands, &
        average=.true.,normalized=.true.)
    call check(maxval(abs(integrated%power-direct%power))<1.0e-13_dp, &
        'single-taper exact band agreement')

    ! Fixed taper weights must be applied to whole spectral matrices,
    ! including the imaginary parts of their off-diagonal entries.
    out = varp_expected_multitaper_spectrum(a,sigma,taper,frequency, &
        taper_weights=weights,normalized=.true.)
    call check(out%status == varp_spectrum_success,'weighted bank accepted')
    allocate(expected,mold=out%power)
    expected = (0.0_dp,0.0_dp)
    do k = 1, 3
        one = varp_expected_taper_spectrum(a,sigma,taper(:,k),frequency,normalized=.true.)
        expected = expected + weights(k)/sum(weights)*one%power
    end do
    call check(maxval(abs(out%power-expected))<1.0e-12_dp, &
        'energy-normalized weighted spectral matrix')
    call check(abs(aimag(out%power(3,1,4)))>1.0e-6_dp, &
        'phase-sensitive off-diagonal spectra')
    call check(maxval(abs(out%lag0_covariance-one%lag0_covariance))<1.0e-13_dp, &
        'stationary covariance unaffected by taper bank')

    integrated = varp_expected_band_multitaper_spectrum(a,sigma,taper,bands, &
        taper_weights=weights,normalized=.true.)
    call check(integrated%status == varp_spectrum_success,'exact weighted bands accepted')
    call check(maxval(abs(integrated%power(1,:,:)+integrated%power(2,:,:) &
        -integrated%power(3,:,:)))<2.0e-12_dp,'exact adjacent-band additivity')
    call check(maxval(abs(integrated%power(3,:,:)-out%lag0_covariance)) &
        <2.0e-12_dp,'normalized full-period energy identity')
    expected_variance = sum(weights(:)*sum(taper*taper,dim=1))/sum(weights)
    direct = varp_expected_band_multitaper_spectrum(a,sigma,taper,bands, &
        taper_weights=weights)
    call check(maxval(abs(direct%power(3,:,:)- &
        expected_variance*out%lag0_covariance))<2.0e-11_dp, &
        'unnormalized full-period energy identity')
    do k = 1, 3
        one = varp_expected_band_taper_spectrum(a,sigma,taper(:,k),bands, &
            average=.true.,normalized=.true.)
        if (k == 1) then
            deallocate(expected)
            allocate(expected,mold=one%power)
            expected = (0.0_dp,0.0_dp)
        end if
        expected = expected + weights(k)/sum(weights)*one%power
    end do
    direct = varp_expected_band_multitaper_spectrum(a,sigma,taper,bands, &
        taper_weights=weights,average=.true.,normalized=.true.)
    call check(maxval(abs(direct%power-expected))<1.0e-12_dp, &
        'complex band-average spectral matrices')
    do q = 1, size(frequency)
        do j = 1, 5
            do i = 1, 5
                call check(abs(out%power(q,i,j)-conjg(out%power(q,j,i))) &
                    <2.0e-11_dp,'weighted spectral matrix Hermitian')
            end do
        end do
    end do

    ! Canonical coherence of a mixture is not the mean of coherences.
    matrix%frequencies = out%frequencies
    matrix%power = out%power
    matrix%status = 0
    reference = canonical_coherence_matrix(matrix,[1,2],[3,4],conditioned_on=[5], &
        ridge=0.01_dp,shrinkage=0.06_dp,return_vectors=.true.)
    mode = varp_expected_multitaper_canonical_coherence(a,sigma,taper,frequency, &
        [1,2],[3,4],conditioned_on=[5],taper_weights=weights,ridge=0.01_dp, &
        shrinkage=0.06_dp,normalized=.true.,return_vectors=.true.)
    call check(mode%status == 0 .and. all(mode%valid),'conditional two-mode result')
    call check(size(mode%squared_modes,2) == 2,'two canonical strengths')
    call check(allocated(mode%coefficients_a),'canonical directions returned')
    call check(maxval(abs(mode%squared_modes-reference%squared_modes))<1.0e-12_dp, &
        'canonical coherence of weighted matrix')
    call check(maxval(abs(mode%coherence-reference%coherence))<1.0e-12_dp, &
        'aggregate canonical coherence of weighted matrix')
    matrix%frequencies = direct%frequencies
    matrix%power = direct%power
    reference = canonical_coherence_matrix(matrix,[1,2],[3,4],conditioned_on=[5])
    mode = varp_expected_band_multitaper_canonical_coherence(a,sigma,taper,bands, &
        [1,2],[3,4],conditioned_on=[5],taper_weights=weights, &
        average=.true.,normalized=.true.)
    call check(mode%status == 0 .and. all(mode%valid),'conditional band coherence')
    call check(maxval(abs(mode%squared_modes-reference%squared_modes))<1.0e-12_dp, &
        'canonical band coherence matches direct fit')
    mode = varp_expected_multitaper_canonical_coherence(a,sigma,taper,frequency, &
        [1,2],[2,3])
    call check(mode%status /= 0,'overlapping groups rejected')

    ! Real DPSS tapers generated by SciFort exercise a standard independent
    ! finite-window bank, including an optional concentration-ratio weight.
    bank = dpss(13,2.5_dp,kmax=3)
    call check(bank%status == dpss_status_success,'DPSS bank constructed')
    direct = varp_expected_multitaper_spectrum(a,sigma,bank%tapers,frequency, &
        taper_weights=bank%ratios,normalized=.true.)
    call check(direct%status == varp_spectrum_success,'DPSS weighted target')
    call check(all(ieee_is_finite(real(direct%power,dp))) .and. &
        all(ieee_is_finite(aimag(direct%power))),'finite DPSS target')

    ! Zero-weight windows do not enter the estimator, including zero-energy
    ! windows. Bad weights and nonstationary models must be rejected.
    taper(:,3) = 0.0_dp
    bad = varp_expected_multitaper_spectrum(a,sigma,taper,frequency, &
        taper_weights=[1.0_dp,0.0_dp,0.0_dp])
    one = varp_expected_taper_spectrum(a,sigma,taper(:,1),frequency)
    call check(bad%status == varp_spectrum_success .and. &
        maxval(abs(bad%power-one%power))<1.0e-12_dp, &
        'zero weights omit unused zero-energy windows')
    bad = varp_expected_multitaper_spectrum(a,sigma,taper,frequency)
    call check(bad%status == varp_spectrum_invalid_input, &
        'selected zero-energy window rejected')
    bad = varp_expected_multitaper_spectrum(a,sigma,taper,frequency, &
        taper_weights=[0.0_dp,0.0_dp,0.0_dp])
    call check(bad%status == varp_spectrum_invalid_input,'all-zero weights rejected')
    bad = varp_expected_multitaper_spectrum(a,sigma,taper,frequency, &
        taper_weights=[-1.0_dp,0.0_dp,0.0_dp])
    call check(bad%status == varp_spectrum_invalid_input,'negative weight rejected')
    bad = varp_expected_multitaper_spectrum(a,sigma,taper,frequency, &
        taper_weights=[1.0_dp])
    call check(bad%status == varp_spectrum_invalid_input,'wrong weight length rejected')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    bad = varp_expected_multitaper_spectrum(a,sigma,taper,frequency, &
        taper_weights=[nan_value,1.0_dp,0.0_dp])
    call check(bad%status == varp_spectrum_invalid_input,'nonfinite weights rejected')
    bad = varp_expected_band_multitaper_spectrum(a,sigma,taper,bands(:,1:0))
    call check(bad%status == varp_spectrum_invalid_input,'zero bands rejected')
    bad = varp_expected_multitaper_spectrum(a,sigma,taper(:,1:0),frequency)
    call check(bad%status == varp_spectrum_invalid_input,'zero tapers rejected')
    a(1,1,1) = 1.6_dp
    bad = varp_expected_multitaper_spectrum(a,sigma,taper(:,1:2),frequency)
    call check(bad%status == varp_spectrum_nonstationary,'unstable VAR(p) rejected')

    print '(a)', 'PASS exact fixed-weight VAR(p) multitaper spectra and canonical bands'
contains
    subroutine check(condition,message)
        logical, intent(in) :: condition !! true means the test passed
        character(*), intent(in) :: message !! failure description
        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine check

end program test_signal_varp_multitaper_target
