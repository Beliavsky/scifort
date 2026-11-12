! SPDX-License-Identifier: MIT
program test_signal_varp_multitaper_public
    use scifort_signal, only : dp, dpss, dpss_result, varp_taper_spectrum_result, &
        canonical_coherence_result, varp_expected_multitaper_spectrum, &
        varp_expected_band_multitaper_spectrum, &
        varp_expected_multitaper_canonical_coherence, &
        varp_expected_band_multitaper_canonical_coherence, &
        varp_spectrum_success
    implicit none
    real(dp) :: a(5,5,2), sigma(5,5), bands(2,2)
    type(dpss_result) :: tapers
    type(varp_taper_spectrum_result) :: point, integrated
    type(canonical_coherence_result) :: canonical_point, canonical_band
    integer :: i

    a = 0.0_dp
    sigma = 0.0_dp
    do i = 1, 5
        a(i,i,1) = 0.35_dp
        a(i,i,2) = -0.14_dp
        sigma(i,i) = 1.0_dp
    end do
    sigma(1,3) = 0.5_dp
    sigma(3,1) = sigma(1,3)
    sigma(2,4) = -0.3_dp
    sigma(4,2) = sigma(2,4)
    tapers = dpss(14,2.5_dp,kmax=3)
    if (tapers%status /= 0) error stop 'DPSS bank construction failed'
    bands(:,1) = [-0.25_dp,0.0_dp]
    bands(:,2) = [0.0_dp,0.25_dp]
    point = varp_expected_multitaper_spectrum(a,sigma,tapers%tapers, &
        [-0.1_dp,0.0_dp,0.3_dp],taper_weights=tapers%ratios,normalized=.true.)
    integrated = varp_expected_band_multitaper_spectrum(a,sigma,tapers%tapers, &
        bands,taper_weights=tapers%ratios,normalized=.true.,average=.true.)
    if (point%status /= varp_spectrum_success .or. &
        integrated%status /= varp_spectrum_success) error stop 'multitaper spectral target'
    canonical_point = varp_expected_multitaper_canonical_coherence( &
        a,sigma,tapers%tapers,[-0.1_dp,0.0_dp,0.3_dp], &
        [1,2],[3,4],conditioned_on=[5],taper_weights=tapers%ratios, &
        normalized=.true.,return_vectors=.true.)
    canonical_band = varp_expected_band_multitaper_canonical_coherence( &
        a,sigma,tapers%tapers,bands,[1,2],[3,4],conditioned_on=[5], &
        taper_weights=tapers%ratios,normalized=.true.,average=.true.)
    if (canonical_point%status /= 0 .or. canonical_band%status /= 0) &
        error stop 'multitaper conditional coherence'
    if (.not. all(canonical_point%valid) .or. .not. all(canonical_band%valid)) &
        error stop 'multitaper canonical validity'
    if (.not. allocated(canonical_point%coefficients_a)) error stop 'canonical directions'
    if (size(canonical_point%squared_modes,2) /= 2) error stop 'canonical mode count'
    print '(a)', 'PASS public VAR(p) multitaper spectral and conditional canonical API'
end program test_signal_varp_multitaper_public
