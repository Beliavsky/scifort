! SPDX-License-Identifier: MIT
program test_signal_varp_population_public
    use scifort_signal, only : dp, varp_expected_taper_spectrum, &
        varp_expected_band_taper_spectrum, varp_expected_canonical_coherence, &
        varp_expected_band_canonical_coherence, varp_taper_spectrum_result, &
        varp_spectrum_success, canonical_coherence_result, canonical_coherence_matrix, &
        multitaper_matrix_result
    implicit none
    type(varp_taper_spectrum_result) :: spectrum, band_spectrum, invalid
    type(canonical_coherence_result) :: point, integrated, reference, ref_band
    type(multitaper_matrix_result) :: s
    real(dp) :: a(5,5,2), sigma(5,5), taper(12), freq(3), bands(2,2)
    integer :: k

    a = 0.0_dp
    do k = 1, 5
        a(k,k,1) = 0.3_dp
        a(k,k,2) = -0.13_dp
        sigma(k,k) = 1.0_dp
    end do
    sigma = 0.0_dp
    do k = 1, 5
        sigma(k,k) = 1.0_dp
    end do
    sigma(1,3) = 0.42_dp
    sigma(3,1) = sigma(1,3)
    sigma(2,4) = -0.3_dp
    sigma(4,2) = sigma(2,4)
    a(1,4,1) = 0.18_dp
    a(3,2,2) = 0.12_dp
    taper = [(0.5_dp-0.5_dp*cos(2.0_dp*acos(-1.0_dp)*real(k-1,dp)/12.0_dp),k=1,12)]
    freq = [0.0_dp,0.17_dp,0.39_dp]
    bands(:,1) = [0.0_dp,0.24_dp]
    bands(:,2) = [-0.2_dp,0.3_dp]

    spectrum = varp_expected_taper_spectrum(a,sigma,taper,freq,normalized=.true.)
    call check(spectrum%status == varp_spectrum_success,'point spectrum produced')
    point = varp_expected_canonical_coherence(a,sigma,taper,freq,[1,2],[3,4], &
        conditioned_on=[5],normalized=.true.,return_vectors=.true.)
    call check(point%status == 0,'conditional canonical frequency target')
    call check(size(point%squared_modes,2) == 2,'two canonical modes')
    call check(all(point%valid),'all point modes valid')
    call check(allocated(point%coefficients_a),'mode coefficients opt-in')
    s%frequencies = spectrum%frequencies
    s%power = spectrum%power
    s%status = 0
    reference = canonical_coherence_matrix(s,[1,2],[3,4],conditioned_on=[5])
    call check(maxval(abs(reference%squared_modes-point%squared_modes))<1.0e-12_dp, &
        'point wrapper matches existing canonical solver')
    call check(maxval(abs(reference%coherence-point%coherence))<1.0e-12_dp, &
        'point wrapper aggregate matches existing canonical solver')

    band_spectrum = varp_expected_band_taper_spectrum(a,sigma,taper,bands, &
        average=.true.,normalized=.true.)
    call check(band_spectrum%status == varp_spectrum_success,'band population target')
    integrated = varp_expected_band_canonical_coherence(a,sigma,taper,bands, &
        [1,2],[3,4],conditioned_on=[5],average=.true.,normalized=.true.)
    call check(integrated%status == 0,'conditional canonical band target')
    s%frequencies = band_spectrum%frequencies
    s%power = band_spectrum%power
    s%status = 0
    ref_band = canonical_coherence_matrix(s,[1,2],[3,4],conditioned_on=[5])
    call check(all(integrated%valid),'valid band canonical targets')
    call check(maxval(abs(ref_band%squared_modes-integrated%squared_modes))<1.0e-12_dp, &
        'exact-band wrapper matches direct conditioned matrix')
    call check(maxval(abs(ref_band%coherence-integrated%coherence))<1.0e-12_dp, &
        'integrated group coherence preserved')
    call check(maxval(abs(integrated%coherence-point%coherence(1)))>1.0e-5_dp, &
        'integrate matrices before computing coherence')

    point = varp_expected_canonical_coherence(a,sigma,taper,freq,[1],[2,3], &
        conditioned_on=[4,5],ridge=0.02_dp,shrinkage=0.05_dp)
    call check(point%status == 0 .and. all(point%valid), &
        'two conditioning channels with explicit ridge/shrinkage')
    call check(size(point%squared_modes,2) == 1,'rectangular group modes')
    point = varp_expected_canonical_coherence(a,sigma,taper,freq,[1,2],[3,4])
    call check(point%status == 0 .and. all(point%valid), &
        'unconditional population canonical coherence')
    point = varp_expected_canonical_coherence(a,sigma,taper,freq,[1,2],[2,3])
    call check(point%status /= 0,'overlapping groups rejected')

    a(1,1,1) = 1.5_dp
    invalid = varp_expected_taper_spectrum(a,sigma,taper,freq)
    call check(invalid%status /= varp_spectrum_success,'unstable process rejected')
    print '(a)', 'PASS public VAR(p) conditional canonical spectral and band targets'

contains
    subroutine check(condition,message)
        logical, intent(in) :: condition !! required assertion
        character(*), intent(in) :: message !! name of failed assertion
        if (.not. condition) then
            print *, 'FAIL ',message
            error stop 1
        end if
    end subroutine check
end program test_signal_varp_population_public
