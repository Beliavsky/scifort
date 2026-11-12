! SPDX-License-Identifier: MIT
program test_signal_band_coherence_public
    use scifort_signal, only : dp, multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result, multiple_coherence_result, &
        multiple_coherence_tf_result, partial_coherence_result, partial_coherence_tf_result, &
        integrate_spectral_bands, integrate_tf_spectral_bands, &
        band_multiple_coherence_matrix, band_multiple_coherence_tf_matrix, &
        band_partial_coherence_matrix, band_partial_coherence_tf_matrix, &
        multitaper_band_multiple_coherence, multitaper_tf_band_multiple_coherence, &
        multitaper_band_partial_coherence, multitaper_tf_band_partial_coherence, &
        partial_coherence_success
    implicit none
    type(multitaper_matrix_result) :: spectral, integrated
    type(multitaper_spectral_matrix_spectrogram_result) :: spectral_tf, integrated_tf
    type(multiple_coherence_result) :: result, direct
    type(multiple_coherence_tf_result) :: result_tf
    type(partial_coherence_result) :: partial
    type(partial_coherence_tf_result) :: partial_tf
    real(dp) :: bands(2,1), x(80,2)
    complex(dp) :: z(80,2)
    integer :: i

    bands(:,1) = [0.05_dp,0.35_dp]
    allocate(spectral%frequencies(3),spectral%power(3,2,2))
    spectral%status = 0
    spectral%frequencies = [0.0_dp,0.20_dp,0.5_dp]
    spectral%power = (0.0_dp,0.0_dp)
    do i = 1, 3
        spectral%power(i,1,1) = (2.0_dp,0.0_dp)
        spectral%power(i,2,2) = (1.0_dp,0.0_dp)
        spectral%power(i,1,2) = cmplx(0.5_dp,0.2_dp,dp)
        spectral%power(i,2,1) = conjg(spectral%power(i,1,2))
    end do
    integrated = integrate_spectral_bands(spectral,bands)
    if (integrated%status /= 0) error stop 'public band integration'
    if (abs(real(integrated%power(1,1,1),dp)-0.6_dp) > 2.0e-14_dp) &
        error stop 'integrated power'
    result = band_multiple_coherence_matrix(spectral,bands)
    if (result%status /= partial_coherence_success) error stop 'public band multiple coherence'
    if (abs(result%coherence(1,1)-0.145_dp) > 2.0e-13_dp) error stop 'multiple coherence value'
    partial = band_partial_coherence_matrix(spectral,bands,return_coherency=.true.)
    if (partial%status /= partial_coherence_success) error stop 'public band partial coherence'

    allocate(spectral_tf%frequencies(3),spectral_tf%times(2),spectral_tf%power(3,2,2,2))
    spectral_tf%status = 0
    spectral_tf%frequencies = spectral%frequencies
    spectral_tf%times = [1.0_dp,2.0_dp]
    spectral_tf%power(:,1,:,:) = spectral%power
    spectral_tf%power(:,2,:,:) = 2.0_dp*spectral%power
    integrated_tf = integrate_tf_spectral_bands(spectral_tf,bands)
    if (integrated_tf%status /= 0) error stop 'public TF integration'
    result_tf = band_multiple_coherence_tf_matrix(spectral_tf,bands)
    partial_tf = band_partial_coherence_tf_matrix(spectral_tf,bands)
    if (result_tf%status /= partial_coherence_success .or. &
        partial_tf%status /= partial_coherence_success) error stop 'public TF band coherence'

    do i = 1, 80
        x(i,1) = sin(0.13_dp*real(i,dp))
        x(i,2) = x(i,1)+0.2_dp*cos(0.21_dp*real(i,dp))
        z(i,1) = cmplx(x(i,1),0.15_dp*x(i,2),dp)
        z(i,2) = cmplx(x(i,2),0.10_dp*x(i,1),dp)
    end do
    direct = multitaper_band_multiple_coherence(x,3.0_dp,bands,kmax=5,ridge=0.03_dp)
    if (direct%status /= partial_coherence_success) error stop 'real direct band API'
    direct = multitaper_band_multiple_coherence(z,3.0_dp,bands,kmax=5,ridge=0.03_dp)
    if (direct%status /= partial_coherence_success) error stop 'complex direct band API'
    result_tf = multitaper_tf_band_multiple_coherence(x,3.0_dp,40,bands, &
        kmax=5,ridge=0.03_dp)
    if (result_tf%status /= partial_coherence_success) error stop 'real TF direct band API'
    result_tf = multitaper_tf_band_multiple_coherence(z,3.0_dp,40,bands, &
        kmax=5,ridge=0.03_dp)
    if (result_tf%status /= partial_coherence_success) error stop 'complex TF direct band API'
    partial = multitaper_band_partial_coherence(x,3.0_dp,bands,kmax=5,ridge=0.03_dp)
    if (partial%status /= partial_coherence_success) error stop 'real direct partial bands'
    partial = multitaper_band_partial_coherence(z,3.0_dp,bands,kmax=5,ridge=0.03_dp)
    if (partial%status /= partial_coherence_success) error stop 'complex direct partial bands'
    partial_tf = multitaper_tf_band_partial_coherence(x,3.0_dp,40,bands, &
        kmax=5,ridge=0.03_dp)
    if (partial_tf%status /= partial_coherence_success) error stop 'real direct TF partial bands'
    partial_tf = multitaper_tf_band_partial_coherence(z,3.0_dp,40,bands, &
        kmax=5,ridge=0.03_dp)
    if (partial_tf%status /= partial_coherence_success) error stop 'complex direct TF partial bands'
    print *, 'PASS public spectral band integration / partial / multiple coherence'
end program test_signal_band_coherence_public
