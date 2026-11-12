! SPDX-License-Identifier: MIT
! Top-level re-export, generic resolution and numeric consistency smoke test.
program test_signal_canonical_coherence_public
    use scifort_signal, only : dp, multitaper_matrix_result, &
        canonical_coherence_result, canonical_coherence_tf_result, &
        canonical_coherence_matrix, canonical_coherence_tf_matrix, &
        band_canonical_coherence_matrix, band_canonical_coherence_tf_matrix, &
        multitaper_canonical_coherence, multitaper_tf_canonical_coherence, &
        multitaper_band_canonical_coherence, multitaper_tf_band_canonical_coherence
    use scifort_signal_multitaper, only : multitaper_spectral_matrix_spectrogram_result
    implicit none

    type(multitaper_matrix_result) :: spectra
    type(multitaper_spectral_matrix_spectrogram_result) :: spectrogram
    type(canonical_coherence_result) :: out, band
    type(canonical_coherence_tf_result) :: moving, band_moving
    real(dp) :: x(96,4), bands(2,1)
    complex(dp) :: z(96,4)
    integer :: i, j

    do i = 1, 96
        do j = 1, 4
            x(i,j) = sin(0.07_dp*i*j)+0.3_dp*cos(0.13_dp*i+0.2_dp*j)
            z(i,j) = cmplx(x(i,j),0.18_dp*cos(0.19_dp*i*j),dp)
        end do
    end do
    out = multitaper_canonical_coherence(x,3.0_dp,[1,2],[3,4], &
        kmax=6,ridge=0.05_dp,shrinkage=0.2_dp)
    if (out%status /= 0 .or. .not. all(out%valid)) error stop 'public real canonical'
    if (size(out%squared_modes,2) /= 2) error stop 'public canonical mode count'
    out = multitaper_canonical_coherence(z,3.0_dp,[1],[3,4],conditioned_on=[2], &
        kmax=6,weighting='adaptive',cv_shrinkage=.true.,ridge=0.05_dp)
    if (out%status /= 0 .or. .not. all(out%valid)) error stop 'public complex canonical'
    moving = multitaper_tf_canonical_coherence(z,3.0_dp,64,[1,2],[3,4], &
        kmax=6,ridge=0.05_dp,shrinkage=0.2_dp)
    if (moving%status /= 0 .or. .not. all(moving%valid)) error stop 'public TF canonical'
    allocate(spectra%frequencies(3),spectra%power(3,4,4))
    spectra%status = 0
    spectra%frequencies = [0.0_dp,0.2_dp,0.5_dp]
    spectra%power = (0.0_dp,0.0_dp)
    do i = 1, 4
        spectra%power(:,i,i) = cmplx(1.0_dp,0.0_dp,dp)
    end do
    out = canonical_coherence_matrix(spectra,[1,2],[3,4])
    if (out%status /= 0 .or. maxval(abs(out%squared_modes)) > 1.0e-13_dp) &
        error stop 'public independent canonical'
    bands(:,1) = [0.05_dp,0.35_dp]
    band = band_canonical_coherence_matrix(spectra,bands,[1,2],[3,4])
    if (band%status /= 0) error stop 'public band canonical'
    band = multitaper_band_canonical_coherence(x,3.0_dp,bands,[1,2],[3,4], &
        kmax=6,ridge=0.05_dp,shrinkage=0.2_dp)
    if (band%status /= 0 .or. .not. all(band%valid)) error stop 'public real band canonical'
    band = multitaper_band_canonical_coherence(z,3.0_dp,bands,[1,2],[3,4], &
        kmax=6,ridge=0.05_dp,shrinkage=0.2_dp)
    if (band%status /= 0 .or. .not. all(band%valid)) error stop 'public complex band canonical'
    allocate(spectrogram%frequencies(3),spectrogram%times(2),spectrogram%power(3,2,4,4))
    spectrogram%status = 0
    spectrogram%frequencies = spectra%frequencies
    spectrogram%times = [1.0_dp,2.0_dp]
    do i = 1, 2
        spectrogram%power(:,i,:,:) = spectra%power
    end do
    moving = canonical_coherence_tf_matrix(spectrogram,[1,2],[3,4])
    if (moving%status /= 0) error stop 'public TF matrix canonical'
    band_moving = band_canonical_coherence_tf_matrix(spectrogram,bands,[1,2],[3,4])
    if (band_moving%status /= 0) error stop 'public TF band matrix canonical'
    band_moving = multitaper_tf_band_canonical_coherence(x,3.0_dp,64,bands,[1,2],[3,4], &
        kmax=6,ridge=0.05_dp,shrinkage=0.2_dp)
    if (band_moving%status /= 0 .or. .not. all(band_moving%valid)) &
        error stop 'public TF band real canonical'
    band_moving = multitaper_tf_band_canonical_coherence(z,3.0_dp,64,bands,[1,2],[3,4], &
        kmax=6,ridge=0.05_dp,shrinkage=0.2_dp)
    if (band_moving%status /= 0 .or. .not. all(band_moving%valid)) &
        error stop 'public TF band complex canonical'
    print *, 'PASS canonical coherence public API'
end program test_signal_canonical_coherence_public
