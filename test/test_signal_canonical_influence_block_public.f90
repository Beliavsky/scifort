! SPDX-License-Identifier: MIT
program test_signal_canonical_influence_block_public
    use scifort_signal, only : dp, canonical_block_length_result, &
        select_canonical_influence_block_length, canonical_moving_block_bootstrap_result, &
        canonical_moving_block_bootstrap, multitaper_canonical_moving_block_bootstrap, &
        multitaper_spectral_matrix_spectrogram, multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_block_length_result) :: pilot
    type(canonical_moving_block_bootstrap_result) :: automatic, fixed, complex_direct
    real(dp) :: x(384,3)
    complex(dp) :: z(384,3)
    integer :: i, j

    do i = 1, 384
        do j = 1, 3
            x(i,j) = sin(0.071_dp*i*j)+0.24_dp*cos(0.043_dp*i+0.16_dp*j)
            z(i,j) = cmplx(x(i,j),0.18_dp*sin(0.051_dp*i*j),dp)
        end do
    end do
    spectra = multitaper_spectral_matrix_spectrogram(x,3.0_dp,24, &
        noverlap=0,kmax=4,weighting='eigen')
    pilot = select_canonical_influence_block_length(spectra,[1],[2], &
        mode_index=1,frequency_index=2,ridge=0.05_dp)
    if (pilot%status /= 0 .or. pilot%block_length < 1 .or. &
        pilot%mode_index /= 1 .or. .not. allocated(pilot%influence)) &
        error stop 'public influence pilot'
    automatic = canonical_moving_block_bootstrap(spectra,[1],[2],n_resamples=25, &
        seed=1122,ridge=0.05_dp,selection_method='influence', &
        selection_mode=1,selection_frequency=2)
    fixed = canonical_moving_block_bootstrap(spectra,[1],[2],pilot%block_length, &
        n_resamples=25,seed=1122,ridge=0.05_dp)
    if (automatic%status /= 0 .or. fixed%status /= 0) &
        error stop 'public influence bootstrap'
    if (automatic%selection_method /= 'influence' .or. &
        automatic%selection_frequency /= 2 .or. &
        any(automatic%coherence_lower95 /= fixed%coherence_lower95)) &
        error stop 'public influence metadata/reproducibility'
    complex_direct = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,24,[1],[2], &
        n_resamples=22,seed=1234,kmax=4,weighting='adaptive',ridge=0.05_dp, &
        bands=reshape([0.05_dp,0.20_dp],[2,1]),selection_method='influence')
    if (complex_direct%status /= 0 .or. .not. complex_direct%automatic_block_length .or. &
        .not. allocated(complex_direct%selection_influence)) &
        error stop 'public complex adaptive influence bootstrap'
    print '(a)', 'PASS public influence-targeted canonical bootstrap interface'
end program test_signal_canonical_influence_block_public
