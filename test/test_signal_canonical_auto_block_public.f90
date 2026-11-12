! SPDX-License-Identifier: MIT
program test_signal_canonical_auto_block_public
    use scifort_signal, only : dp, canonical_block_length_result, select_canonical_block_length, &
        canonical_moving_block_bootstrap_result, canonical_moving_block_bootstrap, &
        multitaper_canonical_moving_block_bootstrap, multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_block_length_result) :: choice
    type(canonical_moving_block_bootstrap_result) :: a, b, c
    real(dp) :: x(320,3)
    complex(dp) :: z(320,3)
    integer :: i, j

    do i = 1, size(x,1)
        do j = 1, 3
            x(i,j) = sin(0.071_dp*i*j)+0.24_dp*cos(0.043_dp*i+0.16_dp*j)
            z(i,j) = cmplx(x(i,j),0.18_dp*sin(0.051_dp*i*j),dp)
        end do
    end do
    spectra = multitaper_spectral_matrix_spectrogram(x,3.0_dp,16, &
        noverlap=0,kmax=4,weighting='eigen')
    choice = select_canonical_block_length(spectra)
    if (choice%status /= 0 .or. choice%block_length < 1 .or. &
        choice%maximum_block_length < choice%block_length) error stop 'public selector'
    a = canonical_moving_block_bootstrap(spectra,[1],[2],n_resamples=25,seed=1122, &
        ridge=0.05_dp)
    b = canonical_moving_block_bootstrap(spectra,[1],[2],choice%block_length, &
        n_resamples=25,seed=1122,ridge=0.05_dp)
    if (a%status /= 0 .or. b%status /= 0 .or. &
        .not. a%automatic_block_length .or. b%automatic_block_length) &
        error stop 'public auto metadata'
    if (any(a%coherence_lower95 /= b%coherence_lower95) .or. &
        any(a%coherence_stderr /= b%coherence_stderr)) &
        error stop 'public auto/fixed reproducibility'
    c = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,16,[1],[2], &
        n_resamples=24,seed=1122,kmax=4,weighting='adaptive',ridge=0.05_dp, &
        bands=reshape([0.05_dp,0.20_dp],[2,1]))
    if (c%status /= 0 .or. .not. c%automatic_block_length .or. &
        size(c%frequencies) /= 1) error stop 'public complex auto-band API'
    print '(a)', 'PASS canonical auto-block selector umbrella API'
end program test_signal_canonical_auto_block_public
