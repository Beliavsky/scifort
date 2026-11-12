! SPDX-License-Identifier: MIT
! Verify the fixed-taper exact-band API is visible through scifort_signal.
program test_signal_fixed_taper_spectrogram_public
    use scifort_signal, only : dp, dpss, dpss_result, &
        fixed_taper_spectral_matrix_spectrogram, &
        canonical_moving_block_bootstrap, canonical_moving_block_bootstrap_result, &
        multitaper_spectral_matrix_spectrogram_result, spectral_status_success
    implicit none
    type(dpss_result) :: bank
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_moving_block_bootstrap_result) :: result
    real(dp) :: x(240,5), bands(2,2)
    integer :: i, j

    bank = dpss(20,2.5_dp,kmax=3,sym=.false.)
    if (bank%status /= 0) error stop 'DPSS bank'
    do j = 1, 5
        do i = 1, 240
            x(i,j) = sin(0.17_dp*real(i*j,dp))+0.34_dp*cos(0.1_dp*real(i+j,dp))
        end do
    end do
    bands(:,1) = [0.0_dp,0.125_dp]
    bands(:,2) = [0.125_dp,0.5_dp]
    spectra = fixed_taper_spectral_matrix_spectrogram( &
        x,bank%tapers,[0.1375_dp],bands=bands,taper_weights=bank%ratios)
    if (spectra%status /= spectral_status_success) error stop 'fixed taper spectrum'
    if (any(shape(spectra%power) /= [3,12,5,5])) error stop 'wrong output shape'
    if (spectra%effective_tapers < 1.0_dp .or. &
        spectra%effective_tapers > 3.0_dp) error stop 'effective taper count'
    result = canonical_moving_block_bootstrap(spectra,[1,2],[3,4], &
        conditioned_on=[5],n_resamples=25,block_length=2,seed=321,circular=.true.)
    if (.not. allocated(result%squared_modes)) error stop 'bootstrap output'
    if (any(shape(result%squared_modes) /= [3,2])) error stop 'canonical root dimensions'
    print '(a)', 'PASS public fixed-taper exact-band CSD and canonical bootstrap API'
end program test_signal_fixed_taper_spectrogram_public
