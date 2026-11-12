! SPDX-License-Identifier: MIT
! Verify moving-block canonical inference through umbrella scifort_signal API.
program test_signal_canonical_moving_bootstrap_public
    use scifort_signal, only : dp, canonical_moving_block_bootstrap, &
        multitaper_canonical_moving_block_bootstrap, canonical_moving_block_bootstrap_result, &
        multitaper_spectral_matrix_spectrogram, multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_moving_block_bootstrap_result) :: a, b, c, d
    real(dp) :: x(256,5), bands(2,2)
    complex(dp) :: z(256,5)
    integer :: i, j

    do i = 1, 256
        do j = 1, 5
            x(i,j) = sin(0.097_dp*i*j)+0.17_dp*cos(0.23_dp*i+0.08_dp*j) &
                +0.07_dp*sin(0.39_dp*i*i+0.20_dp*j)
            z(i,j) = cmplx(x(i,j),0.19_dp*cos(0.078_dp*i*j+0.15_dp*j),dp)
        end do
    end do
    bands(:,1) = [0.04_dp,0.18_dp]
    bands(:,2) = [0.23_dp,0.43_dp]
    a = multitaper_canonical_moving_block_bootstrap(x,3.0_dp,16,[1,2],[3,4],3, &
        n_resamples=26,conditioned_on=[5],kmax=5,weighting='eigen', &
        ridge=0.08_dp,shrinkage=0.15_dp,seed=5431)
    if (a%status /= 0 .or. a%n_windows /= 16 .or. a%blocks_per_resample /= 6 .or. &
        a%block_length /= 3 .or. a%n_resamples /= 26) error stop 'public real MBB'
    b = multitaper_canonical_moving_block_bootstrap(z,3.0_dp,16,[1,2],[3,4],4, &
        n_resamples=25,conditioned_on=[5],kmax=5,weighting='adaptive', &
        ridge=0.09_dp,shrinkage=0.12_dp,seed=5432,bands=bands,circular=.true.)
    if (b%status /= 0 .or. .not. b%circular .or. size(b%frequencies) /= 2) &
        error stop 'public complex circular band MBB'
    spectra = multitaper_spectral_matrix_spectrogram(z,3.0_dp,16, &
        noverlap=0,kmax=5,weighting='adaptive')
    c = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],4,n_resamples=25, &
        conditioned_on=[5],ridge=0.09_dp,shrinkage=0.12_dp, &
        seed=5432,bands=bands,circular=.true.)
    if (c%status /= 0 .or. maxval(abs(c%squared_modes-b%squared_modes)) > 1.0e-12_dp) &
        error stop 'public spectra and direct complex MBB mismatch'
    if (maxval(abs(c%stderr-b%stderr)) > 1.0e-12_dp) &
        error stop 'public complex MBB deterministic inference mismatch'
    d = canonical_moving_block_bootstrap(spectra,[1,2],[3,4],3,conditioned_on=[5], &
        ridge=0.09_dp,shrinkage=0.12_dp,seed=5432,bands=bands)
    if (d%status /= 0 .or. d%n_resamples /= 399 .or. d%circular) &
        error stop 'public default bootstrap repetitions and ordinary geometry'
    print '(a)', 'PASS canonical moving-block bootstrap umbrella public API'
end program test_signal_canonical_moving_bootstrap_public
