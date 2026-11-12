! SPDX-License-Identifier: MIT
! Verify the block extension is exported from the umbrella scifort_signal module.
program test_signal_canonical_block_jackknife_public
    use scifort_signal, only : dp, canonical_block_jackknife, &
        multitaper_canonical_block_jackknife, canonical_window_jackknife_result, &
        canonical_window_jackknife, multitaper_spectral_matrix_spectrogram, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_spectral_matrix_spectrogram_result) :: spectra
    type(canonical_window_jackknife_result) :: a, b, c, d
    real(dp) :: x(256,5), bands(2,2)
    complex(dp) :: z(256,5)
    integer :: i, j

    do i = 1, 256
        do j = 1, 5
            x(i,j) = sin(0.097_dp*i*j)+0.17_dp*cos(0.23_dp*i+0.08_dp*j) &
                + 0.07_dp*sin(0.39_dp*i*i+0.20_dp*j)
            z(i,j) = cmplx(x(i,j),0.19_dp*cos(0.078_dp*i*j+0.15_dp*j),dp)
        end do
    end do
    bands(:,1) = [0.04_dp,0.18_dp]
    bands(:,2) = [0.23_dp,0.43_dp]
    a = multitaper_canonical_block_jackknife(x,3.0_dp,16,[1,2],[3,4],2, &
        conditioned_on=[5],kmax=5,weighting='eigen',ridge=0.08_dp,shrinkage=0.15_dp)
    if (a%status /= 0 .or. a%n_windows /= 16 .or. a%n_blocks /= 8 .or. &
        a%block_length /= 2) error stop 'public real block API'
    b = multitaper_canonical_block_jackknife(z,3.0_dp,16,[1,2],[3,4],4, &
        conditioned_on=[5],kmax=5,weighting='adaptive',ridge=0.09_dp, &
        shrinkage=0.12_dp,bands=bands)
    if (b%status /= 0 .or. b%n_blocks /= 4 .or. size(b%frequencies) /= 2) &
        error stop 'public complex band block API'
    spectra = multitaper_spectral_matrix_spectrogram(z,3.0_dp,16, &
        noverlap=0,kmax=5,weighting='adaptive')
    c = canonical_block_jackknife(spectra,[1,2],[3,4],4,conditioned_on=[5], &
        ridge=0.09_dp,shrinkage=0.12_dp,bands=bands)
    if (c%status /= 0) error stop 'public block matrix API'
    if (maxval(abs(c%squared_modes-b%squared_modes)) > 1.0e-12_dp) &
        error stop 'matrix and wrapped block estimates disagree'
    d = canonical_window_jackknife(spectra,[1,2],[3,4],conditioned_on=[5], &
        ridge=0.09_dp,shrinkage=0.12_dp,block_length=4,bands=bands)
    if (d%status /= 0) error stop 'public window API with explicit block length'
    if (maxval(abs(c%stderr-d%stderr)) > 1.0e-12_dp) &
        error stop 'block alias and explicit block keyword differ'
    print '(a)', 'PASS canonical block jackknife umbrella public API'
end program test_signal_canonical_block_jackknife_public
