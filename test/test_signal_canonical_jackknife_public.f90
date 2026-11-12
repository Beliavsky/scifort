! SPDX-License-Identifier: MIT
! Umbrella API regression: canonical uncertainty from independent windows.
program test_signal_canonical_jackknife_public
    use scifort_signal, only : dp, multitaper_canonical_window_jackknife, &
        canonical_window_jackknife, canonical_window_jackknife_result, &
        multitaper_spectral_matrix_spectrogram, multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(canonical_window_jackknife_result) :: result, direct, band
    type(multitaper_spectral_matrix_spectrogram_result) :: matrix
    real(dp) :: x(192,5), windows(2,2)
    complex(dp) :: z(192,5)
    integer :: i, j

    do i = 1, 192
        do j = 1, 5
            x(i,j) = sin(0.085_dp*i*j)+0.18_dp*cos(0.23_dp*i+0.07_dp*j) &
                + 0.04_dp*sin(0.47_dp*i*i+0.31_dp*j)
            z(i,j) = cmplx(x(i,j),0.17_dp*cos(0.071_dp*i*j+0.17_dp*j),dp)
        end do
    end do
    windows(:,1) = [0.04_dp,0.19_dp]
    windows(:,2) = [0.23_dp,0.37_dp]
    result = multitaper_canonical_window_jackknife(x,3.0_dp,24,[1,2],[3,4], &
        conditioned_on=[5],kmax=5,weighting='eigen',ridge=0.1_dp,shrinkage=0.12_dp)
    if (result%status /= 0 .or. result%n_windows /= 8) error stop 'public real multitaper'
    if (.not. all(result%coherence_valid)) error stop 'public real estimates'
    result = multitaper_canonical_window_jackknife(z,3.0_dp,24,[1,2],[3,4], &
        conditioned_on=[5],kmax=5,weighting='adaptive',ridge=0.1_dp,shrinkage=0.12_dp, &
        bands=windows)
    if (result%status /= 0 .or. size(result%frequencies) /= 2) error stop 'public complex bands'
    matrix = multitaper_spectral_matrix_spectrogram(z,3.0_dp,24, &
        noverlap=0,kmax=5,weighting='adaptive')
    direct = canonical_window_jackknife(matrix,[1,2],[3,4], &
        conditioned_on=[5],ridge=0.1_dp,shrinkage=0.12_dp,bands=windows)
    if (direct%status /= 0) error stop 'public matrix windows'
    if (maxval(abs(direct%squared_modes-result%squared_modes)) > 1.0e-12_dp) &
        error stop 'direct and wrapped estimator disagree'
    band = multitaper_canonical_window_jackknife(x,3.0_dp,24,[1],[2], &
        conditioned_on=[3,4,5],kmax=5,weighting='equal',ridge=0.15_dp,bands=windows)
    if (band%status /= 0 .or. size(band%squared_modes,2) /= 1) error stop 'public singleton'
    print '(a)', 'PASS canonical window jackknife public API'
end program test_signal_canonical_jackknife_public
