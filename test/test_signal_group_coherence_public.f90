! SPDX-License-Identifier: MIT
! Public scifort_signal re-export and generic-resolution smoke test.
program test_signal_group_coherence_public
    use scifort_signal, only : dp, group_coherence_matrix, group_coherence_tf_matrix, &
        multitaper_group_coherence, multitaper_tf_group_coherence, &
        band_group_coherence_matrix, multitaper_band_group_coherence, &
        group_coherence_result, group_coherence_tf_result, multitaper_matrix_result
    implicit none
    type(group_coherence_result) :: result, direct, band
    type(group_coherence_tf_result) :: times
    type(multitaper_matrix_result) :: spectra
    real(dp) :: x(96,3), bounds(2,1)
    complex(dp) :: z(96,3)
    integer :: i

    do i = 1, 96
        x(i,1) = cos(0.09_dp*i)+0.1_dp*sin(0.21_dp*i)
        x(i,2) = sin(0.11_dp*i)+0.15_dp*cos(0.14_dp*i)
        x(i,3) = x(i,1)+0.2_dp*sin(0.27_dp*i)
    end do
    z = cmplx(x,0.25_dp*x,dp)
    result = multitaper_group_coherence(x,3.0_dp,[1],[2],conditioned_on=[3], &
        kmax=6,ridge=0.02_dp,shrinkage=0.2_dp)
    if (result%status /= 0 .or. .not. all(result%valid)) error stop 'public real group API'
    direct = multitaper_group_coherence(z,3.0_dp,[1],[2],conditioned_on=[3], &
        kmax=6,ridge=0.02_dp,shrinkage=0.2_dp)
    if (direct%status /= 0 .or. .not. all(direct%valid)) error stop 'public complex group API'
    times = multitaper_tf_group_coherence(z,3.0_dp,64,[1],[2],conditioned_on=[3], &
        kmax=6,ridge=0.02_dp,shrinkage=0.2_dp)
    if (times%status /= 0 .or. .not. all(times%valid)) error stop 'public time-frequency group API'
    allocate(spectra%frequencies(2),spectra%power(2,3,3))
    spectra%frequencies = [0.0_dp,0.5_dp]
    spectra%power = (0.0_dp,0.0_dp)
    do i = 1, 3
        spectra%power(:,i,i) = cmplx(1.0_dp,0.0_dp,dp)
    end do
    spectra%status = 0
    direct = group_coherence_matrix(spectra,[1],[2],conditioned_on=[3])
    if (direct%status /= 0) error stop 'public direct group API'
    if (maxval(abs(direct%coherence)) > 1.0e-12_dp) error stop 'independent channels'
    bounds(:,1) = [0.1_dp,0.4_dp]
    band = band_group_coherence_matrix(spectra,bounds,[1],[2],conditioned_on=[3])
    if (band%status /= 0) error stop 'public integrated group API'
    band = multitaper_band_group_coherence(x,3.0_dp,bounds,[1],[2], &
        conditioned_on=[3],kmax=6,ridge=0.02_dp,shrinkage=0.2_dp)
    if (band%status /= 0 .or. .not. all(band%valid)) error stop 'public multitaper band API'
    print *, 'PASS group coherence public API'
end program test_signal_group_coherence_public
