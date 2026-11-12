! SPDX-License-Identifier: MIT
! Top-level scifort_signal export and generic resolution for optional vectors.
program test_signal_canonical_vectors_public
    use scifort_signal, only : dp, multitaper_matrix_result, canonical_coherence_result, &
        canonical_coherence_tf_result, canonical_coherence_matrix, &
        multitaper_canonical_coherence, multitaper_tf_canonical_coherence, &
        multitaper_band_canonical_coherence, multitaper_tf_band_canonical_coherence
    implicit none
    type(multitaper_matrix_result) :: spectrum
    type(canonical_coherence_result) :: result, band
    type(canonical_coherence_tf_result) :: moving
    real(dp) :: x(128,4), bands(2,1)
    complex(dp) :: z(128,4)
    integer :: i, j

    do i = 1, 128
        do j = 1, 4
            x(i,j) = sin(0.09_dp*i*j)+0.1_dp*cos(0.11_dp*i+0.2_dp*j)
            z(i,j) = cmplx(x(i,j),0.3_dp*cos(0.13_dp*i*j),dp)
        end do
    end do
    result = multitaper_canonical_coherence(x,3.0_dp,[1,2],[3,4],kmax=5, &
        ridge=0.05_dp,return_vectors=.true.)
    if (result%status /= 0 .or. .not. all(result%valid)) &
        error stop 'public real coefficient vectors'
    if (.not. allocated(result%coefficients_a)) error stop 'public A allocation'
    if (size(result%coefficients_a,2) /= 2) error stop 'public real A shape'
    result = multitaper_canonical_coherence(z,3.0_dp,[1],[3,4],conditioned_on=[2], &
        kmax=5,weighting='adaptive',cv_shrinkage=.true.,ridge=0.05_dp, &
        return_vectors=.true.)
    if (result%status /= 0 .or. .not. all(result%valid)) &
        error stop 'public adaptive complex vectors'
    if (size(result%coefficients_b,2) /= 2) error stop 'public complex B shape'
    moving = multitaper_tf_canonical_coherence(z,3.0_dp,64,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    if (moving%status /= 0 .or. .not. all(moving%valid)) &
        error stop 'public moving complex vectors'
    if (.not. allocated(moving%coefficients_b)) error stop 'public moving allocation'
    bands(:,1) = [0.10_dp,0.35_dp]
    band = multitaper_band_canonical_coherence(x,3.0_dp,bands,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    if (band%status /= 0 .or. .not. all(band%valid)) error stop 'public band vectors'
    moving = multitaper_tf_band_canonical_coherence(z,3.0_dp,64,bands,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    if (moving%status /= 0 .or. .not. all(moving%valid)) &
        error stop 'public moving band vectors'
    allocate(spectrum%frequencies(1),spectrum%power(1,4,4))
    spectrum%status = 0
    spectrum%frequencies = [0.0_dp]
    spectrum%power = (0.0_dp,0.0_dp)
    do i = 1, 4
        spectrum%power(1,i,i) = (1.0_dp,0.0_dp)
    end do
    result = canonical_coherence_matrix(spectrum,[1,2],[3,4])
    if (allocated(result%coefficients_a)) error stop 'public default allocation changed'
    result = canonical_coherence_matrix(spectrum,[1,2],[3,4],return_vectors=.true.)
    if (result%status /= 0 .or. any(result%directions_stable)) &
        error stop 'public zero-root directions'
    print *, 'PASS canonical vectors public API'
end program test_signal_canonical_vectors_public
