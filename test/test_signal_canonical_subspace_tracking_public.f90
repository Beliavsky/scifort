! SPDX-License-Identifier: MIT
! Public umbrella module compiles and exposes canonical subspace tracking.
program test_signal_canonical_subspace_tracking_public
    use scifort_signal, only : dp, canonical_coherence_result, canonical_coherence_tf_result, &
        canonical_subspace_tracks_result, canonical_subspace_tracks_tf_result, &
        track_canonical_subspaces, canonical_subspace_success, canonical_subspace_invalid_input, &
        canonical_coherence_matrix, canonical_coherence_tf_matrix, &
        multitaper_matrix_result, multitaper_spectral_matrix_spectrogram_result
    implicit none
    type(multitaper_matrix_result) :: spec
    type(multitaper_spectral_matrix_spectrogram_result) :: tf_spec
    type(canonical_coherence_result) :: modes
    type(canonical_coherence_tf_result) :: tf_modes
    type(canonical_subspace_tracks_result) :: result
    type(canonical_subspace_tracks_tf_result) :: tf_result
    integer :: f, i
    real(dp) :: strength

    allocate(spec%frequencies(5),spec%power(5,4,4))
    spec%frequencies = [0.0_dp,0.1_dp,0.2_dp,0.3_dp,0.4_dp]
    spec%power = (0.0_dp,0.0_dp)
    spec%status = 0
    do f = 1, 5
        do i = 1, 4
            spec%power(f,i,i) = (1.0_dp,0.0_dp)
        end do
        strength = 0.7_dp+0.001_dp*f
        spec%power(f,1,3) = cmplx(strength,0.0_dp,dp)
        spec%power(f,3,1) = conjg(spec%power(f,1,3))
        spec%power(f,2,4) = cmplx(strength+1.0e-8_dp,0.0_dp,dp)
        spec%power(f,4,2) = conjg(spec%power(f,2,4))
    end do
    modes = canonical_coherence_matrix(spec,[1,2],[3,4],return_vectors=.true.)
    if (modes%status /= 0 .or. .not. all(modes%valid)) error stop 1
    result = track_canonical_subspaces(modes,gap_tolerance=0.001_dp)
    if (result%status /= canonical_subspace_success) error stop 2
    if (.not. all(result%block_size == 2)) error stop 3
    if (.not. all(result%matched(2:,:))) error stop 4
    if (.not. all(result%similarity(2:,:) > 0.99_dp)) error stop 5
    if (.not. all(result%track_id == 1)) error stop 6
    allocate(tf_spec%frequencies(5),tf_spec%times(2),tf_spec%power(5,2,4,4))
    tf_spec%frequencies = spec%frequencies
    tf_spec%times = [1.0_dp,2.0_dp]
    tf_spec%power(:,1,:,:) = spec%power
    tf_spec%power(:,2,:,:) = 2.0_dp*spec%power
    tf_spec%status = 0
    tf_modes = canonical_coherence_tf_matrix(tf_spec,[1,2],[3,4],return_vectors=.true.)
    if (tf_modes%status /= 0) error stop 7
    tf_result = track_canonical_subspaces(tf_modes)
    if (tf_result%status /= canonical_subspace_success) error stop 8
    if (.not. all(tf_result%matched(2:,:,:))) error stop 9
    if (.not. all(shape(tf_result%basis_a) == [5,2,2,2])) error stop 10
    print *, 'PASS public canonical subspace tracking'
end program test_signal_canonical_subspace_tracking_public
