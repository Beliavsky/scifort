! SPDX-License-Identifier: MIT
! Umbrella API: requested conditional spectra feed covariance-aware tracking.
program test_signal_canonical_covmetric_tracking_public
    use scifort_signal, only : dp, canonical_coherence_result, canonical_coherence_tf_result, &
        canonical_subspace_tracks_result, canonical_subspace_tracks_tf_result, &
        multitaper_canonical_coherence, multitaper_tf_canonical_coherence, &
        multitaper_band_canonical_coherence, multitaper_tf_band_canonical_coherence, &
        track_canonical_subspaces, canonical_subspace_success
    implicit none
    real(dp) :: x(128,5), bands(2,2)
    complex(dp) :: z(128,5)
    type(canonical_coherence_result) :: stat
    type(canonical_coherence_tf_result) :: tf
    type(canonical_subspace_tracks_result) :: tracks
    type(canonical_subspace_tracks_tf_result) :: tf_tracks
    integer :: i, j

    do j = 1, 5
        do i = 1, 128
            x(i,j) = sin(0.065_dp*i*j)+0.27_dp*cos(0.13_dp*i+0.17_dp*j)
            z(i,j) = cmplx(x(i,j),0.32_dp*cos(0.19_dp*i*j),dp)
        end do
    end do
    stat = multitaper_canonical_coherence(x,3.0_dp,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,shrinkage=0.3_dp, &
        return_vectors=.true.,return_covariances=.true.)
    if (stat%status /= 0 .or. .not. all(stat%valid)) error stop 'metric real canonical'
    if (.not. allocated(stat%covariance_a)) error stop 'metric real covariance output'
    tracks = track_canonical_subspaces(stat,covariance_metric=.true.,min_similarity=0.4_dp)
    if (tracks%status /= canonical_subspace_success) error stop 'metric real tracker'
    if (.not. tracks%covariance_metric_used) error stop 'metric selected'
    stat = multitaper_canonical_coherence(z,3.0_dp,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,weighting='adaptive',cv_shrinkage=.true., &
        ridge=0.1_dp,return_vectors=.true.,return_covariances=.true.)
    if (stat%status /= 0 .or. .not. all(stat%valid)) error stop 'metric complex canonical'
    tracks = track_canonical_subspaces(stat,covariance_metric=.true.)
    if (tracks%status /= canonical_subspace_success) error stop 'metric complex tracker'
    bands(:,1) = [0.01_dp,0.1_dp]
    bands(:,2) = [0.15_dp,0.25_dp]
    stat = multitaper_band_canonical_coherence(x,3.0_dp,bands,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,return_vectors=.true., &
        return_covariances=.true.)
    if (stat%status /= 0 .or. .not. all(stat%valid)) error stop 'metric real band'
    tracks = track_canonical_subspaces(stat,covariance_metric=.true.)
    if (tracks%status /= canonical_subspace_success) error stop 'metric band tracker'
    stat = multitaper_band_canonical_coherence(z,3.0_dp,bands,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,return_vectors=.true., &
        return_covariances=.true.)
    if (stat%status /= 0 .or. .not. all(stat%valid)) error stop 'metric complex band'
    tf = multitaper_tf_canonical_coherence(z,3.0_dp,64,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,shrinkage=0.3_dp, &
        return_vectors=.true.,return_covariances=.true.)
    if (tf%status /= 0 .or. .not. all(tf%valid)) error stop 'metric complex TF'
    tf_tracks = track_canonical_subspaces(tf,covariance_metric=.true.)
    if (tf_tracks%status /= canonical_subspace_success) error stop 'metric TF tracker'
    if (.not. tf_tracks%covariance_metric_used) error stop 'metric TF selected'
    tf = multitaper_tf_canonical_coherence(x,3.0_dp,64,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,shrinkage=0.3_dp, &
        return_vectors=.true.,return_covariances=.true.)
    if (tf%status /= 0 .or. .not. all(tf%valid)) error stop 'metric real TF'
    tf = multitaper_tf_band_canonical_coherence(z,3.0_dp,64,bands,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,return_vectors=.true., &
        return_covariances=.true.)
    if (tf%status /= 0 .or. .not. all(tf%valid)) error stop 'metric complex TF band'
    tf_tracks = track_canonical_subspaces(tf,covariance_metric=.true.)
    if (tf_tracks%status /= canonical_subspace_success) error stop 'metric TF band tracker'
    tf = multitaper_tf_band_canonical_coherence(x,3.0_dp,64,bands,[1,2],[3,4], &
        conditioned_on=[5],kmax=6,ridge=0.1_dp,return_vectors=.true., &
        return_covariances=.true.)
    if (tf%status /= 0 .or. .not. all(tf%valid)) error stop 'metric real TF band'
    print *, 'PASS covariance-metric canonical tracking public API'
end program test_signal_canonical_covmetric_tracking_public
