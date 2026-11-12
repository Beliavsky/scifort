! SPDX-License-Identifier: MIT
! Ensure frequency tracking is exported and generic resolves for static and TF.
program test_signal_canonical_tracking_public
    use scifort_signal, only : dp, multitaper_canonical_coherence, &
        multitaper_tf_canonical_coherence, canonical_coherence_result, &
        canonical_coherence_tf_result, canonical_mode_tracks_result, &
        canonical_mode_tracks_tf_result, track_canonical_modes, &
        canonical_tracking_success, canonical_tracking_invalid_input
    implicit none
    type(canonical_coherence_result) :: canonical
    type(canonical_coherence_tf_result) :: moving
    type(canonical_mode_tracks_result) :: tracks
    type(canonical_mode_tracks_tf_result) :: tf_tracks
    real(dp) :: x(128,4)
    complex(dp) :: z(128,4)
    integer :: i, j

    do j = 1, 4
        do i = 1, 128
            x(i,j) = sin(0.07_dp*i*j) + 0.4_dp*cos(0.13_dp*i+0.3_dp*j)
            z(i,j) = cmplx(x(i,j),0.18_dp*cos(0.11_dp*i*j),dp)
        end do
    end do
    canonical = multitaper_canonical_coherence(x,3.0_dp,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp,return_vectors=.true.)
    if (canonical%status /= 0) error stop 'public static input failed'
    tracks = track_canonical_modes(canonical)
    if (tracks%status /= canonical_tracking_success) error stop 'public static tracking failed'
    if (size(tracks%source_mode,1) /= size(canonical%frequencies)) &
        error stop 'public static frequency axis'
    if (.not. all(tracks%valid)) error stop 'public static validity'
    canonical = multitaper_canonical_coherence(z,3.0_dp,[1,2],[3,4], &
        kmax=5,weighting='eigen',ridge=0.1_dp,return_vectors=.true.)
    tracks = track_canonical_modes(canonical,min_overlap=0.4_dp)
    if (tracks%status /= canonical_tracking_success) error stop 'public complex tracking failed'
    moving = multitaper_tf_canonical_coherence(z,3.0_dp,64,[1,2],[3,4], &
        kmax=5,ridge=0.08_dp,return_vectors=.true.)
    tf_tracks = track_canonical_modes(moving)
    if (tf_tracks%status /= canonical_tracking_success) error stop 'public TF tracking failed'
    if (any(shape(tf_tracks%source_mode) /= shape(moving%squared_modes))) &
        error stop 'public TF mode shapes'
    if (size(tf_tracks%times) /= size(moving%times)) &
        error stop 'public TF time coordinates'
    canonical = multitaper_canonical_coherence(x,3.0_dp,[1,2],[3,4], &
        kmax=5,ridge=0.05_dp)
    tracks = track_canonical_modes(canonical)
    if (tracks%status /= canonical_tracking_invalid_input) &
        error stop 'public missing vectors validation'
    print *, 'PASS public canonical mode tracking'
end program test_signal_canonical_tracking_public
