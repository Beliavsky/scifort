! SPDX-License-Identifier: MIT
! Regression for optimal frequency matching, common complex phase, and resets.
program test_signal_canonical_tracking
    use scifort_kinds, only : dp
    use scifort_signal_multitaper, only : multitaper_matrix_result, &
        multitaper_spectral_matrix_spectrogram_result
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result, canonical_coherence_matrix, canonical_coherence_tf_matrix
    use scifort_signal_canonical_tracking, only : track_canonical_modes, &
        canonical_mode_tracks_result, canonical_mode_tracks_tf_result, &
        canonical_tracking_success, canonical_tracking_invalid_input
    implicit none
    type(multitaper_matrix_result) :: spec
    type(multitaper_spectral_matrix_spectrogram_result) :: tf_spec
    type(canonical_coherence_result) :: can, alt, no_vectors, difficult
    type(canonical_coherence_tf_result) :: tf_can
    type(canonical_mode_tracks_result) :: tracks
    type(canonical_mode_tracks_tf_result) :: tf_tracks
    real(dp) :: first(6), second(6), angle
    complex(dp) :: cross, pair_phase
    integer :: i, f, t, k

    first = [0.92_dp,0.82_dp,0.62_dp,0.49_dp,0.35_dp,0.22_dp]
    second = [0.18_dp,0.28_dp,0.68_dp,0.78_dp,0.86_dp,0.91_dp]
    allocate(spec%frequencies(6),spec%power(6,4,4))
    spec%status = 0
    spec%frequencies = [(real(i-1,dp),i=1,6)]
    spec%power = (0.0_dp,0.0_dp)
    do f = 1, 6
        do i = 1, 4
            spec%power(f,i,i) = (1.0_dp,0.0_dp)
        end do
        angle = 0.15_dp*f
        cross = first(f)*cmplx(cos(angle),sin(angle),dp)
        spec%power(f,1,3) = cross
        spec%power(f,3,1) = conjg(cross)
        cross = second(f)*cmplx(cos(-0.12_dp*f),sin(-0.12_dp*f),dp)
        spec%power(f,2,4) = cross
        spec%power(f,4,2) = conjg(cross)
    end do
    can = canonical_coherence_matrix(spec,[1,2],[3,4],return_vectors=.true.)
    call require(can%status == 0 .and. all(can%valid), 'valid canonical input')
    call require(all(can%directions_stable), 'isolated singular modes')
    no_vectors = canonical_coherence_matrix(spec,[1,2],[3,4])
    tracks = track_canonical_modes(no_vectors)
    call require(tracks%status == canonical_tracking_invalid_input, 'vectors required')
    tracks = track_canonical_modes(can,min_overlap=0.7_dp)
    call require(tracks%status == canonical_tracking_success, 'tracking status')
    call require(all(tracks%source_mode(1,:) == [1,2]), 'initial ordering')
    call require(all(tracks%source_mode(2,:) == [1,2]), 'no switch before crossing')
    call require(all(tracks%source_mode(3,:) == [2,1]), 'modes exchanged at crossing')
    call require(all(tracks%source_mode(6,:) == [2,1]), 'preserve labels after crossing')
    call require(all(tracks%matched(2:,:)), 'match all subsequent bins')
    call require(tracks%segment_start(1) .and. .not. any(tracks%segment_start(2:)), &
        'first bin starts segment')
    call require(maxval(abs(tracks%squared_modes(:,1)-first**2)) < 3.0e-13_dp, &
        'first coherent mode after crossing')
    call require(maxval(abs(tracks%squared_modes(:,2)-second**2)) < 3.0e-13_dp, &
        'second coherent mode after crossing')
    call require(all(tracks%overlap(2:,:) > 0.98_dp), 'cross-frequency normalized overlap')
    do f = 2, 6
        do k = 1, 2
            pair_phase = sum(conjg(tracks%coefficients_a(f-1,:,k))* &
                tracks%coefficients_a(f,:,k)) + &
                sum(conjg(tracks%coefficients_b(f-1,:,k))* &
                tracks%coefficients_b(f,:,k))
            call require(abs(aimag(pair_phase)) < 1.0e-12_dp, 'paired phase continuity')
            call require(real(pair_phase,dp) > 0.0_dp, 'positive overlap gauge')
            call require(abs(sum(conjg(tracks%coefficients_a(f,:,k))* &
                matmul(spec%power(f,1:2,3:4),tracks%coefficients_b(f,:,k))) - &
                sqrt(tracks%squared_modes(f,k))) < 1.0e-12_dp, &
                'shared phase preserves positive cross spectrum')
        end do
    end do
    call require(all(can%squared_modes(:,1) >= can%squared_modes(:,2)), &
        'original descending API unmodified')

    ! Wrapping FFT indices must split independent contiguous segments.
    alt = can
    alt%frequencies = [0.0_dp,1.0_dp,2.0_dp,-3.0_dp,-2.0_dp,-1.0_dp]
    tracks = track_canonical_modes(alt)
    call require(tracks%segment_start(4), 'wrapped frequency reset')
    call require(.not. any(tracks%matched(4,:)), 'no match across wrap')
    call require(all(tracks%matched(5,:)), 'restart after wrap')

    ! An ambiguous mode, invalid bin, or insufficient similarity breaks labels.
    alt = can
    alt%directions_stable(3,2) = .false.
    tracks = track_canonical_modes(alt)
    call require(tracks%segment_start(3) .and. tracks%segment_start(4), &
        'unstable root breaks both links')
    call require(.not. any(tracks%matched(3:4,:)), 'no match adjacent to repeated root')
    call require(all(tracks%matched(5,:)), 'resume after unstable root')
    alt = can
    alt%valid(3) = .false.
    tracks = track_canonical_modes(alt)
    call require(all(tracks%source_mode(3,:) == 0), 'invalid bin source sentinel')
    call require(tracks%segment_start(4), 'resume following invalid bin')
    alt = can
    alt%coefficients_a(3,:,:) = cmplx(0.0_dp,0.0_dp,dp)
    tracks = track_canonical_modes(alt)
    call require(.not. any(tracks%matched(3:4,:)), 'zero coefficient reset')
    tracks = track_canonical_modes(can,min_overlap=1.01_dp)
    call require(tracks%status == canonical_tracking_invalid_input, 'reject threshold >1')
    tracks = track_canonical_modes(can,min_overlap=-0.01_dp)
    call require(tracks%status == canonical_tracking_invalid_input, 'reject negative threshold')
    tracks = track_canonical_modes(can,min_overlap=1.0_dp)
    call require(tracks%status == canonical_tracking_success, 'valid threshold endpoint')

    ! In three dimensions an assignment chosen greedily row-by-row fails.
    ! The globally optimal one-to-one matching is the identity permutation,
    ! even though row 1 individually prefers column 3.
    allocate(difficult%frequencies(2),difficult%valid(2), &
        difficult%squared_modes(2,3),difficult%directions_stable(2,3), &
        difficult%coefficients_a(2,3,3),difficult%coefficients_b(2,3,3))
    difficult%frequencies = [0.0_dp,1.0_dp]
    difficult%valid = .true.
    difficult%directions_stable = .true.
    difficult%squared_modes(1,:) = [0.9_dp,0.7_dp,0.5_dp]
    difficult%squared_modes(2,:) = [0.85_dp,0.65_dp,0.45_dp]
    difficult%coefficients_a = (0.0_dp,0.0_dp)
    do i = 1, 3
        difficult%coefficients_a(1,i,i) = (1.0_dp,0.0_dp)
    end do
    difficult%coefficients_a(2,:,:) = cmplx(reshape([ &
        -0.61305609981128018_dp, 0.69175755382320481_dp, &
        -0.38161984384032921_dp,-0.4707460503896756_dp, &
        -0.70777386291222899_dp,-0.52673932359452935_dp, &
        -0.63447645703107247_dp,-0.14327472110193437_dp, &
         0.75954722023483201_dp], [3,3]), 0.0_dp,dp)
    difficult%coefficients_b = difficult%coefficients_a
    tracks = track_canonical_modes(difficult,min_overlap=0.6_dp)
    call require(all(tracks%source_mode(2,:) == [1,2,3]), 'optimal global assignment')
    call require(all(tracks%matched(2,:)), 'optimal assignment meets minimum overlap')
    call require(maxval(abs(tracks%overlap(2,:) - &
        [0.61305609981128018_dp,0.70777386291222899_dp,0.75954722023483201_dp])) &
        < 1.0e-12_dp, 'global matching scores')
    tracks = track_canonical_modes(difficult,min_overlap=0.9_dp)
    call require(tracks%segment_start(2) .and. .not. any(tracks%matched(2,:)), &
        'weak overlap resets entire mode group')

    allocate(tf_spec%frequencies(6),tf_spec%times(2),tf_spec%power(6,2,4,4))
    tf_spec%frequencies = spec%frequencies
    tf_spec%times = [0.0_dp,1.0_dp]
    tf_spec%status = 0
    do t = 1, 2
        tf_spec%power(:,t,:,:) = real(t,dp)*spec%power
    end do
    tf_can = canonical_coherence_tf_matrix(tf_spec,[1,2],[3,4],return_vectors=.true.)
    call require(tf_can%status == 0, 'time frequency estimate valid')
    tf_tracks = track_canonical_modes(tf_can)
    call require(tf_tracks%status == canonical_tracking_success, 'TF tracking status')
    call require(all(shape(tf_tracks%source_mode) == [6,2,2]), 'TF tracking shape')
    call require(all(tf_tracks%source_mode(:,1,:) == tf_tracks%source_mode(:,2,:)), &
        'independent windows have same tracked root labels')
    call require(all(tf_tracks%matched(2:,:,:)), 'time-frequency trajectories')
    tf_can%directions_stable(3,2,:) = .false.
    tf_tracks = track_canonical_modes(tf_can)
    call require(tf_tracks%segment_start(3,2) .and. tf_tracks%segment_start(4,2), &
        'second window restarts on ambiguity')
    call require(all(tf_tracks%matched(3,1,:)), 'first window unaffected')
    print *, 'PASS canonical frequency-mode tracking'

contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! requirement to check
        character(*), intent(in) :: label !! description of assertion
        if (condition) return
        print *, 'FAIL:',label
        error stop 1
    end subroutine require
end program test_signal_canonical_tracking
