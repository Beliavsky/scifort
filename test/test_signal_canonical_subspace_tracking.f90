! SPDX-License-Identifier: MIT
! Near-degenerate subspace tracking, unitary alignment and split/merge resets.
program test_signal_canonical_subspace_tracking
    use scifort_kinds, only : dp
    use scifort_signal_canonical_coherence, only : canonical_coherence_result, &
        canonical_coherence_tf_result
    use scifort_signal_canonical_subspace_tracking, only : track_canonical_subspaces, &
        canonical_subspace_tracks_result, canonical_subspace_tracks_tf_result, &
        canonical_subspace_success, canonical_subspace_invalid_input
    implicit none
    type(canonical_coherence_result) :: can, broken, missing
    type(canonical_coherence_tf_result) :: tf_can
    type(canonical_subspace_tracks_result) :: out
    type(canonical_subspace_tracks_tf_result) :: tf_out
    complex(dp) :: rot(2,2), phase, overlap(2,2), eye(2,2)
    real(dp) :: theta
    integer :: f, i, t

    theta = 0.43_dp
    phase = cmplx(cos(0.7_dp),sin(0.7_dp),dp)
    rot(1,1) = cmplx(cos(theta),0.0_dp,dp)*phase
    rot(2,1) = cmplx(sin(theta),0.0_dp,dp)*phase
    rot(1,2) = cmplx(-sin(theta),0.0_dp,dp)*conjg(phase)
    rot(2,2) = cmplx(cos(theta),0.0_dp,dp)*conjg(phase)
    eye = (0.0_dp,0.0_dp)
    eye(1,1) = (1.0_dp,0.0_dp)
    eye(2,2) = (1.0_dp,0.0_dp)
    call require(maxval(abs(matmul(transpose(conjg(rot)),rot)-eye)) < 1.0e-14_dp, &
        'test rotation unitary')

    allocate(can%frequencies(6),can%valid(6),can%squared_modes(6,3), &
        can%coefficients_a(6,3,3),can%coefficients_b(6,3,3), &
        can%directions_stable(6,3))
    can%frequencies = [0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    can%valid = .true.
    can%status = 0
    can%squared_modes(1,:) = [0.800000002_dp,0.8_dp,0.3_dp]
    can%squared_modes(2,:) = [0.800000003_dp,0.8_dp,0.3_dp]
    can%squared_modes(3,:) = [0.88_dp,0.72_dp,0.3_dp]
    can%squared_modes(4,:) = [0.800000002_dp,0.8_dp,0.3_dp]
    can%squared_modes(5,:) = can%squared_modes(4,:)
    can%squared_modes(6,:) = can%squared_modes(4,:)
    can%coefficients_a = (0.0_dp,0.0_dp)
    do f = 1, 6
        do i = 1, 3
            can%coefficients_a(f,i,i) = (1.0_dp,0.0_dp)
        end do
    end do
    can%coefficients_a(2,1:2,1:2) = rot
    can%coefficients_a(3,1:2,1:2) = rot
    can%coefficients_a(4,1:2,1:2) = rot
    can%coefficients_a(5,1:2,1:2) = rot
    can%coefficients_a(6,1:2,1:2) = rot
    can%coefficients_b = can%coefficients_a
    can%directions_stable = .false. ! old tracker would refuse every bin

    out = track_canonical_subspaces(can,gap_tolerance=0.001_dp,min_similarity=0.8_dp)
    call require(out%status == canonical_subspace_success, 'new API succeeds with tied roots')
    call require(all(out%block_size(1,:) == [2,2,1]), 'near-tie group dimension')
    call require(all(out%block_first(1,:) == [1,1,3]), 'block starts at strength rank')
    call require(all(out%track_id(1,:) == [1,1,2]), 'initial subspace IDs')
    call require(all(out%track_id(2,:) == [1,1,2]), 'degenerate subspace linked')
    call require(all(out%matched(2,:)), 'matched despite unstable individual modes')
    call require(maxval(abs(out%basis_a(2,1:2,1:2)-eye)) < 1.0e-12_dp, &
        'Procrustes undoes arbitrary complex unitary rotation')
    call require(maxval(abs(out%basis_b(2,1:2,1:2)-eye)) < 1.0e-12_dp, &
        'A and B use the same rotation')
    call require(minval(out%principal_cosines(2,:)) > 1.0_dp-1.0e-12_dp, &
        'identical subspaces have zero principal angles')
    call require(all(out%track_id(3,1:2) /= 1) .and. &
        all(out%segment_start(3,1:2)), 'splitting block creates new labels')
    call require(out%track_id(3,3) == 2 .and. out%matched(3,3), &
        'unaffected singleton keeps label across split')
    call require(out%track_id(4,1) == out%track_id(4,2), 'merged roots form block')
    call require(all(out%segment_start(4,1:2)), 'new label after merge')
    call require(out%track_id(4,3) == 2, 'unaffected mode survives merge')
    call require(all(out%matched(5:,:)), 'subspace links resume after merge')
    do f = 1, 6
        overlap = matmul(transpose(conjg(out%basis_a(f,1:2,1:2))), &
            out%basis_a(f,1:2,1:2))
        call require(maxval(abs(overlap-eye)) < 1.0e-12_dp, &
            'rotated subspace basis remains orthonormal')
        call require(maxval(abs(out%basis_a(f,:,:)-out%basis_b(f,:,:))) < 1.0e-14_dp, &
            'paired coefficient mixing preserved')
    end do
    call require(all(out%squared_modes == can%squared_modes), &
        'original canonical strengths preserved without reinterpretation')

    broken = can
    broken%frequencies(5) = -2.0_dp
    out = track_canonical_subspaces(broken)
    call require(all(out%segment_start(5,:)), 'FFT wrap creates new labels')
    call require(all(out%matched(6,:)), 'restart follows positive step after wrap')
    broken = can
    broken%valid(5) = .false.
    out = track_canonical_subspaces(broken)
    call require(all(out%track_id(5,:) == 0), 'invalid cell has no label')
    call require(all(out%segment_start(6,:)), 'invalid bin interrupts continuity')
    broken = can
    broken%squared_modes(2,3) = 0.0_dp
    out = track_canonical_subspaces(broken)
    call require(.not. out%matched(2,3), 'zero singular direction not linked')
    call require(.not. out%matched(3,3), 'zero singular predecessor not linked')
    call require(all(out%matched(2,1:2)), 'nonzero subspace unaffected by zero mode')
    missing = can
    deallocate(missing%coefficients_b)
    out = track_canonical_subspaces(missing)
    call require(out%status == canonical_subspace_invalid_input, 'missing vectors rejected')
    out = track_canonical_subspaces(can,gap_tolerance=-0.01_dp)
    call require(out%status == canonical_subspace_invalid_input, 'negative gap rejected')
    out = track_canonical_subspaces(can,min_similarity=1.01_dp)
    call require(out%status == canonical_subspace_invalid_input, 'similarity >1 rejected')
    out = track_canonical_subspaces(can,gap_tolerance=0.0_dp)
    call require(all(out%block_size == 1), 'zero threshold requires exactly equal roots')

    allocate(tf_can%frequencies(6),tf_can%times(2),tf_can%valid(6,2), &
        tf_can%squared_modes(6,2,3),tf_can%coefficients_a(6,2,3,3), &
        tf_can%coefficients_b(6,2,3,3))
    tf_can%frequencies = can%frequencies
    tf_can%times = [1.0_dp,2.0_dp]
    tf_can%valid = .true.
    tf_can%status = 0
    do t = 1, 2
        tf_can%squared_modes(:,t,:) = can%squared_modes
        tf_can%coefficients_a(:,t,:,:) = can%coefficients_a
        tf_can%coefficients_b(:,t,:,:) = can%coefficients_b
    end do
    tf_can%valid(3,2) = .false.
    tf_out = track_canonical_subspaces(tf_can)
    call require(tf_out%status == canonical_subspace_success, 'TF API succeeds')
    call require(all(tf_out%matched(2,:,:)), 'TF independent windows matched')
    call require(all(tf_out%track_id(3,2,:) == 0), 'invalid TF window has no labels')
    call require(all(tf_out%segment_start(4,2,:)), 'invalid TF window resets next bin')
    call require(tf_out%matched(3,1,3), 'other time window independent')
    call require(all(shape(tf_out%principal_cosines) == [6,2,3]), 'TF output shapes')
    print *, 'PASS canonical subspace tracking'

contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! assertion to evaluate
        character(*), intent(in) :: label !! failure diagnostic
        if (condition) return
        print *, 'FAIL: ',label
        error stop 1
    end subroutine require
end program test_signal_canonical_subspace_tracking
