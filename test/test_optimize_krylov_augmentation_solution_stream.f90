! SPDX-License-Identifier: MIT
! Append-only solution-space pruning must agree with span reconstruction.
program test_optimize_krylov_augmentation_solution_stream
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_direction_pruning_stream, retain_lgmres_direction_stream, &
        retain_lgmres_direction
    implicit none

    integer, parameter :: n = 96
    integer, parameter :: capacity = 64
    integer, parameter :: proposals = 60
    type(lgmres_pruning_workspace) :: workspace
    real(dp) :: old_vectors(n, capacity)
    real(dp) :: new_vectors(n, capacity)
    real(dp) :: candidates(n, proposals)
    real(dp) :: vec(n)
    real(dp) :: threshold
    integer :: i, j, trial, old_count, new_count, failures
    logical :: accepted_old, accepted_new

    failures = 0
    do trial = 1, 9
        threshold = 1.0e-8_dp
        if (trial == 2 .or. trial == 7) threshold = 1.0e-13_dp
        if (trial == 3 .or. trial == 8) threshold = 1.0e-18_dp
        do j = 1, proposals
            do i = 1, n
                candidates(i, j) = sin(real(i*7 + j*13, dp)) + &
                    cos(real(i*j + 5*j, dp))
            end do
            if (mod(j, 5) == 0) candidates(:, j) = candidates(:, j - 1) + &
                1.0e-10_dp * candidates(:, j - 2)
            if (mod(j, 6) == 0) candidates(:, j) = candidates(:, j - 1)
            if (mod(j, 11) == 0) candidates(:, j) = 0.0_dp
            if (trial >= 4 .and. mod(j, 7) == 0) &
                candidates(:, j) = candidates(:, j) * (huge(1.0_dp) / 16.0_dp)
            if (trial >= 6 .and. mod(j, 13) == 0) &
                candidates(:, j) = candidates(:, j) * (tiny(1.0_dp) / 16.0_dp)
        end do
        old_vectors = 0.0_dp
        new_vectors = 0.0_dp
        old_count = 0
        new_count = 0
        call reset_lgmres_direction_pruning_stream(workspace)
        do j = 1, proposals
            vec = candidates(:, j)
            call retain_lgmres_direction(old_vectors, old_count, vec, &
                threshold, accepted_old)
            call retain_lgmres_direction_stream(new_vectors, new_count, vec, &
                threshold, accepted_new, workspace)
            call check(accepted_new .eqv. accepted_old, 'append retention flag')
            call check(new_count == old_count, 'append retained count')
            call check(all(new_vectors == old_vectors), 'bitwise retained direction values')
        end do
        ! A changed candidate order must not use an earlier stream basis.
        old_vectors = 0.0_dp
        new_vectors = 0.0_dp
        old_count = 0
        new_count = 0
        call reset_lgmres_direction_pruning_stream(workspace)
        do j = proposals, 1, -1
            vec = candidates(:, j)
            call retain_lgmres_direction(old_vectors, old_count, vec, &
                threshold, accepted_old)
            call retain_lgmres_direction_stream(new_vectors, new_count, vec, &
                threshold, accepted_new, workspace)
            call check(accepted_new .eqv. accepted_old, 'restarted retention flag')
            call check(all(new_vectors == old_vectors), 'restarted exact vectors')
        end do
    end do

    ! Directions closer than the 32*epsilon floor can be retained for a very
    ! small tolerance but must not enlarge the numerical orthogonal span.
    old_vectors = 0.0_dp
    new_vectors = 0.0_dp
    old_count = 0
    new_count = 0
    call reset_lgmres_direction_pruning_stream(workspace)
    do j = 1, 4
        vec = 0.0_dp
        vec(1) = 1.0_dp
        if (j >= 2) vec(2) = real(j - 1, dp) * 1.0e-15_dp
        call retain_lgmres_direction(old_vectors, old_count, vec, &
            1.0e-18_dp, accepted_old)
        call retain_lgmres_direction_stream(new_vectors, new_count, vec, &
            1.0e-18_dp, accepted_new, workspace)
        call check(accepted_old .eqv. accepted_new, 'sub-floor candidate')
        call check(all(new_vectors == old_vectors), 'sub-floor stored values')
    end do

    ! A full append-only stream never performs FIFO eviction. The caller
    ! switches to the original retention routine once capacity is reached.
    old_count = 0
    new_count = 0
    old_vectors = 0.0_dp
    new_vectors = 0.0_dp
    call reset_lgmres_direction_pruning_stream(workspace)
    do j = 1, capacity
        vec = 0.0_dp
        vec(j) = 1.0_dp
        call retain_lgmres_direction(old_vectors, old_count, vec, &
            1.0e-8_dp, accepted_old)
        call retain_lgmres_direction_stream(new_vectors, new_count, vec, &
            1.0e-8_dp, accepted_new, workspace)
        call check(accepted_old .eqv. accepted_new, 'fill capacity')
    end do
    vec = 0.0_dp
    vec(capacity + 1) = 1.0_dp
    call retain_lgmres_direction_stream(new_vectors, new_count, vec, &
        1.0e-8_dp, accepted_new, workspace)
    call check(.not. accepted_new .and. new_count == capacity, 'full stream does not evict')
    call retain_lgmres_direction(new_vectors, new_count, vec, &
        1.0e-8_dp, accepted_new)
    call check(accepted_new, 'legacy FIFO continuation')
    call check(all(new_vectors(:, 1) == old_vectors(:, 2)), 'FIFO removes oldest')

    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_solution_stream: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_solution_stream: PASS'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition !! required test property
        character(len=*), intent(in) :: description !! check label
        if (.not. condition) then
            failures = failures + 1
            if (failures <= 20) print '(a)', 'FAIL: '//description
        end if
    end subroutine check
end program test_optimize_krylov_augmentation_solution_stream
