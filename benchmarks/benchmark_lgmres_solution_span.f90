! SPDX-License-Identifier: MIT
! Standalone timing of append-only LGMRES solution-space pruning.
program benchmark_lgmres_solution_span
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        retain_lgmres_direction, retain_lgmres_direction_stream, &
        reset_lgmres_direction_pruning_stream
    implicit none

    integer, parameter :: n = 192
    integer, parameter :: capacity = 112
    integer, parameter :: proposals = 96
    integer, parameter :: repeats = 12
    real(dp), parameter :: tolerance = 1.0e-8_dp
    real(dp) :: vectors(n, proposals)
    real(dp) :: old_history(n, capacity)
    real(dp) :: new_history(n, capacity)
    real(dp) :: begin_time, end_time, old_seconds, new_seconds
    type(lgmres_pruning_workspace) :: workspace, legacy_work
    integer :: rep, i, j, old_count, new_count
    logical :: old_added, new_added

    do j = 1, proposals
        do i = 1, n
            vectors(i, j) = sin(real(i * 13 + 7 * j, dp)) + &
                cos(real(i * j + 3 * j, dp))
        end do
        if (mod(j, 7) == 0) vectors(:, j) = vectors(:, j - 1)
        if (mod(j, 11) == 0) vectors(:, j) = vectors(:, j - 1) + &
            1.0e-10_dp * vectors(:, j - 2)
        if (mod(j, 13) == 0) vectors(:, j) = 0.0_dp
    end do
    call cpu_time(begin_time)
    do rep = 1, repeats
        old_history = 0.0_dp
        old_count = 0
        do j = 1, proposals
            call retain_lgmres_direction(old_history, old_count, &
                vectors(:, j), tolerance, old_added, workspace=legacy_work)
        end do
    end do
    call cpu_time(end_time)
    old_seconds = end_time - begin_time

    call cpu_time(begin_time)
    do rep = 1, repeats
        new_history = 0.0_dp
        new_count = 0
        call reset_lgmres_direction_pruning_stream(workspace)
        do j = 1, proposals
            call retain_lgmres_direction_stream(new_history, new_count, &
                vectors(:, j), tolerance, new_added, workspace)
        end do
    end do
    call cpu_time(end_time)
    new_seconds = end_time - begin_time

    if (old_count /= new_count .or. .not. all(old_history == new_history)) then
        print '(a)', 'benchmark_lgmres_solution_span: ERROR: streams differ'
        stop 1
    end if
    print '(a,i0)', 'length = ', n
    print '(a,i0)', 'proposals = ', proposals
    print '(a,i0)', 'retained = ', old_count
    print '(a,i0)', 'repetitions = ', repeats
    print '(a,es14.6)', 'reconstruction CPU seconds = ', old_seconds
    print '(a,es14.6)', 'incremental CPU seconds = ', new_seconds
    print '(a,f10.2)', 'speedup = ', old_seconds / new_seconds
end program benchmark_lgmres_solution_span
