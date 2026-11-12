! SPDX-License-Identifier: MIT
! CPU-time comparison of repeated full-history LGMRES pruning decisions.
program benchmark_lgmres_fifo_pruning
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        retain_lgmres_direction, retain_lgmres_direction_stream, &
        retain_lgmres_direction_fifo, reset_lgmres_direction_pruning_stream
    implicit none

    integer, parameter :: n = 192
    integer, parameter :: capacity = 48
    integer, parameter :: proposals = 200
    integer, parameter :: repeats = 12
    real(dp), parameter :: tolerance = 1.0e-8_dp
    real(dp) :: vectors(n, proposals)
    real(dp) :: old_history(n, capacity)
    real(dp) :: new_history(n, capacity)
    real(dp) :: t0, t1, old_time, new_time
    type(lgmres_pruning_workspace) :: old_work
    type(lgmres_pruning_workspace) :: new_work
    integer :: rep, i, j, old_count, new_count
    logical :: old_added, new_added

    do j = 1, proposals
        do i = 1, n
            vectors(i, j) = sin(real(13*i*i + 7*i*j + 11*j*j, dp)) + &
                cos(real(3*i*i + 29*i*j + 5*j*j, dp))
        end do
        if (mod(j, 17) == 0) vectors(:, j) = vectors(:, j - 1)
    end do

    call cpu_time(t0)
    do rep = 1, repeats
        old_history = 0.0_dp
        old_count = 0
        do j = 1, proposals
            call retain_lgmres_direction(old_history, old_count, &
                vectors(:, j), tolerance, old_added, workspace=old_work)
        end do
    end do
    call cpu_time(t1)
    old_time = t1 - t0

    call cpu_time(t0)
    do rep = 1, repeats
        new_history = 0.0_dp
        new_count = 0
        call reset_lgmres_direction_pruning_stream(new_work)
        do j = 1, proposals
            if (new_count < capacity) then
                call retain_lgmres_direction_stream(new_history, new_count, &
                    vectors(:, j), tolerance, new_added, new_work)
            else
                call retain_lgmres_direction_fifo(new_history, new_count, &
                    vectors(:, j), tolerance, new_added, new_work)
            end if
        end do
    end do
    call cpu_time(t1)
    new_time = t1 - t0

    if (old_count /= new_count .or. .not. all(old_history == new_history)) then
        print '(a)', 'benchmark_lgmres_fifo_pruning: ERROR: final histories differ'
        stop 1
    end if
    print '(a,i0)', 'vector length = ', n
    print '(a,i0)', 'history capacity = ', capacity
    print '(a,i0)', 'proposals = ', proposals
    print '(a,i0)', 'retained = ', old_count
    print '(a,i0)', 'repetitions = ', repeats
    print '(a,es14.6)', 'legacy CPU seconds = ', old_time
    print '(a,es14.6)', 'FIFO QR CPU seconds = ', new_time
    print '(a,f10.2)', 'speedup = ', old_time / new_time
end program benchmark_lgmres_fifo_pruning
