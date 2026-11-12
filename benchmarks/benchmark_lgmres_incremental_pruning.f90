! SPDX-License-Identifier: MIT
! Time repeated filtering of one current-Jacobian LGMRES augmentation history.
! Standalone: compile with src/kinds.f90 and src/optimize/krylov_augmentation.f90.
program benchmark_lgmres_incremental_pruning
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_image_pruning_stream, retain_lgmres_image_direction_stream, &
        retain_lgmres_image_direction
    implicit none

    integer, parameter :: n = 192
    integer, parameter :: history = 96
    integer, parameter :: runs = 8
    real(dp) :: directions(n, history)
    real(dp) :: images(n, history)
    real(dp) :: retained(n, history)
    real(dp) :: old_images(n, history)
    real(dp) :: basis(n, history)
    real(dp) :: t0, old_elapsed, new_elapsed
    integer :: i, j, rep, count, old_count
    logical :: added
    type(lgmres_pruning_workspace) :: work

    directions = 0.0_dp
    images = 0.0_dp
    do j = 1, history
        do i = 1, n
            directions(i, j) = sin(0.137_dp * real(i + 17 * j, dp))
            images(i, j) = 0.001_dp * cos(0.3_dp * real(i + j, dp))
        end do
        images(j, j) = images(j, j) + 1.0_dp
        if (mod(j, 8) == 0) images(:, j) = images(:, j - 1) + 1.0e-12_dp * images(:, j)
    end do

    call cpu_time(t0)
    old_count = 0
    do rep = 1, runs
        retained = 0.0_dp
        old_images = 0.0_dp
        old_count = 0
        do j = 1, history
            call retain_lgmres_image_direction(retained, old_count, directions(:, j), &
                old_images, images(:, j), 1.0e-8_dp, added, workspace=work, &
                basis_work=basis)
            if (added) old_images(:, old_count) = images(:, j)
        end do
    end do
    call cpu_time(old_elapsed)
    old_elapsed = old_elapsed - t0

    call cpu_time(t0)
    do rep = 1, runs
        retained = 0.0_dp
        count = 0
        call reset_lgmres_image_pruning_stream(work)
        do j = 1, history
            call retain_lgmres_image_direction_stream(retained, count, directions(:, j), &
                images(:, j), 1.0e-8_dp, added, work, basis_work=basis)
        end do
    end do
    call cpu_time(new_elapsed)
    new_elapsed = new_elapsed - t0
    if (count /= old_count) error stop 'stream and original disagree on retained count'
    print '(a,i0,a,i0,a,i0)', 'n=', n, ', proposed=', history, ', retained=', count
    print '(a,f10.5)', 'rebuild basis (s): ', old_elapsed
    print '(a,f10.5)', 'incremental  (s): ', new_elapsed
    if (new_elapsed > 0.0_dp) print '(a,f9.2)', 'speedup: ', old_elapsed / new_elapsed
end program benchmark_lgmres_incremental_pruning
