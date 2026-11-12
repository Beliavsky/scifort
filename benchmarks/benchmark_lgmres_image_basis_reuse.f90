! SPDX-License-Identifier: MIT
! Time current-Jacobian image pruning plus one candidate insertion.
! This excludes Jacobian-vector products and Arnoldi iteration.
program benchmark_lgmres_image_basis_reuse
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_image_pruning_stream, retain_lgmres_image_direction_stream, &
        retain_lgmres_image_direction, retain_lgmres_image_candidate_cached
    implicit none
    integer, parameter :: n = 192
    integer, parameter :: capacity = 48
    integer, parameter :: repetitions = 250
    real(dp) :: directions(n, capacity + 1)
    real(dp) :: images(n, capacity + 1)
    real(dp) :: saved(n, capacity)
    real(dp) :: old_images(n, capacity)
    real(dp) :: unused_arnoldi(n, capacity)
    real(dp) :: before
    real(dp) :: after
    real(dp) :: old_seconds
    real(dp) :: cached_seconds
    real(dp) :: tolerance
    type(lgmres_pruning_workspace) :: workspace
    integer :: mode
    integer :: trial
    integer :: j
    integer :: i
    integer :: count
    logical :: added
    logical :: valid
    logical :: automatic

    do j = 1, capacity + 1
        do i = 1, n
            directions(i, j) = sin(real(5*i + 11*j, dp)) + cos(real(i*j, dp))
            images(i, j) = sin(real(13*i + 7*j, dp)) + cos(real(i*j + j, dp))
        end do
    end do
    do mode = 1, 2
        automatic = mode == 2
        tolerance = 1.0e-9_dp
        if (automatic) tolerance = 0.0_dp
        call cpu_time(before)
        do trial = 1, repetitions
            count = 0
            saved = 0.0_dp
            old_images = 0.0_dp
            call reset_lgmres_image_pruning_stream(workspace)
            do j = 1, capacity
                call retain_lgmres_image_direction_stream(saved, count, &
                    directions(:, j), images(:, j), tolerance, added, workspace, &
                    auto_threshold=automatic, basis_work=unused_arnoldi)
                if (added) old_images(:, count) = images(:, j)
            end do
            call retain_lgmres_image_direction(saved, count, directions(:, capacity+1), &
                old_images, images(:, capacity+1), tolerance, added, &
                auto_threshold=automatic, workspace=workspace, &
                basis_work=unused_arnoldi)
        end do
        call cpu_time(after)
        old_seconds = after - before
        call cpu_time(before)
        do trial = 1, repetitions
            count = 0
            saved = 0.0_dp
            call reset_lgmres_image_pruning_stream(workspace)
            do j = 1, capacity
                call retain_lgmres_image_direction_stream(saved, count, &
                    directions(:, j), images(:, j), tolerance, added, workspace, &
                    auto_threshold=automatic)
            end do
            call retain_lgmres_image_candidate_cached(saved, count, &
                directions(:, capacity+1), images(:, capacity+1), tolerance, &
                added, workspace, valid, auto_threshold=automatic)
            if (.not. valid) error stop 'invalid cached candidate'
        end do
        call cpu_time(after)
        cached_seconds = after - before
        if (automatic) then
            print '(a)', 'mode=automatic'
        else
            print '(a)', 'mode=fixed'
        end if
        print '(a,i0)', 'repetitions=', repetitions
        print '(a,f10.6)', 'original_seconds=', old_seconds
        print '(a,f10.6)', 'cached_seconds=', cached_seconds
        print '(a,f10.3)', 'speedup=', old_seconds / max(cached_seconds, tiny(1.0_dp))
        print '(a,i0)', 'retained=', count
    end do
end program benchmark_lgmres_image_basis_reuse
