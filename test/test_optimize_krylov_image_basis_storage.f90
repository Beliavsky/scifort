! SPDX-License-Identifier: MIT
! Compare image-basis storage in its own buffer and in an existing image array.
program test_optimize_krylov_image_basis_storage
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_image_pruning_stream, retain_lgmres_image_direction_stream, &
        retain_lgmres_image_candidate_cached
    implicit none
    integer, parameter :: n = 32
    integer, parameter :: k = 9
    real(dp) :: directions(n, k+2), products(n, k+2)
    real(dp) :: retained_old(n, k), retained_new(n, k)
    real(dp) :: shared_images(n, k)
    real(dp) :: independent_image(n)
    type(lgmres_pruning_workspace) :: original_workspace, shared_workspace
    real(dp) :: tolerance
    integer :: i, j, mode, trial, count_old, count_new, failures
    logical :: added_old, added_new, valid_old, valid_new, automatic

    failures = 0
    do j = 1, k+2
        do i = 1, n
            directions(i, j) = sin(real(7*i + 11*j, dp)) + &
                cos(0.13_dp * real(i*j, dp))
            products(i, j) = sin(real(17*i + 3*j, dp)) + &
                cos(0.09_dp * real(i*j + 2*j, dp))
        end do
    end do
    ! Exercise both full-rank FIFO history and a dependent retained image.
    independent_image = products(:, 5)
    do mode = 1, 2
        automatic = mode == 2
        tolerance = 1.0e-8_dp
        if (automatic) tolerance = 0.0_dp
        do trial = 1, 3
            products(:, 5) = independent_image
            if (trial == 2) products(:, 5) = products(:, 1) + 0.5_dp * products(:, 2)
            retained_old = 0.0_dp
            retained_new = 0.0_dp
            shared_images = 0.0_dp
            count_old = 0
            count_new = 0
            call reset_lgmres_image_pruning_stream(original_workspace)
            call reset_lgmres_image_pruning_stream(shared_workspace)
            do j = 1, k
                call retain_lgmres_image_direction_stream(retained_old, count_old, &
                    directions(:, j), products(:, j), tolerance, added_old, &
                    original_workspace, auto_threshold=automatic)
                call retain_lgmres_image_direction_stream(retained_new, count_new, &
                    directions(:, j), products(:, j), tolerance, added_new, &
                    shared_workspace, auto_threshold=automatic, basis_work=shared_images)
                call check(added_old .eqv. added_new, 'stream acceptance differs')
                call check(count_old == count_new, 'stream count differs')
                call check(all(retained_old == retained_new), 'retained columns differ')
            end do
            if (trial /= 2) call check(count_old == k, 'full-rank history expected')
            do j = 1, 2
                if (j == 1) then
                    ! Dependent candidate should be rejected.
                    directions(:, k+1) = directions(:, 1) + directions(:, 2)
                    products(:, k+1) = products(:, 1) + products(:, 2)
                end if
                call retain_lgmres_image_candidate_cached(retained_old, count_old, &
                    directions(:, k+j), products(:, k+j), tolerance, added_old, &
                    original_workspace, valid_old, auto_threshold=automatic)
                call retain_lgmres_image_candidate_cached(retained_new, count_new, &
                    directions(:, k+j), products(:, k+j), tolerance, added_new, &
                    shared_workspace, valid_new, auto_threshold=automatic, &
                    basis_work=shared_images)
                call check(valid_old .and. valid_new, 'candidate cache invalid')
                call check(added_old .eqv. added_new, 'candidate acceptance differs')
                call check(count_old == count_new, 'candidate count differs')
                call check(all(retained_old == retained_new), 'candidate retained vectors differ')
                ! Cache deliberately becomes invalid after the first candidate.
                if (j == 1) then
                    call reset_lgmres_image_pruning_stream(original_workspace)
                    call reset_lgmres_image_pruning_stream(shared_workspace)
                    count_old = 0
                    count_new = 0
                    retained_old = 0.0_dp
                    retained_new = 0.0_dp
                    shared_images = 0.0_dp
                    do i = 1, k
                        call retain_lgmres_image_direction_stream(retained_old, count_old, &
                            directions(:, i), products(:, i), tolerance, added_old, &
                            original_workspace, auto_threshold=automatic)
                        call retain_lgmres_image_direction_stream(retained_new, count_new, &
                            directions(:, i), products(:, i), tolerance, added_new, &
                            shared_workspace, auto_threshold=automatic, basis_work=shared_images)
                    end do
                end if
            end do
        end do
    end do
    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_krylov_image_basis_storage: FAIL ', failures
        error stop 1
    end if
    print '(a)', 'test_optimize_krylov_image_basis_storage: PASS'
contains
    subroutine check(condition, msg)
        logical, intent(in) :: condition !! requirement to verify
        character(len=*), intent(in) :: msg !! failure label
        if (.not. condition) then
            failures = failures + 1
            print '(a)', 'FAIL: '//msg
        end if
    end subroutine check
end program test_optimize_krylov_image_basis_storage
