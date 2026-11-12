! SPDX-License-Identifier: MIT
! The incremental current-Jacobian image-span filter must match reconstruction.
program test_optimize_krylov_augmentation_stream
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_image_pruning_stream, retain_lgmres_image_direction_stream, &
        retain_lgmres_image_direction
    implicit none

    integer, parameter :: n = 48
    integer, parameter :: capacity = 40
    real(dp) :: directions(n, capacity)
    real(dp) :: jac_images(n, capacity)
    real(dp) :: old_vectors(n, capacity)
    real(dp) :: old_images(n, capacity)
    real(dp) :: streamed_vectors(n, capacity)
    real(dp) :: basis(n, capacity)
    real(dp) :: tolerance
    type(lgmres_pruning_workspace) :: workspace
    integer :: count_old, count_new, i, j, trial, failures
    logical :: old_added, new_added, automatic

    failures = 0
    do trial = 1, 8
        automatic = mod(trial, 2) == 0
        tolerance = 1.0e-8_dp
        if (automatic) tolerance = 0.0_dp
        directions = 0.0_dp
        jac_images = 0.0_dp
        do j = 1, capacity
            do i = 1, n
                directions(i, j) = sin(real(7*i+13*j, dp)) + cos(real(i*j+5, dp))
                jac_images(i, j) = sin(real(3*i+11*j, dp)) + cos(real(2*i*j+17, dp))
            end do
            if (mod(j, 5) == 0) then
                jac_images(:, j) = jac_images(:, j - 1) + 1.0e-10_dp * &
                    jac_images(:, max(1, j - 2))
            end if
            if (mod(j, 7) == 0) jac_images(:, j) = 0.0_dp
            if (mod(j, 9) == 0) directions(:, j) = 0.0_dp
        end do
        if (trial >= 5) then
            jac_images(:, 1) = jac_images(:, 1) * (huge(1.0_dp) / 16.0_dp)
            jac_images(:, 2) = jac_images(:, 2) * (tiny(1.0_dp) / 16.0_dp)
        end if
        old_vectors = 0.0_dp
        streamed_vectors = 0.0_dp
        old_images = 0.0_dp
        count_old = 0
        count_new = 0
        basis = 0.0_dp
        call reset_lgmres_image_pruning_stream(workspace)
        do j = 1, capacity
            call retain_lgmres_image_direction(old_vectors, count_old, directions(:, j), &
                old_images, jac_images(:, j), tolerance, old_added, auto_threshold=automatic)
            if (old_added) old_images(:, count_old) = jac_images(:, j)
            if (mod(trial, 4) <= 1) then
                call retain_lgmres_image_direction_stream(streamed_vectors, count_new, &
                    directions(:, j), jac_images(:, j), tolerance, new_added, &
                    workspace, auto_threshold=automatic)
            else
                call retain_lgmres_image_direction_stream(streamed_vectors, count_new, &
                    directions(:, j), jac_images(:, j), tolerance, new_added, &
                    workspace, auto_threshold=automatic, basis_work=basis)
            end if
            call check(old_added .eqv. new_added, 'same retained decision')
            call check(count_old == count_new, 'same retained count')
            call check(all(old_vectors == streamed_vectors), 'bitwise identical retained directions')
        end do
        ! Restart with the same workspace but a changed Jacobian. No old
        ! orthogonal basis is permitted to survive the reset.
        call reset_lgmres_image_pruning_stream(workspace)
        old_vectors = 0.0_dp
        streamed_vectors = 0.0_dp
        old_images = 0.0_dp
        count_old = 0
        count_new = 0
        do j = capacity, 1, -1
            call retain_lgmres_image_direction(old_vectors, count_old, directions(:, j), &
                old_images, jac_images(:, j), tolerance, old_added, auto_threshold=automatic)
            if (old_added) old_images(:, count_old) = jac_images(:, j)
            call retain_lgmres_image_direction_stream(streamed_vectors, count_new, &
                directions(:, j), jac_images(:, j), tolerance, new_added, &
                workspace, auto_threshold=automatic)
            call check(old_added .eqv. new_added, 'reset accepts same reverse order')
            call check(all(old_vectors == streamed_vectors), 'reset keeps identical stored directions')
        end do
    end do
    ! Below the 32*epsilon span floor, a direction may be retained without
    ! enlarging the numerical basis. The incremental path must treat the
    ! following candidates just as a fresh reconstruction would.
    old_vectors = 0.0_dp
    streamed_vectors = 0.0_dp
    old_images = 0.0_dp
    count_old = 0
    count_new = 0
    call reset_lgmres_image_pruning_stream(workspace)
    do j = 1, 3
        directions(:, j) = 0.0_dp
        jac_images(:, j) = 0.0_dp
        directions(j, j) = 1.0_dp
        jac_images(1, j) = 1.0_dp
    end do
    jac_images(2, 2) = 1.0e-15_dp
    jac_images(2, 3) = 1.0e-14_dp
    do j = 1, 3
        call retain_lgmres_image_direction(old_vectors, count_old, directions(:, j), &
            old_images, jac_images(:, j), 1.0e-18_dp, old_added)
        if (old_added) old_images(:, count_old) = jac_images(:, j)
        call retain_lgmres_image_direction_stream(streamed_vectors, count_new, &
            directions(:, j), jac_images(:, j), 1.0e-18_dp, new_added, &
            workspace, stream_limit=3)
        call check(old_added .eqv. new_added, 'sub-floor image accepts same decision')
        call check(all(old_vectors == streamed_vectors), 'sub-floor image retains identical vectors')
    end do
    ! The stream only appends; it must not silently evict an old direction.
    call retain_lgmres_image_direction_stream(streamed_vectors, count_new, &
        directions(:, 4), jac_images(:, 4), 1.0e-18_dp, new_added, &
        workspace, stream_limit=count_new)
    call check(.not. new_added, 'full append-only stream rejects FIFO insertion')
    if (failures /= 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_stream: FAIL ', failures
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_stream: PASS'
contains
    subroutine check(condition, description)
        logical, intent(in) :: condition !! required equivalence
        character(len=*), intent(in) :: description !! regression name
        if (.not. condition) then
            failures = failures + 1
            if (failures <= 12) print '(a)', 'FAIL: '//description
        end if
    end subroutine check
end program test_optimize_krylov_augmentation_stream
