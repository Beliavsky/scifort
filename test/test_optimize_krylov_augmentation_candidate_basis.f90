! SPDX-License-Identifier: MIT
! Compare candidate screening against reconstruction of the current-Jacobian span.
program test_optimize_krylov_augmentation_candidate_basis
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : lgmres_pruning_workspace, &
        reset_lgmres_image_pruning_stream, retain_lgmres_image_direction_stream, &
        retain_lgmres_image_candidate_cached, retain_lgmres_image_direction
    implicit none

    integer, parameter :: n = 48
    integer, parameter :: capacity = 12
    real(dp) :: directions(n, capacity + 1)
    real(dp) :: image_vectors(n, capacity + 1)
    real(dp) :: legacy_vectors(n, capacity)
    real(dp) :: cached_vectors(n, capacity)
    real(dp) :: legacy_images(n, capacity)
    real(dp) :: candidate(n)
    real(dp) :: candidate_image(n)
    real(dp) :: tolerance
    real(dp) :: factor
    type(lgmres_pruning_workspace) :: workspace
    integer :: trial
    integer :: j
    integer :: i
    integer :: limit
    integer :: count_old
    integer :: count_new
    integer :: fails
    logical :: old_added
    logical :: new_added
    logical :: valid
    logical :: automatic

    fails = 0
    do trial = 1, 32
        automatic = mod(trial, 2) == 0
        tolerance = 1.0e-8_dp
        if (automatic) tolerance = 0.0_dp
        limit = mod(trial - 1, capacity) + 1
        legacy_vectors = 0.0_dp
        cached_vectors = 0.0_dp
        legacy_images = 0.0_dp
        count_old = 0
        count_new = 0
        call reset_lgmres_image_pruning_stream(workspace)
        do j = 1, capacity + 1
            do i = 1, n
                directions(i, j) = sin(real(7*i + 11*j + trial, dp)) + &
                    cos(real(5*i*j + trial, dp))
                image_vectors(i, j) = sin(real(3*i + 13*j + trial, dp)) + &
                    cos(real(9*i*j + trial, dp))
            end do
            if (mod(j, 4) == 0) then
                image_vectors(:, j) = image_vectors(:, j - 1) + &
                    1.0e-10_dp * image_vectors(:, max(1, j - 2))
            end if
        end do
        ! Exercise image normalization at both ends of the exponent range.
        if (mod(trial, 3) == 0) image_vectors(:, 1) = &
            image_vectors(:, 1) * (huge(1.0_dp) / 64.0_dp)
        if (mod(trial, 5) == 0) image_vectors(:, 2) = &
            image_vectors(:, 2) * (tiny(1.0_dp) / 64.0_dp)
        do j = 1, limit
            call retain_lgmres_image_direction(legacy_vectors, count_old, &
                directions(:, j), legacy_images, image_vectors(:, j), &
                tolerance, old_added, auto_threshold=automatic)
            if (old_added) legacy_images(:, count_old) = image_vectors(:, j)
            call retain_lgmres_image_direction_stream(cached_vectors, count_new, &
                directions(:, j), image_vectors(:, j), tolerance, new_added, &
                workspace, auto_threshold=automatic, stream_limit=limit)
            call check(old_added .eqv. new_added, 'pre-filter agrees')
            call check(count_old == count_new, 'pre-filter count agrees')
            call check(all(legacy_vectors == cached_vectors), 'pre-filter directions agree')
        end do
        candidate = directions(:, capacity + 1)
        candidate_image = image_vectors(:, capacity + 1)
        if (mod(trial, 7) == 0) candidate_image = image_vectors(:, 1)
        if (mod(trial, 9) == 0) candidate_image = 0.0_dp
        if (mod(trial, 11) == 0) candidate = 0.0_dp
        if (mod(trial, 13) == 0) then
            factor = huge(1.0_dp) / 128.0_dp
            candidate = candidate * factor
        end if
        call retain_lgmres_image_direction(legacy_vectors, count_old, &
            candidate, legacy_images, candidate_image, tolerance, old_added, &
            auto_threshold=automatic)
        call retain_lgmres_image_candidate_cached(cached_vectors, count_new, &
            candidate, candidate_image, tolerance, new_added, workspace, valid, &
            auto_threshold=automatic)
        call check(valid, 'current basis valid')
        call check(old_added .eqv. new_added, 'candidate decision agrees')
        call check(count_old == count_new, 'candidate count agrees')
        call check(all(legacy_vectors == cached_vectors), 'candidate directions agree')
        call retain_lgmres_image_candidate_cached(cached_vectors, count_new, &
            candidate, candidate_image, tolerance, new_added, workspace, valid, &
            auto_threshold=automatic)
        call check(.not. valid, 'one-shot basis cannot be reused twice')
    end do
    if (fails /= 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_candidate_basis: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_candidate_basis: PASS'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected outcome
        character(len=*), intent(in) :: label !! scenario name
        if (.not. condition) then
            fails = fails + 1
            if (fails <= 12) print '(a)', 'FAIL: '//label
        end if
    end subroutine check
end program test_optimize_krylov_augmentation_candidate_basis
