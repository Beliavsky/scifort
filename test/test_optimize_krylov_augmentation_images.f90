! SPDX-License-Identifier: MIT
! Direct tests for rank filtering in the current Jacobian-image space.
program test_optimize_krylov_augmentation_images
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : retain_lgmres_image_direction
    implicit none

    real(dp) :: stored(4, 3)
    real(dp) :: old(4, 3)
    real(dp) :: images(4, 3)
    real(dp) :: direction(4)
    real(dp) :: nan_value
    integer :: count
    integer :: fails
    logical :: added

    fails = 0
    stored = 0.0_dp
    images = 0.0_dp
    count = 0
    call retain_lgmres_image_direction(stored, count, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        images, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(added .and. count == 1 .and. stored(1, 1) == 1.0_dp, &
        'first image retained')
    images(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    old = stored

    ! Solution directions e1 and e2 are orthogonal, but their images under
    ! [1, 1, 0, 0; 0, 1e-10, 0, 0; ...] are nearly parallel.
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images, [1.0_dp, 1.0e-10_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(.not. added .and. count == 1 .and. all(stored == old), &
        'solution-independent but image-dependent direction rejected')

    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images, [1.0_dp, 1.0e-4_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(added .and. count == 2 .and. stored(2, 2) == 1.0_dp, &
        'resolvable image direction accepted and original solution vector saved')
    images(:, 2) = [1.0_dp, 1.0e-4_dp, 0.0_dp, 0.0_dp]
    old = stored
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], &
        images, [2.0_dp, 1.0e-4_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(.not. added .and. all(stored == old), &
        'image within two-column span rejected')

    call retain_lgmres_image_direction(stored, count, [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], &
        images, [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(added .and. count == 3 .and. stored(3, 3) == 1.0_dp, &
        'third independent image retained')
    images(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    old = stored
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 0.0_dp, 0.0_dp, 8.0_dp], &
        images, [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], 1.0e-8_dp, added)
    call check(added .and. count == 3 .and. all(stored(:, 1) == old(:, 2)) .and. &
        all(stored(:, 2) == old(:, 3)) .and. stored(4, 3) == 1.0_dp, &
        'FIFO eviction preserves the solution directions')

    ! Image from a changed Jacobian can make an old independent pair
    ! dependent. The caller rebuilds retained history at each Newton point.
    stored = 0.0_dp
    images = 0.0_dp
    count = 0
    call retain_lgmres_image_direction(stored, count, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        images, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    images(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(added .and. count == 2, 'directions initially independent')
    images(:, 2) = [1.0_dp, 1.0e-10_dp, 0.0_dp, 0.0_dp]
    old = stored
    stored = 0.0_dp
    count = 0
    call retain_lgmres_image_direction(stored, count, old(:, 1), images, &
        images(:, 1), 1.0e-8_dp, added)
    call retain_lgmres_image_direction(stored, count, old(:, 2), images, &
        images(:, 2), 1.0e-8_dp, added)
    call check(.not. added .and. count == 1 .and. stored(1, 1) == 1.0_dp, &
        'changed Jacobian permits stale-image pruning')

    ! Neither zero Jacobian images nor invalid inputs may corrupt history.
    old = stored
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(.not. added .and. all(stored == old), 'null Jacobian image rejected')
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    direction = [0.0_dp, nan_value, 0.0_dp, 0.0_dp]
    call retain_lgmres_image_direction(stored, count, direction, images, &
        [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(.not. added .and. all(stored == old), 'nonfinite direction rejected')
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images, [0.0_dp, nan_value, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(.not. added .and. all(stored == old), 'nonfinite image rejected')
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], -1.0_dp, added)
    call check(.not. added .and. all(stored == old), 'invalid tolerance rejected')
    call retain_lgmres_image_direction(stored, count, [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        images(1:3, :), [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(.not. added .and. all(stored == old), 'wrong image shape rejected')

    ! Scaling both the directions and their images avoids overflow/underflow.
    stored = 0.0_dp
    images = 0.0_dp
    count = 0
    call retain_lgmres_image_direction(stored, count, &
        [huge(1.0_dp) / 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        images, [huge(1.0_dp) / 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(added .and. stored(1, 1) == 1.0_dp, 'huge direction/image scale')
    images(:, 1) = [huge(1.0_dp) / 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call retain_lgmres_image_direction(stored, count, &
        [0.0_dp, tiny(1.0_dp) / 4.0_dp, 0.0_dp, 0.0_dp], &
        images, [0.0_dp, tiny(1.0_dp) / 4.0_dp, 0.0_dp, 0.0_dp], 1.0e-8_dp, added)
    call check(added .and. count == 2 .and. stored(2, 2) == 1.0_dp, &
        'subnormal direction/image scale')

    if (fails > 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_images: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_images: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected outcome
        character(len=*), intent(in) :: label !! test case

        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

end program test_optimize_krylov_augmentation_images
