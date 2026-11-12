! SPDX-License-Identifier: MIT
! Automatic LGMRES Jacobian-image pruning and numerical dependence.
program test_optimize_krylov_augmentation_auto
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : retain_lgmres_image_direction
    implicit none

    real(dp) :: stored(4, 3)
    real(dp) :: old(4, 3)
    real(dp) :: images(4, 3)
    real(dp), parameter :: e1(4) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: e2(4) = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    real(dp), parameter :: e3(4) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: e4(4) = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
    integer :: count
    integer :: fails
    logical :: added

    fails = 0
    stored = 0.0_dp
    images = 0.0_dp
    count = 0

    call retain_lgmres_image_direction(stored, count, e1, images, e1, &
        0.0_dp, added, auto_threshold=.true.)
    call check(added .and. count == 1 .and. all(stored(:, 1) == e1), &
        'first direction accepted with automatic threshold')
    images(:, 1) = e1

    old = stored
    call retain_lgmres_image_direction(stored, count, e2, images, &
        e1 + 1.0e-10_dp*e2, 0.0_dp, added, auto_threshold=.true.)
    call check(.not. added .and. count == 1 .and. all(stored == old), &
        'automatically rejects a nearly dependent image')

    ! This separation is larger than the projection-dependent automatic
    ! tolerance. A fixed 1e-5 cutoff would reject it, but auto retains it.
    call retain_lgmres_image_direction(stored, count, e2, images, &
        e1 + 1.0e-6_dp*e2, 0.0_dp, added, auto_threshold=.true.)
    call check(added .and. count == 2 .and. all(stored(:, 2) == e2), &
        'automatically retains a distinguishable image')
    images(:, 2) = e1 + 1.0e-6_dp*e2
    old = stored
    call retain_lgmres_image_direction(stored, count, e3, images, &
        e2, 0.0_dp, added, auto_threshold=.true.)
    call check(.not. added .and. all(stored == old), &
        'rejects dependence on a combination of stored images')

    call retain_lgmres_image_direction(stored, count, e3, images, &
        e3, 0.0_dp, added, auto_threshold=.true.)
    call check(added .and. count == 3 .and. all(stored(:, 3) == e3), &
        'accepts orthogonal image')
    images(:, 3) = e3
    call retain_lgmres_image_direction(stored, count, e4, images, &
        e4, 0.0_dp, added, auto_threshold=.true.)
    call check(added .and. count == 3 .and. all(stored(:, 1) == e2) .and. &
        all(stored(:, 2) == e3) .and. all(stored(:, 3) == e4), &
        'independent new image evicts oldest direction')

    ! The fixed-tolerance path remains entirely independent of automatic
    ! filtering, and a caller may explicitly disable automatic filtering.
    stored = 0.0_dp
    images = 0.0_dp
    count = 1
    stored(:, 1) = e1
    images(:, 1) = e1
    call retain_lgmres_image_direction(stored, count, e2, images, &
        e1 + 1.0e-10_dp*e2, 1.0e-12_dp, added, auto_threshold=.false.)
    call check(added .and. count == 2, 'fixed tolerance retains smaller independent component')
    old = stored
    call retain_lgmres_image_direction(stored, count, e3, images, &
        e3, 0.0_dp, added, auto_threshold=.false.)
    call check(.not. added .and. all(old == stored), 'fixed zero disables the helper')
    call retain_lgmres_image_direction(stored, count, e3, images, &
        e3, 1.0e-8_dp, added, auto_threshold=.true.)
    call check(.not. added .and. all(old == stored), 'automatic mode rejects positive fixed tolerance')

    ! Nonfinite images must never enter the retained history. Scaling must
    ! also handle subnormal and huge but finite direction/image components.
    call retain_lgmres_image_direction(stored, count, e3, images, &
        [ieee_value(1.0_dp, ieee_quiet_nan), 0.0_dp, 0.0_dp, 0.0_dp], &
        0.0_dp, added, auto_threshold=.true.)
    call check(.not. added .and. all(old == stored), 'nonfinite image rejected')
    stored = 0.0_dp
    images = 0.0_dp
    count = 0
    call retain_lgmres_image_direction(stored, count, &
        [huge(1.0_dp)/4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], images, &
        [huge(1.0_dp)/4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        0.0_dp, added, auto_threshold=.true.)
    call check(added .and. count == 1, 'huge finite direction accepted')
    images(:, 1) = e1
    call retain_lgmres_image_direction(stored, count, &
        [0.0_dp, tiny(1.0_dp)/8.0_dp, 0.0_dp, 0.0_dp], images, &
        [0.0_dp, tiny(1.0_dp)/8.0_dp, 0.0_dp, 0.0_dp], &
        0.0_dp, added, auto_threshold=.true.)
    call check(added .and. count == 2 .and. stored(2, 2) == 1.0_dp, &
        'subnormal finite image accepted')

    if (fails /= 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation_auto: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation_auto: PASS'
contains
    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected test outcome
        character(len=*), intent(in) :: label !! description

        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check
end program test_optimize_krylov_augmentation_auto
