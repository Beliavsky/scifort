! SPDX-License-Identifier: BSD-3-Clause
! Test opt-in LGMRES augmentation-space redundancy filtering.
program test_optimize_krylov_augmentation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_optimize_krylov_augmentation, only : retain_lgmres_direction
    implicit none

    real(dp) :: retained(4, 3)
    real(dp) :: old(4, 3)
    real(dp) :: direction(4)
    real(dp) :: nan_value
    integer :: count
    integer :: fails
    logical :: added

    fails = 0
    retained = 0.0_dp
    count = 0
    call retain_lgmres_direction(retained, count, [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        1.0e-8_dp, added)
    call check(added .and. count == 1 .and. retained(1, 1) == 1.0_dp, &
        'first direction retained and normalized')
    call retain_lgmres_direction(retained, count, [0.0_dp, 2.0_dp, 0.0_dp, 0.0_dp], &
        1.0e-8_dp, added)
    call check(added .and. count == 2 .and. retained(2, 2) == 1.0_dp, &
        'independent second direction retained')
    old = retained

    ! A new direction can be a linear combination even when no individual
    ! retained vector is close to parallel (pairwise cosine about 0.707).
    call retain_lgmres_direction(retained, count, [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], &
        1.0e-8_dp, added)
    call check(.not. added .and. count == 2 .and. all(retained == old), &
        'two-vector linear combination is rejected')
    call retain_lgmres_direction(retained, count, [1.0_dp, 0.0_dp, 1.0e-10_dp, 0.0_dp], &
        1.0e-8_dp, added)
    call check(.not. added .and. count == 2 .and. all(retained == old), &
        'near-dependent direction rejected')
    call retain_lgmres_direction(retained, count, [1.0_dp, 0.0_dp, 1.0e-4_dp, 0.0_dp], &
        1.0e-8_dp, added)
    call check(added .and. count == 3 .and. retained(3, 3) > 0.0_dp, &
        'resolvable independent component retained')

    old = retained
    ! Capacity is FIFO, and stored solution directions are never reorthogonalized
    ! in place: they must remain the original augmentation directions.
    call retain_lgmres_direction(retained, count, [0.0_dp, 0.0_dp, 0.0_dp, 8.0_dp], &
        1.0e-8_dp, added)
    call check(added .and. count == 3 .and. retained(2, 1) == 1.0_dp .and. &
        retained(4, 3) == 1.0_dp .and. all(retained(:, 2) == old(:, 3)), &
        'FIFO discards oldest without changing retained span vectors')

    old = retained
    nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
    direction = 0.0_dp
    call retain_lgmres_direction(retained, count, direction, 1.0e-8_dp, added)
    call check(.not. added .and. count == 3 .and. all(retained == old), &
        'zero direction does not change history')
    direction(1) = nan_value
    call retain_lgmres_direction(retained, count, direction, 1.0e-8_dp, added)
    call check(.not. added .and. count == 3 .and. all(retained == old), &
        'nonfinite direction does not change history')
    call retain_lgmres_direction(retained, count, [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
        -1.0_dp, added)
    call check(.not. added .and. count == 3 .and. all(retained == old), &
        'invalid tolerance leaves history unchanged')

    ! Normalization must not overflow on a huge vector or underflow on a
    ! subnormal one. Reinitialize to leave empty independent directions.
    retained = 0.0_dp
    count = 0
    direction = [huge(1.0_dp) / 2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    call retain_lgmres_direction(retained, count, direction, 1.0e-8_dp, added)
    call check(added .and. count == 1 .and. retained(1, 1) == 1.0_dp, &
        'huge vector normalization')
    direction = [0.0_dp, tiny(1.0_dp) / 4.0_dp, 0.0_dp, 0.0_dp]
    call retain_lgmres_direction(retained, count, direction, 1.0e-8_dp, added)
    call check(added .and. count == 2 .and. retained(2, 2) == 1.0_dp, &
        'subnormal vector normalization')

    if (fails > 0) then
        print '(a,i0)', 'test_optimize_krylov_augmentation: FAIL ', fails
        stop 1
    end if
    print '(a)', 'test_optimize_krylov_augmentation: PASS'

contains

    subroutine check(condition, label)
        logical, intent(in) :: condition !! expected result
        character(len=*), intent(in) :: label !! case being tested

        if (.not. condition) then
            fails = fails + 1
            print '(a)', 'FAIL: '//label
        end if
    end subroutine check

end program test_optimize_krylov_augmentation
