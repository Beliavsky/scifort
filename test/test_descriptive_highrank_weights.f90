! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0.

program test_descriptive_highrank_weights
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, hmean_keepdims, &
        pmean, pmean_keepdims
    implicit none

    real(dp), parameter :: tol = 5.0e-12_dp
    integer :: i
    real(dp) :: nan
    real(dp) :: scalar_weight
    real(dp) :: weights1(4)
    real(dp) :: weights2(4, 3)
    real(dp) :: weights3(1, 1, 3)
    real(dp) :: weights4(2, 1, 1, 3)
    real(dp) :: bad_weights2(2, 2)
    real(dp) :: bad_weights5(1, 1, 1, 1, 1)
    real(dp) :: x4(2, 1, 1, 3)
    real(dp) :: x5(2, 1, 1, 1, 1)
    real(dp) :: x6(2, 1, 1, 1, 1, 1)
    real(dp) :: x7(2, 1, 1, 1, 1, 1, 1)
    real(dp) :: x8(2, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x9(2, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x10(2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x11(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x12(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x13(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x14(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights5(2, 1, 1, 1, 1)
    real(dp) :: weights6(2, 1, 1, 1, 1, 1)
    real(dp) :: weights7(2, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights8(2, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights9(2, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights10(2, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights11(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights12(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights13(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights14(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: x15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: empty15(0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: weights15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: nan_weights15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: r3(:, :, :)
    real(dp), allocatable :: k4(:, :, :, :)
    real(dp), allocatable :: r14(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: k15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

    scalar_weight = 2.0_dp
    x5 = reshape([1.0_dp, 4.0_dp], shape(x5))
    x6 = reshape([1.0_dp, 4.0_dp], shape(x6))
    x7 = reshape([1.0_dp, 4.0_dp], shape(x7))
    x8 = reshape([1.0_dp, 4.0_dp], shape(x8))
    x9 = reshape([1.0_dp, 4.0_dp], shape(x9))
    x10 = reshape([1.0_dp, 4.0_dp], shape(x10))
    x11 = reshape([1.0_dp, 4.0_dp], shape(x11))
    x12 = reshape([1.0_dp, 4.0_dp], shape(x12))
    x13 = reshape([1.0_dp, 4.0_dp], shape(x13))
    x14 = reshape([1.0_dp, 4.0_dp], shape(x14))
    weights5 = 1.0_dp
    weights6 = 1.0_dp
    weights7 = 1.0_dp
    weights8 = 1.0_dp
    weights9 = 1.0_dp
    weights10 = 1.0_dp
    weights11 = 1.0_dp
    weights12 = 1.0_dp
    weights13 = 1.0_dp
    weights14 = 1.0_dp

    call check_close('rank5 weighted gmean', sum(gmean(x5, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank5 weighted hmean', sum(hmean(x5, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank5 weighted pmean', &
        sum(pmean(x5, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank5 full-rank weight branch', &
        sum(gmean(x5, weights=weights5)), 2.0_dp, tol)
    call check_close('rank6 weighted gmean', sum(gmean(x6, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank6 weighted hmean', sum(hmean(x6, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank6 weighted pmean', &
        sum(pmean(x6, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank6 full-rank weight branch', &
        sum(gmean(x6, weights=weights6)), 2.0_dp, tol)
    call check_close('rank7 weighted gmean', sum(gmean(x7, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank7 weighted hmean', sum(hmean(x7, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank7 weighted pmean', &
        sum(pmean(x7, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank7 full-rank weight branch', &
        sum(gmean(x7, weights=weights7)), 2.0_dp, tol)
    call check_close('rank8 weighted gmean', sum(gmean(x8, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank8 weighted hmean', sum(hmean(x8, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank8 weighted pmean', &
        sum(pmean(x8, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank8 full-rank weight branch', &
        sum(gmean(x8, weights=weights8)), 2.0_dp, tol)
    call check_close('rank9 weighted gmean', sum(gmean(x9, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank9 weighted hmean', sum(hmean(x9, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank9 weighted pmean', &
        sum(pmean(x9, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank9 full-rank weight branch', &
        sum(gmean(x9, weights=weights9)), 2.0_dp, tol)
    call check_close('rank10 weighted gmean', sum(gmean(x10, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank10 weighted hmean', sum(hmean(x10, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank10 weighted pmean', &
        sum(pmean(x10, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank10 full-rank weight branch', &
        sum(gmean(x10, weights=weights10)), 2.0_dp, tol)
    call check_close('rank11 weighted gmean', sum(gmean(x11, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank11 weighted hmean', sum(hmean(x11, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank11 weighted pmean', &
        sum(pmean(x11, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank11 full-rank weight branch', &
        sum(gmean(x11, weights=weights11)), 2.0_dp, tol)
    call check_close('rank12 weighted gmean', sum(gmean(x12, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank12 weighted hmean', sum(hmean(x12, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank12 weighted pmean', &
        sum(pmean(x12, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank12 full-rank weight branch', &
        sum(gmean(x12, weights=weights12)), 2.0_dp, tol)
    call check_close('rank13 weighted gmean', sum(gmean(x13, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank13 weighted hmean', sum(hmean(x13, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank13 weighted pmean', &
        sum(pmean(x13, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank13 full-rank weight branch', &
        sum(gmean(x13, weights=weights13)), 2.0_dp, tol)
    call check_close('rank14 weighted gmean', sum(gmean(x14, weights=scalar_weight)), 2.0_dp, tol)
    call check_close('rank14 weighted hmean', sum(hmean(x14, weights=scalar_weight)), 1.6_dp, tol)
    call check_close('rank14 weighted pmean', &
        sum(pmean(x14, 2.0_dp, weights=scalar_weight)), sqrt(8.5_dp), tol)
    call check_close('rank14 full-rank weight branch', &
        sum(gmean(x14, weights=weights14)), 2.0_dp, tol)

    x4 = reshape([(real(i, dp), i = 1, 6)], shape(x4))
    weights2 = reshape([(real(i, dp), i = 1, 12)], shape(weights2))
    weights3 = 1.0_dp
    weights4 = 1.0_dp
    r3 = gmean(x4, axis=1, weights=weights2)
    call check_true('rank4 joint broadcast shape', all(shape(r3) == [1, 4, 3]))
    call check_close('rank4 joint broadcast gmean first', r3(1, 1, 1), &
        1.4142135623730950_dp, tol)
    call check_close('rank4 joint broadcast gmean last', r3(1, 4, 3), &
        5.477225575051661_dp, tol)
    r3 = hmean(x4, axis=1, weights=weights2)
    call check_close('rank4 joint broadcast hmean', r3(1, 2, 2), &
        3.4285714285714284_dp, tol)
    r3 = pmean(x4, 2.0_dp, axis=1, weights=weights2)
    call check_close('rank4 joint broadcast pmean', r3(1, 3, 3), &
        5.522680508593631_dp, tol)
    k4 = gmean_keepdims(x4, axis=1, weights=weights2)
    call check_true('rank4 weighted keepdims shape', all(shape(k4) == [1, 1, 4, 3]))
    call check_close('rank4 weighted keepdims value', k4(1, 1, 4, 3), &
        5.477225575051661_dp, tol)
    r3 = gmean(x4, axis=1, weights=weights3)
    call check_true('rank4 rank3-weight shape', all(shape(r3) == [1, 1, 3]))
    call check_close('rank4 rank3-weight branch', r3(1, 1, 1), &
        1.4142135623730950_dp, tol)
    r3 = gmean(x4, axis=1, weights=weights4)
    call check_close('rank4 full-rank weight branch', r3(1, 1, 3), &
        5.477225575051661_dp, tol)

    bad_weights2 = 1.0_dp
    r3 = gmean(x4, axis=1, weights=bad_weights2)
    call check_true('rank4 invalid broadcast shape', all(shape(r3) == [1, 1, 3]))
    call check_true('rank4 invalid broadcast NaN', all(ieee_is_nan(r3)))
    bad_weights5 = 1.0_dp
    r3 = gmean(x4, axis=1, weights=bad_weights5)
    call check_true('rank4 excessive weight rank NaN', all(ieee_is_nan(r3)))

    x15 = reshape([(real(i, dp), i = 1, 6)], shape(x15))
    weights1 = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    r14 = gmean(x15, axis=2, weights=weights1)
    call check_true('rank15 trailing expansion shape', &
        all(shape(r14) == [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4]))
    call check_close('rank15 trailing expansion first', &
        r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 2.4662120743304703_dp, tol)
    call check_close('rank15 trailing expansion last', &
        r14(2,1,1,1,1,1,1,1,1,1,1,1,1,4), 3.634241185664279_dp, tol)

    weights15 = reshape([(real(i, dp), i = 1, 6)], shape(weights15))
    r14 = gmean(x15, axis=2, weights=weights15)
    call check_close('rank15 full-rank gmean first', &
        r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 3.5266071794790452_dp, tol)
    call check_close('rank15 full-rank gmean second', &
        r14(2,1,1,1,1,1,1,1,1,1,1,1,1,1), 4.364494543886886_dp, tol)
    r14 = hmean(x15, axis=2, weights=weights15)
    call check_close('rank15 full-rank hmean first', &
        r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 3.0_dp, tol)
    r14 = pmean(x15, 2.0_dp, axis=2, weights=weights15)
    call check_close('rank15 full-rank pmean second', &
        r14(2,1,1,1,1,1,1,1,1,1,1,1,1,1), 4.898979485566356_dp, tol)
    k15 = gmean_keepdims(x15, axis=2, weights=weights15)
    call check_true('rank15 weighted gmean keepdims shape', &
        all(shape(k15) == [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
    call check_close('rank15 weighted gmean keepdims value', &
        k15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), 3.5266071794790452_dp, tol)
    k15 = hmean_keepdims(x15, axis=2, weights=weights15)
    call check_close('rank15 weighted hmean keepdims', &
        k15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), 3.0_dp, tol)
    k15 = pmean_keepdims(x15, 2.0_dp, axis=2, weights=weights15)
    call check_close('rank15 weighted pmean keepdims', &
        k15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1), 4.898979485566356_dp, tol)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    nan_weights15 = weights15
    nan_weights15(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = nan
    r14 = gmean(x15, axis=2, weights=nan_weights15, nan_policy='omit')
    call check_close('rank15 NaN weight omit', &
        r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 3.823622456658650_dp, tol)
    r14 = hmean(x15, axis=2, weights=nan_weights15, nan_policy='omit')
    call check_close('rank15 NaN weight hmean omit', &
        r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 3.0_dp, tol)
    r14 = pmean(x15, 2.0_dp, axis=2, weights=nan_weights15, nan_policy='omit')
    call check_close('rank15 NaN weight pmean omit', &
        r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), 4.582575694955840_dp, tol)
    r14 = gmean(x15, axis=2, weights=nan_weights15, nan_policy='propagate')
    call check_true('rank15 NaN weight propagate', &
        ieee_is_nan(r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1)))

    r14 = gmean(empty15, axis=1, weights=scalar_weight)
    call check_true('rank15 weighted empty reduced shape', &
        all(shape(r14) == [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
    call check_true('rank15 weighted empty reduced NaN', all(ieee_is_nan(r14)))

    call check_true('rank15 weighted invalid axis', &
        size(gmean(x15, axis=16, weights=scalar_weight)) == 0)

    print *, 'descriptive high-rank weighted mean tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance

        if (ieee_is_nan(actual)) then
            print *, 'FAIL ', name
            print *, ' actual is NaN'
            error stop 1
        end if
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_descriptive_highrank_weights
