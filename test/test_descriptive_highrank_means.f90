! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference values were generated independently with SciPy 1.17.0.

program test_descriptive_highrank_means
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_stats, only : dp, gmean, gmean_keepdims, hmean, hmean_keepdims, &
        pmean, pmean_keepdims
    implicit none

    real(dp), parameter :: tol = 5.0e-12_dp
    integer :: i
    real(dp) :: nan
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
    real(dp) :: x15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: xn15(2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp) :: empty15(0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    real(dp), allocatable :: r14(:, :, :, :, :, :, :, :, :, :, :, :, :, :)
    real(dp), allocatable :: k15(:, :, :, :, :, :, :, :, :, :, :, :, :, :, :)

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

    call check_close('rank5 gmean', sum(gmean(x5)), 2.0_dp, tol)
    call check_close('rank5 hmean', sum(hmean(x5)), 1.6_dp, tol)
    call check_close('rank5 pmean', sum(pmean(x5, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank6 gmean', sum(gmean(x6)), 2.0_dp, tol)
    call check_close('rank6 hmean', sum(hmean(x6)), 1.6_dp, tol)
    call check_close('rank6 pmean', sum(pmean(x6, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank7 gmean', sum(gmean(x7)), 2.0_dp, tol)
    call check_close('rank7 hmean', sum(hmean(x7)), 1.6_dp, tol)
    call check_close('rank7 pmean', sum(pmean(x7, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank8 gmean', sum(gmean(x8)), 2.0_dp, tol)
    call check_close('rank8 hmean', sum(hmean(x8)), 1.6_dp, tol)
    call check_close('rank8 pmean', sum(pmean(x8, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank9 gmean', sum(gmean(x9)), 2.0_dp, tol)
    call check_close('rank9 hmean', sum(hmean(x9)), 1.6_dp, tol)
    call check_close('rank9 pmean', sum(pmean(x9, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank10 gmean', sum(gmean(x10)), 2.0_dp, tol)
    call check_close('rank10 hmean', sum(hmean(x10)), 1.6_dp, tol)
    call check_close('rank10 pmean', sum(pmean(x10, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank11 gmean', sum(gmean(x11)), 2.0_dp, tol)
    call check_close('rank11 hmean', sum(hmean(x11)), 1.6_dp, tol)
    call check_close('rank11 pmean', sum(pmean(x11, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank12 gmean', sum(gmean(x12)), 2.0_dp, tol)
    call check_close('rank12 hmean', sum(hmean(x12)), 1.6_dp, tol)
    call check_close('rank12 pmean', sum(pmean(x12, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank13 gmean', sum(gmean(x13)), 2.0_dp, tol)
    call check_close('rank13 hmean', sum(hmean(x13)), 1.6_dp, tol)
    call check_close('rank13 pmean', sum(pmean(x13, 2.0_dp)), sqrt(8.5_dp), tol)
    call check_close('rank14 gmean', sum(gmean(x14)), 2.0_dp, tol)
    call check_close('rank14 hmean', sum(hmean(x14)), 1.6_dp, tol)
    call check_close('rank14 pmean', sum(pmean(x14, 2.0_dp)), sqrt(8.5_dp), tol)

    x15 = reshape([(real(i, dp), i = 1, 6)], shape(x15))
    r14 = gmean(x15, axis=2)
    call check_true('rank15 gmean axis2 shape', &
        all(shape(r14) == [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
    call check_close('rank15 gmean axis2 first', r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        2.4662120743304703_dp, tol)
    call check_close('rank15 gmean axis2 second', r14(2,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        3.634241185664279_dp, tol)

    r14 = hmean(x15, axis=2)
    call check_close('rank15 hmean axis2 first', r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        1.9565217391304350_dp, tol)
    call check_close('rank15 hmean axis2 second', r14(2,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        3.2727272727272729_dp, tol)

    r14 = pmean(x15, 2.0_dp, axis=2)
    call check_close('rank15 pmean axis2 first', r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        3.4156502553198660_dp, tol)
    call check_close('rank15 pmean axis2 second', r14(2,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        4.3204937989385740_dp, tol)

    k15 = gmean_keepdims(x15, axis=2)
    call check_true('rank15 keepdims shape', &
        all(shape(k15) == [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
    call check_close('rank15 keepdims first', k15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        2.4662120743304703_dp, tol)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    xn15 = x15
    xn15(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) = nan
    r14 = gmean(xn15, axis=2, nan_policy='omit')
    call check_close('rank15 gmean omit', r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        2.2360679774997898_dp, tol)
    r14 = hmean(xn15, axis=2, nan_policy='omit')
    call check_close('rank15 hmean omit', r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        1.6666666666666667_dp, tol)
    r14 = pmean(xn15, 2.0_dp, axis=2, nan_policy='omit')
    call check_close('rank15 pmean omit', r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1), &
        3.6055512754639891_dp, tol)
    r14 = gmean(xn15, axis=2, nan_policy='propagate')
    call check_true('rank15 gmean propagate NaN', &
        ieee_is_nan(r14(1,1,1,1,1,1,1,1,1,1,1,1,1,1)))

    r14 = gmean(empty15, axis=1)
    call check_true('rank15 empty reduced shape', &
        all(shape(r14) == [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
    call check_true('rank15 empty reduced NaN', all(ieee_is_nan(r14)))
    r14 = gmean(empty15, axis=2)
    call check_true('rank15 empty retained zero extent', size(r14) == 0)
    k15 = gmean_keepdims(empty15, axis=1)
    call check_true('rank15 empty keepdims shape', &
        all(shape(k15) == [1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
    call check_true('rank15 empty keepdims NaN', all(ieee_is_nan(k15)))

    call check_true('rank15 invalid axis', size(gmean(x15, axis=16)) == 0)

    print *, 'descriptive high-rank mean tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance

        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
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

end program test_descriptive_highrank_means
