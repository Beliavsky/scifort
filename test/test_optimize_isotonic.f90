! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Reference output arrays from scipy.optimize.isotonic_regression in SciPy 1.17.
program test_optimize_isotonic
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_isotonic, only : isotonic_regression, isotonic_result, &
        isotonic_status_success, isotonic_status_invalid_input, &
        isotonic_status_numerical_failure
    implicit none

    type(isotonic_result) :: res
    real(dp) :: huge_value
    real(dp) :: nan_value

    res = isotonic_regression([3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 4.0_dp])
    call check(res, [2.0_dp, 2.0_dp, 2.0_dp, 4.5_dp, 4.5_dp], &
        [3.0_dp, 2.0_dp], [1, 4, 6])

    res = isotonic_regression([3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 4.0_dp], &
        weights=[1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 1.0_dp])
    call check(res, [1.5_dp, 1.5_dp, 2.0_dp, 4.8_dp, 4.8_dp], &
        [4.0_dp, 2.0_dp, 5.0_dp], [1, 3, 4, 6])

    res = isotonic_regression([3.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 4.0_dp], &
        increasing=.false.)
    call check(res, [3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp], [5.0_dp], [1, 6])

    res = isotonic_regression([1.0_dp, 1.0_dp, 1.0_dp])
    call check(res, [1.0_dp, 1.0_dp, 1.0_dp], [3.0_dp], [1, 4])

    res = isotonic_regression([1.5_dp, 1.0_dp, 4.0_dp, 6.0_dp, 5.7_dp, 5.0_dp, &
        7.8_dp, 9.0_dp, 7.5_dp, 9.5_dp, 9.0_dp])
    call check(res, [1.25_dp, 1.25_dp, 4.0_dp, &
        5.566666666666667_dp, 5.566666666666667_dp, 5.566666666666667_dp, &
        7.8_dp, 8.25_dp, 8.25_dp, 9.25_dp, 9.25_dp], &
        [2.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [1, 3, 4, 7, 8, 10, 12])

    res = isotonic_regression([10.0_dp, -2.0_dp, 3.0_dp, 6.0_dp, -3.0_dp, 4.0_dp], &
        weights=[1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 5.0_dp, 2.0_dp], &
        increasing=.false.)
    call check(res, [10.0_dp, 2.833333333333333_dp, 2.833333333333333_dp, &
        2.833333333333333_dp, -1.0_dp, -1.0_dp], &
        [1.0_dp, 6.0_dp, 7.0_dp], [1, 2, 5, 7])

    res = isotonic_regression([real(dp)::])
    call check(res, [real(dp)::], [real(dp)::], [1])
    res = isotonic_regression([5.0_dp])
    call check(res, [5.0_dp], [1.0_dp], [1, 2])

    huge_value = huge(1.0_dp) / 8.0_dp
    res = isotonic_regression([huge_value, -huge_value])
    call check(res, [0.0_dp, 0.0_dp], [2.0_dp], [1, 3])

    res = isotonic_regression([2.0_dp, 1.0_dp], weights=[1.0_dp])
    call require(res%status == isotonic_status_invalid_input, 'weight length validation')
    res = isotonic_regression([2.0_dp, 1.0_dp], weights=[1.0_dp, 0.0_dp])
    call require(res%status == isotonic_status_invalid_input, 'zero weight validation')
    res = isotonic_regression([2.0_dp, 1.0_dp], weights=[1.0_dp, -1.0_dp])
    call require(res%status == isotonic_status_invalid_input, 'negative weight validation')
    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    res = isotonic_regression([2.0_dp, nan_value])
    call require(res%status == isotonic_status_invalid_input, 'NaN value validation')
    res = isotonic_regression([2.0_dp, 1.0_dp], weights=[huge_value * 8.0_dp, huge_value * 8.0_dp])
    call require(res%status == isotonic_status_numerical_failure, 'weight overflow guard')

    print '(a)', 'test_optimize_isotonic: PASS'

contains

    subroutine check(res, values, pool_weights, boundaries)
        type(isotonic_result), intent(in) :: res !! fitted result to check
        real(dp), intent(in) :: values(:) !! SciPy reference fitted values
        real(dp), intent(in) :: pool_weights(:) !! reference pooled case weights
        integer, intent(in) :: boundaries(:) !! reference one-based block starts and end

        call require(res%status == isotonic_status_success, 'fit status')
        call require(size(res%x) == size(values), 'fitted length')
        call require(size(res%weights) == size(pool_weights), 'pool length')
        call require(size(res%blocks) == size(boundaries), 'block length')
        call require(all(res%blocks == boundaries), 'block boundaries')
        call require(all(abs(res%x - values) < 8.0e-14_dp * &
            max(1.0_dp, abs(values))), 'fitted values')
        call require(all(abs(res%weights - pool_weights) < 8.0e-14_dp * &
            max(1.0_dp, abs(pool_weights))), 'pooled weights')
    end subroutine check

    subroutine require(condition, message)
        logical, intent(in) :: condition !! assertion condition
        character(*), intent(in) :: message !! description on failure
        if (.not. condition) then
            print '(a,a)', 'test_optimize_isotonic: ', message
            error stop 1
        end if
    end subroutine require

end program test_optimize_isotonic
