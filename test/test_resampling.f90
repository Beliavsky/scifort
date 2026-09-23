! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module test_resampling_statistics
    use scifort_stats, only : dp, mean
    implicit none
    private
    public :: difference_of_means
    public :: sample_mean

contains

    pure function sample_mean(x) result(value)
        real(dp), intent(in) :: x(:) !! observations in one bootstrap resample
        real(dp) :: value

        value = mean(x)
    end function sample_mean

    pure function difference_of_means(x, y) result(value)
        real(dp), intent(in) :: x(:) !! observations assigned to the first group
        real(dp), intent(in) :: y(:) !! observations assigned to the second group
        real(dp) :: value

        value = mean(x) - mean(y)
    end function difference_of_means

end module test_resampling_statistics

program test_resampling
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_random, only : rng_get_state, rng_seed, rng_state
    use scifort_stats, only : bootstrap_percentile, bootstrap_result, dp, &
        permutation_test, permutation_test_result
    use test_resampling_statistics, only : difference_of_means, sample_mean
    implicit none

    real(dp), parameter :: bootstrap_data(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    real(dp), parameter :: permutation_x(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: permutation_y(3) = [7.0_dp, 8.0_dp, 9.0_dp]

    integer(int64) :: after(6)
    integer(int64) :: before(6)
    type(bootstrap_result) :: bootstrap
    type(permutation_test_result) :: permutation
    type(rng_state) :: state
    type(rng_state) :: state_repeat

    call rng_seed(state, 1234)
    bootstrap = bootstrap_percentile(&
        bootstrap_data, sample_mean, state, 2000, 0.90_dp)
    call check_close('bootstrap statistic', bootstrap%statistic, 3.0_dp, 0.0_dp)
    call check_close('bootstrap standard error', bootstrap%standard_error, &
        0.63248435413780368_dp, 3.0e-15_dp)
    call check_close('bootstrap lower', bootstrap%confidence_low, 2.0_dp, 0.0_dp)
    call check_close('bootstrap upper', bootstrap%confidence_high, 4.0_dp, 0.0_dp)
    call check_true('bootstrap draw count', bootstrap%n_resamples == 2000)

    call rng_seed(state_repeat, 1234)
    bootstrap = bootstrap_percentile(&
        bootstrap_data, sample_mean, state_repeat, 2000, 0.90_dp)
    call rng_get_state(state, before)
    call rng_get_state(state_repeat, after)
    call check_true('bootstrap reproducible state', all(before == after))

    call rng_seed(state, 2468)
    permutation = permutation_test(&
        permutation_x, permutation_y, difference_of_means, state, 5000)
    call check_close('permutation statistic', permutation%statistic, -6.0_dp, 0.0_dp)
    call check_close('permutation two-sided pvalue', permutation%pvalue, 0.1_dp, 0.0_dp)
    call check_true('permutation exact flag', permutation%exact)
    call check_true('permutation draw count', permutation%n_resamples == 20)

    call rng_seed(state, 2468)
    permutation = permutation_test(&
        permutation_x, permutation_y, difference_of_means, state, 5000, 'less')
    call check_close('permutation less pvalue', permutation%pvalue, 0.05_dp, 0.0_dp)
    call rng_seed(state, 2468)
    permutation = permutation_test(&
        permutation_x, permutation_y, difference_of_means, state, 5000, 'greater')
    call check_close('permutation greater pvalue', permutation%pvalue, 1.0_dp, 0.0_dp)

    call rng_seed(state, 112233)
    call rng_get_state(state, before)
    bootstrap = bootstrap_percentile(&
        bootstrap_data, sample_mean, state, 1, 0.95_dp)
    call rng_get_state(state, after)
    call check_true('invalid bootstrap preserves rng', all(before == after))
    call check_true('invalid bootstrap result', ieee_is_nan(bootstrap%statistic))

    call rng_get_state(state, before)
    permutation = permutation_test(&
        permutation_x, permutation_y, difference_of_means, state, 100, 'bad')
    call rng_get_state(state, after)
    call check_true('invalid permutation preserves rng', all(before == after))
    call check_true('invalid permutation result', ieee_is_nan(permutation%pvalue))

    print *, 'resampling tests passed'

contains

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_resampling
