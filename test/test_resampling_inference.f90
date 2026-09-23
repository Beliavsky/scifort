! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module test_resampling_inference_callbacks
    use scifort_random, only : rng_state
    use scifort_stats, only : dp, mean
    implicit none
    private

    public :: difference_of_means
    public :: dot_product_statistic
    public :: fixed_null_one
    public :: fixed_null_two_x
    public :: fixed_null_two_y
    public :: fixed_power_one
    public :: fixed_power_two_x
    public :: fixed_power_two_y
    public :: mean_statistic
    public :: one_sample_test_pvalue
    public :: two_sample_test_pvalue

contains

    pure function mean_statistic(sample) result(value)
        real(dp), intent(in) :: sample(:) !! observations whose mean is returned
        real(dp) :: value

        value = mean(sample)
    end function mean_statistic

    pure function difference_of_means(sample_x, sample_y) result(value)
        real(dp), intent(in) :: sample_x(:) !! first sample
        real(dp), intent(in) :: sample_y(:) !! second sample
        real(dp) :: value

        value = mean(sample_x) - mean(sample_y)
    end function difference_of_means

    pure function dot_product_statistic(sample_x, sample_y) result(value)
        real(dp), intent(in) :: sample_x(:) !! first paired sample
        real(dp), intent(in) :: sample_y(:) !! second paired sample
        real(dp) :: value

        value = dot_product(sample_x, sample_y)
    end function dot_product_statistic

    subroutine fixed_null_one(state, sample)
        type(rng_state), intent(inout) :: state !! explicit state, intentionally unchanged
        real(dp), intent(out) :: sample(:) !! deterministic null sample

        sample = 0.0_dp
    end subroutine fixed_null_one

    subroutine fixed_null_two_x(state, sample)
        type(rng_state), intent(inout) :: state !! explicit state, intentionally unchanged
        real(dp), intent(out) :: sample(:) !! deterministic first null sample

        sample = 0.0_dp
    end subroutine fixed_null_two_x

    subroutine fixed_null_two_y(state, sample)
        type(rng_state), intent(inout) :: state !! explicit state, intentionally unchanged
        real(dp), intent(out) :: sample(:) !! deterministic second null sample

        sample = 0.0_dp
    end subroutine fixed_null_two_y

    subroutine fixed_power_one(state, sample)
        type(rng_state), intent(inout) :: state !! explicit state, intentionally unchanged
        real(dp), intent(out) :: sample(:) !! deterministic alternative sample

        sample = 1.0_dp
    end subroutine fixed_power_one

    subroutine fixed_power_two_x(state, sample)
        type(rng_state), intent(inout) :: state !! explicit state, intentionally unchanged
        real(dp), intent(out) :: sample(:) !! deterministic first alternative sample

        sample = 1.0_dp
    end subroutine fixed_power_two_x

    subroutine fixed_power_two_y(state, sample)
        type(rng_state), intent(inout) :: state !! explicit state, intentionally unchanged
        real(dp), intent(out) :: sample(:) !! deterministic second alternative sample

        sample = 0.0_dp
    end subroutine fixed_power_two_y

    pure function one_sample_test_pvalue(sample) result(value)
        real(dp), intent(in) :: sample(:) !! simulated sample
        real(dp) :: value

        if (mean(sample) > 0.5_dp) then
            value = 0.01_dp
        else
            value = 0.50_dp
        end if
    end function one_sample_test_pvalue

    pure function two_sample_test_pvalue(sample_x, sample_y) result(value)
        real(dp), intent(in) :: sample_x(:) !! simulated first sample
        real(dp), intent(in) :: sample_y(:) !! simulated second sample
        real(dp) :: value

        if (mean(sample_x) > mean(sample_y)) then
            value = 0.01_dp
        else
            value = 0.50_dp
        end if
    end function two_sample_test_pvalue

end module test_resampling_inference_callbacks

program test_resampling_inference
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use resampling_inference_reference, only : ref_independent_pvalue, &
        ref_monte_greater, ref_pairings_pvalue, ref_power_one, ref_power_two, &
        ref_samples_pvalue
    use scifort_random, only : rng_get_state, rng_seed, rng_state
    use scifort_stats, only : bootstrap, bootstrap_result, dp, monte_carlo_test, &
        monte_carlo_test_result, permutation_test, permutation_test_result, power, power_result
    use test_resampling_inference_callbacks, only : difference_of_means, &
        dot_product_statistic, fixed_null_one, fixed_null_two_x, fixed_null_two_y, &
        fixed_power_one, fixed_power_two_x, fixed_power_two_y, mean_statistic, &
        one_sample_test_pvalue, two_sample_test_pvalue
    implicit none

    integer(int64) :: after(6)
    integer(int64) :: before(6)
    real(dp), parameter :: bootstrap_x(6) = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 13.0_dp]
    real(dp), parameter :: exact_x(2) = [1.0_dp, 2.0_dp]
    real(dp), parameter :: exact_y(2) = [3.0_dp, 4.0_dp]
    real(dp), parameter :: observed_one(3) = [1.0_dp, 2.0_dp, 3.0_dp]
    real(dp), parameter :: observed_two_x(2) = [2.0_dp, 3.0_dp]
    real(dp), parameter :: observed_two_y(2) = [0.0_dp, 1.0_dp]
    type(bootstrap_result) :: basic
    type(bootstrap_result) :: bca
    type(bootstrap_result) :: less_interval
    type(bootstrap_result) :: percentile
    type(monte_carlo_test_result) :: monte
    type(permutation_test_result) :: permutation
    type(power_result) :: power_value
    type(rng_state) :: state

    call rng_seed(state, 20260922)
    percentile = bootstrap(bootstrap_x, mean_statistic, state, 1000, 0.90_dp, &
        'two-sided', 'percentile')
    call check_close('bootstrap percentile statistic', percentile%statistic, &
        5.3333333333333339_dp, 2.0e-15_dp)
    call check_close('bootstrap percentile se', percentile%standard_error, &
        1.6849682840095592_dp, 2.0e-14_dp)
    call check_close('bootstrap percentile low', percentile%confidence_low, &
        2.6666666666666665_dp, 2.0e-14_dp)
    call check_close('bootstrap percentile high', percentile%confidence_high, &
        8.3333333333333321_dp, 2.0e-14_dp)
    call check_true('bootstrap percentile draws', percentile%n_resamples == 1000)
    call check_close('bootstrap first draw', percentile%bootstrap_distribution(1), &
        4.6666666666666670_dp, 2.0e-15_dp)

    call rng_seed(state, 20260922)
    basic = bootstrap(bootstrap_x, mean_statistic, state, 1000, 0.90_dp, &
        'two-sided', 'basic')
    call check_close('bootstrap basic low', basic%confidence_low, &
        2.3333333333333357_dp, 2.0e-14_dp)
    call check_close('bootstrap basic high', basic%confidence_high, &
        8.0000000000000018_dp, 2.0e-14_dp)

    call rng_seed(state, 20260922)
    bca = bootstrap(bootstrap_x, mean_statistic, state, 1000, 0.90_dp, &
        'two-sided', 'bca')
    call check_close('bootstrap bca low', bca%confidence_low, 3.0_dp, 2.0e-14_dp)
    call check_close('bootstrap bca high', bca%confidence_high, &
        8.8333333333333321_dp, 2.0e-14_dp)

    call rng_seed(state, 20260922)
    less_interval = bootstrap(bootstrap_x, mean_statistic, state, 1000, 0.95_dp, &
        'less', 'percentile')
    call check_true('bootstrap less lower infinite', less_interval%confidence_low < -huge(1.0_dp))
    call check_close('bootstrap less matches two-sided upper', less_interval%confidence_high, &
        percentile%confidence_high, 2.0e-14_dp)

    call rng_seed(state, 1234)
    call rng_get_state(state, before)
    permutation = permutation_test(exact_x, exact_y, difference_of_means, state, 100, &
        'two-sided', 'independent')
    call rng_get_state(state, after)
    call check_true('independent exact flag', permutation%exact)
    call check_true('independent exact count', permutation%n_resamples == 6)
    call check_close('independent exact pvalue', permutation%pvalue, ref_independent_pvalue, 2.0e-16_dp)
    call check_true('exact permutation preserves rng', all(before == after))

    permutation = permutation_test(exact_x, exact_y, difference_of_means, state, 100, &
        'two-sided', 'samples')
    call check_true('samples exact flag', permutation%exact)
    call check_true('samples exact count', permutation%n_resamples == 4)
    call check_close('samples exact pvalue', permutation%pvalue, ref_samples_pvalue, 0.0_dp)

    permutation = permutation_test(exact_x, exact_y, dot_product_statistic, state, 100, &
        'two-sided', 'pairings')
    call check_true('pairings exact flag', permutation%exact)
    call check_true('pairings exact count', permutation%n_resamples == 4)
    call check_close('pairings exact pvalue', permutation%pvalue, ref_pairings_pvalue, 0.0_dp)

    call rng_seed(state, 4567)
    monte = monte_carlo_test(observed_one, fixed_null_one, mean_statistic, state, 9, 'greater')
    call check_close('monte one statistic', monte%statistic, 2.0_dp, 0.0_dp)
    call check_close('monte one greater', monte%pvalue, ref_monte_greater, 2.0e-16_dp)
    call check_true('monte one count', monte%n_resamples == 9)
    call check_true('monte one null distribution', all(monte%null_distribution == 0.0_dp))

    monte = monte_carlo_test(observed_one, fixed_null_one, mean_statistic, state, 9, 'two-sided')
    call check_close('monte one two-sided', monte%pvalue, 0.2_dp, 0.0_dp)

    monte = monte_carlo_test(observed_two_x, observed_two_y, fixed_null_two_x, &
        fixed_null_two_y, difference_of_means, state, 9, 'greater')
    call check_close('monte two statistic', monte%statistic, 2.0_dp, 0.0_dp)
    call check_close('monte two greater', monte%pvalue, 0.1_dp, 0.0_dp)

    call rng_seed(state, 7890)
    power_value = power(fixed_power_one, one_sample_test_pvalue, 8, state, 0.05_dp, 25)
    call check_close('power one', power_value%power, ref_power_one, 0.0_dp)
    call check_true('power one count', power_value%n_resamples == 25)
    call check_true('power one pvalues', all(power_value%pvalues == 0.01_dp))

    power_value = power(fixed_power_one, one_sample_test_pvalue, 8, state, 0.005_dp, 25)
    call check_close('power one strict alpha', power_value%power, 0.0_dp, 0.0_dp)

    power_value = power(fixed_power_two_x, fixed_power_two_y, two_sample_test_pvalue, &
        8, 9, state, 0.05_dp, 25)
    call check_close('power two', power_value%power, ref_power_two, 0.0_dp)
    call check_true('power two pvalues', all(power_value%pvalues == 0.01_dp))

    call rng_get_state(state, before)
    power_value = power(fixed_power_one, one_sample_test_pvalue, 0, state, 0.05_dp, 25)
    call rng_get_state(state, after)
    call check_true('invalid power preserves rng', all(before == after))
    call check_true('invalid power nan', ieee_is_nan(power_value%power))

    print *, 'resampling inference tests passed'

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

end program test_resampling_inference
