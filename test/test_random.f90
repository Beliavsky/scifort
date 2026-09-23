! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_random
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_random
    use scifort_stats
    implicit none

    integer :: failures

    failures = 0
    call test_reference_sequence(failures)
    call test_seed_state(failures)
    call test_uniform_properties(failures)
    call test_scalar_variates(failures)
    call test_array_variates(failures)
    call test_invalid_variates(failures)

    if (failures == 0) then
        print '(a)', 'test_random: PASS'
    else
        print '(a,1x,i0)', 'test_random: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_reference_sequence(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        integer(int64) :: values(6)
        real(dp) :: actual
        type(rng_state) :: state
        real(dp), parameter :: expected(12) = [ &
            0.12902716284369731_dp, &
            0.31409372149341552_dp, &
            0.22514783653411208_dp, &
            0.48840551616554950_dp, &
            0.13814694833531760_dp, &
            0.58469108367820122_dp, &
            0.33150948947956504_dp, &
            0.61974650500018813_dp, &
            0.30372054926599013_dp, &
            0.98176509359140851_dp, &
            0.76550250030179801_dp, &
            0.08176988692848719_dp ]

        do i = 1, size(expected)
            actual = rng_uniform(state)
            call check_close('rng reference sequence', actual, expected(i), 2.0e-16_dp, failures)
        end do
        call rng_get_state(state, values)
        call check_true('rng reference state', all(values == [ &
            1503475344_int64, &
            4076668063_int64, &
            4132706289_int64, &
            723801936_int64, &
            3730956535_int64, &
            3197591836_int64 ]), failures)
    end subroutine test_reference_sequence

    subroutine test_seed_state(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        integer :: status
        integer(int64) :: before(6)
        integer(int64) :: invalid(6)
        integer(int64) :: saved(6)
        real(dp) :: a(16)
        real(dp) :: b(16)
        real(dp) :: x1
        real(dp) :: x2
        type(rng_state) :: state1
        type(rng_state) :: state2

        call rng_seed(state1, 12345)
        call rng_seed(state2, 12345_int64)
        call rng_fill_uniform(state1, a)
        do i = 1, size(b)
            b(i) = rng_uniform(state2)
        end do
        call check_true('rng scalar/int64 seed parity', all(a == b), failures)

        state2 = rng_state()
        do i = 1, size(b)
            b(i) = rng_uniform(state2)
        end do
        call check_true('rng explicit/default seed parity', all(a == b), failures)

        call rng_get_state(state1, saved)
        x1 = rng_uniform(state1)
        call rng_seed(state2, 1)
        call rng_set_state(state2, saved, status)
        call check_true('rng set-state status', status == rng_status_ok, failures)
        x2 = rng_uniform(state2)
        call check_close('rng state restore', x1, x2, 0.0_dp, failures)

        call rng_get_state(state2, before)
        invalid = 0_int64
        call rng_set_state(state2, invalid, status)
        call check_true('rng invalid-state status', status == rng_status_invalid_state, failures)
        call rng_get_state(state2, saved)
        call check_true('rng invalid state unchanged', all(saved == before), failures)
    end subroutine test_seed_state

    subroutine test_uniform_properties(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        real(dp) :: mean
        real(dp) :: second
        real(dp) :: u
        type(rng_state) :: state

        call rng_seed(state, 97531)
        mean = 0.0_dp
        second = 0.0_dp
        do i = 1, 100000
            u = rng_uniform(state)
            call check_true('rng open interval', u > 0.0_dp .and. u < 1.0_dp, failures)
            mean = mean + u
            second = second + u * u
        end do
        mean = mean / 100000.0_dp
        second = second / 100000.0_dp - mean * mean
        call check_close('rng mean smoke', mean, 0.5_dp, 5.0e-3_dp, failures)
        call check_close('rng variance smoke', second, 1.0_dp / 12.0_dp, 3.0e-3_dp, failures)
    end subroutine test_uniform_properties

    subroutine test_scalar_variates(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: actual
        real(dp) :: expected
        type(rng_state) :: state
        type(rng_state) :: reference

        call rng_seed(state, 24680)
        call rng_seed(reference, 24680)

        actual = normal_rvs(state, -0.2_dp, 1.3_dp)
        expected = normal_ppf(rng_uniform(reference), -0.2_dp, 1.3_dp)
        call check_close('rvs normal', actual, expected, 0.0_dp, failures)
        actual = uniform_rvs(state, 0.5_dp, 2.0_dp)
        expected = uniform_ppf(rng_uniform(reference), 0.5_dp, 2.0_dp)
        call check_close('rvs uniform', actual, expected, 0.0_dp, failures)
        actual = exponential_rvs(state, 0.1_dp, 2.3_dp)
        expected = exponential_ppf(rng_uniform(reference), 0.1_dp, 2.3_dp)
        call check_close('rvs exponential', actual, expected, 0.0_dp, failures)
        actual = laplace_rvs(state, -0.3_dp, 1.7_dp)
        expected = laplace_ppf(rng_uniform(reference), -0.3_dp, 1.7_dp)
        call check_close('rvs laplace', actual, expected, 0.0_dp, failures)
        actual = logistic_rvs(state, -0.3_dp, 1.7_dp)
        expected = logistic_ppf(rng_uniform(reference), -0.3_dp, 1.7_dp)
        call check_close('rvs logistic', actual, expected, 0.0_dp, failures)
        actual = cauchy_rvs(state, -0.3_dp, 1.7_dp)
        expected = cauchy_ppf(rng_uniform(reference), -0.3_dp, 1.7_dp)
        call check_close('rvs cauchy', actual, expected, 0.0_dp, failures)
        actual = rayleigh_rvs(state, 0.2_dp, 0.9_dp)
        expected = rayleigh_ppf(rng_uniform(reference), 0.2_dp, 0.9_dp)
        call check_close('rvs rayleigh', actual, expected, 0.0_dp, failures)
        actual = gamma_rvs(state, 2.5_dp, 0.2_dp, 1.3_dp)
        expected = gamma_ppf(rng_uniform(reference), 2.5_dp, 0.2_dp, 1.3_dp)
        call check_close('rvs gamma', actual, expected, 0.0_dp, failures)
        actual = chi2_rvs(state, 3.5_dp, 0.1_dp, 1.2_dp)
        expected = chi2_ppf(rng_uniform(reference), 3.5_dp, 0.1_dp, 1.2_dp)
        call check_close('rvs chi2', actual, expected, 0.0_dp, failures)
        actual = t_rvs(state, 7.0_dp, -0.1_dp, 1.4_dp)
        expected = t_ppf(rng_uniform(reference), 7.0_dp, -0.1_dp, 1.4_dp)
        call check_close('rvs t', actual, expected, 0.0_dp, failures)
        actual = lognormal_rvs(state, 0.8_dp, 0.2_dp, 1.1_dp)
        expected = lognormal_ppf(rng_uniform(reference), 0.8_dp, 0.2_dp, 1.1_dp)
        call check_close('rvs lognormal', actual, expected, 0.0_dp, failures)
        actual = weibull_rvs(state, 1.7_dp, 0.2_dp, 0.9_dp)
        expected = weibull_ppf(rng_uniform(reference), 1.7_dp, 0.2_dp, 0.9_dp)
        call check_close('rvs weibull', actual, expected, 0.0_dp, failures)
        actual = pareto_rvs(state, 2.2_dp, 0.1_dp, 1.1_dp)
        expected = pareto_ppf(rng_uniform(reference), 2.2_dp, 0.1_dp, 1.1_dp)
        call check_close('rvs pareto', actual, expected, 0.0_dp, failures)
        actual = beta_rvs(state, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp)
        expected = beta_ppf(rng_uniform(reference), 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp)
        call check_close('rvs beta', actual, expected, 0.0_dp, failures)
        actual = f_rvs(state, 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp)
        expected = f_ppf(rng_uniform(reference), 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp)
        call check_close('rvs f', actual, expected, 0.0_dp, failures)
        actual = bernoulli_rvs(state, 0.3_dp, 2.0_dp)
        expected = bernoulli_ppf(rng_uniform(reference), 0.3_dp, 2.0_dp)
        call check_close('rvs bernoulli', actual, expected, 0.0_dp, failures)
        actual = poisson_rvs(state, 2.5_dp, 1.0_dp)
        expected = poisson_ppf(rng_uniform(reference), 2.5_dp, 1.0_dp)
        call check_close('rvs poisson', actual, expected, 0.0_dp, failures)
        actual = geometric_rvs(state, 0.3_dp, 2.0_dp)
        expected = geometric_ppf(rng_uniform(reference), 0.3_dp, 2.0_dp)
        call check_close('rvs geometric', actual, expected, 0.0_dp, failures)
        actual = binomial_rvs(state, 7.0_dp, 0.4_dp, 1.0_dp)
        expected = binomial_ppf(rng_uniform(reference), 7.0_dp, 0.4_dp, 1.0_dp)
        call check_close('rvs binomial', actual, expected, 0.0_dp, failures)
        actual = negative_binomial_rvs(state, 2.5_dp, 0.4_dp, 1.0_dp)
        expected = negative_binomial_ppf(rng_uniform(reference), 2.5_dp, 0.4_dp, 1.0_dp)
        call check_close('rvs negative_binomial', actual, expected, 0.0_dp, failures)
    end subroutine test_scalar_variates

    subroutine test_array_variates(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: i
        real(dp) :: actual(3)
        real(dp) :: expected(3)
        type(rng_state) :: state
        type(rng_state) :: reference

        call normal_rvs_array(state, actual, -0.2_dp, 1.3_dp)
        do i = 1, size(expected)
            expected(i) = normal_rvs(reference, -0.2_dp, 1.3_dp)
        end do
        call check_true('rvs array normal', all(actual == expected), failures)
        call uniform_rvs_array(state, actual, 0.5_dp, 2.0_dp)
        do i = 1, size(expected)
            expected(i) = uniform_rvs(reference, 0.5_dp, 2.0_dp)
        end do
        call check_true('rvs array uniform', all(actual == expected), failures)
        call exponential_rvs_array(state, actual, 0.1_dp, 2.3_dp)
        do i = 1, size(expected)
            expected(i) = exponential_rvs(reference, 0.1_dp, 2.3_dp)
        end do
        call check_true('rvs array exponential', all(actual == expected), failures)
        call laplace_rvs_array(state, actual, -0.3_dp, 1.7_dp)
        do i = 1, size(expected)
            expected(i) = laplace_rvs(reference, -0.3_dp, 1.7_dp)
        end do
        call check_true('rvs array laplace', all(actual == expected), failures)
        call logistic_rvs_array(state, actual, -0.3_dp, 1.7_dp)
        do i = 1, size(expected)
            expected(i) = logistic_rvs(reference, -0.3_dp, 1.7_dp)
        end do
        call check_true('rvs array logistic', all(actual == expected), failures)
        call cauchy_rvs_array(state, actual, -0.3_dp, 1.7_dp)
        do i = 1, size(expected)
            expected(i) = cauchy_rvs(reference, -0.3_dp, 1.7_dp)
        end do
        call check_true('rvs array cauchy', all(actual == expected), failures)
        call rayleigh_rvs_array(state, actual, 0.2_dp, 0.9_dp)
        do i = 1, size(expected)
            expected(i) = rayleigh_rvs(reference, 0.2_dp, 0.9_dp)
        end do
        call check_true('rvs array rayleigh', all(actual == expected), failures)
        call gamma_rvs_array(state, actual, 2.5_dp, 0.2_dp, 1.3_dp)
        do i = 1, size(expected)
            expected(i) = gamma_rvs(reference, 2.5_dp, 0.2_dp, 1.3_dp)
        end do
        call check_true('rvs array gamma', all(actual == expected), failures)
        call chi2_rvs_array(state, actual, 3.5_dp, 0.1_dp, 1.2_dp)
        do i = 1, size(expected)
            expected(i) = chi2_rvs(reference, 3.5_dp, 0.1_dp, 1.2_dp)
        end do
        call check_true('rvs array chi2', all(actual == expected), failures)
        call t_rvs_array(state, actual, 7.0_dp, -0.1_dp, 1.4_dp)
        do i = 1, size(expected)
            expected(i) = t_rvs(reference, 7.0_dp, -0.1_dp, 1.4_dp)
        end do
        call check_true('rvs array t', all(actual == expected), failures)
        call lognormal_rvs_array(state, actual, 0.8_dp, 0.2_dp, 1.1_dp)
        do i = 1, size(expected)
            expected(i) = lognormal_rvs(reference, 0.8_dp, 0.2_dp, 1.1_dp)
        end do
        call check_true('rvs array lognormal', all(actual == expected), failures)
        call weibull_rvs_array(state, actual, 1.7_dp, 0.2_dp, 0.9_dp)
        do i = 1, size(expected)
            expected(i) = weibull_rvs(reference, 1.7_dp, 0.2_dp, 0.9_dp)
        end do
        call check_true('rvs array weibull', all(actual == expected), failures)
        call pareto_rvs_array(state, actual, 2.2_dp, 0.1_dp, 1.1_dp)
        do i = 1, size(expected)
            expected(i) = pareto_rvs(reference, 2.2_dp, 0.1_dp, 1.1_dp)
        end do
        call check_true('rvs array pareto', all(actual == expected), failures)
        call beta_rvs_array(state, actual, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp)
        do i = 1, size(expected)
            expected(i) = beta_rvs(reference, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp)
        end do
        call check_true('rvs array beta', all(actual == expected), failures)
        call f_rvs_array(state, actual, 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp)
        do i = 1, size(expected)
            expected(i) = f_rvs(reference, 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp)
        end do
        call check_true('rvs array f', all(actual == expected), failures)
        call bernoulli_rvs_array(state, actual, 0.3_dp, 2.0_dp)
        do i = 1, size(expected)
            expected(i) = bernoulli_rvs(reference, 0.3_dp, 2.0_dp)
        end do
        call check_true('rvs array bernoulli', all(actual == expected), failures)
        call poisson_rvs_array(state, actual, 2.5_dp, 1.0_dp)
        do i = 1, size(expected)
            expected(i) = poisson_rvs(reference, 2.5_dp, 1.0_dp)
        end do
        call check_true('rvs array poisson', all(actual == expected), failures)
        call geometric_rvs_array(state, actual, 0.3_dp, 2.0_dp)
        do i = 1, size(expected)
            expected(i) = geometric_rvs(reference, 0.3_dp, 2.0_dp)
        end do
        call check_true('rvs array geometric', all(actual == expected), failures)
        call binomial_rvs_array(state, actual, 7.0_dp, 0.4_dp, 1.0_dp)
        do i = 1, size(expected)
            expected(i) = binomial_rvs(reference, 7.0_dp, 0.4_dp, 1.0_dp)
        end do
        call check_true('rvs array binomial', all(actual == expected), failures)
        call negative_binomial_rvs_array(state, actual, 2.5_dp, 0.4_dp, 1.0_dp)
        do i = 1, size(expected)
            expected(i) = negative_binomial_rvs(reference, 2.5_dp, 0.4_dp, 1.0_dp)
        end do
        call check_true('rvs array negative_binomial', all(actual == expected), failures)
    end subroutine test_array_variates

    subroutine test_invalid_variates(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: samples(4)
        real(dp) :: u1
        real(dp) :: u2
        real(dp) :: x
        type(rng_state) :: state
        type(rng_state) :: reference

        call rng_seed(state, 86420)
        call rng_seed(reference, 86420)
        x = gamma_rvs(state, -1.0_dp)
        call check_true('invalid scalar rvs NaN', ieee_is_nan(x), failures)
        u1 = rng_uniform(state)
        u2 = rng_uniform(reference)
        call check_close('invalid scalar rvs state', u1, u2, 0.0_dp, failures)

        call rng_seed(state, 13579)
        call rng_seed(reference, 13579)
        call binomial_rvs_array(state, samples, 4.5_dp, 0.3_dp)
        call check_true('invalid array rvs NaN', all(ieee_is_nan(samples)), failures)
        u1 = rng_uniform(state)
        u2 = rng_uniform(reference)
        call check_close('invalid array rvs state', u1, u2, 0.0_dp, failures)
    end subroutine test_invalid_variates

    subroutine check_close(name, actual, expected, atol, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: atol !! permitted absolute error
        integer, intent(inout) :: failures !! running count of failed checks

        if (ieee_is_nan(actual) .or. abs(actual - expected) > atol) then
            print '(a,2(1x,es24.16))', trim(name), actual, expected
            failures = failures + 1
        end if
    end subroutine check_close

    subroutine check_true(name, condition, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        logical, intent(in) :: condition !! condition that must be true
        integer, intent(inout) :: failures !! running count of failed checks

        if (.not. condition) then
            print '(a)', trim(name)
            failures = failures + 1
        end if
    end subroutine check_true

end program test_random
