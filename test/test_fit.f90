! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_fit
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_fit
    use scifort_kinds, only : dp
    implicit none

    integer, parameter :: binary(8) = [0, 1, 1, 0, 1, 0, 1, 1]
    integer, parameter :: binom_data(8) = [2, 3, 4, 5, 3, 6, 4, 5]
    integer, parameter :: counts0(8) = [0, 1, 2, 3, 4, 2, 1, 5]
    integer, parameter :: counts1(8) = [1, 2, 3, 4, 5, 2, 1, 6]
    real(dp), parameter :: x(7) = [-1.2_dp, -0.1_dp, 0.2_dp, 0.6_dp, &
        1.4_dp, 2.0_dp, 0.8_dp]
    real(dp), parameter :: xp(7) = [0.2_dp, 0.5_dp, 0.9_dp, 1.3_dp, &
        1.8_dp, 2.6_dp, 3.7_dp]
    real(dp), parameter :: xb(7) = [0.08_dp, 0.19_dp, 0.31_dp, 0.46_dp, &
        0.62_dp, 0.77_dp, 0.91_dp]
    type(fit_result) :: fit
    integer :: failures
    real(dp) :: expected
    real(dp) :: mean_x
    real(dp) :: scale_x

    failures = 0

    mean_x = sum(x) / real(size(x), dp)
    scale_x = sqrt(sum((x - mean_x)**2) / real(size(x), dp))
    call normal_fit(x, [-3.0_dp, 0.1_dp], [3.0_dp, 4.0_dp], fit, guess=[0.0_dp, 1.0_dp])
    call check_true('normal fit success', fit%success, failures)
    call check_close('normal fit loc', fit%params(1), mean_x, 2.0e-6_dp, failures)
    call check_close('normal fit scale', fit%params(2), scale_x, 2.0e-6_dp, failures)

    expected = sum(xp) / real(size(xp), dp)
    call exponential_fit(xp, [0.0_dp, 0.1_dp], [0.0_dp, 6.0_dp], fit, &
        guess=[0.0_dp, 1.0_dp])
    call check_true('exponential fit success', fit%success, failures)
    call check_equal('exponential fixed loc', fit%params(1), 0.0_dp, failures)
    call check_close('exponential fit scale', fit%params(2), expected, 2.0e-6_dp, failures)

    expected = sum(real(binary, dp)) / real(size(binary), dp)
    call bernoulli_fit(binary, [0.01_dp, 0.0_dp], [0.99_dp, 0.0_dp], fit, &
        guess=[0.5_dp, 0.0_dp])
    call check_true('bernoulli fit success', fit%success, failures)
    call check_close('bernoulli fit p', fit%params(1), expected, 2.0e-6_dp, failures)
    call check_equal('bernoulli fixed loc', fit%params(2), 0.0_dp, failures)

    expected = sum(real(counts0, dp)) / real(size(counts0), dp)
    call poisson_fit(counts0, [0.05_dp, 0.0_dp], [8.0_dp, 0.0_dp], fit, &
        guess=[2.0_dp, 0.0_dp])
    call check_true('poisson fit success', fit%success, failures)
    call check_close('poisson fit mu', fit%params(1), expected, 2.0e-6_dp, failures)

    expected = real(size(counts1), dp) / sum(real(counts1, dp))
    call geometric_fit(counts1, [0.02_dp, 0.0_dp], [0.98_dp, 0.0_dp], fit, &
        guess=[0.4_dp, 0.0_dp])
    call check_true('geometric fit success', fit%success, failures)
    call check_close('geometric fit p', fit%params(1), expected, 2.0e-6_dp, failures)

    expected = sum(real(binom_data, dp)) / (8.0_dp * real(size(binom_data), dp))
    call binomial_fit(binom_data, [8.0_dp, 0.02_dp, 0.0_dp], &
        [8.0_dp, 0.98_dp, 0.0_dp], fit, guess=[8.0_dp, 0.5_dp, 0.0_dp])
    call check_true('binomial fit success', fit%success, failures)
    call check_equal('binomial fixed n', fit%params(1), 8.0_dp, failures)
    call check_close('binomial fit p', fit%params(2), expected, 2.0e-6_dp, failures)
    call check_equal('binomial fixed loc', fit%params(3), 0.0_dp, failures)

    call binomial_fit(binom_data, [6.0_dp, 0.02_dp, 0.0_dp], &
        [12.0_dp, 0.98_dp, 0.0_dp], fit, guess=[9.0_dp, 0.5_dp, 0.0_dp])
    call check_true('binomial free n success', fit%success, failures)
    call check_true('binomial n integral', fit%params(1) == anint(fit%params(1)), failures)

    expected = 3.0_dp / (3.0_dp + sum(real(counts0, dp)) / real(size(counts0), dp))
    call negative_binomial_fit(counts0, [3.0_dp, 0.02_dp, 0.0_dp], &
        [3.0_dp, 0.98_dp, 0.0_dp], fit, guess=[3.0_dp, 0.5_dp, 0.0_dp])
    call check_true('negative binomial fit success', fit%success, failures)
    call check_equal('negative binomial fixed n', fit%params(1), 3.0_dp, failures)
    call check_close('negative binomial fit p', fit%params(2), expected, 2.0e-6_dp, failures)

    call smoke_fixed_continuous(failures)
    call smoke_fixed_discrete(failures)

    call normal_fit(x, [-3.0_dp, -1.0_dp], [3.0_dp, -1.0_dp], fit)
    call check_true('invalid fixed scale fails', .not. fit%success, failures)
    call check_true('invalid fixed scale status', &
        fit%status == fit_status_no_finite_objective, failures)
    call check_true('invalid fixed scale nnlf nan', ieee_is_nan(fit%nllf), failures)

    if (failures /= 0) then
        print '(a,i0)', 'fit test failures: ', failures
        error stop 1
    end if
    print '(a)', 'fit tests passed'

contains


    subroutine smoke_fixed_continuous(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        type(fit_result) :: r

        call uniform_fit(xb, [0.0_dp, 1.0_dp], [0.0_dp, 1.0_dp], r)
        call check_true('uniform fixed smoke', r%success, failures)
        call laplace_fit(x, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], r)
        call check_true('laplace fixed smoke', r%success, failures)
        call logistic_fit(x, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], r)
        call check_true('logistic fixed smoke', r%success, failures)
        call cauchy_fit(x, [0.1_dp, 1.4_dp], [0.1_dp, 1.4_dp], r)
        call check_true('cauchy fixed smoke', r%success, failures)
        call rayleigh_fit(xp, [0.0_dp, 1.4_dp], [0.0_dp, 1.4_dp], r)
        call check_true('rayleigh fixed smoke', r%success, failures)
        call gamma_fit(xp, [2.0_dp, 0.0_dp, 1.3_dp], [2.0_dp, 0.0_dp, 1.3_dp], r)
        call check_true('gamma fixed smoke', r%success, failures)
        call chi2_fit(xp, [4.0_dp, 0.0_dp, 1.0_dp], [4.0_dp, 0.0_dp, 1.0_dp], r)
        call check_true('chi2 fixed smoke', r%success, failures)
        call t_fit(x, [5.0_dp, 0.0_dp, 1.2_dp], [5.0_dp, 0.0_dp, 1.2_dp], r)
        call check_true('t fixed smoke', r%success, failures)
        call lognormal_fit(xp, [0.8_dp, 0.0_dp, 1.0_dp], [0.8_dp, 0.0_dp, 1.0_dp], r)
        call check_true('lognormal fixed smoke', r%success, failures)
        call weibull_fit(xp, [1.5_dp, 0.0_dp, 1.2_dp], [1.5_dp, 0.0_dp, 1.2_dp], r)
        call check_true('weibull fixed smoke', r%success, failures)
        call pareto_fit(1.1_dp + 2.0_dp * xp, [2.0_dp, 0.0_dp, 1.0_dp], &
            [2.0_dp, 0.0_dp, 1.0_dp], r)
        call check_true('pareto fixed smoke', r%success, failures)
        call beta_fit(xb, [2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp], &
            [2.0_dp, 3.0_dp, 0.0_dp, 1.0_dp], r)
        call check_true('beta fixed smoke', r%success, failures)
        call f_fit(xp, [4.0_dp, 7.0_dp, 0.0_dp, 1.0_dp], &
            [4.0_dp, 7.0_dp, 0.0_dp, 1.0_dp], r)
        call check_true('f fixed smoke', r%success, failures)
    end subroutine smoke_fixed_continuous

    subroutine smoke_fixed_discrete(failures)
        integer, intent(inout) :: failures !! running count of failed checks
        type(fit_result) :: r

        call bernoulli_fit(real(binary, dp), [0.4_dp, 0.0_dp], [0.4_dp, 0.0_dp], r)
        call check_true('bernoulli real fixed smoke', r%success, failures)
        call poisson_fit(real(counts0, dp), [2.0_dp, 0.0_dp], [2.0_dp, 0.0_dp], r)
        call check_true('poisson real fixed smoke', r%success, failures)
        call geometric_fit(real(counts1, dp), [0.4_dp, 0.0_dp], [0.4_dp, 0.0_dp], r)
        call check_true('geometric real fixed smoke', r%success, failures)
        call binomial_fit(real(binom_data, dp), [8.0_dp, 0.4_dp, 0.0_dp], &
            [8.0_dp, 0.4_dp, 0.0_dp], r)
        call check_true('binomial real fixed smoke', r%success, failures)
        call negative_binomial_fit(real(counts0, dp), [3.0_dp, 0.5_dp, 0.0_dp], &
            [3.0_dp, 0.5_dp, 0.0_dp], r)
        call check_true('negative binomial real fixed smoke', r%success, failures)
    end subroutine smoke_fixed_discrete

    subroutine check_close(name, actual, expected, tolerance, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! relative tolerance with unit absolute floor
        integer, intent(inout) :: failures !! running count of failed checks

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print '(a,2(1x,es24.16))', trim(name), actual, expected
            failures = failures + 1
        end if
    end subroutine check_close

    subroutine check_equal(name, actual, expected, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! required exact value
        integer, intent(inout) :: failures !! running count of failed checks

        if (actual /= expected) then
            print '(a,2(1x,es24.16))', trim(name), actual, expected
            failures = failures + 1
        end if
    end subroutine check_equal

    subroutine check_true(name, condition, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        logical, intent(in) :: condition !! condition that must be true
        integer, intent(inout) :: failures !! running count of failed checks

        if (.not. condition) then
            print '(a)', trim(name)
            failures = failures + 1
        end if
    end subroutine check_true

end program test_fit
