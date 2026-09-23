! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_likelihood
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_stats
    implicit none

    integer :: failures

    failures = 0
    call test_continuous(failures)
    call test_discrete(failures)
    call test_edges(failures)

    if (failures == 0) then
        print '(a)', 'test_likelihood: PASS'
    else
        print '(a,1x,i0)', 'test_likelihood: FAIL', failures
        error stop 1
    end if

contains

    subroutine test_continuous(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: x(3)
        real(dp) :: expected

        x = [-1.0_dp, 0.2_dp, 1.4_dp]
        expected = sum(normal_logpdf(x, -0.2_dp, 1.3_dp))
        call check_equal('normal loglike', &
            normal_loglikelihood(x, -0.2_dp, 1.3_dp), expected, failures)
        call check_equal('normal nnlf', normal_nnlf(x, -0.2_dp, 1.3_dp), -expected, failures)

        x = [0.2_dp, 0.5_dp, 0.8_dp]
        expected = sum(uniform_logpdf(x, 0.1_dp, 1.2_dp))
        call check_equal('uniform loglike', &
            uniform_loglikelihood(x, 0.1_dp, 1.2_dp), expected, failures)
        call check_equal('exponential loglike', exponential_loglikelihood(x, 0.1_dp, 1.2_dp), &
            sum(exponential_logpdf(x, 0.1_dp, 1.2_dp)), failures)
        call check_equal('laplace loglike', laplace_loglikelihood(x, 0.1_dp, 1.2_dp), &
            sum(laplace_logpdf(x, 0.1_dp, 1.2_dp)), failures)
        call check_equal('logistic loglike', logistic_loglikelihood(x, 0.1_dp, 1.2_dp), &
            sum(logistic_logpdf(x, 0.1_dp, 1.2_dp)), failures)
        call check_equal('cauchy loglike', cauchy_loglikelihood(x, 0.1_dp, 1.2_dp), &
            sum(cauchy_logpdf(x, 0.1_dp, 1.2_dp)), failures)
        call check_equal('rayleigh loglike', rayleigh_loglikelihood(x, 0.1_dp, 1.2_dp), &
            sum(rayleigh_logpdf(x, 0.1_dp, 1.2_dp)), failures)
        call check_equal('gamma loglike', gamma_loglikelihood(x, 2.5_dp, 0.1_dp, 1.2_dp), &
            sum(gamma_logpdf(x, 2.5_dp, 0.1_dp, 1.2_dp)), failures)
        call check_equal('chi2 loglike', chi2_loglikelihood(x, 3.5_dp, 0.1_dp, 1.2_dp), &
            sum(chi2_logpdf(x, 3.5_dp, 0.1_dp, 1.2_dp)), failures)

        x = [-1.0_dp, 0.2_dp, 1.4_dp]
        call check_equal('t loglike', t_loglikelihood(x, 7.0_dp, -0.1_dp, 1.4_dp), &
            sum(t_logpdf(x, 7.0_dp, -0.1_dp, 1.4_dp)), failures)

        x = [0.3_dp, 1.0_dp, 2.0_dp]
        call check_equal('lognormal loglike', lognormal_loglikelihood(x, 0.8_dp, 0.1_dp, 1.1_dp), &
            sum(lognormal_logpdf(x, 0.8_dp, 0.1_dp, 1.1_dp)), failures)
        call check_equal('weibull loglike', weibull_loglikelihood(x, 1.7_dp, 0.1_dp, 1.1_dp), &
            sum(weibull_logpdf(x, 1.7_dp, 0.1_dp, 1.1_dp)), failures)

        x = [1.3_dp, 2.0_dp, 4.0_dp]
        call check_equal('pareto loglike', pareto_loglikelihood(x, 2.2_dp, 0.1_dp, 1.1_dp), &
            sum(pareto_logpdf(x, 2.2_dp, 0.1_dp, 1.1_dp)), failures)

        x = [0.2_dp, 0.5_dp, 0.8_dp]
        call check_equal('beta loglike', beta_loglikelihood(x, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp), &
            sum(beta_logpdf(x, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp)), failures)
        call check_equal('f loglike', f_loglikelihood(x, 5.0_dp, 9.0_dp, 0.1_dp, 1.1_dp), &
            sum(f_logpdf(x, 5.0_dp, 9.0_dp, 0.1_dp, 1.1_dp)), failures)

        expected = gamma_loglikelihood(x, 2.5_dp, 0.1_dp, 1.2_dp)
        call check_equal('gamma nnlf', gamma_nnlf(x, 2.5_dp, 0.1_dp, 1.2_dp), -expected, failures)
        expected = negative_binomial_loglikelihood([0, 1, 4], 2.5_dp, 0.4_dp)
        call check_equal('negative binomial nnlf', &
            negative_binomial_nnlf([0, 1, 4], 2.5_dp, 0.4_dp), &
            -expected, failures)
    end subroutine test_continuous

    subroutine test_discrete(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        integer :: k(3)
        real(dp) :: kr(3)
        real(dp) :: expected

        k = [2, 3, 3]
        kr = real(k, dp)
        expected = sum(bernoulli_logpmf(k, 0.3_dp, 2.0_dp))
        call check_equal('bernoulli integer loglike', &
            bernoulli_loglikelihood(k, 0.3_dp, 2.0_dp), expected, failures)
        call check_equal('bernoulli real loglike', bernoulli_loglikelihood(kr, 0.3_dp, 2.0_dp), &
            sum(bernoulli_logpmf(kr, 0.3_dp, 2.0_dp)), failures)

        k = [1, 3, 5]
        kr = real(k, dp)
        call check_equal('poisson integer loglike', poisson_loglikelihood(k, 2.5_dp, 1.0_dp), &
            sum(poisson_logpmf(k, 2.5_dp, 1.0_dp)), failures)
        call check_equal('poisson real loglike', poisson_loglikelihood(kr, 2.5_dp, 1.0_dp), &
            sum(poisson_logpmf(kr, 2.5_dp, 1.0_dp)), failures)

        k = [3, 4, 7]
        kr = real(k, dp)
        call check_equal('geometric integer loglike', geometric_loglikelihood(k, 0.3_dp, 2.0_dp), &
            sum(geometric_logpmf(k, 0.3_dp, 2.0_dp)), failures)
        call check_equal('geometric real loglike', geometric_loglikelihood(kr, 0.3_dp, 2.0_dp), &
            sum(geometric_logpmf(kr, 0.3_dp, 2.0_dp)), failures)

        k = [1, 3, 5]
        kr = real(k, dp)
        call check_equal('binomial integer loglike', binomial_loglikelihood(k, 7, 0.4_dp, 1.0_dp), &
            sum(binomial_logpmf(k, 7, 0.4_dp, 1.0_dp)), failures)
        call check_equal('binomial real loglike', &
            binomial_loglikelihood(kr, 7.0_dp, 0.4_dp, 1.0_dp), &
            sum(binomial_logpmf(kr, 7.0_dp, 0.4_dp, 1.0_dp)), failures)

        k = [1, 2, 5]
        kr = real(k, dp)
        call check_equal('negative binomial integer loglike', &
            negative_binomial_loglikelihood(k, 2.5_dp, 0.4_dp, 1.0_dp), &
            sum(negative_binomial_logpmf(k, 2.5_dp, 0.4_dp, 1.0_dp)), failures)
        call check_equal('negative binomial real loglike', &
            negative_binomial_loglikelihood(kr, 2.5_dp, 0.4_dp, 1.0_dp), &
            sum(negative_binomial_logpmf(kr, 2.5_dp, 0.4_dp, 1.0_dp)), failures)
    end subroutine test_discrete

    subroutine test_edges(failures)
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: empty(0)
        real(dp) :: x(2)

        call check_equal('empty loglike', normal_loglikelihood(empty), 0.0_dp, failures)
        call check_equal('empty nnlf', normal_nnlf(empty), -0.0_dp, failures)

        x = [0.0_dp, 1.0_dp]
        call check_true('invalid likelihood NaN', &
            ieee_is_nan(normal_loglikelihood(x, 0.0_dp, -1.0_dp)), failures)
        call check_true('invalid nnlf NaN', ieee_is_nan(gamma_nnlf(x, -1.0_dp)), failures)
    end subroutine test_edges

    subroutine check_equal(name, actual, expected, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        integer, intent(inout) :: failures !! running count of failed checks

        if (ieee_is_nan(actual) .neqv. ieee_is_nan(expected)) then
            print '(a,2(1x,es24.16))', trim(name), actual, expected
            failures = failures + 1
        else if (.not. ieee_is_nan(actual) .and. actual /= expected) then
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

end program test_likelihood
