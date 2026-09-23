! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_score
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_digamma, only : digamma_positive
    use scifort_kinds, only : dp
    use scifort_likelihood
    use scifort_score
    implicit none

    integer :: failures
    integer, parameter :: binary(6) = [0, 1, 1, 0, 1, 0]
    integer, parameter :: counts0(6) = [0, 1, 2, 3, 4, 2]
    integer, parameter :: counts1(6) = [1, 2, 3, 4, 5, 2]
    real(dp), parameter :: x(5) = [-1.4_dp, -0.2_dp, 0.7_dp, 1.9_dp, 3.2_dp]
    real(dp), parameter :: xp(5) = [0.35_dp, 0.8_dp, 1.4_dp, 2.6_dp, 4.1_dp]
    real(dp), parameter :: xb(5) = [-0.08_dp, 0.24_dp, 0.70_dp, 1.14_dp, 1.50_dp]
    real(dp), parameter :: xpareto(5) = [1.24_dp, 1.75_dp, 2.70_dp, 4.60_dp, 8.40_dp]
    real(dp), parameter :: tol = 4.0e-6_dp
    real(dp) :: h
    real(dp) :: s2(2)
    real(dp) :: s3(3)
    real(dp) :: s4(4)

    failures = 0

    call check_close('digamma 0.1', digamma_positive(0.1_dp), &
        -1.04237549404110759e1_dp, 3.0e-15_dp, failures)
    call check_close('digamma 1', digamma_positive(1.0_dp), &
        -5.77215664901532866e-1_dp, 3.0e-15_dp, failures)
    call check_close('digamma 2.5', digamma_positive(2.5_dp), &
        7.03156640645243192e-1_dp, 3.0e-15_dp, failures)
    call check_close('digamma 10', digamma_positive(10.0_dp), &
        2.25175258906672093_dp, 3.0e-15_dp, failures)
    call check_close('digamma 100', digamma_positive(100.0_dp), &
        4.60016185273808809_dp, 3.0e-15_dp, failures)
    call check_true('digamma invalid', ieee_is_nan(digamma_positive(0.0_dp)), failures)

    s2 = normal_score(x, 0.3_dp, 1.7_dp)
    h = 1.0e-6_dp
    call check_close('normal loc', s2(1), &
        centered(normal_loglikelihood(x, 0.3_dp + h, 1.7_dp), &
        normal_loglikelihood(x, 0.3_dp - h, 1.7_dp), h), tol, failures)
    call check_close('normal scale', s2(2), &
        centered(normal_loglikelihood(x, 0.3_dp, 1.7_dp + h), &
        normal_loglikelihood(x, 0.3_dp, 1.7_dp - h), h), tol, failures)

    s2 = uniform_score(xb, -0.3_dp, 2.0_dp)
    call check_close('uniform loc', s2(1), &
        centered(uniform_loglikelihood(xb, -0.3_dp + h, 2.0_dp), &
        uniform_loglikelihood(xb, -0.3_dp - h, 2.0_dp), h), tol, failures)
    call check_close('uniform scale', s2(2), &
        centered(uniform_loglikelihood(xb, -0.3_dp, 2.0_dp + h), &
        uniform_loglikelihood(xb, -0.3_dp, 2.0_dp - h), h), tol, failures)

    s2 = exponential_score(xp, -0.2_dp, 1.4_dp)
    call check_close('exponential loc', s2(1), &
        centered(exponential_loglikelihood(xp, -0.2_dp + h, 1.4_dp), &
        exponential_loglikelihood(xp, -0.2_dp - h, 1.4_dp), h), tol, failures)
    call check_close('exponential scale', s2(2), &
        centered(exponential_loglikelihood(xp, -0.2_dp, 1.4_dp + h), &
        exponential_loglikelihood(xp, -0.2_dp, 1.4_dp - h), h), tol, failures)

    s2 = laplace_score(x, 0.15_dp, 1.3_dp)
    call check_close('laplace loc', s2(1), &
        centered(laplace_loglikelihood(x, 0.15_dp + h, 1.3_dp), &
        laplace_loglikelihood(x, 0.15_dp - h, 1.3_dp), h), tol, failures)
    call check_close('laplace scale', s2(2), &
        centered(laplace_loglikelihood(x, 0.15_dp, 1.3_dp + h), &
        laplace_loglikelihood(x, 0.15_dp, 1.3_dp - h), h), tol, failures)

    s2 = logistic_score(x, 0.2_dp, 1.1_dp)
    call check_close('logistic loc', s2(1), &
        centered(logistic_loglikelihood(x, 0.2_dp + h, 1.1_dp), &
        logistic_loglikelihood(x, 0.2_dp - h, 1.1_dp), h), tol, failures)
    call check_close('logistic scale', s2(2), &
        centered(logistic_loglikelihood(x, 0.2_dp, 1.1_dp + h), &
        logistic_loglikelihood(x, 0.2_dp, 1.1_dp - h), h), tol, failures)

    s2 = cauchy_score(x, -0.1_dp, 1.6_dp)
    call check_close('cauchy loc', s2(1), &
        centered(cauchy_loglikelihood(x, -0.1_dp + h, 1.6_dp), &
        cauchy_loglikelihood(x, -0.1_dp - h, 1.6_dp), h), tol, failures)
    call check_close('cauchy scale', s2(2), &
        centered(cauchy_loglikelihood(x, -0.1_dp, 1.6_dp + h), &
        cauchy_loglikelihood(x, -0.1_dp, 1.6_dp - h), h), tol, failures)

    s2 = rayleigh_score(xp, -0.2_dp, 1.3_dp)
    call check_close('rayleigh loc', s2(1), &
        centered(rayleigh_loglikelihood(xp, -0.2_dp + h, 1.3_dp), &
        rayleigh_loglikelihood(xp, -0.2_dp - h, 1.3_dp), h), tol, failures)
    call check_close('rayleigh scale', s2(2), &
        centered(rayleigh_loglikelihood(xp, -0.2_dp, 1.3_dp + h), &
        rayleigh_loglikelihood(xp, -0.2_dp, 1.3_dp - h), h), tol, failures)

    s3 = gamma_score(xp, 2.4_dp, -0.2_dp, 1.3_dp)
    call check_close('gamma a', s3(1), &
        centered(gamma_loglikelihood(xp, 2.4_dp + h, -0.2_dp, 1.3_dp), &
        gamma_loglikelihood(xp, 2.4_dp - h, -0.2_dp, 1.3_dp), h), tol, failures)
    call check_close('gamma loc', s3(2), &
        centered(gamma_loglikelihood(xp, 2.4_dp, -0.2_dp + h, 1.3_dp), &
        gamma_loglikelihood(xp, 2.4_dp, -0.2_dp - h, 1.3_dp), h), tol, failures)
    call check_close('gamma scale', s3(3), &
        centered(gamma_loglikelihood(xp, 2.4_dp, -0.2_dp, 1.3_dp + h), &
        gamma_loglikelihood(xp, 2.4_dp, -0.2_dp, 1.3_dp - h), h), tol, failures)

    s3 = chi2_score(xp, 4.8_dp, -0.2_dp, 1.3_dp)
    call check_close('chi2 df', s3(1), &
        centered(chi2_loglikelihood(xp, 4.8_dp + h, -0.2_dp, 1.3_dp), &
        chi2_loglikelihood(xp, 4.8_dp - h, -0.2_dp, 1.3_dp), h), tol, failures)
    call check_close('chi2 loc', s3(2), &
        centered(chi2_loglikelihood(xp, 4.8_dp, -0.2_dp + h, 1.3_dp), &
        chi2_loglikelihood(xp, 4.8_dp, -0.2_dp - h, 1.3_dp), h), tol, failures)
    call check_close('chi2 scale', s3(3), &
        centered(chi2_loglikelihood(xp, 4.8_dp, -0.2_dp, 1.3_dp + h), &
        chi2_loglikelihood(xp, 4.8_dp, -0.2_dp, 1.3_dp - h), h), tol, failures)

    s3 = t_score(x, 5.5_dp, 0.2_dp, 1.3_dp)
    call check_close('t df', s3(1), &
        centered(t_loglikelihood(x, 5.5_dp + h, 0.2_dp, 1.3_dp), &
        t_loglikelihood(x, 5.5_dp - h, 0.2_dp, 1.3_dp), h), tol, failures)
    call check_close('t loc', s3(2), &
        centered(t_loglikelihood(x, 5.5_dp, 0.2_dp + h, 1.3_dp), &
        t_loglikelihood(x, 5.5_dp, 0.2_dp - h, 1.3_dp), h), tol, failures)
    call check_close('t scale', s3(3), &
        centered(t_loglikelihood(x, 5.5_dp, 0.2_dp, 1.3_dp + h), &
        t_loglikelihood(x, 5.5_dp, 0.2_dp, 1.3_dp - h), h), tol, failures)

    s3 = lognormal_score(xp, 0.8_dp, -0.2_dp, 1.4_dp)
    call check_close('lognormal s', s3(1), &
        centered(lognormal_loglikelihood(xp, 0.8_dp + h, -0.2_dp, 1.4_dp), &
        lognormal_loglikelihood(xp, 0.8_dp - h, -0.2_dp, 1.4_dp), h), tol, failures)
    call check_close('lognormal loc', s3(2), &
        centered(lognormal_loglikelihood(xp, 0.8_dp, -0.2_dp + h, 1.4_dp), &
        lognormal_loglikelihood(xp, 0.8_dp, -0.2_dp - h, 1.4_dp), h), tol, failures)
    call check_close('lognormal scale', s3(3), &
        centered(lognormal_loglikelihood(xp, 0.8_dp, -0.2_dp, 1.4_dp + h), &
        lognormal_loglikelihood(xp, 0.8_dp, -0.2_dp, 1.4_dp - h), h), tol, failures)

    s3 = weibull_score(xp, 1.7_dp, -0.2_dp, 1.5_dp)
    call check_close('weibull c', s3(1), &
        centered(weibull_loglikelihood(xp, 1.7_dp + h, -0.2_dp, 1.5_dp), &
        weibull_loglikelihood(xp, 1.7_dp - h, -0.2_dp, 1.5_dp), h), tol, failures)
    call check_close('weibull loc', s3(2), &
        centered(weibull_loglikelihood(xp, 1.7_dp, -0.2_dp + h, 1.5_dp), &
        weibull_loglikelihood(xp, 1.7_dp, -0.2_dp - h, 1.5_dp), h), tol, failures)
    call check_close('weibull scale', s3(3), &
        centered(weibull_loglikelihood(xp, 1.7_dp, -0.2_dp, 1.5_dp + h), &
        weibull_loglikelihood(xp, 1.7_dp, -0.2_dp, 1.5_dp - h), h), tol, failures)

    s3 = pareto_score(xpareto, 2.3_dp, -0.4_dp, 1.2_dp)
    call check_close('pareto b', s3(1), &
        centered(pareto_loglikelihood(xpareto, 2.3_dp + h, -0.4_dp, 1.2_dp), &
        pareto_loglikelihood(xpareto, 2.3_dp - h, -0.4_dp, 1.2_dp), h), tol, failures)
    call check_close('pareto loc', s3(2), &
        centered(pareto_loglikelihood(xpareto, 2.3_dp, -0.4_dp + h, 1.2_dp), &
        pareto_loglikelihood(xpareto, 2.3_dp, -0.4_dp - h, 1.2_dp), h), tol, failures)
    call check_close('pareto scale', s3(3), &
        centered(pareto_loglikelihood(xpareto, 2.3_dp, -0.4_dp, 1.2_dp + h), &
        pareto_loglikelihood(xpareto, 2.3_dp, -0.4_dp, 1.2_dp - h), h), tol, failures)

    s4 = beta_score(xb, 2.2_dp, 3.1_dp, -0.3_dp, 2.0_dp)
    call check_close('beta a', s4(1), &
        centered(beta_loglikelihood(xb, 2.2_dp + h, 3.1_dp, -0.3_dp, 2.0_dp), &
        beta_loglikelihood(xb, 2.2_dp - h, 3.1_dp, -0.3_dp, 2.0_dp), h), tol, failures)
    call check_close('beta b', s4(2), &
        centered(beta_loglikelihood(xb, 2.2_dp, 3.1_dp + h, -0.3_dp, 2.0_dp), &
        beta_loglikelihood(xb, 2.2_dp, 3.1_dp - h, -0.3_dp, 2.0_dp), h), tol, failures)
    call check_close('beta loc', s4(3), &
        centered(beta_loglikelihood(xb, 2.2_dp, 3.1_dp, -0.3_dp + h, 2.0_dp), &
        beta_loglikelihood(xb, 2.2_dp, 3.1_dp, -0.3_dp - h, 2.0_dp), h), tol, failures)
    call check_close('beta scale', s4(4), &
        centered(beta_loglikelihood(xb, 2.2_dp, 3.1_dp, -0.3_dp, 2.0_dp + h), &
        beta_loglikelihood(xb, 2.2_dp, 3.1_dp, -0.3_dp, 2.0_dp - h), h), tol, failures)

    s4 = f_score(xp, 4.2_dp, 7.5_dp, -0.3_dp, 1.7_dp)
    call check_close('f dfn', s4(1), &
        centered(f_loglikelihood(xp, 4.2_dp + h, 7.5_dp, -0.3_dp, 1.7_dp), &
        f_loglikelihood(xp, 4.2_dp - h, 7.5_dp, -0.3_dp, 1.7_dp), h), tol, failures)
    call check_close('f dfd', s4(2), &
        centered(f_loglikelihood(xp, 4.2_dp, 7.5_dp + h, -0.3_dp, 1.7_dp), &
        f_loglikelihood(xp, 4.2_dp, 7.5_dp - h, -0.3_dp, 1.7_dp), h), tol, failures)
    call check_close('f loc', s4(3), &
        centered(f_loglikelihood(xp, 4.2_dp, 7.5_dp, -0.3_dp + h, 1.7_dp), &
        f_loglikelihood(xp, 4.2_dp, 7.5_dp, -0.3_dp - h, 1.7_dp), h), tol, failures)
    call check_close('f scale', s4(4), &
        centered(f_loglikelihood(xp, 4.2_dp, 7.5_dp, -0.3_dp, 1.7_dp + h), &
        f_loglikelihood(xp, 4.2_dp, 7.5_dp, -0.3_dp, 1.7_dp - h), h), tol, failures)

    call check_close('bernoulli p', first(bernoulli_score(binary, 0.37_dp)), &
        centered(bernoulli_loglikelihood(binary, 0.37_dp + h), &
        bernoulli_loglikelihood(binary, 0.37_dp - h), h), tol, failures)
    call check_close('poisson mu', first(poisson_score(counts0, 2.4_dp)), &
        centered(poisson_loglikelihood(counts0, 2.4_dp + h), &
        poisson_loglikelihood(counts0, 2.4_dp - h), h), tol, failures)
    call check_close('geometric p', first(geometric_score(counts1, 0.31_dp)), &
        centered(geometric_loglikelihood(counts1, 0.31_dp + h), &
        geometric_loglikelihood(counts1, 0.31_dp - h), h), tol, failures)
    call check_close('binomial p', first(binomial_score(counts0, 5, 0.43_dp)), &
        centered(binomial_loglikelihood(counts0, 5, 0.43_dp + h), &
        binomial_loglikelihood(counts0, 5, 0.43_dp - h), h), tol, failures)
    s2 = negative_binomial_score(counts0, 2.7_dp, 0.46_dp)
    call check_close('negative binomial n', s2(1), &
        centered(negative_binomial_loglikelihood(counts0, 2.7_dp + h, 0.46_dp), &
        negative_binomial_loglikelihood(counts0, 2.7_dp - h, 0.46_dp), h), tol, failures)
    call check_close('negative binomial p', s2(2), &
        centered(negative_binomial_loglikelihood(counts0, 2.7_dp, 0.46_dp + h), &
        negative_binomial_loglikelihood(counts0, 2.7_dp, 0.46_dp - h), h), tol, failures)

    call check_close('bernoulli real/int', first(bernoulli_score(real(binary, dp), 0.37_dp)), &
        first(bernoulli_score(binary, 0.37_dp)), 0.0_dp, failures)
    call check_close('poisson real/int', first(poisson_score(real(counts0, dp), 2.4_dp)), &
        first(poisson_score(counts0, 2.4_dp)), 0.0_dp, failures)
    call check_close('geometric real/int', first(geometric_score(real(counts1, dp), 0.31_dp)), &
        first(geometric_score(counts1, 0.31_dp)), 0.0_dp, failures)

    call check_true('bernoulli p boundary', &
        ieee_is_nan(first(bernoulli_score(binary, 0.0_dp))), failures)
    call check_true('poisson mu boundary', &
        ieee_is_nan(first(poisson_score(counts0, 0.0_dp))), failures)
    call check_true('uniform support boundary', &
        all(ieee_is_nan(uniform_score([0.0_dp, 0.5_dp], 0.0_dp, 1.0_dp))), failures)
    s2 = laplace_score([0.0_dp, 1.0_dp], 0.0_dp, 1.0_dp)
    call check_true('laplace loc kink', ieee_is_nan(s2(1)), failures)
    call check_true('laplace scale finite at kink', .not. ieee_is_nan(s2(2)), failures)
    s3 = pareto_score([1.0_dp, 2.0_dp], 2.0_dp, 0.0_dp, 1.0_dp)
    call check_true('pareto shape at boundary', .not. ieee_is_nan(s3(1)), failures)
    call check_true('pareto loc boundary', ieee_is_nan(s3(2)), failures)
    call check_true('pareto scale boundary', ieee_is_nan(s3(3)), failures)

    if (failures /= 0) then
        print '(a,i0)', 'score test failures: ', failures
        error stop 1
    end if
    print '(a)', 'score tests passed'

contains

    pure function first(values) result(value)
        real(dp), intent(in) :: values(:) !! nonempty vector whose first value is requested
        real(dp) :: value

        value = values(1)
    end function first

    pure function centered(plus, minus, step) result(value)
        real(dp), intent(in) :: plus !! objective at parameter + step
        real(dp), intent(in) :: minus !! objective at parameter - step
        real(dp), intent(in) :: step !! positive finite-difference step
        real(dp) :: value

        value = (plus - minus) / (2.0_dp * step)
    end function centered

    subroutine check_close(name, actual, expected, rtol, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: rtol !! relative tolerance with unit absolute floor
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: scale

        scale = max(1.0_dp, abs(expected))
        if (ieee_is_nan(actual) .or. abs(actual - expected) > rtol * scale) then
            print '(a,2(1x,es24.16),1x,es12.4)', trim(name), actual, expected, &
                abs(actual - expected)
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

end program test_score
