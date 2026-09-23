! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_c_api
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use, intrinsic :: iso_c_binding, only : c_double
    use scifort_c_api, only : &
        scifort_normal_pdf_f64, &
        scifort_normal_cdf_f64, &
        scifort_normal_ppf_f64, &
        scifort_uniform_pdf_f64, &
        scifort_uniform_cdf_f64, &
        scifort_uniform_ppf_f64, &
        scifort_exponential_pdf_f64, &
        scifort_exponential_cdf_f64, &
        scifort_exponential_ppf_f64, &
        scifort_laplace_pdf_f64, &
        scifort_laplace_cdf_f64, &
        scifort_laplace_ppf_f64, &
        scifort_logistic_pdf_f64, &
        scifort_logistic_cdf_f64, &
        scifort_logistic_ppf_f64, &
        scifort_cauchy_pdf_f64, &
        scifort_cauchy_cdf_f64, &
        scifort_cauchy_ppf_f64, &
        scifort_rayleigh_pdf_f64, &
        scifort_rayleigh_cdf_f64, &
        scifort_rayleigh_ppf_f64, &
        scifort_gamma_pdf_f64, &
        scifort_gamma_cdf_f64, &
        scifort_gamma_ppf_f64, &
        scifort_chi2_pdf_f64, &
        scifort_chi2_cdf_f64, &
        scifort_chi2_ppf_f64, &
        scifort_t_pdf_f64, &
        scifort_t_cdf_f64, &
        scifort_t_ppf_f64, &
        scifort_lognormal_pdf_f64, &
        scifort_lognormal_cdf_f64, &
        scifort_lognormal_ppf_f64, &
        scifort_weibull_pdf_f64, &
        scifort_weibull_cdf_f64, &
        scifort_weibull_ppf_f64, &
        scifort_pareto_pdf_f64, &
        scifort_pareto_cdf_f64, &
        scifort_pareto_ppf_f64, &
        scifort_beta_pdf_f64, &
        scifort_beta_cdf_f64, &
        scifort_beta_ppf_f64, &
        scifort_f_pdf_f64, &
        scifort_f_cdf_f64, &
        scifort_f_ppf_f64, &
        scifort_bernoulli_pmf_f64, &
        scifort_bernoulli_cdf_f64, &
        scifort_bernoulli_ppf_f64, &
        scifort_poisson_pmf_f64, &
        scifort_poisson_cdf_f64, &
        scifort_poisson_ppf_f64, &
        scifort_geometric_pmf_f64, &
        scifort_geometric_cdf_f64, &
        scifort_geometric_ppf_f64, &
        scifort_binomial_pmf_f64, &
        scifort_binomial_cdf_f64, &
        scifort_binomial_ppf_f64, &
        scifort_negative_binomial_pmf_f64, &
        scifort_negative_binomial_cdf_f64, &
        scifort_negative_binomial_ppf_f64
    use scifort_kinds, only : dp
    use scifort_stats
    implicit none

    integer :: failures

    failures = 0

    call check_close('C API normal pdf', &
        real(scifort_normal_pdf_f64(0.7_c_double, -0.2_c_double, 1.3_c_double), dp), &
        normal_pdf(0.7_dp, -0.2_dp, 1.3_dp), failures)
    call check_close('C API normal cdf', &
        real(scifort_normal_cdf_f64(0.7_c_double, -0.2_c_double, 1.3_c_double), dp), &
        normal_cdf(0.7_dp, -0.2_dp, 1.3_dp), failures)
    call check_close('C API normal ppf', &
        real(scifort_normal_ppf_f64(0.73_c_double, -0.2_c_double, 1.3_c_double), dp), &
        normal_ppf(0.73_dp, -0.2_dp, 1.3_dp), failures)
    call check_close('C API uniform pdf', &
        real(scifort_uniform_pdf_f64(1.2_c_double, 0.5_c_double, 2.0_c_double), dp), &
        uniform_pdf(1.2_dp, 0.5_dp, 2.0_dp), failures)
    call check_close('C API uniform cdf', &
        real(scifort_uniform_cdf_f64(1.2_c_double, 0.5_c_double, 2.0_c_double), dp), &
        uniform_cdf(1.2_dp, 0.5_dp, 2.0_dp), failures)
    call check_close('C API uniform ppf', &
        real(scifort_uniform_ppf_f64(0.30_c_double, 0.5_c_double, 2.0_c_double), dp), &
        uniform_ppf(0.30_dp, 0.5_dp, 2.0_dp), failures)
    call check_close('C API exponential pdf', &
        real(scifort_exponential_pdf_f64(1.2_c_double, 0.1_c_double, 2.3_c_double), dp), &
        exponential_pdf(1.2_dp, 0.1_dp, 2.3_dp), failures)
    call check_close('C API exponential cdf', &
        real(scifort_exponential_cdf_f64(1.2_c_double, 0.1_c_double, 2.3_c_double), dp), &
        exponential_cdf(1.2_dp, 0.1_dp, 2.3_dp), failures)
    call check_close('C API exponential ppf', &
        real(scifort_exponential_ppf_f64(0.40_c_double, 0.1_c_double, 2.3_c_double), dp), &
        exponential_ppf(0.40_dp, 0.1_dp, 2.3_dp), failures)
    call check_close('C API laplace pdf', &
        real(scifort_laplace_pdf_f64(0.8_c_double, -0.3_c_double, 1.7_c_double), dp), &
        laplace_pdf(0.8_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API laplace cdf', &
        real(scifort_laplace_cdf_f64(0.8_c_double, -0.3_c_double, 1.7_c_double), dp), &
        laplace_cdf(0.8_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API laplace ppf', &
        real(scifort_laplace_ppf_f64(0.65_c_double, -0.3_c_double, 1.7_c_double), dp), &
        laplace_ppf(0.65_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API logistic pdf', &
        real(scifort_logistic_pdf_f64(0.8_c_double, -0.3_c_double, 1.7_c_double), dp), &
        logistic_pdf(0.8_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API logistic cdf', &
        real(scifort_logistic_cdf_f64(0.8_c_double, -0.3_c_double, 1.7_c_double), dp), &
        logistic_cdf(0.8_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API logistic ppf', &
        real(scifort_logistic_ppf_f64(0.65_c_double, -0.3_c_double, 1.7_c_double), dp), &
        logistic_ppf(0.65_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API cauchy pdf', &
        real(scifort_cauchy_pdf_f64(0.8_c_double, -0.3_c_double, 1.7_c_double), dp), &
        cauchy_pdf(0.8_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API cauchy cdf', &
        real(scifort_cauchy_cdf_f64(0.8_c_double, -0.3_c_double, 1.7_c_double), dp), &
        cauchy_cdf(0.8_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API cauchy ppf', &
        real(scifort_cauchy_ppf_f64(0.65_c_double, -0.3_c_double, 1.7_c_double), dp), &
        cauchy_ppf(0.65_dp, -0.3_dp, 1.7_dp), failures)
    call check_close('C API rayleigh pdf', &
        real(scifort_rayleigh_pdf_f64(1.4_c_double, 0.2_c_double, 0.9_c_double), dp), &
        rayleigh_pdf(1.4_dp, 0.2_dp, 0.9_dp), failures)
    call check_close('C API rayleigh cdf', &
        real(scifort_rayleigh_cdf_f64(1.4_c_double, 0.2_c_double, 0.9_c_double), dp), &
        rayleigh_cdf(1.4_dp, 0.2_dp, 0.9_dp), failures)
    call check_close('C API rayleigh ppf', &
        real(scifort_rayleigh_ppf_f64(0.55_c_double, 0.2_c_double, 0.9_c_double), dp), &
        rayleigh_ppf(0.55_dp, 0.2_dp, 0.9_dp), failures)
    call check_close('C API gamma pdf', &
        real(scifort_gamma_pdf_f64(3.2_c_double, 2.5_c_double, 0.2_c_double, 1.3_c_double), dp), &
        gamma_pdf(3.2_dp, 2.5_dp, 0.2_dp, 1.3_dp), failures)
    call check_close('C API gamma cdf', &
        real(scifort_gamma_cdf_f64(3.2_c_double, 2.5_c_double, 0.2_c_double, 1.3_c_double), dp), &
        gamma_cdf(3.2_dp, 2.5_dp, 0.2_dp, 1.3_dp), failures)
    call check_close('C API gamma ppf', &
        real(scifort_gamma_ppf_f64(0.60_c_double, 2.5_c_double, 0.2_c_double, 1.3_c_double), dp), &
        gamma_ppf(0.60_dp, 2.5_dp, 0.2_dp, 1.3_dp), failures)
    call check_close('C API chi2 pdf', &
        real(scifort_chi2_pdf_f64(4.0_c_double, 3.5_c_double, 0.1_c_double, 1.2_c_double), dp), &
        chi2_pdf(4.0_dp, 3.5_dp, 0.1_dp, 1.2_dp), failures)
    call check_close('C API chi2 cdf', &
        real(scifort_chi2_cdf_f64(4.0_c_double, 3.5_c_double, 0.1_c_double, 1.2_c_double), dp), &
        chi2_cdf(4.0_dp, 3.5_dp, 0.1_dp, 1.2_dp), failures)
    call check_close('C API chi2 ppf', &
        real(scifort_chi2_ppf_f64(0.70_c_double, 3.5_c_double, 0.1_c_double, 1.2_c_double), dp), &
        chi2_ppf(0.70_dp, 3.5_dp, 0.1_dp, 1.2_dp), failures)
    call check_close('C API t pdf', &
        real(scifort_t_pdf_f64(0.5_c_double, 7.0_c_double, -0.1_c_double, 1.4_c_double), dp), &
        t_pdf(0.5_dp, 7.0_dp, -0.1_dp, 1.4_dp), failures)
    call check_close('C API t cdf', &
        real(scifort_t_cdf_f64(0.5_c_double, 7.0_c_double, -0.1_c_double, 1.4_c_double), dp), &
        t_cdf(0.5_dp, 7.0_dp, -0.1_dp, 1.4_dp), failures)
    call check_close('C API t ppf', &
        real(scifort_t_ppf_f64(0.80_c_double, 7.0_c_double, -0.1_c_double, 1.4_c_double), dp), &
        t_ppf(0.80_dp, 7.0_dp, -0.1_dp, 1.4_dp), failures)
    call check_close('C API lognormal pdf', &
        real(scifort_lognormal_pdf_f64( &
            2.0_c_double, 0.8_c_double, 0.2_c_double, 1.1_c_double), dp), &
        lognormal_pdf(2.0_dp, 0.8_dp, 0.2_dp, 1.1_dp), failures)
    call check_close('C API lognormal cdf', &
        real(scifort_lognormal_cdf_f64( &
            2.0_c_double, 0.8_c_double, 0.2_c_double, 1.1_c_double), dp), &
        lognormal_cdf(2.0_dp, 0.8_dp, 0.2_dp, 1.1_dp), failures)
    call check_close('C API lognormal ppf', &
        real(scifort_lognormal_ppf_f64( &
            0.60_c_double, 0.8_c_double, 0.2_c_double, 1.1_c_double), dp), &
        lognormal_ppf(0.60_dp, 0.8_dp, 0.2_dp, 1.1_dp), failures)
    call check_close('C API weibull pdf', &
        real(scifort_weibull_pdf_f64(1.4_c_double, 1.7_c_double, 0.2_c_double, 0.9_c_double), dp), &
        weibull_pdf(1.4_dp, 1.7_dp, 0.2_dp, 0.9_dp), failures)
    call check_close('C API weibull cdf', &
        real(scifort_weibull_cdf_f64(1.4_c_double, 1.7_c_double, 0.2_c_double, 0.9_c_double), dp), &
        weibull_cdf(1.4_dp, 1.7_dp, 0.2_dp, 0.9_dp), failures)
    call check_close('C API weibull ppf', &
        real(scifort_weibull_ppf_f64( &
            0.40_c_double, 1.7_c_double, 0.2_c_double, 0.9_c_double), dp), &
        weibull_ppf(0.40_dp, 1.7_dp, 0.2_dp, 0.9_dp), failures)
    call check_close('C API pareto pdf', &
        real(scifort_pareto_pdf_f64(3.0_c_double, 2.2_c_double, 0.1_c_double, 1.1_c_double), dp), &
        pareto_pdf(3.0_dp, 2.2_dp, 0.1_dp, 1.1_dp), failures)
    call check_close('C API pareto cdf', &
        real(scifort_pareto_cdf_f64(3.0_c_double, 2.2_c_double, 0.1_c_double, 1.1_c_double), dp), &
        pareto_cdf(3.0_dp, 2.2_dp, 0.1_dp, 1.1_dp), failures)
    call check_close('C API pareto ppf', &
        real(scifort_pareto_ppf_f64(0.70_c_double, 2.2_c_double, 0.1_c_double, 1.1_c_double), dp), &
        pareto_ppf(0.70_dp, 2.2_dp, 0.1_dp, 1.1_dp), failures)
    call check_close('C API beta pdf', &
        real(scifort_beta_pdf_f64( &
            0.6_c_double, 2.5_c_double, 3.1_c_double, 0.1_c_double, 0.9_c_double), dp), &
        beta_pdf(0.6_dp, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp), failures)
    call check_close('C API beta cdf', &
        real(scifort_beta_cdf_f64( &
            0.6_c_double, 2.5_c_double, 3.1_c_double, 0.1_c_double, 0.9_c_double), dp), &
        beta_cdf(0.6_dp, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp), failures)
    call check_close('C API beta ppf', &
        real(scifort_beta_ppf_f64( &
            0.55_c_double, 2.5_c_double, 3.1_c_double, 0.1_c_double, 0.9_c_double), dp), &
        beta_ppf(0.55_dp, 2.5_dp, 3.1_dp, 0.1_dp, 0.9_dp), failures)
    call check_close('C API f pdf', &
        real(scifort_f_pdf_f64( &
            1.3_c_double, 5.0_c_double, 9.0_c_double, 0.2_c_double, 1.1_c_double), dp), &
        f_pdf(1.3_dp, 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp), failures)
    call check_close('C API f cdf', &
        real(scifort_f_cdf_f64( &
            1.3_c_double, 5.0_c_double, 9.0_c_double, 0.2_c_double, 1.1_c_double), dp), &
        f_cdf(1.3_dp, 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp), failures)
    call check_close('C API f ppf', &
        real(scifort_f_ppf_f64( &
            0.60_c_double, 5.0_c_double, 9.0_c_double, 0.2_c_double, 1.1_c_double), dp), &
        f_ppf(0.60_dp, 5.0_dp, 9.0_dp, 0.2_dp, 1.1_dp), failures)
    call check_close('C API bernoulli pmf', &
        real(scifort_bernoulli_pmf_f64(3.0_c_double, 0.3_c_double, 2.0_c_double), dp), &
        bernoulli_pmf(3.0_dp, 0.3_dp, 2.0_dp), failures)
    call check_close('C API bernoulli cdf', &
        real(scifort_bernoulli_cdf_f64(3.0_c_double, 0.3_c_double, 2.0_c_double), dp), &
        bernoulli_cdf(3.0_dp, 0.3_dp, 2.0_dp), failures)
    call check_close('C API bernoulli ppf', &
        real(scifort_bernoulli_ppf_f64(0.60_c_double, 0.3_c_double, 2.0_c_double), dp), &
        bernoulli_ppf(0.60_dp, 0.3_dp, 2.0_dp), failures)
    call check_close('C API poisson pmf', &
        real(scifort_poisson_pmf_f64(5.0_c_double, 2.5_c_double, 1.0_c_double), dp), &
        poisson_pmf(5.0_dp, 2.5_dp, 1.0_dp), failures)
    call check_close('C API poisson cdf', &
        real(scifort_poisson_cdf_f64(5.0_c_double, 2.5_c_double, 1.0_c_double), dp), &
        poisson_cdf(5.0_dp, 2.5_dp, 1.0_dp), failures)
    call check_close('C API poisson ppf', &
        real(scifort_poisson_ppf_f64(0.60_c_double, 2.5_c_double, 1.0_c_double), dp), &
        poisson_ppf(0.60_dp, 2.5_dp, 1.0_dp), failures)
    call check_close('C API geometric pmf', &
        real(scifort_geometric_pmf_f64(5.0_c_double, 0.3_c_double, 2.0_c_double), dp), &
        geometric_pmf(5.0_dp, 0.3_dp, 2.0_dp), failures)
    call check_close('C API geometric cdf', &
        real(scifort_geometric_cdf_f64(5.0_c_double, 0.3_c_double, 2.0_c_double), dp), &
        geometric_cdf(5.0_dp, 0.3_dp, 2.0_dp), failures)
    call check_close('C API geometric ppf', &
        real(scifort_geometric_ppf_f64(0.60_c_double, 0.3_c_double, 2.0_c_double), dp), &
        geometric_ppf(0.60_dp, 0.3_dp, 2.0_dp), failures)
    call check_close('C API binomial pmf', &
        real(scifort_binomial_pmf_f64( &
            5.0_c_double, 7.0_c_double, 0.4_c_double, 1.0_c_double), dp), &
        binomial_pmf(5.0_dp, 7.0_dp, 0.4_dp, 1.0_dp), failures)
    call check_close('C API binomial cdf', &
        real(scifort_binomial_cdf_f64( &
            5.0_c_double, 7.0_c_double, 0.4_c_double, 1.0_c_double), dp), &
        binomial_cdf(5.0_dp, 7.0_dp, 0.4_dp, 1.0_dp), failures)
    call check_close('C API binomial ppf', &
        real(scifort_binomial_ppf_f64( &
            0.60_c_double, 7.0_c_double, 0.4_c_double, 1.0_c_double), dp), &
        binomial_ppf(0.60_dp, 7.0_dp, 0.4_dp, 1.0_dp), failures)
    call check_close('C API negative_binomial pmf', &
        real(scifort_negative_binomial_pmf_f64( &
            6.0_c_double, 2.5_c_double, 0.4_c_double, 1.0_c_double), dp), &
        negative_binomial_pmf(6.0_dp, 2.5_dp, 0.4_dp, 1.0_dp), failures)
    call check_close('C API negative_binomial cdf', &
        real(scifort_negative_binomial_cdf_f64( &
            6.0_c_double, 2.5_c_double, 0.4_c_double, 1.0_c_double), dp), &
        negative_binomial_cdf(6.0_dp, 2.5_dp, 0.4_dp, 1.0_dp), failures)
    call check_close('C API negative_binomial ppf', &
        real(scifort_negative_binomial_ppf_f64( &
            0.60_c_double, 2.5_c_double, 0.4_c_double, 1.0_c_double), dp), &
        negative_binomial_ppf(0.60_dp, 2.5_dp, 0.4_dp, 1.0_dp), failures)

    call check_true('C API invalid gamma shape', ieee_is_nan(real( &
        scifort_gamma_pdf_f64(1.0_c_double, -1.0_c_double, 0.0_c_double, &
        1.0_c_double), dp)), failures)
    call check_true('C API invalid binomial probability', ieee_is_nan(real( &
        scifort_binomial_pmf_f64(1.0_c_double, 4.0_c_double, 1.5_c_double, &
        0.0_c_double), dp)), failures)

    if (failures == 0) then
        print '(a)', 'test_c_api: PASS'
    else
        print '(a,1x,i0)', 'test_c_api: FAIL', failures
        error stop 1
    end if

contains

    subroutine check_close(name, actual, expected, failures)
        character(len=*), intent(in) :: name !! label printed when the check fails
        real(dp), intent(in) :: actual !! value returned by the C ABI wrapper
        real(dp), intent(in) :: expected !! value returned by the native Fortran API
        integer, intent(inout) :: failures !! running count of failed checks

        real(dp) :: tol

        tol = 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(expected))
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tol) then
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

end program test_c_api
