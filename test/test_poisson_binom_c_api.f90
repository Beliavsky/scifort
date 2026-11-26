! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_poisson_binom_c_api
    use iso_c_binding, only : c_double, c_int
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use scifort_poisson_binom
    use scifort_c_api_generated_poisson_binom
    implicit none
    real(c_double) :: p(2), value, expected, argument, logf, derivative(2)
    integer(c_int) :: status
    p = [.2_c_double, .8_c_double]
    argument = 1.5_c_double
    expected = poisson_binom_pmf(argument,p,.5_c_double)
    call scifort_poisson_binom_pmf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = 1.5_c_double
    expected = poisson_binom_logpmf(argument,p,.5_c_double)
    call scifort_poisson_binom_logpmf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = 1.5_c_double
    expected = poisson_binom_cdf(argument,p,.5_c_double)
    call scifort_poisson_binom_cdf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = 1.5_c_double
    expected = poisson_binom_sf(argument,p,.5_c_double)
    call scifort_poisson_binom_sf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = 1.5_c_double
    expected = poisson_binom_logcdf(argument,p,.5_c_double)
    call scifort_poisson_binom_logcdf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = 1.5_c_double
    expected = poisson_binom_logsf(argument,p,.5_c_double)
    call scifort_poisson_binom_logsf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = .6_c_double
    expected = poisson_binom_ppf(argument,p,.5_c_double)
    call scifort_poisson_binom_ppf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    argument = .6_c_double
    expected = poisson_binom_isf(argument,p,.5_c_double)
    call scifort_poisson_binom_isf_explicit_f64(2_c_int,argument,p,.5_c_double,value,status)
    if (status /= 0 .or. .not. ieee_is_finite(value)) stop 1
    if (abs(value-expected) > 1.0e-13_c_double) stop 2
    if (abs(poisson_binom_pmf(1,p)-poisson_binom_pmf(1.0_c_double,p)) > 1.0e-14_c_double) stop 10
    call poisson_binom_logpmf_derivative_p(1.0_c_double,p,logf,derivative)
    if (abs(logf-log(.68_c_double)) > 1.0e-14_c_double) stop 11
    if (any(abs(derivative-[-.6_c_double,.6_c_double]/.68_c_double) > 1.0e-13_c_double)) stop 12
    call poisson_binom_logpmf_derivative_p(0.0_c_double,p(:1),logf,derivative(:1))
    if (abs(derivative(1)+1.25_c_double) > 1.0e-13_c_double) stop 13
    if (poisson_binom_cdf(1.0e100_c_double,p) /= 1.0_c_double) stop 3
    if (poisson_binom_cdf(-1.0e100_c_double,p) /= 0.0_c_double) stop 4
    if (poisson_binom_cdf(-.2_c_double,p) /= 0.0_c_double) stop 5
    call poisson_binom_logpmf_derivative_p(1.0e100_c_double,p,logf,derivative)
    if (ieee_is_finite(logf) .or. ieee_is_nan(logf) .or. logf > 0.0_c_double) stop 6
    if (.not. all(ieee_is_nan(derivative))) stop 7
    call poisson_binom_logpmf_derivative_p(-1.0e100_c_double,p,logf,derivative)
    if (ieee_is_finite(logf) .or. ieee_is_nan(logf) .or. logf > 0.0_c_double) stop 8
    if (.not. all(ieee_is_nan(derivative))) stop 9
    print '(a)', 'Poisson-binomial native/C and conversion regressions passed.'
end program test_poisson_binom_c_api
