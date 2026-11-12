! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_empirical_distribution
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : positive_infinity, quiet_nan
    use scifort_stats, only : confidence_interval, ecdf, ecdf_confidence_interval, &
        ecdf_result, normal_cdf, normal_pdf, normal_ppf, order_statistic, &
        order_statistic_distribution, quantile_test, quantile_test_result
    implicit none

    type(ecdf_confidence_interval) :: eci
    type(ecdf_result) :: eres
    type(order_statistic_distribution) :: os
    type(confidence_interval) :: qci
    type(quantile_test_result) :: qres
    real(dp) :: broken(5)
    real(dp) :: evaluated(5)
    real(dp) :: right(5)
    real(dp) :: sample(5)
    real(dp), allocatable :: x(:)

    eres = ecdf([real(dp) ::])
    call check_scalar(eres%cdf%evaluate(0.0_dp), 0.0_dp, 0.0_dp, 'empty cdf finite')
    call check_scalar(eres%sf%evaluate(0.0_dp), 1.0_dp, 0.0_dp, 'empty sf finite')
    call check_scalar(eres%cdf%evaluate(positive_infinity(0.0_dp)), 1.0_dp, 0.0_dp, &
        'empty cdf positive infinity')

    sample = [6.23_dp, 5.58_dp, 7.06_dp, 6.42_dp, 5.20_dp]
    eres = ecdf(sample)
    call check_vector(eres%cdf%quantiles, [5.20_dp, 5.58_dp, 6.23_dp, 6.42_dp, 7.06_dp], &
        2.0e-14_dp, 'uncensored quantiles')
    call check_vector(eres%cdf%probabilities, [0.2_dp, 0.4_dp, 0.6_dp, 0.8_dp, 1.0_dp], &
        2.0e-14_dp, 'uncensored cdf')
    evaluated = eres%cdf%evaluate([5.0_dp, 5.2_dp, 6.0_dp, 7.5_dp, quiet_nan(0.0_dp)])
    call check_vector(evaluated, [0.0_dp, 0.2_dp, 0.4_dp, 1.0_dp, 1.0_dp], &
        2.0e-14_dp, 'cdf evaluate')

    eci = eres%cdf%confidence_interval(0.9_dp, 'linear')
    call check_vector(eci%low%probabilities(1:4), &
        [0.0_dp, 0.03963063_dp, 0.23963063_dp, 0.50575964_dp], 2.0e-7_dp, &
        'linear Greenwood lower')
    if (.not. ieee_is_nan(eci%low%probabilities(5))) error stop 'linear Greenwood endpoint'

    eci = eres%cdf%confidence_interval(0.9_dp, 'log-log')
    call check_vector(eci%low%probabilities(1:4), &
        [0.04202057_dp, 0.14583471_dp, 0.29020956_dp, 0.47543251_dp], 2.0e-7_dp, &
        'log-log Greenwood lower')

    broken = [77.0_dp, 47.0_dp, 81.0_dp, 56.0_dp, 80.0_dp]
    right = [62.0_dp, 60.0_dp, 43.0_dp, 71.0_dp, 37.0_dp]
    eres = ecdf(broken, right)
    call check_vector(eres%sf%quantiles, &
        [37.0_dp, 43.0_dp, 47.0_dp, 56.0_dp, 60.0_dp, 62.0_dp, 71.0_dp, &
         77.0_dp, 80.0_dp, 81.0_dp], 2.0e-14_dp, 'Kaplan-Meier quantiles')
    call check_vector(eres%sf%probabilities, &
        [1.0_dp, 1.0_dp, 0.875_dp, 0.75_dp, 0.75_dp, 0.75_dp, 0.75_dp, &
         0.5_dp, 0.25_dp, 0.0_dp], 2.0e-14_dp, 'Kaplan-Meier survival')
    eci = eres%sf%confidence_interval(0.95_dp, 'linear')
    call check_vector(eci%low%probabilities(1:9), &
        [1.0_dp, 1.0_dp, 0.64582770_dp, 0.44994302_dp, 0.44994302_dp, &
         0.44994302_dp, 0.44994302_dp, 0.05270146_dp, 0.0_dp], 3.0e-7_dp, &
        'Kaplan-Meier Greenwood lower')
    if (.not. ieee_is_nan(eci%low%probabilities(10))) error stop 'Kaplan-Meier endpoint'

    x = [0.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
    qres = quantile_test(x, q=1.0_dp, p=0.5_dp, alternative='two-sided')
    if (qres%statistic /= 3 .or. qres%statistic_type /= 1) error stop 'quantile two-sided statistic'
    call check_scalar(qres%pvalue, 1.0_dp, 2.0e-14_dp, 'quantile two-sided pvalue')
    qci = qres%confidence_interval(0.9_dp)
    call check_scalar(qci%low, 0.0_dp, 2.0e-14_dp, 'quantile CI low')
    call check_scalar(qci%high, 5.0_dp, 2.0e-14_dp, 'quantile CI high')

    qres = quantile_test(x, q=1.0_dp, p=0.5_dp, alternative='less')
    if (qres%statistic /= 2 .or. qres%statistic_type /= 2) error stop 'quantile less statistic'
    call check_scalar(qres%pvalue, 0.9375_dp, 2.0e-14_dp, 'quantile less pvalue')
    qci = qres%confidence_interval(0.9_dp)
    call check_scalar(qci%high, 4.0_dp, 2.0e-14_dp, 'quantile less CI')

    qres = quantile_test(x, q=1.0_dp, p=quiet_nan(0.0_dp), alternative='greater')
    if (.not. ieee_is_nan(qres%pvalue)) error stop 'quantile NaN p propagation'

    os = order_statistic(base_cdf, base_pdf, base_ppf, 2, 5)
    call check_scalar(os%cdf(0.0_dp), 0.8125_dp, 2.0e-13_dp, 'order statistic cdf')
    call check_scalar(os%pdf(0.0_dp), 0.49867785050179086_dp, 2.0e-13_dp, &
        'order statistic pdf')
    call check_scalar(os%ppf(0.5_dp), -0.48507895332644424_dp, 4.0e-8_dp, &
        'order statistic ppf')
    call check_scalar(os%isf(0.5_dp), -0.48507895332644424_dp, 4.0e-8_dp, &
        'order statistic isf')

    print *, 'empirical-distribution and quantile-inference tests passed'

contains

    pure function base_cdf(value, parameters) result(probability)
        real(dp), intent(in) :: value
        real(dp), intent(in) :: parameters(:)
        real(dp) :: probability

        if (size(parameters) /= 0) error stop 'unexpected normal callback parameters'
        probability = normal_cdf(value)
    end function base_cdf

    pure function base_pdf(value, parameters) result(probability)
        real(dp), intent(in) :: value
        real(dp), intent(in) :: parameters(:)
        real(dp) :: probability

        if (size(parameters) /= 0) error stop 'unexpected normal callback parameters'
        probability = normal_pdf(value)
    end function base_pdf

    pure function base_ppf(probability, parameters) result(value)
        real(dp), intent(in) :: probability
        real(dp), intent(in) :: parameters(:)
        real(dp) :: value

        if (size(parameters) /= 0) error stop 'unexpected normal callback parameters'
        value = normal_ppf(probability)
    end function base_ppf

    subroutine check_scalar(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual
        real(dp), intent(in) :: expected
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_scalar

    subroutine check_vector(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual(:)
        real(dp), intent(in) :: expected(:)
        real(dp), intent(in) :: tolerance
        character(len=*), intent(in) :: label

        if (size(actual) /= size(expected)) then
            print *, 'FAIL ', trim(label), ': size mismatch'
            error stop 1
        end if
        if (any(abs(actual - expected) > tolerance)) then
            print *, 'FAIL ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_vector

end program test_empirical_distribution
