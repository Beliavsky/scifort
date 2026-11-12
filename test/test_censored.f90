! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_censored
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
    use scifort_stats, only : dp, censored_data, make_censored_data, censored_data_valid, &
        censored_size, normal_censored_loglikelihood, exponential_censored_loglikelihood, &
        weibull_censored_loglikelihood, normal_loglikelihood, exponential_loglikelihood, &
        weibull_loglikelihood, normal_fit, exponential_fit, weibull_fit, fit_result, &
        exponential_fit_right_censored, fit_status_no_finite_mle, &
        fit_status_invalid_input, fit_status_no_finite_objective, ecdf, ecdf_result, &
        logrank, logrank_result, censored_loglikelihood, normal_logpdf, normal_logcdf, normal_logsf
    use scifort_stats, only : lognormal_censored_loglikelihood, lognormal_censored_nnlf, &
        lognormal_loglikelihood, lognormal_fit
    use scifort_stats, only : gamma_censored_loglikelihood, gamma_censored_nnlf, &
        gamma_loglikelihood, gamma_fit, gamma_logcdf, gamma_logsf
    use scifort_stats, only : beta_censored_loglikelihood, beta_censored_nnlf, &
        beta_loglikelihood, beta_fit, beta_logcdf, beta_logsf, uniform_loglikelihood
    use scifort_stats, only : logistic_censored_loglikelihood, logistic_censored_nnlf, &
        logistic_loglikelihood, logistic_fit, logistic_logcdf, logistic_logsf
    use scifort_math, only : quiet_nan, negative_infinity, positive_infinity
    implicit none

    type(censored_data) :: data, other, empty
    type(fit_result) :: fit
    type(ecdf_result) :: empirical
    type(logrank_result) :: rank_result
    real(dp) :: inf, nan, interval(2,2), value, base
    integer :: failures, checks, status

    failures = 0
    checks = 0
    inf = positive_infinity(0.0_dp)
    nan = quiet_nan(0.0_dp)
    interval = reshape([1.1_dp, 4.0_dp, 1.4_dp, 4.4_dp], [2,2])
    data = make_censored_data(uncensored=[0.5_dp, 1.2_dp, 2.5_dp], &
        left=[0.3_dp], right=[3.5_dp], interval=interval, status=status)
    call check(status == 0 .and. censored_data_valid(data), 'mixed-data validation')
    call check(censored_size(data) == 7, 'observation count')
    call check(data%n_observations() == 7, 'existing count method includes new groups')
    empirical = ecdf(data)
    call check(all(ieee_is_nan(empirical%cdf%probabilities)), &
        'ECDF rejects unsupported censoring without stopping')
    rank_result = logrank(data, [1.0_dp, 2.0_dp])
    call check(ieee_is_nan(rank_result%pvalue), 'logrank rejects mixed censoring')
    rank_result = logrank([1.0_dp, 2.0_dp], data)
    call check(ieee_is_nan(rank_result%pvalue), 'logrank reversed mixed censoring')
    rank_result = logrank(data, data)
    call check(ieee_is_nan(rank_result%pvalue), 'logrank two mixed samples')

    ! mpmath 1.3.0, 100 digits; tools/generate_censored_reference.py.
    call close(normal_censored_loglikelihood(data, 0.4_dp, 1.3_dp), &
        -19.500020468907458_dp, 2.0e-12_dp, 'normal mixed likelihood')
    call close(censored_loglikelihood(data, [0.4_dp, 1.3_dp], &
        custom_pdf, custom_cdf, custom_sf), -19.500020468907458_dp, 2.0e-12_dp, &
        'user callback likelihood and generic log interval')
    call close(exponential_censored_loglikelihood(data, 0.0_dp, 2.3_dp), &
        -14.101622283408165_dp, 2.0e-12_dp, 'exponential mixed likelihood')
    call close(weibull_censored_loglikelihood(data, 1.7_dp, 0.0_dp, 2.3_dp), &
        -15.3162874574001_dp, 2.0e-12_dp, 'Weibull mixed likelihood')
    call close(weibull_censored_loglikelihood(data, 1.0_dp, 0.0_dp, 2.3_dp), &
        exponential_censored_loglikelihood(data, 0.0_dp, 2.3_dp), &
        2.0e-12_dp, 'shape-one Weibull identity')

    call normal_fit(data, [-4.0_dp, 0.1_dp], [6.0_dp, 6.0_dp], fit, &
        guess=[1.0_dp, 2.0_dp], tolerance=1.0e-9_dp)
    call check(fit%success, 'normal mixed fit')
    ! SciPy 1.15.3 fit(CensoredData(...)); compare parameters, not optimizer paths.
    call close(fit%params(1), 1.9216994688285323_dp, 3.0e-4_dp, 'normal fitted mean')
    call close(fit%params(2), 1.8362526868728573_dp, 3.0e-4_dp, 'normal fitted scale')
    call exponential_fit(data, [0.0_dp, 0.1_dp], [0.0_dp, 10.0_dp], fit)
    call check(fit%success, 'exponential mixed fit')
    call close(fit%params(2), 2.214550781250003_dp, 2.0e-4_dp, 'exponential mixed scale')
    call weibull_fit(data, [0.2_dp, 0.0_dp, 0.1_dp], [6.0_dp, 0.0_dp, 10.0_dp], fit, &
        guess=[1.5_dp, 0.0_dp, 2.0_dp], tolerance=1.0e-9_dp)
    call check(fit%success, 'Weibull mixed fit')
    call close(fit%params(1), 1.0192387961160727_dp, 3.0e-4_dp, 'Weibull fitted shape')
    call close(fit%params(3), 2.22130765498472_dp, 3.0e-4_dp, 'Weibull fitted scale')

    data = make_censored_data(uncensored=[-1.0_dp, 0.2_dp, 1.4_dp])
    call close(normal_censored_loglikelihood(data, 0.4_dp, 1.3_dp), &
        normal_loglikelihood(data%uncensored, 0.4_dp, 1.3_dp), 1.0e-14_dp, &
        'uncensored normal agreement')
    data = make_censored_data(uncensored=[0.2_dp, 1.2_dp, 3.4_dp])
    call close(exponential_censored_loglikelihood(data), &
        exponential_loglikelihood(data%uncensored), 1.0e-14_dp, 'uncensored exponential')
    call close(weibull_censored_loglikelihood(data, 1.7_dp), &
        weibull_loglikelihood(data%uncensored, 1.7_dp), 1.0e-14_dp, 'uncensored Weibull')

    other = make_censored_data(interval=reshape([0.2_dp, 1.2_dp, 3.4_dp, &
        0.2_dp, 1.2_dp, 3.4_dp], [3,2]))
    call close(normal_censored_loglikelihood(other), normal_censored_loglikelihood(data), &
        1.0e-14_dp, 'equal interval endpoints represent exact observations')
    data = make_censored_data(left=[-2.0_dp], right=[3.0_dp])
    other = make_censored_data(interval=reshape([-inf, 3.0_dp, -2.0_dp, inf], [2,2]))
    call close(normal_censored_loglikelihood(data), normal_censored_loglikelihood(other), &
        1.0e-14_dp, 'infinite interval endpoints match left and right groups')

    ! mpmath references use exact binary64 endpoints, at up to 700 digits.
    call set_interval(0.0_dp, 1.0e-20_dp)
    call close(normal_censored_loglikelihood(data), -46.97064039308559_dp, &
        2.0e-13_dp, 'normal interval with equal rounded CDFs')
    call set_interval(40.0_dp, 40.01_dp)
    call close(normal_censored_loglikelihood(data), -805.7174659453684_dp, &
        3.0e-11_dp, 'normal extreme-tail interval')
    call set_interval(1000.0_dp, 1000.0_dp + 1.0e-10_dp)
    call close(exponential_censored_loglikelihood(data), -1023.0254068566554_dp, &
        3.0e-12_dp, 'exponential narrow extreme-tail interval')
    call set_interval(0.0_dp, 1.0e-300_dp)
    call close(exponential_censored_loglikelihood(data, scale=1.0e300_dp), &
        -1381.5510557964274_dp, 2.0e-12_dp, 'interval with underflowed normalized width')
    call set_interval(100.0_dp, 100.0_dp + 1.0e-8_dp)
    call close(weibull_censored_loglikelihood(data, 2.0_dp), -10013.122365005156_dp, &
        3.0e-10_dp, 'Weibull narrow extreme-tail interval')
    call set_interval(1.0_dp, nearest(1.0_dp, 1.0_dp))
    call close(weibull_censored_loglikelihood(data, 0.01_dp), -41.64882357510525_dp, &
        3.0e-13_dp, 'Weibull adjacent-float interval')
    call set_interval(-1.0_dp, 0.0_dp)
    call check(exponential_censored_loglikelihood(data) == negative_infinity(0.0_dp), &
        'impossible support interval')
    call set_interval(-1.0_dp, 1.0_dp)
    value = exponential_censored_loglikelihood(data)
    call set_interval(0.0_dp, 1.0_dp)
    call close(value, exponential_censored_loglikelihood(data), 1.0e-15_dp, &
        'support-clipped interval')
    base = normal_censored_loglikelihood(data)
    call set_interval(0.0_dp, 2.0_dp)
    call check(normal_censored_loglikelihood(data) > base, 'nested interval monotonicity')
    call set_interval(2.0_dp, 3.0_dp)
    call check(weibull_censored_loglikelihood(data, huge(1.0_dp)) == &
        negative_infinity(0.0_dp), 'huge finite Weibull shape does not overflow internally')

    data = make_censored_data(uncensored=[13.0_dp, 22.0_dp, 17.0_dp, 15.0_dp], &
        right=[20.0_dp, 18.0_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(fit%success, 'analytic right-censored fit')
    call close(fit%params(2), 26.25_dp, 1.0e-14_dp, 'total exposure divided by failures')
    call close(fit%nllf, 4.0_dp*(log(26.25_dp) + 1.0_dp), 2.0e-14_dp, &
        'analytic fit likelihood')
    call exponential_fit(data, [0.0_dp, 0.1_dp], [0.0_dp, 100.0_dp], fit)
    call check(fit%success, 'numerical right-censored fit')
    call close(fit%params(2), 26.25_dp, 3.0e-5_dp, 'analytic and numerical fit agree')
    data%uncensored = data%uncensored + 5.0_dp
    data%right = data%right + 5.0_dp
    call exponential_fit_right_censored(data, 5.0_dp, fit)
    call close(fit%params(2), 26.25_dp, 1.0e-14_dp, 'fixed-location translation')
    data = make_censored_data(uncensored=[1.0e308_dp, 1.0e308_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(fit%success, 'finite MLE despite overflow-prone total exposure')
    call close(fit%params(2)/1.0e308_dp, 1.0_dp, 1.0e-14_dp, 'scaled exposure sum')
    data = make_censored_data(uncensored=[nearest(0.0_dp, 1.0_dp), &
        nearest(0.0_dp, 1.0_dp)])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(fit%success, 'subnormal scale estimate')
    call check(fit%params(2) == nearest(0.0_dp, 1.0_dp), 'no premature exposure underflow')
    data = make_censored_data(uncensored=[1.0_dp], right=[-2.0_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call close(fit%params(2), 1.0_dp, 1.0e-15_dp, 'right censor below support adds no exposure')

    data = make_censored_data(right=[1.0_dp, 2.0_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(.not. fit%success .and. fit%status == fit_status_no_finite_mle, &
        'all right-censored has no finite identifiable MLE')
    data = make_censored_data(uncensored=[0.0_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(fit%status == fit_status_no_finite_mle, 'zero-exposure boundary optimum')
    data = make_censored_data(uncensored=[-1.0_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(fit%status == fit_status_no_finite_objective, 'exact value below support')
    data = make_censored_data(left=[1.0_dp])
    call exponential_fit_right_censored(data, 0.0_dp, fit)
    call check(fit%status == fit_status_invalid_input, 'analytic fit rejects left censoring')

    call close(normal_censored_loglikelihood(empty), 0.0_dp, 0.0_dp, 'empty likelihood')
    call normal_fit(empty, [-1.0_dp, 0.1_dp], [1.0_dp, 2.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'empty fit rejected')
    call check(ieee_is_nan(normal_censored_loglikelihood(empty, scale=0.0_dp)), &
        'invalid parameters rejected even with no observations')
    call check(ieee_is_nan(weibull_censored_loglikelihood(empty, nan)), 'NaN shape')
    call check(ieee_is_nan(exponential_censored_loglikelihood(empty, loc=inf)), 'infinite loc')
    data = make_censored_data(uncensored=[nan], status=status)
    call check(status /= 0, 'NaN exact value rejected')
    call check(ieee_is_nan(normal_censored_loglikelihood(data)), 'invalid data yields NaN')
    data = make_censored_data(right=[inf])
    call check(.not. censored_data_valid(data), 'infinite explicit censor limit rejected')
    call set_interval(-inf, inf)
    call check(.not. censored_data_valid(data), 'doubly unbounded interval rejected')
    call set_interval(2.0_dp, 1.0_dp)
    call check(.not. censored_data_valid(data), 'reversed interval rejected')
    call set_interval(nan, 1.0_dp)
    call check(.not. censored_data_valid(data), 'NaN interval rejected')
    data = make_censored_data(interval=reshape([1.0_dp], [1,1]))
    call check(.not. censored_data_valid(data), 'invalid interval shape rejected')

    data = make_censored_data(uncensored=[0.5_dp, 1.2_dp, 2.5_dp], &
        left=[0.3_dp], right=[3.5_dp], interval=interval)
    call close(lognormal_censored_loglikelihood(data, 0.8_dp, scale=2.3_dp), &
        -16.811522795614245_dp, 3.0e-12_dp, 'lognormal mixed reference')
    call close(lognormal_censored_nnlf(data, 0.8_dp, scale=2.3_dp), &
        16.811522795614245_dp, 3.0e-12_dp, 'lognormal negative likelihood')
    call lognormal_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], &
        fit, guess=[1.0_dp, 0.0_dp, 1.0_dp], tolerance=1.0e-9_dp)
    call check(fit%success, 'lognormal mixed fit')
    call close(fit%params(1), 1.2154461568899246_dp, 3.0e-4_dp, 'lognormal fitted shape')
    call close(fit%params(3), 1.373258023610082_dp, 3.0e-4_dp, 'lognormal fitted scale')
    call close(fit%nllf, lognormal_censored_nnlf(data, fit%params(1), &
        fit%params(2), fit%params(3)), 1.0e-13_dp, 'lognormal fit objective')
    call lognormal_fit(data, [0.1_dp, -1.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'lognormal free location rejected')
    call lognormal_fit(data, [0.1_dp, nan, 0.1_dp], [5.0_dp, nan, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'lognormal NaN location bounds rejected')
    call lognormal_fit(empty, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'lognormal empty fit rejected')
    data = make_censored_data(uncensored=exp([-1.0_dp, 0.0_dp, 1.0_dp]))
    call close(lognormal_censored_loglikelihood(data, 0.8_dp), &
        lognormal_loglikelihood(data%uncensored, 0.8_dp), 2.0e-14_dp, &
        'lognormal uncensored array agreement')
    call lognormal_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], [3.0_dp, 0.0_dp, 5.0_dp], &
        fit, tolerance=1.0e-9_dp)
    call check(fit%success, 'lognormal exact-only fit')
    call close(fit%params(1), sqrt(2.0_dp/3.0_dp), 2.0e-7_dp, 'analytic exact shape')
    call close(fit%params(3), 1.0_dp, 2.0e-7_dp, 'analytic exact scale')
    call set_interval(1.0e100_dp, nearest(1.0e100_dp, 1.0_dp))
    call close(lognormal_censored_loglikelihood(data, 0.8_dp, scale=1.0e100_dp), &
        -36.87309372451035_dp, 3.0e-13_dp, 'lognormal adjacent endpoints')
    call set_interval(exp(40.0_dp), exp(40.01_dp))
    call close(lognormal_censored_loglikelihood(data, 1.0_dp), &
        -805.7174659453684_dp, 4.0e-11_dp, 'lognormal extreme tail interval')
    data = make_censored_data(left=[-1.0_dp])
    call check(lognormal_censored_loglikelihood(data, 1.0_dp) == -inf, &
        'lognormal left limit below support')
    data = make_censored_data(right=[0.0_dp])
    call close(lognormal_censored_loglikelihood(data, 1.0_dp), 0.0_dp, 0.0_dp, &
        'lognormal right limit at support')
    data = make_censored_data(uncensored=[0.0_dp])
    call check(lognormal_censored_loglikelihood(data, 1.0_dp) == -inf, &
        'lognormal exact support endpoint')
    call set_interval(-2.0_dp, -1.0_dp)
    call check(lognormal_censored_loglikelihood(data, 1.0_dp) == -inf, &
        'lognormal interval below support')
    call set_interval(-2.0_dp, 1.0_dp)
    call close(lognormal_censored_loglikelihood(data, 1.0_dp), log(0.5_dp), &
        2.0e-15_dp, 'lognormal interval crossing support')
    data = make_censored_data(left=[0.5_dp], right=[2.0_dp])
    other = make_censored_data(interval=reshape([-inf, 2.0_dp, 0.5_dp, inf], [2,2]))
    call close(lognormal_censored_loglikelihood(data, 0.8_dp), &
        lognormal_censored_loglikelihood(other, 0.8_dp), 1.0e-14_dp, &
        'lognormal infinite interval endpoints')
    data = make_censored_data(uncensored=[1.0_dp])
    call set_interval(1.0_dp, 1.0_dp)
    other = make_censored_data(uncensored=[1.0_dp])
    call close(lognormal_censored_loglikelihood(data, 0.8_dp), &
        lognormal_censored_loglikelihood(other, 0.8_dp), 1.0e-14_dp, &
        'lognormal equal endpoints exact density')
    call close(lognormal_censored_loglikelihood(empty, 1.0_dp), 0.0_dp, 0.0_dp, &
        'lognormal empty likelihood')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, 0.0_dp)), &
        'lognormal invalid zero shape')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, -1.0_dp)), &
        'lognormal invalid negative shape')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, nan)), &
        'lognormal NaN shape')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, inf)), &
        'lognormal infinite shape')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, 1.0_dp, loc=inf)), &
        'lognormal infinite location')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, 1.0_dp, scale=0.0_dp)), &
        'lognormal invalid scale')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, 1.0_dp, scale=nan)), &
        'lognormal NaN scale')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, 1.0_dp, loc=nan)), &
        'lognormal NaN location')
    call check(ieee_is_nan(lognormal_censored_loglikelihood(empty, 1.0_dp, scale=inf)), &
        'lognormal infinite scale')
    data = make_censored_data(uncensored=[exp(1.0_dp)])
    call check(lognormal_censored_loglikelihood(data, tiny(1.0_dp)) == -inf, &
        'lognormal tiny shape unrepresentable exact log density')
    call check(lognormal_censored_loglikelihood(data, nearest(0.0_dp, 1.0_dp)) == -inf, &
        'lognormal subnormal shape standardization guard')
    data = make_censored_data(left=[exp(-1.0_dp)])
    call check(lognormal_censored_loglikelihood(data, tiny(1.0_dp)) == -inf, &
        'lognormal tiny shape extreme left tail')
    data = make_censored_data(right=[exp(1.0_dp)])
    call check(lognormal_censored_loglikelihood(data, tiny(1.0_dp)) == -inf, &
        'lognormal tiny shape extreme right tail')
    data = make_censored_data(uncensored=[1.0_dp])
    call close(lognormal_censored_loglikelihood(data, tiny(1.0_dp)), &
        -log(tiny(1.0_dp))-0.91893853320467274178_dp, 2.0e-13_dp, &
        'lognormal tiny shape finite central log density')
    data = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(lognormal_censored_loglikelihood(data, 1.0_dp)), &
        'lognormal invalid observation')

    data = make_censored_data(uncensored=[0.5_dp, 1.2_dp, 2.5_dp], &
        left=[0.3_dp], right=[3.5_dp], interval=interval)
    call close(gamma_censored_loglikelihood(data, 0.6_dp, scale=2.3_dp), &
        -15.04060700971328_dp, 3.0e-12_dp, 'gamma mixed reference')
    call close(gamma_censored_nnlf(data, 0.6_dp, scale=2.3_dp), &
        15.04060700971328_dp, 3.0e-12_dp, 'gamma negative log likelihood')
    call close(gamma_censored_loglikelihood(data, 1.0_dp, scale=2.3_dp), &
        exponential_censored_loglikelihood(data, scale=2.3_dp), 0.0_dp, &
        'gamma shape-one exponential identity')
    call gamma_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], &
        fit, guess=[1.0_dp, 0.0_dp, 2.0_dp], tolerance=1.0e-9_dp)
    call check(fit%success, 'gamma mixed fit')
    call close(fit%params(1), 0.9927593350620091_dp, 3.0e-4_dp, 'gamma fitted shape')
    call close(fit%params(3), 2.232420387730313_dp, 3.0e-4_dp, 'gamma fitted scale')
    call close(fit%nllf, gamma_censored_nnlf(data, fit%params(1), fit%params(2), &
        fit%params(3)), 1.0e-13_dp, 'gamma fit objective')
    call gamma_fit(data, [0.1_dp, -1.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'gamma free location rejected')
    call gamma_fit(data, [0.1_dp, nan, 0.1_dp], [5.0_dp, nan, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'gamma NaN location bounds rejected')
    call gamma_fit(empty, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'gamma empty fit rejected')
    data = make_censored_data(uncensored=[0.2_dp, 1.2_dp, 3.4_dp])
    call close(gamma_censored_loglikelihood(data, 2.5_dp, loc=-0.2_dp, scale=1.3_dp), &
        gamma_loglikelihood(data%uncensored, 2.5_dp, loc=-0.2_dp, scale=1.3_dp), &
        3.0e-14_dp, 'gamma uncensored array agreement')
    call set_interval(0.0_dp, 1.0e-300_dp)
    call close(gamma_censored_loglikelihood(data, 0.4_dp, scale=1.0e300_dp), &
        -552.5008094043986_dp, 3.0e-12_dp, 'gamma boundary and scaled underflow')
    call set_interval(1.0_dp, nearest(1.0_dp, 1.0_dp))
    call close(gamma_censored_loglikelihood(data, 0.4_dp), &
        -37.84033120681894_dp, 3.0e-13_dp, 'gamma adjacent endpoints')
    call set_interval(1000.0_dp, 1000.0_dp+1.0e-10_dp)
    call close(gamma_censored_loglikelihood(data, 3.5_dp), &
        -1006.956992261547_dp, 5.0e-11_dp, 'gamma narrow extreme tail')
    call set_interval(-2.0_dp, 1.0_dp)
    call close(gamma_censored_loglikelihood(data, 0.4_dp), gamma_logcdf(1.0_dp, 0.4_dp), &
        1.0e-14_dp, 'gamma interval crossing support')
    call set_interval(-2.0_dp, 0.0_dp)
    call check(gamma_censored_loglikelihood(data, 0.4_dp) == -inf, &
        'gamma interval ending at support')
    data = make_censored_data(right=[-1.0_dp, 0.0_dp])
    call close(gamma_censored_loglikelihood(data, 0.4_dp), 0.0_dp, 0.0_dp, &
        'gamma right limits below and at support')
    data = make_censored_data(left=[0.0_dp])
    call check(gamma_censored_loglikelihood(data, 0.4_dp) == -inf, &
        'gamma left limit at support')
    data = make_censored_data(uncensored=[0.0_dp])
    call check(gamma_censored_loglikelihood(data, 0.4_dp) == inf, &
        'gamma infinite density at support for small shape')
    call check(gamma_censored_loglikelihood(data, 2.0_dp) == -inf, &
        'gamma zero density at support for shape above one')
    call close(gamma_censored_loglikelihood(data, 1.0_dp, scale=2.0_dp), &
        -log(2.0_dp), 1.0e-15_dp, 'gamma shape-one support density')
    call gamma_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'gamma exact support endpoint fit rejected')
    data = make_censored_data(uncensored=[0.0_dp], left=[-1.0_dp])
    call check(ieee_is_nan(gamma_censored_loglikelihood(data, 0.4_dp)), &
        'gamma opposing infinite contributions yield NaN')
    data = make_censored_data(uncensored=[-1.0_dp, 0.0_dp])
    call check(ieee_is_nan(gamma_censored_loglikelihood(data, 0.4_dp)), &
        'gamma opposing infinities in reverse order')
    data = make_censored_data(left=[0.5_dp], right=[2.0_dp])
    other = make_censored_data(interval=reshape([-inf, 2.0_dp, 0.5_dp, inf], [2,2]))
    call close(gamma_censored_loglikelihood(data, 0.4_dp), &
        gamma_censored_loglikelihood(other, 0.4_dp), 1.0e-14_dp, &
        'gamma infinite interval endpoint agreement')
    call set_interval(1.0_dp, 1.0_dp)
    other = make_censored_data(uncensored=[1.0_dp])
    call close(gamma_censored_loglikelihood(data, 0.4_dp), &
        gamma_censored_loglikelihood(other, 0.4_dp), 1.0e-14_dp, 'gamma equal endpoints')
    call close(gamma_censored_loglikelihood(empty, 0.4_dp), 0.0_dp, 0.0_dp, &
        'gamma empty likelihood')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, 0.0_dp)), 'gamma zero shape')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, -1.0_dp)), 'gamma negative shape')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, nan)), 'gamma NaN shape')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, inf)), 'gamma infinite shape')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, 1.0_dp, scale=nan)), &
        'gamma NaN scale')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, 1.0_dp, scale=0.0_dp)), &
        'gamma zero scale')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, 1.0_dp, scale=inf)), &
        'gamma infinite scale')
    call check(ieee_is_nan(gamma_censored_loglikelihood(empty, 1.0_dp, loc=nan)), &
        'gamma NaN location')
    data = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(gamma_censored_loglikelihood(data, 0.4_dp)), 'gamma invalid data')

    data = make_censored_data(right=[1.0e-300_dp])
    call close(gamma_censored_loglikelihood(data, 1.0e-20_dp, scale=1.0e300_dp), &
        -38.82115764935087_dp, 3.0e-13_dp, 'gamma small-shape scaled-underflow SF')
    data = make_censored_data(uncensored=[1.0e-300_dp])
    call close(gamma_censored_loglikelihood(data, 0.4_dp, scale=1.0e300_dp), &
        137.35842776194093_dp, 5.0e-13_dp, 'gamma scaled-underflow exact density')
    call set_interval(100.0_dp, 101.0_dp)
    call close(gamma_censored_loglikelihood(data, 2.5_dp), &
        -93.82934250728414_dp, 4.0e-12_dp, 'gamma wider upper-tail interval')
    call set_interval(0.1_dp, 0.2_dp)
    call close(gamma_censored_loglikelihood(data, 0.4_dp), &
        -2.088544002566658_dp, 4.0e-13_dp, 'gamma wider lower-tail interval')
    base = gamma_censored_loglikelihood(data, 0.4_dp)
    call set_interval(0.1_dp, 0.3_dp)
    value = gamma_censored_loglikelihood(data, 0.4_dp)
    call check(value >= base .and. value <= 0.0_dp, 'gamma interval mass monotonicity')
    data = make_censored_data(left=[0.1_dp])
    base = gamma_censored_loglikelihood(data, 0.4_dp)
    data = make_censored_data(left=[1.0_dp])
    value = gamma_censored_loglikelihood(data, 0.4_dp)
    call check(value >= base .and. value <= 0.0_dp, 'gamma lower tail monotonicity')
    data = make_censored_data(right=[0.1_dp])
    base = gamma_censored_loglikelihood(data, 0.4_dp)
    data = make_censored_data(right=[1.0_dp])
    value = gamma_censored_loglikelihood(data, 0.4_dp)
    call check(value <= base .and. value <= 0.0_dp, 'gamma upper tail monotonicity')
    data = make_censored_data(uncensored=[2.0_dp, 4.0_dp, 6.0_dp])
    call gamma_fit(data, [2.0_dp, 0.0_dp, 0.1_dp], [2.0_dp, 0.0_dp, 8.0_dp], fit, &
        tolerance=1.0e-9_dp)
    call check(fit%success, 'gamma fixed-shape exact-only fit')
    call close(fit%params(3), 2.0_dp, 1.0e-7_dp, 'gamma analytic fixed-shape scale')
    call set_interval(0.0_dp, 0.0_dp)
    call gamma_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_invalid_input, 'gamma exact interval at endpoint rejected')
    data = make_censored_data(uncensored=[-1.0_dp])
    call gamma_fit(data, [0.1_dp, 0.0_dp, 0.1_dp], [5.0_dp, 0.0_dp, 8.0_dp], fit)
    call check(fit%status == fit_status_no_finite_objective, 'gamma impossible exact sample fit')

    data = make_censored_data(left=[1.0e-300_dp])
    call check(gamma_censored_loglikelihood(data, 1.0e308_dp, scale=1.0e300_dp) == -inf, &
        'gamma huge-shape scaled-underflow lower tail saturates safely')
    data = make_censored_data(right=[1.0e-300_dp])
    call close(gamma_censored_loglikelihood(data, 1.0e308_dp, scale=1.0e300_dp), &
        0.0_dp, 0.0_dp, 'gamma huge-shape scaled-underflow upper tail')

    data = make_censored_data(uncensored=[0.1_dp,0.3_dp,0.6_dp], left=[0.05_dp], &
        right=[0.8_dp], interval=reshape([0.2_dp,0.7_dp,0.25_dp,0.75_dp],[2,2]))
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        -10.367442718990157_dp, 4.0e-12_dp, 'beta mixed reference')
    call close(beta_censored_nnlf(data,0.4_dp,0.6_dp), &
        10.367442718990157_dp, 4.0e-12_dp, 'beta negative log likelihood')
    call close(beta_censored_loglikelihood(data,1.0_dp,1.0_dp), &
        log(0.05_dp)+log(0.2_dp)+2.0_dp*log(0.05_dp), &
        2.0e-14_dp, 'beta uniform mixed-data identity')
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp], &
        fit,guess=[1.0_dp,1.0_dp,0.0_dp,1.0_dp],tolerance=1.0e-9_dp)
    call check(fit%success,'beta mixed fit')
    call close(fit%params(1),0.640291446923525_dp,3.0e-4_dp,'beta fitted first shape')
    call close(fit%params(2),0.9270186194778561_dp,3.0e-4_dp,'beta fitted second shape')
    call close(fit%nllf,beta_censored_nnlf(data,fit%params(1),fit%params(2), &
        fit%params(3),fit%params(4)),1.0e-13_dp,'beta fit objective')
    call beta_fit(data,[0.1_dp,0.1_dp,-1.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta free location rejected')
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,0.5_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta free scale rejected')
    call beta_fit(data,[0.1_dp,0.1_dp,nan,1.0_dp],[5.0_dp,5.0_dp,nan,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta NaN location bounds rejected')
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,nan],[5.0_dp,5.0_dp,0.0_dp,nan],fit)
    call check(fit%status == fit_status_invalid_input,'beta NaN scale bounds rejected')
    call beta_fit(empty,[0.1_dp,0.1_dp,0.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta empty fit rejected')
    data = make_censored_data(uncensored=[0.2_dp,0.4_dp,0.7_dp])
    call close(beta_censored_loglikelihood(data,2.5_dp,3.5_dp,loc=-0.5_dp,scale=2.0_dp), &
        beta_loglikelihood(data%uncensored,2.5_dp,3.5_dp,loc=-0.5_dp,scale=2.0_dp), &
        4.0e-14_dp,'beta uncensored array agreement')
    call close(beta_censored_loglikelihood(data,1.0_dp,1.0_dp,loc=-0.5_dp,scale=2.0_dp), &
        uniform_loglikelihood(data%uncensored,loc=-0.5_dp,scale=2.0_dp), &
        1.0e-14_dp,'beta uniform exact-data identity')
    call set_interval(0.5_dp,nearest(0.5_dp,1.0_dp))
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        -37.23856506488818_dp,4.0e-13_dp,'beta adjacent interior endpoints')
    call set_interval(0.0_dp,1.0e-300_dp)
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        -276.58883210318237_dp,3.0e-12_dp,'beta tiny interval touching lower edge')
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp,scale=1.0e300_dp), &
        -552.8990432624679_dp,4.0e-12_dp,'beta standardized lower distance underflows')
    call set_interval(nearest(1.0_dp,-1.0_dp),1.0_dp)
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        -22.726166393811287_dp,4.0e-13_dp,'beta interval touching upper edge')
    call set_interval(0.99_dp,0.999_dp)
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        -3.7335165416256575_dp,4.0e-13_dp,'beta upper tail interval')
    data = make_censored_data(right=[1.0e300_dp])
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp, &
        loc=1.0e-300_dp,scale=1.0e300_dp),-829.6147195298614_dp,5.0e-12_dp, &
        'beta compensated upper edge retains positive upper tail')
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp, &
        loc=-1.0e-300_dp,scale=1.0e300_dp) == -inf, &
        'beta compensated upper edge recognizes outside point')
    data = make_censored_data(left=[0.2_dp],right=[0.8_dp])
    other = make_censored_data(interval=reshape([-inf,0.8_dp,0.2_dp,inf],[2,2]))
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        beta_censored_loglikelihood(other,0.4_dp,0.6_dp), &
        1.0e-14_dp,'beta infinite interval endpoint agreement')
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        beta_logcdf(0.2_dp,0.4_dp,0.6_dp)+beta_logsf(0.8_dp,0.4_dp,0.6_dp), &
        3.0e-14_dp,'beta ordinary log-tail agreement')
    call set_interval(0.4_dp,0.4_dp)
    other = make_censored_data(uncensored=[0.4_dp])
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        beta_censored_loglikelihood(other,0.4_dp,0.6_dp), &
        1.0e-14_dp,'beta equal endpoints represent exact observations')
    call set_interval(-2.0_dp,2.0_dp)
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp),0.0_dp,0.0_dp, &
        'beta interval covers support')
    call set_interval(-2.0_dp,0.0_dp)
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp) == -inf, &
        'beta interval outside lower support')
    call set_interval(1.0_dp,2.0_dp)
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp) == -inf, &
        'beta interval outside upper support')
    data = make_censored_data(uncensored=[0.0_dp])
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp) == inf, &
        'beta infinite lower endpoint density')
    call check(beta_censored_loglikelihood(data,2.0_dp,0.6_dp) == -inf, &
        'beta zero lower endpoint density')
    call close(beta_censored_loglikelihood(data,1.0_dp,2.0_dp),log(2.0_dp), &
        1.0e-15_dp,'beta finite lower endpoint density')
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta exact lower endpoint fit rejected')
    data = make_censored_data(uncensored=[1.0_dp])
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp) == inf, &
        'beta infinite upper endpoint density')
    call check(beta_censored_loglikelihood(data,0.4_dp,2.0_dp) == -inf, &
        'beta zero upper endpoint density')
    call close(beta_censored_loglikelihood(data,2.0_dp,1.0_dp),log(2.0_dp), &
        1.0e-15_dp,'beta finite upper endpoint density')
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta exact upper endpoint fit rejected')
    call set_interval(1.0_dp,1.0_dp)
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'beta exact interval endpoint fit rejected')
    data = make_censored_data(uncensored=[0.0_dp,1.0_dp])
    call check(ieee_is_nan(beta_censored_loglikelihood(data,0.4_dp,2.0_dp)), &
        'beta opposing infinite densities yield NaN')
    data = make_censored_data(uncensored=[-1.0_dp])
    call beta_fit(data,[0.1_dp,0.1_dp,0.0_dp,1.0_dp],[5.0_dp,5.0_dp,0.0_dp,1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'beta impossible exact sample fit')
    data = make_censored_data(right=[0.0_dp],left=[1.0_dp])
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp),0.0_dp,0.0_dp, &
        'beta certain left and right limits at endpoints')
    data = make_censored_data(left=[0.0_dp])
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp) == -inf, &
        'beta impossible left limit at lower edge')
    data = make_censored_data(right=[1.0_dp])
    call check(beta_censored_loglikelihood(data,0.4_dp,0.6_dp) == -inf, &
        'beta impossible right limit at upper edge')
    call set_interval(0.1_dp,0.2_dp)
    base = beta_censored_loglikelihood(data,0.4_dp,0.6_dp)
    call set_interval(0.1_dp,0.3_dp)
    value = beta_censored_loglikelihood(data,0.4_dp,0.6_dp)
    call check(value >= base .and. value <= 0.0_dp,'beta interval monotonicity and nonnegativity')
    data = make_censored_data(uncensored=[0.2_dp,0.6_dp],left=[0.1_dp],right=[0.9_dp])
    other = make_censored_data(uncensored=[0.8_dp,0.4_dp],left=[0.1_dp],right=[0.9_dp])
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp), &
        beta_censored_loglikelihood(other,0.6_dp,0.4_dp),3.0e-13_dp,'beta reflection identity')
    call close(beta_censored_loglikelihood(empty,0.4_dp,0.6_dp),0.0_dp,0.0_dp, &
        'beta empty likelihood')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.0_dp,0.6_dp)), 'beta zero first shape')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.4_dp,-1.0_dp)), 'beta negative second shape')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,nan,0.6_dp)), 'beta NaN first shape')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.4_dp,inf)), 'beta infinite second shape')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.4_dp,0.6_dp,loc=nan)), 'beta NaN location')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.4_dp,0.6_dp,scale=nan)), 'beta NaN scale')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.4_dp,0.6_dp,scale=0.0_dp)), 'beta zero scale')
    call check(ieee_is_nan(beta_censored_loglikelihood(empty,0.4_dp,0.6_dp, &
        loc=huge(1.0_dp),scale=huge(1.0_dp))), 'beta unrepresentable upper edge')
    data = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(beta_censored_loglikelihood(data,0.4_dp,0.6_dp)), 'beta invalid sample')

    data = make_censored_data(right=[1.0e-300_dp])
    call close(beta_censored_loglikelihood(data,1.0e-20_dp,0.6_dp,scale=1.0e300_dp), &
        -38.82004266769975_dp,4.0e-13_dp,'beta tiny shape scaled-underflow upper tail')
    data = make_censored_data(uncensored=[1.0e-300_dp])
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp,scale=1.0e300_dp), &
        136.9601939038717_dp,5.0e-13_dp,'beta scaled-underflow exact density')
    data = make_censored_data(uncensored=[0.5_dp])
    call close(beta_censored_loglikelihood(data,1000.0_dp,1000.0_dp), &
        3.574534877131522_dp,3.0e-12_dp,'beta large-shape density at mode')
    call set_interval(0.5_dp,nearest(0.5_dp,1.0_dp))
    call close(beta_censored_loglikelihood(data,1000.0_dp,1000.0_dp), &
        -33.16226569254558_dp,4.0e-12_dp,'beta large-shape adjacent interval')
    data = make_censored_data(right=[-1.0e-300_dp])
    call close(beta_censored_loglikelihood(data,0.4_dp,0.6_dp, &
        loc=-1.0e300_dp,scale=1.0e300_dp),-829.6147195298614_dp,5.0e-12_dp, &
        'beta independently normalized upper distance underflows')
    data = make_censored_data(left=[1.0e-300_dp])
    call check(ieee_is_nan(beta_censored_loglikelihood(data,0.4_dp,1.0e308_dp, &
        scale=1.0e24_dp)),'beta unsupported underflow expansion returns NaN safely')

    data = make_censored_data(right=[1.0e-300_dp])
    call close(beta_censored_loglikelihood(data,1.0e-20_dp,1.0e308_dp,scale=1.0e300_dp), &
        -39.541774479640864_dp,4.0e-13_dp,'beta tiny/huge shapes retain small complement')

    call close(beta_censored_loglikelihood(data,1.0e-20_dp,1.0e6_dp,scale=1.0e300_dp), &
        -38.831212206840114_dp,4.0e-13_dp,'beta gamma-ratio powers do not overflow')

    ! Logistic references: independent mpmath formulas and SciPy CensoredData fit.
    data = make_censored_data(uncensored=[0.5_dp,1.2_dp,2.5_dp],left=[0.3_dp], &
        right=[3.5_dp],interval=reshape([1.1_dp,4.0_dp,1.4_dp,4.4_dp],[2,2]))
    call close(logistic_censored_loglikelihood(data,0.4_dp,1.3_dp), &
        -15.999635359473038_dp,2.0e-12_dp,'logistic mixed reference')
    call close(logistic_censored_nnlf(data,0.4_dp,1.3_dp), &
        -logistic_censored_loglikelihood(data,0.4_dp,1.3_dp),0.0_dp,'logistic NNLF')
    call logistic_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit,tolerance=1.0e-10_dp)
    call check(fit%success,'logistic censored fit converges')
    call close(fit%params(1),1.847711699978416_dp,3.0e-4_dp,'logistic SciPy fit center')
    call close(fit%params(2),1.12816063694188_dp,3.0e-4_dp,'logistic SciPy fit scale')
    call close(fit%nllf,logistic_censored_nnlf(data,fit%params(1),fit%params(2)), &
        1.0e-13_dp,'logistic fit objective')
    call logistic_fit(data,[0.4_dp,1.3_dp],[0.4_dp,1.3_dp],fit)
    call check(fit%success,'logistic both parameters fixed')
    call close(fit%nllf,logistic_censored_nnlf(data,0.4_dp,1.3_dp), &
        0.0_dp,'logistic fixed objective')
    call logistic_fit(data,[0.4_dp,0.05_dp],[0.4_dp,10.0_dp],fit)
    call check(fit%success .and. fit%params(1) == 0.4_dp,'logistic fixed location')
    other = make_censored_data(uncensored=[-0.5_dp,-1.2_dp,-2.5_dp],right=[-0.3_dp], &
        left=[-3.5_dp],interval=reshape([-1.4_dp,-4.4_dp,-1.1_dp,-4.0_dp],[2,2]))
    call close(logistic_censored_loglikelihood(data,0.4_dp,1.3_dp), &
        logistic_censored_loglikelihood(other,-0.4_dp,1.3_dp),2.0e-13_dp,'logistic reflection')
    data = make_censored_data(uncensored=[-2.0_dp,0.0_dp,0.5_dp,4.0_dp])
    call close(logistic_censored_loglikelihood(data,0.4_dp,1.3_dp), &
        logistic_loglikelihood(data%uncensored,0.4_dp,1.3_dp),2.0e-14_dp, &
        'logistic exact array agreement')
    call logistic_fit(data%uncensored,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    base = fit%nllf
    call logistic_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call close(fit%nllf,base,1.0e-12_dp,'logistic exact fit overload agreement')
    call set_interval(-1.0_dp,1.0_dp)
    base = logistic_censored_loglikelihood(data)
    call close(base,log(tanh(0.5_dp)),2.0e-15_dp,'logistic symmetric mass identity')
    call set_interval(-2.0_dp,2.0_dp)
    value = logistic_censored_loglikelihood(data)
    call check(value >= base .and. value <= 0.0_dp,'logistic nested interval masses')
    call set_interval(0.0_dp,nearest(0.0_dp,1.0_dp))
    call close(logistic_censored_loglikelihood(data), &
        log(nearest(0.0_dp,1.0_dp))-log(4.0_dp),2.0e-13_dp,'logistic subnormal width')
    call set_interval(1.0e-300_dp,2.0e-300_dp)
    call close(logistic_censored_loglikelihood(data,scale=1.0e300_dp), &
        log(1.0e-300_dp)-log(1.0e300_dp)-log(4.0_dp),3.0e-13_dp, &
        'logistic standardized width underflows')
    call set_interval(1000.0_dp,nearest(1000.0_dp,1.0_dp))
    call close(logistic_censored_loglikelihood(data),-1029.8053287640778_dp, &
        3.0e-13_dp,'logistic adjacent upper tail')
    call set_interval(-1000.0_dp,nearest(-1000.0_dp,1.0_dp))
    call close(logistic_censored_loglikelihood(data),-1029.8053287640778_dp, &
        3.0e-13_dp,'logistic adjacent lower tail')
    call set_interval(-1.0e308_dp,1.0e308_dp)
    call close(logistic_censored_loglikelihood(data),0.0_dp,0.0_dp,'logistic width overflow')
    data = make_censored_data(uncensored=[1.0e308_dp])
    call close(logistic_censored_loglikelihood(data,loc=-1.0e308_dp,scale=1.0e308_dp), &
        -2.0_dp-2.0_dp*log(1.0_dp+exp(-2.0_dp))-log(1.0e308_dp), &
        2.0e-13_dp,'logistic subtraction overflow with finite standardized value')
    call check(logistic_censored_loglikelihood(data,scale=1.0e-300_dp) == -inf, &
        'logistic unrepresentable log density')
    data = make_censored_data(uncensored=[nearest(huge(1.0_dp),-1.0_dp)])
    call close(logistic_censored_loglikelihood(data, &
        loc=-1.5_dp*spacing(huge(1.0_dp)),scale=huge(1.0_dp)), &
        -1.0_dp-2.0_dp*log(1.0_dp+exp(-1.0_dp))-log(huge(1.0_dp)), &
        2.0e-13_dp,'logistic rounded subtraction-overflow threshold')
    data = make_censored_data(uncensored=[0.75_dp*huge(1.0_dp)])
    call close(logistic_censored_loglikelihood(data,scale=0.75_dp), &
        -huge(1.0_dp),0.0_dp,'logistic largest finite standardized magnitude')
    data = make_censored_data(uncensored=[0.0_dp])
    call close(logistic_censored_loglikelihood(data,scale=nearest(0.0_dp,1.0_dp)), &
        -log(4.0_dp)-log(nearest(0.0_dp,1.0_dp)),2.0e-13_dp, &
        'logistic zero displacement with smallest positive scale')
    data = make_censored_data(left=[-1.0e308_dp])
    call check(logistic_censored_loglikelihood(data,scale=1.0e-300_dp) == -inf, &
        'logistic unrepresentable lower log tail')
    data = make_censored_data(right=[1.0e308_dp])
    call check(logistic_censored_loglikelihood(data,scale=1.0e-300_dp) == -inf, &
        'logistic unrepresentable upper log tail')
    data = make_censored_data(left=[-1000.0_dp],right=[1000.0_dp])
    call close(logistic_censored_loglikelihood(data),-2000.0_dp,0.0_dp,'logistic extreme tails')
    call set_interval(-inf,0.4_dp)
    call close(logistic_censored_loglikelihood(data,0.4_dp,1.3_dp), &
        logistic_logcdf(0.4_dp,0.4_dp,1.3_dp),0.0_dp,'logistic infinite lower endpoint')
    call set_interval(0.4_dp,inf)
    call close(logistic_censored_loglikelihood(data,0.4_dp,1.3_dp), &
        logistic_logsf(0.4_dp,0.4_dp,1.3_dp),0.0_dp,'logistic infinite upper endpoint')
    call set_interval(0.4_dp,0.4_dp)
    other = make_censored_data(uncensored=[0.4_dp])
    call close(logistic_censored_loglikelihood(data,0.4_dp,1.3_dp), &
        logistic_censored_loglikelihood(other,0.4_dp,1.3_dp),0.0_dp,'logistic equal endpoints exact')
    call close(logistic_censored_loglikelihood(empty),0.0_dp,0.0_dp,'logistic empty likelihood')
    call check(ieee_is_nan(logistic_censored_loglikelihood(empty,loc=nan)),'logistic NaN center')
    call check(ieee_is_nan(logistic_censored_loglikelihood(empty,loc=inf)),'logistic infinite center')
    call check(ieee_is_nan(logistic_censored_loglikelihood(empty,scale=nan)),'logistic NaN scale')
    call check(ieee_is_nan(logistic_censored_loglikelihood(empty,scale=inf)),'logistic infinite scale')
    call check(ieee_is_nan(logistic_censored_loglikelihood(empty,scale=0.0_dp)),'logistic zero scale')
    call check(ieee_is_nan(logistic_censored_loglikelihood(empty,scale=-1.0_dp)),'logistic negative scale')
    data = make_censored_data(uncensored=[nan])
    call check(ieee_is_nan(logistic_censored_loglikelihood(data)),'logistic invalid exact observation')
    call logistic_fit(data,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'logistic invalid sample fit')
    call logistic_fit(empty,[-5.0_dp,0.05_dp],[8.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'logistic empty sample fit')
    data = make_censored_data(uncensored=[0.0_dp])
    call logistic_fit(data,[1.0_dp,0.05_dp],[-1.0_dp,10.0_dp],fit)
    call check(fit%status == fit_status_invalid_input,'logistic reversed bounds')
    call logistic_fit(data,[0.0_dp,-1.0_dp],[0.0_dp,-1.0_dp],fit)
    call check(fit%status == fit_status_no_finite_objective,'logistic invalid fixed scale')

    if (failures /= 0) error stop 1
    print '(a,i0,a)', 'Censored likelihood and fitting: ', checks, ' checks passed'

contains

    pure function custom_pdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! normal mean and standard deviation
        real(dp) :: value
        value = normal_logpdf(x, params(1), params(2))
    end function custom_pdf

    pure function custom_cdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! normal mean and standard deviation
        real(dp) :: value
        value = normal_logcdf(x, params(1), params(2))
    end function custom_cdf

    pure function custom_sf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! normal mean and standard deviation
        real(dp) :: value
        value = normal_logsf(x, params(1), params(2))
    end function custom_sf

    subroutine set_interval(a, b)
        implicit none
        real(dp), intent(in) :: a !! test lower endpoint
        real(dp), intent(in) :: b !! test upper endpoint
        data = make_censored_data(interval=reshape([a,b], [1,2]))
    end subroutine set_interval

    subroutine check(ok, label)
        implicit none
        logical, intent(in) :: ok !! expected invariant
        character(*), intent(in) :: label !! failure description
        checks = checks + 1
        if (ok) return
        failures = failures + 1
        print '(a)', 'FAIL: '//label
    end subroutine check

    subroutine close(actual, expected, atol, label)
        implicit none
        real(dp), intent(in) :: actual !! computed finite value
        real(dp), intent(in) :: expected !! reference finite value
        real(dp), intent(in) :: atol !! regime-specific absolute tolerance
        character(*), intent(in) :: label !! failure description
        call check(ieee_is_finite(actual), label//' is finite')
        if (.not. ieee_is_finite(actual)) return
        call check(abs(actual - expected) <= atol, label)
        if (abs(actual - expected) > atol) print *, actual, expected
    end subroutine close

end program test_censored
