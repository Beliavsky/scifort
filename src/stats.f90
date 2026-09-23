! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_stats
    use scifort_arcsine, only : arcsine_cdf, arcsine_isf, arcsine_logcdf, &
        arcsine_logpdf, arcsine_logsf, arcsine_pdf, arcsine_ppf, arcsine_sf
    use scifort_halfnorm, only : halfnorm_cdf, halfnorm_isf, halfnorm_logcdf, &
        halfnorm_logpdf, halfnorm_logsf, halfnorm_pdf, halfnorm_ppf, halfnorm_sf
    use scifort_halfcauchy, only : halfcauchy_cdf, halfcauchy_isf, halfcauchy_logcdf, &
        halfcauchy_logpdf, halfcauchy_logsf, halfcauchy_pdf, halfcauchy_ppf, halfcauchy_sf
    use scifort_lomax, only : lomax_cdf, lomax_isf, lomax_logcdf, lomax_logpdf, &
        lomax_logsf, lomax_pdf, lomax_ppf, lomax_sf
    use scifort_chi, only : chi_cdf, chi_isf, chi_logcdf, chi_logpdf, &
        chi_logsf, chi_pdf, chi_ppf, chi_sf
    use scifort_maxwell, only : maxwell_cdf, maxwell_isf, maxwell_logcdf, &
        maxwell_logpdf, maxwell_logsf, maxwell_pdf, maxwell_ppf, maxwell_sf
    use scifort_cosine, only : cosine_cdf, cosine_isf, cosine_logcdf, cosine_logpdf, &
        cosine_logsf, cosine_pdf, cosine_ppf, cosine_sf
    use scifort_semicircular, only : semicircular_cdf, semicircular_isf, &
        semicircular_logcdf, semicircular_logpdf, semicircular_logsf, &
        semicircular_pdf, semicircular_ppf, semicircular_sf
    use scifort_anglit, only : anglit_cdf, anglit_isf, anglit_logcdf, anglit_logpdf, &
        anglit_logsf, anglit_pdf, anglit_ppf, anglit_sf
    use scifort_moyal, only : moyal_cdf, moyal_isf, moyal_logcdf, moyal_logpdf, &
        moyal_logsf, moyal_pdf, moyal_ppf, moyal_sf
    use scifort_hypsecant, only : hypsecant_cdf, hypsecant_isf, hypsecant_logcdf, &
        hypsecant_logpdf, hypsecant_logsf, hypsecant_pdf, hypsecant_ppf, hypsecant_sf
    use scifort_halflogistic, only : halflogistic_cdf, halflogistic_isf, &
        halflogistic_logcdf, halflogistic_logpdf, halflogistic_logsf, &
        halflogistic_pdf, halflogistic_ppf, halflogistic_sf
    use scifort_invgamma, only : invgamma_cdf, invgamma_isf, invgamma_logcdf, &
        invgamma_logpdf, invgamma_logsf, invgamma_pdf, invgamma_ppf, invgamma_sf
    use scifort_invgauss, only : invgauss_cdf, invgauss_isf, invgauss_logcdf, &
        invgauss_logpdf, invgauss_logsf, invgauss_pdf, invgauss_ppf, invgauss_sf
    use scifort_levy, only : levy_cdf, levy_isf, levy_logcdf, levy_logpdf, &
        levy_logsf, levy_pdf, levy_ppf, levy_sf
    use scifort_loglaplace, only : loglaplace_cdf, loglaplace_isf, loglaplace_logcdf, &
        loglaplace_logpdf, loglaplace_logsf, loglaplace_pdf, loglaplace_ppf, loglaplace_sf
    use scifort_bradford, only : bradford_cdf, bradford_isf, bradford_logcdf, &
        bradford_logpdf, bradford_logsf, bradford_pdf, bradford_ppf, bradford_sf
    use scifort_truncexpon, only : truncexpon_cdf, truncexpon_isf, truncexpon_logcdf, &
        truncexpon_logpdf, truncexpon_logsf, truncexpon_pdf, truncexpon_ppf, truncexpon_sf
    use scifort_fisk, only : fisk_cdf, fisk_isf, fisk_logcdf, fisk_logpdf, &
        fisk_logsf, fisk_pdf, fisk_ppf, fisk_sf
    use scifort_dweibull, only : dweibull_cdf, dweibull_isf, dweibull_logcdf, &
        dweibull_logpdf, dweibull_logsf, dweibull_pdf, dweibull_ppf, dweibull_sf
    use scifort_alpha, only : alpha_cdf, alpha_isf, alpha_logcdf, alpha_logpdf, &
        alpha_logsf, alpha_pdf, alpha_ppf, alpha_sf
    use scifort_fatiguelife, only : fatiguelife_cdf, fatiguelife_isf, &
        fatiguelife_logcdf, fatiguelife_logpdf, fatiguelife_logsf, fatiguelife_pdf, &
        fatiguelife_ppf, fatiguelife_sf
    use scifort_genlogistic, only : genlogistic_cdf, genlogistic_isf, &
        genlogistic_logcdf, genlogistic_logpdf, genlogistic_logsf, genlogistic_pdf, &
        genlogistic_ppf, genlogistic_sf
    use scifort_gennorm, only : gennorm_cdf, gennorm_isf, gennorm_logcdf, gennorm_logpdf, &
        gennorm_logsf, gennorm_pdf, gennorm_ppf, gennorm_sf
    use scifort_nakagami, only : nakagami_cdf, nakagami_isf, nakagami_logcdf, &
        nakagami_logpdf, nakagami_logsf, nakagami_pdf, nakagami_ppf, nakagami_sf
    use scifort_powernorm, only : powernorm_cdf, powernorm_isf, powernorm_logcdf, &
        powernorm_logpdf, powernorm_logsf, powernorm_pdf, powernorm_ppf, powernorm_sf
    use scifort_loggamma, only : loggamma_cdf, loggamma_isf, loggamma_logcdf, &
        loggamma_logpdf, loggamma_logsf, loggamma_pdf, loggamma_ppf, loggamma_sf
    use scifort_wald, only : wald_cdf, wald_isf, wald_logcdf, wald_logpdf, &
        wald_logsf, wald_pdf, wald_ppf, wald_sf
    use scifort_gompertz, only : gompertz_cdf, gompertz_isf, gompertz_logcdf, &
        gompertz_logpdf, gompertz_logsf, gompertz_pdf, gompertz_ppf, gompertz_sf
    use scifort_invweibull, only : invweibull_cdf, invweibull_isf, invweibull_logcdf, &
        invweibull_logpdf, invweibull_logsf, invweibull_pdf, invweibull_ppf, invweibull_sf
    use scifort_betaprime, only : betaprime_cdf, betaprime_isf, betaprime_logcdf, &
        betaprime_logpdf, betaprime_logsf, betaprime_pdf, betaprime_ppf, betaprime_sf
    use scifort_burr12, only : burr12_cdf, burr12_isf, burr12_logcdf, burr12_logpdf, &
        burr12_logsf, burr12_pdf, burr12_ppf, burr12_sf
    use scifort_genhalflogistic, only : genhalflogistic_cdf, genhalflogistic_isf, &
        genhalflogistic_logcdf, genhalflogistic_logpdf, genhalflogistic_logsf, &
        genhalflogistic_pdf, genhalflogistic_ppf, genhalflogistic_sf
    use scifort_exponpow, only : exponpow_cdf, exponpow_isf, exponpow_logcdf, &
        exponpow_logpdf, exponpow_logsf, exponpow_pdf, exponpow_ppf, exponpow_sf
    use scifort_exponweib, only : exponweib_cdf, exponweib_isf, exponweib_logcdf, &
        exponweib_logpdf, exponweib_logsf, exponweib_pdf, exponweib_ppf, exponweib_sf
    use scifort_powerlognorm, only : powerlognorm_cdf, powerlognorm_isf, &
        powerlognorm_logcdf, powerlognorm_logpdf, powerlognorm_logsf, powerlognorm_pdf, &
        powerlognorm_ppf, powerlognorm_sf
    use scifort_levy_l, only : levy_l_cdf, levy_l_isf, levy_l_logcdf, levy_l_logpdf, &
        levy_l_logsf, levy_l_pdf, levy_l_ppf, levy_l_sf
    use scifort_weibull_max, only : weibull_max_cdf, weibull_max_isf, &
        weibull_max_logcdf, weibull_max_logpdf, weibull_max_logsf, weibull_max_pdf, &
        weibull_max_ppf, weibull_max_sf
    use scifort_rdist, only : rdist_cdf, rdist_isf, rdist_logcdf, rdist_logpdf, &
        rdist_logsf, rdist_pdf, rdist_ppf, rdist_sf
    use scifort_skewcauchy, only : skewcauchy_cdf, skewcauchy_isf, &
        skewcauchy_logcdf, skewcauchy_logpdf, skewcauchy_logsf, skewcauchy_pdf, &
        skewcauchy_ppf, skewcauchy_sf
    use scifort_dgamma, only : dgamma_cdf, dgamma_isf, dgamma_logcdf, dgamma_logpdf, &
        dgamma_logsf, dgamma_pdf, dgamma_ppf, dgamma_sf
    use scifort_laplace_asymmetric, only : laplace_asymmetric_cdf, laplace_asymmetric_isf, &
        laplace_asymmetric_logcdf, laplace_asymmetric_logpdf, laplace_asymmetric_logsf, &
        laplace_asymmetric_pdf, laplace_asymmetric_ppf, laplace_asymmetric_sf
    use scifort_truncnorm, only : truncnorm_cdf, truncnorm_isf, truncnorm_logcdf, &
        truncnorm_logpdf, truncnorm_logsf, truncnorm_pdf, truncnorm_ppf, truncnorm_sf
    use scifort_loguniform, only : loguniform_cdf, loguniform_isf, loguniform_logcdf, &
        loguniform_logpdf, loguniform_logsf, loguniform_pdf, loguniform_ppf, loguniform_sf
    use scifort_bernoulli, only : bernoulli_cdf, bernoulli_isf, bernoulli_logcdf, &
        bernoulli_logpmf, bernoulli_logsf, bernoulli_pmf, bernoulli_ppf, bernoulli_sf
    use scifort_beta, only : beta_cdf, beta_isf, beta_logcdf, beta_logpdf, &
        beta_logsf, beta_pdf, beta_ppf, beta_sf
    use scifort_binomial, only : binomial_cdf, binomial_isf, binomial_logcdf, binomial_logpmf, &
        binomial_logsf, binomial_pmf, binomial_ppf, binomial_sf
    use scifort_cauchy, only : cauchy_cdf, cauchy_isf, cauchy_logcdf, &
        cauchy_logpdf, cauchy_logsf, cauchy_pdf, cauchy_ppf, cauchy_sf
    use scifort_descriptive, only : central_moment, covariance, mean, median, &
        pearson_correlation, quantile, rankdata, standard_deviation, variance
    use scifort_chi2, only : chi2_cdf, chi2_isf, chi2_logcdf, chi2_logpdf, &
        chi2_logsf, chi2_pdf, chi2_ppf, chi2_sf
    use scifort_exponential, only : exponential_cdf, exponential_isf, &
        exponential_logcdf, exponential_logpdf, exponential_logsf, &
        exponential_pdf, exponential_ppf, exponential_sf
    use scifort_fit, only : fit_result, fit_status_success, fit_status_max_iter, &
        fit_status_invalid_input, fit_status_no_finite_objective, normal_fit, &
        uniform_fit, exponential_fit, laplace_fit, logistic_fit, cauchy_fit, &
        rayleigh_fit, gamma_fit, chi2_fit, t_fit, lognormal_fit, weibull_fit, &
        pareto_fit, beta_fit, f_fit, bernoulli_fit, poisson_fit, geometric_fit, &
        binomial_fit, negative_binomial_fit, gumbel_r_fit, gumbel_l_fit, &
        powerlaw_fit, triang_fit, genpareto_fit, arcsine_fit, halfnorm_fit, &
        halfcauchy_fit, lomax_fit, chi_fit, maxwell_fit, cosine_fit, semicircular_fit, &
        anglit_fit, moyal_fit, hypsecant_fit, halflogistic_fit, &
        invgamma_fit, invgauss_fit, levy_fit, loglaplace_fit, &
        bradford_fit, truncexpon_fit, fisk_fit, dweibull_fit, alpha_fit, &
        fatiguelife_fit, genlogistic_fit, gennorm_fit, nakagami_fit, powernorm_fit, &
        loggamma_fit, wald_fit, gompertz_fit, invweibull_fit, betaprime_fit, burr12_fit, &
        genhalflogistic_fit, exponpow_fit, exponweib_fit, powerlognorm_fit, &
        levy_l_fit, weibull_max_fit, rdist_fit, skewcauchy_fit, dgamma_fit, &
        laplace_asymmetric_fit, truncnorm_fit, loguniform_fit
    use scifort_f_distribution, only : f_cdf, f_isf, f_logcdf, f_logpdf, &
        f_logsf, f_pdf, f_ppf, f_sf
    use scifort_gamma, only : gamma_cdf, gamma_isf, gamma_logcdf, gamma_logpdf, &
        gamma_logsf, gamma_pdf, gamma_ppf, gamma_sf
    use scifort_gumbel_r, only : gumbel_r_cdf, gumbel_r_isf, gumbel_r_logcdf, &
        gumbel_r_logpdf, gumbel_r_logsf, gumbel_r_pdf, gumbel_r_ppf, gumbel_r_sf
    use scifort_gumbel_l, only : gumbel_l_cdf, gumbel_l_isf, gumbel_l_logcdf, &
        gumbel_l_logpdf, gumbel_l_logsf, gumbel_l_pdf, gumbel_l_ppf, gumbel_l_sf
    use scifort_powerlaw, only : powerlaw_cdf, powerlaw_isf, powerlaw_logcdf, &
        powerlaw_logpdf, powerlaw_logsf, powerlaw_pdf, powerlaw_ppf, powerlaw_sf
    use scifort_triang, only : triang_cdf, triang_isf, triang_logcdf, triang_logpdf, &
        triang_logsf, triang_pdf, triang_ppf, triang_sf
    use scifort_genpareto, only : genpareto_cdf, genpareto_isf, genpareto_logcdf, &
        genpareto_logpdf, genpareto_logsf, genpareto_pdf, genpareto_ppf, genpareto_sf
    use scifort_geometric, only : geometric_cdf, geometric_isf, geometric_logcdf, &
        geometric_logpmf, geometric_logsf, geometric_pmf, geometric_ppf, geometric_sf
    use scifort_hypothesis, only : correlation_test_result, mannwhitneyu, &
        mannwhitneyu_result, pearsonr, spearmanr, ttest_1samp, ttest_ind, &
        ttest_rel, ttest_result
    use scifort_hypothesis_extended, only : bartlett, bartlett_result, &
        chi2_contingency, chi2_contingency_result, chisquare, contingency_expected_freq, &
        fisher_exact, fisher_exact_result, fligner, fligner_result, friedmanchisquare, &
        friedmanchisquare_result, f_oneway, f_oneway_result, kruskal, kruskal_result, &
        levene, levene_result, make_sample_group, power_divergence, power_divergence_named, &
        power_divergence_result, ranksums, ranksums_result, sample_group, wilcoxon, wilcoxon_result
    use scifort_goodness_of_fit, only : anderson, anderson_result, cramervonmises, &
        cramervonmises_2samp, jarque_bera, ks_1samp, ks_2samp, kstest, kstest_result, &
        kurtosistest, normaltest, shapiro, shapiro_result, significance_result, skewtest
    use scifort_nonparametric_extended, only : anderson_ksamp, anderson_ksamp_result, &
        ansari, epps_singleton_2samp, median_test, median_test_named, median_test_result, mood
    use scifort_association_extended, only : association_result, brunnermunzel, kendalltau, &
        linregress, linregress_result, page_trend_result, page_trend_test, pointbiserialr, &
        siegelslopes, siegelslopes_result, theilslopes, theilslopes_result
    use scifort_contingency_meta, only : association, barnard_exact, binomtest, &
        binomtest_result, boschloo_exact, combine_pvalues, combined_pvalue_result, &
        confidence_interval, exact_2x2_result, false_discovery_control, margins, odds_ratio, &
        odds_ratio_result, relative_risk, relative_risk_result
    use scifort_multiple_comparisons, only : alexandergovern, alexandergovern_result, &
        dunnett, dunnett_result, matrix_confidence_interval, poisson_means_test, &
        poisson_means_test_result, tukey_hsd, tukey_hsd_result, vector_confidence_interval
    use scifort_kinds, only : dp
    use scifort_laplace, only : laplace_cdf, laplace_isf, laplace_logcdf, &
        laplace_logpdf, laplace_logsf, laplace_pdf, laplace_ppf, laplace_sf
    use scifort_likelihood, only : normal_loglikelihood, normal_nnlf, &
        uniform_loglikelihood, uniform_nnlf, exponential_loglikelihood, exponential_nnlf, &
        laplace_loglikelihood, laplace_nnlf, logistic_loglikelihood, logistic_nnlf, &
        cauchy_loglikelihood, cauchy_nnlf, rayleigh_loglikelihood, rayleigh_nnlf, &
        gamma_loglikelihood, gamma_nnlf, chi2_loglikelihood, chi2_nnlf, &
        t_loglikelihood, t_nnlf, lognormal_loglikelihood, lognormal_nnlf, &
        weibull_loglikelihood, weibull_nnlf, pareto_loglikelihood, pareto_nnlf, &
        beta_loglikelihood, beta_nnlf, f_loglikelihood, f_nnlf, &
        bernoulli_loglikelihood, bernoulli_nnlf, poisson_loglikelihood, poisson_nnlf, &
        geometric_loglikelihood, geometric_nnlf, binomial_loglikelihood, binomial_nnlf, &
        negative_binomial_loglikelihood, negative_binomial_nnlf, &
        gumbel_r_loglikelihood, gumbel_r_nnlf, gumbel_l_loglikelihood, gumbel_l_nnlf, &
        powerlaw_loglikelihood, powerlaw_nnlf, triang_loglikelihood, triang_nnlf, &
        genpareto_loglikelihood, genpareto_nnlf, arcsine_loglikelihood, arcsine_nnlf, &
        halfnorm_loglikelihood, halfnorm_nnlf, halfcauchy_loglikelihood, &
        halfcauchy_nnlf, lomax_loglikelihood, lomax_nnlf, &
        chi_loglikelihood, chi_nnlf, maxwell_loglikelihood, maxwell_nnlf, &
        cosine_loglikelihood, cosine_nnlf, semicircular_loglikelihood, semicircular_nnlf, &
        anglit_loglikelihood, anglit_nnlf, moyal_loglikelihood, moyal_nnlf, &
        hypsecant_loglikelihood, hypsecant_nnlf, &
        halflogistic_loglikelihood, halflogistic_nnlf, &
        invgamma_loglikelihood, invgamma_nnlf, invgauss_loglikelihood, invgauss_nnlf, &
        levy_loglikelihood, levy_nnlf, loglaplace_loglikelihood, loglaplace_nnlf, &
        bradford_loglikelihood, bradford_nnlf, &
        truncexpon_loglikelihood, truncexpon_nnlf, &
        fisk_loglikelihood, fisk_nnlf, dweibull_loglikelihood, dweibull_nnlf, &
        alpha_loglikelihood, alpha_nnlf, fatiguelife_loglikelihood, fatiguelife_nnlf, &
        genlogistic_loglikelihood, genlogistic_nnlf, gennorm_loglikelihood, gennorm_nnlf, &
        nakagami_loglikelihood, nakagami_nnlf, powernorm_loglikelihood, powernorm_nnlf, &
        loggamma_loglikelihood, loggamma_nnlf, wald_loglikelihood, wald_nnlf, &
        gompertz_loglikelihood, gompertz_nnlf, invweibull_loglikelihood, &
        invweibull_nnlf, betaprime_loglikelihood, betaprime_nnlf, &
        burr12_loglikelihood, burr12_nnlf, genhalflogistic_loglikelihood, &
        genhalflogistic_nnlf, exponpow_loglikelihood, exponpow_nnlf, &
        exponweib_loglikelihood, exponweib_nnlf, powerlognorm_loglikelihood, powerlognorm_nnlf, &
        levy_l_loglikelihood, levy_l_nnlf, weibull_max_loglikelihood, weibull_max_nnlf, &
        rdist_loglikelihood, rdist_nnlf, skewcauchy_loglikelihood, skewcauchy_nnlf, &
        dgamma_loglikelihood, dgamma_nnlf, laplace_asymmetric_loglikelihood, &
        laplace_asymmetric_nnlf, truncnorm_loglikelihood, truncnorm_nnlf, &
        loguniform_loglikelihood, loguniform_nnlf
    use scifort_logistic, only : logistic_cdf, logistic_isf, logistic_logcdf, &
        logistic_logpdf, logistic_logsf, logistic_pdf, logistic_ppf, logistic_sf
    use scifort_lognormal, only : lognormal_cdf, lognormal_isf, lognormal_logcdf, &
        lognormal_logpdf, lognormal_logsf, lognormal_pdf, lognormal_ppf, lognormal_sf
    use scifort_negative_binomial, only : negative_binomial_cdf, negative_binomial_isf, &
        negative_binomial_logcdf, negative_binomial_logpmf, negative_binomial_logsf, &
        negative_binomial_pmf, negative_binomial_ppf, negative_binomial_sf
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logpdf, normal_logsf, normal_pdf, normal_ppf, normal_sf
    use scifort_pareto, only : pareto_cdf, pareto_isf, pareto_logcdf, &
        pareto_logpdf, pareto_logsf, pareto_pdf, pareto_ppf, pareto_sf
    use scifort_poisson, only : poisson_cdf, poisson_isf, poisson_logcdf, poisson_logpmf, &
        poisson_logsf, poisson_pmf, poisson_ppf, poisson_sf
    use scifort_resampling, only : bootstrap, bootstrap_percentile, bootstrap_result, &
        monte_carlo_test, monte_carlo_test_result, permutation_test, permutation_test_result, &
        power, power_result
    use scifort_random_variates, only : normal_rvs, normal_rvs_array, uniform_rvs, &
        uniform_rvs_array, exponential_rvs, exponential_rvs_array, laplace_rvs, &
        laplace_rvs_array, logistic_rvs, logistic_rvs_array, cauchy_rvs, cauchy_rvs_array, &
        rayleigh_rvs, rayleigh_rvs_array, gamma_rvs, gamma_rvs_array, chi2_rvs, chi2_rvs_array, &
        t_rvs, t_rvs_array, lognormal_rvs, lognormal_rvs_array, weibull_rvs, weibull_rvs_array, &
        pareto_rvs, pareto_rvs_array, beta_rvs, beta_rvs_array, f_rvs, f_rvs_array, &
        bernoulli_rvs, bernoulli_rvs_array, poisson_rvs, poisson_rvs_array, geometric_rvs, &
        geometric_rvs_array, binomial_rvs, binomial_rvs_array, negative_binomial_rvs, &
        negative_binomial_rvs_array, gumbel_r_rvs, gumbel_r_rvs_array, &
        gumbel_l_rvs, gumbel_l_rvs_array, powerlaw_rvs, powerlaw_rvs_array, &
        triang_rvs, triang_rvs_array, genpareto_rvs, genpareto_rvs_array, &
        arcsine_rvs, arcsine_rvs_array, halfnorm_rvs, halfnorm_rvs_array, &
        halfcauchy_rvs, halfcauchy_rvs_array, lomax_rvs, lomax_rvs_array, &
        chi_rvs, chi_rvs_array, maxwell_rvs, maxwell_rvs_array, &
        cosine_rvs, cosine_rvs_array, semicircular_rvs, semicircular_rvs_array, &
        anglit_rvs, anglit_rvs_array, moyal_rvs, moyal_rvs_array, &
        hypsecant_rvs, hypsecant_rvs_array, halflogistic_rvs, halflogistic_rvs_array, &
        invgamma_rvs, invgamma_rvs_array, invgauss_rvs, invgauss_rvs_array, &
        levy_rvs, levy_rvs_array, loglaplace_rvs, loglaplace_rvs_array, &
        bradford_rvs, bradford_rvs_array, truncexpon_rvs, truncexpon_rvs_array, &
        fisk_rvs, fisk_rvs_array, dweibull_rvs, dweibull_rvs_array, &
        alpha_rvs, alpha_rvs_array, fatiguelife_rvs, fatiguelife_rvs_array, &
        genlogistic_rvs, genlogistic_rvs_array, gennorm_rvs, gennorm_rvs_array, &
        nakagami_rvs, nakagami_rvs_array, powernorm_rvs, powernorm_rvs_array, &
        loggamma_rvs, loggamma_rvs_array, wald_rvs, wald_rvs_array, &
        gompertz_rvs, gompertz_rvs_array, invweibull_rvs, invweibull_rvs_array, &
        betaprime_rvs, betaprime_rvs_array, burr12_rvs, burr12_rvs_array, &
        genhalflogistic_rvs, genhalflogistic_rvs_array, exponpow_rvs, exponpow_rvs_array, &
        exponweib_rvs, exponweib_rvs_array, powerlognorm_rvs, powerlognorm_rvs_array, &
        levy_l_rvs, levy_l_rvs_array, weibull_max_rvs, weibull_max_rvs_array, &
        rdist_rvs, rdist_rvs_array, skewcauchy_rvs, skewcauchy_rvs_array, &
        dgamma_rvs, dgamma_rvs_array, laplace_asymmetric_rvs, &
        laplace_asymmetric_rvs_array, truncnorm_rvs, truncnorm_rvs_array, &
        loguniform_rvs, loguniform_rvs_array
    use scifort_score, only : normal_score, uniform_score, exponential_score, &
        laplace_score, logistic_score, cauchy_score, rayleigh_score, gamma_score, &
        chi2_score, t_score, lognormal_score, weibull_score, pareto_score, beta_score, &
        f_score, bernoulli_score, poisson_score, geometric_score, binomial_score, &
        negative_binomial_score, gumbel_r_score, gumbel_l_score, powerlaw_score, &
        triang_score, genpareto_score, arcsine_score, halfnorm_score, &
        halfcauchy_score, lomax_score, chi_score, maxwell_score, &
        cosine_score, semicircular_score, anglit_score, moyal_score, &
        hypsecant_score, halflogistic_score, invgamma_score, invgauss_score, &
        levy_score, loglaplace_score, bradford_score, truncexpon_score, &
        fisk_score, dweibull_score, alpha_score, fatiguelife_score, &
        genlogistic_score, gennorm_score, nakagami_score, powernorm_score, &
        loggamma_score, wald_score, gompertz_score, invweibull_score, &
        betaprime_score, burr12_score, genhalflogistic_score, exponpow_score, &
        exponweib_score, powerlognorm_score, levy_l_score, weibull_max_score, &
        rdist_score, skewcauchy_score, dgamma_score, laplace_asymmetric_score, &
        truncnorm_score, loguniform_score
    use scifort_rayleigh, only : rayleigh_cdf, rayleigh_isf, rayleigh_logcdf, &
        rayleigh_logpdf, rayleigh_logsf, rayleigh_pdf, rayleigh_ppf, rayleigh_sf
    use scifort_student_t, only : t_cdf, t_isf, t_logcdf, t_logpdf, &
        t_logsf, t_pdf, t_ppf, t_sf
    use scifort_uniform, only : uniform_cdf, uniform_isf, uniform_logcdf, &
        uniform_logpdf, uniform_logsf, uniform_pdf, uniform_ppf, uniform_sf
    use scifort_weibull, only : weibull_cdf, weibull_isf, weibull_logcdf, &
        weibull_logpdf, weibull_logsf, weibull_pdf, weibull_ppf, weibull_sf
    use scifort_foldnorm, only : foldnorm_cdf, foldnorm_isf, foldnorm_logcdf, &
        foldnorm_logpdf, foldnorm_logsf, foldnorm_pdf, foldnorm_ppf, foldnorm_sf
    use scifort_foldcauchy, only : foldcauchy_cdf, foldcauchy_isf, foldcauchy_logcdf, &
        foldcauchy_logpdf, foldcauchy_logsf, foldcauchy_pdf, foldcauchy_ppf, foldcauchy_sf
    use scifort_recipinvgauss, only : recipinvgauss_cdf, recipinvgauss_isf, &
        recipinvgauss_logcdf, recipinvgauss_logpdf, recipinvgauss_logsf, &
        recipinvgauss_pdf, recipinvgauss_ppf, recipinvgauss_sf
    use scifort_truncpareto, only : truncpareto_cdf, truncpareto_isf, truncpareto_logcdf, &
        truncpareto_logpdf, truncpareto_logsf, truncpareto_pdf, truncpareto_ppf, truncpareto_sf
    use scifort_likelihood, only : foldnorm_loglikelihood, foldnorm_nnlf, &
        foldcauchy_loglikelihood, foldcauchy_nnlf, recipinvgauss_loglikelihood, &
        recipinvgauss_nnlf, truncpareto_loglikelihood, truncpareto_nnlf
    use scifort_score, only : foldnorm_score, foldcauchy_score, recipinvgauss_score, truncpareto_score
    use scifort_fit, only : foldnorm_fit, foldcauchy_fit, recipinvgauss_fit, truncpareto_fit
    use scifort_random_variates, only : foldnorm_rvs, foldnorm_rvs_array, &
        foldcauchy_rvs, foldcauchy_rvs_array, recipinvgauss_rvs, recipinvgauss_rvs_array, &
        truncpareto_rvs, truncpareto_rvs_array
    use scifort_exponnorm, only : exponnorm_cdf, exponnorm_isf, exponnorm_logcdf, &
        exponnorm_logpdf, exponnorm_logsf, exponnorm_pdf, exponnorm_ppf, exponnorm_sf
    use scifort_johnsonsb, only : johnsonsb_cdf, johnsonsb_isf, johnsonsb_logcdf, &
        johnsonsb_logpdf, johnsonsb_logsf, johnsonsb_pdf, johnsonsb_ppf, johnsonsb_sf
    use scifort_johnsonsu, only : johnsonsu_cdf, johnsonsu_isf, johnsonsu_logcdf, &
        johnsonsu_logpdf, johnsonsu_logsf, johnsonsu_pdf, johnsonsu_ppf, johnsonsu_sf
    use scifort_trapezoid, only : trapezoid_cdf, trapezoid_isf, trapezoid_logcdf, &
        trapezoid_logpdf, trapezoid_logsf, trapezoid_pdf, trapezoid_ppf, trapezoid_sf
    use scifort_likelihood, only : exponnorm_loglikelihood, exponnorm_nnlf, &
        johnsonsb_loglikelihood, johnsonsb_nnlf, johnsonsu_loglikelihood, johnsonsu_nnlf, &
        trapezoid_loglikelihood, trapezoid_nnlf
    use scifort_score, only : exponnorm_score, johnsonsb_score, johnsonsu_score, trapezoid_score
    use scifort_fit, only : exponnorm_fit, johnsonsb_fit, johnsonsu_fit, trapezoid_fit
    use scifort_random_variates, only : exponnorm_rvs, exponnorm_rvs_array, &
        johnsonsb_rvs, johnsonsb_rvs_array, johnsonsu_rvs, johnsonsu_rvs_array, &
        trapezoid_rvs, trapezoid_rvs_array
    use scifort_burr, only : burr_cdf, burr_isf, burr_logcdf, burr_logpdf, &
        burr_logsf, burr_pdf, burr_ppf, burr_sf
    use scifort_mielke, only : mielke_cdf, mielke_isf, mielke_logcdf, mielke_logpdf, &
        mielke_logsf, mielke_pdf, mielke_ppf, mielke_sf
    use scifort_gibrat, only : gibrat_cdf, gibrat_isf, gibrat_logcdf, gibrat_logpdf, &
        gibrat_logsf, gibrat_pdf, gibrat_ppf, gibrat_sf
    use scifort_wrapcauchy, only : wrapcauchy_cdf, wrapcauchy_isf, wrapcauchy_logcdf, &
        wrapcauchy_logpdf, wrapcauchy_logsf, wrapcauchy_pdf, wrapcauchy_ppf, wrapcauchy_sf
    use scifort_likelihood, only : burr_loglikelihood, burr_nnlf, mielke_loglikelihood, &
        mielke_nnlf, gibrat_loglikelihood, gibrat_nnlf, wrapcauchy_loglikelihood, &
        wrapcauchy_nnlf
    use scifort_score, only : burr_score, mielke_score, gibrat_score, wrapcauchy_score
    use scifort_fit, only : burr_fit, mielke_fit, gibrat_fit, wrapcauchy_fit
    use scifort_random_variates, only : burr_rvs, burr_rvs_array, mielke_rvs, &
        mielke_rvs_array, gibrat_rvs, gibrat_rvs_array, wrapcauchy_rvs, &
        wrapcauchy_rvs_array
    use scifort_genextreme, only : genextreme_cdf, genextreme_isf, genextreme_logcdf, &
        genextreme_logpdf, genextreme_logsf, genextreme_pdf, genextreme_ppf, genextreme_sf
    use scifort_kappa3, only : kappa3_cdf, kappa3_isf, kappa3_logcdf, kappa3_logpdf, &
        kappa3_logsf, kappa3_pdf, kappa3_ppf, kappa3_sf
    use scifort_kappa4, only : kappa4_cdf, kappa4_isf, kappa4_logcdf, kappa4_logpdf, &
        kappa4_logsf, kappa4_pdf, kappa4_ppf, kappa4_sf
    use scifort_truncweibull_min, only : truncweibull_min_cdf, truncweibull_min_isf, &
        truncweibull_min_logcdf, truncweibull_min_logpdf, truncweibull_min_logsf, &
        truncweibull_min_pdf, truncweibull_min_ppf, truncweibull_min_sf
    use scifort_likelihood, only : genextreme_loglikelihood, genextreme_nnlf, &
        kappa3_loglikelihood, kappa3_nnlf, kappa4_loglikelihood, kappa4_nnlf, &
        truncweibull_min_loglikelihood, truncweibull_min_nnlf
    use scifort_score, only : genextreme_score, kappa3_score, kappa4_score, truncweibull_min_score
    use scifort_fit, only : genextreme_fit, kappa3_fit, kappa4_fit, truncweibull_min_fit
    use scifort_random_variates, only : genextreme_rvs, genextreme_rvs_array, &
        kappa3_rvs, kappa3_rvs_array, kappa4_rvs, kappa4_rvs_array, &
        truncweibull_min_rvs, truncweibull_min_rvs_array
    use scifort_gengamma, only : gengamma_cdf, gengamma_isf, gengamma_logcdf, &
        gengamma_logpdf, gengamma_logsf, gengamma_pdf, gengamma_ppf, gengamma_sf
    use scifort_halfgennorm, only : halfgennorm_cdf, halfgennorm_isf, halfgennorm_logcdf, &
        halfgennorm_logpdf, halfgennorm_logsf, halfgennorm_pdf, halfgennorm_ppf, halfgennorm_sf
    use scifort_argus, only : argus_cdf, argus_isf, argus_logcdf, argus_logpdf, &
        argus_logsf, argus_pdf, argus_ppf, argus_sf
    use scifort_erlang, only : erlang_cdf, erlang_isf, erlang_logcdf, erlang_logpdf, &
        erlang_logsf, erlang_pdf, erlang_ppf, erlang_sf
    use scifort_likelihood, only : gengamma_loglikelihood, gengamma_nnlf, &
        halfgennorm_loglikelihood, halfgennorm_nnlf, argus_loglikelihood, argus_nnlf, &
        erlang_loglikelihood, erlang_nnlf
    use scifort_score, only : gengamma_score, halfgennorm_score, argus_score, erlang_score
    use scifort_fit, only : gengamma_fit, halfgennorm_fit, argus_fit, erlang_fit
    use scifort_random_variates, only : gengamma_rvs, gengamma_rvs_array, &
        halfgennorm_rvs, halfgennorm_rvs_array, argus_rvs, argus_rvs_array, &
        erlang_rvs, erlang_rvs_array
    use scifort_crystalball, only : crystalball_cdf, crystalball_isf, &
        crystalball_logcdf, crystalball_logpdf, crystalball_logsf, crystalball_pdf, &
        crystalball_ppf, crystalball_sf
    use scifort_jf_skew_t, only : jf_skew_t_cdf, jf_skew_t_isf, jf_skew_t_logcdf, &
        jf_skew_t_logpdf, jf_skew_t_logsf, jf_skew_t_pdf, jf_skew_t_ppf, jf_skew_t_sf
    use scifort_pearson3, only : pearson3_cdf, pearson3_isf, pearson3_logcdf, &
        pearson3_logpdf, pearson3_logsf, pearson3_pdf, pearson3_ppf, pearson3_sf
    use scifort_rel_breitwigner, only : rel_breitwigner_cdf, rel_breitwigner_isf, &
        rel_breitwigner_logcdf, rel_breitwigner_logpdf, rel_breitwigner_logsf, &
        rel_breitwigner_pdf, rel_breitwigner_ppf, rel_breitwigner_sf
    use scifort_likelihood, only : crystalball_loglikelihood, crystalball_nnlf, &
        jf_skew_t_loglikelihood, jf_skew_t_nnlf, pearson3_loglikelihood, pearson3_nnlf, &
        rel_breitwigner_loglikelihood, rel_breitwigner_nnlf
    use scifort_score, only : crystalball_score, jf_skew_t_score, pearson3_score, &
        rel_breitwigner_score
    use scifort_fit, only : crystalball_fit, jf_skew_t_fit, pearson3_fit, &
        rel_breitwigner_fit
    use scifort_random_variates, only : crystalball_rvs, crystalball_rvs_array, &
        jf_skew_t_rvs, jf_skew_t_rvs_array, pearson3_rvs, pearson3_rvs_array, &
        rel_breitwigner_rvs, rel_breitwigner_rvs_array
    use scifort_genexpon, only : genexpon_cdf, genexpon_isf, genexpon_logcdf, &
        genexpon_logpdf, genexpon_logsf, genexpon_pdf, genexpon_ppf, genexpon_sf
    use scifort_skewnorm, only : skewnorm_cdf, skewnorm_isf, skewnorm_logcdf, &
        skewnorm_logpdf, skewnorm_logsf, skewnorm_pdf, skewnorm_ppf, skewnorm_sf
    use scifort_tukeylambda, only : tukeylambda_cdf, tukeylambda_isf, &
        tukeylambda_logcdf, tukeylambda_logpdf, tukeylambda_logsf, tukeylambda_pdf, &
        tukeylambda_ppf, tukeylambda_sf
    use scifort_rice, only : rice_cdf, rice_isf, rice_logcdf, rice_logpdf, &
        rice_logsf, rice_pdf, rice_ppf, rice_sf
    use scifort_likelihood, only : genexpon_loglikelihood, genexpon_nnlf, &
        skewnorm_loglikelihood, skewnorm_nnlf, tukeylambda_loglikelihood, &
        tukeylambda_nnlf, rice_loglikelihood, rice_nnlf
    use scifort_score, only : genexpon_score, skewnorm_score, tukeylambda_score, rice_score
    use scifort_fit, only : genexpon_fit, skewnorm_fit, tukeylambda_fit, rice_fit
    use scifort_random_variates, only : genexpon_rvs, genexpon_rvs_array, &
        skewnorm_rvs, skewnorm_rvs_array, tukeylambda_rvs, tukeylambda_rvs_array, &
        rice_rvs, rice_rvs_array
    use scifort_dpareto_lognorm, only : dpareto_lognorm_cdf, dpareto_lognorm_isf, &
        dpareto_lognorm_logcdf, dpareto_lognorm_logpdf, dpareto_lognorm_logsf, &
        dpareto_lognorm_pdf, dpareto_lognorm_ppf, dpareto_lognorm_sf
    use scifort_vonmises, only : vonmises_cdf, vonmises_isf, vonmises_logcdf, &
        vonmises_logpdf, vonmises_logsf, vonmises_pdf, vonmises_ppf, vonmises_sf
    use scifort_vonmises_line, only : vonmises_line_cdf, vonmises_line_isf, &
        vonmises_line_logcdf, vonmises_line_logpdf, vonmises_line_logsf, &
        vonmises_line_pdf, vonmises_line_ppf, vonmises_line_sf
    use scifort_kstwobign, only : kstwobign_cdf, kstwobign_isf, kstwobign_logcdf, &
        kstwobign_logpdf, kstwobign_logsf, kstwobign_pdf, kstwobign_ppf, kstwobign_sf
    use scifort_likelihood, only : dpareto_lognorm_loglikelihood, dpareto_lognorm_nnlf, &
        vonmises_loglikelihood, vonmises_nnlf, vonmises_line_loglikelihood, &
        vonmises_line_nnlf, kstwobign_loglikelihood, kstwobign_nnlf
    use scifort_score, only : dpareto_lognorm_score, vonmises_score, &
        vonmises_line_score, kstwobign_score
    use scifort_fit, only : dpareto_lognorm_fit, vonmises_fit, vonmises_line_fit, &
        kstwobign_fit
    use scifort_random_variates, only : dpareto_lognorm_rvs, dpareto_lognorm_rvs_array, &
        vonmises_rvs, vonmises_rvs_array, vonmises_line_rvs, vonmises_line_rvs_array, &
        kstwobign_rvs, kstwobign_rvs_array
    use scifort_irwinhall, only : irwinhall_cdf, irwinhall_isf, irwinhall_logcdf, &
        irwinhall_logpdf, irwinhall_logsf, irwinhall_pdf, irwinhall_ppf, irwinhall_sf
    use scifort_ksone, only : ksone_cdf, ksone_isf, ksone_logcdf, ksone_logpdf, &
        ksone_logsf, ksone_pdf, ksone_ppf, ksone_sf
    use scifort_kstwo, only : kstwo_cdf, kstwo_isf, kstwo_logcdf, kstwo_logpdf, &
        kstwo_logsf, kstwo_pdf, kstwo_ppf, kstwo_sf
    use scifort_ncx2, only : ncx2_cdf, ncx2_isf, ncx2_logcdf, ncx2_logpdf, ncx2_logsf, &
        ncx2_pdf, ncx2_ppf, ncx2_sf
    use scifort_ncf, only : ncf_cdf, ncf_isf, ncf_logcdf, ncf_logpdf, ncf_logsf, ncf_pdf, &
        ncf_ppf, ncf_sf
    use scifort_likelihood, only : irwinhall_loglikelihood, irwinhall_nnlf, &
        ksone_loglikelihood, ksone_nnlf, kstwo_loglikelihood, kstwo_nnlf, &
        ncx2_loglikelihood, ncx2_nnlf, ncf_loglikelihood, ncf_nnlf
    use scifort_score, only : irwinhall_score, ksone_score, kstwo_score, ncx2_score, ncf_score
    use scifort_fit, only : irwinhall_fit, ksone_fit, kstwo_fit, ncx2_fit, ncf_fit
    use scifort_random_variates, only : irwinhall_rvs, irwinhall_rvs_array, &
        ksone_rvs, ksone_rvs_array, kstwo_rvs, kstwo_rvs_array, &
        ncx2_rvs, ncx2_rvs_array, ncf_rvs, ncf_rvs_array
    use scifort_randint, only : randint_cdf, randint_isf, randint_logcdf, randint_logpmf, &
        randint_logsf, randint_pmf, randint_ppf, randint_sf
    use scifort_planck, only : planck_cdf, planck_isf, planck_logcdf, planck_logpmf, &
        planck_logsf, planck_pmf, planck_ppf, planck_sf
    use scifort_dlaplace, only : dlaplace_cdf, dlaplace_isf, dlaplace_logcdf, dlaplace_logpmf, &
        dlaplace_logsf, dlaplace_pmf, dlaplace_ppf, dlaplace_sf
    use scifort_logser, only : logser_cdf, logser_isf, logser_logcdf, logser_logpmf, &
        logser_logsf, logser_pmf, logser_ppf, logser_sf
    use scifort_likelihood, only : randint_loglikelihood, randint_nnlf, planck_loglikelihood, &
        planck_nnlf, dlaplace_loglikelihood, dlaplace_nnlf, logser_loglikelihood, logser_nnlf
    use scifort_score, only : randint_score, planck_score, dlaplace_score, logser_score
    use scifort_fit, only : randint_fit, planck_fit, dlaplace_fit, logser_fit
    use scifort_random_variates, only : randint_rvs, randint_rvs_array, planck_rvs, &
        planck_rvs_array, dlaplace_rvs, dlaplace_rvs_array, logser_rvs, logser_rvs_array
    use scifort_betabinom, only : betabinom_cdf, betabinom_isf, betabinom_logcdf, &
        betabinom_logpmf, betabinom_logsf, betabinom_pmf, betabinom_ppf, betabinom_sf
    use scifort_hypergeom, only : hypergeom_cdf, hypergeom_isf, hypergeom_logcdf, &
        hypergeom_logpmf, hypergeom_logsf, hypergeom_pmf, hypergeom_ppf, hypergeom_sf
    use scifort_nhypergeom, only : nhypergeom_cdf, nhypergeom_isf, nhypergeom_logcdf, &
        nhypergeom_logpmf, nhypergeom_logsf, nhypergeom_pmf, nhypergeom_ppf, nhypergeom_sf
    use scifort_boltzmann, only : boltzmann_cdf, boltzmann_isf, boltzmann_logcdf, &
        boltzmann_logpmf, boltzmann_logsf, boltzmann_pmf, boltzmann_ppf, boltzmann_sf
    use scifort_betanbinom, only : betanbinom_cdf, betanbinom_isf, betanbinom_logcdf, &
        betanbinom_logpmf, betanbinom_logsf, betanbinom_pmf, betanbinom_ppf, betanbinom_sf
    use scifort_yulesimon, only : yulesimon_cdf, yulesimon_isf, yulesimon_logcdf, &
        yulesimon_logpmf, yulesimon_logsf, yulesimon_pmf, yulesimon_ppf, yulesimon_sf
    use scifort_zipf, only : zipf_cdf, zipf_isf, zipf_logcdf, zipf_logpmf, &
        zipf_logsf, zipf_pmf, zipf_ppf, zipf_sf
    use scifort_zipfian, only : zipfian_cdf, zipfian_isf, zipfian_logcdf, &
        zipfian_logpmf, zipfian_logsf, zipfian_pmf, zipfian_ppf, zipfian_sf
    use scifort_likelihood, only : betabinom_loglikelihood, betabinom_nnlf, &
        hypergeom_loglikelihood, hypergeom_nnlf, nhypergeom_loglikelihood, nhypergeom_nnlf, &
        boltzmann_loglikelihood, boltzmann_nnlf, betanbinom_loglikelihood, betanbinom_nnlf, &
        yulesimon_loglikelihood, yulesimon_nnlf, zipf_loglikelihood, zipf_nnlf, &
        zipfian_loglikelihood, zipfian_nnlf
    use scifort_score, only : betabinom_score, hypergeom_score, nhypergeom_score, boltzmann_score, &
        betanbinom_score, yulesimon_score, zipf_score, zipfian_score
    use scifort_fit, only : betabinom_fit, hypergeom_fit, nhypergeom_fit, boltzmann_fit, &
        betanbinom_fit, yulesimon_fit, zipf_fit, zipfian_fit
    use scifort_random_variates, only : betabinom_rvs, betabinom_rvs_array, &
        hypergeom_rvs, hypergeom_rvs_array, nhypergeom_rvs, nhypergeom_rvs_array, &
        boltzmann_rvs, boltzmann_rvs_array, betanbinom_rvs, betanbinom_rvs_array, &
        yulesimon_rvs, yulesimon_rvs_array, zipf_rvs, zipf_rvs_array, &
        zipfian_rvs, zipfian_rvs_array
    use scifort_geninvgauss, only : geninvgauss_cdf, geninvgauss_isf, geninvgauss_logcdf, &
        geninvgauss_logpdf, geninvgauss_logsf, geninvgauss_pdf, geninvgauss_ppf, geninvgauss_sf
    use scifort_norminvgauss, only : norminvgauss_cdf, norminvgauss_isf, norminvgauss_logcdf, &
        norminvgauss_logpdf, norminvgauss_logsf, norminvgauss_pdf, norminvgauss_ppf, norminvgauss_sf
    use scifort_skellam, only : skellam_cdf, skellam_isf, skellam_logcdf, skellam_logpmf, &
        skellam_logsf, skellam_pmf, skellam_ppf, skellam_sf
    use scifort_likelihood, only : geninvgauss_loglikelihood, geninvgauss_nnlf, &
        norminvgauss_loglikelihood, norminvgauss_nnlf, skellam_loglikelihood, skellam_nnlf
    use scifort_score, only : geninvgauss_score, norminvgauss_score, skellam_score
    use scifort_fit, only : geninvgauss_fit, norminvgauss_fit, skellam_fit
    use scifort_random_variates, only : geninvgauss_rvs, geninvgauss_rvs_array, &
        norminvgauss_rvs, norminvgauss_rvs_array, skellam_rvs, skellam_rvs_array
    use scifort_genhyperbolic, only : genhyperbolic_cdf, genhyperbolic_isf, genhyperbolic_logcdf, &
        genhyperbolic_logpdf, genhyperbolic_logsf, genhyperbolic_pdf, genhyperbolic_ppf, genhyperbolic_sf
    use scifort_nchypergeom_fisher, only : nchypergeom_fisher_cdf, nchypergeom_fisher_isf, &
        nchypergeom_fisher_logcdf, nchypergeom_fisher_logpmf, nchypergeom_fisher_logsf, &
        nchypergeom_fisher_pmf, nchypergeom_fisher_ppf, nchypergeom_fisher_sf
    use scifort_likelihood, only : genhyperbolic_loglikelihood, genhyperbolic_nnlf, &
        nchypergeom_fisher_loglikelihood, nchypergeom_fisher_nnlf
    use scifort_score, only : genhyperbolic_score, nchypergeom_fisher_score
    use scifort_fit, only : genhyperbolic_fit, nchypergeom_fisher_fit
    use scifort_random_variates, only : genhyperbolic_rvs, genhyperbolic_rvs_array, &
        nchypergeom_fisher_rvs, nchypergeom_fisher_rvs_array
    use scifort_nct, only : nct_cdf, nct_isf, nct_logcdf, nct_logpdf, &
        nct_logsf, nct_pdf, nct_ppf, nct_sf
    use scifort_gausshyper, only : gausshyper_cdf, gausshyper_isf, gausshyper_logcdf, &
        gausshyper_logpdf, gausshyper_logsf, gausshyper_pdf, gausshyper_ppf, gausshyper_sf
    use scifort_likelihood, only : nct_loglikelihood, nct_nnlf, gausshyper_loglikelihood, gausshyper_nnlf
    use scifort_score, only : nct_score, gausshyper_score
    use scifort_fit, only : nct_fit, gausshyper_fit
    use scifort_random_variates, only : nct_rvs, nct_rvs_array, gausshyper_rvs, gausshyper_rvs_array
    use scifort_landau, only : landau_cdf, landau_isf, landau_logcdf, landau_logpdf, &
        landau_logsf, landau_pdf, landau_ppf, landau_sf
    use scifort_likelihood, only : landau_loglikelihood, landau_nnlf
    use scifort_score, only : landau_score
    use scifort_fit, only : landau_fit
    use scifort_random_variates, only : landau_rvs, landau_rvs_array
    use scifort_nchypergeom_wallenius, only : nchypergeom_wallenius_cdf, nchypergeom_wallenius_isf, &
        nchypergeom_wallenius_logcdf, nchypergeom_wallenius_logpmf, nchypergeom_wallenius_logsf, &
        nchypergeom_wallenius_pmf, nchypergeom_wallenius_ppf, nchypergeom_wallenius_sf
    use scifort_poisson_binom, only : poisson_binom_cdf, poisson_binom_isf, poisson_binom_logcdf, &
        poisson_binom_logpmf, poisson_binom_logsf, poisson_binom_pmf, poisson_binom_ppf, poisson_binom_sf
    use scifort_likelihood, only : nchypergeom_wallenius_loglikelihood, nchypergeom_wallenius_nnlf, &
        poisson_binom_loglikelihood, poisson_binom_nnlf
    use scifort_score, only : nchypergeom_wallenius_score, poisson_binom_score
    use scifort_fit, only : nchypergeom_wallenius_fit, poisson_binom_fit
    use scifort_random_variates, only : nchypergeom_wallenius_rvs, nchypergeom_wallenius_rvs_array, &
        poisson_binom_rvs, poisson_binom_rvs_array
    use scifort_levy_stable, only : levy_stable_cdf, levy_stable_isf, &
        levy_stable_logcdf, levy_stable_logpdf, levy_stable_logsf, levy_stable_pdf, &
        levy_stable_ppf, levy_stable_sf
    use scifort_studentized_range, only : studentized_range_cdf, studentized_range_isf, &
        studentized_range_logcdf, studentized_range_logpdf, studentized_range_logsf, &
        studentized_range_pdf, studentized_range_ppf, studentized_range_sf
    use scifort_likelihood, only : levy_stable_loglikelihood, levy_stable_nnlf, &
        studentized_range_loglikelihood, studentized_range_nnlf
    use scifort_score, only : levy_stable_score, studentized_range_score
    use scifort_fit, only : levy_stable_fit, studentized_range_fit
    use scifort_random_variates, only : levy_stable_rvs, levy_stable_rvs_array, &
        studentized_range_rvs, studentized_range_rvs_array
    use scifort_multivariate_normal, only : multivariate_normal_cdf, multivariate_normal_entropy, &
        multivariate_normal_fit, multivariate_normal_logcdf, multivariate_normal_logpdf, &
        multivariate_normal_marginal, multivariate_normal_pdf, multivariate_normal_rvs, &
        multivariate_normal_rvs_array, multivariate_status_success, multivariate_status_invalid_shape, &
        multivariate_status_invalid_parameter, multivariate_status_linalg_failure
    use scifort_gaussian_kde, only : gaussian_kde, gaussian_kde_init, &
        gaussian_kde_status_success, gaussian_kde_status_invalid_shape, &
        gaussian_kde_status_invalid_parameter, gaussian_kde_status_singular, &
        gaussian_kde_status_linalg_failure
    use scifort_qmc, only : qmc_discrepancy, qmc_geometric_discrepancy, qmc_halton, &
        qmc_halton_init, qmc_latin_hypercube, qmc_latin_hypercube_init, qmc_scale, &
        qmc_sobol, qmc_sobol_init, qmc_multinomial, qmc_multinomial_init, &
        qmc_multivariate_normal, qmc_multivariate_normal_init, qmc_poisson_disk, &
        qmc_poisson_disk_init, qmc_to_integers, qmc_update_discrepancy, qmc_van_der_corput, &
        qmc_status_success, qmc_status_invalid_shape, qmc_status_invalid_parameter, &
        qmc_status_out_of_bounds
    use scifort_multivariate_t, only : multivariate_t_cdf, multivariate_t_entropy, &
        multivariate_t_logpdf, multivariate_t_marginal, multivariate_t_pdf, multivariate_t_rvs, &
        multivariate_t_rvs_array, multivariate_t_status_success, multivariate_t_status_invalid_shape, &
        multivariate_t_status_invalid_parameter, multivariate_t_status_linalg_failure
    use scifort_dirichlet, only : dirichlet_cov, dirichlet_entropy, dirichlet_logpdf, &
        dirichlet_mean, dirichlet_pdf, dirichlet_rvs, dirichlet_rvs_array, dirichlet_var, &
        dirichlet_status_success, dirichlet_status_invalid_parameter, dirichlet_status_invalid_shape
    use scifort_multinomial, only : multinomial_cov, multinomial_entropy, multinomial_logpmf, &
        multinomial_mean, multinomial_pmf, multinomial_rvs, multinomial_rvs_array, &
        multinomial_status_success, multinomial_status_invalid_parameter, multinomial_status_invalid_shape
    use scifort_dirichlet_multinomial, only : dirichlet_multinomial_cov, &
        dirichlet_multinomial_logpmf, dirichlet_multinomial_mean, dirichlet_multinomial_pmf, &
        dirichlet_multinomial_var
    use scifort_multivariate_hypergeom, only : multivariate_hypergeom_cov, &
        multivariate_hypergeom_logpmf, multivariate_hypergeom_mean, multivariate_hypergeom_pmf, &
        multivariate_hypergeom_rvs, multivariate_hypergeom_rvs_array, multivariate_hypergeom_var, &
        multivariate_hypergeom_status_success, multivariate_hypergeom_status_invalid_parameter, &
        multivariate_hypergeom_status_invalid_shape
    use scifort_normal_inverse_gamma, only : normal_inverse_gamma_logpdf, &
        normal_inverse_gamma_mean, normal_inverse_gamma_pdf, normal_inverse_gamma_rvs, &
        normal_inverse_gamma_rvs_array, normal_inverse_gamma_var
    use scifort_matrix_normal, only : matrix_normal_entropy, matrix_normal_logpdf, matrix_normal_pdf, &
        matrix_normal_rvs, matrix_normal_rvs_array, matrix_normal_status_success, &
        matrix_normal_status_invalid_shape, matrix_normal_status_invalid_parameter, &
        matrix_normal_status_linalg_failure
    use scifort_wishart, only : wishart_entropy, wishart_logpdf, wishart_mean, wishart_mode, &
        wishart_pdf, wishart_rvs, wishart_rvs_array, wishart_var, wishart_status_success, &
        wishart_status_invalid_shape, wishart_status_invalid_parameter, wishart_status_linalg_failure
    use scifort_invwishart, only : invwishart_entropy, invwishart_logpdf, invwishart_mean, &
        invwishart_mode, invwishart_pdf, invwishart_rvs, invwishart_rvs_array, invwishart_var, &
        invwishart_status_success, invwishart_status_invalid_shape, invwishart_status_invalid_parameter, &
        invwishart_status_linalg_failure
    use scifort_matrix_t, only : matrix_t_logpdf, matrix_t_pdf, matrix_t_rvs, matrix_t_rvs_array, &
        matrix_t_status_success, matrix_t_status_invalid_shape, matrix_t_status_invalid_parameter, &
        matrix_t_status_linalg_failure
    use scifort_uniform_direction, only : uniform_direction_rvs, uniform_direction_rvs_array, &
        uniform_direction_status_success, uniform_direction_status_invalid_shape, &
        uniform_direction_status_numerical_failure
    use scifort_vonmises_fisher, only : vonmises_fisher_entropy, vonmises_fisher_fit, &
        vonmises_fisher_logpdf, vonmises_fisher_pdf, vonmises_fisher_rvs, &
        vonmises_fisher_rvs_array, vonmises_fisher_status_success, &
        vonmises_fisher_status_invalid_shape, vonmises_fisher_status_invalid_parameter, &
        vonmises_fisher_status_numerical_failure
    use scifort_ortho_group, only : ortho_group_rvs, ortho_group_rvs_array, &
        ortho_group_status_success, ortho_group_status_invalid_shape, &
        ortho_group_status_numerical_failure
    use scifort_special_ortho_group, only : special_ortho_group_rvs, &
        special_ortho_group_rvs_array, special_ortho_group_status_success, &
        special_ortho_group_status_invalid_shape, special_ortho_group_status_numerical_failure
    use scifort_unitary_group, only : unitary_group_rvs, unitary_group_rvs_array, &
        unitary_group_status_success, unitary_group_status_invalid_shape, &
        unitary_group_status_numerical_failure
    use scifort_random_correlation, only : random_correlation_rvs, &
        random_correlation_status_success, random_correlation_status_invalid_shape, &
        random_correlation_status_invalid_parameter, random_correlation_status_numerical_failure
    use scifort_random_table, only : random_table_logpmf, random_table_mean, random_table_pmf, &
        random_table_rvs, random_table_rvs_array, random_table_status_success, &
        random_table_status_invalid_shape, random_table_status_invalid_parameter, &
        random_table_status_numerical_failure
    implicit none
    private

    public :: exponnorm_cdf, exponnorm_isf, exponnorm_logcdf, exponnorm_logpdf
    public :: exponnorm_logsf, exponnorm_pdf, exponnorm_ppf, exponnorm_sf
    public :: johnsonsb_cdf, johnsonsb_isf, johnsonsb_logcdf, johnsonsb_logpdf
    public :: johnsonsb_logsf, johnsonsb_pdf, johnsonsb_ppf, johnsonsb_sf
    public :: johnsonsu_cdf, johnsonsu_isf, johnsonsu_logcdf, johnsonsu_logpdf
    public :: johnsonsu_logsf, johnsonsu_pdf, johnsonsu_ppf, johnsonsu_sf
    public :: trapezoid_cdf, trapezoid_isf, trapezoid_logcdf, trapezoid_logpdf
    public :: trapezoid_logsf, trapezoid_pdf, trapezoid_ppf, trapezoid_sf
    public :: exponnorm_loglikelihood, exponnorm_nnlf, exponnorm_score, exponnorm_fit
    public :: johnsonsb_loglikelihood, johnsonsb_nnlf, johnsonsb_score, johnsonsb_fit
    public :: johnsonsu_loglikelihood, johnsonsu_nnlf, johnsonsu_score, johnsonsu_fit
    public :: trapezoid_loglikelihood, trapezoid_nnlf, trapezoid_score, trapezoid_fit
    public :: exponnorm_rvs, exponnorm_rvs_array
    public :: johnsonsb_rvs, johnsonsb_rvs_array
    public :: johnsonsu_rvs, johnsonsu_rvs_array
    public :: trapezoid_rvs, trapezoid_rvs_array
    public :: burr_cdf, burr_isf, burr_logcdf, burr_logpdf, burr_logsf, burr_pdf
    public :: burr_ppf, burr_sf, burr_loglikelihood, burr_nnlf, burr_score, burr_fit
    public :: burr_rvs, burr_rvs_array
    public :: mielke_cdf, mielke_isf, mielke_logcdf, mielke_logpdf, mielke_logsf
    public :: mielke_pdf, mielke_ppf, mielke_sf, mielke_loglikelihood, mielke_nnlf
    public :: mielke_score, mielke_fit, mielke_rvs, mielke_rvs_array
    public :: gibrat_cdf, gibrat_isf, gibrat_logcdf, gibrat_logpdf, gibrat_logsf
    public :: gibrat_pdf, gibrat_ppf, gibrat_sf, gibrat_loglikelihood, gibrat_nnlf
    public :: gibrat_score, gibrat_fit, gibrat_rvs, gibrat_rvs_array
    public :: wrapcauchy_cdf, wrapcauchy_isf, wrapcauchy_logcdf, wrapcauchy_logpdf
    public :: wrapcauchy_logsf, wrapcauchy_pdf, wrapcauchy_ppf, wrapcauchy_sf
    public :: wrapcauchy_loglikelihood, wrapcauchy_nnlf, wrapcauchy_score, wrapcauchy_fit
    public :: wrapcauchy_rvs, wrapcauchy_rvs_array
    public :: genextreme_cdf, genextreme_isf, genextreme_logcdf, genextreme_logpdf
    public :: genextreme_logsf, genextreme_pdf, genextreme_ppf, genextreme_sf
    public :: genextreme_loglikelihood, genextreme_nnlf, genextreme_score, genextreme_fit
    public :: genextreme_rvs, genextreme_rvs_array
    public :: kappa3_cdf, kappa3_isf, kappa3_logcdf, kappa3_logpdf, kappa3_logsf
    public :: kappa3_pdf, kappa3_ppf, kappa3_sf, kappa3_loglikelihood, kappa3_nnlf
    public :: kappa3_score, kappa3_fit, kappa3_rvs, kappa3_rvs_array
    public :: kappa4_cdf, kappa4_isf, kappa4_logcdf, kappa4_logpdf, kappa4_logsf
    public :: kappa4_pdf, kappa4_ppf, kappa4_sf, kappa4_loglikelihood, kappa4_nnlf
    public :: kappa4_score, kappa4_fit, kappa4_rvs, kappa4_rvs_array
    public :: truncweibull_min_cdf, truncweibull_min_isf, truncweibull_min_logcdf
    public :: truncweibull_min_logpdf, truncweibull_min_logsf, truncweibull_min_pdf
    public :: truncweibull_min_ppf, truncweibull_min_sf, truncweibull_min_loglikelihood
    public :: truncweibull_min_nnlf, truncweibull_min_score, truncweibull_min_fit
    public :: truncweibull_min_rvs, truncweibull_min_rvs_array

    public :: gengamma_cdf, gengamma_isf, gengamma_logcdf, gengamma_logpdf
    public :: gengamma_logsf, gengamma_pdf, gengamma_ppf, gengamma_sf
    public :: gengamma_loglikelihood, gengamma_nnlf, gengamma_score, gengamma_fit
    public :: gengamma_rvs, gengamma_rvs_array
    public :: halfgennorm_cdf, halfgennorm_isf, halfgennorm_logcdf, halfgennorm_logpdf
    public :: halfgennorm_logsf, halfgennorm_pdf, halfgennorm_ppf, halfgennorm_sf
    public :: halfgennorm_loglikelihood, halfgennorm_nnlf, halfgennorm_score, halfgennorm_fit
    public :: halfgennorm_rvs, halfgennorm_rvs_array
    public :: argus_cdf, argus_isf, argus_logcdf, argus_logpdf, argus_logsf
    public :: argus_pdf, argus_ppf, argus_sf, argus_loglikelihood, argus_nnlf
    public :: argus_score, argus_fit, argus_rvs, argus_rvs_array
    public :: erlang_cdf, erlang_isf, erlang_logcdf, erlang_logpdf, erlang_logsf
    public :: erlang_pdf, erlang_ppf, erlang_sf, erlang_loglikelihood, erlang_nnlf
    public :: erlang_score, erlang_fit, erlang_rvs, erlang_rvs_array

    public :: crystalball_cdf, crystalball_isf, crystalball_logcdf, crystalball_logpdf
    public :: crystalball_logsf, crystalball_pdf, crystalball_ppf, crystalball_sf
    public :: crystalball_loglikelihood, crystalball_nnlf, crystalball_score, crystalball_fit
    public :: crystalball_rvs, crystalball_rvs_array
    public :: jf_skew_t_cdf, jf_skew_t_isf, jf_skew_t_logcdf, jf_skew_t_logpdf
    public :: jf_skew_t_logsf, jf_skew_t_pdf, jf_skew_t_ppf, jf_skew_t_sf
    public :: jf_skew_t_loglikelihood, jf_skew_t_nnlf, jf_skew_t_score, jf_skew_t_fit
    public :: jf_skew_t_rvs, jf_skew_t_rvs_array
    public :: pearson3_cdf, pearson3_isf, pearson3_logcdf, pearson3_logpdf
    public :: pearson3_logsf, pearson3_pdf, pearson3_ppf, pearson3_sf
    public :: pearson3_loglikelihood, pearson3_nnlf, pearson3_score, pearson3_fit
    public :: pearson3_rvs, pearson3_rvs_array
    public :: rel_breitwigner_cdf, rel_breitwigner_isf, rel_breitwigner_logcdf
    public :: rel_breitwigner_logpdf, rel_breitwigner_logsf, rel_breitwigner_pdf
    public :: rel_breitwigner_ppf, rel_breitwigner_sf, rel_breitwigner_loglikelihood
    public :: rel_breitwigner_nnlf, rel_breitwigner_score, rel_breitwigner_fit
    public :: rel_breitwigner_rvs, rel_breitwigner_rvs_array
    public :: genexpon_cdf, genexpon_isf, genexpon_logcdf, genexpon_logpdf
    public :: genexpon_logsf, genexpon_pdf, genexpon_ppf, genexpon_sf
    public :: genexpon_loglikelihood, genexpon_nnlf, genexpon_score, genexpon_fit
    public :: genexpon_rvs, genexpon_rvs_array
    public :: skewnorm_cdf, skewnorm_isf, skewnorm_logcdf, skewnorm_logpdf
    public :: skewnorm_logsf, skewnorm_pdf, skewnorm_ppf, skewnorm_sf
    public :: skewnorm_loglikelihood, skewnorm_nnlf, skewnorm_score, skewnorm_fit
    public :: skewnorm_rvs, skewnorm_rvs_array
    public :: tukeylambda_cdf, tukeylambda_isf, tukeylambda_logcdf, tukeylambda_logpdf
    public :: tukeylambda_logsf, tukeylambda_pdf, tukeylambda_ppf, tukeylambda_sf
    public :: tukeylambda_loglikelihood, tukeylambda_nnlf, tukeylambda_score, tukeylambda_fit
    public :: tukeylambda_rvs, tukeylambda_rvs_array
    public :: rice_cdf, rice_isf, rice_logcdf, rice_logpdf, rice_logsf, rice_pdf
    public :: rice_ppf, rice_sf, rice_loglikelihood, rice_nnlf, rice_score, rice_fit
    public :: rice_rvs, rice_rvs_array
    public :: dpareto_lognorm_cdf, dpareto_lognorm_isf, dpareto_lognorm_logcdf
    public :: dpareto_lognorm_logpdf, dpareto_lognorm_logsf, dpareto_lognorm_pdf
    public :: dpareto_lognorm_ppf, dpareto_lognorm_sf, dpareto_lognorm_loglikelihood
    public :: dpareto_lognorm_nnlf, dpareto_lognorm_score, dpareto_lognorm_fit
    public :: dpareto_lognorm_rvs, dpareto_lognorm_rvs_array
    public :: vonmises_cdf, vonmises_isf, vonmises_logcdf, vonmises_logpdf
    public :: vonmises_logsf, vonmises_pdf, vonmises_ppf, vonmises_sf
    public :: vonmises_loglikelihood, vonmises_nnlf, vonmises_score, vonmises_fit
    public :: vonmises_rvs, vonmises_rvs_array
    public :: vonmises_line_cdf, vonmises_line_isf, vonmises_line_logcdf
    public :: vonmises_line_logpdf, vonmises_line_logsf, vonmises_line_pdf
    public :: vonmises_line_ppf, vonmises_line_sf, vonmises_line_loglikelihood
    public :: vonmises_line_nnlf, vonmises_line_score, vonmises_line_fit
    public :: vonmises_line_rvs, vonmises_line_rvs_array
    public :: kstwobign_cdf, kstwobign_isf, kstwobign_logcdf, kstwobign_logpdf
    public :: kstwobign_logsf, kstwobign_pdf, kstwobign_ppf, kstwobign_sf
    public :: kstwobign_loglikelihood, kstwobign_nnlf, kstwobign_score, kstwobign_fit
    public :: kstwobign_rvs, kstwobign_rvs_array
    public :: irwinhall_cdf, irwinhall_isf, irwinhall_logcdf, irwinhall_logpdf
    public :: irwinhall_logsf, irwinhall_pdf, irwinhall_ppf, irwinhall_sf
    public :: irwinhall_loglikelihood, irwinhall_nnlf, irwinhall_score, irwinhall_fit
    public :: irwinhall_rvs, irwinhall_rvs_array
    public :: ksone_cdf, ksone_isf, ksone_logcdf, ksone_logpdf, ksone_logsf, ksone_pdf
    public :: ksone_ppf, ksone_sf, ksone_loglikelihood, ksone_nnlf, ksone_score, ksone_fit
    public :: ksone_rvs, ksone_rvs_array
    public :: kstwo_cdf, kstwo_isf, kstwo_logcdf, kstwo_logpdf, kstwo_logsf, kstwo_pdf
    public :: kstwo_ppf, kstwo_sf, kstwo_loglikelihood, kstwo_nnlf, kstwo_score, kstwo_fit
    public :: kstwo_rvs, kstwo_rvs_array
    public :: ncx2_cdf, ncx2_isf, ncx2_logcdf, ncx2_logpdf, ncx2_logsf, ncx2_pdf
    public :: ncx2_ppf, ncx2_sf, ncx2_loglikelihood, ncx2_nnlf, ncx2_score, ncx2_fit
    public :: ncx2_rvs, ncx2_rvs_array
    public :: ncf_cdf, ncf_isf, ncf_logcdf, ncf_logpdf, ncf_logsf, ncf_pdf, ncf_ppf, ncf_sf
    public :: ncf_loglikelihood, ncf_nnlf, ncf_score, ncf_fit, ncf_rvs, ncf_rvs_array
    public :: randint_cdf, randint_isf, randint_logcdf, randint_logpmf, randint_logsf
    public :: randint_pmf, randint_ppf, randint_sf, randint_loglikelihood, randint_nnlf
    public :: randint_score, randint_fit, randint_rvs, randint_rvs_array
    public :: planck_cdf, planck_isf, planck_logcdf, planck_logpmf, planck_logsf
    public :: planck_pmf, planck_ppf, planck_sf, planck_loglikelihood, planck_nnlf
    public :: planck_score, planck_fit, planck_rvs, planck_rvs_array
    public :: dlaplace_cdf, dlaplace_isf, dlaplace_logcdf, dlaplace_logpmf, dlaplace_logsf
    public :: dlaplace_pmf, dlaplace_ppf, dlaplace_sf, dlaplace_loglikelihood, dlaplace_nnlf
    public :: dlaplace_score, dlaplace_fit, dlaplace_rvs, dlaplace_rvs_array
    public :: logser_cdf, logser_isf, logser_logcdf, logser_logpmf, logser_logsf
    public :: logser_pmf, logser_ppf, logser_sf, logser_loglikelihood, logser_nnlf
    public :: logser_score, logser_fit, logser_rvs, logser_rvs_array
    public :: betabinom_cdf, betabinom_isf, betabinom_logcdf, betabinom_logpmf
    public :: betabinom_logsf, betabinom_pmf, betabinom_ppf, betabinom_sf
    public :: betabinom_loglikelihood, betabinom_nnlf, betabinom_score, betabinom_fit
    public :: betabinom_rvs, betabinom_rvs_array
    public :: hypergeom_cdf, hypergeom_isf, hypergeom_logcdf, hypergeom_logpmf
    public :: hypergeom_logsf, hypergeom_pmf, hypergeom_ppf, hypergeom_sf
    public :: hypergeom_loglikelihood, hypergeom_nnlf, hypergeom_score, hypergeom_fit
    public :: hypergeom_rvs, hypergeom_rvs_array
    public :: nhypergeom_cdf, nhypergeom_isf, nhypergeom_logcdf, nhypergeom_logpmf
    public :: nhypergeom_logsf, nhypergeom_pmf, nhypergeom_ppf, nhypergeom_sf
    public :: nhypergeom_loglikelihood, nhypergeom_nnlf, nhypergeom_score, nhypergeom_fit
    public :: nhypergeom_rvs, nhypergeom_rvs_array
    public :: boltzmann_cdf, boltzmann_isf, boltzmann_logcdf, boltzmann_logpmf
    public :: boltzmann_logsf, boltzmann_pmf, boltzmann_ppf, boltzmann_sf
    public :: boltzmann_loglikelihood, boltzmann_nnlf, boltzmann_score, boltzmann_fit
    public :: boltzmann_rvs, boltzmann_rvs_array
    public :: betanbinom_cdf, betanbinom_isf, betanbinom_logcdf, betanbinom_logpmf
    public :: betanbinom_logsf, betanbinom_pmf, betanbinom_ppf, betanbinom_sf
    public :: betanbinom_loglikelihood, betanbinom_nnlf, betanbinom_score, betanbinom_fit
    public :: betanbinom_rvs, betanbinom_rvs_array
    public :: yulesimon_cdf, yulesimon_isf, yulesimon_logcdf, yulesimon_logpmf
    public :: yulesimon_logsf, yulesimon_pmf, yulesimon_ppf, yulesimon_sf
    public :: yulesimon_loglikelihood, yulesimon_nnlf, yulesimon_score, yulesimon_fit
    public :: yulesimon_rvs, yulesimon_rvs_array
    public :: zipf_cdf, zipf_isf, zipf_logcdf, zipf_logpmf, zipf_logsf, zipf_pmf, zipf_ppf, zipf_sf
    public :: zipf_loglikelihood, zipf_nnlf, zipf_score, zipf_fit, zipf_rvs, zipf_rvs_array
    public :: zipfian_cdf, zipfian_isf, zipfian_logcdf, zipfian_logpmf
    public :: zipfian_logsf, zipfian_pmf, zipfian_ppf, zipfian_sf
    public :: zipfian_loglikelihood, zipfian_nnlf, zipfian_score, zipfian_fit
    public :: zipfian_rvs, zipfian_rvs_array
    public :: geninvgauss_cdf, geninvgauss_isf, geninvgauss_logcdf, geninvgauss_logpdf
    public :: geninvgauss_logsf, geninvgauss_pdf, geninvgauss_ppf, geninvgauss_sf
    public :: geninvgauss_loglikelihood, geninvgauss_nnlf, geninvgauss_score, geninvgauss_fit
    public :: geninvgauss_rvs, geninvgauss_rvs_array
    public :: norminvgauss_cdf, norminvgauss_isf, norminvgauss_logcdf, norminvgauss_logpdf
    public :: norminvgauss_logsf, norminvgauss_pdf, norminvgauss_ppf, norminvgauss_sf
    public :: norminvgauss_loglikelihood, norminvgauss_nnlf, norminvgauss_score, norminvgauss_fit
    public :: norminvgauss_rvs, norminvgauss_rvs_array
    public :: skellam_cdf, skellam_isf, skellam_logcdf, skellam_logpmf, skellam_logsf
    public :: skellam_pmf, skellam_ppf, skellam_sf, skellam_loglikelihood, skellam_nnlf
    public :: skellam_score, skellam_fit, skellam_rvs, skellam_rvs_array
    public :: genhyperbolic_cdf, genhyperbolic_isf, genhyperbolic_logcdf, genhyperbolic_logpdf
    public :: genhyperbolic_logsf, genhyperbolic_pdf, genhyperbolic_ppf, genhyperbolic_sf
    public :: genhyperbolic_loglikelihood, genhyperbolic_nnlf, genhyperbolic_score, genhyperbolic_fit
    public :: genhyperbolic_rvs, genhyperbolic_rvs_array
    public :: nchypergeom_fisher_cdf, nchypergeom_fisher_isf, nchypergeom_fisher_logcdf
    public :: nchypergeom_fisher_logpmf, nchypergeom_fisher_logsf, nchypergeom_fisher_pmf
    public :: nchypergeom_fisher_ppf, nchypergeom_fisher_sf, nchypergeom_fisher_loglikelihood
    public :: nchypergeom_fisher_nnlf, nchypergeom_fisher_score, nchypergeom_fisher_fit
    public :: nchypergeom_fisher_rvs, nchypergeom_fisher_rvs_array
    public :: nct_cdf, nct_isf, nct_logcdf, nct_logpdf, nct_logsf, nct_pdf, nct_ppf, nct_sf
    public :: nct_loglikelihood, nct_nnlf, nct_score, nct_fit, nct_rvs, nct_rvs_array
    public :: gausshyper_cdf, gausshyper_isf, gausshyper_logcdf, gausshyper_logpdf
    public :: gausshyper_logsf, gausshyper_pdf, gausshyper_ppf, gausshyper_sf
    public :: gausshyper_loglikelihood, gausshyper_nnlf, gausshyper_score, gausshyper_fit
    public :: gausshyper_rvs, gausshyper_rvs_array
    public :: landau_cdf, landau_isf, landau_logcdf, landau_logpdf
    public :: landau_logsf, landau_pdf, landau_ppf, landau_sf
    public :: landau_loglikelihood, landau_nnlf, landau_score, landau_fit
    public :: landau_rvs, landau_rvs_array
    public :: nchypergeom_wallenius_cdf, nchypergeom_wallenius_isf, nchypergeom_wallenius_logcdf
    public :: nchypergeom_wallenius_logpmf, nchypergeom_wallenius_logsf, nchypergeom_wallenius_pmf
    public :: nchypergeom_wallenius_ppf, nchypergeom_wallenius_sf
    public :: nchypergeom_wallenius_loglikelihood, nchypergeom_wallenius_nnlf
    public :: nchypergeom_wallenius_score, nchypergeom_wallenius_fit
    public :: nchypergeom_wallenius_rvs, nchypergeom_wallenius_rvs_array
    public :: poisson_binom_cdf, poisson_binom_isf, poisson_binom_logcdf, poisson_binom_logpmf
    public :: poisson_binom_logsf, poisson_binom_pmf, poisson_binom_ppf, poisson_binom_sf
    public :: poisson_binom_loglikelihood, poisson_binom_nnlf, poisson_binom_score, poisson_binom_fit
    public :: poisson_binom_rvs, poisson_binom_rvs_array

    public :: bernoulli_cdf
    public :: bernoulli_isf
    public :: bernoulli_logcdf
    public :: bernoulli_logpmf
    public :: bernoulli_logsf
    public :: bernoulli_pmf
    public :: bernoulli_ppf
    public :: bernoulli_sf
    public :: beta_cdf
    public :: beta_isf
    public :: beta_logcdf
    public :: beta_logpdf
    public :: beta_logsf
    public :: beta_pdf
    public :: beta_ppf
    public :: beta_sf
    public :: binomial_cdf
    public :: binomial_isf
    public :: binomial_logcdf
    public :: binomial_logpmf
    public :: binomial_logsf
    public :: binomial_pmf
    public :: binomial_ppf
    public :: binomial_sf
    public :: cauchy_cdf
    public :: cauchy_isf
    public :: cauchy_logcdf
    public :: cauchy_logpdf
    public :: cauchy_logsf
    public :: cauchy_pdf
    public :: cauchy_ppf
    public :: cauchy_sf
    public :: central_moment
    public :: covariance
    public :: mean
    public :: median
    public :: pearson_correlation
    public :: quantile
    public :: rankdata
    public :: standard_deviation
    public :: variance
    public :: chi2_cdf
    public :: chi2_isf
    public :: chi2_logcdf
    public :: chi2_logpdf
    public :: chi2_logsf
    public :: chi2_pdf
    public :: chi2_ppf
    public :: chi2_sf
    public :: dp
    public :: fit_result
    public :: fit_status_success
    public :: fit_status_max_iter
    public :: fit_status_invalid_input
    public :: fit_status_no_finite_objective
    public :: normal_fit
    public :: uniform_fit
    public :: exponential_fit
    public :: laplace_fit
    public :: logistic_fit
    public :: cauchy_fit
    public :: rayleigh_fit
    public :: gamma_fit
    public :: chi2_fit
    public :: t_fit
    public :: lognormal_fit
    public :: weibull_fit
    public :: pareto_fit
    public :: beta_fit
    public :: f_fit
    public :: bernoulli_fit
    public :: poisson_fit
    public :: geometric_fit
    public :: binomial_fit
    public :: negative_binomial_fit
    public :: gumbel_r_fit
    public :: gumbel_l_fit
    public :: powerlaw_fit
    public :: triang_fit
    public :: genpareto_fit
    public :: arcsine_fit
    public :: halfnorm_fit
    public :: halfcauchy_fit
    public :: lomax_fit
    public :: chi_fit
    public :: maxwell_fit
    public :: cosine_fit
    public :: semicircular_fit
    public :: anglit_fit
    public :: moyal_fit
    public :: hypsecant_fit
    public :: halflogistic_fit
    public :: invgamma_fit
    public :: invgauss_fit
    public :: levy_fit
    public :: loglaplace_fit
    public :: normal_score
    public :: uniform_score
    public :: exponential_score
    public :: laplace_score
    public :: logistic_score
    public :: cauchy_score
    public :: rayleigh_score
    public :: gamma_score
    public :: chi2_score
    public :: t_score
    public :: lognormal_score
    public :: weibull_score
    public :: pareto_score
    public :: beta_score
    public :: f_score
    public :: bernoulli_score
    public :: poisson_score
    public :: geometric_score
    public :: binomial_score
    public :: negative_binomial_score
    public :: gumbel_r_score
    public :: gumbel_l_score
    public :: powerlaw_score
    public :: triang_score
    public :: genpareto_score
    public :: arcsine_score
    public :: halfnorm_score
    public :: halfcauchy_score
    public :: lomax_score
    public :: chi_score
    public :: maxwell_score
    public :: cosine_score
    public :: semicircular_score
    public :: anglit_score
    public :: moyal_score
    public :: hypsecant_score
    public :: halflogistic_score
    public :: invgamma_score
    public :: invgauss_score
    public :: levy_score
    public :: loglaplace_score
    public :: normal_loglikelihood
    public :: normal_nnlf
    public :: uniform_loglikelihood
    public :: uniform_nnlf
    public :: exponential_loglikelihood
    public :: exponential_nnlf
    public :: laplace_loglikelihood
    public :: laplace_nnlf
    public :: logistic_loglikelihood
    public :: logistic_nnlf
    public :: cauchy_loglikelihood
    public :: cauchy_nnlf
    public :: rayleigh_loglikelihood
    public :: rayleigh_nnlf
    public :: gamma_loglikelihood
    public :: gamma_nnlf
    public :: chi2_loglikelihood
    public :: chi2_nnlf
    public :: t_loglikelihood
    public :: t_nnlf
    public :: lognormal_loglikelihood
    public :: lognormal_nnlf
    public :: weibull_loglikelihood
    public :: weibull_nnlf
    public :: pareto_loglikelihood
    public :: pareto_nnlf
    public :: beta_loglikelihood
    public :: beta_nnlf
    public :: f_loglikelihood
    public :: f_nnlf
    public :: bernoulli_loglikelihood
    public :: bernoulli_nnlf
    public :: poisson_loglikelihood
    public :: poisson_nnlf
    public :: geometric_loglikelihood
    public :: geometric_nnlf
    public :: binomial_loglikelihood
    public :: binomial_nnlf
    public :: negative_binomial_loglikelihood
    public :: negative_binomial_nnlf
    public :: gumbel_r_loglikelihood
    public :: gumbel_r_nnlf
    public :: gumbel_l_loglikelihood
    public :: gumbel_l_nnlf
    public :: powerlaw_loglikelihood
    public :: powerlaw_nnlf
    public :: triang_loglikelihood
    public :: triang_nnlf
    public :: genpareto_loglikelihood
    public :: genpareto_nnlf
    public :: arcsine_loglikelihood
    public :: arcsine_nnlf
    public :: halfnorm_loglikelihood
    public :: halfnorm_nnlf
    public :: halfcauchy_loglikelihood
    public :: halfcauchy_nnlf
    public :: lomax_loglikelihood
    public :: lomax_nnlf
    public :: chi_loglikelihood
    public :: chi_nnlf
    public :: maxwell_loglikelihood
    public :: maxwell_nnlf
    public :: cosine_loglikelihood
    public :: cosine_nnlf
    public :: semicircular_loglikelihood
    public :: semicircular_nnlf
    public :: anglit_loglikelihood
    public :: anglit_nnlf
    public :: moyal_loglikelihood
    public :: moyal_nnlf
    public :: hypsecant_loglikelihood
    public :: hypsecant_nnlf
    public :: halflogistic_loglikelihood
    public :: halflogistic_nnlf
    public :: invgamma_loglikelihood
    public :: invgamma_nnlf
    public :: invgauss_loglikelihood
    public :: invgauss_nnlf
    public :: levy_loglikelihood
    public :: levy_nnlf
    public :: loglaplace_loglikelihood
    public :: loglaplace_nnlf
    public :: invgamma_cdf
    public :: invgamma_isf
    public :: invgamma_logcdf
    public :: invgamma_logpdf
    public :: invgamma_logsf
    public :: invgamma_pdf
    public :: invgamma_ppf
    public :: invgamma_sf
    public :: invgauss_cdf
    public :: invgauss_isf
    public :: invgauss_logcdf
    public :: invgauss_logpdf
    public :: invgauss_logsf
    public :: invgauss_pdf
    public :: invgauss_ppf
    public :: invgauss_sf
    public :: levy_cdf
    public :: levy_isf
    public :: levy_logcdf
    public :: levy_logpdf
    public :: levy_logsf
    public :: levy_pdf
    public :: levy_ppf
    public :: levy_sf
    public :: loglaplace_cdf
    public :: loglaplace_isf
    public :: loglaplace_logcdf
    public :: loglaplace_logpdf
    public :: loglaplace_logsf
    public :: loglaplace_pdf
    public :: loglaplace_ppf
    public :: loglaplace_sf
    public :: exponential_cdf
    public :: exponential_isf
    public :: exponential_logcdf
    public :: exponential_logpdf
    public :: exponential_logsf
    public :: exponential_pdf
    public :: exponential_ppf
    public :: exponential_sf
    public :: f_cdf
    public :: f_isf
    public :: f_logcdf
    public :: f_logpdf
    public :: f_logsf
    public :: f_pdf
    public :: f_ppf
    public :: f_sf
    public :: gamma_cdf
    public :: gamma_isf
    public :: gamma_logcdf
    public :: gamma_logpdf
    public :: gamma_logsf
    public :: gamma_pdf
    public :: gamma_ppf
    public :: gamma_sf
    public :: geometric_cdf
    public :: geometric_isf
    public :: geometric_logcdf
    public :: geometric_logpmf
    public :: geometric_logsf
    public :: geometric_pmf
    public :: geometric_ppf
    public :: geometric_sf
    public :: correlation_test_result
    public :: mannwhitneyu
    public :: mannwhitneyu_result
    public :: pearsonr
    public :: spearmanr
    public :: ttest_1samp
    public :: ttest_ind
    public :: ttest_rel
    public :: ttest_result
    public :: bartlett
    public :: bartlett_result
    public :: chi2_contingency
    public :: chi2_contingency_result
    public :: chisquare
    public :: contingency_expected_freq
    public :: fisher_exact
    public :: fisher_exact_result
    public :: fligner
    public :: fligner_result
    public :: friedmanchisquare
    public :: friedmanchisquare_result
    public :: f_oneway
    public :: f_oneway_result
    public :: kruskal
    public :: kruskal_result
    public :: levene
    public :: levene_result
    public :: make_sample_group
    public :: power_divergence
    public :: power_divergence_named
    public :: power_divergence_result
    public :: ranksums
    public :: ranksums_result
    public :: sample_group
    public :: wilcoxon
    public :: wilcoxon_result
    public :: anderson, anderson_result
    public :: cramervonmises, cramervonmises_2samp
    public :: jarque_bera, ks_1samp, ks_2samp, kstest, kstest_result
    public :: kurtosistest, normaltest, shapiro, shapiro_result
    public :: anderson_ksamp, anderson_ksamp_result, ansari, epps_singleton_2samp
    public :: median_test, median_test_named, median_test_result, mood
    public :: alexandergovern, alexandergovern_result
    public :: dunnett, dunnett_result
    public :: matrix_confidence_interval, vector_confidence_interval
    public :: poisson_means_test, poisson_means_test_result
    public :: tukey_hsd, tukey_hsd_result
    public :: association_result, brunnermunzel, kendalltau
    public :: linregress, linregress_result, page_trend_result, page_trend_test
    public :: pointbiserialr, siegelslopes, siegelslopes_result
    public :: theilslopes, theilslopes_result
    public :: association, barnard_exact, binomtest, binomtest_result, boschloo_exact
    public :: combine_pvalues, combined_pvalue_result, confidence_interval, exact_2x2_result
    public :: false_discovery_control, margins, odds_ratio, odds_ratio_result
    public :: relative_risk, relative_risk_result
    public :: significance_result, skewtest
    public :: laplace_cdf
    public :: laplace_isf
    public :: laplace_logcdf
    public :: laplace_logpdf
    public :: laplace_logsf
    public :: laplace_pdf
    public :: laplace_ppf
    public :: laplace_sf
    public :: logistic_cdf
    public :: logistic_isf
    public :: logistic_logcdf
    public :: logistic_logpdf
    public :: logistic_logsf
    public :: logistic_pdf
    public :: logistic_ppf
    public :: logistic_sf
    public :: lognormal_cdf
    public :: lognormal_isf
    public :: lognormal_logcdf
    public :: lognormal_logpdf
    public :: lognormal_logsf
    public :: lognormal_pdf
    public :: lognormal_ppf
    public :: lognormal_sf
    public :: negative_binomial_cdf
    public :: negative_binomial_isf
    public :: negative_binomial_logcdf
    public :: negative_binomial_logpmf
    public :: negative_binomial_logsf
    public :: negative_binomial_pmf
    public :: negative_binomial_ppf
    public :: negative_binomial_sf
    public :: normal_rvs
    public :: normal_rvs_array
    public :: uniform_rvs
    public :: uniform_rvs_array
    public :: exponential_rvs
    public :: exponential_rvs_array
    public :: laplace_rvs
    public :: laplace_rvs_array
    public :: logistic_rvs
    public :: logistic_rvs_array
    public :: cauchy_rvs
    public :: cauchy_rvs_array
    public :: rayleigh_rvs
    public :: rayleigh_rvs_array
    public :: gamma_rvs
    public :: gamma_rvs_array
    public :: chi2_rvs
    public :: chi2_rvs_array
    public :: t_rvs
    public :: t_rvs_array
    public :: lognormal_rvs
    public :: lognormal_rvs_array
    public :: weibull_rvs
    public :: weibull_rvs_array
    public :: pareto_rvs
    public :: pareto_rvs_array
    public :: beta_rvs
    public :: beta_rvs_array
    public :: f_rvs
    public :: f_rvs_array
    public :: bernoulli_rvs
    public :: bernoulli_rvs_array
    public :: poisson_rvs
    public :: poisson_rvs_array
    public :: geometric_rvs
    public :: geometric_rvs_array
    public :: binomial_rvs
    public :: binomial_rvs_array
    public :: negative_binomial_rvs
    public :: negative_binomial_rvs_array
    public :: gumbel_r_rvs
    public :: gumbel_r_rvs_array
    public :: gumbel_l_rvs
    public :: gumbel_l_rvs_array
    public :: powerlaw_rvs
    public :: powerlaw_rvs_array
    public :: triang_rvs
    public :: triang_rvs_array
    public :: genpareto_rvs
    public :: genpareto_rvs_array
    public :: arcsine_rvs
    public :: arcsine_rvs_array
    public :: halfnorm_rvs
    public :: halfnorm_rvs_array
    public :: halfcauchy_rvs
    public :: halfcauchy_rvs_array
    public :: lomax_rvs
    public :: lomax_rvs_array
    public :: chi_rvs
    public :: chi_rvs_array
    public :: maxwell_rvs
    public :: maxwell_rvs_array
    public :: cosine_rvs
    public :: cosine_rvs_array
    public :: semicircular_rvs
    public :: semicircular_rvs_array
    public :: anglit_rvs
    public :: anglit_rvs_array
    public :: moyal_rvs
    public :: moyal_rvs_array
    public :: hypsecant_rvs
    public :: hypsecant_rvs_array
    public :: halflogistic_rvs
    public :: halflogistic_rvs_array
    public :: invgamma_rvs
    public :: invgamma_rvs_array
    public :: invgauss_rvs
    public :: invgauss_rvs_array
    public :: levy_rvs
    public :: levy_rvs_array
    public :: loglaplace_rvs
    public :: loglaplace_rvs_array
    public :: gumbel_r_cdf
    public :: gumbel_r_isf
    public :: gumbel_r_logcdf
    public :: gumbel_r_logpdf
    public :: gumbel_r_logsf
    public :: gumbel_r_pdf
    public :: gumbel_r_ppf
    public :: gumbel_r_sf
    public :: gumbel_l_cdf
    public :: gumbel_l_isf
    public :: gumbel_l_logcdf
    public :: gumbel_l_logpdf
    public :: gumbel_l_logsf
    public :: gumbel_l_pdf
    public :: gumbel_l_ppf
    public :: gumbel_l_sf
    public :: powerlaw_cdf
    public :: powerlaw_isf
    public :: powerlaw_logcdf
    public :: powerlaw_logpdf
    public :: powerlaw_logsf
    public :: powerlaw_pdf
    public :: powerlaw_ppf
    public :: powerlaw_sf
    public :: triang_cdf
    public :: triang_isf
    public :: triang_logcdf
    public :: triang_logpdf
    public :: triang_logsf
    public :: triang_pdf
    public :: triang_ppf
    public :: triang_sf
    public :: genpareto_cdf
    public :: genpareto_isf
    public :: genpareto_logcdf
    public :: genpareto_logpdf
    public :: genpareto_logsf
    public :: genpareto_pdf
    public :: genpareto_ppf
    public :: genpareto_sf
    public :: arcsine_cdf
    public :: arcsine_isf
    public :: arcsine_logcdf
    public :: arcsine_logpdf
    public :: arcsine_logsf
    public :: arcsine_pdf
    public :: arcsine_ppf
    public :: arcsine_sf
    public :: halfnorm_cdf
    public :: halfnorm_isf
    public :: halfnorm_logcdf
    public :: halfnorm_logpdf
    public :: halfnorm_logsf
    public :: halfnorm_pdf
    public :: halfnorm_ppf
    public :: halfnorm_sf
    public :: halfcauchy_cdf
    public :: halfcauchy_isf
    public :: halfcauchy_logcdf
    public :: halfcauchy_logpdf
    public :: halfcauchy_logsf
    public :: halfcauchy_pdf
    public :: halfcauchy_ppf
    public :: halfcauchy_sf
    public :: lomax_cdf
    public :: lomax_isf
    public :: lomax_logcdf
    public :: lomax_logpdf
    public :: lomax_logsf
    public :: lomax_pdf
    public :: lomax_ppf
    public :: lomax_sf
    public :: chi_cdf
    public :: chi_isf
    public :: chi_logcdf
    public :: chi_logpdf
    public :: chi_logsf
    public :: chi_pdf
    public :: chi_ppf
    public :: chi_sf
    public :: maxwell_cdf
    public :: maxwell_isf
    public :: maxwell_logcdf
    public :: maxwell_logpdf
    public :: maxwell_logsf
    public :: maxwell_pdf
    public :: maxwell_ppf
    public :: maxwell_sf
    public :: cosine_cdf
    public :: cosine_isf
    public :: cosine_logcdf
    public :: cosine_logpdf
    public :: cosine_logsf
    public :: cosine_pdf
    public :: cosine_ppf
    public :: cosine_sf
    public :: semicircular_cdf
    public :: semicircular_isf
    public :: semicircular_logcdf
    public :: semicircular_logpdf
    public :: semicircular_logsf
    public :: semicircular_pdf
    public :: semicircular_ppf
    public :: semicircular_sf
    public :: anglit_cdf
    public :: anglit_isf
    public :: anglit_logcdf
    public :: anglit_logpdf
    public :: anglit_logsf
    public :: anglit_pdf
    public :: anglit_ppf
    public :: anglit_sf
    public :: moyal_cdf
    public :: moyal_isf
    public :: moyal_logcdf
    public :: moyal_logpdf
    public :: moyal_logsf
    public :: moyal_pdf
    public :: moyal_ppf
    public :: moyal_sf
    public :: hypsecant_cdf
    public :: hypsecant_isf
    public :: hypsecant_logcdf
    public :: hypsecant_logpdf
    public :: hypsecant_logsf
    public :: hypsecant_pdf
    public :: hypsecant_ppf
    public :: hypsecant_sf
    public :: halflogistic_cdf
    public :: halflogistic_isf
    public :: halflogistic_logcdf
    public :: halflogistic_logpdf
    public :: halflogistic_logsf
    public :: halflogistic_pdf
    public :: halflogistic_ppf
    public :: halflogistic_sf
    public :: normal_cdf
    public :: normal_isf
    public :: normal_logcdf
    public :: normal_logpdf
    public :: normal_logsf
    public :: normal_pdf
    public :: normal_ppf
    public :: normal_sf
    public :: pareto_cdf
    public :: pareto_isf
    public :: pareto_logcdf
    public :: pareto_logpdf
    public :: pareto_logsf
    public :: pareto_pdf
    public :: pareto_ppf
    public :: pareto_sf
    public :: poisson_cdf
    public :: poisson_isf
    public :: poisson_logcdf
    public :: poisson_logpmf
    public :: poisson_logsf
    public :: poisson_pmf
    public :: poisson_ppf
    public :: poisson_sf
    public :: rayleigh_cdf
    public :: rayleigh_isf
    public :: rayleigh_logcdf
    public :: rayleigh_logpdf
    public :: rayleigh_logsf
    public :: rayleigh_pdf
    public :: rayleigh_ppf
    public :: rayleigh_sf
    public :: bootstrap
    public :: bootstrap_percentile
    public :: bootstrap_result
    public :: monte_carlo_test
    public :: monte_carlo_test_result
    public :: permutation_test
    public :: permutation_test_result
    public :: power
    public :: power_result
    public :: t_cdf
    public :: t_isf
    public :: t_logcdf
    public :: t_logpdf
    public :: t_logsf
    public :: t_pdf
    public :: t_ppf
    public :: t_sf
    public :: uniform_cdf
    public :: uniform_isf
    public :: uniform_logcdf
    public :: uniform_logpdf
    public :: uniform_logsf
    public :: uniform_pdf
    public :: uniform_ppf
    public :: uniform_sf
    public :: weibull_cdf
    public :: weibull_isf
    public :: weibull_logcdf
    public :: weibull_logpdf
    public :: weibull_logsf
    public :: weibull_pdf
    public :: weibull_ppf
    public :: weibull_sf

    public :: bradford_cdf
    public :: bradford_isf
    public :: bradford_logcdf
    public :: bradford_logpdf
    public :: bradford_logsf
    public :: bradford_pdf
    public :: bradford_ppf
    public :: bradford_sf
    public :: bradford_loglikelihood
    public :: bradford_nnlf
    public :: bradford_score
    public :: bradford_fit
    public :: bradford_rvs
    public :: bradford_rvs_array
    public :: truncexpon_cdf
    public :: truncexpon_isf
    public :: truncexpon_logcdf
    public :: truncexpon_logpdf
    public :: truncexpon_logsf
    public :: truncexpon_pdf
    public :: truncexpon_ppf
    public :: truncexpon_sf
    public :: truncexpon_loglikelihood
    public :: truncexpon_nnlf
    public :: truncexpon_score
    public :: truncexpon_fit
    public :: truncexpon_rvs
    public :: truncexpon_rvs_array
    public :: fisk_cdf
    public :: fisk_isf
    public :: fisk_logcdf
    public :: fisk_logpdf
    public :: fisk_logsf
    public :: fisk_pdf
    public :: fisk_ppf
    public :: fisk_sf
    public :: fisk_loglikelihood
    public :: fisk_nnlf
    public :: fisk_score
    public :: fisk_fit
    public :: fisk_rvs
    public :: fisk_rvs_array
    public :: dweibull_cdf
    public :: dweibull_isf
    public :: dweibull_logcdf
    public :: dweibull_logpdf
    public :: dweibull_logsf
    public :: dweibull_pdf
    public :: dweibull_ppf
    public :: dweibull_sf
    public :: dweibull_loglikelihood
    public :: dweibull_nnlf
    public :: dweibull_score
    public :: dweibull_fit
    public :: dweibull_rvs
    public :: dweibull_rvs_array

    public :: alpha_cdf
    public :: alpha_isf
    public :: alpha_logcdf
    public :: alpha_logpdf
    public :: alpha_logsf
    public :: alpha_pdf
    public :: alpha_ppf
    public :: alpha_sf
    public :: alpha_loglikelihood
    public :: alpha_nnlf
    public :: alpha_score
    public :: alpha_fit
    public :: alpha_rvs
    public :: alpha_rvs_array
    public :: fatiguelife_cdf
    public :: fatiguelife_isf
    public :: fatiguelife_logcdf
    public :: fatiguelife_logpdf
    public :: fatiguelife_logsf
    public :: fatiguelife_pdf
    public :: fatiguelife_ppf
    public :: fatiguelife_sf
    public :: fatiguelife_loglikelihood
    public :: fatiguelife_nnlf
    public :: fatiguelife_score
    public :: fatiguelife_fit
    public :: fatiguelife_rvs
    public :: fatiguelife_rvs_array
    public :: genlogistic_cdf
    public :: genlogistic_isf
    public :: genlogistic_logcdf
    public :: genlogistic_logpdf
    public :: genlogistic_logsf
    public :: genlogistic_pdf
    public :: genlogistic_ppf
    public :: genlogistic_sf
    public :: genlogistic_loglikelihood
    public :: genlogistic_nnlf
    public :: genlogistic_score
    public :: genlogistic_fit
    public :: genlogistic_rvs
    public :: genlogistic_rvs_array
    public :: gennorm_cdf
    public :: gennorm_isf
    public :: gennorm_logcdf
    public :: gennorm_logpdf
    public :: gennorm_logsf
    public :: gennorm_pdf
    public :: gennorm_ppf
    public :: gennorm_sf
    public :: gennorm_loglikelihood
    public :: gennorm_nnlf
    public :: gennorm_score
    public :: gennorm_fit
    public :: gennorm_rvs
    public :: gennorm_rvs_array
    public :: nakagami_cdf
    public :: nakagami_isf
    public :: nakagami_logcdf
    public :: nakagami_logpdf
    public :: nakagami_logsf
    public :: nakagami_pdf
    public :: nakagami_ppf
    public :: nakagami_sf
    public :: nakagami_loglikelihood
    public :: nakagami_nnlf
    public :: nakagami_score
    public :: nakagami_fit
    public :: nakagami_rvs
    public :: nakagami_rvs_array
    public :: powernorm_cdf
    public :: powernorm_isf
    public :: powernorm_logcdf
    public :: powernorm_logpdf
    public :: powernorm_logsf
    public :: powernorm_pdf
    public :: powernorm_ppf
    public :: powernorm_sf
    public :: powernorm_loglikelihood
    public :: powernorm_nnlf
    public :: powernorm_score
    public :: powernorm_fit
    public :: powernorm_rvs
    public :: powernorm_rvs_array
    public :: loggamma_cdf
    public :: loggamma_isf
    public :: loggamma_logcdf
    public :: loggamma_logpdf
    public :: loggamma_logsf
    public :: loggamma_pdf
    public :: loggamma_ppf
    public :: loggamma_sf
    public :: loggamma_loglikelihood
    public :: loggamma_nnlf
    public :: loggamma_score
    public :: loggamma_fit
    public :: loggamma_rvs
    public :: loggamma_rvs_array
    public :: wald_cdf
    public :: wald_isf
    public :: wald_logcdf
    public :: wald_logpdf
    public :: wald_logsf
    public :: wald_pdf
    public :: wald_ppf
    public :: wald_sf
    public :: wald_loglikelihood
    public :: wald_nnlf
    public :: wald_score
    public :: wald_fit
    public :: wald_rvs
    public :: wald_rvs_array
    public :: gompertz_cdf
    public :: gompertz_isf
    public :: gompertz_logcdf
    public :: gompertz_logpdf
    public :: gompertz_logsf
    public :: gompertz_pdf
    public :: gompertz_ppf
    public :: gompertz_sf
    public :: gompertz_loglikelihood
    public :: gompertz_nnlf
    public :: gompertz_score
    public :: gompertz_fit
    public :: gompertz_rvs
    public :: gompertz_rvs_array
    public :: invweibull_cdf
    public :: invweibull_isf
    public :: invweibull_logcdf
    public :: invweibull_logpdf
    public :: invweibull_logsf
    public :: invweibull_pdf
    public :: invweibull_ppf
    public :: invweibull_sf
    public :: invweibull_loglikelihood
    public :: invweibull_nnlf
    public :: invweibull_score
    public :: invweibull_fit
    public :: invweibull_rvs
    public :: invweibull_rvs_array
    public :: betaprime_cdf
    public :: betaprime_isf
    public :: betaprime_logcdf
    public :: betaprime_logpdf
    public :: betaprime_logsf
    public :: betaprime_pdf
    public :: betaprime_ppf
    public :: betaprime_sf
    public :: betaprime_loglikelihood
    public :: betaprime_nnlf
    public :: betaprime_score
    public :: betaprime_fit
    public :: betaprime_rvs
    public :: betaprime_rvs_array
    public :: burr12_cdf
    public :: burr12_isf
    public :: burr12_logcdf
    public :: burr12_logpdf
    public :: burr12_logsf
    public :: burr12_pdf
    public :: burr12_ppf
    public :: burr12_sf
    public :: burr12_loglikelihood
    public :: burr12_nnlf
    public :: burr12_score
    public :: burr12_fit
    public :: burr12_rvs
    public :: burr12_rvs_array
    public :: genhalflogistic_cdf
    public :: genhalflogistic_isf
    public :: genhalflogistic_logcdf
    public :: genhalflogistic_logpdf
    public :: genhalflogistic_logsf
    public :: genhalflogistic_pdf
    public :: genhalflogistic_ppf
    public :: genhalflogistic_sf
    public :: genhalflogistic_loglikelihood
    public :: genhalflogistic_nnlf
    public :: genhalflogistic_score
    public :: genhalflogistic_fit
    public :: genhalflogistic_rvs
    public :: genhalflogistic_rvs_array
    public :: exponpow_cdf
    public :: exponpow_isf
    public :: exponpow_logcdf
    public :: exponpow_logpdf
    public :: exponpow_logsf
    public :: exponpow_pdf
    public :: exponpow_ppf
    public :: exponpow_sf
    public :: exponpow_loglikelihood
    public :: exponpow_nnlf
    public :: exponpow_score
    public :: exponpow_fit
    public :: exponpow_rvs
    public :: exponpow_rvs_array
    public :: exponweib_cdf
    public :: exponweib_isf
    public :: exponweib_logcdf
    public :: exponweib_logpdf
    public :: exponweib_logsf
    public :: exponweib_pdf
    public :: exponweib_ppf
    public :: exponweib_sf
    public :: exponweib_loglikelihood
    public :: exponweib_nnlf
    public :: exponweib_score
    public :: exponweib_fit
    public :: exponweib_rvs
    public :: exponweib_rvs_array
    public :: powerlognorm_cdf
    public :: powerlognorm_isf
    public :: powerlognorm_logcdf
    public :: powerlognorm_logpdf
    public :: powerlognorm_logsf
    public :: powerlognorm_pdf
    public :: powerlognorm_ppf
    public :: powerlognorm_sf
    public :: powerlognorm_loglikelihood
    public :: powerlognorm_nnlf
    public :: powerlognorm_score
    public :: powerlognorm_fit
    public :: powerlognorm_rvs
    public :: powerlognorm_rvs_array
    public :: levy_l_cdf
    public :: levy_l_isf
    public :: levy_l_logcdf
    public :: levy_l_logpdf
    public :: levy_l_logsf
    public :: levy_l_pdf
    public :: levy_l_ppf
    public :: levy_l_sf
    public :: levy_l_loglikelihood
    public :: levy_l_nnlf
    public :: levy_l_score
    public :: levy_l_fit
    public :: levy_l_rvs
    public :: levy_l_rvs_array
    public :: weibull_max_cdf
    public :: weibull_max_isf
    public :: weibull_max_logcdf
    public :: weibull_max_logpdf
    public :: weibull_max_logsf
    public :: weibull_max_pdf
    public :: weibull_max_ppf
    public :: weibull_max_sf
    public :: weibull_max_loglikelihood
    public :: weibull_max_nnlf
    public :: weibull_max_score
    public :: weibull_max_fit
    public :: weibull_max_rvs
    public :: weibull_max_rvs_array
    public :: rdist_cdf
    public :: rdist_isf
    public :: rdist_logcdf
    public :: rdist_logpdf
    public :: rdist_logsf
    public :: rdist_pdf
    public :: rdist_ppf
    public :: rdist_sf
    public :: rdist_loglikelihood
    public :: rdist_nnlf
    public :: rdist_score
    public :: rdist_fit
    public :: rdist_rvs
    public :: rdist_rvs_array
    public :: skewcauchy_cdf
    public :: skewcauchy_isf
    public :: skewcauchy_logcdf
    public :: skewcauchy_logpdf
    public :: skewcauchy_logsf
    public :: skewcauchy_pdf
    public :: skewcauchy_ppf
    public :: skewcauchy_sf
    public :: skewcauchy_loglikelihood
    public :: skewcauchy_nnlf
    public :: skewcauchy_score
    public :: skewcauchy_fit
    public :: skewcauchy_rvs
    public :: skewcauchy_rvs_array
    public :: dgamma_cdf
    public :: dgamma_isf
    public :: dgamma_logcdf
    public :: dgamma_logpdf
    public :: dgamma_logsf
    public :: dgamma_pdf
    public :: dgamma_ppf
    public :: dgamma_sf
    public :: dgamma_loglikelihood
    public :: dgamma_nnlf
    public :: dgamma_score
    public :: dgamma_fit
    public :: dgamma_rvs
    public :: dgamma_rvs_array
    public :: laplace_asymmetric_cdf
    public :: laplace_asymmetric_isf
    public :: laplace_asymmetric_logcdf
    public :: laplace_asymmetric_logpdf
    public :: laplace_asymmetric_logsf
    public :: laplace_asymmetric_pdf
    public :: laplace_asymmetric_ppf
    public :: laplace_asymmetric_sf
    public :: laplace_asymmetric_loglikelihood
    public :: laplace_asymmetric_nnlf
    public :: laplace_asymmetric_score
    public :: laplace_asymmetric_fit
    public :: laplace_asymmetric_rvs
    public :: laplace_asymmetric_rvs_array
    public :: truncnorm_cdf
    public :: truncnorm_isf
    public :: truncnorm_logcdf
    public :: truncnorm_logpdf
    public :: truncnorm_logsf
    public :: truncnorm_pdf
    public :: truncnorm_ppf
    public :: truncnorm_sf
    public :: truncnorm_loglikelihood
    public :: truncnorm_nnlf
    public :: truncnorm_score
    public :: truncnorm_fit
    public :: truncnorm_rvs
    public :: truncnorm_rvs_array
    public :: loguniform_cdf
    public :: loguniform_isf
    public :: loguniform_logcdf
    public :: loguniform_logpdf
    public :: loguniform_logsf
    public :: loguniform_pdf
    public :: loguniform_ppf
    public :: loguniform_sf
    public :: loguniform_loglikelihood
    public :: loguniform_nnlf
    public :: loguniform_score
    public :: loguniform_fit
    public :: loguniform_rvs
    public :: loguniform_rvs_array
    public :: foldnorm_cdf, foldnorm_isf, foldnorm_logcdf, foldnorm_logpdf
    public :: foldnorm_logsf, foldnorm_pdf, foldnorm_ppf, foldnorm_sf
    public :: foldnorm_loglikelihood, foldnorm_nnlf, foldnorm_score, foldnorm_fit
    public :: foldnorm_rvs, foldnorm_rvs_array
    public :: foldcauchy_cdf, foldcauchy_isf, foldcauchy_logcdf, foldcauchy_logpdf
    public :: foldcauchy_logsf, foldcauchy_pdf, foldcauchy_ppf, foldcauchy_sf
    public :: foldcauchy_loglikelihood, foldcauchy_nnlf, foldcauchy_score, foldcauchy_fit
    public :: foldcauchy_rvs, foldcauchy_rvs_array
    public :: recipinvgauss_cdf, recipinvgauss_isf, recipinvgauss_logcdf
    public :: recipinvgauss_logpdf, recipinvgauss_logsf, recipinvgauss_pdf
    public :: recipinvgauss_ppf, recipinvgauss_sf, recipinvgauss_loglikelihood
    public :: recipinvgauss_nnlf, recipinvgauss_score, recipinvgauss_fit
    public :: recipinvgauss_rvs, recipinvgauss_rvs_array
    public :: truncpareto_cdf, truncpareto_isf, truncpareto_logcdf, truncpareto_logpdf
    public :: truncpareto_logsf, truncpareto_pdf, truncpareto_ppf, truncpareto_sf
    public :: truncpareto_loglikelihood, truncpareto_nnlf, truncpareto_score, truncpareto_fit
    public :: truncpareto_rvs, truncpareto_rvs_array
    public :: levy_stable_cdf, levy_stable_isf, levy_stable_logcdf, levy_stable_logpdf
    public :: levy_stable_logsf, levy_stable_pdf, levy_stable_ppf, levy_stable_sf
    public :: levy_stable_loglikelihood, levy_stable_nnlf, levy_stable_score, levy_stable_fit
    public :: levy_stable_rvs, levy_stable_rvs_array
    public :: studentized_range_cdf, studentized_range_isf, studentized_range_logcdf
    public :: studentized_range_logpdf, studentized_range_logsf, studentized_range_pdf
    public :: studentized_range_ppf, studentized_range_sf
    public :: studentized_range_loglikelihood, studentized_range_nnlf
    public :: studentized_range_score, studentized_range_fit
    public :: studentized_range_rvs, studentized_range_rvs_array
    public :: multivariate_normal_cdf, multivariate_normal_entropy, multivariate_normal_fit
    public :: multivariate_normal_logcdf, multivariate_normal_logpdf, multivariate_normal_marginal
    public :: multivariate_normal_pdf, multivariate_normal_rvs, multivariate_normal_rvs_array
    public :: multivariate_status_success, multivariate_status_invalid_shape
    public :: multivariate_status_invalid_parameter, multivariate_status_linalg_failure
    public :: gaussian_kde, gaussian_kde_init
    public :: gaussian_kde_status_success, gaussian_kde_status_invalid_shape
    public :: gaussian_kde_status_invalid_parameter, gaussian_kde_status_singular
    public :: gaussian_kde_status_linalg_failure
    public :: qmc_discrepancy, qmc_geometric_discrepancy, qmc_halton, qmc_halton_init
    public :: qmc_latin_hypercube, qmc_latin_hypercube_init, qmc_scale, qmc_to_integers
    public :: qmc_sobol, qmc_sobol_init, qmc_multinomial, qmc_multinomial_init
    public :: qmc_multivariate_normal, qmc_multivariate_normal_init
    public :: qmc_poisson_disk, qmc_poisson_disk_init
    public :: qmc_update_discrepancy, qmc_van_der_corput
    public :: qmc_status_success, qmc_status_invalid_shape, qmc_status_invalid_parameter
    public :: qmc_status_out_of_bounds
    public :: multivariate_t_cdf, multivariate_t_entropy, multivariate_t_logpdf, multivariate_t_marginal
    public :: multivariate_t_pdf, multivariate_t_rvs, multivariate_t_rvs_array
    public :: multivariate_t_status_success, multivariate_t_status_invalid_shape
    public :: multivariate_t_status_invalid_parameter, multivariate_t_status_linalg_failure
    public :: dirichlet_cov, dirichlet_entropy, dirichlet_logpdf, dirichlet_mean
    public :: dirichlet_pdf, dirichlet_rvs, dirichlet_rvs_array, dirichlet_var
    public :: dirichlet_status_success, dirichlet_status_invalid_parameter, dirichlet_status_invalid_shape
    public :: multinomial_cov, multinomial_entropy, multinomial_logpmf, multinomial_mean
    public :: multinomial_pmf, multinomial_rvs, multinomial_rvs_array
    public :: multinomial_status_success, multinomial_status_invalid_parameter, multinomial_status_invalid_shape
    public :: dirichlet_multinomial_cov, dirichlet_multinomial_logpmf
    public :: dirichlet_multinomial_mean, dirichlet_multinomial_pmf, dirichlet_multinomial_var
    public :: multivariate_hypergeom_cov, multivariate_hypergeom_logpmf, multivariate_hypergeom_mean
    public :: multivariate_hypergeom_pmf, multivariate_hypergeom_rvs, multivariate_hypergeom_rvs_array
    public :: multivariate_hypergeom_var, multivariate_hypergeom_status_success
    public :: multivariate_hypergeom_status_invalid_parameter, multivariate_hypergeom_status_invalid_shape
    public :: normal_inverse_gamma_logpdf, normal_inverse_gamma_mean, normal_inverse_gamma_pdf
    public :: normal_inverse_gamma_rvs, normal_inverse_gamma_rvs_array, normal_inverse_gamma_var
    public :: matrix_normal_entropy, matrix_normal_logpdf, matrix_normal_pdf, matrix_normal_rvs
    public :: matrix_normal_rvs_array, matrix_normal_status_success, matrix_normal_status_invalid_shape
    public :: matrix_normal_status_invalid_parameter, matrix_normal_status_linalg_failure
    public :: wishart_entropy, wishart_logpdf, wishart_mean, wishart_mode, wishart_pdf
    public :: wishart_rvs, wishart_rvs_array, wishart_var, wishart_status_success
    public :: wishart_status_invalid_shape, wishart_status_invalid_parameter, wishart_status_linalg_failure
    public :: invwishart_entropy, invwishart_logpdf, invwishart_mean, invwishart_mode, invwishart_pdf
    public :: invwishart_rvs, invwishart_rvs_array, invwishart_var, invwishart_status_success
    public :: invwishart_status_invalid_shape, invwishart_status_invalid_parameter
    public :: invwishart_status_linalg_failure
    public :: matrix_t_logpdf, matrix_t_pdf, matrix_t_rvs, matrix_t_rvs_array
    public :: matrix_t_status_success, matrix_t_status_invalid_shape, matrix_t_status_invalid_parameter
    public :: matrix_t_status_linalg_failure
    public :: uniform_direction_rvs, uniform_direction_rvs_array
    public :: uniform_direction_status_success, uniform_direction_status_invalid_shape
    public :: uniform_direction_status_numerical_failure
    public :: vonmises_fisher_entropy, vonmises_fisher_fit, vonmises_fisher_logpdf
    public :: vonmises_fisher_pdf, vonmises_fisher_rvs, vonmises_fisher_rvs_array
    public :: vonmises_fisher_status_success, vonmises_fisher_status_invalid_shape
    public :: vonmises_fisher_status_invalid_parameter, vonmises_fisher_status_numerical_failure
    public :: ortho_group_rvs, ortho_group_rvs_array, ortho_group_status_success
    public :: ortho_group_status_invalid_shape, ortho_group_status_numerical_failure
    public :: special_ortho_group_rvs, special_ortho_group_rvs_array
    public :: special_ortho_group_status_success, special_ortho_group_status_invalid_shape
    public :: special_ortho_group_status_numerical_failure
    public :: unitary_group_rvs, unitary_group_rvs_array, unitary_group_status_success
    public :: unitary_group_status_invalid_shape, unitary_group_status_numerical_failure
    public :: random_correlation_rvs, random_correlation_status_success
    public :: random_correlation_status_invalid_shape, random_correlation_status_invalid_parameter
    public :: random_correlation_status_numerical_failure
    public :: random_table_logpmf, random_table_mean, random_table_pmf, random_table_rvs
    public :: random_table_rvs_array, random_table_status_success, random_table_status_invalid_shape
    public :: random_table_status_invalid_parameter, random_table_status_numerical_failure
end module scifort_stats
