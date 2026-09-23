! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_size_t
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_normal, only : normal_pdf, normal_cdf, normal_ppf
    use scifort_uniform, only : uniform_pdf, uniform_cdf, uniform_ppf
    use scifort_exponential, only : exponential_pdf, exponential_cdf, exponential_ppf
    use scifort_laplace, only : laplace_pdf, laplace_cdf, laplace_ppf
    use scifort_logistic, only : logistic_pdf, logistic_cdf, logistic_ppf
    use scifort_cauchy, only : cauchy_pdf, cauchy_cdf, cauchy_ppf
    use scifort_rayleigh, only : rayleigh_pdf, rayleigh_cdf, rayleigh_ppf
    use scifort_gamma, only : gamma_pdf, gamma_cdf, gamma_ppf
    use scifort_chi2, only : chi2_pdf, chi2_cdf, chi2_ppf
    use scifort_student_t, only : t_pdf, t_cdf, t_ppf
    use scifort_lognormal, only : lognormal_pdf, lognormal_cdf, lognormal_ppf
    use scifort_weibull, only : weibull_pdf, weibull_cdf, weibull_ppf
    use scifort_pareto, only : pareto_pdf, pareto_cdf, pareto_ppf
    use scifort_beta, only : beta_pdf, beta_cdf, beta_ppf
    use scifort_f_distribution, only : f_pdf, f_cdf, f_ppf
    use scifort_gumbel_r, only : gumbel_r_pdf, gumbel_r_cdf, gumbel_r_ppf
    use scifort_gumbel_l, only : gumbel_l_pdf, gumbel_l_cdf, gumbel_l_ppf
    use scifort_powerlaw, only : powerlaw_pdf, powerlaw_cdf, powerlaw_ppf
    use scifort_triang, only : triang_pdf, triang_cdf, triang_ppf
    use scifort_genpareto, only : genpareto_pdf, genpareto_cdf, genpareto_ppf
    use scifort_arcsine, only : arcsine_pdf, arcsine_cdf, arcsine_ppf
    use scifort_halfnorm, only : halfnorm_pdf, halfnorm_cdf, halfnorm_ppf
    use scifort_halfcauchy, only : halfcauchy_pdf, halfcauchy_cdf, halfcauchy_ppf
    use scifort_lomax, only : lomax_pdf, lomax_cdf, lomax_ppf
    use scifort_chi, only : chi_pdf, chi_cdf, chi_ppf
    use scifort_maxwell, only : maxwell_pdf, maxwell_cdf, maxwell_ppf
    use scifort_cosine, only : cosine_pdf, cosine_cdf, cosine_ppf
    use scifort_semicircular, only : semicircular_pdf, semicircular_cdf, semicircular_ppf
    use scifort_anglit, only : anglit_pdf, anglit_cdf, anglit_ppf
    use scifort_moyal, only : moyal_pdf, moyal_cdf, moyal_ppf
    use scifort_hypsecant, only : hypsecant_pdf, hypsecant_cdf, hypsecant_ppf
    use scifort_halflogistic, only : halflogistic_pdf, halflogistic_cdf, halflogistic_ppf
    use scifort_invgamma, only : invgamma_pdf, invgamma_cdf, invgamma_ppf
    use scifort_invgauss, only : invgauss_pdf, invgauss_cdf, invgauss_ppf
    use scifort_levy, only : levy_pdf, levy_cdf, levy_ppf
    use scifort_loglaplace, only : loglaplace_pdf, loglaplace_cdf, loglaplace_ppf
    use scifort_bradford, only : bradford_pdf, bradford_cdf, bradford_ppf
    use scifort_truncexpon, only : truncexpon_pdf, truncexpon_cdf, truncexpon_ppf
    use scifort_fisk, only : fisk_pdf, fisk_cdf, fisk_ppf
    use scifort_dweibull, only : dweibull_pdf, dweibull_cdf, dweibull_ppf
    use scifort_alpha, only : alpha_pdf, alpha_cdf, alpha_ppf
    use scifort_fatiguelife, only : fatiguelife_pdf, fatiguelife_cdf, fatiguelife_ppf
    use scifort_genlogistic, only : genlogistic_pdf, genlogistic_cdf, genlogistic_ppf
    use scifort_gennorm, only : gennorm_pdf, gennorm_cdf, gennorm_ppf
    use scifort_nakagami, only : nakagami_pdf, nakagami_cdf, nakagami_ppf
    use scifort_powernorm, only : powernorm_pdf, powernorm_cdf, powernorm_ppf
    use scifort_loggamma, only : loggamma_pdf, loggamma_cdf, loggamma_ppf
    use scifort_wald, only : wald_pdf, wald_cdf, wald_ppf
    use scifort_gompertz, only : gompertz_pdf, gompertz_cdf, gompertz_ppf
    use scifort_invweibull, only : invweibull_pdf, invweibull_cdf, invweibull_ppf
    use scifort_betaprime, only : betaprime_pdf, betaprime_cdf, betaprime_ppf
    use scifort_burr12, only : burr12_pdf, burr12_cdf, burr12_ppf
    use scifort_genhalflogistic, only : genhalflogistic_pdf, genhalflogistic_cdf, &
        genhalflogistic_ppf
    use scifort_exponpow, only : exponpow_pdf, exponpow_cdf, exponpow_ppf
    use scifort_exponweib, only : exponweib_pdf, exponweib_cdf, exponweib_ppf
    use scifort_powerlognorm, only : powerlognorm_pdf, powerlognorm_cdf, powerlognorm_ppf
    use scifort_levy_l, only : levy_l_pdf, levy_l_cdf, levy_l_ppf
    use scifort_weibull_max, only : weibull_max_pdf, weibull_max_cdf, weibull_max_ppf
    use scifort_rdist, only : rdist_pdf, rdist_cdf, rdist_ppf
    use scifort_skewcauchy, only : skewcauchy_pdf, skewcauchy_cdf, skewcauchy_ppf
    use scifort_dgamma, only : dgamma_pdf, dgamma_cdf, dgamma_ppf
    use scifort_laplace_asymmetric, only : laplace_asymmetric_pdf, &
        laplace_asymmetric_cdf, laplace_asymmetric_ppf
    use scifort_truncnorm, only : truncnorm_pdf, truncnorm_cdf, truncnorm_ppf
    use scifort_loguniform, only : loguniform_pdf, loguniform_cdf, loguniform_ppf
    use scifort_foldnorm, only : foldnorm_pdf, foldnorm_cdf, foldnorm_ppf
    use scifort_foldcauchy, only : foldcauchy_pdf, foldcauchy_cdf, foldcauchy_ppf
    use scifort_recipinvgauss, only : recipinvgauss_pdf, recipinvgauss_cdf, recipinvgauss_ppf
    use scifort_truncpareto, only : truncpareto_pdf, truncpareto_cdf, truncpareto_ppf
    use scifort_exponnorm, only : exponnorm_pdf, exponnorm_cdf, exponnorm_ppf
    use scifort_johnsonsb, only : johnsonsb_pdf, johnsonsb_cdf, johnsonsb_ppf
    use scifort_johnsonsu, only : johnsonsu_pdf, johnsonsu_cdf, johnsonsu_ppf
    use scifort_trapezoid, only : trapezoid_pdf, trapezoid_cdf, trapezoid_ppf
    use scifort_burr, only : burr_pdf, burr_cdf, burr_ppf
    use scifort_mielke, only : mielke_pdf, mielke_cdf, mielke_ppf
    use scifort_gibrat, only : gibrat_pdf, gibrat_cdf, gibrat_ppf
    use scifort_wrapcauchy, only : wrapcauchy_pdf, wrapcauchy_cdf, wrapcauchy_ppf
    use scifort_genextreme, only : genextreme_pdf, genextreme_cdf, genextreme_ppf
    use scifort_kappa3, only : kappa3_pdf, kappa3_cdf, kappa3_ppf
    use scifort_kappa4, only : kappa4_pdf, kappa4_cdf, kappa4_ppf
    use scifort_truncweibull_min, only : truncweibull_min_pdf, truncweibull_min_cdf, &
        truncweibull_min_ppf
    use scifort_gengamma, only : gengamma_pdf, gengamma_cdf, gengamma_ppf
    use scifort_halfgennorm, only : halfgennorm_pdf, halfgennorm_cdf, halfgennorm_ppf
    use scifort_argus, only : argus_pdf, argus_cdf, argus_ppf
    use scifort_erlang, only : erlang_pdf, erlang_cdf, erlang_ppf
    use scifort_crystalball, only : crystalball_pdf, crystalball_cdf, crystalball_ppf
    use scifort_jf_skew_t, only : jf_skew_t_pdf, jf_skew_t_cdf, jf_skew_t_ppf
    use scifort_pearson3, only : pearson3_pdf, pearson3_cdf, pearson3_ppf
    use scifort_rel_breitwigner, only : rel_breitwigner_pdf, rel_breitwigner_cdf, &
        rel_breitwigner_ppf
    use scifort_genexpon, only : genexpon_pdf, genexpon_cdf, genexpon_ppf
    use scifort_skewnorm, only : skewnorm_pdf, skewnorm_cdf, skewnorm_ppf
    use scifort_tukeylambda, only : tukeylambda_pdf, tukeylambda_cdf, tukeylambda_ppf
    use scifort_rice, only : rice_pdf, rice_cdf, rice_ppf
    use scifort_dpareto_lognorm, only : dpareto_lognorm_pdf, dpareto_lognorm_cdf, &
        dpareto_lognorm_ppf
    use scifort_vonmises, only : vonmises_pdf, vonmises_cdf, vonmises_ppf
    use scifort_vonmises_line, only : vonmises_line_pdf, vonmises_line_cdf, vonmises_line_ppf
    use scifort_kstwobign, only : kstwobign_pdf, kstwobign_cdf, kstwobign_ppf
    use scifort_irwinhall, only : irwinhall_pdf, irwinhall_cdf, irwinhall_ppf
    use scifort_ksone, only : ksone_pdf, ksone_cdf, ksone_ppf
    use scifort_kstwo, only : kstwo_pdf, kstwo_cdf, kstwo_ppf
    use scifort_levy_stable, only : levy_stable_pdf, levy_stable_cdf, levy_stable_ppf
    use scifort_studentized_range, only : studentized_range_pdf, studentized_range_cdf, &
        studentized_range_ppf
    use scifort_ncx2, only : ncx2_pdf, ncx2_cdf, ncx2_ppf
    use scifort_ncf, only : ncf_pdf, ncf_cdf, ncf_ppf
    use scifort_randint, only : randint_pmf, randint_cdf, randint_ppf
    use scifort_planck, only : planck_pmf, planck_cdf, planck_ppf
    use scifort_dlaplace, only : dlaplace_pmf, dlaplace_cdf, dlaplace_ppf
    use scifort_logser, only : logser_pmf, logser_cdf, logser_ppf
    use scifort_betabinom, only : betabinom_pmf, betabinom_cdf, betabinom_ppf
    use scifort_hypergeom, only : hypergeom_pmf, hypergeom_cdf, hypergeom_ppf
    use scifort_nhypergeom, only : nhypergeom_pmf, nhypergeom_cdf, nhypergeom_ppf
    use scifort_boltzmann, only : boltzmann_pmf, boltzmann_cdf, boltzmann_ppf
    use scifort_betanbinom, only : betanbinom_pmf, betanbinom_cdf, betanbinom_ppf
    use scifort_yulesimon, only : yulesimon_pmf, yulesimon_cdf, yulesimon_ppf
    use scifort_zipf, only : zipf_pmf, zipf_cdf, zipf_ppf
    use scifort_zipfian, only : zipfian_pmf, zipfian_cdf, zipfian_ppf
    use scifort_geninvgauss, only : geninvgauss_pdf, geninvgauss_cdf, geninvgauss_ppf
    use scifort_norminvgauss, only : norminvgauss_pdf, norminvgauss_cdf, norminvgauss_ppf
    use scifort_skellam, only : skellam_pmf, skellam_cdf, skellam_ppf
    use scifort_genhyperbolic, only : genhyperbolic_pdf, genhyperbolic_cdf, genhyperbolic_ppf
    use scifort_nchypergeom_fisher, only : nchypergeom_fisher_pmf, nchypergeom_fisher_cdf, &
        nchypergeom_fisher_ppf
    use scifort_nchypergeom_wallenius, only : nchypergeom_wallenius_pmf, &
        nchypergeom_wallenius_cdf, nchypergeom_wallenius_ppf
    use scifort_poisson_binom, only : poisson_binom_pmf, poisson_binom_cdf, poisson_binom_ppf
    use scifort_nct, only : nct_pdf, nct_cdf, nct_ppf
    use scifort_gausshyper, only : gausshyper_pdf, gausshyper_cdf, gausshyper_ppf
    use scifort_landau, only : landau_pdf, landau_cdf, landau_ppf
    use scifort_bernoulli, only : bernoulli_pmf, bernoulli_cdf, bernoulli_ppf
    use scifort_poisson, only : poisson_pmf, poisson_cdf, poisson_ppf
    use scifort_geometric, only : geometric_pmf, geometric_cdf, geometric_ppf
    use scifort_binomial, only : binomial_pmf, binomial_cdf, binomial_ppf
    use scifort_negative_binomial, only : negative_binomial_cdf, &
        negative_binomial_pmf, negative_binomial_ppf
    use scifort_kinds, only : dp
    implicit none
    private

    integer(c_int), parameter :: status_ok = 0_c_int
    integer(c_int), parameter :: status_invalid_argument = 1_c_int


    public :: scifort_exponnorm_pdf_f64, scifort_exponnorm_cdf_f64, scifort_exponnorm_ppf_f64
    public :: scifort_johnsonsb_pdf_f64, scifort_johnsonsb_cdf_f64, scifort_johnsonsb_ppf_f64
    public :: scifort_johnsonsu_pdf_f64, scifort_johnsonsu_cdf_f64, scifort_johnsonsu_ppf_f64
    public :: scifort_trapezoid_pdf_f64, scifort_trapezoid_cdf_f64, scifort_trapezoid_ppf_f64
    public :: scifort_burr_pdf_f64, scifort_burr_cdf_f64, scifort_burr_ppf_f64
    public :: scifort_mielke_pdf_f64, scifort_mielke_cdf_f64, scifort_mielke_ppf_f64
    public :: scifort_gibrat_pdf_f64, scifort_gibrat_cdf_f64, scifort_gibrat_ppf_f64
    public :: scifort_wrapcauchy_pdf_f64, scifort_wrapcauchy_cdf_f64
    public :: scifort_wrapcauchy_ppf_f64
    public :: scifort_genextreme_pdf_f64, scifort_genextreme_cdf_f64, scifort_genextreme_ppf_f64
    public :: scifort_kappa3_pdf_f64, scifort_kappa3_cdf_f64, scifort_kappa3_ppf_f64
    public :: scifort_kappa4_pdf_f64, scifort_kappa4_cdf_f64, scifort_kappa4_ppf_f64
    public :: scifort_truncweibull_min_pdf_f64, scifort_truncweibull_min_cdf_f64
    public :: scifort_truncweibull_min_ppf_f64
    public :: scifort_gengamma_pdf_f64, scifort_gengamma_cdf_f64, scifort_gengamma_ppf_f64
    public :: scifort_halfgennorm_pdf_f64, scifort_halfgennorm_cdf_f64
    public :: scifort_halfgennorm_ppf_f64
    public :: scifort_argus_pdf_f64, scifort_argus_cdf_f64, scifort_argus_ppf_f64
    public :: scifort_erlang_pdf_f64, scifort_erlang_cdf_f64, scifort_erlang_ppf_f64
    public :: scifort_crystalball_pdf_f64, scifort_crystalball_cdf_f64
    public :: scifort_crystalball_ppf_f64
    public :: scifort_jf_skew_t_pdf_f64, scifort_jf_skew_t_cdf_f64, scifort_jf_skew_t_ppf_f64
    public :: scifort_pearson3_pdf_f64, scifort_pearson3_cdf_f64, scifort_pearson3_ppf_f64
    public :: scifort_rel_breitwigner_pdf_f64, scifort_rel_breitwigner_cdf_f64
    public :: scifort_rel_breitwigner_ppf_f64
    public :: scifort_genexpon_pdf_f64, scifort_genexpon_cdf_f64, scifort_genexpon_ppf_f64
    public :: scifort_skewnorm_pdf_f64, scifort_skewnorm_cdf_f64, scifort_skewnorm_ppf_f64
    public :: scifort_tukeylambda_pdf_f64, scifort_tukeylambda_cdf_f64
    public :: scifort_tukeylambda_ppf_f64
    public :: scifort_rice_pdf_f64, scifort_rice_cdf_f64, scifort_rice_ppf_f64
    public :: scifort_dpareto_lognorm_pdf_f64, scifort_dpareto_lognorm_cdf_f64
    public :: scifort_dpareto_lognorm_ppf_f64
    public :: scifort_vonmises_pdf_f64, scifort_vonmises_cdf_f64, scifort_vonmises_ppf_f64
    public :: scifort_vonmises_line_pdf_f64, scifort_vonmises_line_cdf_f64
    public :: scifort_vonmises_line_ppf_f64
    public :: scifort_kstwobign_pdf_f64, scifort_kstwobign_cdf_f64, scifort_kstwobign_ppf_f64
    public :: scifort_irwinhall_pdf_f64, scifort_irwinhall_cdf_f64, scifort_irwinhall_ppf_f64
    public :: scifort_ksone_pdf_f64, scifort_ksone_cdf_f64, scifort_ksone_ppf_f64
    public :: scifort_kstwo_pdf_f64, scifort_kstwo_cdf_f64, scifort_kstwo_ppf_f64
    public :: scifort_levy_stable_pdf_f64, scifort_levy_stable_cdf_f64
    public :: scifort_levy_stable_ppf_f64
    public :: scifort_studentized_range_pdf_f64, scifort_studentized_range_cdf_f64
    public :: scifort_studentized_range_ppf_f64
    public :: scifort_ncx2_pdf_f64, scifort_ncx2_cdf_f64, scifort_ncx2_ppf_f64
    public :: scifort_ncf_pdf_f64, scifort_ncf_cdf_f64, scifort_ncf_ppf_f64
    public :: scifort_randint_pmf_f64, scifort_randint_cdf_f64, scifort_randint_ppf_f64
    public :: scifort_planck_pmf_f64, scifort_planck_cdf_f64, scifort_planck_ppf_f64
    public :: scifort_dlaplace_pmf_f64, scifort_dlaplace_cdf_f64, scifort_dlaplace_ppf_f64
    public :: scifort_logser_pmf_f64, scifort_logser_cdf_f64, scifort_logser_ppf_f64
    public :: scifort_betabinom_pmf_f64, scifort_betabinom_cdf_f64, scifort_betabinom_ppf_f64
    public :: scifort_hypergeom_pmf_f64, scifort_hypergeom_cdf_f64, scifort_hypergeom_ppf_f64
    public :: scifort_nhypergeom_pmf_f64, scifort_nhypergeom_cdf_f64, scifort_nhypergeom_ppf_f64
    public :: scifort_boltzmann_pmf_f64, scifort_boltzmann_cdf_f64, scifort_boltzmann_ppf_f64
    public :: scifort_betanbinom_pmf_f64, scifort_betanbinom_cdf_f64, scifort_betanbinom_ppf_f64
    public :: scifort_yulesimon_pmf_f64, scifort_yulesimon_cdf_f64, scifort_yulesimon_ppf_f64
    public :: scifort_zipf_pmf_f64, scifort_zipf_cdf_f64, scifort_zipf_ppf_f64
    public :: scifort_zipfian_pmf_f64, scifort_zipfian_cdf_f64, scifort_zipfian_ppf_f64
    public :: scifort_geninvgauss_pdf_f64, scifort_geninvgauss_cdf_f64, scifort_geninvgauss_ppf_f64
    public :: scifort_norminvgauss_pdf_f64, scifort_norminvgauss_cdf_f64, scifort_norminvgauss_ppf_f64
    public :: scifort_skellam_pmf_f64, scifort_skellam_cdf_f64, scifort_skellam_ppf_f64
    public :: scifort_genhyperbolic_pdf_f64, scifort_genhyperbolic_cdf_f64, scifort_genhyperbolic_ppf_f64
    public :: scifort_nchypergeom_fisher_pmf_f64, scifort_nchypergeom_fisher_cdf_f64
    public :: scifort_nchypergeom_fisher_ppf_f64
    public :: scifort_nct_pdf_f64, scifort_nct_cdf_f64, scifort_nct_ppf_f64
    public :: scifort_gausshyper_pdf_f64, scifort_gausshyper_cdf_f64, scifort_gausshyper_ppf_f64
    public :: scifort_landau_pdf_f64, scifort_landau_cdf_f64, scifort_landau_ppf_f64
    public :: scifort_nchypergeom_wallenius_pmf_f64, scifort_nchypergeom_wallenius_cdf_f64
    public :: scifort_nchypergeom_wallenius_ppf_f64
    public :: scifort_poisson_binom_pmf_f64, scifort_poisson_binom_cdf_f64
    public :: scifort_poisson_binom_ppf_f64

    public :: scifort_bernoulli_cdf_f64
    public :: scifort_bernoulli_pmf_f64
    public :: scifort_bernoulli_ppf_f64
    public :: scifort_beta_cdf_f64
    public :: scifort_beta_pdf_f64
    public :: scifort_beta_ppf_f64
    public :: scifort_binomial_cdf_f64
    public :: scifort_binomial_pmf_f64
    public :: scifort_binomial_ppf_f64
    public :: scifort_cauchy_cdf_f64
    public :: scifort_cauchy_pdf_f64
    public :: scifort_cauchy_ppf_f64
    public :: scifort_chi2_cdf_f64
    public :: scifort_chi2_pdf_f64
    public :: scifort_chi2_ppf_f64
    public :: scifort_exponential_cdf_f64
    public :: scifort_exponential_pdf_f64
    public :: scifort_exponential_ppf_f64
    public :: scifort_f_cdf_f64
    public :: scifort_f_pdf_f64
    public :: scifort_f_ppf_f64
    public :: scifort_gumbel_r_pdf_f64
    public :: scifort_gumbel_r_cdf_f64
    public :: scifort_gumbel_r_ppf_f64
    public :: scifort_gumbel_l_pdf_f64
    public :: scifort_gumbel_l_cdf_f64
    public :: scifort_gumbel_l_ppf_f64
    public :: scifort_powerlaw_pdf_f64
    public :: scifort_powerlaw_cdf_f64
    public :: scifort_powerlaw_ppf_f64
    public :: scifort_triang_pdf_f64
    public :: scifort_triang_cdf_f64
    public :: scifort_triang_ppf_f64
    public :: scifort_genpareto_pdf_f64
    public :: scifort_genpareto_cdf_f64
    public :: scifort_genpareto_ppf_f64
    public :: scifort_arcsine_pdf_f64
    public :: scifort_arcsine_cdf_f64
    public :: scifort_arcsine_ppf_f64
    public :: scifort_halfnorm_pdf_f64
    public :: scifort_halfnorm_cdf_f64
    public :: scifort_halfnorm_ppf_f64
    public :: scifort_halfcauchy_pdf_f64
    public :: scifort_halfcauchy_cdf_f64
    public :: scifort_halfcauchy_ppf_f64
    public :: scifort_lomax_pdf_f64
    public :: scifort_lomax_cdf_f64
    public :: scifort_lomax_ppf_f64
    public :: scifort_chi_pdf_f64
    public :: scifort_chi_cdf_f64
    public :: scifort_chi_ppf_f64
    public :: scifort_maxwell_pdf_f64
    public :: scifort_maxwell_cdf_f64
    public :: scifort_maxwell_ppf_f64
    public :: scifort_cosine_pdf_f64
    public :: scifort_cosine_cdf_f64
    public :: scifort_cosine_ppf_f64
    public :: scifort_semicircular_pdf_f64
    public :: scifort_semicircular_cdf_f64
    public :: scifort_semicircular_ppf_f64
    public :: scifort_anglit_pdf_f64
    public :: scifort_anglit_cdf_f64
    public :: scifort_anglit_ppf_f64
    public :: scifort_moyal_pdf_f64
    public :: scifort_moyal_cdf_f64
    public :: scifort_moyal_ppf_f64
    public :: scifort_hypsecant_pdf_f64
    public :: scifort_hypsecant_cdf_f64
    public :: scifort_hypsecant_ppf_f64
    public :: scifort_halflogistic_pdf_f64
    public :: scifort_halflogistic_cdf_f64
    public :: scifort_halflogistic_ppf_f64
    public :: scifort_invgamma_pdf_f64
    public :: scifort_invgamma_cdf_f64
    public :: scifort_invgamma_ppf_f64
    public :: scifort_invgauss_pdf_f64
    public :: scifort_invgauss_cdf_f64
    public :: scifort_invgauss_ppf_f64
    public :: scifort_levy_pdf_f64
    public :: scifort_levy_cdf_f64
    public :: scifort_levy_ppf_f64
    public :: scifort_loglaplace_pdf_f64
    public :: scifort_loglaplace_cdf_f64
    public :: scifort_loglaplace_ppf_f64
    public :: scifort_bradford_pdf_f64
    public :: scifort_bradford_cdf_f64
    public :: scifort_bradford_ppf_f64
    public :: scifort_truncexpon_pdf_f64
    public :: scifort_truncexpon_cdf_f64
    public :: scifort_truncexpon_ppf_f64
    public :: scifort_fisk_pdf_f64
    public :: scifort_fisk_cdf_f64
    public :: scifort_fisk_ppf_f64
    public :: scifort_dweibull_pdf_f64
    public :: scifort_dweibull_cdf_f64
    public :: scifort_dweibull_ppf_f64
    public :: scifort_alpha_pdf_f64
    public :: scifort_alpha_cdf_f64
    public :: scifort_alpha_ppf_f64
    public :: scifort_fatiguelife_pdf_f64
    public :: scifort_fatiguelife_cdf_f64
    public :: scifort_fatiguelife_ppf_f64
    public :: scifort_genlogistic_pdf_f64
    public :: scifort_genlogistic_cdf_f64
    public :: scifort_genlogistic_ppf_f64
    public :: scifort_gennorm_pdf_f64
    public :: scifort_gennorm_cdf_f64
    public :: scifort_gennorm_ppf_f64
    public :: scifort_nakagami_pdf_f64
    public :: scifort_nakagami_cdf_f64
    public :: scifort_nakagami_ppf_f64
    public :: scifort_powernorm_pdf_f64
    public :: scifort_powernorm_cdf_f64
    public :: scifort_powernorm_ppf_f64
    public :: scifort_loggamma_pdf_f64
    public :: scifort_loggamma_cdf_f64
    public :: scifort_loggamma_ppf_f64
    public :: scifort_wald_pdf_f64
    public :: scifort_wald_cdf_f64
    public :: scifort_wald_ppf_f64
    public :: scifort_gompertz_pdf_f64
    public :: scifort_gompertz_cdf_f64
    public :: scifort_gompertz_ppf_f64
    public :: scifort_invweibull_pdf_f64
    public :: scifort_invweibull_cdf_f64
    public :: scifort_invweibull_ppf_f64
    public :: scifort_betaprime_pdf_f64
    public :: scifort_betaprime_cdf_f64
    public :: scifort_betaprime_ppf_f64
    public :: scifort_burr12_pdf_f64
    public :: scifort_burr12_cdf_f64
    public :: scifort_burr12_ppf_f64
    public :: scifort_genhalflogistic_pdf_f64
    public :: scifort_genhalflogistic_cdf_f64
    public :: scifort_genhalflogistic_ppf_f64
    public :: scifort_exponpow_pdf_f64
    public :: scifort_exponpow_cdf_f64
    public :: scifort_exponpow_ppf_f64
    public :: scifort_exponweib_pdf_f64
    public :: scifort_exponweib_cdf_f64
    public :: scifort_exponweib_ppf_f64
    public :: scifort_powerlognorm_pdf_f64
    public :: scifort_powerlognorm_cdf_f64
    public :: scifort_powerlognorm_ppf_f64
    public :: scifort_gamma_cdf_f64
    public :: scifort_gamma_pdf_f64
    public :: scifort_gamma_ppf_f64
    public :: scifort_geometric_cdf_f64
    public :: scifort_geometric_pmf_f64
    public :: scifort_geometric_ppf_f64
    public :: scifort_laplace_cdf_f64
    public :: scifort_laplace_pdf_f64
    public :: scifort_laplace_ppf_f64
    public :: scifort_logistic_cdf_f64
    public :: scifort_logistic_pdf_f64
    public :: scifort_logistic_ppf_f64
    public :: scifort_lognormal_cdf_f64
    public :: scifort_lognormal_pdf_f64
    public :: scifort_lognormal_ppf_f64
    public :: scifort_negative_binomial_cdf_f64
    public :: scifort_negative_binomial_pmf_f64
    public :: scifort_negative_binomial_ppf_f64
    public :: scifort_normal_cdf_f64
    public :: scifort_normal_cdf_vec_f64
    public :: scifort_normal_pdf_f64
    public :: scifort_normal_pdf_vec_f64
    public :: scifort_normal_ppf_f64
    public :: scifort_normal_ppf_vec_f64
    public :: scifort_pareto_cdf_f64
    public :: scifort_pareto_pdf_f64
    public :: scifort_pareto_ppf_f64
    public :: scifort_poisson_cdf_f64
    public :: scifort_poisson_pmf_f64
    public :: scifort_poisson_ppf_f64
    public :: scifort_rayleigh_cdf_f64
    public :: scifort_rayleigh_pdf_f64
    public :: scifort_rayleigh_ppf_f64
    public :: scifort_t_cdf_f64
    public :: scifort_t_pdf_f64
    public :: scifort_t_ppf_f64
    public :: scifort_uniform_cdf_f64
    public :: scifort_uniform_pdf_f64
    public :: scifort_uniform_ppf_f64
    public :: scifort_version_major
    public :: scifort_version_minor
    public :: scifort_version_patch
    public :: scifort_weibull_cdf_f64
    public :: scifort_weibull_pdf_f64
    public :: scifort_weibull_ppf_f64
    public :: scifort_levy_l_pdf_f64
    public :: scifort_levy_l_cdf_f64
    public :: scifort_levy_l_ppf_f64
    public :: scifort_weibull_max_pdf_f64
    public :: scifort_weibull_max_cdf_f64
    public :: scifort_weibull_max_ppf_f64
    public :: scifort_rdist_pdf_f64
    public :: scifort_rdist_cdf_f64
    public :: scifort_rdist_ppf_f64
    public :: scifort_skewcauchy_pdf_f64
    public :: scifort_skewcauchy_cdf_f64
    public :: scifort_skewcauchy_ppf_f64
    public :: scifort_dgamma_pdf_f64
    public :: scifort_dgamma_cdf_f64
    public :: scifort_dgamma_ppf_f64
    public :: scifort_laplace_asymmetric_pdf_f64
    public :: scifort_laplace_asymmetric_cdf_f64
    public :: scifort_laplace_asymmetric_ppf_f64
    public :: scifort_truncnorm_pdf_f64
    public :: scifort_truncnorm_cdf_f64
    public :: scifort_truncnorm_ppf_f64
    public :: scifort_loguniform_pdf_f64
    public :: scifort_loguniform_cdf_f64
    public :: scifort_loguniform_ppf_f64

    public :: scifort_foldnorm_pdf_f64, scifort_foldnorm_cdf_f64, scifort_foldnorm_ppf_f64
    public :: scifort_foldcauchy_pdf_f64, scifort_foldcauchy_cdf_f64, scifort_foldcauchy_ppf_f64
    public :: scifort_recipinvgauss_pdf_f64, scifort_recipinvgauss_cdf_f64
    public :: scifort_recipinvgauss_ppf_f64
    public :: scifort_truncpareto_pdf_f64, scifort_truncpareto_cdf_f64, scifort_truncpareto_ppf_f64

contains

    function scifort_version_major() result(version) bind(c, name="scifort_version_major")
        integer(c_int) :: version

        version = 0_c_int
    end function scifort_version_major

    function scifort_version_minor() result(version) bind(c, name="scifort_version_minor")
        integer(c_int) :: version

        version = 1_c_int
    end function scifort_version_minor

    function scifort_version_patch() result(version) bind(c, name="scifort_version_patch")
        integer(c_int) :: version

        version = 0_c_int
    end function scifort_version_patch

    function scifort_normal_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_normal_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: y

        y = real(normal_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_normal_pdf_f64

    function scifort_normal_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_normal_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: y

        y = real(normal_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_normal_cdf_f64

    function scifort_normal_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_normal_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: x

        x = real(normal_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_normal_ppf_f64

    subroutine scifort_normal_pdf_vec_f64(n, x, loc, scale, y, status) &
            bind(c, name="scifort_normal_pdf_vec_f64")
        integer(c_size_t), value, intent(in) :: n !! number of elements of x and y
        real(c_double), intent(in) :: x(*) !! points of evaluation
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double), intent(out) :: y(*) !! densities; unchanged if status is nonzero
        integer(c_int), intent(out) :: status !! 0 on success, 1 for an invalid loc or scale

        integer(c_size_t) :: i

        if (.not. valid_c_loc_scale(loc, scale)) then
            status = status_invalid_argument
            return
        end if

        do i = 1_c_size_t, n
            y(i) = real(normal_pdf(real(x(i), dp), real(loc, dp), &
                real(scale, dp)), c_double)
        end do
        status = status_ok
    end subroutine scifort_normal_pdf_vec_f64

    subroutine scifort_normal_cdf_vec_f64(n, x, loc, scale, y, status) &
            bind(c, name="scifort_normal_cdf_vec_f64")
        integer(c_size_t), value, intent(in) :: n !! number of elements of x and y
        real(c_double), intent(in) :: x(*) !! points of evaluation
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double), intent(out) :: y(*) !! lower-tail probabilities; unchanged if status /= 0
        integer(c_int), intent(out) :: status !! 0 on success, 1 for an invalid loc or scale

        integer(c_size_t) :: i

        if (.not. valid_c_loc_scale(loc, scale)) then
            status = status_invalid_argument
            return
        end if

        do i = 1_c_size_t, n
            y(i) = real(normal_cdf(real(x(i), dp), real(loc, dp), &
                real(scale, dp)), c_double)
        end do
        status = status_ok
    end subroutine scifort_normal_cdf_vec_f64

    subroutine scifort_normal_ppf_vec_f64(n, p, loc, scale, x, status) &
            bind(c, name="scifort_normal_ppf_vec_f64")
        integer(c_size_t), value, intent(in) :: n !! number of elements of p and x
        real(c_double), intent(in) :: p(*) !! lower-tail probabilities
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double), intent(out) :: x(*) !! quantiles; unchanged if status is nonzero
        integer(c_int), intent(out) :: status !! 0 on success, 1 for an invalid loc or scale

        integer(c_size_t) :: i

        if (.not. valid_c_loc_scale(loc, scale)) then
            status = status_invalid_argument
            return
        end if

        do i = 1_c_size_t, n
            x(i) = real(normal_ppf(real(p(i), dp), real(loc, dp), &
                real(scale, dp)), c_double)
        end do
        status = status_ok
    end subroutine scifort_normal_ppf_vec_f64


    function scifort_uniform_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_uniform_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(uniform_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_uniform_pdf_f64

    function scifort_uniform_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_uniform_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(uniform_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_uniform_cdf_f64

    function scifort_uniform_ppf_f64(probability, loc, scale) result(x) &
            bind(c, name="scifort_uniform_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(uniform_ppf(real(probability, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_uniform_ppf_f64

    function scifort_exponential_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_exponential_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(exponential_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponential_pdf_f64

    function scifort_exponential_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_exponential_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(exponential_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponential_cdf_f64

    function scifort_exponential_ppf_f64(probability, loc, scale) result(x) &
            bind(c, name="scifort_exponential_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(exponential_ppf(real(probability, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponential_ppf_f64

    function scifort_laplace_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_laplace_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(laplace_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_laplace_pdf_f64

    function scifort_laplace_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_laplace_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(laplace_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_laplace_cdf_f64

    function scifort_laplace_ppf_f64(probability, loc, scale) result(x) &
            bind(c, name="scifort_laplace_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(laplace_ppf(real(probability, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_laplace_ppf_f64

    function scifort_logistic_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_logistic_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(logistic_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_logistic_pdf_f64

    function scifort_logistic_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_logistic_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(logistic_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_logistic_cdf_f64

    function scifort_logistic_ppf_f64(probability, loc, scale) result(x) &
            bind(c, name="scifort_logistic_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(logistic_ppf(real(probability, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_logistic_ppf_f64

    function scifort_cauchy_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_cauchy_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(cauchy_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_cauchy_pdf_f64

    function scifort_cauchy_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_cauchy_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(cauchy_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_cauchy_cdf_f64

    function scifort_cauchy_ppf_f64(probability, loc, scale) result(x) &
            bind(c, name="scifort_cauchy_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(cauchy_ppf(real(probability, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_cauchy_ppf_f64

    function scifort_rayleigh_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_rayleigh_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(rayleigh_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rayleigh_pdf_f64

    function scifort_rayleigh_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_rayleigh_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(rayleigh_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rayleigh_cdf_f64

    function scifort_rayleigh_ppf_f64(probability, loc, scale) result(x) &
            bind(c, name="scifort_rayleigh_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(rayleigh_ppf(real(probability, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rayleigh_ppf_f64

    function scifort_gamma_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_gamma_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: a !! first positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(gamma_pdf(real(x, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gamma_pdf_f64

    function scifort_gamma_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_gamma_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: a !! first positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(gamma_cdf(real(x, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gamma_cdf_f64

    function scifort_gamma_ppf_f64(probability, a, loc, scale) result(x) &
            bind(c, name="scifort_gamma_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: a !! first positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(gamma_ppf(real(probability, dp), real(a, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_gamma_ppf_f64

    function scifort_chi2_pdf_f64(x, df, loc, scale) result(y) &
            bind(c, name="scifort_chi2_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(chi2_pdf(real(x, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_chi2_pdf_f64

    function scifort_chi2_cdf_f64(x, df, loc, scale) result(y) &
            bind(c, name="scifort_chi2_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(chi2_cdf(real(x, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_chi2_cdf_f64

    function scifort_chi2_ppf_f64(probability, df, loc, scale) result(x) &
            bind(c, name="scifort_chi2_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(chi2_ppf(real(probability, dp), real(df, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_chi2_ppf_f64

    function scifort_t_pdf_f64(x, df, loc, scale) result(y) &
            bind(c, name="scifort_t_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(t_pdf(real(x, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_t_pdf_f64

    function scifort_t_cdf_f64(x, df, loc, scale) result(y) &
            bind(c, name="scifort_t_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(t_cdf(real(x, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_t_cdf_f64

    function scifort_t_ppf_f64(probability, df, loc, scale) result(x) &
            bind(c, name="scifort_t_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(t_ppf(real(probability, dp), real(df, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_t_ppf_f64

    function scifort_lognormal_pdf_f64(x, s, loc, scale) result(y) &
            bind(c, name="scifort_lognormal_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: s !! positive lognormal shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(lognormal_pdf(real(x, dp), real(s, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_lognormal_pdf_f64

    function scifort_lognormal_cdf_f64(x, s, loc, scale) result(y) &
            bind(c, name="scifort_lognormal_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: s !! positive lognormal shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(lognormal_cdf(real(x, dp), real(s, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_lognormal_cdf_f64

    function scifort_lognormal_ppf_f64(probability, s, loc, scale) result(x) &
            bind(c, name="scifort_lognormal_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: s !! positive lognormal shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(lognormal_ppf(real(probability, dp), real(s, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_lognormal_ppf_f64

    function scifort_weibull_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_weibull_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: c !! positive Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(weibull_pdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_weibull_pdf_f64

    function scifort_weibull_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_weibull_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: c !! positive Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(weibull_cdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_weibull_cdf_f64

    function scifort_weibull_ppf_f64(probability, c, loc, scale) result(x) &
            bind(c, name="scifort_weibull_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(weibull_ppf(real(probability, dp), real(c, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_weibull_ppf_f64

    function scifort_pareto_pdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_pareto_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: b !! second positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(pareto_pdf(real(x, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_pareto_pdf_f64

    function scifort_pareto_cdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_pareto_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: b !! second positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(pareto_cdf(real(x, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_pareto_cdf_f64

    function scifort_pareto_ppf_f64(probability, b, loc, scale) result(x) &
            bind(c, name="scifort_pareto_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: b !! second positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(pareto_ppf(real(probability, dp), real(b, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_pareto_ppf_f64

    function scifort_beta_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_beta_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: a !! first positive shape parameter
        real(c_double), value, intent(in) :: b !! second positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(beta_pdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_beta_pdf_f64

    function scifort_beta_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_beta_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: a !! first positive shape parameter
        real(c_double), value, intent(in) :: b !! second positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(beta_cdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_beta_cdf_f64

    function scifort_beta_ppf_f64(probability, a, b, loc, scale) result(x) &
            bind(c, name="scifort_beta_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: a !! first positive shape parameter
        real(c_double), value, intent(in) :: b !! second positive shape parameter
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(beta_ppf(real(probability, dp), real(a, dp), real(b, dp),  &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_beta_ppf_f64

    function scifort_f_pdf_f64(x, dfn, dfd, loc, scale) result(y) &
            bind(c, name="scifort_f_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: dfn !! positive numerator degrees of freedom
        real(c_double), value, intent(in) :: dfd !! positive denominator degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(f_pdf(real(x, dp), real(dfn, dp), real(dfd, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_f_pdf_f64

    function scifort_f_cdf_f64(x, dfn, dfd, loc, scale) result(y) &
            bind(c, name="scifort_f_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the function is evaluated
        real(c_double), value, intent(in) :: dfn !! positive numerator degrees of freedom
        real(c_double), value, intent(in) :: dfd !! positive denominator degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(f_cdf(real(x, dp), real(dfn, dp), real(dfd, dp), real(loc, dp),  &
            real(scale, dp)), c_double)
    end function scifort_f_cdf_f64

    function scifort_f_ppf_f64(probability, dfn, dfd, loc, scale) result(x) &
            bind(c, name="scifort_f_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: dfn !! positive numerator degrees of freedom
        real(c_double), value, intent(in) :: dfd !! positive denominator degrees of freedom
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(f_ppf(real(probability, dp), real(dfn, dp), real(dfd, dp),  &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_f_ppf_f64

    function scifort_gumbel_r_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_gumbel_r_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(gumbel_r_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gumbel_r_pdf_f64

    function scifort_gumbel_r_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_gumbel_r_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(gumbel_r_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gumbel_r_cdf_f64

    function scifort_gumbel_r_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_gumbel_r_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(gumbel_r_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gumbel_r_ppf_f64

    function scifort_gumbel_l_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_gumbel_l_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(gumbel_l_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gumbel_l_pdf_f64

    function scifort_gumbel_l_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_gumbel_l_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(gumbel_l_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gumbel_l_cdf_f64

    function scifort_gumbel_l_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_gumbel_l_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(gumbel_l_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gumbel_l_ppf_f64

    function scifort_powerlaw_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_powerlaw_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: a !! finite positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! density

        y = real(powerlaw_pdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powerlaw_pdf_f64

    function scifort_powerlaw_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_powerlaw_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: a !! finite positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! lower-tail probability

        y = real(powerlaw_cdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powerlaw_cdf_f64

    function scifort_powerlaw_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_powerlaw_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: a !! finite positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: x !! quantile

        x = real(powerlaw_ppf(real(p, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powerlaw_ppf_f64

    function scifort_triang_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_triang_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! standardized mode in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! density

        y = real(triang_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_triang_pdf_f64

    function scifort_triang_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_triang_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: c !! standardized mode in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! lower-tail probability

        y = real(triang_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_triang_cdf_f64

    function scifort_triang_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_triang_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! standardized mode in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: x !! quantile

        x = real(triang_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_triang_ppf_f64

    function scifort_genpareto_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genpareto_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! finite shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(genpareto_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genpareto_pdf_f64

    function scifort_genpareto_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genpareto_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(c_double), value, intent(in) :: c !! finite shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(genpareto_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genpareto_cdf_f64

    function scifort_genpareto_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_genpareto_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! finite shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(genpareto_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genpareto_ppf_f64

    function scifort_arcsine_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_arcsine_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! density

        y = real(arcsine_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_arcsine_pdf_f64

    function scifort_arcsine_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_arcsine_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! lower-tail probability

        y = real(arcsine_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_arcsine_cdf_f64

    function scifort_arcsine_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_arcsine_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: x !! quantile

        x = real(arcsine_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_arcsine_ppf_f64

    function scifort_halfnorm_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_halfnorm_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(halfnorm_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halfnorm_pdf_f64

    function scifort_halfnorm_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_halfnorm_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(halfnorm_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halfnorm_cdf_f64

    function scifort_halfnorm_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_halfnorm_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(halfnorm_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halfnorm_ppf_f64

    function scifort_halfcauchy_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_halfcauchy_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(halfcauchy_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halfcauchy_pdf_f64

    function scifort_halfcauchy_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_halfcauchy_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(halfcauchy_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halfcauchy_cdf_f64

    function scifort_halfcauchy_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_halfcauchy_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(halfcauchy_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halfcauchy_ppf_f64

    function scifort_lomax_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_lomax_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(lomax_pdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_lomax_pdf_f64

    function scifort_lomax_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_lomax_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(lomax_cdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_lomax_cdf_f64

    function scifort_lomax_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_lomax_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(lomax_ppf(real(p, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_lomax_ppf_f64

    function scifort_chi_pdf_f64(x, df, loc, scale) result(y) &
            bind(c, name="scifort_chi_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(chi_pdf(real(x, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_chi_pdf_f64

    function scifort_chi_cdf_f64(x, df, loc, scale) result(y) &
            bind(c, name="scifort_chi_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(chi_cdf(real(x, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_chi_cdf_f64

    function scifort_chi_ppf_f64(p, df, loc, scale) result(x) &
            bind(c, name="scifort_chi_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: df !! positive degrees of freedom
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(chi_ppf(real(p, dp), real(df, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_chi_ppf_f64

    function scifort_maxwell_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_maxwell_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(maxwell_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_maxwell_pdf_f64

    function scifort_maxwell_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_maxwell_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(maxwell_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_maxwell_cdf_f64

    function scifort_maxwell_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_maxwell_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(maxwell_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_maxwell_ppf_f64

    function scifort_cosine_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_cosine_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(cosine_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_cosine_pdf_f64

    function scifort_cosine_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_cosine_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(cosine_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_cosine_cdf_f64

    function scifort_cosine_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_cosine_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(cosine_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_cosine_ppf_f64

    function scifort_semicircular_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_semicircular_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density

        y = real(semicircular_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_semicircular_pdf_f64

    function scifort_semicircular_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_semicircular_cdf_f64")
        real(c_double), value, intent(in) :: x !! point x of the lower-tail probability
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(semicircular_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_semicircular_cdf_f64

    function scifort_semicircular_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_semicircular_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(semicircular_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_semicircular_ppf_f64

    function scifort_bernoulli_pmf_f64(k, p, loc) result(y) &
            bind(c, name="scifort_bernoulli_pmf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! probability mass

        y = real(bernoulli_pmf(real(k, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_bernoulli_pmf_f64

    function scifort_bernoulli_cdf_f64(k, p, loc) result(y) &
            bind(c, name="scifort_bernoulli_cdf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! lower-tail probability

        y = real(bernoulli_cdf(real(k, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_bernoulli_cdf_f64

    function scifort_bernoulli_ppf_f64(probability, p, loc) result(x) &
            bind(c, name="scifort_bernoulli_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: x !! quantile

        x = real(bernoulli_ppf(real(probability, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_bernoulli_ppf_f64

    function scifort_poisson_pmf_f64(k, mu, loc) result(y) &
            bind(c, name="scifort_poisson_pmf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: mu !! nonnegative Poisson mean
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! probability mass

        y = real(poisson_pmf(real(k, dp), real(mu, dp), real(loc, dp)), c_double)
    end function scifort_poisson_pmf_f64

    function scifort_poisson_cdf_f64(k, mu, loc) result(y) &
            bind(c, name="scifort_poisson_cdf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: mu !! nonnegative Poisson mean
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! lower-tail probability

        y = real(poisson_cdf(real(k, dp), real(mu, dp), real(loc, dp)), c_double)
    end function scifort_poisson_cdf_f64

    function scifort_poisson_ppf_f64(probability, mu, loc) result(x) &
            bind(c, name="scifort_poisson_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: mu !! nonnegative Poisson mean
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: x !! quantile

        x = real(poisson_ppf(real(probability, dp), real(mu, dp), real(loc, dp)), c_double)
    end function scifort_poisson_ppf_f64

    function scifort_geometric_pmf_f64(k, p, loc) result(y) &
            bind(c, name="scifort_geometric_pmf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! probability mass

        y = real(geometric_pmf(real(k, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_geometric_pmf_f64

    function scifort_geometric_cdf_f64(k, p, loc) result(y) &
            bind(c, name="scifort_geometric_cdf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! lower-tail probability

        y = real(geometric_cdf(real(k, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_geometric_cdf_f64

    function scifort_geometric_ppf_f64(probability, p, loc) result(x) &
            bind(c, name="scifort_geometric_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: x !! quantile

        x = real(geometric_ppf(real(probability, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_geometric_ppf_f64

    function scifort_binomial_pmf_f64(k, n, p, loc) result(y) &
            bind(c, name="scifort_binomial_pmf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: n !! trial count or negative-binomial success shape
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! probability mass

        y = real(binomial_pmf(real(k, dp), real(n, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_binomial_pmf_f64

    function scifort_binomial_cdf_f64(k, n, p, loc) result(y) &
            bind(c, name="scifort_binomial_cdf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: n !! trial count or negative-binomial success shape
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! lower-tail probability

        y = real(binomial_cdf(real(k, dp), real(n, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_binomial_cdf_f64

    function scifort_binomial_ppf_f64(probability, n, p, loc) result(x) &
            bind(c, name="scifort_binomial_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: n !! trial count or negative-binomial success shape
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: x !! quantile

        x = real(binomial_ppf(real(probability, dp), real(n, dp), real(p, dp),  &
            real(loc, dp)), c_double)
    end function scifort_binomial_ppf_f64

    function scifort_negative_binomial_pmf_f64(k, n, p, loc) result(y) &
            bind(c, name="scifort_negative_binomial_pmf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: n !! trial count or negative-binomial success shape
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! probability mass

        y = real(negative_binomial_pmf(real(k, dp), real(n, dp), real(p, dp),  &
            real(loc, dp)), c_double)
    end function scifort_negative_binomial_pmf_f64

    function scifort_negative_binomial_cdf_f64(k, n, p, loc) result(y) &
            bind(c, name="scifort_negative_binomial_cdf_f64")
        real(c_double), value, intent(in) :: k !! count or count threshold in the shifted support
        real(c_double), value, intent(in) :: n !! trial count or negative-binomial success shape
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: y !! lower-tail probability

        y = real(negative_binomial_cdf(real(k, dp), real(n, dp), real(p, dp),  &
            real(loc, dp)), c_double)
    end function scifort_negative_binomial_cdf_f64

    function scifort_negative_binomial_ppf_f64(probability, n, p, loc) result(x) &
            bind(c, name="scifort_negative_binomial_ppf_f64")
        real(c_double), value, intent(in) :: probability !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: n !! trial count or negative-binomial success shape
        real(c_double), value, intent(in) :: p !! success probability
        real(c_double), value, intent(in) :: loc !! location or support shift
        real(c_double) :: x !! quantile

        x = real(negative_binomial_ppf(real(probability, dp), real(n, dp),  &
            real(p, dp), real(loc, dp)), c_double)
    end function scifort_negative_binomial_ppf_f64
    function scifort_anglit_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_anglit_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(anglit_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_anglit_pdf_f64

    function scifort_anglit_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_anglit_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(anglit_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_anglit_cdf_f64

    function scifort_anglit_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_anglit_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(anglit_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_anglit_ppf_f64

    function scifort_moyal_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_moyal_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(moyal_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_moyal_pdf_f64

    function scifort_moyal_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_moyal_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(moyal_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_moyal_cdf_f64

    function scifort_moyal_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_moyal_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(moyal_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_moyal_ppf_f64

    function scifort_hypsecant_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_hypsecant_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(hypsecant_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_hypsecant_pdf_f64

    function scifort_hypsecant_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_hypsecant_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(hypsecant_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_hypsecant_cdf_f64

    function scifort_hypsecant_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_hypsecant_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(hypsecant_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_hypsecant_ppf_f64

    function scifort_halflogistic_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_halflogistic_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(halflogistic_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halflogistic_pdf_f64

    function scifort_halflogistic_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_halflogistic_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(halflogistic_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halflogistic_cdf_f64

    function scifort_halflogistic_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_halflogistic_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(halflogistic_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_halflogistic_ppf_f64

    function scifort_invgamma_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_invgamma_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: a !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(invgamma_pdf(real(x, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_invgamma_pdf_f64

    function scifort_invgamma_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_invgamma_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: a !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(invgamma_cdf(real(x, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_invgamma_cdf_f64

    function scifort_invgamma_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_invgamma_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: a !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(invgamma_ppf(real(p, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_invgamma_ppf_f64

    function scifort_invgauss_pdf_f64(x, mu_shape, loc, scale) result(y) &
            bind(c, name="scifort_invgauss_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: mu_shape !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(invgauss_pdf(real(x, dp), real(mu_shape, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_invgauss_pdf_f64

    function scifort_invgauss_cdf_f64(x, mu_shape, loc, scale) result(y) &
            bind(c, name="scifort_invgauss_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: mu_shape !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(invgauss_cdf(real(x, dp), real(mu_shape, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_invgauss_cdf_f64

    function scifort_invgauss_ppf_f64(p, mu_shape, loc, scale) result(x) &
            bind(c, name="scifort_invgauss_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: mu_shape !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(invgauss_ppf(real(p, dp), real(mu_shape, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_invgauss_ppf_f64

    function scifort_levy_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_levy_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(levy_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_levy_pdf_f64

    function scifort_levy_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_levy_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(levy_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_levy_cdf_f64

    function scifort_levy_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_levy_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(levy_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_levy_ppf_f64

    function scifort_loglaplace_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_loglaplace_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(loglaplace_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loglaplace_pdf_f64

    function scifort_loglaplace_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_loglaplace_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(loglaplace_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loglaplace_cdf_f64

    function scifort_loglaplace_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_loglaplace_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(loglaplace_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loglaplace_ppf_f64


    function scifort_bradford_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_bradford_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! probability density

        y = real(bradford_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_bradford_pdf_f64

    function scifort_bradford_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_bradford_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! lower-tail probability

        y = real(bradford_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_bradford_cdf_f64

    function scifort_bradford_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_bradford_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: x !! quantile

        x = real(bradford_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_bradford_ppf_f64

    function scifort_truncexpon_pdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_truncexpon_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: b !! positive standardized upper endpoint
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(truncexpon_pdf(real(x, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncexpon_pdf_f64

    function scifort_truncexpon_cdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_truncexpon_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: b !! positive standardized upper endpoint
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(truncexpon_cdf(real(x, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncexpon_cdf_f64

    function scifort_truncexpon_ppf_f64(p, b, loc, scale) result(x) &
            bind(c, name="scifort_truncexpon_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: b !! positive standardized upper endpoint
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(truncexpon_ppf(real(p, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncexpon_ppf_f64

    function scifort_fisk_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_fisk_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(fisk_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_fisk_pdf_f64

    function scifort_fisk_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_fisk_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(fisk_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_fisk_cdf_f64

    function scifort_fisk_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_fisk_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(fisk_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_fisk_ppf_f64

    function scifort_dweibull_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_dweibull_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density

        y = real(dweibull_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_dweibull_pdf_f64

    function scifort_dweibull_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_dweibull_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability

        y = real(dweibull_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_dweibull_cdf_f64

    function scifort_dweibull_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_dweibull_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile

        x = real(dweibull_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_dweibull_ppf_f64

    function scifort_alpha_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_alpha_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: a !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(alpha_pdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_alpha_pdf_f64

    function scifort_alpha_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_alpha_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: a !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(alpha_cdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_alpha_cdf_f64

    function scifort_alpha_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_alpha_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: a !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(alpha_ppf(real(p, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_alpha_ppf_f64

    function scifort_fatiguelife_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_fatiguelife_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(fatiguelife_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_fatiguelife_pdf_f64

    function scifort_fatiguelife_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_fatiguelife_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(fatiguelife_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_fatiguelife_cdf_f64

    function scifort_fatiguelife_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_fatiguelife_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(fatiguelife_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_fatiguelife_ppf_f64

    function scifort_genlogistic_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genlogistic_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(genlogistic_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genlogistic_pdf_f64

    function scifort_genlogistic_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genlogistic_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(genlogistic_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genlogistic_cdf_f64

    function scifort_genlogistic_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_genlogistic_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(genlogistic_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genlogistic_ppf_f64

    function scifort_gennorm_pdf_f64(x, beta, loc, scale) result(y) &
            bind(c, name="scifort_gennorm_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: beta !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(gennorm_pdf(real(x, dp), real(beta, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gennorm_pdf_f64

    function scifort_gennorm_cdf_f64(x, beta, loc, scale) result(y) &
            bind(c, name="scifort_gennorm_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: beta !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(gennorm_cdf(real(x, dp), real(beta, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gennorm_cdf_f64

    function scifort_gennorm_ppf_f64(p, beta, loc, scale) result(x) &
            bind(c, name="scifort_gennorm_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: beta !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(gennorm_ppf(real(p, dp), real(beta, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gennorm_ppf_f64

    pure logical function valid_c_loc_scale(loc, scale) result(valid)
        real(c_double), intent(in) :: loc !! mean to check
        real(c_double), intent(in) :: scale !! standard deviation to check

        valid = ieee_is_finite(loc) .and. ieee_is_finite(scale) .and. &
            scale > 0.0_c_double
    end function valid_c_loc_scale

    function scifort_nakagami_pdf_f64(x, nu, loc, scale) result(y) &
            bind(c, name="scifort_nakagami_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: nu !! positive Nakagami shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(nakagami_pdf(real(x, dp), real(nu, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_nakagami_pdf_f64

    function scifort_nakagami_cdf_f64(x, nu, loc, scale) result(y) &
            bind(c, name="scifort_nakagami_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: nu !! positive Nakagami shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(nakagami_cdf(real(x, dp), real(nu, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_nakagami_cdf_f64

    function scifort_nakagami_ppf_f64(p, nu, loc, scale) result(x) &
            bind(c, name="scifort_nakagami_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: nu !! positive Nakagami shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(nakagami_ppf(real(p, dp), real(nu, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_nakagami_ppf_f64

    function scifort_powernorm_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_powernorm_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive power-normal shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(powernorm_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powernorm_pdf_f64

    function scifort_powernorm_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_powernorm_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive power-normal shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(powernorm_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powernorm_cdf_f64

    function scifort_powernorm_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_powernorm_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive power-normal shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(powernorm_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powernorm_ppf_f64

    function scifort_loggamma_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_loggamma_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive log-gamma shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(loggamma_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loggamma_pdf_f64

    function scifort_loggamma_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_loggamma_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive log-gamma shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(loggamma_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loggamma_cdf_f64

    function scifort_loggamma_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_loggamma_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive log-gamma shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(loggamma_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loggamma_ppf_f64

    function scifort_wald_pdf_f64(x, loc, scale) result(y) bind(c, name="scifort_wald_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(wald_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_wald_pdf_f64

    function scifort_wald_cdf_f64(x, loc, scale) result(y) bind(c, name="scifort_wald_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(wald_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_wald_cdf_f64

    function scifort_wald_ppf_f64(p, loc, scale) result(x) bind(c, name="scifort_wald_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(wald_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_wald_ppf_f64


    function scifort_gompertz_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_gompertz_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive Gompertz shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(gompertz_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gompertz_pdf_f64

    function scifort_gompertz_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_gompertz_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive Gompertz shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(gompertz_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gompertz_cdf_f64

    function scifort_gompertz_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_gompertz_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive Gompertz shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(gompertz_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gompertz_ppf_f64

    function scifort_invweibull_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_invweibull_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive inverse-Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(invweibull_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_invweibull_pdf_f64

    function scifort_invweibull_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_invweibull_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive inverse-Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(invweibull_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_invweibull_cdf_f64

    function scifort_invweibull_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_invweibull_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive inverse-Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(invweibull_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_invweibull_ppf_f64

    function scifort_betaprime_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_betaprime_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: a !! positive first beta-prime shape parameter
        real(c_double), value, intent(in) :: b !! positive second beta-prime shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(betaprime_pdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_betaprime_pdf_f64

    function scifort_betaprime_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_betaprime_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: a !! positive first beta-prime shape parameter
        real(c_double), value, intent(in) :: b !! positive second beta-prime shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(betaprime_cdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_betaprime_cdf_f64

    function scifort_betaprime_ppf_f64(p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_betaprime_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: a !! positive first beta-prime shape parameter
        real(c_double), value, intent(in) :: b !! positive second beta-prime shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(betaprime_ppf(real(p, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_betaprime_ppf_f64

    function scifort_burr12_pdf_f64(x, c, d, loc, scale) result(y) &
            bind(c, name="scifort_burr12_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive first Burr-XII shape parameter
        real(c_double), value, intent(in) :: d !! positive second Burr-XII shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(burr12_pdf(real(x, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_burr12_pdf_f64

    function scifort_burr12_cdf_f64(x, c, d, loc, scale) result(y) &
            bind(c, name="scifort_burr12_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive first Burr-XII shape parameter
        real(c_double), value, intent(in) :: d !! positive second Burr-XII shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(burr12_cdf(real(x, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_burr12_cdf_f64

    function scifort_burr12_ppf_f64(p, c, d, loc, scale) result(x) &
            bind(c, name="scifort_burr12_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0, 1]
        real(c_double), value, intent(in) :: c !! positive first Burr-XII shape parameter
        real(c_double), value, intent(in) :: d !! positive second Burr-XII shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(burr12_ppf(real(p, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_burr12_ppf_f64

    function scifort_genhalflogistic_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genhalflogistic_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(genhalflogistic_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genhalflogistic_pdf_f64

    function scifort_genhalflogistic_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genhalflogistic_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(genhalflogistic_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genhalflogistic_cdf_f64

    function scifort_genhalflogistic_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_genhalflogistic_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: c !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(genhalflogistic_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genhalflogistic_ppf_f64

    function scifort_exponpow_pdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_exponpow_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: b !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(exponpow_pdf(real(x, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponpow_pdf_f64

    function scifort_exponpow_cdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_exponpow_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: b !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(exponpow_cdf(real(x, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponpow_cdf_f64

    function scifort_exponpow_ppf_f64(p, b, loc, scale) result(x) &
            bind(c, name="scifort_exponpow_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: b !! positive shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(exponpow_ppf(real(p, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponpow_ppf_f64

    function scifort_exponweib_pdf_f64(x, a, c, loc, scale) result(y) &
            bind(c, name="scifort_exponweib_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: a !! positive exponentiation shape parameter
        real(c_double), value, intent(in) :: c !! positive Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(exponweib_pdf(real(x, dp), real(a, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_exponweib_pdf_f64

    function scifort_exponweib_cdf_f64(x, a, c, loc, scale) result(y) &
            bind(c, name="scifort_exponweib_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: a !! positive exponentiation shape parameter
        real(c_double), value, intent(in) :: c !! positive Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(exponweib_cdf(real(x, dp), real(a, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_exponweib_cdf_f64

    function scifort_exponweib_ppf_f64(p, a, c, loc, scale) result(x) &
            bind(c, name="scifort_exponweib_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: a !! positive exponentiation shape parameter
        real(c_double), value, intent(in) :: c !! positive Weibull shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(exponweib_ppf(real(p, dp), real(a, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_exponweib_ppf_f64

    function scifort_powerlognorm_pdf_f64(x, c, s, loc, scale) result(y) &
            bind(c, name="scifort_powerlognorm_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive power shape parameter
        real(c_double), value, intent(in) :: s !! positive lognormal shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(powerlognorm_pdf(real(x, dp), real(c, dp), real(s, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powerlognorm_pdf_f64

    function scifort_powerlognorm_cdf_f64(x, c, s, loc, scale) result(y) &
            bind(c, name="scifort_powerlognorm_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive power shape parameter
        real(c_double), value, intent(in) :: s !! positive lognormal shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(powerlognorm_cdf(real(x, dp), real(c, dp), real(s, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powerlognorm_cdf_f64

    function scifort_powerlognorm_ppf_f64(p, c, s, loc, scale) result(x) &
            bind(c, name="scifort_powerlognorm_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: c !! positive power shape parameter
        real(c_double), value, intent(in) :: s !! positive lognormal shape parameter
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(powerlognorm_ppf(real(p, dp), real(c, dp), real(s, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_powerlognorm_ppf_f64

    function scifort_levy_l_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_levy_l_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! upper support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(levy_l_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_levy_l_pdf_f64

    function scifort_levy_l_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_levy_l_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! upper support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(levy_l_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_levy_l_cdf_f64

    function scifort_levy_l_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_levy_l_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: loc !! upper support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(levy_l_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_levy_l_ppf_f64

    function scifort_weibull_max_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_weibull_max_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive Weibull-maximum shape parameter
        real(c_double), value, intent(in) :: loc !! upper support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(weibull_max_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_weibull_max_pdf_f64

    function scifort_weibull_max_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_weibull_max_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive Weibull-maximum shape parameter
        real(c_double), value, intent(in) :: loc !! upper support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(weibull_max_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_weibull_max_cdf_f64

    function scifort_weibull_max_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_weibull_max_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: c !! positive Weibull-maximum shape parameter
        real(c_double), value, intent(in) :: loc !! upper support endpoint
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(weibull_max_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_weibull_max_ppf_f64

    function scifort_rdist_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_rdist_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: c !! positive R-distribution shape parameter
        real(c_double), value, intent(in) :: loc !! center of support
        real(c_double), value, intent(in) :: scale !! positive half-width
        real(c_double) :: y !! probability density
        y = real(rdist_pdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rdist_pdf_f64

    function scifort_rdist_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_rdist_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: c !! positive R-distribution shape parameter
        real(c_double), value, intent(in) :: loc !! center of support
        real(c_double), value, intent(in) :: scale !! positive half-width
        real(c_double) :: y !! lower-tail probability
        y = real(rdist_cdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rdist_cdf_f64

    function scifort_rdist_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_rdist_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: c !! positive R-distribution shape parameter
        real(c_double), value, intent(in) :: loc !! center of support
        real(c_double), value, intent(in) :: scale !! positive half-width
        real(c_double) :: x !! quantile
        x = real(rdist_ppf(real(p, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rdist_ppf_f64

    function scifort_skewcauchy_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_skewcauchy_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! probability density
        y = real(skewcauchy_pdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_skewcauchy_pdf_f64

    function scifort_skewcauchy_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_skewcauchy_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(skewcauchy_cdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_skewcauchy_cdf_f64

    function scifort_skewcauchy_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_skewcauchy_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(skewcauchy_ppf(real(p, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_skewcauchy_ppf_f64


    function scifort_dgamma_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_dgamma_pdf_f64")
        real(c_double), value :: x, a, loc, scale
        real(c_double) :: y
        y = real(dgamma_pdf(real(x, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_dgamma_pdf_f64

    function scifort_dgamma_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_dgamma_cdf_f64")
        real(c_double), value :: x, a, loc, scale
        real(c_double) :: y
        y = real(dgamma_cdf(real(x, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_dgamma_cdf_f64

    function scifort_dgamma_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_dgamma_ppf_f64")
        real(c_double), value :: p, a, loc, scale
        real(c_double) :: x
        x = real(dgamma_ppf(real(p, dp), real(a, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_dgamma_ppf_f64

    function scifort_laplace_asymmetric_pdf_f64(x, kappa, loc, scale) result(y) &
            bind(c, name="scifort_laplace_asymmetric_pdf_f64")
        real(c_double), value :: x, kappa, loc, scale
        real(c_double) :: y
        y = real(laplace_asymmetric_pdf(real(x, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_laplace_asymmetric_pdf_f64

    function scifort_laplace_asymmetric_cdf_f64(x, kappa, loc, scale) result(y) &
            bind(c, name="scifort_laplace_asymmetric_cdf_f64")
        real(c_double), value :: x, kappa, loc, scale
        real(c_double) :: y
        y = real(laplace_asymmetric_cdf(real(x, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_laplace_asymmetric_cdf_f64

    function scifort_laplace_asymmetric_ppf_f64(p, kappa, loc, scale) result(x) &
            bind(c, name="scifort_laplace_asymmetric_ppf_f64")
        real(c_double), value :: p, kappa, loc, scale
        real(c_double) :: x
        x = real(laplace_asymmetric_ppf(real(p, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_laplace_asymmetric_ppf_f64

    function scifort_truncnorm_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_truncnorm_pdf_f64")
        real(c_double), value :: x, a, b, loc, scale
        real(c_double) :: y
        y = real(truncnorm_pdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncnorm_pdf_f64

    function scifort_truncnorm_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_truncnorm_cdf_f64")
        real(c_double), value :: x, a, b, loc, scale
        real(c_double) :: y
        y = real(truncnorm_cdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncnorm_cdf_f64

    function scifort_truncnorm_ppf_f64(p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_truncnorm_ppf_f64")
        real(c_double), value :: p, a, b, loc, scale
        real(c_double) :: x
        x = real(truncnorm_ppf(real(p, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncnorm_ppf_f64

    function scifort_loguniform_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_loguniform_pdf_f64")
        real(c_double), value :: x, a, b, loc, scale
        real(c_double) :: y
        y = real(loguniform_pdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loguniform_pdf_f64

    function scifort_loguniform_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_loguniform_cdf_f64")
        real(c_double), value :: x, a, b, loc, scale
        real(c_double) :: y
        y = real(loguniform_cdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loguniform_cdf_f64

    function scifort_loguniform_ppf_f64(p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_loguniform_ppf_f64")
        real(c_double), value :: p, a, b, loc, scale
        real(c_double) :: x
        x = real(loguniform_ppf(real(p, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_loguniform_ppf_f64

    function scifort_foldnorm_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_foldnorm_pdf_f64")
        real(c_double), value :: x, c, loc, scale
        real(c_double) :: y
        y = real(foldnorm_pdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_foldnorm_pdf_f64

    function scifort_foldnorm_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_foldnorm_cdf_f64")
        real(c_double), value :: x, c, loc, scale
        real(c_double) :: y
        y = real(foldnorm_cdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_foldnorm_cdf_f64

    function scifort_foldnorm_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_foldnorm_ppf_f64")
        real(c_double), value :: p, c, loc, scale
        real(c_double) :: x
        x = real(foldnorm_ppf(real(p, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_foldnorm_ppf_f64

    function scifort_foldcauchy_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_foldcauchy_pdf_f64")
        real(c_double), value :: x, c, loc, scale
        real(c_double) :: y
        y = real(foldcauchy_pdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_foldcauchy_pdf_f64

    function scifort_foldcauchy_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_foldcauchy_cdf_f64")
        real(c_double), value :: x, c, loc, scale
        real(c_double) :: y
        y = real(foldcauchy_cdf(real(x, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_foldcauchy_cdf_f64

    function scifort_foldcauchy_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_foldcauchy_ppf_f64")
        real(c_double), value :: p, c, loc, scale
        real(c_double) :: x
        x = real(foldcauchy_ppf(real(p, dp), real(c, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_foldcauchy_ppf_f64

    function scifort_recipinvgauss_pdf_f64(x, mu_shape, loc, scale) result(y) &
            bind(c, name="scifort_recipinvgauss_pdf_f64")
        real(c_double), value :: x, mu_shape, loc, scale
        real(c_double) :: y
        y = real(recipinvgauss_pdf(real(x, dp), real(mu_shape, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_recipinvgauss_pdf_f64

    function scifort_recipinvgauss_cdf_f64(x, mu_shape, loc, scale) result(y) &
            bind(c, name="scifort_recipinvgauss_cdf_f64")
        real(c_double), value :: x, mu_shape, loc, scale
        real(c_double) :: y
        y = real(recipinvgauss_cdf(real(x, dp), real(mu_shape, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_recipinvgauss_cdf_f64

    function scifort_recipinvgauss_ppf_f64(p, mu_shape, loc, scale) result(x) &
            bind(c, name="scifort_recipinvgauss_ppf_f64")
        real(c_double), value :: p, mu_shape, loc, scale
        real(c_double) :: x
        x = real(recipinvgauss_ppf(real(p, dp), real(mu_shape, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_recipinvgauss_ppf_f64

    function scifort_truncpareto_pdf_f64(x, b, c, loc, scale) result(y) &
            bind(c, name="scifort_truncpareto_pdf_f64")
        real(c_double), value :: x, b, c, loc, scale
        real(c_double) :: y
        y = real(truncpareto_pdf(real(x, dp), real(b, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncpareto_pdf_f64

    function scifort_truncpareto_cdf_f64(x, b, c, loc, scale) result(y) &
            bind(c, name="scifort_truncpareto_cdf_f64")
        real(c_double), value :: x, b, c, loc, scale
        real(c_double) :: y
        y = real(truncpareto_cdf(real(x, dp), real(b, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncpareto_cdf_f64

    function scifort_truncpareto_ppf_f64(p, b, c, loc, scale) result(x) &
            bind(c, name="scifort_truncpareto_ppf_f64")
        real(c_double), value :: p, b, c, loc, scale
        real(c_double) :: x
        x = real(truncpareto_ppf(real(p, dp), real(b, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_truncpareto_ppf_f64

    function scifort_exponnorm_pdf_f64(x, k, loc, scale) result(y) &
            bind(c, name="scifort_exponnorm_pdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: k !! positive exponential-to-normal scale ratio
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density
        y = real(exponnorm_pdf(real(x, dp), real(k, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponnorm_pdf_f64

    function scifort_exponnorm_cdf_f64(x, k, loc, scale) result(y) &
            bind(c, name="scifort_exponnorm_cdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: k !! positive exponential-to-normal scale ratio
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! lower-tail probability
        y = real(exponnorm_cdf(real(x, dp), real(k, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponnorm_cdf_f64

    function scifort_exponnorm_ppf_f64(p, k, loc, scale) result(x) &
            bind(c, name="scifort_exponnorm_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: k !! positive exponential-to-normal scale ratio
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(exponnorm_ppf(real(p, dp), real(k, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_exponnorm_ppf_f64

    function scifort_johnsonsb_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_johnsonsb_pdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: a !! finite first Johnson shape parameter
        real(c_double), value, intent(in) :: b !! positive second Johnson shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density or lower-tail probability
        y = real(johnsonsb_pdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_johnsonsb_pdf_f64

    function scifort_johnsonsb_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_johnsonsb_cdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: a !! finite first Johnson shape parameter
        real(c_double), value, intent(in) :: b !! positive second Johnson shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density or lower-tail probability
        y = real(johnsonsb_cdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_johnsonsb_cdf_f64

    function scifort_johnsonsb_ppf_f64(p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_johnsonsb_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: a !! finite first Johnson shape parameter
        real(c_double), value, intent(in) :: b !! positive second Johnson shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(johnsonsb_ppf(real(p, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_johnsonsb_ppf_f64

    function scifort_johnsonsu_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_johnsonsu_pdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: a !! finite first Johnson shape parameter
        real(c_double), value, intent(in) :: b !! positive second Johnson shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density or lower-tail probability
        y = real(johnsonsu_pdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_johnsonsu_pdf_f64

    function scifort_johnsonsu_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_johnsonsu_cdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: a !! finite first Johnson shape parameter
        real(c_double), value, intent(in) :: b !! positive second Johnson shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y !! density or lower-tail probability
        y = real(johnsonsu_cdf(real(x, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_johnsonsu_cdf_f64

    function scifort_johnsonsu_ppf_f64(p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_johnsonsu_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: a !! finite first Johnson shape parameter
        real(c_double), value, intent(in) :: b !! positive second Johnson shape parameter
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x !! quantile
        x = real(johnsonsu_ppf(real(p, dp), real(a, dp), real(b, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_johnsonsu_ppf_f64

    function scifort_trapezoid_pdf_f64(x, c, d, loc, scale) result(y) &
            bind(c, name="scifort_trapezoid_pdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: c !! standardized left plateau edge in [0,1]
        real(c_double), value, intent(in) :: d !! standardized right plateau edge in [c,1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! density or lower-tail probability
        y = real(trapezoid_pdf(real(x, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_trapezoid_pdf_f64

    function scifort_trapezoid_cdf_f64(x, c, d, loc, scale) result(y) &
            bind(c, name="scifort_trapezoid_cdf_f64")
        real(c_double), value, intent(in) :: x !! evaluation point
        real(c_double), value, intent(in) :: c !! standardized left plateau edge in [0,1]
        real(c_double), value, intent(in) :: d !! standardized right plateau edge in [c,1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: y !! density or lower-tail probability
        y = real(trapezoid_cdf(real(x, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_trapezoid_cdf_f64

    function scifort_trapezoid_ppf_f64(p, c, d, loc, scale) result(x) &
            bind(c, name="scifort_trapezoid_ppf_f64")
        real(c_double), value, intent(in) :: p !! lower-tail probability in [0,1]
        real(c_double), value, intent(in) :: c !! standardized left plateau edge in [0,1]
        real(c_double), value, intent(in) :: d !! standardized right plateau edge in [c,1]
        real(c_double), value, intent(in) :: loc !! lower support endpoint
        real(c_double), value, intent(in) :: scale !! positive support width
        real(c_double) :: x !! quantile
        x = real(trapezoid_ppf(real(p, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_trapezoid_ppf_f64

    function scifort_burr_pdf_f64(x, c, d, loc, scale) result(y) &
            bind(c, name="scifort_burr_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! positive first Burr shape
        real(c_double), value :: d !! positive second Burr shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(burr_pdf(real(x, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_burr_pdf_f64

    function scifort_burr_cdf_f64(x, c, d, loc, scale) result(y) &
            bind(c, name="scifort_burr_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! positive first Burr shape
        real(c_double), value :: d !! positive second Burr shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(burr_cdf(real(x, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_burr_cdf_f64

    function scifort_burr_ppf_f64(p, c, d, loc, scale) result(x) &
            bind(c, name="scifort_burr_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: c !! positive first Burr shape
        real(c_double), value :: d !! positive second Burr shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(burr_ppf(real(p, dp), real(c, dp), real(d, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_burr_ppf_f64

    function scifort_mielke_pdf_f64(x, k, s, loc, scale) result(y) &
            bind(c, name="scifort_mielke_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: k !! positive first beta-kappa shape
        real(c_double), value :: s !! positive second beta-kappa shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(mielke_pdf(real(x, dp), real(k, dp), real(s, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_mielke_pdf_f64

    function scifort_mielke_cdf_f64(x, k, s, loc, scale) result(y) &
            bind(c, name="scifort_mielke_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: k !! positive first beta-kappa shape
        real(c_double), value :: s !! positive second beta-kappa shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(mielke_cdf(real(x, dp), real(k, dp), real(s, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_mielke_cdf_f64

    function scifort_mielke_ppf_f64(p, k, s, loc, scale) result(x) &
            bind(c, name="scifort_mielke_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: k !! positive first beta-kappa shape
        real(c_double), value :: s !! positive second beta-kappa shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(mielke_ppf(real(p, dp), real(k, dp), real(s, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_mielke_ppf_f64

    function scifort_gibrat_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_gibrat_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(gibrat_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gibrat_pdf_f64

    function scifort_gibrat_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_gibrat_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(gibrat_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gibrat_cdf_f64

    function scifort_gibrat_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_gibrat_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(gibrat_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_gibrat_ppf_f64

    function scifort_wrapcauchy_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_wrapcauchy_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! concentration parameter in (0,1)
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive angular scale
        real(c_double) :: y
        y = real(wrapcauchy_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_wrapcauchy_pdf_f64

    function scifort_wrapcauchy_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_wrapcauchy_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! concentration parameter in (0,1)
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive angular scale
        real(c_double) :: y
        y = real(wrapcauchy_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_wrapcauchy_cdf_f64

    function scifort_wrapcauchy_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_wrapcauchy_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: c !! concentration parameter in (0,1)
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive angular scale
        real(c_double) :: x
        x = real(wrapcauchy_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_wrapcauchy_ppf_f64

    function scifort_genextreme_pdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genextreme_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! finite generalized-extreme-value shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(genextreme_pdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genextreme_pdf_f64

    function scifort_genextreme_cdf_f64(x, c, loc, scale) result(y) &
            bind(c, name="scifort_genextreme_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! finite generalized-extreme-value shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(genextreme_cdf(real(x, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genextreme_cdf_f64

    function scifort_genextreme_ppf_f64(p, c, loc, scale) result(x) &
            bind(c, name="scifort_genextreme_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: c !! finite generalized-extreme-value shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(genextreme_ppf(real(p, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_genextreme_ppf_f64

    function scifort_kappa3_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_kappa3_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive kappa-3 shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(kappa3_pdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_kappa3_pdf_f64

    function scifort_kappa3_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_kappa3_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive kappa-3 shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(kappa3_cdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_kappa3_cdf_f64

    function scifort_kappa3_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_kappa3_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive kappa-3 shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(kappa3_ppf(real(p, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_kappa3_ppf_f64

    function scifort_kappa4_pdf_f64(x, h, k, loc, scale) result(y) &
            bind(c, name="scifort_kappa4_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: h !! finite first kappa-4 shape
        real(c_double), value :: k !! finite second kappa-4 shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(kappa4_pdf(real(x, dp), real(h, dp), real(k, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_kappa4_pdf_f64

    function scifort_kappa4_cdf_f64(x, h, k, loc, scale) result(y) &
            bind(c, name="scifort_kappa4_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: h !! finite first kappa-4 shape
        real(c_double), value :: k !! finite second kappa-4 shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(kappa4_cdf(real(x, dp), real(h, dp), real(k, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_kappa4_cdf_f64

    function scifort_kappa4_ppf_f64(p, h, k, loc, scale) result(x) &
            bind(c, name="scifort_kappa4_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: h !! finite first kappa-4 shape
        real(c_double), value :: k !! finite second kappa-4 shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(kappa4_ppf(real(p, dp), real(h, dp), real(k, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_kappa4_ppf_f64

    function scifort_truncweibull_min_pdf_f64(x, c, a, b, loc, scale) result(y) &
            bind(c, name="scifort_truncweibull_min_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! positive Weibull shape
        real(c_double), value :: a !! standardized lower truncation point
        real(c_double), value :: b !! standardized upper truncation point
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(truncweibull_min_pdf(real(x, dp), real(c, dp), real(a, dp), &
            real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_truncweibull_min_pdf_f64

    function scifort_truncweibull_min_cdf_f64(x, c, a, b, loc, scale) result(y) &
            bind(c, name="scifort_truncweibull_min_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: c !! positive Weibull shape
        real(c_double), value :: a !! standardized lower truncation point
        real(c_double), value :: b !! standardized upper truncation point
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(truncweibull_min_cdf(real(x, dp), real(c, dp), real(a, dp), &
            real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_truncweibull_min_cdf_f64

    function scifort_truncweibull_min_ppf_f64(p, c, a, b, loc, scale) result(x) &
            bind(c, name="scifort_truncweibull_min_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: c !! positive Weibull shape
        real(c_double), value :: a !! standardized lower truncation point
        real(c_double), value :: b !! standardized upper truncation point
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(truncweibull_min_ppf(real(p, dp), real(c, dp), real(a, dp), &
            real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_truncweibull_min_ppf_f64

    function scifort_gengamma_pdf_f64(x, a, c, loc, scale) result(y) &
            bind(c, name="scifort_gengamma_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive gamma-shape parameter
        real(c_double), value :: c !! finite nonzero power parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(gengamma_pdf(real(x, dp), real(a, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gengamma_pdf_f64

    function scifort_gengamma_cdf_f64(x, a, c, loc, scale) result(y) &
            bind(c, name="scifort_gengamma_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive gamma-shape parameter
        real(c_double), value :: c !! finite nonzero power parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(gengamma_cdf(real(x, dp), real(a, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gengamma_cdf_f64

    function scifort_gengamma_ppf_f64(p, a, c, loc, scale) result(x) &
            bind(c, name="scifort_gengamma_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive gamma-shape parameter
        real(c_double), value :: c !! finite nonzero power parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(gengamma_ppf(real(p, dp), real(a, dp), real(c, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_gengamma_ppf_f64

    function scifort_halfgennorm_pdf_f64(x, beta, loc, scale) result(y) &
            bind(c, name="scifort_halfgennorm_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: beta !! positive generalized-normal shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(halfgennorm_pdf(real(x, dp), real(beta, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_halfgennorm_pdf_f64

    function scifort_halfgennorm_cdf_f64(x, beta, loc, scale) result(y) &
            bind(c, name="scifort_halfgennorm_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: beta !! positive generalized-normal shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(halfgennorm_cdf(real(x, dp), real(beta, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_halfgennorm_cdf_f64

    function scifort_halfgennorm_ppf_f64(p, beta, loc, scale) result(x) &
            bind(c, name="scifort_halfgennorm_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: beta !! positive generalized-normal shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(halfgennorm_ppf(real(p, dp), real(beta, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_halfgennorm_ppf_f64

    function scifort_argus_pdf_f64(x, chi, loc, scale) result(y) &
            bind(c, name="scifort_argus_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: chi !! positive ARGUS shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive support width
        real(c_double) :: y
        y = real(argus_pdf(real(x, dp), real(chi, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_argus_pdf_f64

    function scifort_argus_cdf_f64(x, chi, loc, scale) result(y) &
            bind(c, name="scifort_argus_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: chi !! positive ARGUS shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive support width
        real(c_double) :: y
        y = real(argus_cdf(real(x, dp), real(chi, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_argus_cdf_f64

    function scifort_argus_ppf_f64(p, chi, loc, scale) result(x) &
            bind(c, name="scifort_argus_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: chi !! positive ARGUS shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive support width
        real(c_double) :: x
        x = real(argus_ppf(real(p, dp), real(chi, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_argus_ppf_f64

    function scifort_erlang_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_erlang_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive Erlang/gamma shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(erlang_pdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_erlang_pdf_f64

    function scifort_erlang_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_erlang_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive Erlang/gamma shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(erlang_cdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_erlang_cdf_f64

    function scifort_erlang_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_erlang_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive Erlang/gamma shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(erlang_ppf(real(p, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_erlang_ppf_f64


    function scifort_crystalball_pdf_f64(x, beta, m, loc, scale) result(y) &
            bind(c, name="scifort_crystalball_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: beta !! positive tail-transition magnitude
        real(c_double), value :: m !! power-law exponent, > 1
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(crystalball_pdf(real(x, dp), real(beta, dp), real(m, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_crystalball_pdf_f64

    function scifort_crystalball_cdf_f64(x, beta, m, loc, scale) result(y) &
            bind(c, name="scifort_crystalball_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: beta !! positive tail-transition magnitude
        real(c_double), value :: m !! power-law exponent, > 1
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(crystalball_cdf(real(x, dp), real(beta, dp), real(m, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_crystalball_cdf_f64

    function scifort_crystalball_ppf_f64(p, beta, m, loc, scale) result(x) &
            bind(c, name="scifort_crystalball_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: beta !! positive tail-transition magnitude
        real(c_double), value :: m !! power-law exponent, > 1
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(crystalball_ppf(real(p, dp), real(beta, dp), real(m, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_crystalball_ppf_f64

    function scifort_jf_skew_t_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_jf_skew_t_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive left-shape parameter
        real(c_double), value :: b !! positive right-shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(jf_skew_t_pdf(real(x, dp), real(a, dp), real(b, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_jf_skew_t_pdf_f64

    function scifort_jf_skew_t_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_jf_skew_t_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive left-shape parameter
        real(c_double), value :: b !! positive right-shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(jf_skew_t_cdf(real(x, dp), real(a, dp), real(b, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_jf_skew_t_cdf_f64

    function scifort_jf_skew_t_ppf_f64(p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_jf_skew_t_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive left-shape parameter
        real(c_double), value :: b !! positive right-shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(jf_skew_t_ppf(real(p, dp), real(a, dp), real(b, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_jf_skew_t_ppf_f64

    function scifort_pearson3_pdf_f64(x, skew, loc, scale) result(y) &
            bind(c, name="scifort_pearson3_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: skew !! finite Pearson III skewness parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(pearson3_pdf(real(x, dp), real(skew, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_pearson3_pdf_f64

    function scifort_pearson3_cdf_f64(x, skew, loc, scale) result(y) &
            bind(c, name="scifort_pearson3_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: skew !! finite Pearson III skewness parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(pearson3_cdf(real(x, dp), real(skew, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_pearson3_cdf_f64

    function scifort_pearson3_ppf_f64(p, skew, loc, scale) result(x) &
            bind(c, name="scifort_pearson3_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: skew !! finite Pearson III skewness parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(pearson3_ppf(real(p, dp), real(skew, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_pearson3_ppf_f64

    function scifort_rel_breitwigner_pdf_f64(x, rho, loc, scale) result(y) &
            bind(c, name="scifort_rel_breitwigner_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: rho !! positive resonance-to-width ratio
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(rel_breitwigner_pdf(real(x, dp), real(rho, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_rel_breitwigner_pdf_f64

    function scifort_rel_breitwigner_cdf_f64(x, rho, loc, scale) result(y) &
            bind(c, name="scifort_rel_breitwigner_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: rho !! positive resonance-to-width ratio
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(rel_breitwigner_cdf(real(x, dp), real(rho, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_rel_breitwigner_cdf_f64

    function scifort_rel_breitwigner_ppf_f64(p, rho, loc, scale) result(x) &
            bind(c, name="scifort_rel_breitwigner_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: rho !! positive resonance-to-width ratio
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(rel_breitwigner_ppf(real(p, dp), real(rho, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_rel_breitwigner_ppf_f64

    function scifort_genexpon_pdf_f64(x, a, b, c, loc, scale) result(y) &
            bind(c, name="scifort_genexpon_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive first rate parameter
        real(c_double), value :: b !! positive second rate parameter
        real(c_double), value :: c !! positive decay parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(genexpon_pdf(real(x, dp), real(a, dp), real(b, dp), real(c, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_genexpon_pdf_f64

    function scifort_genexpon_cdf_f64(x, a, b, c, loc, scale) result(y) &
            bind(c, name="scifort_genexpon_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive first rate parameter
        real(c_double), value :: b !! positive second rate parameter
        real(c_double), value :: c !! positive decay parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(genexpon_cdf(real(x, dp), real(a, dp), real(b, dp), real(c, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_genexpon_cdf_f64

    function scifort_genexpon_ppf_f64(p, a, b, c, loc, scale) result(x) &
            bind(c, name="scifort_genexpon_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive first rate parameter
        real(c_double), value :: b !! positive second rate parameter
        real(c_double), value :: c !! positive decay parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(genexpon_ppf(real(p, dp), real(a, dp), real(b, dp), real(c, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_genexpon_ppf_f64

    function scifort_skewnorm_pdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_skewnorm_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! finite skewness parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(skewnorm_pdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_skewnorm_pdf_f64

    function scifort_skewnorm_cdf_f64(x, a, loc, scale) result(y) &
            bind(c, name="scifort_skewnorm_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! finite skewness parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(skewnorm_cdf(real(x, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_skewnorm_cdf_f64

    function scifort_skewnorm_ppf_f64(p, a, loc, scale) result(x) &
            bind(c, name="scifort_skewnorm_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! finite skewness parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(skewnorm_ppf(real(p, dp), real(a, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_skewnorm_ppf_f64

    function scifort_tukeylambda_pdf_f64(x, lam, loc, scale) result(y) &
            bind(c, name="scifort_tukeylambda_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: lam !! finite Tukey lambda shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(tukeylambda_pdf(real(x, dp), real(lam, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_tukeylambda_pdf_f64

    function scifort_tukeylambda_cdf_f64(x, lam, loc, scale) result(y) &
            bind(c, name="scifort_tukeylambda_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: lam !! finite Tukey lambda shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(tukeylambda_cdf(real(x, dp), real(lam, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_tukeylambda_cdf_f64

    function scifort_tukeylambda_ppf_f64(p, lam, loc, scale) result(x) &
            bind(c, name="scifort_tukeylambda_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: lam !! finite Tukey lambda shape
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(tukeylambda_ppf(real(p, dp), real(lam, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_tukeylambda_ppf_f64

    function scifort_rice_pdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_rice_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: b !! nonnegative Rice shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(rice_pdf(real(x, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rice_pdf_f64

    function scifort_rice_cdf_f64(x, b, loc, scale) result(y) &
            bind(c, name="scifort_rice_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: b !! nonnegative Rice shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(rice_cdf(real(x, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rice_cdf_f64

    function scifort_rice_ppf_f64(p, b, loc, scale) result(x) &
            bind(c, name="scifort_rice_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: b !! nonnegative Rice shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(rice_ppf(real(p, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_rice_ppf_f64


    function scifort_dpareto_lognorm_pdf_f64(x, u, s, a, b, loc, scale) result(y) &
            bind(c, name="scifort_dpareto_lognorm_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: u !! log-location shape
        real(c_double), value :: s !! positive log-scale shape
        real(c_double), value :: a !! positive right-tail shape
        real(c_double), value :: b !! positive left-tail shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive outer scale
        real(c_double) :: y
        y = real(dpareto_lognorm_pdf(real(x, dp), real(u, dp), real(s, dp), &
            real(a, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_dpareto_lognorm_pdf_f64

    function scifort_dpareto_lognorm_cdf_f64(x, u, s, a, b, loc, scale) result(y) &
            bind(c, name="scifort_dpareto_lognorm_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: u !! log-location shape
        real(c_double), value :: s !! positive log-scale shape
        real(c_double), value :: a !! positive right-tail shape
        real(c_double), value :: b !! positive left-tail shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive outer scale
        real(c_double) :: y
        y = real(dpareto_lognorm_cdf(real(x, dp), real(u, dp), real(s, dp), &
            real(a, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_dpareto_lognorm_cdf_f64

    function scifort_dpareto_lognorm_ppf_f64(p, u, s, a, b, loc, scale) result(x) &
            bind(c, name="scifort_dpareto_lognorm_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: u !! log-location shape
        real(c_double), value :: s !! positive log-scale shape
        real(c_double), value :: a !! positive right-tail shape
        real(c_double), value :: b !! positive left-tail shape
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive outer scale
        real(c_double) :: x
        x = real(dpareto_lognorm_ppf(real(p, dp), real(u, dp), real(s, dp), &
            real(a, dp), real(b, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_dpareto_lognorm_ppf_f64

    function scifort_vonmises_pdf_f64(x, kappa, loc, scale) result(y) &
            bind(c, name="scifort_vonmises_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: kappa !! nonnegative concentration
        real(c_double), value :: loc !! circular location
        real(c_double), value :: scale !! positive angular scale
        real(c_double) :: y
        y = real(vonmises_pdf(real(x, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_vonmises_pdf_f64

    function scifort_vonmises_cdf_f64(x, kappa, loc, scale) result(y) &
            bind(c, name="scifort_vonmises_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: kappa !! nonnegative concentration
        real(c_double), value :: loc !! circular location
        real(c_double), value :: scale !! positive angular scale
        real(c_double) :: y
        y = real(vonmises_cdf(real(x, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_vonmises_cdf_f64

    function scifort_vonmises_ppf_f64(p, kappa, loc, scale) result(x) &
            bind(c, name="scifort_vonmises_ppf_f64")
        real(c_double), value :: p !! central-period lower-tail probability in [0,1]
        real(c_double), value :: kappa !! nonnegative concentration
        real(c_double), value :: loc !! circular location
        real(c_double), value :: scale !! positive angular scale
        real(c_double) :: x
        x = real(vonmises_ppf(real(p, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_vonmises_ppf_f64

    function scifort_vonmises_line_pdf_f64(x, kappa, loc, scale) result(y) &
            bind(c, name="scifort_vonmises_line_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: kappa !! nonnegative concentration
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(vonmises_line_pdf(real(x, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_vonmises_line_pdf_f64

    function scifort_vonmises_line_cdf_f64(x, kappa, loc, scale) result(y) &
            bind(c, name="scifort_vonmises_line_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: kappa !! nonnegative concentration
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(vonmises_line_cdf(real(x, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_vonmises_line_cdf_f64

    function scifort_vonmises_line_ppf_f64(p, kappa, loc, scale) result(x) &
            bind(c, name="scifort_vonmises_line_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: kappa !! nonnegative concentration
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(vonmises_line_ppf(real(p, dp), real(kappa, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_vonmises_line_ppf_f64

    function scifort_kstwobign_pdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_kstwobign_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(kstwobign_pdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_kstwobign_pdf_f64

    function scifort_kstwobign_cdf_f64(x, loc, scale) result(y) &
            bind(c, name="scifort_kstwobign_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: y
        y = real(kstwobign_cdf(real(x, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_kstwobign_cdf_f64

    function scifort_kstwobign_ppf_f64(p, loc, scale) result(x) &
            bind(c, name="scifort_kstwobign_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale
        real(c_double) :: x
        x = real(kstwobign_ppf(real(p, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_kstwobign_ppf_f64


    function scifort_irwinhall_pdf_f64(x, n, loc, scale) result(y) &
            bind(c, name="scifort_irwinhall_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: n !! positive integer number of summed uniforms
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(irwinhall_pdf(real(x, dp), real(n, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_irwinhall_pdf_f64

    function scifort_irwinhall_cdf_f64(x, n, loc, scale) result(y) &
            bind(c, name="scifort_irwinhall_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: n !! positive integer number of summed uniforms
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(irwinhall_cdf(real(x, dp), real(n, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_irwinhall_cdf_f64

    function scifort_irwinhall_ppf_f64(p, n, loc, scale) result(x) &
            bind(c, name="scifort_irwinhall_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: n !! positive integer number of summed uniforms
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(irwinhall_ppf(real(p, dp), real(n, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_irwinhall_ppf_f64

    function scifort_ksone_pdf_f64(x, n, loc, scale) result(y) &
            bind(c, name="scifort_ksone_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: n !! positive integer sample-size shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(ksone_pdf(real(x, dp), real(n, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_ksone_pdf_f64

    function scifort_ksone_cdf_f64(x, n, loc, scale) result(y) &
            bind(c, name="scifort_ksone_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: n !! positive integer sample-size shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(ksone_cdf(real(x, dp), real(n, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_ksone_cdf_f64

    function scifort_ksone_ppf_f64(p, n, loc, scale) result(x) &
            bind(c, name="scifort_ksone_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: n !! positive integer sample-size shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(ksone_ppf(real(p, dp), real(n, dp), real(loc, dp), &
            real(scale, dp)), c_double)
    end function scifort_ksone_ppf_f64

    function scifort_kstwo_pdf_f64(x, n, loc, scale) result(y) &
            bind(c, name="scifort_kstwo_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: n !! positive integer sample size
        real(c_double), value :: loc !! lower support location
        real(c_double), value :: scale !! positive support scale
        real(c_double) :: y
        y = real(kstwo_pdf(real(x, dp), real(n, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_kstwo_pdf_f64

    function scifort_kstwo_cdf_f64(x, n, loc, scale) result(y) &
            bind(c, name="scifort_kstwo_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: n !! positive integer sample size
        real(c_double), value :: loc !! lower support location
        real(c_double), value :: scale !! positive support scale
        real(c_double) :: y
        y = real(kstwo_cdf(real(x, dp), real(n, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_kstwo_cdf_f64

    function scifort_kstwo_ppf_f64(p, n, loc, scale) result(x) &
            bind(c, name="scifort_kstwo_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: n !! positive integer sample size
        real(c_double), value :: loc !! lower support location
        real(c_double), value :: scale !! positive support scale
        real(c_double) :: x
        x = real(kstwo_ppf(real(p, dp), real(n, dp), real(loc, dp), real(scale, dp)), c_double)
    end function scifort_kstwo_ppf_f64


    function scifort_levy_stable_pdf_f64(x, alpha, beta, loc, scale) result(y) &
            bind(c, name="scifort_levy_stable_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: alpha !! stability shape parameter in (0,2]
        real(c_double), value :: beta !! skewness shape parameter in [-1,1]
        real(c_double), value :: loc !! S1 location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(levy_stable_pdf(real(x,dp),real(alpha,dp),real(beta,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_levy_stable_pdf_f64

    function scifort_levy_stable_cdf_f64(x, alpha, beta, loc, scale) result(y) &
            bind(c, name="scifort_levy_stable_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: alpha !! stability shape parameter in (0,2]
        real(c_double), value :: beta !! skewness shape parameter in [-1,1]
        real(c_double), value :: loc !! S1 location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(levy_stable_cdf(real(x,dp),real(alpha,dp),real(beta,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_levy_stable_cdf_f64

    function scifort_levy_stable_ppf_f64(p, alpha, beta, loc, scale) result(x) &
            bind(c, name="scifort_levy_stable_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: alpha !! stability shape parameter in (0,2]
        real(c_double), value :: beta !! skewness shape parameter in [-1,1]
        real(c_double), value :: loc !! S1 location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(levy_stable_ppf(real(p,dp),real(alpha,dp),real(beta,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_levy_stable_ppf_f64

    function scifort_studentized_range_pdf_f64(x, k, df, loc, scale) result(y) &
            bind(c, name="scifort_studentized_range_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: k !! number-of-means shape parameter, > 1
        real(c_double), value :: df !! degrees of freedom, > 0
        real(c_double), value :: loc !! lower support location
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(studentized_range_pdf(real(x,dp),real(k,dp),real(df,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_studentized_range_pdf_f64

    function scifort_studentized_range_cdf_f64(x, k, df, loc, scale) result(y) &
            bind(c, name="scifort_studentized_range_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: k !! number-of-means shape parameter, > 1
        real(c_double), value :: df !! degrees of freedom, > 0
        real(c_double), value :: loc !! lower support location
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(studentized_range_cdf(real(x,dp),real(k,dp),real(df,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_studentized_range_cdf_f64

    function scifort_studentized_range_ppf_f64(p, k, df, loc, scale) result(x) &
            bind(c, name="scifort_studentized_range_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: k !! number-of-means shape parameter, > 1
        real(c_double), value :: df !! degrees of freedom, > 0
        real(c_double), value :: loc !! lower support location
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(studentized_range_ppf(real(p,dp),real(k,dp),real(df,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_studentized_range_ppf_f64

    function scifort_ncx2_pdf_f64(x, df, nc, loc, scale) result(y) &
            bind(c, name="scifort_ncx2_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: df !! positive degrees of freedom
        real(c_double), value :: nc !! nonnegative noncentrality parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(ncx2_pdf(real(x, dp), real(df, dp), real(nc, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_ncx2_pdf_f64

    function scifort_ncx2_cdf_f64(x, df, nc, loc, scale) result(y) &
            bind(c, name="scifort_ncx2_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: df !! positive degrees of freedom
        real(c_double), value :: nc !! nonnegative noncentrality parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(ncx2_cdf(real(x, dp), real(df, dp), real(nc, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_ncx2_cdf_f64

    function scifort_ncx2_ppf_f64(p, df, nc, loc, scale) result(x) &
            bind(c, name="scifort_ncx2_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: df !! positive degrees of freedom
        real(c_double), value :: nc !! nonnegative noncentrality parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(ncx2_ppf(real(p, dp), real(df, dp), real(nc, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_ncx2_ppf_f64

    function scifort_ncf_pdf_f64(x, dfn, dfd, nc, loc, scale) result(y) &
            bind(c, name="scifort_ncf_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: dfn !! positive numerator degrees of freedom
        real(c_double), value :: dfd !! positive denominator degrees of freedom
        real(c_double), value :: nc !! nonnegative noncentrality parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(ncf_pdf(real(x, dp), real(dfn, dp), real(dfd, dp), real(nc, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_ncf_pdf_f64

    function scifort_ncf_cdf_f64(x, dfn, dfd, nc, loc, scale) result(y) &
            bind(c, name="scifort_ncf_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: dfn !! positive numerator degrees of freedom
        real(c_double), value :: dfd !! positive denominator degrees of freedom
        real(c_double), value :: nc !! nonnegative noncentrality parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(ncf_cdf(real(x, dp), real(dfn, dp), real(dfd, dp), real(nc, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_ncf_cdf_f64

    function scifort_ncf_ppf_f64(p, dfn, dfd, nc, loc, scale) result(x) &
            bind(c, name="scifort_ncf_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: dfn !! positive numerator degrees of freedom
        real(c_double), value :: dfd !! positive denominator degrees of freedom
        real(c_double), value :: nc !! nonnegative noncentrality parameter
        real(c_double), value :: loc !! lower support endpoint
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(ncf_ppf(real(p, dp), real(dfn, dp), real(dfd, dp), real(nc, dp), &
            real(loc, dp), real(scale, dp)), c_double)
    end function scifort_ncf_ppf_f64

    function scifort_randint_pmf_f64(k, low, high, loc) result(y) &
            bind(c, name="scifort_randint_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: low !! integer lower bound, inclusive
        real(c_double), value :: high !! integer upper bound, exclusive
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(randint_pmf(real(k, dp), real(low, dp), real(high, dp), &
            real(loc, dp)), c_double)
    end function scifort_randint_pmf_f64

    function scifort_randint_cdf_f64(k, low, high, loc) result(y) &
            bind(c, name="scifort_randint_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: low !! integer lower bound, inclusive
        real(c_double), value :: high !! integer upper bound, exclusive
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(randint_cdf(real(k, dp), real(low, dp), real(high, dp), &
            real(loc, dp)), c_double)
    end function scifort_randint_cdf_f64

    function scifort_randint_ppf_f64(p, low, high, loc) result(k) &
            bind(c, name="scifort_randint_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: low !! integer lower bound, inclusive
        real(c_double), value :: high !! integer upper bound, exclusive
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(randint_ppf(real(p, dp), real(low, dp), real(high, dp), &
            real(loc, dp)), c_double)
    end function scifort_randint_ppf_f64

    function scifort_planck_pmf_f64(k, lambda, loc) result(y) &
            bind(c, name="scifort_planck_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: lambda !! positive exponential rate
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(planck_pmf(real(k, dp), real(lambda, dp), real(loc, dp)), c_double)
    end function scifort_planck_pmf_f64

    function scifort_planck_cdf_f64(k, lambda, loc) result(y) &
            bind(c, name="scifort_planck_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: lambda !! positive exponential rate
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(planck_cdf(real(k, dp), real(lambda, dp), real(loc, dp)), c_double)
    end function scifort_planck_cdf_f64

    function scifort_planck_ppf_f64(p, lambda, loc) result(k) &
            bind(c, name="scifort_planck_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: lambda !! positive exponential rate
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(planck_ppf(real(p, dp), real(lambda, dp), real(loc, dp)), c_double)
    end function scifort_planck_ppf_f64

    function scifort_dlaplace_pmf_f64(k, a, loc) result(y) &
            bind(c, name="scifort_dlaplace_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: a !! positive decay parameter
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(dlaplace_pmf(real(k, dp), real(a, dp), real(loc, dp)), c_double)
    end function scifort_dlaplace_pmf_f64

    function scifort_dlaplace_cdf_f64(k, a, loc) result(y) &
            bind(c, name="scifort_dlaplace_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: a !! positive decay parameter
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(dlaplace_cdf(real(k, dp), real(a, dp), real(loc, dp)), c_double)
    end function scifort_dlaplace_cdf_f64

    function scifort_dlaplace_ppf_f64(p, a, loc) result(k) &
            bind(c, name="scifort_dlaplace_ppf_f64")
        real(c_double), value :: p !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive decay parameter
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(dlaplace_ppf(real(p, dp), real(a, dp), real(loc, dp)), c_double)
    end function scifort_dlaplace_ppf_f64

    function scifort_logser_pmf_f64(k, p, loc) result(y) &
            bind(c, name="scifort_logser_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: p !! shape probability in (0,1)
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(logser_pmf(real(k, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_logser_pmf_f64

    function scifort_logser_cdf_f64(k, p, loc) result(y) &
            bind(c, name="scifort_logser_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: p !! shape probability in (0,1)
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(logser_cdf(real(k, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_logser_cdf_f64

    function scifort_logser_ppf_f64(q, p, loc) result(k) &
            bind(c, name="scifort_logser_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: p !! shape probability in (0,1)
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(logser_ppf(real(q, dp), real(p, dp), real(loc, dp)), c_double)
    end function scifort_logser_ppf_f64


    function scifort_betabinom_pmf_f64(k, n, a, b, loc) result(y) &
            bind(c, name="scifort_betabinom_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: n !! nonnegative integer number of trials
        real(c_double), value :: a !! positive first beta shape
        real(c_double), value :: b !! positive second beta shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(betabinom_pmf(real(k, dp), real(n, dp), real(a, dp), real(b, dp), real(loc, dp)), c_double)
    end function scifort_betabinom_pmf_f64

    function scifort_betabinom_cdf_f64(k, n, a, b, loc) result(y) &
            bind(c, name="scifort_betabinom_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: n !! nonnegative integer number of trials
        real(c_double), value :: a !! positive first beta shape
        real(c_double), value :: b !! positive second beta shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(betabinom_cdf(real(k, dp), real(n, dp), real(a, dp), real(b, dp), real(loc, dp)), c_double)
    end function scifort_betabinom_cdf_f64

    function scifort_betabinom_ppf_f64(q, n, a, b, loc) result(k) &
            bind(c, name="scifort_betabinom_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: n !! nonnegative integer number of trials
        real(c_double), value :: a !! positive first beta shape
        real(c_double), value :: b !! positive second beta shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(betabinom_ppf(real(q, dp), real(n, dp), real(a, dp), real(b, dp), real(loc, dp)), c_double)
    end function scifort_betabinom_ppf_f64

    function scifort_hypergeom_pmf_f64(k, m, n, draws, loc) result(y) &
            bind(c, name="scifort_hypergeom_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer Type-I count
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(hypergeom_pmf(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), real(loc, dp)), c_double)
    end function scifort_hypergeom_pmf_f64

    function scifort_hypergeom_cdf_f64(k, m, n, draws, loc) result(y) &
            bind(c, name="scifort_hypergeom_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer Type-I count
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(hypergeom_cdf(real(k, dp), real(m, dp), real(n, dp), real(draws, dp), real(loc, dp)), c_double)
    end function scifort_hypergeom_cdf_f64

    function scifort_hypergeom_ppf_f64(q, m, n, draws, loc) result(k) &
            bind(c, name="scifort_hypergeom_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer Type-I count
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(hypergeom_ppf(real(q, dp), real(m, dp), real(n, dp), real(draws, dp), real(loc, dp)), c_double)
    end function scifort_hypergeom_ppf_f64

    function scifort_nhypergeom_pmf_f64(k, m, n, r, loc) result(y) &
            bind(c, name="scifort_nhypergeom_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: m !! nonnegative integer population size
        real(c_double), value :: n !! integer success count
        real(c_double), value :: r !! integer failures required
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(nhypergeom_pmf(real(k, dp), real(m, dp), real(n, dp), real(r, dp), real(loc, dp)), c_double)
    end function scifort_nhypergeom_pmf_f64

    function scifort_nhypergeom_cdf_f64(k, m, n, r, loc) result(y) &
            bind(c, name="scifort_nhypergeom_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: m !! nonnegative integer population size
        real(c_double), value :: n !! integer success count
        real(c_double), value :: r !! integer failures required
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(nhypergeom_cdf(real(k, dp), real(m, dp), real(n, dp), real(r, dp), real(loc, dp)), c_double)
    end function scifort_nhypergeom_cdf_f64

    function scifort_nhypergeom_ppf_f64(q, m, n, r, loc) result(k) &
            bind(c, name="scifort_nhypergeom_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: m !! nonnegative integer population size
        real(c_double), value :: n !! integer success count
        real(c_double), value :: r !! integer failures required
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(nhypergeom_ppf(real(q, dp), real(m, dp), real(n, dp), real(r, dp), real(loc, dp)), c_double)
    end function scifort_nhypergeom_ppf_f64

    function scifort_boltzmann_pmf_f64(k, lambda, n, loc) result(y) &
            bind(c, name="scifort_boltzmann_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: lambda !! positive exponential rate
        real(c_double), value :: n !! positive integer support size
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(boltzmann_pmf(real(k, dp), real(lambda, dp), real(n, dp), real(loc, dp)), c_double)
    end function scifort_boltzmann_pmf_f64

    function scifort_boltzmann_cdf_f64(k, lambda, n, loc) result(y) &
            bind(c, name="scifort_boltzmann_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: lambda !! positive exponential rate
        real(c_double), value :: n !! positive integer support size
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(boltzmann_cdf(real(k, dp), real(lambda, dp), real(n, dp), real(loc, dp)), c_double)
    end function scifort_boltzmann_cdf_f64

    function scifort_boltzmann_ppf_f64(q, lambda, n, loc) result(k) &
            bind(c, name="scifort_boltzmann_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: lambda !! positive exponential rate
        real(c_double), value :: n !! positive integer support size
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(boltzmann_ppf(real(q, dp), real(lambda, dp), real(n, dp), real(loc, dp)), c_double)
    end function scifort_boltzmann_ppf_f64


    function scifort_betanbinom_pmf_f64(k, n, a, b, loc) result(y) &
            bind(c, name="scifort_betanbinom_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: n !! positive integer success count
        real(c_double), value :: a !! positive first beta shape
        real(c_double), value :: b !! positive second beta shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(betanbinom_pmf(real(k, dp), real(n, dp), real(a, dp), real(b, dp), real(loc, dp)), c_double)
    end function scifort_betanbinom_pmf_f64

    function scifort_betanbinom_cdf_f64(k, n, a, b, loc) result(y) &
            bind(c, name="scifort_betanbinom_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: n !! positive integer success count
        real(c_double), value :: a !! positive first beta shape
        real(c_double), value :: b !! positive second beta shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(betanbinom_cdf(real(k, dp), real(n, dp), real(a, dp), real(b, dp), real(loc, dp)), c_double)
    end function scifort_betanbinom_cdf_f64

    function scifort_betanbinom_ppf_f64(q, n, a, b, loc) result(k) &
            bind(c, name="scifort_betanbinom_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: n !! positive integer success count
        real(c_double), value :: a !! positive first beta shape
        real(c_double), value :: b !! positive second beta shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(betanbinom_ppf(real(q, dp), real(n, dp), real(a, dp), real(b, dp), real(loc, dp)), c_double)
    end function scifort_betanbinom_ppf_f64

    function scifort_yulesimon_pmf_f64(k, alpha, loc) result(y) &
            bind(c, name="scifort_yulesimon_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: alpha !! positive Yule-Simon shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(yulesimon_pmf(real(k, dp), real(alpha, dp), real(loc, dp)), c_double)
    end function scifort_yulesimon_pmf_f64

    function scifort_yulesimon_cdf_f64(k, alpha, loc) result(y) &
            bind(c, name="scifort_yulesimon_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: alpha !! positive Yule-Simon shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(yulesimon_cdf(real(k, dp), real(alpha, dp), real(loc, dp)), c_double)
    end function scifort_yulesimon_cdf_f64

    function scifort_yulesimon_ppf_f64(q, alpha, loc) result(k) &
            bind(c, name="scifort_yulesimon_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: alpha !! positive Yule-Simon shape
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(yulesimon_ppf(real(q, dp), real(alpha, dp), real(loc, dp)), c_double)
    end function scifort_yulesimon_ppf_f64

    function scifort_zipf_pmf_f64(k, a, loc) result(y) bind(c, name="scifort_zipf_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: a !! power exponent, strictly greater than one
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(zipf_pmf(real(k, dp), real(a, dp), real(loc, dp)), c_double)
    end function scifort_zipf_pmf_f64

    function scifort_zipf_cdf_f64(k, a, loc) result(y) bind(c, name="scifort_zipf_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: a !! power exponent, strictly greater than one
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(zipf_cdf(real(k, dp), real(a, dp), real(loc, dp)), c_double)
    end function scifort_zipf_cdf_f64

    function scifort_zipf_ppf_f64(q, a, loc) result(k) bind(c, name="scifort_zipf_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: a !! power exponent, strictly greater than one
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(zipf_ppf(real(q, dp), real(a, dp), real(loc, dp)), c_double)
    end function scifort_zipf_ppf_f64

    function scifort_zipfian_pmf_f64(k, a, n, loc) result(y) bind(c, name="scifort_zipfian_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: a !! nonnegative power exponent
        real(c_double), value :: n !! positive integer upper support bound
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(zipfian_pmf(real(k, dp), real(a, dp), real(n, dp), real(loc, dp)), c_double)
    end function scifort_zipfian_pmf_f64

    function scifort_zipfian_cdf_f64(k, a, n, loc) result(y) bind(c, name="scifort_zipfian_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: a !! nonnegative power exponent
        real(c_double), value :: n !! positive integer upper support bound
        real(c_double), value :: loc !! real support shift
        real(c_double) :: y
        y = real(zipfian_cdf(real(k, dp), real(a, dp), real(n, dp), real(loc, dp)), c_double)
    end function scifort_zipfian_cdf_f64

    function scifort_zipfian_ppf_f64(q, a, n, loc) result(k) bind(c, name="scifort_zipfian_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: a !! nonnegative power exponent
        real(c_double), value :: n !! positive integer upper support bound
        real(c_double), value :: loc !! real support shift
        real(c_double) :: k
        k = real(zipfian_ppf(real(q, dp), real(a, dp), real(n, dp), real(loc, dp)), c_double)
    end function scifort_zipfian_ppf_f64

    function scifort_geninvgauss_pdf_f64(x, p, b, loc, scale) result(y) &
            bind(c, name="scifort_geninvgauss_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: p !! real shape parameter
        real(c_double), value :: b !! strictly positive shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(geninvgauss_pdf(real(x,dp),real(p,dp),real(b,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_geninvgauss_pdf_f64

    function scifort_geninvgauss_cdf_f64(x, p, b, loc, scale) result(y) &
            bind(c, name="scifort_geninvgauss_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: p !! real shape parameter
        real(c_double), value :: b !! strictly positive shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(geninvgauss_cdf(real(x,dp),real(p,dp),real(b,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_geninvgauss_cdf_f64

    function scifort_geninvgauss_ppf_f64(q, p, b, loc, scale) result(x) &
            bind(c, name="scifort_geninvgauss_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: p !! real shape parameter
        real(c_double), value :: b !! strictly positive shape parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(geninvgauss_ppf(real(q,dp),real(p,dp),real(b,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_geninvgauss_ppf_f64

    function scifort_norminvgauss_pdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_norminvgauss_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive tail shape
        real(c_double), value :: b !! skew shape with abs(b) < a
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(norminvgauss_pdf(real(x,dp),real(a,dp),real(b,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_norminvgauss_pdf_f64

    function scifort_norminvgauss_cdf_f64(x, a, b, loc, scale) result(y) &
            bind(c, name="scifort_norminvgauss_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive tail shape
        real(c_double), value :: b !! skew shape with abs(b) < a
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(norminvgauss_cdf(real(x,dp),real(a,dp),real(b,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_norminvgauss_cdf_f64

    function scifort_norminvgauss_ppf_f64(q, a, b, loc, scale) result(x) &
            bind(c, name="scifort_norminvgauss_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive tail shape
        real(c_double), value :: b !! skew shape with abs(b) < a
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(norminvgauss_ppf(real(q,dp),real(a,dp),real(b,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_norminvgauss_ppf_f64

    function scifort_skellam_pmf_f64(k, mu1, mu2, loc) result(y) &
            bind(c, name="scifort_skellam_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: mu1 !! positive first Poisson mean
        real(c_double), value :: mu2 !! positive second Poisson mean
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        y = real(skellam_pmf(real(k,dp),real(mu1,dp),real(mu2,dp),real(loc,dp)),c_double)
    end function scifort_skellam_pmf_f64

    function scifort_skellam_cdf_f64(k, mu1, mu2, loc) result(y) &
            bind(c, name="scifort_skellam_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: mu1 !! positive first Poisson mean
        real(c_double), value :: mu2 !! positive second Poisson mean
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        y = real(skellam_cdf(real(k,dp),real(mu1,dp),real(mu2,dp),real(loc,dp)),c_double)
    end function scifort_skellam_cdf_f64

    function scifort_skellam_ppf_f64(q, mu1, mu2, loc) result(k) &
            bind(c, name="scifort_skellam_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: mu1 !! positive first Poisson mean
        real(c_double), value :: mu2 !! positive second Poisson mean
        real(c_double), value :: loc !! support shift
        real(c_double) :: k
        k = real(skellam_ppf(real(q,dp),real(mu1,dp),real(mu2,dp),real(loc,dp)),c_double)
    end function scifort_skellam_ppf_f64


    function scifort_genhyperbolic_pdf_f64(x, p, a, b, loc, scale) result(y) &
            bind(c, name="scifort_genhyperbolic_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: p !! real tail shape parameter
        real(c_double), value :: a !! positive shape parameter
        real(c_double), value :: b !! skew shape constrained by a and p
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(genhyperbolic_pdf(real(x,dp),real(p,dp),real(a,dp),real(b,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_genhyperbolic_pdf_f64

    function scifort_genhyperbolic_cdf_f64(x, p, a, b, loc, scale) result(y) &
            bind(c, name="scifort_genhyperbolic_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: p !! real tail shape parameter
        real(c_double), value :: a !! positive shape parameter
        real(c_double), value :: b !! skew shape constrained by a and p
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(genhyperbolic_cdf(real(x,dp),real(p,dp),real(a,dp),real(b,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_genhyperbolic_cdf_f64

    function scifort_genhyperbolic_ppf_f64(q, p, a, b, loc, scale) result(x) &
            bind(c, name="scifort_genhyperbolic_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: p !! real tail shape parameter
        real(c_double), value :: a !! positive shape parameter
        real(c_double), value :: b !! skew shape constrained by a and p
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(genhyperbolic_ppf(real(q,dp),real(p,dp),real(a,dp),real(b,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_genhyperbolic_ppf_f64

    function scifort_nchypergeom_fisher_pmf_f64(k, m, n, draws, odds, loc) result(y) &
            bind(c, name="scifort_nchypergeom_fisher_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer number of Type-I objects
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: odds !! strictly positive Type-I odds ratio
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        y = real(nchypergeom_fisher_pmf(real(k,dp),real(m,dp),real(n,dp),real(draws,dp), &
            real(odds,dp),real(loc,dp)),c_double)
    end function scifort_nchypergeom_fisher_pmf_f64

    function scifort_nchypergeom_fisher_cdf_f64(k, m, n, draws, odds, loc) result(y) &
            bind(c, name="scifort_nchypergeom_fisher_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer number of Type-I objects
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: odds !! strictly positive Type-I odds ratio
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        y = real(nchypergeom_fisher_cdf(real(k,dp),real(m,dp),real(n,dp),real(draws,dp), &
            real(odds,dp),real(loc,dp)),c_double)
    end function scifort_nchypergeom_fisher_cdf_f64

    function scifort_nchypergeom_fisher_ppf_f64(q, m, n, draws, odds, loc) result(k) &
            bind(c, name="scifort_nchypergeom_fisher_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer number of Type-I objects
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: odds !! strictly positive Type-I odds ratio
        real(c_double), value :: loc !! support shift
        real(c_double) :: k
        k = real(nchypergeom_fisher_ppf(real(q,dp),real(m,dp),real(n,dp),real(draws,dp), &
            real(odds,dp),real(loc,dp)),c_double)
    end function scifort_nchypergeom_fisher_ppf_f64

    function scifort_nct_pdf_f64(x, df, nc, loc, scale) result(y) &
            bind(c, name="scifort_nct_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: df !! positive degrees of freedom
        real(c_double), value :: nc !! noncentrality parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(nct_pdf(real(x,dp),real(df,dp),real(nc,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_nct_pdf_f64

    function scifort_nct_cdf_f64(x, df, nc, loc, scale) result(y) &
            bind(c, name="scifort_nct_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: df !! positive degrees of freedom
        real(c_double), value :: nc !! noncentrality parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(nct_cdf(real(x,dp),real(df,dp),real(nc,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_nct_cdf_f64

    function scifort_nct_ppf_f64(q, df, nc, loc, scale) result(x) &
            bind(c, name="scifort_nct_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: df !! positive degrees of freedom
        real(c_double), value :: nc !! noncentrality parameter
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(nct_ppf(real(q,dp),real(df,dp),real(nc,dp),real(loc,dp),real(scale,dp)),c_double)
    end function scifort_nct_ppf_f64

    function scifort_gausshyper_pdf_f64(x, a, b, c, z, loc, scale) result(y) &
            bind(c, name="scifort_gausshyper_pdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive left beta shape
        real(c_double), value :: b !! positive right beta shape
        real(c_double), value :: c !! real hypergeometric exponent
        real(c_double), value :: z !! tilt parameter greater than -1
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(gausshyper_pdf(real(x,dp),real(a,dp),real(b,dp),real(c,dp),real(z,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_gausshyper_pdf_f64

    function scifort_gausshyper_cdf_f64(x, a, b, c, z, loc, scale) result(y) &
            bind(c, name="scifort_gausshyper_cdf_f64")
        real(c_double), value :: x !! evaluation point
        real(c_double), value :: a !! positive left beta shape
        real(c_double), value :: b !! positive right beta shape
        real(c_double), value :: c !! real hypergeometric exponent
        real(c_double), value :: z !! tilt parameter greater than -1
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: y
        y = real(gausshyper_cdf(real(x,dp),real(a,dp),real(b,dp),real(c,dp),real(z,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_gausshyper_cdf_f64

    function scifort_gausshyper_ppf_f64(q, a, b, c, z, loc, scale) result(x) &
            bind(c, name="scifort_gausshyper_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: a !! positive left beta shape
        real(c_double), value :: b !! positive right beta shape
        real(c_double), value :: c !! real hypergeometric exponent
        real(c_double), value :: z !! tilt parameter greater than -1
        real(c_double), value :: loc !! location parameter
        real(c_double), value :: scale !! positive scale parameter
        real(c_double) :: x
        x = real(gausshyper_ppf(real(q,dp),real(a,dp),real(b,dp),real(c,dp),real(z,dp), &
            real(loc,dp),real(scale,dp)),c_double)
    end function scifort_gausshyper_ppf_f64

    function scifort_landau_pdf_f64(x, loc, scale) result(y) bind(c, name="scifort_landau_pdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the density is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y

        y = real(landau_pdf(real(x,dp), real(loc,dp), real(scale,dp)), c_double)
    end function scifort_landau_pdf_f64

    function scifort_landau_cdf_f64(x, loc, scale) result(y) bind(c, name="scifort_landau_cdf_f64")
        real(c_double), value, intent(in) :: x !! point at which the CDF is evaluated
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: y

        y = real(landau_cdf(real(x,dp), real(loc,dp), real(scale,dp)), c_double)
    end function scifort_landau_cdf_f64

    function scifort_landau_ppf_f64(q, loc, scale) result(x) bind(c, name="scifort_landau_ppf_f64")
        real(c_double), value, intent(in) :: q !! lower-tail probability
        real(c_double), value, intent(in) :: loc !! location parameter
        real(c_double), value, intent(in) :: scale !! positive scale parameter
        real(c_double) :: x

        x = real(landau_ppf(real(q,dp), real(loc,dp), real(scale,dp)), c_double)
    end function scifort_landau_ppf_f64

    function scifort_nchypergeom_wallenius_pmf_f64(k,m,n,draws,odds,loc) result(y) &
            bind(c,name="scifort_nchypergeom_wallenius_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer number of Type-I objects
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: odds !! strictly positive Type-I odds ratio
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        y=real(nchypergeom_wallenius_pmf(real(k,dp),real(m,dp),real(n,dp),real(draws,dp), &
            real(odds,dp),real(loc,dp)),c_double)
    end function scifort_nchypergeom_wallenius_pmf_f64

    function scifort_nchypergeom_wallenius_cdf_f64(k,m,n,draws,odds,loc) result(y) &
            bind(c,name="scifort_nchypergeom_wallenius_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer number of Type-I objects
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: odds !! strictly positive Type-I odds ratio
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        y=real(nchypergeom_wallenius_cdf(real(k,dp),real(m,dp),real(n,dp),real(draws,dp), &
            real(odds,dp),real(loc,dp)),c_double)
    end function scifort_nchypergeom_wallenius_cdf_f64

    function scifort_nchypergeom_wallenius_ppf_f64(q,m,n,draws,odds,loc) result(k) &
            bind(c,name="scifort_nchypergeom_wallenius_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), value :: m !! positive integer population size
        real(c_double), value :: n !! integer number of Type-I objects
        real(c_double), value :: draws !! integer sample size
        real(c_double), value :: odds !! strictly positive Type-I odds ratio
        real(c_double), value :: loc !! support shift
        real(c_double) :: k
        k=real(nchypergeom_wallenius_ppf(real(q,dp),real(m,dp),real(n,dp),real(draws,dp), &
            real(odds,dp),real(loc,dp)),c_double)
    end function scifort_nchypergeom_wallenius_ppf_f64

    function scifort_poisson_binom_pmf_f64(k,p,n,loc) result(y) bind(c,name="scifort_poisson_binom_pmf_f64")
        real(c_double), value :: k !! lattice point
        real(c_double), intent(in) :: p(*) !! Bernoulli success probabilities
        integer(c_size_t), value :: n !! number of Bernoulli probabilities
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        real(dp), allocatable :: pp(:)
        integer :: nn
        nn=int(n); allocate(pp(nn)); pp=real(p(1:nn),dp)
        y=real(poisson_binom_pmf(real(k,dp),pp,real(loc,dp)),c_double)
    end function scifort_poisson_binom_pmf_f64

    function scifort_poisson_binom_cdf_f64(k,p,n,loc) result(y) bind(c,name="scifort_poisson_binom_cdf_f64")
        real(c_double), value :: k !! evaluation point
        real(c_double), intent(in) :: p(*) !! Bernoulli success probabilities
        integer(c_size_t), value :: n !! number of Bernoulli probabilities
        real(c_double), value :: loc !! support shift
        real(c_double) :: y
        real(dp), allocatable :: pp(:)
        integer :: nn
        nn=int(n); allocate(pp(nn)); pp=real(p(1:nn),dp)
        y=real(poisson_binom_cdf(real(k,dp),pp,real(loc,dp)),c_double)
    end function scifort_poisson_binom_cdf_f64

    function scifort_poisson_binom_ppf_f64(q,p,n,loc) result(k) bind(c,name="scifort_poisson_binom_ppf_f64")
        real(c_double), value :: q !! lower-tail probability in [0,1]
        real(c_double), intent(in) :: p(*) !! Bernoulli success probabilities
        integer(c_size_t), value :: n !! number of Bernoulli probabilities
        real(c_double), value :: loc !! support shift
        real(c_double) :: k
        real(dp), allocatable :: pp(:)
        integer :: nn
        nn=int(n); allocate(pp(nn)); pp=real(p(1:nn),dp)
        k=real(poisson_binom_ppf(real(q,dp),pp,real(loc,dp)),c_double)
    end function scifort_poisson_binom_ppf_f64

end module scifort_c_api
