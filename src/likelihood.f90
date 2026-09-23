! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_likelihood
    use scifort_kinds, only : dp
    use scifort_normal, only : normal_logpdf
    use scifort_uniform, only : uniform_logpdf
    use scifort_exponential, only : exponential_logpdf
    use scifort_laplace, only : laplace_logpdf
    use scifort_logistic, only : logistic_logpdf
    use scifort_cauchy, only : cauchy_logpdf
    use scifort_rayleigh, only : rayleigh_logpdf
    use scifort_gamma, only : gamma_logpdf
    use scifort_chi2, only : chi2_logpdf
    use scifort_student_t, only : t_logpdf
    use scifort_lognormal, only : lognormal_logpdf
    use scifort_weibull, only : weibull_logpdf
    use scifort_pareto, only : pareto_logpdf
    use scifort_beta, only : beta_logpdf
    use scifort_f_distribution, only : f_logpdf
    use scifort_gumbel_r, only : gumbel_r_logpdf
    use scifort_gumbel_l, only : gumbel_l_logpdf
    use scifort_powerlaw, only : powerlaw_logpdf
    use scifort_triang, only : triang_logpdf
    use scifort_genpareto, only : genpareto_logpdf
    use scifort_arcsine, only : arcsine_logpdf
    use scifort_halfnorm, only : halfnorm_logpdf
    use scifort_halfcauchy, only : halfcauchy_logpdf
    use scifort_lomax, only : lomax_logpdf
    use scifort_chi, only : chi_logpdf
    use scifort_maxwell, only : maxwell_logpdf
    use scifort_cosine, only : cosine_logpdf
    use scifort_semicircular, only : semicircular_logpdf
    use scifort_anglit, only : anglit_logpdf
    use scifort_moyal, only : moyal_logpdf
    use scifort_landau, only : landau_logpdf
    use scifort_hypsecant, only : hypsecant_logpdf
    use scifort_halflogistic, only : halflogistic_logpdf
    use scifort_invgamma, only : invgamma_logpdf
    use scifort_invgauss, only : invgauss_logpdf
    use scifort_levy, only : levy_logpdf
    use scifort_loglaplace, only : loglaplace_logpdf
    use scifort_bradford, only : bradford_logpdf
    use scifort_truncexpon, only : truncexpon_logpdf
    use scifort_fisk, only : fisk_logpdf
    use scifort_dweibull, only : dweibull_logpdf
    use scifort_alpha, only : alpha_logpdf
    use scifort_fatiguelife, only : fatiguelife_logpdf
    use scifort_genlogistic, only : genlogistic_logpdf
    use scifort_gennorm, only : gennorm_logpdf
    use scifort_nakagami, only : nakagami_logpdf
    use scifort_powernorm, only : powernorm_logpdf
    use scifort_loggamma, only : loggamma_logpdf
    use scifort_wald, only : wald_logpdf
    use scifort_gompertz, only : gompertz_logpdf
    use scifort_invweibull, only : invweibull_logpdf
    use scifort_betaprime, only : betaprime_logpdf
    use scifort_burr12, only : burr12_logpdf
    use scifort_genhalflogistic, only : genhalflogistic_logpdf
    use scifort_exponpow, only : exponpow_logpdf
    use scifort_exponweib, only : exponweib_logpdf
    use scifort_powerlognorm, only : powerlognorm_logpdf
    use scifort_levy_l, only : levy_l_logpdf
    use scifort_weibull_max, only : weibull_max_logpdf
    use scifort_rdist, only : rdist_logpdf
    use scifort_skewcauchy, only : skewcauchy_logpdf
    use scifort_dgamma, only : dgamma_logpdf
    use scifort_laplace_asymmetric, only : laplace_asymmetric_logpdf
    use scifort_truncnorm, only : truncnorm_logpdf
    use scifort_loguniform, only : loguniform_logpdf
    use scifort_foldnorm, only : foldnorm_logpdf
    use scifort_foldcauchy, only : foldcauchy_logpdf
    use scifort_recipinvgauss, only : recipinvgauss_logpdf
    use scifort_truncpareto, only : truncpareto_logpdf
    use scifort_exponnorm, only : exponnorm_logpdf
    use scifort_johnsonsb, only : johnsonsb_logpdf
    use scifort_johnsonsu, only : johnsonsu_logpdf
    use scifort_trapezoid, only : trapezoid_logpdf
    use scifort_burr, only : burr_logpdf
    use scifort_mielke, only : mielke_logpdf
    use scifort_gibrat, only : gibrat_logpdf
    use scifort_wrapcauchy, only : wrapcauchy_logpdf
    use scifort_genextreme, only : genextreme_logpdf
    use scifort_kappa3, only : kappa3_logpdf
    use scifort_kappa4, only : kappa4_logpdf
    use scifort_truncweibull_min, only : truncweibull_min_logpdf
    use scifort_gengamma, only : gengamma_logpdf
    use scifort_halfgennorm, only : halfgennorm_logpdf
    use scifort_argus, only : argus_logpdf
    use scifort_erlang, only : erlang_logpdf
    use scifort_crystalball, only : crystalball_logpdf
    use scifort_jf_skew_t, only : jf_skew_t_logpdf
    use scifort_pearson3, only : pearson3_logpdf
    use scifort_rel_breitwigner, only : rel_breitwigner_logpdf
    use scifort_genexpon, only : genexpon_logpdf
    use scifort_skewnorm, only : skewnorm_logpdf
    use scifort_tukeylambda, only : tukeylambda_logpdf
    use scifort_rice, only : rice_logpdf
    use scifort_dpareto_lognorm, only : dpareto_lognorm_logpdf
    use scifort_vonmises, only : vonmises_logpdf
    use scifort_vonmises_line, only : vonmises_line_logpdf
    use scifort_kstwobign, only : kstwobign_logpdf
    use scifort_irwinhall, only : irwinhall_logpdf
    use scifort_ksone, only : ksone_logpdf
    use scifort_kstwo, only : kstwo_logpdf
    use scifort_levy_stable, only : levy_stable_logpdf
    use scifort_studentized_range, only : studentized_range_logpdf
    use scifort_ncx2, only : ncx2_logpdf
    use scifort_ncf, only : ncf_logpdf
    use scifort_bernoulli, only : bernoulli_logpmf
    use scifort_poisson, only : poisson_logpmf
    use scifort_geometric, only : geometric_logpmf
    use scifort_binomial, only : binomial_logpmf
    use scifort_negative_binomial, only : negative_binomial_logpmf
    use scifort_randint, only : randint_logpmf
    use scifort_planck, only : planck_logpmf
    use scifort_dlaplace, only : dlaplace_logpmf
    use scifort_logser, only : logser_logpmf
    use scifort_betabinom, only : betabinom_logpmf
    use scifort_hypergeom, only : hypergeom_logpmf
    use scifort_nhypergeom, only : nhypergeom_logpmf
    use scifort_boltzmann, only : boltzmann_logpmf
    use scifort_betanbinom, only : betanbinom_logpmf
    use scifort_yulesimon, only : yulesimon_logpmf
    use scifort_zipf, only : zipf_logpmf
    use scifort_zipfian, only : zipfian_logpmf
    use scifort_geninvgauss, only : geninvgauss_logpdf
    use scifort_norminvgauss, only : norminvgauss_logpdf
    use scifort_skellam, only : skellam_logpmf
    use scifort_genhyperbolic, only : genhyperbolic_logpdf
    use scifort_nchypergeom_fisher, only : nchypergeom_fisher_logpmf
    use scifort_nchypergeom_wallenius, only : nchypergeom_wallenius_logpmf
    use scifort_poisson_binom, only : poisson_binom_logpmf
    use scifort_nct, only : nct_logpdf
    use scifort_gausshyper, only : gausshyper_logpdf
    implicit none
    private

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
    public :: landau_loglikelihood
    public :: landau_nnlf
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
    public :: bradford_loglikelihood
    public :: bradford_nnlf
    public :: truncexpon_loglikelihood
    public :: truncexpon_nnlf
    public :: fisk_loglikelihood
    public :: fisk_nnlf
    public :: dweibull_loglikelihood
    public :: dweibull_nnlf
    public :: alpha_loglikelihood
    public :: alpha_nnlf
    public :: fatiguelife_loglikelihood
    public :: fatiguelife_nnlf
    public :: genlogistic_loglikelihood
    public :: genlogistic_nnlf
    public :: gennorm_loglikelihood
    public :: gennorm_nnlf
    public :: nakagami_loglikelihood
    public :: nakagami_nnlf
    public :: powernorm_loglikelihood
    public :: powernorm_nnlf
    public :: loggamma_loglikelihood
    public :: loggamma_nnlf
    public :: wald_loglikelihood
    public :: wald_nnlf
    public :: gompertz_loglikelihood
    public :: gompertz_nnlf
    public :: invweibull_loglikelihood
    public :: invweibull_nnlf
    public :: betaprime_loglikelihood
    public :: betaprime_nnlf
    public :: burr12_loglikelihood
    public :: burr12_nnlf
    public :: genhalflogistic_loglikelihood
    public :: genhalflogistic_nnlf
    public :: exponpow_loglikelihood
    public :: exponpow_nnlf
    public :: exponweib_loglikelihood
    public :: exponweib_nnlf
    public :: powerlognorm_loglikelihood
    public :: powerlognorm_nnlf
    public :: levy_l_loglikelihood
    public :: levy_l_nnlf
    public :: weibull_max_loglikelihood
    public :: weibull_max_nnlf
    public :: rdist_loglikelihood
    public :: rdist_nnlf
    public :: skewcauchy_loglikelihood
    public :: skewcauchy_nnlf
    public :: dgamma_loglikelihood
    public :: dgamma_nnlf
    public :: laplace_asymmetric_loglikelihood
    public :: laplace_asymmetric_nnlf
    public :: truncnorm_loglikelihood
    public :: truncnorm_nnlf
    public :: loguniform_loglikelihood
    public :: loguniform_nnlf
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

    interface bernoulli_loglikelihood
        module procedure bernoulli_loglikelihood_real
        module procedure bernoulli_loglikelihood_int
    end interface bernoulli_loglikelihood

    interface bernoulli_nnlf
        module procedure bernoulli_nnlf_real
        module procedure bernoulli_nnlf_int
    end interface bernoulli_nnlf

    interface poisson_loglikelihood
        module procedure poisson_loglikelihood_real
        module procedure poisson_loglikelihood_int
    end interface poisson_loglikelihood

    interface poisson_nnlf
        module procedure poisson_nnlf_real
        module procedure poisson_nnlf_int
    end interface poisson_nnlf

    interface geometric_loglikelihood
        module procedure geometric_loglikelihood_real
        module procedure geometric_loglikelihood_int
    end interface geometric_loglikelihood

    interface geometric_nnlf
        module procedure geometric_nnlf_real
        module procedure geometric_nnlf_int
    end interface geometric_nnlf

    interface binomial_loglikelihood
        module procedure binomial_loglikelihood_real
        module procedure binomial_loglikelihood_int
    end interface binomial_loglikelihood

    interface binomial_nnlf
        module procedure binomial_nnlf_real
        module procedure binomial_nnlf_int
    end interface binomial_nnlf

    interface negative_binomial_loglikelihood
        module procedure negative_binomial_loglikelihood_real
        module procedure negative_binomial_loglikelihood_int
    end interface negative_binomial_loglikelihood

    interface negative_binomial_nnlf
        module procedure negative_binomial_nnlf_real
        module procedure negative_binomial_nnlf_int
    end interface negative_binomial_nnlf

    interface randint_loglikelihood
        module procedure randint_loglikelihood_real
        module procedure randint_loglikelihood_int
    end interface randint_loglikelihood

    interface randint_nnlf
        module procedure randint_nnlf_real
        module procedure randint_nnlf_int
    end interface randint_nnlf

    interface planck_loglikelihood
        module procedure planck_loglikelihood_real
        module procedure planck_loglikelihood_int
    end interface planck_loglikelihood

    interface planck_nnlf
        module procedure planck_nnlf_real
        module procedure planck_nnlf_int
    end interface planck_nnlf

    interface dlaplace_loglikelihood
        module procedure dlaplace_loglikelihood_real
        module procedure dlaplace_loglikelihood_int
    end interface dlaplace_loglikelihood

    interface dlaplace_nnlf
        module procedure dlaplace_nnlf_real
        module procedure dlaplace_nnlf_int
    end interface dlaplace_nnlf

    interface logser_loglikelihood
        module procedure logser_loglikelihood_real
        module procedure logser_loglikelihood_int
    end interface logser_loglikelihood

    interface logser_nnlf
        module procedure logser_nnlf_real
        module procedure logser_nnlf_int
    end interface logser_nnlf

    interface betabinom_loglikelihood
        module procedure betabinom_loglikelihood_real
        module procedure betabinom_loglikelihood_int
    end interface betabinom_loglikelihood

    interface betabinom_nnlf
        module procedure betabinom_nnlf_real
        module procedure betabinom_nnlf_int
    end interface betabinom_nnlf

    interface hypergeom_loglikelihood
        module procedure hypergeom_loglikelihood_real
        module procedure hypergeom_loglikelihood_int
    end interface hypergeom_loglikelihood

    interface hypergeom_nnlf
        module procedure hypergeom_nnlf_real
        module procedure hypergeom_nnlf_int
    end interface hypergeom_nnlf

    interface nhypergeom_loglikelihood
        module procedure nhypergeom_loglikelihood_real
        module procedure nhypergeom_loglikelihood_int
    end interface nhypergeom_loglikelihood

    interface nhypergeom_nnlf
        module procedure nhypergeom_nnlf_real
        module procedure nhypergeom_nnlf_int
    end interface nhypergeom_nnlf

    interface boltzmann_loglikelihood
        module procedure boltzmann_loglikelihood_real
        module procedure boltzmann_loglikelihood_int
    end interface boltzmann_loglikelihood

    interface boltzmann_nnlf
        module procedure boltzmann_nnlf_real
        module procedure boltzmann_nnlf_int
    end interface boltzmann_nnlf

    interface betanbinom_loglikelihood
        module procedure betanbinom_loglikelihood_real
        module procedure betanbinom_loglikelihood_int
    end interface betanbinom_loglikelihood

    interface betanbinom_nnlf
        module procedure betanbinom_nnlf_real
        module procedure betanbinom_nnlf_int
    end interface betanbinom_nnlf

    interface yulesimon_loglikelihood
        module procedure yulesimon_loglikelihood_real
        module procedure yulesimon_loglikelihood_int
    end interface yulesimon_loglikelihood

    interface yulesimon_nnlf
        module procedure yulesimon_nnlf_real
        module procedure yulesimon_nnlf_int
    end interface yulesimon_nnlf

    interface zipf_loglikelihood
        module procedure zipf_loglikelihood_real
        module procedure zipf_loglikelihood_int
    end interface zipf_loglikelihood

    interface zipf_nnlf
        module procedure zipf_nnlf_real
        module procedure zipf_nnlf_int
    end interface zipf_nnlf

    interface zipfian_loglikelihood
        module procedure zipfian_loglikelihood_real
        module procedure zipfian_loglikelihood_int
    end interface zipfian_loglikelihood

    interface zipfian_nnlf
        module procedure zipfian_nnlf_real
        module procedure zipfian_nnlf_int
    end interface zipfian_nnlf

    public :: foldnorm_loglikelihood, foldnorm_nnlf
    public :: foldcauchy_loglikelihood, foldcauchy_nnlf
    public :: recipinvgauss_loglikelihood, recipinvgauss_nnlf
    public :: truncpareto_loglikelihood, truncpareto_nnlf
    public :: exponnorm_loglikelihood, exponnorm_nnlf
    public :: johnsonsb_loglikelihood, johnsonsb_nnlf
    public :: johnsonsu_loglikelihood, johnsonsu_nnlf
    public :: trapezoid_loglikelihood, trapezoid_nnlf
    public :: burr_loglikelihood, burr_nnlf
    public :: mielke_loglikelihood, mielke_nnlf
    public :: gibrat_loglikelihood, gibrat_nnlf
    public :: wrapcauchy_loglikelihood, wrapcauchy_nnlf
    public :: genextreme_loglikelihood, genextreme_nnlf
    public :: kappa3_loglikelihood, kappa3_nnlf
    public :: kappa4_loglikelihood, kappa4_nnlf
    public :: truncweibull_min_loglikelihood, truncweibull_min_nnlf
    public :: gengamma_loglikelihood, gengamma_nnlf
    public :: halfgennorm_loglikelihood, halfgennorm_nnlf
    public :: argus_loglikelihood, argus_nnlf
    public :: erlang_loglikelihood, erlang_nnlf
    public :: crystalball_loglikelihood, crystalball_nnlf
    public :: jf_skew_t_loglikelihood, jf_skew_t_nnlf
    public :: pearson3_loglikelihood, pearson3_nnlf
    public :: rel_breitwigner_loglikelihood, rel_breitwigner_nnlf
    public :: genexpon_loglikelihood, genexpon_nnlf
    public :: skewnorm_loglikelihood, skewnorm_nnlf
    public :: tukeylambda_loglikelihood, tukeylambda_nnlf
    public :: rice_loglikelihood, rice_nnlf
    public :: dpareto_lognorm_loglikelihood, dpareto_lognorm_nnlf
    public :: vonmises_loglikelihood, vonmises_nnlf
    public :: vonmises_line_loglikelihood, vonmises_line_nnlf
    public :: kstwobign_loglikelihood, kstwobign_nnlf
    public :: irwinhall_loglikelihood, irwinhall_nnlf
    public :: ksone_loglikelihood, ksone_nnlf
    public :: kstwo_loglikelihood, kstwo_nnlf
    public :: levy_stable_loglikelihood, levy_stable_nnlf
    public :: studentized_range_loglikelihood, studentized_range_nnlf
    public :: ncx2_loglikelihood, ncx2_nnlf
    public :: ncf_loglikelihood, ncf_nnlf
    public :: randint_loglikelihood, randint_nnlf
    public :: planck_loglikelihood, planck_nnlf
    public :: dlaplace_loglikelihood, dlaplace_nnlf
    public :: logser_loglikelihood, logser_nnlf
    public :: betabinom_loglikelihood, betabinom_nnlf
    public :: hypergeom_loglikelihood, hypergeom_nnlf
    public :: nhypergeom_loglikelihood, nhypergeom_nnlf
    public :: boltzmann_loglikelihood, boltzmann_nnlf
    public :: betanbinom_loglikelihood, betanbinom_nnlf
    public :: yulesimon_loglikelihood, yulesimon_nnlf
    public :: zipf_loglikelihood, zipf_nnlf
    public :: zipfian_loglikelihood, zipfian_nnlf
    public :: geninvgauss_loglikelihood, geninvgauss_nnlf
    public :: norminvgauss_loglikelihood, norminvgauss_nnlf
    public :: skellam_loglikelihood, skellam_nnlf
    public :: genhyperbolic_loglikelihood, genhyperbolic_nnlf
    public :: nchypergeom_fisher_loglikelihood, nchypergeom_fisher_nnlf
    public :: nct_loglikelihood, nct_nnlf
    public :: gausshyper_loglikelihood, gausshyper_nnlf
    public :: nchypergeom_wallenius_loglikelihood, nchypergeom_wallenius_nnlf
    public :: poisson_binom_loglikelihood, poisson_binom_nnlf

    interface skellam_loglikelihood
        module procedure skellam_loglikelihood_real
        module procedure skellam_loglikelihood_int
    end interface skellam_loglikelihood

    interface skellam_nnlf
        module procedure skellam_nnlf_real
        module procedure skellam_nnlf_int
    end interface skellam_nnlf

    interface nchypergeom_fisher_loglikelihood
        module procedure nchypergeom_fisher_loglikelihood_real
        module procedure nchypergeom_fisher_loglikelihood_int
    end interface nchypergeom_fisher_loglikelihood

    interface nchypergeom_fisher_nnlf
        module procedure nchypergeom_fisher_nnlf_real
        module procedure nchypergeom_fisher_nnlf_int
    end interface nchypergeom_fisher_nnlf

    interface nchypergeom_wallenius_loglikelihood
        module procedure nchypergeom_wallenius_loglikelihood_real
        module procedure nchypergeom_wallenius_loglikelihood_int
    end interface nchypergeom_wallenius_loglikelihood

    interface nchypergeom_wallenius_nnlf
        module procedure nchypergeom_wallenius_nnlf_real
        module procedure nchypergeom_wallenius_nnlf_int
    end interface nchypergeom_wallenius_nnlf

    interface poisson_binom_loglikelihood
        module procedure poisson_binom_loglikelihood_real
        module procedure poisson_binom_loglikelihood_int
    end interface poisson_binom_loglikelihood

    interface poisson_binom_nnlf
        module procedure poisson_binom_nnlf_real
        module procedure poisson_binom_nnlf_int
    end interface poisson_binom_nnlf

contains

    pure function normal_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(normal_logpdf(data, loc, scale))
    end function normal_loglikelihood

    pure function normal_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -normal_loglikelihood(data, loc, scale)
    end function normal_nnlf

    pure function uniform_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(uniform_logpdf(data, loc, scale))
    end function uniform_loglikelihood

    pure function uniform_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -uniform_loglikelihood(data, loc, scale)
    end function uniform_nnlf

    pure function exponential_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(exponential_logpdf(data, loc, scale))
    end function exponential_loglikelihood

    pure function exponential_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -exponential_loglikelihood(data, loc, scale)
    end function exponential_nnlf

    pure function laplace_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(laplace_logpdf(data, loc, scale))
    end function laplace_loglikelihood

    pure function laplace_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -laplace_loglikelihood(data, loc, scale)
    end function laplace_nnlf

    pure function logistic_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(logistic_logpdf(data, loc, scale))
    end function logistic_loglikelihood

    pure function logistic_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -logistic_loglikelihood(data, loc, scale)
    end function logistic_nnlf

    pure function cauchy_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(cauchy_logpdf(data, loc, scale))
    end function cauchy_loglikelihood

    pure function cauchy_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -cauchy_loglikelihood(data, loc, scale)
    end function cauchy_nnlf

    pure function rayleigh_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(rayleigh_logpdf(data, loc, scale))
    end function rayleigh_loglikelihood

    pure function rayleigh_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -rayleigh_loglikelihood(data, loc, scale)
    end function rayleigh_nnlf

    pure function gamma_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(gamma_logpdf(data, a, loc, scale))
    end function gamma_loglikelihood

    pure function gamma_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -gamma_loglikelihood(data, a, loc, scale)
    end function gamma_nnlf

    pure function chi2_loglikelihood(data, df, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(chi2_logpdf(data, df, loc, scale))
    end function chi2_loglikelihood

    pure function chi2_nnlf(data, df, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -chi2_loglikelihood(data, df, loc, scale)
    end function chi2_nnlf

    pure function t_loglikelihood(data, df, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(t_logpdf(data, df, loc, scale))
    end function t_loglikelihood

    pure function t_nnlf(data, df, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -t_loglikelihood(data, df, loc, scale)
    end function t_nnlf

    pure function lognormal_loglikelihood(data, s, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: s !! log-space shape, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(lognormal_logpdf(data, s, loc, scale))
    end function lognormal_loglikelihood

    pure function lognormal_nnlf(data, s, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: s !! log-space shape, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -lognormal_loglikelihood(data, s, loc, scale)
    end function lognormal_nnlf

    pure function weibull_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(weibull_logpdf(data, c, loc, scale))
    end function weibull_loglikelihood

    pure function weibull_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -weibull_loglikelihood(data, c, loc, scale)
    end function weibull_nnlf

    pure function pareto_loglikelihood(data, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(pareto_logpdf(data, b, loc, scale))
    end function pareto_loglikelihood

    pure function pareto_nnlf(data, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -pareto_loglikelihood(data, b, loc, scale)
    end function pareto_nnlf

    pure function beta_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! first shape parameter, > 0
        real(dp), intent(in) :: b !! second shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(beta_logpdf(data, a, b, loc, scale))
    end function beta_loglikelihood

    pure function beta_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! first shape parameter, > 0
        real(dp), intent(in) :: b !! second shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -beta_loglikelihood(data, a, b, loc, scale)
    end function beta_nnlf

    pure function f_loglikelihood(data, dfn, dfd, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(f_logpdf(data, dfn, dfd, loc, scale))
    end function f_loglikelihood

    pure function f_nnlf(data, dfn, dfd, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -f_loglikelihood(data, dfn, dfd, loc, scale)
    end function f_nnlf

    pure function gumbel_r_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(gumbel_r_logpdf(data, loc, scale))
    end function gumbel_r_loglikelihood

    pure function gumbel_r_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -gumbel_r_loglikelihood(data, loc, scale)
    end function gumbel_r_nnlf

    pure function gumbel_l_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(gumbel_l_logpdf(data, loc, scale))
    end function gumbel_l_loglikelihood

    pure function gumbel_l_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -gumbel_l_loglikelihood(data, loc, scale)
    end function gumbel_l_nnlf

    pure function powerlaw_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: loglike

        loglike = sum(powerlaw_logpdf(data, a, loc, scale))
    end function powerlaw_loglikelihood

    pure function powerlaw_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: value

        value = -powerlaw_loglikelihood(data, a, loc, scale)
    end function powerlaw_nnlf

    pure function triang_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: loglike

        loglike = sum(triang_logpdf(data, c, loc, scale))
    end function triang_loglikelihood

    pure function triang_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: value

        value = -triang_loglikelihood(data, c, loc, scale)
    end function triang_nnlf

    pure function genpareto_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(genpareto_logpdf(data, c, loc, scale))
    end function genpareto_loglikelihood

    pure function genpareto_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -genpareto_loglikelihood(data, c, loc, scale)
    end function genpareto_nnlf

    pure function arcsine_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: loglike

        loglike = sum(arcsine_logpdf(data, loc, scale))
    end function arcsine_loglikelihood

    pure function arcsine_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: value

        value = -arcsine_loglikelihood(data, loc, scale)
    end function arcsine_nnlf

    pure function halfnorm_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(halfnorm_logpdf(data, loc, scale))
    end function halfnorm_loglikelihood

    pure function halfnorm_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -halfnorm_loglikelihood(data, loc, scale)
    end function halfnorm_nnlf

    pure function halfcauchy_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(halfcauchy_logpdf(data, loc, scale))
    end function halfcauchy_loglikelihood

    pure function halfcauchy_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -halfcauchy_loglikelihood(data, loc, scale)
    end function halfcauchy_nnlf

    pure function lomax_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(lomax_logpdf(data, c, loc, scale))
    end function lomax_loglikelihood

    pure function lomax_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -lomax_loglikelihood(data, c, loc, scale)
    end function lomax_nnlf

    pure function chi_loglikelihood(data, df, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(chi_logpdf(data, df, loc, scale))
    end function chi_loglikelihood

    pure function chi_nnlf(data, df, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -chi_loglikelihood(data, df, loc, scale)
    end function chi_nnlf

    pure function maxwell_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(maxwell_logpdf(data, loc, scale))
    end function maxwell_loglikelihood

    pure function maxwell_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -maxwell_loglikelihood(data, loc, scale)
    end function maxwell_nnlf

    pure function cosine_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(cosine_logpdf(data, loc, scale))
    end function cosine_loglikelihood

    pure function cosine_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -cosine_loglikelihood(data, loc, scale)
    end function cosine_nnlf

    pure function semicircular_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike

        loglike = sum(semicircular_logpdf(data, loc, scale))
    end function semicircular_loglikelihood

    pure function semicircular_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value

        value = -semicircular_loglikelihood(data, loc, scale)
    end function semicircular_nnlf

    pure function bernoulli_loglikelihood_real(data, p, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(bernoulli_logpmf(data, p, loc))
    end function bernoulli_loglikelihood_real

    pure function bernoulli_nnlf_real(data, p, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -bernoulli_loglikelihood_real(data, p, loc)
    end function bernoulli_nnlf_real

    pure function bernoulli_loglikelihood_int(data, p, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(bernoulli_logpmf(data, p, loc))
    end function bernoulli_loglikelihood_int

    pure function bernoulli_nnlf_int(data, p, loc) result(value)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -bernoulli_loglikelihood_int(data, p, loc)
    end function bernoulli_nnlf_int

    pure function poisson_loglikelihood_real(data, mu, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(poisson_logpmf(data, mu, loc))
    end function poisson_loglikelihood_real

    pure function poisson_nnlf_real(data, mu, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -poisson_loglikelihood_real(data, mu, loc)
    end function poisson_nnlf_real

    pure function poisson_loglikelihood_int(data, mu, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(poisson_logpmf(data, mu, loc))
    end function poisson_loglikelihood_int

    pure function poisson_nnlf_int(data, mu, loc) result(value)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: mu !! mean, >= 0
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -poisson_loglikelihood_int(data, mu, loc)
    end function poisson_nnlf_int

    pure function geometric_loglikelihood_real(data, p, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(geometric_logpmf(data, p, loc))
    end function geometric_loglikelihood_real

    pure function geometric_nnlf_real(data, p, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -geometric_loglikelihood_real(data, p, loc)
    end function geometric_nnlf_real

    pure function geometric_loglikelihood_int(data, p, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(geometric_logpmf(data, p, loc))
    end function geometric_loglikelihood_int

    pure function geometric_nnlf_int(data, p, loc) result(value)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -geometric_loglikelihood_int(data, p, loc)
    end function geometric_nnlf_int

    pure function binomial_loglikelihood_real(data, n, p, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer value
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(binomial_logpmf(data, n, p, loc))
    end function binomial_loglikelihood_real

    pure function binomial_nnlf_real(data, n, p, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer value
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -binomial_loglikelihood_real(data, n, p, loc)
    end function binomial_nnlf_real

    pure function binomial_loglikelihood_int(data, n, p, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent count observations
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: loglike

        loglike = sum(binomial_logpmf(data, n, p, loc))
    end function binomial_loglikelihood_int

    pure function binomial_nnlf_int(data, n, p, loc) result(value)
        integer, intent(in) :: data(:) !! independent count observations
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability in [0, 1]
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: value

        value = -binomial_loglikelihood_int(data, n, p, loc)
    end function binomial_nnlf_int

    pure function negative_binomial_loglikelihood_real(data, n, p, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike

        loglike = sum(negative_binomial_logpmf(data, n, p, loc))
    end function negative_binomial_loglikelihood_real

    pure function negative_binomial_nnlf_real(data, n, p, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value

        value = -negative_binomial_loglikelihood_real(data, n, p, loc)
    end function negative_binomial_nnlf_real

    pure function negative_binomial_loglikelihood_int(data, n, p, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike

        loglike = sum(negative_binomial_logpmf(data, n, p, loc))
    end function negative_binomial_loglikelihood_int

    pure function negative_binomial_nnlf_int(data, n, p, loc) result(value)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! required successes, finite and > 0
        real(dp), intent(in) :: p !! success probability in (0, 1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value

        value = -negative_binomial_loglikelihood_int(data, n, p, loc)
    end function negative_binomial_nnlf_int

    pure function anglit_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(anglit_logpdf(data, loc, scale))
    end function anglit_loglikelihood

    pure function anglit_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -anglit_loglikelihood(data, loc, scale)
    end function anglit_nnlf

    pure function moyal_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(moyal_logpdf(data, loc, scale))
    end function moyal_loglikelihood

    pure function moyal_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -moyal_loglikelihood(data, loc, scale)
    end function moyal_nnlf

    pure function landau_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(landau_logpdf(data, loc, scale))
    end function landau_loglikelihood

    pure function landau_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -landau_loglikelihood(data, loc, scale)
    end function landau_nnlf


    pure function hypsecant_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(hypsecant_logpdf(data, loc, scale))
    end function hypsecant_loglikelihood

    pure function hypsecant_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -hypsecant_loglikelihood(data, loc, scale)
    end function hypsecant_nnlf

    pure function halflogistic_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(halflogistic_logpdf(data, loc, scale))
    end function halflogistic_loglikelihood

    pure function halflogistic_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -halflogistic_loglikelihood(data, loc, scale)
    end function halflogistic_nnlf

    pure function invgamma_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive inverse-gamma shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(invgamma_logpdf(data, a, loc, scale))
    end function invgamma_loglikelihood

    pure function invgamma_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive inverse-gamma shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -invgamma_loglikelihood(data, a, loc, scale)
    end function invgamma_nnlf

    pure function invgauss_loglikelihood(data, mu_shape, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(invgauss_logpdf(data, mu_shape, loc, scale))
    end function invgauss_loglikelihood

    pure function invgauss_nnlf(data, mu_shape, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -invgauss_loglikelihood(data, mu_shape, loc, scale)
    end function invgauss_nnlf

    pure function levy_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(levy_logpdf(data, loc, scale))
    end function levy_loglikelihood

    pure function levy_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -levy_loglikelihood(data, loc, scale)
    end function levy_nnlf

    pure function loglaplace_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive log-Laplace shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(loglaplace_logpdf(data, c, loc, scale))
    end function loglaplace_loglikelihood

    pure function loglaplace_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive log-Laplace shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -loglaplace_loglikelihood(data, c, loc, scale)
    end function loglaplace_nnlf


    pure function bradford_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Bradford shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: loglike

        loglike = sum(bradford_logpdf(data, c, loc, scale))
    end function bradford_loglikelihood

    pure function bradford_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Bradford shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: value

        value = -bradford_loglikelihood(data, c, loc, scale)
    end function bradford_nnlf

    pure function truncexpon_loglikelihood(data, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! positive standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(truncexpon_logpdf(data, b, loc, scale))
    end function truncexpon_loglikelihood

    pure function truncexpon_nnlf(data, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! positive standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -truncexpon_loglikelihood(data, b, loc, scale)
    end function truncexpon_nnlf

    pure function fisk_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Fisk shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(fisk_logpdf(data, c, loc, scale))
    end function fisk_loglikelihood

    pure function fisk_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Fisk shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -fisk_loglikelihood(data, c, loc, scale)
    end function fisk_nnlf

    pure function dweibull_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive double-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(dweibull_logpdf(data, c, loc, scale))
    end function dweibull_loglikelihood

    pure function dweibull_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive double-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -dweibull_loglikelihood(data, c, loc, scale)
    end function dweibull_nnlf

    pure function alpha_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive alpha shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(alpha_logpdf(data, a, loc, scale))
    end function alpha_loglikelihood

    pure function alpha_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive alpha shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -alpha_loglikelihood(data, a, loc, scale)
    end function alpha_nnlf

    pure function fatiguelife_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive fatigue-life shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(fatiguelife_logpdf(data, c, loc, scale))
    end function fatiguelife_loglikelihood

    pure function fatiguelife_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive fatigue-life shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -fatiguelife_loglikelihood(data, c, loc, scale)
    end function fatiguelife_nnlf

    pure function genlogistic_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive generalized-logistic shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(genlogistic_logpdf(data, c, loc, scale))
    end function genlogistic_loglikelihood

    pure function genlogistic_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive generalized-logistic shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -genlogistic_loglikelihood(data, c, loc, scale)
    end function genlogistic_nnlf

    pure function gennorm_loglikelihood(data, beta, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: beta !! positive generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(gennorm_logpdf(data, beta, loc, scale))
    end function gennorm_loglikelihood

    pure function gennorm_nnlf(data, beta, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: beta !! positive generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -gennorm_loglikelihood(data, beta, loc, scale)
    end function gennorm_nnlf

    pure function nakagami_loglikelihood(data, nu, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(nakagami_logpdf(data, nu, loc, scale))
    end function nakagami_loglikelihood

    pure function nakagami_nnlf(data, nu, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -nakagami_loglikelihood(data, nu, loc, scale)
    end function nakagami_nnlf

    pure function powernorm_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive power-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(powernorm_logpdf(data, c, loc, scale))
    end function powernorm_loglikelihood

    pure function powernorm_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive power-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -powernorm_loglikelihood(data, c, loc, scale)
    end function powernorm_nnlf

    pure function loggamma_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(loggamma_logpdf(data, c, loc, scale))
    end function loggamma_loglikelihood

    pure function loggamma_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -loggamma_loglikelihood(data, c, loc, scale)
    end function loggamma_nnlf

    pure function wald_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(wald_logpdf(data, loc, scale))
    end function wald_loglikelihood

    pure function wald_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -wald_loglikelihood(data, loc, scale)
    end function wald_nnlf


    pure function gompertz_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Gompertz shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(gompertz_logpdf(data, c, loc, scale))
    end function gompertz_loglikelihood

    pure function gompertz_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Gompertz shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -gompertz_loglikelihood(data, c, loc, scale)
    end function gompertz_nnlf

    pure function invweibull_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(invweibull_logpdf(data, c, loc, scale))
    end function invweibull_loglikelihood

    pure function invweibull_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -invweibull_loglikelihood(data, c, loc, scale)
    end function invweibull_nnlf

    pure function betaprime_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive first beta-prime shape parameter
        real(dp), intent(in) :: b !! positive second beta-prime shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(betaprime_logpdf(data, a, b, loc, scale))
    end function betaprime_loglikelihood

    pure function betaprime_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive first beta-prime shape parameter
        real(dp), intent(in) :: b !! positive second beta-prime shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -betaprime_loglikelihood(data, a, b, loc, scale)
    end function betaprime_nnlf

    pure function burr12_loglikelihood(data, c, d, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive first Burr-XII shape parameter
        real(dp), intent(in) :: d !! positive second Burr-XII shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(burr12_logpdf(data, c, d, loc, scale))
    end function burr12_loglikelihood

    pure function burr12_nnlf(data, c, d, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive first Burr-XII shape parameter
        real(dp), intent(in) :: d !! positive second Burr-XII shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -burr12_loglikelihood(data, c, d, loc, scale)
    end function burr12_nnlf

    pure function genhalflogistic_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive generalized half-logistic shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(genhalflogistic_logpdf(data, c, loc, scale))
    end function genhalflogistic_loglikelihood

    pure function genhalflogistic_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive generalized half-logistic shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -genhalflogistic_loglikelihood(data, c, loc, scale)
    end function genhalflogistic_nnlf

    pure function exponpow_loglikelihood(data, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! positive exponential-power shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(exponpow_logpdf(data, b, loc, scale))
    end function exponpow_loglikelihood

    pure function exponpow_nnlf(data, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! positive exponential-power shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -exponpow_loglikelihood(data, b, loc, scale)
    end function exponpow_nnlf

    pure function exponweib_loglikelihood(data, a, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive exponentiation shape parameter
        real(dp), intent(in) :: c !! positive Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(exponweib_logpdf(data, a, c, loc, scale))
    end function exponweib_loglikelihood

    pure function exponweib_nnlf(data, a, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive exponentiation shape parameter
        real(dp), intent(in) :: c !! positive Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -exponweib_loglikelihood(data, a, c, loc, scale)
    end function exponweib_nnlf

    pure function powerlognorm_loglikelihood(data, c, s, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive power shape parameter
        real(dp), intent(in) :: s !! positive lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(powerlognorm_logpdf(data, c, s, loc, scale))
    end function powerlognorm_loglikelihood

    pure function powerlognorm_nnlf(data, c, s, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive power shape parameter
        real(dp), intent(in) :: s !! positive lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -powerlognorm_loglikelihood(data, c, s, loc, scale)
    end function powerlognorm_nnlf

    pure function levy_l_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(levy_l_logpdf(data, loc, scale))
    end function levy_l_loglikelihood

    pure function levy_l_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -levy_l_loglikelihood(data, loc, scale)
    end function levy_l_nnlf

    pure function weibull_max_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Weibull-maximum shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(weibull_max_logpdf(data, c, loc, scale))
    end function weibull_max_loglikelihood

    pure function weibull_max_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive Weibull-maximum shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -weibull_max_loglikelihood(data, c, loc, scale)
    end function weibull_max_nnlf

    pure function rdist_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive R-distribution shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: loglike
        loglike = sum(rdist_logpdf(data, c, loc, scale))
    end function rdist_loglikelihood

    pure function rdist_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive R-distribution shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: value
        value = -rdist_loglikelihood(data, c, loc, scale)
    end function rdist_nnlf

    pure function skewcauchy_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(skewcauchy_logpdf(data, a, loc, scale))
    end function skewcauchy_loglikelihood

    pure function skewcauchy_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -skewcauchy_loglikelihood(data, a, loc, scale)
    end function skewcauchy_nnlf


    pure function dgamma_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(dgamma_logpdf(data, a, loc, scale))
    end function dgamma_loglikelihood

    pure function dgamma_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -dgamma_loglikelihood(data, a, loc, scale)
    end function dgamma_nnlf

    pure function laplace_asymmetric_loglikelihood(data, kappa, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(laplace_asymmetric_logpdf(data, kappa, loc, scale))
    end function laplace_asymmetric_loglikelihood

    pure function laplace_asymmetric_nnlf(data, kappa, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -laplace_asymmetric_loglikelihood(data, kappa, loc, scale)
    end function laplace_asymmetric_nnlf

    pure function truncnorm_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: loglike
        loglike = sum(truncnorm_logpdf(data, a, b, loc, scale))
    end function truncnorm_loglikelihood

    pure function truncnorm_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: value
        value = -truncnorm_loglikelihood(data, a, b, loc, scale)
    end function truncnorm_nnlf

    pure function loguniform_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(loguniform_logpdf(data, a, b, loc, scale))
    end function loguniform_loglikelihood

    pure function loguniform_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -loguniform_loglikelihood(data, a, b, loc, scale)
    end function loguniform_nnlf

    pure function foldnorm_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(foldnorm_logpdf(data, c, loc, scale))
    end function foldnorm_loglikelihood

    pure function foldnorm_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -foldnorm_loglikelihood(data, c, loc, scale)
    end function foldnorm_nnlf

    pure function foldcauchy_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(foldcauchy_logpdf(data, c, loc, scale))
    end function foldcauchy_loglikelihood

    pure function foldcauchy_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -foldcauchy_loglikelihood(data, c, loc, scale)
    end function foldcauchy_nnlf

    pure function recipinvgauss_loglikelihood(data, mu_shape, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: mu_shape !! positive reciprocal inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(recipinvgauss_logpdf(data, mu_shape, loc, scale))
    end function recipinvgauss_loglikelihood

    pure function recipinvgauss_nnlf(data, mu_shape, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: mu_shape !! positive reciprocal inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -recipinvgauss_loglikelihood(data, mu_shape, loc, scale)
    end function recipinvgauss_nnlf

    pure function truncpareto_loglikelihood(data, b, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(truncpareto_logpdf(data, b, c, loc, scale))
    end function truncpareto_loglikelihood

    pure function truncpareto_nnlf(data, b, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -truncpareto_loglikelihood(data, b, c, loc, scale)
    end function truncpareto_nnlf

    pure function exponnorm_loglikelihood(data, k, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: k !! positive exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: loglike
        loglike = sum(exponnorm_logpdf(data, k, loc, scale))
    end function exponnorm_loglikelihood

    pure function exponnorm_nnlf(data, k, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: k !! positive exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: value
        value = -exponnorm_loglikelihood(data, k, loc, scale)
    end function exponnorm_nnlf

    pure function johnsonsb_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: loglike
        loglike = sum(johnsonsb_logpdf(data, a, b, loc, scale))
    end function johnsonsb_loglikelihood

    pure function johnsonsb_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: value
        value = -johnsonsb_loglikelihood(data, a, b, loc, scale)
    end function johnsonsb_nnlf

    pure function johnsonsu_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(johnsonsu_logpdf(data, a, b, loc, scale))
    end function johnsonsu_loglikelihood

    pure function johnsonsu_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -johnsonsu_loglikelihood(data, a, b, loc, scale)
    end function johnsonsu_nnlf

    pure function trapezoid_loglikelihood(data, c, d, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! standardized left plateau edge in [0,1]
        real(dp), intent(in) :: d !! standardized right plateau edge in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: loglike
        loglike = sum(trapezoid_logpdf(data, c, d, loc, scale))
    end function trapezoid_loglikelihood

    pure function trapezoid_nnlf(data, c, d, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! standardized left plateau edge in [0,1]
        real(dp), intent(in) :: d !! standardized right plateau edge in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: value
        value = -trapezoid_loglikelihood(data, c, d, loc, scale)
    end function trapezoid_nnlf

    pure function burr_loglikelihood(data, c, d, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(burr_logpdf(data, c, d, loc, scale))
    end function burr_loglikelihood

    pure function burr_nnlf(data, c, d, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -burr_loglikelihood(data, c, d, loc, scale)
    end function burr_nnlf

    pure function mielke_loglikelihood(data, k, s, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: k !! positive finite beta-kappa first shape parameter
        real(dp), intent(in) :: s !! positive finite beta-kappa second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(mielke_logpdf(data, k, s, loc, scale))
    end function mielke_loglikelihood

    pure function mielke_nnlf(data, k, s, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: k !! positive finite beta-kappa first shape parameter
        real(dp), intent(in) :: s !! positive finite beta-kappa second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -mielke_loglikelihood(data, k, s, loc, scale)
    end function mielke_nnlf

    pure function gibrat_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(gibrat_logpdf(data, loc, scale))
    end function gibrat_loglikelihood

    pure function gibrat_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -gibrat_loglikelihood(data, loc, scale)
    end function gibrat_nnlf

    pure function wrapcauchy_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: loglike
        loglike = sum(wrapcauchy_logpdf(data, c, loc, scale))
    end function wrapcauchy_loglikelihood

    pure function wrapcauchy_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: value
        value = -wrapcauchy_loglikelihood(data, c, loc, scale)
    end function wrapcauchy_nnlf

    pure function genextreme_loglikelihood(data, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(genextreme_logpdf(data, c, loc, scale))
    end function genextreme_loglikelihood

    pure function genextreme_nnlf(data, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -genextreme_loglikelihood(data, c, loc, scale)
    end function genextreme_nnlf

    pure function kappa3_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(kappa3_logpdf(data, a, loc, scale))
    end function kappa3_loglikelihood

    pure function kappa3_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -kappa3_loglikelihood(data, a, loc, scale)
    end function kappa3_nnlf

    pure function kappa4_loglikelihood(data, h, k, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(kappa4_logpdf(data, h, k, loc, scale))
    end function kappa4_loglikelihood

    pure function kappa4_nnlf(data, h, k, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -kappa4_loglikelihood(data, h, k, loc, scale)
    end function kappa4_nnlf

    pure function truncweibull_min_loglikelihood(data, c, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(truncweibull_min_logpdf(data, c, a, b, loc, scale))
    end function truncweibull_min_loglikelihood

    pure function truncweibull_min_nnlf(data, c, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -truncweibull_min_loglikelihood(data, c, a, b, loc, scale)
    end function truncweibull_min_nnlf

    pure function gengamma_loglikelihood(data, a, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(gengamma_logpdf(data, a, c, loc, scale))
    end function gengamma_loglikelihood

    pure function gengamma_nnlf(data, a, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -gengamma_loglikelihood(data, a, c, loc, scale)
    end function gengamma_nnlf

    pure function halfgennorm_loglikelihood(data, beta, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(halfgennorm_logpdf(data, beta, loc, scale))
    end function halfgennorm_loglikelihood

    pure function halfgennorm_nnlf(data, beta, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -halfgennorm_loglikelihood(data, beta, loc, scale)
    end function halfgennorm_nnlf

    pure function argus_loglikelihood(data, chi, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: loglike
        loglike = sum(argus_logpdf(data, chi, loc, scale))
    end function argus_loglikelihood

    pure function argus_nnlf(data, chi, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: value
        value = -argus_loglikelihood(data, chi, loc, scale)
    end function argus_nnlf

    pure function erlang_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(erlang_logpdf(data, a, loc, scale))
    end function erlang_loglikelihood

    pure function erlang_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -erlang_loglikelihood(data, a, loc, scale)
    end function erlang_nnlf

    pure function crystalball_loglikelihood(data, beta, m, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(crystalball_logpdf(data, beta, m, loc, scale))
    end function crystalball_loglikelihood

    pure function crystalball_nnlf(data, beta, m, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -crystalball_loglikelihood(data, beta, m, loc, scale)
    end function crystalball_nnlf

    pure function jf_skew_t_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(jf_skew_t_logpdf(data, a, b, loc, scale))
    end function jf_skew_t_loglikelihood

    pure function jf_skew_t_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -jf_skew_t_loglikelihood(data, a, b, loc, scale)
    end function jf_skew_t_nnlf

    pure function pearson3_loglikelihood(data, skew, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(pearson3_logpdf(data, skew, loc, scale))
    end function pearson3_loglikelihood

    pure function pearson3_nnlf(data, skew, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -pearson3_loglikelihood(data, skew, loc, scale)
    end function pearson3_nnlf

    pure function rel_breitwigner_loglikelihood(data, rho, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(rel_breitwigner_logpdf(data, rho, loc, scale))
    end function rel_breitwigner_loglikelihood

    pure function rel_breitwigner_nnlf(data, rho, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -rel_breitwigner_loglikelihood(data, rho, loc, scale)
    end function rel_breitwigner_nnlf

    pure function genexpon_loglikelihood(data, a, b, c, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(genexpon_logpdf(data, a, b, c, loc, scale))
    end function genexpon_loglikelihood

    pure function genexpon_nnlf(data, a, b, c, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -genexpon_loglikelihood(data, a, b, c, loc, scale)
    end function genexpon_nnlf

    pure function skewnorm_loglikelihood(data, a, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(skewnorm_logpdf(data, a, loc, scale))
    end function skewnorm_loglikelihood

    pure function skewnorm_nnlf(data, a, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: a !! finite skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -skewnorm_loglikelihood(data, a, loc, scale)
    end function skewnorm_nnlf

    pure function tukeylambda_loglikelihood(data, lam, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(tukeylambda_logpdf(data, lam, loc, scale))
    end function tukeylambda_loglikelihood

    pure function tukeylambda_nnlf(data, lam, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -tukeylambda_loglikelihood(data, lam, loc, scale)
    end function tukeylambda_nnlf

    pure function rice_loglikelihood(data, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(rice_logpdf(data, b, loc, scale))
    end function rice_loglikelihood

    pure function rice_nnlf(data, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -rice_loglikelihood(data, b, loc, scale)
    end function rice_nnlf


    pure function dpareto_lognorm_loglikelihood(data, u, s, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: loglike
        loglike = sum(dpareto_lognorm_logpdf(data, u, s, a, b, loc, scale))
    end function dpareto_lognorm_loglikelihood

    pure function dpareto_lognorm_nnlf(data, u, s, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent real-valued observations
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: value
        value = -dpareto_lognorm_loglikelihood(data, u, s, a, b, loc, scale)
    end function dpareto_lognorm_nnlf

    pure function vonmises_loglikelihood(data, kappa, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent angular observations
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: loglike
        loglike = sum(vonmises_logpdf(data, kappa, loc, scale))
    end function vonmises_loglikelihood

    pure function vonmises_nnlf(data, kappa, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent angular observations
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: value
        value = -vonmises_loglikelihood(data, kappa, loc, scale)
    end function vonmises_nnlf

    pure function vonmises_line_loglikelihood(data, kappa, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent bounded angular observations
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! support center (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: loglike
        loglike = sum(vonmises_line_logpdf(data, kappa, loc, scale))
    end function vonmises_line_loglikelihood

    pure function vonmises_line_nnlf(data, kappa, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent bounded angular observations
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! support center (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: value
        value = -vonmises_line_loglikelihood(data, kappa, loc, scale)
    end function vonmises_line_nnlf

    pure function kstwobign_loglikelihood(data, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent positive observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(kstwobign_logpdf(data, loc, scale))
    end function kstwobign_loglikelihood

    pure function kstwobign_nnlf(data, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent positive observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -kstwobign_loglikelihood(data, loc, scale)
    end function kstwobign_nnlf


    pure function irwinhall_loglikelihood(data, n, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: loglike
        loglike = sum(irwinhall_logpdf(data, n, loc, scale))
    end function irwinhall_loglikelihood

    pure function irwinhall_nnlf(data, n, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: value
        value = -irwinhall_loglikelihood(data, n, loc, scale)
    end function irwinhall_nnlf

    pure function ksone_loglikelihood(data, n, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: loglike
        loglike = sum(ksone_logpdf(data, n, loc, scale))
    end function ksone_loglikelihood

    pure function ksone_nnlf(data, n, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: value
        value = -ksone_loglikelihood(data, n, loc, scale)
    end function ksone_nnlf

    pure function kstwo_loglikelihood(data, n, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent finite-sample two-sided KS observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: loglike
        loglike = sum(kstwo_logpdf(data, n, loc, scale))
    end function kstwo_loglikelihood

    pure function kstwo_nnlf(data, n, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent finite-sample two-sided KS observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: value
        value = -kstwo_loglikelihood(data, n, loc, scale)
    end function kstwo_nnlf

    pure function ncx2_loglikelihood(data, df, nc, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(ncx2_logpdf(data, df, nc, loc, scale))
    end function ncx2_loglikelihood

    pure function ncx2_nnlf(data, df, nc, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -ncx2_loglikelihood(data, df, nc, loc, scale)
    end function ncx2_nnlf

    pure function ncf_loglikelihood(data, dfn, dfd, nc, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(ncf_logpdf(data, dfn, dfd, nc, loc, scale))
    end function ncf_loglikelihood

    pure function ncf_nnlf(data, dfn, dfd, nc, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -ncf_loglikelihood(data, dfn, dfd, nc, loc, scale)
    end function ncf_nnlf


    pure function randint_loglikelihood_real(data, low, high, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(randint_logpmf(data, low, high, loc))
    end function randint_loglikelihood_real

    pure function randint_nnlf_real(data, low, high, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -randint_loglikelihood_real(data, low, high, loc)
    end function randint_nnlf_real

    pure function randint_loglikelihood_int(data, low, high, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: low !! integer lower shape bound, inclusive
        integer, intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(randint_logpmf(data, low, high, loc))
    end function randint_loglikelihood_int

    pure function randint_nnlf_int(data, low, high, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: low !! integer lower shape bound, inclusive
        integer, intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -randint_loglikelihood_int(data, low, high, loc)
    end function randint_nnlf_int

    pure function planck_loglikelihood_real(data, lambda, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(planck_logpmf(data, lambda, loc))
    end function planck_loglikelihood_real

    pure function planck_nnlf_real(data, lambda, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -planck_loglikelihood_real(data, lambda, loc)
    end function planck_nnlf_real

    pure function planck_loglikelihood_int(data, lambda, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(planck_logpmf(data, lambda, loc))
    end function planck_loglikelihood_int

    pure function planck_nnlf_int(data, lambda, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -planck_loglikelihood_int(data, lambda, loc)
    end function planck_nnlf_int

    pure function dlaplace_loglikelihood_real(data, a, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(dlaplace_logpmf(data, a, loc))
    end function dlaplace_loglikelihood_real

    pure function dlaplace_nnlf_real(data, a, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -dlaplace_loglikelihood_real(data, a, loc)
    end function dlaplace_nnlf_real

    pure function dlaplace_loglikelihood_int(data, a, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(dlaplace_logpmf(data, a, loc))
    end function dlaplace_loglikelihood_int

    pure function dlaplace_nnlf_int(data, a, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -dlaplace_loglikelihood_int(data, a, loc)
    end function dlaplace_nnlf_int

    pure function logser_loglikelihood_real(data, p, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(logser_logpmf(data, p, loc))
    end function logser_loglikelihood_real

    pure function logser_nnlf_real(data, p, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -logser_loglikelihood_real(data, p, loc)
    end function logser_nnlf_real

    pure function logser_loglikelihood_int(data, p, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(logser_logpmf(data, p, loc))
    end function logser_loglikelihood_int

    pure function logser_nnlf_int(data, p, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -logser_loglikelihood_int(data, p, loc)
    end function logser_nnlf_int


    pure function betabinom_loglikelihood_real(data, n, a, b, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(betabinom_logpmf(data, n, a, b, loc))
    end function betabinom_loglikelihood_real

    pure function betabinom_nnlf_real(data, n, a, b, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -betabinom_loglikelihood_real(data, n, a, b, loc)
    end function betabinom_nnlf_real

    pure function betabinom_loglikelihood_int(data, n, a, b, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(betabinom_logpmf(data, n, a, b, loc))
    end function betabinom_loglikelihood_int

    pure function betabinom_nnlf_int(data, n, a, b, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -betabinom_loglikelihood_int(data, n, a, b, loc)
    end function betabinom_nnlf_int

    pure function hypergeom_loglikelihood_real(data, m, n, draws, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer Type-I count
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(hypergeom_logpmf(data, m, n, draws, loc))
    end function hypergeom_loglikelihood_real

    pure function hypergeom_nnlf_real(data, m, n, draws, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer Type-I count
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -hypergeom_loglikelihood_real(data, m, n, draws, loc)
    end function hypergeom_nnlf_real

    pure function hypergeom_loglikelihood_int(data, m, n, draws, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! Type-I count
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(hypergeom_logpmf(data, m, n, draws, loc))
    end function hypergeom_loglikelihood_int

    pure function hypergeom_nnlf_int(data, m, n, draws, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! Type-I count
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -hypergeom_loglikelihood_int(data, m, n, draws, loc)
    end function hypergeom_nnlf_int

    pure function nhypergeom_loglikelihood_real(data, m, n, r, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer success count
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(nhypergeom_logpmf(data, m, n, r, loc))
    end function nhypergeom_loglikelihood_real

    pure function nhypergeom_nnlf_real(data, m, n, r, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer success count
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -nhypergeom_loglikelihood_real(data, m, n, r, loc)
    end function nhypergeom_nnlf_real

    pure function nhypergeom_loglikelihood_int(data, m, n, r, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! success count
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(nhypergeom_logpmf(data, m, n, r, loc))
    end function nhypergeom_loglikelihood_int

    pure function nhypergeom_nnlf_int(data, m, n, r, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! success count
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -nhypergeom_loglikelihood_int(data, m, n, r, loc)
    end function nhypergeom_nnlf_int

    pure function boltzmann_loglikelihood_real(data, lambda, n, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(boltzmann_logpmf(data, lambda, n, loc))
    end function boltzmann_loglikelihood_real

    pure function boltzmann_nnlf_real(data, lambda, n, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -boltzmann_loglikelihood_real(data, lambda, n, loc)
    end function boltzmann_nnlf_real

    pure function boltzmann_loglikelihood_int(data, lambda, n, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(boltzmann_logpmf(data, lambda, n, loc))
    end function boltzmann_loglikelihood_int

    pure function boltzmann_nnlf_int(data, lambda, n, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -boltzmann_loglikelihood_int(data, lambda, n, loc)
    end function boltzmann_nnlf_int


    pure function betanbinom_loglikelihood_real(data, n, a, b, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(betanbinom_logpmf(data, n, a, b, loc))
    end function betanbinom_loglikelihood_real

    pure function betanbinom_nnlf_real(data, n, a, b, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -betanbinom_loglikelihood_real(data, n, a, b, loc)
    end function betanbinom_nnlf_real

    pure function betanbinom_loglikelihood_int(data, n, a, b, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(betanbinom_logpmf(data, n, a, b, loc))
    end function betanbinom_loglikelihood_int

    pure function betanbinom_nnlf_int(data, n, a, b, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -betanbinom_loglikelihood_int(data, n, a, b, loc)
    end function betanbinom_nnlf_int

    pure function yulesimon_loglikelihood_real(data, alpha, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(yulesimon_logpmf(data, alpha, loc))
    end function yulesimon_loglikelihood_real

    pure function yulesimon_nnlf_real(data, alpha, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -yulesimon_loglikelihood_real(data, alpha, loc)
    end function yulesimon_nnlf_real

    pure function yulesimon_loglikelihood_int(data, alpha, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(yulesimon_logpmf(data, alpha, loc))
    end function yulesimon_loglikelihood_int

    pure function yulesimon_nnlf_int(data, alpha, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -yulesimon_loglikelihood_int(data, alpha, loc)
    end function yulesimon_nnlf_int

    pure function zipf_loglikelihood_real(data, a, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(zipf_logpmf(data, a, loc))
    end function zipf_loglikelihood_real

    pure function zipf_nnlf_real(data, a, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -zipf_loglikelihood_real(data, a, loc)
    end function zipf_nnlf_real

    pure function zipf_loglikelihood_int(data, a, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(zipf_logpmf(data, a, loc))
    end function zipf_loglikelihood_int

    pure function zipf_nnlf_int(data, a, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -zipf_loglikelihood_int(data, a, loc)
    end function zipf_nnlf_int

    pure function zipfian_loglikelihood_real(data, a, n, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(zipfian_logpmf(data, a, n, loc))
    end function zipfian_loglikelihood_real

    pure function zipfian_nnlf_real(data, a, n, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -zipfian_loglikelihood_real(data, a, n, loc)
    end function zipfian_nnlf_real

    pure function zipfian_loglikelihood_int(data, a, n, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(zipfian_logpmf(data, a, n, loc))
    end function zipfian_loglikelihood_int

    pure function zipfian_nnlf_int(data, a, n, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -zipfian_loglikelihood_int(data, a, n, loc)
    end function zipfian_nnlf_int

    pure function geninvgauss_loglikelihood(data, p, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: p !! real generalized-inverse-Gaussian shape
        real(dp), intent(in) :: b !! strictly positive generalized-inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(geninvgauss_logpdf(data, p, b, loc, scale))
    end function geninvgauss_loglikelihood

    pure function geninvgauss_nnlf(data, p, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: p !! real generalized-inverse-Gaussian shape
        real(dp), intent(in) :: b !! strictly positive generalized-inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -geninvgauss_loglikelihood(data, p, b, loc, scale)
    end function geninvgauss_nnlf

    pure function norminvgauss_loglikelihood(data, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive normal-inverse-Gaussian tail shape
        real(dp), intent(in) :: b !! normal-inverse-Gaussian skew shape with abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(norminvgauss_logpdf(data, a, b, loc, scale))
    end function norminvgauss_loglikelihood

    pure function norminvgauss_nnlf(data, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive normal-inverse-Gaussian tail shape
        real(dp), intent(in) :: b !! normal-inverse-Gaussian skew shape with abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -norminvgauss_loglikelihood(data, a, b, loc, scale)
    end function norminvgauss_nnlf

    pure function skellam_loglikelihood_real(data, mu1, mu2, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: loglike
        loglike = sum(skellam_logpmf(data, mu1, mu2, loc))
    end function skellam_loglikelihood_real

    pure function skellam_nnlf_real(data, mu1, mu2, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: value
        value = -skellam_loglikelihood_real(data, mu1, mu2, loc)
    end function skellam_nnlf_real

    pure function skellam_loglikelihood_int(data, mu1, mu2, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: loglike
        loglike = sum(skellam_logpmf(data, mu1, mu2, loc))
    end function skellam_loglikelihood_int

    pure function skellam_nnlf_int(data, mu1, mu2, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: value
        value = -skellam_loglikelihood_int(data, mu1, mu2, loc)
    end function skellam_nnlf_int


    pure function genhyperbolic_loglikelihood(data, p, a, b, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: p !! real generalized-hyperbolic tail shape
        real(dp), intent(in) :: a !! positive generalized-hyperbolic shape
        real(dp), intent(in) :: b !! generalized-hyperbolic skew shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(genhyperbolic_logpdf(data, p, a, b, loc, scale))
    end function genhyperbolic_loglikelihood

    pure function genhyperbolic_nnlf(data, p, a, b, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: p !! real generalized-hyperbolic tail shape
        real(dp), intent(in) :: a !! positive generalized-hyperbolic shape
        real(dp), intent(in) :: b !! generalized-hyperbolic skew shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -genhyperbolic_loglikelihood(data, p, a, b, loc, scale)
    end function genhyperbolic_nnlf

    pure function nchypergeom_fisher_loglikelihood_real(data, m, n, draws, odds, loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(nchypergeom_fisher_logpmf(data, m, n, draws, odds, loc))
    end function nchypergeom_fisher_loglikelihood_real

    pure function nchypergeom_fisher_nnlf_real(data, m, n, draws, odds, loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -nchypergeom_fisher_loglikelihood_real(data, m, n, draws, odds, loc)
    end function nchypergeom_fisher_nnlf_real

    pure function nchypergeom_fisher_loglikelihood_int(data, m, n, draws, odds, loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        loglike = sum(nchypergeom_fisher_logpmf(data, m, n, draws, odds, loc))
    end function nchypergeom_fisher_loglikelihood_int

    pure function nchypergeom_fisher_nnlf_int(data, m, n, draws, odds, loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value = -nchypergeom_fisher_loglikelihood_int(data, m, n, draws, odds, loc)
    end function nchypergeom_fisher_nnlf_int

    pure function nct_loglikelihood(data, df, nc, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! positive noncentral-t degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(nct_logpdf(data, df, nc, loc, scale))
    end function nct_loglikelihood

    pure function nct_nnlf(data, df, nc, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! positive noncentral-t degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -nct_loglikelihood(data, df, nc, loc, scale)
    end function nct_nnlf

    pure function gausshyper_loglikelihood(data, a, b, c, z, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike
        loglike = sum(gausshyper_logpdf(data, a, b, c, z, loc, scale))
    end function gausshyper_loglikelihood

    pure function gausshyper_nnlf(data, a, b, c, z, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value
        value = -gausshyper_loglikelihood(data, a, b, c, z, loc, scale)
    end function gausshyper_nnlf

    pure function nchypergeom_wallenius_loglikelihood_real(data,m,n,draws,odds,loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        integer :: i
        loglike=0.0_dp
        do i=1,size(data)
            loglike=loglike+nchypergeom_wallenius_logpmf(data(i),m,n,draws,odds,loc)
        end do
    end function nchypergeom_wallenius_loglikelihood_real

    pure function nchypergeom_wallenius_nnlf_real(data,m,n,draws,odds,loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value=-nchypergeom_wallenius_loglikelihood_real(data,m,n,draws,odds,loc)
    end function nchypergeom_wallenius_nnlf_real

    pure function nchypergeom_wallenius_loglikelihood_int(data,m,n,draws,odds,loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        integer :: i
        loglike=0.0_dp
        do i=1,size(data)
            loglike=loglike+nchypergeom_wallenius_logpmf(data(i),m,n,draws,odds,loc)
        end do
    end function nchypergeom_wallenius_loglikelihood_int

    pure function nchypergeom_wallenius_nnlf_int(data,m,n,draws,odds,loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value=-nchypergeom_wallenius_loglikelihood_int(data,m,n,draws,odds,loc)
    end function nchypergeom_wallenius_nnlf_int

    pure function poisson_binom_loglikelihood_real(data,p,loc) result(loglike)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        integer :: i
        loglike=0.0_dp
        do i=1,size(data); loglike=loglike+poisson_binom_logpmf(data(i),p,loc); end do
    end function poisson_binom_loglikelihood_real

    pure function poisson_binom_nnlf_real(data,p,loc) result(value)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value=-poisson_binom_loglikelihood_real(data,p,loc)
    end function poisson_binom_nnlf_real

    pure function poisson_binom_loglikelihood_int(data,p,loc) result(loglike)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: loglike
        integer :: i
        loglike=0.0_dp
        do i=1,size(data); loglike=loglike+poisson_binom_logpmf(data(i),p,loc); end do
    end function poisson_binom_loglikelihood_int

    pure function poisson_binom_nnlf_int(data,p,loc) result(value)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: value
        value=-poisson_binom_loglikelihood_int(data,p,loc)
    end function poisson_binom_nnlf_int

    pure function levy_stable_loglikelihood(data, alpha, beta, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent Levy-stable observations
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(levy_stable_logpdf(data, alpha, beta, loc, scale))
    end function levy_stable_loglikelihood

    pure function levy_stable_nnlf(data, alpha, beta, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent Levy-stable observations
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -levy_stable_loglikelihood(data, alpha, beta, loc, scale)
    end function levy_stable_nnlf

    pure function studentized_range_loglikelihood(data, k, df, loc, scale) result(loglike)
        real(dp), intent(in) :: data(:) !! independent studentized-range observations
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: loglike

        loglike = sum(studentized_range_logpdf(data, k, df, loc, scale))
    end function studentized_range_loglikelihood

    pure function studentized_range_nnlf(data, k, df, loc, scale) result(value)
        real(dp), intent(in) :: data(:) !! independent studentized-range observations
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: value

        value = -studentized_range_loglikelihood(data, k, df, loc, scale)
    end function studentized_range_nnlf

end module scifort_likelihood
