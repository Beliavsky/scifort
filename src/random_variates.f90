! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_random_variates
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_uniform
    use scifort_normal, only : normal_ppf
    use scifort_uniform, only : uniform_ppf
    use scifort_exponential, only : exponential_ppf
    use scifort_laplace, only : laplace_ppf
    use scifort_logistic, only : logistic_ppf
    use scifort_cauchy, only : cauchy_ppf
    use scifort_rayleigh, only : rayleigh_ppf
    use scifort_gamma, only : gamma_ppf
    use scifort_chi2, only : chi2_ppf
    use scifort_student_t, only : t_ppf
    use scifort_lognormal, only : lognormal_ppf
    use scifort_weibull, only : weibull_ppf
    use scifort_pareto, only : pareto_ppf
    use scifort_beta, only : beta_ppf
    use scifort_f_distribution, only : f_ppf
    use scifort_gumbel_r, only : gumbel_r_ppf
    use scifort_gumbel_l, only : gumbel_l_ppf
    use scifort_powerlaw, only : powerlaw_ppf
    use scifort_triang, only : triang_ppf
    use scifort_genpareto, only : genpareto_ppf
    use scifort_arcsine, only : arcsine_ppf
    use scifort_halfnorm, only : halfnorm_ppf
    use scifort_halfcauchy, only : halfcauchy_ppf
    use scifort_lomax, only : lomax_ppf
    use scifort_chi, only : chi_ppf
    use scifort_maxwell, only : maxwell_ppf
    use scifort_cosine, only : cosine_ppf
    use scifort_semicircular, only : semicircular_ppf
    use scifort_anglit, only : anglit_ppf
    use scifort_moyal, only : moyal_ppf
    use scifort_hypsecant, only : hypsecant_ppf
    use scifort_halflogistic, only : halflogistic_ppf
    use scifort_invgamma, only : invgamma_ppf
    use scifort_invgauss, only : invgauss_ppf
    use scifort_levy, only : levy_ppf
    use scifort_loglaplace, only : loglaplace_ppf
    use scifort_bradford, only : bradford_ppf
    use scifort_truncexpon, only : truncexpon_ppf
    use scifort_fisk, only : fisk_ppf
    use scifort_dweibull, only : dweibull_ppf
    use scifort_alpha, only : alpha_ppf
    use scifort_fatiguelife, only : fatiguelife_ppf
    use scifort_genlogistic, only : genlogistic_ppf
    use scifort_gennorm, only : gennorm_ppf
    use scifort_nakagami, only : nakagami_ppf
    use scifort_powernorm, only : powernorm_ppf
    use scifort_loggamma, only : loggamma_ppf
    use scifort_wald, only : wald_ppf
    use scifort_gompertz, only : gompertz_ppf
    use scifort_invweibull, only : invweibull_ppf
    use scifort_betaprime, only : betaprime_ppf
    use scifort_burr12, only : burr12_ppf
    use scifort_genhalflogistic, only : genhalflogistic_ppf
    use scifort_exponpow, only : exponpow_ppf
    use scifort_exponweib, only : exponweib_ppf
    use scifort_powerlognorm, only : powerlognorm_ppf
    use scifort_levy_l, only : levy_l_ppf
    use scifort_weibull_max, only : weibull_max_ppf
    use scifort_rdist, only : rdist_ppf
    use scifort_skewcauchy, only : skewcauchy_ppf
    use scifort_dgamma, only : dgamma_ppf
    use scifort_laplace_asymmetric, only : laplace_asymmetric_ppf
    use scifort_truncnorm, only : truncnorm_ppf
    use scifort_loguniform, only : loguniform_ppf
    use scifort_foldnorm, only : foldnorm_ppf
    use scifort_foldcauchy, only : foldcauchy_ppf
    use scifort_recipinvgauss, only : recipinvgauss_ppf
    use scifort_truncpareto, only : truncpareto_ppf
    use scifort_exponnorm, only : exponnorm_ppf
    use scifort_johnsonsb, only : johnsonsb_ppf
    use scifort_johnsonsu, only : johnsonsu_ppf
    use scifort_trapezoid, only : trapezoid_ppf
    use scifort_burr, only : burr_ppf
    use scifort_mielke, only : mielke_ppf
    use scifort_gibrat, only : gibrat_ppf
    use scifort_wrapcauchy, only : wrapcauchy_ppf
    use scifort_genextreme, only : genextreme_ppf
    use scifort_kappa3, only : kappa3_ppf
    use scifort_kappa4, only : kappa4_ppf
    use scifort_truncweibull_min, only : truncweibull_min_ppf
    use scifort_gengamma, only : gengamma_ppf
    use scifort_halfgennorm, only : halfgennorm_ppf
    use scifort_argus, only : argus_ppf
    use scifort_erlang, only : erlang_ppf
    use scifort_crystalball, only : crystalball_ppf
    use scifort_jf_skew_t, only : jf_skew_t_ppf
    use scifort_pearson3, only : pearson3_ppf
    use scifort_rel_breitwigner, only : rel_breitwigner_ppf
    use scifort_genexpon, only : genexpon_ppf
    use scifort_skewnorm, only : skewnorm_ppf
    use scifort_tukeylambda, only : tukeylambda_ppf
    use scifort_rice, only : rice_ppf
    use scifort_dpareto_lognorm, only : dpareto_lognorm_ppf
    use scifort_vonmises, only : vonmises_ppf
    use scifort_vonmises_line, only : vonmises_line_ppf
    use scifort_kstwobign, only : kstwobign_ppf
    use scifort_irwinhall, only : irwinhall_ppf
    use scifort_ksone, only : ksone_ppf
    use scifort_kstwo, only : kstwo_ppf
    use scifort_levy_stable, only : levy_stable_ppf
    use scifort_studentized_range, only : studentized_range_ppf
    use scifort_ncx2, only : ncx2_ppf
    use scifort_ncf, only : ncf_ppf
    use scifort_bernoulli, only : bernoulli_ppf
    use scifort_poisson, only : poisson_ppf
    use scifort_geometric, only : geometric_ppf
    use scifort_binomial, only : binomial_ppf
    use scifort_negative_binomial, only : negative_binomial_ppf
    use scifort_randint, only : randint_ppf
    use scifort_planck, only : planck_ppf
    use scifort_dlaplace, only : dlaplace_ppf
    use scifort_logser, only : logser_ppf
    use scifort_betabinom, only : betabinom_ppf
    use scifort_hypergeom, only : hypergeom_ppf
    use scifort_nhypergeom, only : nhypergeom_ppf
    use scifort_boltzmann, only : boltzmann_ppf
    use scifort_betanbinom, only : betanbinom_ppf
    use scifort_yulesimon, only : yulesimon_ppf
    use scifort_zipf, only : zipf_ppf
    use scifort_zipfian, only : zipfian_ppf
    use scifort_geninvgauss, only : geninvgauss_ppf
    use scifort_norminvgauss, only : norminvgauss_ppf
    use scifort_skellam, only : skellam_ppf
    use scifort_genhyperbolic, only : genhyperbolic_ppf
    use scifort_nchypergeom_fisher, only : nchypergeom_fisher_ppf
    use scifort_nchypergeom_wallenius, only : nchypergeom_wallenius_ppf
    use scifort_poisson_binom, only : poisson_binom_ppf
    use scifort_nct, only : nct_ppf
    use scifort_gausshyper, only : gausshyper_ppf
    use scifort_landau, only : landau_ppf
    implicit none
    private

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
    public :: bradford_rvs
    public :: bradford_rvs_array
    public :: truncexpon_rvs
    public :: truncexpon_rvs_array
    public :: fisk_rvs
    public :: fisk_rvs_array
    public :: dweibull_rvs
    public :: dweibull_rvs_array
    public :: alpha_rvs
    public :: alpha_rvs_array
    public :: fatiguelife_rvs
    public :: fatiguelife_rvs_array
    public :: genlogistic_rvs
    public :: genlogistic_rvs_array
    public :: gennorm_rvs
    public :: gennorm_rvs_array
    public :: nakagami_rvs
    public :: nakagami_rvs_array
    public :: powernorm_rvs
    public :: powernorm_rvs_array
    public :: loggamma_rvs
    public :: loggamma_rvs_array
    public :: wald_rvs
    public :: wald_rvs_array
    public :: gompertz_rvs
    public :: gompertz_rvs_array
    public :: invweibull_rvs
    public :: invweibull_rvs_array
    public :: betaprime_rvs
    public :: betaprime_rvs_array
    public :: burr12_rvs
    public :: burr12_rvs_array
    public :: genhalflogistic_rvs
    public :: genhalflogistic_rvs_array
    public :: exponpow_rvs
    public :: exponpow_rvs_array
    public :: exponweib_rvs
    public :: exponweib_rvs_array
    public :: powerlognorm_rvs
    public :: powerlognorm_rvs_array
    public :: levy_l_rvs
    public :: levy_l_rvs_array
    public :: weibull_max_rvs
    public :: weibull_max_rvs_array
    public :: rdist_rvs
    public :: rdist_rvs_array
    public :: skewcauchy_rvs
    public :: skewcauchy_rvs_array
    public :: dgamma_rvs
    public :: dgamma_rvs_array
    public :: laplace_asymmetric_rvs
    public :: laplace_asymmetric_rvs_array
    public :: truncnorm_rvs
    public :: truncnorm_rvs_array
    public :: loguniform_rvs
    public :: loguniform_rvs_array
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

    public :: foldnorm_rvs, foldnorm_rvs_array
    public :: foldcauchy_rvs, foldcauchy_rvs_array
    public :: recipinvgauss_rvs, recipinvgauss_rvs_array
    public :: truncpareto_rvs, truncpareto_rvs_array
    public :: exponnorm_rvs, exponnorm_rvs_array
    public :: johnsonsb_rvs, johnsonsb_rvs_array
    public :: johnsonsu_rvs, johnsonsu_rvs_array
    public :: trapezoid_rvs, trapezoid_rvs_array
    public :: burr_rvs, burr_rvs_array
    public :: mielke_rvs, mielke_rvs_array
    public :: gibrat_rvs, gibrat_rvs_array
    public :: wrapcauchy_rvs, wrapcauchy_rvs_array
    public :: genextreme_rvs, genextreme_rvs_array
    public :: kappa3_rvs, kappa3_rvs_array
    public :: kappa4_rvs, kappa4_rvs_array
    public :: truncweibull_min_rvs, truncweibull_min_rvs_array
    public :: gengamma_rvs, gengamma_rvs_array
    public :: halfgennorm_rvs, halfgennorm_rvs_array
    public :: argus_rvs, argus_rvs_array
    public :: erlang_rvs, erlang_rvs_array
    public :: crystalball_rvs, crystalball_rvs_array
    public :: jf_skew_t_rvs, jf_skew_t_rvs_array
    public :: pearson3_rvs, pearson3_rvs_array
    public :: rel_breitwigner_rvs, rel_breitwigner_rvs_array
    public :: genexpon_rvs, genexpon_rvs_array
    public :: skewnorm_rvs, skewnorm_rvs_array
    public :: tukeylambda_rvs, tukeylambda_rvs_array
    public :: rice_rvs, rice_rvs_array
    public :: dpareto_lognorm_rvs, dpareto_lognorm_rvs_array
    public :: vonmises_rvs, vonmises_rvs_array
    public :: vonmises_line_rvs, vonmises_line_rvs_array
    public :: kstwobign_rvs, kstwobign_rvs_array
    public :: irwinhall_rvs, irwinhall_rvs_array
    public :: ksone_rvs, ksone_rvs_array
    public :: kstwo_rvs, kstwo_rvs_array
    public :: levy_stable_rvs, levy_stable_rvs_array
    public :: studentized_range_rvs, studentized_range_rvs_array
    public :: ncx2_rvs, ncx2_rvs_array
    public :: ncf_rvs, ncf_rvs_array
    public :: randint_rvs, randint_rvs_array
    public :: planck_rvs, planck_rvs_array
    public :: dlaplace_rvs, dlaplace_rvs_array
    public :: logser_rvs, logser_rvs_array
    public :: betabinom_rvs, betabinom_rvs_array
    public :: hypergeom_rvs, hypergeom_rvs_array
    public :: nhypergeom_rvs, nhypergeom_rvs_array
    public :: boltzmann_rvs, boltzmann_rvs_array
    public :: betanbinom_rvs, betanbinom_rvs_array
    public :: yulesimon_rvs, yulesimon_rvs_array
    public :: zipf_rvs, zipf_rvs_array
    public :: zipfian_rvs, zipfian_rvs_array
    public :: geninvgauss_rvs, geninvgauss_rvs_array
    public :: norminvgauss_rvs, norminvgauss_rvs_array
    public :: skellam_rvs, skellam_rvs_array
    public :: genhyperbolic_rvs, genhyperbolic_rvs_array
    public :: nchypergeom_fisher_rvs, nchypergeom_fisher_rvs_array
    public :: nct_rvs, nct_rvs_array
    public :: gausshyper_rvs, gausshyper_rvs_array
    public :: landau_rvs, landau_rvs_array
    public :: nchypergeom_wallenius_rvs, nchypergeom_wallenius_rvs_array
    public :: poisson_binom_rvs, poisson_binom_rvs_array

contains

    function normal_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = normal_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = normal_ppf(rng_uniform(state), mu, sigma)
    end function normal_rvs

    subroutine normal_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = normal_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = normal_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine normal_rvs_array

    function uniform_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = uniform_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = uniform_ppf(rng_uniform(state), mu, sigma)
    end function uniform_rvs

    subroutine uniform_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = uniform_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = uniform_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine uniform_rvs_array

    function exponential_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponential_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = exponential_ppf(rng_uniform(state), mu, sigma)
    end function exponential_rvs

    subroutine exponential_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponential_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = exponential_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine exponential_rvs_array

    function laplace_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = laplace_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = laplace_ppf(rng_uniform(state), mu, sigma)
    end function laplace_rvs

    subroutine laplace_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = laplace_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = laplace_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine laplace_rvs_array

    function logistic_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = logistic_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = logistic_ppf(rng_uniform(state), mu, sigma)
    end function logistic_rvs

    subroutine logistic_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = logistic_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = logistic_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine logistic_rvs_array

    function cauchy_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = cauchy_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = cauchy_ppf(rng_uniform(state), mu, sigma)
    end function cauchy_rvs

    subroutine cauchy_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = cauchy_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = cauchy_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine cauchy_rvs_array

    function rayleigh_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = rayleigh_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = rayleigh_ppf(rng_uniform(state), mu, sigma)
    end function rayleigh_rvs

    subroutine rayleigh_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = rayleigh_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = rayleigh_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine rayleigh_rvs_array

    function gamma_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = gamma_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gamma_ppf(rng_uniform(state), a, mu, sigma)
    end function gamma_rvs

    subroutine gamma_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = gamma_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = gamma_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine gamma_rvs_array

    function chi2_rvs(state, df, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = chi2_ppf(0.5_dp, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = chi2_ppf(rng_uniform(state), df, mu, sigma)
    end function chi2_rvs

    subroutine chi2_rvs_array(state, samples, df, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = chi2_ppf(0.5_dp, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = chi2_ppf(rng_uniform(state), df, mu, sigma)
        end do
    end subroutine chi2_rvs_array

    function t_rvs(state, df, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = t_ppf(0.5_dp, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = t_ppf(rng_uniform(state), df, mu, sigma)
    end function t_rvs

    subroutine t_rvs_array(state, samples, df, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: df !! degrees of freedom
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = t_ppf(0.5_dp, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = t_ppf(rng_uniform(state), df, mu, sigma)
        end do
    end subroutine t_rvs_array

    function lognormal_rvs(state, s, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: s !! positive lognormal shape
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = lognormal_ppf(0.5_dp, s, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = lognormal_ppf(rng_uniform(state), s, mu, sigma)
    end function lognormal_rvs

    subroutine lognormal_rvs_array(state, samples, s, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: s !! positive lognormal shape
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = lognormal_ppf(0.5_dp, s, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = lognormal_ppf(rng_uniform(state), s, mu, sigma)
        end do
    end subroutine lognormal_rvs_array

    function weibull_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive Weibull shape
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = weibull_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = weibull_ppf(rng_uniform(state), c, mu, sigma)
    end function weibull_rvs

    subroutine weibull_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive Weibull shape
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = weibull_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = weibull_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine weibull_rvs_array

    function pareto_rvs(state, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: b !! positive second shape parameter
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = pareto_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = pareto_ppf(rng_uniform(state), b, mu, sigma)
    end function pareto_rvs

    subroutine pareto_rvs_array(state, samples, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: b !! positive second shape parameter
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = pareto_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = pareto_ppf(rng_uniform(state), b, mu, sigma)
        end do
    end subroutine pareto_rvs_array

    function beta_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! positive second shape parameter
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = beta_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = beta_ppf(rng_uniform(state), a, b, mu, sigma)
    end function beta_rvs

    subroutine beta_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive shape parameter
        real(dp), intent(in) :: b !! positive second shape parameter
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = beta_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = beta_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine beta_rvs_array

    function f_rvs(state, dfn, dfd, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: dfn !! numerator degrees of freedom
        real(dp), intent(in) :: dfd !! denominator degrees of freedom
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = f_ppf(0.5_dp, dfn, dfd, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = f_ppf(rng_uniform(state), dfn, dfd, mu, sigma)
    end function f_rvs

    subroutine f_rvs_array(state, samples, dfn, dfd, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: dfn !! numerator degrees of freedom
        real(dp), intent(in) :: dfd !! denominator degrees of freedom
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = f_ppf(0.5_dp, dfn, dfd, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = f_ppf(rng_uniform(state), dfn, dfd, mu, sigma)
        end do
    end subroutine f_rvs_array

    function gumbel_r_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = gumbel_r_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gumbel_r_ppf(rng_uniform(state), mu, sigma)
    end function gumbel_r_rvs

    subroutine gumbel_r_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = gumbel_r_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = gumbel_r_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine gumbel_r_rvs_array

    function gumbel_l_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = gumbel_l_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gumbel_l_ppf(rng_uniform(state), mu, sigma)
    end function gumbel_l_rvs

    subroutine gumbel_l_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = gumbel_l_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = gumbel_l_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine gumbel_l_rvs_array

    function powerlaw_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = powerlaw_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = powerlaw_ppf(rng_uniform(state), a, mu, sigma)
    end function powerlaw_rvs

    subroutine powerlaw_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = powerlaw_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = powerlaw_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine powerlaw_rvs_array

    function triang_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = triang_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = triang_ppf(rng_uniform(state), c, mu, sigma)
    end function triang_rvs

    subroutine triang_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! standardized mode in [0, 1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = triang_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = triang_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine triang_rvs_array

    function genpareto_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = genpareto_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = genpareto_ppf(rng_uniform(state), c, mu, sigma)
    end function genpareto_rvs

    subroutine genpareto_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = genpareto_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = genpareto_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine genpareto_rvs_array

    function arcsine_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = arcsine_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = arcsine_ppf(rng_uniform(state), mu, sigma)
    end function arcsine_rvs

    subroutine arcsine_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = arcsine_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = arcsine_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine arcsine_rvs_array

    function halfnorm_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = halfnorm_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = halfnorm_ppf(rng_uniform(state), mu, sigma)
    end function halfnorm_rvs

    subroutine halfnorm_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = halfnorm_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = halfnorm_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine halfnorm_rvs_array

    function halfcauchy_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = halfcauchy_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = halfcauchy_ppf(rng_uniform(state), mu, sigma)
    end function halfcauchy_rvs

    subroutine halfcauchy_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = halfcauchy_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = halfcauchy_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine halfcauchy_rvs_array

    function lomax_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = lomax_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = lomax_ppf(rng_uniform(state), c, mu, sigma)
    end function lomax_rvs

    subroutine lomax_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = lomax_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = lomax_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine lomax_rvs_array

    function chi_rvs(state, df, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = chi_ppf(0.5_dp, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = chi_ppf(rng_uniform(state), df, mu, sigma)
    end function chi_rvs

    subroutine chi_rvs_array(state, samples, df, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = chi_ppf(0.5_dp, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = chi_ppf(rng_uniform(state), df, mu, sigma)
        end do
    end subroutine chi_rvs_array

    function maxwell_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = maxwell_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = maxwell_ppf(rng_uniform(state), mu, sigma)
    end function maxwell_rvs

    subroutine maxwell_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = maxwell_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = maxwell_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine maxwell_rvs_array

    function cosine_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = cosine_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = cosine_ppf(rng_uniform(state), mu, sigma)
    end function cosine_rvs

    subroutine cosine_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = cosine_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = cosine_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine cosine_rvs_array

    function semicircular_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = semicircular_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = semicircular_ppf(rng_uniform(state), mu, sigma)
    end function semicircular_rvs

    subroutine semicircular_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = semicircular_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = semicircular_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine semicircular_rvs_array

    function anglit_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = anglit_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = anglit_ppf(rng_uniform(state), mu, sigma)
    end function anglit_rvs

    subroutine anglit_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = anglit_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = anglit_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine anglit_rvs_array

    function moyal_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = moyal_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = moyal_ppf(rng_uniform(state), mu, sigma)
    end function moyal_rvs

    subroutine moyal_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = moyal_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = moyal_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine moyal_rvs_array

    function bernoulli_rvs(state, p, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp) :: x

        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = bernoulli_ppf(0.5_dp, p, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = bernoulli_ppf(rng_uniform(state), p, shift)
    end function bernoulli_rvs

    subroutine bernoulli_rvs_array(state, samples, p, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)

        integer :: i
        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = bernoulli_ppf(0.5_dp, p, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = bernoulli_ppf(rng_uniform(state), p, shift)
        end do
    end subroutine bernoulli_rvs_array

    function poisson_rvs(state, mu, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: mu !! Poisson mean
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp) :: x

        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = poisson_ppf(0.5_dp, mu, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = poisson_ppf(rng_uniform(state), mu, shift)
    end function poisson_rvs

    subroutine poisson_rvs_array(state, samples, mu, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: mu !! Poisson mean
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)

        integer :: i
        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = poisson_ppf(0.5_dp, mu, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = poisson_ppf(rng_uniform(state), mu, shift)
        end do
    end subroutine poisson_rvs_array

    function geometric_rvs(state, p, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp) :: x

        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = geometric_ppf(0.5_dp, p, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = geometric_ppf(rng_uniform(state), p, shift)
    end function geometric_rvs

    subroutine geometric_rvs_array(state, samples, p, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)

        integer :: i
        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = geometric_ppf(0.5_dp, p, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = geometric_ppf(rng_uniform(state), p, shift)
        end do
    end subroutine geometric_rvs_array

    function binomial_rvs(state, n, p, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: n !! trial count or success shape
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp) :: x

        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = binomial_ppf(0.5_dp, n, p, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = binomial_ppf(rng_uniform(state), n, p, shift)
    end function binomial_rvs

    subroutine binomial_rvs_array(state, samples, n, p, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: n !! trial count or success shape
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)

        integer :: i
        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = binomial_ppf(0.5_dp, n, p, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = binomial_ppf(rng_uniform(state), n, p, shift)
        end do
    end subroutine binomial_rvs_array

    function negative_binomial_rvs(state, n, p, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: n !! trial count or success shape
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)
        real(dp) :: x

        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = negative_binomial_ppf(0.5_dp, n, p, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = negative_binomial_ppf(rng_uniform(state), n, p, shift)
    end function negative_binomial_rvs

    subroutine negative_binomial_rvs_array(state, samples, n, p, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: n !! trial count or success shape
        real(dp), intent(in) :: p !! success probability
        real(dp), intent(in), optional :: loc !! location or support shift (default 0)

        integer :: i
        real(dp) :: shift
        real(dp) :: probe

        shift = get_shift(loc)
        probe = negative_binomial_ppf(0.5_dp, n, p, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if

        do i = 1, size(samples)
            samples(i) = negative_binomial_ppf(rng_uniform(state), n, p, shift)
        end do
    end subroutine negative_binomial_rvs_array

    function hypsecant_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = hypsecant_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = hypsecant_ppf(rng_uniform(state), mu, sigma)
    end function hypsecant_rvs

    subroutine hypsecant_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = hypsecant_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = hypsecant_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine hypsecant_rvs_array

    function halflogistic_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = halflogistic_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = halflogistic_ppf(rng_uniform(state), mu, sigma)
    end function halflogistic_rvs

    subroutine halflogistic_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = halflogistic_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = halflogistic_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine halflogistic_rvs_array

    function invgamma_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive inverse-gamma shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = invgamma_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = invgamma_ppf(rng_uniform(state), a, mu, sigma)
    end function invgamma_rvs

    subroutine invgamma_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive inverse-gamma shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = invgamma_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = invgamma_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine invgamma_rvs_array

    function invgauss_rvs(state, mu_shape, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = invgauss_ppf(0.5_dp, mu_shape, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = invgauss_ppf(rng_uniform(state), mu_shape, mu, sigma)
    end function invgauss_rvs

    subroutine invgauss_rvs_array(state, samples, mu_shape, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = invgauss_ppf(0.5_dp, mu_shape, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = invgauss_ppf(rng_uniform(state), mu_shape, mu, sigma)
        end do
    end subroutine invgauss_rvs_array

    function levy_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = levy_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = levy_ppf(rng_uniform(state), mu, sigma)
    end function levy_rvs

    subroutine levy_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = levy_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = levy_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine levy_rvs_array

    function loglaplace_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive log-Laplace shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = loglaplace_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = loglaplace_ppf(rng_uniform(state), c, mu, sigma)
    end function loglaplace_rvs

    subroutine loglaplace_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive log-Laplace shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = loglaplace_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = loglaplace_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine loglaplace_rvs_array


    function bradford_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive Bradford shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = bradford_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = bradford_ppf(rng_uniform(state), c, mu, sigma)
    end function bradford_rvs

    subroutine bradford_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive Bradford shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = bradford_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = bradford_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine bradford_rvs_array

    function truncexpon_rvs(state, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: b !! positive standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncexpon_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = truncexpon_ppf(rng_uniform(state), b, mu, sigma)
    end function truncexpon_rvs

    subroutine truncexpon_rvs_array(state, samples, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: b !! positive standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncexpon_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = truncexpon_ppf(rng_uniform(state), b, mu, sigma)
        end do
    end subroutine truncexpon_rvs_array

    function fisk_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive Fisk shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = fisk_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = fisk_ppf(rng_uniform(state), c, mu, sigma)
    end function fisk_rvs

    subroutine fisk_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive Fisk shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = fisk_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = fisk_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine fisk_rvs_array

    function dweibull_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive double-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = dweibull_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = dweibull_ppf(rng_uniform(state), c, mu, sigma)
    end function dweibull_rvs

    subroutine dweibull_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive double-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)

        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = dweibull_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = dweibull_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine dweibull_rvs_array

    function alpha_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive alpha shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = alpha_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = alpha_ppf(rng_uniform(state), a, mu, sigma)
    end function alpha_rvs

    subroutine alpha_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive alpha shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = alpha_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = alpha_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine alpha_rvs_array

    function fatiguelife_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive fatigue-life shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = fatiguelife_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = fatiguelife_ppf(rng_uniform(state), c, mu, sigma)
    end function fatiguelife_rvs

    subroutine fatiguelife_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive fatigue-life shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = fatiguelife_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = fatiguelife_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine fatiguelife_rvs_array

    function genlogistic_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive generalized-logistic shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genlogistic_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = genlogistic_ppf(rng_uniform(state), c, mu, sigma)
    end function genlogistic_rvs

    subroutine genlogistic_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive generalized-logistic shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genlogistic_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = genlogistic_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine genlogistic_rvs_array

    function gennorm_rvs(state, beta, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: beta !! positive generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gennorm_ppf(0.5_dp, beta, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gennorm_ppf(rng_uniform(state), beta, mu, sigma)
    end function gennorm_rvs

    subroutine gennorm_rvs_array(state, samples, beta, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: beta !! positive generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gennorm_ppf(0.5_dp, beta, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = gennorm_ppf(rng_uniform(state), beta, mu, sigma)
        end do
    end subroutine gennorm_rvs_array

    function nakagami_rvs(state, nu, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = nakagami_ppf(0.5_dp, nu, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = nakagami_ppf(rng_uniform(state), nu, mu, sigma)
    end function nakagami_rvs

    subroutine nakagami_rvs_array(state, samples, nu, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = nakagami_ppf(0.5_dp, nu, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = nakagami_ppf(rng_uniform(state), nu, mu, sigma)
        end do
    end subroutine nakagami_rvs_array

    function powernorm_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive power-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = powernorm_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = powernorm_ppf(rng_uniform(state), c, mu, sigma)
    end function powernorm_rvs

    subroutine powernorm_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive power-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = powernorm_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = powernorm_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine powernorm_rvs_array

    function loggamma_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = loggamma_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = loggamma_ppf(rng_uniform(state), c, mu, sigma)
    end function loggamma_rvs

    subroutine loggamma_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = loggamma_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = loggamma_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine loggamma_rvs_array

    function wald_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = wald_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = wald_ppf(rng_uniform(state), mu, sigma)
    end function wald_rvs

    subroutine wald_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = wald_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = wald_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine wald_rvs_array


    function gompertz_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive Gompertz shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gompertz_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gompertz_ppf(rng_uniform(state), c, mu, sigma)
    end function gompertz_rvs

    subroutine gompertz_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive Gompertz shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gompertz_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = gompertz_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine gompertz_rvs_array

    function invweibull_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = invweibull_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = invweibull_ppf(rng_uniform(state), c, mu, sigma)
    end function invweibull_rvs

    subroutine invweibull_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = invweibull_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = invweibull_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine invweibull_rvs_array

    function betaprime_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive first beta-prime shape parameter
        real(dp), intent(in) :: b !! positive second beta-prime shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = betaprime_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = betaprime_ppf(rng_uniform(state), a, b, mu, sigma)
    end function betaprime_rvs

    subroutine betaprime_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive first beta-prime shape parameter
        real(dp), intent(in) :: b !! positive second beta-prime shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = betaprime_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = betaprime_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine betaprime_rvs_array

    function burr12_rvs(state, c, d, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive first Burr-XII shape parameter
        real(dp), intent(in) :: d !! positive second Burr-XII shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = burr12_ppf(0.5_dp, c, d, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = burr12_ppf(rng_uniform(state), c, d, mu, sigma)
    end function burr12_rvs

    subroutine burr12_rvs_array(state, samples, c, d, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive first Burr-XII shape parameter
        real(dp), intent(in) :: d !! positive second Burr-XII shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = burr12_ppf(0.5_dp, c, d, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = burr12_ppf(rng_uniform(state), c, d, mu, sigma)
        end do
    end subroutine burr12_rvs_array


    function genhalflogistic_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive generalized half-logistic shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genhalflogistic_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = genhalflogistic_ppf(rng_uniform(state), c, mu, sigma)
    end function genhalflogistic_rvs

    subroutine genhalflogistic_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive generalized half-logistic shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genhalflogistic_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = genhalflogistic_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine genhalflogistic_rvs_array

    function exponpow_rvs(state, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: b !! positive exponential-power shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponpow_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = exponpow_ppf(rng_uniform(state), b, mu, sigma)
    end function exponpow_rvs

    subroutine exponpow_rvs_array(state, samples, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: b !! positive exponential-power shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponpow_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = exponpow_ppf(rng_uniform(state), b, mu, sigma)
        end do
    end subroutine exponpow_rvs_array

    function exponweib_rvs(state, a, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive exponentiation shape parameter
        real(dp), intent(in) :: c !! positive Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponweib_ppf(0.5_dp, a, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = exponweib_ppf(rng_uniform(state), a, c, mu, sigma)
    end function exponweib_rvs

    subroutine exponweib_rvs_array(state, samples, a, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive exponentiation shape parameter
        real(dp), intent(in) :: c !! positive Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponweib_ppf(0.5_dp, a, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = exponweib_ppf(rng_uniform(state), a, c, mu, sigma)
        end do
    end subroutine exponweib_rvs_array

    function powerlognorm_rvs(state, c, s, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive power shape parameter
        real(dp), intent(in) :: s !! positive lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = powerlognorm_ppf(0.5_dp, c, s, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = powerlognorm_ppf(rng_uniform(state), c, s, mu, sigma)
    end function powerlognorm_rvs

    subroutine powerlognorm_rvs_array(state, samples, c, s, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive power shape parameter
        real(dp), intent(in) :: s !! positive lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = powerlognorm_ppf(0.5_dp, c, s, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = powerlognorm_ppf(rng_uniform(state), c, s, mu, sigma)
        end do
    end subroutine powerlognorm_rvs_array

    function levy_l_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = levy_l_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = levy_l_ppf(rng_uniform(state), mu, sigma)
    end function levy_l_rvs

    subroutine levy_l_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = levy_l_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = levy_l_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine levy_l_rvs_array

    function weibull_max_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive Weibull-maximum shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = weibull_max_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = weibull_max_ppf(rng_uniform(state), c, mu, sigma)
    end function weibull_max_rvs

    subroutine weibull_max_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive Weibull-maximum shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = weibull_max_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = weibull_max_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine weibull_max_rvs_array

    function rdist_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive R-distribution shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = rdist_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = rdist_ppf(rng_uniform(state), c, mu, sigma)
    end function rdist_rvs

    subroutine rdist_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive R-distribution shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = rdist_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = rdist_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine rdist_rvs_array

    function skewcauchy_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = skewcauchy_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = skewcauchy_ppf(rng_uniform(state), a, mu, sigma)
    end function skewcauchy_rvs

    subroutine skewcauchy_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = skewcauchy_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = skewcauchy_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine skewcauchy_rvs_array


    function dgamma_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = dgamma_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = dgamma_ppf(rng_uniform(state), a, mu, sigma)
    end function dgamma_rvs

    subroutine dgamma_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = dgamma_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = dgamma_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine dgamma_rvs_array

    function laplace_asymmetric_rvs(state, kappa, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = laplace_asymmetric_ppf(0.5_dp, kappa, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = laplace_asymmetric_ppf(rng_uniform(state), kappa, mu, sigma)
    end function laplace_asymmetric_rvs

    subroutine laplace_asymmetric_rvs_array(state, samples, kappa, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = laplace_asymmetric_ppf(0.5_dp, kappa, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = laplace_asymmetric_ppf(rng_uniform(state), kappa, mu, sigma)
        end do
    end subroutine laplace_asymmetric_rvs_array

    function truncnorm_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncnorm_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = truncnorm_ppf(rng_uniform(state), a, b, mu, sigma)
    end function truncnorm_rvs

    subroutine truncnorm_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncnorm_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = truncnorm_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine truncnorm_rvs_array

    function loguniform_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = loguniform_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = loguniform_ppf(rng_uniform(state), a, b, mu, sigma)
    end function loguniform_rvs

    subroutine loguniform_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = loguniform_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = loguniform_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine loguniform_rvs_array


    function foldnorm_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = foldnorm_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = foldnorm_ppf(rng_uniform(state), c, mu, sigma)
    end function foldnorm_rvs

    subroutine foldnorm_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = foldnorm_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = foldnorm_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine foldnorm_rvs_array

    function foldcauchy_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = foldcauchy_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = foldcauchy_ppf(rng_uniform(state), c, mu, sigma)
    end function foldcauchy_rvs

    subroutine foldcauchy_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = foldcauchy_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = foldcauchy_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine foldcauchy_rvs_array

    function recipinvgauss_rvs(state, mu_shape, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: mu_shape !! positive reciprocal inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = recipinvgauss_ppf(0.5_dp, mu_shape, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = recipinvgauss_ppf(rng_uniform(state), mu_shape, mu, sigma)
    end function recipinvgauss_rvs

    subroutine recipinvgauss_rvs_array(state, samples, mu_shape, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: mu_shape !! positive reciprocal inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = recipinvgauss_ppf(0.5_dp, mu_shape, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = recipinvgauss_ppf(rng_uniform(state), mu_shape, mu, sigma)
        end do
    end subroutine recipinvgauss_rvs_array

    function truncpareto_rvs(state, b, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncpareto_ppf(0.5_dp, b, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = truncpareto_ppf(rng_uniform(state), b, c, mu, sigma)
    end function truncpareto_rvs

    subroutine truncpareto_rvs_array(state, samples, b, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncpareto_ppf(0.5_dp, b, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = truncpareto_ppf(rng_uniform(state), b, c, mu, sigma)
        end do
    end subroutine truncpareto_rvs_array


    function exponnorm_rvs(state, k, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: k !! positive exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponnorm_ppf(0.5_dp, k, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = exponnorm_ppf(rng_uniform(state), k, mu, sigma)
    end function exponnorm_rvs

    subroutine exponnorm_rvs_array(state, samples, k, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: k !! positive exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = exponnorm_ppf(0.5_dp, k, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = exponnorm_ppf(rng_uniform(state), k, mu, sigma)
        end do
    end subroutine exponnorm_rvs_array

    function johnsonsb_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = johnsonsb_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = johnsonsb_ppf(rng_uniform(state), a, b, mu, sigma)
    end function johnsonsb_rvs

    subroutine johnsonsb_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = johnsonsb_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = johnsonsb_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine johnsonsb_rvs_array

    function johnsonsu_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = johnsonsu_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = johnsonsu_ppf(rng_uniform(state), a, b, mu, sigma)
    end function johnsonsu_rvs

    subroutine johnsonsu_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = johnsonsu_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = johnsonsu_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine johnsonsu_rvs_array

    function trapezoid_rvs(state, c, d, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! standardized left plateau edge in [0,1]
        real(dp), intent(in) :: d !! standardized right plateau edge in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: x, mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = trapezoid_ppf(0.5_dp, c, d, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = trapezoid_ppf(rng_uniform(state), c, d, mu, sigma)
    end function trapezoid_rvs

    subroutine trapezoid_rvs_array(state, samples, c, d, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! standardized left plateau edge in [0,1]
        real(dp), intent(in) :: d !! standardized right plateau edge in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = trapezoid_ppf(0.5_dp, c, d, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = trapezoid_ppf(rng_uniform(state), c, d, mu, sigma)
        end do
    end subroutine trapezoid_rvs_array

    function burr_rvs(state, c, d, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = burr_ppf(0.5_dp, c, d, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = burr_ppf(rng_uniform(state), c, d, mu, sigma)
    end function burr_rvs

    subroutine burr_rvs_array(state, samples, c, d, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = burr_ppf(0.5_dp, c, d, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = burr_ppf(rng_uniform(state), c, d, mu, sigma)
        end do
    end subroutine burr_rvs_array

    function mielke_rvs(state, k, s, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: k !! positive finite beta-kappa first shape parameter
        real(dp), intent(in) :: s !! positive finite beta-kappa second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = mielke_ppf(0.5_dp, k, s, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = mielke_ppf(rng_uniform(state), k, s, mu, sigma)
    end function mielke_rvs

    subroutine mielke_rvs_array(state, samples, k, s, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: k !! positive finite beta-kappa first shape parameter
        real(dp), intent(in) :: s !! positive finite beta-kappa second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = mielke_ppf(0.5_dp, k, s, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = mielke_ppf(rng_uniform(state), k, s, mu, sigma)
        end do
    end subroutine mielke_rvs_array

    function gibrat_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gibrat_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gibrat_ppf(rng_uniform(state), mu, sigma)
    end function gibrat_rvs

    subroutine gibrat_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gibrat_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = gibrat_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine gibrat_rvs_array

    function wrapcauchy_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = wrapcauchy_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = wrapcauchy_ppf(rng_uniform(state), c, mu, sigma)
    end function wrapcauchy_rvs

    subroutine wrapcauchy_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = wrapcauchy_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = wrapcauchy_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine wrapcauchy_rvs_array

    function genextreme_rvs(state, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genextreme_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = genextreme_ppf(rng_uniform(state), c, mu, sigma)
    end function genextreme_rvs

    subroutine genextreme_rvs_array(state, samples, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genextreme_ppf(0.5_dp, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = genextreme_ppf(rng_uniform(state), c, mu, sigma)
        end do
    end subroutine genextreme_rvs_array

    function kappa3_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive finite kappa-3 shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = kappa3_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = kappa3_ppf(rng_uniform(state), a, mu, sigma)
    end function kappa3_rvs

    subroutine kappa3_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive finite kappa-3 shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = kappa3_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = kappa3_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine kappa3_rvs_array

    function kappa4_rvs(state, h, k, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: h !! finite first kappa-4 shape parameter
        real(dp), intent(in) :: k !! finite second kappa-4 shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = kappa4_ppf(0.5_dp, h, k, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = kappa4_ppf(rng_uniform(state), h, k, mu, sigma)
    end function kappa4_rvs

    subroutine kappa4_rvs_array(state, samples, h, k, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: h !! finite first kappa-4 shape parameter
        real(dp), intent(in) :: k !! finite second kappa-4 shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = kappa4_ppf(0.5_dp, h, k, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = kappa4_ppf(rng_uniform(state), h, k, mu, sigma)
        end do
    end subroutine kappa4_rvs_array

    function truncweibull_min_rvs(state, c, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncweibull_min_ppf(0.5_dp, c, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = truncweibull_min_ppf(rng_uniform(state), c, a, b, mu, sigma)
    end function truncweibull_min_rvs

    subroutine truncweibull_min_rvs_array(state, samples, c, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = truncweibull_min_ppf(0.5_dp, c, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = truncweibull_min_ppf(rng_uniform(state), c, a, b, mu, sigma)
        end do
    end subroutine truncweibull_min_rvs_array


    function gengamma_rvs(state, a, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gengamma_ppf(0.5_dp, a, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = gengamma_ppf(rng_uniform(state), a, c, mu, sigma)
    end function gengamma_rvs

    subroutine gengamma_rvs_array(state, samples, a, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = gengamma_ppf(0.5_dp, a, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = gengamma_ppf(rng_uniform(state), a, c, mu, sigma)
        end do
    end subroutine gengamma_rvs_array

    function halfgennorm_rvs(state, beta, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: beta !! positive finite generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = halfgennorm_ppf(0.5_dp, beta, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = halfgennorm_ppf(rng_uniform(state), beta, mu, sigma)
    end function halfgennorm_rvs

    subroutine halfgennorm_rvs_array(state, samples, beta, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: beta !! positive finite generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = halfgennorm_ppf(0.5_dp, beta, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = halfgennorm_ppf(rng_uniform(state), beta, mu, sigma)
        end do
    end subroutine halfgennorm_rvs_array

    function argus_rvs(state, chi, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = argus_ppf(0.5_dp, chi, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = argus_ppf(rng_uniform(state), chi, mu, sigma)
    end function argus_rvs

    subroutine argus_rvs_array(state, samples, chi, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = argus_ppf(0.5_dp, chi, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = argus_ppf(rng_uniform(state), chi, mu, sigma)
        end do
    end subroutine argus_rvs_array

    function erlang_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = erlang_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = erlang_ppf(rng_uniform(state), a, mu, sigma)
    end function erlang_rvs

    subroutine erlang_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = erlang_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = erlang_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine erlang_rvs_array


    function crystalball_rvs(state, beta, m, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = crystalball_ppf(0.5_dp, beta, m, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = crystalball_ppf(rng_uniform(state), beta, m, mu, sigma)
    end function crystalball_rvs

    subroutine crystalball_rvs_array(state, samples, beta, m, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = crystalball_ppf(0.5_dp, beta, m, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = crystalball_ppf(rng_uniform(state), beta, m, mu, sigma)
        end do
    end subroutine crystalball_rvs_array

    function jf_skew_t_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = jf_skew_t_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = jf_skew_t_ppf(rng_uniform(state), a, b, mu, sigma)
    end function jf_skew_t_rvs

    subroutine jf_skew_t_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = jf_skew_t_ppf(0.5_dp, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = jf_skew_t_ppf(rng_uniform(state), a, b, mu, sigma)
        end do
    end subroutine jf_skew_t_rvs_array

    function pearson3_rvs(state, skew, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = pearson3_ppf(0.5_dp, skew, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = pearson3_ppf(rng_uniform(state), skew, mu, sigma)
    end function pearson3_rvs

    subroutine pearson3_rvs_array(state, samples, skew, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = pearson3_ppf(0.5_dp, skew, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = pearson3_ppf(rng_uniform(state), skew, mu, sigma)
        end do
    end subroutine pearson3_rvs_array

    function rel_breitwigner_rvs(state, rho, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = rel_breitwigner_ppf(0.5_dp, rho, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = rel_breitwigner_ppf(rng_uniform(state), rho, mu, sigma)
    end function rel_breitwigner_rvs

    subroutine rel_breitwigner_rvs_array(state, samples, rho, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = rel_breitwigner_ppf(0.5_dp, rho, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = rel_breitwigner_ppf(rng_uniform(state), rho, mu, sigma)
        end do
    end subroutine rel_breitwigner_rvs_array


    function genexpon_rvs(state, a, b, c, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genexpon_ppf(0.5_dp, a, b, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = genexpon_ppf(rng_uniform(state), a, b, c, mu, sigma)
    end function genexpon_rvs

    subroutine genexpon_rvs_array(state, samples, a, b, c, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = genexpon_ppf(0.5_dp, a, b, c, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = genexpon_ppf(rng_uniform(state), a, b, c, mu, sigma)
        end do
    end subroutine genexpon_rvs_array

    function skewnorm_rvs(state, a, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! finite skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = skewnorm_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = skewnorm_ppf(rng_uniform(state), a, mu, sigma)
    end function skewnorm_rvs

    subroutine skewnorm_rvs_array(state, samples, a, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! finite skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = skewnorm_ppf(0.5_dp, a, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = skewnorm_ppf(rng_uniform(state), a, mu, sigma)
        end do
    end subroutine skewnorm_rvs_array

    function tukeylambda_rvs(state, lam, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = tukeylambda_ppf(0.5_dp, lam, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = tukeylambda_ppf(rng_uniform(state), lam, mu, sigma)
    end function tukeylambda_rvs

    subroutine tukeylambda_rvs_array(state, samples, lam, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = tukeylambda_ppf(0.5_dp, lam, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = tukeylambda_ppf(rng_uniform(state), lam, mu, sigma)
        end do
    end subroutine tukeylambda_rvs_array

    function rice_rvs(state, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = rice_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = rice_ppf(rng_uniform(state), b, mu, sigma)
    end function rice_rvs

    subroutine rice_rvs_array(state, samples, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = rice_ppf(0.5_dp, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = rice_ppf(rng_uniform(state), b, mu, sigma)
        end do
    end subroutine rice_rvs_array


    function dpareto_lognorm_rvs(state, u, s, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = dpareto_lognorm_ppf(0.5_dp, u, s, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = dpareto_lognorm_ppf(rng_uniform(state), u, s, a, b, mu, sigma)
    end function dpareto_lognorm_rvs

    subroutine dpareto_lognorm_rvs_array(state, samples, u, s, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = dpareto_lognorm_ppf(0.5_dp, u, s, a, b, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = dpareto_lognorm_ppf(rng_uniform(state), u, s, a, b, mu, sigma)
        end do
    end subroutine dpareto_lognorm_rvs_array

    function vonmises_rvs(state, kappa, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = vonmises_ppf(0.5_dp, kappa, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = vonmises_ppf(rng_uniform(state), kappa, mu, sigma)
    end function vonmises_rvs

    subroutine vonmises_rvs_array(state, samples, kappa, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = vonmises_ppf(0.5_dp, kappa, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = vonmises_ppf(rng_uniform(state), kappa, mu, sigma)
        end do
    end subroutine vonmises_rvs_array

    function vonmises_line_rvs(state, kappa, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! support center (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = vonmises_line_ppf(0.5_dp, kappa, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = vonmises_line_ppf(rng_uniform(state), kappa, mu, sigma)
    end function vonmises_line_rvs

    subroutine vonmises_line_rvs_array(state, samples, kappa, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! support center (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = vonmises_line_ppf(0.5_dp, kappa, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = vonmises_line_ppf(rng_uniform(state), kappa, mu, sigma)
        end do
    end subroutine vonmises_line_rvs_array

    function kstwobign_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = kstwobign_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = kstwobign_ppf(rng_uniform(state), mu, sigma)
    end function kstwobign_rvs

    subroutine kstwobign_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = kstwobign_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = kstwobign_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine kstwobign_rvs_array

    function irwinhall_rvs(state, n, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = irwinhall_ppf(0.5_dp, n, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = irwinhall_ppf(rng_uniform(state), n, mu, sigma)
    end function irwinhall_rvs

    subroutine irwinhall_rvs_array(state, samples, n, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = irwinhall_ppf(0.5_dp, n, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = irwinhall_ppf(rng_uniform(state), n, mu, sigma)
        end do
    end subroutine irwinhall_rvs_array

    function ksone_rvs(state, n, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: n !! positive integer sample-size shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = ksone_ppf(0.5_dp, n, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = ksone_ppf(rng_uniform(state), n, mu, sigma)
    end function ksone_rvs

    subroutine ksone_rvs_array(state, samples, n, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: n !! positive integer sample-size shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = ksone_ppf(0.5_dp, n, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = ksone_ppf(rng_uniform(state), n, mu, sigma)
        end do
    end subroutine ksone_rvs_array

    function kstwo_rvs(state, n, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced on success
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: x, mu, sigma, probe

        call get_loc_scale(loc, scale, mu, sigma)
        probe = kstwo_ppf(0.5_dp, n, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = kstwo_ppf(rng_uniform(state), n, mu, sigma)
    end function kstwo_rvs

    subroutine kstwo_rvs_array(state, samples, n, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced per generated sample
        real(dp), intent(out) :: samples(:) !! generated finite-sample two-sided KS observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: mu, sigma, probe
        integer :: i

        call get_loc_scale(loc, scale, mu, sigma)
        probe = kstwo_ppf(0.5_dp, n, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = kstwo_ppf(rng_uniform(state), n, mu, sigma)
        end do
    end subroutine kstwo_rvs_array

    function ncx2_rvs(state, df, nc, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = ncx2_ppf(0.5_dp, df, nc, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = ncx2_ppf(rng_uniform(state), df, nc, mu, sigma)
    end function ncx2_rvs

    subroutine ncx2_rvs_array(state, samples, df, nc, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = ncx2_ppf(0.5_dp, df, nc, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = ncx2_ppf(rng_uniform(state), df, nc, mu, sigma)
        end do
    end subroutine ncx2_rvs_array

    function ncf_rvs(state, dfn, dfd, nc, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = ncf_ppf(0.5_dp, dfn, dfd, nc, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = ncf_ppf(rng_uniform(state), dfn, dfd, nc, mu, sigma)
    end function ncf_rvs

    subroutine ncf_rvs_array(state, samples, dfn, dfd, nc, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu, probe, sigma
        call get_loc_scale(loc, scale, mu, sigma)
        probe = ncf_ppf(0.5_dp, dfn, dfd, nc, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = ncf_ppf(rng_uniform(state), dfn, dfd, nc, mu, sigma)
        end do
    end subroutine ncf_rvs_array

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! requested location, absent for zero
        real(dp), intent(in), optional :: scale !! requested scale, absent for one
        real(dp), intent(out) :: mu !! resolved location
        real(dp), intent(out) :: sigma !! resolved scale

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure function get_shift(loc) result(shift)
        real(dp), intent(in), optional :: loc !! requested support shift, absent for zero
        real(dp) :: shift

        shift = 0.0_dp
        if (present(loc)) shift = loc
    end function get_shift


    function randint_rvs(state, low, high, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = randint_ppf(0.5_dp, low, high, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = randint_ppf(rng_uniform(state), low, high, shift)
    end function randint_rvs

    subroutine randint_rvs_array(state, samples, low, high, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = randint_ppf(0.5_dp, low, high, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = randint_ppf(rng_uniform(state), low, high, shift)
        end do
    end subroutine randint_rvs_array

    function planck_rvs(state, lambda, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = planck_ppf(0.5_dp, lambda, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = planck_ppf(rng_uniform(state), lambda, shift)
    end function planck_rvs

    subroutine planck_rvs_array(state, samples, lambda, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = planck_ppf(0.5_dp, lambda, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = planck_ppf(rng_uniform(state), lambda, shift)
        end do
    end subroutine planck_rvs_array

    function dlaplace_rvs(state, a, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = dlaplace_ppf(0.5_dp, a, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = dlaplace_ppf(rng_uniform(state), a, shift)
    end function dlaplace_rvs

    subroutine dlaplace_rvs_array(state, samples, a, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = dlaplace_ppf(0.5_dp, a, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = dlaplace_ppf(rng_uniform(state), a, shift)
        end do
    end subroutine dlaplace_rvs_array

    function logser_rvs(state, p, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = logser_ppf(0.5_dp, p, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = logser_ppf(rng_uniform(state), p, shift)
    end function logser_rvs

    subroutine logser_rvs_array(state, samples, p, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = logser_ppf(0.5_dp, p, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = logser_ppf(rng_uniform(state), p, shift)
        end do
    end subroutine logser_rvs_array


    function betabinom_rvs(state, n, a, b, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = betabinom_ppf(0.5_dp, n, a, b, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = betabinom_ppf(rng_uniform(state), n, a, b, shift)
    end function betabinom_rvs

    subroutine betabinom_rvs_array(state, samples, n, a, b, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = betabinom_ppf(0.5_dp, n, a, b, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = betabinom_ppf(rng_uniform(state), n, a, b, shift)
        end do
    end subroutine betabinom_rvs_array

    function hypergeom_rvs(state, m, n, draws, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer Type-I count
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = hypergeom_ppf(0.5_dp, m, n, draws, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = hypergeom_ppf(rng_uniform(state), m, n, draws, shift)
    end function hypergeom_rvs

    subroutine hypergeom_rvs_array(state, samples, m, n, draws, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer Type-I count
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = hypergeom_ppf(0.5_dp, m, n, draws, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = hypergeom_ppf(rng_uniform(state), m, n, draws, shift)
        end do
    end subroutine hypergeom_rvs_array

    function nhypergeom_rvs(state, m, n, r, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer success count
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = nhypergeom_ppf(0.5_dp, m, n, r, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = nhypergeom_ppf(rng_uniform(state), m, n, r, shift)
    end function nhypergeom_rvs

    subroutine nhypergeom_rvs_array(state, samples, m, n, r, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer success count
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = nhypergeom_ppf(0.5_dp, m, n, r, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = nhypergeom_ppf(rng_uniform(state), m, n, r, shift)
        end do
    end subroutine nhypergeom_rvs_array

    function boltzmann_rvs(state, lambda, n, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = boltzmann_ppf(0.5_dp, lambda, n, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = boltzmann_ppf(rng_uniform(state), lambda, n, shift)
    end function boltzmann_rvs

    subroutine boltzmann_rvs_array(state, samples, lambda, n, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = boltzmann_ppf(0.5_dp, lambda, n, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = boltzmann_ppf(rng_uniform(state), lambda, n, shift)
        end do
    end subroutine boltzmann_rvs_array


    function betanbinom_rvs(state, n, a, b, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = betanbinom_ppf(0.5_dp, n, a, b, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = betanbinom_ppf(rng_uniform(state), n, a, b, shift)
    end function betanbinom_rvs

    subroutine betanbinom_rvs_array(state, samples, n, a, b, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = betanbinom_ppf(0.5_dp, n, a, b, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = betanbinom_ppf(rng_uniform(state), n, a, b, shift)
        end do
    end subroutine betanbinom_rvs_array

    function yulesimon_rvs(state, alpha, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = yulesimon_ppf(0.5_dp, alpha, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = yulesimon_ppf(rng_uniform(state), alpha, shift)
    end function yulesimon_rvs

    subroutine yulesimon_rvs_array(state, samples, alpha, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = yulesimon_ppf(0.5_dp, alpha, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = yulesimon_ppf(rng_uniform(state), alpha, shift)
        end do
    end subroutine yulesimon_rvs_array

    function zipf_rvs(state, a, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = zipf_ppf(0.5_dp, a, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = zipf_ppf(rng_uniform(state), a, shift)
    end function zipf_rvs

    subroutine zipf_rvs_array(state, samples, a, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = zipf_ppf(0.5_dp, a, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = zipf_ppf(rng_uniform(state), a, shift)
        end do
    end subroutine zipf_rvs_array

    function zipfian_rvs(state, a, n, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = zipfian_ppf(0.5_dp, a, n, shift)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = zipfian_ppf(rng_uniform(state), a, n, shift)
    end function zipfian_rvs

    subroutine zipfian_rvs_array(state, samples, a, n, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        real(dp) :: shift, probe
        shift = get_shift(loc)
        probe = zipfian_ppf(0.5_dp, a, n, shift)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = zipfian_ppf(rng_uniform(state), a, n, shift)
        end do
    end subroutine zipfian_rvs_array

    function geninvgauss_rvs(state, p, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: p !! real generalized-inverse-Gaussian shape
        real(dp), intent(in) :: b !! strictly positive generalized-inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, probe
        call get_loc_scale(loc, scale, mu, sigma)
        probe=geninvgauss_ppf(0.5_dp,p,b,mu,sigma)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=geninvgauss_ppf(rng_uniform(state),p,b,mu,sigma)
    end function geninvgauss_rvs

    subroutine geninvgauss_rvs_array(state, samples, p, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: p !! real generalized-inverse-Gaussian shape
        real(dp), intent(in) :: b !! strictly positive generalized-inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu,sigma,probe
        integer :: i
        call get_loc_scale(loc,scale,mu,sigma)
        probe=geninvgauss_ppf(0.5_dp,p,b,mu,sigma)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples); samples(i)=geninvgauss_ppf(rng_uniform(state),p,b,mu,sigma); end do
    end subroutine geninvgauss_rvs_array

    function norminvgauss_rvs(state, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive normal-inverse-Gaussian tail shape
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x,mu,sigma,probe
        call get_loc_scale(loc,scale,mu,sigma)
        probe=norminvgauss_ppf(0.5_dp,a,b,mu,sigma)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=norminvgauss_ppf(rng_uniform(state),a,b,mu,sigma)
    end function norminvgauss_rvs

    subroutine norminvgauss_rvs_array(state, samples, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive normal-inverse-Gaussian tail shape
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu,sigma,probe
        integer :: i
        call get_loc_scale(loc,scale,mu,sigma)
        probe=norminvgauss_ppf(0.5_dp,a,b,mu,sigma)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples); samples(i)=norminvgauss_ppf(rng_uniform(state),a,b,mu,sigma); end do
    end subroutine norminvgauss_rvs_array

    function skellam_rvs(state, mu1, mu2, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: x,shift,probe
        shift=get_shift(loc); probe=skellam_ppf(0.5_dp,mu1,mu2,shift)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=skellam_ppf(rng_uniform(state),mu1,mu2,shift)
    end function skellam_rvs

    subroutine skellam_rvs_array(state, samples, mu1, mu2, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: shift,probe
        integer :: i
        shift=get_shift(loc); probe=skellam_ppf(0.5_dp,mu1,mu2,shift)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples); samples(i)=skellam_ppf(rng_uniform(state),mu1,mu2,shift); end do
    end subroutine skellam_rvs_array


    function genhyperbolic_rvs(state, p, a, b, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: p !! real generalized-hyperbolic tail shape
        real(dp), intent(in) :: a !! positive generalized-hyperbolic shape
        real(dp), intent(in) :: b !! generalized-hyperbolic skew shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, probe
        call get_loc_scale(loc,scale,mu,sigma)
        probe=genhyperbolic_ppf(0.5_dp,p,a,b,mu,sigma)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=genhyperbolic_ppf(rng_uniform(state),p,a,b,mu,sigma)
    end function genhyperbolic_rvs

    subroutine genhyperbolic_rvs_array(state, samples, p, a, b, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: p !! real generalized-hyperbolic tail shape
        real(dp), intent(in) :: a !! positive generalized-hyperbolic shape
        real(dp), intent(in) :: b !! generalized-hyperbolic skew shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu, sigma, probe
        integer :: i
        call get_loc_scale(loc,scale,mu,sigma)
        probe=genhyperbolic_ppf(0.5_dp,p,a,b,mu,sigma)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples); samples(i)=genhyperbolic_ppf(rng_uniform(state),p,a,b,mu,sigma); end do
    end subroutine genhyperbolic_rvs_array

    function nchypergeom_fisher_rvs(state, m, n, draws, odds, loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x, shift, probe
        shift=get_shift(loc); probe=nchypergeom_fisher_ppf(0.5_dp,m,n,draws,odds,shift)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=nchypergeom_fisher_ppf(rng_uniform(state),m,n,draws,odds,shift)
    end function nchypergeom_fisher_rvs

    subroutine nchypergeom_fisher_rvs_array(state, samples, m, n, draws, odds, loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: shift, probe
        integer :: i
        shift=get_shift(loc); probe=nchypergeom_fisher_ppf(0.5_dp,m,n,draws,odds,shift)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples)
            samples(i)=nchypergeom_fisher_ppf(rng_uniform(state),m,n,draws,odds,shift)
        end do
    end subroutine nchypergeom_fisher_rvs_array

    function nct_rvs(state, df, nc, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: df !! positive noncentral-t degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, probe
        call get_loc_scale(loc,scale,mu,sigma)
        probe=nct_ppf(0.5_dp,df,nc,mu,sigma)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=nct_ppf(rng_uniform(state),df,nc,mu,sigma)
    end function nct_rvs

    subroutine nct_rvs_array(state, samples, df, nc, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: df !! positive noncentral-t degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu, sigma, probe
        integer :: i
        call get_loc_scale(loc,scale,mu,sigma)
        probe=nct_ppf(0.5_dp,df,nc,mu,sigma)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples); samples(i)=nct_ppf(rng_uniform(state),df,nc,mu,sigma); end do
    end subroutine nct_rvs_array

    function gausshyper_rvs(state, a, b, c, z, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, probe
        call get_loc_scale(loc,scale,mu,sigma)
        probe=gausshyper_ppf(0.5_dp,a,b,c,z,mu,sigma)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        x=gausshyper_ppf(rng_uniform(state),a,b,c,z,mu,sigma)
    end function gausshyper_rvs

    subroutine gausshyper_rvs_array(state, samples, a, b, c, z, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: z !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu, sigma, probe
        integer :: i
        call get_loc_scale(loc,scale,mu,sigma)
        probe=gausshyper_ppf(0.5_dp,a,b,c,z,mu,sigma)
        if (ieee_is_nan(probe)) then; samples=probe; return; end if
        do i=1,size(samples); samples(i)=gausshyper_ppf(rng_uniform(state),a,b,c,z,mu,sigma); end do
    end subroutine gausshyper_rvs_array

    function landau_rvs(state, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = landau_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = landau_ppf(rng_uniform(state), mu, sigma)
    end function landau_rvs

    subroutine landau_rvs_array(state, samples, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        integer :: i
        real(dp) :: mu
        real(dp) :: probe
        real(dp) :: sigma

        call get_loc_scale(loc, scale, mu, sigma)
        probe = landau_ppf(0.5_dp, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = landau_ppf(rng_uniform(state), mu, sigma)
        end do
    end subroutine landau_rvs_array

    function nchypergeom_wallenius_rvs(state,m,n,draws,odds,loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by sequential draws
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x, shift, probe, q, den
        integer :: j, aleft, bleft, count
        shift=get_shift(loc); probe=nchypergeom_wallenius_ppf(0.5_dp,m,n,draws,odds,shift)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        aleft=int(n); bleft=int(m-n); count=0
        do j=1,int(draws)
            den=odds*real(aleft,dp)+real(bleft,dp)
            if (aleft<=0) then
                bleft=bleft-1
            else if (bleft<=0) then
                count=count+1; aleft=aleft-1
            else
                q=odds*real(aleft,dp)/den
                if (rng_uniform(state)<q) then
                    count=count+1; aleft=aleft-1
                else
                    bleft=bleft-1
                end if
            end if
        end do
        x=shift+real(count,dp)
    end function nchypergeom_wallenius_rvs

    subroutine nchypergeom_wallenius_rvs_array(state,samples,m,n,draws,odds,loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by sequential draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        do i=1,size(samples)
            samples(i)=nchypergeom_wallenius_rvs(state,m,n,draws,odds,loc)
        end do
    end subroutine nchypergeom_wallenius_rvs_array

    function poisson_binom_rvs(state,p,loc) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by Bernoulli draws
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: x, shift, probe
        integer :: i, count
        shift=get_shift(loc); probe=poisson_binom_ppf(0.5_dp,p,shift)
        if (ieee_is_nan(probe)) then; x=probe; return; end if
        count=0
        do i=1,size(p)
            if (rng_uniform(state)<p(i)) count=count+1
        end do
        x=shift+real(count,dp)
    end function poisson_binom_rvs

    subroutine poisson_binom_rvs_array(state,samples,p,loc)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by Bernoulli draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities in [0,1]
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        integer :: i
        do i=1,size(samples); samples(i)=poisson_binom_rvs(state,p,loc); end do
    end subroutine poisson_binom_rvs_array

    function levy_stable_rvs(state, alpha, beta, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, probe

        call get_loc_scale(loc, scale, mu, sigma)
        probe = levy_stable_ppf(0.5_dp, alpha, beta, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = levy_stable_ppf(rng_uniform(state), alpha, beta, mu, sigma)
    end function levy_stable_rvs

    subroutine levy_stable_rvs_array(state, samples, alpha, beta, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu, sigma, probe
        integer :: i

        call get_loc_scale(loc, scale, mu, sigma)
        probe = levy_stable_ppf(0.5_dp, alpha, beta, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = levy_stable_ppf(rng_uniform(state), alpha, beta, mu, sigma)
        end do
    end subroutine levy_stable_rvs_array

    function studentized_range_rvs(state, k, df, loc, scale) result(x)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by one valid draw
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, mu, sigma, probe

        call get_loc_scale(loc, scale, mu, sigma)
        probe = studentized_range_ppf(0.5_dp, k, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            x = probe
            return
        end if
        x = studentized_range_ppf(rng_uniform(state), k, df, mu, sigma)
    end function studentized_range_rvs

    subroutine studentized_range_rvs_array(state, samples, k, df, loc, scale)
        type(rng_state), intent(inout) :: state !! explicit RNG state advanced by valid draws
        real(dp), intent(out) :: samples(:) !! generated random variates
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: mu, sigma, probe
        integer :: i

        call get_loc_scale(loc, scale, mu, sigma)
        probe = studentized_range_ppf(0.5_dp, k, df, mu, sigma)
        if (ieee_is_nan(probe)) then
            samples = probe
            return
        end if
        do i = 1, size(samples)
            samples(i) = studentized_range_ppf(rng_uniform(state), k, df, mu, sigma)
        end do
    end subroutine studentized_range_rvs_array

end module scifort_random_variates
