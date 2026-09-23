! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Score (log-likelihood gradient) functions for the implemented probability
! families. Continuous-family result ordering follows shape parameters first,
! then loc and scale. Discrete families score only differentiable distribution
! parameters: integer support shifts and binomial n are intentionally omitted.
!
! Moving-support families have nonregular location/scale derivatives at support
! boundaries. The affected score components return NaN there rather than
! reporting a misleading one-sided derivative.

module scifort_score
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_inv_sqrt_two_pi, scifort_log_sqrt_two_pi, &
        scifort_log_two, scifort_pi
    use scifort_incomplete_gamma, only : gammainc
    use scifort_digamma, only : digamma_positive
    use scifort_kinds, only : dp
    use scifort_likelihood, only : normal_loglikelihood, uniform_loglikelihood, &
        exponential_loglikelihood, laplace_loglikelihood, logistic_loglikelihood, &
        cauchy_loglikelihood, rayleigh_loglikelihood, gamma_loglikelihood, &
        chi2_loglikelihood, t_loglikelihood, lognormal_loglikelihood, &
        weibull_loglikelihood, pareto_loglikelihood, beta_loglikelihood, &
        f_loglikelihood, gumbel_r_loglikelihood, gumbel_l_loglikelihood, &
        powerlaw_loglikelihood, triang_loglikelihood, genpareto_loglikelihood, &
        arcsine_loglikelihood, halfnorm_loglikelihood, halfcauchy_loglikelihood, &
        lomax_loglikelihood, chi_loglikelihood, maxwell_loglikelihood, &
        cosine_loglikelihood, semicircular_loglikelihood, &
        anglit_loglikelihood, moyal_loglikelihood, &
        hypsecant_loglikelihood, halflogistic_loglikelihood, &
        invgamma_loglikelihood, invgauss_loglikelihood, levy_loglikelihood, &
        loglaplace_loglikelihood, bradford_loglikelihood, truncexpon_loglikelihood, &
        fisk_loglikelihood, dweibull_loglikelihood, alpha_loglikelihood, &
        fatiguelife_loglikelihood, genlogistic_loglikelihood, gennorm_loglikelihood, &
        nakagami_loglikelihood, powernorm_loglikelihood, loggamma_loglikelihood, &
        wald_loglikelihood, gompertz_loglikelihood, invweibull_loglikelihood, &
        betaprime_loglikelihood, burr12_loglikelihood, genhalflogistic_loglikelihood, &
        exponpow_loglikelihood, exponweib_loglikelihood, powerlognorm_loglikelihood, &
        levy_l_loglikelihood, weibull_max_loglikelihood, rdist_loglikelihood, &
        skewcauchy_loglikelihood, dgamma_loglikelihood, &
        laplace_asymmetric_loglikelihood, truncnorm_loglikelihood, &
        loguniform_loglikelihood, foldnorm_loglikelihood, &
        foldcauchy_loglikelihood, recipinvgauss_loglikelihood, &
        truncpareto_loglikelihood, exponnorm_loglikelihood, &
        johnsonsb_loglikelihood, johnsonsu_loglikelihood, trapezoid_loglikelihood, &
        burr_loglikelihood, mielke_loglikelihood, gibrat_loglikelihood, &
        wrapcauchy_loglikelihood, genextreme_loglikelihood, kappa3_loglikelihood, &
        kappa4_loglikelihood, truncweibull_min_loglikelihood, &
        gengamma_loglikelihood, halfgennorm_loglikelihood, argus_loglikelihood, &
        erlang_loglikelihood, crystalball_loglikelihood, jf_skew_t_loglikelihood, &
        pearson3_loglikelihood, rel_breitwigner_loglikelihood, &
        genexpon_loglikelihood, skewnorm_loglikelihood, tukeylambda_loglikelihood, &
        rice_loglikelihood, dpareto_lognorm_loglikelihood, &
        vonmises_loglikelihood, vonmises_line_loglikelihood, kstwobign_loglikelihood, &
        irwinhall_loglikelihood, ksone_loglikelihood, ncx2_loglikelihood, &
        ncf_loglikelihood, bernoulli_loglikelihood, &
        poisson_loglikelihood, &
        geometric_loglikelihood, binomial_loglikelihood, &
        negative_binomial_loglikelihood, randint_loglikelihood, planck_loglikelihood, &
        dlaplace_loglikelihood, logser_loglikelihood, &
        betabinom_loglikelihood, hypergeom_loglikelihood, nhypergeom_loglikelihood, &
        boltzmann_loglikelihood, betanbinom_loglikelihood, yulesimon_loglikelihood, &
        zipf_loglikelihood, zipfian_loglikelihood, geninvgauss_loglikelihood, &
        norminvgauss_loglikelihood, skellam_loglikelihood, genhyperbolic_loglikelihood, &
        nchypergeom_fisher_loglikelihood, nct_loglikelihood, gausshyper_loglikelihood, &
        landau_loglikelihood, nchypergeom_wallenius_loglikelihood, poisson_binom_loglikelihood, &
        kstwo_loglikelihood, levy_stable_loglikelihood, &
        studentized_range_loglikelihood
    use scifort_fisk, only : fisk_cdf
    use scifort_truncnorm, only : truncnorm_pdf
    use scifort_normal, only : normal_cdf, normal_logcdf, normal_logpdf, normal_logsf
    use scifort_math, only : expm1_safe, log1p_safe, log1pexp, quiet_nan
    use scifort_tukeylambda, only : tukeylambda_cdf
    use scifort_rice, only : rice_i1_i0_ratio
    use scifort_vonmises, only : vonmises_i1_i0_ratio
    use scifort_kstwobign, only : kstwobign_logpdf_derivative
    use scifort_irwinhall, only : irwinhall_logpdf_derivative
    use scifort_ksone, only : ksone_logpdf_derivative
    use scifort_kstwo, only : kstwo_logpdf_derivative
    use scifort_ncx2, only : ncx2_logpdf_derivatives
    use scifort_ncf, only : ncf_logpdf_derivatives
    use scifort_zeta, only : hurwitz_zeta, hurwitz_zeta_derivative
    use scifort_geninvgauss, only : geninvgauss_logpdf_derivatives
    use scifort_norminvgauss, only : norminvgauss_logpdf_derivatives
    use scifort_skellam, only : skellam_logpmf_derivatives
    use scifort_genhyperbolic, only : genhyperbolic_logpdf_derivatives
    use scifort_nchypergeom_fisher, only : nchypergeom_fisher_logpmf_derivative_odds
    use scifort_nchypergeom_wallenius, only : nchypergeom_wallenius_logpmf_derivative_odds
    use scifort_poisson_binom, only : poisson_binom_logpmf_derivative_p
    use scifort_nct, only : nct_logpdf_derivatives
    use scifort_gausshyper, only : gausshyper_logpdf_derivatives
    use scifort_landau, only : landau_logpdf_derivatives
    implicit none
    private

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
    public :: landau_score
    public :: hypsecant_score
    public :: halflogistic_score
    public :: invgamma_score
    public :: invgauss_score
    public :: levy_score
    public :: loglaplace_score
    public :: bradford_score
    public :: truncexpon_score
    public :: fisk_score
    public :: dweibull_score
    public :: alpha_score
    public :: fatiguelife_score
    public :: genlogistic_score
    public :: gennorm_score
    public :: nakagami_score
    public :: powernorm_score
    public :: loggamma_score
    public :: wald_score
    public :: gompertz_score
    public :: invweibull_score
    public :: betaprime_score
    public :: burr12_score
    public :: genhalflogistic_score
    public :: exponpow_score
    public :: exponweib_score
    public :: powerlognorm_score
    public :: levy_l_score
    public :: weibull_max_score
    public :: rdist_score
    public :: skewcauchy_score
    public :: dgamma_score
    public :: laplace_asymmetric_score
    public :: truncnorm_score
    public :: loguniform_score
    public :: bernoulli_score
    public :: poisson_score
    public :: geometric_score
    public :: binomial_score
    public :: negative_binomial_score
    public :: foldnorm_score
    public :: foldcauchy_score
    public :: recipinvgauss_score
    public :: truncpareto_score
    public :: exponnorm_score
    public :: johnsonsb_score
    public :: johnsonsu_score
    public :: trapezoid_score
    public :: burr_score
    public :: mielke_score
    public :: gibrat_score
    public :: wrapcauchy_score
    public :: genextreme_score
    public :: kappa3_score
    public :: kappa4_score
    public :: truncweibull_min_score
    public :: gengamma_score
    public :: halfgennorm_score
    public :: argus_score
    public :: erlang_score
    public :: crystalball_score
    public :: jf_skew_t_score
    public :: pearson3_score
    public :: rel_breitwigner_score
    public :: genexpon_score
    public :: skewnorm_score
    public :: tukeylambda_score
    public :: rice_score
    public :: dpareto_lognorm_score
    public :: vonmises_score
    public :: vonmises_line_score
    public :: kstwobign_score
    public :: irwinhall_score
    public :: ksone_score
    public :: kstwo_score
    public :: levy_stable_score
    public :: studentized_range_score
    public :: ncx2_score
    public :: ncf_score
    public :: randint_score
    public :: planck_score
    public :: dlaplace_score
    public :: logser_score
    public :: betabinom_score
    public :: hypergeom_score
    public :: nhypergeom_score
    public :: boltzmann_score
    public :: betanbinom_score
    public :: yulesimon_score
    public :: zipf_score
    public :: zipfian_score
    public :: geninvgauss_score
    public :: norminvgauss_score
    public :: skellam_score
    public :: genhyperbolic_score
    public :: nchypergeom_fisher_score
    public :: nct_score
    public :: gausshyper_score
    public :: nchypergeom_wallenius_score
    public :: poisson_binom_score

    interface bernoulli_score
        module procedure bernoulli_score_real
        module procedure bernoulli_score_int
    end interface bernoulli_score

    interface poisson_score
        module procedure poisson_score_real
        module procedure poisson_score_int
    end interface poisson_score

    interface geometric_score
        module procedure geometric_score_real
        module procedure geometric_score_int
    end interface geometric_score

    interface binomial_score
        module procedure binomial_score_real
        module procedure binomial_score_int
    end interface binomial_score

    interface negative_binomial_score
        module procedure negative_binomial_score_real
        module procedure negative_binomial_score_int
    end interface negative_binomial_score

    interface randint_score
        module procedure randint_score_real
        module procedure randint_score_int
    end interface randint_score

    interface planck_score
        module procedure planck_score_real
        module procedure planck_score_int
    end interface planck_score

    interface dlaplace_score
        module procedure dlaplace_score_real
        module procedure dlaplace_score_int
    end interface dlaplace_score

    interface logser_score
        module procedure logser_score_real
        module procedure logser_score_int
    end interface logser_score

    interface betabinom_score
        module procedure betabinom_score_real
        module procedure betabinom_score_int
    end interface betabinom_score

    interface hypergeom_score
        module procedure hypergeom_score_real
        module procedure hypergeom_score_int
    end interface hypergeom_score

    interface nhypergeom_score
        module procedure nhypergeom_score_real
        module procedure nhypergeom_score_int
    end interface nhypergeom_score

    interface boltzmann_score
        module procedure boltzmann_score_real
        module procedure boltzmann_score_int
    end interface boltzmann_score

    interface betanbinom_score
        module procedure betanbinom_score_real
        module procedure betanbinom_score_int
    end interface betanbinom_score

    interface yulesimon_score
        module procedure yulesimon_score_real
        module procedure yulesimon_score_int
    end interface yulesimon_score

    interface zipf_score
        module procedure zipf_score_real
        module procedure zipf_score_int
    end interface zipf_score

    interface zipfian_score
        module procedure zipfian_score_real
        module procedure zipfian_score_int
    end interface zipfian_score

    interface skellam_score
        module procedure skellam_score_real
        module procedure skellam_score_int
    end interface skellam_score

    interface nchypergeom_fisher_score
        module procedure nchypergeom_fisher_score_real
        module procedure nchypergeom_fisher_score_int
    end interface nchypergeom_fisher_score

    interface nchypergeom_wallenius_score
        module procedure nchypergeom_wallenius_score_real
        module procedure nchypergeom_wallenius_score_int
    end interface nchypergeom_wallenius_score

    interface poisson_binom_score
        module procedure poisson_binom_score_real
        module procedure poisson_binom_score_int
    end interface poisson_binom_score

contains

    pure function normal_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(normal_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(z) / sigma
        score(2) = sum(z * z - 1.0_dp) / sigma
    end function normal_score

    pure function uniform_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma

        if (.not. valid_loglike(uniform_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        if (any(data <= mu) .or. any(data >= mu + sigma)) then
            score = nan_value()
            return
        end if
        score(1) = 0.0_dp
        score(2) = -real(size(data), dp) / sigma
    end function uniform_score

    pure function exponential_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! inverse rate, > 0 (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(exponential_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = real(size(data), dp) / sigma
        score(2) = sum(z - 1.0_dp) / sigma
        if (any(data == mu)) score(1) = nan_value()
    end function exponential_score

    pure function laplace_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! mean absolute deviation, > 0 (default 1)
        real(dp) :: score(2)

        integer :: i
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(laplace_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = 0.0_dp
        do i = 1, size(data)
            if (data(i) > mu) score(1) = score(1) + 1.0_dp / sigma
            if (data(i) < mu) score(1) = score(1) - 1.0_dp / sigma
        end do
        if (any(data == mu)) score(1) = nan_value()
        score(2) = sum(abs(z) - 1.0_dp) / sigma
    end function laplace_score

    pure function logistic_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: g(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(logistic_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        g = tanh(0.5_dp * z)
        score(1) = sum(g) / sigma
        score(2) = sum(z * g - 1.0_dp) / sigma
    end function logistic_score

    pure function cauchy_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: denom(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(cauchy_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        denom = 1.0_dp + z * z
        score(1) = sum(2.0_dp * z / denom) / sigma
        score(2) = sum((z * z - 1.0_dp) / denom) / sigma
    end function cauchy_score

    pure function rayleigh_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(rayleigh_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum(z - 1.0_dp / z) / sigma
        score(2) = sum(z * z - 2.0_dp) / sigma
    end function rayleigh_score

    pure function gamma_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: psi_a
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(gamma_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        psi_a = digamma_positive(a)
        score(1) = sum(log(z) - psi_a)
        score(2) = sum(1.0_dp - (a - 1.0_dp) / z) / sigma
        score(3) = sum(z - a) / sigma
    end function gamma_score

    pure function chi2_score(data, df, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: a
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(chi2_loglikelihood(data, df, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        a = 0.5_dp * df
        score(1) = 0.5_dp * sum(log(z) - scifort_log_two - digamma_positive(a))
        score(2) = sum(0.5_dp - (a - 1.0_dp) / z) / sigma
        score(3) = sum(0.5_dp * z - a) / sigma
    end function chi2_score

    pure function t_score(data, df, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: denom(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))
        real(dp) :: z2(size(data))

        if (.not. valid_loglike(t_loglikelihood(data, df, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        z2 = z * z
        denom = df + z2
        score(1) = 0.5_dp * real(size(data), dp) * ( &
            digamma_positive(0.5_dp * (df + 1.0_dp)) - &
            digamma_positive(0.5_dp * df) - 1.0_dp / df) - &
            0.5_dp * sum(log(1.0_dp + z2 / df)) + &
            0.5_dp * (df + 1.0_dp) / df * sum(z2 / denom)
        score(2) = (df + 1.0_dp) * sum(z / denom) / sigma
        score(3) = sum(-1.0_dp + (df + 1.0_dp) * z2 / denom) / sigma
    end function t_score

    pure function lognormal_score(data, s, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: s !! log-space shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: logz(size(data))
        real(dp) :: mu
        real(dp) :: s2
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(lognormal_loglikelihood(data, s, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        logz = log(z)
        s2 = s * s
        score(1) = sum(-1.0_dp / s + logz * logz / (s * s2))
        score(2) = sum((1.0_dp + logz / s2) / z) / sigma
        score(3) = sum(logz / s2) / sigma
    end function lognormal_score

    pure function weibull_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: logz(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))
        real(dp) :: zc(size(data))

        if (.not. valid_loglike(weibull_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        logz = log(z)
        zc = exp(c * logz)
        score(1) = sum(1.0_dp / c + logz * (1.0_dp - zc))
        score(2) = sum(c * zc / z - (c - 1.0_dp) / z) / sigma
        score(3) = c * sum(zc - 1.0_dp) / sigma
    end function weibull_score

    pure function pareto_score(data, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! shape parameter, > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(pareto_loglikelihood(data, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(1.0_dp / b - log(z))
        score(2) = (b + 1.0_dp) * sum(1.0_dp / z) / sigma
        score(3) = real(size(data), dp) * b / sigma
        if (any(z == 1.0_dp)) then
            score(2) = nan_value()
            score(3) = nan_value()
        end if
    end function pareto_score

    pure function beta_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! first shape parameter, > 0
        real(dp), intent(in) :: b !! second shape parameter, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: score(4)

        real(dp) :: mu
        real(dp) :: psi_ab
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(beta_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp) .or. any(z >= 1.0_dp)) then
            score = nan_value()
            return
        end if
        psi_ab = digamma_positive(a + b)
        score(1) = sum(log(z) - digamma_positive(a) + psi_ab)
        score(2) = sum(log(1.0_dp - z) - digamma_positive(b) + psi_ab)
        score(3) = sum((b - 1.0_dp) / (1.0_dp - z) - (a - 1.0_dp) / z) / sigma
        score(4) = sum(-a + z * (b - 1.0_dp) / (1.0_dp - z)) / sigma
    end function beta_score

    pure function f_score(data, dfn, dfd, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: dfn !! numerator degrees of freedom, > 0
        real(dp), intent(in) :: dfd !! denominator degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(4)

        real(dp) :: a
        real(dp) :: ab
        real(dp) :: b
        real(dp) :: denom(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(f_loglikelihood(data, dfn, dfd, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        a = 0.5_dp * dfn
        b = 0.5_dp * dfd
        ab = a + b
        denom = dfd + dfn * z
        score(1) = 0.5_dp * sum(log(dfn) + 1.0_dp + log(z) - log(denom) - &
            digamma_positive(a) + digamma_positive(ab)) - ab * sum(z / denom)
        score(2) = 0.5_dp * sum(log(dfd) + 1.0_dp - log(denom) - &
            digamma_positive(b) + digamma_positive(ab)) - ab * sum(1.0_dp / denom)
        score(3) = sum(ab * dfn / denom - (a - 1.0_dp) / z) / sigma
        score(4) = sum(-a + ab * dfn * z / denom) / sigma
    end function f_score

    pure function gumbel_r_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: e(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(gumbel_r_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        e = exp(-z)
        score(1) = sum(1.0_dp - e) / sigma
        score(2) = sum(-1.0_dp + z * (1.0_dp - e)) / sigma
    end function gumbel_r_score

    pure function gumbel_l_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: e(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(gumbel_l_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        e = exp(z)
        score(1) = sum(e - 1.0_dp) / sigma
        score(2) = sum(-1.0_dp - z + z * e) / sigma
    end function gumbel_l_score

    pure function powerlaw_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! finite positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(powerlaw_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp) .or. any(z >= 1.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = real(size(data), dp) / a + sum(log(z))
        score(2) = -(a - 1.0_dp) * sum(1.0_dp / z) / sigma
        score(3) = -a * real(size(data), dp) / sigma
    end function powerlaw_score

    pure function triang_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! standardized mode, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: score(3)

        integer :: i
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(triang_loglikelihood(data, c, loc, scale)) .or. &
            c <= 0.0_dp .or. c >= 1.0_dp) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp) .or. any(z >= 1.0_dp) .or. any(z == c)) then
            score = nan_value()
            return
        end if

        score = 0.0_dp
        do i = 1, size(data)
            if (z(i) < c) then
                score(1) = score(1) - 1.0_dp / c
                score(2) = score(2) - 1.0_dp / (sigma * z(i))
                score(3) = score(3) - 2.0_dp / sigma
            else
                score(1) = score(1) + 1.0_dp / (1.0_dp - c)
                score(2) = score(2) + 1.0_dp / (sigma * (1.0_dp - z(i)))
                score(3) = score(3) + &
                    (2.0_dp * z(i) - 1.0_dp) / (sigma * (1.0_dp - z(i)))
            end if
        end do
    end function triang_score

    pure function genpareto_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! finite shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        integer :: i
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t(size(data))
        real(dp) :: z(size(data))

        if (.not. valid_loglike(genpareto_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        t = 1.0_dp + c * z
        if (any(z <= 0.0_dp) .or. any(t <= 0.0_dp)) then
            score = nan_value()
            return
        end if

        score = 0.0_dp
        do i = 1, size(data)
            score(1) = score(1) + genpareto_shape_term(c, z(i))
        end do
        score(2) = sum((c + 1.0_dp) / t) / sigma
        score(3) = sum((z - 1.0_dp) / t) / sigma
    end function genpareto_score

    pure function arcsine_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(arcsine_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp) .or. any(z >= 1.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum((1.0_dp - 2.0_dp * z) / (2.0_dp * z * (1.0_dp - z))) / sigma
        score(2) = -sum(1.0_dp / (1.0_dp - z)) / (2.0_dp * sigma)
    end function arcsine_score

    pure function halfnorm_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(halfnorm_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(z) / sigma
        if (any(data == mu)) score(1) = nan_value()
        score(2) = sum(z * z - 1.0_dp) / sigma
    end function halfnorm_score

    pure function halfcauchy_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(halfcauchy_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(2.0_dp * z / (1.0_dp + z * z)) / sigma
        if (any(data == mu)) score(1) = nan_value()
        score(2) = sum((z * z - 1.0_dp) / (1.0_dp + z * z)) / sigma
    end function halfcauchy_score

    pure function lomax_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(lomax_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = real(size(data), dp) / c - sum(log1p_safe(z))
        score(2) = sum((c + 1.0_dp) / (1.0_dp + z)) / sigma
        if (any(data == mu)) score(2) = nan_value()
        score(3) = sum((c * z - 1.0_dp) / (1.0_dp + z)) / sigma
    end function lomax_score

    pure function chi_score(data, df, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! degrees of freedom, finite and > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(chi_loglikelihood(data, df, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum(log(z) - 0.5_dp * scifort_log_two - &
            0.5_dp * digamma_positive(0.5_dp * df))
        score(2) = sum(z - (df - 1.0_dp) / z) / sigma
        score(3) = sum(z * z - df) / sigma
    end function chi_score

    pure function maxwell_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(maxwell_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum(z - 2.0_dp / z) / sigma
        score(2) = sum(z * z - 3.0_dp) / sigma
    end function maxwell_score

    pure function cosine_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t(size(data))
        real(dp) :: z(size(data))

        if (.not. valid_loglike(cosine_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        t = tan(0.5_dp * z)
        score(1) = sum(t) / sigma
        score(2) = sum(z * t - 1.0_dp) / sigma
    end function cosine_score

    pure function semicircular_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: denominator(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(semicircular_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        denominator = 1.0_dp - z * z
        score(1) = sum(z / denominator) / sigma
        score(2) = sum(z * z / denominator - 1.0_dp) / sigma
    end function semicircular_score

    pure function anglit_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t(size(data))
        real(dp) :: z(size(data))

        if (.not. valid_loglike(anglit_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        t = tan(2.0_dp * z)
        score(1) = 2.0_dp * sum(t) / sigma
        score(2) = sum(2.0_dp * z * t - 1.0_dp) / sigma
    end function anglit_score

    pure function moyal_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: e(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(moyal_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        e = exp(-z)
        score(1) = 0.5_dp * sum(1.0_dp - e) / sigma
        score(2) = sum(-1.0_dp + 0.5_dp * z * (1.0_dp - e)) / sigma
    end function moyal_score

    pure function landau_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(2)

        integer :: i
        real(dp) :: dz
        real(dp) :: logf
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        if (.not. valid_loglike(landau_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            call landau_logpdf_derivatives(z, logf, dz)
            if (.not. valid_loglike(logf)) then
                score = nan_value()
                return
            end if
            score(1) = score(1) - dz / sigma
            score(2) = score(2) - (1.0_dp + z * dz) / sigma
        end do
    end function landau_score

    pure function hypsecant_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(hypsecant_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(tanh(z)) / sigma
        score(2) = sum(z * tanh(z) - 1.0_dp) / sigma
    end function hypsecant_score

    pure function halflogistic_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(halflogistic_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(tanh(0.5_dp * z)) / sigma
        if (any(data == mu)) score(1) = nan_value()
        score(2) = sum(z * tanh(0.5_dp * z) - 1.0_dp) / sigma
    end function halflogistic_score

    pure elemental function genpareto_shape_term(c, z) result(value)
        real(dp), intent(in) :: c !! finite generalized-Pareto shape parameter
        real(dp), intent(in) :: z !! positive standardized observation in the support interior
        real(dp) :: value

        if (c == 0.0_dp) then
            value = 0.5_dp * z * z - z
        else if (abs(c) <= 1.0e-5_dp .and. abs(c * z) <= 1.0e-3_dp) then
            value = 0.5_dp * z * z - z + &
                c * (z * z - (2.0_dp / 3.0_dp) * z**3) + &
                c * c * (0.75_dp * z**4 - z**3)
        else
            value = log1p_safe(c * z) / (c * c) - &
                (1.0_dp + 1.0_dp / c) * z / (1.0_dp + c * z)
        end if
    end function genpareto_shape_term

    pure function bernoulli_score_real(data, p, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(bernoulli_loglikelihood(data, p, loc)) .or. &
            p <= 0.0_dp .or. p >= 1.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = sum(count / p - (1.0_dp - count) / (1.0_dp - p))
    end function bernoulli_score_real

    pure function bernoulli_score_int(data, p, loc) result(score)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        score = bernoulli_score_real(real(data, dp), p, loc)
    end function bernoulli_score_int

    pure function poisson_score_real(data, mu, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: mu !! mean, strictly positive for a finite score
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(poisson_loglikelihood(data, mu, loc)) .or. mu <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = sum(count / mu - 1.0_dp)
    end function poisson_score_real

    pure function poisson_score_int(data, mu, loc) result(score)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: mu !! mean, strictly positive for a finite score
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        score = poisson_score_real(real(data, dp), mu, loc)
    end function poisson_score_int

    pure function geometric_score_real(data, p, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(geometric_loglikelihood(data, p, loc)) .or. &
            p <= 0.0_dp .or. p >= 1.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = sum(1.0_dp / p - (count - 1.0_dp) / (1.0_dp - p))
    end function geometric_score_real

    pure function geometric_score_int(data, p, loc) result(score)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        score = geometric_score_real(real(data, dp), p, loc)
    end function geometric_score_int

    pure function binomial_score_real(data, n, p, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! number of trials, a nonnegative integer value
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(binomial_loglikelihood(data, n, p, loc)) .or. &
            p <= 0.0_dp .or. p >= 1.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = sum(count / p - (n - count) / (1.0_dp - p))
    end function binomial_score_real

    pure function binomial_score_int(data, n, p, loc) result(score)
        integer, intent(in) :: data(:) !! independent count observations
        integer, intent(in) :: n !! number of trials, >= 0
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(1)

        score = binomial_score_real(real(data, dp), real(n, dp), p, loc)
    end function binomial_score_int

    pure function negative_binomial_score_real(data, n, p, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! required successes, > 0
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(2)

        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(negative_binomial_loglikelihood(data, n, p, loc)) .or. &
            p <= 0.0_dp .or. p >= 1.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = sum(digamma_positive(n + count) - digamma_positive(n) + log(p))
        score(2) = sum(n / p - count / (1.0_dp - p))
    end function negative_binomial_score_real

    pure function negative_binomial_score_int(data, n, p, loc) result(score)
        integer, intent(in) :: data(:) !! independent count observations
        real(dp), intent(in) :: n !! required successes, > 0
        real(dp), intent(in) :: p !! success probability, strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! integer support shift (default 0)
        real(dp) :: score(2)

        score = negative_binomial_score_real(real(data, dp), n, p, loc)
    end function negative_binomial_score_int

    pure function invgamma_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive inverse-gamma shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(invgamma_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(-digamma_positive(a) - log(z))
        score(2) = sum((a + 1.0_dp) / z - 1.0_dp / (z * z)) / sigma
        score(3) = sum(a - 1.0_dp / z) / sigma
    end function invgamma_score

    pure function invgauss_score(data, mu_shape, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: mu_shape !! positive inverse-Gaussian shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(invgauss_loglikelihood(data, mu_shape, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(z - mu_shape) / (mu_shape * mu_shape * mu_shape)
        score(2) = sum(1.5_dp / z + 0.5_dp / (mu_shape * mu_shape) - &
            0.5_dp / (z * z)) / sigma
        score(3) = sum(0.5_dp + 0.5_dp * z / (mu_shape * mu_shape) - &
            0.5_dp / z) / sigma
    end function invgauss_score

    pure function levy_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(2)

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(levy_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(1.5_dp / z - 0.5_dp / (z * z)) / sigma
        score(2) = sum(0.5_dp - 0.5_dp / z) / sigma
    end function levy_score

    pure function loglaplace_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive log-Laplace shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)

        integer :: i
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(loglaplace_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum(1.0_dp / c - abs(log(z)))
        score(2:3) = 0.0_dp
        do i = 1, size(z)
            if (z(i) < 1.0_dp) then
                score(2) = score(2) + (1.0_dp - c) / (z(i) * sigma)
                score(3) = score(3) - c / sigma
            else if (z(i) > 1.0_dp) then
                score(2) = score(2) + (1.0_dp + c) / (z(i) * sigma)
                score(3) = score(3) + c / sigma
            else
                score(2:3) = nan_value()
                return
            end if
        end do
    end function loglaplace_score


    pure function bradford_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive Bradford shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: score(3)

        real(dp) :: denom(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(bradford_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp) .or. any(z >= 1.0_dp)) then
            score = nan_value()
            return
        end if
        denom = 1.0_dp + c * z
        score(1) = sum(1.0_dp / c - 1.0_dp / ((1.0_dp + c) * log1p_safe(c)) - &
            z / denom)
        score(2) = sum(c / denom) / sigma
        score(3) = -sum(1.0_dp / denom) / sigma
    end function bradford_score

    pure function truncexpon_score(data, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! positive standardized upper endpoint
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)

        real(dp) :: denom
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(truncexpon_loglikelihood(data, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp) .or. any(z >= b)) then
            score = nan_value()
            return
        end if
        denom = -expm1_safe(-b)
        score(1) = -real(size(data), dp) * exp(-b) / denom
        score(2) = real(size(data), dp) / sigma
        score(3) = sum(z - 1.0_dp) / sigma
    end function truncexpon_score

    pure function fisk_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive Fisk shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)

        integer :: i
        real(dp) :: f
        real(dp) :: lz
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z

        if (.not. valid_loglike(fisk_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            lz = log(z)
            f = fisk_cdf(data(i), c, mu, sigma)
            score(1) = score(1) + 1.0_dp / c + lz * (1.0_dp - 2.0_dp * f)
            score(2) = score(2) + (2.0_dp * c * f - (c - 1.0_dp)) / (sigma * z)
            score(3) = score(3) + c * (2.0_dp * f - 1.0_dp) / sigma
        end do
    end function fisk_score

    pure function dweibull_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive double-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)

        integer :: i
        real(dp) :: logr
        real(dp) :: mu
        real(dp) :: power
        real(dp) :: r
        real(dp) :: sigma
        real(dp) :: z

        if (.not. valid_loglike(dweibull_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            r = abs(z)
            if (r <= 0.0_dp) then
                score = nan_value()
                return
            end if
            logr = log(r)
            power = exp(c * logr)
            score(1) = score(1) + 1.0_dp / c + logr * (1.0_dp - power)
            score(2) = score(2) + sign(1.0_dp, z) * &
                (c * power - (c - 1.0_dp)) / (sigma * r)
            score(3) = score(3) + c * (power - 1.0_dp) / sigma
        end do
    end function dweibull_score

    pure function alpha_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive alpha shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mills
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w(size(data))
        real(dp) :: z(size(data))
        if (.not. valid_loglike(alpha_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        w = a - 1.0_dp / z
        mills = exp(normal_logpdf(a) - normal_logcdf(a))
        score(1) = -sum(w) - real(size(data), dp) * mills
        score(2) = sum(2.0_dp / z + w / (z * z)) / sigma
        score(3) = sum(1.0_dp + w / z) / sigma
    end function alpha_score

    pure function fatiguelife_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive fatigue-life shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: g(size(data))
        real(dp) :: h(size(data))
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: z(size(data))
        if (.not. valid_loglike(fatiguelife_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        h = z - 2.0_dp + 1.0_dp / z
        g = 1.0_dp / (z + 1.0_dp) - 1.5_dp / z - &
            0.5_dp * (1.0_dp - 1.0_dp / (z * z)) / (c * c)
        score(1) = sum(-1.0_dp / c + h / (c * c * c))
        score(2) = -sum(g) / sigma
        score(3) = sum(-1.0_dp - z * g) / sigma
    end function fatiguelife_score

    pure function genlogistic_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive generalized-logistic shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: g(size(data))
        real(dp) :: l(size(data))
        real(dp) :: mu
        real(dp) :: r(size(data))
        real(dp) :: sigma
        real(dp) :: z(size(data))
        if (.not. valid_loglike(genlogistic_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        l = log1p_exp_neg(z)
        r = exp(-log1p_exp_pos(z))
        g = -1.0_dp + (c + 1.0_dp) * r
        score(1) = real(size(data), dp) / c - sum(l)
        score(2) = -sum(g) / sigma
        score(3) = sum(-1.0_dp - z * g) / sigma
    end function genlogistic_score

    pure function gennorm_score(data, beta, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: beta !! positive generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: az
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: z
        if (.not. valid_loglike(gennorm_loglikelihood(data, beta, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        score(1) = real(size(data), dp) * (1.0_dp / beta + &
            digamma_positive(1.0_dp / beta) / (beta * beta))
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            az = abs(z)
            if (az == 0.0_dp) then
                if (beta <= 1.0_dp) then
                    score = nan_value()
                    return
                end if
                score(3) = score(3) - 1.0_dp / sigma
            else
                t = exp(beta * log(az))
                score(1) = score(1) - t * log(az)
                score(2) = score(2) + beta * sign(1.0_dp, z) * &
                    exp((beta - 1.0_dp) * log(az)) / sigma
                score(3) = score(3) + (-1.0_dp + beta * t) / sigma
            end if
        end do
    end function gennorm_score

    pure function nakagami_score(data, nu, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: nu !! positive Nakagami shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: z(size(data))
        if (.not. valid_loglike(nakagami_loglikelihood(data, nu, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum(1.0_dp + log(nu) - digamma_positive(nu) + &
            2.0_dp * log(z) - z * z)
        score(2) = sum(2.0_dp * nu * z - (2.0_dp * nu - 1.0_dp) / z) / sigma
        score(3) = sum(2.0_dp * nu * (z * z - 1.0_dp)) / sigma
    end function nakagami_score

    pure function powernorm_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive power-normal shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: z(size(data)), mills(size(data))
        if (.not. valid_loglike(powernorm_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        mills = exp(normal_logpdf(z) - normal_logcdf(-z))
        score(1) = sum(1.0_dp / c + normal_logcdf(-z))
        score(2) = sum(z + (c - 1.0_dp) * mills) / sigma
        score(3) = sum(-1.0_dp + z * z + (c - 1.0_dp) * z * mills) / sigma
    end function powernorm_score

    pure function loggamma_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive log-gamma shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: z(size(data)), ez(size(data))
        if (.not. valid_loglike(loggamma_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        ez = exp(z)
        score(1) = sum(z - digamma_positive(c))
        score(2) = sum(ez - c) / sigma
        score(3) = sum(-1.0_dp - c * z + z * ez) / sigma
    end function loggamma_score

    pure function wald_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(2)
        real(dp) :: mu, sigma
        real(dp) :: z(size(data))
        if (.not. valid_loglike(wald_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        score(1) = sum(1.5_dp / z + 0.5_dp - 0.5_dp / (z * z)) / sigma
        score(2) = sum(0.5_dp + 0.5_dp * z - 0.5_dp / z) / sigma
    end function wald_score


    pure function gompertz_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive Gompertz shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: ez(size(data)), z(size(data))
        if (.not. valid_loglike(gompertz_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        ez = exp(z)
        score(1) = sum(1.0_dp / c - expm1_safe(z))
        score(2) = sum(c * ez - 1.0_dp) / sigma
        score(3) = sum(-1.0_dp - z + c * z * ez) / sigma
    end function gompertz_score

    pure function invweibull_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive inverse-Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: lz(size(data)), t(size(data)), z(size(data))
        if (.not. valid_loglike(invweibull_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        lz = log(z)
        t = exp(-c * lz)
        score(1) = sum(1.0_dp / c + (t - 1.0_dp) * lz)
        score(2) = sum((c + 1.0_dp - c * t) / z) / sigma
        score(3) = sum(c * (1.0_dp - t)) / sigma
    end function invweibull_score

    pure function betaprime_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive first beta-prime shape parameter
        real(dp), intent(in) :: b !! positive second beta-prime shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma
        real(dp) :: lz(size(data)), l1(size(data)), lprime(size(data)), z(size(data))
        if (.not. valid_loglike(betaprime_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        lz = log(z)
        l1 = log1p_safe(z)
        lprime = (a - 1.0_dp) / z - (a + b) / (1.0_dp + z)
        score(1) = sum(lz - l1 - digamma_positive(a) + digamma_positive(a + b))
        score(2) = sum(-l1 - digamma_positive(b) + digamma_positive(a + b))
        score(3) = -sum(lprime) / sigma
        score(4) = -sum(1.0_dp + z * lprime) / sigma
    end function betaprime_score

    pure function burr12_score(data, c, d, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive first Burr-XII shape parameter
        real(dp), intent(in) :: d !! positive second Burr-XII shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: mu, sigma
        real(dp) :: lz(size(data)), lp(size(data)), r(size(data)), s(size(data))
        real(dp) :: lprime(size(data)), z(size(data)), e
        if (.not. valid_loglike(burr12_loglikelihood(data, c, d, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        lz = log(z)
        s = c * lz
        do i = 1, size(data)
            lp(i) = log1pexp(s(i))
            if (s(i) >= 0.0_dp) then
                r(i) = 1.0_dp / (1.0_dp + exp(-s(i)))
            else
                e = exp(s(i))
                r(i) = e / (1.0_dp + e)
            end if
        end do
        lprime = (c - 1.0_dp - c * (d + 1.0_dp) * r) / z
        score(1) = sum(1.0_dp / c + lz * (1.0_dp - (d + 1.0_dp) * r))
        score(2) = sum(1.0_dp / d - lp)
        score(3) = -sum(lprime) / sigma
        score(4) = -sum(1.0_dp + z * lprime) / sigma
    end function burr12_score

    pure function genhalflogistic_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive generalized half-logistic shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: mu, sigma
        real(dp) :: gc(size(data)), gz(size(data)), lt, ratio, t, u, vc, z
        if (.not. valid_loglike(genhalflogistic_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp / c) then
                score = nan_value()
                return
            end if
            t = 1.0_dp - c * z
            lt = log(t)
            u = exp(lt / c)
            ratio = u / (1.0_dp + u)
            gz(i) = (c - 1.0_dp + 2.0_dp * ratio) / t
            vc = -z / (c * t) - lt / (c * c)
            gc(i) = -lt / (c * c) - (1.0_dp / c - 1.0_dp) * z / t - 2.0_dp * ratio * vc
        end do
        score(1) = sum(gc)
        score(2) = -sum(gz) / sigma
        score(3) = -sum(1.0_dp + ((data - mu) / sigma) * gz) / sigma
    end function genhalflogistic_score

    pure function exponpow_score(data, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! positive exponential-power shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: e(size(data)), lz(size(data)), xb(size(data)), z(size(data)), gz(size(data))
        if (.not. valid_loglike(exponpow_loglikelihood(data, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        lz = log(z)
        xb = exp(b * lz)
        e = exp(xb)
        gz = (b - 1.0_dp + b * xb * (1.0_dp - e)) / z
        score(1) = sum(1.0_dp / b + lz + xb * lz * (1.0_dp - e))
        score(2) = -sum(gz) / sigma
        score(3) = -sum(1.0_dp + z * gz) / sigma
    end function exponpow_score

    pure function exponweib_score(data, a, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive exponentiation shape parameter
        real(dp), intent(in) :: c !! positive Weibull shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma
        real(dp) :: lz(size(data)), q(size(data)), r(size(data)), t(size(data))
        real(dp) :: z(size(data)), gz(size(data)), ratio(size(data))
        if (.not. valid_loglike(exponweib_loglikelihood(data, a, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        lz = log(z)
        t = exp(c * lz)
        q = exp(-t)
        r = -expm1_safe(-t)
        ratio = q / r
        gz = ((a - 1.0_dp) * c * t * ratio - c * t + c - 1.0_dp) / z
        score(1) = sum(1.0_dp / a + log(r))
        score(2) = sum(1.0_dp / c + lz + t * lz * ((a - 1.0_dp) * ratio - 1.0_dp))
        score(3) = -sum(gz) / sigma
        score(4) = -sum(1.0_dp + z * gz) / sigma
    end function exponweib_score

    pure function powerlognorm_score(data, c, s, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive power shape parameter
        real(dp), intent(in) :: s !! positive lognormal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma
        real(dp) :: lz(size(data)), w(size(data)), z(size(data)), ratio(size(data)), gz(size(data))
        if (.not. valid_loglike(powerlognorm_loglikelihood(data, c, s, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        lz = log(z)
        w = lz / s
        ratio = exp(normal_logpdf(w) - normal_logcdf(-w))
        gz = (-1.0_dp + (-w - (c - 1.0_dp) * ratio) / s) / z
        score(1) = sum(1.0_dp / c + normal_logcdf(-w))
        score(2) = sum((-1.0_dp + w * w + (c - 1.0_dp) * w * ratio) / s)
        score(3) = -sum(gz) / sigma
        score(4) = -sum(1.0_dp + z * gz) / sigma
    end function powerlognorm_score

    pure function levy_l_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(2)
        real(dp) :: mu, sigma
        real(dp) :: y(size(data)), z(size(data))

        if (.not. valid_loglike(levy_l_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        y = -z
        if (any(y <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        score(1) = sum(-1.5_dp / y + 0.5_dp / (y * y)) / sigma
        score(2) = sum(0.5_dp - 0.5_dp / y) / sigma
    end function levy_l_score

    pure function weibull_max_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive Weibull-maximum shape parameter
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: logy(size(data)), mu, sigma, y(size(data)), yc(size(data)), z(size(data))

        if (.not. valid_loglike(weibull_max_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        y = -z
        if (any(y <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        logy = log(y)
        yc = exp(c * logy)
        score(1) = sum(1.0_dp / c + logy * (1.0_dp - yc))
        score(2) = -sum(c * yc / y - (c - 1.0_dp) / y) / sigma
        score(3) = c * sum(yc - 1.0_dp) / sigma
    end function weibull_max_score

    pure function rdist_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive R-distribution shape parameter
        real(dp), intent(in), optional :: loc !! center of support (default 0)
        real(dp), intent(in), optional :: scale !! positive half-width (default 1)
        real(dp) :: score(3)
        real(dp) :: denom(size(data)), mu, sigma, z(size(data))

        if (.not. valid_loglike(rdist_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(abs(z) >= 1.0_dp)) then
            score = nan_value()
            return
        end if
        denom = 1.0_dp - z * z
        score(1) = 0.5_dp * sum(log(denom) - digamma_positive(0.5_dp * c) + &
            digamma_positive(0.5_dp * (c + 1.0_dp)))
        score(2) = sum((c - 2.0_dp) * z / denom) / sigma
        score(3) = sum(((c - 1.0_dp) * z * z - 1.0_dp) / denom) / sigma
    end function rdist_score

    pure function skewcauchy_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! skewness parameter strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: h, invq, mu, q, r, sigma, sgn, z

        if (.not. valid_loglike(skewcauchy_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z < 0.0_dp) then
                sgn = -1.0_dp
                h = 1.0_dp - a
            else if (z > 0.0_dp) then
                sgn = 1.0_dp
                h = 1.0_dp + a
            else
                sgn = 0.0_dp
                h = 1.0_dp
            end if
            q = z / h
            if (abs(q) <= 1.0_dp) then
                r = q * q / (1.0_dp + q * q)
                score(2) = score(2) + 2.0_dp * q / (h * (1.0_dp + q * q) * sigma)
            else
                invq = 1.0_dp / q
                r = 1.0_dp / (1.0_dp + invq * invq)
                score(2) = score(2) + 2.0_dp * invq / (h * (1.0_dp + invq * invq) * sigma)
            end if
            score(1) = score(1) + 2.0_dp * sgn * r / h
            score(3) = score(3) + (2.0_dp * r - 1.0_dp) / sigma
        end do
    end function skewcauchy_score


    pure function dgamma_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive finite shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        real(dp) :: mu, sigma
        real(dp) :: r(size(data)), sgn(size(data)), z(size(data))

        if (.not. valid_loglike(dgamma_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        r = abs(z)
        if (any(r <= 0.0_dp)) then
            score = nan_value()
            return
        end if
        sgn = merge(1.0_dp, -1.0_dp, z > 0.0_dp)
        score(1) = sum(log(r) - digamma_positive(a))
        score(2) = sum(sgn * (1.0_dp - (a - 1.0_dp) / r)) / sigma
        score(3) = sum(r - a) / sigma
    end function dgamma_score

    pure function laplace_asymmetric_score(data, kappa, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: kappa !! positive finite asymmetry parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: base, mu, sigma, z

        if (.not. valid_loglike(laplace_asymmetric_loglikelihood(data, kappa, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        base = (1.0_dp - kappa * kappa) / (kappa * (1.0_dp + kappa * kappa))
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z >= 0.0_dp) then
                score(1) = score(1) + base - z
                score(2) = score(2) + kappa / sigma
                score(3) = score(3) + (kappa * z - 1.0_dp) / sigma
            else
                score(1) = score(1) + base - z / (kappa * kappa)
                score(2) = score(2) - 1.0_dp / (kappa * sigma)
                score(3) = score(3) + (-z / kappa - 1.0_dp) / sigma
            end if
        end do
        if (any(data == mu)) score(2) = nan_value()
    end function laplace_asymmetric_score

    pure function truncnorm_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! finite standardized lower truncation point
        real(dp), intent(in) :: b !! finite standardized upper truncation point greater than a
        real(dp), intent(in), optional :: loc !! parent-normal location (default 0)
        real(dp), intent(in), optional :: scale !! parent-normal standard deviation (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, pa, pb, sigma
        real(dp) :: z(size(data))

        if (.not. valid_loglike(truncnorm_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= a) .or. any(z >= b)) then
            score = nan_value()
            return
        end if
        pa = truncnorm_pdf(a, a, b)
        pb = truncnorm_pdf(b, a, b)
        score(1) = real(size(data), dp) * pa
        score(2) = -real(size(data), dp) * pb
        score(3) = sum(z) / sigma
        score(4) = sum(z * z - 1.0_dp) / sigma
    end function truncnorm_score

    pure function loguniform_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive finite standardized lower endpoint
        real(dp), intent(in) :: b !! finite standardized upper endpoint greater than a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma, width
        real(dp) :: z(size(data))

        if (.not. valid_loglike(loguniform_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        z = (data - mu) / sigma
        if (any(z <= a) .or. any(z >= b)) then
            score = nan_value()
            return
        end if
        width = log(b) - log(a)
        score(1) = real(size(data), dp) / (a * width)
        score(2) = -real(size(data), dp) / (b * width)
        score(3) = sum(1.0_dp / z) / sigma
        score(4) = 0.0_dp
    end function loguniform_score

    pure elemental function log1p_exp_neg(z) result(value)
        real(dp), intent(in) :: z !! real argument
        real(dp) :: value
        if (z >= 0.0_dp) then
            value = log1p_safe(exp(-z))
        else
            value = -z + log1p_safe(exp(z))
        end if
    end function log1p_exp_neg

    pure elemental function log1p_exp_pos(z) result(value)
        real(dp), intent(in) :: z !! real argument
        real(dp) :: value
        if (z >= 0.0_dp) then
            value = z + log1p_safe(exp(-z))
        else
            value = log1p_safe(exp(z))
        end if
    end function log1p_exp_pos

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = optional_loc(loc)
        sigma = 1.0_dp
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

    pure function optional_loc(loc) result(value)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp) :: value

        value = 0.0_dp
        if (present(loc)) value = loc
    end function optional_loc

    pure function valid_loglike(value) result(valid)
        real(dp), intent(in) :: value !! summed log-likelihood to validate
        logical :: valid

        valid = ieee_is_finite(value)
    end function valid_loglike

    pure function nan_value() result(value)
        real(dp) :: value

        value = quiet_nan(0.0_dp)
    end function nan_value

    pure function foldnorm_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: fc, fz, mu, r, sigma, wa, wb, z
        if (.not. valid_loglike(foldnorm_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            r = exp(-2.0_dp * z * c)
            wa = 1.0_dp / (1.0_dp + r)
            wb = r * wa
            fc = wa * (z - c) - wb * (z + c)
            fz = -wa * (z - c) - wb * (z + c)
            score(1) = score(1) + fc
            score(2) = score(2) - fz / sigma
            score(3) = score(3) - (1.0_dp + z * fz) / sigma
        end do
    end function foldnorm_score

    pure function foldcauchy_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: a, b, fc, fz, mu, ru, rv, sigma, u, v, wa, wb, z
        if (.not. valid_loglike(foldcauchy_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            u = z - c
            v = z + c
            a = inv_one_plus_square_score(u)
            b = inv_one_plus_square_score(v)
            wa = a / (a + b)
            wb = b / (a + b)
            ru = x_over_one_plus_square(u)
            rv = x_over_one_plus_square(v)
            fc = 2.0_dp * (wa * ru - wb * rv)
            fz = -2.0_dp * (wa * ru + wb * rv)
            score(1) = score(1) + fc
            score(2) = score(2) - fz / sigma
            score(3) = score(3) - (1.0_dp + z * fz) / sigma
        end do
    end function foldcauchy_score

    pure function recipinvgauss_score(data, mu_shape, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: mu_shape !! positive reciprocal inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: fm, fz, mu, sigma, z
        if (.not. valid_loglike(recipinvgauss_loglikelihood(data, mu_shape, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            fm = 1.0_dp / (mu_shape**3 * z) - 1.0_dp / (mu_shape * mu_shape)
            fz = 0.5_dp / (mu_shape * mu_shape * z * z) - 0.5_dp / z - 0.5_dp
            score(1) = score(1) + fm
            score(2) = score(2) - fz / sigma
            score(3) = score(3) - (1.0_dp + z * fz) / sigma
        end do
    end function recipinvgauss_score

    pure function truncpareto_score(data, b, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: fz, invden, lnc, mu, sigma, t, z
        if (.not. valid_loglike(truncpareto_loglikelihood(data, b, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        if (b == 0.0_dp .or. c <= 1.0_dp) then
            score = nan_value()
            return
        end if
        lnc = log(c)
        t = b * lnc
        invden = inverse_expm1_score(t)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 1.0_dp .or. z >= c) then
                score = nan_value()
                return
            end if
            fz = -(b + 1.0_dp) / z
            score(1) = score(1) + 1.0_dp / b - lnc * invden - log(z)
            score(2) = score(2) - b * invden / c
            score(3) = score(3) - fz / sigma
            score(4) = score(4) - (1.0_dp + z * fz) / sigma
        end do
    end function truncpareto_score


    pure function exponnorm_score(data, k, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: k !! positive exponential-to-normal scale ratio
        real(dp), intent(in), optional :: loc !! normal-component location (default 0)
        real(dp), intent(in), optional :: scale !! normal-component scale (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: fz, invk, mu, ratio, sigma, w, z
        if (.not. valid_loglike(exponnorm_loglikelihood(data, k, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        invk = 1.0_dp / k
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            w = z - invk
            ratio = exp(normal_logpdf(w) - normal_logcdf(w))
            score(1) = score(1) + invk * invk * (z + ratio - invk) - invk
            fz = -invk + ratio
            score(2) = score(2) - fz / sigma
            score(3) = score(3) - (1.0_dp + z * fz) / sigma
        end do
    end function exponnorm_score

    pure function johnsonsb_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! finite first Johnson SB shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SB shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: fz, mu, q, sigma, t, z
        if (.not. valid_loglike(johnsonsb_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp) then
                score = nan_value()
                return
            end if
            t = log(z) - log1p_safe(-z)
            q = a + b * t
            score(1) = score(1) - q
            score(2) = score(2) + 1.0_dp / b - q * t
            fz = (2.0_dp * z - 1.0_dp - q * b) / (z * (1.0_dp - z))
            score(3) = score(3) - fz / sigma
            score(4) = score(4) - (1.0_dp + z * fz) / sigma
        end do
    end function johnsonsb_score

    pure function johnsonsu_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! finite first Johnson SU shape parameter
        real(dp), intent(in) :: b !! positive second Johnson SU shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: fz, h, mu, q, sigma, t, z
        if (.not. valid_loglike(johnsonsu_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            h = hypot(z, 1.0_dp)
            t = asinh(z)
            q = a + b * t
            score(1) = score(1) - q
            score(2) = score(2) + 1.0_dp / b - q * t
            fz = -x_over_one_plus_square(z) - q * b / h
            score(3) = score(3) - fz / sigma
            score(4) = score(4) - (1.0_dp + z * fz) / sigma
        end do
    end function johnsonsu_score

    pure function trapezoid_score(data, c, d, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! standardized left plateau edge in [0,1]
        real(dp), intent(in) :: d !! standardized right plateau edge in [c,1]
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive support width (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: dc, dd, den, fz, mu, sigma, z
        if (.not. valid_loglike(trapezoid_loglikelihood(data, c, d, loc, scale))) then
            score = nan_value()
            return
        end if
        if (c <= 0.0_dp .or. d >= 1.0_dp .or. c >= d) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        den = 1.0_dp + d - c
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp .or. z == c .or. z == d) then
                score = nan_value()
                return
            end if
            dc = 1.0_dp / den
            dd = -1.0_dp / den
            if (z < c) then
                dc = dc - 1.0_dp / c
                fz = 1.0_dp / z
            else if (z < d) then
                fz = 0.0_dp
            else
                dd = dd + 1.0_dp / (1.0_dp - d)
                fz = -1.0_dp / (1.0_dp - z)
            end if
            score(1) = score(1) + dc
            score(2) = score(2) + dd
            score(3) = score(3) - fz / sigma
            score(4) = score(4) - (1.0_dp + z * fz) / sigma
        end do
    end function trapezoid_score

    pure elemental function inv_one_plus_square_score(x) result(y)
        real(dp), intent(in) :: x !! real argument
        real(dp) :: y, a, r
        a = abs(x)
        if (a <= 1.0_dp) then
            y = 1.0_dp / (1.0_dp + x * x)
        else
            r = 1.0_dp / a
            y = r * r / (1.0_dp + r * r)
        end if
    end function inv_one_plus_square_score

    pure elemental function x_over_one_plus_square(x) result(y)
        real(dp), intent(in) :: x !! real argument
        real(dp) :: y, a, r
        a = abs(x)
        if (a <= 1.0_dp) then
            y = x / (1.0_dp + x * x)
        else
            r = 1.0_dp / x
            y = r / (1.0_dp + r * r)
        end if
    end function x_over_one_plus_square

    pure elemental function inverse_expm1_score(x) result(y)
        real(dp), intent(in) :: x !! argument of 1/expm1(x)
        real(dp) :: y
        if (x > log(huge(1.0_dp))) then
            y = 0.0_dp
        else
            y = 1.0_dp / expm1_safe(x)
        end if
    end function inverse_expm1_score

    pure function burr_score(data, c, d, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive finite first Burr shape parameter
        real(dp), intent(in) :: d !! positive finite second Burr shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: fz, l, lz, mu, r, sigma, t, z, zfz

        if (.not. valid_loglike(burr_loglikelihood(data, c, d, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            lz = log(z)
            t = c * lz
            r = logistic01(-t)
            l = log1pexp(-t)
            score(1) = score(1) + 1.0_dp / c - lz + (d + 1.0_dp) * lz * r
            score(2) = score(2) + 1.0_dp / d - l
            zfz = -(c + 1.0_dp) + (d + 1.0_dp) * c * r
            fz = zfz / z
            score(3) = score(3) - fz / sigma
            score(4) = score(4) - (1.0_dp + zfz) / sigma
        end do
    end function burr_score

    pure function mielke_score(data, k, s, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: k !! positive finite beta-kappa first shape parameter
        real(dp), intent(in) :: s !! positive finite beta-kappa second shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: fz, l, lz, mu, sigma, t, u, z, zfz

        if (.not. valid_loglike(mielke_loglikelihood(data, k, s, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            lz = log(z)
            t = s * lz
            u = logistic01(t)
            l = log1pexp(t)
            score(1) = score(1) + 1.0_dp / k + lz - l / s
            score(2) = score(2) + k * l / (s * s) - &
                (1.0_dp + k / s) * lz * u
            zfz = (k - 1.0_dp) - (s + k) * u
            fz = zfz / z
            score(3) = score(3) - fz / sigma
            score(4) = score(4) - (1.0_dp + zfz) / sigma
        end do
    end function mielke_score

    pure function gibrat_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(2)
        real(dp) :: base(3)

        base = lognormal_score(data, 1.0_dp, loc, scale)
        score = base(2:3)
    end function gibrat_score

    pure function wrapcauchy_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! wrapped-Cauchy concentration parameter in (0,1)
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: den, fz, mu, sigma, sh, z, zfz

        if (.not. valid_loglike(wrapcauchy_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= 2.0_dp * scifort_pi) then
                score = nan_value()
                return
            end if
            sh = sin(0.5_dp * z)
            den = (1.0_dp - c) * (1.0_dp - c) + 4.0_dp * c * sh * sh
            score(1) = score(1) - 2.0_dp * c / ((1.0_dp - c) * (1.0_dp + c)) - &
                2.0_dp * (c - cos(z)) / den
            fz = -2.0_dp * c * sin(z) / den
            zfz = z * fz
            score(2) = score(2) - fz / sigma
            score(3) = score(3) - (1.0_dp + zfz) / sigma
        end do
    end function wrapcauchy_score


    pure function genextreme_score(data, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! finite generalized-extreme-value shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: ac, dz, logt, logu, mu, sigma, t, u, z

        if (.not. valid_loglike(genextreme_loglikelihood(data, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (c == 0.0_dp) then
                t = exp(-z)
                score(1) = score(1) + z - 0.5_dp * z * z * (1.0_dp - t)
                dz = t - 1.0_dp
            else
                u = 1.0_dp - c * z
                if (u <= 0.0_dp) then
                    score = nan_value()
                    return
                end if
                logu = log1p_safe(-c * z)
                logt = logu / c
                t = exp(logt)
                ac = (-c * z / u - logu) / (c * c)
                score(1) = score(1) + (1.0_dp - t) * ac + z / u
                dz = (t - 1.0_dp + c) / u
            end if
            score(2) = score(2) - dz / sigma
            score(3) = score(3) - (1.0_dp + z * dz) / sigma
        end do
    end function genextreme_score

    pure function kappa3_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive finite kappa shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: den, dz, logden, logz, mu, sigma, z, za

        if (.not. valid_loglike(kappa3_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            logz = log(z)
            za = exp(a * logz)
            den = a + za
            logden = log(den)
            score(1) = score(1) + 1.0_dp / a + logden / (a * a) - &
                (1.0_dp + 1.0_dp / a) * (1.0_dp + za * logz) / den
            dz = -(a + 1.0_dp) * za / (z * den)
            score(2) = score(2) - dz / sigma
            score(3) = score(3) - (1.0_dp + z * dz) / sigma
        end do
    end function kappa3_score

    pure function kappa4_score(data, h, k, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: h !! finite first kappa shape parameter
        real(dp), intent(in) :: k !! finite second kappa shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: ak, base, dz, logbase, logt, logv, mu, sigma, t, v, z

        if (.not. valid_loglike(kappa4_loglikelihood(data, h, k, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (k == 0.0_dp) then
                base = 1.0_dp
                logbase = 0.0_dp
                logt = -z
                ak = -0.5_dp * z * z
            else
                base = 1.0_dp - k * z
                if (base <= 0.0_dp) then
                    score = nan_value()
                    return
                end if
                logbase = log1p_safe(-k * z)
                logt = logbase / k
                ak = (-k * z / base - logbase) / (k * k)
            end if
            t = exp(logt)
            v = 1.0_dp - h * t
            if (v <= 0.0_dp) then
                score = nan_value()
                return
            end if
            if (h == 0.0_dp) then
                score(1) = score(1) + t - 0.5_dp * t * t
                logv = 0.0_dp
            else
                logv = log1p_safe(-h * t)
                score(1) = score(1) - logv / (h * h) - &
                    (1.0_dp / h - 1.0_dp) * t / v
            end if
            score(2) = score(2) + z / base + ak * (1.0_dp - t) / v
            dz = (k - 1.0_dp) / base + (1.0_dp - h) * t / (base * v)
            score(3) = score(3) - dz / sigma
            score(4) = score(4) - (1.0_dp + z * dz) / sigma
        end do
    end function kappa4_score

    pure function truncweibull_min_score(data, c, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: c !! positive finite Weibull shape parameter
        real(dp), intent(in) :: a !! standardized lower truncation point, >= 0
        real(dp), intent(in) :: b !! standardized upper truncation point, > a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(5)
        integer :: i
        real(dp) :: ac, alog, apowm1, bc, blog, bpowm1, dcden, den, dz
        real(dp) :: ea, eb, mu, sigma, z, zc, logz

        if (.not. valid_loglike(truncweibull_min_loglikelihood(data, c, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        ac = a ** c
        bc = b ** c
        ea = exp(-ac)
        eb = exp(-bc)
        den = ea - eb
        if (a == 0.0_dp) then
            alog = 0.0_dp
            if (c > 1.0_dp) then
                apowm1 = 0.0_dp
            else if (c == 1.0_dp) then
                apowm1 = 1.0_dp
            else
                score = nan_value()
                return
            end if
        else
            alog = log(a)
            apowm1 = a ** (c - 1.0_dp)
        end if
        blog = log(b)
        bpowm1 = b ** (c - 1.0_dp)
        dcden = -ac * alog * ea + bc * blog * eb
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= a .or. z >= b .or. z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            logz = log(z)
            zc = z ** c
            score(1) = score(1) + 1.0_dp / c + logz - zc * logz - dcden / den
            score(2) = score(2) + c * apowm1 * ea / den
            score(3) = score(3) - c * bpowm1 * eb / den
            dz = ((c - 1.0_dp) / z - c * zc / z)
            score(4) = score(4) - dz / sigma
            score(5) = score(5) - (1.0_dp + z * dz) / sigma
        end do
    end function truncweibull_min_score

    pure elemental function logistic01(x) result(y)
        real(dp), intent(in) :: x !! real logistic transform argument
        real(dp) :: y, e
        if (x >= 0.0_dp) then
            e = exp(-x)
            y = 1.0_dp / (1.0_dp + e)
        else
            e = exp(x)
            y = e / (1.0_dp + e)
        end if
    end function logistic01

    pure function gengamma_score(data, a, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive finite gamma-shape parameter
        real(dp), intent(in) :: c !! finite nonzero power parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: dlogdz, logz, mu, psi_a, sigma, t, z

        if (.not. valid_loglike(gengamma_loglikelihood(data, a, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        psi_a = digamma_positive(a)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            logz = log(z)
            t = exp(c * logz)
            score(1) = score(1) + c * logz - psi_a
            score(2) = score(2) + 1.0_dp / c + a * logz - t * logz
            dlogdz = (c * a - 1.0_dp - c * t) / z
            score(3) = score(3) - dlogdz / sigma
            score(4) = score(4) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function gengamma_score

    pure function halfgennorm_score(data, beta, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: mu, psi_b, sigma, t, z, logz, dlogdz

        if (.not. valid_loglike(halfgennorm_loglikelihood(data, beta, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        psi_b = digamma_positive(1.0_dp / beta)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            logz = log(z)
            t = exp(beta * logz)
            score(1) = score(1) + 1.0_dp / beta + psi_b / (beta * beta) - t * logz
            dlogdz = -beta * t / z
            score(2) = score(2) - dlogdz / sigma
            score(3) = score(3) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function halfgennorm_score

    pure function argus_score(data, chi, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: chi !! positive finite ARGUS shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! support width, > 0 (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: dlogdz, dlogphi, mu, phi, sigma, z, z2

        if (.not. valid_loglike(argus_loglikelihood(data, chi, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        phi = 0.5_dp * gammainc(1.5_dp, 0.5_dp * chi * chi)
        if (.not. (phi > 0.0_dp)) then
            score = nan_value()
            return
        end if
        dlogphi = chi * chi * scifort_inv_sqrt_two_pi * exp(-0.5_dp * chi * chi) / phi
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp) then
                score = nan_value()
                return
            end if
            z2 = z * z
            score(1) = score(1) + 3.0_dp / chi - dlogphi - chi * (1.0_dp - z2)
            dlogdz = 1.0_dp / z - z / (1.0_dp - z2) + chi * chi * z
            score(2) = score(2) - dlogdz / sigma
            score(3) = score(3) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function argus_score

    pure function erlang_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive shape; conventionally integer-valued
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        if (.not. valid_loglike(erlang_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        score = gamma_score(data, a, loc, scale)
    end function erlang_score


    pure function crystalball_score(data, beta, m, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: beta !! positive tail-transition magnitude
        real(dp), intent(in) :: m !! power-law exponent, > 1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: d, dlogdz, dlogn_beta, dlogn_m, dmd_beta, dmd_m
        real(dp) :: e, gauss, mu, mden, sigma, tail, z

        if (.not. valid_loglike(crystalball_loglikelihood(data, beta, m, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        e = exp(-0.5_dp * beta * beta)
        tail = e * m / (beta * (m - 1.0_dp))
        gauss = exp(scifort_log_sqrt_two_pi) * normal_cdf(beta)
        mden = tail + gauss
        dmd_beta = tail * (-beta - 1.0_dp / beta) + e
        dmd_m = -e / (beta * (m - 1.0_dp) ** 2)
        dlogn_beta = -dmd_beta / mden
        dlogn_m = -dmd_m / mden
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            score(1) = score(1) + dlogn_beta
            score(2) = score(2) + dlogn_m
            if (z > -beta) then
                dlogdz = -z
            else
                d = m / beta - beta - z
                score(1) = score(1) - m / beta - beta + &
                    m * (m / (beta * beta) + 1.0_dp) / d
                score(2) = score(2) + log(m / beta) + 1.0_dp - log(d) - &
                    m / (beta * d)
                dlogdz = m / d
            end if
            score(3) = score(3) - dlogdz / sigma
            score(4) = score(4) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function crystalball_score

    pure function jf_skew_t_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive left-shape parameter
        real(dp), intent(in) :: b !! positive right-shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: dlogdz, drds, drdz, h, mu, psi_a, psi_b, psi_s
        real(dp) :: r, rterm, s, sigma, z

        if (.not. valid_loglike(jf_skew_t_loglikelihood(data, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        s = a + b
        psi_a = digamma_positive(a)
        psi_b = digamma_positive(b)
        psi_s = digamma_positive(s)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            h = s + z * z
            r = z / sqrt(h)
            drds = -0.5_dp * r / h
            drdz = s / (h * sqrt(h))
            rterm = drds * ((a + 0.5_dp) / (1.0_dp + r) - &
                (b + 0.5_dp) / (1.0_dp - r))
            score(1) = score(1) + log1p_safe(r) + rterm - scifort_log_two - &
                psi_a + psi_s - 0.5_dp / s
            score(2) = score(2) + log1p_safe(-r) + rterm - scifort_log_two - &
                psi_b + psi_s - 0.5_dp / s
            dlogdz = drdz * ((a + 0.5_dp) / (1.0_dp + r) - &
                (b + 0.5_dp) / (1.0_dp - r))
            score(3) = score(3) - dlogdz / sigma
            score(4) = score(4) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function jf_skew_t_score

    pure function pearson3_score(data, skew, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: skew !! finite Pearson III skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: alpha, beta, dalpha, dlogdz, dt, mu, psi_a, sigma, t, z

        if (.not. valid_loglike(pearson3_loglikelihood(data, skew, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        if (abs(skew) < 1.6e-5_dp) then
            do i = 1, size(data)
                z = (data(i) - mu) / sigma
                dlogdz = -z
                score(2) = score(2) - dlogdz / sigma
                score(3) = score(3) - (1.0_dp + z * dlogdz) / sigma
            end do
            return
        end if
        beta = 2.0_dp / skew
        alpha = beta * beta
        psi_a = digamma_positive(alpha)
        dalpha = -2.0_dp * alpha / skew
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            t = alpha + beta * z
            if (t <= 0.0_dp) then
                score = nan_value()
                return
            end if
            dt = -(alpha + t) / skew
            score(1) = score(1) - 1.0_dp / skew + &
                dalpha * (log(t) - psi_a) + dt * ((alpha - 1.0_dp) / t - 1.0_dp)
            dlogdz = beta * ((alpha - 1.0_dp) / t - 1.0_dp)
            score(2) = score(2) - dlogdz / sigma
            score(3) = score(3) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function pearson3_score

    pure function rel_breitwigner_score(data, rho, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: rho !! positive resonance-to-width ratio
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: dlogc, dlogdz, dq, mu, q, sigma, u, z

        if (.not. valid_loglike(rel_breitwigner_loglikelihood(data, rho, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        u = sqrt(1.0_dp + 1.0_dp / (rho * rho))
        dlogc = -(u + 2.0_dp) / (2.0_dp * u * u * (1.0_dp + u) * rho ** 3)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            q = z * z / rho - rho
            dq = -z * z / (rho * rho) - 1.0_dp
            score(1) = score(1) + dlogc - 2.0_dp * q * dq / (1.0_dp + q * q)
            dlogdz = -4.0_dp * z * q / (rho * (1.0_dp + q * q))
            score(2) = score(2) - dlogdz / sigma
            score(3) = score(3) - (1.0_dp + z * dlogdz) / sigma
        end do
    end function rel_breitwigner_score

    pure function genexpon_score(data, a, b, c, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive first rate parameter
        real(dp), intent(in) :: b !! positive second rate parameter
        real(dp), intent(in) :: c !! positive decay parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(5)
        integer :: i
        real(dp) :: e, gz, h, mu, sigma, u, z

        if (.not. valid_loglike(genexpon_loglikelihood(data, a, b, c, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            e = exp(-c * z)
            u = 1.0_dp - e
            h = a + b * u
            gz = b * c * e / h - (a + b) + b * e
            score(1) = score(1) + 1.0_dp / h - z
            score(2) = score(2) + u / h - z + u / c
            score(3) = score(3) + b * z * e / h + b * (c * z * e - u) / (c * c)
            score(4) = score(4) - gz / sigma
            score(5) = score(5) - (1.0_dp + z * gz) / sigma
        end do
    end function genexpon_score

    pure function skewnorm_score(data, a, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! finite skewness parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: gz, mu, r, sigma, u, z

        if (.not. valid_loglike(skewnorm_loglikelihood(data, a, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            u = a * z
            r = exp(normal_logpdf(u) - normal_logcdf(u))
            gz = -z + a * r
            score(1) = score(1) + z * r
            score(2) = score(2) - gz / sigma
            score(3) = score(3) - (1.0_dp + z * gz) / sigma
        end do
    end function skewnorm_score

    pure function tukeylambda_score(data, lam, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: lam !! finite Tukey lambda shape parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: d, dl, dpv, lp, lq, mu, p, q, qlam, sigma, u, v, z, gz

        if (.not. valid_loglike(tukeylambda_loglikelihood(data, lam, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            p = tukeylambda_cdf(data(i), lam, mu, sigma)
            if (p <= 0.0_dp .or. p >= 1.0_dp) then
                score = nan_value()
                return
            end if
            q = 1.0_dp - p
            lp = log(p)
            lq = log(q)
            u = exp((lam - 1.0_dp) * lp)
            v = exp((lam - 1.0_dp) * lq)
            d = u + v
            dpv = (lam - 1.0_dp) * (u / p - v / q)
            dl = u * lp + v * lq
            if (abs(lam) < 1.0e-6_dp) then
                qlam = 0.5_dp * (lp * lp - lq * lq) + &
                    lam * (lp**3 - lq**3) / 3.0_dp
            else
                qlam = (lam * (exp(lam * lp) * lp - exp(lam * lq) * lq) - &
                    (exp(lam * lp) - exp(lam * lq))) / (lam * lam)
            end if
            gz = -dpv / (d * d)
            score(1) = score(1) - dl / d + dpv * qlam / (d * d)
            score(2) = score(2) - gz / sigma
            score(3) = score(3) - (1.0_dp + z * gz) / sigma
        end do
    end function tukeylambda_score

    pure function rice_score(data, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: b !! nonnegative finite Rice shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: gz, mu, r, sigma, z

        if (.not. valid_loglike(rice_loglikelihood(data, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            r = rice_i1_i0_ratio(z * b)
            gz = 1.0_dp / z - z + b * r
            score(1) = score(1) + z * r - b
            score(2) = score(2) - gz / sigma
            score(3) = score(3) - (1.0_dp + z * gz) / sigma
        end do
    end function rice_score


    pure function dpareto_lognorm_score(data, u, s, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent positive observations
        real(dp), intent(in) :: u !! finite lognormal location shape parameter
        real(dp), intent(in) :: s !! positive lognormal scale shape parameter
        real(dp), intent(in) :: a !! positive upper-Pareto shape parameter
        real(dp), intent(in) :: b !! positive lower-Pareto shape parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive outer scale (default 1)
        real(dp) :: score(6)
        integer :: i
        real(dp) :: gt, lr1, lr2, lrs, mu, q1, q2, sigma, w1, w2, x1, x2, y, z

        if (.not. valid_loglike(dpareto_lognorm_loglikelihood(data, u, s, a, b, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            y = (data(i) - mu) / sigma
            if (y <= 0.0_dp) then
                score = nan_value()
                return
            end if
            z = (log(y) - u) / s
            x1 = a * s - z
            x2 = b * s + z
            lr1 = normal_logsf(x1) - normal_logpdf(x1)
            lr2 = normal_logsf(x2) - normal_logpdf(x2)
            lrs = max(lr1, lr2) + log(exp(lr1 - max(lr1, lr2)) + exp(lr2 - max(lr1, lr2)))
            w1 = exp(lr1 - lrs)
            w2 = exp(lr2 - lrs)
            q1 = w1 * (x1 - exp(-lr1))
            q2 = w2 * (x2 - exp(-lr2))
            score(1) = score(1) + z / s + (q1 - q2) / s
            score(2) = score(2) + z * z / s + q1 * (a + z / s) + q2 * (b - z / s)
            score(3) = score(3) + 1.0_dp / a - 1.0_dp / (a + b) + s * q1
            score(4) = score(4) + 1.0_dp / b - 1.0_dp / (a + b) + s * q2
            gt = -1.0_dp - z / s + (q2 - q1) / s
            score(5) = score(5) - gt / (sigma * y)
            score(6) = score(6) - (1.0_dp + gt) / sigma
        end do
    end function dpareto_lognorm_score

    pure function vonmises_score(data, kappa, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent angular observations
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! circular location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: gz, mu, ratio, sigma, z

        if (.not. valid_loglike(vonmises_loglikelihood(data, kappa, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        ratio = vonmises_i1_i0_ratio(kappa)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            gz = -kappa * sin(z)
            score(1) = score(1) + cos(z) - ratio
            score(2) = score(2) - gz / sigma
            score(3) = score(3) - (1.0_dp + z * gz) / sigma
        end do
    end function vonmises_score

    pure function vonmises_line_score(data, kappa, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent bounded angular observations
        real(dp), intent(in) :: kappa !! nonnegative concentration parameter
        real(dp), intent(in), optional :: loc !! support center (default 0)
        real(dp), intent(in), optional :: scale !! positive angular scale (default 1)
        real(dp) :: score(3)
        integer :: i
        real(dp) :: gz, mu, ratio, sigma, z

        if (.not. valid_loglike(vonmises_line_loglikelihood(data, kappa, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        ratio = vonmises_i1_i0_ratio(kappa)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (abs(z) >= scifort_pi) then
                score = nan_value()
                return
            end if
            gz = -kappa * sin(z)
            score(1) = score(1) + cos(z) - ratio
            score(2) = score(2) - gz / sigma
            score(3) = score(3) - (1.0_dp + z * gz) / sigma
        end do
    end function vonmises_line_score

    pure function kstwobign_score(data, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent positive observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(2)
        integer :: i
        real(dp) :: gz, mu, sigma, z

        if (.not. valid_loglike(kstwobign_loglikelihood(data, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            gz = kstwobign_logpdf_derivative(z)
            score(1) = score(1) - gz / sigma
            score(2) = score(2) - (1.0_dp + z * gz) / sigma
        end do
    end function kstwobign_score


    pure function irwinhall_score(data, n, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: n !! positive integer number of summed uniforms
        real(dp), intent(in), optional :: loc !! support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale (default 1)
        real(dp) :: score(2)
        integer :: i
        real(dp) :: gz, mu, sigma, z

        if (.not. valid_loglike(irwinhall_loglikelihood(data, n, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= n) then
                score = nan_value()
                return
            end if
            gz = irwinhall_logpdf_derivative(z, n)
            if (.not. ieee_is_finite(gz)) then
                score = nan_value()
                return
            end if
            score(1) = score(1) - gz / sigma
            score(2) = score(2) - (1.0_dp + z * gz) / sigma
        end do
    end function irwinhall_score

    pure function ksone_score(data, n, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: score(2)
        integer :: i
        real(dp) :: gz, mu, sigma, z

        if (.not. valid_loglike(ksone_loglikelihood(data, n, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp .or. z >= 1.0_dp) then
                score = nan_value()
                return
            end if
            gz = ksone_logpdf_derivative(z, n)
            if (.not. ieee_is_finite(gz)) then
                score = nan_value()
                return
            end if
            score(1) = score(1) - gz / sigma
            score(2) = score(2) - (1.0_dp + z * gz) / sigma
        end do
    end function ksone_score

    pure function kstwo_score(data, n, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent finite-sample two-sided KS observations
        real(dp), intent(in) :: n !! positive integer sample size
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive support scale (default 1)
        real(dp) :: score(2)
        integer :: i
        real(dp) :: gz, mu, sigma, z

        if (.not. valid_loglike(kstwo_loglikelihood(data, n, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.5_dp / n .or. z >= 1.0_dp) then
                score = nan_value()
                return
            end if
            gz = kstwo_logpdf_derivative(z, n)
            if (.not. ieee_is_finite(gz)) then
                score = nan_value()
                return
            end if
            score(1) = score(1) - gz / sigma
            score(2) = score(2) - (1.0_dp + z * gz) / sigma
        end do
    end function kstwo_score

    pure function ncx2_score(data, df, nc, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! positive degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        integer :: i
        real(dp) :: mu, sigma, z, logf, ddf, dnc, gz

        if (.not. valid_loglike(ncx2_loglikelihood(data, df, nc, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            call ncx2_logpdf_derivatives(z, df, nc, logf, ddf, dnc, gz)
            if (.not. ieee_is_finite(logf) .or. .not. ieee_is_finite(gz)) then
                score = nan_value()
                return
            end if
            score(1) = score(1) + ddf
            score(2) = score(2) + dnc
            score(3) = score(3) - gz / sigma
            score(4) = score(4) - (1.0_dp + z * gz) / sigma
        end do
    end function ncx2_score

    pure function ncf_score(data, dfn, dfd, nc, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: dfn !! positive numerator degrees of freedom
        real(dp), intent(in) :: dfd !! positive denominator degrees of freedom
        real(dp), intent(in) :: nc !! nonnegative noncentrality parameter
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(5)
        integer :: i
        real(dp) :: mu, sigma, z, logf, dd1, dd2, dnc, gz

        if (.not. valid_loglike(ncf_loglikelihood(data, dfn, dfd, nc, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            if (z <= 0.0_dp) then
                score = nan_value()
                return
            end if
            call ncf_logpdf_derivatives(z, dfn, dfd, nc, logf, dd1, dd2, dnc, gz)
            if (.not. ieee_is_finite(logf) .or. .not. ieee_is_finite(gz)) then
                score = nan_value()
                return
            end if
            score(1) = score(1) + dd1
            score(2) = score(2) + dd2
            score(3) = score(3) + dnc
            score(4) = score(4) - gz / sigma
            score(5) = score(5) - (1.0_dp + z * gz) / sigma
        end do
    end function ncf_score


    pure function randint_score_real(data, low, high, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: low !! integer lower shape bound, inclusive
        real(dp), intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(0)

        score = [real(dp) ::]
        ! randint has no differentiable shape parameter. Evaluating the
        ! likelihood here still validates the family/observation combination.
        if (.not. valid_loglike(randint_loglikelihood(data, low, high, loc))) return
    end function randint_score_real

    pure function randint_score_int(data, low, high, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: low !! integer lower shape bound, inclusive
        integer, intent(in) :: high !! integer upper shape bound, exclusive
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(0)
        score = randint_score_real(real(data, dp), real(low, dp), real(high, dp), loc)
    end function randint_score_int

    pure function planck_score_real(data, lambda, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(planck_loglikelihood(data, lambda, loc)) .or. lambda <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = real(size(data), dp) / expm1_safe(lambda) - sum(count)
    end function planck_score_real

    pure function planck_score_int(data, lambda, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = planck_score_real(real(data, dp), lambda, loc)
    end function planck_score_int

    pure function dlaplace_score_real(data, a, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(dlaplace_loglikelihood(data, a, loc)) .or. a <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = real(size(data), dp) / sinh(a) - sum(abs(count))
    end function dlaplace_score_real

    pure function dlaplace_score_int(data, a, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! positive discrete-Laplace rate
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = dlaplace_score_real(real(data, dp), a, loc)
    end function dlaplace_score_int

    pure function logser_score_real(data, p, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset, norm

        if (.not. valid_loglike(logser_loglikelihood(data, p, loc)) .or. &
            p <= 0.0_dp .or. p >= 1.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        norm = -log1p_safe(-p)
        score(1) = sum(count) / p - real(size(data), dp) / ((1.0_dp - p) * norm)
    end function logser_score_real

    pure function logser_score_int(data, p, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: p !! probability parameter strictly between 0 and 1
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = logser_score_real(real(data, dp), p, loc)
    end function logser_score_int


    pure function betabinom_score_real(data, n, a, b, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: n !! nonnegative integer number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(2)
        real(dp) :: count(size(data))
        real(dp) :: common
        real(dp) :: offset

        if (.not. valid_loglike(betabinom_loglikelihood(data, n, a, b, loc)) .or. &
            a <= 0.0_dp .or. b <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        common = -real(size(data), dp) * digamma_positive(n + a + b) + &
            real(size(data), dp) * digamma_positive(a + b)
        score(1) = sum(digamma_positive(count + a)) + common - &
            real(size(data), dp) * digamma_positive(a)
        score(2) = sum(digamma_positive(n - count + b)) + common - &
            real(size(data), dp) * digamma_positive(b)
    end function betabinom_score_real

    pure function betabinom_score_int(data, n, a, b, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: n !! nonnegative number of trials
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(2)
        score = betabinom_score_real(real(data, dp), real(n, dp), a, b, loc)
    end function betabinom_score_int

    pure function hypergeom_score_real(data, m, n, draws, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer Type-I count
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(0)
        score = [real(dp) ::]
        if (.not. valid_loglike(hypergeom_loglikelihood(data, m, n, draws, loc))) return
    end function hypergeom_score_real

    pure function hypergeom_score_int(data, m, n, draws, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! Type-I count
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(0)
        score = hypergeom_score_real(real(data, dp), real(m, dp), real(n, dp), real(draws, dp), loc)
    end function hypergeom_score_int

    pure function nhypergeom_score_real(data, m, n, r, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! nonnegative integer population size
        real(dp), intent(in) :: n !! integer success count
        real(dp), intent(in) :: r !! integer failures required
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(0)
        score = [real(dp) ::]
        if (.not. valid_loglike(nhypergeom_loglikelihood(data, m, n, r, loc))) return
    end function nhypergeom_score_real

    pure function nhypergeom_score_int(data, m, n, r, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! population size
        integer, intent(in) :: n !! success count
        integer, intent(in) :: r !! failures required
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(0)
        score = nhypergeom_score_real(real(data, dp), real(m, dp), real(n, dp), real(r, dp), loc)
    end function nhypergeom_score_int

    pure function boltzmann_score_real(data, lambda, n, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        real(dp), intent(in) :: n !! positive integer support size
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(boltzmann_loglikelihood(data, lambda, n, loc)) .or. lambda <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = real(size(data), dp) / expm1_safe(lambda) - sum(count) - &
            real(size(data), dp) * n / expm1_safe(lambda * n)
    end function boltzmann_score_real

    pure function boltzmann_score_int(data, lambda, n, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: lambda !! positive exponential rate
        integer, intent(in) :: n !! positive support size
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = boltzmann_score_real(real(data, dp), lambda, real(n, dp), loc)
    end function boltzmann_score_int


    pure function betanbinom_score_real(data, n, a, b, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(2)
        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(betanbinom_loglikelihood(data, n, a, b, loc)) .or. &
            a <= 0.0_dp .or. b <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = real(size(data), dp) * (digamma_positive(a + n) - digamma_positive(a) + &
            digamma_positive(a + b)) - sum(digamma_positive(a + b + n + count))
        score(2) = sum(digamma_positive(b + count) - digamma_positive(a + b + n + count)) - &
            real(size(data), dp) * (digamma_positive(b) - digamma_positive(a + b))
    end function betanbinom_score_real

    pure function betanbinom_score_int(data, n, a, b, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: n !! positive integer success count
        real(dp), intent(in) :: a !! positive first beta shape
        real(dp), intent(in) :: b !! positive second beta shape
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(2)
        score = betanbinom_score_real(real(data, dp), real(n, dp), a, b, loc)
    end function betanbinom_score_int

    pure function yulesimon_score_real(data, alpha, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset

        if (.not. valid_loglike(yulesimon_loglikelihood(data, alpha, loc)) .or. alpha <= 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        score(1) = real(size(data), dp) * (1.0_dp / alpha + digamma_positive(alpha + 1.0_dp)) - &
            sum(digamma_positive(count + alpha + 1.0_dp))
    end function yulesimon_score_real

    pure function yulesimon_score_int(data, alpha, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: alpha !! positive Yule-Simon shape
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = yulesimon_score_real(real(data, dp), alpha, loc)
    end function yulesimon_score_int

    pure function zipf_score_real(data, a, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset, z, dz

        if (.not. valid_loglike(zipf_loglikelihood(data, a, loc)) .or. a <= 1.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        z = hurwitz_zeta(a, 1.0_dp)
        dz = hurwitz_zeta_derivative(a, 1.0_dp)
        score(1) = -sum(log(count)) - real(size(data), dp) * dz / z
    end function zipf_score_real

    pure function zipf_score_int(data, a, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! power exponent, strictly greater than one
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = zipf_score_real(real(data, dp), a, loc)
    end function zipf_score_int

    pure function zipfian_score_real(data, a, n, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: a !! nonnegative power exponent
        real(dp), intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        real(dp) :: count(size(data))
        real(dp) :: offset, h, hlog, weight
        integer :: j

        if (.not. valid_loglike(zipfian_loglikelihood(data, a, n, loc)) .or. a < 0.0_dp) then
            score = nan_value()
            return
        end if
        offset = optional_loc(loc)
        count = data - offset
        h = 0.0_dp
        hlog = 0.0_dp
        do j = 1, int(n)
            weight = exp(-a * log(real(j, dp)))
            h = h + weight
            hlog = hlog + weight * log(real(j, dp))
        end do
        score(1) = -sum(log(count)) + real(size(data), dp) * hlog / h
    end function zipfian_score_real

    pure function zipfian_score_int(data, a, n, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: a !! nonnegative power exponent
        integer, intent(in) :: n !! positive integer upper support bound
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(1)
        score = zipfian_score_real(real(data, dp), a, real(n, dp), loc)
    end function zipfian_score_int

    pure function geninvgauss_score(data, p, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: p !! real generalized-inverse-Gaussian shape
        real(dp), intent(in) :: b !! strictly positive generalized-inverse-Gaussian shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma, z, logf, dp_shape, db_shape, dz
        integer :: i
        if (.not. valid_loglike(geninvgauss_loglikelihood(data,p,b,loc,scale))) then
            score = nan_value(); return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i=1,size(data)
            z=(data(i)-mu)/sigma
            call geninvgauss_logpdf_derivatives(z,p,b,logf,dp_shape,db_shape,dz)
            score(1)=score(1)+dp_shape
            score(2)=score(2)+db_shape
            score(3)=score(3)-dz/sigma
            score(4)=score(4)-(1.0_dp+z*dz)/sigma
        end do
    end function geninvgauss_score

    pure function norminvgauss_score(data, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive normal-inverse-Gaussian tail shape
        real(dp), intent(in) :: b !! skew shape satisfying abs(b) < a
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma, z, logf, da, db, dz
        integer :: i
        if (.not. valid_loglike(norminvgauss_loglikelihood(data,a,b,loc,scale))) then
            score = nan_value(); return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i=1,size(data)
            z=(data(i)-mu)/sigma
            call norminvgauss_logpdf_derivatives(z,a,b,logf,da,db,dz)
            score(1)=score(1)+da
            score(2)=score(2)+db
            score(3)=score(3)-dz/sigma
            score(4)=score(4)-(1.0_dp+z*dz)/sigma
        end do
    end function norminvgauss_score

    pure function skellam_score_real(data, mu1, mu2, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(2)
        real(dp) :: shift, count, logf, d1, d2
        integer :: i
        if (.not. valid_loglike(skellam_loglikelihood(data,mu1,mu2,loc))) then
            score=nan_value(); return
        end if
        shift=optional_loc(loc); score=0.0_dp
        do i=1,size(data)
            count=data(i)-shift
            call skellam_logpmf_derivatives(count,mu1,mu2,logf,d1,d2)
            score(1)=score(1)+d1
            score(2)=score(2)+d2
        end do
    end function skellam_score_real

    pure function skellam_score_int(data, mu1, mu2, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: mu1 !! positive first Poisson mean
        real(dp), intent(in) :: mu2 !! positive second Poisson mean
        real(dp), intent(in), optional :: loc !! lattice shift (default 0)
        real(dp) :: score(2)
        score=skellam_score_real(real(data,dp),mu1,mu2,loc)
    end function skellam_score_int


    pure function genhyperbolic_score(data, p, a, b, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: p !! real generalized-hyperbolic tail shape
        real(dp), intent(in) :: a !! positive generalized-hyperbolic shape
        real(dp), intent(in) :: b !! generalized-hyperbolic skew shape
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(5)
        real(dp) :: mu, sigma, z, logf, dp_shape, da, db, dz
        integer :: i
        if (.not. valid_loglike(genhyperbolic_loglikelihood(data,p,a,b,loc,scale))) then
            score = nan_value(); return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            call genhyperbolic_logpdf_derivatives(z,p,a,b,logf,dp_shape,da,db,dz)
            if (.not. ieee_is_finite(dp_shape) .or. .not. ieee_is_finite(da) .or. &
                .not. ieee_is_finite(db) .or. .not. ieee_is_finite(dz)) then
                score = nan_value(); return
            end if
            score(1) = score(1) + dp_shape
            score(2) = score(2) + da
            score(3) = score(3) + db
            score(4) = score(4) - dz / sigma
            score(5) = score(5) - (1.0_dp + z * dz) / sigma
        end do
    end function genhyperbolic_score

    pure function nchypergeom_fisher_score_real(data, m, n, draws, odds, loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: score(1)
        real(dp) :: shift, count, logf, dodds
        integer :: i
        if (.not. valid_loglike(nchypergeom_fisher_loglikelihood(data,m,n,draws,odds,loc))) then
            score = nan_value(); return
        end if
        shift = optional_loc(loc); score = 0.0_dp
        do i = 1, size(data)
            count = data(i) - shift
            call nchypergeom_fisher_logpmf_derivative_odds(count,m,n,draws,odds,logf,dodds)
            if (.not. ieee_is_finite(dodds)) then
                score = nan_value(); return
            end if
            score(1) = score(1) + dodds
        end do
    end function nchypergeom_fisher_score_real

    pure function nchypergeom_fisher_score_int(data, m, n, draws, odds, loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: score(1)
        score = nchypergeom_fisher_score_real(real(data,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_fisher_score_int

    pure function nct_score(data, df, nc, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: df !! positive noncentral-t degrees of freedom
        real(dp), intent(in) :: nc !! finite noncentrality parameter
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma, z, logf, ddf, dnc, dz
        integer :: i
        if (.not. valid_loglike(nct_loglikelihood(data,df,nc,loc,scale))) then
            score = nan_value(); return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            z = (data(i) - mu) / sigma
            call nct_logpdf_derivatives(z,df,nc,logf,ddf,dnc,dz)
            if (.not. ieee_is_finite(ddf) .or. .not. ieee_is_finite(dnc) .or. &
                .not. ieee_is_finite(dz)) then
                score = nan_value(); return
            end if
            score(1) = score(1) + ddf
            score(2) = score(2) + dnc
            score(3) = score(3) - dz / sigma
            score(4) = score(4) - (1.0_dp + z * dz) / sigma
        end do
    end function nct_score

    pure function gausshyper_score(data, a, b, c, zshape, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent observations
        real(dp), intent(in) :: a !! positive first beta-shape parameter
        real(dp), intent(in) :: b !! positive second beta-shape parameter
        real(dp), intent(in) :: c !! finite hypergeometric tilt exponent
        real(dp), intent(in) :: zshape !! finite tilt parameter greater than -1
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(6)
        real(dp) :: mu, sigma, q, logf, da, db, dc, dzshape, dq
        integer :: i
        if (.not. valid_loglike(gausshyper_loglikelihood(data,a,b,c,zshape,loc,scale))) then
            score = nan_value(); return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        score = 0.0_dp
        do i = 1, size(data)
            q = (data(i) - mu) / sigma
            call gausshyper_logpdf_derivatives(q,a,b,c,zshape,logf,da,db,dc,dzshape,dq)
            if (.not. ieee_is_finite(da) .or. .not. ieee_is_finite(db) .or. &
                .not. ieee_is_finite(dc) .or. .not. ieee_is_finite(dzshape) .or. &
                .not. ieee_is_finite(dq)) then
                score = nan_value(); return
            end if
            score(1) = score(1) + da
            score(2) = score(2) + db
            score(3) = score(3) + dc
            score(4) = score(4) + dzshape
            score(5) = score(5) - dq / sigma
            score(6) = score(6) - (1.0_dp + q * dq) / sigma
        end do
    end function gausshyper_score

    pure function nchypergeom_wallenius_score_real(data,m,n,draws,odds,loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: m !! positive integer population size
        real(dp), intent(in) :: n !! integer number of Type-I objects
        real(dp), intent(in) :: draws !! integer sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: score(1)
        real(dp) :: shift, count, logf, dodds
        integer :: i
        if (.not.valid_loglike(nchypergeom_wallenius_loglikelihood(data,m,n,draws,odds,loc))) then
            score=nan_value(); return
        end if
        shift=optional_loc(loc); score=0.0_dp
        do i=1,size(data)
            count=data(i)-shift
            call nchypergeom_wallenius_logpmf_derivative_odds(count,m,n,draws,odds,logf,dodds)
            if (.not.ieee_is_finite(dodds)) then; score=nan_value(); return; end if
            score(1)=score(1)+dodds
        end do
    end function nchypergeom_wallenius_score_real

    pure function nchypergeom_wallenius_score_int(data,m,n,draws,odds,loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        integer, intent(in) :: m !! positive population size
        integer, intent(in) :: n !! number of Type-I objects
        integer, intent(in) :: draws !! sample size
        real(dp), intent(in) :: odds !! strictly positive Type-I odds ratio
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: score(1)
        score=nchypergeom_wallenius_score_real(real(data,dp),real(m,dp),real(n,dp),real(draws,dp),odds,loc)
    end function nchypergeom_wallenius_score_int

    pure function poisson_binom_score_real(data,p,loc) result(score)
        real(dp), intent(in) :: data(:) !! independent lattice observations
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: score(size(p))
        real(dp) :: tmp(size(p)), logf, shift
        integer :: i
        if (.not.valid_loglike(poisson_binom_loglikelihood(data,p,loc))) then
            score=nan_value(); return
        end if
        shift=optional_loc(loc); score=0.0_dp
        do i=1,size(data)
            call poisson_binom_logpmf_derivative_p(data(i)-shift,p,logf,tmp)
            if (any(.not.ieee_is_finite(tmp))) then; score=nan_value(); return; end if
            score=score+tmp
        end do
    end function poisson_binom_score_real

    pure function poisson_binom_score_int(data,p,loc) result(score)
        integer, intent(in) :: data(:) !! independent integer observations
        real(dp), intent(in) :: p(:) !! Bernoulli success probabilities
        real(dp), intent(in), optional :: loc !! support shift (default 0)
        real(dp) :: score(size(p))
        score=poisson_binom_score_real(real(data,dp),p,loc)
    end function poisson_binom_score_int

    pure function levy_stable_score(data, alpha, beta, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent Levy-stable observations
        real(dp), intent(in) :: alpha !! stability shape parameter in (0,2]
        real(dp), intent(in) :: beta !! skewness shape parameter in [-1,1]
        real(dp), intent(in), optional :: loc !! S1 location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma, ha, hb, hl, hs, fp, fm, f0

        if (.not. valid_loglike(levy_stable_loglikelihood(data, alpha, beta, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        f0 = levy_stable_loglikelihood(data, alpha, beta, mu, sigma)

        ha = 2.0e-5_dp * max(1.0_dp, abs(alpha))
        if (alpha - ha > 0.0_dp .and. alpha + ha <= 2.0_dp .and. &
            .not. ((alpha - ha < 1.0_dp) .and. (alpha + ha > 1.0_dp))) then
            fp = levy_stable_loglikelihood(data, alpha + ha, beta, mu, sigma)
            fm = levy_stable_loglikelihood(data, alpha - ha, beta, mu, sigma)
            score(1) = (fp - fm) / (2.0_dp * ha)
        else if (alpha + ha <= 2.0_dp) then
            fp = levy_stable_loglikelihood(data, alpha + ha, beta, mu, sigma)
            score(1) = (fp - f0) / ha
        else if (alpha - ha > 0.0_dp) then
            fm = levy_stable_loglikelihood(data, alpha - ha, beta, mu, sigma)
            score(1) = (f0 - fm) / ha
        else
            score = nan_value()
            return
        end if

        hb = 2.0e-5_dp
        if (beta - hb >= -1.0_dp .and. beta + hb <= 1.0_dp) then
            fp = levy_stable_loglikelihood(data, alpha, beta + hb, mu, sigma)
            fm = levy_stable_loglikelihood(data, alpha, beta - hb, mu, sigma)
            score(2) = (fp - fm) / (2.0_dp * hb)
        else if (beta + hb <= 1.0_dp) then
            fp = levy_stable_loglikelihood(data, alpha, beta + hb, mu, sigma)
            score(2) = (fp - f0) / hb
        else
            fm = levy_stable_loglikelihood(data, alpha, beta - hb, mu, sigma)
            score(2) = (f0 - fm) / hb
        end if

        hl = 2.0e-6_dp * max(1.0_dp, sigma)
        fp = levy_stable_loglikelihood(data, alpha, beta, mu + hl, sigma)
        fm = levy_stable_loglikelihood(data, alpha, beta, mu - hl, sigma)
        score(3) = (fp - fm) / (2.0_dp * hl)

        hs = 2.0e-6_dp * sigma
        fp = levy_stable_loglikelihood(data, alpha, beta, mu, sigma + hs)
        fm = levy_stable_loglikelihood(data, alpha, beta, mu, sigma - hs)
        score(4) = (fp - fm) / (2.0_dp * hs)
        if (any(.not. ieee_is_finite(score))) score = nan_value()
    end function levy_stable_score

    pure function studentized_range_score(data, k, df, loc, scale) result(score)
        real(dp), intent(in) :: data(:) !! independent studentized-range observations
        real(dp), intent(in) :: k !! number-of-means shape parameter, > 1
        real(dp), intent(in) :: df !! degrees of freedom, > 0
        real(dp), intent(in), optional :: loc !! lower support location (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: score(4)
        real(dp) :: mu, sigma, hk, hd, hl, hs, fp, fm, f0

        if (.not. valid_loglike(studentized_range_loglikelihood(data, k, df, loc, scale))) then
            score = nan_value()
            return
        end if
        call get_loc_scale(loc, scale, mu, sigma)
        f0 = studentized_range_loglikelihood(data, k, df, mu, sigma)

        hk = 2.0e-5_dp * max(1.0_dp, abs(k))
        if (k - hk > 1.0_dp) then
            fp = studentized_range_loglikelihood(data, k + hk, df, mu, sigma)
            fm = studentized_range_loglikelihood(data, k - hk, df, mu, sigma)
            score(1) = (fp - fm) / (2.0_dp * hk)
        else
            fp = studentized_range_loglikelihood(data, k + hk, df, mu, sigma)
            score(1) = (fp - f0) / hk
        end if

        hd = 2.0e-5_dp * max(1.0_dp, abs(df))
        if (df - hd > 0.0_dp) then
            fp = studentized_range_loglikelihood(data, k, df + hd, mu, sigma)
            fm = studentized_range_loglikelihood(data, k, df - hd, mu, sigma)
            score(2) = (fp - fm) / (2.0_dp * hd)
        else
            fp = studentized_range_loglikelihood(data, k, df + hd, mu, sigma)
            score(2) = (fp - f0) / hd
        end if

        hl = 2.0e-6_dp * max(1.0_dp, sigma)
        fp = studentized_range_loglikelihood(data, k, df, mu + hl, sigma)
        fm = studentized_range_loglikelihood(data, k, df, mu - hl, sigma)
        score(3) = (fp - fm) / (2.0_dp * hl)

        hs = 2.0e-6_dp * sigma
        fp = studentized_range_loglikelihood(data, k, df, mu, sigma + hs)
        fm = studentized_range_loglikelihood(data, k, df, mu, sigma - hs)
        score(4) = (fp - fm) / (2.0_dp * hs)
        if (any(.not. ieee_is_finite(score))) score = nan_value()
    end function studentized_range_score

end module scifort_score
