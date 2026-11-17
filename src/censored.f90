! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent implementation of the censored likelihood specified in
! https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html
module scifort_censored
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_log_sqrt_two_pi, scifort_log_two, scifort_pi
    use scifort_empirical_distribution, only : censored_data
    use scifort_math, only : quiet_nan, negative_infinity, positive_infinity, valid_loc_scale, &
        expm1_safe, log1p_safe, log1pexp
    use scifort_normal, only : normal_logpdf, normal_logcdf, normal_logsf
    use scifort_exponential, only : exponential_logpdf, exponential_logcdf, exponential_logsf
    use scifort_log_gamma, only : log_gamma_one_plus, stirling_remainder, log_beta, &
        stirling_threshold, log_gamma_difference, log_gamma_ratio_scaled
    use scifort_incomplete_beta, only : incomplete_beta_xy, log_beta_kernel
    use scifort_gamma, only : gamma_logpdf
    use scifort_exponential_integrals, only : exp1
    use scifort_incomplete_gamma, only : log_gammainc, log_gammaincc, log_gamma_kernel
    use scifort_weibull, only : weibull_logpdf, weibull_logcdf, weibull_logsf
    implicit none
    private

    public :: censored_data
    public :: chi2_censored_loglikelihood, chi2_censored_nnlf
    public :: chi_censored_loglikelihood, chi_censored_nnlf
    public :: nakagami_censored_loglikelihood, nakagami_censored_nnlf
    public :: gennorm_censored_loglikelihood, gennorm_censored_nnlf
    public :: halfgennorm_censored_loglikelihood, halfgennorm_censored_nnlf
    public :: rdist_censored_loglikelihood, rdist_censored_nnlf
    public :: foldnorm_censored_loglikelihood, foldnorm_censored_nnlf
    public :: foldcauchy_censored_loglikelihood, foldcauchy_censored_nnlf
    public :: genpareto_censored_loglikelihood, genpareto_censored_nnlf
    public :: loguniform_censored_loglikelihood, loguniform_censored_nnlf
    public :: truncpareto_censored_loglikelihood, truncpareto_censored_nnlf
    public :: genhalflogistic_censored_loglikelihood, genhalflogistic_censored_nnlf
    public :: truncexpon_censored_loglikelihood, truncexpon_censored_nnlf
    public :: johnsonsb_censored_loglikelihood, johnsonsb_censored_nnlf
    public :: kappa3_censored_loglikelihood, kappa3_censored_nnlf
    public :: cosine_censored_loglikelihood, cosine_censored_nnlf
    public :: anglit_censored_loglikelihood, anglit_censored_nnlf
    public :: semicircular_censored_loglikelihood, semicircular_censored_nnlf
    public :: fatiguelife_censored_loglikelihood, fatiguelife_censored_nnlf
    public :: johnsonsu_censored_loglikelihood, johnsonsu_censored_nnlf
    public :: genexpon_censored_loglikelihood, genexpon_censored_nnlf
    public :: exponpow_censored_loglikelihood, exponpow_censored_nnlf
    public :: exponweib_censored_loglikelihood, exponweib_censored_nnlf
    public :: mielke_censored_loglikelihood, mielke_censored_nnlf
    public :: powernorm_censored_loglikelihood, powernorm_censored_nnlf
    public :: powerlognorm_censored_loglikelihood, powerlognorm_censored_nnlf
    public :: skewcauchy_censored_loglikelihood, skewcauchy_censored_nnlf
    public :: invgamma_censored_loglikelihood, invgamma_censored_nnlf
    public :: betaprime_censored_loglikelihood, betaprime_censored_nnlf
    public :: maxwell_censored_loglikelihood, maxwell_censored_nnlf
    public :: dgamma_censored_loglikelihood, dgamma_censored_nnlf
    public :: dweibull_censored_loglikelihood, dweibull_censored_nnlf
    public :: loggamma_censored_loglikelihood, loggamma_censored_nnlf
    public :: bradford_censored_loglikelihood, bradford_censored_nnlf
    public :: triang_censored_loglikelihood, triang_censored_nnlf
    public :: trapezoid_censored_loglikelihood, trapezoid_censored_nnlf
    public :: moyal_censored_loglikelihood, moyal_censored_nnlf
    public :: levy_censored_loglikelihood, levy_censored_nnlf
    public :: levy_l_censored_loglikelihood, levy_l_censored_nnlf
    public :: genlogistic_censored_loglikelihood, genlogistic_censored_nnlf
    public :: gompertz_censored_loglikelihood, gompertz_censored_nnlf
    public :: invweibull_censored_loglikelihood, invweibull_censored_nnlf
    public :: powerlaw_censored_loglikelihood, powerlaw_censored_nnlf
    public :: arcsine_censored_loglikelihood, arcsine_censored_nnlf
    public :: uniform_censored_loglikelihood, uniform_censored_nnlf
    public :: burr12_censored_loglikelihood, burr12_censored_nnlf
    public :: burr_censored_loglikelihood, burr_censored_nnlf
    public :: fisk_censored_loglikelihood, fisk_censored_nnlf
    public :: lomax_censored_loglikelihood, lomax_censored_nnlf
    public :: pareto_censored_loglikelihood, pareto_censored_nnlf
    public :: loglaplace_censored_loglikelihood, loglaplace_censored_nnlf

    integer, parameter, public :: censored_status_success = 0
    integer, parameter, public :: censored_status_invalid_input = 1

    abstract interface
        pure function censored_log_function(x, params) result(value)
            import dp
            implicit none
            real(dp), intent(in) :: x !! evaluation point, possibly infinite
            real(dp), intent(in) :: params(:) !! finite distribution parameters
            real(dp) :: value
        end function censored_log_function
        pure function censored_log_interval(a, b, params) result(value)
            import dp
            implicit none
            real(dp), intent(in) :: a !! lower endpoint, strictly less than b
            real(dp), intent(in) :: b !! upper endpoint, possibly infinite
            real(dp), intent(in) :: params(:) !! finite distribution parameters
            real(dp) :: value
        end function censored_log_interval
    end interface

    public :: make_censored_data, censored_data_valid, censored_size, censored_loglikelihood
    public :: censored_log_function, censored_log_interval
    public :: normal_censored_loglikelihood, normal_censored_nnlf
    public :: exponential_censored_loglikelihood, exponential_censored_nnlf
    public :: weibull_censored_loglikelihood, weibull_censored_nnlf
    public :: lognormal_censored_loglikelihood, lognormal_censored_nnlf
    public :: gibrat_censored_loglikelihood, gibrat_censored_nnlf
    public :: gamma_censored_loglikelihood, gamma_censored_nnlf
    public :: erlang_censored_loglikelihood, erlang_censored_nnlf
    public :: beta_censored_loglikelihood, beta_censored_nnlf
    public :: rayleigh_censored_loglikelihood, rayleigh_censored_nnlf
    public :: weibull_max_censored_loglikelihood, weibull_max_censored_nnlf
    public :: halfcauchy_censored_loglikelihood, halfcauchy_censored_nnlf
    public :: halfnorm_censored_loglikelihood, halfnorm_censored_nnlf
    public :: hypsecant_censored_loglikelihood, hypsecant_censored_nnlf
    public :: laplace_asymmetric_censored_loglikelihood, laplace_asymmetric_censored_nnlf
    public :: halflogistic_censored_loglikelihood, halflogistic_censored_nnlf
    public :: cauchy_censored_loglikelihood, cauchy_censored_nnlf
    public :: laplace_censored_loglikelihood, laplace_censored_nnlf
    public :: logistic_censored_loglikelihood, logistic_censored_nnlf
    public :: gumbel_r_censored_loglikelihood, gumbel_r_censored_nnlf
    public :: gumbel_l_censored_loglikelihood, gumbel_l_censored_nnlf

contains

    ! Gamma-power laws: Y = rate*((x-loc)/scale)**power ~ Gamma(a,1).
    ! Derived from the public density specifications; see CODE_PROVENANCE.md.
    pure function gamma_power_logy(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value
        value = params(5)*lognormal_logz(x,params(1:3))+params(4)
    end function gamma_power_logy

    pure function gamma_power_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(1) < 1.0_dp/params(5)) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp/params(5)) value = log(params(5)) &
                +params(1)*params(4)-log_gamma(params(1))-log(params(3))
            return
        end if
        logy = gamma_power_logy(x,params)
        value = loggamma_lpdf(logy,[params(1),0.0_dp,1.0_dp]) &
            +log(params(5))-log_distance(x,params(2))
    end function gamma_power_lpdf

    pure function gamma_power_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp) then
            value = log_difference(0.0_dp,gamma_power_lsf(x,params))
            return
        end if
        value = loggamma_lcdf(gamma_power_logy(x,params),[params(1),0.0_dp,1.0_dp])
    end function gamma_power_lcdf

    pure function gamma_power_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value, logy
        real(dp), parameter :: euler_gamma = 0.577215664901532860606512090082402431_dp
        value = 0.0_dp
        if (x <= params(2)) return
        logy = gamma_power_logy(x,params)
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp .and. logy <= log(1.5_dp)) then
            ! Q(a,Y)=a*E1(Y)*(1+O(a*(1+abs(log(Y))))). Keep log(a).
            if (logy < log(tiny(1.0_dp))) then
                value = log(params(1))+log(-logy-euler_gamma)
            else
                value = log(params(1))+log(exp1(exp(logy)))
            end if
        else
            value = loggamma_lsf(logy,[params(1),0.0_dp,1.0_dp])
        end if
    end function gamma_power_lsf

    pure function gamma_power_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! gamma shape, loc, scale, log rate, power
        real(dp) :: value, logwidth, logq, logy, q, t(3), d(3), r(3)
        real(dp) :: first(3), second(3), correction, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = gamma_power_lcdf(b,params)
            return
        end if
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(2))
        logy = gamma_power_logy(a,params)
        ! Bound relative width and both terms of the log-density variation.
        if (logq+max(0.0_dp,log(params(5))+log(params(1)), &
                log(params(5))+logy) <= log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            d = log1p_safe(t*q)
            first = t
            second = params(5)*t
            if (q > 0.0_dp) then
                first = d/q
                second = expm1_safe(params(5)*d)/q
            end if
            r = exp(log(params(1))+log(params(5))+logq)*first-d &
                -exp(logy+logq)*second
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = gamma_power_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,gamma_power_lcdf,gamma_power_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(gamma_power_lcdf(a,params),gamma_power_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function gamma_power_interval

    pure function chi2_censored_loglikelihood(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(df)) return
        if (df <= 0.0_dp) return
        params = [0.5_dp*df,0.0_dp,1.0_dp,-scifort_log_two,1.0_dp]
        if (params(1) > 1.0e6_dp) return
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gamma_power_lpdf,gamma_power_lcdf, &
            gamma_power_lsf,gamma_power_interval)
    end function chi2_censored_loglikelihood

    pure function chi2_censored_nnlf(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -chi2_censored_loglikelihood(data,df,loc,scale)
    end function chi2_censored_nnlf

    pure function chi_censored_loglikelihood(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(df)) return
        if (df <= 0.0_dp) return
        params = [0.5_dp*df,0.0_dp,1.0_dp,-scifort_log_two,2.0_dp]
        if (params(1) > 1.0e6_dp) return
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gamma_power_lpdf,gamma_power_lcdf, &
            gamma_power_lsf,gamma_power_interval)
    end function chi_censored_loglikelihood

    pure function chi_censored_nnlf(data, df, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: df !! degrees of freedom in (0,2e6], with positive fl(df/2)
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -chi_censored_loglikelihood(data,df,loc,scale)
    end function chi_censored_nnlf

    pure function nakagami_censored_loglikelihood(data, nu, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: nu !! Nakagami shape in (0,1e6]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(nu)) return
        if (nu <= 0.0_dp) return
        params = [nu,0.0_dp,1.0_dp,log(nu),2.0_dp]
        if (params(1) > 1.0e6_dp) return
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        if (.not. all(ieee_is_finite(params))) return
        if (params(1) <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gamma_power_lpdf,gamma_power_lcdf, &
            gamma_power_lsf,gamma_power_interval)
    end function nakagami_censored_loglikelihood

    pure function nakagami_censored_nnlf(data, nu, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: nu !! Nakagami shape in (0,1e6]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -nakagami_censored_loglikelihood(data,nu,loc,scale)
    end function nakagami_censored_nnlf


    pure function rdist_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, left, right, z, w, lz, lw
        if (params(1) > epsilon(1.0_dp)/4096.0_dp) then
            value = beta_lpdf(x,params)
            return
        end if
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        value = positive_infinity(0.0_dp)
        if (left == 0.0_dp .or. right == 0.0_dp) return
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        value = log(params(1))-scifort_log_two-log(params(4))-lz-lw
    end function rdist_lpdf

    pure subroutine rdist_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: left, right
        call beta_distances(x,params,left,right)
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (left <= 0.0_dp) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        if (params(1) <= epsilon(1.0_dp)/4096.0_dp .or. left == right) then
            lp = -scifort_log_two
            lq = lp
        else
            call beta_logtails(x,params,lp,lq)
        end if
    end subroutine rdist_tails

    pure function rdist_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, unused
        call rdist_tails(x,params,value,unused)
    end function rdist_lcdf

    pure function rdist_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, unused
        call rdist_tails(x,params,unused,value)
    end function rdist_lsf

    pure function rdist_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! equal half-shapes, rounded lower edge, full width
        real(dp) :: value, left, right, lower, upper, span
        if (params(1) > epsilon(1.0_dp)/4096.0_dp) then
            value = beta_interval(a,b,params)
            return
        end if
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = rdist_lsf(a,params)
            return
        end if
        upper = log(right)
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = rdist_lcdf(b,params)
            return
        end if
        lower = log_distance(a,params(3))
        span = log_distance(b,a)
        value = log(params(1))-scifort_log_two &
            +bounded_logadd(burr_logsoftplus(span-lower),burr_logsoftplus(span-upper))
    end function rdist_interval

    pure function generalized_normal_norm(a) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive finite reciprocal power
        real(dp) :: value, product
        if (a < stirling_threshold) then
            value = -log_gamma_one_plus(a)
        else
            product = bounded_cutoff_product(a,log(a)-1.0_dp)
            value = negative_infinity(0.0_dp)
            if (.not. ieee_is_finite(product)) return
            value = -product-0.5_dp*log(a)-scifort_log_sqrt_two_pi-stirling_remainder(a)
        end if
    end function generalized_normal_norm

    pure function generalized_normal_power(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, lz
        value = negative_infinity(0.0_dp)
        if (x == params(2)) return
        call kappa3_power(x,params,lz,value)
    end function generalized_normal_power

    pure subroutine generalized_normal_gamma(ly, params, lp, lq)
        implicit none
        real(dp), intent(in) :: ly !! log powered coordinate, possibly infinite
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp), intent(out) :: lp !! log probability between location and positive observation
        real(dp), intent(out) :: lq !! log positive half-law survival probability
        real(dp) :: y
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (ly == negative_infinity(0.0_dp)) return
        if (params(6) /= 0.0_dp) return
        if (ly > log(huge(1.0_dp))) then
            lp = 0.0_dp
            lq = negative_infinity(0.0_dp)
        else if (ly < log(tiny(1.0_dp))) then
            lp = gamma_small_kernel(params(4),ly)
            lq = log_difference(0.0_dp,lp)
        else
            y = exp(ly)
            if (params(4) > 1.0_dp) then
                if (params(4) > huge(1.0_dp)/max(1.0_dp,log(params(4)))) then
                    ! Here y is near one: the lower-gamma series correction is O(1/a).
                    lp = params(5)+params(4)*ly-y
                    lq = log_difference(0.0_dp,lp)
                    return
                end if
            end if
            lp = log_gammainc(params(4),y)
            lq = log_gammaincc(params(4),y)
        end if
    end subroutine generalized_normal_gamma

    pure subroutine generalized_normal_point(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above the half-law location
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, normalizer, flag
        real(dp), intent(out) :: lp !! log positive half-law CDF
        real(dp), intent(out) :: lq !! log positive half-law SF
        real(dp) :: lz, ly
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (x == params(2) .or. params(6) /= 0.0_dp) return
        call kappa3_power(x,params,lz,ly)
        if (ly < log(tiny(1.0_dp))) then
            ! (z**beta)**(1/beta) = z, even if beta*log(z) overflows.
            lp = lz+params(5)
            lq = log_difference(0.0_dp,lp)
        else
            call generalized_normal_gamma(ly,params,lp,lq)
        end if
    end subroutine generalized_normal_point

    pure function halfgennorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (params(6) /= 0.0_dp) return
        value = params(5)-log(params(3))+negative_hazard(generalized_normal_power(x,params))
    end function halfgennorm_lpdf

    pure function gennorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, reflected(6)
        reflected = params
        reflected(2) = -params(2)
        if (x < params(2)) then
            value = halfgennorm_lpdf(-x,reflected)-scifort_log_two
        else
            value = halfgennorm_lpdf(x,params)-scifort_log_two
        end if
    end function gennorm_lpdf

    pure subroutine generalized_normal_tails(x, params, half, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        logical, intent(in) :: half !! restrict support to the positive half
        real(dp), intent(out) :: lp !! log lower-tail probability
        real(dp), intent(out) :: lq !! log upper-tail probability
        real(dp) :: first, second, reflected(6)
        if (half) then
            lp = negative_infinity(0.0_dp)
            lq = 0.0_dp
            if (x <= params(2)) return
            call generalized_normal_point(x,params,lp,lq)
        else
            reflected = params
            reflected(2) = -params(2)
            if (x < params(2)) then
                call generalized_normal_point(-x,reflected,first,second)
                lp = second-scifort_log_two
                lq = log1p_safe(-exp(lp))
            else
                call generalized_normal_point(x,params,first,second)
                lq = second-scifort_log_two
                lp = log1p_safe(-exp(lq))
            end if
        end if
    end subroutine generalized_normal_tails

    pure function halfgennorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, lo, la, lb, lr, ld, lp, lq, rp, rq, t(3), r(3), correction
        integer :: j
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lo = max(a,params(2))
        lb = generalized_normal_power(b,params)
        if (lb < log(epsilon(1.0_dp))-2.0_dp) then
            value = params(5)-log(params(3))+log_distance(b,lo)
            if (params(6) /= 0.0_dp) value = negative_infinity(0.0_dp)
            return
        end if
        call generalized_normal_point(b,params,rp,rq)
        if (lo == params(2)) then
            value = rp
            return
        end if
        la = generalized_normal_power(lo,params)
        call generalized_normal_point(lo,params,lp,lq)
        if (la > log(huge(1.0_dp))) return
        lr = log_distance(b,lo)-log_distance(lo,params(2))
        ld = la+log_expm1_from_log(log(params(1))+burr_logsoftplus(lr))
        if (lr <= log(0.001_dp) .and. ld <= log(0.001_dp) &
                .and. log(params(1))+burr_logsoftplus(lr) <= log(0.001_dp)) then
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            do j = 1, 3
                ld = la+log_expm1_from_log(log(params(1))+burr_logsoftplus(log(t(j))+lr))
                r(j) = negative_hazard(ld)
            end do
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = halfgennorm_lpdf(lo,params)+log_distance(b,lo)+log1p_safe(correction)
        else if (rp <= lq) then
            value = log_difference(rp,lp)
        else
            value = log_difference(lq,rq)
        end if
        if (value > -scifort_log_two) value = log1p_safe(-exp(bounded_logadd(lp,rq)))
    end function halfgennorm_interval

    pure function gennorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, reflected(6), left, right, unused
        reflected = params
        reflected(2) = -params(2)
        if (b <= params(2)) then
            value = halfgennorm_interval(-b,-a,reflected)-scifort_log_two
        else if (a >= params(2)) then
            value = halfgennorm_interval(a,b,params)-scifort_log_two
        else
            call generalized_normal_point(-a,reflected,left,unused)
            call generalized_normal_point(b,params,right,unused)
            value = bounded_logadd(left,right)-scifort_log_two
        end if
        if (value > -scifort_log_two) then
            call generalized_normal_tails(a,params,.false.,left,unused)
            call generalized_normal_tails(b,params,.false.,unused,right)
            value = log1p_safe(-exp(bounded_logadd(left,right)))
        end if
    end function gennorm_interval

    pure function gennorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.false.,value,unused)
    end function gennorm_lcdf

    pure function gennorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.false.,unused,value)
    end function gennorm_lsf

    pure function halfgennorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.true.,value,unused)
    end function halfgennorm_lcdf

    pure function halfgennorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! power, location, scale, reciprocal power, log normalizer, overflow flag
        real(dp) :: value, unused
        call generalized_normal_tails(x,params,.true.,unused,value)
    end function halfgennorm_lsf


    pure function gennorm_censored_loglikelihood(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(6)
        params = [beta,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (beta <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        params(4) = location_scale_z(1.0_dp,[0.0_dp,beta])
        if (.not. ieee_is_finite(params(4))) return
        params(5) = generalized_normal_norm(params(4))
        if (.not. ieee_is_finite(params(5))) then
            params(5) = 0.0_dp
            params(6) = 1.0_dp
        end if
        value = censored_loglikelihood(data,params,gennorm_lpdf,gennorm_lcdf, &
            gennorm_lsf,gennorm_interval)
    end function gennorm_censored_loglikelihood

    pure function gennorm_censored_nnlf(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -gennorm_censored_loglikelihood(data, beta, loc, scale)
    end function gennorm_censored_nnlf

    pure function halfgennorm_censored_loglikelihood(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(6)
        params = [beta,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (beta <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        params(4) = location_scale_z(1.0_dp,[0.0_dp,beta])
        if (.not. ieee_is_finite(params(4))) return
        params(5) = generalized_normal_norm(params(4))
        if (.not. ieee_is_finite(params(5))) then
            params(5) = 0.0_dp
            params(6) = 1.0_dp
        end if
        value = censored_loglikelihood(data,params,halfgennorm_lpdf,halfgennorm_lcdf, &
            halfgennorm_lsf,halfgennorm_interval)
    end function halfgennorm_censored_loglikelihood

    pure function halfgennorm_censored_nnlf(data, beta, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: beta !! positive finite generalized-normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -halfgennorm_censored_loglikelihood(data, beta, loc, scale)
    end function halfgennorm_censored_nnlf

    pure function rdist_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite symmetric-beta shape, at most 1e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4), center, width
        params = [c,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (c > 1.0e6_dp) return
        if (params(3) <= 0.0_dp) return
        if (0.5_dp*c <= 0.0_dp) return
        if (params(3) > 0.5_dp*huge(1.0_dp)) return
        center = params(2)
        width = 2.0_dp*params(3)
        if (center < 0.0_dp) then
            if (params(3) > huge(1.0_dp)+center) return
        end if
        params = [0.5_dp*c,0.5_dp*c,center-params(3),width]
        if (.not. beta_parameters_valid(params)) return
        if (params(3) >= center .or. params(3)+width <= center) return
        value = censored_loglikelihood(data,params,rdist_lpdf,rdist_lcdf, &
            rdist_lsf,rdist_interval)
    end function rdist_censored_loglikelihood

    pure function rdist_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite symmetric-beta shape, at most 1e6
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -rdist_censored_loglikelihood(data, c, loc, scale)
    end function rdist_censored_nnlf


    pure subroutine folded_coordinates(x, params, minus, plus, lm, lp, lz)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp), intent(out) :: minus !! standardized distance minus offset
        real(dp), intent(out) :: plus !! standardized distance plus offset
        real(dp), intent(out) :: lm !! log absolute minus coordinate
        real(dp), intent(out) :: lp !! log plus coordinate
        real(dp), intent(out) :: lz !! log standardized distance
        real(dp) :: z, half, lc
        lz = negative_infinity(0.0_dp)
        if (x > params(2)) lz = log_distance(x,params(2))-log(params(3))
        lc = bounded_log(params(1))
        lp = bounded_logadd(lz,lc)
        z = location_scale_z(x,params(2:3))
        plus = finite_sum(z,params(1))
        minus = z-params(1)
        if (.not. ieee_is_finite(z)) then
            half = location_scale_z(0.5_dp*x,[0.5_dp*params(2),params(3)])
            minus = half-0.5_dp*params(1)
            minus = finite_sum(minus,minus)
        end if
        if (ieee_is_finite(minus)) then
            lm = bounded_log(abs(minus))
        else
            lm = lz+log1p_safe(-exp(lc-lz))
        end if
    end subroutine folded_coordinates

    pure function folded_normal_mass(a, b, lw) result(value)
        implicit none
        real(dp), intent(in) :: a !! standardized lower endpoint, possibly infinite
        real(dp), intent(in) :: b !! standardized upper endpoint, possibly coincident after rounding
        real(dp), intent(in) :: lw !! log exact standardized positive interval width
        real(dp) :: value, omitted
        value = normal_small_interval(a,lw)
        if (.not. ieee_is_nan(value)) return
        if (a >= 0.0_dp) then
            value = log_difference(normal_logsf(a),normal_logsf(b))
        else
            value = log_difference(normal_logcdf(b),normal_logcdf(a))
            if (b > 0.0_dp) then
                omitted = exp(normal_logcdf(a))+exp(normal_logsf(b))
                if (omitted <= 0.5_dp) value = log1p_safe(-omitted)
            end if
        end if
    end function folded_normal_mass

    pure function folded_cauchy_mass(a, b, la, lb, lw) result(value)
        implicit none
        real(dp), intent(in) :: a !! standardized lower endpoint
        real(dp), intent(in) :: b !! standardized upper endpoint
        real(dp), intent(in) :: la !! log absolute exact lower coordinate
        real(dp), intent(in) :: lb !! log absolute exact upper coordinate
        real(dp), intent(in) :: lw !! log exact positive standardized interval width
        real(dp) :: value, omitted
        if (a >= 0.0_dp .or. b <= 0.0_dp) then
            value = cauchy_logatan_exp(lw-log1pexp(la+lb))-log(scifort_pi)
        else if (min(la,lb) >= 0.0_dp) then
            omitted = exp(cauchy_log_small_tail(la))+exp(cauchy_log_small_tail(lb))
            value = log1p_safe(-omitted)
        else
            value = bounded_logadd(cauchy_logatan_exp(la),cauchy_logatan_exp(lb)) &
                -log(scifort_pi)
        end if
    end function folded_cauchy_mass

    pure function foldnorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, a, b, la, lb, lz
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        call folded_coordinates(x,params,a,b,la,lb,lz)
        value = bounded_logadd(finite_normal_lpdf(a),finite_normal_lpdf(b))-log(params(3))
    end function foldnorm_lpdf

    pure function foldcauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, a, b, la, lb, lz
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        call folded_coordinates(x,params,a,b,la,lb,lz)
        value = bounded_logadd(-log1pexp(2.0_dp*la),-log1pexp(2.0_dp*lb)) &
            -log(scifort_pi)-log(params(3))
    end function foldcauchy_lpdf

    pure subroutine folded_tails(x, params, normal, lcdf, lsf)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        logical, intent(in) :: normal !! normal rather than Cauchy base law
        real(dp), intent(out) :: lcdf !! log lower-tail probability
        real(dp), intent(out) :: lsf !! log upper-tail probability
        real(dp) :: a, b, la, lb, lz, first
        lcdf = negative_infinity(0.0_dp)
        lsf = 0.0_dp
        if (x <= params(2)) return
        call folded_coordinates(x,params,a,b,la,lb,lz)
        if (normal) then
            lcdf = folded_normal_mass(-b,a,scifort_log_two+lz)
            lsf = bounded_logadd(normal_logsf(a),normal_logsf(b))
        else
            lcdf = folded_cauchy_mass(-b,a,lb,la,scifort_log_two+lz)
            first = cauchy_log_small_tail(la)
            if (a < 0.0_dp) first = log1p_safe(-exp(first))
            lsf = bounded_logadd(first,cauchy_log_small_tail(lb))
        end if
        if (lcdf <= lsf) then
            lcdf = min(lcdf,0.0_dp)
            lsf = log1p_safe(-exp(lcdf))
        else
            lsf = min(lsf,0.0_dp)
            lcdf = log1p_safe(-exp(lsf))
        end if
    end subroutine folded_tails

    pure function folded_interval(a, b, params, normal) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint greater than a
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        logical, intent(in) :: normal !! normal rather than Cauchy base law
        real(dp) :: value, lo, am, ap, bm, bp, lam, lap, lbm, lbp, lz, lw, omitted, unused
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lo = max(a,params(2))
        call folded_coordinates(lo,params,am,ap,lam,lap,lz)
        call folded_coordinates(b,params,bm,bp,lbm,lbp,lz)
        lw = log_distance(b,lo)-log(params(3))
        if (normal) then
            value = bounded_logadd(folded_normal_mass(am,bm,lw),folded_normal_mass(ap,bp,lw))
        else
            value = bounded_logadd(folded_cauchy_mass(am,bm,lam,lbm,lw), &
                folded_cauchy_mass(ap,bp,lap,lbp,lw))
        end if
        if (value > -scifort_log_two) then
            call folded_tails(lo,params,normal,omitted,unused)
            call folded_tails(b,params,normal,unused,lz)
            value = log1p_safe(-(exp(omitted)+exp(lz)))
        end if
    end function folded_interval

    pure function genpareto_logh(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! interior point strictly above location
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, lu, lw, nh
        if (params(1) < 0.0_dp) then
            call genhalflogistic_coordinates(x,[-params(1),params(2),params(4)],lu,lw,nh,value)
        else
            value = log_distance(x,params(2))-log(params(3))
            if (params(1) > 0.0_dp) value = &
                burr_logsoftplus(log(params(1))+value)-log(params(1))
        end if
    end function genpareto_logh

    pure function genpareto_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, left, right, lu, lw, nh, lh
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (params(1) < 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            if (right < 0.0_dp) return
            if (right == 0.0_dp) then
                if (params(1) == -1.0_dp) value = -log(params(4))
                if (params(1) < -1.0_dp) value = positive_infinity(0.0_dp)
                return
            end if
            value = -log(-params(1))-log(params(4))
            if (left == 0.0_dp) return
            call genhalflogistic_coordinates(x,[-params(1),params(2),params(4)],lu,lw,nh,lh)
            if (params(1) <= -0.5_dp) then
                value = value+((1.0_dp+params(1))/(-params(1)))*lw
            else
                value = value+nh-lw
            end if
        else
            value = -log(params(3))
            if (x == params(2)) return
            value = value+negative_hazard(genpareto_logh(x,params))
            if (params(1) > 0.0_dp) value = value &
                -log1pexp(log(params(1))+log_distance(x,params(2))-log(params(3)))
        end if
    end function genpareto_lpdf

    pure function genpareto_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, left, right
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        if (params(1) < 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            value = 0.0_dp
            if (right <= 0.0_dp) return
        end if
        value = log_hazard_mass(genpareto_logh(x,params))
    end function genpareto_lcdf

    pure function genpareto_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, left, right
        value = 0.0_dp
        if (x <= params(2)) return
        if (params(1) < 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            value = negative_infinity(0.0_dp)
            if (right <= 0.0_dp) return
        end if
        value = negative_hazard(genpareto_logh(x,params))
    end function genpareto_lsf

    pure function genpareto_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! shape, location, scale, bounded width if negative
        real(dp) :: value, lo, left, right, delta, denom
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lo = max(a,params(2))
        if (params(1) < 0.0_dp) then
            call beta_distances(b,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
            if (right <= 0.0_dp) then
                value = genpareto_lsf(lo,params)
                return
            end if
            delta = burr_logsoftplus(log_distance(b,lo)-log(right))-log(-params(1))
        else if (params(1) == 0.0_dp) then
            delta = log_distance(b,lo)-log(params(3))
        else
            denom = log(params(3))
            if (lo > params(2)) denom = bounded_logadd(denom, &
                log(params(1))+log_distance(lo,params(2)))
            delta = burr_logsoftplus(log(params(1))+log_distance(b,lo)-denom)-log(params(1))
        end if
        value = genpareto_lsf(lo,params)+log_hazard_mass(delta)
    end function genpareto_interval

    pure function foldnorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.true.,value,unused)
    end function foldnorm_lcdf

    pure function foldnorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.true.,unused,value)
    end function foldnorm_lsf

    pure function foldnorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value
        value = folded_interval(a,b,params,.true.)
    end function foldnorm_interval

    pure function foldcauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.false.,value,unused)
    end function foldcauchy_lcdf

    pure function foldcauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value, unused
        call folded_tails(x,params,.false.,unused,value)
    end function foldcauchy_lsf

    pure function foldcauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! nonnegative offset, location, positive scale
        real(dp) :: value
        value = folded_interval(a,b,params,.false.)
    end function foldcauchy_interval


    pure function foldnorm_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c < 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,foldnorm_lpdf,foldnorm_lcdf, &
            foldnorm_lsf,foldnorm_interval)
    end function foldnorm_censored_loglikelihood

    pure function foldnorm_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -foldnorm_censored_loglikelihood(data, c, loc, scale)
    end function foldnorm_censored_nnlf

    pure function foldcauchy_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c < 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,foldcauchy_lpdf,foldcauchy_lcdf, &
            foldcauchy_lsf,foldcauchy_interval)
    end function foldcauchy_censored_loglikelihood

    pure function foldcauchy_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! nonnegative finite folding offset
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -foldcauchy_censored_loglikelihood(data, c, loc, scale)
    end function foldcauchy_censored_nnlf

    pure function genpareto_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite generalized Pareto shape, including zero
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [c,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(3) <= 0.0_dp) return
        if (c < 0.0_dp) then
            if (-c < 1.0_dp) then
                if (params(3) > huge(1.0_dp)*(-c)) return
            end if
            params(4) = params(3)/(-c)
            if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(2),params(4)])) return
            if (params(2)+params(4) <= params(2)) return
        end if
        value = censored_loglikelihood(data,params,genpareto_lpdf,genpareto_lcdf, &
            genpareto_lsf,genpareto_interval)
    end function genpareto_censored_loglikelihood

    pure function genpareto_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite generalized Pareto shape, including zero
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genpareto_censored_loglikelihood(data, c, loc, scale)
    end function genpareto_censored_nnlf


    pure function bounded_cutoff_product(a, b) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive finite cutoff multiplier
        real(dp), intent(in) :: b !! positive finite scale
        real(dp) :: value, half
        value = quiet_nan(0.0_dp)
        if (a > 1.0_dp) then
            if (b > huge(1.0_dp)/a) return
            if (b <= (0.5_dp*huge(1.0_dp))/a) then
                value = a*b
                return
            end if
            ! The quotient guard can round upward. Test the half product before doubling.
            half = (0.5_dp*a)*b
            if (half > 0.5_dp*huge(1.0_dp)) return
            value = 2.0_dp*half
        else
            value = a*b
        end if
    end function bounded_cutoff_product

    pure function bounded_power_valid(params) result(valid)
        implicit none
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        logical :: valid
        valid = .false.
        if (params(3) <= 0.0_dp .or. params(4) <= params(3)) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(2),params(4)])) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)+params(4) <= params(2)+params(3)) return
        valid = .true.
    end function bounded_power_valid

    pure subroutine bounded_power_gaps(x, params, left, right)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp), intent(out) :: left !! signed gap from lower cutoff
        real(dp), intent(out) :: right !! signed gap to upper cutoff
        real(dp) :: unused
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],unused,right)
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,left)
        left = -left
    end subroutine bounded_power_gaps

    pure function bounded_power_range(params) result(value)
        implicit none
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value
        value = burr_logsoftplus(log(params(4)-params(3))-log(params(3)))
    end function bounded_power_range

    pure function bounded_power_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, left, right, lr, lv, anchor, lk, lb, beta, v2, term, total
        integer :: j
        value = negative_infinity(0.0_dp)
        call bounded_power_gaps(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        lr = bounded_power_range(params)
        value = -log_distance(x,params(2))-lr
        if (params(1) == 0.0_dp) return
        lk = log(abs(params(1)))
        if (params(1) > 0.0_dp) then
            lv = burr_logsoftplus(bounded_log(left)-log(params(3)))
        else
            lv = burr_logsoftplus(bounded_log(right)-log_distance(x,params(2)))
        end if
        anchor = negative_hazard(lk+lv)
        lb = lk+lr
        if (lb <= 0.0_dp) then
            beta = exp(lb)
            v2 = (0.5_dp*beta)**2
            term = 1.0_dp
            total = 0.0_dp
            do j = 1, 10
                term = term*v2/real((2*j)*(2*j+1),dp)
                total = total+term
            end do
            value = value+0.5_dp*beta-log1p_safe(total)+anchor
        else
            value = lk-log_distance(x,params(2))+anchor-log_hazard_mass(lb)
        end if
    end function bounded_power_lpdf

    pure subroutine bounded_power_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: left, right, lv, lw, lr, lk, norm
        call bounded_power_gaps(x,params,left,right)
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (left <= 0.0_dp) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        lr = bounded_power_range(params)
        lv = burr_logsoftplus(log(left)-log(params(3)))
        lw = burr_logsoftplus(log(right)-log_distance(x,params(2)))
        if (params(1) == 0.0_dp) then
            lp = lv-lr
            lq = lw-lr
        else
            lk = log(abs(params(1)))
            norm = log_hazard_mass(lk+lr)
            lp = log_hazard_mass(lk+lv)-norm
            lq = log_hazard_mass(lk+lw)-norm
            if (params(1) > 0.0_dp) then
                lq = negative_hazard(lk+lv)+lq
            else
                lp = negative_hazard(lk+lw)+lp
            end if
        end if
        if (lp <= lq) then
            lp = min(lp,0.0_dp)
            lq = log1p_safe(-exp(lp))
        else
            lq = min(lq,0.0_dp)
            lp = log1p_safe(-exp(lq))
        end if
    end subroutine bounded_power_tails

    pure function bounded_power_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, unused
        call bounded_power_tails(x,params,value,unused)
    end function bounded_power_lcdf

    pure function bounded_power_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, unused
        call bounded_power_tails(x,params,unused,value)
    end function bounded_power_lsf

    pure function bounded_power_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! power, location, positive physical cutoffs
        real(dp) :: value, left, right, omitted, span, lk, anchor
        value = negative_infinity(0.0_dp)
        call bounded_power_gaps(b,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = bounded_power_lsf(a,params)
            return
        end if
        call bounded_power_gaps(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = bounded_power_lcdf(b,params)
            return
        end if
        omitted = exp(bounded_power_lcdf(a,params))+exp(bounded_power_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        span = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(2)))
        if (params(1) == 0.0_dp) then
            value = span-bounded_power_range(params)
        else
            lk = log(abs(params(1)))
            if (params(1) > 0.0_dp) then
                anchor = burr_logsoftplus(log(left)-log(params(3)))
            else
                call bounded_power_gaps(b,params,left,right)
                anchor = burr_logsoftplus(log(right)-log_distance(b,params(2)))
            end if
            value = negative_hazard(lk+anchor)+log_hazard_mass(lk+span) &
                -log_hazard_mass(lk+bounded_power_range(params))
        end if
    end function bounded_power_interval

    pure subroutine genhalflogistic_coordinates(x, params, lu, lw, nh, lh)
        implicit none
        real(dp), intent(in) :: x !! interior observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp), intent(out) :: lu !! log lower coordinate
        real(dp), intent(out) :: lw !! log upper coordinate
        real(dp), intent(out) :: nh !! negative transformed hazard
        real(dp), intent(out) :: lh !! log transformed hazard
        real(dp) :: left, right, u, w, bp(4)
        bp = [1.0_dp,1.0_dp,params(2),params(3)]
        call beta_distances(x,bp,left,right)
        call beta_coordinates(x,bp,left,right,u,w,lu,lw)
        if (lu < log(epsilon(1.0_dp))-2.0_dp) then
            lh = lu-log(params(1))
        else
            lh = log(-lw)-log(params(1))
        end if
        nh = negative_hazard(lh)
    end subroutine genhalflogistic_coordinates

    pure function genhalflogistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, left, right, lu, lw, nh, lh
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        value = negative_infinity(0.0_dp)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp) then
            value = -scifort_log_two-log(params(1))-log(params(3))
        else if (right == 0.0_dp) then
            if (params(1) == 1.0_dp) value = scifort_log_two-log(params(3))
            if (params(1) > 1.0_dp) value = positive_infinity(0.0_dp)
        else
            call genhalflogistic_coordinates(x,params,lu,lw,nh,lh)
            value = scifort_log_two-log(params(1))-log(params(3))-lw+nh &
                -2.0_dp*log1pexp(nh)
        end if
    end function genhalflogistic_lpdf

    pure subroutine genhalflogistic_tails(x, params, lp, lq)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp), intent(out) :: lp !! log lower probability
        real(dp), intent(out) :: lq !! log upper probability
        real(dp) :: left, right, lu, lw, nh, lh
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        lp = negative_infinity(0.0_dp)
        lq = 0.0_dp
        if (left <= 0.0_dp) return
        lp = 0.0_dp
        lq = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        call genhalflogistic_coordinates(x,params,lu,lw,nh,lh)
        lp = log_hazard_mass(lh)-log1pexp(nh)
        lq = scifort_log_two+nh-log1pexp(nh)
        if (lp <= lq) then
            lq = log1p_safe(-exp(lp))
        else
            lp = log1p_safe(-exp(lq))
        end if
    end subroutine genhalflogistic_tails

    pure function genhalflogistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, unused
        call genhalflogistic_tails(x,params,value,unused)
    end function genhalflogistic_lcdf

    pure function genhalflogistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, unused
        call genhalflogistic_tails(x,params,unused,value)
    end function genhalflogistic_lsf

    pure function genhalflogistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, rounded physical width
        real(dp) :: value, left, right, omitted, delta, lu, lw, nha, nhb, lh
        value = negative_infinity(0.0_dp)
        call beta_distances(b,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = genhalflogistic_lsf(a,params)
            return
        end if
        delta = burr_logsoftplus(log_distance(b,a)-log(right))-log(params(1))
        call beta_distances(a,[1.0_dp,1.0_dp,params(2),params(3)],left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = genhalflogistic_lcdf(b,params)
            return
        end if
        omitted = exp(genhalflogistic_lcdf(a,params))+exp(genhalflogistic_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call genhalflogistic_coordinates(a,params,lu,lw,nha,lh)
        call genhalflogistic_coordinates(b,params,lu,lw,nhb,lh)
        value = scifort_log_two+nha+log_hazard_mass(delta) &
            -log1pexp(nha)-log1pexp(nhb)
    end function genhalflogistic_interval


    pure function loguniform_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite lower shape, less than b
        real(dp), intent(in) :: b !! finite upper shape, greater than a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,b) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        if (a >= b) return
        if (b > 1.0_dp) then
            if (params(4) > huge(1.0_dp)/b) return
        end if
        params = [0.0_dp,params(3),bounded_cutoff_product(a,params(4)), &
            bounded_cutoff_product(b,params(4))]
        if (.not. all(ieee_is_finite(params))) return
        if (.not. bounded_power_valid(params)) return
        value = censored_loglikelihood(data,params,bounded_power_lpdf,bounded_power_lcdf, &
            bounded_power_lsf,bounded_power_interval)
    end function loguniform_censored_loglikelihood

    pure function loguniform_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite lower shape, less than b
        real(dp), intent(in) :: b !! finite upper shape, greater than a
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -loguniform_censored_loglikelihood(data, a, b, loc, scale)
    end function loguniform_censored_nnlf

    pure function truncpareto_censored_loglikelihood(data, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite nonzero power, either sign
        real(dp), intent(in) :: c !! finite upper cutoff, greater than one
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [b,c,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b == 0.0_dp .or. c <= 1.0_dp) return
        if (params(4) <= 0.0_dp) return
        if (c > 1.0_dp) then
            if (params(4) > huge(1.0_dp)/c) return
        end if
        params = [b,params(3),params(4),bounded_cutoff_product(c,params(4))]
        if (.not. all(ieee_is_finite(params))) return
        if (.not. bounded_power_valid(params)) return
        value = censored_loglikelihood(data,params,bounded_power_lpdf,bounded_power_lcdf, &
            bounded_power_lsf,bounded_power_interval)
    end function truncpareto_censored_loglikelihood

    pure function truncpareto_censored_nnlf(data, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite nonzero power, either sign
        real(dp), intent(in) :: c !! finite upper cutoff, greater than one
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -truncpareto_censored_loglikelihood(data, b, c, loc, scale)
    end function truncpareto_censored_nnlf

    pure function genhalflogistic_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        if (c < 1.0_dp) then
            if (params(3) > huge(1.0_dp)*c) return
        end if
        params(3) = params(3)/c
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(2),params(3)])) return
        if (params(2)+params(3) <= params(2)) return
        value = censored_loglikelihood(data,params,genhalflogistic_lpdf,genhalflogistic_lcdf, &
            genhalflogistic_lsf,genhalflogistic_interval)
    end function genhalflogistic_censored_loglikelihood

    pure function genhalflogistic_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genhalflogistic_censored_loglikelihood(data, c, loc, scale)
    end function genhalflogistic_censored_nnlf


    ! Independently integrated public density specifications; see CODE_PROVENANCE.md.
    pure subroutine truncexpon_gaps(x, params, left, right)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp), intent(out) :: left !! signed distance above lower support edge
        real(dp), intent(out) :: right !! compensated distance below upper support edge
        call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(4)],left,right)
    end subroutine truncexpon_gaps

    pure function truncexpon_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, left, right, cutoff, logcutoff, term, series, v2, gap, unused
        integer :: k
        call truncexpon_gaps(x,params,left,right)
        value = negative_infinity(0.0_dp)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        logcutoff = log(params(4))-log(params(3))
        if (logcutoff <= 0.0_dp) then
            ! exp(-z)/(1-exp(-cutoff)) = exp(cutoff/2-z)
            ! / (cutoff*sinh(cutoff/2)/(cutoff/2)). Retain the nearly uniform log density.
            cutoff = exp(logcutoff)
            v2 = (0.5_dp*cutoff)**2
            term = 1.0_dp
            series = 0.0_dp
            do k = 1, 10
                term = term*v2/real((2*k)*(2*k+1),dp)
                series = series+term
            end do
            if (params(4) >= 2.0_dp*tiny(1.0_dp)) then
                call beta_distances(x,[1.0_dp,1.0_dp,params(2),0.5_dp*params(4)],unused,gap)
                value = gap/params(3)
            else
                value = 0.5_dp*cutoff-location_scale_z(x,params(2:3))
            end if
            value = value-log(params(4))-log1p_safe(series)
            return
        end if
        value = -location_scale_z(x,params(2:3))-log(params(3)) &
            -log_hazard_mass(logcutoff)
    end function truncexpon_lpdf

    pure subroutine truncexpon_tails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: left, right, normalization
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        call truncexpon_gaps(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logp = 0.0_dp
            logq = negative_infinity(0.0_dp)
            return
        end if
        normalization = log_hazard_mass(log(params(4))-log(params(3)))
        logp = log_hazard_mass(log_distance(x,params(2))-log(params(3)))-normalization
        logq = -location_scale_z(x,params(2:3)) &
            +log_hazard_mass(log(right)-log(params(3)))-normalization
        if (logp <= logq) then
            logp = min(0.0_dp,logp)
            logq = log1p_safe(-exp(logp))
        else
            logq = min(0.0_dp,logq)
            logp = log1p_safe(-exp(logq))
        end if
    end subroutine truncexpon_tails

    pure function truncexpon_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, unused
        call truncexpon_tails(x,params,value,unused)
    end function truncexpon_lcdf

    pure function truncexpon_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, unused
        call truncexpon_tails(x,params,unused,value)
    end function truncexpon_lsf

    pure function truncexpon_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, rate scale, rounded width
        real(dp) :: value, left, right, omitted
        call truncexpon_gaps(b,params,left,right)
        value = negative_infinity(0.0_dp)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = truncexpon_lsf(a,params)
            return
        end if
        call truncexpon_gaps(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = truncexpon_lcdf(b,params)
            return
        end if
        omitted = exp(truncexpon_lcdf(a,params))+exp(truncexpon_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
        else
            value = -location_scale_z(a,params(2:3)) &
                +log_hazard_mass(log_distance(b,a)-log(params(3))) &
                -log_hazard_mass(log(params(4))-log(params(3)))
        end if
    end function truncexpon_interval

    pure function johnsonsb_affine(a, b, t) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite first Johnson shape
        real(dp), intent(in) :: b !! positive finite second Johnson shape
        real(dp), intent(in) :: t !! finite logit coordinate
        real(dp) :: value, halfsum
        if (abs(t) > 1.0_dp) then
            if (b > huge(1.0_dp)/abs(t)) then
                value = sign(positive_infinity(0.0_dp),t)
                if (0.5_dp*b > huge(1.0_dp)/abs(t)) return
                halfsum = finite_sum(0.5_dp*a,(0.5_dp*b)*t)
                if (abs(halfsum) > 0.5_dp*huge(1.0_dp)) return
                value = 2.0_dp*halfsum
                return
            end if
        end if
        value = finite_sum(a,b*t)
    end function johnsonsb_affine

    pure subroutine johnsonsb_coordinates(x, params, q, lz, lw)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly inside support
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp), intent(out) :: q !! transformed normal coordinate, possibly infinite
        real(dp), intent(out) :: lz !! log standardized distance above lower edge
        real(dp), intent(out) :: lw !! log standardized distance below upper edge
        real(dp) :: left, right, z, w, t, gap, unused, v, product
        call beta_distances(x,params,left,right)
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        t = lz-lw
        if (abs(t) < 0.1_dp .and. params(4) < 2.0_dp*tiny(1.0_dp)) then
            ! Subnormal widths can have an unrepresentable half-step midpoint.
            v = (left-right)/params(4)
            t = log1p_safe(v)-log1p_safe(-v)
            q = johnsonsb_affine(params(1),params(2),t)
            return
        end if
        if (abs(t) < 0.1_dp .and. 0.5_dp*params(4) > 0.0_dp) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(3),0.5_dp*params(4)],unused,gap)
            if (gap == 0.0_dp) then
                q = params(1)
                return
            end if
            v = (-2.0_dp*gap)/params(4)
            if (abs(v) < sqrt(epsilon(1.0_dp))) then
                product = signed_exp(log(params(2))+log(abs(gap))-log(params(4)) &
                    +2.0_dp*scifort_log_two,-sign(1.0_dp,gap))
                q = finite_sum(params(1),product)
                return
            end if
            t = log1p_safe(v)-log1p_safe(-v)
        end if
        q = johnsonsb_affine(params(1),params(2),t)
    end subroutine johnsonsb_coordinates

    pure function johnsonsb_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, left, right, q, lz, lw
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (min(left,right) <= 0.0_dp) return
        call johnsonsb_coordinates(x,params,q,lz,lw)
        value = finite_normal_lpdf(q)+log(params(2))-log(params(4))-lz-lw
    end function johnsonsb_lpdf

    pure subroutine johnsonsb_tails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: left, right, q, lz, lw
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logp = 0.0_dp
            logq = negative_infinity(0.0_dp)
            return
        end if
        call johnsonsb_coordinates(x,params,q,lz,lw)
        logp = normal_logcdf(q)
        logq = normal_logsf(q)
    end subroutine johnsonsb_tails

    pure function johnsonsb_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, unused
        call johnsonsb_tails(x,params,value,unused)
    end function johnsonsb_lcdf

    pure function johnsonsb_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, unused
        call johnsonsb_tails(x,params,unused,value)
    end function johnsonsb_lsf

    pure function johnsonsb_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! finite a, positive b, lower edge, positive width
        real(dp) :: value, left, right, q, lz, lw, span, first, second, omitted
        call beta_distances(b,params,left,right)
        value = negative_infinity(0.0_dp)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            value = johnsonsb_lsf(a,params)
            return
        end if
        second = log(right)
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = johnsonsb_lcdf(b,params)
            return
        end if
        span = log_distance(b,a)
        first = burr_logsoftplus(span-log_distance(a,params(3)))
        second = burr_logsoftplus(span-second)
        call johnsonsb_coordinates(a,params,q,lz,lw)
        value = normal_small_interval(q,log(params(2))+bounded_logadd(first,second))
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,johnsonsb_lcdf,johnsonsb_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(johnsonsb_lcdf(a,params),johnsonsb_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function johnsonsb_interval

    pure subroutine kappa3_power(x, params, lz, power)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above lower support edge
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp), intent(out) :: lz !! log standardized distance, possibly rounded to zero
        real(dp), intent(out) :: power !! shape times log standardized distance
        real(dp) :: z, gap, unused, logr
        z = location_scale_z(x,params(2:3))
        if (z >= 0.5_dp .and. z <= 2.0_dp) then
            if (params(2) > 0.0_dp) then
                if (params(3) > huge(1.0_dp)-params(2)) then
                    ! Reorder scale-(x-loc) when loc+scale is unrepresentable.
                    call beta_distances(-params(2),[1.0_dp,1.0_dp,-x,params(3)],unused,gap)
                else
                    call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,gap)
                end if
            else
                call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,gap)
            end if
            lz = 0.0_dp
            power = 0.0_dp
            if (gap == 0.0_dp) return
            logr = log(abs(gap))-log(params(3))
            if (logr < log(epsilon(1.0_dp))-2.0_dp) then
                lz = -gap/params(3)
                power = signed_exp(log(params(1))+logr,-sign(1.0_dp,gap))
                return
            end if
            lz = log1p_safe(-gap/params(3))
        else
            lz = log_distance(x,params(2))-log(params(3))
        end if
        power = sign(positive_infinity(0.0_dp),lz)
        if (params(1) > 1.0_dp) then
            if (abs(lz) > huge(1.0_dp)/params(1)) return
        end if
        power = params(1)*lz
    end subroutine kappa3_power

    pure subroutine kappa3_hazard(x, params, q, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above lower support edge
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp), intent(out) :: q !! log(a)-a*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log reversed cumulative hazard
        real(dp), intent(out) :: nh !! negative reversed cumulative hazard
        real(dp) :: power, lz
        call kappa3_power(x,params,lz,power)
        q = -power
        if (ieee_is_finite(power)) q = finite_sum(log(params(1)),-power)
        if (q >= 40.0_dp) then
            nh = lz-log(params(1))/params(1)-log1pexp(-q)/params(1)
            logh = log(-nh)
        else
            logh = burr_logsoftplus(q)-log(params(1))
            nh = negative_hazard(logh)
        end if
    end subroutine kappa3_hazard

    pure function kappa3_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, q, logh, nh, la
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        la = log(params(1))
        if (x == params(2)) then
            value = -log(params(3))
            if (la < 0.0_dp) value = signed_exp(log(-la)-la,1.0_dp)-log(params(3))
            if (la > 0.0_dp) value = -la/params(1)-log(params(3))
            return
        end if
        call kappa3_hazard(x,params,q,logh,nh)
        if (q >= 40.0_dp) then
            value = -la/params(1)-log(params(3)) &
                -(1.0_dp+1.0_dp/params(1))*log1pexp(-q)
        else
            value = nh-log_distance(x,params(2))-log1pexp(-q)
        end if
    end function kappa3_lpdf

    pure function kappa3_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, q, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call kappa3_hazard(x,params,q,logh,value)
    end function kappa3_lcdf

    pure function kappa3_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, q, logh, nh
        value = 0.0_dp
        if (x <= params(2)) return
        call kappa3_hazard(x,params,q,logh,nh)
        value = hazard_complement(logh,nh)
    end function kappa3_lsf

    pure function kappa3_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, qa, qb, logh, nh, span, delta, diff, u, omitted, d
        value = kappa3_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(kappa3_lcdf(a,params))+exp(kappa3_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call kappa3_hazard(a,params,qa,logh,nh)
        call kappa3_hazard(b,params,qb,logh,value)
        span = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(2)))
        if (qb >= 40.0_dp) then
            diff = span
        else
            delta = log(params(1))+span
            if (delta < log(0.5_dp)) then
                u = log_expm1_from_log(delta)
                diff = burr_logsoftplus(-log1pexp(-qb)+u)-log(params(1))
            else if (qb >= 0.0_dp) then
                diff = delta
                if (delta < log(huge(1.0_dp))) then
                    d = exp(delta)
                    diff = log(d+log1pexp(-qa)-log1pexp(-qb))
                end if
                diff = diff-log(params(1))
            else if (qa <= 0.0_dp) then
                u = -log1pexp(-qa)+log_hazard_mass(delta)
                diff = u
                if (u >= log(epsilon(1.0_dp))-2.0_dp) &
                    diff = log(-log1p_safe(-exp(u)))
                diff = diff-log(params(1))
            else if (ieee_is_finite(qa)) then
                diff = log(log1pexp(qa)-log1pexp(qb))-log(params(1))
            else
                diff = log(-lognormal_logz(a,params))
            end if
        end if
        value = value+log_hazard_mass(diff)
    end function kappa3_interval


    pure function truncexpon_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite standardized cutoff
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [b,0.0_dp,1.0_dp,0.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        if (b > 1.0_dp) then
            if (params(3) > huge(1.0_dp)/b) return
        end if
        params(4) = b*params(3)
        if (params(4) <= 0.0_dp) return
        if (params(2) > 0.0_dp) then
            if (params(4) > huge(1.0_dp)-params(2)) return
        end if
        if (params(2)+params(4) <= params(2)) return
        value = censored_loglikelihood(data,params,truncexpon_lpdf,truncexpon_lcdf, &
            truncexpon_lsf,truncexpon_interval)
    end function truncexpon_censored_loglikelihood

    pure function truncexpon_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite standardized cutoff
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -truncexpon_censored_loglikelihood(data, b, loc, scale)
    end function truncexpon_censored_nnlf

    pure function johnsonsb_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite normal-logit shift
        real(dp), intent(in) :: b !! positive finite normal-logit scale
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite support width, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(3),params(4)])) return
        if (params(3)+params(4) <= params(3)) return
        value = censored_loglikelihood(data,params,johnsonsb_lpdf,johnsonsb_lcdf, &
            johnsonsb_lsf,johnsonsb_interval)
    end function johnsonsb_censored_loglikelihood

    pure function johnsonsb_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite normal-logit shift
        real(dp), intent(in) :: b !! positive finite normal-logit scale
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite support width, default one
        real(dp) :: value
        value = -johnsonsb_censored_loglikelihood(data, a, b, loc, scale)
    end function johnsonsb_censored_nnlf

    pure function kappa3_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite kappa exponent
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (a <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,kappa3_lpdf,kappa3_lcdf, &
            kappa3_lsf,kappa3_interval)
    end function kappa3_censored_loglikelihood

    pure function kappa3_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite kappa exponent
        real(dp), intent(in), optional :: loc !! finite lower support edge, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -kappa3_censored_loglikelihood(data, a, loc, scale)
    end function kappa3_censored_nnlf


    pure function cosine_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3), sigma
        params = [1.0_dp,0.0_dp,1.0_dp]
        sigma = 1.0_dp
        if (present(loc)) params(2) = loc
        if (present(scale)) sigma = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(2),sigma)) return
        if (sigma > huge(1.0_dp)/3.1415926535897931_dp) return
        params(3) = 3.1415926535897931_dp*sigma
        if (params(3) > huge(1.0_dp)-abs(params(2))) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)-params(3) >= params(2)) return
        value = censored_loglikelihood(data,params,compact_lpdf,compact_lcdf, &
            compact_lsf,compact_interval)
    end function cosine_censored_loglikelihood

    pure function cosine_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -cosine_censored_loglikelihood(data,loc,scale)
    end function cosine_censored_nnlf

    pure function anglit_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3), sigma
        params = [2.0_dp,0.0_dp,1.0_dp]
        sigma = 1.0_dp
        if (present(loc)) params(2) = loc
        if (present(scale)) sigma = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(2),sigma)) return
        params(3) = 0.78539816339744828_dp*sigma
        if (params(3) > huge(1.0_dp)-abs(params(2))) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)-params(3) >= params(2)) return
        value = censored_loglikelihood(data,params,compact_lpdf,compact_lcdf, &
            compact_lsf,compact_interval)
    end function anglit_censored_loglikelihood

    pure function anglit_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -anglit_censored_loglikelihood(data,loc,scale)
    end function anglit_censored_nnlf

    pure function semicircular_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3), sigma
        params = [3.0_dp,0.0_dp,1.0_dp]
        sigma = 1.0_dp
        if (present(loc)) params(2) = loc
        if (present(scale)) sigma = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(2),sigma)) return
        params(3) = 1_dp*sigma
        if (params(3) > huge(1.0_dp)-abs(params(2))) return
        if (params(2)+params(3) <= params(2)) return
        if (params(2)-params(3) >= params(2)) return
        value = censored_loglikelihood(data,params,compact_lpdf,compact_lcdf, &
            compact_lsf,compact_interval)
    end function semicircular_censored_loglikelihood

    pure function semicircular_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite support center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -semicircular_censored_loglikelihood(data,loc,scale)
    end function semicircular_censored_nnlf

    ! Independent integration of the SciPy public density specifications.
    ! The rounded half-width is the numerical support contract; see the guide.
    pure subroutine compact_gap(x, params, gap, logq)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp), intent(out) :: gap !! signed distance below the nearer support edge
        real(dp), intent(out) :: logq !! log distance to the nearer edge divided by half-width
        real(dp) :: unused
        if (x >= params(2)) then
            call beta_distances(x,[1.0_dp,1.0_dp,params(2),params(3)],unused,gap)
        else
            call beta_distances(-x,[1.0_dp,1.0_dp,-params(2),params(3)],unused,gap)
        end if
        logq = negative_infinity(0.0_dp)
        if (gap > 0.0_dp) logq = min(0.0_dp,log(gap)-log(params(3)))
    end subroutine compact_gap

    pure function compact_logsin(logx) result(value)
        implicit none
        real(dp), intent(in) :: logx !! log angle in (0, pi/2]
        real(dp) :: value
        value = logx
        ! The omitted relative correction is x**2/6 < epsilon/6.
        if (logx > 0.5_dp*log(epsilon(1.0_dp))) value = log(sin(exp(logx)))
    end function compact_logsin

    pure function compact_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, gap, lq
        call compact_gap(x,params,gap,lq)
        value = negative_infinity(0.0_dp)
        if (gap <= 0.0_dp) return
        value = compact_density(lq,params)
    end function compact_lpdf

    pure function compact_density(lq, params) result(value)
        implicit none
        real(dp), intent(in) :: lq !! log nearer endpoint gap divided by half-width
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value
        select case (int(params(1)))
        case (1)
            value = 2.0_dp*compact_logsin(lq+log(scifort_pi/2.0_dp))
        case (2)
            value = log(scifort_pi/4.0_dp)+compact_logsin(lq+log(scifort_pi/2.0_dp))
        case (3)
            value = log(2.0_dp/scifort_pi)+0.5_dp*(lq+log(2.0_dp-exp(lq)))
        end select
        value = value-log(params(3))
    end function compact_density

    pure function compact_small_tail(logq, family) result(value)
        implicit none
        real(dp), intent(in) :: logq !! log distance to the nearer endpoint in half-width units
        integer, intent(in) :: family !! 1 cosine, 2 anglit, 3 semicircular
        real(dp) :: value, q, term, total, v
        integer :: k
        q = exp(logq)
        select case (family)
        case (1)
            if (q <= 0.5_dp) then
                ! Integrate sin(pi*q/2)**2: q-sin(pi*q)/pi, without subtraction.
                v = (scifort_pi*q)**2
                term = 1.0_dp
                total = term
                do k = 2, 16
                    term = -term*v/real((2*k)*(2*k+1),dp)
                    total = total+term
                end do
                value = 3.0_dp*logq+2.0_dp*log(scifort_pi)-log(12.0_dp)+log(total)
            else
                value = log((q-sin(scifort_pi*q)/scifort_pi)/2.0_dp)
            end if
        case (2)
            value = 2.0_dp*compact_logsin(logq+log(scifort_pi/4.0_dp))
        case (3)
            if (q <= 0.25_dp) then
                ! Integrate sqrt(q*(2-q)) using the binomial series for sqrt(1-q/2).
                term = 1.0_dp
                total = term
                do k = 1, 24
                    term = term*(real(k,dp)-1.5_dp)*q/(2.0_dp*real(k,dp))
                    total = total+term*3.0_dp/real(2*k+3,dp)
                end do
                value = 1.5_dp*logq+log(4.0_dp*sqrt(2.0_dp)/(3.0_dp*scifort_pi)) &
                    +log(total)
            else
                value = log((acos(1.0_dp-q)-(1.0_dp-q)*sqrt(q*(2.0_dp-q)))/scifort_pi)
            end if
        end select
    end function compact_small_tail

    pure subroutine compact_tails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp), intent(out) :: logp !! log lower-tail probability
        real(dp), intent(out) :: logq !! log upper-tail probability
        real(dp) :: gap, lq, small, large
        call compact_gap(x,params,gap,lq)
        small = negative_infinity(0.0_dp)
        large = 0.0_dp
        if (gap > 0.0_dp) then
            small = min(-scifort_log_two,compact_small_tail(lq,int(params(1))))
            large = log1p_safe(-exp(small))
        end if
        logp = small
        logq = large
        if (x >= params(2)) then
            logp = large
            logq = small
        end if
    end subroutine compact_tails

    pure function compact_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, unused
        call compact_tails(x,params,value,unused)
    end function compact_lcdf

    pure function compact_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, unused
        call compact_tails(x,params,unused,value)
    end function compact_lsf

    pure function compact_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint greater than a
        real(dp), intent(in) :: params(:) !! family number, center, rounded half-width
        real(dp) :: value, ga, gb, la, lb, lw, anchor, correction, rho, lq, t(3), r(3)
        integer :: k
        call compact_gap(a,params,ga,la)
        call compact_gap(b,params,gb,lb)
        lw = log_distance(b,a)
        if (min(ga,gb) > 0.0_dp .and. lw < log(huge(1.0_dp))-1.0_dp) then
            if (lw-min(la,lb)-log(params(3)) <= log(0.001_dp)) then
                anchor = compact_lpdf(a,params)
                t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                    0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
                do k = 1, 3
                    ! Integrate in gap coordinates: physical quadrature nodes can
                    ! round to the original endpoints for adjacent observations.
                    if (a >= params(2)) then
                        rho = exp(lw-log(params(3))-la)
                        lq = la+log1p_safe(-t(k)*rho)
                    else if (b <= params(2)) then
                        rho = exp(lw-log(params(3))-la)
                        lq = min(0.0_dp,la+log1p_safe(t(k)*rho))
                    else
                        rho = exp(lw-log(params(3)))
                        lq = log1p_safe(-abs(expm1_safe(la)+t(k)*rho))
                    end if
                    r(k) = compact_density(lq,params)-anchor
                end do
                correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                    +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
                value = anchor+lw+log1p_safe(correction)
                return
            end if
        end if
        value = tail_interval(a,b,params,compact_lcdf,compact_lsf)
        if (value > -scifort_log_two) then
            value = log1p_safe(-exp(bounded_logadd(compact_lcdf(a,params),compact_lsf(b,params))))
        end if
    end function compact_interval


    pure function finite_normal_lpdf(z) result(value)
        implicit none
        real(dp), intent(in) :: z !! normal coordinate, possibly infinite
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        value = (-0.5_dp*z)*z-scifort_log_sqrt_two_pi
    end function finite_normal_lpdf

    pure function signed_exp(logabs, direction) result(value)
        implicit none
        real(dp), intent(in) :: logabs !! log nonnegative magnitude
        real(dp), intent(in) :: direction !! sign to apply, zero treated as positive
        real(dp) :: value
        value = sign(positive_infinity(0.0_dp),direction)
        if (logabs <= log(huge(1.0_dp))) value = sign(exp(logabs),direction)
    end function signed_exp

    pure subroutine fatiguelife_coordinates(x, params, w, v, logjac)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp), intent(out) :: w !! transformed normal coordinate
        real(dp), intent(out) :: v !! log standardized distance
        real(dp), intent(out) :: logjac !! log derivative of the normal coordinate
        real(dp) :: logw
        v = lognormal_logz(x,params)
        logjac = 0.5_dp*abs(v)+log1p_safe(exp(-abs(v)))-scifort_log_two &
            -log(params(1))-log_distance(x,params(2))
        w = 0.0_dp
        if (v == 0.0_dp) return
        if (abs(v) < 0.5_dp) then
            logw = log(abs(expm1_safe(v)))-0.5_dp*v-log(params(1))
        else
            logw = 0.5_dp*abs(v)+log1p_safe(-exp(-abs(v)))-log(params(1))
        end if
        w = signed_exp(logw,v)
    end subroutine fatiguelife_coordinates

    pure function fatiguelife_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, w, v, logjac
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call fatiguelife_coordinates(x,params,w,v,logjac)
        value = finite_normal_lpdf(w)+logjac
    end function fatiguelife_lpdf

    pure function fatiguelife_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, w, v, logjac
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call fatiguelife_coordinates(x,params,w,v,logjac)
        value = normal_logcdf(w)
    end function fatiguelife_lcdf

    pure function fatiguelife_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, w, v, logjac
        value = 0.0_dp
        if (x <= params(2)) return
        call fatiguelife_coordinates(x,params,w,v,logjac)
        value = normal_logsf(w)
    end function fatiguelife_lsf

    pure function fatiguelife_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, wa, wb, va, vb, logjac, logspan, logwidth, omitted
        value = fatiguelife_lcdf(b,params)
        if (a <= params(2)) return
        call fatiguelife_coordinates(a,params,wa,va,logjac)
        call fatiguelife_coordinates(b,params,wb,vb,logjac)
        logspan = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(2)))
        logwidth = log_expm1_from_log(logspan-scifort_log_two) &
            +bounded_logadd(0.5_dp*va,-0.5_dp*vb)-log(params(1))
        value = normal_small_interval(wa,logwidth)
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,fatiguelife_lcdf,fatiguelife_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(fatiguelife_lcdf(a,params),fatiguelife_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function fatiguelife_interval

    pure function finite_sum(a, b) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite first summand
        real(dp), intent(in) :: b !! finite second summand
        real(dp) :: value
        if (a > 0.0_dp .and. b > 0.0_dp) then
            if (b > huge(1.0_dp)-a) then
                value = positive_infinity(0.0_dp)
                return
            end if
        else if (a < 0.0_dp .and. b < 0.0_dp) then
            if (b < -huge(1.0_dp)-a) then
                value = negative_infinity(0.0_dp)
                return
            end if
        end if
        value = a+b
    end function finite_sum

    pure subroutine johnsonsu_coordinates(x, params, q, t, loghypot)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp), intent(out) :: q !! transformed normal coordinate
        real(dp), intent(out) :: t !! inverse hyperbolic sine of standardized point
        real(dp), intent(out) :: loghypot !! log sqrt(1+z*z)
        real(dp) :: z, logz, product, halfsum
        logz = cauchy_log_distance(x,params(3:4))
        z = location_scale_z(x,params(3:4))
        loghypot = 0.5_dp*log1pexp(2.0_dp*logz)
        if (logz < log(epsilon(1.0_dp))-2.0_dp) then
            t = z
            product = signed_exp(log(params(2))+logz,sign(1.0_dp,z))
            if (x < params(3)) product = -abs(product)
            q = finite_sum(params(1),product)
            return
        else if (logz > log(1.0e8_dp)) then
            t = sign(logz+scifort_log_two,z)
        else
            t = asinh(z)
        end if
        if (abs(t) > 1.0_dp) then
            if (params(2) > huge(1.0_dp)/abs(t)) then
                ! Halving before the product allows an opposing finite a
                ! to rescue a product between huge() and 2*huge().
                q = sign(positive_infinity(0.0_dp),t)
                if (0.5_dp*params(2) > huge(1.0_dp)/abs(t)) return
                halfsum = finite_sum(0.5_dp*params(1),(0.5_dp*params(2))*t)
                if (abs(halfsum) > 0.5_dp*huge(1.0_dp)) return
                q = 2.0_dp*halfsum
                return
            end if
        end if
        q = finite_sum(params(1),params(2)*t)
    end subroutine johnsonsu_coordinates

    pure function johnsonsu_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, q, t, loghypot
        call johnsonsu_coordinates(x,params,q,t,loghypot)
        value = finite_normal_lpdf(q)+log(params(2))-log(params(4))-loghypot
    end function johnsonsu_lpdf

    pure function johnsonsu_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, q, t, loghypot
        call johnsonsu_coordinates(x,params,q,t,loghypot)
        value = normal_logcdf(q)
    end function johnsonsu_lcdf

    pure function johnsonsu_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, q, t, loghypot
        call johnsonsu_coordinates(x,params,q,t,loghypot)
        value = normal_logsf(q)
    end function johnsonsu_lsf

    pure function johnsonsu_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite a, positive b, location, scale
        real(dp) :: value, qa, qb, ta, tb, la, lb, logr, logwidth, r, omitted
        call johnsonsu_coordinates(a,params,qa,ta,la)
        call johnsonsu_coordinates(b,params,qb,tb,lb)
        ! tanh((asinh(zb)-asinh(za))/2)=(zb-za)/(hypot(1,za)+hypot(1,zb)).
        logr = log_distance(b,a)-log(params(4))-bounded_logadd(la,lb)
        if (logr < log(epsilon(1.0_dp))-2.0_dp) then
            logwidth = scifort_log_two+logr
        else if (logr <= -scifort_log_two) then
            r = exp(logr)
            logwidth = log(log1p_safe(2.0_dp*r/(1.0_dp-r)))
        else
            logwidth = log(tb-ta)
        end if
        logwidth = log(params(2))+logwidth
        value = normal_small_interval(qa,logwidth)
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,johnsonsu_lcdf,johnsonsu_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(johnsonsu_lcdf(a,params),johnsonsu_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function johnsonsu_interval

    pure function genexpon_log_remainder(logv) result(value)
        implicit none
        real(dp), intent(in) :: logv !! log nonnegative v in v+expm1(-v)
        real(dp) :: value, v, term, series
        integer :: k
        if (logv < log(epsilon(1.0_dp))-2.0_dp) then
            value = 2.0_dp*logv-scifort_log_two
        else if (logv <= log(0.5_dp)) then
            v = exp(logv)
            term = 1.0_dp
            series = 1.0_dp
            do k = 3, 22
                term = -term*v/real(k,dp)
                series = series+term
            end do
            value = 2.0_dp*logv-scifort_log_two+log(series)
        else if (logv > log(40.0_dp)) then
            value = logv+log1p_safe(-exp(-logv))
        else
            v = exp(logv)
            value = log(v+expm1_safe(-v))
        end if
    end function genexpon_log_remainder

    pure function genexpon_logh(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(4)) return
        logz = lognormal_logz(x,[1.0_dp,params(4),params(5)])
        value = bounded_logadd(log(params(1))+logz,log(params(2))-log(params(3)) &
            +genexpon_log_remainder(log(params(3))+logz))
    end function genexpon_logh

    pure function genexpon_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value, logv
        value = negative_infinity(0.0_dp)
        if (x < params(4)) return
        value = log(params(1))-log(params(5))
        if (x == params(4)) return
        logv = log(params(3))+lognormal_logz(x,[1.0_dp,params(4),params(5)])
        value = bounded_logadd(log(params(1)),log(params(2))+log_hazard_mass(logv)) &
            -log(params(5))+negative_hazard(genexpon_logh(x,params))
    end function genexpon_lpdf

    pure function genexpon_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value
        value = log_hazard_mass(genexpon_logh(x,params))
    end function genexpon_lcdf

    pure function genexpon_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value
        value = negative_hazard(genexpon_logh(x,params))
    end function genexpon_lsf

    pure function genexpon_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive a, b, c, location, scale
        real(dp) :: value, lower, logdelta, logv, logdv, first, second, third, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(4)) return
        lower = max(a,params(4))
        omitted = exp(genexpon_lcdf(lower,params))+exp(genexpon_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        logdelta = log_distance(b,lower)-log(params(5))
        logv = negative_infinity(0.0_dp)
        if (lower > params(4)) then
            logv = log(params(3))+lognormal_logz(lower,[1.0_dp,params(4),params(5)])
        end if
        logdv = log(params(3))+logdelta
        first = log(params(1))+logdelta
        second = log(params(2))+logdelta+log_hazard_mass(logv)
        third = log(params(2))-log(params(3))+negative_hazard(logv) &
            +genexpon_log_remainder(logdv)
        value = negative_hazard(genexpon_logh(lower,params)) &
            +log_hazard_mass(bounded_logadd(first,bounded_logadd(second,third)))
    end function genexpon_interval


    pure function fatiguelife_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite fatigue-life shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,fatiguelife_lpdf,fatiguelife_lcdf, &
            fatiguelife_lsf,fatiguelife_interval)
    end function fatiguelife_censored_loglikelihood

    pure function fatiguelife_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite fatigue-life shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -fatiguelife_censored_loglikelihood(data, c, loc, scale)
    end function fatiguelife_censored_nnlf

    pure function johnsonsu_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first Johnson shape
        real(dp), intent(in) :: b !! positive finite Johnson scale shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,johnsonsu_lpdf,johnsonsu_lcdf, &
            johnsonsu_lsf,johnsonsu_interval)
    end function johnsonsu_censored_loglikelihood

    pure function johnsonsu_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first Johnson shape
        real(dp), intent(in) :: b !! positive finite Johnson scale shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -johnsonsu_censored_loglikelihood(data, a, b, loc, scale)
    end function johnsonsu_censored_nnlf

    pure function genexpon_censored_loglikelihood(data, a, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite baseline rate
        real(dp), intent(in) :: b !! positive finite additional rate
        real(dp), intent(in) :: c !! positive finite rate-transition parameter
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(5)
        params = [a,b,c,0.0_dp,1.0_dp]
        if (present(loc)) params(4) = loc
        if (present(scale)) params(5) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,b,c) <= 0.0_dp) return
        if (params(5) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,genexpon_lpdf,genexpon_lcdf, &
            genexpon_lsf,genexpon_interval)
    end function genexpon_censored_loglikelihood

    pure function genexpon_censored_nnlf(data, a, b, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite baseline rate
        real(dp), intent(in) :: b !! positive finite additional rate
        real(dp), intent(in) :: c !! positive finite rate-transition parameter
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -genexpon_censored_loglikelihood(data, a, b, c, loc, scale)
    end function genexpon_censored_nnlf


    pure function weibull_log_increment(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint strictly above location
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive Weibull shape, location, scale
        real(dp) :: value, logspan, logdelta
        ! exp(logt(b))-exp(logt(a)), preserving the original endpoint gap.
        logspan = log(params(1))+burr_logsoftplus( &
            log_distance(b,a)-log_distance(a,params(2)))
        if (logspan > log(50.0_dp)) then
            value = weibull_log_hazard(b,params)
        else
            logdelta = log_expm1_from_log(logspan)
            value = weibull_log_hazard(a,params)+logdelta
        end if
    end function weibull_log_increment

    pure function exponpow_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value, logt, t
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(3))
            return
        end if
        logt = weibull_log_hazard(x,params)
        if (logt > log(log(huge(1.0_dp)))) return
        t = exp(logt)
        value = log(params(1))-log_distance(x,params(2))+logt+t-expm1_safe(t)
    end function exponpow_lpdf

    pure function exponpow_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        value = log_hazard_mass(log_expm1_from_log(weibull_log_hazard(x,params)))
    end function exponpow_lcdf

    pure function exponpow_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value, logt
        value = 0.0_dp
        if (x <= params(2)) return
        logt = weibull_log_hazard(x,params)
        value = negative_hazard(log_expm1_from_log(logt))
        if (logt >= log(tiny(1.0_dp)) .and. logt <= log(log(huge(1.0_dp)))) &
            value = -expm1_safe(exp(logt))
    end function exponpow_lsf

    pure function exponpow_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive power shape, location, scale
        real(dp) :: value, logt, logdelta, omitted
        value = exponpow_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(exponpow_lcdf(a,params))+exp(exponpow_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        value = exponpow_lsf(a,params)
        if (.not. ieee_is_finite(value)) return
        logt = weibull_log_hazard(a,params)
        logdelta = weibull_log_increment(a,b,params)
        value = value+log_hazard_mass(exp(logt)+log_expm1_from_log(logdelta))
    end function exponpow_interval

    pure subroutine exponweib_hazard(x, params, logt, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp), intent(out) :: logt !! c*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log(-a*log(1-exp(-z**c)))
        real(dp), intent(out) :: nh !! negative reversed cumulative hazard
        real(dp) :: v, base, a
        a = params(1)
        logt = weibull_log_hazard(x,params(2:4))
        if (logt < log(epsilon(1.0_dp))-2.0_dp) then
            ! For tiny z**c, log(F_weibull) = c*log(z) to relative roundoff.
            ! Form a*c*abs(log(z)) without requiring c*log(z) to be finite.
            v = lognormal_logz(x,params(2:4))
            logh = log(a)+log(params(2))+log(-v)
            nh = negative_hazard(logh)
            if (ieee_is_finite(logt)) then
                if (a <= 1.0_dp) then
                    nh = a*logt
                else if (abs(logt) <= huge(1.0_dp)/a) then
                    nh = a*logt
                end if
            end if
        else if (logt > log(36.0_dp)) then
            ! -log(1-exp(-t)) / exp(-t) differs from one by < 2e-16.
            logh = log(a)+negative_hazard(logt)
            nh = negative_hazard(logh)
        else
            base = -log_hazard_mass(logt)
            logh = log(a)+log(base)
            nh = negative_hazard(logh)
            if (base <= 1.0_dp) then
                nh = -a*base
            else if (a <= huge(1.0_dp)/base) then
                nh = -a*base
            end if
        end if
    end subroutine exponweib_hazard

    pure function exponweib_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logt, logh, nh, product
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            if (params(1) > 1.0_dp) then
                if (params(2) > huge(1.0_dp)/params(1)) return
            end if
            product = params(1)*params(2)
            if (product < 1.0_dp) value = positive_infinity(0.0_dp)
            if (product == 1.0_dp) value = -log(params(4))
            return
        end if
        call exponweib_hazard(x,params,logt,logh,nh)
        if (.not. ieee_is_finite(nh)) return
        if (logt > log(huge(1.0_dp))) return
        value = log(params(1))+log(params(2))-log_distance(x,params(3))+nh
        if (logt >= log(epsilon(1.0_dp))-2.0_dp) &
            value = value+logt-log_expm1_from_log(logt)
    end function exponweib_lpdf

    pure function exponweib_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logt, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call exponweib_hazard(x,params,logt,logh,value)
    end function exponweib_lcdf

    pure function exponweib_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logt, logh, nh
        value = 0.0_dp
        if (x <= params(3)) return
        call exponweib_hazard(x,params,logt,logh,nh)
        value = hazard_complement(logh,nh)
    end function exponweib_lsf

    pure function exponweib_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive a, c, location, scale
        real(dp) :: value, logta, logtb, logha, loghb, nha, nhb, logdelta, logu, omitted
        value = exponweib_lcdf(b,params)
        if (a <= params(3)) return
        omitted = exp(exponweib_lcdf(a,params))+exp(exponweib_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call exponweib_hazard(a,params,logta,logha,nha)
        call exponweib_hazard(b,params,logtb,loghb,nhb)
        if (logtb < log(epsilon(1.0_dp))-2.0_dp) then
            logdelta = log(params(1))+log(params(2))+burr_logsoftplus( &
                log_distance(b,a)-log_distance(a,params(3)))
        else if (logta == negative_infinity(0.0_dp)) then
            logdelta = logha
        else
            value = negative_infinity(0.0_dp)
            if (logta > log(huge(1.0_dp))) return
            logdelta = weibull_log_increment(a,b,params(2:4))
            logu = -exp(logta)+log_hazard_mass(logdelta)-log_hazard_mass(logta)
            logdelta = log(params(1))+burr_logsoftplus(logu)
        end if
        value = nhb+log_hazard_mass(logdelta)
    end function exponweib_interval

    pure subroutine mielke_hazard(x, params, q, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp), intent(out) :: q !! -s*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log((k/s)*log1p(z**(-s)))
        real(dp), intent(out) :: nh !! negative reversed cumulative hazard
        real(dp) :: v, base, ratio
        q = -weibull_log_hazard(x,params(2:4))
        if (q >= 40.0_dp) then
            v = -lognormal_logz(x,params(2:4))
            logh = log(params(1))+log(v)
            nh = negative_hazard(logh)
            if (params(1) <= 1.0_dp) then
                nh = -params(1)*v
            else if (v <= huge(1.0_dp)/params(1)) then
                nh = -params(1)*v
            end if
        else
            logh = log(params(1))-log(params(2))+burr_logsoftplus(q)
            nh = negative_hazard(logh)
            ! Direct finite ratios/products retain ordinary-tail accuracy.
            if (params(2) < 1.0_dp) then
                if (params(1) > huge(1.0_dp)*params(2)) return
            end if
            ratio = params(1)/params(2)
            base = log1pexp(q)
            if (ratio < tiny(1.0_dp) .or. base < tiny(1.0_dp)) return
            if (ratio > 1.0_dp) then
                if (base > huge(1.0_dp)/ratio) return
            end if
            nh = -ratio*base
        end if
    end subroutine mielke_hazard

    pure function mielke_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, q, logh, nh
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(4))
            return
        end if
        call mielke_hazard(x,params,q,logh,nh)
        value = log(params(1))-log_distance(x,params(3))-log1pexp(-q)+nh
    end function mielke_lpdf

    pure function mielke_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, q, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call mielke_hazard(x,params,q,logh,value)
    end function mielke_lcdf

    pure function mielke_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, q, logh, nh
        value = 0.0_dp
        if (x <= params(3)) return
        call mielke_hazard(x,params,q,logh,nh)
        value = hazard_complement(logh,nh)
    end function mielke_lsf

    pure function mielke_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive k, s, location, scale
        real(dp) :: value, qa, qb, logh, nh, logspan, logdelta, logdiff, logu, delta, omitted
        value = mielke_lcdf(b,params)
        if (a <= params(3)) return
        omitted = exp(mielke_lcdf(a,params))+exp(mielke_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call mielke_hazard(a,params,qa,logh,nh)
        call mielke_hazard(b,params,qb,logh,value)
        logspan = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(3)))
        if (qb >= 40.0_dp) then
            logdiff = log(params(1))+logspan
        else
            logdelta = log(params(2))+logspan
            if (logdelta < log(0.5_dp)) then
                logu = log_expm1_from_log(logdelta)
                logdiff = burr_logsoftplus(-log1pexp(-qb)+logu)
            else if (qb >= 0.0_dp) then
                logdiff = logdelta
                if (logdelta < log(huge(1.0_dp))) then
                    delta = exp(logdelta)
                    logdiff = log(delta+log1pexp(-qa)-log1pexp(-qb))
                end if
            else if (qa <= 0.0_dp) then
                logu = -log1pexp(-qa)+log_hazard_mass(logdelta)
                logdiff = logu
                if (logu >= log(epsilon(1.0_dp))-2.0_dp) &
                    logdiff = log(-log1p_safe(-exp(logu)))
            else if (ieee_is_finite(qa)) then
                logdiff = log(log1pexp(qa)-log1pexp(qb))
            else
                logdiff = log(params(2))+log(-lognormal_logz(a,params(2:4)))
            end if
            logdiff = log(params(1))-log(params(2))+logdiff
        end if
        value = value+log_hazard_mass(logdiff)
    end function mielke_interval


    pure function exponpow_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite power shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [b,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,exponpow_lpdf,exponpow_lcdf, &
            exponpow_lsf,exponpow_interval)
    end function exponpow_censored_loglikelihood

    pure function exponpow_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! positive finite power shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -exponpow_censored_loglikelihood(data, b, loc, scale)
    end function exponpow_censored_nnlf

    pure function exponweib_censored_loglikelihood(data, a, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite exponentiation shape
        real(dp), intent(in) :: c !! positive finite Weibull shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [a,c,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,c) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,exponweib_lpdf,exponweib_lcdf, &
            exponweib_lsf,exponweib_interval)
    end function exponweib_censored_loglikelihood

    pure function exponweib_censored_nnlf(data, a, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! positive finite exponentiation shape
        real(dp), intent(in) :: c !! positive finite Weibull shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -exponweib_censored_loglikelihood(data, a, c, loc, scale)
    end function exponweib_censored_nnlf

    pure function mielke_censored_loglikelihood(data, k, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! positive finite first shape
        real(dp), intent(in) :: s !! positive finite second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [k,s,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(k,s) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,mielke_lpdf,mielke_lcdf, &
            mielke_lsf,mielke_interval)
    end function mielke_censored_loglikelihood

    pure function mielke_censored_nnlf(data, k, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: k !! positive finite first shape
        real(dp), intent(in) :: s !! positive finite second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -mielke_censored_loglikelihood(data, k, s, loc, scale)
    end function mielke_censored_nnlf


    pure function normal_log_mills(logz) result(value)
        implicit none
        real(dp), intent(in) :: logz !! log positive standardized point, z >= 20
        real(dp) :: value, inverse_square, term, series
        integer :: k
        ! Q(z)*z/phi(z), DLMF 7.12.1. At z >= 20 the first omitted
        ! term after twelve corrections is below 2e-20.
        inverse_square = exp(-2.0_dp*logz)
        term = 1.0_dp
        series = 1.0_dp
        do k = 1, 12
            term = -term*real(2*k-1,dp)*inverse_square
            series = series+term
        end do
        value = log(series)
    end function normal_log_mills

    pure function normal_log_hazard(z, logz) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point, possibly infinite
        real(dp), intent(in) :: logz !! log absolute standardized point
        real(dp) :: value
        if (z >= 20.0_dp) then
            value = logz-normal_log_mills(logz)
        else if (z < -20.0_dp) then
            value = negative_hazard(2.0_dp*logz-scifort_log_two)-scifort_log_sqrt_two_pi
        else
            value = normal_logpdf(z)-normal_logsf(z)
        end if
    end function normal_log_hazard

    pure function power_normal_logh(z, logz, c) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point, possibly infinite
        real(dp), intent(in) :: logz !! log absolute standardized point
        real(dp), intent(in) :: c !! positive finite power
        real(dp) :: value, correction
        if (z >= 20.0_dp) then
            correction = 2.0_dp*(logz+scifort_log_sqrt_two_pi-normal_log_mills(logz))
            correction = correction*exp(-2.0_dp*logz)
            value = log(c)+2.0_dp*logz-scifort_log_two+log1p_safe(correction)
        else if (z < -8.0_dp) then
            ! -log(1-Phi(z)) / Phi(z) differs from one by < 4e-16.
            value = log(c)+normal_logcdf(z)
        else
            value = log(c)+log(-normal_logsf(z))
        end if
    end function power_normal_logh

    pure subroutine power_coordinates(x, params, z, logz)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp), intent(out) :: z !! standardized point, possibly infinite
        real(dp), intent(out) :: logz !! log absolute standardized point
        z = location_scale_z(x,params(2:3))
        logz = cauchy_log_distance(x,params(2:3))
        if (z /= 0.0_dp .and. ieee_is_finite(z)) logz = log(abs(z))
    end subroutine power_coordinates

    pure function powernorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz
        call power_coordinates(x,params,z,logz)
        value = log(params(1))+normal_log_hazard(z,logz)-log(params(3)) &
            +negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powernorm_lpdf

    pure function powernorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz
        call power_coordinates(x,params,z,logz)
        value = log_hazard_mass(power_normal_logh(z,logz,params(1)))
    end function powernorm_lcdf

    pure function powernorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz
        call power_coordinates(x,params,z,logz)
        value = negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powernorm_lsf

    pure function power_normal_small_interval(z, logz, logwidth, c) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized lower endpoint, possibly infinite
        real(dp), intent(in) :: logz !! log absolute standardized lower endpoint
        real(dp), intent(in) :: logwidth !! log width in normal coordinates
        real(dp), intent(in) :: c !! positive finite power
        real(dp) :: value, logdelta, q, nodes(3), ratios(3), w, hazard, correction
        integer :: i
        value = quiet_nan(0.0_dp)
        if (z >= 20.0_dp) then
            if (logwidth-logz > log(0.001_dp)) return
        else
            if (logwidth+max(0.0_dp,logz) > log(0.001_dp)) return
        end if
        hazard = normal_log_hazard(z,logz)
        if (.not. ieee_is_finite(hazard)) then
            value = negative_infinity(0.0_dp)
            return
        end if
        if (z >= 1.0e8_dp) then
            ! The Mills correction to the integrated hazard is <= 1/z**2.
            q = exp(logwidth-logz)
            correction = log1p_safe(0.5_dp*q)
        else
            ! Integrate the normal hazard, rather than subtract nearly equal
            ! log survivals. Three-point Gauss-Legendre, DLMF 3.5(v).
            nodes = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            do i = 1, 3
                w = z+nodes(i)*exp(logwidth)
                if (w == 0.0_dp) then
                    ratios(i) = normal_log_hazard(w,negative_infinity(0.0_dp))-hazard
                else
                    ratios(i) = normal_log_hazard(w,log(abs(w)))-hazard
                end if
            end do
            correction = log1p_safe((4.0_dp/9.0_dp)*expm1_safe(ratios(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(ratios(2))+expm1_safe(ratios(3))))
        end if
        logdelta = log(c)+logwidth+hazard+correction
        value = negative_hazard(power_normal_logh(z,logz,c))+log_hazard_mass(logdelta)
    end function power_normal_small_interval

    pure function powernorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! power, location, positive scale
        real(dp) :: value, z, logz, omitted
        call power_coordinates(a,params,z,logz)
        value = power_normal_small_interval(z,logz, &
            log_distance(b,a)-log(params(3)),params(1))
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,powernorm_lcdf,powernorm_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(powernorm_lcdf(a,params),powernorm_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function powernorm_interval

    pure subroutine powerlog_coordinates(x, params, z, logz)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp), intent(out) :: z !! normal coordinate, possibly infinite
        real(dp), intent(out) :: logz !! log absolute normal coordinate
        real(dp) :: lz
        lz = lognormal_logz(x,params(2:4))
        z = lognormal_standardized(x,params(2:4))
        logz = negative_infinity(0.0_dp)
        if (lz /= 0.0_dp) logz = log(abs(lz))-log(params(2))
        if (z /= 0.0_dp .and. ieee_is_finite(z)) logz = log(abs(z))
    end subroutine powerlog_coordinates

    pure function powerlognorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call powerlog_coordinates(x,params,z,logz)
        value = log(params(1))+normal_log_hazard(z,logz)-log(params(2)) &
            -log_distance(x,params(3))+negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powerlognorm_lpdf

    pure function powerlognorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call powerlog_coordinates(x,params,z,logz)
        value = log_hazard_mass(power_normal_logh(z,logz,params(1)))
    end function powerlognorm_lcdf

    pure function powerlognorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz
        value = 0.0_dp
        if (x <= params(3)) return
        call powerlog_coordinates(x,params,z,logz)
        value = negative_hazard(power_normal_logh(z,logz,params(1)))
    end function powerlognorm_lsf

    pure function powerlognorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! power, lognormal shape, location, scale
        real(dp) :: value, z, logz, logr, logspan, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(3)) return
        if (a <= params(3)) then
            value = powerlognorm_lcdf(b,params)
            return
        end if
        logr = log_distance(b,a)-log_distance(a,params(3))
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        call powerlog_coordinates(a,params,z,logz)
        value = power_normal_small_interval(z,logz,logspan-log(params(2)),params(1))
        if (ieee_is_nan(value)) then
            value = tail_interval(a,b,params,powerlognorm_lcdf,powerlognorm_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(powerlognorm_lcdf(a,params),powerlognorm_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function powerlognorm_interval

    pure function skewcauchy_logside(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value
        value = log1p_safe(params(1))
        if (x < params(2)) value = log1p_safe(-params(1))
    end function skewcauchy_logside

    pure function skewcauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, logz
        logz = cauchy_log_distance(x,params(2:3))-skewcauchy_logside(x,params)
        value = -log(scifort_pi)-log(params(3))-log1pexp(2.0_dp*logz)
    end function skewcauchy_lpdf

    pure function skewcauchy_outer_tail(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite evaluation point
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, logh, logz
        logh = skewcauchy_logside(x,params)
        logz = cauchy_log_distance(x,params(2:3))-logh
        value = logh-log(scifort_pi)+cauchy_logatan_exp(-logz)
    end function skewcauchy_outer_tail

    pure function skewcauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, center, logh, logz
        center = log1p_safe(-params(1))-scifort_log_two
        if (params(1) < 0.0_dp) center = log1p_safe(-0.5_dp*(1.0_dp+params(1)))
        value = center
        if (x == params(2)) return
        value = skewcauchy_outer_tail(x,params)
        if (x > params(2)) then
            if (value < -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                logh = log1p_safe(params(1))
                logz = cauchy_log_distance(x,params(2:3))-logh
                value = bounded_logadd(center,logh-log(scifort_pi)+cauchy_logatan_exp(logz))
            end if
        end if
    end function skewcauchy_lcdf

    pure function skewcauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, center, logh, logz
        center = log1p_safe(params(1))-scifort_log_two
        if (params(1) > 0.0_dp) center = log1p_safe(-0.5_dp*(1.0_dp-params(1)))
        value = center
        if (x == params(2)) return
        value = skewcauchy_outer_tail(x,params)
        if (x < params(2)) then
            if (value < -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                logh = log1p_safe(-params(1))
                logz = cauchy_log_distance(x,params(2:3))-logh
                value = bounded_logadd(center,logh-log(scifort_pi)+cauchy_logatan_exp(logz))
            end if
        end if
    end function skewcauchy_lsf

    pure function skewcauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! skewness, location, positive scale
        real(dp) :: value, logh, loga, logb, ratio, left, right, omitted
        if (b <= params(2) .or. a >= params(2)) then
            logh = log1p_safe(params(1))
            if (b <= params(2)) logh = log1p_safe(-params(1))
            loga = cauchy_log_distance(a,params(2:3))-logh
            logb = cauchy_log_distance(b,params(2:3))-logh
            ratio = log_distance(b,a)-log(params(3))-logh-log1pexp(loga+logb)
            value = logh-log(scifort_pi)+cauchy_logatan_exp(ratio)
        else
            logh = log1p_safe(-params(1))
            loga = cauchy_log_distance(a,params(2:3))-logh
            left = logh-log(scifort_pi)+cauchy_logatan_exp(loga)
            logh = log1p_safe(params(1))
            logb = cauchy_log_distance(b,params(2:3))-logh
            right = logh-log(scifort_pi)+cauchy_logatan_exp(logb)
            value = min(0.0_dp,bounded_logadd(left,right))
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(skewcauchy_lcdf(a,params),skewcauchy_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function skewcauchy_interval


    pure function powernorm_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,powernorm_lpdf,powernorm_lcdf, &
            powernorm_lsf,powernorm_interval)
    end function powernorm_censored_loglikelihood

    pure function powernorm_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -powernorm_censored_loglikelihood(data, c, loc, scale)
    end function powernorm_censored_nnlf

    pure function powerlognorm_censored_loglikelihood(data, c, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in) :: s !! positive finite lognormal shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(4)
        params = [c,s,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(c,s) <= 0.0_dp) return
        if (params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,powerlognorm_lpdf,powerlognorm_lcdf, &
            powerlognorm_lsf,powerlognorm_interval)
    end function powerlognorm_censored_loglikelihood

    pure function powerlognorm_censored_nnlf(data, c, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! positive finite normal power
        real(dp), intent(in) :: s !! positive finite lognormal shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -powerlognorm_censored_loglikelihood(data, c, s, loc, scale)
    end function powerlognorm_censored_nnlf

    pure function skewcauchy_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! skewness, strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (abs(a) >= 1.0_dp) return
        if (params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,skewcauchy_lpdf,skewcauchy_lcdf, &
            skewcauchy_lsf,skewcauchy_interval)
    end function skewcauchy_censored_loglikelihood

    pure function skewcauchy_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! skewness, strictly between -1 and 1
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -skewcauchy_censored_loglikelihood(data, a, loc, scale)
    end function skewcauchy_censored_nnlf


    pure function invgamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logy = -lognormal_logz(x,params)
        value = loggamma_lpdf(logy,[params(1),0.0_dp,1.0_dp]) &
            -log_distance(x,params(2))
    end function invgamma_lpdf

    pure function invgamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logy = -lognormal_logz(x,params)
        value = loggamma_lsf(logy,[params(1),0.0_dp,1.0_dp])
    end function invgamma_lcdf

    pure function invgamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logy
        value = 0.0_dp
        if (x <= params(2)) return
        logy = -lognormal_logz(x,params)
        value = loggamma_lcdf(logy,[params(1),0.0_dp,1.0_dp])
    end function invgamma_lsf

    pure function invgamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, location, scale
        real(dp) :: value, logwidth, logq, logy, q, t(3), r(3), correction, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = invgamma_lcdf(b,params)
            return
        end if
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(2))
        logy = -lognormal_logz(a,params)
        if (logq+max(log(params(1)+1.0_dp),logy) <= log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = -(params(1)+1.0_dp)*log1p_safe(t*q)+exp(logy+logq)*t/(1.0_dp+t*q)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = invgamma_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,invgamma_lcdf,invgamma_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(invgamma_lcdf(a,params),invgamma_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function invgamma_interval

    pure function maxwell_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        logy = 2.0_dp*lognormal_logz(x,[1.0_dp,params(1),params(2)])-scifort_log_two
        value = loggamma_lpdf(logy,[1.5_dp,0.0_dp,1.0_dp]) &
            -log_distance(x,params(1))+scifort_log_two
    end function maxwell_lpdf

    pure function maxwell_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logy
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        logy = 2.0_dp*lognormal_logz(x,[1.0_dp,params(1),params(2)])-scifort_log_two
        value = loggamma_lcdf(logy,[1.5_dp,0.0_dp,1.0_dp])
    end function maxwell_lcdf

    pure function maxwell_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logy
        value = 0.0_dp
        if (x <= params(1)) return
        logy = 2.0_dp*lognormal_logz(x,[1.0_dp,params(1),params(2)])-scifort_log_two
        value = loggamma_lsf(logy,[1.5_dp,0.0_dp,1.0_dp])
    end function maxwell_lsf

    pure function maxwell_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logwidth, logq, logy, q, t(3), r(3), correction, omitted
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = maxwell_lcdf(b,params)
            return
        end if
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(1))
        logy = 2.0_dp*lognormal_logz(a,[1.0_dp,params(1),params(2)])-scifort_log_two
        if (logq+max(0.0_dp,logy) <= log(0.001_dp)) then
            q = exp(logq)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = 2.0_dp*log1p_safe(t*q)-exp(logy+logq)*t*(2.0_dp+t*q)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = maxwell_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,maxwell_lcdf,maxwell_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(maxwell_lcdf(a,params),maxwell_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function maxwell_interval

    pure subroutine betaprime_coordinates(x, params, u, v, logu, logv)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp), intent(out) :: u !! z/(1+z), possibly rounded to zero or one
        real(dp), intent(out) :: v !! 1/(1+z), retained independently
        real(dp), intent(out) :: logu !! log transformed lower coordinate
        real(dp), intent(out) :: logv !! log transformed upper coordinate
        real(dp) :: logz, e
        logz = lognormal_logz(x,[1.0_dp,params(3),params(4)])
        logu = -log1pexp(-logz)
        logv = -log1pexp(logz)
        if (logz <= 0.0_dp) then
            e = exp(logz)
            u = e/(1.0_dp+e)
            v = 1.0_dp/(1.0_dp+e)
        else
            e = exp(-logz)
            u = 1.0_dp/(1.0_dp+e)
            v = e/(1.0_dp+e)
        end if
    end subroutine betaprime_coordinates

    pure function betaprime_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, u, v, logu, logv
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = log(params(2))-log(params(4))
            return
        end if
        call betaprime_coordinates(x,params,u,v,logu,logv)
        if (min(u,v) > 0.0_dp) then
            value = log_beta_kernel(params(1),params(2),u,v)
        else if (u == 0.0_dp) then
            value = betaprime_endpoint_tail(params(1),params(2),logu)+log(params(1))
        else
            value = betaprime_endpoint_tail(params(2),params(1),logv)+log(params(2))
        end if
        value = value-log_distance(x,params(3))
    end function betaprime_lpdf

    pure function betaprime_endpoint_tail(a, b, logu) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive shape at a coordinate rounded to zero
        real(dp), intent(in) :: b !! positive opposite shape; a+b finite
        real(dp), intent(in) :: logu !! finite log coordinate below the subnormal range
        real(dp) :: value, normalization, factor
        if (a < 1.0_dp) then
            value = beta_small_tail(a,b,logu)
        else
            ! Combine normalization before multiplication if a*logu overflows.
            normalization = beta_small_tail(a,b,0.0_dp)
            factor = logu+normalization/a
            value = negative_infinity(0.0_dp)
            if (a > huge(1.0_dp)/abs(factor)) return
            value = a*factor
        end if
    end function betaprime_endpoint_tail

    pure subroutine betaprime_logtails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp), intent(out) :: logp !! log lower probability
        real(dp), intent(out) :: logq !! log upper probability
        real(dp) :: u, v, logu, logv, p, q
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        if (x <= params(3)) return
        call betaprime_coordinates(x,params,u,v,logu,logv)
        if (u == 0.0_dp) then
            logp = betaprime_endpoint_tail(params(1),params(2),logu)
            logq = log_difference(0.0_dp,logp)
        else if (v == 0.0_dp) then
            logq = betaprime_endpoint_tail(params(2),params(1),logv)
            logp = log_difference(0.0_dp,logq)
        else
            call incomplete_beta_xy(params(1),params(2),u,v,p,q,logp,logq)
        end if
    end subroutine betaprime_logtails

    pure function betaprime_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, logq
        call betaprime_logtails(x,params,value,logq)
    end function betaprime_lcdf

    pure function betaprime_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, logp
        call betaprime_logtails(x,params,logp,value)
    end function betaprime_lsf

    pure function betaprime_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shapes with finite sum, location, scale
        real(dp) :: value, logwidth, logq, logq2, u, v, logu, logv, q, q2
        real(dp) :: t(3), r(3), correction, omitted, shape_sum
        value = negative_infinity(0.0_dp)
        if (b <= params(3)) return
        if (a <= params(3)) then
            value = betaprime_lcdf(b,params)
            return
        end if
        call betaprime_coordinates(a,params,u,v,logu,logv)
        logwidth = log_distance(b,a)
        logq = logwidth-log_distance(a,params(3))
        logq2 = logq+logu
        shape_sum = params(1)+params(2)
        if (logq+log(max(1.0_dp,abs(params(1)-1.0_dp))) <= log(0.001_dp) &
                .and. logq2+log(shape_sum) <= log(0.001_dp)) then
            q = exp(logq)
            q2 = exp(logq2)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = (params(1)-1.0_dp)*log1p_safe(t*q)-shape_sum*log1p_safe(t*q2)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = betaprime_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,betaprime_lcdf,betaprime_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(betaprime_lcdf(a,params),betaprime_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function betaprime_interval


    pure function invgamma_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,params(3)) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,invgamma_lpdf,invgamma_lcdf, &
            invgamma_lsf,invgamma_interval)
    end function invgamma_censored_loglikelihood

    pure function invgamma_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -invgamma_censored_loglikelihood(data, a, loc, scale)
    end function invgamma_censored_nnlf

    pure function betaprime_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in) :: b !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(4)
        params = [a,b,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (min(a,b,params(4)) <= 0.0_dp) return
        if (0.5_dp*a+0.5_dp*b > 0.5_dp*huge(1.0_dp)) return
        value = censored_loglikelihood(data,params,betaprime_lpdf,betaprime_lcdf, &
            betaprime_lsf,betaprime_interval)
    end function betaprime_censored_loglikelihood

    pure function betaprime_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive first shape
        real(dp), intent(in) :: b !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -betaprime_censored_loglikelihood(data, a, b, loc, scale)
    end function betaprime_censored_nnlf

    pure function maxwell_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp,1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,maxwell_lpdf,maxwell_lcdf, &
            maxwell_lsf,maxwell_interval)
    end function maxwell_censored_loglikelihood

    pure function maxwell_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -maxwell_censored_loglikelihood(data, loc, scale)
    end function maxwell_censored_nnlf


    pure function dgamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = gamma_lpdf(y,p)-scifort_log_two
    end function dgamma_lpdf

    pure function dgamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        if (x == params(2)) then
            value = -scifort_log_two
            return
        end if
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = gamma_lsf(y,p)-scifort_log_two
        if (.not. (x < params(2))) value = log1p_safe(-exp(value))
    end function dgamma_lcdf

    pure function dgamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = gamma_lsf(y,p)-scifort_log_two
        if (.not. (x >= params(2))) value = log1p_safe(-exp(value))
    end function dgamma_lsf

    pure function dgamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, p(3), left, right, omitted
        p = [params(1),-params(2),params(3)]
        if (b <= params(2)) then
            value = gamma_interval(-b,-a,p)-scifort_log_two
        else if (a >= params(2)) then
            value = gamma_interval(a,b,params)-scifort_log_two
        else
            left = gamma_lcdf(-a,p)
            right = gamma_lcdf(b,params)
            value = bounded_logadd(left,right)-scifort_log_two
            ! Preserve both excluded tails for nearly full intervals.
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(gamma_lsf(-a,p), &
                    gamma_lsf(b,params))-scifort_log_two
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function dgamma_interval

    pure function dweibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = radial_weibull_lpdf(y,p)-scifort_log_two
    end function dweibull_lpdf

    pure function dweibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        if (x == params(2)) then
            value = -scifort_log_two
            return
        end if
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = radial_weibull_lsf(y,p)-scifort_log_two
        if (.not. (x < params(2))) value = log1p_safe(-exp(value))
    end function dweibull_lcdf

    pure function dweibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, y, p(3)
        y = x
        p = params
        if (x < params(2)) then
            y = -x
            p(2) = -params(2)
        end if
        value = radial_weibull_lsf(y,p)-scifort_log_two
        if (.not. (x >= params(2))) value = log1p_safe(-exp(value))
    end function dweibull_lsf

    pure function dweibull_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, p(3), left, right, omitted
        p = [params(1),-params(2),params(3)]
        if (b <= params(2)) then
            value = weibull_interval(-b,-a,p)-scifort_log_two
        else if (a >= params(2)) then
            value = weibull_interval(a,b,params)-scifort_log_two
        else
            left = radial_weibull_lcdf(-a,p)
            right = radial_weibull_lcdf(b,params)
            value = bounded_logadd(left,right)-scifort_log_two
            ! Preserve both excluded tails for nearly full intervals.
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(radial_weibull_lsf(-a,p), &
                    radial_weibull_lsf(b,params))-scifort_log_two
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function dweibull_interval

    pure function radial_weibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logh
        if (x == params(2)) then
            value = weibull_logpdf(x,params(1),params(2),params(3))
            return
        end if
        logh = weibull_log_hazard(x,params)
        value = negative_infinity(0.0_dp)
        if (.not. ieee_is_finite(logh)) return
        if (logh > log(huge(1.0_dp))) return
        value = log(params(1))+logh-log_distance(x,params(2))-exp(logh)
    end function radial_weibull_lpdf

    pure function radial_weibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x == params(2)) return
        value = log_exponential_mass(weibull_log_hazard(x,params))
    end function radial_weibull_lcdf

    pure function radial_weibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point at or above location
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x == params(2)) return
        value = negative_hazard(weibull_log_hazard(x,params))
    end function radial_weibull_lsf

    pure function loggamma_left_kernel(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point with standardized coordinate below log(tiny)
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, product
        ! Retain c*z even if z itself overflows for a very small shape.
        z = location_scale_z(x,params(2:3))
        if (ieee_is_finite(z)) then
            product = negative_infinity(0.0_dp)
            if (params(1) <= huge(1.0_dp)/abs(z)) product = params(1)*z
        else
            product = negative_hazard(log(params(1))+log_distance(params(2),x)-log(params(3)))
        end if
        if (params(1) < 10.0_dp) then
            value = product-log_gamma_one_plus(params(1))
        else
            z = location_scale_z(x,params(2:3))
            value = negative_infinity(0.0_dp)
            if (ieee_is_finite(z)) value = gamma_small_kernel(params(1),z)
        end if
    end function loggamma_left_kernel

    pure function loggamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, u, factor, term
        integer :: k
        z = location_scale_z(x,params(2:3))
        value = negative_infinity(0.0_dp)
        if (z < log(tiny(1.0_dp))) then
            value = log(params(1))+loggamma_left_kernel(x,params)-log(params(3))
        else if (params(1) < 10.0_dp) then
            if (z > log(huge(1.0_dp))) return
            value = log(params(1))+params(1)*z-exp(z)-log_gamma_one_plus(params(1)) &
                -log(params(3))
        else
            ! Gamma saddle form retains the cancellation at the log-gamma mode.
            u = z-log(params(1))
            if (u > log(huge(1.0_dp))) return
            if (abs(u) < 0.1_dp) then
                term = (0.5_dp*u)*u
                factor = term
                do k = 3, 20
                    term = term*u/real(k,dp)
                    factor = factor+term
                end do
            else
                factor = expm1_safe(u)-u
            end if
            if (factor > 1.0_dp) then
                if (params(1) > huge(1.0_dp)/factor) return
            end if
            value = 0.5_dp*log(params(1))-scifort_log_sqrt_two_pi &
                -params(1)*factor-stirling_remainder(params(1))-log(params(3))
        end if
    end function loggamma_lpdf

    pure function loggamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z
        z = location_scale_z(x,params(2:3))
        if (z < log(tiny(1.0_dp))) then
            value = loggamma_left_kernel(x,params)
        else if (z > log(huge(1.0_dp))) then
            value = 0.0_dp
        else
            value = log_gammainc(params(1),exp(z))
        end if
    end function loggamma_lcdf

    pure function loggamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logp, u
        z = location_scale_z(x,params(2:3))
        if (z < log(tiny(1.0_dp))) then
            logp = loggamma_left_kernel(x,params)
            value = 0.0_dp
            if (ieee_is_finite(logp)) value = log_exponential_mass(bounded_log(-logp))
        else if (z > log(huge(1.0_dp))) then
            ! For exp(z)>huge, the endpoint expansion has negligible relative
            ! correction O(exp(z)/(exp(z)-c)**2). Retain a representable log tail.
            u = z-log(params(1))
            value = loggamma_lpdf(x,params)
            if (ieee_is_finite(value)) value = value+log(params(3)) &
                -log(params(1))-log_expm1_from_log(log(u))
        else
            value = log_gammaincc(params(1),exp(z))
        end if
    end function loggamma_lsf

    pure function loggamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logwidth, logscaled, d, t(3), r(3), correction, omitted
        z = location_scale_z(a,params(2:3))
        logwidth = log_distance(b,a)
        logscaled = logwidth-log(params(3))
        ! Bound c*d and exp(z)*d independently; the density varies by <0.002.
        if (logscaled+max(log(params(1)),z) <= log(0.001_dp)) then
            d = 0.0_dp
            if (z > log(tiny(1.0_dp))) d = exp(logscaled)
            t = [0.5_dp,0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = t*exp(log(params(1))+logscaled)
            if (z > log(tiny(1.0_dp))) then
                ! expm1(t*d)/d remains t when d underflows.
                if (d > 0.0_dp) then
                    r = r-exp(z+logscaled)*(expm1_safe(t*d)/d)
                else
                    r = r-t*exp(z+logscaled)
                end if
            end if
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) &
                +(5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = loggamma_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,loggamma_lcdf,loggamma_lsf)
            if (value > -scifort_log_two) then
                omitted = bounded_logadd(loggamma_lcdf(a,params),loggamma_lsf(b,params))
                value = log1p_safe(-exp(omitted))
            end if
        end if
    end function loggamma_interval


    pure function dgamma_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [a,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (a <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,dgamma_lpdf,dgamma_lcdf, &
            dgamma_lsf,dgamma_interval)
    end function dgamma_censored_loglikelihood

    pure function dgamma_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -dgamma_censored_loglikelihood(data,a,loc,scale)
    end function dgamma_censored_nnlf

    pure function dweibull_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,dweibull_lpdf,dweibull_lcdf, &
            dweibull_lsf,dweibull_interval)
    end function dweibull_censored_loglikelihood

    pure function dweibull_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -dweibull_censored_loglikelihood(data,c,loc,scale)
    end function dweibull_censored_nnlf

    pure function loggamma_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,loggamma_lpdf,loggamma_lcdf, &
            loggamma_lsf,loggamma_interval)
    end function loggamma_censored_loglikelihood

    pure function loggamma_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -loggamma_censored_loglikelihood(data,c,loc,scale)
    end function loggamma_censored_nnlf


    pure function bounded_logadd(a, b) result(value)
        implicit none
        real(dp), intent(in) :: a !! log first nonnegative term, possibly negative infinity
        real(dp), intent(in) :: b !! log second nonnegative term, possibly negative infinity
        real(dp) :: value, big
        big = max(a,b)
        value = big
        if (.not. ieee_is_finite(big)) return
        value = big+log1p_safe(exp(min(a,b)-big))
    end function bounded_logadd

    pure function bounded_log(x) result(value)
        implicit none
        real(dp), intent(in) :: x !! nonnegative dimensionless shape or width
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x > 0.0_dp) value = log(x)
    end function bounded_log

    pure function bounded_left(logz, logw, t) result(left)
        implicit none
        real(dp), intent(in) :: logz !! log coordinate strictly between zero and one
        real(dp), intent(in) :: logw !! log complementary coordinate, retained independently
        real(dp), intent(in) :: t !! breakpoint in [0,1]
        logical :: left
        left = .false.
        if (t <= 0.0_dp) return
        left = .true.
        if (t >= 1.0_dp) return
        if (t <= 0.5_dp) then
            left = logz < log(t)
        else
            left = logw > log1p_safe(-t)
        end if
    end function bounded_left

    pure function bounded_right(logz, logw, t) result(right)
        implicit none
        real(dp), intent(in) :: logz !! log coordinate strictly between zero and one
        real(dp), intent(in) :: logw !! log complementary coordinate, retained independently
        real(dp), intent(in) :: t !! breakpoint in [0,1]
        logical :: right
        right = .false.
        if (t >= 1.0_dp) return
        right = .true.
        if (t <= 0.0_dp) return
        if (t <= 0.5_dp) then
            right = logz > log(t)
        else
            right = logw < log1p_safe(-t)
        end if
    end function bounded_right

    pure function bounded_gap(t, logz, logw) result(value)
        implicit none
        real(dp), intent(in) :: t !! coordinate in (0,1], no smaller than z
        real(dp), intent(in) :: logz !! log lower coordinate, possibly negative infinity
        real(dp), intent(in) :: logw !! log complementary coordinate, retained independently
        real(dp) :: value, large, small
        if (t == 1.0_dp) then
            value = logw
        else
            if (t > 0.5_dp) then
                large = logw
                small = log1p_safe(-t)
            else
                large = log(t)
                small = logz
            end if
            ! Rounded coordinates may straddle a breakpoint by one ulp.
            value = log_difference(large,min(large,small))
        end if
    end function bounded_gap

    pure function bounded_right_gap(t, logz, logw) result(value)
        implicit none
        real(dp), intent(in) :: t !! breakpoint in [0,1), no greater than z
        real(dp), intent(in) :: logz !! log lower coordinate, retained independently
        real(dp), intent(in) :: logw !! log complementary coordinate
        real(dp) :: value
        if (t <= 0.5_dp) then
            value = log_difference(logz,min(logz,bounded_log(t)))
        else
            value = bounded_gap(1.0_dp-t,logw,logz)
        end if
    end function bounded_right_gap

    pure function trapezoid_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, left, right, z, w, lz, lw, c, d
        c = params(1)
        d = params(2)
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp .and. c > 0.0_dp) return
        if (right == 0.0_dp .and. d < 1.0_dp) return
        value = scifort_log_two-log1p_safe(d-c)-log(params(4))
        if (left == 0.0_dp .or. right == 0.0_dp) return
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        if (bounded_left(lz,lw,c)) value = value+lz-log(c)
        if (bounded_right(lz,lw,d)) value = value+lw-log1p_safe(-d)
    end function trapezoid_lpdf

    pure subroutine trapezoid_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp), intent(out) :: logcdf !! log lower-tail probability
        real(dp), intent(out) :: logsf !! log upper-tail probability
        real(dp) :: left, right, z, w, lz, lw, c, d, lc, lr, lh, gap, terms
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logcdf = 0.0_dp
            logsf = negative_infinity(0.0_dp)
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        c = params(1)
        d = params(2)
        lc = bounded_log(c)
        lr = bounded_log(1.0_dp-d)
        lh = scifort_log_two-log1p_safe(d-c)
        if (bounded_left(lz,lw,c)) then
            logcdf = lh-scifort_log_two+2.0_dp*lz-lc
            gap = bounded_gap(c,lz,lw)
            terms = bounded_logadd(lr,scifort_log_two+bounded_log(d-c))
            terms = bounded_logadd(terms,gap+log1p_safe(exp(lz-lc)))
            logsf = lh-scifort_log_two+terms
        else if (bounded_right(lz,lw,d)) then
            logsf = lh-scifort_log_two+2.0_dp*lw-lr
            gap = bounded_right_gap(d,lz,lw)
            terms = bounded_logadd(lc,scifort_log_two+bounded_log(d-c))
            terms = bounded_logadd(terms,gap+log1p_safe(exp(lw-lr)))
            logcdf = lh-scifort_log_two+terms
        else
            logcdf = lh+lz+log1p_safe(-0.5_dp*exp(lc-lz))
            logsf = lh+lw+log1p_safe(-0.5_dp*exp(lr-lw))
        end if
        if (logcdf <= logsf) then
            logsf = log1p_safe(-exp(logcdf))
        else
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine trapezoid_tails

    pure function trapezoid_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left-censoring limit
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, unused
        call trapezoid_tails(x,params,value,unused)
    end function trapezoid_lcdf

    pure function trapezoid_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right-censoring limit
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, unused
        call trapezoid_tails(x,params,unused,value)
    end function trapezoid_lsf

    pure function trapezoid_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! plateau start/end in [0,1], location, width
        real(dp) :: value, left, right, z, w, lza, lwa, lzb, lwb, ld, lc, lr, lh
        real(dp) :: omitted, deficit, gap
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = trapezoid_lcdf(b,params)
            return
        end if
        call beta_coordinates(a,params,left,right,z,w,lza,lwa)
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = trapezoid_lsf(a,params)
            return
        end if
        call beta_coordinates(b,params,left,right,z,w,lzb,lwb)
        omitted = exp(trapezoid_lcdf(a,params))+exp(trapezoid_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        ld = log_distance(b,a)-log(params(4))
        lc = bounded_log(params(1))
        lr = bounded_log(1.0_dp-params(2))
        lh = scifort_log_two-log1p_safe(params(2)-params(1))
        if (.not. bounded_right(lzb,lwb,params(1))) then
            value = lh+ld+bounded_logadd(lza,lzb)-scifort_log_two-lc
        else if (.not. bounded_left(lza,lwa,params(2))) then
            value = lh+ld+bounded_logadd(lwa,lwb)-scifort_log_two-lr
        else
            ! Integrate the constant plateau height minus the two triangular deficits.
            deficit = 0.0_dp
            if (bounded_left(lza,lwa,params(1))) then
                gap = bounded_gap(params(1),lza,lwa)
                deficit = exp(2.0_dp*gap-scifort_log_two-lc-ld)
            end if
            if (bounded_right(lzb,lwb,params(2))) then
                gap = bounded_right_gap(params(2),lzb,lwb)
                deficit = deficit+exp(2.0_dp*gap-scifort_log_two-lr-ld)
            end if
            value = lh+ld+log1p_safe(-deficit)
        end if
    end function trapezoid_interval

    pure function bradford_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, left, right, z, w, lz, lw, lc
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        lc = log(params(1))
        value = -log(params(4))
        if (lc < log(epsilon(1.0_dp))-2.0_dp) return
        lz = negative_infinity(0.0_dp)
        if (right == 0.0_dp) lz = 0.0_dp
        if (left > 0.0_dp .and. right > 0.0_dp) then
            call beta_coordinates(x,params,left,right,z,w,lz,lw)
        end if
        value = lc-log(params(4))-burr_logsoftplus(lc)-log1pexp(lc+lz)
    end function bradford_lpdf

    pure subroutine bradford_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp), intent(out) :: logcdf !! log lower-tail probability
        real(dp), intent(out) :: logsf !! log upper-tail probability
        real(dp) :: left, right, z, w, lz, lw, lc, ln
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logcdf = 0.0_dp
            logsf = negative_infinity(0.0_dp)
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,lz,lw)
        lc = log(params(1))
        if (lc < log(epsilon(1.0_dp))-2.0_dp) then
            logcdf = lz
            logsf = lw
        else
            ln = burr_logsoftplus(lc)
            logcdf = burr_logsoftplus(lc+lz)-ln
            logsf = burr_logsoftplus(lc+lw-log1pexp(lc+lz))-ln
        end if
        if (logcdf <= logsf) then
            logsf = log1p_safe(-exp(logcdf))
        else
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine bradford_tails

    pure function bradford_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left-censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, unused
        call bradford_tails(x,params,value,unused)
    end function bradford_lcdf

    pure function bradford_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right-censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, unused
        call bradford_tails(x,params,unused,value)
    end function bradford_lsf

    pure function bradford_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, unused slot, location, width
        real(dp) :: value, left, right, z, w, lz, lw, lc, ld, omitted
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = bradford_lcdf(b,params)
            return
        end if
        call beta_coordinates(a,params,left,right,z,w,lz,lw)
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = bradford_lsf(a,params)
            return
        end if
        omitted = exp(bradford_lcdf(a,params))+exp(bradford_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        lc = log(params(1))
        ld = log_distance(b,a)-log(params(4))
        value = ld
        if (lc < log(epsilon(1.0_dp))-2.0_dp) return
        value = burr_logsoftplus(lc+ld-log1pexp(lc+lz))-burr_logsoftplus(lc)
    end function bradford_interval

    pure function bradford_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive Bradford shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [c,1.0_dp,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(3),params(4)])) return
        value = censored_loglikelihood(data,params,bradford_lpdf,bradford_lcdf, &
            bradford_lsf,bradford_interval)
    end function bradford_censored_loglikelihood

    pure function bradford_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive Bradford shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -bradford_censored_loglikelihood(data,c,loc,scale)
    end function bradford_censored_nnlf

    pure function triang_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite standardized mode in [0,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = trapezoid_censored_loglikelihood(data,c,c,loc,scale)
    end function triang_censored_loglikelihood

    pure function triang_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite standardized mode in [0,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -triang_censored_loglikelihood(data,c,loc,scale)
    end function triang_censored_nnlf

    pure function trapezoid_censored_loglikelihood(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite plateau start in [0,1]
        real(dp), intent(in) :: d !! finite plateau end in [c,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [c,d,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c < 0.0_dp .or. d < c .or. d > 1.0_dp) return
        if (.not. beta_parameters_valid([1.0_dp,1.0_dp,params(3),params(4)])) return
        value = censored_loglikelihood(data,params,trapezoid_lpdf,trapezoid_lcdf, &
            trapezoid_lsf,trapezoid_interval)
    end function trapezoid_censored_loglikelihood

    pure function trapezoid_censored_nnlf(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite plateau start in [0,1]
        real(dp), intent(in) :: d !! finite plateau end in [c,1]
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -trapezoid_censored_loglikelihood(data,c, d,loc,scale)
    end function trapezoid_censored_nnlf


    pure subroutine normal_transform_tails(logt, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: logt !! log positive half-normal coordinate, possibly infinite
        real(dp), intent(out) :: logcdf !! log half-normal CDF
        real(dp), intent(out) :: logsf !! log half-normal survival probability
        real(dp) :: t
        if (logt <= log(0.01_dp)) then
            logcdf = scifort_log_two+normal_small_interval(0.0_dp,logt)
            logsf = log1p_safe(-exp(logcdf))
        else if (logt <= 0.0_dp) then
            t = exp(logt)
            logcdf = log(erf(t/sqrt(2.0_dp)))
            logsf = log1p_safe(-exp(logcdf))
        else
            logsf = negative_infinity(0.0_dp)
            if (logt <= 0.5_dp*(log(huge(1.0_dp))+scifort_log_two)) then
                logsf = scifort_log_two+normal_logsf(exp(logt))
            end if
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine normal_transform_tails

    pure function normal_transform_exponent(logt) result(value)
        implicit none
        real(dp), intent(in) :: logt !! log positive half-normal coordinate, possibly infinite
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (logt > 0.5_dp*(log(huge(1.0_dp))+scifort_log_two)) return
        value = 0.0_dp
        if (logt < log(tiny(1.0_dp))) return
        value = negative_hazard(2.0_dp*logt-scifort_log_two)
    end function normal_transform_exponent

    pure function normal_transform_interval(loglo, loghi, logdelta) result(value)
        implicit none
        real(dp), intent(in) :: loglo !! log smaller half-normal coordinate
        real(dp), intent(in) :: loghi !! log larger half-normal coordinate
        real(dp), intent(in) :: logdelta !! log positive difference of the two log coordinates
        real(dp) :: value, t, logwidth
        value = negative_infinity(0.0_dp)
        if (loglo > 0.5_dp*(log(huge(1.0_dp))+scifort_log_two)) return
        t = exp(loglo)
        logwidth = loghi+log_hazard_mass(logdelta)
        ! NaN requests a tail-difference fallback when the interval is not narrow.
        value = normal_small_interval(t,logwidth)
        if (.not. ieee_is_nan(value)) value = value+scifort_log_two
    end function normal_transform_interval

    pure function moyal_logt(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, z
        z = location_scale_z(x,params)
        value = -0.5_dp*z
        if (.not. ieee_is_finite(z)) then
            ! The half coordinate can be finite when the full coordinate overflows.
            ! Such overflow requires scale < huge()/2, so doubling scale is safe.
            value = -location_scale_z(x,[params(1),2.0_dp*params(2)])
        end if
    end function moyal_logt

    pure function moyal_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, logt, exponent
        logt = moyal_logt(x,params)
        exponent = normal_transform_exponent(logt)
        value = negative_infinity(0.0_dp)
        if (.not. ieee_is_finite(exponent)) return
        value = logt+exponent-scifort_log_sqrt_two_pi-log(params(2))
    end function moyal_lpdf

    pure function moyal_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, unused
        call normal_transform_tails(moyal_logt(x,params),unused,value)
    end function moyal_lcdf

    pure function moyal_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, unused
        call normal_transform_tails(moyal_logt(x,params),value,unused)
    end function moyal_lsf

    pure function moyal_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite location and positive scale
        real(dp) :: value, omitted, logdelta
        omitted = exp(moyal_lcdf(a,params))+exp(moyal_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        logdelta = log_distance(b,a)-log(params(2))-scifort_log_two
        value = normal_transform_interval(moyal_logt(b,params),moyal_logt(a,params),logdelta)
        if (ieee_is_nan(value)) value = tail_interval(a,b,params,moyal_lcdf,moyal_lsf)
    end function moyal_interval

    pure function levy_logt(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = 0.5_dp*(log(params(2))-log_distance(x,params(1)))
    end function levy_logt

    pure function levy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, logt
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        logt = levy_logt(x,params)
        value = logt-log_distance(x,params(1))-scifort_log_sqrt_two_pi &
            +normal_transform_exponent(logt)
    end function levy_lpdf

    pure function levy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, unused
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        call normal_transform_tails(levy_logt(x,params),unused,value)
    end function levy_lcdf

    pure function levy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, unused
        value = 0.0_dp
        if (x <= params(1)) return
        call normal_transform_tails(levy_logt(x,params),value,unused)
    end function levy_lsf

    pure function levy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, omitted, logdelta
        value = levy_lcdf(b,params)
        if (a <= params(1)) return
        omitted = exp(levy_lcdf(a,params))+exp(levy_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        logdelta = burr_logsoftplus(log_distance(b,a)-log_distance(a,params(1))) &
            -scifort_log_two
        value = normal_transform_interval(levy_logt(b,params),levy_logt(a,params),logdelta)
        if (ieee_is_nan(value)) value = tail_interval(a,b,params,levy_lcdf,levy_lsf)
    end function levy_interval

    pure function levy_l_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_lpdf(-x,[-params(1),params(2)])
    end function levy_l_lpdf

    pure function levy_l_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_lsf(-x,[-params(1),params(2)])
    end function levy_l_lcdf

    pure function levy_l_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_lcdf(-x,[-params(1),params(2)])
    end function levy_l_lsf

    pure function levy_l_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite upper support endpoint and positive scale
        real(dp) :: value
        value = levy_interval(-b,-a,[-params(1),params(2)])
    end function levy_l_interval

    pure function moyal_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, moyal_lpdf, moyal_lcdf, &
            moyal_lsf, moyal_interval)
    end function moyal_censored_loglikelihood

    pure function moyal_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -moyal_censored_loglikelihood(data, loc, scale)
    end function moyal_censored_nnlf

    pure function levy_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, levy_lpdf, levy_lcdf, &
            levy_lsf, levy_interval)
    end function levy_censored_loglikelihood

    pure function levy_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -levy_censored_loglikelihood(data, loc, scale)
    end function levy_censored_nnlf

    pure function levy_l_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite upper support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, levy_l_lpdf, levy_l_lcdf, &
            levy_l_lsf, levy_l_interval)
    end function levy_l_censored_loglikelihood

    pure function levy_l_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite upper support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -levy_l_censored_loglikelihood(data, loc, scale)
    end function levy_l_censored_nnlf


    pure function hazard_complement(logh, nh) result(value)
        implicit none
        real(dp), intent(in) :: logh !! log positive hazard, possibly infinite
        real(dp), intent(in) :: nh !! negative hazard, possibly negative infinity
        real(dp) :: value
        value = log_hazard_mass(logh)
        if (nh <= -scifort_log_two) value = log1p_safe(-exp(nh))
    end function hazard_complement

    pure subroutine genlogistic_hazard(x, params, z, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp), intent(out) :: z !! standardized coordinate, possibly infinite
        real(dp), intent(out) :: logh !! log(c*softplus(-z)), possibly infinite
        real(dp), intent(out) :: nh !! negative c*softplus(-z), possibly infinite
        real(dp) :: base
        z = location_scale_z(x,params(2:3))
        logh = log(params(1))+burr_logsoftplus(-z)
        if (.not. ieee_is_finite(z)) then
            if (z < 0.0_dp) then
                logh = log(params(1))+log_distance(params(2),x)-log(params(3))
            end if
        end if
        nh = negative_hazard(logh)
        if (.not. ieee_is_finite(z)) return
        base = log1pexp(-z)
        if (base < tiny(1.0_dp)) return
        if (base > 1.0_dp) then
            if (params(1) > huge(1.0_dp)/base) return
        end if
        nh = -params(1)*base
    end subroutine genlogistic_hazard

    pure function genlogistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logh, nh
        call genlogistic_hazard(x,params,z,logh,nh)
        value = log(params(1))-log(params(3))-log1pexp(z)+nh
    end function genlogistic_lpdf

    pure function genlogistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logh
        call genlogistic_hazard(x,params,z,logh,value)
    end function genlogistic_lcdf

    pure function genlogistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, z, logh, nh
        call genlogistic_hazard(x,params,z,logh,nh)
        value = hazard_complement(logh,nh)
    end function genlogistic_lsf

    pure function genlogistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, za, zb, logh, anchor, nh, omitted, ld, delta, lu, diff
        omitted = exp(genlogistic_lcdf(a,params))+exp(genlogistic_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call genlogistic_hazard(a,params,za,logh,nh)
        call genlogistic_hazard(b,params,zb,logh,anchor)
        value = anchor
        if (.not. ieee_is_finite(anchor)) return
        ld = log_distance(b,a)-log(params(3))
        ! Difference softplus(-za)-softplus(-zb), retaining the original width.
        if (ld < log(0.5_dp)) then
            lu = ld
            if (ld >= log(epsilon(1.0_dp))-2.0_dp) lu = log(expm1_safe(exp(ld)))
            diff = burr_logsoftplus(-log1pexp(zb)+lu)
        else if (zb <= 0.0_dp) then
            diff = ld
            if (ld < log(huge(1.0_dp))) then
                delta = exp(ld)
                diff = log(delta+log1pexp(za)-log1pexp(zb))
            end if
        else if (za >= 0.0_dp) then
            lu = -log1pexp(za)+log_hazard_mass(ld)
            diff = lu
            if (lu >= log(epsilon(1.0_dp))-2.0_dp) diff = log(-log1p_safe(-exp(lu)))
        else if (ieee_is_finite(za)) then
            diff = log(log1pexp(-za)-log1pexp(-zb))
        else
            diff = log_distance(params(2),a)-log(params(3))
        end if
        value = anchor+log_hazard_mass(log(params(1))+diff)
    end function genlogistic_interval

    pure function log_expm1_from_log(logz) result(value)
        implicit none
        real(dp), intent(in) :: logz !! log positive standardized distance, possibly infinite
        real(dp) :: value, z
        value = logz
        if (logz < log(epsilon(1.0_dp))-2.0_dp) return
        value = positive_infinity(0.0_dp)
        if (logz > log(huge(1.0_dp))) return
        z = exp(logz)
        value = z+log1p_safe(-exp(-z))
        if (z < 0.5_dp) value = log(expm1_safe(z))
    end function log_expm1_from_log

    pure subroutine gompertz_hazard(x, params, z, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp), intent(out) :: z !! nonnegative standardized coordinate, possibly infinite
        real(dp), intent(out) :: logh !! log(c*expm1(z)), possibly infinite
        real(dp), intent(out) :: nh !! negative c*expm1(z), possibly negative infinity
        real(dp) :: logz, base
        z = location_scale_z(x,params(2:3))
        logz = log_distance(x,params(2))-log(params(3))
        logh = log(params(1))+log_expm1_from_log(logz)
        ! Avoid a log/exp round trip for representable z, especially near unit mass.
        if (z >= 0.5_dp) logh = log(params(1))+z+log1p_safe(-exp(-z))
        nh = negative_hazard(logh)
        if (z <= 0.0_dp .or. z > log(huge(1.0_dp))) return
        base = expm1_safe(z)
        if (base < tiny(1.0_dp)) return
        if (base > 1.0_dp) then
            if (params(1) > huge(1.0_dp)/base) return
        end if
        nh = -params(1)*base
    end subroutine gompertz_hazard

    pure function gompertz_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh, nh
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            value = log(params(1))-log(params(3))
            return
        end if
        call gompertz_hazard(x,params,z,logh,nh)
        if (.not. ieee_is_finite(nh)) return
        value = log(params(1))-log(params(3))+z+nh
    end function gompertz_lpdf

    pure function gompertz_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh, nh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        call gompertz_hazard(x,params,z,logh,nh)
        value = hazard_complement(logh,nh)
    end function gompertz_lcdf

    pure function gompertz_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh
        value = 0.0_dp
        if (x <= params(2)) return
        call gompertz_hazard(x,params,z,logh,value)
    end function gompertz_lsf

    pure function gompertz_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z, logh, anchor, omitted, logwidth
        value = gompertz_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(gompertz_lcdf(a,params))+exp(gompertz_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call gompertz_hazard(a,params,z,logh,anchor)
        value = anchor
        if (.not. ieee_is_finite(anchor)) return
        logwidth = log_distance(b,a)-log(params(3))
        value = anchor+log_hazard_mass(log(params(1))+z+log_expm1_from_log(logwidth))
    end function gompertz_interval

    pure function invweibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, logh, nh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logh = -weibull_log_hazard(x,params)
        if (.not. ieee_is_finite(logh)) return
        nh = negative_hazard(logh)
        value = log(params(1))-log_distance(x,params(2))+logh+nh
    end function invweibull_lpdf

    pure function invweibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        value = negative_hazard(-weibull_log_hazard(x,params))
    end function invweibull_lcdf

    pure function invweibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, logh
        value = 0.0_dp
        if (x <= params(2)) return
        logh = -weibull_log_hazard(x,params)
        value = hazard_complement(logh,negative_hazard(logh))
    end function invweibull_lsf

    pure function invweibull_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, omitted, anchor, logdelta, logh
        value = invweibull_lcdf(b,params)
        if (a <= params(2)) return
        omitted = exp(invweibull_lcdf(a,params))+exp(invweibull_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        anchor = value
        if (.not. ieee_is_finite(anchor)) return
        logdelta = log(params(1))+burr_logsoftplus( &
            log_distance(b,a)-log_distance(a,params(2)))
        logh = -weibull_log_hazard(a,params)
        value = anchor+log_hazard_mass(logh+log_hazard_mass(logdelta))
    end function invweibull_interval

    pure function genlogistic_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,genlogistic_lpdf,genlogistic_lcdf, &
            genlogistic_lsf,genlogistic_interval)
    end function genlogistic_censored_loglikelihood

    pure function genlogistic_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -genlogistic_censored_loglikelihood(data,c,loc,scale)
    end function genlogistic_censored_nnlf

    pure function gompertz_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,gompertz_lpdf,gompertz_lcdf, &
            gompertz_lsf,gompertz_interval)
    end function gompertz_censored_loglikelihood

    pure function gompertz_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -gompertz_censored_loglikelihood(data,c,loc,scale)
    end function gompertz_censored_nnlf

    pure function invweibull_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c,0.0_dp,1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,invweibull_lpdf,invweibull_lcdf, &
            invweibull_lsf,invweibull_interval)
    end function invweibull_censored_loglikelihood

    pure function invweibull_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -invweibull_censored_loglikelihood(data,c,loc,scale)
    end function invweibull_censored_nnlf


    pure function uniform_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = powerlaw_censored_loglikelihood(data,1.0_dp,loc,scale)
    end function uniform_censored_loglikelihood

    pure function uniform_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -uniform_censored_loglikelihood(data,loc,scale)
    end function uniform_censored_nnlf

    pure function arcsine_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [0.5_dp,0.5_dp,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. beta_parameters_valid(params)) return
        value = censored_loglikelihood(data,params,beta_lpdf,beta_lcdf,beta_lsf,arcsine_interval)
    end function arcsine_censored_loglikelihood

    pure function arcsine_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -arcsine_censored_loglikelihood(data,loc,scale)
    end function arcsine_censored_nnlf

    pure function arcsine_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! half beta shapes, finite endpoint, positive width
        real(dp) :: value, omitted
        omitted = exp(beta_lcdf(a,params))+exp(beta_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
        else
            value = beta_interval(a,b,params)
        end if
    end function arcsine_interval

    pure function powerlaw_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive power-function shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value, params(4)
        params = [a,1.0_dp,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. beta_parameters_valid(params)) return
        value = censored_loglikelihood(data,params,powerlaw_lpdf,powerlaw_lcdf, &
            powerlaw_lsf,powerlaw_interval)
    end function powerlaw_censored_loglikelihood

    pure function powerlaw_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive power-function shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive support width, default one
        real(dp) :: value
        value = -powerlaw_censored_loglikelihood(data,a,loc,scale)
    end function powerlaw_censored_nnlf

    pure function powerlaw_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, z, w, logz, logw, factor
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -log(params(4))
            return
        end if
        if (right == 0.0_dp .or. params(1) == 1.0_dp) then
            value = log(params(1))-log(params(4))
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        factor = params(1)-1.0_dp
        if (factor > 1.0_dp) then
            if (-logz > huge(1.0_dp)/factor) return
        end if
        value = log(params(1))-log(params(4))+factor*logz
    end function powerlaw_lpdf

    pure subroutine powerlaw_hazard(x, params, nh, logh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly inside support
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp), intent(out) :: nh !! log CDF, possibly negative infinity
        real(dp), intent(out) :: logh !! log of negative log CDF
        real(dp) :: left, right, z, w, logz, logw
        call beta_distances(x,params,left,right)
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        if (w < tiny(1.0_dp)) then
            ! -log(z)=w+O(w*w); retain log(w) when division rounds to zero
            ! or to a subnormal with poor relative precision.
            logh = log(params(1))+logw
            nh = negative_hazard(logh)
            return
        end if
        logh = log(params(1))+log(-logz)
        nh = negative_infinity(0.0_dp)
        if (params(1) > 1.0_dp) then
            if (-logz > huge(1.0_dp)/params(1)) return
        end if
        nh = params(1)*logz
    end subroutine powerlaw_hazard

    pure function powerlaw_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, logh
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        value = 0.0_dp
        if (right <= 0.0_dp) return
        call powerlaw_hazard(x,params,value,logh)
    end function powerlaw_lcdf

    pure function powerlaw_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, nh, logh
        value = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        value = negative_infinity(0.0_dp)
        if (right <= 0.0_dp) return
        call powerlaw_hazard(x,params,nh,logh)
        if (nh <= -scifort_log_two) then
            value = log1p_safe(-exp(nh))
        else
            value = log_hazard_mass(logh)
        end if
    end function powerlaw_lsf

    pure function powerlaw_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, unit beta shape, location, width
        real(dp) :: value, left, right, logr, logdelta, omitted
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = powerlaw_lcdf(b,params)
            return
        end if
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = powerlaw_lsf(a,params)
            return
        end if
        omitted = exp(powerlaw_lcdf(a,params))+exp(powerlaw_lsf(b,params))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
        else if (params(1) == 1.0_dp) then
            value = log_distance(b,a)-log(params(4))
        else
            logr = log_distance(b,a)-log_distance(a,params(3))
            logdelta = log(params(1))+burr_logsoftplus(logr)
            value = powerlaw_lcdf(b,params)+log_hazard_mass(logdelta)
        end if
    end function powerlaw_interval


    pure function fisk_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = burr12_censored_loglikelihood(data,c,1.0_dp,loc,scale)
    end function fisk_censored_loglikelihood

    pure function fisk_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -fisk_censored_loglikelihood(data,c,loc,scale)
    end function fisk_censored_nnlf

    pure function burr_censored_loglikelihood(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(4)
        params = [c,d,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. d <= 0.0_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,burr_lpdf,burr_lcdf, &
            burr_lsf,burr_interval)
    end function burr_censored_loglikelihood

    pure function burr_censored_nnlf(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -burr_censored_loglikelihood(data,c,d,loc,scale)
    end function burr_censored_nnlf

    pure function burr12_censored_loglikelihood(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(4)
        params = [c,d,0.0_dp,1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. d <= 0.0_dp .or. params(4) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,burr12_lpdf,burr12_lcdf, &
            burr12_lsf,burr12_interval)
    end function burr12_censored_loglikelihood

    pure function burr12_censored_nnlf(data, c, d, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive first shape
        real(dp), intent(in) :: d !! finite positive second shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -burr12_censored_loglikelihood(data,c,d,loc,scale)
    end function burr12_censored_nnlf

    pure function burr_logsoftplus(x) result(value)
        implicit none
        real(dp), intent(in) :: x !! log argument of log1p(exp(x)), possibly infinite
        real(dp) :: value
        ! For exp(x) below roundoff, log(log1p(exp(x))) = x + O(exp(x)).
        value = x
        if (x >= log(epsilon(1.0_dp))-2.0_dp) value = log(log1pexp(x))
    end function burr_logsoftplus

    pure subroutine burr_hazard(x, params, reverse, q, logh, nh)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        real(dp), intent(out) :: q !! signed c*log(z), possibly infinite
        real(dp), intent(out) :: logh !! log(d*log1p(exp(q))), possibly infinite
        real(dp), intent(out) :: nh !! negative d*log1p(exp(q)), possibly infinite
        real(dp) :: v, logbase, base, product
        v = lognormal_logz(x,[params(1),params(3),params(4)])
        q = weibull_log_hazard(x,[params(1),params(3),params(4)])
        if (reverse) q = -q
        if (ieee_is_finite(q)) then
            logbase = burr_logsoftplus(q)
        else if (q > 0.0_dp) then
            ! d can compensate c even when c*log(z) exceeds huge().
            logbase = log(params(1))+log(abs(v))
        else
            logbase = negative_infinity(0.0_dp)
        end if
        logh = log(params(2))+logbase
        nh = negative_hazard(logh)
        if (.not. ieee_is_finite(q)) then
            if (q > 0.0_dp) then
                ! Here c is large and abs(log(z)) > 1. Group c*d first.
                if (params(2) <= huge(1.0_dp)/params(1)) then
                    product = params(1)*params(2)
                    logh = log(product)+log(abs(v))
                    nh = negative_infinity(0.0_dp)
                    if (product <= huge(1.0_dp)/abs(v)) nh = -product*abs(v)
                end if
            end if
            return
        end if
        base = log1pexp(q)
        if (base < tiny(1.0_dp)) return
        if (base > 1.0_dp) then
            if (params(2) > huge(1.0_dp)/base) return
        end if
        ! Direct multiplication avoids amplifying log/exp rounding near unit mass.
        nh = -params(2)*base
    end subroutine burr_hazard

    pure function burr_density(x, params, reverse) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        real(dp) :: value, exponent, q, logh, nh
        value = negative_infinity(0.0_dp)
        if (x < params(3)) return
        if (x == params(3)) then
            exponent = params(1)
            if (reverse) then
                if (params(1) > 1.0_dp) then
                    if (params(2) > huge(1.0_dp)/params(1)) return
                end if
                exponent = params(1)*params(2)
            end if
            if (exponent < 1.0_dp) value = positive_infinity(0.0_dp)
            if (exponent == 1.0_dp) then
                value = log(params(1))+log(params(2))-log(params(4))
            end if
            return
        end if
        call burr_hazard(x,params,reverse,q,logh,nh)
        value = log(params(1))+log(params(2))-log_distance(x,params(3)) &
            -log1pexp(-q)+nh
    end function burr_density

    pure function burr_tail(x, params, reverse, complement) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        logical, intent(in) :: complement !! true to complement the direct exp(-hazard) tail
        real(dp) :: value, q, logh, nh
        value = 0.0_dp
        if (reverse .neqv. complement) value = negative_infinity(0.0_dp)
        if (x <= params(3)) return
        call burr_hazard(x,params,reverse,q,logh,nh)
        value = nh
        if (complement) then
            if (nh <= -scifort_log_two) then
                value = log1p_safe(-exp(nh))
            else
                value = log_hazard_mass(logh)
            end if
        end if
    end function burr_tail

    pure function burr_mass(a, b, params, reverse) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        logical, intent(in) :: reverse !! true for Type III, false for Type XII
        real(dp) :: value, qa, qb, qlo, qhi, la, lb, logdelta, logr, span, delta
        real(dp) :: logdiff, logu, omitted, anchor, unused, nh
        value = burr_tail(b,params,reverse,.not. reverse)
        if (a <= params(3)) return
        omitted = exp(burr_tail(a,params,reverse,.not. reverse)) &
            +exp(burr_tail(b,params,reverse,reverse))
        if (omitted <= 0.5_dp) then
            value = log1p_safe(-omitted)
            return
        end if
        call burr_hazard(a,params,reverse,qa,unused,nh)
        anchor = nh
        call burr_hazard(b,params,reverse,qb,unused,nh)
        if (reverse) anchor = nh
        qlo = min(qa,qb)
        qhi = max(qa,qb)
        logr = log_distance(b,a)-log_distance(a,params(3))
        logdelta = log(params(1))+burr_logsoftplus(logr)
        if (logdelta < log(0.5_dp)) then
            ! softplus(qhi)-softplus(qlo)=softplus(log(sigmoid(qlo)*expm1(delta))).
            logu = logdelta
            if (logdelta >= log(epsilon(1.0_dp))-2.0_dp) then
                delta = exp(logdelta)
                logu = log(expm1_safe(delta))
            end if
            logdiff = burr_logsoftplus(-log1pexp(-qlo)+logu)
        else if (qlo >= 0.0_dp) then
            logdiff = logdelta
            if (logdelta < log(huge(1.0_dp))) then
                delta = exp(logdelta)
                span = delta+log1pexp(-qhi)-log1pexp(-qlo)
                logdiff = log(span)
            end if
        else if (qhi <= 0.0_dp) then
            ! 1-exp(-difference)=sigmoid(qhi)*(1-exp(-delta)).
            logu = -log1pexp(-qhi)+log_hazard_mass(logdelta)
            logdiff = logu
            if (logu >= log(epsilon(1.0_dp))-2.0_dp) then
                logdiff = log(-log1p_safe(-exp(logu)))
            end if
        else if (ieee_is_finite(qhi)) then
            logdiff = log(log1pexp(qhi)-log1pexp(qlo))
        else
            ! Crossing interval with an unrepresentably large positive endpoint.
            la = abs(lognormal_logz(a,[params(1),params(3),params(4)]))
            lb = abs(lognormal_logz(b,[params(1),params(3),params(4)]))
            span = lb
            if (reverse) span = la
            logdiff = log(params(1))+log(span)
        end if
        value = anchor+log_hazard_mass(log(params(2))+logdiff)
    end function burr_mass

    pure function burr_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_density(x,params,.true.)
    end function burr_lpdf

    pure function burr_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.true.,.false.)
    end function burr_lcdf

    pure function burr_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.true.,.true.)
    end function burr_lsf

    pure function burr_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_mass(a,b,params,.true.)
    end function burr_interval

    pure function burr12_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_density(x,params,.false.)
    end function burr12_lpdf

    pure function burr12_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.false.,.true.)
    end function burr12_lcdf

    pure function burr12_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_tail(x,params,.false.,.false.)
    end function burr12_lsf

    pure function burr12_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, greater than a
        real(dp), intent(in) :: params(:) !! c, d, finite endpoint, positive scale
        real(dp) :: value
        value = burr_mass(a,b,params,.false.)
    end function burr12_interval


    pure function lomax_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,lomax_lpdf,lomax_lcdf, &
            lomax_lsf,lomax_interval)
    end function lomax_censored_loglikelihood

    pure function lomax_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -lomax_censored_loglikelihood(data,c,loc,scale)
    end function lomax_censored_nnlf

    pure function lomax_base(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z
        z = location_scale_z(x,params(2:3))
        if (ieee_is_finite(z)) then
            value = log1p_safe(z)
        else
            value = log1pexp(lognormal_logz(x,params))
        end if
    end function lomax_base

    pure function lomax_negative_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, base
        base = lomax_base(x,params)
        if (base < tiny(1.0_dp)) then
            ! Subnormal displacements lose relative precision; shape can amplify it.
            value = negative_hazard(log(params(1))+lognormal_logz(x,params))
        else
            value = negative_infinity(0.0_dp)
            if (base > 1.0_dp) then
                if (params(1) > huge(1.0_dp)/base) return
            end if
            value = -params(1)*base
        end if
    end function lomax_negative_hazard

    pure function lomax_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        value = log(params(1))-log(params(3))
        if (x == params(2)) return
        value = value-lomax_base(x,params)+lomax_negative_hazard(x,params)
    end function lomax_lpdf

    pure function lomax_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, h, base, logh
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        h = lomax_negative_hazard(x,params)
        if (h <= -scifort_log_two) then
            value = log1p_safe(-exp(h))
        else
            base = lomax_base(x,params)
            if (base < tiny(1.0_dp)) then
                logh = log(params(1))+lognormal_logz(x,params)
            else
                logh = log(params(1))+log(base)
            end if
            value = log_hazard_mass(logh)
        end if
    end function lomax_lcdf

    pure function lomax_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(2)) return
        value = lomax_negative_hazard(x,params)
    end function lomax_lsf

    pure function lomax_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, logs, logd, logdenom, logr, logspan
        value = lomax_lcdf(b,params)
        if (a <= params(2)) return
        ! Denominator scale+a-loc is a positive sum which may exceed huge().
        logs = log(params(3))
        logd = log_distance(a,params(2))
        logdenom = max(logs,logd)+log1pexp(-abs(logs-logd))
        logr = log_distance(b,a)-logdenom
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        value = lomax_lsf(a,params)+log_hazard_mass(log(params(1))+logspan)
    end function lomax_interval

    pure function pareto_censored_loglikelihood(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location shift, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [b, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (b <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,pareto_lpdf,pareto_lcdf, &
            pareto_lsf,pareto_interval)
    end function pareto_censored_loglikelihood

    pure function pareto_censored_nnlf(data, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: b !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite location shift, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -pareto_censored_loglikelihood(data,b,loc,scale)
    end function pareto_censored_nnlf

    pure function pareto_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logz = lognormal_logz(x,params)
        if (logz < 0.0_dp) return
        value = log(params(1))-log_distance(x,params(2))-weibull_log_hazard(x,params)
    end function pareto_lpdf

    pure function pareto_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logz = lognormal_logz(x,params)
        if (logz <= 0.0_dp) return
        ! Keep small hazards in log form, but avoid a log/exp round trip for
        ! large hazards: its error would be amplified in a nearly unit CDF.
        value = weibull_log_hazard(x,params)
        if (value >= scifort_log_two) then
            value = log1p_safe(-exp(-value))
        else
            value = log_hazard_mass(log(params(1))+log(logz))
        end if
    end function pareto_lcdf

    pure function pareto_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(2)) return
        value = -max(0.0_dp,weibull_log_hazard(x,params))
    end function pareto_lsf

    pure function pareto_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite location, positive scale
        real(dp) :: value, logr, logspan
        ! Clip intervals crossing the support boundary without forming loc+scale.
        value = pareto_lcdf(b,params)
        if (a <= params(2)) return
        if (lognormal_logz(a,params) <= 0.0_dp) return
        logr = log_distance(b,a)-log_distance(a,params(2))
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        value = pareto_lsf(a,params)+log_hazard_mass(log(params(1))+logspan)
    end function pareto_interval

    pure function loglaplace_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value, params(3)
        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (c <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data,params,loglaplace_lpdf,loglaplace_lcdf, &
            loglaplace_lsf,loglaplace_interval)
    end function loglaplace_censored_loglikelihood

    pure function loglaplace_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite positive shape
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! finite positive scale, default one
        real(dp) :: value
        value = -loglaplace_censored_loglikelihood(data,c,loc,scale)
    end function loglaplace_censored_nnlf

    pure function loglaplace_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            if (params(1) < 1.0_dp) value = positive_infinity(0.0_dp)
            if (params(1) == 1.0_dp) value = -scifort_log_two-log(params(3))
            return
        end if
        value = log(params(1))-scifort_log_two-log_distance(x,params(2)) &
            -abs(weibull_log_hazard(x,params))
    end function loglaplace_lpdf

    pure function loglaplace_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite left censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        z = weibull_log_hazard(x,params)
        if (z <= 0.0_dp) then
            value = z-scifort_log_two
        else
            value = log1p_safe(-0.5_dp*exp(-z))
        end if
    end function loglaplace_lcdf

    pure function loglaplace_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite right censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, z
        value = 0.0_dp
        if (x <= params(2)) return
        z = weibull_log_hazard(x,params)
        if (z >= 0.0_dp) then
            value = -z-scifort_log_two
        else
            value = log1p_safe(-0.5_dp*exp(z))
        end if
    end function loglaplace_lsf

    pure function loglaplace_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite endpoint, positive scale
        real(dp) :: value, za, zb, logr, logspan, omitted, left, right, large, small
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = loglaplace_lcdf(b,params)
            return
        end if
        za = lognormal_logz(a,params)
        zb = lognormal_logz(b,params)
        ! Work with the log of c*log((b-loc)/(a-loc)); retain adjacent widths.
        logr = log_distance(b,a)-log_distance(a,params(2))
        logspan = logr
        if (logr >= log(epsilon(1.0_dp))-2.0_dp) logspan = log(log1pexp(logr))
        if (zb <= 0.0_dp) then
            value = loglaplace_lcdf(b,params)+log_hazard_mass(log(params(1))+logspan)
        else if (za >= 0.0_dp) then
            value = loglaplace_lsf(a,params)+log_hazard_mass(log(params(1))+logspan)
        else
            omitted = exp(loglaplace_lcdf(a,params))+exp(loglaplace_lsf(b,params))
            if (omitted <= 0.5_dp) then
                value = log1p_safe(-omitted)
            else
                left = -scifort_log_two+log_hazard_mass(log(params(1))+log(-za))
                right = -scifort_log_two+log_hazard_mass(log(params(1))+log(zb))
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function loglaplace_interval


    pure function laplace_asymmetric_censored_loglikelihood(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! finite positive asymmetry parameter
        real(dp), intent(in), optional :: loc !! finite cusp location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(3)
        params = [kappa, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (kappa <= 0.0_dp .or. params(3) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, laplace_asymmetric_lpdf, &
            laplace_asymmetric_lcdf, laplace_asymmetric_lsf, laplace_asymmetric_interval)
    end function laplace_asymmetric_censored_loglikelihood

    pure function laplace_asymmetric_censored_nnlf(data, kappa, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: kappa !! finite positive asymmetry parameter
        real(dp), intent(in), optional :: loc !! finite cusp location, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -laplace_asymmetric_censored_loglikelihood(data, kappa, loc, scale)
    end function laplace_asymmetric_censored_nnlf

    pure function laplace_asymmetric_log_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point on either side of the cusp
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logs
        value = negative_infinity(0.0_dp)
        if (x == params(2)) return
        t = log(params(1))
        logs = log(params(3))
        value = log_distance(max(x,params(2)),min(x,params(2)))
        if (x > params(2)) then
            value = value+(t-logs)
        else
            value = value-(t+logs)
        end if
    end function laplace_asymmetric_log_hazard

    pure function laplace_asymmetric_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t
        t = abs(log(params(1)))
        value = negative_hazard(laplace_asymmetric_log_hazard(x,params)) &
            -t-log1p_safe(exp(-t)**2)-log(params(3))
    end function laplace_asymmetric_lpdf

    pure function laplace_asymmetric_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logh, left, right, large, small
        t = log(params(1))
        logh = laplace_asymmetric_log_hazard(x,params)
        left = -log1pexp(-2.0_dp*t)
        right = -log1pexp(2.0_dp*t)
        if (x <= params(2)) then
            value = left+negative_hazard(logh)
        else
            value = right+negative_hazard(logh)
            if (value <= -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                ! F = left weight + right weight*(1-exp(-H)); retain underflowed F.
                right = right+log_hazard_mass(logh)
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function laplace_asymmetric_lcdf

    pure function laplace_asymmetric_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logh, left, right, large, small
        t = log(params(1))
        logh = laplace_asymmetric_log_hazard(x,params)
        left = -log1pexp(-2.0_dp*t)
        right = -log1pexp(2.0_dp*t)
        if (x >= params(2)) then
            value = right+negative_hazard(logh)
        else
            value = left+negative_hazard(logh)
            if (value <= -scifort_log_two) then
                value = log1p_safe(-exp(value))
            else
                ! S = right weight + left weight*(1-exp(-H)); retain underflowed S.
                left = left+log_hazard_mass(logh)
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function laplace_asymmetric_lsf

    pure function laplace_asymmetric_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! positive shape, finite cusp, positive scale
        real(dp) :: value, t, logs, left, right, large, small, omitted
        t = log(params(1))
        logs = log(params(3))
        if (b <= params(2)) then
            value = laplace_asymmetric_lcdf(b,params) &
                + log_hazard_mass(log_distance(b,a)-(logs+t))
        else if (a >= params(2)) then
            value = laplace_asymmetric_lsf(a,params) &
                + log_hazard_mass(log_distance(b,a)+(t-logs))
        else
            omitted = exp(laplace_asymmetric_lcdf(a,params)) &
                + exp(laplace_asymmetric_lsf(b,params))
            if (omitted <= 0.5_dp) then
                value = log1p_safe(-omitted)
            else
                ! Add positive masses from each side rather than subtracting near-unit tails.
                left = -log1pexp(-2.0_dp*t) &
                    + log_hazard_mass(log_distance(params(2),a)-(logs+t))
                right = -log1pexp(2.0_dp*t) &
                    + log_hazard_mass(log_distance(b,params(2))+(t-logs))
                large = max(left,right)
                small = min(left,right)
                value = large+log1p_safe(exp(small-large))
            end if
        end if
    end function laplace_asymmetric_interval

    pure function hypsecant_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, hypsecant_lpdf, hypsecant_lcdf, &
            hypsecant_lsf, hypsecant_interval)
    end function hypsecant_censored_loglikelihood

    pure function hypsecant_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite center, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -hypsecant_censored_loglikelihood(data, loc, scale)
    end function hypsecant_censored_nnlf

    pure function hypsecant_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z, t
        z = abs(location_scale_z(x,params))
        t = exp(-z)
        value = -z-log1p_safe(t*t)+(scifort_log_two-log(scifort_pi))-log(params(2))
    end function hypsecant_lpdf

    pure function hypsecant_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x,params)
        value = (scifort_log_two-log(scifort_pi))+cauchy_logatan_exp(-abs(z))
        if (z > 0.0_dp) value = log1p_safe(-exp(value))
    end function hypsecant_lcdf

    pure function hypsecant_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x,params)
        value = (scifort_log_two-log(scifort_pi))+cauchy_logatan_exp(-abs(z))
        if (z < 0.0_dp) value = log1p_safe(-exp(value))
    end function hypsecant_lsf

    pure function hypsecant_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, za, zb, logmass, logratio, omitted, large, small
        value = negative_infinity(0.0_dp)
        za = location_scale_z(a,params)
        zb = location_scale_z(b,params)
        logmass = log_exponential_mass(log_distance(b,a)-log(params(2)))
        if (a >= params(1) .or. b <= params(1)) then
            ! Reflect a negative interval; 0 <= za <= zb thereafter.
            if (b <= params(1)) then
                large = -za
                za = -zb
                zb = large
            end if
            if (.not. ieee_is_finite(za)) return
            ! atan(exp(-za))-atan(exp(-zb)) has a positive difference formula.
            logratio = -za+logmass-log1p_safe(exp(-za)*exp(-zb))
        else
            omitted = exp(hypsecant_lcdf(a,params))+exp(hypsecant_lsf(b,params))
            if (omitted <= 0.5_dp) then
                value = log1p_safe(-omitted)
                return
            end if
            ! tan(angle)=(1-exp(-(zb-za)))/(exp(za)+exp(-zb)).
            large = max(za,-zb)
            small = min(za,-zb)
            logratio = logmass-(large+log1p_safe(exp(small-large)))
        end if
        value = (scifort_log_two-log(scifort_pi))+cauchy_logatan_exp(logratio)
    end function hypsecant_interval

    pure function halfnorm_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite normal scale, default one
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, halfnorm_lpdf, halfnorm_lcdf, &
            halfnorm_lsf, halfnorm_interval)
    end function halfnorm_censored_loglikelihood

    pure function halfnorm_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite normal scale, default one
        real(dp) :: value
        value = -halfnorm_censored_loglikelihood(data, loc, scale)
    end function halfnorm_censored_nnlf

    pure function halfnorm_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, z
        value = negative_infinity(0.0_dp)
        if (x < params(1)) return
        z = location_scale_z(x, params)
        if (z > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        value = (-0.5_dp*z)*z + (scifort_log_two-scifort_log_sqrt_two_pi) &
            - log(params(2))
    end function halfnorm_lpdf

    pure function halfnorm_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, other
        call halfnorm_tails(x, params, value, other)
    end function halfnorm_lcdf

    pure function halfnorm_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, other
        call halfnorm_tails(x, params, other, value)
    end function halfnorm_lsf

    pure subroutine halfnorm_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp), intent(out) :: logcdf !! log lower-tail probability, -infinity at/below support
        real(dp), intent(out) :: logsf !! log upper-tail probability, zero at/below support
        real(dp) :: logz, z
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        if (x <= params(1)) return
        logz = log_distance(x, params(1))-log(params(2))
        if (logz <= 0.0_dp) then
            if (logz <= log(0.01_dp)) then
                ! F(x)=2*P(loc < N(loc,scale) <= x); retain log widths below tiny.
                logcdf = scifort_log_two + normal_small_interval(0.0_dp, logz)
            else
                z = location_scale_z(x, params)
                logcdf = log(erf(z/sqrt(2.0_dp)))
            end if
            logsf = log1p_safe(-exp(logcdf))
        else
            z = location_scale_z(x, params)
            logsf = scifort_log_two + normal_logsf(z)
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine halfnorm_tails

    pure function halfnorm_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, greater than a
        real(dp), intent(in) :: params(:) !! finite support endpoint and positive finite scale
        real(dp) :: value, lower_sf, upper_sf, z, logscaled
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = halfnorm_lcdf(b, params)
            return
        end if
        lower_sf = halfnorm_lsf(a, params)
        if (.not. ieee_is_finite(lower_sf)) return
        upper_sf = halfnorm_lsf(b, params)
        if (upper_sf-lower_sf <= -scifort_log_two) then
            ! Preserve tiny excluded probabilities for nearly full-support intervals.
            value = lower_sf + log1p_safe(-exp(upper_sf-lower_sf))
        else
            z = location_scale_z(a, params)
            logscaled = log_distance(b,a)-log(params(2))
            value = normal_small_interval(z,logscaled)
            if (ieee_is_nan(value)) then
                value = tail_interval(a,b,params,halfnorm_lcdf,halfnorm_lsf)
            else
                value = value+scifort_log_two
            end if
        end if
    end function halfnorm_interval

    pure function gibrat_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = lognormal_censored_loglikelihood(data, 1.0_dp, loc, scale)
    end function gibrat_censored_loglikelihood

    pure function gibrat_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -gibrat_censored_loglikelihood(data, loc, scale)
    end function gibrat_censored_nnlf

    ! Erlang is gamma, including the existing positive noninteger-shape convention.
    pure function erlang_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape, conventionally an integer
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = gamma_censored_loglikelihood(data, a, loc, scale)
    end function erlang_censored_loglikelihood

    pure function erlang_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite positive shape, conventionally an integer
        real(dp), intent(in), optional :: loc !! finite lower support endpoint, default zero
        real(dp), intent(in), optional :: scale !! positive finite scale, default one
        real(dp) :: value
        value = -erlang_censored_loglikelihood(data, a, loc, scale)
    end function erlang_censored_nnlf

    pure function rayleigh_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, rayleigh_lpdf, rayleigh_lcdf, &
            rayleigh_lsf, rayleigh_interval)
    end function rayleigh_censored_loglikelihood

    pure function rayleigh_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -rayleigh_censored_loglikelihood(data, loc, scale)
    end function rayleigh_censored_nnlf

    pure function weibull_max_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value, params(3)

        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(c)) return
        if (c <= 0.0_dp .or. .not. valid_loc_scale(params(2), params(3))) return
        value = censored_loglikelihood(data, params, weibull_max_lpdf, weibull_max_lcdf, &
            weibull_max_lsf, weibull_max_interval)
    end function weibull_max_censored_loglikelihood

    pure function weibull_max_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! upper support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value
        value = -weibull_max_censored_loglikelihood(data, c, loc, scale)
    end function weibull_max_censored_nnlf

    pure function negative_hazard(logh) result(value)
        implicit none
        real(dp), intent(in) :: logh !! log nonnegative cumulative hazard, possibly infinite
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (logh > log(huge(1.0_dp))) return
        value = -exp(logh)
    end function negative_hazard

    pure function log_hazard_mass(logh) result(value)
        implicit none
        real(dp), intent(in) :: logh !! log positive cumulative-hazard increment, possibly infinite
        real(dp) :: value
        if (logh <= log(scifort_log_two)) then
            value = log_exponential_mass(logh)
        else
            value = log1p_safe(-exp(negative_hazard(logh)))
        end if
    end function log_hazard_mass

    pure function rayleigh_log_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above the support endpoint
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = 2.0_dp*(log_distance(x,params(1))-log(params(2))) - scifort_log_two
    end function rayleigh_log_hazard

    pure function rayleigh_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        value = negative_hazard(rayleigh_log_hazard(x,params))
        value = value + log_distance(x,params(1)) - 2.0_dp*log(params(2))
    end function rayleigh_lpdf

    pure function rayleigh_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(1)) return
        value = log_hazard_mass(rayleigh_log_hazard(x,params))
    end function rayleigh_lcdf

    pure function rayleigh_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(1)) return
        value = negative_hazard(rayleigh_log_hazard(x,params))
    end function rayleigh_lsf

    pure function rayleigh_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive scale
        real(dp) :: value, loga, logb, logdelta
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = rayleigh_lcdf(b,params)
        else
            value = rayleigh_lsf(a,params)
            if (.not. ieee_is_finite(value)) return
            loga = log_distance(a,params(1))
            logb = log_distance(b,params(1))
            ! Delta H = (b-a)*((b-loc)+(a-loc))/(2*scale**2).
            logdelta = log_distance(b,a)+logb+log1p_safe(exp(loga-logb)) &
                -scifort_log_two-2.0_dp*log(params(2))
            value = value + log_hazard_mass(logdelta)
        end if
    end function rayleigh_interval

    pure function weibull_max_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value, logh
        value = negative_infinity(0.0_dp)
        if (x > params(2)) return
        if (x == params(2)) then
            value = weibull_logpdf(0.0_dp,params(1),scale=params(3))
            return
        end if
        logh = weibull_log_hazard(-x,[params(1),-params(2),params(3)])
        if (logh > log(huge(1.0_dp))) return
        value = log(params(1))-log_distance(params(2),x)+logh-exp(logh)
    end function weibull_max_lpdf

    pure function weibull_max_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value
        value = 0.0_dp
        if (x >= params(2)) return
        value = negative_hazard(weibull_log_hazard(-x,[params(1),-params(2),params(3)]))
    end function weibull_max_lcdf

    pure function weibull_max_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x >= params(2)) return
        value = log_hazard_mass(weibull_log_hazard(-x,[params(1),-params(2),params(3)]))
    end function weibull_max_lsf

    pure function weibull_max_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, finite upper support endpoint, scale
        real(dp) :: value, lower_cdf, upper_cdf
        value = negative_infinity(0.0_dp)
        if (a >= params(2)) return
        if (b >= params(2)) then
            value = weibull_max_lsf(a,params)
        else
            upper_cdf = weibull_max_lcdf(b,params)
            if (.not. ieee_is_finite(upper_cdf)) return
            lower_cdf = weibull_max_lcdf(a,params)
            if (lower_cdf-upper_cdf <= -scifort_log_two) then
                value = upper_cdf+log1p_safe(-exp(lower_cdf-upper_cdf))
            else
                value = weibull_interval(-b,-a,[params(1),-params(2),params(3)])
            end if
        end if
    end function weibull_max_interval


    pure function halfcauchy_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, halfcauchy_lpdf, halfcauchy_lcdf, &
            halfcauchy_lsf, halfcauchy_interval)
    end function halfcauchy_censored_loglikelihood

    pure function halfcauchy_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -halfcauchy_censored_loglikelihood(data, loc, scale)
    end function halfcauchy_censored_nnlf

    pure function halfcauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(1)) return
        value = scifort_log_two + cauchy_lpdf(x, params)
    end function halfcauchy_lpdf

    pure function halfcauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halfcauchy_tails(x, params, value, other)
    end function halfcauchy_lcdf

    pure function halfcauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halfcauchy_tails(x, params, other, value)
    end function halfcauchy_lsf

    pure function halfcauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, lower_sf, upper_sf
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = halfcauchy_lcdf(b, params)
        else
            lower_sf = halfcauchy_lsf(a, params)
            if (.not. ieee_is_finite(lower_sf)) return
            upper_sf = halfcauchy_lsf(b, params)
            if (upper_sf - lower_sf <= -scifort_log_two) then
                ! Wide intervals retain tiny omitted probabilities near unit mass.
                value = lower_sf + log1p_safe(-exp(upper_sf - lower_sf))
            else
                value = scifort_log_two + cauchy_interval(a, b, params)
            end if
        end if
    end function halfcauchy_interval

    pure function halflogistic_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, halflogistic_lpdf, halflogistic_lcdf, &
            halflogistic_lsf, halflogistic_interval)
    end function halflogistic_censored_loglikelihood

    pure function halflogistic_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -halflogistic_censored_loglikelihood(data, loc, scale)
    end function halflogistic_censored_nnlf

    pure function halflogistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x < params(1)) return
        if (x == params(1)) then
            value = -scifort_log_two - log(params(2))
        else
            value = scifort_log_two + logistic_lpdf(x, params)
        end if
    end function halflogistic_lpdf

    pure function halflogistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halflogistic_tails(x, params, value, other)
    end function halflogistic_lcdf

    pure function halflogistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, other
        call halflogistic_tails(x, params, other, value)
    end function halflogistic_lsf

    pure function halflogistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp) :: value, lower_sf, upper_sf
        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        if (a <= params(1)) then
            value = halflogistic_lcdf(b, params)
        else
            lower_sf = halflogistic_lsf(a, params)
            if (.not. ieee_is_finite(lower_sf)) return
            upper_sf = halflogistic_lsf(b, params)
            if (upper_sf - lower_sf <= -scifort_log_two) then
                ! Wide intervals retain tiny omitted probabilities near unit mass.
                value = lower_sf + log1p_safe(-exp(upper_sf - lower_sf))
            else
                value = scifort_log_two + logistic_interval(a, b, params)
            end if
        end if
    end function halflogistic_interval

    pure subroutine halfcauchy_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp), intent(out) :: logcdf !! log lower-tail probability, -infinity at/below support
        real(dp), intent(out) :: logsf !! log upper-tail probability, zero at/below support
        real(dp) :: logz
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        if (x <= params(1)) return
        logz = cauchy_log_distance(x, params)
        if (logz <= 0.0_dp) then
            logcdf = scifort_log_two + cauchy_logatan_exp(logz) - log(scifort_pi)
            logsf = log1p_safe(-exp(logcdf))
        else
            logsf = scifort_log_two + cauchy_log_small_tail(logz)
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine halfcauchy_tails

    pure subroutine halflogistic_tails(x, params, logcdf, logsf)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite lower support endpoint and positive finite scale
        real(dp), intent(out) :: logcdf !! log lower-tail probability, -infinity at/below support
        real(dp), intent(out) :: logsf !! log upper-tail probability, zero at/below support
        real(dp) :: logz
        logcdf = negative_infinity(0.0_dp)
        logsf = 0.0_dp
        if (x <= params(1)) return
        logz = cauchy_log_distance(x, params)
        if (logz <= 0.0_dp) then
            ! F(z) = (1-exp(-z))/(1+exp(-z)); retain log(z) if z underflows.
            logcdf = log_exponential_mass(logz) - log1pexp(-exp(logz))
            logsf = log1p_safe(-exp(logcdf))
        else
            logsf = scifort_log_two + logistic_lsf(x, params)
            logcdf = log1p_safe(-exp(logsf))
        end if
    end subroutine halflogistic_tails


    pure function cauchy_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite median (default 0)
        real(dp), intent(in), optional :: scale !! positive finite half width at half maximum (default 1)
        real(dp) :: value, params(2)
        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, cauchy_lpdf, cauchy_lcdf, &
            cauchy_lsf, cauchy_interval)
    end function cauchy_censored_loglikelihood

    pure function cauchy_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! finite median (default 0)
        real(dp), intent(in), optional :: scale !! positive finite half width at half maximum (default 1)
        real(dp) :: value
        value = -cauchy_censored_loglikelihood(data, loc, scale)
    end function cauchy_censored_nnlf

    pure function cauchy_log_distance(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x == params(1)) return
        value = log_distance(max(x, params(1)), min(x, params(1))) - log(params(2))
    end function cauchy_log_distance

    pure function cauchy_logatan_exp(logx) result(value)
        implicit none
        real(dp), intent(in) :: logx !! log positive argument; negative infinity means zero
        real(dp) :: value
        if (logx < 0.5_dp * log(epsilon(1.0_dp))) then
            ! atan(t)/t differs from 1 by at most t**2/3 for 0 <= t <= 1.
            value = logx
        else if (logx <= 0.0_dp) then
            value = log(atan(exp(logx)))
        else
            ! Reciprocal identity keeps EXP bounded and subtraction well conditioned.
            value = log(0.5_dp * scifort_pi - atan(exp(-logx)))
        end if
    end function cauchy_logatan_exp

    pure function cauchy_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value, logz
        logz = cauchy_log_distance(x, params)
        value = -log(scifort_pi) - log(params(2)) - log1pexp(2.0_dp * logz)
    end function cauchy_lpdf

    pure function cauchy_log_small_tail(logz) result(value)
        implicit none
        real(dp), intent(in) :: logz !! log absolute standardized displacement, possibly -infinity
        real(dp) :: value
        if (logz >= 0.0_dp) then
            value = cauchy_logatan_exp(-logz) - log(scifort_pi)
        else
            value = log(0.5_dp - atan(exp(logz)) / scifort_pi)
        end if
    end function cauchy_log_small_tail

    pure function cauchy_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value
        value = cauchy_log_small_tail(cauchy_log_distance(x, params))
        if (x > params(1)) value = log1p_safe(-exp(value))
    end function cauchy_lcdf

    pure function cauchy_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value
        value = cauchy_log_small_tail(cauchy_log_distance(x, params))
        if (x < params(1)) value = log1p_safe(-exp(value))
    end function cauchy_lsf

    pure function cauchy_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite median and positive finite scale
        real(dp) :: value, loga, logb, logratio, left, right, large, small
        loga = cauchy_log_distance(a, params)
        logb = cauchy_log_distance(b, params)
        if (b <= params(1) .or. a >= params(1)) then
            ! Same-side arctangent difference has a positive denominator:
            ! tan(angle) = ((b-a)/scale) / (1 + za*zb).
            logratio = log_distance(b, a) - log(params(2)) - log1pexp(loga + logb)
            value = cauchy_logatan_exp(logratio) - log(scifort_pi)
        else if (min(loga, logb) >= 0.0_dp) then
            ! Wide crossing interval: retain even tiny missing tail probabilities.
            left = exp(cauchy_log_small_tail(loga))
            right = exp(cauchy_log_small_tail(logb))
            value = log1p_safe(-(left + right))
        else
            ! Narrow crossing interval: add positive angles, never subtract CDFs.
            left = cauchy_logatan_exp(loga)
            right = cauchy_logatan_exp(logb)
            large = max(left, right)
            small = min(left, right)
            value = min(0.0_dp, large + log1p_safe(exp(small - large)) - log(scifort_pi))
        end if
    end function cauchy_interval


    pure function laplace_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite mean absolute deviation (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, laplace_lpdf, laplace_lcdf, &
            laplace_lsf, laplace_interval)
    end function laplace_censored_loglikelihood

    pure function laplace_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite mean absolute deviation (default 1)
        real(dp) :: value
        value = -laplace_censored_loglikelihood(data, loc, scale)
    end function laplace_censored_nnlf

    pure function laplace_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        value = -abs(location_scale_z(x, params)) - scifort_log_two - log(params(2))
    end function laplace_lpdf

    pure function laplace_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x, params)
        if (z <= 0.0_dp) then
            value = z - scifort_log_two
        else
            value = log1p_safe(-0.5_dp * exp(-z))
        end if
    end function laplace_lcdf

    pure function laplace_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, z
        z = location_scale_z(x, params)
        if (z >= 0.0_dp) then
            value = -z - scifort_log_two
        else
            value = log1p_safe(-0.5_dp * exp(z))
        end if
    end function laplace_lsf

    pure function laplace_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, left, right, large, small, logscale
        logscale = log(params(2))
        if (b <= params(1)) then
            ! Left-side exponential mass, anchored at the endpoint nearest the center.
            value = laplace_lcdf(b, params) + &
                log_exponential_mass(log_distance(b, a) - logscale)
        else if (a >= params(1)) then
            value = laplace_lsf(a, params) + &
                log_exponential_mass(log_distance(b, a) - logscale)
        else
            ! Split at the cusp; log-add the two positive half-masses.
            left = log_exponential_mass(log_distance(params(1), a) - logscale)
            right = log_exponential_mass(log_distance(b, params(1)) - logscale)
            large = max(left, right)
            small = min(left, right)
            value = large + log1p_safe(0.5_dp * expm1_safe(small - large))
        end if
    end function laplace_interval


    pure function gumbel_r_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(3)

        params = [0.0_dp, 1.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, gumbel_lpdf, gumbel_lcdf, &
            gumbel_lsf, gumbel_interval)
    end function gumbel_r_censored_loglikelihood

    pure function gumbel_r_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -gumbel_r_censored_loglikelihood(data, loc, scale)
    end function gumbel_r_censored_nnlf

    pure function gumbel_l_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(3)

        params = [0.0_dp, 1.0_dp, -1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, gumbel_lpdf, gumbel_lcdf, &
            gumbel_lsf, gumbel_interval)
    end function gumbel_l_censored_loglikelihood

    pure function gumbel_l_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! location (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -gumbel_l_censored_loglikelihood(data, loc, scale)
    end function gumbel_l_censored_nnlf

    pure function gumbel_right_logcdf(z) result(value)
        implicit none
        real(dp), intent(in) :: z !! standardized point for the right-skewed family
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (-z > log(huge(1.0_dp))) return
        value = -exp(-z)
    end function gumbel_right_logcdf

    pure function gumbel_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! exact observation
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, z
        z = params(3) * location_scale_z(x, params)
        value = gumbel_right_logcdf(z)
        if (.not. ieee_is_finite(value)) return
        value = value - z - log(params(2))
    end function gumbel_lpdf

    pure function gumbel_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring limit
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, z
        z = params(3) * location_scale_z(x, params)
        if (params(3) > 0.0_dp) then
            value = gumbel_right_logcdf(z)
        else
            value = log_exponential_mass(-z)
        end if
    end function gumbel_lcdf

    pure function gumbel_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring limit
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, z
        z = params(3) * location_scale_z(x, params)
        if (params(3) > 0.0_dp) then
            value = log_exponential_mass(-z)
        else
            value = gumbel_right_logcdf(z)
        end if
    end function gumbel_lsf

    pure function gumbel_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! location, positive scale, direction (+1 or -1)
        real(dp) :: value, za, zb, logwidth, logincrement
        if (params(3) > 0.0_dp) then
            za = location_scale_z(a, params)
            zb = location_scale_z(b, params)
        else
            ! Reflection interchanges the interval endpoints and the two tails.
            za = -location_scale_z(b, params)
            zb = -location_scale_z(a, params)
        end if
        logwidth = log_distance(b, a) - log(params(2))
        ! exp(-za)-exp(-zb) = exp(-za)*(1-exp(-(b-a)/scale)).
        logincrement = -za + log_exponential_mass(logwidth)
        value = gumbel_right_logcdf(zb) + log_exponential_mass(logincrement)
    end function gumbel_interval

    pure function logistic_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (params(2) <= 0.0_dp) return
        value = censored_loglikelihood(data, params, logistic_lpdf, logistic_lcdf, &
            logistic_lsf, logistic_interval)
    end function logistic_censored_loglikelihood

    pure function logistic_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! center (default 0)
        real(dp), intent(in), optional :: scale !! positive finite scale (default 1)
        real(dp) :: value
        value = -logistic_censored_loglikelihood(data, loc, scale)
    end function logistic_censored_nnlf

    pure function location_scale_z(x, params) result(z)
        implicit none
        real(dp), intent(in) :: x !! evaluation point, possibly infinite
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: z, delta, factor, ratio
        integer :: power

        if (.not. ieee_is_finite(x)) then
            z = x
            return
        end if
        factor = 1.0_dp
        if ((x > 0.0_dp .and. params(1) < 0.0_dp) .or. &
            (x < 0.0_dp .and. params(1) > 0.0_dp)) then
            if (abs(x) >= huge(x) - abs(params(1))) factor = 2.0_dp
        end if
        ! Halve near subtraction overflow; retain nearby differences elsewhere.
        delta = x / factor - params(1) / factor
        if (delta == 0.0_dp) then
            z = delta
            return
        end if
        z = sign(positive_infinity(0.0_dp), delta)
        ! Divide bounded significands, then scale by a checked binary exponent.
        ! A rounded huge*scale threshold can otherwise miss division overflow.
        ratio = fraction(delta) / fraction(params(2))
        power = exponent(delta) - exponent(params(2)) + exponent(ratio)
        if (factor == 2.0_dp) power = power + 1
        if (power > maxexponent(z)) return
        z = scale(fraction(ratio), power)
    end function location_scale_z

    pure function logistic_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! exact observation
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value, a
        a = abs(location_scale_z(x, params))
        value = -a - 2.0_dp * log1p_safe(exp(-a)) - log(params(2))
    end function logistic_lpdf

    pure function logistic_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! upper censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        value = -log1pexp(-location_scale_z(x, params))
    end function logistic_lcdf

    pure function logistic_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! lower censoring limit
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        value = -log1pexp(location_scale_z(x, params))
    end function logistic_lsf

    pure function logistic_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! finite center and positive finite scale
        real(dp) :: value
        ! F(b)-F(a) = F(b)*S(a)*(1-exp(-(b-a)/scale)).
        value = logistic_lcdf(b, params) + logistic_lsf(a, params) + &
            log_exponential_mass(log_distance(b, a) - log(params(2)))
    end function logistic_interval

    function make_censored_data(uncensored, left, right, interval, status) result(data)
        implicit none
        real(dp), intent(in), optional :: uncensored(:) !! finite exact observations
        real(dp), intent(in), optional :: left(:) !! finite upper censoring limits
        real(dp), intent(in), optional :: right(:) !! finite lower censoring limits
        real(dp), intent(in), optional :: interval(:,:) !! (n,2) lower/upper limits
        integer, intent(out), optional :: status !! 0 valid, 1 invalid; invalid data is retained
        type(censored_data) :: data

        if (present(uncensored)) data%uncensored = uncensored
        if (present(left)) data%left = left
        if (present(right)) data%right = right
        if (present(interval)) data%interval = interval
        if (present(status)) then
            status = censored_status_invalid_input
            if (censored_data_valid(data)) status = censored_status_success
        end if
    end function make_censored_data

    pure function censored_data_valid(data) result(valid)
        implicit none
        type(censored_data), intent(in) :: data !! observations; unallocated groups are empty
        logical :: valid
        real(dp) :: a, b
        integer :: i

        valid = .false.
        if (allocated(data%uncensored)) then
            if (.not. all(ieee_is_finite(data%uncensored))) return
        end if
        if (allocated(data%left)) then
            if (.not. all(ieee_is_finite(data%left))) return
        end if
        if (allocated(data%right)) then
            if (.not. all(ieee_is_finite(data%right))) return
        end if
        if (allocated(data%interval)) then
            if (size(data%interval, 2) /= 2) return
            do i = 1, size(data%interval, 1)
                a = data%interval(i,1)
                b = data%interval(i,2)
                if (ieee_is_nan(a) .or. ieee_is_nan(b)) return
                if (a > b) return
                if (.not. ieee_is_finite(a) .and. .not. ieee_is_finite(b)) return
            end do
        end if
        valid = .true.
    end function censored_data_valid

    pure function censored_size(data) result(n)
        implicit none
        type(censored_data), intent(in) :: data !! observations; unallocated groups are empty
        integer :: n
        n = 0
        if (allocated(data%uncensored)) n = n + size(data%uncensored)
        if (allocated(data%left)) n = n + size(data%left)
        if (allocated(data%right)) n = n + size(data%right)
        if (allocated(data%interval)) n = n + size(data%interval, 1)
    end function censored_size

    pure function censored_loglikelihood(data, params, logpdf, logcdf, logsf, &
            loginterval) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: params(:) !! finite parameters; callbacks validate their domain
        procedure(censored_log_function) :: logpdf !! log density callback
        procedure(censored_log_function) :: logcdf !! stable log lower-tail callback
        procedure(censored_log_function) :: logsf !! stable log upper-tail callback
        procedure(censored_log_interval), optional :: loginterval !! optional stable interval mass
        real(dp) :: value
        real(dp) :: a, b, term
        integer :: i

        value = quiet_nan(0.0_dp)
        if (.not. censored_data_valid(data)) return
        if (.not. all(ieee_is_finite(params))) return
        value = 0.0_dp
        if (allocated(data%uncensored)) then
            do i = 1, size(data%uncensored)
                value = add_log_term(value, logpdf(data%uncensored(i), params))
            end do
        end if
        if (allocated(data%left)) then
            do i = 1, size(data%left)
                value = add_log_term(value, logcdf(data%left(i), params))
            end do
        end if
        if (allocated(data%right)) then
            do i = 1, size(data%right)
                value = add_log_term(value, logsf(data%right(i), params))
            end do
        end if
        if (allocated(data%interval)) then
            do i = 1, size(data%interval, 1)
                a = data%interval(i,1)
                b = data%interval(i,2)
                if (a == b) then
                    term = logpdf(a, params)
                else if (.not. ieee_is_finite(a)) then
                    term = logcdf(b, params)
                else if (.not. ieee_is_finite(b)) then
                    term = logsf(a, params)
                else if (present(loginterval)) then
                    term = loginterval(a, b, params)
                else
                    term = tail_interval(a, b, params, logcdf, logsf)
                end if
                value = add_log_term(value, term)
            end do
        end if
    end function censored_loglikelihood

    pure function add_log_term(total, term) result(value)
        implicit none
        real(dp), intent(in) :: total !! accumulated log likelihood, possibly infinite
        real(dp), intent(in) :: term !! next contribution, possibly infinite
        real(dp) :: value
        value = quiet_nan(0.0_dp)
        if (ieee_is_nan(total) .or. ieee_is_nan(term)) return
        if (.not. ieee_is_finite(total) .and. .not. ieee_is_finite(term)) then
            if (total /= term) return
        end if
        value = total + term
    end function add_log_term

    pure function log_difference(large, small) result(value)
        implicit none
        real(dp), intent(in) :: large !! log of larger nonnegative probability
        real(dp), intent(in) :: small !! log of smaller nonnegative probability
        real(dp) :: value

        value = quiet_nan(large)
        if (ieee_is_nan(large) .or. ieee_is_nan(small)) return
        if (small > large) return
        if (large == small) then
            value = negative_infinity(large)
        else
            value = large + log(-expm1_safe(small - large))
        end if
    end function log_difference

    pure function tail_interval(a, b, params, logcdf, logsf) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! distribution parameters
        procedure(censored_log_function) :: logcdf !! stable log lower-tail callback
        procedure(censored_log_function) :: logsf !! stable log upper-tail callback
        real(dp) :: value, upper_cdf, lower_sf

        upper_cdf = logcdf(b, params)
        lower_sf = logsf(a, params)
        if (upper_cdf <= lower_sf) then
            value = log_difference(upper_cdf, logcdf(a, params))
        else
            value = log_difference(lower_sf, logsf(b, params))
        end if
    end function tail_interval

    pure function log_distance(b, a) result(value)
        implicit none
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp) :: value, big, small

        if (a < 0.0_dp .and. b > 0.0_dp) then
            big = max(abs(a), abs(b))
            small = min(abs(a), abs(b))
            value = log(big) + log1p_safe(small / big)
        else
            value = log(b - a)
        end if
    end function log_distance

    pure function log_exponential_mass(logwidth) result(value)
        implicit none
        real(dp), intent(in) :: logwidth !! log positive cumulative-hazard increment
        real(dp) :: value

        if (logwidth < log(epsilon(1.0_dp)) - 2.0_dp) then
            value = logwidth
        else if (logwidth > log(huge(1.0_dp))) then
            value = 0.0_dp
        else
            value = log(-expm1_safe(-exp(logwidth)))
        end if
    end function log_exponential_mass

    pure function exponential_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! location and positive scale
        real(dp) :: value, lower

        value = negative_infinity(0.0_dp)
        if (b <= params(1)) return
        lower = max(a, params(1))
        value = exponential_logsf(lower, params(1), params(2)) + &
            log_exponential_mass(log_distance(b, lower) - log(params(2)))
    end function exponential_interval

    pure function normal_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! mean and positive standard deviation
        real(dp) :: value, z, logwidth, logscaled

        z = (a - params(1)) / params(2)
        logwidth = log_distance(b, a)
        logscaled = logwidth - log(params(2))
        value = normal_small_interval(z, logscaled)
        if (ieee_is_nan(value)) value = tail_interval(a, b, params, normal_lcdf, normal_lsf)
    end function normal_interval

    pure function normal_small_interval(lower_z, logscaled) result(value)
        implicit none
        real(dp), intent(in) :: lower_z !! standardized lower endpoint
        real(dp), intent(in) :: logscaled !! log standardized interval width
        real(dp) :: value, z, halfwidth, offset, correction

        ! Three-point Gauss-Legendre density integral; NaN requests tail fallback.
        ! Nodes 0,+/-sqrt(3/5), weights 8/9,5/9,5/9 (DLMF 3.5(v)).
        value = quiet_nan(0.0_dp)
        z = lower_z
        if (ieee_is_finite(z)) then
            if (logscaled + log(max(1.0_dp, abs(z))) <= log(0.01_dp)) then
                if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) then
                    value = negative_infinity(0.0_dp)
                    return
                end if
                halfwidth = 0.5_dp*exp(logscaled)
                z = z + halfwidth
                offset = halfwidth*sqrt(3.0_dp/5.0_dp)
                correction = (5.0_dp/18.0_dp)*( &
                    expm1_safe(-z*offset - 0.5_dp*offset**2) + &
                    expm1_safe(z*offset - 0.5_dp*offset**2))
                value = (-0.5_dp*z)*z - scifort_log_sqrt_two_pi + &
                    logscaled + log1p_safe(correction)
                return
            end if
        end if
    end function normal_small_interval

    pure function lognormal_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logr, logspan, z

        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = lognormal_lcdf(b, params)
            return
        end if
        ! log(b-loc)-log(a-loc) loses adjacent endpoint widths. Instead use
        ! log1p((b-a)/(a-loc)), retaining its logarithm if the ratio underflows.
        logr = log_distance(b, a) - log_distance(a, params(2))
        if (logr < log(epsilon(1.0_dp)) - 2.0_dp) then
            logspan = logr
        else
            logspan = log(log1pexp(logr))
        end if
        z = lognormal_standardized(a, params)
        value = normal_small_interval(z, logspan - log(params(1)))
        if (ieee_is_nan(value)) then
            value = tail_interval(a, b, params, lognormal_lcdf, lognormal_lsf)
        end if
    end function lognormal_interval

    pure function lognormal_standardized(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz
        logz = lognormal_logz(x, params)
        if (params(1) < 1.0_dp) then
            if (abs(logz) > huge(1.0_dp)*params(1)) then
                value = sign(positive_infinity(0.0_dp), logz)
                return
            end if
        end if
        value = logz/params(1)
    end function lognormal_standardized

    pure function lognormal_logz(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value
        ! Reuse the stable log distance/scale calculation, with unit exponent.
        value = weibull_log_hazard(x, [1.0_dp, params(2), params(3)])
    end function lognormal_logz

    pure function weibull_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, lower, loglo, loghi, logratio, logr, delta, logdelta

        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        lower = max(a, params(2))
        loghi = weibull_log_hazard(b, params)
        if (lower == params(2)) then
            value = log_exponential_mass(loghi)
            return
        end if
        loglo = weibull_log_hazard(lower, params)
        if (loglo > log(huge(1.0_dp))) return
        logratio = log_distance(b, lower) - log_distance(lower, params(2))
        if (logratio < log(epsilon(1.0_dp)) - 2.0_dp) then
            logr = log(params(1)) + logratio
        else
            logr = log(params(1)) + log(log1pexp(logratio))
        end if
        ! Here the omitted factor exp(-delta) is below exp(-50) < 2e-22.
        if (logr > log(50.0_dp)) then
            logdelta = loghi
        else if (logr < log(epsilon(1.0_dp)) - 2.0_dp) then
            logdelta = loglo + logr
        else
            delta = exp(logr)
            logdelta = loglo + log(expm1_safe(delta))
        end if
        value = -exp(loglo) + log_exponential_mass(logdelta)
    end function weibull_interval

    pure function weibull_log_hazard(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz, distance
        logical :: finite_distance
        logz = log_distance(x, params(2)) - log(params(3))
        ! Avoid subtracting nearly equal large logarithms when the normalized
        ! distance is near one: a large shape would amplify that lost precision.
        finite_distance = .true.
        if (params(2) < 0.0_dp .and. x > 0.0_dp) then
            finite_distance = x <= huge(1.0_dp) + params(2)
        end if
        if (finite_distance) then
            distance = x - params(2)
            if (0.5_dp*distance <= params(3) .and. distance >= 0.5_dp*params(3)) then
                logz = log1p_safe((distance - params(3))/params(3))
            end if
        end if
        if (params(1) > 1.0_dp) then
            if (abs(logz) > huge(1.0_dp)/params(1)) then
                value = sign(positive_infinity(0.0_dp), logz)
                return
            end if
        end if
        value = params(1)*logz
    end function weibull_log_hazard

    pure function normal_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(1), params(2))) return
        value = censored_loglikelihood(data, params, normal_lpdf, normal_lcdf, &
            normal_lsf, normal_interval)
    end function normal_censored_loglikelihood

    pure function exponential_censored_loglikelihood(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! inverse rate, > 0 (default 1)
        real(dp) :: value, params(2)

        params = [0.0_dp, 1.0_dp]
        if (present(loc)) params(1) = loc
        if (present(scale)) params(2) = scale
        value = quiet_nan(0.0_dp)
        if (.not. valid_loc_scale(params(1), params(2))) return
        value = censored_loglikelihood(data, params, exponential_lpdf, exponential_lcdf, &
            exponential_lsf, exponential_interval)
    end function exponential_censored_loglikelihood

    pure function weibull_censored_loglikelihood(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value, params(3)

        params = [c, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. ieee_is_finite(c)) return
        if (c <= 0.0_dp .or. .not. valid_loc_scale(params(2), params(3))) return
        value = censored_loglikelihood(data, params, weibull_lpdf, weibull_lcdf, &
            weibull_lsf, weibull_interval)
    end function weibull_censored_loglikelihood

    pure function normal_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! mean (default 0)
        real(dp), intent(in), optional :: scale !! standard deviation, > 0 (default 1)
        real(dp) :: value
        value = -normal_censored_loglikelihood(data, loc, scale)
    end function normal_censored_nnlf

    pure function exponential_censored_nnlf(data, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! inverse rate, > 0 (default 1)
        real(dp) :: value
        value = -exponential_censored_loglikelihood(data, loc, scale)
    end function exponential_censored_nnlf

    pure function weibull_censored_nnlf(data, c, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: c !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! scale, > 0 (default 1)
        real(dp) :: value
        value = -weibull_censored_loglikelihood(data, c, loc, scale)
    end function weibull_censored_nnlf

    pure function normal_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! mean and standard deviation
        real(dp) :: value
        value = normal_logpdf(x, params(1), params(2))
    end function normal_lpdf

    pure function normal_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! mean and standard deviation
        real(dp) :: value
        value = normal_logcdf(x, params(1), params(2))
    end function normal_lcdf

    pure function normal_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! mean and standard deviation
        real(dp) :: value
        value = normal_logsf(x, params(1), params(2))
    end function normal_lsf

    pure function exponential_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! location and scale
        real(dp) :: value
        value = exponential_logpdf(x, params(1), params(2))
    end function exponential_lpdf

    pure function exponential_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! location and scale
        real(dp) :: value
        value = exponential_logcdf(x, params(1), params(2))
    end function exponential_lcdf

    pure function exponential_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! location and scale
        real(dp) :: value
        value = exponential_logsf(x, params(1), params(2))
    end function exponential_lsf

    pure function weibull_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! shape, location, scale
        real(dp) :: value
        value = weibull_logpdf(x, params(1), params(2), params(3))
    end function weibull_lpdf

    pure function weibull_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! shape, location, scale
        real(dp) :: value
        value = weibull_logcdf(x, params(1), params(2), params(3))
    end function weibull_lcdf

    pure function weibull_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! evaluation point
        real(dp), intent(in) :: params(:) !! shape, location, scale
        real(dp) :: value
        value = weibull_logsf(x, params(1), params(2), params(3))
    end function weibull_lsf

    pure function lognormal_censored_loglikelihood(data, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: s !! finite log standard deviation, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite scale, > 0 (default 1)
        real(dp) :: value, params(3)

        params = [s, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (s <= 0.0_dp .or. .not. valid_loc_scale(params(2), params(3))) return
        value = censored_loglikelihood(data, params, lognormal_lpdf, lognormal_lcdf, &
            lognormal_lsf, lognormal_interval)
    end function lognormal_censored_loglikelihood

    pure function lognormal_censored_nnlf(data, s, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: s !! finite log standard deviation, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite scale, > 0 (default 1)
        real(dp) :: value
        value = -lognormal_censored_loglikelihood(data, s, loc, scale)
    end function lognormal_censored_nnlf

    pure function lognormal_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, z
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        z = lognormal_standardized(x, params)
        if (abs(z) > sqrt(huge(1.0_dp))*sqrt(2.0_dp)) return
        value = normal_logpdf(z) - log(params(1)) - log_distance(x, params(2))
    end function lognormal_lpdf

    pure function lognormal_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        value = normal_logcdf(lognormal_standardized(x, params))
    end function lognormal_lcdf

    pure function lognormal_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value
        value = 0.0_dp
        if (x <= params(2)) return
        value = normal_logsf(lognormal_standardized(x, params))
    end function lognormal_lsf

    pure function gamma_censored_loglikelihood(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite inverse rate, > 0 (default 1)
        real(dp) :: value, params(3)
        params = [a, 0.0_dp, 1.0_dp]
        if (present(loc)) params(2) = loc
        if (present(scale)) params(3) = scale
        value = quiet_nan(0.0_dp)
        if (.not. all(ieee_is_finite(params))) return
        if (a <= 0.0_dp .or. params(3) <= 0.0_dp) return
        if (a == 1.0_dp) then
            value = exponential_censored_loglikelihood(data, params(2), params(3))
        else
            value = censored_loglikelihood(data, params, gamma_lpdf, gamma_lcdf, &
                gamma_lsf, gamma_interval)
        end if
    end function gamma_censored_loglikelihood

    pure function gamma_censored_nnlf(data, a, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite inverse rate, > 0 (default 1)
        real(dp) :: value
        value = -gamma_censored_loglikelihood(data, a, loc, scale)
    end function gamma_censored_nnlf

    pure function gamma_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz, distance_log, z
        value = negative_infinity(0.0_dp)
        if (x < params(2)) return
        if (x == params(2)) then
            value = gamma_logpdf(x, params(1), params(2), params(3))
            return
        end if
        distance_log = log_distance(x, params(2))
        logz = distance_log - log(params(3))
        if (logz > log(huge(1.0_dp))) return
        if (logz < log(tiny(1.0_dp))) then
            ! The omitted exp(-z) correction is below the smallest normal real.
            value = log(params(1))+gamma_small_kernel(params(1),logz)-distance_log
        else
            z = gamma_scaled(x,params,logz)
            if (params(1) < 10.0_dp) then
                value = (params(1)-1.0_dp)*logz-z-log_gamma(params(1))-log(params(3))
            else
                value = log(params(1))+log_gamma_kernel(params(1),z)-distance_log
            end if
        end if
    end function gamma_lpdf

    pure function gamma_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz
        value = negative_infinity(0.0_dp)
        if (x <= params(2)) return
        logz = log_distance(x, params(2))-log(params(3))
        if (logz > log(huge(1.0_dp))) then
            value = 0.0_dp
        else if (logz < log(tiny(1.0_dp))) then
            ! P(a,z) = z**a/Gamma(a+1) * (1 + O(z)); DLMF 8.7.
            value = gamma_small_kernel(params(1),logz)
        else
            value = log_gammainc(params(1), gamma_scaled(x,params,logz))
        end if
    end function gamma_lcdf

    pure function gamma_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logz, logp
        value = 0.0_dp
        if (x <= params(2)) return
        logz = log_distance(x, params(2))-log(params(3))
        if (logz > log(huge(1.0_dp))) then
            value = negative_infinity(0.0_dp)
        else if (logz < log(tiny(1.0_dp))) then
            logp = gamma_lcdf(x, params)
            value = log_exponential_mass(log(-logp))
        else
            value = log_gammaincc(params(1), gamma_scaled(x,params,logz))
        end if
    end function gamma_lsf

    pure function gamma_small_kernel(a, logz) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite positive shape
        real(dp), intent(in) :: logz !! log standardized argument below log(tiny)
        real(dp) :: value, factor
        if (a < 10.0_dp) then
            value = a*logz-log_gamma_one_plus(a)
        else
            ! Stirling form avoids overflow in Gamma(a+1) and in a*log(z).
            factor = -logz+log(a)-1.0_dp
            value = negative_infinity(0.0_dp)
            if (a > huge(1.0_dp)/factor) return
            value = -a*factor-0.5_dp*log(a)-scifort_log_sqrt_two_pi-stirling_remainder(a)
        end if
    end function gamma_small_kernel

    pure function gamma_scaled(x, params, logz) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly above location
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp), intent(in) :: logz !! log standardized distance within normal real range
        real(dp) :: value
        ! Direct division retains near-mode accuracy for large shapes. Use the
        ! logarithm only if the unscaled distance or ratio could overflow.
        value = exp(logz)
        if (logz > log(huge(1.0_dp))-1.0_dp) return
        if (params(2) < 0.0_dp .and. x > 0.0_dp) then
            if (x > huge(1.0_dp)+params(2)) return
        end if
        value = (x-params(2))/params(3)
    end function gamma_scaled

    pure function gamma_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper endpoint, > a
        real(dp), intent(in) :: params(:) !! positive shape, location, positive scale
        real(dp) :: value, logwidth, logrelative, logscaled, q, d, correction, t(3), r(3)
        value = negative_infinity(0.0_dp)
        if (b <= params(2)) return
        if (a <= params(2)) then
            value = gamma_lcdf(b, params)
            return
        end if
        logwidth = log_distance(b,a)
        logrelative = logwidth-log_distance(a,params(2))
        logscaled = logwidth-log(params(3))
        ! Bound both relative endpoint separation and density variation.
        ! Three-point Gauss-Legendre is exact through degree five (DLMF 3.5).
        if (logrelative+log(max(1.0_dp,abs(params(1)-1.0_dp))) <= log(0.001_dp) &
                .and. logscaled <= log(0.001_dp)) then
            q = exp(logrelative)
            d = exp(logscaled)
            t = [0.5_dp, 0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = (params(1)-1.0_dp)*log1p_safe(t*q)-t*d
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) + &
                (5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = gamma_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,gamma_lcdf,gamma_lsf)
        end if
    end function gamma_interval

    pure function beta_censored_loglikelihood(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first shape, > 0
        real(dp), intent(in) :: b !! finite second shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite support width, > 0 (default 1)
        real(dp) :: value, params(4)
        params = [a, b, 0.0_dp, 1.0_dp]
        if (present(loc)) params(3) = loc
        if (present(scale)) params(4) = scale
        value = quiet_nan(0.0_dp)
        if (.not. beta_parameters_valid(params)) return
        value = censored_loglikelihood(data,params,beta_lpdf,beta_lcdf,beta_lsf,beta_interval)
    end function beta_censored_loglikelihood

    pure function beta_censored_nnlf(data, a, b, loc, scale) result(value)
        implicit none
        type(censored_data), intent(in) :: data !! exact and censored observations
        real(dp), intent(in) :: a !! finite first shape, > 0
        real(dp), intent(in) :: b !! finite second shape, > 0
        real(dp), intent(in), optional :: loc !! lower support endpoint (default 0)
        real(dp), intent(in), optional :: scale !! finite support width, > 0 (default 1)
        real(dp) :: value
        value = -beta_censored_loglikelihood(data,a,b,loc,scale)
    end function beta_censored_nnlf

    pure function beta_parameters_valid(params) result(valid)
        implicit none
        real(dp), intent(in) :: params(:) !! first/second shapes, location, support width
        logical :: valid
        valid = .false.
        if (.not. all(ieee_is_finite(params))) return
        if (min(params(1),params(2),params(4)) <= 0.0_dp) return
        ! This initial interface requires a finite rounded upper support edge.
        if (params(3) > 0.0_dp) then
            if (params(4) > huge(1.0_dp)-params(3)) return
        end if
        valid = .true.
    end function beta_parameters_valid

    pure subroutine beta_distances(x, params, left, right)
        implicit none
        real(dp), intent(in) :: x !! finite observation or censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp), intent(out) :: left !! distance above lower edge; negative outside support
        real(dp), intent(out) :: right !! distance below upper edge; negative outside support
        real(dp) :: upper, correction
        ! FastTwoSum with operands ordered by magnitude (Ogita et al. 2005).
        upper = params(3)+params(4)
        if (abs(params(3)) >= params(4)) then
            correction = (params(3)-upper)+params(4)
        else
            correction = (params(4)-upper)+params(3)
        end if
        left = 0.0_dp
        right = params(4)
        if (x == params(3)) return
        if (x < params(3)) then
            left = -1.0_dp
            return
        end if
        left = positive_infinity(0.0_dp)
        if (params(3) < 0.0_dp .and. x > 0.0_dp) then
            if (x <= huge(1.0_dp)+params(3)) left = x-params(3)
        else
            left = x-params(3)
        end if
        if (x > 0.0_dp .and. upper < 0.0_dp) then
            if (x > huge(1.0_dp)+upper) then
                right = -1.0_dp
                return
            end if
        else if (x < 0.0_dp .and. upper > 0.0_dp) then
            if (upper > huge(1.0_dp)+x) then
                right = huge(1.0_dp)
                return
            end if
        end if
        right = (upper-x)+correction
    end subroutine beta_distances

    pure subroutine beta_coordinates(x, params, left, right, z, w, logz, logw)
        implicit none
        real(dp), intent(in) :: x !! finite point strictly inside support
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp), intent(in) :: left !! positive distance above lower edge, possibly infinite
        real(dp), intent(in) :: right !! positive distance below upper edge
        real(dp), intent(out) :: z !! standardized lower distance, possibly rounded to 0 or 1
        real(dp), intent(out) :: w !! standardized upper distance, formed independently of z
        real(dp), intent(out) :: logz !! finite log standardized lower distance
        real(dp), intent(out) :: logw !! finite log standardized upper distance
        if (ieee_is_finite(left)) then
            z = left/params(4)
        else
            z = x/params(4)-params(3)/params(4)
        end if
        w = right/params(4)
        logz = log_distance(x,params(3))-log(params(4))
        logw = log(right)-log(params(4))
        if (w < 0.5_dp) logz = log1p_safe(-w)
        if (z < 0.5_dp) logw = log1p_safe(-z)
    end subroutine beta_coordinates

    pure function beta_lpdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite exact observation
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, left, right, z, w, logz, logw, shape, other
        value = negative_infinity(0.0_dp)
        call beta_distances(x,params,left,right)
        if (left < 0.0_dp .or. right < 0.0_dp) return
        if (left == 0.0_dp .or. right == 0.0_dp) then
            shape = params(1)
            other = params(2)
            if (right == 0.0_dp) then
                shape = params(2)
                other = params(1)
            end if
            if (shape < 1.0_dp) value = positive_infinity(0.0_dp)
            if (shape == 1.0_dp) value = log(other)-log(params(4))
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        if (min(params(1),params(2)) >= stirling_threshold .and. min(z,w) > 0.0_dp) then
            value = log_beta_kernel(params(1),params(2),z,w)-logz-logw-log(params(4))
        else
            value = (params(1)-1.0_dp)*logz+(params(2)-1.0_dp)*logw - &
                log_beta(params(1),params(2))-log(params(4))
        end if
    end function beta_lpdf

    pure subroutine beta_logtails(x, params, logp, logq)
        implicit none
        real(dp), intent(in) :: x !! finite censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp), intent(out) :: logp !! stable log lower probability
        real(dp), intent(out) :: logq !! stable log upper probability
        real(dp) :: left, right, z, w, logz, logw, p, q
        logp = negative_infinity(0.0_dp)
        logq = 0.0_dp
        call beta_distances(x,params,left,right)
        if (left <= 0.0_dp) return
        if (right <= 0.0_dp) then
            logp = 0.0_dp
            logq = negative_infinity(0.0_dp)
            return
        end if
        call beta_coordinates(x,params,left,right,z,w,logz,logw)
        if (z == 0.0_dp) then
            logp = quiet_nan(0.0_dp)
            logq = logp
            if (log(params(2))+logz > log(epsilon(1.0_dp))-2.0_dp) return
            logp = beta_small_tail(params(1),params(2),logz)
            logq = log_difference(0.0_dp,logp)
        else if (w == 0.0_dp) then
            logp = quiet_nan(0.0_dp)
            logq = logp
            if (log(params(1))+logw > log(epsilon(1.0_dp))-2.0_dp) return
            logq = beta_small_tail(params(2),params(1),logw)
            logp = log_difference(0.0_dp,logq)
        else
            call incomplete_beta_xy(params(1),params(2),z,w,p,q,logp,logq)
        end if
    end subroutine beta_logtails

    pure function beta_small_tail(a, b, logz) result(value)
        implicit none
        real(dp), intent(in) :: a !! positive shape at the nearby endpoint
        real(dp), intent(in) :: b !! positive shape at the far endpoint
        real(dp), intent(in) :: logz !! log coordinate rounded to zero in ordinary arithmetic
        real(dp) :: value, normalization
        ! The gamma-difference form retains the O(a) normalization when a is tiny.
        if (a < sqrt(epsilon(1.0_dp))*min(b,1.0_dp)) then
            if (b > 1.0_dp/epsilon(1.0_dp)) then
                ! For a < sqrt(epsilon), the gamma-ratio correction is O(a/b),
                ! below roundoff in this tail; avoid powers of enormous b.
                normalization = a*log(b)-log_gamma_one_plus(a)
            else if (b >= stirling_threshold) then
                normalization = a*log(b)+log_gamma_ratio_scaled(b,a)-log_gamma_one_plus(a)
            else
                normalization = log_gamma_difference(b,a)-log_gamma_one_plus(a)
            end if
        else
            normalization = -log(a)-log_beta(a,b)
        end if
        value = a*logz+normalization
    end function beta_small_tail

    pure function beta_lcdf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite upper censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, logq
        call beta_logtails(x,params,value,logq)
    end function beta_lcdf

    pure function beta_lsf(x, params) result(value)
        implicit none
        real(dp), intent(in) :: x !! finite lower censoring limit
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, logp
        call beta_logtails(x,params,logp,value)
    end function beta_lsf

    pure function beta_interval(a, b, params) result(value)
        implicit none
        real(dp), intent(in) :: a !! finite lower interval endpoint
        real(dp), intent(in) :: b !! finite upper interval endpoint, > a
        real(dp), intent(in) :: params(:) !! valid first/second shapes, location, width
        real(dp) :: value, left, right, logwidth, logql, logqr, ql, qr, t(3), r(3), correction
        value = negative_infinity(0.0_dp)
        call beta_distances(b,params,left,right)
        if (left <= 0.0_dp) return
        call beta_distances(a,params,left,right)
        if (right <= 0.0_dp) return
        if (left <= 0.0_dp) then
            value = beta_lcdf(b,params)
            return
        end if
        call beta_distances(b,params,left,right)
        if (right <= 0.0_dp) then
            value = beta_lsf(a,params)
            return
        end if
        call beta_distances(a,params,left,right)
        logwidth = log_distance(b,a)
        logql = logwidth-log_distance(a,params(3))
        logqr = logwidth-log(right)
        if (logql+log(max(1.0_dp,abs(params(1)-1.0_dp))) <= log(0.001_dp) .and. &
                logqr+log(max(1.0_dp,abs(params(2)-1.0_dp))) <= log(0.001_dp)) then
            ql = exp(logql)
            qr = exp(logqr)
            t = [0.5_dp, 0.5_dp-sqrt(3.0_dp/5.0_dp)/2.0_dp, &
                0.5_dp+sqrt(3.0_dp/5.0_dp)/2.0_dp]
            r = (params(1)-1.0_dp)*log1p_safe(t*ql) + &
                (params(2)-1.0_dp)*log1p_safe(-t*qr)
            correction = (4.0_dp/9.0_dp)*expm1_safe(r(1)) + &
                (5.0_dp/18.0_dp)*(expm1_safe(r(2))+expm1_safe(r(3)))
            value = beta_lpdf(a,params)+logwidth+log1p_safe(correction)
        else
            value = tail_interval(a,b,params,beta_lcdf,beta_lsf)
        end if
    end function beta_interval

end module scifort_censored
