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
    use scifort_incomplete_gamma, only : log_gammainc, log_gammaincc, log_gamma_kernel
    use scifort_weibull, only : weibull_logpdf, weibull_logcdf, weibull_logsf
    implicit none
    private

    public :: censored_data
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
