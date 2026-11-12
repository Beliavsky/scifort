! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Independent implementation of the censored likelihood specified in
! https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.fit.html
module scifort_censored
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_log_sqrt_two_pi
    use scifort_empirical_distribution, only : censored_data
    use scifort_math, only : quiet_nan, negative_infinity, positive_infinity, valid_loc_scale, &
        expm1_safe, log1p_safe, log1pexp
    use scifort_normal, only : normal_logpdf, normal_logcdf, normal_logsf
    use scifort_exponential, only : exponential_logpdf, exponential_logcdf, exponential_logsf
    use scifort_log_gamma, only : log_gamma_one_plus, stirling_remainder
    use scifort_gamma, only : gamma_logpdf
    use scifort_incomplete_gamma, only : log_gammainc, log_gammaincc, log_gamma_kernel
    use scifort_weibull, only : weibull_logpdf, weibull_logcdf, weibull_logsf
    implicit none
    private

    public :: censored_data

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
    public :: gamma_censored_loglikelihood, gamma_censored_nnlf

contains

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

end module scifort_censored
