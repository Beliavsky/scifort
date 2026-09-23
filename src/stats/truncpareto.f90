! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Truncated Pareto distribution, matching scipy.stats.truncpareto.
! Standardized support is 1 <= z <= c, with b /= 0 and c > 1.

module scifort_truncpareto
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: truncpareto_cdf, truncpareto_isf, truncpareto_logcdf, truncpareto_logpdf
    public :: truncpareto_logsf, truncpareto_pdf, truncpareto_ppf, truncpareto_sf

contains

    pure elemental function truncpareto_pdf(x, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logy
        logy = truncpareto_logpdf(x, b, c, loc, scale)
        if (ieee_is_nan(logy)) then
            y = logy
        else if (logy == negative_infinity(logy)) then
            y = 0.0_dp
        else
            y = exp(logy)
        end if
    end function truncpareto_pdf

    pure elemental function truncpareto_logpdf(x, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z < 1.0_dp .or. z > c .or. .not. ieee_is_finite(z)) then
                y = negative_infinity(x)
            else
                y = log(abs(b)) - log_abs_normalizer(b, c) - (b + 1.0_dp) * log(z) - log(sigma)
            end if
        end if
    end function truncpareto_logpdf

    pure elemental function truncpareto_cdf(x, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logp, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 1.0_dp) then
                y = 0.0_dp
            else if (z >= c) then
                y = 1.0_dp
            else
                logp = log_abs_one_minus_power(b, z) - log_abs_normalizer(b, c)
                y = exp(logp)
            end if
        end if
    end function truncpareto_cdf

    pure elemental function truncpareto_sf(x, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, logq
        logq = truncpareto_logsf(x, b, c, loc, scale)
        if (ieee_is_nan(logq)) then
            y = logq
        else if (logq == negative_infinity(logq)) then
            y = 0.0_dp
        else if (logq >= 0.0_dp) then
            y = 1.0_dp
        else
            y = exp(logq)
        end if
    end function truncpareto_sf

    pure elemental function truncpareto_logcdf(x, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 1.0_dp) then
                y = negative_infinity(x)
            else if (z >= c) then
                y = 0.0_dp
            else
                y = log_abs_one_minus_power(b, z) - log_abs_normalizer(b, c)
            end if
        end if
    end function truncpareto_logcdf

    pure elemental function truncpareto_logsf(x, b, c, loc, scale) result(y)
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: y, delta, lognum, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            y = quiet_nan(x)
        else if (ieee_is_nan(x)) then
            y = quiet_nan(x)
        else
            z = (x - mu) / sigma
            if (z <= 1.0_dp) then
                y = 0.0_dp
            else if (z >= c) then
                y = negative_infinity(x)
            else
                delta = log(c) - log(z)
                if (b > 0.0_dp) then
                    lognum = -b * log(z) + log_one_minus_exp(-b * delta)
                else
                    lognum = -b * log(c) + log_one_minus_exp(b * delta)
                end if
                y = lognum - log_abs_normalizer(b, c)
            end if
        end if
    end function truncpareto_logsf

    pure elemental function truncpareto_ppf(p, b, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, logq, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu + sigma
        else if (p == 1.0_dp) then
            x = affine_bounded(mu, sigma, c, p)
        else
            logq = logaddexp(log1p_safe(-p), log(p) - b * log(c))
            z = exp(-logq / b)
            x = affine_bounded(mu, sigma, z, p)
        end if
    end function truncpareto_ppf

    pure elemental function truncpareto_isf(p, b, c, loc, scale) result(x)
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite standardized upper endpoint greater than one
        real(dp), intent(in), optional :: loc !! location parameter (default 0)
        real(dp), intent(in), optional :: scale !! positive scale parameter (default 1)
        real(dp) :: x, logq, mu, sigma, z
        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. valid_shapes(b, c) .or. .not. valid_loc_scale(mu, sigma)) then
            x = quiet_nan(p)
        else if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) then
            x = quiet_nan(p)
        else if (p == 1.0_dp) then
            x = mu + sigma
        else if (p == 0.0_dp) then
            x = affine_bounded(mu, sigma, c, p)
        else
            logq = logaddexp(log(p), log1p_safe(-p) - b * log(c))
            z = exp(-logq / b)
            x = affine_bounded(mu, sigma, z, p)
        end if
    end function truncpareto_isf

    pure elemental function log_abs_normalizer(b, c) result(y)
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: c !! finite upper endpoint greater than one
        real(dp) :: y
        y = log_abs_one_minus_exp(-b * log(c))
    end function log_abs_normalizer

    pure elemental function log_abs_one_minus_power(b, z) result(y)
        real(dp), intent(in) :: b !! finite nonzero Pareto exponent
        real(dp), intent(in) :: z !! standardized point greater than one
        real(dp) :: y
        y = log_abs_one_minus_exp(-b * log(z))
    end function log_abs_one_minus_power

    pure elemental function log_abs_one_minus_exp(t) result(y)
        real(dp), intent(in) :: t !! nonzero exponent in log(abs(1-exp(t)))
        real(dp) :: y
        if (t < 0.0_dp) then
            y = log(-expm1_safe(t))
        else
            y = t + log1p_safe(-exp(-t))
        end if
    end function log_abs_one_minus_exp

    pure elemental function log_one_minus_exp(t) result(y)
        real(dp), intent(in) :: t !! negative exponent in log(1-exp(t))
        real(dp) :: y
        if (t < -0.6931471805599453_dp) then
            y = log1p_safe(-exp(t))
        else
            y = log(-expm1_safe(t))
        end if
    end function log_one_minus_exp

    pure elemental function logaddexp(a, b) result(y)
        real(dp), intent(in) :: a !! first logarithm
        real(dp), intent(in) :: b !! second logarithm
        real(dp) :: y, m
        m = max(a, b)
        y = m + log(exp(a - m) + exp(b - m))
    end function logaddexp

    pure elemental logical function valid_shapes(b, c) result(ok)
        real(dp), intent(in) :: b !! candidate Pareto exponent
        real(dp), intent(in) :: c !! candidate standardized upper endpoint
        ok = ieee_is_finite(b) .and. b /= 0.0_dp .and. ieee_is_finite(c) .and. c > 1.0_dp
    end function valid_shapes

    pure elemental function affine_bounded(mu, sigma, z, seed) result(x)
        real(dp), intent(in) :: mu !! finite location parameter
        real(dp), intent(in) :: sigma !! positive finite scale parameter
        real(dp), intent(in) :: z !! finite standardized point
        real(dp), intent(in) :: seed !! NaN construction seed if overflow occurs
        real(dp) :: x
        if (z > (huge(1.0_dp) - abs(mu)) / sigma) then
            x = quiet_nan(seed)
        else
            x = mu + sigma * z
        end if
    end function affine_bounded

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when absent
        real(dp), intent(out) :: sigma !! scale, 1 when absent
        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_truncpareto
