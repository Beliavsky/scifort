! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Weibull distribution (minimum extreme value) with shape c > 0, location loc, and scale > 0:
! with z = (x - loc) / scale >= 0, f(x) = (c / scale) * z**(c - 1) * exp(-z**c).
! This matches the parameterization of scipy.stats.weibull_min.

module scifort_weibull
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: weibull_cdf
    public :: weibull_isf
    public :: weibull_logcdf
    public :: weibull_logpdf
    public :: weibull_logsf
    public :: weibull_pdf
    public :: weibull_ppf
    public :: weibull_sf

contains

    pure elemental function weibull_pdf(x, c, loc, scale) result(y)
        ! Evaluates the probability density function of the Weibull distribution.
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: lnp
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 0.0_dp) then
            y = 0.0_dp
        else if (z == 0.0_dp) then
            if (c < 1.0_dp) then
                y = positive_infinity(z)
            else if (c == 1.0_dp) then
                y = 1.0_dp / sigma
            else
                y = 0.0_dp
            end if
        else
            w = c * log(z)
            if (w > 709.0_dp) then
                y = 0.0_dp
            else
                t = exp(w)
                if (t > 709.0_dp) then
                    y = 0.0_dp
                else
                    lnp = log(c) - log(sigma) + (c - 1.0_dp) * log(z) - t
                    if (lnp < -709.0_dp) then
                        y = 0.0_dp
                    else if (lnp > 709.0_dp) then
                        y = positive_infinity(z)
                    else
                        y = exp(lnp)
                    end if
                end if
            end if
        end if
    end function weibull_pdf

    pure elemental function weibull_logpdf(x, c, loc, scale) result(y)
        ! Evaluates the log probability density function of the Weibull distribution.
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 0.0_dp) then
            y = negative_infinity(x)
        else if (z == 0.0_dp) then
            if (c < 1.0_dp) then
                y = positive_infinity(z)
            else if (c == 1.0_dp) then
                y = -log(sigma)
            else
                y = negative_infinity(z)
            end if
        else
            w = c * log(z)
            if (w > 709.0_dp) then
                y = negative_infinity(x)
            else
                t = exp(w)
                y = log(c) - log(sigma) + (c - 1.0_dp) * log(z) - t
            end if
        end if
    end function weibull_logpdf

    pure elemental function weibull_cdf(x, c, loc, scale) result(y)
        ! Evaluates the cumulative distribution function of the Weibull distribution.
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            w = c * log(z)
            if (w > 709.0_dp) then
                y = 1.0_dp
            else if (w < -709.0_dp) then
                y = exp(w)
            else
                t = exp(w)
                y = -expm1_safe(-t)
            end if
        end if
    end function weibull_cdf

    pure elemental function weibull_sf(x, c, loc, scale) result(y)
        ! Evaluates the survival function P(X > x) of the Weibull distribution.
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 1.0_dp
        else
            w = c * log(z)
            if (w > 709.0_dp) then
                y = 0.0_dp
            else
                t = exp(w)
                if (t > 709.0_dp) then
                    y = 0.0_dp
                else
                    y = exp(-t)
                end if
            end if
        end if
    end function weibull_sf

    pure elemental function weibull_logcdf(x, c, loc, scale) result(y)
        ! Evaluates the logarithm of the cumulative distribution function of the Weibull distribution.
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: t
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = negative_infinity(x)
        else
            w = c * log(z)
            if (w > 709.0_dp) then
                y = 0.0_dp
            else if (w < -709.0_dp) then
                y = w
            else
                t = exp(w)
                if (t < scifort_log_two) then
                    y = log(-expm1_safe(-t))
                else if (t > 709.0_dp) then
                    y = 0.0_dp
                else
                    y = log1p_safe(-exp(-t))
                end if
            end if
        end if
    end function weibull_logcdf

    pure elemental function weibull_logsf(x, c, loc, scale) result(y)
        ! Evaluates the logarithm of the survival function of the Weibull distribution.
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 0.0_dp) then
            y = 0.0_dp
        else
            w = c * log(z)
            if (w > 709.0_dp) then
                y = negative_infinity(x)
            else
                y = -exp(w)
            end if
        end if
    end function weibull_logsf

    pure elemental function weibull_ppf(p, c, loc, scale) result(x)
        ! Evaluates the percent point function (quantile) of the Weibull distribution.
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: v
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
            return
        end if

        if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            v = -log1p_safe(-p)
            z = exp(log(v) / c)
            x = mu + sigma * z
        end if
    end function weibull_ppf

    pure elemental function weibull_isf(p, c, loc, scale) result(x)
        ! Evaluates the inverse survival function of the Weibull distribution.
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: c !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! lower end of the support (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: v
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(c))) then
            x = quiet_nan(p)
            return
        end if

        if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = mu
        else
            v = -log(p)
            z = exp(log(v) / c)
            x = mu + sigma * z
        end if
    end function weibull_isf

    pure elemental logical function valid_shape(c) result(valid)
        ! Returns true if the shape parameter is finite and strictly positive.
        real(dp), intent(in) :: c !! shape parameter to check

        valid = ieee_is_finite(c) .and. c > 0.0_dp
    end function valid_shape

    pure subroutine get_loc_scale(loc, scale, mu, sigma)
        ! Extracts location and scale arguments, substituting default values if absent.
        real(dp), intent(in), optional :: loc !! location argument, if present
        real(dp), intent(in), optional :: scale !! scale argument, if present
        real(dp), intent(out) :: mu !! location, 0 when loc is absent
        real(dp), intent(out) :: sigma !! scale, 1 when scale is absent

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(loc)) mu = loc
        if (present(scale)) sigma = scale
    end subroutine get_loc_scale

end module scifort_weibull
