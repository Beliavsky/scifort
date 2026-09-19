! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Pareto distribution with shape b > 0, location loc, and scale > 0:
! with z = (x - loc) / scale >= 1, f(x) = (b / scale) / z**(b + 1).
! This matches the parameterization of scipy.stats.pareto.

module scifort_pareto
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_constants, only : scifort_log_two
    use scifort_kinds, only : dp
    use scifort_math, only : expm1_safe, log1p_safe, negative_infinity, &
        positive_infinity, quiet_nan, valid_loc_scale
    implicit none
    private

    public :: pareto_cdf
    public :: pareto_isf
    public :: pareto_logcdf
    public :: pareto_logpdf
    public :: pareto_logsf
    public :: pareto_pdf
    public :: pareto_ppf
    public :: pareto_sf

contains

    pure elemental function pareto_pdf(x, b, loc, scale) result(y)
        ! Evaluates the probability density function of the Pareto distribution.
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 1.0_dp) then
            y = 0.0_dp
        else if (z == 1.0_dp) then
            y = b / sigma
        else
            w = log(z)
            y = (b / sigma) * exp(-(b + 1.0_dp) * w)
        end if
    end function pareto_pdf

    pure elemental function pareto_logpdf(x, b, loc, scale) result(y)
        ! Evaluates the log probability density function of the Pareto distribution.
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 1.0_dp) then
            y = negative_infinity(x)
        else
            w = log(z)
            y = log(b) - log(sigma) - (b + 1.0_dp) * w
        end if
    end function pareto_logpdf

    pure elemental function pareto_cdf(x, b, loc, scale) result(y)
        ! Evaluates the cumulative distribution function of the Pareto distribution.
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: u
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 1.0_dp) then
            y = 0.0_dp
        else
            w = log(z)
            u = b * w
            if (u > 709.0_dp) then
                y = 1.0_dp
            else
                y = -expm1_safe(-u)
            end if
        end if
    end function pareto_cdf

    pure elemental function pareto_sf(x, b, loc, scale) result(y)
        ! Evaluates the survival function P(X > x) of the Pareto distribution.
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: u
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 1.0_dp) then
            y = 1.0_dp
        else
            w = log(z)
            u = b * w
            if (u > 709.0_dp) then
                y = 0.0_dp
            else
                y = exp(-u)
            end if
        end if
    end function pareto_sf

    pure elemental function pareto_logcdf(x, b, loc, scale) result(y)
        ! Evaluates the logarithm of the cumulative distribution function of the Pareto distribution.
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: u
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z <= 1.0_dp) then
            y = negative_infinity(x)
        else
            w = log(z)
            u = b * w
            if (u < scifort_log_two) then
                y = log(-expm1_safe(-u))
            else if (u > 709.0_dp) then
                y = 0.0_dp
            else
                y = log1p_safe(-exp(-u))
            end if
        end if
    end function pareto_logcdf

    pure elemental function pareto_logsf(x, b, loc, scale) result(y)
        ! Evaluates the logarithm of the survival function of the Pareto distribution.
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            y = quiet_nan(x)
            return
        end if

        if (ieee_is_nan(x)) then
            y = quiet_nan(x)
            return
        end if

        z = (x - mu) / sigma
        if (z < 1.0_dp) then
            y = 0.0_dp
        else
            w = log(z)
            y = -b * w
        end if
    end function pareto_logsf

    pure elemental function pareto_ppf(p, b, loc, scale) result(x)
        ! Evaluates the percent point function (quantile) of the Pareto distribution.
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: v
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            x = quiet_nan(p)
            return
        end if

        if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = mu + sigma
        else if (p == 1.0_dp) then
            x = positive_infinity(p)
        else
            v = -log1p_safe(-p) / b
            if (v > 709.0_dp) then
                x = positive_infinity(p)
            else
                z = exp(v)
                x = mu + sigma * z
            end if
        end if
    end function pareto_ppf

    pure elemental function pareto_isf(p, b, loc, scale) result(x)
        ! Evaluates the inverse survival function of the Pareto distribution.
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: b !! shape parameter, finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: v
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(b))) then
            x = quiet_nan(p)
            return
        end if

        if (ieee_is_nan(p) .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = quiet_nan(p)
        else if (p == 0.0_dp) then
            x = positive_infinity(p)
        else if (p == 1.0_dp) then
            x = mu + sigma
        else
            v = -log(p) / b
            if (v > 709.0_dp) then
                x = positive_infinity(p)
            else
                z = exp(v)
                x = mu + sigma * z
            end if
        end if
    end function pareto_isf

    pure elemental logical function valid_shape(b) result(valid)
        ! Returns true if the shape parameter is finite and strictly positive.
        real(dp), intent(in) :: b !! shape parameter to check

        valid = ieee_is_finite(b) .and. b > 0.0_dp
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

end module scifort_pareto
