! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

! Lognormal distribution with shape s > 0, location loc, and scale > 0:
! with z = (x - loc) / scale, X is distributed such that log(z) / s is standard normal.
! This matches the parameterization of scipy.stats.lognorm.

module scifort_lognormal
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_math, only : negative_infinity, positive_infinity, quiet_nan, &
        valid_loc_scale
    use scifort_normal, only : normal_cdf, normal_isf, normal_logcdf, &
        normal_logpdf, normal_logsf, normal_pdf, normal_ppf, normal_sf
    implicit none
    private

    public :: lognormal_cdf
    public :: lognormal_isf
    public :: lognormal_logcdf
    public :: lognormal_logpdf
    public :: lognormal_logsf
    public :: lognormal_pdf
    public :: lognormal_ppf
    public :: lognormal_sf

contains

    pure elemental function lognormal_pdf(x, s, loc, scale) result(y)
        ! Evaluates the probability density function of the lognormal distribution.
        real(dp), intent(in) :: x !! point at which the density is evaluated
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: denom
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = log(z) / s
            denom = s * sigma * z
            if (denom == 0.0_dp) then
                y = 0.0_dp
            else
                y = normal_pdf(w) / denom
            end if
        end if
    end function lognormal_pdf

    pure elemental function lognormal_logpdf(x, s, loc, scale) result(y)
        ! Evaluates the log probability density function of the lognormal distribution.
        real(dp), intent(in) :: x !! point at which the log density is evaluated
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: lnp
        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = log(z) / s
            lnp = normal_logpdf(w)
            if (.not. ieee_is_finite(lnp)) then
                y = lnp
            else
                y = lnp - log(s) - log(sigma) - log(z)
            end if
        end if
    end function lognormal_logpdf

    pure elemental function lognormal_cdf(x, s, loc, scale) result(y)
        ! Evaluates the cumulative distribution function of the lognormal distribution.
        real(dp), intent(in) :: x !! point x of the lower-tail probability P(X <= x)
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = log(z) / s
            y = normal_cdf(w)
        end if
    end function lognormal_cdf

    pure elemental function lognormal_sf(x, s, loc, scale) result(y)
        ! Evaluates the survival function P(X > x) of the lognormal distribution.
        real(dp), intent(in) :: x !! point x of the upper-tail probability P(X > x)
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = log(z) / s
            y = normal_sf(w)
        end if
    end function lognormal_sf

    pure elemental function lognormal_logcdf(x, s, loc, scale) result(y)
        ! Evaluates the logarithm of the cumulative distribution function of the lognormal distribution.
        real(dp), intent(in) :: x !! point x of log(P(X <= x))
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = log(z) / s
            y = normal_logcdf(w)
        end if
    end function lognormal_logcdf

    pure elemental function lognormal_logsf(x, s, loc, scale) result(y)
        ! Evaluates the logarithm of the survival function of the lognormal distribution.
        real(dp), intent(in) :: x !! point x of log(P(X > x))
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: y

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = log(z) / s
            y = normal_logsf(w)
        end if
    end function lognormal_logsf

    pure elemental function lognormal_ppf(p, s, loc, scale) result(x)
        ! Evaluates the percent point function (quantile) of the lognormal distribution.
        real(dp), intent(in) :: p !! lower-tail probability in [0, 1]
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: u
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = normal_ppf(p)
            u = s * w
            if (u > 709.0_dp) then
                x = positive_infinity(p)
            else
                z = exp(u)
                x = mu + sigma * z
            end if
        end if
    end function lognormal_ppf

    pure elemental function lognormal_isf(p, s, loc, scale) result(x)
        ! Evaluates the inverse survival function of the lognormal distribution.
        real(dp), intent(in) :: p !! upper-tail probability in [0, 1]
        real(dp), intent(in) :: s !! shape parameter (standard deviation of log), finite and > 0
        real(dp), intent(in), optional :: loc !! location shift (default 0)
        real(dp), intent(in), optional :: scale !! scale parameter, finite and > 0 (default 1)
        real(dp) :: x

        real(dp) :: mu
        real(dp) :: sigma
        real(dp) :: u
        real(dp) :: w
        real(dp) :: z

        call get_loc_scale(loc, scale, mu, sigma)
        if (.not. (valid_loc_scale(mu, sigma) .and. valid_shape(s))) then
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
            w = normal_isf(p)
            u = s * w
            if (u > 709.0_dp) then
                x = positive_infinity(p)
            else
                z = exp(u)
                x = mu + sigma * z
            end if
        end if
    end function lognormal_isf

    pure elemental logical function valid_shape(s) result(valid)
        ! Returns true if the shape parameter is finite and strictly positive.
        real(dp), intent(in) :: s !! shape parameter to check

        valid = ieee_is_finite(s) .and. s > 0.0_dp
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

end module scifort_lognormal
